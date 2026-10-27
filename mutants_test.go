package holemath

import (
	"bytes"
	"context"
	"flag"
	"fmt"
	"go/ast"
	"go/parser"
	"go/printer"
	"go/token"
	"os"
	"os/exec"
	"path/filepath"
	"regexp"
	"sort"
	"strconv"
	"strings"
	"sync"
	"testing"
	"time"
)

// The mutation run. It changes the package's source one small thing at
// a time, an operator or a number or a guard, and runs the whole suite
// against each change. A change the suite does not notice is a survivor,
// and a survivor is a line of code no test is holding to anything.
//
// It is slow, so it runs only when asked:
//
//	GOWORK=off go test -run TestMutants -mutants -timeout 3h
//
// It runs two mutants at a time unless -mutants.workers says otherwise,
// and at two it takes about six minutes for some seven hundred. With
// -mutants.report it writes what became of each, and with
// -mutants.again and that report it runs only the ones that survived,
// which is the quick way to try a new test against them.
//
// A survivor is answered in one of two ways. If the change makes the
// code wrong, a test is written that fails on it; those tests are in
// mutations_test.go, each named for what it kills. If the change makes
// no difference to anything the package promises, it is listed in
// equivalentMutants with the reason, and the run lets it by.
var (
	runMutants     = flag.Bool("mutants", false, "run the mutation tests, which take minutes")
	mutantsOnly    = flag.String("mutants.only", "", "mutate only where file:function matches this regular expression")
	mutantsWithout = flag.String("mutants.without", "", "leave these test files, comma separated, out of the suite the mutants are run against")
	mutantsReport  = flag.String("mutants.report", "", "write every mutant and what became of it to this file")
	mutantsWorkers = flag.Int("mutants.workers", 2, "run this many mutants at once, each on a copy of the package of its own")
	mutantsAgain   = flag.String("mutants.again", "", "run only the mutants this earlier report says survived")
	mutantsList    = flag.Bool("mutants.list", false, "write the report of every mutant's name and run none of them")
)

// mutantTimeout is how long a mutant's suite may run before the mutant
// is taken as killed. The suite takes about two seconds unmutated; a
// mutant that makes a loop endless or a march a thousand times finer
// has been noticed.
const mutantTimeout = 20 * time.Second

// mutant is one change to one file: where, what it was, what it became,
// and the changed file's source.
type mutant struct {
	file, function string
	line           int
	before, after  string
	source         []byte
}

// name is how a mutant is known in equivalentMutants and in the report:
// by its file, its function and the text changed, not by its line,
// which moves whenever anything above it does.
func (m mutant) name() string {
	return fmt.Sprintf("%s %s: %s => %s", m.file, m.function, m.before, m.after)
}

// shape is a name with the formatting taken out: no whitespace, and no
// comma before a closing bracket. A statement wrapped over more lines, or
// given a trailing comma, prints differently and is the same statement, so
// equivalentMutants is looked up by shape and not by the text as printed.
func shape(name string) string {
	return strings.NewReplacer(",)", ")", ",}", "}").Replace(strings.Join(strings.Fields(name), ""))
}

// TestShapeIgnoresFormatting holds the reason shape exists: the same
// statement, wrapped or given a trailing comma, is the same mutant.
func TestShapeIgnoresFormatting(t *testing.T) {
	one := "geodesic.go f: x = g(a, b, c) + d => x = g(a, b, c) - d"
	wrapped := "geodesic.go f: x = g( a, b, c, ) +\n\td => x = g( a, b, c, ) -\n\td"
	if shape(one) != shape(wrapped) {
		t.Errorf("one statement has two shapes: %q and %q", shape(one), shape(wrapped))
	}
	if shape(one) == shape("geodesic.go f: x = g(a, b, c) + d => x = g(a, b, c) * d") {
		t.Errorf("two different changes have one shape")
	}
}

// The operators swapped. An arithmetic operator has one partner; a
// comparison has two, its edge moved and itself turned about.
var (
	arithmeticSwap = map[token.Token]token.Token{
		token.ADD: token.SUB, token.SUB: token.ADD, token.MUL: token.QUO, token.QUO: token.MUL,
		token.LAND: token.LOR, token.LOR: token.LAND,
		token.ADD_ASSIGN: token.SUB_ASSIGN, token.SUB_ASSIGN: token.ADD_ASSIGN,
		token.MUL_ASSIGN: token.QUO_ASSIGN, token.QUO_ASSIGN: token.MUL_ASSIGN,
		token.INC: token.DEC, token.DEC: token.INC,
	}
	comparisonSwap = map[token.Token][]token.Token{
		token.LSS: {token.LEQ, token.GEQ}, token.LEQ: {token.LSS, token.GTR},
		token.GTR: {token.GEQ, token.LEQ}, token.GEQ: {token.GTR, token.LSS},
		token.EQL: {token.NEQ}, token.NEQ: {token.EQL},
	}
	callSwap = map[string]string{"Min": "Max", "Max": "Min"}
)

// mutantsOf is every mutant of one source file.
func mutantsOf(t *testing.T, path string) []mutant {
	fileSet := token.NewFileSet()
	file, err := parser.ParseFile(fileSet, path, nil, parser.ParseComments)
	if err != nil {
		t.Fatalf("parsing %s: %v", path, err)
	}
	render := func(node any) string {
		var out bytes.Buffer
		if err := printer.Fprint(&out, fileSet, node); err != nil {
			t.Fatalf("printing %s: %v", path, err)
		}
		return out.String()
	}
	oneLine := func(node ast.Node) string {
		// a constant is named without the comment over it, which is not
		// what was changed and would make the name change with the prose
		if spec, ok := node.(*ast.ValueSpec); ok {
			doc := spec.Doc
			spec.Doc = nil
			defer func() { spec.Doc = doc }()
		}
		return strings.Join(strings.Fields(render(node)), " ")
	}
	var mutants []mutant
	// record makes one mutant: change is applied, the file and the
	// statement around the change are printed, and undo puts it back.
	record := func(function string, around, at ast.Node, change, undo func()) {
		before := oneLine(around)
		change()
		mutants = append(mutants, mutant{
			file: filepath.Base(path), function: function, line: fileSet.Position(at.Pos()).Line,
			before: before, after: oneLine(around), source: []byte(render(file)),
		})
		undo()
	}
	for _, declaration := range file.Decls {
		function := "package"
		if decl, ok := declaration.(*ast.FuncDecl); ok {
			function = decl.Name.Name
		}
		// around is the statement or declaration a node sits in, which is
		// what the report shows: the nearest of its parents that is one
		var parents []ast.Node
		around := func() ast.Node {
			for i := len(parents) - 1; i >= 0; i-- {
				switch parents[i].(type) {
				case *ast.AssignStmt, *ast.ReturnStmt, *ast.ExprStmt, *ast.IncDecStmt, *ast.ValueSpec, *ast.DeclStmt:
					return parents[i]
				case *ast.IfStmt:
					return parents[i].(*ast.IfStmt).Cond
				case *ast.ForStmt:
					return parents[i].(*ast.ForStmt).Cond
				case *ast.RangeStmt:
					return parents[i].(*ast.RangeStmt).X
				}
			}
			return parents[len(parents)-1]
		}
		ast.Inspect(declaration, func(node ast.Node) bool {
			if node == nil {
				parents = parents[:len(parents)-1]
				return true
			}
			parents = append(parents, node)
			switch node := node.(type) {
			case *ast.ArrayType:
				// an array's length is not arithmetic: leave it be
				parents = parents[:len(parents)-1]
				return false
			case *ast.BinaryExpr:
				original := node.Op
				swaps := comparisonSwap[original]
				if swap, ok := arithmeticSwap[original]; ok {
					swaps = []token.Token{swap}
				}
				for _, swap := range swaps {
					record(function, around(), node, func() { node.Op = swap }, func() { node.Op = original })
				}
			case *ast.UnaryExpr:
				if node.Op == token.SUB {
					record(function, around(), node, func() { node.Op = token.ADD }, func() { node.Op = token.SUB })
				}
			case *ast.AssignStmt:
				if swap, ok := arithmeticSwap[node.Tok]; ok {
					original := node.Tok
					record(function, node, node, func() { node.Tok = swap }, func() { node.Tok = original })
				}
			case *ast.IncDecStmt:
				original := node.Tok
				record(function, node, node, func() { node.Tok = arithmeticSwap[original] }, func() { node.Tok = original })
			case *ast.BasicLit:
				original := node.Value
				var changed string
				switch node.Kind {
				case token.INT:
					value, err := strconv.ParseInt(original, 0, 64)
					if err != nil {
						return true
					}
					changed = strconv.FormatInt(value+1, 10)
				case token.FLOAT:
					value, err := strconv.ParseFloat(original, 64)
					if err != nil {
						return true
					}
					changed = strconv.FormatFloat(2*value+1, 'g', -1, 64)
				default:
					return true
				}
				record(function, around(), node, func() { node.Value = changed }, func() { node.Value = original })
			case *ast.CallExpr:
				if selector, ok := node.Fun.(*ast.SelectorExpr); ok {
					if swap, ok := callSwap[selector.Sel.Name]; ok {
						original := selector.Sel.Name
						record(function, around(), node, func() { selector.Sel.Name = swap }, func() { selector.Sel.Name = original })
					}
				}
			case *ast.IfStmt:
				// the guard turned about, and the guarded statements taken
				// away. The guard is shown by a stand-in that holds it, since
				// what is changed is which expression the if holds.
				condition := node.Cond
				shown := &ast.ParenExpr{X: condition}
				record(function, shown, node, func() {
					node.Cond = &ast.UnaryExpr{Op: token.NOT, X: &ast.ParenExpr{X: condition}}
					shown.X = node.Cond
				}, func() {
					node.Cond = condition
					shown.X = condition
				})
				body := node.Body.List
				record(function, node, node, func() { node.Body.List = nil }, func() { node.Body.List = body })
			}
			return true
		})
	}
	return mutants
}

// TestMutants is the mutation run: see the top of this file.
func TestMutants(t *testing.T) {
	if !*runMutants {
		t.Skip("the mutation run takes minutes: ask for it with -mutants")
	}
	sources, err := filepath.Glob("*.go")
	if err != nil {
		t.Fatal(err)
	}
	without := map[string]bool{}
	for name := range strings.SplitSeq(*mutantsWithout, ",") {
		without[name] = true
	}
	var only *regexp.Regexp
	if *mutantsOnly != "" {
		only = regexp.MustCompile(*mutantsOnly)
	}
	// again is the mutants an earlier report says survived: after tests
	// are added and no source is changed, they are the only ones whose
	// fate can have changed
	var again map[string]bool
	if *mutantsAgain != "" {
		report, err := os.ReadFile(*mutantsAgain)
		if err != nil {
			t.Fatal(err)
		}
		again = map[string]bool{}
		for line := range strings.SplitSeq(string(report), "\n") {
			if name, survived := strings.CutPrefix(line, "survived"); survived {
				name, _, _ = strings.Cut(strings.TrimSpace(name), "\t")
				// the report puts file:line before the name; the name starts at the next space
				_, name, _ = strings.Cut(name, " ")
				again[name] = true
			}
		}
	}
	var mutants []mutant
	var suite []string
	for _, source := range sources {
		switch {
		case without[source]:
		case strings.HasSuffix(source, "_test.go"):
			suite = append(suite, source)
		default:
			suite = append(suite, source)
			for _, m := range mutantsOf(t, source) {
				if (only == nil || only.MatchString(m.file+":"+m.function)) && (again == nil || again[m.name()]) {
					mutants = append(mutants, m)
				}
			}
		}
	}
	suite = append(suite, "go.mod")
	t.Logf("%d mutants", len(mutants))
	if *mutantsList {
		var names []string
		for _, m := range mutants {
			names = append(names, fmt.Sprintf("listed     %s:%d %s", m.file, m.line, m.name()))
		}
		if err := os.WriteFile(*mutantsReport, []byte(strings.Join(names, "\n")+"\n"), 0o644); err != nil {
			t.Fatal(err)
		}
		return
	}

	equivalentShapes := map[string]string{}
	for name, reason := range equivalentMutants {
		if other, twice := equivalentShapes[shape(name)]; twice {
			t.Fatalf("equivalentMutants lists one mutant twice, as %q and again: %s", other, name)
		}
		equivalentShapes[shape(name)] = reason
	}

	// each worker has a copy of the package of its own, and writes one
	// mutated file over the original at a time
	workers := max(1, *mutantsWorkers)
	outcomes := make([]string, len(mutants))
	failing := make([]string, len(mutants))
	next := make(chan int)
	var group sync.WaitGroup
	for worker := range workers {
		directory := filepath.Join(t.TempDir(), fmt.Sprintf("worker%d", worker), "libholemath")
		if err := os.MkdirAll(directory, 0o755); err != nil {
			t.Fatal(err)
		}
		for _, name := range suite {
			content, err := os.ReadFile(name)
			if err != nil {
				t.Fatal(err)
			}
			if err := os.WriteFile(filepath.Join(directory, name), content, 0o644); err != nil {
				t.Fatal(err)
			}
		}
		group.Go(func() {
			for index := range next {
				outcomes[index], failing[index] = runMutant(directory, mutants[index])
			}
		})
	}
	// tally reads what has become of the mutants so far and writes the
	// report; it is asked as the run goes, so a run that is stopped has
	// still said what it found, and once more at the end
	tally := func(upTo int) (counts map[string]int, survivors []string) {
		counts = map[string]int{}
		var report []string
		for index, m := range mutants[:upTo] {
			outcome := outcomes[index]
			if outcome == "survived" && equivalentShapes[shape(m.name())] != "" {
				outcome = "equivalent"
			}
			counts[outcome]++
			report = append(report, fmt.Sprintf("%-10s %s:%d %s\t%s", outcome, m.file, m.line, m.name(), failing[index]))
			if outcome == "survived" {
				survivors = append(survivors, fmt.Sprintf("%s:%d\t%s", m.file, m.line, m.name()))
			}
		}
		if *mutantsReport != "" {
			if err := os.WriteFile(*mutantsReport, []byte(strings.Join(report, "\n")+"\n"), 0o644); err != nil {
				t.Fatal(err)
			}
		}
		return counts, survivors
	}
	const sayEvery = 50
	for index := range mutants {
		next <- index
		// the mutants before the workers' own are done by now
		if done := index - workers; done > 0 && done%sayEvery == 0 {
			counts, _ := tally(done)
			fmt.Fprintf(os.Stderr, "mutants: %d of %d done: killed %d, survived %d\n", done, len(mutants), counts["killed"], counts["survived"])
		}
	}
	close(next)
	group.Wait()
	counts, survivors := tally(len(mutants))
	t.Logf("killed %d, survived %d, equivalent %d, would not build %d", counts["killed"], counts["survived"], counts["equivalent"], counts["unbuildable"])
	for _, survivor := range survivors {
		t.Errorf("survived: %s", survivor)
	}
}

// TestEquivalentMutantsAreCurrent is the check that the list cannot go
// stale in silence, made cheap: it reads the source and does not run the
// suite, so a name the source no longer has fails in an ordinary test
// run and not at the end of a six minute one.
func TestEquivalentMutantsAreCurrent(t *testing.T) {
	sources, err := filepath.Glob("*.go")
	if err != nil {
		t.Fatal(err)
	}
	seen := map[string]bool{}
	for _, source := range sources {
		if !strings.HasSuffix(source, "_test.go") {
			for _, m := range mutantsOf(t, source) {
				seen[shape(m.name())] = true
			}
		}
	}
	var stale []string
	for name := range equivalentMutants {
		if !seen[shape(name)] {
			stale = append(stale, name)
		}
	}
	sort.Strings(stale)
	for _, name := range stale {
		t.Errorf("equivalentMutants lists a mutant the source no longer has: %s", name)
	}
}

// runMutant writes one mutant into a worker's copy of the package, runs
// the suite there, puts the original back, and says what became of the
// mutant and which test failed on it first.
func runMutant(directory string, m mutant) (outcome, failing string) {
	path := filepath.Join(directory, m.file)
	original, err := os.ReadFile(path)
	if err != nil {
		return "broken: " + err.Error(), ""
	}
	defer os.WriteFile(path, original, 0o644)
	if err := os.WriteFile(path, m.source, 0o644); err != nil {
		return "broken: " + err.Error(), ""
	}
	// the build is not the suite's to time, so it has as long again
	bound, cancel := context.WithTimeout(context.Background(), 2*mutantTimeout)
	defer cancel()
	// the first test to fail has killed the mutant, and the rest need not
	// run; the source is vetted where it is written, not here
	command := exec.CommandContext(bound, "go", "test", "-count=1", "-failfast", "-vet=off", "-timeout", mutantTimeout.String(), ".")
	command.Dir = directory
	command.Env = append(os.Environ(), "GOWORK=off", "GOFLAGS=")
	output, err := command.CombinedOutput()
	if err == nil {
		return "survived", ""
	}
	if bytes.Contains(output, []byte("[build failed]")) || bytes.Contains(output, []byte("[setup failed]")) {
		return "unbuildable", ""
	}
	var tests []string
	for _, match := range regexp.MustCompile(`(?m)^\s*--- FAIL: (\S+)`).FindAllSubmatch(output, -1) {
		tests = append(tests, string(match[1]))
	}
	if bytes.Contains(output, []byte("panic: test timed out")) || bound.Err() != nil {
		tests = append(tests, "(timed out)")
	} else if bytes.Contains(output, []byte("panic:")) {
		tests = append(tests, "(panicked)")
	}
	return "killed", strings.Join(tests, ",")
}

// equivalentMutants is the changes the run lets by, each with why: the
// change makes no difference to anything the package promises, so no
// test could tell it from the original and none should try. A name here
// that the source no longer has fails the run, so the list cannot go
// stale in silence.
var equivalentMutants = map[string]string{
	`frame.go FrameOf: seeds := []Vec4{{0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}} => seeds := []Vec4{{0, 2, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}}`:                                                                                                             scaledSeed,
	`frame.go FrameOf: seeds := []Vec4{{0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}} => seeds := []Vec4{{0, 1, 0, 0}, {0, 0, 2, 0}, {0, 0, 0, 1}}`:                                                                                                             scaledSeed,
	`frame.go FrameOf: seeds := []Vec4{{0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}} => seeds := []Vec4{{0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 2}}`:                                                                                                             scaledSeed,
	`frame.go FrameOf: seeds := []Vec4{{0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}} => seeds := []Vec4{{0, 1, 0, 0}, {0, 1, 1, 0}, {0, 0, 0, 1}}`:                                                                                                             earlierAxisInSeed,
	`frame.go FrameOf: seeds := []Vec4{{0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}} => seeds := []Vec4{{0, 1, 0, 0}, {0, 0, 1, 0}, {0, 1, 0, 1}}`:                                                                                                             earlierAxisInSeed,
	`frame.go FrameOf: seeds := []Vec4{{0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}} => seeds := []Vec4{{0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 1, 1}}`:                                                                                                             earlierAxisInSeed,
	`frame.go FrameOf: seeds := []Vec4{{0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}} => seeds := []Vec4{{0, 1, 0, 0}, {0, 0, 1, 0}, {1, 0, 0, 1}}`:                                                                                                             lastSeed,
	`frame.go Advance: reach < singularityMargin*h.M || radius < innerHorizon => reach <= singularityMargin*h.M || radius < innerHorizon`:                                                                                                                  onTheEdge,
	`frame.go Advance: reach < singularityMargin*h.M || radius < innerHorizon => reach < singularityMargin*h.M || radius <= innerHorizon`:                                                                                                                  onTheEdge,
	`frame.go Advance: stepSize := math.Min(remainingProperTime, h.Pace(reach)/2) => stepSize := math.Min(remainingProperTime, h.Pace(reach)/3)`:                                                                                                           aChoice,
	`frame.go Advance: return r, h.R(r.Space()) < h.Horizon() => return r, h.R(r.Space()) <= h.Horizon()`:                                                                                                                                                  onTheEdge,
	`frame.go package: classifyFarFactor = 4 => classifyFarFactor = 5`:                                                                                                                                                                                     aChoice,
	`frame.go package: classifyFarFloor = 40 => classifyFarFloor = 41`:                                                                                                                                                                                     aChoice,
	`frame.go package: classifySteps = 5000 => classifySteps = 5001`:                                                                                                                                                                                       aChoice,
	`geodesic.go package: hoverRatio = 1000 => hoverRatio = 1001`:                                                                                                                                                                                          aChoice,
	`geodesic.go package: turnedBack = 1e-3 => turnedBack = 1.002`:                                                                                                                                                                                         onlySooner,
	`geodesic.go noReturn: next := y - (y*y*y+linear*y+constant)/(3*y*y+linear) => next := y - (y*y*y+linear*y+constant)/(4*y*y+linear)`:                                                                                                                   wrongSlope,
	`geodesic.go goingInward: return radiusSquared*(point.x*velocity.x+point.y*velocity.y)+(radiusSquared+h.A*h.A)*point.z*velocity.z < 0 => return radiusSquared*(point.x*velocity.x+point.y*velocity.y)+(radiusSquared+h.A*h.A)*point.z*velocity.z <= 0`: onTheEdge,
	`geodesic.go runsAway: return -timeVelocity > hoverRatio*math.Abs(timeMomentum) => return -timeVelocity >= hoverRatio*math.Abs(timeMomentum)`:                                                                                                          onTheEdge,
	`geodesic.go March: radius > farRadius => radius >= farRadius`:                                                                                                                                                                                         onTheEdge,
	`geodesic.go pacedStep: notComingBack := radius < limits.noReturn && h.goingInward(point, velocity, radius) => notComingBack := radius <= limits.noReturn && h.goingInward(point, velocity, radius)`:                                                   onTheEdge,
	`geodesic.go pacedStep: hovering := radius < limits.horizon && runsAway(velocity.t, momentum.t) => hovering := radius <= limits.horizon && runsAway(velocity.t, momentum.t)`:                                                                           onTheEdge,
	`geodesic.go ClosestApproach: i < steps && radius <= leastRadius*(1+turnedBack) => i < steps && radius < leastRadius*(1+turnedBack)`:                                                                                                                   onTheEdge,
	`geodesic.go ClosestApproach: i < steps && radius <= leastRadius*(1+turnedBack) => i < steps && radius <= leastRadius*(2+turnedBack)`:                                                                                                                  onlySooner,
	`geodesic.go ClosestApproach: if point, momentum, givenUp = h.pacedStep(point, momentum, radius, limits); givenUp { break } => if point, momentum, givenUp = h.pacedStep(point, momentum, radius, limits); givenUp { }`:                                givenUpStaysPut,
	`gpu.go SetCamera: math.IsNaN(velocity[0]) => math.IsNaN(velocity[1])`:                                                                                                                                                                                 numbersTogether,
}

// The reasons a change makes no difference.
const (
	// onTheEdge: a comparison of two measured numbers, changed only in
	// which way it falls when they are exactly equal. They are equal
	// only by arrangement, and which way that falls is no promise.
	onTheEdge = "the two differ only when a float lands exactly on the edge"
	// aChoice: a tolerance, changed to another that would do as well.
	// The tests hold what the march and the fall answer, not the number
	// they were tuned with.
	aChoice = "a tolerance that is a choice: this value would serve as well"
	// onlySooner: a tolerance whose only work is to stop a march that
	// already has its answer.
	onlySooner = "stops the march sooner or later, with the same answer"
	// scaledSeed: gram and schmidt scales every axis to unit length.
	scaledSeed = "a seed twice as long gives the same axis once it is scaled to unit length"
	// earlierAxisInSeed: gram and schmidt takes out of each seed what
	// lies along the axes before it, and x, or x and y, is all in those.
	earlierAxisInSeed = "what is added to the seed lies along the axes already made, and is taken out again"
	// lastSeed: the last axis is whatever is left square to the other
	// three; any seed with some z in it finds the same one.
	lastSeed = "the last axis is fixed by the other three, whatever seed finds it"
	// wrongSlope: newton's method with the slope too steep takes shorter
	// steps toward the same root and stops only on it.
	wrongSlope = "a wrong slope makes newton's method slower and not wrong: it stops only at a root"
	// givenUpStaysPut: a ray the march has given up is not stepped, so
	// asking again changes nothing but the time taken.
	givenUpStaysPut = "a ray given up is not moved, so going on asks the same question to the same answer"
	// numbersTogether: orbiting scales the time, x and y parts of a
	// velocity by one root, so all three are numbers or none is.
	numbersTogether = "the velocity's parts are all numbers or all not, so either one says"
)
