package holemath

import (
	"math"
	"testing"
)

// The mutation tests. Each one is here because the mutation run, in
// mutants_test.go, found a change to the source that every other test
// let by; each says what that change was and fails on it. They add
// nothing to the coverage, which was whole without them: they are the
// second set, the one that tests the tests.
//
// What they found the rest of the suite blind to, in kind: a hole of
// any mass but one; a ray's time; an observer moving three ways at
// once; the last step of a budget; the corners where a radius is
// nought; and a photon whose energy is not about one.

// A hole has one length, its mass. Make the mass four times as great
// and the spin with it, and every length the package answers is four
// times as great, and every angle, ratio, fate and momentum as it was.
// Nearly every other test is of a hole of mass one, where a mass left
// out, or divided by where it should multiply, cannot be seen: the
// horizon's M*M as M/M, the photon sphere's 3*M as 3/M, the ISCO's A/M
// as A*M. Four is a power of two, so the scaled answers have the same
// digits and are held closely.
func TestMutationEverythingScalesWithTheMass(t *testing.T) {
	const scale = 4.0
	near := func(got, want float64) bool {
		return math.Abs(got-want) <= 1e-11*math.Max(1, math.Abs(want)) || (math.IsNaN(got) && math.IsNaN(want))
	}
	scaled := func(x Vec3) Vec3 { return Vec3{scale * x[0], scale * x[1], scale * x[2]} }
	scaledRay := func(r Ray) Ray {
		return Ray{X: Vec4{scale * r.X[0], scale * r.X[1], scale * r.X[2], scale * r.X[3]}, P: r.P}
	}
	sameRay := func(got, want Ray) bool {
		for k := range 4 {
			if !near(got.X[k], want.X[k]) || !near(got.P[k], want.P[k]) {
				return false
			}
		}
		return true
	}
	for _, spin := range []float64{0, 0.6, 0.95} {
		small, large := Hole{M: 1, A: spin}, Hole{M: scale, A: scale * spin}
		for name, radius := range map[string]func(Hole) float64{
			"Horizon":         Hole.Horizon,
			"innerHorizon":    Hole.innerHorizon,
			"PhotonSphere":    Hole.PhotonSphere,
			"Shadow":          Hole.Shadow,
			"ISCO prograde":   func(h Hole) float64 { return h.ISCO(true) },
			"ISCO retrograde": func(h Hole) float64 { return h.ISCO(false) },
		} {
			if got, want := radius(large), scale*radius(small); !near(got, want) || math.IsNaN(got) {
				t.Errorf("spin %v: %s of the larger hole is %v, want %v", spin, name, got, want)
			}
		}
		if got, want := large.Pace(scale*7), scale*small.Pace(7); !near(got, want) {
			t.Errorf("spin %v: the pace is %v, want %v", spin, got, want)
		}
		if got, want := large.Falling(scale*10), scaledRay(small.Falling(10)); !sameRay(got, want) {
			t.Errorf("spin %v: Falling from ten of the larger hole's masses is %v, want %v", spin, got, want)
		}
		// the metric, the observers and what a ray keeps, on and off the
		// equator and near the axis
		for _, position := range []Vec3{{7, 2, 3}, {12, -5, 0}, {0.3, 0.2, 9}, {2.6, 1, 0.4}} {
			if got, want := large.R(scaled(position)), scale*small.R(position); !near(got, want) {
				t.Errorf("spin %v at %v: R is %v, want %v", spin, position, got, want)
			}
			if got, want := large.ringDistance(scaled(position)), scale*small.ringDistance(position); !near(got, want) {
				t.Errorf("spin %v at %v: the ring is %v away, want %v", spin, position, got, want)
			}
			gotMetric, wantMetric := large.Metric(scaled(position)), small.Metric(position)
			gotInverse, wantInverse := large.Inverse(scaled(position)), small.Inverse(position)
			for i := range 4 {
				for j := range 4 {
					if !near(gotMetric[i][j], wantMetric[i][j]) || !near(gotInverse[i][j], wantInverse[i][j]) {
						t.Fatalf("spin %v at %v: the metric or its inverse changed with the mass at %d,%d", spin, position, i, j)
					}
				}
			}
			for name, observer := range map[string]func(Hole, Vec3) Vec4{"Static": Hole.Static, "Circular": Hole.Circular} {
				got, want := observer(large, scaled(position)), observer(small, position)
				for k := range 4 {
					if !near(got[k], want[k]) {
						t.Errorf("spin %v at %v: %s is %v round the larger hole and %v round the smaller", spin, position, name, got, want)
						break
					}
				}
			}
			body, inside := small.Advance(small.Dropped(position), 30)
			largeBody, largeInside := large.Advance(large.Dropped(scaled(position)), scale*30)
			if inside != largeInside || !(sameRay(largeBody, scaledRay(body)) || math.IsNaN(body.P[0])) {
				t.Errorf("spin %v from %v: a body falls to %v, inside %v, round the larger hole and %v, inside %v, round the smaller", spin, position, largeBody, largeInside, body, inside)
			}
			if math.IsNaN(small.Static(position)[0]) {
				continue // inside the ergosurface: no still frame to send rays from
			}
			frame, largeFrame := small.FrameOf(position, small.Static(position)), large.FrameOf(scaled(position), large.Static(scaled(position)))
			for axis := range 4 {
				for k := range 4 {
					if !near(largeFrame[axis][k], frame[axis][k]) || math.IsNaN(frame[axis][k]) {
						t.Fatalf("spin %v at %v: the frame changed with the mass", spin, position)
					}
				}
			}
			direction := Vec3{-0.8, 0.36, 0.48}
			rays := map[string][2]Ray{
				"light leaving":  {small.Null(position, frame, direction), large.Null(scaled(position), largeFrame, direction)},
				"light arriving": {small.Arriving(position, frame, direction), large.Arriving(scaled(position), largeFrame, direction)},
				"a body dropped": {small.Dropped(position), large.Dropped(scaled(position))},
			}
			for name, pair := range rays {
				ray, largeRay := pair[0], pair[1]
				if !sameRay(largeRay, scaledRay(ray)) {
					t.Errorf("spin %v at %v, %s: the larger hole's ray is %v, want %v", spin, position, name, largeRay, scaledRay(ray))
					continue
				}
				checks := map[string][2]float64{
					"H":               {large.H(largeRay), small.H(ray)},
					"E":               {large.E(largeRay), small.E(ray)},
					"Lz":              {large.Lz(largeRay), scale * small.Lz(ray)},
					"Carter":          {large.Carter(largeRay), scale * scale * small.Carter(ray)},
					"ImpactParameter": {large.ImpactParameter(largeRay), scale * small.ImpactParameter(ray)},
					"PaceIn":          {large.PaceIn(largeRay), scale * small.PaceIn(ray)},
					"ClosestApproach": {large.ClosestApproach(largeRay, 3000), scale * small.ClosestApproach(ray, 3000)},
				}
				for what, pair := range checks {
					if !near(pair[0], pair[1]) || math.IsNaN(pair[0]) {
						t.Errorf("spin %v at %v, %s: %s is %v round the larger hole, want %v", spin, position, name, what, pair[0], pair[1])
					}
				}
				if large.Hovering(largeRay) != small.Hovering(ray) || large.ClassifyRay(largeRay) != small.ClassifyRay(ray) {
					t.Errorf("spin %v at %v, %s: hovering or the class changed with the mass", spin, position, name)
				}
				if got, want := large.Step(largeRay, scale*0.05), scaledRay(small.Step(ray, 0.05)); !sameRay(got, want) {
					t.Errorf("spin %v at %v, %s: a step round the larger hole ends at %v, want %v", spin, position, name, got, want)
				}
				end, fate := small.March(ray, 500, 3000)
				largeEnd, largeFate := large.March(largeRay, scale*500, 3000)
				if fate != largeFate || !sameRay(largeEnd, scaledRay(end)) {
					t.Errorf("spin %v at %v, %s: the march ends %v at %v round the larger hole and %v at %v round the smaller", spin, position, name, largeFate, largeEnd, fate, end)
				}
			}
		}
	}
}

// boosted is the four-velocity of something moving at a velocity past
// the observer with a frame, the velocity given along the frame's axes
// as a part of light's.
func boosted(frame Frame, velocity Vec3) Vec4 {
	gamma := 1 / math.Sqrt(1-velocity[0]*velocity[0]-velocity[1]*velocity[1]-velocity[2]*velocity[2])
	var fourVelocity Vec4
	for k := range 4 {
		fourVelocity[k] = gamma * (frame[0][k] + velocity[0]*frame[1][k] + velocity[1]*frame[2][k] + velocity[2]*frame[3][k])
	}
	return fourVelocity
}

// FrameOf makes its axes from the coordinate axes by gram and schmidt,
// x first: the first axis lies in the plane of the observer's velocity
// and the x axis, the second in the space of those and the y axis, and
// the four are orthonormal and keep the coordinates' handedness. An
// observer held still cannot tell a seed of x from a seed of x plus t,
// which is the mutant: its velocity is the t axis. One moving three
// ways at once can.
func TestMutationTheFrameIsBuiltXFirst(t *testing.T) {
	hole := Hole{M: 1, A: 0.6}
	position := Vec3{9, 2, 3}
	observer := boosted(hole.FrameOf(position, hole.Static(position)), Vec3{0.3, 0.2, 0.4})
	frame := hole.FrameOf(position, observer)
	if frame[0] != observer {
		t.Fatalf("the frame's time axis is %v, not the observer's velocity %v", frame[0], observer)
	}
	flat := [4]Vec4{{-1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}}
	for i := range 4 {
		for j := range 4 {
			if got := Dot(hole.Lower(position, frame[i]), frame[j]); !(math.Abs(got-flat[i][j]) <= 1e-12) {
				t.Errorf("axes %d and %d have inner product %v, want %v", i, j, got, flat[i][j])
			}
		}
	}
	// the first axis is some of x and some of the velocity, and nothing
	// else: what it has along y says how much velocity, and its parts
	// along z and t are that much velocity and no more
	ofVelocity := frame[1][2] / observer[2]
	if !(math.Abs(frame[1][3]-ofVelocity*observer[3]) <= 1e-12) || !(math.Abs(frame[1][0]-ofVelocity*observer[0]) <= 1e-12) {
		t.Errorf("the first axis %v is not in the plane of the velocity %v and the x axis", frame[1], observer)
	}
	if !(frame[1][1]-ofVelocity*observer[1] > 0) {
		t.Errorf("the first axis %v points against x", frame[1])
	}
	// the second is some of x, y and the velocity: its part along z says
	// how much velocity, and its part along t is that much and no more
	ofVelocity = frame[2][3] / observer[3]
	if !(math.Abs(frame[2][0]-ofVelocity*observer[0]) <= 1e-12) {
		t.Errorf("the second axis %v is not in the space of the velocity, x and y", frame[2])
	}
	if !(frame[2][2]-ofVelocity*observer[2] > 0) {
		t.Errorf("the second axis %v points against y", frame[2])
	}
	// the third is what is left, on the side of z that keeps the hand
	if !(determinant4(frame) > 0) {
		t.Errorf("the frame has turned its hand about: its determinant is %v", determinant4(frame))
	}
}

// determinant4 is the determinant of four vectors' components, by
// cofactors along the first.
func determinant4(rows [4]Vec4) float64 {
	minor := func(skip int) float64 {
		var m [3][3]float64
		for i := 1; i < 4; i++ {
			column := 0
			for j := range 4 {
				if j != skip {
					m[i-1][column] = rows[i][j]
					column++
				}
			}
		}
		return m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1]) - m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0]) + m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0])
	}
	return rows[0][0]*minor(0) - rows[0][1]*minor(1) + rows[0][2]*minor(2) - rows[0][3]*minor(3)
}

// Every ray the package makes starts at time nought. The mutant started
// Null's and Dropped's at one, and nothing looked.
func TestMutationRaysStartAtTimeNought(t *testing.T) {
	hole := Hole{M: 1, A: 0.6}
	position := Vec3{9, 2, 3}
	frame := hole.FrameOf(position, hole.Static(position))
	rays := map[string]Ray{
		"Null":     hole.Null(position, frame, Vec3{0.6, 0, 0.8}),
		"Arriving": hole.Arriving(position, frame, Vec3{0.6, 0, 0.8}),
		"Dropped":  hole.Dropped(position),
		"Falling":  hole.Falling(10),
	}
	for name, ray := range rays {
		if ray.X[0] != 0 {
			t.Errorf("%s starts its ray at time %v", name, ray.X[0])
		}
	}
	if got := hole.Null(position, frame, Vec3{0.6, 0, 0.8}).Space(); got != position {
		t.Errorf("Null starts its ray at %v, not %v", got, position)
	}
}

// A ray carries its time, and nothing tested it: the mutants ran the
// time backward, divided it by the step where it should multiply, and
// weighed the step's four rates wrongly for the time alone. In these
// coordinates light falling straight in keeps t plus r exactly, and
// light going straight out of a still hole keeps t less r less four M
// times the logarithm of r less two M. An arriving ray's time runs
// back, by the light's travel time.
func TestMutationARayKeepsItsTime(t *testing.T) {
	hole := Hole{M: 1}
	outer, inner := Vec3{20, 0, 0}, Vec3{5, 0, 0}
	outerFrame, innerFrame := hole.FrameOf(outer, hole.Static(outer)), hole.FrameOf(inner, hole.Static(inner))

	fallingIn := hole.Null(outer, outerFrame, Vec3{-1, 0, 0})
	for hole.R(fallingIn.Space()) > 5 {
		fallingIn = hole.Step(fallingIn, hole.PaceIn(fallingIn))
		if sum := fallingIn.X[0] + fallingIn.X[1]; !(math.Abs(sum-20) <= 1e-12) {
			t.Fatalf("light falling in from 20 has t plus r of %v at radius %v", sum, fallingIn.X[1])
		}
	}
	if !(fallingIn.X[0] > 14) {
		t.Errorf("light took a time of %v to fall from 20 to 5", fallingIn.X[0])
	}

	kept := func(r Ray) float64 { return r.X[0] - r.X[1] - 4*math.Log(r.X[1]-2) }
	goingOut := hole.Null(inner, innerFrame, Vec3{1, 0, 0})
	atStart := kept(goingOut)
	for hole.R(goingOut.Space()) < 20 {
		goingOut = hole.Step(goingOut, hole.PaceIn(goingOut))
	}
	if !(math.Abs(kept(goingOut)-atStart) <= 1e-5) {
		t.Errorf("light going out from 5 to %v changed t less r less four log(r less two) by %v", goingOut.X[1], kept(goingOut)-atStart)
	}

	// seen from 5 looking out: the arriving ray runs back up the path of
	// light that fell in, so its time falls as its radius rises
	lookingOut, fate := hole.March(hole.Arriving(inner, innerFrame, Vec3{1, 0, 0}), 20, 10000)
	if fate != Escaped {
		t.Fatalf("an arriving ray looking out had fate %v", fate)
	}
	if sum := lookingOut.X[0] + lookingOut.X[1]; !(math.Abs(sum-5) <= 1e-12) {
		t.Errorf("an arriving ray looking out from 5 has t plus r of %v at radius %v", sum, lookingOut.X[1])
	}
}

// The time velocity that sizes a step and says a ray hovers is the time
// part of Velocity, worked a shorter way. The mutants turned signs
// about in it and crossed its terms, and since only its size against
// the energy was ever used, nothing saw.
func TestMutationTimeVelocityIsVelocitysTimePart(t *testing.T) {
	for _, hole := range []Hole{{M: 1}, {M: 1, A: 0.6}, {M: 2.5, A: 2}} {
		for _, position := range []Vec3{{7, 2, 3}, {-3, 4, -6}, {2.4, 0.3, 0.2}, {0, 0, 9}} {
			for _, momentum := range []Vec4{{-1, 0.3, 0.8, -0.5}, {0.7, -1.1, 0.2, 0.4}, {-2, -0.6, 0.1, 1.3}} {
				ray := Ray{X: Vec4{0, position[0], position[1], position[2]}, P: momentum}
				if got, want := hole.timeVelocity(ray), hole.Velocity(ray)[0]; !(math.Abs(got-want) <= 1e-12*math.Max(1, math.Abs(want))) {
					t.Errorf("hole %v at %v, momentum %v: the time velocity is %v and Velocity's time part is %v", hole, position, momentum, got, want)
				}
			}
		}
	}
}

// Where the kernel has no value it answers flat space, and it knows
// where that is by two tests: a radius of nought, and a denominator of
// nought. The mutants asked for both at once, or asked for one in place
// of nought. Each test stands alone at one of three points: on the disc
// inside the ring a hair off the equator, where the radius rounds to
// nought and the denominator does not; so near the centre of a still
// hole that the denominator underflows and the radius does not; and at
// a radius of exactly one, where neither is nought and the answer is
// not flat.
func TestMutationWhereTheKernelIsFlat(t *testing.T) {
	flat := Vec4{1, 0, 0, 0}
	momentum := Vec4{-1, 0.3, -0.2, 0.7}
	flatVelocity := Vec4{1, 0.3, -0.2, 0.7}

	spinning := Hole{M: 1, A: 0.5}
	onTheDisc := Vec3{0.1, 0, 1e-10}
	if spinning.R(onTheDisc) != 0 {
		t.Fatalf("the radius a hair off the disc is %v, not the nought this test is for", spinning.R(onTheDisc))
	}
	still := Hole{M: 1}
	nearTheCentre := Vec3{1e-90, 0, 0}
	if still.R(nearTheCentre) == 0 {
		t.Fatalf("the radius near the centre is nought, and this test is for one that is not")
	}
	for name, at := range map[string]struct {
		hole  Hole
		point Vec3
	}{
		"on the disc, a hair off the equator": {spinning, onTheDisc},
		"near the centre of a still hole":     {still, nearTheCentre},
	} {
		if scalar, vector := at.hole.kernel(at.point); scalar != 0 || vector != flat {
			t.Errorf("%s the kernel is %v, %v, not flat space's", name, scalar, vector)
		}
		velocity, momentumChange := at.hole.derive(Ray{X: Vec4{0, at.point[0], at.point[1], at.point[2]}, P: momentum})
		if velocity != flatVelocity || momentumChange != (Vec4{}) {
			t.Errorf("%s hamilton's equations give %v and %v, not flat space's", name, velocity, momentumChange)
		}
	}

	atOne := Vec3{1, 0, 0}
	if scalar, vector := still.kernel(atOne); scalar != 2 || vector != (Vec4{1, 1, 0, 0}) {
		t.Errorf("at a radius of one the kernel is %v, %v, want 2 and the x axis", scalar, vector)
	}
	if _, momentumChange := still.derive(Ray{X: Vec4{0, 1, 0, 0}, P: momentum}); momentumChange == (Vec4{}) {
		t.Errorf("at a radius of one the momentum does not change, as if space were flat there")
	}
}

// A budget of steps is that many steps. The mutants began the count at
// one, and let ClosestApproach's run to one more; with budgets in the
// thousands and rays that stop of themselves, nothing saw. March with a
// budget of one takes one step of PaceIn's size and is stuck, and
// ClosestApproach with a budget of n answers the least radius of the
// start and the n points after it.
func TestMutationABudgetIsThatManySteps(t *testing.T) {
	hole := Hole{M: 1, A: 0.6}
	position := Vec3{20, 3, 2}
	frame := hole.FrameOf(position, hole.Static(position))
	ray := hole.Null(position, frame, Vec3{-0.96, 0.2, 0.2})
	walked := []Ray{ray}
	for range 3 {
		last := walked[len(walked)-1]
		walked = append(walked, hole.Step(last, hole.PaceIn(last)))
	}
	// held to a part in a million million and not to the last digit: the
	// compiler fuses a multiply and an add differently in two places
	sameRay := func(got, want Ray) bool {
		for k := range 4 {
			if !(math.Abs(got.X[k]-want.X[k]) <= 1e-12*math.Max(1, math.Abs(want.X[k]))) || !(math.Abs(got.P[k]-want.P[k]) <= 1e-12) {
				return false
			}
		}
		return true
	}
	for budget := range 4 {
		if got, fate := hole.March(ray, 1000, budget); fate != Stuck || !sameRay(got, walked[budget]) {
			t.Errorf("a march of %d steps ended %v at %v; that many steps by hand end at %v", budget, fate, got, walked[budget])
		}
		if got, want := hole.ClosestApproach(ray, budget), hole.R(walked[budget].Space()); !(math.Abs(got-want) <= 1e-12*want) {
			t.Errorf("the closest approach in %d steps is %v; the ray, falling, stands at %v after them", budget, got, want)
		}
	}
	// matter on a circle never turns back and is never captured: only
	// the budget stops it, and a count that runs the wrong way never does
	orbit := Vec3{10, 0, 0}
	still := Hole{M: 1}
	circling := Ray{X: Vec4{0, 10, 0, 0}, P: still.Lower(orbit, still.Circular(orbit))}
	if got := still.ClosestApproach(circling, 200); !(math.Abs(got-10) <= 1e-6) {
		t.Errorf("matter on a circle at 10 came nearest at %v in 200 steps", got)
	}
}

// The angle between two directions was tested only through
// DeflectionAngle, on rays bent a few hundredths of a radian in planes
// through an axis, where most of a cross product is nought: seventeen
// mutants of it lived. Here it is held to the arccosine on directions
// with no nought in them, and to its ends.
func TestMutationAngleBetween(t *testing.T) {
	pairs := [][2]Vec3{
		{{1, 2, 3}, {-2, 0.5, 4}},
		{{0.3, -0.7, 0.2}, {0.9, 0.4, -1.1}},
		{{5, 1, -2}, {-1, 3, 0.5}},
		{{2, -3, 1}, {1, 4, 6}},
	}
	for _, pair := range pairs {
		a, b := pair[0], pair[1]
		dot := a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
		lengths := math.Sqrt((a[0]*a[0] + a[1]*a[1] + a[2]*a[2]) * (b[0]*b[0] + b[1]*b[1] + b[2]*b[2]))
		if got, want := angleBetween(a, b), math.Acos(dot/lengths); !(math.Abs(got-want) <= 1e-12) {
			t.Errorf("the angle between %v and %v is %v, want %v", a, b, got, want)
		}
		if got, want := angleBetween(b, a), angleBetween(a, b); got != want {
			t.Errorf("the angle from %v to %v is %v and back is %v", a, b, want, got)
		}
	}
	ends := map[string]struct {
		a, b Vec3
		want float64
	}{
		"a direction and itself":       {Vec3{1, 2, 3}, Vec3{1, 2, 3}, 0},
		"a direction and its opposite": {Vec3{1, 2, 3}, Vec3{-2, -4, -6}, math.Pi},
		"two at right angles":          {Vec3{1, 2, 3}, Vec3{3, 0, -1}, math.Pi / 2},
		"two a millionth apart":        {Vec3{1, 0, 0}, Vec3{1, 1e-6, 0}, 1e-6},
	}
	for name, end := range ends {
		if got := angleBetween(end.a, end.b); !(math.Abs(got-end.want) <= 1e-15+1e-12*end.want) {
			t.Errorf("the angle between %s is %v, want %v", name, got, end.want)
		}
	}
}

// The innermost stable orbit is where r squared, less six M r, plus or
// less eight A times the root of M r, less three A squared, is nought:
// plus for prograde. The tests of ISCO were at nought spin and at full,
// with tolerances of a hundredth, and five mutants of the formula
// between them lived, three by answering not-a-number to a test that
// let it by.
func TestMutationISCOSolvesItsEquation(t *testing.T) {
	for _, hole := range []Hole{{M: 1}, {M: 1, A: 0.3}, {M: 1, A: 0.5}, {M: 1, A: 0.9}, {M: 1, A: 0.99}, {M: 3, A: 1.5}} {
		for _, prograde := range []bool{true, false} {
			sense := -1.0
			if prograde {
				sense = 1
			}
			radius := hole.ISCO(prograde)
			residue := radius*radius - 6*hole.M*radius + sense*8*hole.A*math.Sqrt(hole.M*radius) - 3*hole.A*hole.A
			if !(math.Abs(residue) <= 1e-9*hole.M*hole.M) {
				t.Errorf("hole %v, prograde %v: the innermost stable orbit at %v leaves its equation %v from nought", hole, prograde, radius, residue)
			}
			if !(radius > hole.Horizon()) {
				t.Errorf("hole %v, prograde %v: the innermost stable orbit at %v is not outside the horizon at %v", hole, prograde, radius, hole.Horizon())
			}
		}
		if hole.A > 0 && !(hole.ISCO(true) < 6*hole.M && hole.ISCO(false) > 6*hole.M) {
			t.Errorf("hole %v: prograde %v and retrograde %v are not either side of six M", hole, hole.ISCO(true), hole.ISCO(false))
		}
	}
}

// The distance to the ring is the hypotenuse of the distance from the
// ring in the plane and the height above it. The mutants crossed x, y
// and z in the first of those, and every fall tested came in along the
// x axis or down the z axis, where the crossed terms are nought.
func TestMutationRingDistance(t *testing.T) {
	hole := Hole{M: 1, A: 0.6}
	for _, point := range []Vec3{{3, 4, 2}, {0.3, -0.4, 0.1}, {-5, 12, -7}} {
		want := math.Hypot(math.Hypot(point[0], point[1])-0.6, point[2])
		if got := hole.ringDistance(point); !(math.Abs(got-want) <= 1e-14*want) {
			t.Errorf("the ring is %v from %v, want %v", got, point, want)
		}
	}
	if got := hole.ringDistance(Vec3{0.6, 0, 0}); got != 0 {
		t.Errorf("a point on the ring is %v from it", got)
	}
}

// thrown is a body thrown straight up along the x axis of a still hole
// from a radius, just fast enough to come to rest at a height.
func thrown(hole Hole, from, height float64) Ray {
	position := Vec3{from, 0, 0}
	// a body at rest at the height has energy root(1 - 2M/height) a
	// unit of mass, and keeps it; the still observer at the start sees
	// that as a gamma
	gamma := math.Sqrt(1-2*hole.M/height) / math.Sqrt(1-2*hole.M/from)
	speed := math.Sqrt(1 - 1/(gamma*gamma))
	velocity := boosted(hole.FrameOf(position, hole.Static(position)), Vec3{speed, 0, 0})
	return Ray{X: Vec4{0, from, 0, 0}, P: hole.Lower(position, velocity)}
}

// ClassifyRay's far radius has a floor of forty M for a ray that starts
// close in. The mutant divided by the mass where it should multiply,
// and at a mass of one that is the same floor. A body thrown from three
// M to twenty M and no further has not escaped: it falls back. Round a
// hole of mass two the mutant's floor is ten M, and calls it escaped.
func TestMutationClassifyingARayThatStartsCloseIn(t *testing.T) {
	for _, hole := range []Hole{{M: 1}, {M: 2}} {
		body := thrown(hole, 3*hole.M, 20*hole.M)
		if !(math.Abs(hole.H(body)+0.5) <= 1e-12) {
			t.Fatalf("mass %v: the thrown body is not a body: H is %v", hole.M, hole.H(body))
		}
		if got := hole.ClassifyRay(body); got != RayPlunges {
			t.Errorf("mass %v: a body thrown from three M to twenty M and back is classed %v, want plunges", hole.M, got)
		}
		// thrown hard enough to reach fifty M, it is past forty and classed
		// as escaping, as the doc says a bound body that goes far is
		if got := hole.ClassifyRay(thrown(hole, 3*hole.M, 50*hole.M)); got != RayEscapes {
			t.Errorf("mass %v: a body thrown from three M to fifty M is classed %v, want escapes", hole.M, got)
		}
	}
}

// A photon's energy is its own business: light of a thousand times the
// frequency goes the same way. Hovering is a ratio of the time velocity
// to the energy and does not change when the momentum is scaled, and
// the step PaceIn gives is as many times shorter as the momentum is
// greater, so it covers the same ground. The mutant divided the
// hovering ratio by the energy where it should multiply; every ray
// tested had an energy of about one, where that is nearly the same.
func TestMutationAPhotonsEnergyDoesNotMoveItsEdge(t *testing.T) {
	hole := Hole{M: 1}
	position := Vec3{6, 0, 0}
	frame := hole.FrameOf(position, hole.Static(position))
	ray := hole.Arriving(position, frame, Vec3{-1, 0, 0})
	var beforeTheEdge Ray
	for !hole.Hovering(ray) {
		beforeTheEdge = ray
		ray = hole.Step(ray, hole.PaceIn(ray))
	}
	for _, scale := range []float64{1e-3, 0.5, 3, 1e3} {
		scaled := func(r Ray) Ray {
			return Ray{X: r.X, P: Vec4{scale * r.P[0], scale * r.P[1], scale * r.P[2], scale * r.P[3]}}
		}
		if !hole.Hovering(scaled(ray)) {
			t.Errorf("a hovering ray stopped hovering when its energy was scaled by %v", scale)
		}
		if hole.Hovering(scaled(beforeTheEdge)) {
			t.Errorf("a ray a step short of hovering began to hover when its energy was scaled by %v", scale)
		}
		if got, want := scale*hole.PaceIn(scaled(beforeTheEdge)), hole.PaceIn(beforeTheEdge); !(math.Abs(got-want) <= 1e-12*want) {
			t.Errorf("a ray near the edge steps %v; with its energy scaled by %v it covers %v in a step", want, scale, got)
		}
	}
}
