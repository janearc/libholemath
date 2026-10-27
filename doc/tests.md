# running the tests

the tests are the spec, and they cover:

- the horizon, the photon sphere, the shadow at root 27 M, the iscos,
  1.75 seconds of arc at the sun, H, E, Lz and the carter constant kept,
  step halving converging, the subrings thinning by e to the minus pi
- a camera reaching the horizon, the arriving ray's round trip, a pixel
  at the shadow coming back black, bending the same in any plane, the
  fall from ten M taking pi over two root five hundred
- matter on a circle staying on it, light of any energy going the same
  way, an eye inside the horizon still seeing out, and a ray given up at
  its photon orbit having the fate it would have had at the horizon.

## running the tests

    GOWORK=off go test ./...

the go.work above this directory does not list this module, and without
the gowork variable set to off, go will not run anything here.

    GOWORK=off go test -run '^$' -bench . ./...

the benchmarks: the radius, one evaluation of hamilton's equations, one
step, a whole march of a ray that escapes and of one that is captured.

    GOWORK=off go test -run '^$' -fuzz FuzzRadii -fuzztime 1m

the fuzzing. fuzz_test.go holds four targets: FuzzRadii, FuzzMetric,
FuzzObserver and FuzzMarch. each holds the package to something true of
every hole and every point, and go's fuzzer looks for the one that is not.

run as part of the ordinary tests they take their seeds only and cost
nothing. `-fuzz` takes one target at a time.

a failing input is written to testdata/fuzz/<target>/. keep it there once
the bug is fixed, and it runs as a seed from then on.

    GOWORK=off go test -run Property -rapid.checks 1000 ./...

the property tests. property_test.go states eighteen things true of
every hole in a range and every point outside its horizon, and rapid draws
the cases and shrinks one that fails to the smallest that still does.
the default is a hundred draws each; ask for more with `-rapid.checks`.

    GOWORK=off go test -run TestMutants -mutants -timeout 3h

the mutation run. it changes the source one operator, number or guard at
a time and runs the suite against each change; a change nothing notices
is a line nothing is testing.

it takes about six minutes at two mutants at a time, so it only runs when
asked. mutants_test.go says how to read what it prints, and
mutations_test.go is the tests it asked for.
