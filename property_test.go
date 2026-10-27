package holemath

import (
	"math"
	"testing"

	"pgregory.net/rapid"
)

// The property tests. Each states something true of every hole and every
// point in its range, and rapid draws many of them, and when one breaks
// it shrinks the case to the smallest one that still breaks it. The
// fuzz targets in fuzz_test.go hunt the same kind of claim with coverage
// guidance; these say what the claim is, and keep the draws they found.
//
// Run with the ordinary tests:
//
//	GOWORK=off go test -run Property ./...
//
// and for more draws, -rapid.checks 1000.

// genHole draws a hole that Valid says is one: a mass from half to five,
// and a spin from nought to just under 0.95 of it.
func genHole(t *rapid.T) Hole {
	mass := rapid.Float64Range(0.5, 5).Draw(t, "mass")
	spin := rapid.Float64Range(0, 0.95).Draw(t, "spin fraction") * mass
	return Hole{M: mass, A: spin}
}

// genOutside draws a point well outside the ergosurface and the horizon,
// at three to forty M from the centre, where a still observer exists and
// no step of a test comes near the hole.
func genOutside(t *rapid.T, h Hole) Vec3 {
	reach := rapid.Float64Range(3*h.M, 40*h.M).Draw(t, "distance")
	cosTheta := rapid.Float64Range(-1, 1).Draw(t, "cos theta")
	phi := rapid.Float64Range(0, 2*math.Pi).Draw(t, "phi")
	sinTheta := math.Sqrt(1 - cosTheta*cosTheta)
	return Vec3{
		reach * sinTheta * math.Cos(phi),
		reach * sinTheta * math.Sin(phi),
		reach * cosTheta,
	}
}

// genDirection draws a unit vector in the observer's space.
func genDirection(t *rapid.T) Vec3 {
	cosTheta := rapid.Float64Range(-1, 1).Draw(t, "direction cos theta")
	phi := rapid.Float64Range(0, 2*math.Pi).Draw(t, "direction phi")
	sinTheta := math.Sqrt(1 - cosTheta*cosTheta)
	return Vec3{sinTheta * math.Cos(phi), sinTheta * math.Sin(phi), cosTheta}
}

// ---- the hole and its radii ----

func TestPropertyValidHolesAreValid(t *testing.T) {
	rapid.Check(t, func(t *rapid.T) {
		hole := genHole(t)
		if !hole.Valid() {
			t.Fatalf("%+v is drawn as a hole and Valid says not", hole)
		}
	})
}

func TestPropertyHorizonsOrdered(t *testing.T) {
	rapid.Check(t, func(t *rapid.T) {
		hole := genHole(t)
		if !(hole.innerHorizon() <= hole.M && hole.M <= hole.Horizon() && hole.Horizon() <= 2*hole.M) {
			t.Fatalf("%+v: horizons %v, %v, %v are not in order", hole, hole.innerHorizon(), hole.M, hole.Horizon())
		}
	})
}

func TestPropertyRadiusIsSymmetric(t *testing.T) {
	// the radius depends on x and y through their squares, and on z
	// through its square, so flipping any sign leaves it where it was
	rapid.Check(t, func(t *rapid.T) {
		hole := genHole(t)
		x := genOutside(t, hole)
		flipped := Vec3{-x[0], -x[1], -x[2]}
		if !within(hole.R(x), hole.R(flipped), 1e-12) {
			t.Fatalf("%+v: R at %v is %v but at %v it is %v", hole, x, hole.R(x), flipped, hole.R(flipped))
		}
	})
}

func TestPropertyRadiusIsTheDistanceWithNoSpin(t *testing.T) {
	rapid.Check(t, func(t *rapid.T) {
		mass := rapid.Float64Range(0.5, 5).Draw(t, "mass")
		hole := Hole{M: mass}
		x := genOutside(t, hole)
		distance := math.Sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2])
		if !within(hole.R(x), distance, 1e-12) {
			t.Fatalf("%+v at %v: R is %v, the distance is %v", hole, x, hole.R(x), distance)
		}
	})
}

// the radii are nested as spin grows: a faster hole has its horizons and
// its orbits in at least as close, and the prograde orbit only moves in
func TestPropertyOrbitsMoveInWithSpin(t *testing.T) {
	rapid.Check(t, func(t *rapid.T) {
		mass := rapid.Float64Range(0.5, 5).Draw(t, "mass")
		low := rapid.Float64Range(0, 0.9).Draw(t, "low spin fraction")
		high := rapid.Float64Range(low, 0.95).Draw(t, "high spin fraction")
		slow, fast := Hole{M: mass, A: low * mass}, Hole{M: mass, A: high * mass}
		slack := 1e-9 * mass
		if fast.PhotonOrbit(true) > slow.PhotonOrbit(true)+slack {
			t.Fatalf("prograde photon orbit moves out with spin: %v at %v, %v at %v",
				slow.PhotonOrbit(true), low, fast.PhotonOrbit(true), high)
		}
		if fast.ISCO(true) > slow.ISCO(true)+slack {
			t.Fatalf("prograde ISCO moves out with spin: %v at %v, %v at %v",
				slow.ISCO(true), low, fast.ISCO(true), high)
		}
		if fast.Horizon() > slow.Horizon()+slack {
			t.Fatalf("horizon moves out with spin: %v at %v, %v at %v",
				slow.Horizon(), low, fast.Horizon(), high)
		}
	})
}

func TestPropertyRetrogradeOrbitsAreOutside(t *testing.T) {
	rapid.Check(t, func(t *rapid.T) {
		hole := genHole(t)
		slack := 1e-9 * hole.M
		if hole.PhotonOrbit(false) < hole.PhotonOrbit(true)-slack {
			t.Fatalf("%+v: retrograde photon orbit %v is inside the prograde %v", hole, hole.PhotonOrbit(false), hole.PhotonOrbit(true))
		}
		if hole.ISCO(false) < hole.ISCO(true)-slack {
			t.Fatalf("%+v: retrograde ISCO %v is inside the prograde %v", hole, hole.ISCO(false), hole.ISCO(true))
		}
		// and each lies on its side of the still hole's six M, to nine M at most
		if hole.ISCO(true) > 6*hole.M+slack || hole.ISCO(false) < 6*hole.M-slack || hole.ISCO(false) > 9*hole.M+slack {
			t.Fatalf("%+v: ISCO %v prograde and %v retrograde are not about 6M, %v", hole, hole.ISCO(true), hole.ISCO(false), 6*hole.M)
		}
	})
}

// ---- the metric ----

func TestPropertyMetricTimesInverseIsIdentity(t *testing.T) {
	rapid.Check(t, func(t *rapid.T) {
		hole := genHole(t)
		point := genOutside(t, hole)
		scalar, _ := hole.kernel(point)
		tolerance := 1e-9 * (1 + scalar*scalar)
		metric, inverse := hole.Metric(point), hole.Inverse(point)
		for i := range 4 {
			for j := range 4 {
				var sum float64
				for k := range 4 {
					sum += metric[i][k] * inverse[k][j]
				}
				want := 0.0
				if i == j {
					want = 1
				}
				if !within(sum, want, tolerance) {
					t.Fatalf("%+v at %v: metric times inverse at %d,%d is %v", hole, point, i, j, sum)
				}
			}
		}
	})
}

func TestPropertyRaiseUndoesLower(t *testing.T) {
	rapid.Check(t, func(t *rapid.T) {
		hole := genHole(t)
		point := genOutside(t, hole)
		var vector Vec4
		for k := range 4 {
			vector[k] = rapid.Float64Range(-10, 10).Draw(t, "component")
		}
		scalar, _ := hole.kernel(point)
		tolerance := 1e-9 * (1 + scalar*scalar)
		back := hole.Raise(point, hole.Lower(point, vector))
		for k := range 4 {
			if !within(back[k], vector[k], tolerance) {
				t.Fatalf("%+v at %v: Raise of Lower of %v is %v", hole, point, vector, back)
			}
		}
	})
}

// the kernel's vector is null for every point, so the light cone is the
// same shape everywhere: a vector is null when its lowered form has
// nothing in it along itself beyond the flat part
func TestPropertyKernelVectorIsNull(t *testing.T) {
	rapid.Check(t, func(t *rapid.T) {
		hole := genHole(t)
		point := genOutside(t, hole)
		_, null := hole.kernel(point)
		flat := -null[0]*null[0] + null[1]*null[1] + null[2]*null[2] + null[3]*null[3]
		if !within(flat, 0, 1e-9) {
			t.Fatalf("%+v at %v: kernel vector %v is not null: %v", hole, point, null, flat)
		}
	})
}

// ---- observers and rays ----

func TestPropertyStaticObserverIsUnitTimelike(t *testing.T) {
	rapid.Check(t, func(t *rapid.T) {
		hole := genHole(t)
		point := genOutside(t, hole)
		u := hole.Static(point)
		if math.IsNaN(u[0]) {
			t.Fatalf("%+v at %v: no still observer, which the draw should not reach", hole, point)
		}
		if !isUnitTimelike(hole, point, u) {
			t.Fatalf("%+v at %v: static observer %v is not unit timelike", hole, point, u)
		}
	})
}

func TestPropertyNullRayHasUnitEnergyAndNoLength(t *testing.T) {
	rapid.Check(t, func(t *rapid.T) {
		hole := genHole(t)
		point := genOutside(t, hole)
		u := hole.Static(point)
		frame := hole.FrameOf(point, u)
		ray := hole.Null(point, frame, genDirection(t))
		if energy := Dot(ray.P, frame[0]); !within(energy, -1, 1e-8) {
			t.Fatalf("%+v at %v: a photon from a still observer has energy %v, want -1", hole, point, energy)
		}
		if h := hole.H(ray); !within(h, 0, 1e-8) {
			t.Fatalf("%+v at %v: a photon is not null: H is %v", hole, point, h)
		}
	})
}

// a null ray stays null as it is marched, and the quantities a geodesic
// keeps stay where they were: energy, angular momentum and the carter
// constant, over a few steps near where a camera would look
func TestPropertyMarchingKeepsWhatAGeodesicKeeps(t *testing.T) {
	rapid.Check(t, func(t *rapid.T) {
		hole := genHole(t)
		point := genOutside(t, hole)
		frame := hole.FrameOf(point, hole.Static(point))
		ray := hole.Null(point, frame, genDirection(t))
		energy, angular, carter := -ray.P[0], hole.Lz(ray), hole.Carter(ray)
		start := hole.R(ray.Space())
		for range 5 {
			ray = hole.Step(ray, 0.02*start)
			if !within(hole.H(ray), 0, 1e-6) {
				t.Fatalf("%+v: a photon stopped being null, H is %v", hole, hole.H(ray))
			}
		}
		if !within(-ray.P[0], energy, 1e-6) {
			t.Fatalf("%+v: the energy went from %v to %v", hole, energy, -ray.P[0])
		}
		if !within(hole.Lz(ray), angular, 1e-5) {
			t.Fatalf("%+v: the angular momentum went from %v to %v", hole, angular, hole.Lz(ray))
		}
		if !within(hole.Carter(ray), carter, 1e-5) {
			t.Fatalf("%+v: the carter constant went from %v to %v", hole, carter, hole.Carter(ray))
		}
	})
}

// hamilton's function is minus a half for a body of unit mass, and
// stays there: the scale of H is not seen by a null ray, which has none
func TestPropertyABodyHasHamiltonianMinusAHalf(t *testing.T) {
	rapid.Check(t, func(t *rapid.T) {
		hole := genHole(t)
		body := hole.Dropped(genOutside(t, hole))
		if h := hole.H(body); !within(h, -0.5, 1e-8) {
			t.Fatalf("%+v: a body let go has H %v, want -0.5", hole, h)
		}
		for range 5 {
			body = hole.Step(body, 0.02*hole.R(body.Space()))
		}
		if h := hole.H(body); !within(h, -0.5, 1e-6) {
			t.Fatalf("%+v: a body stepped has H %v, want -0.5", hole, h)
		}
	})
}

func TestPropertyPaceIsPositiveAndFinite(t *testing.T) {
	rapid.Check(t, func(t *rapid.T) {
		hole := genHole(t)
		point := genOutside(t, hole)
		frame := hole.FrameOf(point, hole.Static(point))
		ray := hole.Null(point, frame, genDirection(t))
		if step := hole.PaceIn(ray); !(step > 0) || math.IsInf(step, 0) {
			t.Fatalf("%+v at %v: the step is %v", hole, point, step)
		}
		if step := hole.Pace(hole.R(point)); !(step > 0) || math.IsInf(step, 0) {
			t.Fatalf("%+v at %v: Pace is %v", hole, point, step)
		}
	})
}

// the angle between two directions is symmetric, nought for a direction
// with itself, and the same under any positive scaling of either
func TestPropertyAngleBetween(t *testing.T) {
	rapid.Check(t, func(t *rapid.T) {
		a := Vec3(genDirection(t))
		b := Vec3(genDirection(t))
		scale := rapid.Float64Range(0.1, 10).Draw(t, "scale")
		ab := angleBetween(a, b)
		if !within(ab, angleBetween(b, a), 1e-12) {
			t.Fatalf("the angle from %v to %v is %v, the other way %v", a, b, ab, angleBetween(b, a))
		}
		if ab < 0 || ab > math.Pi {
			t.Fatalf("the angle from %v to %v is %v, out of [0, pi]", a, b, ab)
		}
		if got := angleBetween(a, a); !within(got, 0, 1e-12) {
			t.Fatalf("the angle from %v to itself is %v", a, got)
		}
		scaled := Vec3{a[0] * scale, a[1] * scale, a[2] * scale}
		if !within(angleBetween(scaled, b), ab, 1e-9) {
			t.Fatalf("scaling %v by %v changed the angle to %v from %v to %v", a, scale, b, ab, angleBetween(scaled, b))
		}
	})
}

// ---- the march ----

// a ray at nought spin, from far away, with impact parameter b, is
// captured below the shadow and escapes above it: the shadow is root
// twenty seven M, and the test leaves a margin of five hundredths of it
func TestPropertyShadowAtNoughtSpin(t *testing.T) {
	rapid.Check(t, func(t *rapid.T) {
		mass := rapid.Float64Range(0.5, 5).Draw(t, "mass")
		hole := Hole{M: mass}
		impact := rapid.Float64Range(0.05, 2).Draw(t, "impact in shadows") * hole.Shadow()
		distance := 100 * mass
		ray := Ray{
			X: Vec4{0, -distance, impact, 0},
			P: Vec4{-1, 1, 0, 0},
		}
		if !within(hole.ImpactParameter(ray), impact, 1e-9) {
			t.Fatalf("%+v: the impact parameter is %v, drawn as %v", hole, hole.ImpactParameter(ray), impact)
		}
		fate := hole.ClassifyRay(ray)
		switch {
		case impact < 0.95*hole.Shadow() && fate != RayPlunges:
			t.Fatalf("%+v: impact %v is inside the shadow and the ray %v", hole, impact, fate)
		case impact > 1.05*hole.Shadow() && fate != RayEscapes:
			t.Fatalf("%+v: impact %v is outside the shadow and the ray %v", hole, impact, fate)
		}
	})
}

// a captured ray at nought spin is given up inside three M, which is
// where every answer lies at that spin, and a ray that escapes ends past
// the far radius it was marched to
func TestPropertyMarchFatesAtNoughtSpin(t *testing.T) {
	rapid.Check(t, func(t *rapid.T) {
		mass := rapid.Float64Range(0.5, 5).Draw(t, "mass")
		hole := Hole{M: mass}
		point := genOutside(t, hole)
		frame := hole.FrameOf(point, hole.Static(point))
		ray := hole.Null(point, frame, genDirection(t))
		far := 4 * hole.R(point)
		end, fate := hole.March(ray, far, 5000)
		switch fate {
		case Captured:
			if !(hole.R(end.Space()) < 3*mass*(1+1e-6)) {
				t.Fatalf("%+v: captured at radius %v, outside the photon sphere %v", hole, hole.R(end.Space()), 3*mass)
			}
		case Escaped:
			if !(hole.R(end.Space()) > far) {
				t.Fatalf("%+v: escaped at radius %v, inside the far radius %v", hole, hole.R(end.Space()), far)
			}
		}
	})
}

// the closest approach is no further out than the ray starts and no
// nearer the centre than the horizon of a still hole allows for light
// that has come from outside
func TestPropertyClosestApproachIsBetweenTheHorizonAndTheStart(t *testing.T) {
	rapid.Check(t, func(t *rapid.T) {
		hole := genHole(t)
		point := genOutside(t, hole)
		frame := hole.FrameOf(point, hole.Static(point))
		ray := hole.Null(point, frame, genDirection(t))
		start := hole.R(point)
		closest := hole.ClosestApproach(ray, 2000)
		if math.IsNaN(closest) {
			t.Fatalf("%+v at %v: the closest approach is not a number", hole, point)
		}
		if closest > start*(1+1e-12) {
			t.Fatalf("%+v at %v: the closest approach %v is past the start %v", hole, point, closest, start)
		}
		if closest < hole.innerHorizon()*(1-1e-9) {
			t.Fatalf("%+v at %v: the closest approach %v is inside the inner horizon %v", hole, point, closest, hole.innerHorizon())
		}
	})
}
