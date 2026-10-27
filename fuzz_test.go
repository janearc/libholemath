package holemath

import (
	"math"
	"testing"
)

// The fuzz targets. Each holds the package to something that is true of
// every hole and every point, not of the few the other tests pick, and
// lets go's fuzzer look for the one that is not.
//
// Run as tests they take their seeds only, and cost nothing. To search:
//
//	GOWORK=off go test -run '^$' -fuzz FuzzRadii -fuzztime 1m
//
// One target at a time; -fuzz takes one. A failing input is written to
// testdata/fuzz/<target>/ and MUST be kept there once the bug it found is
// fixed, so it is run as a seed from then on.

// fit folds any number into [low, high). A number already there is itself,
// so a seed reads as what it is; one that is not a number, or is infinite,
// is low. The fuzzer is free to try the whole of a float and the target is
// only given the part of it that means something.
func fit(v, low, high float64) float64 {
	if math.IsNaN(v) || math.IsInf(v, 0) {
		return low
	}
	return low + math.Mod(math.Abs(v), high-low)
}

// within is whether two numbers agree to a relative tolerance, with an
// absolute one of the same size for numbers near nought.
func within(a, b, tolerance float64) bool {
	return math.Abs(a-b) <= tolerance*(1+math.Max(math.Abs(a), math.Abs(b)))
}

// fuzzHole is a hole that Valid agrees is one, from two numbers of any
// size: a mass from a hundredth to a hundred, and a spin from nought to
// just under it.
func fuzzHole(mass, spinFraction float64) Hole {
	m := fit(mass, 0.01, 100)
	return Hole{M: m, A: fit(spinFraction, 0, 1) * m}
}

// FuzzRadii holds the radii of a hole in their order: the horizons round
// the mass, the photon orbit outside the horizon, the innermost stable
// orbit outside the photon orbit, and the retrograde of each outside the
// prograde, between the still hole's value and the extreme.
func FuzzRadii(f *testing.F) {
	f.Add(1.0, 0.0)
	f.Add(1.0, 0.5)
	f.Add(1.0, 0.998)
	f.Add(7.3, 0.999999)
	f.Fuzz(func(t *testing.T, mass, spinFraction float64) {
		hole := fuzzHole(mass, spinFraction)
		if !hole.Valid() {
			t.Fatalf("%+v is made to be a hole and Valid says it is not", hole)
		}
		m := hole.M
		slack := 1e-9 * m
		pro, retro := hole.PhotonOrbit(true), hole.PhotonOrbit(false)
		isco, iscoRetro := hole.ISCO(true), hole.ISCO(false)
		for name, below := range map[string][2]float64{
			"inner horizon under M":        {hole.innerHorizon(), m},
			"M under the outer horizon":    {m, hole.Horizon()},
			"outer horizon under 2M":       {hole.Horizon(), 2 * m},
			"horizon under prograde orbit": {hole.Horizon(), pro},
			"prograde orbit under 3M":      {pro, 3 * m},
			"3M under retrograde orbit":    {3 * m, retro},
			"retrograde orbit under 4M":    {retro, 4 * m},
			"prograde orbit under isco":    {pro, isco},
			"prograde isco under 6M":       {isco, 6 * m},
			"6M under retrograde isco":     {6 * m, iscoRetro},
			"retrograde isco under 9M":     {iscoRetro, 9 * m},
			"retrograde orbit under isco":  {retro, iscoRetro},
		} {
			if !(below[0] <= below[1]+slack) {
				t.Errorf("%+v: %s fails: %v is over %v", hole, name, below[0], below[1])
			}
		}
	})
}

// FuzzMetric holds the metric to being the metric: its inverse is its
// inverse, Lower is the metric on a vector, Raise undoes Lower, and the
// kernel's vector is null. At any point not within a tenth of the mass of
// the ring singularity, which is where the numbers stop agreeing.
func FuzzMetric(f *testing.F) {
	f.Add(1.0, 0.0, 10.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0)
	f.Add(1.0, 0.9, 3.0, 4.0, 5.0, 1.0, 2.0, -1.0, 0.5)
	f.Add(2.0, 0.5, 0.5, 0.2, 0.1, 0.0, 1.0, 1.0, 1.0)
	f.Fuzz(func(t *testing.T, mass, spinFraction, x, y, z, v0, v1, v2, v3 float64) {
		hole := fuzzHole(mass, spinFraction)
		reach := 50 * hole.M
		point := Vec3{
			fit(x, -reach, reach), fit(y, -reach, reach), fit(z, -reach, reach),
		}
		if hole.R(point) < 0.1*hole.M {
			t.Skip("within a tenth of the mass of the singularity")
		}
		vector := Vec4{
			fit(v0, -10, 10), fit(v1, -10, 10), fit(v2, -10, 10), fit(v3, -10, 10),
		}
		scalar, null := hole.kernel(point)
		// the numbers lose digits as the square of the scalar
		tolerance := 1e-9 * (1 + scalar*scalar)

		if flat := -null[0]*null[0] + null[1]*null[1] + null[2]*null[2] + null[3]*null[3]; !within(flat, 0, tolerance) {
			t.Errorf("%+v at %v: the kernel vector is not null: %v", hole, point, flat)
		}
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
					t.Errorf("%+v at %v: metric times inverse at %d,%d is %v", hole, point, i, j, sum)
				}
			}
		}
		lowered := hole.Lower(point, vector)
		byMetric := mul(metric, vector)
		for k := range 4 {
			if !within(lowered[k], byMetric[k], tolerance) {
				t.Errorf("%+v at %v: Lower is %v and the metric on it is %v", hole, point, lowered, byMetric)
				break
			}
		}
		back := hole.Raise(point, lowered)
		for k := range 4 {
			if !within(back[k], vector[k], tolerance) {
				t.Errorf("%+v at %v: Raise of Lower of %v is %v", hole, point, vector, back)
				break
			}
		}
	})
}

// FuzzObserver holds an observer's frame to being a frame, and the light
// it sees to coming from where Seen says. For a still observer or one
// going round with the matter of a disc, wherever there is one: the four
// axes are square to one another and to unit length, the ray Arriving
// builds from a direction is null with energy one, and Seen of it is that
// direction.
func FuzzObserver(f *testing.F) {
	f.Add(1.0, 0.0, 10.0, 0.0, 0.0, false, 1.0, 0.5)
	f.Add(1.0, 0.9, 6.0, 3.0, 2.0, true, 2.0, 4.0)
	f.Add(1.0, 0.5, 5.0, -4.0, 1.0, false, 0.3, 5.5)
	f.Fuzz(func(t *testing.T, mass, spinFraction, x, y, z float64, circling bool, polar, azimuth float64) {
		hole := fuzzHole(mass, spinFraction)
		reach := 50 * hole.M
		point := Vec3{
			fit(x, -reach, reach), fit(y, -reach, reach), fit(z, -reach, reach),
		}
		velocity := hole.Static(point)
		if circling {
			velocity = hole.Circular(point)
		}
		if math.IsNaN(velocity[0]) {
			t.Skip("nothing goes so there")
		}
		// near the ergosurface the velocity runs off and so do the digits
		scalar, _ := hole.kernel(point)
		speed := math.Max(math.Abs(velocity[0]), 1)
		if speed > 1e3 || scalar > 1e3 {
			t.Skip("too near the edge of where there is a frame")
		}
		tolerance := 1e-8 * speed * speed * (1 + scalar*scalar)

		frame := hole.FrameOf(point, velocity)
		for a := range 4 {
			for b := range 4 {
				got := Dot(hole.Lower(point, frame[a]), frame[b])
				want := 0.0
				switch {
				case a != b:
				case a == 0:
					want = -1
				default:
					want = 1
				}
				if !within(got, want, tolerance) {
					t.Fatalf("%+v at %v: frame axes %d and %d have product %v, want %v", hole, point, a, b, got, want)
				}
			}
		}

		theta, phi := fit(polar, 0, math.Pi), fit(azimuth, 0, 2*math.Pi)
		direction := Vec3{
			math.Sin(theta) * math.Cos(phi), math.Sin(theta) * math.Sin(phi), math.Cos(theta),
		}
		ray := hole.Arriving(point, frame, direction)
		if energy := Dot(ray.P, frame[0]); !within(energy, 1, tolerance) {
			t.Errorf("%+v at %v: the arriving ray has energy %v to its observer", hole, point, energy)
		}
		if h := hole.H(ray); !within(h, 0, tolerance) {
			t.Errorf("%+v at %v: the arriving ray is not null: H is %v", hole, point, h)
		}
		seen := hole.Seen(frame, ray)
		for k := range 3 {
			if !within(seen[k], direction[k], tolerance) {
				t.Errorf("%+v at %v: Seen is %v of a ray from %v", hole, point, seen, direction)
				break
			}
		}
	})
}

// FuzzMarch holds a march to ending, and to ending as it says. A ray of
// any momentum from anywhere, not light and not a body, any number of
// steps: March comes back, with a fate it was given; it is the same march
// twice; a ray it says escaped is past the far radius where it stopped;
// no steps leave the ray where it was; and ClosestApproach is never
// further than where the ray began.
func FuzzMarch(f *testing.F) {
	f.Add(1.0, 0.0, 30.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 200.0, 2.0)
	f.Add(1.0, 0.9, 20.0, 2.0, 1.0, -1.0, 0.1, 0.3, 0.1, 200.0, 2.0)
	f.Add(1.0, 0.5, 6.0, 0.0, 0.5, 0.0, 0.1, -1.0, 0.0, 300.0, 1.0)
	f.Add(1.0, 0.0, 2.5, 0.0, 0.0, -1.0, -1.0, 0.0, 0.0, 50.0, 0.5)
	f.Fuzz(func(t *testing.T, mass, spinFraction, x, y, z, p0, p1, p2, p3, budget, far float64) {
		hole := fuzzHole(mass, spinFraction)
		reach := 100 * hole.M
		ray := Ray{
			X: Vec4{0, fit(x, -reach, reach), fit(y, -reach, reach), fit(z, -reach, reach)},
			P: Vec4{fit(p0, -5, 5), fit(p1, -5, 5), fit(p2, -5, 5), fit(p3, -5, 5)},
		}
		steps := int(fit(budget, 0, 400))
		// the far radius is further out than the start, so the ray has
		// somewhere to go
		start := hole.R(ray.Space())
		farRadius := start * (1 + fit(far, 0.05, 3))

		end, fate := hole.March(ray, farRadius, steps)
		if fate != Escaped && fate != Captured && fate != Stuck {
			t.Fatalf("%+v %+v: March gave a fate that is none of the three: %d", hole, ray, fate)
		}
		again, fateAgain := hole.March(ray, farRadius, steps)
		if fate != fateAgain || !sameRay(end, again) {
			t.Errorf("%+v %+v: March twice gave %v %d and then %v %d", hole, ray, end, fate, again, fateAgain)
		}
		if fate == Escaped {
			if radius := hole.R(end.Space()); !(radius > farRadius) {
				t.Errorf("%+v %+v: escaped, and stopped at radius %v inside the far radius %v", hole, ray, radius, farRadius)
			}
		}
		if steps == 0 && (fate != Stuck || !sameRay(end, ray)) {
			t.Errorf("%+v %+v: no steps moved the ray or gave it a fate: %v %d", hole, ray, end, fate)
		}

		if closest := hole.ClosestApproach(ray, steps); !math.IsNaN(closest) && closest > start {
			t.Errorf("%+v %+v: the closest approach %v is further than the start %v", hole, ray, closest, start)
		}
	})
}

// sameRay is whether two rays are the same to the last bit, a number
// that is not a number being the same as another.
func sameRay(a, b Ray) bool {
	for k := range 4 {
		if math.Float64bits(a.X[k]) != math.Float64bits(b.X[k]) && !(math.IsNaN(a.X[k]) && math.IsNaN(b.X[k])) {
			return false
		}
		if math.Float64bits(a.P[k]) != math.Float64bits(b.P[k]) && !(math.IsNaN(a.P[k]) && math.IsNaN(b.P[k])) {
			return false
		}
	}
	return true
}
