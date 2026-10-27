package holemath

import (
	"math"
	"testing"
	"time"
)

// The exact numbers relativity hands over: the horizon, the photon
// sphere, the shadow, the innermost stable orbits.
func TestTheExactRadii(t *testing.T) {
	schwarzschildHole := Hole{M: 1}
	if schwarzschildHole.Horizon() != 2 || schwarzschildHole.PhotonSphere() != 3 || !(math.Abs(schwarzschildHole.Shadow()-math.Sqrt(27)) <= 1e-12) {
		t.Errorf("schwarzschild: horizon %v, photon sphere %v, shadow %v", schwarzschildHole.Horizon(), schwarzschildHole.PhotonSphere(), schwarzschildHole.Shadow())
	}
	if isco := schwarzschildHole.ISCO(true); !(math.Abs(isco-6) <= 1e-9) {
		t.Errorf("the isco at nought spin is %v", isco)
	}
	kerrHole := Hole{M: 1, A: 0.999999}
	if horizon := kerrHole.Horizon(); !(math.Abs(horizon-1) <= 0.002) {
		t.Errorf("the horizon at full spin is %v", horizon)
	}
	if prograde, retrograde := kerrHole.ISCO(true), kerrHole.ISCO(false); !(math.Abs(prograde-1) <= 0.02) || !(math.Abs(retrograde-9) <= 0.01) {
		t.Errorf("the iscos at full spin are %v and %v", prograde, retrograde)
	}
	halfSpinHole := Hole{M: 1, A: 0.5}
	if horizon := halfSpinHole.Horizon(); !(math.Abs(horizon-(1+math.Sqrt(0.75))) <= 1e-12) {
		t.Errorf("the horizon at half spin is %v", horizon)
	}
}

// A hole the package is right for has a mass above nought and a spin
// from nought to just under it. Outside that the radii are not a
// hole's: at a spin of 1.2 M the innermost orbit comes out inside the
// horizon that is not there.
func TestValid(t *testing.T) {
	valid := []Hole{{M: 1}, {M: 1, A: 0.6}, {M: 2.5, A: 2}, {M: 1, A: 0.999999}, {M: 1e-9}, {M: math.MaxFloat64}}
	for _, hole := range valid {
		if !hole.Valid() {
			t.Errorf("hole %v should be valid", hole)
		}
	}
	invalid := map[string]Hole{
		"no mass":               {},
		"mass below nought":     {M: -1},
		"spin at the mass":      {M: 1, A: 1},
		"spin above the mass":   {M: 1, A: 1.2},
		"spin below nought":     {M: 1, A: -0.5},
		"a mass not a number":   {M: math.NaN()},
		"a spin not a number":   {M: 1, A: math.NaN()},
		"spin and no mass":      {A: 0.5},
		"an infinite mass":      {M: math.Inf(1)},
		"a spin above, no sign": {M: 2, A: 2.0000001},
	}
	for name, hole := range invalid {
		if hole.Valid() {
			t.Errorf("%s: hole %v should not be valid", name, hole)
		}
	}
	tooFast := Hole{M: 1, A: 1.2}
	if !(tooFast.ISCO(true) < tooFast.Horizon()) {
		t.Errorf("a hole spun past its mass has its innermost orbit at %v outside its horizon at %v; the doc of Valid says otherwise", tooFast.ISCO(true), tooFast.Horizon())
	}
}

// The metric and its inverse are inverses, and flat far away.
func TestTheInverseIsTheInverse(t *testing.T) {
	hole := Hole{M: 1, A: 0.7}
	for _, position := range []Vec3{{5, 2, 1}, {3, -4, 0.5}, {-2, 1, 6}} {
		metric, inverseMetric := hole.Metric(position), hole.Inverse(position)
		for i := range 4 {
			for j := range 4 {
				sum := 0.0
				for k := range 4 {
					sum += metric[i][k] * inverseMetric[k][j]
				}
				want := 0.0
				if i == j {
					want = 1
				}
				if !(math.Abs(sum-want) <= 1e-9) {
					t.Fatalf("at %v the metric times its inverse has %v at %d,%d", position, sum, i, j)
				}
			}
		}
	}
	farMetric := hole.Metric(Vec3{1e6, 0, 0})
	if !(math.Abs(farMetric[0][0]+1) <= 1e-5) || !(math.Abs(farMetric[1][1]-1) <= 1e-5) || !(math.Abs(farMetric[0][1]) <= 1e-5) {
		t.Errorf("far away the metric is %v", farMetric)
	}
}

// A null ray keeps H at nought and its energy and angular momentum,
// and halving the step changes where it ends by little.
func TestAGeodesicKeepsWhatItShould(t *testing.T) {
	hole := Hole{M: 1, A: 0.6}
	position := Vec3{-30, 6, 2}
	frame := hole.FrameOf(position, hole.Static(position))
	ray := hole.Null(position, frame, Vec3{1, 0, 0})
	if hamiltonianValue := hole.H(ray); !(math.Abs(hamiltonianValue) <= 1e-10) {
		t.Fatalf("at the start H is %v", hamiltonianValue)
	}
	energy, angularMomentum := hole.E(ray), hole.Lz(ray)
	endRay, fate := hole.March(ray, 60, 20000)
	if fate != Escaped {
		t.Fatalf("the ray's fate is %v", fate)
	}
	if hamiltonianValue := hole.H(endRay); !(math.Abs(hamiltonianValue) <= 1e-6) {
		t.Errorf("at the end H is %v", hamiltonianValue)
	}
	if !(math.Abs(hole.E(endRay)-energy) <= 1e-9) || !(math.Abs(hole.Lz(endRay)-angularMomentum) <= 1e-6*math.Abs(angularMomentum)) {
		t.Errorf("energy went %v to %v, angular momentum %v to %v", energy, hole.E(endRay), angularMomentum, hole.Lz(endRay))
	}
	// the same march at half the pace lands within a small part of the
	// bend of the first
	fineStepRay := ray
	for i := 0; i < 400000; i++ {
		radius := hole.R(fineStepRay.Space())
		if radius > 60 {
			break
		}
		fineStepRay = hole.Step(fineStepRay, hole.Pace(radius)/2)
	}
	velocityEnd, velocityFine := hole.Velocity(endRay), hole.Velocity(fineStepRay)
	separationRatio := math.Hypot(math.Hypot(velocityEnd[1]-velocityFine[1], velocityEnd[2]-velocityFine[2]), velocityEnd[3]-velocityFine[3]) / math.Hypot(math.Hypot(velocityEnd[1], velocityEnd[2]), velocityEnd[3])
	if !(separationRatio <= 1e-3) {
		t.Errorf("halving the step turns the ray by %v", separationRatio)
	}
}

// deflect is the angle a ray from far away is bent by passing a hole
// with no spin at an impact parameter, and its fate.
func deflect(hole Hole, impactParameter float64) (float64, Fate) {
	// from far enough that the hole is a point, and a little further
	// still for the far radius, so the start is inside it
	distance := math.Max(2000, 20*impactParameter)
	position := Vec3{-distance, impactParameter, 0}
	frame := hole.FrameOf(position, hole.Static(position))
	ray := hole.Null(position, frame, Vec3{1, 0, 0})
	endRay, fate := hole.March(ray, 1.2*distance, 200000)
	velocity := hole.Velocity(endRay)
	return math.Atan2(velocity[2], velocity[1]), fate
}

// A ray grazing the sun bends by 1.75 seconds of arc, and a ray at the
// shadow's edge is just captured, just escapes a hair further out.
func TestTheSunBendsLightAndTheShadowCaptures(t *testing.T) {
	sunHole := Hole{M: 1}
	// the sun's radius in units of its mass, 1.477 km
	angle, fate := deflect(sunHole, 696000/1.477)
	if fate != Escaped {
		t.Fatalf("a ray grazing the sun was %v", fate)
	}
	arcsec := math.Abs(angle) / (math.Pi / 180 / 3600)
	if !(math.Abs(arcsec-1.75) <= 0.02) {
		t.Errorf("a ray grazing the sun bends %.3f seconds of arc", arcsec)
	}
	if _, fate := deflect(sunHole, 5.15); fate != Captured {
		t.Errorf("inside the shadow's edge a ray was %v", fate)
	}
	if _, fate := deflect(sunHole, 5.25); fate != Escaped {
		t.Errorf("outside the shadow's edge a ray was %v", fate)
	}
}

// The bending of light is the same in any plane through a still hole,
// and at a hundred M it is the weak field's: four M over b, and fifteen
// pi over four times M over b squared after it.
func TestDeflectionInAnyPlane(t *testing.T) {
	hole := Hole{M: 1}
	impactParameter := 100.0
	want := 4/impactParameter + 15*math.Pi/4/(impactParameter*impactParameter)
	starts := map[string]Vec3{
		"the xy plane":   {-1e4, impactParameter, 0},
		"the xz plane":   {-1e4, 0, impactParameter},
		"a tilted plane": {-1e4, 0.6 * impactParameter, 0.8 * impactParameter},
	}
	var angles []float64
	for plane, start := range starts {
		frame := hole.FrameOf(start, hole.Static(start))
		angle := hole.DeflectionAngle(hole.Null(start, frame, Vec3{1, 0, 0}), 1e5, 2000)
		if !(math.Abs(angle-want) < 0.003*want) {
			t.Errorf("in %s the ray bends by %v; the weak field says %v", plane, angle, want)
		}
		angles = append(angles, angle)
	}
	for _, angle := range angles {
		if !(math.Abs(angle-angles[0]) < 1e-6) {
			t.Errorf("the bending differs between planes: %v", angles)
		}
	}
}

// A ray going straight out is not bent, and a ray that falls in has no
// deflection at all: the answer is not a number, and not the nought
// that would say it went straight.
func TestDeflectionOfRaysThatAreNotBent(t *testing.T) {
	hole := Hole{M: 1}
	position := Vec3{100, 0, 0}
	frame := hole.FrameOf(position, hole.Static(position))
	if angle := hole.DeflectionAngle(hole.Null(position, frame, Vec3{1, 0, 0}), 1000, 2000); angle != 0 {
		t.Errorf("a ray going straight out bends by %v", angle)
	}
	if angle := hole.DeflectionAngle(hole.Null(position, frame, Vec3{-1, 0, 0}), 1000, 2000); !math.IsNaN(angle) {
		t.Errorf("a ray that falls in has a deflection of %v", angle)
	}
}

// A photon let go sideways at three M circles the hole for a while.
func TestThePhotonSphereHoldsALightFor(t *testing.T) {
	hole := Hole{M: 1}
	position := Vec3{3, 0, 0}
	frame := hole.FrameOf(position, hole.Static(position))
	ray := hole.Null(position, frame, Vec3{0, 1, 0})
	turned := 0.0
	lastAngle := math.Atan2(position[1], position[0])
	for i := 0; i < 4000; i++ {
		ray = hole.Step(ray, 0.005)
		spacePoint := ray.Space()
		radius := hole.R(spacePoint)
		if !(math.Abs(radius-3) <= 0.06) {
			break
		}
		currentAngle := math.Atan2(spacePoint[1], spacePoint[0])
		angleChange := currentAngle - lastAngle
		for angleChange > math.Pi {
			angleChange -= 2 * math.Pi
		}
		for angleChange < -math.Pi {
			angleChange += 2 * math.Pi
		}
		turned += angleChange
		lastAngle = currentAngle
	}
	if !(turned >= 2*math.Pi) {
		t.Errorf("the photon turned %.2f before leaving three M", turned)
	}
}

// Light from matter going round the hole is redder from the receding
// side and bluer from the approaching side as a far observer sees it,
// and gravity reddens light from near the hole even at rest. The far
// observer is still, so only the momentum's time part reaches it, and
// a ray keeps that: the momentum at the emitter serves at both ends.
func TestRedshiftRunsTheRightWay(t *testing.T) {
	hole := Hole{M: 1}
	farObserver := Vec4{1, 0, 0, 0}
	photon := Vec4{-1, 0, 0, 0} // a photon of unit energy at infinity, heading some way
	position := Vec3{8, 0, 0}
	circularVelocity := hole.Circular(position)
	// seen from plus y, the matter at plus x moves toward the observer;
	// a photon on its way there has its momentum along plus y
	towardPhoton := Vec4{-1, 0, 1, 0}
	awayPhoton := Vec4{-1, 0, -1, 0}
	if blueShift, redShift := Shift(towardPhoton, farObserver, towardPhoton, circularVelocity), Shift(awayPhoton, farObserver, awayPhoton, circularVelocity); !(blueShift > redShift) {
		t.Errorf("approaching %v, receding %v", blueShift, redShift)
	}
	staticVelocity := hole.Static(position)
	if redshiftValue := Shift(photon, farObserver, photon, staticVelocity); !(redshiftValue < 1 && redshiftValue > 0.8) {
		t.Errorf("light from rest at eight M arrives at %v of its energy", redshiftValue)
	}
}

// The shift between a still eye and matter going round the hole takes
// the photon's momentum at each end, and is the number the two things
// the ray keeps say it is: the eye's clock rate times E, over the
// matter's times E less its turning rate times Lz. The eye's momentum
// used at both ends, which one-momentum forms invite, is wrong by a
// part in a hundred here.
func TestShiftFromMatterOnACircle(t *testing.T) {
	hole := Hole{M: 1, A: 0.6}
	eye := Vec3{3, -40, 6}
	frame := hole.FrameOf(eye, hole.Static(eye))
	ray := hole.Arriving(eye, frame, Vec3{0.12, 1, -0.15})
	momentumAtEye := ray.P
	steps := 0
	for ray.X[3] > 0 && steps < 5000 {
		ray = hole.Step(ray, hole.PaceIn(ray))
		steps++
	}
	disc := ray.Space()
	if radius := hole.R(disc); !(ray.X[3] <= 0 && radius >= hole.ISCO(true) && radius <= 20) {
		t.Fatalf("the ray did not land on the disc: at %v, radius %v, after %d steps", disc, radius, steps)
	}
	matter := hole.Circular(disc)
	turningRate := matter[2] / matter[0] / disc[0]
	want := hole.Static(eye)[0] * hole.E(ray) / (matter[0] * (hole.E(ray) - turningRate*hole.Lz(ray)))
	got := Shift(momentumAtEye, hole.Static(eye), ray.P, matter)
	if !(math.Abs(got-want) <= 1e-6*want) {
		t.Errorf("the shift from matter on a circle is %v; E and Lz say %v", got, want)
	}
	if oneMomentum := Shift(momentumAtEye, hole.Static(eye), momentumAtEye, matter); math.Abs(oneMomentum-want) < 1e-3*want {
		t.Errorf("the eye's momentum at both ends gives %v, as good as the right %v, so this ray tests nothing", oneMomentum, want)
	}
}

// A falling camera from rest at ten M reaches the horizon in finite
// proper time, its frame staying orthonormal on the way.
func TestACameraFallsIn(t *testing.T) {
	hole := Hole{M: 1}
	camera := hole.Falling(10)
	if hamiltonianValue := hole.H(camera); !(math.Abs(hamiltonianValue+0.5) <= 1e-9) {
		t.Fatalf("a body of unit mass has H %v", hamiltonianValue)
	}
	insideHorizon := false
	for properTime := 0.0; properTime < 60 && !insideHorizon; properTime += 1 {
		camera, insideHorizon = hole.Advance(camera, 1)
		position := camera.Space()
		if hole.R(position) > hole.Horizon()*1.05 {
			frame := hole.FrameOf(position, hole.Raise(position, camera.P))
			metric := hole.Metric(position)
			if normValue := Dot(mul(metric, frame[1]), frame[1]); !(math.Abs(normValue-1) <= 1e-6) {
				t.Fatalf("the frame's first axis has length %v at %v", normValue, position)
			}
			if orthogonalityValue := Dot(mul(metric, frame[0]), frame[2]); !(math.Abs(orthogonalityValue) <= 1e-6) {
				t.Fatalf("the frame's time and second axes are not square at %v: %v", position, orthogonalityValue)
			}
		}
	}
	if !insideHorizon {
		t.Error("the camera did not reach the horizon in sixty M of proper time")
	}
}

// Falling drops its body at the radius it is given, on the x axis: round
// a spinning hole that is not the same number as x.
func TestFallingStartsAtTheGivenRadius(t *testing.T) {
	for _, spin := range []float64{0, 0.6, 0.99} {
		hole := Hole{M: 1, A: spin}
		body := hole.Falling(10)
		if radius := hole.R(body.Space()); !(math.Abs(radius-10) <= 1e-12) {
			t.Errorf("spin %v: Falling(10) starts at radius %v", spin, radius)
		}
		if body.X[2] != 0 || body.X[3] != 0 || !(body.X[1] > 0) {
			t.Errorf("spin %v: Falling(10) starts at %v, off the x axis", spin, body.Space())
		}
		if body != hole.Dropped(body.Space()) {
			t.Errorf("spin %v: Falling is not Dropped on the x axis", spin)
		}
	}
}

// A body is carried no further than the chart goes, however much proper
// time it is given: it stops inside, still a body of unit mass, and
// nothing comes out the far side. Round a still hole that is at the
// singularity; round a spinning one at the ring or the inner horizon,
// from the equator, from off it and down the axis.
func TestAdvanceStopsWhereTheChartEnds(t *testing.T) {
	starts := []Vec3{{10, 0, 0}, {6, 0, 5}, {0, 0, 10}, {8, 3, 0.5}}
	for _, spin := range []float64{0, 0.1, 0.6, 0.99} {
		hole := Hole{M: 1, A: spin}
		for _, start := range starts {
			for _, properTime := range []float64{60, 400} {
				body, inside := hole.Advance(hole.Dropped(start), properTime)
				radius := hole.R(body.Space())
				if !inside || !(radius < hole.Horizon()) {
					t.Errorf("spin %v, from %v, after %v of proper time: inside %v at radius %v", spin, start, properTime, inside, radius)
				}
				if drift := math.Abs(hole.H(body) + 0.5); !(drift < 1e-4) {
					t.Errorf("spin %v, from %v, after %v of proper time: H is %v from minus a half, so this is no longer a body", spin, start, properTime, drift)
				}
				if spin == 0 && !(radius <= 2*singularityMargin) {
					t.Errorf("round a still hole, from %v, the body stopped at radius %v, short of the singularity", start, radius)
				}
				if spin == 0.6 && !(radius < hole.innerHorizon() || hole.ringDistance(body.Space()) < singularityMargin) {
					t.Errorf("spin %v, from %v: the body stopped at radius %v, outside the inner horizon and off the ring", spin, start, radius)
				}
			}
		}
	}
}

// From rest at ten M the singularity is pi over two times the root of
// five hundred away in proper time, 35.12: the body is still falling at
// 35 and stopped at 36, and at 30 it has not reached the horizon.
func TestTheFallTakesTheTimeItShould(t *testing.T) {
	hole := Hole{M: 1}
	timeToSingularity := math.Pi / 2 * math.Sqrt(10*10*10/2.0)
	before, _ := hole.Advance(hole.Falling(10), timeToSingularity-0.1)
	after, _ := hole.Advance(hole.Falling(10), timeToSingularity+0.1)
	if radius := hole.R(before.Space()); !(radius > 0.2 && radius < 0.6) {
		t.Errorf("a tenth short of the singularity the body is at radius %v", radius)
	}
	if radius := hole.R(after.Space()); !(radius <= 2*singularityMargin) {
		t.Errorf("a tenth past the singularity the body is at radius %v", radius)
	}
	if body, inside := hole.Advance(hole.Falling(10), 30); inside || !(hole.R(body.Space()) >= 4) {
		t.Errorf("after 30 of proper time: inside %v at radius %v", inside, hole.R(body.Space()))
	}
}

// turning is how far round the hole a ray at an impact parameter
// turns before it escapes, in radians, for a hole with no spin.
func turning(hole Hole, impactParameter float64) float64 {
	position := Vec3{-2000, impactParameter, 0}
	frame := hole.FrameOf(position, hole.Static(position))
	ray := hole.Null(position, frame, Vec3{1, 0, 0})
	turned := 0.0
	lastAngle := math.Atan2(position[1], position[0])
	for i := 0; i < 400000; i++ {
		spacePoint := ray.Space()
		radius := hole.R(spacePoint)
		if radius > 2400 || radius < hole.Horizon()*1.001 {
			break
		}
		ray = hole.Step(ray, math.Min(hole.Pace(radius), 0.05))
		newPosition := ray.Space()
		currentAngle := math.Atan2(newPosition[1], newPosition[0])
		angleChange := currentAngle - lastAngle
		for angleChange > math.Pi {
			angleChange -= 2 * math.Pi
		}
		for angleChange < -math.Pi {
			angleChange += 2 * math.Pi
		}
		turned += angleChange
		lastAngle = currentAngle
	}
	return math.Abs(turned)
}

// The photon ring's subrings thin by e to the minus pi: the impact
// parameter at which a ray turns one more half circle before escaping
// lies closer to the shadow's edge by that factor each time.
func TestTheSubringsThinByEToTheMinusPi(t *testing.T) {
	hole := Hole{M: 1}
	shadowEdge := hole.Shadow()
	// the impact parameter at which the turning first reaches a given
	// angle, by bisection between the edge and a wide pass
	findParameterAtAngle := func(angle float64) float64 {
		low, high := shadowEdge+1e-7, shadowEdge+3
		for i := 0; i < 50; i++ {
			middle := (low + high) / 2
			if turning(hole, middle) > angle {
				low = middle
			} else {
				high = middle
			}
		}
		return (low + high) / 2
	}
	// a ray passing straight is turned by nought; one that goes half
	// round the hole by pi more; the orders are each pi further
	parameter1, parameter2 := findParameterAtAngle(1.5*math.Pi), findParameterAtAngle(2.5*math.Pi)
	ratio := (parameter2 - shadowEdge) / (parameter1 - shadowEdge)
	if !(math.Abs(ratio-math.Exp(-math.Pi)) <= 0.3*math.Exp(-math.Pi)) {
		t.Errorf("the second subring lies at %v of the first's distance from the edge; e to the minus pi is %v", ratio, math.Exp(-math.Pi))
	}
}

// mul is a matrix times a vector, for checking the metric the long way.
func mul(matrix Mat4, vector Vec4) Vec4 {
	var result Vec4
	for i := range 4 {
		for j := range 4 {
			result[i] += matrix[i][j] * vector[j]
		}
	}
	return result
}

// deriveSlow is hamilton's equations the long way, as the check on
// derive: the velocity by the whole inverse metric on the momentum, and
// the momentum's change by central differences of that metric.
func (h Hole) deriveSlow(r Ray) (velocity, momentumChange Vec4) {
	spacePoint := r.Space()
	velocity = mul(h.Inverse(spacePoint), r.P)
	epsilon := 1e-5 * math.Max(1, h.R(spacePoint))
	for i := range 3 {
		positionUp, positionDown := spacePoint, spacePoint
		positionUp[i] += epsilon
		positionDown[i] -= epsilon
		inverseMetricUp, inverseMetricDown := h.Inverse(positionUp), h.Inverse(positionDown)
		sum := 0.0
		for a := range 4 {
			for b := range 4 {
				sum += (inverseMetricUp[a][b] - inverseMetricDown[a][b]) * r.P[a] * r.P[b]
			}
		}
		momentumChange[i+1] = -0.5 * sum / (2 * epsilon)
	}
	return velocity, momentumChange
}

// The exact derivative of the metric agrees with central differences,
// with and without spin, at points on and off the equator, inside the
// ergosphere and out of it, for a momentum of light's kind and one that
// is no geodesic's at all: hamilton's equations do not ask.
func TestTheExactGradientMatchesTheDifferences(t *testing.T) {
	for _, hole := range []Hole{{M: 1}, {M: 1, A: 0.7}, {M: 2.5, A: 2}} {
		for _, position := range []Vec3{{5, 2, 1}, {3, -4, 0.5}, {-2, 1, 6}, {4, 0, 0}, {0, 0, 7}} {
			for _, momentum := range []Vec4{{-1, 0.3, 0.8, -0.5}, {0.7, -1.1, 0.2, 0.4}} {
				ray := Ray{X: Vec4{0, position[0], position[1], position[2]}, P: momentum}
				velocityExact, momentumChangeExact := hole.derive(ray)
				velocitySlow, momentumChangeSlow := hole.deriveSlow(ray)
				for i := range 4 {
					if !(math.Abs(velocityExact[i]-velocitySlow[i]) <= 1e-9) {
						t.Fatalf("hole %v at %v: velocity %v exact, %v by differences", hole, position, velocityExact, velocitySlow)
					}
					if !(math.Abs(momentumChangeExact[i]-momentumChangeSlow[i]) <= 1e-6*math.Max(1, math.Abs(momentumChangeSlow[i]))) {
						t.Fatalf("hole %v at %v: momentum change %v exact, %v by differences", hole, position, momentumChangeExact, momentumChangeSlow)
					}
				}
			}
		}
	}
}

// An arriving ray measures positive energies, says approaching gas is
// bluer and a static emitter deep in the well redder, and in a
// spinning hole runs the path a leaving ray runs in the hole spun the
// other way, which is time's reverse.
func TestTheArrivingRay(t *testing.T) {
	hole := Hole{M: 1}
	position := Vec3{20, 0, 0}
	frame := hole.FrameOf(position, hole.Static(position))
	ray := hole.Arriving(position, frame, Vec3{-1, 0.1, 0})
	if energy := Dot(ray.P, hole.Static(position)); !(energy > 0) {
		t.Fatalf("the arriving ray's energy is %v", energy)
	}
	// aimed at the hole, it goes inward
	if velocity := hole.Velocity(ray); !(velocity[1] < 0) || !(velocity[0] < 0) {
		t.Fatalf("the arriving ray heads %v: it should go inward and backward in time", velocity)
	}
	velocityParameter := 0.4
	gamma := 1 / math.Sqrt(1-velocityParameter*velocityParameter)
	flatHole := Hole{M: 1e-9}
	flatFrame := flatHole.FrameOf(position, flatHole.Static(position))
	flatRay := flatHole.Arriving(position, flatFrame, Vec3{1, 0, 0})
	approachingVelocity := Vec4{gamma, -gamma * velocityParameter, 0, 0}
	if shiftValue := Shift(flatRay.P, flatHole.Static(position), flatRay.P, approachingVelocity); !(math.Abs(shiftValue-1/(gamma*(1-velocityParameter))) <= 1e-6) {
		t.Errorf("approaching gas shifts by %v, not %v", shiftValue, 1/(gamma*(1-velocityParameter)))
	}
	// a static emitter at 4M seen by a static eye at 20M: the shift is
	// the ratio of the roots of one minus two over r
	incomingRay := hole.Arriving(position, frame, Vec3{-1, 0, 0})
	momentumAtEye := incomingRay.P
	for hole.R(incomingRay.Space()) > 4 {
		incomingRay = hole.Step(incomingRay, 0.01)
	}
	emitterPosition := incomingRay.Space()
	expectedShift := math.Sqrt(1-2/hole.R(emitterPosition)) / math.Sqrt(1-2/20.0)
	if shiftValue := Shift(momentumAtEye, hole.Static(position), incomingRay.P, hole.Static(emitterPosition)); !(math.Abs(shiftValue-expectedShift) <= 2e-3) {
		t.Errorf("light from 4M to 20M shifts by %v, not %v", shiftValue, expectedShift)
	}
	// the round trip: the arriving ray marched in, then the photon it
	// stands for, its momentum negated, marched back out, comes home
	// with its momentum the arriving ray's negated, and the observer
	// sees it come from direction
	spinningHole := Hole{M: 1, A: 0.9}
	spinningFrame := spinningHole.FrameOf(position, spinningHole.Static(position))
	direction := Vec3{-1, 0.25, 0.1}
	arrivingRay := spinningHole.Arriving(position, spinningFrame, direction)
	home := arrivingRay
	for i := 0; i < 400; i++ {
		arrivingRay = spinningHole.Step(arrivingRay, 0.05)
	}
	if !(spinningHole.R(arrivingRay.Space()) <= 19) {
		t.Fatalf("the arriving ray got nowhere: at %v", arrivingRay.Space())
	}
	photon := Ray{X: arrivingRay.X, P: Vec4{-arrivingRay.P[0], -arrivingRay.P[1], -arrivingRay.P[2], -arrivingRay.P[3]}}
	for i := 0; i < 400; i++ {
		photon = spinningHole.Step(photon, 0.05)
	}
	for i := range 3 {
		if !(math.Abs(photon.Space()[i]-position[i]) <= 1e-4) {
			t.Fatalf("the photon came back to %v, not %v", photon.Space(), position)
		}
		if !(math.Abs(photon.P[i+1]+home.P[i+1]) <= 1e-4) {
			t.Fatalf("the photon came back with momentum %v, the arriving ray left with %v", photon.P, home.P)
		}
	}
	observedDirection := spinningHole.Seen(spinningFrame, photon)
	directionNorm := math.Sqrt(direction[0]*direction[0] + direction[1]*direction[1] + direction[2]*direction[2])
	for i := range 3 {
		if !(math.Abs(observedDirection[i]-direction[i]/directionNorm) <= 1e-3) {
			t.Fatalf("the observer sees the photon come from %v, not %v", observedDirection, direction)
		}
	}
	// round a still hole the arriving ray and the leaving ray run the
	// same path; round a spinning one they part
	for _, spin := range []float64{0, 0.9} {
		stillHole := Hole{M: 1, A: spin}
		stillFrame := stillHole.FrameOf(position, stillHole.Static(position))
		incomingRay, leavingRay := stillHole.Arriving(position, stillFrame, direction), stillHole.Null(position, stillFrame, direction)
		for i := 0; i < 400; i++ {
			incomingRay, leavingRay = stillHole.Step(incomingRay, 0.05), stillHole.Step(leavingRay, 0.05)
		}
		maxSeparation := 0.0
		for i := range 3 {
			maxSeparation = math.Max(maxSeparation, math.Abs(incomingRay.Space()[i]-leavingRay.Space()[i]))
		}
		if spin == 0 && !(maxSeparation <= 1e-6) {
			t.Errorf("round a still hole the rays part by %v", maxSeparation)
		}
		if spin != 0 && !(maxSeparation >= 1e-2) {
			t.Errorf("round a spinning hole the rays part by only %v", maxSeparation)
		}
	}
}

// An arriving ray aimed at the hole, marched at PaceIn, is found
// hovering at the horizon in a few hundred steps with its hamiltonian
// still nought to a percent of the energy squared, which is Pace's own
// drift close in; from inside the horizon, the rays that fell in cross
// back out and the rest hover.
func TestAnArrivingRayLandsOnTheHorizon(t *testing.T) {
	hole := Hole{M: 1}
	position := Vec3{0, -40, 0}
	frame := hole.FrameOf(position, hole.Static(position))
	for _, impactParameter := range []float64{0, 4} {
		ray := hole.Arriving(position, frame, Vec3{impactParameter / 40, 1, 0})
		steps := 0
		for !hole.Hovering(ray) && steps < 4000 {
			ray = hole.Step(ray, hole.PaceIn(ray))
			steps++
		}
		if radius := hole.R(ray.Space()); !hole.Hovering(ray) || radius > hole.Horizon()*1.02 {
			t.Errorf("b=%v: the ray is at %v after %d steps", impactParameter, radius, steps)
		}
		if hamiltonianValue := math.Abs(hole.H(ray)); !(hamiltonianValue <= 1e-2) {
			t.Errorf("b=%v: H is %v at the horizon", impactParameter, hamiltonianValue)
		}
		if steps > 1000 {
			t.Errorf("b=%v: %d steps to the horizon", impactParameter, steps)
		}
	}
	// a camera that fell to one and a half masses, looking about
	camera := hole.Falling(20)
	for hole.R(camera.Space()) > 1.5 {
		camera, _ = hole.Advance(camera, 0.1)
	}
	cameraPosition := camera.Space()
	cameraFrame := hole.FrameOf(cameraPosition, hole.Velocity(camera))
	raysCrossOut, raysHover := 0, 0
	for _, direction := range []Vec3{{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0.7, 0.7, 0}, {-0.7, 0.7, 0}} {
		ray := hole.Arriving(cameraPosition, cameraFrame, direction)
		for i := 0; i < 4000 && !hole.Hovering(ray) && hole.R(ray.Space()) < 5; i++ {
			ray = hole.Step(ray, hole.PaceIn(ray))
		}
		switch {
		case hole.R(ray.Space()) >= 5:
			raysCrossOut++
		case hole.Hovering(ray):
			raysHover++
		default:
			t.Errorf("from inside, the ray toward %v is at %v, neither out nor hovering", direction, hole.R(ray.Space()))
		}
	}
	if raysCrossOut == 0 || raysHover == 0 {
		t.Errorf("from inside the horizon %d rays cross out and %d hover", raysCrossOut, raysHover)
	}
}

// March is what a renderer calls for a pixel: an arriving ray aimed at
// the hole is captured, not carried through the horizon and out the far
// side, and one aimed away escapes, from an eye near or far.
func TestMarchOfAnArrivingRayAtTheShadow(t *testing.T) {
	hole := Hole{M: 1}
	for _, eyeRadius := range []float64{2.5, 6, 40} {
		eye := Vec3{eyeRadius, 0, 0}
		frame := hole.FrameOf(eye, hole.Static(eye))
		if _, fate := hole.March(hole.Arriving(eye, frame, Vec3{-1, 0, 0}), 1000, 20000); fate != Captured {
			t.Errorf("eye at %v, looking at the hole: fate %v, want Captured", eyeRadius, fate)
		}
		if _, fate := hole.March(hole.Arriving(eye, frame, Vec3{1, 0, 0}), 1000, 20000); fate != Escaped {
			t.Errorf("eye at %v, looking away: fate %v, want Escaped", eyeRadius, fate)
		}
	}
	// seventeen degrees off the hole from six M is still inside the shadow
	eye := Vec3{6, 0, 0}
	frame := hole.FrameOf(eye, hole.Static(eye))
	offAxis := Vec3{-math.Cos(17 * math.Pi / 180), math.Sin(17 * math.Pi / 180), 0}
	if _, fate := hole.March(hole.Arriving(eye, frame, offAxis), 1000, 20000); fate != Captured {
		t.Errorf("seventeen degrees off the hole from six M: fate %v, want Captured", fate)
	}
}

// The circular photon orbits in the equator are where r squared, less
// three M r, plus or less twice A times the root of M r, is nought:
// plus for prograde. Three M at nought spin, M and four M at full; and
// the prograde one is the edge of where matter can go round at kepler's
// rate at all.
func TestPhotonOrbit(t *testing.T) {
	still := Hole{M: 1}
	if prograde, retrograde := still.PhotonOrbit(true), still.PhotonOrbit(false); prograde != 3 || retrograde != 3 || prograde != still.PhotonSphere() {
		t.Errorf("at nought spin the photon orbits are at %v and %v", prograde, retrograde)
	}
	full := Hole{M: 1, A: 0.999999}
	if prograde, retrograde := full.PhotonOrbit(true), full.PhotonOrbit(false); !(math.Abs(prograde-1) <= 0.01) || !(math.Abs(retrograde-4) <= 0.001) {
		t.Errorf("at full spin the photon orbits are at %v and %v", prograde, retrograde)
	}
	for _, hole := range []Hole{{M: 1, A: 0.3}, {M: 1, A: 0.6}, {M: 1, A: 0.9}, {M: 1, A: 0.99}, {M: 3, A: 1.5}} {
		for _, prograde := range []bool{true, false} {
			sense := -1.0
			if prograde {
				sense = 1
			}
			radius := hole.PhotonOrbit(prograde)
			if residue := radius*radius - 3*hole.M*radius + sense*2*hole.A*math.Sqrt(hole.M*radius); !(math.Abs(residue) <= 1e-12*hole.M*hole.M) {
				t.Errorf("hole %v, prograde %v: the photon orbit at %v leaves its equation %v from nought", hole, prograde, radius, residue)
			}
		}
		if !(hole.Horizon() < hole.PhotonOrbit(true) && hole.PhotonOrbit(true) < 3*hole.M && 3*hole.M < hole.PhotonOrbit(false) && hole.PhotonOrbit(true) < hole.ISCO(true)) {
			t.Errorf("hole %v: horizon %v, photon orbits %v and %v, innermost stable orbit %v are out of order", hole, hole.Horizon(), hole.PhotonOrbit(true), hole.PhotonOrbit(false), hole.ISCO(true))
		}
		justOutside, justInside := onEquator(hole, 1.001*hole.PhotonOrbit(true)), onEquator(hole, 0.999*hole.PhotonOrbit(true))
		if !isUnitTimelike(hole, justOutside, hole.Circular(justOutside)) || !math.IsNaN(hole.Circular(justInside)[0]) {
			t.Errorf("hole %v: matter can go round just inside the photon orbit at %v, or cannot just outside", hole, hole.PhotonOrbit(true))
		}
	}
}

// marchedToTheEnd is a ray's fate the long way: stepped by hand at
// PaceIn until it is far, or hovering, or light sent forward has gone
// inside the horizon, with no giving up before. It is the check
// on March, which gives a ray up as soon as it cannot come back.
func marchedToTheEnd(hole Hole, ray Ray, farRadius float64, steps int) Fate {
	for range steps {
		radius := hole.R(ray.Space())
		if radius > farRadius {
			return Escaped
		}
		if hole.Hovering(ray) || (hole.timeVelocity(ray) > 0 && radius < hole.Horizon()) {
			return Captured
		}
		ray = hole.Step(ray, hole.PaceIn(ray))
	}
	return Stuck
}

// A ray that is inside the photon orbit belonging to its own angular
// momentum, and still going inward, is not coming back, and March gives
// it up there and does not creep it down to the horizon. That changes
// no fate: looking every way at once, at four spins, from far, from
// near, from over the pole and from just outside the innermost photon
// orbit, March and the long way agree on every ray.
func TestGivingARayUpEarlyChangesNoFate(t *testing.T) {
	// directions spread evenly over the whole sky, on a spiral
	const count = 240
	var sky []Vec3
	for i := range count {
		height := 1 - 2*(float64(i)+0.5)/count
		turn := float64(i) * math.Pi * (3 - math.Sqrt(5))
		across := math.Sqrt(1 - height*height)
		sky = append(sky, Vec3{across * math.Cos(turn), across * math.Sin(turn), height})
	}
	for _, spin := range []float64{0, 0.5, 0.9, 0.99} {
		hole := Hole{M: 1, A: spin}
		for _, eye := range []Vec3{{12 * math.Cos(0.17), 0, 12 * math.Sin(0.17)}, {0.5, 0.3, 7}, {-6, 3, -2}, onEquator(hole, math.Max(2.02, 1.02*hole.PhotonOrbit(true)))} {
			frame := hole.FrameOf(eye, hole.Static(eye))
			fates := map[Fate]int{}
			for _, direction := range sky {
				ray := hole.Arriving(eye, frame, direction)
				_, fate := hole.March(ray, 200, 20000)
				fates[fate]++
				if want := marchedToTheEnd(hole, ray, 200, 20000); fate != want {
					t.Errorf("spin %v, eye at %v, looking toward %v: March says %v and the long way says %v", spin, eye, direction, fate, want)
				}
			}
			if fates[Captured] < 3 || fates[Escaped] < 3 || fates[Stuck] != 0 {
				t.Errorf("spin %v, eye at %v: the fates were %v, which does not have both kinds in it", spin, eye, fates)
			}
		}
	}
}

// Going inward is losing radius: for points and velocities with no
// nought in them, on and off the equator and near the axis, at three
// spins, it is what the radius does a short way along the velocity.
func TestGoingInwardIsLosingRadius(t *testing.T) {
	points := []Vec3{{7, 2, 3}, {-3, 4, -6}, {2.4, 0.3, 0.2}, {0.2, -0.1, 9}, {1.3, 0.4, -0.9}, {12, -5, 0.01}}
	velocities := []Vec3{{1, 0.3, -0.2}, {-0.4, 1, 0.7}, {-0.6, -0.5, 0.9}, {0.2, -0.9, -1}, {-1, 0.1, 0.05}, {0.05, 0.02, -1}, {0.05, 0.02, 1}}
	inward, outward := 0, 0
	for _, hole := range []Hole{{M: 1}, {M: 1, A: 0.6}, {M: 1, A: 0.99}} {
		for _, point := range points {
			for _, velocity := range velocities {
				radius := hole.R(point)
				const short = 1e-6
				along := Vec3{point[0] + short*velocity[0], point[1] + short*velocity[1], point[2] + short*velocity[2]}
				change := hole.R(along) - radius
				if math.Abs(change) < 1e-9 {
					continue // across the radius, near enough: either answer will do
				}
				got := hole.goingInward(quad{0, point[0], point[1], point[2]}, quad{1, velocity[0], velocity[1], velocity[2]}, radius)
				if got != (change < 0) {
					t.Errorf("hole %v at %v moving %v: going inward is %v, and the radius changes by %v", hole, point, velocity, got, change)
				}
				if got {
					inward++
				} else {
					outward++
				}
			}
		}
	}
	if inward < 20 || outward < 20 {
		t.Errorf("%d cases going inward and %d outward, which tests one kind", inward, outward)
	}
}

// A captured ray is handed back where it was given up, and for light
// seen from far away that is inside its photon orbit and well outside
// the horizon: three M round a still hole, where the march stops and
// goes no nearer.
func TestACapturedRayIsGivenUpAtItsPhotonOrbit(t *testing.T) {
	hole := Hole{M: 1}
	eye := Vec3{25, 0, 0}
	frame := hole.FrameOf(eye, hole.Static(eye))
	for _, direction := range []Vec3{{-1, 0, 0}, {-0.99, 0.1, 0.1}, {-0.98, -0.15, 0.1}} {
		end, fate := hole.March(hole.Arriving(eye, frame, direction), 200, 20000)
		if radius := hole.R(end.Space()); fate != Captured || !(radius < 3 && radius > 2.7) {
			t.Errorf("looking toward %v: fate %v, given up at radius %v, want captured just inside three M", direction, fate, radius)
		}
	}
	// the radius a ray cannot come back from is its own: round a spinning
	// hole it is the prograde photon orbit for light going round with
	// the hole as hard as light can, the retrograde one for light going
	// as hard against it, and between them for light over the poles
	spinning := Hole{M: 1, A: 0.9}
	turning := func(radius float64) float64 { // Lz over E of the spherical photon orbit at a radius
		return -(radius*radius*radius - 3*radius*radius + 0.81*radius + 0.81) / (0.9 * (radius - 1))
	}
	for _, radius := range []float64{spinning.PhotonOrbit(true), 2, 2.6, 3.2, spinning.PhotonOrbit(false)} {
		// any ray with that much angular momentum about the axis to its energy
		ray := Ray{X: Vec4{0, 10, 0, 0}, P: Vec4{-1, 0.2, turning(radius) / 10, 0.1}}
		if got := spinning.noReturn(quadOf(ray.X), quadOf(ray.P)); !(math.Abs(got-radius) <= 1e-9) {
			t.Errorf("light with Lz over E of %v cannot come back from inside %v; the photon orbit with that much is at %v", turning(radius), got, radius)
		}
	}
	// light with more angular momentum than any photon orbit has, or with
	// no energy, has no orbit of its own: the innermost serves
	for name, momentum := range map[string]Vec4{"too much angular momentum": {-1, 0, 3, 0}, "no energy": {0, 0.2, 0.3, 0.1}} {
		if got := spinning.noReturn(quad{0, 10, 0, 0}, quadOf(momentum)); got != spinning.PhotonOrbit(true) {
			t.Errorf("light with %s cannot come back from inside %v, want the innermost photon orbit at %v", name, got, spinning.PhotonOrbit(true))
		}
	}
}

// A photon's energy is its own business: light of a thousand times the
// frequency, or a thousandth, goes the same way, and a march takes it
// there by the same steps, to the same place and the same fate, with
// its momentum that many times as great. A camera moving fast or
// sitting deep makes rays whose energy is far from one.
func TestAPhotonsEnergyDoesNotChangeItsPath(t *testing.T) {
	hole := Hole{M: 1, A: 0.6}
	eye := Vec3{20, 3, 2}
	frame := hole.FrameOf(eye, hole.Static(eye))
	directions := map[string]Vec3{
		"at the hole":   {-0.985, -0.14, -0.1},
		"past the hole": {-0.94, 0.3, 0.16},
		"away from it":  {0.8, 0.36, 0.48},
	}
	fates := map[Fate]bool{}
	for name, direction := range directions {
		ray := hole.Arriving(eye, frame, direction)
		end, fate := hole.March(ray, 500, 20000)
		fates[fate] = true
		for _, scale := range []float64{1e-3, 0.3, 3, 1e3} {
			scaled := Ray{X: ray.X, P: Vec4{scale * ray.P[0], scale * ray.P[1], scale * ray.P[2], scale * ray.P[3]}}
			scaledEnd, scaledFate := hole.March(scaled, 500, 20000)
			if scaledFate != fate {
				t.Errorf("light looking %s had fate %v, and fate %v with its energy scaled by %v", name, fate, scaledFate, scale)
				continue
			}
			for k := range 4 {
				if !(math.Abs(scaledEnd.X[k]-end.X[k]) <= 1e-9*math.Max(1, math.Abs(end.X[k]))) || !(math.Abs(scaledEnd.P[k]-scale*end.P[k]) <= 1e-9*scale*math.Max(1, math.Abs(end.P[k]))) {
					t.Errorf("light looking %s ended at %v with momentum %v; with its energy scaled by %v it ended at %v with %v", name, end.X, end.P, scale, scaledEnd.X, scaledEnd.P)
					break
				}
			}
		}
	}
	if !fates[Captured] || !fates[Escaped] {
		t.Errorf("the rays' fates were %v, which leaves one kind untested", fates)
	}
}

// A camera that has fallen through the horizon still sees out: the
// light that fell in after it. March carries such a ray back out across
// the horizon to the sky, gives up the rays that hover, and agrees
// pixel for pixel with a march made by hand at PaceIn; light sent
// forward from inside is captured where it stands.
func TestMarchFromInsideTheHorizon(t *testing.T) {
	for _, spin := range []float64{0, 0.6} {
		hole := Hole{M: 1, A: spin}
		camera := hole.Falling(20)
		for hole.R(camera.Space()) > 0.75*hole.Horizon() {
			camera, _ = hole.Advance(camera, 0.1)
		}
		position := camera.Space()
		frame := hole.FrameOf(position, hole.Velocity(camera))
		sky, black := 0, 0
		for _, direction := range []Vec3{{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0.7, 0.7, 0}, {-0.7, 0.7, 0}, {0.6, 0, 0.8}, {0, -0.6, -0.8}} {
			arriving := hole.Arriving(position, frame, direction)
			byHand := arriving
			for i := 0; i < 20000 && !hole.Hovering(byHand) && hole.R(byHand.Space()) < 50; i++ {
				byHand = hole.Step(byHand, hole.PaceIn(byHand))
			}
			want := Captured
			if !hole.Hovering(byHand) {
				want = Escaped
			}
			if _, fate := hole.March(arriving, 50, 20000); fate != want {
				t.Errorf("spin %v, from inside at %v, looking toward %v: March says %v and a march by hand says %v", spin, position, direction, fate, want)
			}
			if want == Escaped {
				sky++
			} else {
				black++
			}
			if _, fate := hole.March(hole.Null(position, frame, direction), 50, 20000); fate != Captured {
				t.Errorf("spin %v: light sent toward %v from inside the horizon had fate %v", spin, direction, fate)
			}
		}
		if sky == 0 || black == 0 {
			t.Errorf("spin %v: from inside, %d rays reach the sky and %d are black, which tests one kind only", spin, sky, black)
		}
	}
}

// Hovering is a thing arriving rays do. Light leaving outward from just
// outside the horizon has its time run fast too, forward, and it
// escapes, from a thousandth of the horizon's radius outside it and
// from a hundredth of that: March must not call it captured. Light sent
// inward from the same places, and light sent any way from inside, is.
func TestMarchOfALeavingRayFromTheEdge(t *testing.T) {
	hole := Hole{M: 1}
	for _, radius := range []float64{2.003, 2.0005, 2.00002} {
		position := Vec3{radius, 0, 0}
		frame := hole.FrameOf(position, hole.Static(position))
		leaving := hole.Null(position, frame, Vec3{1, 0, 0})
		if hole.Hovering(leaving) {
			t.Errorf("a leaving ray at %v is not hovering", radius)
		}
		if _, fate := hole.March(leaving, 100, 20000); fate != Escaped {
			t.Errorf("light leaving outward from %v: fate %v, want Escaped", radius, fate)
		}
		if _, fate := hole.March(hole.Null(position, frame, Vec3{-1, 0, 0}), 100, 20000); fate != Captured {
			t.Errorf("light sent inward from %v: fate %v, want Captured", radius, fate)
		}
	}
	inside := Vec3{1.5, 0, 0}
	for _, momentum := range []Vec4{{-1, 0.5, 0, 0}, {-1, -0.5, 0.3, 0}, {-1, 0, 0, 0.4}} {
		ray := Ray{X: Vec4{0, inside[0], inside[1], inside[2]}, P: momentum}
		if !(hole.Velocity(ray)[0] > 0) {
			t.Fatalf("the momentum %v does not run forward in time, and this test is for one that does", momentum)
		}
		if _, fate := hole.March(ray, 100, 20000); fate != Captured {
			t.Errorf("something sent forward from inside the horizon with momentum %v: fate %v, want Captured", momentum, fate)
		}
	}
}

// A ray that is not a number is stuck at once, not after the whole
// budget of steps: it is handed back as it came, its momentum unstepped.
func TestMarchDoesNotStepARayThatIsNotANumber(t *testing.T) {
	hole := Hole{M: 1}
	ray := Ray{X: Vec4{0, math.NaN(), 0, 0}, P: Vec4{-1, 1, 0, 0}}
	stopped, fate := hole.March(ray, 100, 10)
	if fate != Stuck {
		t.Errorf("fate %v, want Stuck", fate)
	}
	if stopped.P != ray.P {
		t.Errorf("the march stepped a ray that is not a number: momentum %v", stopped.P)
	}
}

// onEquator is the point on the x axis at a kerr-schild radius.
func onEquator(hole Hole, radius float64) Vec3 {
	return Vec3{math.Sqrt(radius*radius + hole.A*hole.A), 0, 0}
}

// isUnitTimelike is whether a four-velocity has length minus one in the
// metric at a point.
func isUnitTimelike(hole Hole, position Vec3, velocity Vec4) bool {
	return math.Abs(Dot(hole.Lower(position, velocity), velocity)+1) < 1e-9
}

// An observer can be held still outside the ergosurface and nowhere
// inside it: that is two M on the equator of a hole of any spin, and
// the horizon itself on the axis. Where nothing is still the answer is
// not a number.
func TestStaticStopsAtTheErgosurface(t *testing.T) {
	hole := Hole{M: 1, A: 0.9}
	for _, radius := range []float64{2.01, 2.5, 40} {
		position := onEquator(hole, radius)
		if static := hole.Static(position); !isUnitTimelike(hole, position, static) || static[1] != 0 || static[2] != 0 || static[3] != 0 {
			t.Errorf("on the equator at %v a still observer has velocity %v", radius, static)
		}
	}
	for _, radius := range []float64{1.99, 1.8, 1.5} {
		if static := hole.Static(onEquator(hole, radius)); !math.IsNaN(static[0]) {
			t.Errorf("on the equator at %v, inside the ergosurface and outside the horizon at %v, a still observer has velocity %v", radius, hole.Horizon(), static)
		}
	}
	onAxis := Vec3{0, 0, 1.5}
	if static := hole.Static(onAxis); !isUnitTimelike(hole, onAxis, static) {
		t.Errorf("on the axis at 1.5, outside the horizon, a still observer has velocity %v", static)
	}
	still := Hole{M: 1}
	if static := still.Static(Vec3{2.01, 0, 0}); !isUnitTimelike(still, Vec3{2.01, 0, 0}, static) {
		t.Errorf("round a still hole, just outside the horizon, a still observer has velocity %v", static)
	}
	if static := still.Static(Vec3{1.99, 0, 0}); !math.IsNaN(static[0]) {
		t.Errorf("round a still hole, inside the horizon, a still observer has velocity %v", static)
	}
}

// Matter on a circle goes round at kepler's rate with the hole's
// correction, and stays on its circle: with spin and without, a quarter
// of the way round it is at the radius it started at. Inside the
// ergosphere, where nothing is still, it still goes round.
func TestCircularIsACircle(t *testing.T) {
	for _, hole := range []Hole{{M: 1}, {M: 1, A: 0.6}, {M: 1, A: 0.9}, {M: 2, A: 1.2}} {
		for _, radiusInMasses := range []float64{1.8, 4, 10} {
			radius := radiusInMasses * hole.M
			if hole.A < 0.9*hole.M && radiusInMasses < 4 {
				continue // inside this hole's photon orbit
			}
			position := onEquator(hole, radius)
			circular := hole.Circular(position)
			if !isUnitTimelike(hole, position, circular) {
				t.Fatalf("hole %v at %v: circular velocity %v is not of unit length", hole, radius, circular)
			}
			rate := circular[2] / circular[0] / position[0]
			if want := math.Sqrt(hole.M) / (math.Pow(radius, 1.5) + hole.A*math.Sqrt(hole.M)); !(math.Abs(rate-want) <= 1e-12) {
				t.Errorf("hole %v at %v: matter goes round at %v; kepler's rate is %v", hole, radius, rate, want)
			}
			body := Ray{X: Vec4{0, position[0], position[1], position[2]}, P: hole.Lower(position, circular)}
			quarterTurnInSteps := int(math.Pi / 2 / rate / circular[0] / (hole.Pace(radius) / 8))
			for range quarterTurnInSteps {
				body = hole.Step(body, hole.Pace(radius)/8)
			}
			if drift := math.Abs(hole.R(body.Space())/radius - 1); !(drift <= 1e-5) {
				t.Errorf("hole %v: matter on a circle at %v is at %v a quarter of the way round", hole, radius, hole.R(body.Space()))
			}
			if angle := math.Atan2(body.X[2], body.X[1]); !(math.Abs(angle-math.Pi/2) <= 0.01) {
				t.Errorf("hole %v at %v: after a quarter turn's steps the matter has turned %v", hole, radius, angle)
			}
		}
	}
}

// Inside the photon orbit nothing goes round at kepler's rate, and the
// answer is not a number: not the still observer's velocity, which is
// some other observer's and would shade a disc there wrongly in
// silence. The points here are inside the prograde photon orbit, at
// 3 M and 2.35 M, and outside the ergosurface, so a still observer
// does exist at each.
func TestCircularInsideThePhotonOrbit(t *testing.T) {
	cases := []struct {
		hole   Hole
		radius float64
	}{
		{Hole{M: 1}, 2.5},
		{Hole{M: 1, A: 0.5}, 2.2},
	}
	for _, c := range cases {
		position := onEquator(c.hole, c.radius)
		if static := c.hole.Static(position); !isUnitTimelike(c.hole, position, static) {
			t.Fatalf("hole %v at %v: no still observer, so this tests nothing", c.hole, c.radius)
		}
		if circular := c.hole.Circular(position); !math.IsNaN(circular[0]) {
			t.Errorf("hole %v: at %v, inside the photon orbit, matter on a circle has velocity %v", c.hole, c.radius, circular)
		}
	}
}

// Kernel returns a safe default when radius is zero.
func TestKernelAtOrigin(t *testing.T) {
	// At the origin (0,0,0), R returns 0 and kernel should handle it safely
	kerrHole := Hole{M: 1, A: 0.5}
	originPosition := Vec3{0.0, 0.0, 0.0}
	scalar, vector := kerrHole.kernel(originPosition)
	// kernel handles radius==0 by returning 0 scalar and identity vector
	if scalar != 0 {
		t.Errorf("at origin, scalar should be 0, got %v", scalar)
	}
	if vector != (Vec4{1, 0, 0, 0}) {
		t.Errorf("at origin, vector should be {1,0,0,0}, got %v", vector)
	}
}

// Where the kernel has no value, hamilton's equations are flat
// space's: the velocity is the flat momentum and the momentum does not
// change.
func TestDeriveAtOrigin(t *testing.T) {
	kerrHole := Hole{M: 1, A: 0.5}
	ray := Ray{P: Vec4{-1, 0.3, -0.2, 0.7}}
	velocity, momentumChange := kerrHole.derive(ray)
	if velocity != (Vec4{1, 0.3, -0.2, 0.7}) {
		t.Errorf("at the origin the velocity is %v", velocity)
	}
	if momentumChange != (Vec4{}) {
		t.Errorf("at the origin the momentum changes by %v", momentumChange)
	}
}

// March with zero steps budget returns Stuck immediately.
func TestMarchWithZeroSteps(t *testing.T) {
	hole := Hole{M: 1}
	position := Vec3{20, 0, 0}
	frame := hole.FrameOf(position, hole.Static(position))
	ray := hole.Null(position, frame, Vec3{1, 0, 0})
	// March with 0 steps should return immediately with Stuck
	resultRay, fate := hole.March(ray, 60, 0)
	if fate != Stuck {
		t.Errorf("march with 0 steps should return Stuck, got %v", fate)
	}
	// The ray should be unchanged
	for i := range 4 {
		if resultRay.X[i] != ray.X[i] || resultRay.P[i] != ray.P[i] {
			t.Errorf("march with 0 steps should not change the ray")
		}
	}
}

// timelike is a body of unit mass at a point, moving at half the speed
// of light past the observer held still there.
func timelike(hole Hole, position Vec3) Ray {
	frame := hole.FrameOf(position, hole.Static(position))
	speedAlongX, speedAlongY := 0.4, 0.3
	gamma := 1 / math.Sqrt(1-speedAlongX*speedAlongX-speedAlongY*speedAlongY)
	var velocity Vec4
	for k := range 4 {
		velocity[k] = gamma * (frame[0][k] + speedAlongX*frame[1][k] + speedAlongY*frame[2][k])
	}
	return Ray{X: Vec4{0, position[0], position[1], position[2]}, P: hole.Lower(position, velocity)}
}

// carterByAngles is the carter constant as the textbooks write it, in
// boyer and lindquist's angle, with the angle taken and the division by
// its sine done: the check on Carter, which does neither.
func carterByAngles(hole Hole, ray Ray) float64 {
	position, momentum := ray.Space(), ray.P
	radius := hole.R(position)
	sinTheta := math.Sqrt((position[0]*position[0] + position[1]*position[1]) / (radius*radius + hole.A*hole.A))
	cosTheta := position[2] / radius
	thetaMomentum := cosTheta/sinTheta*(position[0]*momentum[1]+position[1]*momentum[2]) - radius*sinTheta*momentum[3]
	energy, angularMomentum, massSquared := hole.E(ray), hole.Lz(ray), -2*hole.H(ray)
	return thetaMomentum*thetaMomentum + cosTheta*cosTheta*(hole.A*hole.A*(massSquared-energy*energy)+angularMomentum*angularMomentum/(sinTheta*sinTheta))
}

// The carter constant is kept: along light and along a body, with and
// without spin, off the equator, on a pass inside six M and on the way
// out, it stays where it started while the radius changes many fold.
func TestCarterIsKept(t *testing.T) {
	position := Vec3{30, 8, 12}
	for _, spin := range []float64{0, 0.6, 0.9} {
		hole := Hole{M: 1, A: spin}
		frame := hole.FrameOf(position, hole.Static(position))
		rays := map[string]Ray{
			"light going out":  hole.Arriving(position, frame, Vec3{0.9, -0.3, -0.2}),
			"light passing by": hole.Null(position, frame, Vec3{-0.93, -0.33, -0.17}),
			"a body":           timelike(hole, position),
		}
		for name, ray := range rays {
			carterAtStart := hole.Carter(ray)
			if !(carterAtStart > 1) {
				t.Fatalf("spin %v, %s: carter is %v at the start, off the equator", spin, name, carterAtStart)
			}
			leastRadius, mostRadius := hole.R(ray.Space()), hole.R(ray.Space())
			for step := 0; step < 600; step++ {
				ray = hole.Step(ray, hole.Pace(hole.R(ray.Space())))
				radius := hole.R(ray.Space())
				leastRadius, mostRadius = math.Min(leastRadius, radius), math.Max(mostRadius, radius)
				if carter := hole.Carter(ray); !(math.Abs(carter-carterAtStart) <= 1e-7*carterAtStart) {
					t.Fatalf("spin %v, %s: carter went from %v to %v by step %d, at radius %v", spin, name, carterAtStart, carter, step, radius)
				}
			}
			if !(mostRadius >= 3*leastRadius) {
				t.Errorf("spin %v, %s: the radius only ran from %v to %v, which tests little", spin, name, leastRadius, mostRadius)
			}
		}
	}
}

// Carter is the textbook's constant: the same number as the form with
// the angle taken, at points off the axis, for light and for a body.
func TestCarterIsTheTextbookForm(t *testing.T) {
	for _, hole := range []Hole{{M: 1}, {M: 1, A: 0.6}, {M: 2.5, A: 2}} {
		for _, position := range []Vec3{{30, 8, 12}, {5, 2, 1}, {-3, 4, -6}, {0.5, 0.2, 7}} {
			frame := hole.FrameOf(position, hole.Static(position))
			for name, ray := range map[string]Ray{
				"light":  hole.Null(position, frame, Vec3{0.3, 0.8, -0.5}),
				"a body": timelike(hole, position),
			} {
				got, want := hole.Carter(ray), carterByAngles(hole, ray)
				if !(math.Abs(got-want) <= 1e-9*math.Max(1, math.Abs(want))) {
					t.Errorf("hole %v at %v, %s: carter %v, the textbook form %v", hole, position, name, got, want)
				}
			}
		}
	}
}

// At nought spin the carter constant is the square of the angular
// momentum less the square of its part about z.
func TestCarterAtNoughtSpin(t *testing.T) {
	hole := Hole{M: 1}
	ray := Ray{X: Vec4{0, 7, -3, 4}, P: Vec4{-1, 0.3, -0.2, 0.7}}
	position := ray.Space()
	angularMomentum := Vec3{
		position[1]*ray.P[3] - position[2]*ray.P[2],
		position[2]*ray.P[1] - position[0]*ray.P[3],
		position[0]*ray.P[2] - position[1]*ray.P[1],
	}
	want := angularMomentum[0]*angularMomentum[0] + angularMomentum[1]*angularMomentum[1]
	if got := hole.Carter(ray); !(math.Abs(got-want) <= 1e-12*want) {
		t.Errorf("carter at nought spin is %v; the angular momentum squared less its z part squared is %v", got, want)
	}
	if !(math.Abs(angularMomentum[2]-hole.Lz(ray)) <= 1e-12) {
		t.Errorf("Lz is %v, the angular momentum about z is %v", hole.Lz(ray), angularMomentum[2])
	}
}

// The spin axis is a point like any other: on it the carter constant is
// the limit the textbook form takes beside it, and a ray lying in the
// equator has none.
func TestCarterOnTheAxisAndTheEquator(t *testing.T) {
	hole := Hole{M: 1, A: 0.6}
	momentum := Vec4{-1, 0.3, -0.2, 0.7}
	onAxis := Ray{X: Vec4{0, 0, 0, 20}, P: momentum}
	radius := hole.R(onAxis.Space())
	limit := (radius*radius+hole.A*hole.A)*(momentum[1]*momentum[1]+momentum[2]*momentum[2]) + hole.A*hole.A*(-2*hole.H(onAxis)-momentum[0]*momentum[0])
	if got := hole.Carter(onAxis); !(math.Abs(got-limit) <= 1e-12*math.Abs(limit)) {
		t.Errorf("on the axis carter is %v, the limit is %v", got, limit)
	}
	besideAxis := Ray{X: Vec4{0, 1e-6, 2e-6, 20}, P: momentum}
	if on, beside := hole.Carter(onAxis), carterByAngles(hole, besideAxis); !(math.Abs(on-beside) <= 1e-5*math.Abs(on)) {
		t.Errorf("on the axis carter is %v, beside it the textbook form is %v", on, beside)
	}
	position := Vec3{8, 3, 0}
	frame := hole.FrameOf(position, hole.Static(position))
	inEquator := hole.Null(position, frame, Vec3{0.6, 0.8, 0})
	if got := hole.Carter(inEquator); got != 0 {
		t.Errorf("a ray lying in the equator has carter %v", got)
	}
}

// Seen undoes Arriving and Null: light built to arrive from a direction
// is seen coming from it, whether it is carried as the arriving ray or
// as the photon itself, by an observer held still and by one going
// round the hole.
func TestSeenIsWhereTheLightComesFrom(t *testing.T) {
	hole := Hole{M: 1, A: 0.6}
	position := Vec3{20, 0, 3}
	direction := Vec3{0.6, 0, 0.8}
	for name, observer := range map[string]Vec4{"still": hole.Static(position), "going round": hole.Circular(position)} {
		frame := hole.FrameOf(position, observer)
		arriving := hole.Arriving(position, frame, direction)
		photon := hole.Null(position, frame, Vec3{-direction[0], -direction[1], -direction[2]})
		for kind, ray := range map[string]Ray{"the arriving ray": arriving, "the photon": photon} {
			seen := hole.Seen(frame, ray)
			for i := range 3 {
				if !(math.Abs(seen[i]-direction[i]) <= 1e-12) {
					t.Errorf("an observer %s sees %s come from %v; it was built to come from %v", name, kind, seen, direction)
					break
				}
			}
		}
	}
}

// A ray with no momentum across the observer's axes is seen from no
// direction: nought, and not the nought over nought that normalising
// it would give.
func TestSeenOfARayWithNoDirection(t *testing.T) {
	hole := Hole{M: 1}
	position := Vec3{20, 0, 0}
	frame := hole.FrameOf(position, hole.Static(position))
	if seen := hole.Seen(frame, Ray{X: Vec4{0, 20, 0, 0}}); seen != (Vec3{}) {
		t.Errorf("a ray with no momentum is seen from %v", seen)
	}
}

// The impact parameter is how far from the hole a ray from far away
// would pass were it not bent, in whatever plane it lies. Round a still
// hole the ray keeps it, and it is captured when that is less than the
// shadow and not otherwise.
func TestImpactParameterInAnyPlane(t *testing.T) {
	hole := Hole{M: 1}
	planes := map[string]Vec3{"the xy plane": {0, 1, 0}, "the xz plane": {0, 0, 1}, "a tilted plane": {0, 0.6, 0.8}}
	for plane, across := range planes {
		for _, offset := range []float64{5.1, 5.3, 30} {
			start := Vec3{-2000, offset * across[1], offset * across[2]}
			frame := hole.FrameOf(start, hole.Static(start))
			ray := hole.Null(start, frame, Vec3{1, 0, 0})
			impactParameter := hole.ImpactParameter(ray)
			if !(math.Abs(impactParameter-offset) <= 0.005*offset) {
				t.Errorf("in %s a ray aimed to pass at %v has impact parameter %v", plane, offset, impactParameter)
			}
			end, fate := hole.March(ray, 2400, 200000)
			if kept := hole.ImpactParameter(end); !(math.Abs(kept-impactParameter) <= 1e-6*impactParameter) {
				t.Errorf("in %s the impact parameter went from %v to %v along the ray", plane, impactParameter, kept)
			}
			if captured, insideShadow := fate == Captured, impactParameter < hole.Shadow(); captured != insideShadow {
				t.Errorf("in %s a ray of impact parameter %v, shadow %v, had fate %v", plane, impactParameter, hole.Shadow(), fate)
			}
		}
	}
}

// Round a spinning hole the impact parameter of a ray from far away is
// still its distance from the hole's centre on the sky, and a ray down
// the spin axis has none.
func TestImpactParameterWithSpin(t *testing.T) {
	hole := Hole{M: 1, A: 0.9}
	start := Vec3{-2000, 30, 40}
	frame := hole.FrameOf(start, hole.Static(start))
	if impactParameter := hole.ImpactParameter(hole.Null(start, frame, Vec3{1, 0, 0})); !(math.Abs(impactParameter-50) <= 0.005*50) {
		t.Errorf("a ray aimed to pass fifty M from a spinning hole has impact parameter %v", impactParameter)
	}
	onAxis := Vec3{0, 0, 2000}
	axisFrame := hole.FrameOf(onAxis, hole.Static(onAxis))
	if impactParameter := hole.ImpactParameter(hole.Null(onAxis, axisFrame, Vec3{0, 0, -1})); !(math.Abs(impactParameter) < 1e-9) {
		t.Errorf("a ray down the spin axis has impact parameter %v", impactParameter)
	}
}

// Light with no energy has no finite impact parameter: the answer is
// infinite, and not the nought that would say it was aimed at the hole.
func TestImpactParameterOfARayWithNoEnergy(t *testing.T) {
	hole := Hole{M: 1}
	if impactParameter := hole.ImpactParameter(Ray{X: Vec4{0, 10, 0, 0}, P: Vec4{0, 0, 1, 0}}); !math.IsInf(impactParameter, 1) {
		t.Errorf("a ray with no energy has impact parameter %v", impactParameter)
	}
}

// Light passing a still hole comes nearest at the radius its impact
// parameter fixes, b squared being that radius cubed over the radius
// less two M: from just outside the photon sphere to far out. A ray on
// its way out is nearest where it starts, and a ray that falls in comes
// nearer than the photon orbit, where it is given up.
func TestClosestApproach(t *testing.T) {
	hole := Hole{M: 1}
	for _, nearest := range []float64{3.2, 4, 6, 10, 40} {
		impactParameter := math.Sqrt(nearest * nearest * nearest / (nearest - 2))
		start := Vec3{-2000, 0.6 * impactParameter, 0.8 * impactParameter}
		frame := hole.FrameOf(start, hole.Static(start))
		ray := hole.Null(start, frame, Vec3{1, 0, 0})
		// the start is not at infinity, so ask the ray what its impact
		// parameter came out as and hold the answer to that
		got := hole.ClosestApproach(ray, 100000)
		if kept := math.Sqrt(got * got * got / (got - 2)); !(math.Abs(kept-hole.ImpactParameter(ray)) <= 1e-3*kept) {
			t.Errorf("a ray of impact parameter %v came nearest at %v, which belongs to an impact parameter of %v", hole.ImpactParameter(ray), got, kept)
		}
		if !(math.Abs(got-nearest) <= 0.01*nearest) {
			t.Errorf("a ray aimed to come nearest at %v came nearest at %v", nearest, got)
		}
	}
	position := Vec3{20, 0, 0}
	frame := hole.FrameOf(position, hole.Static(position))
	if got := hole.ClosestApproach(hole.Null(position, frame, Vec3{1, 0, 0}), 10000); got != 20 {
		t.Errorf("a ray going straight out from 20 came nearest at %v", got)
	}
	for name, ray := range map[string]Ray{
		"light falling in":            hole.Null(position, frame, Vec3{-1, 0, 0}),
		"an arriving ray at the hole": hole.Arriving(position, frame, Vec3{-1, 0, 0}),
	} {
		if got := hole.ClosestApproach(ray, 10000); !(got < hole.PhotonOrbit(true)) {
			t.Errorf("%s came nearest at %v, outside the photon orbit it cannot come back from", name, got)
		}
	}
	if got := hole.ClosestApproach(hole.Null(position, frame, Vec3{-1, 0, 0}), 0); got != 20 {
		t.Errorf("with no steps to take, the nearest a ray from 20 came is %v", got)
	}
	// a ray that is not a number: the answer is not one either, and it
	// comes at once, where spending this budget would take ten seconds
	notANumber := Ray{X: Vec4{0, math.NaN(), 0, 0}, P: Vec4{-1, 1, 0, 0}}
	begun := time.Now()
	if got := hole.ClosestApproach(notANumber, 100_000_000); !math.IsNaN(got) {
		t.Errorf("a ray that is not a number came nearest at %v", got)
	}
	if spent := time.Since(begun); spent > time.Second {
		t.Errorf("a ray that is not a number took %v to give up", spent)
	}
}

// A ray is classed by where it goes and not by how many steps it was
// watched for: from a hundred M, light aimed at the hole plunges, light
// aimed to miss it or away from it escapes, and from close in the same;
// matter going round on a circle circulates.
func TestClassifyRay(t *testing.T) {
	hole := Hole{M: 1}
	far, near := Vec3{100, 0, 0}, Vec3{2.5, 0, 0}
	farFrame, nearFrame := hole.FrameOf(far, hole.Static(far)), hole.FrameOf(near, hole.Static(near))
	orbit := Vec3{10, 0, 0}
	cases := []struct {
		name string
		ray  Ray
		want RayClassification
	}{
		{"from a hundred M, straight in", hole.Null(far, farFrame, Vec3{-1, 0, 0}), RayPlunges},
		{"from a hundred M, in but missing by ten M", hole.Null(far, farFrame, Vec3{-0.995, math.Sqrt(1 - 0.995*0.995), 0}), RayEscapes},
		{"from a hundred M, straight out", hole.Null(far, farFrame, Vec3{1, 0, 0}), RayEscapes},
		{"from two and a half M, straight in", hole.Null(near, nearFrame, Vec3{-1, 0, 0}), RayPlunges},
		{"from two and a half M, straight out", hole.Null(near, nearFrame, Vec3{1, 0, 0}), RayEscapes},
		{"matter on a circle at ten M", Ray{X: Vec4{0, 10, 0, 0}, P: hole.Lower(orbit, hole.Circular(orbit))}, RayCirculates},
	}
	for _, c := range cases {
		if got := hole.ClassifyRay(c.ray); got != c.want {
			t.Errorf("%s: classed %v, want %v", c.name, got, c.want)
		}
	}
}

// The components of a covector on an observer's frame are along the
// four-velocity first and then along the frame's x, y and z axes: the
// frame's own axes, lowered, come out as the flat metric's rows, and
// light sent along a direction as minus one and that direction.
func TestToObserverFrame(t *testing.T) {
	hole := Hole{M: 1, A: 0.6}
	position := Vec3{10, 2, 3}
	frame := hole.FrameOf(position, hole.Circular(position))
	flat := [4]Vec4{{-1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}}
	for axis := range 4 {
		components := hole.ToObserverFrame(frame, hole.Lower(position, frame[axis]))
		for k := range 4 {
			if !(math.Abs(components[k]-flat[axis][k]) <= 1e-12) {
				t.Errorf("axis %d of the frame has components %v on the frame, want %v", axis, components, flat[axis])
				break
			}
		}
	}
	direction := Vec3{0.48, -0.6, 0.64}
	components := hole.ToObserverFrame(frame, hole.Null(position, frame, direction).P)
	want := Vec4{-1, direction[0], direction[1], direction[2]}
	for k := range 4 {
		if !(math.Abs(components[k]-want[k]) <= 1e-12) {
			t.Errorf("light sent along %v has components %v on the frame, want %v", direction, components, want)
			break
		}
	}
}
