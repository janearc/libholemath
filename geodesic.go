package holemath

import "math"

// Ray is a geodesic's state: where it is in spacetime and its covariant
// momentum. For light the momentum is null; for a body it is timelike
// with unit mass.
type Ray struct {
	X Vec4
	P Vec4
}

// Fate is how a march ended.
type Fate int

// The fates.
const (
	Escaped  Fate = iota // past the far radius, and straight from there
	Captured             // not coming back: see March
	// the step budget ran out, or the ray stopped being a number
	Stuck
)

// The march's tolerances. They are choices, not physics, and each is
// made once, here.
const (
	// paceFraction is Pace's step, in the affine parameter, as a part
	// of the radius it is taken at.
	paceFraction = 0.02
	// timeFraction is the time a march's step spans, PaceIn's, as a
	// part of the radius it is taken at. At three hundredths a march
	// costs what it did stepping a fiftieth of the radius in the affine
	// parameter, and its worst ray lands several times nearer.
	timeFraction = 0.03
	// hoverRatio is how many times its energy an arriving ray's time
	// velocity is when the march gives it up as come from the edge.
	hoverRatio = 1000
	// turnedBack is how far above its least radius a ray has risen, as
	// a part of that radius, when ClosestApproach takes it to be on its
	// way out and stops.
	turnedBack = 1e-3
)

// limits is the radii a march decides a ray's fate by, worked out once
// for the ray before it is stepped.
type limits struct {
	horizon float64 // inside it an arriving ray may be found hovering
	// inside it a ray still going inward is not coming back
	noReturn float64
}

// limitsOf is the limits for a ray.
func (h Hole) limitsOf(point, momentum quad) limits {
	return limits{
		horizon:  h.Horizon(),
		noReturn: h.noReturn(point, momentum),
	}
}

// noReturn is the radius inside which a ray that is still going inward
// is never coming back out: the radius of the photon orbit that has the
// ray's own angular momentum about the axis for its energy.
//
// Light can circle a spinning hole at any radius between the prograde
// and retrograde orbits of the equator, on an orbit that leaves the
// equator, and at each radius only with one ratio of Lz to E:
//
//	Lz/E = -(r^3 - 3 M r^2 + A^2 r + A^2 M) / (A (r - M))
//
// by bardeen. A ray keeps that ratio, and whatever else it has, the places
// it can turn about lie either side of the orbit with its ratio: outside it
// if it came from outside, and then it is not inside going inward; inside it
// if it came from the hole, and then going inward it is past its turn.
//
// So inside that radius and going inward there is no turn ahead. The
// radius is M + y, with y the one root above nought of the cubic the
// equation makes,
//
//	y^3 + (A^2 + A Lz/E - 3 M^2) y - 2 M (M^2 - A^2) = 0
//
// which is below nought at nought and convex above it, so newton's
// method from above the root falls onto it and does not pass, and it is
// run until it stops falling.
//
// The innermost photon orbit, the prograde one, is inside all the
// others and bounds the answer below. A ray with no energy has no
// ratio, and the innermost serves for it. At nought spin every answer
// is 3M.
func (h Hole) noReturn(point, momentum quad) float64 {
	innermost := h.PhotonOrbit(true)
	turning := (point.x*momentum.y - point.y*momentum.x) / -momentum.t
	if math.IsInf(turning, 0) || math.IsNaN(turning) {
		return innermost
	}
	linear := h.A*h.A + h.A*turning - 3*h.M*h.M
	constant := -2 * h.M * (h.M*h.M - h.A*h.A)
	// the root is below the root of the linear term's size plus the cube
	// root of the constant's
	y := math.Sqrt(math.Abs(linear)) + math.Cbrt(-constant)
	for {
		next := y - (y*y*y+linear*y+constant)/(3*y*y+linear)
		if !(next < y) {
			break
		}
		y = next
	}
	return math.Max(innermost, h.M+y)
}

// goingInward is whether a point moving at a velocity is losing radius:
// the sign of the radius's gradient on the velocity, the gradient being
// along (x r^2, y r^2, z (r^2 + A^2)).
func (h Hole) goingInward(point, velocity quad, radius float64) bool {
	radiusSquared := radius * radius
	return radiusSquared*(point.x*velocity.x+point.y*velocity.y)+
		(radiusSquared+h.A*h.A)*point.z*velocity.z < 0
}

// Space is the ray's point in space.
func (r Ray) Space() Vec3 { return space(r.X) }

// space is the part in space of a point or a vector in spacetime.
func space(v Vec4) Vec3 { return Vec3{v[1], v[2], v[3]} }

// H is hamilton's function, half the inverse metric on the momentum:
// nought along a null ray and minus a half along a body of unit mass,
// kept along any geodesic.
func (h Hole) H(r Ray) float64 {
	return 0.5 * Dot(r.P, h.Velocity(r))
}

// E is the energy a geodesic keeps, minus the time momentum.
func (h Hole) E(r Ray) float64 { return -r.P[0] }

// Lz is the angular momentum about the spin axis a geodesic keeps.
func (h Hole) Lz(r Ray) float64 { return r.X[1]*r.P[2] - r.X[2]*r.P[1] }

// Carter is the carter constant, which a geodesic keeps beside E and
// Lz. In the angle theta that the kerr-schild coordinates share with
// boyer and lindquist's, z = R cos theta, it is
//
//	Q = p_theta^2 + cos^2 theta (A^2 (mu^2 - E^2) + Lz^2 / sin^2 theta)
//
// with mu^2 minus twice H: nought for light, one for a body of unit
// mass. At nought spin it is the square of the angular momentum less
// Lz squared. On the disc inside the ring, where R is nought, it is not
// a number.
func (h Hole) Carter(r Ray) float64 {
	x := r.Space()
	radius := h.R(x)
	cosTheta := x[2] / radius
	energy, massSquared := h.E(r), -2*h.H(r)
	return h.acrossAxisSquared(x, r.P, radius) +
		cosTheta*cosTheta*h.A*h.A*(massSquared-energy*energy)
}

// acrossAxisSquared is the square of the angular momentum across the
// spin axis, for a momentum at a point of known radius: Lx squared plus
// Ly squared at nought spin, and with spin the first and last terms of
// Carter,
//
//	p_theta^2 + cos^2 theta Lz^2 / sin^2 theta
//
// It is worked as the squared length of one vector in the xy plane, so
// no angle is taken, nothing is divided by sin theta, and the spin axis
// is a point like any other. With spin a geodesic does not keep it;
// Carter is the part it keeps.
func (h Hole) acrossAxisSquared(x Vec3, p Vec4, radius float64) float64 {
	radiusSpinSum := radius*radius + h.A*h.A
	cosTheta := x[2] / radius
	acrossX := cosTheta*radiusSpinSum*p[1] - radius*p[3]*x[0]
	acrossY := cosTheta*radiusSpinSum*p[2] - radius*p[3]*x[1]
	return (acrossX*acrossX + acrossY*acrossY) / radiusSpinSum
}

// Velocity is the contravariant velocity along the geodesic, dx by the
// affine parameter: the inverse metric on the momentum.
func (h Hole) Velocity(r Ray) Vec4 { return h.Raise(r.Space(), r.P) }

// quad is four numbers as the march holds them, a point or a momentum or
// the rate of either.
//
// It is a struct of four where Vec4 is an array of four because the
// compiler keeps a struct of up to four numbers in registers, in a
// function and between functions, and sends an array through memory; with
// arrays that traffic was more than half the cost of a step.
type quad struct{ t, x, y, z float64 }

// quadOf is a Vec4 as a quad.
func quadOf(v Vec4) quad { return quad{v[0], v[1], v[2], v[3]} }

// vec4 is a quad as a Vec4.
func (q quad) vec4() Vec4 { return Vec4{q.t, q.x, q.y, q.z} }

// moved is a quad carried along a rate by a step.
func (q quad) moved(rate quad, stepSize float64) quad {
	return quad{
		q.t + rate.t*stepSize, q.x + rate.x*stepSize,
		q.y + rate.y*stepSize, q.z + rate.z*stepSize,
	}
}

// rungeKutta is the four rates of a runge-kutta step weighed together:
// six times the step's mean rate.
func rungeKutta(k1, k2, k3, k4 quad) quad {
	return quad{
		k1.t + 2*(k2.t+k3.t) + k4.t,
		k1.x + 2*(k2.x+k3.x) + k4.x,
		k1.y + 2*(k2.y+k3.y) + k4.y,
		k1.z + 2*(k2.z+k3.z) + k4.z,
	}
}

// derive is hamilton's equations: the velocity is the inverse metric on
// the momentum, and the momentum's change is minus half the gradient of
// the inverse metric on the momentum twice. Time is a symmetry, so the
// time momentum does not change and the energy is kept.
func (h Hole) derive(r Ray) (velocity, momentumChange Vec4) {
	point, momentum := quadOf(r.X), quadOf(r.P)
	pointRate, momentumRate := h.rateAt(
		point, momentum, h.radius(point.x, point.y, point.z),
	)
	return pointRate.vec4(), momentumRate.vec4()
}

// rateAt is hamilton's equations for a point and a momentum whose
// radius is already in hand.
//
// The inverse metric is flat minus scalar * vector * vector, so with
// reach the vector on the momentum, the hamiltonian is the flat one
// less half of scalar * reach squared, and
//
//	velocity        = flat momentum - scalar * reach * vector
//	momentum change = reach * (reach/2 * grad scalar + scalar * grad reach)
//
// Both gradients are a multiple of the radius's gradient plus a short
// vector of their own, so the whole of it is three divisions and no
// loop, written out in a line as a shader would have it. The kernel is
// worked here again and not called for, to share those divisions.
func (h Hole) rateAt(
	point, momentum quad, radius float64,
) (velocity, momentumChange quad) {
	spin := h.A
	radiusSquared, spinSquared := radius*radius, spin*spin
	radiusSpinSum := radiusSquared + spinSquared
	// r^4 + A^2 z^2, which is nought only on the ring and the disc
	denominator := radiusSquared*radiusSquared + spinSquared*point.z*point.z
	if denominator == 0 || radius == 0 {
		// flat space, as kernel answers there
		return quad{
			-momentum.t, momentum.x, momentum.y, momentum.z,
		}, quad{}
	}
	overRadius, overRadiusSpinSum, overDenominator :=
		1/radius, 1/radiusSpinSum, 1/denominator

	scalar := 2 * h.M * radius * radiusSquared * overDenominator
	vectorX := (radius*point.x + spin*point.y) * overRadiusSpinSum
	vectorY := (radius*point.y - spin*point.x) * overRadiusSpinSum
	vectorZ := point.z * overRadius
	reachInPlane := vectorX*momentum.x + vectorY*momentum.y
	reachAlongAxis := vectorZ * momentum.z
	reach := reachInPlane + reachAlongAxis - momentum.t
	pull := scalar * reach
	velocity = quad{
		pull - momentum.t, momentum.x - pull*vectorX,
		momentum.y - pull*vectorY, momentum.z - pull*vectorZ,
	}

	// the radius's gradient, from the quartic that defines it:
	// (x r^3, y r^3, z r (r^2 + A^2)) / (r^4 + A^2 z^2)
	radiusGradInPlane := radius * radiusSquared * overDenominator
	radiusGradZ := point.z * radius * radiusSpinSum * overDenominator

	// grad scalar / scalar = scalarAlongRadius * grad r -
	// scalarAlongAxis * zhat
	scalarAlongRadius := (3*spinSquared*point.z*point.z -
		radiusSquared*radiusSquared) * overRadius * overDenominator
	scalarAlongAxis := 2 * spinSquared * point.z * overDenominator
	// grad reach = reachAlongRadius * grad r + the vector's own turn with
	// the point
	reachAlongRadius := (point.x*momentum.x+point.y*momentum.y-
		2*radius*reachInPlane)*overRadiusSpinSum -
		reachAlongAxis*overRadius
	alongRadius := 0.5*reach*scalarAlongRadius + reachAlongRadius

	momentumChange = quad{
		0,
		pull * (alongRadius*point.x*radiusGradInPlane +
			(radius*momentum.x-spin*momentum.y)*overRadiusSpinSum),
		pull * (alongRadius*point.y*radiusGradInPlane +
			(radius*momentum.y+spin*momentum.x)*overRadiusSpinSum),
		pull * (alongRadius*radiusGradZ + momentum.z*overRadius -
			0.5*reach*scalarAlongAxis),
	}
	return velocity, momentumChange
}

// rateOf is rateAt for a point whose radius is not in hand.
func (h Hole) rateOf(point, momentum quad) (velocity, momentumChange quad) {
	return h.rateAt(point, momentum, h.radius(point.x, point.y, point.z))
}

// Step is one runge-kutta step of the affine parameter.
func (h Hole) Step(r Ray, stepSize float64) Ray {
	point, momentum := quadOf(r.X), quadOf(r.P)
	velocity, momentumChange := h.rateOf(point, momentum)
	point, momentum = h.stepFrom(
		point, momentum, velocity, momentumChange, stepSize,
	)
	return Ray{X: point.vec4(), P: momentum.vec4()}
}

// stepFrom is a runge-kutta step for a point and a momentum whose rates
// where they stand are already in hand: the first of the step's four,
// which March has worked out to size the step by.
func (h Hole) stepFrom(
	point, momentum, velocity1, momentumChange1 quad, stepSize float64,
) (quad, quad) {
	velocity2, momentumChange2 := h.rateOf(
		point.moved(velocity1, stepSize/2),
		momentum.moved(momentumChange1, stepSize/2),
	)
	velocity3, momentumChange3 := h.rateOf(
		point.moved(velocity2, stepSize/2),
		momentum.moved(momentumChange2, stepSize/2),
	)
	velocity4, momentumChange4 := h.rateOf(
		point.moved(velocity3, stepSize),
		momentum.moved(momentumChange3, stepSize),
	)
	pointStep := rungeKutta(velocity1, velocity2, velocity3, velocity4)
	momentumStep := rungeKutta(
		momentumChange1, momentumChange2,
		momentumChange3, momentumChange4,
	)
	return point.moved(pointStep, stepSize/6),
		momentum.moved(momentumStep, stepSize/6)
}

// Pace is a step's size at a radius: a fiftieth of the radius, so the
// steps are fine near the hole, where the bending is, and coarse far
// from it, where space is flat and a long step is exact.
//
// The step is in the affine parameter, which is about distance for a
// ray from far away whose energy is one, and proper time for a body. It
// is not distance for light of any other energy: that light goes as
// many times too far in a step as its energy is from one.
//
// PaceIn is the step for a ray, and is what March takes.
func (h Hole) Pace(r float64) float64 {
	return paceFraction * r
}

// PaceIn is the step for a ray where it stands: the step in which its
// time moves on, or back, by three hundredths of its radius.
//
// Far from the hole that is three hundredths of the radius in distance
// too, whatever the light's energy, since a photon's energy has no say in
// where it goes.
//
// Near the hole's edge it is less, and for an arriving ray that came
// from the edge it shrinks without end: that light left the edge in the
// far past, its time runs away as the ray comes to it, and so does its
// momentum, which a step sized by distance would overshoot.
//
// A ray crossing the horizon the other way, light that fell in, keeps its
// time in hand and crosses at a full step, as the chart lets.
func (h Hole) PaceIn(r Ray) float64 {
	return timeFraction * h.R(r.Space()) / math.Abs(h.timeVelocity(r))
}

// Hovering is whether an arriving ray has come from the hole's edge:
// its time runs away backward, a thousand times its energy's, and the
// march would never get there. Such a ray is black.
//
// Light leaving the edge outward has its time run as fast forward, and
// is not hovering: it gets out.
func (h Hole) Hovering(r Ray) bool {
	return runsAway(h.timeVelocity(r), r.P[0])
}

// timeVelocity is the time part of Velocity alone, dt by the affine
// parameter: above nought for a ray marched forward in time and below
// for an arriving one, and never nought, anywhere, since the chart's
// time is a time inside the horizon too.
func (h Hole) timeVelocity(r Ray) float64 {
	scalar, vector := h.kernel(r.Space())
	return -r.P[0] + scalar*(-r.P[0]+
		vector[1]*r.P[1]+vector[2]*r.P[2]+vector[3]*r.P[3])
}

// runsAway is whether a ray's time runs backward hoverRatio times
// faster than its energy.
func runsAway(timeVelocity, timeMomentum float64) bool {
	return -timeVelocity > hoverRatio*math.Abs(timeMomentum)
}

// March carries a ray until it is past far from the hole, captured, or
// out of steps, and says which; the ray handed back is where it
// stopped. It marches a leaving ray forward and an arriving ray back
// alike.
//
// A ray is captured as soon as it cannot come back: when it is inside
// the photon orbit that has its own angular momentum and is still going
// inward.
//
// For light seen from outside that is well short of the horizon, three M
// round a still hole, and the ray is handed back from there; that pixel is
// black.
//
// Light sent forward from inside the horizon is always going inward, and
// is captured where it starts. Light sent outward from just outside the
// horizon is not captured, however near: it gets out.
//
// An arriving ray is not captured for being inside the horizon: an eye
// that has fallen through still sees the light that fell in after it,
// and March carries that light back out to where it came from. Inside
// the horizon an arriving ray is captured when it is found hovering.
//
// Each step is a runge-kutta step of PaceIn's size, and a step that
// would cross the far radius is taken anyway, so the ray may end a
// little past it. A ray that is not a number is stuck at once.
func (h Hole) March(r Ray, farRadius float64, steps int) (Ray, Fate) {
	point, momentum := quadOf(r.X), quadOf(r.P)
	limits := h.limitsOf(point, momentum)
	fate := Stuck
	for i := 0; i < steps; i++ {
		radius := h.radius(point.x, point.y, point.z)
		if math.IsNaN(radius) {
			break
		}
		if radius > farRadius {
			fate = Escaped
			break
		}
		var givenUp bool
		if point, momentum, givenUp = h.pacedStep(
			point, momentum, radius, limits,
		); givenUp {
			fate = Captured
			break
		}
	}
	return Ray{X: point.vec4(), P: momentum.vec4()}, fate
}

// pacedStep is one step of a march, from a point whose radius is in
// hand: a runge-kutta step of PaceIn's size, or no step at all and true
// if the march gives the ray up there as captured.
//
// The time velocity that sizes the step and says whether the ray hovers
// is the time part of the step's first rate, and the velocity that says
// whether it is going inward is the rest of it, so the rate is worked
// out once here and the step is given it.
func (h Hole) pacedStep(
	point, momentum quad, radius float64, limits limits,
) (quad, quad, bool) {
	velocity, momentumChange := h.rateAt(point, momentum, radius)
	notComingBack := radius < limits.noReturn &&
		h.goingInward(point, velocity, radius)
	hovering := radius < limits.horizon && runsAway(velocity.t, momentum.t)
	if notComingBack || hovering {
		return point, momentum, true
	}
	point, momentum = h.stepFrom(
		point, momentum, velocity, momentumChange,
		timeFraction*radius/math.Abs(velocity.t),
	)
	return point, momentum, false
}

// ImpactParameter is how far from the hole's centre a ray of light from
// far away would pass were it not bent: its angular momentum over its
// energy, the whole angular momentum, across the spin axis and about
// it, so the answer is the same in whatever plane the ray lies.
//
// Round a still hole the ray keeps it, and is captured when it is less
// than Shadow.
//
// Round a spinning one it is the ray's distance from the centre on the sky
// of whoever is far along it, right in the limit of far away; there is no
// one radius to hold it against there, since capture is decided by Lz over
// E and Carter over E squared apart.
//
// A ray with no energy has none, and the answer is infinite.
func (h Hole) ImpactParameter(r Ray) float64 {
	x := r.Space()
	angularMomentum := h.Lz(r)
	squared := h.acrossAxisSquared(x, r.P, h.R(x)) +
		angularMomentum*angularMomentum
	return math.Sqrt(squared) / math.Abs(h.E(r))
}

// ClosestApproach is the least radius a ray comes to: it is marched, as
// March marches it, until it has turned back outward, been captured, or
// taken all its steps, and the least radius it stood at, the last
// included, is the answer.
//
// A ray on its way out from the start answers the radius it started at,
// and a captured ray the radius it was given up at, which is inside its
// photon orbit and says no more than that.
//
// The answer is the least of the radii stepped at, so it is over by up
// to the square of half a step over twice the radius: half a part in
// ten thousand. For a ray that is not a number it is not a number, at
// once.
func (h Hole) ClosestApproach(r Ray, steps int) float64 {
	point, momentum := quadOf(r.X), quadOf(r.P)
	limits := h.limitsOf(point, momentum)
	radius := h.radius(point.x, point.y, point.z)
	leastRadius := radius
	// go on only while the ray has not turned back; a radius that is not
	// a number fails the comparison, and that stops the march too
	for i := 0; i < steps && radius <= leastRadius*(1+turnedBack); i++ {
		var givenUp bool
		if point, momentum, givenUp = h.pacedStep(
			point, momentum, radius, limits,
		); givenUp {
			break
		}
		radius = h.radius(point.x, point.y, point.z)
		leastRadius = math.Min(leastRadius, radius)
	}
	return leastRadius
}

// DeflectionAngle is the angle between the way a ray sets out and the
// way it leaves past the far radius, in radians from nought to pi: how
// far the hole bent it, in whatever plane it was bent.
//
// The two ways are the ray's velocity in space, which is where it is
// going only where space is flat, so the start and the far radius MUST
// both be far from the hole.
//
// A ray that goes more than half round reads as the angle between its two
// ways, not as its turning. A ray that does not escape has no deflection,
// and the answer is not a number.
func (h Hole) DeflectionAngle(r Ray, farRadius float64, steps int) float64 {
	end, fate := h.March(r, farRadius, steps)
	if fate != Escaped {
		return math.NaN()
	}
	return angleBetween(space(h.Velocity(r)), space(h.Velocity(end)))
}

// angleBetween is the angle between two directions in space, from
// nought to pi, by the arctangent of the cross product's length over
// the dot product: it keeps its digits at small angles, where the
// arccosine of the dot product alone loses them.
func angleBetween(a, b Vec3) float64 {
	cross := Vec3{
		a[1]*b[2] - a[2]*b[1],
		a[2]*b[0] - a[0]*b[2],
		a[0]*b[1] - a[1]*b[0],
	}
	crossLength := math.Sqrt(
		cross[0]*cross[0] + cross[1]*cross[1] + cross[2]*cross[2],
	)
	return math.Atan2(crossLength, a[0]*b[0]+a[1]*b[1]+a[2]*b[2])
}
