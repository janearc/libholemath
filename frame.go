package holemath

import "math"

// Static is the four-velocity of an observer held still at a point: the
// time axis scaled to unit length.
//
// Nothing is still inside the ergosurface, where the hole's spin drags
// everything round with it: on the equator that is inside two M at any
// spin, on the axis it is the horizon, and at nought spin it is the
// horizon everywhere. There the answer is not a number.
func (h Hole) Static(x Vec3) Vec4 {
	return h.orbiting(x, 0)
}

// Circular is the four-velocity of matter going round the hole on a
// circle at a point, prograde, as a thin disc's does: kepler's rate
// with the hole's correction, which is the same rate by the clock at
// infinity in these coordinates as in boyer and lindquist's.
//
// Inside the prograde photon orbit nothing goes round at that rate, and
// the answer is not a number. It can be had inside the ergosphere,
// where Static cannot.
func (h Hole) Circular(x Vec3) Vec4 {
	radius := h.R(x)
	sqrtMass := math.Sqrt(h.M)
	return h.orbiting(x, sqrtMass/(radius*math.Sqrt(radius)+h.A*sqrtMass))
}

// orbiting is the four-velocity of an observer going round the spin
// axis at a point, at an angular velocity by the clock at infinity: the
// time axis and the turn about z in that proportion, scaled to unit
// length.
//
// Where nothing can go round at that rate, the answer is not a number.
func (h Hole) orbiting(x Vec3, angularVelocity float64) Vec4 {
	velocity := Vec4{1, -angularVelocity * x[1], angularVelocity * x[0], 0}
	// the root of a negative is not a number, which is the answer wanted
	scale := math.Sqrt(-1 / Dot(h.Lower(x, velocity), velocity))
	return Vec4{scale, scale * velocity[1], scale * velocity[2], 0}
}

// Frame is an observer's orthonormal tetrad at a point: the observer's
// four-velocity first, then three spatial axes square to it and to one
// another in the metric, made from the coordinate axes by gram and
// schmidt, x first.
type Frame [4]Vec4

// FrameOf is the tetrad of an observer with four-velocity u at x.
func (h Hole) FrameOf(x Vec3, u Vec4) Frame {
	scalar, null := h.kernel(x)
	innerProduct := func(a, b Vec4) float64 {
		return Dot(lowered(scalar, null, a), b)
	}
	var frame Frame
	frame[0] = u
	seeds := []Vec4{{0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}}
	for i, seed := range seeds {
		vector := seed
		// take out the parts along the axes so far; the time axis is
		// timelike, so its part comes in with the other sign
		vector = sub(vector, frame[0], -innerProduct(vector, frame[0]))
		for j := 1; j <= i; j++ {
			vector = sub(
				vector, frame[j],
				innerProduct(vector, frame[j]),
			)
		}
		norm := math.Sqrt(innerProduct(vector, vector))
		for k := range 4 {
			vector[k] /= norm
		}
		frame[i+1] = vector
	}
	return frame
}

// sub subtracts b times scale from a.
func sub(vec, subtract Vec4, scale float64) Vec4 {
	for i := range 4 {
		vec[i] -= subtract[i] * scale
	}
	return vec
}

// Null is a ray of light leaving x in the direction direction as the
// observer with tetrad frame sees it: direction is a direction in the
// observer's own space, its length not used.
//
// The ray's momentum is the observer's time axis plus the unit direction
// along the spatial axes, lowered, so that it is null.
//
// The ray starts at time nought, and its energy as that observer measures
// it is one.
func (h Hole) Null(x Vec3, frame Frame, direction Vec3) Ray {
	length := math.Sqrt(
		direction[0]*direction[0] + direction[1]*direction[1] +
			direction[2]*direction[2],
	)
	if length > 0 {
		direction = Vec3{
			direction[0] / length, direction[1] / length,
			direction[2] / length,
		}
	}
	var rayVector Vec4
	for k := range 4 {
		rayVector[k] = frame[0][k] +
			direction[0]*frame[1][k] +
			direction[1]*frame[2][k] +
			direction[2]*frame[3][k]
	}
	return Ray{X: Vec4{0, x[0], x[1], x[2]}, P: h.Lower(x, rayVector)}
}

// Seen is the direction an observer with tetrad frame sees a ray's light
// come from, a unit direction in the observer's own space, which is
// where the light appears on that observer's sky after aberration.
//
// It is the momentum's parts along the spatial axes, normalised, and
// turned to face the way the light came: about, for a photon's own
// momentum, and not at all for an arriving ray's, which Arriving has
// negated already.
//
// The momentum's part along the observer's time axis says which it is,
// below nought for a photon and above for an arriving ray, so Seen of the
// ray Arriving built from a direction is that direction.
func (h Hole) Seen(frame Frame, r Ray) Vec3 {
	components := h.ToObserverFrame(frame, r.P)
	direction := Vec3{components[1], components[2], components[3]}
	magnitude := math.Sqrt(
		direction[0]*direction[0] + direction[1]*direction[1] +
			direction[2]*direction[2],
	)
	if magnitude == 0 {
		return direction
	}
	// the length to divide by carries the facing as its sign
	signedLength := math.Copysign(magnitude, components[0])
	for i := range 3 {
		direction[i] /= signedLength
	}
	return direction
}

// Falling is a body let go from rest at a radius on the x axis, as a
// camera falling into the hole: its geodesic from there, timelike with
// unit mass, marched by Advance.
//
// The radius is R, the kerr-schild one; round a spinning hole that point's
// x is the root of R squared plus A squared.
func (h Hole) Falling(initialRadius float64) Ray {
	return h.Dropped(
		Vec3{math.Sqrt(initialRadius*initialRadius + h.A*h.A), 0, 0},
	)
}

// singularityMargin is how near the singularity Advance carries a body
// before it stops, as a part of M. It is a choice: the pace goes to
// nought there, and a step taken any nearer lands in the flat space
// the kernel falls back to at nought radius, on the far side.
const singularityMargin = 0.01

// ringDistance is how far a point is from the singularity, the ring of
// radius A in the equator. At nought spin the ring is the centre and
// the distance is R.
func (h Hole) ringDistance(x Vec3) float64 {
	fromRingInPlane := math.Sqrt(x[0]*x[0]+x[1]*x[1]) - h.A
	return math.Sqrt(fromRingInPlane*fromRingInPlane + x[2]*x[2])
}

// Advance carries a body's geodesic on by a proper time and says
// whether it is inside the horizon. Its steps are half the pace at the
// nearer of its radius and its distance from the ring: close to the
// ring the hole changes over that distance, which is much the shorter.
//
// Inside means inside the horizon itself, to the last digit: the margin
// a march of light stops at is for rays and is not used here.
//
// A body is carried no further than the chart goes.
//
// It stops, inside, with the rest of its proper time unspent, within
// singularityMargin of the ring or of the disc the ring spans, and at the
// inner horizon of a spinning hole: past that one a body off the equator
// turns back outward, and its way out is not in these coordinates.
func (h Hole) Advance(r Ray, properTime float64) (Ray, bool) {
	innerHorizon := h.innerHorizon()
	remainingProperTime := properTime
	for remainingProperTime > 0 {
		x := r.Space()
		radius := h.R(x)
		reach := math.Min(radius, h.ringDistance(x))
		if reach < singularityMargin*h.M || radius < innerHorizon {
			return r, true
		}
		stepSize := math.Min(remainingProperTime, h.Pace(reach)/2)
		r = h.Step(r, stepSize)
		remainingProperTime -= stepSize
	}
	return r, h.R(r.Space()) < h.Horizon()
}

// Arriving is the ray of light arriving at x from the direction direction as
// the observer with tetrad frame sees it, set up to be marched forward in
// its parameter, which runs backward in time: a picture is traced from the
// eye.
//
// Its momentum is the arriving photon's negated, so the energies an
// observer measures along it, Dot(P, u), come out positive. Round a
// still hole it runs the path of the ray of Null leaving toward direction;
// round a spinning one the two part, since the dragging turns them.
func (h Hole) Arriving(x Vec3, frame Frame, direction Vec3) Ray {
	ray := h.Null(
		x, frame, Vec3{-direction[0], -direction[1], -direction[2]},
	)
	for k := range 4 {
		ray.P[k] = -ray.P[k]
	}
	return ray
}

// Shift is the ratio of a photon's energy at the eye to its energy at
// the emitter, along a ray: its momentum at the eye on the observer's
// four-velocity there, over its momentum at the emitter on the emitter's
// four-velocity there.
//
// Over one the light arrives bluer, under one redder; doppler beaming and
// the hole's own redshift are both in it.
//
// The momentum is not the same at the two ends: a ray keeps only its
// time part. One momentum serves at both ends only where that end's
// observer is still, with nothing but a time part to its velocity.
func Shift(momentumAtEye, observer, momentumAtEmitter, emitter Vec4) float64 {
	return Dot(momentumAtEye, observer) / Dot(momentumAtEmitter, emitter)
}

// Dropped is a body let go from rest at a point, at time nought, as a
// camera is that falls in from anywhere; Falling is the same on the x
// axis.
func (h Hole) Dropped(x Vec3) Ray {
	return Ray{X: Vec4{0, x[0], x[1], x[2]}, P: h.Lower(x, h.Static(x))}
}

// RayClassification is where a ray or a body ends up, as ClassifyRay
// says it.
type RayClassification int

// The classifications.
const (
	RayEscapes RayClassification = iota // got far away
	RayPlunges                          // was captured
	// did neither in the steps it was given
	RayCirculates
)

// ClassifyRay's reach. They are choices, each made once, here.
const (
	// classifyFarFactor is how many times its starting radius out a ray
	// has got when it is classed as escaping.
	classifyFarFactor = 4
	// classifyFarFloor is the least radius that counts as far, in M,
	// for a ray that starts close in.
	classifyFarFloor = 40
	// classifySteps is the march's budget. Light from a hundred M takes
	// about a hundred and thirty steps to fall in; matter on a circle
	// at ten M goes round seven times in this many.
	classifySteps = 5000
)

// classOf is the classification each fate of a march stands for.
var classOf = [...]RayClassification{
	Escaped:  RayEscapes,
	Captured: RayPlunges,
	Stuck:    RayCirculates,
}

// ClassifyRay says whether a ray or a body escapes, plunges into the hole
// or circulates, by marching it and seeing where it goes.
//
// It escapes on getting classifyFarFactor times its starting radius out, or
// classifyFarFloor M if that is further; it circulates if it has done
// neither after classifySteps steps.
//
// A bound body whose orbit reaches past that radius is classed as
// escaping, and a ray that is not a number as circulating.
func (h Hole) ClassifyRay(r Ray) RayClassification {
	farRadius := math.Max(
		classifyFarFactor*h.R(r.Space()), classifyFarFloor*h.M,
	)
	_, fate := h.March(r, farRadius, classifySteps)
	return classOf[fate]
}

// ToObserverFrame is the components of a covariant vector on an
// observer's tetrad: along the four-velocity first, then along the
// frame's x, y and z axes. A contravariant vector is made covariant
// with Lower first.
//
// For a photon's momentum the first is minus the energy the observer
// measures, and the other three its direction of travel times that
// energy.
func (h Hole) ToObserverFrame(frame Frame, covector Vec4) Vec4 {
	return Vec4{
		Dot(covector, frame[0]),
		Dot(covector, frame[1]),
		Dot(covector, frame[2]),
		Dot(covector, frame[3]),
	}
}
