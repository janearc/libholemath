// Package holemath is the maths of a black hole for a renderer: a kerr
// hole in kerr-schild form, null and timelike geodesics marched by
// hamilton's equations, the quantities a geodesic keeps, the redshift
// between two observers, and the frame of a camera falling in.
//
// It imports nothing but the standard library, so the relativity can
// leave; what it does not know is anything about light or pictures.
//
// One file leans toward a renderer, gpu.go: the hole and a camera as
// the float32 constants a shader takes, and a step that reports what a
// ray keeps, to hold a shader's own steps against. There is no GPU in
// it and nothing of a picture.
//
// Units: G and c are one, and a hole of mass M has its lengths in M.
// The signature is minus, plus, plus, plus; a point is (t, x, y, z)
// and a momentum is covariant.
package holemath

import "math"

// Vec3 is a point or a direction in space.
type Vec3 [3]float64

// Vec4 is a point in spacetime, or a covariant or contravariant vector
// there; which one is said where it is used.
type Vec4 [4]float64

// Mat4 is a metric or its inverse.
type Mat4 [4][4]float64

// Hole is a black hole: its mass M and its spin A, the angular momentum
// per unit mass, from nought for schwarzschild up to just under M. The
// spin axis is z, and the hole turns from x toward y.
//
// Nothing here stops a Hole being made with any two numbers. Valid says
// whether they are a hole's.
type Hole struct {
	M, A float64
}

// Valid is whether the hole is one the other methods are right for: a
// finite mass above nought and a spin from nought to just under it.
//
// A spin at or above the mass is no black hole, since it has no
// horizon; Horizon answers M for it and ISCO a radius inside that. A
// spin below nought is this hole seen from below: turn z about instead,
// since ISCO and Circular take the spin's sense to be the one above.
func (h Hole) Valid() bool {
	// a spin from nought to under the mass puts the mass above nought
	return h.A >= 0 && h.A < h.M && h.M <= math.MaxFloat64
}

// R is the kerr-schild radius at a point: the schwarzschild r at nought
// spin, and the oblate coordinate that reaches the ring singularity
// with spin.
//
// On the disc the ring spans it is nought, and within a hair of that
// disc it rounds to nought as well, being there the difference of two
// numbers that agree to all their digits. No march comes there: a body
// is stopped at the inner horizon and light at the outer.
func (h Hole) R(x Vec3) float64 {
	return h.radius(x[0], x[1], x[2])
}

// radius is R of a point given as three numbers, which is how the march
// holds it.
func (h Hole) radius(x, y, z float64) float64 {
	positionSquaredSum := x*x + y*y + z*z
	spinSquared := h.A * h.A
	discriminant := positionSquaredSum - spinSquared
	return math.Sqrt(0.5 * (discriminant +
		math.Sqrt(discriminant*discriminant+4*spinSquared*z*z)))
}

// kernel is the kerr-schild scalar and the null covector at a
// point: the metric is flat plus scalar * vector * vector.
//
// On the disc inside the ring the radius is nought and the kernel has
// no value; it answers flat space there, which nothing should lean on.
func (h Hole) kernel(x Vec3) (scalar float64, vector Vec4) {
	radius := h.R(x)
	spin := h.A
	radiusSquared, spinSquared := radius*radius, spin*spin
	denominator := radiusSquared*radiusSquared + spinSquared*x[2]*x[2]
	if denominator == 0 || radius == 0 {
		return 0, Vec4{1, 0, 0, 0}
	}
	scalar = 2 * h.M * radius * radiusSquared / denominator
	radiusSpinSum := radiusSquared + spinSquared
	vector = Vec4{
		1,
		(radius*x[0] + spin*x[1]) / radiusSpinSum,
		(radius*x[1] - spin*x[0]) / radiusSpinSum,
		x[2] / radius,
	}
	return scalar, vector
}

// Metric is the covariant metric at a point: flat plus scalar * vector *
// vector.
func (h Hole) Metric(x Vec3) Mat4 {
	scalar, vector := h.kernel(x)
	var metric Mat4
	for i := range 4 {
		for j := range 4 {
			metric[i][j] = scalar * vector[i] * vector[j]
		}
	}
	metric[0][0] -= 1
	metric[1][1]++
	metric[2][2]++
	metric[3][3]++
	return metric
}

// Inverse is the contravariant metric at a point: flat minus scalar *
// vector * vector with vector raised, which the kerr-schild form gives
// exactly.
func (h Hole) Inverse(x Vec3) Mat4 {
	scalar, vector := h.kernel(x)
	raisedVector := Vec4{-vector[0], vector[1], vector[2], vector[3]}
	var metric Mat4
	for i := range 4 {
		for j := range 4 {
			metric[i][j] = -scalar * raisedVector[i] *
				raisedVector[j]
		}
	}
	metric[0][0] -= 1
	metric[1][1]++
	metric[2][2]++
	metric[3][3]++
	return metric
}

// Lower is a contravariant vector made covariant with the metric at x.
func (h Hole) Lower(x Vec3, vector Vec4) Vec4 {
	scalar, null := h.kernel(x)
	return lowered(scalar, null, vector)
}

// Raise is a covariant vector made contravariant with the inverse at x.
func (h Hole) Raise(x Vec3, covector Vec4) Vec4 {
	scalar, null := h.kernel(x)
	return raised(scalar, null, covector)
}

// lowered is a vector made covariant by a kernel already in hand: the
// flat metric on it, plus scalar times the null covector times the null
// covector on it. The metric is a flat one and a rank of one more, and
// this uses it without building its sixteen entries.
//
// The null covector's time part is one, always, and is not multiplied
// by here or in raised.
func lowered(scalar float64, null, vector Vec4) Vec4 {
	pull := scalar * (vector[0] + null[1]*vector[1] +
		null[2]*vector[2] + null[3]*vector[3])
	return Vec4{
		pull - vector[0],
		pull*null[1] + vector[1],
		pull*null[2] + vector[2],
		pull*null[3] + vector[3],
	}
}

// raised is a covector made contravariant by a kernel already in hand:
// the flat inverse on it, less scalar times the null vector times the
// null vector on it, the null vector being the null covector with its
// time part turned about.
func raised(scalar float64, null, covector Vec4) Vec4 {
	pull := scalar * (covector[1]*null[1] + covector[2]*null[2] +
		covector[3]*null[3] - covector[0])
	return Vec4{
		pull - covector[0],
		covector[1] - pull*null[1],
		covector[2] - pull*null[2],
		covector[3] - pull*null[3],
	}
}

// Dot is the inner product of a covariant and a contravariant vector.
func Dot(covector, vector Vec4) float64 {
	return covector[0]*vector[0] + covector[1]*vector[1] +
		covector[2]*vector[2] + covector[3]*vector[3]
}

// Horizon is the outer horizon's radius: 2M at nought spin.
func (h Hole) Horizon() float64 {
	return h.M + h.horizonSplit()
}

// innerHorizon is the inner horizon's radius: nought at nought spin,
// where it is the singularity, and up to M at full spin.
func (h Hole) innerHorizon() float64 {
	return h.M - h.horizonSplit()
}

// horizonSplit is how far the two horizons lie either side of M.
func (h Hole) horizonSplit() float64 {
	return math.Sqrt(math.Max(0, h.M*h.M-h.A*h.A))
}

// PhotonSphere is the radius of the circular photon orbit of a hole
// with no spin, 3M; with spin the orbits split prograde and retrograde
// and this is the schwarzschild value. PhotonOrbit is the two.
func (h Hole) PhotonSphere() float64 { return 3 * h.M }

// PhotonOrbit is the radius of the circular photon orbit in the equator,
// prograde or retrograde: 3M at nought spin, down to M prograde and up
// to 4M retrograde at full spin, by bardeen, press and teukolsky.
//
// Every other photon orbit, the ones that leave the equator, lies between
// the two.
func (h Hole) PhotonOrbit(prograde bool) float64 {
	against := h.A / h.M
	if prograde {
		against = -against
	}
	return 2 * h.M * (1 + math.Cos(2*math.Acos(against)/3))
}

// Shadow is the radius of the shadow a hole with no spin casts on the
// sky of a far observer, root twenty-seven M: the impact parameter at
// which a ray is just captured.
func (h Hole) Shadow() float64 { return math.Sqrt(27) * h.M }

// ISCO is the innermost stable circular orbit, prograde or retrograde:
// 6M at nought spin, down to M prograde and up to 9M retrograde at
// full spin, by bardeen, press and teukolsky.
func (h Hole) ISCO(prograde bool) float64 {
	mass := h.M
	normalizedSpin := h.A / mass
	spinSquared := normalizedSpin * normalizedSpin
	term1 := 1 + math.Cbrt(1-spinSquared)*
		(math.Cbrt(1+normalizedSpin)+math.Cbrt(1-normalizedSpin))
	term2 := math.Sqrt(3*spinSquared + term1*term1)
	discriminant := math.Sqrt((3 - term1) * (3 + term1 + 2*term2))
	if prograde {
		return mass * (3 + term2 - discriminant)
	}
	return mass * (3 + term2 + discriminant)
}
