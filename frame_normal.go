package holemath

import "math"

// Normal is the four-velocity of the observer normal to the chart's
// slices of constant time at a point: it carries no angular momentum, is
// carried round by the dragging and nothing else, and exists everywhere
// outside the horizon, ergosphere included, where Static does not.
//
// It is not still: in this chart it falls inward, slowly far out and
// faster close in, even round a still hole. Where Static exists, it is
// the one to prefer for a thing held in place.
func (h Hole) Normal(x Vec3) Vec4 {
	scalar, null := h.kernel(x)
	scale := 1 / math.Sqrt(1+scalar)
	return Vec4{
		scale * (1 + scalar),
		-scale * scalar * null[1],
		-scale * scalar * null[2],
		-scale * scalar * null[3],
	}
}
