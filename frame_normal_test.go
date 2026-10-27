package holemath

import (
	"math"
	"testing"
)

// Normal is a unit timelike velocity, carries no angular momentum, and
// exists inside the ergosphere where Static does not.
func TestTheNormalObserver(t *testing.T) {
	spinning := Hole{M: 1, A: 0.9}
	inside := Vec3{1.8, 0.2, 0.1}
	if static := spinning.Static(inside); !math.IsNaN(static[0]) {
		t.Fatalf("static inside the ergosphere is %v, not nan", static)
	}
	for _, point := range []Vec3{inside, {30, 4, -2}, {0, 3, 0.5}} {
		velocity := spinning.Normal(point)
		norm := Dot(spinning.Lower(point, velocity), velocity)
		if math.Abs(norm+1) > 1e-9 {
			t.Errorf("normal at %v has norm %v", point, norm)
		}
		covariant := spinning.Lower(point, velocity)
		if angular := point[0]*covariant[2] - point[1]*covariant[1]; math.Abs(angular) > 1e-9 {
			t.Errorf("normal at %v carries angular momentum %v", point, angular)
		}
	}
	still := Hole{M: 1}
	out := Vec3{10, 0, 0}
	if velocity := still.Normal(out); !(velocity[1] < 0) || velocity[2] != 0 {
		t.Errorf("normal round a still hole at %v is %v, and should fall inward", out, velocity)
	}
}

// Null normalises its direction, so any length gives the same null
// ray; a direction of no length is left alone and does not make nan.
func TestNullTakesAnyLength(t *testing.T) {
	hole := Hole{M: 1, A: 0.6}
	position := Vec3{-30, 6, 2}
	frame := hole.FrameOf(position, hole.Static(position))
	unit := hole.Null(position, frame, Vec3{1, 0, 0})
	long := hole.Null(position, frame, Vec3{7, 0, 0})
	for i := range 4 {
		if math.Abs(unit.P[i]-long.P[i]) > 1e-12 {
			t.Errorf("momentum %d: unit %v, long %v", i, unit.P[i], long.P[i])
		}
	}
	if h := hole.H(hole.Null(position, frame, Vec3{3, 4, 12})); math.Abs(h) > 1e-10 {
		t.Errorf("H of a ray of length 13 is %v", h)
	}
	none := hole.Null(position, frame, Vec3{})
	for i := range 4 {
		if math.IsNaN(none.P[i]) {
			t.Errorf("a direction of no length gives nan in momentum %d", i)
		}
	}
}
