package holemath

import "testing"

// benchRay is an arriving ray off the equator of a spinning hole, where
// no term of the kernel or its gradient vanishes.
func benchRay() (Hole, Ray) {
	hole := Hole{M: 1, A: 0.6}
	position := Vec3{12, 5, 3}
	frame := hole.FrameOf(position, hole.Static(position))
	return hole, hole.Arriving(position, frame, Vec3{-0.9, 0.3, 0.2})
}

// The radius: two square roots, taken at every evaluation of the kernel.
func BenchmarkR(b *testing.B) {
	hole, ray := benchRay()
	position := ray.Space()
	for b.Loop() {
		hole.R(position)
	}
}

// One evaluation of hamilton's equations, four of which make a step.
func BenchmarkDerive(b *testing.B) {
	hole, ray := benchRay()
	for b.Loop() {
		hole.derive(ray)
	}
}

// One runge-kutta step, the unit a renderer pays for.
func BenchmarkStep(b *testing.B) {
	hole, ray := benchRay()
	for b.Loop() {
		hole.Step(ray, 0.1)
	}
}

// A whole march of one pixel's ray, from the eye to its fate.
func BenchmarkMarch(b *testing.B) {
	hole, ray := benchRay()
	for b.Loop() {
		hole.March(ray, 1000, 20000)
	}
}

// A whole march of a pixel inside the shadow: a ray that is captured.
func BenchmarkMarchCaptured(b *testing.B) {
	hole, ray := benchRay()
	position := ray.Space()
	frame := hole.FrameOf(position, hole.Static(position))
	captured := hole.Arriving(position, frame, Vec3{-0.9, -0.38, -0.22})
	if _, fate := hole.March(captured, 1000, 20000); fate != Captured {
		b.Fatalf("the ray this benchmark is for had fate %v", fate)
	}
	for b.Loop() {
		hole.March(captured, 1000, 20000)
	}
}

// Hamilton's function, which every check of a ray asks for.
func BenchmarkH(b *testing.B) {
	hole, ray := benchRay()
	for b.Loop() {
		hole.H(ray)
	}
}

// A covector raised, which is what Velocity is.
func BenchmarkRaise(b *testing.B) {
	hole, ray := benchRay()
	position := ray.Space()
	for b.Loop() {
		hole.Raise(position, ray.P)
	}
}

// An observer's frame: once for a camera, and once for every emitter
// whose own directions are wanted.
func BenchmarkFrameOf(b *testing.B) {
	hole, ray := benchRay()
	position := ray.Space()
	observer := hole.Circular(position)
	for b.Loop() {
		hole.FrameOf(position, observer)
	}
}

// The wrapper's step with everything it reports.
func BenchmarkGPUStepRay(b *testing.B) {
	hole, ray := benchRay()
	gpu := NewGPUBlackHole(hole.M, hole.A)
	gpu.SetCamera(ray.Space())
	for b.Loop() {
		stepped := ray
		gpu.StepRay(&stepped, 0.1)
	}
}
