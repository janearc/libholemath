package holemath

import (
	"errors"
	"math"
	"testing"
	"unsafe"
)

// GPU layer is properly initialized.
func TestGPUBlackHoleCreation(t *testing.T) {
	mass, spin := 1.0, 0.6
	gpu := NewGPUBlackHole(mass, spin)

	if gpu.hole.M != mass || gpu.hole.A != spin {
		t.Errorf("GPU wrapper not initialized correctly: M=%v (want %v), A=%v (want %v)",
			gpu.hole.M, mass, gpu.hole.A, spin)
	}
}

// A camera can be held still outside the ergosurface and nowhere inside
// it. SetCamera refuses a point inside, says why, and leaves the camera
// it had.
func TestGPUSetCameraWhereNothingIsStill(t *testing.T) {
	gpu := NewGPUBlackHole(1, 0.9)
	outside := Vec3{20, 0, 0}
	if err := gpu.SetCamera(outside); err != nil {
		t.Fatalf("a camera at %v was refused: %v", outside, err)
	}
	before := *gpu
	inside := Vec3{1.8, 0, 0} // radius 1.56: outside the horizon at 1.44, inside the ergosurface at 2
	err := gpu.SetCamera(inside)
	if !errors.Is(err, ErrNothingStill) {
		t.Fatalf("a camera inside the ergosurface gave %v, want ErrNothingStill", err)
	}
	if *gpu != before {
		t.Errorf("a refused camera was set anyway: %v", gpu.cameraPos)
	}
	ray := gpu.hole.Arriving(outside, gpu.cameraFrame, Vec3{-1, 0.3, 0})
	if step := gpu.StepRay(&ray, 0.1); math.IsNaN(step.Redshift) {
		t.Errorf("after a refused camera the redshift is not a number")
	}
}

// Before a camera is set there is no camera to shift against, and the
// redshift is not a number: not a nought that reads as a shift.
func TestGPURedshiftBeforeSetCamera(t *testing.T) {
	gpu := NewGPUBlackHole(1, 0)
	position := Vec3{20, 0, 0}
	ray := gpu.hole.Null(position, gpu.hole.FrameOf(position, gpu.hole.Static(position)), Vec3{-1, 0.3, 0})
	if step := gpu.StepRay(&ray, 0.1); !math.IsNaN(step.Redshift) {
		t.Errorf("with no camera set the redshift is %v", step.Redshift)
	}
}

// The redshift StepRay reports is the still camera's energy over a
// still emitter's where the ray is: the ratio of the two clocks' rates,
// the root of one less two M over the radius at each.
func TestGPURedshiftBetweenStillObservers(t *testing.T) {
	gpu := NewGPUBlackHole(1, 0)
	position := Vec3{20, 0, 0}
	if err := gpu.SetCamera(position); err != nil {
		t.Fatal(err)
	}
	ray := gpu.hole.Arriving(position, gpu.cameraFrame, Vec3{-1, 0.2, 0})
	var step RayStep
	for range 40 {
		step = gpu.StepRay(&ray, 0.2)
	}
	want := math.Sqrt(1-2/step.Radius) / math.Sqrt(1-2/20.0)
	if !(math.Abs(step.Redshift-want) <= 1e-9) {
		t.Errorf("from radius %v to a camera at 20 the redshift is %v, want %v", step.Radius, step.Redshift, want)
	}
	if !(step.Radius <= 14) {
		t.Errorf("the ray only got to radius %v, which tests little", step.Radius)
	}
}

// The shader's struct is laid out as ShaderParams is: this is the
// layout, and a change here is a change to make there too.
func TestShaderParamsLayout(t *testing.T) {
	var params ShaderParams
	if size := unsafe.Sizeof(params); size != 96 {
		t.Errorf("ShaderParams is %d bytes, and the shader's struct is 96", size)
	}
	offsets := map[string][2]uintptr{
		"BlackHoleMass": {unsafe.Offsetof(params.BlackHoleMass), 0},
		"BlackHoleSpin": {unsafe.Offsetof(params.BlackHoleSpin), 4},
		"Horizon":       {unsafe.Offsetof(params.Horizon), 8},
		"PhotonSphere":  {unsafe.Offsetof(params.PhotonSphere), 12},
		"Shadow":        {unsafe.Offsetof(params.Shadow), 16},
		"CameraPos":     {unsafe.Offsetof(params.CameraPos), 20},
		"CameraFrame":   {unsafe.Offsetof(params.CameraFrame), 32},
	}
	for field, offset := range offsets {
		if offset[0] != offset[1] {
			t.Errorf("%s is at byte %d, and the shader's struct has it at %d", field, offset[0], offset[1])
		}
	}
}

// The horizon is the hole's from the start, before any camera is set:
// the shader is told two M.
func TestGPUHorizonBeforeSetCamera(t *testing.T) {
	gpu := NewGPUBlackHole(1, 0)
	if horizon := gpu.ShaderParams().Horizon; horizon != 2 {
		t.Errorf("before a camera is set, the shader's horizon is %v", horizon)
	}
}

// StepRay steps the ray it is handed, and the step it reports is where
// that ray now is; a batch stepped twice has gone twice as far.
func TestGPUStepRayAdvancesTheRay(t *testing.T) {
	gpu := NewGPUBlackHole(1, 0)
	position := Vec3{20, 0, 0}
	gpu.SetCamera(position)
	ray := gpu.hole.Arriving(position, gpu.cameraFrame, Vec3{-1, 0.1, 0})
	before := ray
	step := gpu.StepRay(&ray, 0.4)
	if ray == before {
		t.Errorf("StepRay left the ray where it was")
	}
	if ray.X != step.Position || ray.P != step.Momentum {
		t.Errorf("StepRay reports %v, %v and left the ray at %v, %v", step.Position, step.Momentum, ray.X, ray.P)
	}
	if want := gpu.hole.Step(before, 0.4); ray != want {
		t.Errorf("StepRay took the ray to %v; one step of the hole's takes it to %v", ray, want)
	}
	rays := []Ray{before, before}
	first := gpu.StepRaysBatch(nil, rays, 0.4)
	firstPosition := first[0].Position
	second := gpu.StepRaysBatch(first, rays, 0.4)
	if firstPosition == second[0].Position {
		t.Errorf("a batch stepped twice reports the same place twice: %v", firstPosition)
	}
	if want := gpu.hole.Step(gpu.hole.Step(before, 0.4), 0.4); rays[1] != want {
		t.Errorf("a batch stepped twice left its ray at %v; two steps take it to %v", rays[1], want)
	}
}

// SetCamera keeps the camera's place and works out its frame.
func TestGPUSetCamera(t *testing.T) {
	gpu := NewGPUBlackHole(1, 0.5)
	position := Vec3{20, 0, 0}

	gpu.SetCamera(position)

	// Camera position should be stored
	if gpu.cameraPos != position {
		t.Errorf("camera position not stored correctly")
	}

	// Frame should be computed
	for i := range 4 {
		for j := range 4 {
			if math.IsNaN(float64(gpu.cameraFrame[i][j])) {
				t.Errorf("camera frame contains NaN")
			}
		}
	}
}

// StepRay returns correctly formatted output.
func TestGPUStepRay(t *testing.T) {
	gpu := NewGPUBlackHole(1, 0.6)
	position := Vec3{20, 0, 0}
	gpu.SetCamera(position)

	frame := gpu.cameraFrame
	ray := gpu.hole.Null(position, frame, Vec3{1, 0, 0})

	step := gpu.StepRay(&ray, 0.1)

	// All outputs should be finite
	if math.IsNaN(step.Redshift) || math.IsInf(step.Redshift, 0) {
		t.Errorf("redshift not finite: %v", step.Redshift)
	}
	if math.IsNaN(step.Energy) || math.IsInf(step.Energy, 0) {
		t.Errorf("energy not finite: %v", step.Energy)
	}
	if math.IsNaN(step.Hamiltonian) || math.IsInf(step.Hamiltonian, 0) {
		t.Errorf("hamiltonian not finite: %v", step.Hamiltonian)
	}
	if math.IsNaN(step.Radius) || step.Radius <= 0 {
		t.Errorf("radius not valid: %v", step.Radius)
	}
}

// ShaderParams exports correct float32 values.
func TestGPUShaderParams(t *testing.T) {
	gpu := NewGPUBlackHole(1, 0.5)
	position := Vec3{20, 0, 0}
	gpu.SetCamera(position)

	params := gpu.ShaderParams()

	// Values should be finite and match hole
	if params.BlackHoleMass != 1.0 || params.BlackHoleSpin != 0.5 {
		t.Errorf("shader params not exported correctly")
	}

	if params.Horizon <= 0 || math.IsNaN(float64(params.Horizon)) {
		t.Errorf("shader params horizon invalid: %v", params.Horizon)
	}

	// Camera frame should be populated
	anyNonZero := false
	for i := range 4 {
		for j := range 4 {
			if params.CameraFrame[i][j] != 0 {
				anyNonZero = true
			}
		}
	}
	if !anyNonZero {
		t.Errorf("shader params camera frame is all zeros")
	}
}

// Multiple SetCamera calls update state correctly.
func TestGPUMultipleCameras(t *testing.T) {
	gpu := NewGPUBlackHole(1, 0.6)

	pos1 := Vec3{20, 0, 0}
	gpu.SetCamera(pos1)
	frame1 := gpu.cameraFrame

	pos2 := Vec3{30, 0, 0}
	gpu.SetCamera(pos2)
	frame2 := gpu.cameraFrame

	// Frames should be different
	same := true
	for i := range 4 {
		for j := range 4 {
			if frame1[i][j] != frame2[i][j] {
				same = false
				break
			}
		}
	}
	if same {
		t.Errorf("different camera positions should produce different frames")
	}
}

// A batch reports on every ray in order, into the slice it is handed
// when that is long enough, so stepping a batch again and again
// allocates nothing; a slice too short is replaced, and no rays is no
// reports.
func TestGPUStepRaysBatch(t *testing.T) {
	gpu := NewGPUBlackHole(1, 0.5)
	position := Vec3{20, 0, 0}
	if err := gpu.SetCamera(position); err != nil {
		t.Fatal(err)
	}
	rays := make([]Ray, 5)
	for i := range rays {
		rays[i] = gpu.hole.Arriving(position, gpu.cameraFrame, Vec3{-1, 0.1 * float64(i), 0})
	}
	want := make([]Ray, len(rays))
	for i := range rays {
		want[i] = gpu.hole.Step(rays[i], 0.1)
	}
	reports := gpu.StepRaysBatch(nil, rays, 0.1)
	if len(reports) != len(rays) {
		t.Fatalf("%d reports for %d rays", len(reports), len(rays))
	}
	for i := range rays {
		if rays[i] != want[i] || reports[i].Position != want[i].X || reports[i].Radius != gpu.hole.R(want[i].Space()) {
			t.Errorf("ray %d: stepped to %v and reported at %v; one step takes it to %v", i, rays[i], reports[i].Position, want[i])
		}
	}
	reused := gpu.StepRaysBatch(reports, rays, 0.1)
	if &reused[0] != &reports[0] {
		t.Errorf("a slice long enough was not reused")
	}
	if allocations := testing.AllocsPerRun(100, func() { reports = gpu.StepRaysBatch(reports, rays, 0.1) }); allocations != 0 {
		t.Errorf("stepping a batch into its own reports allocates %v times", allocations)
	}
	if short := gpu.StepRaysBatch(make([]RayStep, 0, 2), rays, 0.1); len(short) != len(rays) {
		t.Errorf("a slice too short gave %d reports for %d rays", len(short), len(rays))
	}
	if none := gpu.StepRaysBatch(reports, nil, 0.1); len(none) != 0 {
		t.Errorf("no rays gave %d reports", len(none))
	}
}

// Hole hands back the hole the wrapper was made of, a copy: whatever
// is done to it, the wrapper's own is as it was.
func TestGPUHole(t *testing.T) {
	gpu := NewGPUBlackHole(1, 0.6)
	hole := gpu.Hole()
	if hole != (Hole{M: 1, A: 0.6}) {
		t.Errorf("the wrapper's hole is %v", hole)
	}
	hole.A = 0
	if gpu.Hole().A != 0.6 {
		t.Errorf("changing the copy changed the wrapper's hole to %v", gpu.Hole())
	}
}
