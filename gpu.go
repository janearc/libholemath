package holemath

import (
	"errors"
	"fmt"
	"math"
)

// GPUBlackHole is a hole and a camera held still beside it, as a shader
// wants them: ShaderParams is the two as float32 constants, and StepRay
// steps a ray on the processor and reports what it keeps, for holding a
// shader's own steps against.
//
// There is no GPU in it. It steps one ray at a time on one core, and
// its batches are loops. SetCamera changes it, so it is not to be
// shared between goroutines while one of them sets the camera.
type GPUBlackHole struct {
	hole      Hole
	cameraPos Vec3
	// the still observer's at cameraPos; not a number until a camera is set
	cameraVelocity Vec4
	cameraFrame    Frame
}

// RayStep is a ray after a step and what it keeps there: a geodesic
// holds Energy, AngularMomentum, CarterConstant and Hamiltonian fixed,
// so their drift from step to step is the stepping's error.
type RayStep struct {
	// Spacetime position and covariant momentum
	Position Vec4
	Momentum Vec4
	// Redshift is the photon's energy at the camera over its energy as
	// an emitter held still where the ray now is would have it. It is
	// not a number before a camera is set, and where nothing is still.
	Redshift        float64
	Energy          float64
	AngularMomentum float64
	CarterConstant  float64
	Hamiltonian     float64
	// Trajectory info
	Radius float64
}

// ShaderParams is the hole and the camera as constants for a shader, in
// float32, which is what a shader has.
//
// Its layout is a contract with the shader's own struct: five floats, a
// packed three, and a four by four at byte 32, 96 bytes in all. A field
// moved here is a field moved there, and a test holds the offsets.
type ShaderParams struct {
	BlackHoleMass float32
	BlackHoleSpin float32
	Horizon       float32
	PhotonSphere  float32
	Shadow        float32
	CameraPos     [3]float32
	CameraFrame   [4][4]float32
}

// ErrNothingStill is SetCamera's refusal of a point where no observer
// can be held still: inside the ergosurface of a spinning hole, or
// inside the horizon of any.
var ErrNothingStill = errors.New("holemath: nothing can be held still there")

// NewGPUBlackHole is a hole of a mass and a spin with no camera yet.
// SetCamera gives it one.
func NewGPUBlackHole(mass, spin float64) *GPUBlackHole {
	notANumber := math.NaN()
	return &GPUBlackHole{
		hole: Hole{M: mass, A: spin},
		cameraVelocity: Vec4{
			notANumber, notANumber, notANumber, notANumber,
		},
	}
}

// SetCamera puts a camera held still at a point and works out its
// frame, once, for every step and every ShaderParams after it.
//
// Where nothing is still there is no such camera: it answers
// ErrNothingStill and leaves the camera where it was.
func (gb *GPUBlackHole) SetCamera(position Vec3) error {
	velocity := gb.hole.Static(position)
	if math.IsNaN(velocity[0]) {
		return fmt.Errorf(
			"%w: %v, at radius %v of a hole of mass %v and spin %v",
			ErrNothingStill, position, gb.hole.R(position),
			gb.hole.M, gb.hole.A,
		)
	}
	gb.cameraPos, gb.cameraVelocity = position, velocity
	gb.cameraFrame = gb.hole.FrameOf(position, velocity)
	return nil
}

// StepRay takes a ray one runge-kutta step of a size, in place, and
// reports where it now is and what it keeps there.
//
// The report costs more than twice what the step does, so a march that
// does not read it is better made with Hole's March or Step.
func (gb *GPUBlackHole) StepRay(ray *Ray, stepSize float64) RayStep {
	stepped := gb.hole.Step(*ray, stepSize)
	*ray = stepped
	return RayStep{
		Position: stepped.X,
		Momentum: stepped.P,
		Redshift: Shift(
			stepped.P, gb.cameraVelocity,
			stepped.P, gb.hole.Static(stepped.Space()),
		),
		Energy:          gb.hole.E(stepped),
		AngularMomentum: gb.hole.Lz(stepped),
		CarterConstant:  gb.hole.Carter(stepped),
		Hamiltonian:     gb.hole.H(stepped),
		Radius:          gb.hole.R(stepped.Space()),
	}
}

// StepRaysBatch takes each ray of a slice one step of a size, in place,
// and reports on each in the same order.
//
// The reports go into the slice handed in when it is long enough, so a
// caller that keeps handing back what this returns allocates once; nil
// will do for the first call.
func (gb *GPUBlackHole) StepRaysBatch(
	reports []RayStep, rays []Ray, stepSize float64,
) []RayStep {
	if cap(reports) < len(rays) {
		reports = make([]RayStep, len(rays))
	}
	reports = reports[:len(rays)]
	for i := range rays {
		reports[i] = gb.StepRay(&rays[i], stepSize)
	}
	return reports
}

// Hole is the hole itself, for everything a Hole can be asked: its
// metric, its marches, what its rays keep.
func (gb *GPUBlackHole) Hole() Hole {
	return gb.hole
}

// ShaderParams is the hole and the camera as constants for a shader.
// Before a camera is set the camera's part is all nought.
func (gb *GPUBlackHole) ShaderParams() ShaderParams {
	params := ShaderParams{
		BlackHoleMass: float32(gb.hole.M),
		BlackHoleSpin: float32(gb.hole.A),
		Horizon:       float32(gb.hole.Horizon()),
		PhotonSphere:  float32(gb.hole.PhotonSphere()),
		Shadow:        float32(gb.hole.Shadow()),
	}
	for i := range 3 {
		params.CameraPos[i] = float32(gb.cameraPos[i])
	}
	for i := range 4 {
		for j := range 4 {
			params.CameraFrame[i][j] = float32(gb.cameraFrame[i][j])
		}
	}
	return params
}
