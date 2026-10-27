# definitions

## definitions

every type and function in the package, one line each.

types:

```
Vec3               Holds a point or a direction in space.
Vec4               Holds a spacetime point, or a vector there.
Mat4               Holds a metric or its inverse.
Hole               Holds a black hole's mass M and spin A.
Ray                Holds a geodesic's position X and momentum P.
Fate               Names how a march ended.
Frame              Holds an observer's orthonormal tetrad.
RayClassification  Names a ray's kind: escapes, plunges, circulates.
GPUBlackHole       Holds a hole and a still camera, for a shader.
RayStep            Holds a stepped ray and its kept quantities.
ShaderParams       Holds the hole and camera as float32 constants.
```

functions and values:

```
Dot              Contracts a covector with a vector.
Shift            Gives the photon energy ratio, eye over emitter.
NewGPUBlackHole  Makes a GPUBlackHole from a mass and a spin.
ErrNothingStill  Says no observer can be held still at a point.
```

methods of Hole: the metric:

```
Valid         Says whether the mass and spin are a black hole's.
R             Returns the kerr-schild radius at a point.
Metric        Returns the covariant metric at a point.
Inverse       Returns the contravariant metric at a point.
Lower         Makes a vector covariant at a point.
Raise         Makes a covector contravariant at a point.
Horizon       Returns the outer horizon radius.
PhotonSphere  Returns 3M, the photon orbit at nought spin.
PhotonOrbit   Returns the equator's photon orbit, either sense.
Shadow        Returns root 27 M, the shadow at nought spin.
ISCO          Returns the innermost stable circular orbit.
```

methods of Hole and Ray: geodesics:

```
Space            Returns the ray's point in space.
H                Returns hamilton's function for a ray.
E                Returns the energy a ray keeps.
Lz               Returns the angular momentum about the z axis.
Carter           Returns the carter constant of a ray.
Velocity         Returns dx over the affine parameter.
Step             Advances a ray by one runge-kutta step.
Pace             Returns a fiftieth of a radius, a step for a body.
PaceIn           Returns the step for a ray where it stands.
Hovering         Says whether an arriving ray came from the edge.
March            Steps a ray until it escapes, is captured or runs out.
ImpactParameter  Returns the whole angular momentum over E.
ClosestApproach  Returns the least radius on a ray's march.
DeflectionAngle  Returns the angle a ray is bent by, in any plane.
ClassifyRay      Says whether a ray escapes, plunges or circulates.
```

methods of Hole: observers:

```
Static           Returns the velocity of an observer held still.
Circular         Returns the velocity of matter circling prograde.
FrameOf          Returns an observer's tetrad at a point.
Null             Returns the light ray leaving in a frame direction.
Arriving         Returns an arriving ray, to be marched back.
Seen             Returns the direction a ray's light comes from.
Falling          Returns a body dropped at rest on the x axis.
Dropped          Returns a body dropped at rest at a point.
Advance          Carries a body on by a proper time.
ToObserverFrame  Returns a covector's components on a frame.
```

methods of GPUBlackHole:

```
SetCamera        Sets a still camera, or says nothing is still there.
StepRay          Steps a ray in place and reports what it keeps.
StepRaysBatch    Steps each ray in a slice, in place, in order.
Hole             Returns the hole itself.
ShaderParams     Returns the hole and camera as float32 constants.
```

unexported helpers:

```
radius             Returns R of a point given as three numbers.
kernel             Returns the kerr-schild scalar and null covector.
lowered            Lowers a vector with a kernel already in hand.
raised             Raises a covector with a kernel already in hand.
horizonSplit       Returns how far each horizon lies from M.
innerHorizon       Returns the inner horizon radius.
limits             Holds the radii a march decides a ray's fate by.
limitsOf           Works those radii out for a ray.
noReturn           Returns the photon orbit with a ray's own Lz over E.
goingInward        Says whether a moving point is losing radius.
space              Returns the space part of a spacetime vector.
acrossAxisSquared  Returns the squared angular momentum across z.
quad               Holds four numbers where registers can keep them.
quadOf             Makes a quad of a Vec4.
vec4               Makes a Vec4 of a quad.
moved              Carries a quad along a rate by a step.
rungeKutta         Weighs a step's four rates together.
derive             Returns hamilton's equations for a ray, exactly.
rateAt             Does the same for a point of known radius.
rateOf             Does the same for a point of unknown radius.
stepFrom           Takes a step whose first rate is in hand.
timeVelocity       Returns dt over the affine parameter.
runsAway           Says whether a ray's time runs away backward.
pacedStep          Takes one step of a march, or gives the ray up.
angleBetween       Returns the angle between two directions.
orbiting           Returns the velocity of an observer going round z.
sub                Subtracts a scaled vector from a vector.
ringDistance       Returns the distance to the ring singularity.
classOf            Maps a march's fate to a ray's kind.
```
