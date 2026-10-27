# to do

## simple lens planes

Lensing far from any one hole, where a thin plane of mass bends light
as it passes. Galaxies and clusters bend background light into arcs and
Einstein rings.

In microlensing a passing star briefly brightens one behind it.

The usual models are a point mass and a singular isothermal sphere,
plus a sheet of mass for the rest. Each one maps where a ray crosses
the plane to how far it's bent. slayer would use them the way it uses
the Kerr hole today.

## the full Kerr march

slayer marches its own rays through a hole today, in
render/trace/field.go, using holemath's geodesics. The march itself
belongs here.

Given where a ray starts and which way it goes, the march follows the
geodesic until the ray escapes the reach or crosses the horizon, and
says which happened.

When it escapes it says where and how far the light has shifted. slayer
then calls it instead of stepping the geodesic itself.
