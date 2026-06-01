# Ballbox

(Scroll to the bottom for some demo gifs)

A lightweight, single-header 3D physics engine written in C. Because apparently that's a thing that didn't exist and I kinda need one for my other project. It does spheres, boxes, and static SDF (signed distance functions). You should probably try one of the fancy engines/libraries first, with the C bindings and whatnot. When that proves to be too much of a headache, and you feel like all you need is balls and boxes, or balls in a box, then you come back here.

Balls in a box!

![Andy Samberg Gif](https://i.gifer.com/Jm0S.gif)

**Features:**
- ✅ Dynamic spheres and cubes with realistic physics
- ✅ Static spheres and cubes. Which is also realistic, I suppose
- ✅ Compound rigid bodies with multiple sphere/box colliders *NEW*
- ✅ Impulse-based collision resolution
- ✅ Mass, velocity, and angular velocity control
- ✅ Collision callbacks
- ✅ Force and torque application
- ✅ Single header file
- ✅ Data Oriented Design (to the best of my abilities, at least)
- ✅ Static SDFs (Signed Distance Functions. I need it for my other thing) and implicit surfaces
- ✅ Constraint solver for hinges, with simple fixed and simple hinge for now *NEW*

## Quick Start

1. Include the header:

```c
#define BALLBOX_IMPLEMENTATION
#include "ballbox.h"
```

Only define `BALLBOX_IMPLEMENTATION` in one C/C++ file.

2. Create a physics world with `CreatePhysicsWorld()`.
3. Add objects with `AddBox()` and `AddSphere()`.
4. Update the simulation with `PhysicsStep()`.

The first two arguments to `CreatePhysicsWorld()` are collider capacities: max sphere colliders and max box colliders. The third argument is contact/constraint capacity.

```c
#include <stdio.h>

#define BALLBOX_IMPLEMENTATION
#include "ballbox.h"

int main(void) {
    PhysicsWorld* world = CreatePhysicsWorld(32, 32, 128);
    world->gravity = (Vec3){ 0.0f, -9.8f, 0.0f };

    int ball = AddSphere(world, (Vec3){ 0.0f, 5.0f, 0.0f }, 0.5f);
    SetBodyMass(world, ball, 1.0f);
    SetBodyInertiaAuto(world, ball);

    int ground = AddBox(world,
                        (Vec3){ 0.0f, -0.5f, 0.0f },
                        (Vec3){ 8.0f, 0.5f, 8.0f },
                        QuatIdentity());
    SetBodyStatic(world, ground, true);

    for (int i = 0; i < 240; i++) {
        PhysicsStep(world, 1.0f / 60.0f);

        Vec3 pos = world->bodies->positions[ball];
        printf("ball: %.2f %.2f %.2f\n", pos.x, pos.y, pos.z);
    }

    DestroyPhysicsWorld(world);
    return 0;
}
```

`AddSphere()` and `AddBox()` are convenience functions. They create one rigid body and one collider, and return the body index. Prefer `SetBodyMass()`, `SetBodyStatic()`, `AddBodyForce()`, `AddBodyTorque()`, `SetBodyVelocity()`, and `SetBodyAngularVelocity()` when working with that returned handle.

## Compound Bodies

For compound objects, create one rigid body and attach multiple colliders to it. Collision detection runs against the child colliders, but collision response moves the single parent body.

```c
int ship = AddRigidBody(world, (Vec3){ 0.0f, 8.0f, 0.0f }, QuatIdentity());

AddBoxCollider(world, ship,
               (Vec3){ 0.0f, 0.0f, 0.0f },
               QuatIdentity(),
               (Vec3){ 2.0f, 0.25f, 0.5f });

AddBoxCollider(world, ship,
               (Vec3){ -1.5f, 0.0f, 0.0f },
               QuatIdentity(),
               (Vec3){ 0.25f, 0.25f, 1.5f });

AddBoxCollider(world, ship,
               (Vec3){ 1.5f, 0.0f, 0.0f },
               QuatIdentity(),
               (Vec3){ 0.25f, 0.25f, 1.5f });

AddSphereCollider(world, ship, (Vec3){ 0.0f, 0.45f, 0.0f }, 0.35f);

SetBodyMass(world, ship, 50.0f);
SetBodyInertiaAuto(world, ship);
```

Colliders attached to the same body do not collide with each other.

## Forces and Torques

```c
AddBodyForce(world, ship, (Vec3){ 0.0f, 100.0f, 0.0f });
AddBodyTorque(world, ship, (Vec3){ 0.0f, 10.0f, 0.0f });
```

The old collider-specific helpers still exist:

```c
AddSphereForce(world, sphereColliderIndex, force);
AddBoxTorque(world, boxColliderIndex, torque);
```

But new code should usually use the body functions.

## Collision Callbacks

Collision callbacks are registered on colliders, not bodies.

```c
void OnSphereCollision(int selfCollider,
                       int otherColliderType,
                       int otherCollider,
                       CollisionContact* contact) {
    printf("sphere collider %d hit collider %d, penetration %.3f\n",
           selfCollider,
           otherCollider,
           contact->penetration);
}

int body = AddRigidBody(world, (Vec3){ 0.0f, 4.0f, 0.0f }, QuatIdentity());
int sphereCollider = AddSphereCollider(world, body, (Vec3){ 0.0f, 0.0f, 0.0f }, 0.5f);

SetSphereCollisionCallback(world, sphereCollider, OnSphereCollision);
```

The `CollisionContact` itself contains both body and collider IDs:

```c
contact->bodyA_index;
contact->colliderA_index;
contact->bodyB_index;
contact->colliderB_index;
```

## Joints

Joints connect rigid body indices.

### Fixed Joint

Keeps two bodies locked together at an anchor, including relative orientation.

```c
int a = AddBox(world, (Vec3){ -1, 5, 0 }, (Vec3){ 0.5f, 0.5f, 0.5f }, QuatIdentity());
int b = AddBox(world, (Vec3){  1, 5, 0 }, (Vec3){ 0.5f, 0.5f, 0.5f }, QuatIdentity());

AddFixedJoint(world, a, b, (Vec3){ 0.0f, 5.0f, 0.0f });
```

### Hinge Joint

Keeps two anchors together and allows rotation around one axis.

```c
int post = AddBox(world, (Vec3){ 0, 2, 0 }, (Vec3){ 0.1f, 2.0f, 0.1f }, QuatIdentity());
SetBodyStatic(world, post, true);

int door = AddBox(world, (Vec3){ 1.0f, 2, 0 }, (Vec3){ 1.0f, 1.2f, 0.1f }, QuatIdentity());

AddHingeJoint(world,
              post,
              door,
              (Vec3){ 0.0f, 2.0f, 0.0f },
              (Vec3){ 0.0f, 1.0f, 0.0f });
```

Optional hinge motor:

```c
int hinge = AddHingeJoint(world, bodyA, bodyB, anchor, axis);
SetHingeMotor(world, hinge, true, 2.0f, 0.25f);
```

### Spring Joint

This one is a soft force-based joint, not a hard solver constraint. It applies Hooke force plus damping between two anchors.

```c
int a = AddSphere(world, (Vec3){ -1, 6, 0 }, 0.4f);
int b = AddSphere(world, (Vec3){  1, 6, 0 }, 0.4f);

AddSpringJoint(world,
               a,
               b,
               world->bodies->positions[a],
               world->bodies->positions[b],
               1.0f,   // rest length
               35.0f,  // stiffness
               5.0f);  // damping
```

## Static SDF Collider

Signed distance functions can be defined to serve as static colliders for dynamic objects. It is a niche feature, for sure, but it's here if you need it.

```c
#include <math.h>

float WavyGroundSDF(Vec3 point) {
    float wavy_y = sinf(point.x) * sinf(point.z);
    return point.y - wavy_y;
}

PhysicsWorld* world = CreatePhysicsWorld(50, 50, 500);
UseSDF(world, WavyGroundSDF);
```

If you can provide a normal, do it. It avoids finite-difference normal estimation, which costs six SDF evaluations per contact.

```c
Vec3 WavyGroundNormal(Vec3 point) {
    float dx = -cosf(point.x) * sinf(point.z);
    float dz = -sinf(point.x) * cosf(point.z);
    return Vec3Normalize((Vec3){ dx, 1.0f, dz });
}

UseSDFNormal(world, WavyGroundNormal);
```

Box-SDF collision can run in two modes:

```c
SetSDFBoxCollisionMode(world, PHYSICS_SDF_BOX_SIGNED_DISTANCE);
```

That is faster and expects a real signed-distance field with negative interior values.

```c
SetSDFBoxCollisionMode(world, PHYSICS_SDF_BOX_SURFACE_TRACE);
```

That is the default. It traces from the box center toward boundary samples, so it can catch implicit or unsigned surfaces that only approach zero at the surface.

You can tune it:

```c
SetSDFSurfaceTraceSettings(world,
                           PHYSICS_SDF_BOX_TRACE_CORNERS, // or CORNERS_AND_FACES, or FULL
                           6,     // trace steps
                           4,     // refine steps
                           0.05f, // surface epsilon
                           4);    // max contacts per box
```

## Settings

You can modify the simulation through `world->settings`:

```c
world->settings.linear_damping = 0.99f;
world->settings.angular_damping = 0.98f;
world->settings.restitution = 0.0f;
world->settings.friction = 0.1f;
world->settings.solver_iterations = 28;
world->settings.baumgarte_bias = 0.25f;
world->settings.allowed_penetration = 0.02f;
```

There are also helper functions:

```c
SetDamping(world, 0.99f, 0.98f);
SetRestitution(world, 0.0f);
SetCorrectionParameters(world, 0.25f, 0.02f);
SetSolverIterations(world, 28);
SetManifoldSettings(world, 0.01f, -0.5f, 8);
```

The sleep system exists, but I implemented it while trying to fix jitter and then neglected it. I suggest not relying on it too much yet.

## Development Notes

During development, edit:

```text
include/physics/physics_core.h
include/physics/physics_core.c
```

Then regenerate the single-header distribution:

```powershell
python tools\merge_to_ballbox.py
```

If using `ballbox.h` in this repo's test app, make sure `build/premake5.lua` does not compile `include/**.c`, otherwise you will get duplicate symbol linker errors.

# Demo GIFs


![demogif](./physicswip2.gif)
![demogif](./physicswip3.gif)
![demogif](./physicswip_constraints_velcontrol.gif)
The previous gif is testing how the simulation handles direct velocity control and if I can constraint objects to remain together with just forces.

# To do

More demos and more gifs of new features.

# License

MIT License.

# Contributing

Found a bug? Have a feature request? Open an issue or submit a PR. 

Want to help me pay some bills? check out [my Patreon blog](https://www.patreon.com/Diegomakesgames), where I talk about game dev and other stuff.

