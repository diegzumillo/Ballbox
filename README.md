# Ballbox

(Scroll to the bottom for some demo gifs)

> **⚠️ WARNING**  
> I released this thing just now. It worked in my test scene on my computer. Still needs testing.


A lightweight, single-header 3D physics engine written in C. Because apparently that's a thing that didn't exist and I kinda need one for my other project. It does spheres, boxes, and static SDF (signed distance functions). You should probably try one of the fancy engines/libraries first, with the C bindings and whatnot. When that proves to be too much of a headache, and you feel like all you need is balls and boxes, or balls in a box, then you come back here.

Balls in a box!

![Andy Samberg Gif](https://i.gifer.com/Jm0S.gif)

**Features:**
- ✅ Dynamic spheres and cubes with realistic physics
- ✅ Static spheres and cubes. Which is also realistic, I suppose
- ✅ Impulse-based collision resolution
- ✅ Mass, velocity, and angular velocity control
- ✅ Collision callbacks
- ✅ Force and torque application
- ✅ Single header file
- ✅ Data Oriented Design (to the best of my abilities, at least)
- ✅ Static SDFs (Signed Distance Functions. I need it for my other thing)
- ⏳ Constraints, like hinges, springs. So we can use boxes and spheres to build more complex objects.

## Quick Start

1. Include the header: `#include "ballbox.h"`
2. You also need to add `#define BALLBOX_IMPLEMENTATION` in your source file, or before the include
3. Create a physics world with `CreatePhysicsWorld()`
4. Add objects with `AddBox()` and `AddSphere()`
5. Update the simulation with `PhysicsStep()`

As for the arguments, I asked a generic AI to look at my test scene that uses raylib, "hey, make a shorter program that's easier to follow that illustrates ballbox". And it spat this out. It looks solid to me!

```c

  #include <stdio.h>
  #define BALLBOX_IMPLEMENTATION
  #include "ballbox.h"

  int main() {
      // Create a physics world
      PhysicsWorld* world = CreatePhysicsWorld(10, 10, 50);
      world->gravity = (Vec3){0, -9.8f, 0};  // Earth gravity

      // Add a sphere that will fall
      int sphere = AddSphere(world, (Vec3){0, 5, 0}, 1.0f);

      // Add a static box as ground
      int ground = AddBox(world, (Vec3){0, -1, 0}, (Vec3){10, 1, 10}, QuatIdentity());
      SetBoxStatic(world, ground, true);


      // Collision callback test function
      void OnSphereHitsBox(int selfIndex, int otherBodyType, int otherIndex, CollisionContact* contact) {
          if (otherBodyType == 1) { // 1 = box
              printf("Sphere %d hit Box %d! Penetration: %.3f\n", selfIndex, otherIndex, contact->penetration);
          }
      }
    
      // Register collision callback for the sphere
      SetSphereCollisionCallback(world, sphere, OnSphereHitsBox);

      // Simulate for 1000 steps. (Or in a while loop, like any sane game developer. This AI is mad.)
      for (int i = 0; i < 1000; i++) {
          // Apply a sideways force for the first 30 steps
          if (i < 30) {
              AddSphereForce(world, sphere, (Vec3){2.0f, 0, 0});  // Push right
          }

          // You should probably use some kind of deltatime your engine provides here
          PhysicsStep(world, 1.0f/60.0f);  

          // Print sphere position every 10 steps
          if (i % 10 == 0) {
              Vec3 pos = world->spheres->positions[sphere];
              printf("Step %d: Sphere at (%.2f, %.2f, %.2f)\n", i, pos.x, pos.y, pos.z);
          }
      }

      DestroyPhysicsWorld(world);
      return 0;
  }

```

You could also add torque with AddSphereTorque(world, sphere, torque) to make it spin, but I have neglected rotation in the spheres side so I'm not sure how useful that will be. But do use AddBoxForce() and AddBoxTorque() for boxes.

## Static SDF collider

Signed distance functions can be defined to serve as colliders for your dynamic objects. Although a niche feature, for sure, it's here if you need it. It's very simple to use. You define your SDF, then, after you create the world with CreatePhysicsWorld, you call UseSDF() with the world and your function. And that's it.

```c
// Create your own function. You can make an entire scene here, fractals, primitive objects, etc. 
float WavyGroundSDF(Vec3 point) {
    float wavy_y = sinf(point.x) * sinf(point.z);
    return point.y - wavy_y;
}

// Then...
PhysicsWorld* world = CreatePhysicsWorld(50, 50, 500);
UseSDF(world, WavyGroundSDF);
```
The SDF-box collision detection is an approximation. The test is done on the vertices and the centers of the faces. It can look unnatural when the SDF's scale of detail matches the boxes colliding with it. 


## Settings

You can modify the parameters of the simulation using the settings struct:
```c
typedef struct {
     // Damping (most commonly adjusted)
     float linear_damping;      // 0.99f - reduces sliding
     float angular_damping;     // 0.99f - reduces spinning
     
     // Collision response
     float restitution;         // 0.0f - bounciness (0=no bounce, 1=full bounce)
     float friction;            // 0.3f - surface friction coefficient
     int solver_iterations;     // 8 - constraint solver iterations
     
     // Constraint solver settings (Baumgarte stabilization)
     float baumgarte_bias;      // 0.2f - position correction bias factor
     float allowed_penetration; // 0.01f - allowed penetration before correction
     float velocity_threshold;   // 1.0f - relative velocity threshold for restitution
     
     // Sleep system: I implemented while trying to fix jittering and then neglected it. I suggest ignoring it for now.
     float sleep_linear_threshold;   // 0.01f - linear velocity threshold for sleep
     float sleep_angular_threshold;  // 0.01f - angular velocity threshold for sleep
     float sleep_time_required;      // 0.5f - time at low energy before sleeping
     
     // Contact manifold settings (for box-box collisions)
     float manifold_contact_tolerance;   // 0.01f - tolerance for vertex-inside-box detection
     float manifold_penetration_tolerance; // -0.01f - minimum penetration to accept (negative = allow touching)
     int manifold_max_contacts;          // 4 - maximum contact points per manifold
     
     // Numerical tolerances (rarely changed)
     float collision_epsilon;   // 1e-6f - collision detection threshold
     float sdf_normal_epsilon;  // 0.001f - SDF normal estimation step size
     float vector_normalize_epsilon; // 1e-6f - vector normalization threshold
     float quaternion_epsilon;  // 1e-6f - quaternion operations threshold
 } PhysicsSettings;
```

So you can just call `world->settings.linear_damping = 0.9f`. I did make some functions to alter some of these parameters, but it's honestly easier to just modify the struct directly. 

# Demo GIFs


![demogif](./physicswip2.gif)
![demogif](./physicswip3.gif)
![demogif](./physicswip_constraints_velcontrol.gif)
The previous gif is testing how the simulation handles direct velocity control and if I can constraint objects to remain together with just forces.


# License

MIT License.

# Contributing

Found a bug? Have a feature request? Open an issue or submit a PR. 

Want to help me pay some bills? check out [my Patreon blog](https://www.patreon.com/Diegomakesgames), where I talk about game dev and other stuff.

