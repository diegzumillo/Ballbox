# Ballbox


![demogif](./physicswip2.gif)

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
- ✅ Force and torque application
- ✅ Single header file
- ✅ Data Oriented Design (to the best of my abilities, at least)
- ✅ Static SDFs (Signed Distance Functions. I need it for my other thing)
- ⏳ Constraints, like hinges, springs. So we can use boxes and spheres to build more complex objects.

## Quick Start

1. Include the header: `#include "ballbox.h"`
2. Create a physics world with `CreatePhysicsWorld()`
3. Add objects with `AddBox()` and `AddSphere()`
4. Update the simulation with `PhysicsStep()`

As for the arguments, I asked a generic AI to look at my test scene that uses raylib, "hey, make a shorter program that's easier to follow that illustrates ballbox". And it spat this out. It looks solid to me!

```c
#include <stdio.h>
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

      // Simulate for 1000 steps. (Or in a while loop, like any sane game developer. This AI is mad.)
      for (int i = 0; i < 1000; i++) {
          // Apply a sideways force for the first 30 steps
          if (i < 30) {
              AddSphereForce(world, sphere, (Vec3){2.0f, 0, 0});  // Push right
          }

          PhysicsStep(world, 1.0f/60.0f);  // You should probably soem kind of deltatime your engine provides

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

# Static SDF collider

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

# Notes

There are lots of improvements to be made here. Many parameters are still hard coded deep in the code, like restitution, drag, etc. Eventually these will be exposed in a more user-friendly way, and others will be assigned on an object-basis. 

The SDF-box collision detection is an approximation. The test is done on the vertices and the centers of the faces. It will look unnatural when the SDF's scale of detail matches the boxes colliding with it. 


## License

MIT License.

## Contributing

Found a bug? Have a feature request? Open an issue or submit a PR. 

Want to help me pay some bills? check out [my Patreon blog](https://www.patreon.com/Diegomakesgames), where I talk about game dev and other stuff.

