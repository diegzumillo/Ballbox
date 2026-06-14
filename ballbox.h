/*
 * BALLBOX - 3D Physics Engine
 * 
 * A lightweight, data-oriented 3D physics library featuring:
 * - 3D rigid body dynamics (spheres and boxes)
 * - Efficient collision detection and resolution
 * - Impulse-based physics simulation
 * - Structure-of-arrays (SoA) design for cache performance
 * 
 * Usage:
 *   #define BALLBOX_IMPLEMENTATION
 *   #include "ballbox.h"
 * 
 * License: Public Domain / MIT (choose your preference)
 */

#ifndef BALLBOX_H
#define BALLBOX_H

#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>

#define PHYSICS_PI 3.14159265358979323846f
#define PHYSICS_DEG2RAD (PHYSICS_PI / 180.0f)
#define PHYSICS_RAD2DEG (180.0f / PHYSICS_PI)

typedef struct {
    float x, y, z;
} Vec3;

typedef struct {
    float x, y, z, w;
} Quat;

typedef struct {
    Vec3 position;
    float radius;
} Sphere;

typedef struct {
    Vec3 position;
    Vec3 size;     // half-extents
    Quat rotation;
} Box;

typedef struct {
    int bodyA_type;     // 0 = sphere collider, 1 = box collider, 2 = SDF
    int bodyA_index;    // Rigid body index, or -1 for SDF
    int colliderA_index; // Collider index within its collider system, or -1 for SDF
    int bodyB_type;     // 0 = sphere collider, 1 = box collider, 2 = SDF
    int bodyB_index;    // Rigid body index, or -1 for SDF
    int colliderB_index; // Collider index within its collider system, or -1 for SDF
    Vec3 contact_point;
    Vec3 contact_normal; // From A to B
    float penetration;
} CollisionContact;

// Function pointer type for collision callbacks
typedef void (*OnCollision)(int selfIndex, int otherBodyType, int otherIndex, CollisionContact* contact);

typedef struct {
    Vec3* positions;
    Quat* rotations;
    Vec3* velocities;
    Vec3* forces;
    float* masses;
    Vec3* angular_velocities;
    Vec3* torques;
    Vec3* inertias;        // diagonal inertia tensor in local body space
    bool* is_static;
    bool* is_sleeping;
    float* sleep_timers;
    int count;
    int capacity;
} RigidBodySystem;

typedef struct {
    Vec3* positions;
    float* radii;
    Vec3* velocities;
    Vec3* forces;
    float* masses;
    Vec3* angular_velocities;
    Vec3* torques;
    float* inertias;
    bool* is_static;
    bool* is_sleeping;
    float* sleep_timers;
    OnCollision* collision_callbacks;
    int* body_indices;
    Vec3* local_positions;
    int count;
    int capacity;
} SphereSystem;

typedef struct {
    Vec3* positions;
    Vec3* sizes;           // half-extents
    Quat* rotations;
    Vec3* velocities;
    Vec3* forces;
    float* masses;
    Vec3* angular_velocities;
    Vec3* torques;
    Vec3* inertias;        // diagonal inertia tensor
    bool* is_static;
    bool* is_sleeping;
    float* sleep_timers;
    OnCollision* collision_callbacks;
    int* body_indices;
    Vec3* local_positions;
    Quat* local_rotations;
    int count;
    int capacity;
} BoxSystem;

typedef struct {
    CollisionContact* contacts;
    int count;
    int capacity;
} CollisionCollection;

// Contact manifold for box-box collisions (configurable max contact points)
#define MAX_MANIFOLD_CONTACTS 8  // Must be >= manifold_max_contacts setting
typedef struct {
    Vec3 points[MAX_MANIFOLD_CONTACTS];         // Contact points in world space
    float penetrations[MAX_MANIFOLD_CONTACTS];  // Penetration depth for each point
    Vec3 normal;           // Consistent contact normal for entire manifold
    int point_count;       // Number of valid points (0 to manifold_max_contacts)
    
    // Body identification
    int bodyA_type;        // 0 = sphere collider, 1 = box collider, 2 = SDF
    int bodyA_index;       // Rigid body index, or -1 for SDF
    int colliderA_index;   // Collider index within its collider system, or -1 for SDF
    int bodyB_type;        // 0 = sphere collider, 1 = box collider, 2 = SDF
    int bodyB_index;       // Rigid body index, or -1 for SDF
    int colliderB_index;   // Collider index within its collider system, or -1 for SDF
} ContactManifold;

// Contact constraint for Sequential Impulse solver
typedef struct {
    int bodyA_type;     // 0 = sphere collider, 1 = box collider, 2 = SDF
    int bodyA_index;    // Rigid body index, or -1 for SDF
    int colliderA_index; // Collider index within its collider system, or -1 for SDF
    int bodyB_type;     // 0 = sphere collider, 1 = box collider, 2 = SDF
    int bodyB_index;    // Rigid body index, or -1 for SDF
    int colliderB_index; // Collider index within its collider system, or -1 for SDF
    
    Vec3 contact_point;
    Vec3 contact_normal;    // From A to B
    float penetration;
    
    // Pre-computed constraint data
    float normal_mass;      // Effective mass for normal constraint
    float tangent_mass[2];  // Effective mass for friction constraints
    Vec3 tangent[2];        // Friction direction vectors
    
    // Accumulated impulses (for warm starting)
    float normal_impulse;
    float tangent_impulse[2];
    
    // Bias for position correction (Baumgarte stabilization)
    float position_bias;
    
    // Contact properties
    float restitution;
    float friction;
} ContactConstraint;

typedef struct {
    ContactConstraint* constraints;
    int count;
    int capacity;
} ConstraintCollection;

typedef enum {
    PHYSICS_JOINT_FIXED = 0,
    PHYSICS_JOINT_HINGE = 1,
    PHYSICS_JOINT_PISTON = 2
} PhysicsJointType;

typedef struct {
    PhysicsJointType type;
    int bodyA_index;
    int bodyB_index;

    Vec3 local_anchor_a;
    Vec3 local_anchor_b;

    // Hinge axis in each body's local space. Fixed joints use this only for
    // optional debug/introspection.
    Vec3 local_axis_a;
    Vec3 local_axis_b;

    // For fixed joints, preserves B's starting orientation relative to A.
    Quat reference_rotation;

    bool enable_motor;
    float motor_speed;
    float max_motor_impulse;

    // Piston spring settings. stiffness pulls toward target_travel along the
    // piston axis; damping resists relative velocity along that axis.
    float stiffness;
    float damping;

    // Piston joints constrain motion to local_axis_a/local_axis_b, with
    // optional travel limits and a soft target along that axis.
    float min_travel;
    float max_travel;
    float target_travel;
    float limit_restitution;
    bool enable_spring;
} PhysicsJoint;

typedef struct {
    PhysicsJoint* joints;
    int count;
    int capacity;
} JointSystem;

// Function pointer type for user-defined SDF.
//
// Ballbox treats negative values as "inside" and positive values as "outside"
// when PHYSICS_SDF_BOX_SIGNED_DISTANCE is used. PHYSICS_SDF_BOX_SURFACE_TRACE
// can also work with unsigned/implicit fields that only approach zero at the
// surface, which is useful for raymarched content that is not a mathematically
// strict signed-distance field.
typedef float (*SDFFunction)(Vec3 point);

// Optional analytic normal callback for SDF collisions.
//
// If provided through UseSDFNormal(), Ballbox uses this instead of finite
// differencing the SDF. This is usually faster and more stable for complex or
// approximate fields because finite differencing costs six SDF evaluations per
// normal and can become noisy near thin/implicit surfaces.
typedef Vec3 (*SDFNormalFunction)(Vec3 point);

typedef enum {
    // Fast path for true signed-distance fields. Box-SDF contacts are generated
    // by sampling boundary points and accepting samples with negative distance.
    PHYSICS_SDF_BOX_SIGNED_DISTANCE = 0,

    // Robust path for implicit, unsigned, or approximate fields. Box-SDF
    // contacts are generated by tracing from the box center to boundary samples
    // and looking for zero-surface hits or near-surface samples.
    PHYSICS_SDF_BOX_SURFACE_TRACE = 1
} PhysicsSDFBoxCollisionMode;

typedef enum {
    // 8 rays from the box center to its corners. Lowest SDF query count.
    PHYSICS_SDF_BOX_TRACE_CORNERS = 0,

    // 14 rays: corners plus face centers. Better for broad face contacts.
    PHYSICS_SDF_BOX_TRACE_CORNERS_AND_FACES = 1,

    // 26 rays: corners, face centers, and edge centers. Most stable, most expensive.
    PHYSICS_SDF_BOX_TRACE_FULL = 2
} PhysicsSDFBoxTraceSamples;

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
    
    // Sleep system
    float sleep_linear_threshold;   // 0.01f - linear velocity threshold for sleep
    float sleep_angular_threshold;  // 0.01f - angular velocity threshold for sleep
    float sleep_time_required;      // 0.5f - time at low energy before sleeping
    
    // Contact manifold settings (for box-box collisions)
    float manifold_contact_tolerance;   // 0.01f - tolerance for vertex-inside-box detection
    float manifold_penetration_tolerance; // -0.01f - minimum penetration to accept (negative = allow touching)
    int manifold_max_contacts;          // 4 - maximum contact points per manifold

    // Box-SDF collision settings.
    //
    // Default mode is PHYSICS_SDF_BOX_SURFACE_TRACE with corner rays only. That
    // costs more than signed-distance sampling, but it catches surfaces passing
    // through a box even when the field never returns negative values. For true
    // signed SDFs, switch to PHYSICS_SDF_BOX_SIGNED_DISTANCE for better speed.
    PhysicsSDFBoxCollisionMode sdf_box_collision_mode;
    PhysicsSDFBoxTraceSamples sdf_box_trace_samples; // Corners fastest; full adds face/edge centers
    int sdf_box_trace_steps;           // Samples per center-to-boundary ray
    int sdf_box_trace_refine_steps;    // Bisection steps after a sign crossing
    float sdf_box_surface_epsilon;     // Contact shell for implicit/non-distance surfaces
    int sdf_box_max_contacts;          // Maximum SDF contacts emitted per box
    
    // Numerical tolerances (rarely changed)
    float collision_epsilon;   // 1e-6f - collision detection threshold
    float sdf_normal_epsilon;  // 0.001f - SDF normal estimation step size
    float vector_normalize_epsilon; // 1e-6f - vector normalization threshold
    float quaternion_epsilon;  // 1e-6f - quaternion operations threshold
} PhysicsSettings;

typedef struct {
    RigidBodySystem* bodies;
    SphereSystem* spheres;
    BoxSystem* boxes;
    CollisionCollection* collisions;
    ConstraintCollection* constraints;  // Sequential Impulse constraints
    JointSystem* joints;
    Vec3 gravity;
    float dt;
    bool usesSDF;
    SDFFunction SDFCollider;
    SDFNormalFunction SDFNormal;
    PhysicsSettings settings;
} PhysicsWorld;

// Create/destroy physics world
PhysicsWorld* CreatePhysicsWorld(int max_spheres, int max_boxes, int max_collisions);
void DestroyPhysicsWorld(PhysicsWorld* world);

// Rigid body system management.
//
// Ballbox now separates rigid bodies from colliders. A rigid body owns motion
// state: position, rotation, velocity, mass, inertia, forces, and sleep/static
// flags. Colliders own shape data and are attached to a body by index.
//
// Simple convenience path:
//   int ball = AddSphere(world, pos, radius);        // body + one sphere collider
//   SetBodyMass(world, ball, 2.0f);
//   SetBodyInertiaAuto(world, ball);
//
// Compound body path:
//   int ship = AddRigidBody(world, pos, QuatIdentity());
//   AddBoxCollider(world, ship, local_pos_a, local_rot_a, half_extents_a);
//   AddSphereCollider(world, ship, local_pos_b, radius_b);
//   SetBodyMass(world, ship, 50.0f);
//   SetBodyInertiaAuto(world, ship);
//
// Collision detection runs against every collider, but collision response
// applies impulses to the owning rigid body. Colliders attached to the same
// body are ignored during collision collection.
RigidBodySystem* CreateRigidBodySystem(int capacity);
void DestroyRigidBodySystem(RigidBodySystem* system);
int AddRigidBody(PhysicsWorld* world, Vec3 position, Quat rotation);
void SetBodyMass(PhysicsWorld* world, int bodyIndex, float mass);
void SetBodyInertia(PhysicsWorld* world, int bodyIndex, Vec3 inertia);
void SetBodyInertiaAuto(PhysicsWorld* world, int bodyIndex);
void SetBodyStatic(PhysicsWorld* world, int bodyIndex, bool is_static);
void AddBodyForce(PhysicsWorld* world, int bodyIndex, Vec3 force);
void AddBodyTorque(PhysicsWorld* world, int bodyIndex, Vec3 torque);
void SetBodyVelocity(PhysicsWorld* world, int bodyIndex, Vec3 velocity);
void SetBodyAngularVelocity(PhysicsWorld* world, int bodyIndex, Vec3 angular_velocity);
void SyncCollidersFromBodies(PhysicsWorld* world);

// SDF system.
//
// UseSDF() registers a static scene collider. UseSDFNormal() is optional but
// recommended when the scene can provide a cheap analytic/estimated normal.
void UseSDF(PhysicsWorld* world, SDFFunction sdf_function);
void UseSDFNormal(PhysicsWorld* world, SDFNormalFunction sdf_normal_function);

// Physics settings convenience functions
void SetDamping(PhysicsWorld* world, float linear, float angular);
void SetRestitution(PhysicsWorld* world, float restitution);
void SetCorrectionParameters(PhysicsWorld* world, float baumgarte_bias, float allowed_penetration);
void SetSolverIterations(PhysicsWorld* world, int iterations);
void SetManifoldSettings(PhysicsWorld* world, float contact_tolerance, float penetration_tolerance, int max_contacts);

// Select how boxes collide against the registered SDF.
//
// PHYSICS_SDF_BOX_SIGNED_DISTANCE is faster and expects negative interior
// values. PHYSICS_SDF_BOX_SURFACE_TRACE is more expensive but supports
// unsigned/implicit fields by marching from the box center toward boundary
// samples.
void SetSDFBoxCollisionMode(PhysicsWorld* world, PhysicsSDFBoxCollisionMode mode);

// Tune the surface-trace mode.
//
// Approximate SDF query count per box is:
//   ray_count * (trace_steps + possible refine_steps) + contact_count normals
// where ray_count is 8, 14, or 26 depending on the sample mode. Supplying an
// SDF normal callback avoids six extra SDF evaluations per emitted contact.
void SetSDFSurfaceTraceSettings(PhysicsWorld* world, PhysicsSDFBoxTraceSamples samples, int trace_steps, 
                                int refine_steps, float surface_epsilon, int max_contacts);
void ResetPhysicsSettingsToDefaults(PhysicsWorld* world);

// Main simulation step - call once per frame
void PhysicsStep(PhysicsWorld* world, float deltaTime);

void ApplyForces(PhysicsWorld* world);
void IntegrateVelocities(PhysicsWorld* world);
void IntegratePositions(PhysicsWorld* world);
void CollectCollisions(PhysicsWorld* world);
void ResolveCollisions(PhysicsWorld* world);
void CleanupPhysics(PhysicsWorld* world);

// New constraint-based collision resolution
void GenerateContactConstraints(PhysicsWorld* world);
void SolveConstraints(PhysicsWorld* world);
void PreSolveConstraints(PhysicsWorld* world);
void SolveVelocityConstraints(PhysicsWorld* world);
void SolvePositionConstraints(PhysicsWorld* world);

// Constraint system management
ConstraintCollection* CreateConstraintCollection(int capacity);
void DestroyConstraintCollection(ConstraintCollection* constraints);
void ClearConstraints(ConstraintCollection* constraints);

// Joint system management.
//
// Joints connect rigid body indices, not collider indices. They are solved in
// the same iterative phase as contact constraints. Fixed joints constrain both
// anchors and relative orientation. Hinge joints constrain the anchors and keep
// the hinge axes aligned, leaving rotation around the hinge axis free.
// Piston joints constrain lateral anchor drift and relative orientation while
// allowing bounded travel along one axis, with optional spring and limit bounce.
JointSystem* CreateJointSystem(int capacity);
void DestroyJointSystem(JointSystem* system);
int AddFixedJoint(PhysicsWorld* world, int bodyA, int bodyB, Vec3 worldAnchor);
int AddHingeJoint(PhysicsWorld* world, int bodyA, int bodyB, Vec3 worldAnchor, Vec3 worldAxis);
int AddPistonJoint(PhysicsWorld* world, int bodyA, int bodyB, Vec3 worldAnchorA, Vec3 worldAnchorB,
                   Vec3 worldAxis, float min_travel, float max_travel);
void SetHingeMotor(PhysicsWorld* world, int jointIndex, bool enabled, float motor_speed, float max_motor_impulse);
void SetPistonSpring(PhysicsWorld* world, int jointIndex, bool enabled, float target_travel, float stiffness, float damping);
void SetPistonLimitBounce(PhysicsWorld* world, int jointIndex, float restitution);
void SolveJointVelocityConstraints(PhysicsWorld* world);
void SolveJointPositionConstraints(PhysicsWorld* world);

Vec3 GetBodyVelocityAtPoint(PhysicsWorld* world, int bodyType, int bodyIndex, Vec3 point);
float CalculateEffectiveMass(PhysicsWorld* world, CollisionContact* contact);
void ApplyImpulse(PhysicsWorld* world, int bodyType, int bodyIndex, Vec3 impulse, Vec3 contactPoint);
void PositionalCorrection(PhysicsWorld* world, CollisionContact* contact);

Vec3 Vec3Add(Vec3 a, Vec3 b);
Vec3 Vec3Sub(Vec3 a, Vec3 b);
Vec3 Vec3Scale(Vec3 v, float s);
float Vec3Dot(Vec3 a, Vec3 b);
Vec3 Vec3Cross(Vec3 a, Vec3 b);
float Vec3Length(Vec3 v);
Vec3 Vec3Normalize(Vec3 v);

Quat QuatIdentity();
Quat QuatMultiply(Quat q1, Quat q2);
Quat QuatNormalize(Quat q);
Quat QuatConjugate(Quat q);
Quat QuatFromAxisAngle(Vec3 axis, float angle);
Quat QuatFromEuler(float pitch, float yaw, float roll);
void QuatToAxisAngle(Quat q, Vec3* axis, float* angle);
Vec3 QuatRotateVec3(Quat q, Vec3 v);

// Sphere collider system management.
//
// AddSphere() is a convenience wrapper that creates a rigid body and one sphere
// collider, returning the body index. AddSphereCollider() attaches a sphere
// collider to an existing body and returns the sphere collider index.
//
// The SetSphere*/AddSphere* compatibility functions take sphere collider
// indices. Prefer SetBody*/AddBody* when using the body index returned by
// AddSphere() or AddRigidBody().
SphereSystem* CreateSphereSystem(int capacity);
void DestroySphereSystem(SphereSystem* system);
int AddSphere(PhysicsWorld* world, Vec3 position, float radius); // convenience: creates body + collider, returns body index
int AddSphereCollider(PhysicsWorld* world, int bodyIndex, Vec3 localPosition, float radius);
void SetSphereMass(PhysicsWorld* world, int index, float mass); // compatibility: index is a sphere collider index
void SetSphereStatic(PhysicsWorld* world, int index, bool is_static);
void AddSphereForce(PhysicsWorld* world, int index, Vec3 force);
void AddSphereTorque(PhysicsWorld* world, int index, Vec3 torque);
bool CheckSphereCollision(Vec3 pos1, float radius1, Vec3 pos2, float radius2);
bool CheckSphereCollisionWithData(Vec3 pos1, float radius1, Vec3 pos2, float radius2, CollisionContact* contact);
bool CheckSphereSystemCollisions(SphereSystem* system, int sphere_index);

// Box collider system management.
//
// AddBox() is a convenience wrapper that creates a rigid body and one box
// collider, returning the body index. AddBoxCollider() attaches a box collider
// to an existing body and returns the box collider index.
//
// The SetBox*/AddBox* compatibility functions take box collider indices. Prefer
// SetBody*/AddBody* when using the body index returned by AddBox() or
// AddRigidBody().
BoxSystem* CreateBoxSystem(int capacity);
void DestroyBoxSystem(BoxSystem* system);
int AddBox(PhysicsWorld* world, Vec3 position, Vec3 size, Quat rotation); // convenience: creates body + collider, returns body index
int AddBoxCollider(PhysicsWorld* world, int bodyIndex, Vec3 localPosition, Quat localRotation, Vec3 size);
void SetBoxMass(PhysicsWorld* world, int index, float mass); // compatibility: index is a box collider index
void SetBoxStatic(PhysicsWorld* world, int index, bool is_static);
void AddBoxForce(PhysicsWorld* world, int index, Vec3 force);
void AddBoxTorque(PhysicsWorld* world, int index, Vec3 torque);
bool CheckBoxSystemCollisions(BoxSystem* system, int box_index);

// Collision callback registration
void SetSphereCollisionCallback(PhysicsWorld* world, int sphereIndex, OnCollision callback);
void SetBoxCollisionCallback(PhysicsWorld* world, int boxIndex, OnCollision callback);
bool CheckSphereBoxCollision(Vec3 spherePos, float sphereRadius, Vec3 boxPos, Vec3 boxSize, Quat boxRot);
bool CheckSphereBoxCollisionWithData(Vec3 spherePos, float sphereRadius, Vec3 boxPos, Vec3 boxSize, Quat boxRot, CollisionContact* contact);
bool CheckBoxCollision(Vec3 pos1, Vec3 size1, Quat rot1, Vec3 pos2, Vec3 size2, Quat rot2);
bool CheckBoxCollisionWithData(Vec3 pos1, Vec3 size1, Quat rot1, Vec3 pos2, Vec3 size2, Quat rot2, CollisionContact* contact);

// Contact manifold system for box-box collisions
void GenerateBoxBoxManifold(PhysicsWorld* world, int boxIndexA, int boxIndexB, ContactManifold* manifold);
void SelectOptimalContactPoints(Vec3* candidate_points, float* penetrations, int candidate_count, 
                               ContactManifold* manifold, Vec3 normal, int max_contacts);
void AddManifoldToConstraints(PhysicsWorld* world, ContactManifold* manifold);

// SDF functions for static colliders
Vec3 SDF_EstimateNormal(PhysicsWorld* world, Vec3 point);
bool CheckBoxSDFCollision(PhysicsWorld* world, Vec3 boxPos, Vec3 boxSize, Quat boxRot, CollisionContact* contact);
void GenerateBoxSDFContacts(PhysicsWorld* world, int boxIndex);

#ifdef BALLBOX_IMPLEMENTATION

RigidBodySystem* CreateRigidBodySystem(int capacity) {
    RigidBodySystem* system = (RigidBodySystem*)malloc(sizeof(RigidBodySystem));
    if (!system) return NULL;
    
    system->positions = (Vec3*)malloc(sizeof(Vec3) * capacity);
    system->rotations = (Quat*)malloc(sizeof(Quat) * capacity);
    system->velocities = (Vec3*)malloc(sizeof(Vec3) * capacity);
    system->forces = (Vec3*)malloc(sizeof(Vec3) * capacity);
    system->masses = (float*)malloc(sizeof(float) * capacity);
    system->angular_velocities = (Vec3*)malloc(sizeof(Vec3) * capacity);
    system->torques = (Vec3*)malloc(sizeof(Vec3) * capacity);
    system->inertias = (Vec3*)malloc(sizeof(Vec3) * capacity);
    system->is_static = (bool*)malloc(sizeof(bool) * capacity);
    system->is_sleeping = (bool*)malloc(sizeof(bool) * capacity);
    system->sleep_timers = (float*)malloc(sizeof(float) * capacity);
    system->count = 0;
    system->capacity = capacity;
    
    if (!system->positions || !system->rotations || !system->velocities || !system->forces ||
        !system->masses || !system->angular_velocities || !system->torques || !system->inertias ||
        !system->is_static || !system->is_sleeping || !system->sleep_timers) {
        DestroyRigidBodySystem(system);
        return NULL;
    }
    
    for (int i = 0; i < capacity; i++) {
        system->positions[i] = (Vec3){ 0.0f, 0.0f, 0.0f };
        system->rotations[i] = QuatIdentity();
        system->velocities[i] = (Vec3){ 0.0f, 0.0f, 0.0f };
        system->forces[i] = (Vec3){ 0.0f, 0.0f, 0.0f };
        system->masses[i] = 1.0f;
        system->angular_velocities[i] = (Vec3){ 0.0f, 0.0f, 0.0f };
        system->torques[i] = (Vec3){ 0.0f, 0.0f, 0.0f };
        system->inertias[i] = (Vec3){ 1.0f, 1.0f, 1.0f };
        system->is_static[i] = false;
        system->is_sleeping[i] = false;
        system->sleep_timers[i] = 0.0f;
    }
    
    return system;
}

void DestroyRigidBodySystem(RigidBodySystem* system) {
    if (!system) return;
    
    if (system->positions) free(system->positions);
    if (system->rotations) free(system->rotations);
    if (system->velocities) free(system->velocities);
    if (system->forces) free(system->forces);
    if (system->masses) free(system->masses);
    if (system->angular_velocities) free(system->angular_velocities);
    if (system->torques) free(system->torques);
    if (system->inertias) free(system->inertias);
    if (system->is_static) free(system->is_static);
    if (system->is_sleeping) free(system->is_sleeping);
    if (system->sleep_timers) free(system->sleep_timers);
    free(system);
}

SphereSystem* CreateSphereSystem(int capacity) {
    SphereSystem* system = (SphereSystem*)malloc(sizeof(SphereSystem));
    if (!system) return NULL;
    
    system->positions = (Vec3*)malloc(sizeof(Vec3) * capacity);
    system->radii = (float*)malloc(sizeof(float) * capacity);
    system->velocities = (Vec3*)malloc(sizeof(Vec3) * capacity);
    system->forces = (Vec3*)malloc(sizeof(Vec3) * capacity);
    system->masses = (float*)malloc(sizeof(float) * capacity);
    system->angular_velocities = (Vec3*)malloc(sizeof(Vec3) * capacity);
    system->torques = (Vec3*)malloc(sizeof(Vec3) * capacity);
    system->inertias = (float*)malloc(sizeof(float) * capacity);
    system->is_static = (bool*)malloc(sizeof(bool) * capacity);
    system->is_sleeping = (bool*)malloc(sizeof(bool) * capacity);
    system->sleep_timers = (float*)malloc(sizeof(float) * capacity);
    system->collision_callbacks = (OnCollision*)malloc(sizeof(OnCollision) * capacity);
    system->body_indices = (int*)malloc(sizeof(int) * capacity);
    system->local_positions = (Vec3*)malloc(sizeof(Vec3) * capacity);
    system->count = 0;
    system->capacity = capacity;
    
    if (!system->positions || !system->radii || !system->velocities || !system->forces || 
        !system->masses || !system->angular_velocities || !system->torques || !system->inertias || 
        !system->is_static || !system->is_sleeping || !system->sleep_timers || !system->collision_callbacks ||
        !system->body_indices || !system->local_positions) {
        DestroySphereSystem(system);
        return NULL;
    }
    
    for (int i = 0; i < capacity; i++) {
        system->positions[i] = (Vec3){ 0.0f, 0.0f, 0.0f };
        system->radii[i] = 1.0f;
        system->velocities[i] = (Vec3){ 0.0f, 0.0f, 0.0f };
        system->forces[i] = (Vec3){ 0.0f, 0.0f, 0.0f };
        system->masses[i] = 1.0f;
        system->angular_velocities[i] = (Vec3){ 0.0f, 0.0f, 0.0f };
        system->torques[i] = (Vec3){ 0.0f, 0.0f, 0.0f };
        system->inertias[i] = 1.0f;
        system->is_static[i] = false;
        system->is_sleeping[i] = false;
        system->sleep_timers[i] = 0.0f;
        system->collision_callbacks[i] = NULL;
        system->body_indices[i] = -1;
        system->local_positions[i] = (Vec3){ 0.0f, 0.0f, 0.0f };
    }
    
    return system;
}

void DestroySphereSystem(SphereSystem* system) {
    if (!system) return;
    
    if (system->positions) free(system->positions);
    if (system->radii) free(system->radii);
    if (system->velocities) free(system->velocities);
    if (system->forces) free(system->forces);
    if (system->masses) free(system->masses);
    if (system->angular_velocities) free(system->angular_velocities);
    if (system->torques) free(system->torques);
    if (system->inertias) free(system->inertias);
    if (system->is_static) free(system->is_static);
    if (system->is_sleeping) free(system->is_sleeping);
    if (system->sleep_timers) free(system->sleep_timers);
    if (system->collision_callbacks) free(system->collision_callbacks);
    if (system->body_indices) free(system->body_indices);
    if (system->local_positions) free(system->local_positions);
    free(system);
}

int AddRigidBody(PhysicsWorld* world, Vec3 position, Quat rotation) {
    if (!world || !world->bodies) return -1;
    RigidBodySystem* system = world->bodies;
    if (system->count >= system->capacity) return -1;
    
    int index = system->count++;
    system->positions[index] = position;
    system->rotations[index] = QuatNormalize(rotation);
    system->velocities[index] = (Vec3){ 0.0f, 0.0f, 0.0f };
    system->forces[index] = (Vec3){ 0.0f, 0.0f, 0.0f };
    system->masses[index] = 1.0f;
    system->angular_velocities[index] = (Vec3){ 0.0f, 0.0f, 0.0f };
    system->torques[index] = (Vec3){ 0.0f, 0.0f, 0.0f };
    system->inertias[index] = (Vec3){ 1.0f, 1.0f, 1.0f };
    system->is_static[index] = false;
    system->is_sleeping[index] = false;
    system->sleep_timers[index] = 0.0f;
    return index;
}

void SetBodyMass(PhysicsWorld* world, int bodyIndex, float mass) {
    if (!world || !world->bodies || bodyIndex < 0 || bodyIndex >= world->bodies->count || mass <= 0.0f) return;
    world->bodies->masses[bodyIndex] = mass;
}

void SetBodyInertia(PhysicsWorld* world, int bodyIndex, Vec3 inertia) {
    if (!world || !world->bodies || bodyIndex < 0 || bodyIndex >= world->bodies->count) return;
    
    const float min_inertia = 1e-6f;
    world->bodies->inertias[bodyIndex] = (Vec3){
        fmaxf(inertia.x, min_inertia),
        fmaxf(inertia.y, min_inertia),
        fmaxf(inertia.z, min_inertia)
    };
}

void SetBodyStatic(PhysicsWorld* world, int bodyIndex, bool is_static) {
    if (!world || !world->bodies || bodyIndex < 0 || bodyIndex >= world->bodies->count) return;
    
    RigidBodySystem* system = world->bodies;
    system->is_static[bodyIndex] = is_static;
    
    if (is_static) {
        system->velocities[bodyIndex] = (Vec3){ 0.0f, 0.0f, 0.0f };
        system->angular_velocities[bodyIndex] = (Vec3){ 0.0f, 0.0f, 0.0f };
        system->forces[bodyIndex] = (Vec3){ 0.0f, 0.0f, 0.0f };
        system->torques[bodyIndex] = (Vec3){ 0.0f, 0.0f, 0.0f };
    }
}

void AddBodyForce(PhysicsWorld* world, int bodyIndex, Vec3 force) {
    if (!world || !world->bodies || bodyIndex < 0 || bodyIndex >= world->bodies->count) return;
    
    if (world->bodies->is_sleeping[bodyIndex]) {
        world->bodies->is_sleeping[bodyIndex] = false;
        world->bodies->sleep_timers[bodyIndex] = 0.0f;
    }
    
    world->bodies->forces[bodyIndex] = Vec3Add(world->bodies->forces[bodyIndex], force);
}

void AddBodyTorque(PhysicsWorld* world, int bodyIndex, Vec3 torque) {
    if (!world || !world->bodies || bodyIndex < 0 || bodyIndex >= world->bodies->count) return;
    
    if (world->bodies->is_sleeping[bodyIndex]) {
        world->bodies->is_sleeping[bodyIndex] = false;
        world->bodies->sleep_timers[bodyIndex] = 0.0f;
    }
    
    world->bodies->torques[bodyIndex] = Vec3Add(world->bodies->torques[bodyIndex], torque);
}

void SetBodyVelocity(PhysicsWorld* world, int bodyIndex, Vec3 velocity) {
    if (!world || !world->bodies || bodyIndex < 0 || bodyIndex >= world->bodies->count) return;
    world->bodies->velocities[bodyIndex] = velocity;
    world->bodies->is_sleeping[bodyIndex] = false;
    world->bodies->sleep_timers[bodyIndex] = 0.0f;
}

void SetBodyAngularVelocity(PhysicsWorld* world, int bodyIndex, Vec3 angular_velocity) {
    if (!world || !world->bodies || bodyIndex < 0 || bodyIndex >= world->bodies->count) return;
    world->bodies->angular_velocities[bodyIndex] = angular_velocity;
    world->bodies->is_sleeping[bodyIndex] = false;
    world->bodies->sleep_timers[bodyIndex] = 0.0f;
}

void SetBodyInertiaAuto(PhysicsWorld* world, int bodyIndex) {
    if (!world || !world->bodies || bodyIndex < 0 || bodyIndex >= world->bodies->count) return;
    
    float total_volume = 0.0f;
    for (int i = 0; world->spheres && i < world->spheres->count; i++) {
        if (world->spheres->body_indices[i] != bodyIndex) continue;
        float r = world->spheres->radii[i];
        total_volume += (4.0f / 3.0f) * PHYSICS_PI * r * r * r;
    }
    for (int i = 0; world->boxes && i < world->boxes->count; i++) {
        if (world->boxes->body_indices[i] != bodyIndex) continue;
        Vec3 s = world->boxes->sizes[i];
        total_volume += (s.x * 2.0f) * (s.y * 2.0f) * (s.z * 2.0f);
    }
    
    if (total_volume <= 0.0f) return;
    
    float body_mass = world->bodies->masses[bodyIndex];
    Vec3 inertia = { 0.0f, 0.0f, 0.0f };
    
    for (int i = 0; world->spheres && i < world->spheres->count; i++) {
        if (world->spheres->body_indices[i] != bodyIndex) continue;
        float r = world->spheres->radii[i];
        float volume = (4.0f / 3.0f) * PHYSICS_PI * r * r * r;
        float mass = body_mass * (volume / total_volume);
        float local_inertia = (2.0f / 5.0f) * mass * r * r;
        Vec3 offset = world->spheres->local_positions[i];
        float offset_sq = Vec3Dot(offset, offset);
        inertia.x += local_inertia + mass * (offset_sq - offset.x * offset.x);
        inertia.y += local_inertia + mass * (offset_sq - offset.y * offset.y);
        inertia.z += local_inertia + mass * (offset_sq - offset.z * offset.z);
        world->spheres->masses[i] = mass;
        world->spheres->inertias[i] = local_inertia;
    }
    
    for (int i = 0; world->boxes && i < world->boxes->count; i++) {
        if (world->boxes->body_indices[i] != bodyIndex) continue;
        Vec3 s = world->boxes->sizes[i];
        float w = s.x * 2.0f;
        float h = s.y * 2.0f;
        float d = s.z * 2.0f;
        float volume = w * h * d;
        float mass = body_mass * (volume / total_volume);
        Vec3 local_inertia = {
            (mass / 12.0f) * (h*h + d*d),
            (mass / 12.0f) * (w*w + d*d),
            (mass / 12.0f) * (w*w + h*h)
        };
        Vec3 offset = world->boxes->local_positions[i];
        float offset_sq = Vec3Dot(offset, offset);
        inertia.x += local_inertia.x + mass * (offset_sq - offset.x * offset.x);
        inertia.y += local_inertia.y + mass * (offset_sq - offset.y * offset.y);
        inertia.z += local_inertia.z + mass * (offset_sq - offset.z * offset.z);
        world->boxes->masses[i] = mass;
        world->boxes->inertias[i] = local_inertia;
    }
    
    SetBodyInertia(world, bodyIndex, inertia);
}

void SyncCollidersFromBodies(PhysicsWorld* world) {
    if (!world || !world->bodies) return;
    
    if (world->spheres) {
        for (int i = 0; i < world->spheres->count; i++) {
            int bodyIndex = world->spheres->body_indices[i];
            if (bodyIndex < 0 || bodyIndex >= world->bodies->count) continue;
            
            Vec3 local_world = QuatRotateVec3(world->bodies->rotations[bodyIndex], world->spheres->local_positions[i]);
            world->spheres->positions[i] = Vec3Add(world->bodies->positions[bodyIndex], local_world);
            world->spheres->velocities[i] = world->bodies->velocities[bodyIndex];
            world->spheres->forces[i] = world->bodies->forces[bodyIndex];
            world->spheres->masses[i] = world->bodies->masses[bodyIndex];
            world->spheres->angular_velocities[i] = world->bodies->angular_velocities[bodyIndex];
            world->spheres->torques[i] = world->bodies->torques[bodyIndex];
            world->spheres->is_static[i] = world->bodies->is_static[bodyIndex];
            world->spheres->is_sleeping[i] = world->bodies->is_sleeping[bodyIndex];
            world->spheres->sleep_timers[i] = world->bodies->sleep_timers[bodyIndex];
        }
    }
    
    if (world->boxes) {
        for (int i = 0; i < world->boxes->count; i++) {
            int bodyIndex = world->boxes->body_indices[i];
            if (bodyIndex < 0 || bodyIndex >= world->bodies->count) continue;
            
            Vec3 local_world = QuatRotateVec3(world->bodies->rotations[bodyIndex], world->boxes->local_positions[i]);
            world->boxes->positions[i] = Vec3Add(world->bodies->positions[bodyIndex], local_world);
            world->boxes->rotations[i] = QuatMultiply(world->bodies->rotations[bodyIndex], world->boxes->local_rotations[i]);
            world->boxes->velocities[i] = world->bodies->velocities[bodyIndex];
            world->boxes->forces[i] = world->bodies->forces[bodyIndex];
            world->boxes->masses[i] = world->bodies->masses[bodyIndex];
            world->boxes->angular_velocities[i] = world->bodies->angular_velocities[bodyIndex];
            world->boxes->torques[i] = world->bodies->torques[bodyIndex];
            world->boxes->is_static[i] = world->bodies->is_static[bodyIndex];
            world->boxes->is_sleeping[i] = world->bodies->is_sleeping[bodyIndex];
            world->boxes->sleep_timers[i] = world->bodies->sleep_timers[bodyIndex];
        }
    }
}

int AddSphere(PhysicsWorld* world, Vec3 position, float radius) {
    int bodyIndex = AddRigidBody(world, position, QuatIdentity());
    if (bodyIndex < 0) return -1;
    
    int colliderIndex = AddSphereCollider(world, bodyIndex, (Vec3){ 0.0f, 0.0f, 0.0f }, radius);
    if (colliderIndex < 0) {
        world->bodies->count--;
        return -1;
    }
    
    return bodyIndex;
}

int AddSphereCollider(PhysicsWorld* world, int bodyIndex, Vec3 localPosition, float radius) {
    if (!world || !world->spheres) return -1;
    SphereSystem* system = world->spheres;
    if (!world->bodies || bodyIndex < 0 || bodyIndex >= world->bodies->count) return -1;
    if (!system || system->count >= system->capacity) return -1;
    
    int index = system->count;
    system->body_indices[index] = bodyIndex;
    system->local_positions[index] = localPosition;
    system->positions[index] = Vec3Add(world->bodies->positions[bodyIndex], QuatRotateVec3(world->bodies->rotations[bodyIndex], localPosition));
    system->radii[index] = radius;
    system->velocities[index] = (Vec3){ 0.0f, 0.0f, 0.0f };
    system->forces[index] = (Vec3){ 0.0f, 0.0f, 0.0f };
    system->masses[index] = world->bodies->masses[bodyIndex];
    system->angular_velocities[index] = (Vec3){ 0.0f, 0.0f, 0.0f };
    system->torques[index] = (Vec3){ 0.0f, 0.0f, 0.0f };
    system->inertias[index] = (2.0f / 5.0f) * world->bodies->masses[bodyIndex] * radius * radius;
    world->bodies->inertias[bodyIndex] = (Vec3){ system->inertias[index], system->inertias[index], system->inertias[index] };
    system->is_static[index] = world->bodies->is_static[bodyIndex];
    system->is_sleeping[index] = world->bodies->is_sleeping[bodyIndex];
    system->sleep_timers[index] = world->bodies->sleep_timers[bodyIndex];
    system->collision_callbacks[index] = NULL;
    system->count++;
    
    return index;
}

void UpdateSphere(SphereSystem* system, int index, Vec3 position) {
    if (!system || index < 0 || index >= system->count) return;
    
    system->positions[index] = position;
}

void SetSphereMass(PhysicsWorld* world, int index, float mass) {
    if (!world || !world->spheres) return;
    SphereSystem* system = world->spheres;
    if (!system || index < 0 || index >= system->count || mass <= 0.0f) return;
    
    system->masses[index] = mass;
    int bodyIndex = system->body_indices[index];
    SetBodyMass(world, bodyIndex, mass);
    
    float radius = system->radii[index];
    system->inertias[index] = (2.0f / 5.0f) * mass * radius * radius;
    SetBodyInertia(world, bodyIndex, (Vec3){ system->inertias[index], system->inertias[index], system->inertias[index] });
}

void SetSphereStatic(PhysicsWorld* world, int index, bool is_static) {
    if (!world || !world->spheres) return;
    SphereSystem* system = world->spheres;
    if (!system || index < 0 || index >= system->count) return;
    
    system->is_static[index] = is_static;
    SetBodyStatic(world, system->body_indices[index], is_static);
    
    if (is_static) {
        system->velocities[index] = (Vec3){ 0.0f, 0.0f, 0.0f };
        system->angular_velocities[index] = (Vec3){ 0.0f, 0.0f, 0.0f };
        system->forces[index] = (Vec3){ 0.0f, 0.0f, 0.0f };
        system->torques[index] = (Vec3){ 0.0f, 0.0f, 0.0f };
    }
}

void AddSphereForce(PhysicsWorld* world, int index, Vec3 force) {
    if (!world || !world->spheres) return;
    SphereSystem* system = world->spheres;
    if (!system || index < 0 || index >= system->count) return;
    
    AddBodyForce(world, system->body_indices[index], force);
    system->forces[index] = Vec3Add(system->forces[index], force);
}

void AddSphereTorque(PhysicsWorld* world, int index, Vec3 torque) {
    if (!world || !world->spheres) return;
    SphereSystem* system = world->spheres;
    if (!system || index < 0 || index >= system->count) return;
    
    AddBodyTorque(world, system->body_indices[index], torque);
    system->torques[index] = Vec3Add(system->torques[index], torque);
}

void SetSphereCollisionCallback(PhysicsWorld* world, int sphereIndex, OnCollision callback) {
    if (!world || !world->spheres || sphereIndex < 0 || sphereIndex >= world->spheres->count) return;
    world->spheres->collision_callbacks[sphereIndex] = callback;
}

bool CheckSphereCollision(Vec3 pos1, float radius1, Vec3 pos2, float radius2) {
    return CheckSphereCollisionWithData(pos1, radius1, pos2, radius2, NULL);
}

bool CheckSphereCollisionWithData(Vec3 posA, float radiusA, Vec3 posB, float radiusB, CollisionContact* contact) {
    Vec3 n = Vec3Sub(posB, posA);
    
    float distance = Vec3Length(n);
    float combined_radius = radiusA + radiusB;
    
    if (distance < combined_radius) {
        if (contact) {
            if (distance > 1e-6f) {
                contact->contact_normal = Vec3Scale(n, 1.0f / distance);
            } else {
                contact->contact_normal = (Vec3){ 1.0f, 0.0f, 0.0f };
            }
            
            contact->contact_point = Vec3Add(posA, Vec3Scale(contact->contact_normal, radiusA));
            contact->penetration = combined_radius - distance;
        }
        return true;
    }
    return false;
}

bool CheckSphereSystemCollisions(SphereSystem* system, int sphere_index) {
    if (!system || sphere_index < 0 || sphere_index >= system->count) return false;
    
    Vec3 pos1 = system->positions[sphere_index];
    float radius1 = system->radii[sphere_index];
    
    for (int i = 0; i < system->count; i++) {
        if (i == sphere_index) continue;
        
        Vec3 pos2 = system->positions[i];
        float radius2 = system->radii[i];
        
        if (CheckSphereCollision(pos1, radius1, pos2, radius2)) {
            return true;
        }
    }
    
    return false;
}

float GetVec3Component(Vec3 v, int index) {
    switch(index) {
        case 0: return v.x;
        case 1: return v.y; 
        case 2: return v.z;
        default: return 0.0f;
    }
}
Vec3 Vec3Add(Vec3 a, Vec3 b) {
    Vec3 result = { a.x + b.x, a.y + b.y, a.z + b.z };
    return result;
}

Vec3 Vec3Sub(Vec3 a, Vec3 b) {
    Vec3 result = { a.x - b.x, a.y - b.y, a.z - b.z };
    return result;
}

Vec3 Vec3Scale(Vec3 v, float s) {
    Vec3 result = { v.x * s, v.y * s, v.z * s };
    return result;
}

float Vec3Dot(Vec3 a, Vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vec3 Vec3Cross(Vec3 a, Vec3 b) {
    Vec3 result = {
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    };
    return result;
}
float Vec3Length(Vec3 v) {
    return sqrtf(v.x * v.x + v.y * v.y + v.z * v.z);
}
Vec3 Vec3Normalize(Vec3 v) {
    float length = Vec3Length(v);
    if (length < 1e-6f) return (Vec3){0, 0, 0};
    return Vec3Scale(v, 1.0f / length);
}

// Quaternion operations
Quat QuatIdentity() {
    Quat result = { 0.0f, 0.0f, 0.0f, 1.0f };
    return result;
}

Quat QuatMultiply(Quat q1, Quat q2) {
    Quat result;
    result.x = q1.w * q2.x + q1.x * q2.w + q1.y * q2.z - q1.z * q2.y;
    result.y = q1.w * q2.y - q1.x * q2.z + q1.y * q2.w + q1.z * q2.x;
    result.z = q1.w * q2.z + q1.x * q2.y - q1.y * q2.x + q1.z * q2.w;
    result.w = q1.w * q2.w - q1.x * q2.x - q1.y * q2.y - q1.z * q2.z;
    return result;
}

Quat QuatNormalize(Quat q) {
    float length = sqrtf(q.x * q.x + q.y * q.y + q.z * q.z + q.w * q.w);
    if (length == 0.0f) return QuatIdentity();
    
    Quat result;
    result.x = q.x / length;
    result.y = q.y / length;
    result.z = q.z / length;
    result.w = q.w / length;
    return result;
}

Quat QuatFromAxisAngle(Vec3 axis, float angle) {
    Vec3 axisNorm = Vec3Scale(axis, 1.0f / sqrtf(Vec3Dot(axis, axis)));
    float halfAngle = angle * 0.5f;
    float sinHalf = sinf(halfAngle);
    
    Quat result;
    result.x = axisNorm.x * sinHalf;
    result.y = axisNorm.y * sinHalf;
    result.z = axisNorm.z * sinHalf;
    result.w = cosf(halfAngle);
    return result;
}

Quat QuatFromEuler(float pitch, float yaw, float roll) {
    float halfPitch = pitch * 0.5f;
    float halfYaw = yaw * 0.5f;
    float halfRoll = roll * 0.5f;
    
    float cosP = cosf(halfPitch);
    float sinP = sinf(halfPitch);
    float cosY = cosf(halfYaw);
    float sinY = sinf(halfYaw);
    float cosR = cosf(halfRoll);
    float sinR = sinf(halfRoll);
    
    Quat result;
    result.x = sinP * cosY * cosR - cosP * sinY * sinR;
    result.y = cosP * sinY * cosR + sinP * cosY * sinR;
    result.z = cosP * cosY * sinR - sinP * sinY * cosR;
    result.w = cosP * cosY * cosR + sinP * sinY * sinR;
    return result;
}

void QuatToAxisAngle(Quat q, Vec3* axis, float* angle) {
    // Normalize quaternion first
    q = QuatNormalize(q);
    
    // Handle identity quaternion (no rotation)
    if (fabsf(q.w) >= 1.0f) {
        *angle = 0.0f;
        axis->x = 0.0f;
        axis->y = 1.0f;  // Default to Y axis
        axis->z = 0.0f;
        return;
    }
    
    *angle = 2.0f * acosf(q.w);
    float sinHalfAngle = sqrtf(1.0f - q.w * q.w);
    
    if (sinHalfAngle < 1e-6f) {
        // Avoid division by zero
        axis->x = 0.0f;
        axis->y = 1.0f;
        axis->z = 0.0f;
    } else {
        axis->x = q.x / sinHalfAngle;
        axis->y = q.y / sinHalfAngle;
        axis->z = q.z / sinHalfAngle;
    }
}

Quat QuatConjugate(Quat q) {
    Quat conjugate = { -q.x, -q.y, -q.z, q.w };
    return conjugate;
}

Vec3 QuatRotateVec3(Quat q, Vec3 v) {
    // v' = q * v * q^-1 (quaternion rotation)
    Vec3 u = { q.x, q.y, q.z };
    float s = q.w;
    
    Vec3 vprime = Vec3Add(Vec3Add(Vec3Scale(v, s*s - Vec3Dot(u,u)), 
                                  Vec3Scale(u, 2.0f * Vec3Dot(u,v))), 
                          Vec3Scale(Vec3Cross(u,v), 2.0f * s));
    return vprime;
}

// Box system management
BoxSystem* CreateBoxSystem(int capacity) {
    BoxSystem* system = (BoxSystem*)malloc(sizeof(BoxSystem));
    if (!system) return NULL;
    
    system->positions = (Vec3*)malloc(sizeof(Vec3) * capacity);
    system->sizes = (Vec3*)malloc(sizeof(Vec3) * capacity);
    system->rotations = (Quat*)malloc(sizeof(Quat) * capacity);
    system->velocities = (Vec3*)malloc(sizeof(Vec3) * capacity);
    system->forces = (Vec3*)malloc(sizeof(Vec3) * capacity);
    system->masses = (float*)malloc(sizeof(float) * capacity);
    system->angular_velocities = (Vec3*)malloc(sizeof(Vec3) * capacity);
    system->torques = (Vec3*)malloc(sizeof(Vec3) * capacity);
    system->inertias = (Vec3*)malloc(sizeof(Vec3) * capacity);
    system->is_static = (bool*)malloc(sizeof(bool) * capacity);
    system->is_sleeping = (bool*)malloc(sizeof(bool) * capacity);
    system->sleep_timers = (float*)malloc(sizeof(float) * capacity);
    system->collision_callbacks = (OnCollision*)malloc(sizeof(OnCollision) * capacity);
    system->body_indices = (int*)malloc(sizeof(int) * capacity);
    system->local_positions = (Vec3*)malloc(sizeof(Vec3) * capacity);
    system->local_rotations = (Quat*)malloc(sizeof(Quat) * capacity);
    system->count = 0;
    system->capacity = capacity;
    
    if (!system->positions || !system->sizes || !system->rotations || 
        !system->velocities || !system->forces || !system->masses ||
        !system->angular_velocities || !system->torques || !system->inertias || 
        !system->is_static || !system->is_sleeping || !system->sleep_timers || !system->collision_callbacks ||
        !system->body_indices || !system->local_positions || !system->local_rotations) {
        DestroyBoxSystem(system);
        return NULL;
    }
    
    for (int i = 0; i < capacity; i++) {
        system->positions[i] = (Vec3){ 0.0f, 0.0f, 0.0f };
        system->sizes[i] = (Vec3){ 1.0f, 1.0f, 1.0f };
        system->rotations[i] = QuatIdentity();
        system->velocities[i] = (Vec3){ 0.0f, 0.0f, 0.0f };
        system->forces[i] = (Vec3){ 0.0f, 0.0f, 0.0f };
        system->masses[i] = 1.0f;
        system->angular_velocities[i] = (Vec3){ 0.0f, 0.0f, 0.0f };
        system->torques[i] = (Vec3){ 0.0f, 0.0f, 0.0f };
        system->inertias[i] = (Vec3){ 1.0f / 6.0f, 1.0f / 6.0f, 1.0f / 6.0f };
        system->is_static[i] = false;
        system->is_sleeping[i] = false;
        system->sleep_timers[i] = 0.0f;
        system->collision_callbacks[i] = NULL;
        system->body_indices[i] = -1;
        system->local_positions[i] = (Vec3){ 0.0f, 0.0f, 0.0f };
        system->local_rotations[i] = QuatIdentity();
    }
    
    return system;
}

void DestroyBoxSystem(BoxSystem* system) {
    if (!system) return;
    
    if (system->positions) free(system->positions);
    if (system->sizes) free(system->sizes);
    if (system->rotations) free(system->rotations);
    if (system->velocities) free(system->velocities);
    if (system->forces) free(system->forces);
    if (system->masses) free(system->masses);
    if (system->angular_velocities) free(system->angular_velocities);
    if (system->torques) free(system->torques);
    if (system->inertias) free(system->inertias);
    if (system->is_static) free(system->is_static);
    if (system->is_sleeping) free(system->is_sleeping);
    if (system->sleep_timers) free(system->sleep_timers);
    if (system->collision_callbacks) free(system->collision_callbacks);
    if (system->body_indices) free(system->body_indices);
    if (system->local_positions) free(system->local_positions);
    if (system->local_rotations) free(system->local_rotations);
    free(system);
}

int AddBox(PhysicsWorld* world, Vec3 position, Vec3 size, Quat rotation) {
    int bodyIndex = AddRigidBody(world, position, rotation);
    if (bodyIndex < 0) return -1;
    
    int colliderIndex = AddBoxCollider(world, bodyIndex, (Vec3){ 0.0f, 0.0f, 0.0f }, QuatIdentity(), size);
    if (colliderIndex < 0) {
        world->bodies->count--;
        return -1;
    }
    
    return bodyIndex;
}

int AddBoxCollider(PhysicsWorld* world, int bodyIndex, Vec3 localPosition, Quat localRotation, Vec3 size) {
    if (!world || !world->boxes) return -1;
    BoxSystem* system = world->boxes;
    if (!world->bodies || bodyIndex < 0 || bodyIndex >= world->bodies->count) return -1;
    if (!system || system->count >= system->capacity) return -1;
    
    int index = system->count;
    system->body_indices[index] = bodyIndex;
    system->local_positions[index] = localPosition;
    system->local_rotations[index] = QuatNormalize(localRotation);
    system->positions[index] = Vec3Add(world->bodies->positions[bodyIndex], QuatRotateVec3(world->bodies->rotations[bodyIndex], localPosition));
    system->sizes[index] = size;
    system->rotations[index] = QuatMultiply(world->bodies->rotations[bodyIndex], system->local_rotations[index]);
    system->velocities[index] = world->bodies->velocities[bodyIndex];
    system->forces[index] = (Vec3){ 0.0f, 0.0f, 0.0f };
    system->masses[index] = world->bodies->masses[bodyIndex];
    system->angular_velocities[index] = world->bodies->angular_velocities[bodyIndex];
    system->torques[index] = (Vec3){ 0.0f, 0.0f, 0.0f };
    
    float m = world->bodies->masses[bodyIndex];
    float w = size.x * 2.0f;
    float h = size.y * 2.0f;
    float d = size.z * 2.0f;
    system->inertias[index] = (Vec3){
        (m / 12.0f) * (h*h + d*d),
        (m / 12.0f) * (w*w + d*d),
        (m / 12.0f) * (w*w + h*h)
    };
    world->bodies->inertias[bodyIndex] = system->inertias[index];
    
    // Initialize sleeping/static state explicitly
    system->is_static[index] = world->bodies->is_static[bodyIndex];
    system->is_sleeping[index] = world->bodies->is_sleeping[bodyIndex];
    system->sleep_timers[index] = world->bodies->sleep_timers[bodyIndex];
    
    // Ensure forces and torques are completely zeroed (double-check for first frame issues)
    system->forces[index] = (Vec3){ 0.0f, 0.0f, 0.0f };
    system->torques[index] = (Vec3){ 0.0f, 0.0f, 0.0f };
    
    system->count++;
    
    return index;
}

void UpdateBox(BoxSystem* system, int index, Vec3 position, Quat rotation) {
    if (!system || index < 0 || index >= system->count) return;
    
    system->positions[index] = position;
    system->rotations[index] = rotation;
}

void SetBoxMass(PhysicsWorld* world, int index, float mass) {
    if (!world || !world->boxes) return;
    BoxSystem* system = world->boxes;
    if (!system || index < 0 || index >= system->count || mass <= 0.0f) return;
    
    system->masses[index] = mass;
    int bodyIndex = system->body_indices[index];
    SetBodyMass(world, bodyIndex, mass);
    
    Vec3 size = system->sizes[index];
    float w = size.x * 2.0f;
    float h = size.y * 2.0f;
    float d = size.z * 2.0f;
    system->inertias[index] = (Vec3){
        (mass / 12.0f) * (h*h + d*d),
        (mass / 12.0f) * (w*w + d*d),
        (mass / 12.0f) * (w*w + h*h)
    };
    SetBodyInertia(world, bodyIndex, system->inertias[index]);
}

void SetBoxStatic(PhysicsWorld* world, int index, bool is_static) {
    if (!world || !world->boxes) return;
    BoxSystem* system = world->boxes;
    if (!system || index < 0 || index >= system->count) return;
    
    system->is_static[index] = is_static;
    SetBodyStatic(world, system->body_indices[index], is_static);
    
    if (is_static) {
        system->velocities[index] = (Vec3){ 0.0f, 0.0f, 0.0f };
        system->angular_velocities[index] = (Vec3){ 0.0f, 0.0f, 0.0f };
        system->forces[index] = (Vec3){ 0.0f, 0.0f, 0.0f };
        system->torques[index] = (Vec3){ 0.0f, 0.0f, 0.0f };
    }
}

void AddBoxForce(PhysicsWorld* world, int index, Vec3 force) {
    if (!world || !world->boxes) return;
    BoxSystem* system = world->boxes;
    if (!system || index < 0 || index >= system->count) return;
    
    AddBodyForce(world, system->body_indices[index], force);
    system->forces[index] = Vec3Add(system->forces[index], force);
}

void AddBoxTorque(PhysicsWorld* world, int index, Vec3 torque) {
    if (!world || !world->boxes) return;
    BoxSystem* system = world->boxes;
    if (!system || index < 0 || index >= system->count) return;
    
    AddBodyTorque(world, system->body_indices[index], torque);
    system->torques[index] = Vec3Add(system->torques[index], torque);
}

void SetBoxCollisionCallback(PhysicsWorld* world, int boxIndex, OnCollision callback) {
    if (!world || !world->boxes || boxIndex < 0 || boxIndex >= world->boxes->count) return;
    world->boxes->collision_callbacks[boxIndex] = callback;
}

typedef struct {
    Vec3 center;
    Vec3 halfExtents;  // Half-sizes along each axis
    Vec3 axes[3];      // Local x, y, z axes (unit vectors)
} OBB;

OBB CreateOBB(Vec3 position, Vec3 size, Quat rotation) {
    OBB obb;
    obb.center = position;
    obb.halfExtents = size;  // Your size is already half-extents
    obb.axes[0] = QuatRotateVec3(rotation, (Vec3){1, 0, 0});
    obb.axes[1] = QuatRotateVec3(rotation, (Vec3){0, 1, 0});
    obb.axes[2] = QuatRotateVec3(rotation, (Vec3){0, 0, 1});
    return obb;
}

// Get all 8 vertices of the OBB
void GetOBBVertices(const OBB* obb, Vec3 vertices[8]) {
    Vec3 hX = Vec3Scale(obb->axes[0], obb->halfExtents.x);
    Vec3 hY = Vec3Scale(obb->axes[1], obb->halfExtents.y);
    Vec3 hZ = Vec3Scale(obb->axes[2], obb->halfExtents.z);

    vertices[0] = Vec3Add(Vec3Add(Vec3Add(obb->center, hX), hY), hZ);
    vertices[1] = Vec3Add(Vec3Add(Vec3Sub(obb->center, hZ), hX), hY);
    vertices[2] = Vec3Add(Vec3Add(Vec3Add(obb->center, hX), hZ), Vec3Scale(hY, -1));
    vertices[3] = Vec3Add(Vec3Sub(Vec3Add(obb->center, hX), hY), Vec3Scale(hZ, -1));
    vertices[4] = Vec3Add(Vec3Add(Vec3Add(obb->center, hY), hZ), Vec3Scale(hX, -1));
    vertices[5] = Vec3Add(Vec3Sub(Vec3Add(obb->center, hY), hZ), Vec3Scale(hX, -1));
    vertices[6] = Vec3Add(Vec3Sub(Vec3Add(obb->center, hZ), hX), Vec3Scale(hY, -1));
    vertices[7] = Vec3Sub(Vec3Sub(Vec3Sub(obb->center, hX), hY), hZ);
}

// Project OBB onto an axis
void ProjectOBB(const OBB* obb, Vec3 axis, float* min, float* max) {
    float centerProjection = Vec3Dot(obb->center, axis);
    float extent = obb->halfExtents.x * fabsf(Vec3Dot(obb->axes[0], axis)) +
                   obb->halfExtents.y * fabsf(Vec3Dot(obb->axes[1], axis)) +
                   obb->halfExtents.z * fabsf(Vec3Dot(obb->axes[2], axis));
    
    *min = centerProjection - extent;
    *max = centerProjection + extent;
}

// Check overlap on a single axis and return penetration depth
bool OverlapOnAxis(const OBB* obb1, const OBB* obb2, Vec3 axis, float* penetrationDepth) {
    float min1, max1, min2, max2;
    ProjectOBB(obb1, axis, &min1, &max1);
    ProjectOBB(obb2, axis, &min2, &max2);
    
    // Check for separation
    if (max1 < min2 || max2 < min1) {
        return false; // No overlap
    }
    
    // Calculate penetration depth
    float overlap = fminf(max1, max2) - fmaxf(min1, min2);
    *penetrationDepth = overlap;
    
    // Preserve direction - if obb1 center is "before" obb2 center on this axis
    if (min1 < min2) {
        *penetrationDepth *= -1.0f;
    }
    
    return true;
}

// ============================================================================
// CONTACT GENERATION FUNCTIONS
// ============================================================================

// Check if a point is inside the bounds of a face
bool IsPointInFaceBounds(Vec3 point, Vec3 faceCenter, Vec3 u, Vec3 v, float uHalf, float vHalf) {
    Vec3 U = Vec3Normalize(u);
    Vec3 V = Vec3Normalize(v);
    
    Vec3 rel = Vec3Sub(point, faceCenter);
    float uCoord = Vec3Dot(rel, U);
    float vCoord = Vec3Dot(rel, V);
    
    return (fabsf(uCoord) <= uHalf && fabsf(vCoord) <= vHalf);
}

// Signed distance from point to plane
float SignedDistanceToPlane(Vec3 point, Vec3 planeOrigin, Vec3 planeNormal) {
    return Vec3Dot(Vec3Sub(point, planeOrigin), planeNormal);
}

// Vertex-face collision detection - checks all 6 faces of the face OBB
Vec3 VertexFaceCollision(const OBB* vertexOBB, const OBB* faceOBB, float* smallestPenetrationDepth) {
    Vec3 vertices[8];
    GetOBBVertices(vertexOBB, vertices);
    
    Vec3 closestVertex = {0, 0, 0};
    *smallestPenetrationDepth = 999999.0f;
    
    // Check all 6 faces of the face OBB (3 axes × 2 directions each)
    for (int i = 0; i < 3; i++) {
        Vec3 faceNormal = faceOBB->axes[i];
        Vec3 u = faceOBB->axes[(i + 1) % 3];
        Vec3 v = faceOBB->axes[(i + 2) % 3];
        float uHalf = GetVec3Component(faceOBB->halfExtents, (i + 1) % 3);
        float vHalf = GetVec3Component(faceOBB->halfExtents, (i + 2) % 3);
        
        // Check both faces along this axis (positive and negative)
        for (int faceDir = -1; faceDir <= 1; faceDir += 2) {
            Vec3 faceDirNormal = Vec3Scale(faceNormal, (float)faceDir);
            float faceOffset = GetVec3Component(faceOBB->halfExtents, i) * (float)faceDir;
            Vec3 faceCenter = Vec3Add(faceOBB->center, Vec3Scale(faceNormal, faceOffset));
            
            // Check each vertex from vertex OBB against this face
            for (int j = 0; j < 8; j++) {
                float signedDistance = SignedDistanceToPlane(vertices[j], faceCenter, faceDirNormal);
                
                // Only consider vertices that are penetrating (negative signed distance)
                if (signedDistance < 0.0f) {
                    float penetrationDepth = -signedDistance;
                    
                    // Check if vertex is within face bounds
                    if (IsPointInFaceBounds(vertices[j], faceCenter, u, v, uHalf, vHalf)) {
                        if (penetrationDepth < *smallestPenetrationDepth) {
                            *smallestPenetrationDepth = penetrationDepth;
                            closestVertex = vertices[j];
                        }
                    }
                }
            }
        }
    }
    
    return closestVertex;
}

// Squared distance between two line segments
float SquaredDistanceBetweenEdges(Vec3 p1, Vec3 q1, Vec3 p2, Vec3 q2) {
    Vec3 d1 = Vec3Sub(q1, p1);
    Vec3 d2 = Vec3Sub(q2, p2);
    Vec3 r = Vec3Sub(p1, p2);
    
    float a = Vec3Dot(d1, d1);
    float e = Vec3Dot(d2, d2);
    float f = Vec3Dot(d2, r);
    
    float s, t;
    if (a <= 1e-6f && e <= 1e-6f) {
        return Vec3Dot(r, r);
    }
    if (a <= 1e-6f) {
        s = 0.0f;
        t = fmaxf(0.0f, fminf(1.0f, f / e));
    } else {
        float c = Vec3Dot(d1, r);
        if (e <= 1e-6f) {
            t = 0.0f;
            s = fmaxf(0.0f, fminf(1.0f, -c / a));
        } else {
            float b = Vec3Dot(d1, d2);
            float denom = a * e - b * b;
            if (denom != 0.0f) {
                s = fmaxf(0.0f, fminf(1.0f, (b * f - c * e) / denom));
            } else {
                s = 0.0f;
            }
            t = fmaxf(0.0f, fminf(1.0f, (b * s + f) / e));
        }
    }
    
    Vec3 c1 = Vec3Add(p1, Vec3Scale(d1, s));
    Vec3 c2 = Vec3Add(p2, Vec3Scale(d2, t));
    Vec3 diff = Vec3Sub(c1, c2);
    return Vec3Dot(diff, diff);
}

// Find closest point between two edges
Vec3 ClosestPointBetweenEdges(const OBB* obb1, const OBB* obb2) {
    // Get edges for both OBBs (simplified - just use the main axes)
    Vec3 edges1[12], edges2[12];
    Vec3 vertices1[8], vertices2[8];
    
    GetOBBVertices(obb1, vertices1);
    GetOBBVertices(obb2, vertices2);
    
    float minDistSquared = 999999.0f;
    Vec3 closestPoint = {0, 0, 0};
    
    // This is simplified - in the full implementation you'd test all 12 edges of each box
    // For now, just test a few representative edges
    Vec3 edges1_simple[3] = {
        Vec3Sub(vertices1[1], vertices1[0]),  // X-axis edge
        Vec3Sub(vertices1[2], vertices1[0]),  // Y-axis edge  
        Vec3Sub(vertices1[4], vertices1[0])   // Z-axis edge
    };
    
    Vec3 edges2_simple[3] = {
        Vec3Sub(vertices2[1], vertices2[0]),
        Vec3Sub(vertices2[2], vertices2[0]),
        Vec3Sub(vertices2[4], vertices2[0])
    };
    
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            float distSquared = SquaredDistanceBetweenEdges(
                vertices1[0], Vec3Add(vertices1[0], edges1_simple[i]),
                vertices2[0], Vec3Add(vertices2[0], edges2_simple[j])
            );
            
            if (distSquared < minDistSquared) {
                minDistSquared = distSquared;
                // Calculate actual closest point (simplified)
                closestPoint = Vec3Scale(Vec3Add(obb1->center, obb2->center), 0.5f);
            }
        }
    }
    
    return closestPoint;
}

// ============================================================================
// SAT COLLISION DETECTION WITH CONTACT GENERATION
// ============================================================================

bool SATCollision(Vec3 pos1, Vec3 size1, Quat rot1, Vec3 pos2, Vec3 size2, Quat rot2, CollisionContact* contact) {
    OBB obb1 = CreateOBB(pos1, size1, rot1);
    OBB obb2 = CreateOBB(pos2, size2, rot2);
    
    Vec3 testAxes[15];
    float minPenetrationDepth = 999999.0f;
    Vec3 smallestAxis = {0, 0, 0};
    int axisType = -1;
    
    // 3 axes from OBB1 (face normals)
    testAxes[0] = obb1.axes[0];
    testAxes[1] = obb1.axes[1]; 
    testAxes[2] = obb1.axes[2];
    
    // 3 axes from OBB2 (face normals)
    testAxes[3] = obb2.axes[0];
    testAxes[4] = obb2.axes[1];
    testAxes[5] = obb2.axes[2];
    
    // 9 cross product axes (edge-edge)
    int idx = 6;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            Vec3 crossAxis = Vec3Cross(obb1.axes[i], obb2.axes[j]);
            float crossLength = Vec3Length(crossAxis);
            if (crossLength > 1e-6f) {
                testAxes[idx++] = Vec3Scale(crossAxis, 1.0f / crossLength);
            }
        }
    }
    
    // Test for overlap on all axes
    for (int i = 0; i < idx; i++) {
        Vec3 axis = testAxes[i];
        float penetrationDepth = 0.0f;
        
        if (!OverlapOnAxis(&obb1, &obb2, axis, &penetrationDepth)) {
            return false; // Separating axis found, no collision
        }
        
        // Track axis with smallest penetration depth
        if (fabsf(penetrationDepth) < fabsf(minPenetrationDepth)) {
            minPenetrationDepth = penetrationDepth;
            smallestAxis = Vec3Scale(axis, (penetrationDepth < 0) ? -1.0f : 1.0f);
            axisType = i;
        }
    }
    
    // Collision detected - generate contact data
    if (contact) {
        contact->penetration = fabsf(minPenetrationDepth);
        contact->contact_normal = smallestAxis;
        
        // Determine collision type and generate contact point
        if (axisType >= 0 && axisType <= 5) {
            // Vertex-face collision
            float distance1, distance2;
            Vec3 contactPoint1 = VertexFaceCollision(&obb1, &obb2, &distance1);
            Vec3 contactPoint2 = VertexFaceCollision(&obb2, &obb1, &distance2);
            
            contact->contact_point = (distance1 < distance2) ? contactPoint1 : contactPoint2;
        } else if (axisType >= 6 && axisType <= 14) {
            // Edge-edge collision
            contact->contact_point = ClosestPointBetweenEdges(&obb1, &obb2);
        } else {
            // Fallback to center point
            contact->contact_point = Vec3Scale(Vec3Add(pos1, pos2), 0.5f);
        }
    }
    
    return true;
}


bool CheckBoxCollision(Vec3 pos1, Vec3 size1, Quat rot1, Vec3 pos2, Vec3 size2, Quat rot2) {
    return CheckBoxCollisionWithData(pos1, size1, rot1, pos2, size2, rot2, NULL);
}
bool CheckBoxCollisionWithData(Vec3 posA, Vec3 sizeA, Quat rotA, Vec3 posB, Vec3 sizeB, Quat rotB, CollisionContact* contact) {
    // Run SAT collision test
    bool collision = SATCollision(posA, sizeA, rotA, posB, sizeB, rotB, contact);
    
    if (collision && contact) {
        // Randy's convention: normal should point from A to B
        Vec3 centerToCenter = Vec3Sub(posB, posA);
        
        // Check if normal is pointing the wrong way
        if (Vec3Dot(contact->contact_normal, centerToCenter) < 0.0f) {
            contact->contact_normal = Vec3Scale(contact->contact_normal, -1.0f);
        }
    }
    
    return collision;
}

// ============================================================================
// CONSTRAINT SYSTEM MANAGEMENT
// ============================================================================

ConstraintCollection* CreateConstraintCollection(int capacity) {
    ConstraintCollection* constraints = (ConstraintCollection*)malloc(sizeof(ConstraintCollection));
    if (!constraints) return NULL;
    
    constraints->constraints = (ContactConstraint*)malloc(sizeof(ContactConstraint) * capacity);
    if (!constraints->constraints) {
        free(constraints);
        return NULL;
    }
    
    constraints->count = 0;
    constraints->capacity = capacity;
    
    // Initialize constraints to prevent random garbage
    for (int i = 0; i < capacity; i++) {
        ContactConstraint* c = &constraints->constraints[i];
        memset(c, 0, sizeof(ContactConstraint));
    }
    
    return constraints;
}

void DestroyConstraintCollection(ConstraintCollection* constraints) {
    if (!constraints) return;
    
    if (constraints->constraints) free(constraints->constraints);
    free(constraints);
}

void ClearConstraints(ConstraintCollection* constraints) {
    if (!constraints) return;
    constraints->count = 0;
}

Vec3 GetConstraintBodyVelocityAtPoint(PhysicsWorld* world, int bodyType, int bodyIndex, Vec3 point);

JointSystem* CreateJointSystem(int capacity) {
    JointSystem* system = (JointSystem*)malloc(sizeof(JointSystem));
    if (!system) return NULL;

    system->joints = (PhysicsJoint*)malloc(sizeof(PhysicsJoint) * capacity);
    if (!system->joints) {
        free(system);
        return NULL;
    }

    system->count = 0;
    system->capacity = capacity;

    for (int i = 0; i < capacity; i++) {
        memset(&system->joints[i], 0, sizeof(PhysicsJoint));
        system->joints[i].bodyA_index = -1;
        system->joints[i].bodyB_index = -1;
        system->joints[i].reference_rotation = QuatIdentity();
    }

    return system;
}

void DestroyJointSystem(JointSystem* system) {
    if (!system) return;

    if (system->joints) free(system->joints);
    free(system);
}

static bool IsValidBodyIndex(PhysicsWorld* world, int bodyIndex) {
    return world && world->bodies && bodyIndex >= 0 && bodyIndex < world->bodies->count;
}

static Vec3 BodyWorldPoint(PhysicsWorld* world, int bodyIndex, Vec3 localPoint) {
    Vec3 rotated = QuatRotateVec3(world->bodies->rotations[bodyIndex], localPoint);
    return Vec3Add(world->bodies->positions[bodyIndex], rotated);
}

static Vec3 BodyLocalPoint(PhysicsWorld* world, int bodyIndex, Vec3 worldPoint) {
    Vec3 relative = Vec3Sub(worldPoint, world->bodies->positions[bodyIndex]);
    return QuatRotateVec3(QuatConjugate(world->bodies->rotations[bodyIndex]), relative);
}

static Vec3 BodyLocalVector(PhysicsWorld* world, int bodyIndex, Vec3 worldVector) {
    return QuatRotateVec3(QuatConjugate(world->bodies->rotations[bodyIndex]), worldVector);
}

int AddFixedJoint(PhysicsWorld* world, int bodyA, int bodyB, Vec3 worldAnchor) {
    if (!world || !world->joints || world->joints->count >= world->joints->capacity) return -1;
    if (!IsValidBodyIndex(world, bodyA) || !IsValidBodyIndex(world, bodyB) || bodyA == bodyB) return -1;

    int index = world->joints->count++;
    PhysicsJoint* joint = &world->joints->joints[index];
    memset(joint, 0, sizeof(PhysicsJoint));

    joint->type = PHYSICS_JOINT_FIXED;
    joint->bodyA_index = bodyA;
    joint->bodyB_index = bodyB;
    joint->local_anchor_a = BodyLocalPoint(world, bodyA, worldAnchor);
    joint->local_anchor_b = BodyLocalPoint(world, bodyB, worldAnchor);
    joint->local_axis_a = (Vec3){ 1.0f, 0.0f, 0.0f };
    joint->local_axis_b = BodyLocalVector(world, bodyB, QuatRotateVec3(world->bodies->rotations[bodyA], joint->local_axis_a));
    joint->reference_rotation = QuatMultiply(QuatConjugate(world->bodies->rotations[bodyA]), world->bodies->rotations[bodyB]);

    return index;
}

int AddHingeJoint(PhysicsWorld* world, int bodyA, int bodyB, Vec3 worldAnchor, Vec3 worldAxis) {
    if (!world || !world->joints || world->joints->count >= world->joints->capacity) return -1;
    if (!IsValidBodyIndex(world, bodyA) || !IsValidBodyIndex(world, bodyB) || bodyA == bodyB) return -1;

    Vec3 axis = Vec3Normalize(worldAxis);
    if (Vec3Length(axis) < world->settings.vector_normalize_epsilon) return -1;

    int index = world->joints->count++;
    PhysicsJoint* joint = &world->joints->joints[index];
    memset(joint, 0, sizeof(PhysicsJoint));

    joint->type = PHYSICS_JOINT_HINGE;
    joint->bodyA_index = bodyA;
    joint->bodyB_index = bodyB;
    joint->local_anchor_a = BodyLocalPoint(world, bodyA, worldAnchor);
    joint->local_anchor_b = BodyLocalPoint(world, bodyB, worldAnchor);
    joint->local_axis_a = BodyLocalVector(world, bodyA, axis);
    joint->local_axis_b = BodyLocalVector(world, bodyB, axis);
    joint->reference_rotation = QuatMultiply(QuatConjugate(world->bodies->rotations[bodyA]), world->bodies->rotations[bodyB]);
    joint->max_motor_impulse = 0.0f;

    return index;
}

int AddPistonJoint(PhysicsWorld* world, int bodyA, int bodyB, Vec3 worldAnchorA, Vec3 worldAnchorB,
                   Vec3 worldAxis, float min_travel, float max_travel) {
    if (!world || !world->joints || world->joints->count >= world->joints->capacity) return -1;
    if (!IsValidBodyIndex(world, bodyA) || !IsValidBodyIndex(world, bodyB) || bodyA == bodyB) return -1;

    Vec3 axis = Vec3Normalize(worldAxis);
    if (Vec3Length(axis) < world->settings.vector_normalize_epsilon) return -1;

    int index = world->joints->count++;
    PhysicsJoint* joint = &world->joints->joints[index];
    memset(joint, 0, sizeof(PhysicsJoint));

    joint->type = PHYSICS_JOINT_PISTON;
    joint->bodyA_index = bodyA;
    joint->bodyB_index = bodyB;
    joint->local_anchor_a = BodyLocalPoint(world, bodyA, worldAnchorA);
    joint->local_anchor_b = BodyLocalPoint(world, bodyB, worldAnchorB);
    joint->local_axis_a = BodyLocalVector(world, bodyA, axis);
    joint->local_axis_b = BodyLocalVector(world, bodyB, axis);
    joint->reference_rotation = QuatMultiply(QuatConjugate(world->bodies->rotations[bodyA]), world->bodies->rotations[bodyB]);
    joint->min_travel = fminf(min_travel, max_travel);
    joint->max_travel = fmaxf(min_travel, max_travel);
    joint->target_travel = fmaxf(joint->min_travel, fminf(joint->max_travel, Vec3Dot(Vec3Sub(worldAnchorB, worldAnchorA), axis)));
    joint->limit_restitution = 0.0f;
    joint->enable_spring = false;

    return index;
}

void SetHingeMotor(PhysicsWorld* world, int jointIndex, bool enabled, float motor_speed, float max_motor_impulse) {
    if (!world || !world->joints || jointIndex < 0 || jointIndex >= world->joints->count) return;

    PhysicsJoint* joint = &world->joints->joints[jointIndex];
    if (joint->type != PHYSICS_JOINT_HINGE) return;

    joint->enable_motor = enabled;
    joint->motor_speed = motor_speed;
    joint->max_motor_impulse = fmaxf(max_motor_impulse, 0.0f);
}

void SetPistonSpring(PhysicsWorld* world, int jointIndex, bool enabled, float target_travel, float stiffness, float damping) {
    if (!world || !world->joints || jointIndex < 0 || jointIndex >= world->joints->count) return;

    PhysicsJoint* joint = &world->joints->joints[jointIndex];
    if (joint->type != PHYSICS_JOINT_PISTON) return;

    joint->enable_spring = enabled;
    joint->target_travel = fmaxf(joint->min_travel, fminf(joint->max_travel, target_travel));
    joint->stiffness = fmaxf(stiffness, 0.0f);
    joint->damping = fmaxf(damping, 0.0f);
}

void SetPistonLimitBounce(PhysicsWorld* world, int jointIndex, float restitution) {
    if (!world || !world->joints || jointIndex < 0 || jointIndex >= world->joints->count) return;

    PhysicsJoint* joint = &world->joints->joints[jointIndex];
    if (joint->type != PHYSICS_JOINT_PISTON) return;

    joint->limit_restitution = fmaxf(0.0f, fminf(restitution, 1.0f));
}

// ============================================================================
// CONSTRAINT-BASED COLLISION RESOLUTION
// ============================================================================

void GenerateContactConstraints(PhysicsWorld* world) {
    if (!world || !world->collisions || !world->constraints) return;
    
    // Clear previous constraints
    ClearConstraints(world->constraints);
    
    
    // Generate one constraint per collision contact
    for (int i = 0; i < world->collisions->count; i++) {
        if (world->constraints->count >= world->constraints->capacity) break;
        
        CollisionContact* contact = &world->collisions->contacts[i];
        ContactConstraint* constraint = &world->constraints->constraints[world->constraints->count];
        
        // Copy basic contact data
        constraint->bodyA_type = contact->bodyA_type;
        constraint->bodyA_index = contact->bodyA_index;
        constraint->colliderA_index = contact->colliderA_index;
        constraint->bodyB_type = contact->bodyB_type;
        constraint->bodyB_index = contact->bodyB_index;
        constraint->colliderB_index = contact->colliderB_index;
        constraint->contact_point = contact->contact_point;
        constraint->contact_normal = contact->contact_normal;
        constraint->penetration = contact->penetration;
        
        // Set material properties (use world defaults for now)
        constraint->restitution = world->settings.restitution;
        constraint->friction = world->settings.friction;
        
        // Initialize accumulated impulses to zero (warm starting will be added later)
        constraint->normal_impulse = 0.0f;
        constraint->tangent_impulse[0] = 0.0f;
        constraint->tangent_impulse[1] = 0.0f;
        
        
        world->constraints->count++;
    }
}

void SolveConstraints(PhysicsWorld* world) {
    if (!world || !world->constraints) return;
    
    // Generate constraints from collisions
    GenerateContactConstraints(world);
    
    // Pre-solve: compute constraint properties
    PreSolveConstraints(world);
    
    // Solve velocity constraints iteratively
    for (int iter = 0; iter < world->settings.solver_iterations; iter++) {
        SolveVelocityConstraints(world);
        SolveJointVelocityConstraints(world);
    }
}

// ============================================================================
// CONSTRAINT SOLVER HELPER FUNCTIONS
// ============================================================================

// Get body velocity at a point
Vec3 GetConstraintBodyVelocityAtPoint(PhysicsWorld* world, int bodyType, int bodyIndex, Vec3 point) {
    if (bodyType == 2 || !world || !world->bodies) {
        return (Vec3){0, 0, 0};
    }
    if (bodyIndex < 0 || bodyIndex >= world->bodies->count || world->bodies->is_static[bodyIndex]) {
        return (Vec3){0, 0, 0};
    }
    
    Vec3 linear_vel = world->bodies->velocities[bodyIndex];
    Vec3 angular_vel = world->bodies->angular_velocities[bodyIndex];
    Vec3 body_pos = world->bodies->positions[bodyIndex];
    Vec3 r = Vec3Sub(point, body_pos);
    Vec3 angular_component = Vec3Cross(angular_vel, r);
    
    return Vec3Add(linear_vel, angular_component);
}

// Get body inverse mass
float GetConstraintBodyInverseMass(PhysicsWorld* world, int bodyType, int bodyIndex) {
    if (bodyType == 2 || !world || !world->bodies) return 0.0f;
    if (bodyIndex < 0 || bodyIndex >= world->bodies->count || world->bodies->is_static[bodyIndex]) {
        return 0.0f;
    }
    return 1.0f / world->bodies->masses[bodyIndex];
}

// Get body inverse inertia in local body space
Vec3 GetConstraintBodyInverseInertia(PhysicsWorld* world, int bodyType, int bodyIndex) {
    if (bodyType == 2 || !world || !world->bodies) return (Vec3){0, 0, 0};
    if (bodyIndex < 0 || bodyIndex >= world->bodies->count || world->bodies->is_static[bodyIndex]) {
        return (Vec3){0, 0, 0};
    }
    
    Vec3 inertia = world->bodies->inertias[bodyIndex];
    return (Vec3){
        (inertia.x > 1e-6f) ? 1.0f / inertia.x : 0.0f,
        (inertia.y > 1e-6f) ? 1.0f / inertia.y : 0.0f,
        (inertia.z > 1e-6f) ? 1.0f / inertia.z : 0.0f
    };
}

void PreSolveConstraints(PhysicsWorld* world) {
    if (!world || !world->constraints) return;
    
    float dt = world->dt;
    
    for (int i = 0; i < world->constraints->count; i++) {
        ContactConstraint* c = &world->constraints->constraints[i];
        
        // Get body data
        Vec3 posA = (c->bodyA_type == 2 || c->bodyA_index < 0) ? c->contact_point : world->bodies->positions[c->bodyA_index];
        Vec3 posB = (c->bodyB_type == 2 || c->bodyB_index < 0) ? c->contact_point : world->bodies->positions[c->bodyB_index];
        Vec3 rA = Vec3Sub(c->contact_point, posA);
        Vec3 rB = Vec3Sub(c->contact_point, posB);
        
        float invMassA = GetConstraintBodyInverseMass(world, c->bodyA_type, c->bodyA_index);
        float invMassB = GetConstraintBodyInverseMass(world, c->bodyB_type, c->bodyB_index);
        Vec3 invInertiaA = GetConstraintBodyInverseInertia(world, c->bodyA_type, c->bodyA_index);
        Vec3 invInertiaB = GetConstraintBodyInverseInertia(world, c->bodyB_type, c->bodyB_index);
        
        Vec3 normal = c->contact_normal;
        
        // Calculate normal effective mass
        Vec3 rnA = Vec3Cross(rA, normal);
        Vec3 rnB = Vec3Cross(rB, normal);
        float invMassSum = invMassA + invMassB + 
                          Vec3Dot(rnA, (Vec3){rnA.x * invInertiaA.x, rnA.y * invInertiaA.y, rnA.z * invInertiaA.z}) +
                          Vec3Dot(rnB, (Vec3){rnB.x * invInertiaB.x, rnB.y * invInertiaB.y, rnB.z * invInertiaB.z});
        c->normal_mass = (invMassSum > 0.0f) ? 1.0f / invMassSum : 0.0f;
        
        // Generate friction tangent vectors
        Vec3 tangent1, tangent2;
        if (fabsf(normal.x) >= 0.57735f) { // 1/sqrt(3)
            tangent1 = (Vec3){normal.y, -normal.x, 0.0f};
        } else {
            tangent1 = (Vec3){0.0f, normal.z, -normal.y};
        }
        tangent1 = Vec3Normalize(tangent1);
        tangent2 = Vec3Cross(normal, tangent1);
        
        c->tangent[0] = tangent1;
        c->tangent[1] = tangent2;
        
        // Calculate tangent effective masses
        for (int j = 0; j < 2; j++) {
            Vec3 rt_A = Vec3Cross(rA, c->tangent[j]);
            Vec3 rt_B = Vec3Cross(rB, c->tangent[j]);
            float tangent_invMassSum = invMassA + invMassB + 
                                     Vec3Dot(rt_A, (Vec3){rt_A.x * invInertiaA.x, rt_A.y * invInertiaA.y, rt_A.z * invInertiaA.z}) +
                                     Vec3Dot(rt_B, (Vec3){rt_B.x * invInertiaB.x, rt_B.y * invInertiaB.y, rt_B.z * invInertiaB.z});
            c->tangent_mass[j] = (tangent_invMassSum > 0.0f) ? 1.0f / tangent_invMassSum : 0.0f;
        }
        
        // Calculate position bias for Baumgarte stabilization
        float allowedPenetration = world->settings.allowed_penetration;
        float biasFactor = world->settings.baumgarte_bias / dt;
        c->position_bias = biasFactor * fmaxf(0.0f, c->penetration - allowedPenetration);
    }
}

// Apply impulse to a body in the constraint solver
void ApplyConstraintImpulse(PhysicsWorld* world, int bodyType, int bodyIndex, Vec3 linear_impulse, Vec3 angular_impulse) {
    if (bodyType == 2 || !world || !world->bodies) return;
    if (bodyIndex < 0 || bodyIndex >= world->bodies->count || world->bodies->is_static[bodyIndex]) {
        return;
    }
    
    float invMass = 1.0f / world->bodies->masses[bodyIndex];
    Vec3 invInertia = GetConstraintBodyInverseInertia(world, bodyType, bodyIndex);
    
    world->bodies->velocities[bodyIndex] = Vec3Add(world->bodies->velocities[bodyIndex],
                                                  Vec3Scale(linear_impulse, invMass));
    
    Quat rotation = world->bodies->rotations[bodyIndex];
    Vec3 localAngularImpulse = QuatRotateVec3(QuatConjugate(rotation), angular_impulse);
    Vec3 localAngularDelta = (Vec3){
        localAngularImpulse.x * invInertia.x,
        localAngularImpulse.y * invInertia.y,
        localAngularImpulse.z * invInertia.z
    };
    Vec3 worldAngularDelta = QuatRotateVec3(rotation, localAngularDelta);
    
    world->bodies->angular_velocities[bodyIndex] = Vec3Add(world->bodies->angular_velocities[bodyIndex],
                                                          worldAngularDelta);
    
    if (world->bodies->is_sleeping[bodyIndex]) {
        world->bodies->is_sleeping[bodyIndex] = false;
        world->bodies->sleep_timers[bodyIndex] = 0.0f;
    }
}

static Vec3 ApplyBodyInverseInertia(PhysicsWorld* world, int bodyIndex, Vec3 worldVector) {
    if (!IsValidBodyIndex(world, bodyIndex) || world->bodies->is_static[bodyIndex]) {
        return (Vec3){ 0.0f, 0.0f, 0.0f };
    }

    Vec3 invInertia = GetConstraintBodyInverseInertia(world, 0, bodyIndex);
    Quat rotation = world->bodies->rotations[bodyIndex];
    Vec3 localVector = QuatRotateVec3(QuatConjugate(rotation), worldVector);
    Vec3 localResult = {
        localVector.x * invInertia.x,
        localVector.y * invInertia.y,
        localVector.z * invInertia.z
    };
    return QuatRotateVec3(rotation, localResult);
}

static float CalculateJointLinearMass(PhysicsWorld* world, int bodyA, int bodyB, Vec3 rA, Vec3 rB, Vec3 axis) {
    float invMassA = GetConstraintBodyInverseMass(world, 0, bodyA);
    float invMassB = GetConstraintBodyInverseMass(world, 0, bodyB);
    Vec3 rnA = Vec3Cross(rA, axis);
    Vec3 rnB = Vec3Cross(rB, axis);

    float k = invMassA + invMassB +
              Vec3Dot(rnA, ApplyBodyInverseInertia(world, bodyA, rnA)) +
              Vec3Dot(rnB, ApplyBodyInverseInertia(world, bodyB, rnB));

    return (k > 1e-8f) ? 1.0f / k : 0.0f;
}

static float CalculateJointAngularMass(PhysicsWorld* world, int bodyA, int bodyB, Vec3 axis) {
    Vec3 invA = ApplyBodyInverseInertia(world, bodyA, axis);
    Vec3 invB = ApplyBodyInverseInertia(world, bodyB, axis);
    float k = Vec3Dot(axis, Vec3Add(invA, invB));
    return (k > 1e-8f) ? 1.0f / k : 0.0f;
}

static Vec3 GetPistonAxis(PhysicsWorld* world, PhysicsJoint* joint) {
    int bodyA = joint->bodyA_index;
    int bodyB = joint->bodyB_index;
    Vec3 axisA = Vec3Normalize(QuatRotateVec3(world->bodies->rotations[bodyA], joint->local_axis_a));
    Vec3 axisB = Vec3Normalize(QuatRotateVec3(world->bodies->rotations[bodyB], joint->local_axis_b));
    Vec3 axis = Vec3Normalize(Vec3Add(axisA, axisB));
    if (Vec3Length(axis) < world->settings.vector_normalize_epsilon) {
        axis = axisA;
    }
    if (Vec3Length(axis) < world->settings.vector_normalize_epsilon) {
        axis = (Vec3){ 1.0f, 0.0f, 0.0f };
    }
    return axis;
}

static void BuildPerpendicularBasis(Vec3 axis, Vec3* tangentA, Vec3* tangentB) {
    Vec3 reference = fabsf(axis.y) < 0.9f ? (Vec3){ 0.0f, 1.0f, 0.0f } : (Vec3){ 1.0f, 0.0f, 0.0f };
    *tangentA = Vec3Normalize(Vec3Cross(reference, axis));
    *tangentB = Vec3Normalize(Vec3Cross(axis, *tangentA));
}

static void SolveJointLinearVelocityAxis(PhysicsWorld* world, int bodyA, int bodyB,
                                         Vec3 anchorA, Vec3 anchorB, Vec3 axis,
                                         float positionError, float velocityBias) {
    Vec3 rA = Vec3Sub(anchorA, world->bodies->positions[bodyA]);
    Vec3 rB = Vec3Sub(anchorB, world->bodies->positions[bodyB]);
    Vec3 vA = GetConstraintBodyVelocityAtPoint(world, 0, bodyA, anchorA);
    Vec3 vB = GetConstraintBodyVelocityAtPoint(world, 0, bodyB, anchorB);
    float relVel = Vec3Dot(Vec3Sub(vB, vA), axis);
    float bias = (world->settings.baumgarte_bias / world->dt) * positionError + velocityBias;
    float mass = CalculateJointLinearMass(world, bodyA, bodyB, rA, rB, axis);
    float impulseMagnitude = -(relVel + bias) * mass;
    Vec3 impulse = Vec3Scale(axis, impulseMagnitude);

    ApplyConstraintImpulse(world, 0, bodyA, Vec3Scale(impulse, -1.0f), Vec3Cross(rA, Vec3Scale(impulse, -1.0f)));
    ApplyConstraintImpulse(world, 0, bodyB, impulse, Vec3Cross(rB, impulse));
}

static void SolveJointAnchorVelocity(PhysicsWorld* world, PhysicsJoint* joint) {
    int bodyA = joint->bodyA_index;
    int bodyB = joint->bodyB_index;
    Vec3 anchorA = BodyWorldPoint(world, bodyA, joint->local_anchor_a);
    Vec3 anchorB = BodyWorldPoint(world, bodyB, joint->local_anchor_b);
    Vec3 error = Vec3Sub(anchorB, anchorA);
    Vec3 axes[3] = {
        { 1.0f, 0.0f, 0.0f },
        { 0.0f, 1.0f, 0.0f },
        { 0.0f, 0.0f, 1.0f }
    };

    for (int i = 0; i < 3; i++) {
        Vec3 axis = axes[i];
        SolveJointLinearVelocityAxis(world, bodyA, bodyB, anchorA, anchorB, axis, Vec3Dot(error, axis), 0.0f);
    }
}

static void SolvePistonVelocity(PhysicsWorld* world, PhysicsJoint* joint) {
    int bodyA = joint->bodyA_index;
    int bodyB = joint->bodyB_index;
    Vec3 anchorA = BodyWorldPoint(world, bodyA, joint->local_anchor_a);
    Vec3 anchorB = BodyWorldPoint(world, bodyB, joint->local_anchor_b);
    Vec3 delta = Vec3Sub(anchorB, anchorA);
    Vec3 axis = GetPistonAxis(world, joint);
    Vec3 tangentA;
    Vec3 tangentB;
    BuildPerpendicularBasis(axis, &tangentA, &tangentB);

    SolveJointLinearVelocityAxis(world, bodyA, bodyB, anchorA, anchorB, tangentA, Vec3Dot(delta, tangentA), 0.0f);
    SolveJointLinearVelocityAxis(world, bodyA, bodyB, anchorA, anchorB, tangentB, Vec3Dot(delta, tangentB), 0.0f);

    Vec3 vA = GetConstraintBodyVelocityAtPoint(world, 0, bodyA, anchorA);
    Vec3 vB = GetConstraintBodyVelocityAtPoint(world, 0, bodyB, anchorB);
    float relVel = Vec3Dot(Vec3Sub(vB, vA), axis);
    float travel = Vec3Dot(delta, axis);

    if (joint->enable_spring && (joint->stiffness > 0.0f || joint->damping > 0.0f)) {
        Vec3 rA = Vec3Sub(anchorA, world->bodies->positions[bodyA]);
        Vec3 rB = Vec3Sub(anchorB, world->bodies->positions[bodyB]);
        float iterationScale = world->settings.solver_iterations > 0 ? 1.0f / (float)world->settings.solver_iterations : 1.0f;
        float displacement = travel - joint->target_travel;
        float springForce = -(joint->stiffness * displacement + joint->damping * relVel);
        Vec3 impulse = Vec3Scale(axis, springForce * world->dt * iterationScale);

        ApplyConstraintImpulse(world, 0, bodyA, Vec3Scale(impulse, -1.0f), Vec3Cross(rA, Vec3Scale(impulse, -1.0f)));
        ApplyConstraintImpulse(world, 0, bodyB, impulse, Vec3Cross(rB, impulse));
    }

    vA = GetConstraintBodyVelocityAtPoint(world, 0, bodyA, anchorA);
    vB = GetConstraintBodyVelocityAtPoint(world, 0, bodyB, anchorB);
    relVel = Vec3Dot(Vec3Sub(vB, vA), axis);
    travel = Vec3Dot(Vec3Sub(anchorB, anchorA), axis);

    if (travel < joint->min_travel) {
        float velocityBias = (relVel < -world->settings.velocity_threshold) ? joint->limit_restitution * relVel : 0.0f;
        SolveJointLinearVelocityAxis(world, bodyA, bodyB, anchorA, anchorB, axis, travel - joint->min_travel, velocityBias);
    } else if (travel > joint->max_travel) {
        float limitRelVel = -relVel;
        float velocityBias = (limitRelVel < -world->settings.velocity_threshold) ? joint->limit_restitution * limitRelVel : 0.0f;
        SolveJointLinearVelocityAxis(world, bodyA, bodyB, anchorA, anchorB, Vec3Scale(axis, -1.0f),
                                     joint->max_travel - travel, velocityBias);
    }
}

static void SolveJointAngularVelocityVector(PhysicsWorld* world, int bodyA, int bodyB, Vec3 correctionVelocity) {
    Vec3 axes[3] = {
        { 1.0f, 0.0f, 0.0f },
        { 0.0f, 1.0f, 0.0f },
        { 0.0f, 0.0f, 1.0f }
    };

    for (int i = 0; i < 3; i++) {
        Vec3 axis = axes[i];
        float component = Vec3Dot(correctionVelocity, axis);
        float mass = CalculateJointAngularMass(world, bodyA, bodyB, axis);
        float impulseMagnitude = -component * mass;
        Vec3 impulse = Vec3Scale(axis, impulseMagnitude);

        ApplyConstraintImpulse(world, 0, bodyA, (Vec3){ 0.0f, 0.0f, 0.0f }, Vec3Scale(impulse, -1.0f));
        ApplyConstraintImpulse(world, 0, bodyB, (Vec3){ 0.0f, 0.0f, 0.0f }, impulse);
    }
}

static Vec3 GetFixedJointAngularError(PhysicsWorld* world, PhysicsJoint* joint) {
    Quat targetB = QuatMultiply(world->bodies->rotations[joint->bodyA_index], joint->reference_rotation);
    Quat errorQuat = QuatMultiply(world->bodies->rotations[joint->bodyB_index], QuatConjugate(targetB));
    errorQuat = QuatNormalize(errorQuat);

    if (errorQuat.w < 0.0f) {
        errorQuat.x = -errorQuat.x;
        errorQuat.y = -errorQuat.y;
        errorQuat.z = -errorQuat.z;
        errorQuat.w = -errorQuat.w;
    }

    Vec3 axis;
    float angle;
    QuatToAxisAngle(errorQuat, &axis, &angle);
    if (angle > PHYSICS_PI) {
        angle -= 2.0f * PHYSICS_PI;
    }

    return Vec3Scale(axis, angle);
}

static void SolveJointAngularVelocity(PhysicsWorld* world, PhysicsJoint* joint) {
    int bodyA = joint->bodyA_index;
    int bodyB = joint->bodyB_index;
    Vec3 wA = world->bodies->angular_velocities[bodyA];
    Vec3 wB = world->bodies->angular_velocities[bodyB];
    Vec3 relW = Vec3Sub(wB, wA);
    float biasFactor = world->settings.baumgarte_bias / world->dt;

    if (joint->type == PHYSICS_JOINT_FIXED || joint->type == PHYSICS_JOINT_PISTON) {
        Vec3 angularError = GetFixedJointAngularError(world, joint);
        Vec3 correctionVelocity = Vec3Add(relW, Vec3Scale(angularError, biasFactor));
        SolveJointAngularVelocityVector(world, bodyA, bodyB, correctionVelocity);
    } else if (joint->type == PHYSICS_JOINT_HINGE) {
        Vec3 axisA = Vec3Normalize(QuatRotateVec3(world->bodies->rotations[bodyA], joint->local_axis_a));
        Vec3 axisB = Vec3Normalize(QuatRotateVec3(world->bodies->rotations[bodyB], joint->local_axis_b));
        Vec3 hingeAxis = Vec3Normalize(Vec3Add(axisA, axisB));
        if (Vec3Length(hingeAxis) < world->settings.vector_normalize_epsilon) {
            hingeAxis = axisA;
        }

        Vec3 angularError = Vec3Cross(axisA, axisB);
        Vec3 relPerp = Vec3Sub(relW, Vec3Scale(hingeAxis, Vec3Dot(relW, hingeAxis)));
        Vec3 correctionVelocity = Vec3Add(relPerp, Vec3Scale(angularError, biasFactor));
        SolveJointAngularVelocityVector(world, bodyA, bodyB, correctionVelocity);

        if (joint->enable_motor && joint->max_motor_impulse > 0.0f) {
            float relSpeed = Vec3Dot(Vec3Sub(world->bodies->angular_velocities[bodyB],
                                             world->bodies->angular_velocities[bodyA]), hingeAxis);
            float mass = CalculateJointAngularMass(world, bodyA, bodyB, hingeAxis);
            float impulseMagnitude = (joint->motor_speed - relSpeed) * mass;
            impulseMagnitude = fmaxf(-joint->max_motor_impulse, fminf(joint->max_motor_impulse, impulseMagnitude));
            Vec3 impulse = Vec3Scale(hingeAxis, impulseMagnitude);
            ApplyConstraintImpulse(world, 0, bodyA, (Vec3){ 0.0f, 0.0f, 0.0f }, Vec3Scale(impulse, -1.0f));
            ApplyConstraintImpulse(world, 0, bodyB, (Vec3){ 0.0f, 0.0f, 0.0f }, impulse);
        }
    }
}

void SolveJointVelocityConstraints(PhysicsWorld* world) {
    if (!world || !world->joints || !world->bodies) return;

    for (int i = 0; i < world->joints->count; i++) {
        PhysicsJoint* joint = &world->joints->joints[i];
        if (!IsValidBodyIndex(world, joint->bodyA_index) || !IsValidBodyIndex(world, joint->bodyB_index)) continue;

        if (joint->type == PHYSICS_JOINT_FIXED || joint->type == PHYSICS_JOINT_HINGE) {
            SolveJointAnchorVelocity(world, joint);
            SolveJointAngularVelocity(world, joint);
        } else if (joint->type == PHYSICS_JOINT_PISTON) {
            SolvePistonVelocity(world, joint);
            SolveJointAngularVelocity(world, joint);
        }
    }
}

static void ApplyBodyRotationCorrection(PhysicsWorld* world, int bodyIndex, Vec3 delta) {
    if (!IsValidBodyIndex(world, bodyIndex) || world->bodies->is_static[bodyIndex]) return;

    float angle = Vec3Length(delta);
    if (angle < 1e-6f) return;

    Vec3 axis = Vec3Scale(delta, 1.0f / angle);
    Quat dq = QuatFromAxisAngle(axis, angle);
    world->bodies->rotations[bodyIndex] = QuatNormalize(QuatMultiply(dq, world->bodies->rotations[bodyIndex]));
}

static void SolveJointAnchorPosition(PhysicsWorld* world, PhysicsJoint* joint) {
    int bodyA = joint->bodyA_index;
    int bodyB = joint->bodyB_index;
    Vec3 anchorA = BodyWorldPoint(world, bodyA, joint->local_anchor_a);
    Vec3 anchorB = BodyWorldPoint(world, bodyB, joint->local_anchor_b);
    Vec3 error = Vec3Sub(anchorB, anchorA);

    float invMassA = GetConstraintBodyInverseMass(world, 0, bodyA);
    float invMassB = GetConstraintBodyInverseMass(world, 0, bodyB);
    float totalInvMass = invMassA + invMassB;
    if (totalInvMass <= 1e-8f) return;

    Vec3 correction = Vec3Scale(error, 0.6f / totalInvMass);
    if (invMassA > 0.0f) {
        world->bodies->positions[bodyA] = Vec3Add(world->bodies->positions[bodyA], Vec3Scale(correction, invMassA));
    }
    if (invMassB > 0.0f) {
        world->bodies->positions[bodyB] = Vec3Sub(world->bodies->positions[bodyB], Vec3Scale(correction, invMassB));
    }
}

static void SolvePistonPosition(PhysicsWorld* world, PhysicsJoint* joint) {
    int bodyA = joint->bodyA_index;
    int bodyB = joint->bodyB_index;
    Vec3 anchorA = BodyWorldPoint(world, bodyA, joint->local_anchor_a);
    Vec3 anchorB = BodyWorldPoint(world, bodyB, joint->local_anchor_b);
    Vec3 delta = Vec3Sub(anchorB, anchorA);
    Vec3 axis = GetPistonAxis(world, joint);
    float travel = Vec3Dot(delta, axis);
    Vec3 lateralError = Vec3Sub(delta, Vec3Scale(axis, travel));
    Vec3 limitError = { 0.0f, 0.0f, 0.0f };

    if (travel < joint->min_travel) {
        limitError = Vec3Scale(axis, travel - joint->min_travel);
    } else if (travel > joint->max_travel) {
        limitError = Vec3Scale(axis, travel - joint->max_travel);
    }

    Vec3 error = Vec3Add(lateralError, limitError);
    float invMassA = GetConstraintBodyInverseMass(world, 0, bodyA);
    float invMassB = GetConstraintBodyInverseMass(world, 0, bodyB);
    float totalInvMass = invMassA + invMassB;
    if (totalInvMass <= 1e-8f) return;

    Vec3 correction = Vec3Scale(error, 0.6f / totalInvMass);
    if (invMassA > 0.0f) {
        world->bodies->positions[bodyA] = Vec3Add(world->bodies->positions[bodyA], Vec3Scale(correction, invMassA));
    }
    if (invMassB > 0.0f) {
        world->bodies->positions[bodyB] = Vec3Sub(world->bodies->positions[bodyB], Vec3Scale(correction, invMassB));
    }
}

static void SolveJointAngularPosition(PhysicsWorld* world, PhysicsJoint* joint) {
    int bodyA = joint->bodyA_index;
    int bodyB = joint->bodyB_index;
    float invMassA = GetConstraintBodyInverseMass(world, 0, bodyA);
    float invMassB = GetConstraintBodyInverseMass(world, 0, bodyB);
    float totalInvMass = invMassA + invMassB;
    if (totalInvMass <= 1e-8f) return;

    Vec3 error = { 0.0f, 0.0f, 0.0f };
    if (joint->type == PHYSICS_JOINT_FIXED || joint->type == PHYSICS_JOINT_PISTON) {
        error = GetFixedJointAngularError(world, joint);
    } else if (joint->type == PHYSICS_JOINT_HINGE) {
        Vec3 axisA = Vec3Normalize(QuatRotateVec3(world->bodies->rotations[bodyA], joint->local_axis_a));
        Vec3 axisB = Vec3Normalize(QuatRotateVec3(world->bodies->rotations[bodyB], joint->local_axis_b));
        error = Vec3Cross(axisA, axisB);
    }

    Vec3 correction = Vec3Scale(error, 0.4f / totalInvMass);
    ApplyBodyRotationCorrection(world, bodyA, Vec3Scale(correction, invMassA));
    ApplyBodyRotationCorrection(world, bodyB, Vec3Scale(correction, -invMassB));
}

void SolveJointPositionConstraints(PhysicsWorld* world) {
    if (!world || !world->joints || !world->bodies) return;

    for (int i = 0; i < world->joints->count; i++) {
        PhysicsJoint* joint = &world->joints->joints[i];
        if (!IsValidBodyIndex(world, joint->bodyA_index) || !IsValidBodyIndex(world, joint->bodyB_index)) continue;

        if (joint->type == PHYSICS_JOINT_FIXED || joint->type == PHYSICS_JOINT_HINGE) {
            SolveJointAnchorPosition(world, joint);
            SolveJointAngularPosition(world, joint);
        } else if (joint->type == PHYSICS_JOINT_PISTON) {
            SolvePistonPosition(world, joint);
            SolveJointAngularPosition(world, joint);
        }
    }

    SyncCollidersFromBodies(world);
}

void SolveVelocityConstraints(PhysicsWorld* world) {
    if (!world || !world->constraints) return;
    
    
    for (int i = 0; i < world->constraints->count; i++) {
        ContactConstraint* c = &world->constraints->constraints[i];
        
        // Get current body positions and velocities
        Vec3 posA = (c->bodyA_type == 2 || c->bodyA_index < 0) ? c->contact_point : world->bodies->positions[c->bodyA_index];
        Vec3 posB = (c->bodyB_type == 2 || c->bodyB_index < 0) ? c->contact_point : world->bodies->positions[c->bodyB_index];
        
        Vec3 rA = Vec3Sub(c->contact_point, posA);
        Vec3 rB = Vec3Sub(c->contact_point, posB);
        
        Vec3 vA = GetConstraintBodyVelocityAtPoint(world, c->bodyA_type, c->bodyA_index, c->contact_point);
        Vec3 vB = GetConstraintBodyVelocityAtPoint(world, c->bodyB_type, c->bodyB_index, c->contact_point);
        
        // Relative velocity
        Vec3 relativeVel = Vec3Sub(vB, vA);
        
        // === Solve Normal Constraint ===
        float normalVel = Vec3Dot(relativeVel, c->contact_normal);
        
        
        // Compute normal impulse (including restitution and bias)
        float restitution = c->restitution;
        float velocityBias = 0.0f;
        
        // Apply restitution only if relative velocity is above threshold
        if (normalVel < -world->settings.velocity_threshold) {
            velocityBias = -restitution * normalVel;
        }
        
        float jn = c->normal_mass * (-(normalVel + velocityBias) + c->position_bias);
        
        // Clamp accumulated impulse (no pulling apart)
        float oldNormalImpulse = c->normal_impulse;
        c->normal_impulse = fmaxf(0.0f, c->normal_impulse + jn);
        jn = c->normal_impulse - oldNormalImpulse;
        
        // Apply normal impulse
        Vec3 normalImpulse = Vec3Scale(c->contact_normal, jn);
        ApplyConstraintImpulse(world, c->bodyA_type, c->bodyA_index, 
                              Vec3Scale(normalImpulse, -1.0f), Vec3Cross(rA, Vec3Scale(normalImpulse, -1.0f)));
        ApplyConstraintImpulse(world, c->bodyB_type, c->bodyB_index, 
                              normalImpulse, Vec3Cross(rB, normalImpulse));
        
        // === Solve Friction Constraints ===
        if (c->friction > 0.0f) {
            // Recalculate relative velocity after normal impulse
            vA = GetConstraintBodyVelocityAtPoint(world, c->bodyA_type, c->bodyA_index, c->contact_point);
            vB = GetConstraintBodyVelocityAtPoint(world, c->bodyB_type, c->bodyB_index, c->contact_point);
            relativeVel = Vec3Sub(vB, vA);
            
            for (int j = 0; j < 2; j++) {
                float tangentVel = Vec3Dot(relativeVel, c->tangent[j]);
                float jt = c->tangent_mass[j] * (-tangentVel);
                
                // Coulomb friction limit
                float maxFriction = c->friction * c->normal_impulse;
                
                // Clamp tangent impulse
                float oldTangentImpulse = c->tangent_impulse[j];
                c->tangent_impulse[j] = fmaxf(-maxFriction, fminf(maxFriction, c->tangent_impulse[j] + jt));
                jt = c->tangent_impulse[j] - oldTangentImpulse;
                
                // Apply tangent impulse
                Vec3 tangentImpulse = Vec3Scale(c->tangent[j], jt);
                ApplyConstraintImpulse(world, c->bodyA_type, c->bodyA_index, 
                                      Vec3Scale(tangentImpulse, -1.0f), Vec3Cross(rA, Vec3Scale(tangentImpulse, -1.0f)));
                ApplyConstraintImpulse(world, c->bodyB_type, c->bodyB_index, 
                                      tangentImpulse, Vec3Cross(rB, tangentImpulse));
            }
        }
    }
}

void SolvePositionConstraints(PhysicsWorld* world) {
    if (!world || !world->constraints) return;
    
    // Direct position correction without affecting velocities
    // This helps reduce penetration without velocity-based oscillations
    
    for (int i = 0; i < world->constraints->count; i++) {
        ContactConstraint* c = &world->constraints->constraints[i];
        
        if (c->penetration <= world->settings.allowed_penetration) continue;
        
        // Get body masses
        float invMassA = GetConstraintBodyInverseMass(world, c->bodyA_type, c->bodyA_index);
        float invMassB = GetConstraintBodyInverseMass(world, c->bodyB_type, c->bodyB_index);
        float totalInvMass = invMassA + invMassB;
        
        if (totalInvMass < 1e-10f) continue;
        
        // Calculate position correction (smaller factor for stability)
        float correctionAmount = (c->penetration - world->settings.allowed_penetration) * 0.8f;
        Vec3 correction = Vec3Scale(c->contact_normal, correctionAmount / totalInvMass);
        
        // Apply position corrections directly to rigid bodies
        if (c->bodyA_type != 2 && c->bodyA_index >= 0 && invMassA > 0.0f) {
            Vec3 correctionA = Vec3Scale(correction, -invMassA);
            world->bodies->positions[c->bodyA_index] =
                Vec3Add(world->bodies->positions[c->bodyA_index], correctionA);
        }
        
        if (c->bodyB_type != 2 && c->bodyB_index >= 0 && invMassB > 0.0f) {
            Vec3 correctionB = Vec3Scale(correction, invMassB);
            world->bodies->positions[c->bodyB_index] =
                Vec3Add(world->bodies->positions[c->bodyB_index], correctionB);
        }
    }
    
    SyncCollidersFromBodies(world);
}



bool CheckBoxSystemCollisions(BoxSystem* system, int box_index) {
    if (!system || box_index < 0 || box_index >= system->count) return false;
    
    Vec3 pos1 = system->positions[box_index];
    Vec3 size1 = system->sizes[box_index];
    Quat rot1 = system->rotations[box_index];
    
    for (int i = 0; i < system->count; i++) {
        if (i == box_index) continue;
        
        Vec3 pos2 = system->positions[i];
        Vec3 size2 = system->sizes[i];
        Quat rot2 = system->rotations[i];
        
        if (CheckBoxCollision(pos1, size1, rot1, pos2, size2, rot2)) {
            return true;
        }
    }
    
    return false;
}

bool CheckSphereBoxCollision(Vec3 spherePos, float sphereRadius, Vec3 boxPos, Vec3 boxSize, Quat boxRot) {
    return CheckSphereBoxCollisionWithData(spherePos, sphereRadius, boxPos, boxSize, boxRot, NULL);
}

bool CheckSphereBoxCollisionWithData(Vec3 spherePos, float sphereRadius, Vec3 boxPos, Vec3 boxSize, Quat boxRot, CollisionContact* contact) {
    // Transform sphere position to box's local coordinate system
    // This makes the OBB effectively an AABB
    Vec3 sphereRelative = Vec3Sub(spherePos, boxPos);
    
    // Create inverse quaternion (conjugate since it's normalized)
    Quat invBoxRot = { -boxRot.x, -boxRot.y, -boxRot.z, boxRot.w };
    
    // Transform sphere position to box's local space
    Vec3 sphereLocal = QuatRotateVec3(invBoxRot, sphereRelative);
    
    // Find closest point on AABB (box in local space) to sphere center
    Vec3 closestPoint;
    closestPoint.x = fmaxf(-boxSize.x, fminf(sphereLocal.x, boxSize.x));
    closestPoint.y = fmaxf(-boxSize.y, fminf(sphereLocal.y, boxSize.y));
    closestPoint.z = fmaxf(-boxSize.z, fminf(sphereLocal.z, boxSize.z));
    
    // Calculate distance from sphere center to closest point
    Vec3 diff = Vec3Sub(sphereLocal, closestPoint);
    float distSquared = Vec3Dot(diff, diff);
    
    // Check if distance is within sphere radius
    if (distSquared <= (sphereRadius * sphereRadius)) {
        if (contact) {
            float distance = sqrtf(distSquared);
            
            if (distance > 1e-6f) {
                // Normal points from box to sphere (in local space)
                contact->contact_normal = Vec3Scale(diff, 1.0f / distance);
            } else {
                // Sphere center is inside box - use closest axis
                Vec3 absLocal = { fabsf(sphereLocal.x), fabsf(sphereLocal.y), fabsf(sphereLocal.z) };
                if (absLocal.x >= absLocal.y && absLocal.x >= absLocal.z) {
                    contact->contact_normal = sphereLocal.x > 0 ? (Vec3){1,0,0} : (Vec3){-1,0,0};
                } else if (absLocal.y >= absLocal.z) {
                    contact->contact_normal = sphereLocal.y > 0 ? (Vec3){0,1,0} : (Vec3){0,-1,0};
                } else {
                    contact->contact_normal = sphereLocal.z > 0 ? (Vec3){0,0,1} : (Vec3){0,0,-1};
                }
            }
            
            // Transform normal back to world space
            contact->contact_normal = QuatRotateVec3(boxRot, contact->contact_normal);
            
            // Randy's convention: normal points from sphere (A) to box (B)
            Vec3 sphereToBox = Vec3Sub(boxPos, spherePos);
            if (Vec3Dot(contact->contact_normal, sphereToBox) < 0.0f) {
                // Flip normal to point A->B (sphere->box)
                contact->contact_normal = Vec3Scale(contact->contact_normal, -1.0f);
            }
            
            // Contact point calculation (Randy's approach)
            Vec3 actualContact = Vec3Add(spherePos, 
                Vec3Scale(contact->contact_normal, -sphereRadius));
            contact->contact_point = actualContact;
            
            contact->penetration = sphereRadius - distance;
        }
        return true;
    }
    return false;
}

// Physics world management
PhysicsWorld* CreatePhysicsWorld(int max_spheres, int max_boxes, int max_collisions) {
    PhysicsWorld* world = (PhysicsWorld*)malloc(sizeof(PhysicsWorld));
    if (!world) return NULL;
    
    int max_bodies = max_spheres + max_boxes;
    if (max_bodies < 1) max_bodies = 1;
    world->bodies = CreateRigidBodySystem(max_bodies);
    world->spheres = CreateSphereSystem(max_spheres);
    world->boxes = CreateBoxSystem(max_boxes);
    
    // Create collision collection
    world->collisions = (CollisionCollection*)malloc(sizeof(CollisionCollection));
    if (world->collisions) {
        world->collisions->contacts = (CollisionContact*)malloc(sizeof(CollisionContact) * max_collisions);
        world->collisions->count = 0;
        world->collisions->capacity = max_collisions;
        
        // Initialize collision contacts to prevent random garbage
        for (int i = 0; i < max_collisions; i++) {
            world->collisions->contacts[i] = (CollisionContact){0};
        }
    }
    
    // Create constraint collection
    world->constraints = CreateConstraintCollection(max_collisions);
    world->joints = CreateJointSystem(max_bodies);
    
    // Set default values
    world->gravity = (Vec3){ 0.0f, 0.0f, 0.0f }; // No gravity for now
    world->dt = 1.0f / 60.0f;
    world->usesSDF = false;
    world->SDFCollider = NULL;
    world->SDFNormal = NULL;
    
    // Initialize physics settings with defaults
    ResetPhysicsSettingsToDefaults(world);
    
    // Check for allocation failures
    if (!world->bodies || !world->spheres || !world->boxes || !world->collisions || !world->collisions->contacts || !world->constraints || !world->joints) {
        DestroyPhysicsWorld(world);
        return NULL;
    }
    
    return world;
}

void DestroyPhysicsWorld(PhysicsWorld* world) {
    if (!world) return;
    
    if (world->bodies) DestroyRigidBodySystem(world->bodies);
    if (world->spheres) DestroySphereSystem(world->spheres);
    if (world->boxes) DestroyBoxSystem(world->boxes);
    
    if (world->collisions) {
        if (world->collisions->contacts) free(world->collisions->contacts);
        free(world->collisions);
    }
    
    if (world->constraints) DestroyConstraintCollection(world->constraints);
    if (world->joints) DestroyJointSystem(world->joints);
    
    
    free(world);
}

// Register user SDF function and enable SDF collision detection
void UseSDF(PhysicsWorld* world, SDFFunction sdf_function) {
    if (!world) return;
    
    world->SDFCollider = sdf_function;
    world->usesSDF = (sdf_function != NULL);
}

void UseSDFNormal(PhysicsWorld* world, SDFNormalFunction sdf_normal_function) {
    if (!world) return;
    
    world->SDFNormal = sdf_normal_function;
}

// Update sleep states for all rigid bodies
void UpdateSleepStates(PhysicsWorld* world) {
    if (!world || !world->bodies) return;
    
    float dt = world->dt;
    float linear_threshold = world->settings.sleep_linear_threshold;
    float angular_threshold = world->settings.sleep_angular_threshold;
    float time_required = world->settings.sleep_time_required;
    
    // Debug output (once per second)
    static float debug_timer = 0.0f;
    debug_timer += dt;
    if (debug_timer > 1.0f) {
        //printf("Sleep thresholds: linear=%.6f, angular=%.6f, time_required=%.3f\n", 
        //       linear_threshold, angular_threshold, time_required);
        debug_timer = 0.0f;
    }
    
    for (int i = 0; i < world->bodies->count; i++) {
        if (world->bodies->is_static[i]) continue;
        
        Vec3 velocity = world->bodies->velocities[i];
        Vec3 angular_velocity = world->bodies->angular_velocities[i];
        float linear_speed = Vec3Length(velocity);
        float angular_speed = Vec3Length(angular_velocity);
        bool low_energy = (linear_speed < linear_threshold && angular_speed < angular_threshold);
        
        if (low_energy) {
            world->bodies->sleep_timers[i] += dt;
            
            if (world->bodies->sleep_timers[i] >= time_required && !world->bodies->is_sleeping[i]) {
                world->bodies->is_sleeping[i] = true;
                world->bodies->velocities[i] = (Vec3){ 0.0f, 0.0f, 0.0f };
                world->bodies->angular_velocities[i] = (Vec3){ 0.0f, 0.0f, 0.0f };
            }
        } else {
            world->bodies->sleep_timers[i] = 0.0f;
            world->bodies->is_sleeping[i] = false;
        }
    }
    
    SyncCollidersFromBodies(world);
}

// Main physics simulation step - CORRECT Box2D pipeline order
void PhysicsStep(PhysicsWorld* world, float deltaTime) {
    if (!world) return;
    
    // Skip frame if deltaTime is zero or too small (Raylib first frame issue)
    if (deltaTime <= 0.0f || deltaTime < 1e-6f) {
        return;
    }
    
    world->dt = deltaTime;
    SyncCollidersFromBodies(world);
    
    // 1. Apply Forces
    ApplyForces(world);

    // 2. Integrate Velocities ONLY (not positions yet!)
    IntegrateVelocities(world);
    
    // 3. Collision Detection (at current positions)
    CollectCollisions(world);
    
    // 4. Sequential Impulse Constraint Solver (modifies velocities)
    SolveConstraints(world);
    
    // 5. Integrate Positions (after constraint solving)
    IntegratePositions(world);

    // 5b. Joint position stabilization
    SolveJointPositionConstraints(world);
    
    // 6. Update sleep states
    UpdateSleepStates(world);
    
    // 7. Cleanup
    CleanupPhysics(world);
}

// Physics sub-systems implementation
void ApplyForces(PhysicsWorld* world) {
    if (!world || !world->bodies) return;
    
    for (int i = 0; i < world->bodies->count; i++) {
        if (world->bodies->is_static[i]) continue;
        Vec3 gravityForce = Vec3Scale(world->gravity, world->bodies->masses[i]);
        world->bodies->forces[i] = Vec3Add(world->bodies->forces[i], gravityForce);
    }
}

// BOX2D PIPELINE: Integrate velocities first, positions later after constraint solving
void IntegrateVelocities(PhysicsWorld* world) {
    if (!world) return;
    if (world->bodies) {
        float dt = world->dt;
        
        for (int i = 0; i < world->bodies->count; i++) {
            if (world->bodies->is_static[i] || world->bodies->is_sleeping[i]) continue;
            
            Vec3 acceleration = Vec3Scale(world->bodies->forces[i], 1.0f / world->bodies->masses[i]);
            world->bodies->velocities[i] = Vec3Add(world->bodies->velocities[i], Vec3Scale(acceleration, dt));
            
            Quat rotation = world->bodies->rotations[i];
            Vec3 localTorque = QuatRotateVec3(QuatConjugate(rotation), world->bodies->torques[i]);
            Vec3 inertia = world->bodies->inertias[i];
            float min_inertia = 1e-6f;
            Vec3 localAngularAcceleration = {
                (inertia.x > min_inertia) ? localTorque.x / inertia.x : 0.0f,
                (inertia.y > min_inertia) ? localTorque.y / inertia.y : 0.0f,
                (inertia.z > min_inertia) ? localTorque.z / inertia.z : 0.0f
            };
            Vec3 angularAcceleration = QuatRotateVec3(rotation, localAngularAcceleration);
            world->bodies->angular_velocities[i] = Vec3Add(world->bodies->angular_velocities[i], Vec3Scale(angularAcceleration, dt));
        }
        
        float linear_damping = world->settings.linear_damping;
        float angular_damping = world->settings.angular_damping;
        
        for (int i = 0; i < world->bodies->count; i++) {
            if (world->bodies->is_static[i] || world->bodies->is_sleeping[i]) continue;
            world->bodies->velocities[i] = Vec3Scale(world->bodies->velocities[i], linear_damping);
            world->bodies->angular_velocities[i] = Vec3Scale(world->bodies->angular_velocities[i], angular_damping);
            world->bodies->forces[i] = (Vec3){ 0.0f, 0.0f, 0.0f };
            world->bodies->torques[i] = (Vec3){ 0.0f, 0.0f, 0.0f };
        }
        
        SyncCollidersFromBodies(world);
        return;
    }
    
    float dt = world->dt;
    
    // Integrate sphere velocities (linear + angular)
    for (int i = 0; i < world->spheres->count; i++) {
        // Skip static and sleeping objects
        if (world->spheres->is_static[i] || world->spheres->is_sleeping[i]) continue;
        
        // Calculate acceleration: a = F / m
        Vec3 acceleration = Vec3Scale(world->spheres->forces[i], 1.0f / world->spheres->masses[i]);
        
        // Integrate velocity: v += a * dt
        world->spheres->velocities[i] = Vec3Add(world->spheres->velocities[i], 
                                               Vec3Scale(acceleration, dt));
        
        // Angular integration for spheres
        // Calculate angular acceleration: α = τ / I
        Vec3 angular_acceleration = Vec3Scale(world->spheres->torques[i], 1.0f / world->spheres->inertias[i]);
        
        // Integrate angular velocity: ω += α * dt
        world->spheres->angular_velocities[i] = Vec3Add(world->spheres->angular_velocities[i], 
                                                       Vec3Scale(angular_acceleration, dt));
    }
    
    // Integrate box velocities (linear + angular)
    for (int i = 0; i < world->boxes->count; i++) {
        // Skip static and sleeping objects
        if (world->boxes->is_static[i] || world->boxes->is_sleeping[i]) continue;
        
        // Linear integration
        // Calculate acceleration: a = F / m
        Vec3 acceleration = Vec3Scale(world->boxes->forces[i], 1.0f / world->boxes->masses[i]);
        
        // Integrate velocity: v += a * dt
        world->boxes->velocities[i] = Vec3Add(world->boxes->velocities[i], 
                                             Vec3Scale(acceleration, dt));
        
        // Angular integration
        // Calculate angular acceleration in world space: α = I⁻¹ * τ 
        // For diagonal inertia tensor: α = τ / I (component-wise)
        Vec3 angular_acceleration;
        // Prevent division by zero with minimum inertia threshold
        float min_inertia = 1e-6f;
        angular_acceleration.x = (world->boxes->inertias[i].x > min_inertia) ? 
            world->boxes->torques[i].x / world->boxes->inertias[i].x : 0.0f;
        angular_acceleration.y = (world->boxes->inertias[i].y > min_inertia) ? 
            world->boxes->torques[i].y / world->boxes->inertias[i].y : 0.0f;
        angular_acceleration.z = (world->boxes->inertias[i].z > min_inertia) ? 
            world->boxes->torques[i].z / world->boxes->inertias[i].z : 0.0f;
        
        // Integrate angular velocity: ω += α * dt
        world->boxes->angular_velocities[i] = Vec3Add(world->boxes->angular_velocities[i], 
                                                     Vec3Scale(angular_acceleration, dt));
    }
    
    // Apply damping to velocities after integration
    float linear_damping = world->settings.linear_damping;
    float angular_damping = world->settings.angular_damping;
    
    // Apply damping to spheres
    for (int i = 0; i < world->spheres->count; i++) {
        if (world->spheres->is_static[i] || world->spheres->is_sleeping[i]) continue;
        world->spheres->velocities[i] = Vec3Scale(world->spheres->velocities[i], linear_damping);
        world->spheres->angular_velocities[i] = Vec3Scale(world->spheres->angular_velocities[i], angular_damping);
    }
    
    // Apply damping to boxes
    for (int i = 0; i < world->boxes->count; i++) {
        if (world->boxes->is_static[i] || world->boxes->is_sleeping[i]) continue;
        world->boxes->velocities[i] = Vec3Scale(world->boxes->velocities[i], linear_damping);
        world->boxes->angular_velocities[i] = Vec3Scale(world->boxes->angular_velocities[i], angular_damping);
    }
    
    // Clear accumulated forces and torques after integration
    for (int i = 0; i < world->spheres->count; i++) {
        world->spheres->forces[i] = (Vec3){ 0.0f, 0.0f, 0.0f };
        world->spheres->torques[i] = (Vec3){ 0.0f, 0.0f, 0.0f };
    }
    for (int i = 0; i < world->boxes->count; i++) {
        world->boxes->forces[i] = (Vec3){ 0.0f, 0.0f, 0.0f };
        world->boxes->torques[i] = (Vec3){ 0.0f, 0.0f, 0.0f };
    }
}

// BOX2D PIPELINE: Integrate positions after constraint solving
void IntegratePositions(PhysicsWorld* world) {
    if (!world) return;
    if (world->bodies) {
        float dt = world->dt;
        
        for (int i = 0; i < world->bodies->count; i++) {
            if (world->bodies->is_static[i] || world->bodies->is_sleeping[i]) continue;
            
            world->bodies->positions[i] = Vec3Add(world->bodies->positions[i], Vec3Scale(world->bodies->velocities[i], dt));
            
            Vec3 omega = world->bodies->angular_velocities[i];
            Quat omega_quat = { omega.x, omega.y, omega.z, 0.0f };
            Quat current_rotation = world->bodies->rotations[i];
            Quat q_dot = QuatMultiply(omega_quat, current_rotation);
            q_dot.x *= 0.5f;
            q_dot.y *= 0.5f;
            q_dot.z *= 0.5f;
            q_dot.w *= 0.5f;
            
            world->bodies->rotations[i].x = current_rotation.x + q_dot.x * dt;
            world->bodies->rotations[i].y = current_rotation.y + q_dot.y * dt;
            world->bodies->rotations[i].z = current_rotation.z + q_dot.z * dt;
            world->bodies->rotations[i].w = current_rotation.w + q_dot.w * dt;
            world->bodies->rotations[i] = QuatNormalize(world->bodies->rotations[i]);
        }
        
        SyncCollidersFromBodies(world);
        return;
    }
    
    float dt = world->dt;
    
    // Integrate sphere positions
    for (int i = 0; i < world->spheres->count; i++) {
        // Skip static and sleeping objects
        if (world->spheres->is_static[i] || world->spheres->is_sleeping[i]) continue;
        
        // Integrate position: pos += v * dt
        world->spheres->positions[i] = Vec3Add(world->spheres->positions[i], 
                                              Vec3Scale(world->spheres->velocities[i], dt));
        // Note: Spheres don't need orientation integration since they're rotationally symmetric
    }
    
    // Integrate box positions and orientations
    for (int i = 0; i < world->boxes->count; i++) {
        // Skip static and sleeping objects
        if (world->boxes->is_static[i] || world->boxes->is_sleeping[i]) continue;
        
        // Integrate position: pos += v * dt
        world->boxes->positions[i] = Vec3Add(world->boxes->positions[i], 
                                            Vec3Scale(world->boxes->velocities[i], dt));
        
        // Integrate rotation quaternion: q += 0.5 * ω * q * dt
        // Create angular velocity quaternion
        Vec3 omega = world->boxes->angular_velocities[i];
        Quat omega_quat = { omega.x, omega.y, omega.z, 0.0f };
        
        // Current rotation
        Quat current_rotation = world->boxes->rotations[i];
        
        // Calculate quaternion derivative: dq/dt = 0.5 * ω * q
        Quat q_dot = QuatMultiply(omega_quat, current_rotation);
        q_dot.x *= 0.5f;
        q_dot.y *= 0.5f; 
        q_dot.z *= 0.5f;
        q_dot.w *= 0.5f;
        
        // Integrate: q += dq/dt * dt
        world->boxes->rotations[i].x = current_rotation.x + q_dot.x * dt;
        world->boxes->rotations[i].y = current_rotation.y + q_dot.y * dt;
        world->boxes->rotations[i].z = current_rotation.z + q_dot.z * dt;
        world->boxes->rotations[i].w = current_rotation.w + q_dot.w * dt;
        
        // Normalize quaternion to prevent drift
        world->boxes->rotations[i] = QuatNormalize(world->boxes->rotations[i]);
    }
}

void CollectCollisions(PhysicsWorld* world) {
    if (!world || !world->collisions) return;
    
    world->collisions->count = 0;
    
    
    // Sphere-Sphere collisions (A=first sphere, B=second sphere)
    for (int i = 0; i < world->spheres->count; i++) {
        for (int j = i + 1; j < world->spheres->count; j++) {
            if (world->collisions->count >= world->collisions->capacity) break;
            if (world->spheres->body_indices[i] == world->spheres->body_indices[j]) continue;
            
            CollisionContact* contact = &world->collisions->contacts[world->collisions->count];
            if (CheckSphereCollisionWithData(
                world->spheres->positions[i], world->spheres->radii[i],
                world->spheres->positions[j], world->spheres->radii[j], contact)) {
                
                contact->bodyA_type = 0; contact->bodyA_index = world->spheres->body_indices[i]; contact->colliderA_index = i;
                contact->bodyB_type = 0; contact->bodyB_index = world->spheres->body_indices[j]; contact->colliderB_index = j;
                world->collisions->count++;
                
                // Fire callbacks for sphere-sphere collision
                if (world->spheres->collision_callbacks[i]) {
                    world->spheres->collision_callbacks[i](i, 0, j, contact);
                }
                if (world->spheres->collision_callbacks[j]) {
                    world->spheres->collision_callbacks[j](j, 0, i, contact);
                }
            }
        }
    }
    
    // Box-Box collisions using contact manifolds (A=first box, B=second box)
    for (int i = 0; i < world->boxes->count; i++) {
        for (int j = i + 1; j < world->boxes->count; j++) {
            if (world->boxes->body_indices[i] == world->boxes->body_indices[j]) continue;
            // Generate contact manifold for this box pair
            ContactManifold manifold;
            GenerateBoxBoxManifold(world, i, j, &manifold);
            
            // Add manifold contact points to collision collection
            if (manifold.point_count > 0) {
                for (int k = 0; k < manifold.point_count; k++) {
                    if (world->collisions->count >= world->collisions->capacity) break;
                    
                    CollisionContact* contact = &world->collisions->contacts[world->collisions->count];
                    contact->bodyA_type = manifold.bodyA_type;
                    contact->bodyA_index = manifold.bodyA_index;
                    contact->colliderA_index = manifold.colliderA_index;
                    contact->bodyB_type = manifold.bodyB_type;
                    contact->bodyB_index = manifold.bodyB_index;
                    contact->colliderB_index = manifold.colliderB_index;
                    contact->contact_point = manifold.points[k];
                    contact->contact_normal = manifold.normal;  // Consistent normal
                    contact->penetration = manifold.penetrations[k];
                    
                    world->collisions->count++;
                }
                
                // Fire callbacks for box-box collision (once per collision pair)
                if (world->boxes->collision_callbacks[i]) {
                    // Use the first contact point for the callback
                    CollisionContact* first_contact = &world->collisions->contacts[world->collisions->count - manifold.point_count];
                    world->boxes->collision_callbacks[i](i, 1, j, first_contact);
                }
                if (world->boxes->collision_callbacks[j]) {
                    // Use the first contact point for the callback
                    CollisionContact* first_contact = &world->collisions->contacts[world->collisions->count - manifold.point_count];
                    world->boxes->collision_callbacks[j](j, 1, i, first_contact);
                }
            }
        }
    }
    
    // Sphere-Box collisions (A=sphere, B=box)
    for (int i = 0; i < world->spheres->count; i++) {
        for (int j = 0; j < world->boxes->count; j++) {
            if (world->collisions->count >= world->collisions->capacity) break;
            if (world->spheres->body_indices[i] == world->boxes->body_indices[j]) continue;
            
            CollisionContact* contact = &world->collisions->contacts[world->collisions->count];
            if (CheckSphereBoxCollisionWithData(
                world->spheres->positions[i], world->spheres->radii[i],
                world->boxes->positions[j], world->boxes->sizes[j], world->boxes->rotations[j], contact)) {
                
                contact->bodyA_type = 0; contact->bodyA_index = world->spheres->body_indices[i]; contact->colliderA_index = i;  // A = sphere
                contact->bodyB_type = 1; contact->bodyB_index = world->boxes->body_indices[j]; contact->colliderB_index = j;  // B = box
                world->collisions->count++;
                
                // Fire callbacks for sphere-box collision
                if (world->spheres->collision_callbacks[i]) {
                    world->spheres->collision_callbacks[i](i, 1, j, contact);
                }
                if (world->boxes->collision_callbacks[j]) {
                    world->boxes->collision_callbacks[j](j, 0, i, contact);
                }
            }
        }
    }
    
    // Sphere-SDF collisions (A=sphere, B=SDF)
    if (world->usesSDF && world->SDFCollider) {
        for (int i = 0; i < world->spheres->count; i++) {
            if (world->collisions->count >= world->collisions->capacity) break;
            
            Vec3 sphere_center = world->spheres->positions[i];
            float sphere_radius = world->spheres->radii[i];
            
            // Evaluate SDF at sphere center using user function
            float d = world->SDFCollider(sphere_center);
            
            // Check if sphere intersects SDF (d < radius means collision)
            if (d < sphere_radius) {
                CollisionContact* contact = &world->collisions->contacts[world->collisions->count];
                
                // Calculate collision data
                contact->penetration =  sphere_radius-d;
                // SDF gradient points outward from surface, but we need normal from sphere to surface (A to B)
                Vec3 sdf_gradient = SDF_EstimateNormal(world, sphere_center);
                contact->contact_normal = Vec3Scale(sdf_gradient, -1.0f); // Negate to point from sphere to surface
                contact->contact_point = Vec3Sub(sphere_center, Vec3Scale(contact->contact_normal, sphere_radius));
                
                contact->bodyA_type = 0; contact->bodyA_index = world->spheres->body_indices[i]; contact->colliderA_index = i;  // A = sphere
                contact->bodyB_type = 2; contact->bodyB_index = -1; contact->colliderB_index = -1; // B = SDF
                world->collisions->count++;
                
                // Fire callback for sphere-SDF collision
                if (world->spheres->collision_callbacks[i]) {
                    world->spheres->collision_callbacks[i](i, 2, 0, contact);
                }
            }
        }
    }
    
    // Box-SDF collisions with MULTIPLE contact points per box
    if (world->usesSDF && world->SDFCollider) {
        for (int i = 0; i < world->boxes->count; i++) {
            // Generate multiple contact points for each box-SDF collision
            GenerateBoxSDFContacts(world, i);
        }
    }
}

// ============================================================================
// COLLISION RESOLUTION SYSTEM - IMPULSE-BASED
// ============================================================================

// Get velocity at a point on a body (v + ω × r)
Vec3 GetBodyVelocityAtPoint(PhysicsWorld* world, int bodyType, int bodyIndex, Vec3 point) {
    if (world && world->bodies) {
        return GetConstraintBodyVelocityAtPoint(world, bodyType, bodyIndex, point);
    }
    
    Vec3 bodyPos, bodyVel, bodyAngVel;
    //float bodyAngVelLen = Vec3Length(bodyAngVel);
    //if( bodyAngVelLen*bodyAngVelLen < 1e-10f) {  // Make a helper function for square of vec.
    //    return bodyVel;
    //}
    
    if (bodyType == 0) { // Sphere
        bodyPos = world->spheres->positions[bodyIndex];
        bodyVel = world->spheres->velocities[bodyIndex];
        bodyAngVel = world->spheres->angular_velocities[bodyIndex];
    } else if (bodyType == 1) { // Box
        bodyPos = world->boxes->positions[bodyIndex];
        bodyVel = world->boxes->velocities[bodyIndex];
        bodyAngVel = world->boxes->angular_velocities[bodyIndex];
    } else { // SDF (type 2) - static, no velocity
        bodyPos = (Vec3){ 0.0f, 0.0f, 0.0f }; // Not used for static objects
        bodyVel = (Vec3){ 0.0f, 0.0f, 0.0f };
        bodyAngVel = (Vec3){ 0.0f, 0.0f, 0.0f };
    }
    
    // Calculate r = contact_point - body_center
    Vec3 r = Vec3Sub(point, bodyPos);
    
    // Velocity at point = linear_velocity + angular_velocity × r
    Vec3 angularContribution = Vec3Cross(bodyAngVel, r);
    return Vec3Add(bodyVel, angularContribution);
}

// Calculate effective mass for impulse resolution
float CalculateEffectiveMass(PhysicsWorld* world, CollisionContact* contact) {
    if (world && world->bodies) {
        float inv_m1 = GetConstraintBodyInverseMass(world, contact->bodyA_type, contact->bodyA_index);
        float inv_m2 = GetConstraintBodyInverseMass(world, contact->bodyB_type, contact->bodyB_index);
        return (inv_m1 + inv_m2 > 0.0f) ? 1.0f / (inv_m1 + inv_m2) : 0.0f;
    }
    
    Vec3 n = contact->contact_normal;
    Vec3 r1 = Vec3Sub(contact->contact_point, 
                     (contact->bodyA_type == 0) ? world->spheres->positions[contact->bodyA_index] 
                   : (contact->bodyA_type == 1) ? world->boxes->positions[contact->bodyA_index]
                                               : (Vec3){ 0.0f, 0.0f, 0.0f }); // SDF position not used
    Vec3 r2 = Vec3Sub(contact->contact_point, 
                     (contact->bodyB_type == 0) ? world->spheres->positions[contact->bodyB_index] 
                   : (contact->bodyB_type == 1) ? world->boxes->positions[contact->bodyB_index]
                                               : (Vec3){ 0.0f, 0.0f, 0.0f }); // SDF position not used
    
    // Check for static objects and handle infinite mass
    bool is_static1 = (contact->bodyA_type == 0) ? world->spheres->is_static[contact->bodyA_index] 
                    : (contact->bodyA_type == 1) ? world->boxes->is_static[contact->bodyA_index]
                                                 : true; // SDF is always static
    bool is_static2 = (contact->bodyB_type == 0) ? world->spheres->is_static[contact->bodyB_index] 
                    : (contact->bodyB_type == 1) ? world->boxes->is_static[contact->bodyB_index]
                                                 : true; // SDF is always static
    
    // Get masses (use very large mass for static objects)
    //float m1 = is_static1 ? 1e10f : ((contact->bodyA_type == 0) ? world->spheres->masses[contact->bodyA_index] 
    //                                                            : world->boxes->masses[contact->bodyA_index]);
    //float m2 = is_static2 ? 1e10f : ((contact->bodyB_type == 0) ? world->spheres->masses[contact->bodyB_index] 
    //
    float inv_m1 = is_static1 ? 0.0f : 
        1.0f / ((contact->bodyA_type == 0) ? world->spheres->masses[contact->bodyA_index] 
              : (contact->bodyA_type == 1) ? world->boxes->masses[contact->bodyA_index]
                                           : 1.0f); // SDF fallback (not used since static)

    float inv_m2 = is_static2 ? 0.0f : 
        1.0f / ((contact->bodyB_type == 0) ? world->spheres->masses[contact->bodyB_index] 
              : (contact->bodyB_type == 1) ? world->boxes->masses[contact->bodyB_index]
                                           : 1.0f); // SDF fallback (not used since static)
                                                         
    
    // Calculate inertia terms (skip for static objects)
    float inertia_term1 = 0.0f;
    float inertia_term2 = 0.0f;
    
    if (!is_static1) {
        if (contact->bodyA_type == 0) { // Sphere A
            float I1 = world->spheres->inertias[contact->bodyA_index];
            Vec3 r1_cross_n = Vec3Cross(r1, n);
            inertia_term1 = Vec3Dot(r1_cross_n, r1_cross_n) / I1;
        } else { // Box A  
            // Transform r×n to local space for proper inertia calculation
            Quat boxRotation = world->boxes->rotations[contact->bodyA_index];
            Vec3 r1_cross_n_world = Vec3Cross(r1, n);
            Vec3 r1_cross_n_local = QuatRotateVec3(QuatConjugate(boxRotation), r1_cross_n_world);
            
            // Apply with diagonal inertia in local space: (r×n)_local · I⁻¹_local · (r×n)_local
            Vec3 I1 = world->boxes->inertias[contact->bodyA_index];
            inertia_term1 = (r1_cross_n_local.x * r1_cross_n_local.x) / I1.x + 
                           (r1_cross_n_local.y * r1_cross_n_local.y) / I1.y + 
                           (r1_cross_n_local.z * r1_cross_n_local.z) / I1.z;
        }
    }
    
    if (!is_static2) {
        if (contact->bodyB_type == 0) { // Sphere B
            float I2 = world->spheres->inertias[contact->bodyB_index];
            Vec3 r2_cross_n = Vec3Cross(r2, n);
            inertia_term2 = Vec3Dot(r2_cross_n, r2_cross_n) / I2;
        } else { // Box B
            // Transform r×n to local space for proper inertia calculation
            Quat boxRotation = world->boxes->rotations[contact->bodyB_index];
            Vec3 r2_cross_n_world = Vec3Cross(r2, n);
            Vec3 r2_cross_n_local = QuatRotateVec3(QuatConjugate(boxRotation), r2_cross_n_world);
            
            // Apply with diagonal inertia in local space: (r×n)_local · I⁻¹_local · (r×n)_local
            Vec3 I2 = world->boxes->inertias[contact->bodyB_index];
            inertia_term2 = (r2_cross_n_local.x * r2_cross_n_local.x) / I2.x + 
                           (r2_cross_n_local.y * r2_cross_n_local.y) / I2.y + 
                           (r2_cross_n_local.z * r2_cross_n_local.z) / I2.z;
        }
    }
    
    // Effective mass = 1 / (1/m1 + 1/m2 + inertia_terms)
    //return 1.0f / (1.0f/m1 + 1.0f/m2 + inertia_term1 + inertia_term2);
    return 1.0f / (inv_m1 + inv_m2 + inertia_term1 + inertia_term2);
}

// Apply impulse to a body (both linear and angular components)
void ApplyImpulse(PhysicsWorld* world, int bodyType, int bodyIndex, Vec3 impulse, Vec3 contactPoint) {
    if (world && world->bodies) {
        Vec3 bodyPos = (bodyType == 2 || bodyIndex < 0) ? contactPoint : world->bodies->positions[bodyIndex];
        Vec3 r = Vec3Sub(contactPoint, bodyPos);
        ApplyConstraintImpulse(world, bodyType, bodyIndex, impulse, Vec3Cross(r, impulse));
        SyncCollidersFromBodies(world);
        return;
    }
    
    if (bodyType == 0) { // Sphere
        // Skip static objects
        if (world->spheres->is_static[bodyIndex]) return;
        
        // Wake up sleeping objects when impulse is applied
        if (world->spheres->is_sleeping[bodyIndex]) {
            world->spheres->is_sleeping[bodyIndex] = false;
            world->spheres->sleep_timers[bodyIndex] = 0.0f;
        }
        
        // Apply linear impulse: Δv = J / m
        float inv_mass = 1.0f / world->spheres->masses[bodyIndex];
        Vec3 deltaVel = Vec3Scale(impulse, inv_mass);
        world->spheres->velocities[bodyIndex] = Vec3Add(world->spheres->velocities[bodyIndex], deltaVel);
        
        // Apply angular impulse: Δω = (r × J) / I
        Vec3 bodyPos = world->spheres->positions[bodyIndex];
        Vec3 r = Vec3Sub(contactPoint, bodyPos);
        Vec3 torque = Vec3Cross(r, impulse);
        float inv_inertia = 1.0f / world->spheres->inertias[bodyIndex];
        Vec3 deltaAngVel = Vec3Scale(torque, inv_inertia);
        world->spheres->angular_velocities[bodyIndex] = Vec3Add(world->spheres->angular_velocities[bodyIndex], deltaAngVel);
        
    } else if (bodyType == 1) { // Box
        // Skip static objects
        if (world->boxes->is_static[bodyIndex]) return;
        
        // Wake up sleeping objects when impulse is applied
        if (world->boxes->is_sleeping[bodyIndex]) {
            world->boxes->is_sleeping[bodyIndex] = false;
            world->boxes->sleep_timers[bodyIndex] = 0.0f;
        }
        
        // Apply linear impulse: Δv = J / m
        float inv_mass = 1.0f / world->boxes->masses[bodyIndex];
        Vec3 deltaVel = Vec3Scale(impulse, inv_mass);
        world->boxes->velocities[bodyIndex] = Vec3Add(world->boxes->velocities[bodyIndex], deltaVel);
        
        // Apply angular impulse with local space transformation: Δω = R * I⁻¹_local * R^T * (r × J)
        Vec3 bodyPos = world->boxes->positions[bodyIndex];
        Vec3 r = Vec3Sub(contactPoint, bodyPos);
        Vec3 torque = Vec3Cross(r, impulse);
        
        // Get box rotation and transform torque to local space
        Quat boxRotation = world->boxes->rotations[bodyIndex];
        Vec3 torqueLocal = QuatRotateVec3(QuatConjugate(boxRotation), torque);
        
        // Apply with diagonal inertia in local space
        Vec3 inertia = world->boxes->inertias[bodyIndex];
        Vec3 deltaAngVelLocal = {
            torqueLocal.x / inertia.x,
            torqueLocal.y / inertia.y, 
            torqueLocal.z / inertia.z
        };
        
        // Transform angular velocity change back to world space
        Vec3 deltaAngVel = QuatRotateVec3(boxRotation, deltaAngVelLocal);
        world->boxes->angular_velocities[bodyIndex] = Vec3Add(world->boxes->angular_velocities[bodyIndex], deltaAngVel);
    } else { // SDF (type 2) - static, no impulse applied
        return;
    }
}

// OLD FUNCTION REMOVED - Using Baumgarte stabilization in constraint solver instead


// OLD FUNCTION REMOVED - Using constraint-based solver instead


void CleanupPhysics(PhysicsWorld* world) {
    if (!world) return;
    
    // Clear collision list for next frame
    if (world->collisions) {
        world->collisions->count = 0;
    }
    
    // Forces are already cleared in IntegrateMotion()
}


// Estimate normal at a point using finite differences
Vec3 SDF_EstimateNormal(PhysicsWorld* world, Vec3 point) {
    if (world->SDFNormal) {
        Vec3 normal = Vec3Normalize(world->SDFNormal(point));
        if (Vec3Length(normal) > world->settings.vector_normalize_epsilon) {
            return normal;
        }
    }

    float epsilon = world->settings.sdf_normal_epsilon;
    Vec3 eps_x = { epsilon, 0.0f, 0.0f };
    Vec3 eps_y = { 0.0f, epsilon, 0.0f };
    Vec3 eps_z = { 0.0f, 0.0f, epsilon };
    
    // Central difference for more accurate gradient
    float nx = (world->SDFCollider(Vec3Add(point, eps_x)) - world->SDFCollider(Vec3Sub(point, eps_x))) / (2.0f * epsilon);
    float ny = (world->SDFCollider(Vec3Add(point, eps_y)) - world->SDFCollider(Vec3Sub(point, eps_y))) / (2.0f * epsilon);
    float nz = (world->SDFCollider(Vec3Add(point, eps_z)) - world->SDFCollider(Vec3Sub(point, eps_z))) / (2.0f * epsilon);
    
    Vec3 normal = { nx, ny, nz };
    return Vec3Normalize(normal);
}

typedef struct {
    Vec3 point;
    float value;
    float penetration;
} SDFContactCandidate;

static int ClampPhysicsInt(int value, int min_value, int max_value) {
    if (value < min_value) return min_value;
    if (value > max_value) return max_value;
    return value;
}

static int BuildBoxSDFTraceSamples(Vec3 boxPos, Vec3 boxSize, Quat boxRot, 
                                   PhysicsSDFBoxTraceSamples sample_mode, Vec3* samples, int max_samples) {
    Vec3 axes[3] = {
        QuatRotateVec3(boxRot, (Vec3){1.0f, 0.0f, 0.0f}),
        QuatRotateVec3(boxRot, (Vec3){0.0f, 1.0f, 0.0f}),
        QuatRotateVec3(boxRot, (Vec3){0.0f, 0.0f, 1.0f})
    };
    float extents[3] = { boxSize.x, boxSize.y, boxSize.z };
    int count = 0;
    
    for (int ix = -1; ix <= 1; ix++) {
        for (int iy = -1; iy <= 1; iy++) {
            for (int iz = -1; iz <= 1; iz++) {
                int non_zero_axes = (ix != 0) + (iy != 0) + (iz != 0);
                if (non_zero_axes == 0) continue;
                if (sample_mode == PHYSICS_SDF_BOX_TRACE_CORNERS && non_zero_axes != 3) continue;
                if (sample_mode == PHYSICS_SDF_BOX_TRACE_CORNERS_AND_FACES && non_zero_axes == 2) continue;
                if (count >= max_samples) return count;
                
                Vec3 p = boxPos;
                p = Vec3Add(p, Vec3Scale(axes[0], (float)ix * extents[0]));
                p = Vec3Add(p, Vec3Scale(axes[1], (float)iy * extents[1]));
                p = Vec3Add(p, Vec3Scale(axes[2], (float)iz * extents[2]));
                samples[count++] = p;
            }
        }
    }
    
    return count;
}

static bool IsDuplicateSDFContact(SDFContactCandidate* candidates, int candidate_count, Vec3 point, float min_separation) {
    float min_separation_sq = min_separation * min_separation;
    for (int i = 0; i < candidate_count; i++) {
        Vec3 delta = Vec3Sub(candidates[i].point, point);
        if (Vec3Dot(delta, delta) <= min_separation_sq) {
            return true;
        }
    }
    return false;
}

static void InsertSDFContactCandidate(SDFContactCandidate* candidates, int* candidate_count, int max_candidates,
                                      Vec3 point, float value, float penetration, float min_separation) {
    if (!candidates || !candidate_count || max_candidates <= 0 || penetration <= 0.0f) return;
    if (IsDuplicateSDFContact(candidates, *candidate_count, point, min_separation)) return;
    
    int insert_index = *candidate_count;
    if (*candidate_count < max_candidates) {
        (*candidate_count)++;
    } else if (penetration > candidates[max_candidates - 1].penetration) {
        insert_index = max_candidates - 1;
    } else {
        return;
    }
    
    while (insert_index > 0 && candidates[insert_index - 1].penetration < penetration) {
        if (insert_index < max_candidates) {
            candidates[insert_index] = candidates[insert_index - 1];
        }
        insert_index--;
    }
    
    candidates[insert_index].point = point;
    candidates[insert_index].value = value;
    candidates[insert_index].penetration = penetration;
}

static Vec3 RefineSDFSignCrossing(PhysicsWorld* world, Vec3 a, Vec3 b, float value_a, int refine_steps) {
    Vec3 low = a;
    Vec3 high = b;
    float low_value = value_a;
    
    for (int i = 0; i < refine_steps; i++) {
        Vec3 mid = Vec3Scale(Vec3Add(low, high), 0.5f);
        float mid_value = world->SDFCollider(mid);
        
        if ((low_value <= 0.0f && mid_value <= 0.0f) || (low_value >= 0.0f && mid_value >= 0.0f)) {
            low = mid;
            low_value = mid_value;
        } else {
            high = mid;
        }
    }
    
    return Vec3Scale(Vec3Add(low, high), 0.5f);
}

static bool TraceSDFSurfaceOnSegment(PhysicsWorld* world, Vec3 start, float start_value, Vec3 end,
                                     SDFContactCandidate* out_candidate) {
    int trace_steps = ClampPhysicsInt(world->settings.sdf_box_trace_steps, 1, 128);
    int refine_steps = ClampPhysicsInt(world->settings.sdf_box_trace_refine_steps, 0, 32);
    float surface_epsilon = fmaxf(world->settings.sdf_box_surface_epsilon, 0.0f);
    float contact_margin = fmaxf(surface_epsilon, world->settings.allowed_penetration);
    
    Vec3 segment = Vec3Sub(end, start);
    Vec3 previous_point = start;
    float previous_value = start_value;
    Vec3 best_point = start;
    float best_value = start_value;
    float best_abs_value = fabsf(start_value);
    bool found_sign_crossing = false;
    
    for (int step = 1; step <= trace_steps; step++) {
        float t = (float)step / (float)trace_steps;
        Vec3 p = Vec3Add(start, Vec3Scale(segment, t));
        float value = world->SDFCollider(p);
        float abs_value = fabsf(value);
        
        if (abs_value < best_abs_value) {
            best_abs_value = abs_value;
            best_value = value;
            best_point = p;
        }
        
        if ((previous_value <= 0.0f && value >= 0.0f) || (previous_value >= 0.0f && value <= 0.0f)) {
            if (fabsf(previous_value - value) > world->settings.collision_epsilon) {
                best_point = RefineSDFSignCrossing(world, previous_point, p, previous_value, refine_steps);
                best_value = world->SDFCollider(best_point);
                best_abs_value = fabsf(best_value);
                found_sign_crossing = true;
                break;
            }
        }
        
        previous_point = p;
        previous_value = value;
    }
    
    if (!found_sign_crossing && best_value > surface_epsilon && best_abs_value > surface_epsilon) {
        return false;
    }
    
    if (out_candidate) {
        out_candidate->point = best_point;
        out_candidate->value = best_value;
        out_candidate->penetration = contact_margin - fminf(best_value, best_abs_value);
        if (out_candidate->penetration <= 0.0f) {
            out_candidate->penetration = contact_margin;
        }
    }
    
    return true;
}

static int FindBoxSDFSurfaceTraceContacts(PhysicsWorld* world, Vec3 boxPos, Vec3 boxSize, Quat boxRot,
                                          SDFContactCandidate* candidates, int max_candidates) {
    if (!world || !world->SDFCollider || !candidates || max_candidates <= 0) return 0;
    
    Vec3 sample_points[26];
    PhysicsSDFBoxTraceSamples sample_mode = world->settings.sdf_box_trace_samples;
    if (sample_mode < PHYSICS_SDF_BOX_TRACE_CORNERS || sample_mode > PHYSICS_SDF_BOX_TRACE_FULL) {
        sample_mode = PHYSICS_SDF_BOX_TRACE_CORNERS;
    }
    
    int sample_count = BuildBoxSDFTraceSamples(boxPos, boxSize, boxRot, sample_mode, sample_points, 26);
    float center_value = world->SDFCollider(boxPos);
    float surface_epsilon = fmaxf(world->settings.sdf_box_surface_epsilon, 0.0f);
    float contact_margin = fmaxf(surface_epsilon, world->settings.allowed_penetration);
    float min_separation = fmaxf(contact_margin * 0.5f, world->settings.collision_epsilon);
    int candidate_count = 0;
    
    for (int i = 0; i < sample_count; i++) {
        SDFContactCandidate candidate;
        if (!TraceSDFSurfaceOnSegment(world, boxPos, center_value, sample_points[i], &candidate)) {
            continue;
        }
        
        InsertSDFContactCandidate(candidates, &candidate_count, max_candidates,
                                  candidate.point, candidate.value, candidate.penetration, min_separation);
    }
    
    return candidate_count;
}

static int FindBoxSDFSignedDistanceContacts(PhysicsWorld* world, Vec3 boxPos, Vec3 boxSize, Quat boxRot,
                                            SDFContactCandidate* candidates, int max_candidates) {
    if (!world || !world->SDFCollider || !candidates || max_candidates <= 0) return 0;
    
    float circumscribed_radius = Vec3Length(boxSize);
    float center_distance = world->SDFCollider(boxPos);
    
    if (center_distance > circumscribed_radius) {
        return 0;
    }
    
    Vec3 sample_points[26];
    int sample_count = BuildBoxSDFTraceSamples(boxPos, boxSize, boxRot, PHYSICS_SDF_BOX_TRACE_FULL, sample_points, 26);
    float contact_threshold = -fmaxf(world->settings.collision_epsilon, 0.001f);
    float min_separation = fmaxf(world->settings.allowed_penetration * 0.5f, world->settings.collision_epsilon);
    int candidate_count = 0;
    
    for (int i = 0; i < sample_count; i++) {
        float value = world->SDFCollider(sample_points[i]);
        if (value < contact_threshold) {
            InsertSDFContactCandidate(candidates, &candidate_count, max_candidates,
                                      sample_points[i], value, -value, min_separation);
        }
    }
    
    return candidate_count;
}

static Vec3 GetBoxSDFContactNormal(PhysicsWorld* world, Vec3 boxPos, Vec3 contactPoint) {
    Vec3 sdf_gradient = SDF_EstimateNormal(world, contactPoint);
    Vec3 normal = Vec3Scale(sdf_gradient, -1.0f);
    Vec3 center_to_contact = Vec3Sub(contactPoint, boxPos);
    
    if (Vec3Dot(normal, center_to_contact) < 0.0f) {
        normal = Vec3Scale(normal, -1.0f);
    }
    
    return normal;
}

// Box-SDF collision detection using either signed-distance sampling or surface tracing
bool CheckBoxSDFCollision(PhysicsWorld* world, Vec3 boxPos, Vec3 boxSize, Quat boxRot, CollisionContact* contact) {
    if (!world || !world->SDFCollider) return false;
    
    SDFContactCandidate candidates[26];
    int max_contacts = ClampPhysicsInt(world->settings.sdf_box_max_contacts, 1, 26);
    int candidate_count = 0;
    
    if (world->settings.sdf_box_collision_mode == PHYSICS_SDF_BOX_SIGNED_DISTANCE) {
        candidate_count = FindBoxSDFSignedDistanceContacts(world, boxPos, boxSize, boxRot, candidates, max_contacts);
    } else {
        candidate_count = FindBoxSDFSurfaceTraceContacts(world, boxPos, boxSize, boxRot, candidates, max_contacts);
    }
    
    if (candidate_count <= 0) return false;
    
    if (contact) {
        contact->penetration = candidates[0].penetration;
        contact->contact_normal = GetBoxSDFContactNormal(world, boxPos, candidates[0].point);
        contact->contact_point = candidates[0].point;
    }
    
    return true;
}

// Generate multiple contact points for box-SDF collision (for stable stacking)
void GenerateBoxSDFContacts(PhysicsWorld* world, int boxIndex) {
    if (!world || !world->usesSDF || !world->SDFCollider) return;
    if (boxIndex < 0 || boxIndex >= world->boxes->count) return;
    
    Vec3 boxPos = world->boxes->positions[boxIndex];
    Vec3 boxSize = world->boxes->sizes[boxIndex];
    Quat boxRot = world->boxes->rotations[boxIndex];
    
    SDFContactCandidate candidates[26];
    int max_contacts = ClampPhysicsInt(world->settings.sdf_box_max_contacts, 1, 26);
    int candidate_count = 0;
    
    if (world->settings.sdf_box_collision_mode == PHYSICS_SDF_BOX_SIGNED_DISTANCE) {
        candidate_count = FindBoxSDFSignedDistanceContacts(world, boxPos, boxSize, boxRot, candidates, max_contacts);
    } else {
        candidate_count = FindBoxSDFSurfaceTraceContacts(world, boxPos, boxSize, boxRot, candidates, max_contacts);
    }
    
    for (int i = 0; i < candidate_count; i++) {
        if (world->collisions->count >= world->collisions->capacity) break;
        
        CollisionContact* contact = &world->collisions->contacts[world->collisions->count];
        contact->penetration = candidates[i].penetration;
        contact->contact_normal = GetBoxSDFContactNormal(world, boxPos, candidates[i].point);
        contact->contact_point = candidates[i].point;
        contact->bodyA_type = 1; contact->bodyA_index = world->boxes->body_indices[boxIndex]; contact->colliderA_index = boxIndex;
        contact->bodyB_type = 2; contact->bodyB_index = -1; contact->colliderB_index = -1;
        world->collisions->count++;
    }
}

// ============================================================================
// PHYSICS SETTINGS CONVENIENCE FUNCTIONS
// ============================================================================

void ResetPhysicsSettingsToDefaults(PhysicsWorld* world) {
    if (!world) return;
    
    // Damping (working values from testing)
    world->settings.linear_damping = 0.98f;
    world->settings.angular_damping = 0.98f;
    
    // Collision response (working values from testing)
    world->settings.restitution = 0.0f;      // No bounciness for stable stacking
    world->settings.friction = 0.1f;         // Lower friction works better with manifolds
    world->settings.solver_iterations = 25;  // More iterations for accuracy
    
    // Constraint solver settings (Baumgarte stabilization) 
    world->settings.baumgarte_bias = 0.3f;       // High Baumgarte for strong position correction
    world->settings.allowed_penetration = 0.05f; // Larger allowed penetration
    world->settings.velocity_threshold = 0.02f;   // Lower threshold for restitution
    
    // Sleep system
    world->settings.sleep_linear_threshold = 0.1f;
    world->settings.sleep_angular_threshold = 0.1f;
    world->settings.sleep_time_required = 0.05f;   // Sleep faster for stability
    
    // Contact manifold settings (working values from testing)
    world->settings.manifold_contact_tolerance = 0.01f;
    world->settings.manifold_penetration_tolerance = -0.5f;  // Very permissive for touching
    world->settings.manifold_max_contacts = 8;                // More contacts for stability

    // Box-SDF collision settings. Surface trace handles implicit or unsigned
    // fields while keeping the default to 8 center-to-corner rays.
    world->settings.sdf_box_collision_mode = PHYSICS_SDF_BOX_SURFACE_TRACE;
    world->settings.sdf_box_trace_samples = PHYSICS_SDF_BOX_TRACE_CORNERS;
    world->settings.sdf_box_trace_steps = 6;
    world->settings.sdf_box_trace_refine_steps = 4;
    world->settings.sdf_box_surface_epsilon = 0.05f;
    world->settings.sdf_box_max_contacts = 4;
    
    // Numerical tolerances (rarely changed)
    world->settings.collision_epsilon = 1e-6f;
    world->settings.sdf_normal_epsilon = 0.001f;
    world->settings.vector_normalize_epsilon = 1e-6f;
    world->settings.quaternion_epsilon = 1e-6f;
}

void SetDamping(PhysicsWorld* world, float linear, float angular) {
    if (!world) return;
    world->settings.linear_damping = linear;
    world->settings.angular_damping = angular;
}

void SetRestitution(PhysicsWorld* world, float restitution) {
    if (!world) return;
    world->settings.restitution = restitution;
}

void SetCorrectionParameters(PhysicsWorld* world, float baumgarte_bias, float allowed_penetration) {
    if (!world) return;
    world->settings.baumgarte_bias = baumgarte_bias;
    world->settings.allowed_penetration = allowed_penetration;
}

void SetSolverIterations(PhysicsWorld* world, int iterations) {
    if (!world) return;
    world->settings.solver_iterations = iterations;
}

void SetManifoldSettings(PhysicsWorld* world, float contact_tolerance, float penetration_tolerance, int max_contacts) {
    if (!world) return;
    world->settings.manifold_contact_tolerance = contact_tolerance;
    world->settings.manifold_penetration_tolerance = penetration_tolerance;
    world->settings.manifold_max_contacts = (max_contacts > MAX_MANIFOLD_CONTACTS) ? MAX_MANIFOLD_CONTACTS : max_contacts;
}

void SetSDFBoxCollisionMode(PhysicsWorld* world, PhysicsSDFBoxCollisionMode mode) {
    if (!world) return;
    if (mode < PHYSICS_SDF_BOX_SIGNED_DISTANCE || mode > PHYSICS_SDF_BOX_SURFACE_TRACE) return;
    
    world->settings.sdf_box_collision_mode = mode;
}

void SetSDFSurfaceTraceSettings(PhysicsWorld* world, PhysicsSDFBoxTraceSamples samples, int trace_steps, 
                                int refine_steps, float surface_epsilon, int max_contacts) {
    if (!world) return;
    
    if (samples >= PHYSICS_SDF_BOX_TRACE_CORNERS && samples <= PHYSICS_SDF_BOX_TRACE_FULL) {
        world->settings.sdf_box_trace_samples = samples;
    }
    
    world->settings.sdf_box_trace_steps = ClampPhysicsInt(trace_steps, 1, 128);
    world->settings.sdf_box_trace_refine_steps = ClampPhysicsInt(refine_steps, 0, 32);
    world->settings.sdf_box_surface_epsilon = fmaxf(surface_epsilon, 0.0f);
    world->settings.sdf_box_max_contacts = ClampPhysicsInt(max_contacts, 1, 26);
}

// ================================
// Contact Manifold System
// ================================

void GenerateBoxBoxManifold(PhysicsWorld* world, int boxIndexA, int boxIndexB, ContactManifold* manifold) {
    if (!world || !manifold) return;
    if (boxIndexA < 0 || boxIndexA >= world->boxes->count) return;
    if (boxIndexB < 0 || boxIndexB >= world->boxes->count) return;
    if (boxIndexA == boxIndexB) return;
    
    // Initialize manifold
    manifold->point_count = 0;
    manifold->bodyA_type = 1; // box
    manifold->bodyA_index = world->boxes->body_indices[boxIndexA];
    manifold->colliderA_index = boxIndexA;
    manifold->bodyB_type = 1; // box
    manifold->bodyB_index = world->boxes->body_indices[boxIndexB];
    manifold->colliderB_index = boxIndexB;
    
    Vec3 posA = world->boxes->positions[boxIndexA];
    Vec3 sizeA = world->boxes->sizes[boxIndexA];
    Quat rotA = world->boxes->rotations[boxIndexA];
    
    Vec3 posB = world->boxes->positions[boxIndexB];
    Vec3 sizeB = world->boxes->sizes[boxIndexB];
    Quat rotB = world->boxes->rotations[boxIndexB];
    
    // First check if boxes are colliding at all using existing SAT test
    CollisionContact basicContact;
    if (!CheckBoxCollisionWithData(posA, sizeA, rotA, posB, sizeB, rotB, &basicContact)) {
        return; // No collision
    }
    
    // Collect all candidate contact points from multiple sources
    Vec3 candidatePoints[24];  // Increased capacity: 8 vertices + 8 vertices + potential edge contacts
    float candidatePenetrations[24];
    int candidateCount = 0;
    
    // FIRST: Always include the SAT-detected contact point (this is the most accurate)
    candidatePoints[candidateCount] = basicContact.contact_point;
    candidatePenetrations[candidateCount] = basicContact.penetration;
    candidateCount++;
    
    // Get vertices of both boxes
    Vec3 verticesA[8], verticesB[8];
    OBB obbA = CreateOBB(posA, sizeA, rotA);
    OBB obbB = CreateOBB(posB, sizeB, rotB);
    GetOBBVertices(&obbA, verticesA);
    GetOBBVertices(&obbB, verticesB);
    
    // SECOND: Check vertices of box A that are inside box B
    for (int i = 0; i < 8 && candidateCount < 24; i++) {
        Vec3 vertexA = verticesA[i];
        
        // Transform vertex A to box B's local space
        Vec3 relativePos = Vec3Sub(vertexA, posB);
        Quat invRotB = QuatConjugate(rotB);
        Vec3 localPos = QuatRotateVec3(invRotB, relativePos);
        
        // Check if vertex is inside or touching box B using configurable tolerance
        float tolerance = world->settings.manifold_contact_tolerance;
        if (fabsf(localPos.x) <= sizeB.x + tolerance && 
            fabsf(localPos.y) <= sizeB.y + tolerance && 
            fabsf(localPos.z) <= sizeB.z + tolerance) {
            
            // Calculate penetration depth (how far inside the box)
            float penetrationX = sizeB.x - fabsf(localPos.x);
            float penetrationY = sizeB.y - fabsf(localPos.y);
            float penetrationZ = sizeB.z - fabsf(localPos.z);
            float minPenetration = fminf(penetrationX, fminf(penetrationY, penetrationZ));
            
            // Accept penetrations based on configurable threshold
            if (minPenetration > world->settings.manifold_penetration_tolerance) {
                candidatePoints[candidateCount] = vertexA;
                candidatePenetrations[candidateCount] = minPenetration;
                candidateCount++;
            }
        }
    }
    
    // THIRD: Check vertices of box B that are inside box A
    for (int i = 0; i < 8 && candidateCount < 24; i++) {
        Vec3 vertexB = verticesB[i];
        
        // Transform vertex B to box A's local space
        Vec3 relativePos = Vec3Sub(vertexB, posA);
        Quat invRotA = QuatConjugate(rotA);
        Vec3 localPos = QuatRotateVec3(invRotA, relativePos);
        
        // Check if vertex is inside or touching box A using configurable tolerance
        float tolerance = world->settings.manifold_contact_tolerance;
        if (fabsf(localPos.x) <= sizeA.x + tolerance && 
            fabsf(localPos.y) <= sizeA.y + tolerance && 
            fabsf(localPos.z) <= sizeA.z + tolerance) {
            
            // Calculate penetration depth
            float penetrationX = sizeA.x - fabsf(localPos.x);
            float penetrationY = sizeA.y - fabsf(localPos.y);
            float penetrationZ = sizeA.z - fabsf(localPos.z);
            float minPenetration = fminf(penetrationX, fminf(penetrationY, penetrationZ));
            
            // Accept penetrations based on configurable threshold
            if (minPenetration > world->settings.manifold_penetration_tolerance) {
                candidatePoints[candidateCount] = vertexB;
                candidatePenetrations[candidateCount] = minPenetration;
                candidateCount++;
            }
        }
    }
    
    // FOURTH: If we still don't have enough contact points, check for edge-edge intersections
    // This handles cases like two boxes touching edge-to-edge
    if (candidateCount < 3) {
        // Simple approach: Check if any edges of box A are close to any faces of box B
        Vec3 faceNormals[6] = {
            QuatRotateVec3(rotB, (Vec3){1, 0, 0}),   // +X face normal
            QuatRotateVec3(rotB, (Vec3){-1, 0, 0}),  // -X face normal
            QuatRotateVec3(rotB, (Vec3){0, 1, 0}),   // +Y face normal
            QuatRotateVec3(rotB, (Vec3){0, -1, 0}),  // -Y face normal
            QuatRotateVec3(rotB, (Vec3){0, 0, 1}),   // +Z face normal
            QuatRotateVec3(rotB, (Vec3){0, 0, -1})   // -Z face normal
        };
        
        // Sample points along edges of box A and check distance to box B
        for (int edge = 0; edge < 12 && candidateCount < 24; edge++) {
            // Get edge endpoints (simplified - just checking a few key edges)
            Vec3 edgeStart, edgeEnd;
            if (edge < 4) { // Bottom face edges
                edgeStart = verticesA[edge];
                edgeEnd = verticesA[(edge + 1) % 4];
            } else if (edge < 8) { // Top face edges
                edgeStart = verticesA[edge];
                edgeEnd = verticesA[4 + ((edge - 4 + 1) % 4)];
            } else { // Vertical edges
                edgeStart = verticesA[edge - 8];
                edgeEnd = verticesA[edge - 8 + 4];
            }
            
            // Sample middle point of this edge
            Vec3 edgeMidpoint = Vec3Scale(Vec3Add(edgeStart, edgeEnd), 0.5f);
            
            // Check if this point is close to any face of box B
            Vec3 relativePos = Vec3Sub(edgeMidpoint, posB);
            Quat invRotB = QuatConjugate(rotB);
            Vec3 localPos = QuatRotateVec3(invRotB, relativePos);
            
            // Check if point is close to box B surface
            bool nearFace = false;
            float minDistToSurface = 1000.0f;
            
            if (fabsf(localPos.x) <= sizeB.x + world->settings.manifold_contact_tolerance &&
                fabsf(localPos.y) <= sizeB.y + world->settings.manifold_contact_tolerance &&
                fabsf(localPos.z) <= sizeB.z + world->settings.manifold_contact_tolerance) {
                
                float distX = sizeB.x - fabsf(localPos.x);
                float distY = sizeB.y - fabsf(localPos.y);
                float distZ = sizeB.z - fabsf(localPos.z);
                minDistToSurface = fminf(distX, fminf(distY, distZ));
                
                if (minDistToSurface > world->settings.manifold_penetration_tolerance) {
                    nearFace = true;
                }
            }
            
            if (nearFace) {
                candidatePoints[candidateCount] = edgeMidpoint;
                candidatePenetrations[candidateCount] = minDistToSurface;
                candidateCount++;
            }
        }
    }
    
    if (candidateCount == 0) return;
    
    // Use the basic contact normal as the manifold normal (from deepest penetration)
    manifold->normal = basicContact.contact_normal;
    
    // Select optimal contact points from candidates
    SelectOptimalContactPoints(candidatePoints, candidatePenetrations, candidateCount, 
                              manifold, manifold->normal, world->settings.manifold_max_contacts);
}

void SelectOptimalContactPoints(Vec3* candidate_points, float* penetrations, int candidate_count, 
                               ContactManifold* manifold, Vec3 normal, int max_contacts) {
    if (!candidate_points || !penetrations || !manifold || candidate_count <= 0) return;
    
    if (candidate_count == 1) {
        // Only one contact point
        manifold->points[0] = candidate_points[0];
        manifold->penetrations[0] = penetrations[0];
        manifold->point_count = 1;
        return;
    }
    
    if (candidate_count == 2) {
        // Two contact points - use both
        manifold->points[0] = candidate_points[0];
        manifold->penetrations[0] = penetrations[0];
        manifold->points[1] = candidate_points[1];
        manifold->penetrations[1] = penetrations[1];
        manifold->point_count = 2;
        return;
    }
    
    // For 3+ contact points, select up to 4 optimal points
    // Algorithm: Select points that maximize the contact area
    
    // Step 1: Find the deepest contact point (guaranteed to be included)
    int deepestIndex = 0;
    for (int i = 1; i < candidate_count; i++) {
        if (penetrations[i] > penetrations[deepestIndex]) {
            deepestIndex = i;
        }
    }
    
    manifold->points[0] = candidate_points[deepestIndex];
    manifold->penetrations[0] = penetrations[deepestIndex];
    manifold->point_count = 1;
    
    if (candidate_count <= max_contacts) {
        // If we have few enough candidates, use all of them
        int idx = 1;
        for (int i = 0; i < candidate_count; i++) {
            if (i != deepestIndex) {
                manifold->points[idx] = candidate_points[i];
                manifold->penetrations[idx] = penetrations[i];
                idx++;
            }
        }
        manifold->point_count = candidate_count;
        return;
    }
    
    // For more points than allowed, select (max_contacts-1) more points that maximize area
    // Project points onto contact plane (perpendicular to normal)
    Vec3 tangent1, tangent2;
    
    // Generate orthogonal tangent vectors to the normal
    if (fabsf(normal.x) > 0.1f) {
        tangent1 = Vec3Normalize((Vec3){0, 1, 0});
    } else {
        tangent1 = Vec3Normalize((Vec3){1, 0, 0});
    }
    tangent1 = Vec3Normalize(Vec3Sub(tangent1, Vec3Scale(normal, Vec3Dot(tangent1, normal))));
    tangent2 = Vec3Cross(normal, tangent1);
    
    // Find center of contact points
    Vec3 center = {0, 0, 0};
    for (int i = 0; i < candidate_count; i++) {
        center = Vec3Add(center, candidate_points[i]);
    }
    center = Vec3Scale(center, 1.0f / candidate_count);
    
    // Select points furthest from center in different directions
    // (This creates a good spread for stability)
    int remaining_slots = max_contacts - 1;  // Already have deepest point
    
    // Create array to track best distances and indices
    float* maxDists = (float*)malloc(remaining_slots * sizeof(float));
    int* bestIdxs = (int*)malloc(remaining_slots * sizeof(int));
    
    // Initialize arrays
    for (int j = 0; j < remaining_slots; j++) {
        maxDists[j] = -1.0f;
        bestIdxs[j] = -1;
    }
    
    // Find the furthest points from center
    for (int i = 0; i < candidate_count; i++) {
        if (i == deepestIndex) continue;
        
        float dist = Vec3Length(Vec3Sub(candidate_points[i], center));
        
        // Insert this distance in sorted order
        for (int j = 0; j < remaining_slots; j++) {
            if (dist > maxDists[j]) {
                // Shift lower values down
                for (int k = remaining_slots - 1; k > j; k--) {
                    maxDists[k] = maxDists[k-1];
                    bestIdxs[k] = bestIdxs[k-1];
                }
                maxDists[j] = dist;
                bestIdxs[j] = i;
                break;
            }
        }
    }
    
    // Add the selected points
    for (int j = 0; j < remaining_slots; j++) {
        if (bestIdxs[j] >= 0) {
            manifold->points[manifold->point_count] = candidate_points[bestIdxs[j]];
            manifold->penetrations[manifold->point_count] = penetrations[bestIdxs[j]];
            manifold->point_count++;
        }
    }
    
    free(maxDists);
    free(bestIdxs);
}

void AddManifoldToConstraints(PhysicsWorld* world, ContactManifold* manifold) {
    if (!world || !manifold || manifold->point_count == 0) return;
    
    // Add one constraint for each contact point in the manifold
    for (int i = 0; i < manifold->point_count; i++) {
        if (world->constraints->count >= world->constraints->capacity) break;
        
        ContactConstraint* constraint = &world->constraints->constraints[world->constraints->count];
        
        constraint->bodyA_type = manifold->bodyA_type;
        constraint->bodyA_index = manifold->bodyA_index;
        constraint->colliderA_index = manifold->colliderA_index;
        constraint->bodyB_type = manifold->bodyB_type;
        constraint->bodyB_index = manifold->bodyB_index;
        constraint->colliderB_index = manifold->colliderB_index;
        
        constraint->contact_point = manifold->points[i];
        constraint->contact_normal = manifold->normal;  // Consistent normal for all points
        constraint->penetration = manifold->penetrations[i];
        
        // Use physics settings
        constraint->restitution = world->settings.restitution;
        constraint->friction = world->settings.friction;
        
        world->constraints->count++;
    }
}

#endif // BALLBOX_IMPLEMENTATION

#endif // BALLBOX_H
