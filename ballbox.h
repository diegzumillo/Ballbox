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
    int bodyA_type;     // 0 = sphere, 1 = box
    int bodyA_index;
    int bodyB_type;     // 0 = sphere, 1 = box
    int bodyB_index;
    Vec3 contact_point;
    Vec3 contact_normal; // From A to B
    float penetration;
} CollisionContact;

// Function pointer type for collision callbacks
typedef void (*OnCollision)(int selfIndex, int otherBodyType, int otherIndex, CollisionContact* contact);

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
    int bodyA_type;        // 0 = sphere, 1 = box, 2 = SDF
    int bodyA_index;
    int bodyB_type;        // 0 = sphere, 1 = box, 2 = SDF
    int bodyB_index;
} ContactManifold;

// Contact constraint for Sequential Impulse solver
typedef struct {
    int bodyA_type;     // 0 = sphere, 1 = box, 2 = SDF
    int bodyA_index;
    int bodyB_type;     // 0 = sphere, 1 = box, 2 = SDF  
    int bodyB_index;
    
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

// Function pointer type for user-defined SDF
typedef float (*SDFFunction)(Vec3 point);

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
    
    // Numerical tolerances (rarely changed)
    float collision_epsilon;   // 1e-6f - collision detection threshold
    float sdf_normal_epsilon;  // 0.001f - SDF normal estimation step size
    float vector_normalize_epsilon; // 1e-6f - vector normalization threshold
    float quaternion_epsilon;  // 1e-6f - quaternion operations threshold
} PhysicsSettings;

typedef struct {
    SphereSystem* spheres;
    BoxSystem* boxes;
    CollisionCollection* collisions;
    ConstraintCollection* constraints;  // Sequential Impulse constraints
    Vec3 gravity;
    float dt;
    bool usesSDF;
    SDFFunction SDFCollider;
    PhysicsSettings settings;
} PhysicsWorld;

// Create/destroy physics world
PhysicsWorld* CreatePhysicsWorld(int max_spheres, int max_boxes, int max_collisions);
void DestroyPhysicsWorld(PhysicsWorld* world);

// SDF system
void UseSDF(PhysicsWorld* world, SDFFunction sdf_function);

// Physics settings convenience functions
void SetDamping(PhysicsWorld* world, float linear, float angular);
void SetRestitution(PhysicsWorld* world, float restitution);
void SetCorrectionParameters(PhysicsWorld* world, float baumgarte_bias, float allowed_penetration);
void SetSolverIterations(PhysicsWorld* world, int iterations);
void SetManifoldSettings(PhysicsWorld* world, float contact_tolerance, float penetration_tolerance, int max_contacts);
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

// Sphere system management
SphereSystem* CreateSphereSystem(int capacity);
void DestroySphereSystem(SphereSystem* system);
int AddSphere(PhysicsWorld* world, Vec3 position, float radius);
void SetSphereMass(PhysicsWorld* world, int index, float mass);
void SetSphereStatic(PhysicsWorld* world, int index, bool is_static);
void AddSphereForce(PhysicsWorld* world, int index, Vec3 force);
void AddSphereTorque(PhysicsWorld* world, int index, Vec3 torque);
bool CheckSphereCollision(Vec3 pos1, float radius1, Vec3 pos2, float radius2);
bool CheckSphereCollisionWithData(Vec3 pos1, float radius1, Vec3 pos2, float radius2, CollisionContact* contact);
bool CheckSphereSystemCollisions(SphereSystem* system, int sphere_index);

// Box system management  
BoxSystem* CreateBoxSystem(int capacity);
void DestroyBoxSystem(BoxSystem* system);
int AddBox(PhysicsWorld* world, Vec3 position, Vec3 size, Quat rotation); // size = half-extents
void SetBoxMass(PhysicsWorld* world, int index, float mass);
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
    system->count = 0;
    system->capacity = capacity;
    
    if (!system->positions || !system->radii || !system->velocities || !system->forces || 
        !system->masses || !system->angular_velocities || !system->torques || !system->inertias || 
        !system->is_static || !system->is_sleeping || !system->sleep_timers || !system->collision_callbacks) {
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
    free(system);
}

int AddSphere(PhysicsWorld* world, Vec3 position, float radius) {
    if (!world || !world->spheres) return -1;
    SphereSystem* system = world->spheres;
    if (!system || system->count >= system->capacity) return -1;
    
    int index = system->count;
    system->positions[index] = position;
    system->radii[index] = radius;
    system->velocities[index] = (Vec3){ 0.0f, 0.0f, 0.0f };
    system->forces[index] = (Vec3){ 0.0f, 0.0f, 0.0f };
    system->masses[index] = 1.0f;
    system->angular_velocities[index] = (Vec3){ 0.0f, 0.0f, 0.0f };
    system->torques[index] = (Vec3){ 0.0f, 0.0f, 0.0f };
    system->inertias[index] = (2.0f / 5.0f) * system->masses[index] * radius * radius;
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
    
    float radius = system->radii[index];
    system->inertias[index] = (2.0f / 5.0f) * mass * radius * radius;
}

void SetSphereStatic(PhysicsWorld* world, int index, bool is_static) {
    if (!world || !world->spheres) return;
    SphereSystem* system = world->spheres;
    if (!system || index < 0 || index >= system->count) return;
    
    system->is_static[index] = is_static;
    
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
    
    // Wake up sleeping objects when force is applied
    if (system->is_sleeping[index]) {
        system->is_sleeping[index] = false;
        system->sleep_timers[index] = 0.0f;
    }
    
    system->forces[index] = Vec3Add(system->forces[index], force);
}

void AddSphereTorque(PhysicsWorld* world, int index, Vec3 torque) {
    if (!world || !world->spheres) return;
    SphereSystem* system = world->spheres;
    if (!system || index < 0 || index >= system->count) return;
    
    // Wake up sleeping objects when torque is applied
    if (system->is_sleeping[index]) {
        system->is_sleeping[index] = false;
        system->sleep_timers[index] = 0.0f;
    }
    
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
    system->count = 0;
    system->capacity = capacity;
    
    if (!system->positions || !system->sizes || !system->rotations || 
        !system->velocities || !system->forces || !system->masses ||
        !system->angular_velocities || !system->torques || !system->inertias || 
        !system->is_static || !system->is_sleeping || !system->sleep_timers || !system->collision_callbacks) {
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
    free(system);
}

int AddBox(PhysicsWorld* world, Vec3 position, Vec3 size, Quat rotation) {
    if (!world || !world->boxes) return -1;
    BoxSystem* system = world->boxes;
    if (!system || system->count >= system->capacity) return -1;
    
    int index = system->count;
    system->positions[index] = position;
    system->sizes[index] = size;
    system->rotations[index] = rotation;
    system->velocities[index] = (Vec3){ 0.0f, 0.0f, 0.0f };
    system->forces[index] = (Vec3){ 0.0f, 0.0f, 0.0f };
    system->masses[index] = 1.0f;
    system->angular_velocities[index] = (Vec3){ 0.0f, 0.0f, 0.0f };
    system->torques[index] = (Vec3){ 0.0f, 0.0f, 0.0f };
    
    float m = system->masses[index];
    float w = size.x * 2.0f;
    float h = size.y * 2.0f;
    float d = size.z * 2.0f;
    system->inertias[index] = (Vec3){
        (m / 12.0f) * (h*h + d*d),
        (m / 12.0f) * (w*w + d*d),
        (m / 12.0f) * (w*w + h*h)
    };
    
    // Initialize sleeping/static state explicitly
    system->is_static[index] = false;
    system->is_sleeping[index] = false;
    system->sleep_timers[index] = 0.0f;
    
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
    
    Vec3 size = system->sizes[index];
    float w = size.x * 2.0f;
    float h = size.y * 2.0f;
    float d = size.z * 2.0f;
    system->inertias[index] = (Vec3){
        (mass / 12.0f) * (h*h + d*d),
        (mass / 12.0f) * (w*w + d*d),
        (mass / 12.0f) * (w*w + h*h)
    };
}

void SetBoxStatic(PhysicsWorld* world, int index, bool is_static) {
    if (!world || !world->boxes) return;
    BoxSystem* system = world->boxes;
    if (!system || index < 0 || index >= system->count) return;
    
    system->is_static[index] = is_static;
    
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
    
    // Wake up sleeping objects when force is applied
    if (system->is_sleeping[index]) {
        system->is_sleeping[index] = false;
        system->sleep_timers[index] = 0.0f;
    }
    
    system->forces[index] = Vec3Add(system->forces[index], force);
}

void AddBoxTorque(PhysicsWorld* world, int index, Vec3 torque) {
    if (!world || !world->boxes) return;
    BoxSystem* system = world->boxes;
    if (!system || index < 0 || index >= system->count) return;
    
    // Wake up sleeping objects when torque is applied
    if (system->is_sleeping[index]) {
        system->is_sleeping[index] = false;
        system->sleep_timers[index] = 0.0f;
    }
    
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
        constraint->bodyB_type = contact->bodyB_type;
        constraint->bodyB_index = contact->bodyB_index;
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
    }
}

// ============================================================================
// CONSTRAINT SOLVER HELPER FUNCTIONS
// ============================================================================

// Get body velocity at a point
Vec3 GetConstraintBodyVelocityAtPoint(PhysicsWorld* world, int bodyType, int bodyIndex, Vec3 point) {
    if (bodyType == 0) { // Sphere
        if (bodyIndex < 0 || bodyIndex >= world->spheres->count || world->spheres->is_static[bodyIndex]) 
            return (Vec3){0, 0, 0};
        
        Vec3 linear_vel = world->spheres->velocities[bodyIndex];
        Vec3 angular_vel = world->spheres->angular_velocities[bodyIndex];
        Vec3 body_pos = world->spheres->positions[bodyIndex];
        Vec3 r = Vec3Sub(point, body_pos);
        Vec3 angular_component = Vec3Cross(angular_vel, r);
        
        return Vec3Add(linear_vel, angular_component);
        
    } else if (bodyType == 1) { // Box
        if (bodyIndex < 0 || bodyIndex >= world->boxes->count || world->boxes->is_static[bodyIndex]) 
            return (Vec3){0, 0, 0};
        
        Vec3 linear_vel = world->boxes->velocities[bodyIndex];
        Vec3 angular_vel = world->boxes->angular_velocities[bodyIndex];
        Vec3 body_pos = world->boxes->positions[bodyIndex];
        Vec3 r = Vec3Sub(point, body_pos);
        Vec3 angular_component = Vec3Cross(angular_vel, r);
        
        return Vec3Add(linear_vel, angular_component);
        
    } else { // SDF - static
        return (Vec3){0, 0, 0};
    }
}

// Get body inverse mass
float GetConstraintBodyInverseMass(PhysicsWorld* world, int bodyType, int bodyIndex) {
    if (bodyType == 0) { // Sphere
        if (bodyIndex < 0 || bodyIndex >= world->spheres->count || world->spheres->is_static[bodyIndex]) 
            return 0.0f;
        return 1.0f / world->spheres->masses[bodyIndex];
    } else if (bodyType == 1) { // Box
        if (bodyIndex < 0 || bodyIndex >= world->boxes->count || world->boxes->is_static[bodyIndex]) 
            return 0.0f;
        return 1.0f / world->boxes->masses[bodyIndex];
    } else { // SDF - static
        return 0.0f;
    }
}

// Get body inverse inertia (simplified - assumes scalar inertia for spheres)
Vec3 GetConstraintBodyInverseInertia(PhysicsWorld* world, int bodyType, int bodyIndex) {
    if (bodyType == 0) { // Sphere
        if (bodyIndex < 0 || bodyIndex >= world->spheres->count || world->spheres->is_static[bodyIndex]) 
            return (Vec3){0, 0, 0};
        float inv_inertia = 1.0f / world->spheres->inertias[bodyIndex];
        return (Vec3){inv_inertia, inv_inertia, inv_inertia};
    } else if (bodyType == 1) { // Box
        if (bodyIndex < 0 || bodyIndex >= world->boxes->count || world->boxes->is_static[bodyIndex]) 
            return (Vec3){0, 0, 0};
        Vec3 inertia = world->boxes->inertias[bodyIndex];
        return (Vec3){1.0f / inertia.x, 1.0f / inertia.y, 1.0f / inertia.z};
    } else { // SDF - static
        return (Vec3){0, 0, 0};
    }
}

void PreSolveConstraints(PhysicsWorld* world) {
    if (!world || !world->constraints) return;
    
    float dt = world->dt;
    
    for (int i = 0; i < world->constraints->count; i++) {
        ContactConstraint* c = &world->constraints->constraints[i];
        
        // Get body data
        Vec3 rA = Vec3Sub(c->contact_point, 
                         (c->bodyA_type == 0) ? world->spheres->positions[c->bodyA_index] 
                                              : (c->bodyA_type == 1) ? world->boxes->positions[c->bodyA_index]
                                                                     : c->contact_point);
        Vec3 rB = Vec3Sub(c->contact_point, 
                         (c->bodyB_type == 0) ? world->spheres->positions[c->bodyB_index] 
                                              : (c->bodyB_type == 1) ? world->boxes->positions[c->bodyB_index]
                                                                     : c->contact_point);
        
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
    if (bodyType == 0) { // Sphere
        if (bodyIndex < 0 || bodyIndex >= world->spheres->count || world->spheres->is_static[bodyIndex]) 
            return;
        
        float invMass = 1.0f / world->spheres->masses[bodyIndex];
        float invInertia = 1.0f / world->spheres->inertias[bodyIndex];
        
        world->spheres->velocities[bodyIndex] = Vec3Add(world->spheres->velocities[bodyIndex], 
                                                       Vec3Scale(linear_impulse, invMass));
        world->spheres->angular_velocities[bodyIndex] = Vec3Add(world->spheres->angular_velocities[bodyIndex], 
                                                               Vec3Scale(angular_impulse, invInertia));
        
        // Wake up sleeping bodies
        if (world->spheres->is_sleeping[bodyIndex]) {
            world->spheres->is_sleeping[bodyIndex] = false;
            world->spheres->sleep_timers[bodyIndex] = 0.0f;
        }
        
    } else if (bodyType == 1) { // Box
        if (bodyIndex < 0 || bodyIndex >= world->boxes->count || world->boxes->is_static[bodyIndex]) 
            return;
        
        float invMass = 1.0f / world->boxes->masses[bodyIndex];
        Vec3 invInertia = GetConstraintBodyInverseInertia(world, bodyType, bodyIndex);
        
        world->boxes->velocities[bodyIndex] = Vec3Add(world->boxes->velocities[bodyIndex], 
                                                     Vec3Scale(linear_impulse, invMass));
        
        // For boxes, need to transform angular impulse to local space
        Quat rotation = world->boxes->rotations[bodyIndex];
        Vec3 localAngularImpulse = QuatRotateVec3(QuatConjugate(rotation), angular_impulse);
        Vec3 localAngularDelta = (Vec3){
            localAngularImpulse.x * invInertia.x,
            localAngularImpulse.y * invInertia.y,
            localAngularImpulse.z * invInertia.z
        };
        Vec3 worldAngularDelta = QuatRotateVec3(rotation, localAngularDelta);
        
        world->boxes->angular_velocities[bodyIndex] = Vec3Add(world->boxes->angular_velocities[bodyIndex], 
                                                             worldAngularDelta);
        
        // Wake up sleeping bodies
        if (world->boxes->is_sleeping[bodyIndex]) {
            world->boxes->is_sleeping[bodyIndex] = false;
            world->boxes->sleep_timers[bodyIndex] = 0.0f;
        }
    }
    // SDF bodies are static - no impulse applied
}

void SolveVelocityConstraints(PhysicsWorld* world) {
    if (!world || !world->constraints) return;
    
    
    for (int i = 0; i < world->constraints->count; i++) {
        ContactConstraint* c = &world->constraints->constraints[i];
        
        // Get current body positions and velocities
        Vec3 posA = (c->bodyA_type == 0) ? world->spheres->positions[c->bodyA_index] 
                                        : (c->bodyA_type == 1) ? world->boxes->positions[c->bodyA_index]
                                                               : c->contact_point;
        Vec3 posB = (c->bodyB_type == 0) ? world->spheres->positions[c->bodyB_index] 
                                        : (c->bodyB_type == 1) ? world->boxes->positions[c->bodyB_index]
                                                               : c->contact_point;
        
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
        
        // Apply position corrections directly
        if (c->bodyA_type == 0 && invMassA > 0.0f) { // Sphere A
            Vec3 correctionA = Vec3Scale(correction, -invMassA);
            world->spheres->positions[c->bodyA_index] = 
                Vec3Add(world->spheres->positions[c->bodyA_index], correctionA);
        } else if (c->bodyA_type == 1 && invMassA > 0.0f) { // Box A
            Vec3 correctionA = Vec3Scale(correction, -invMassA);
            world->boxes->positions[c->bodyA_index] = 
                Vec3Add(world->boxes->positions[c->bodyA_index], correctionA);
        }
        
        if (c->bodyB_type == 0 && invMassB > 0.0f) { // Sphere B
            Vec3 correctionB = Vec3Scale(correction, invMassB);
            world->spheres->positions[c->bodyB_index] = 
                Vec3Add(world->spheres->positions[c->bodyB_index], correctionB);
        } else if (c->bodyB_type == 1 && invMassB > 0.0f) { // Box B
            Vec3 correctionB = Vec3Scale(correction, invMassB);
            world->boxes->positions[c->bodyB_index] = 
                Vec3Add(world->boxes->positions[c->bodyB_index], correctionB);
        }
        // SDF bodies (type 2) are static - no position correction needed
    }
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
    
    // Set default values
    world->gravity = (Vec3){ 0.0f, 0.0f, 0.0f }; // No gravity for now
    world->dt = 1.0f / 60.0f;
    world->usesSDF = false;
    world->SDFCollider = NULL;
    
    // Initialize physics settings with defaults
    ResetPhysicsSettingsToDefaults(world);
    
    // Check for allocation failures
    if (!world->spheres || !world->boxes || !world->collisions || !world->collisions->contacts || !world->constraints) {
        DestroyPhysicsWorld(world);
        return NULL;
    }
    
    return world;
}

void DestroyPhysicsWorld(PhysicsWorld* world) {
    if (!world) return;
    
    if (world->spheres) DestroySphereSystem(world->spheres);
    if (world->boxes) DestroyBoxSystem(world->boxes);
    
    if (world->collisions) {
        if (world->collisions->contacts) free(world->collisions->contacts);
        free(world->collisions);
    }
    
    if (world->constraints) DestroyConstraintCollection(world->constraints);
    
    
    free(world);
}

// Register user SDF function and enable SDF collision detection
void UseSDF(PhysicsWorld* world, SDFFunction sdf_function) {
    if (!world) return;
    
    world->SDFCollider = sdf_function;
    world->usesSDF = (sdf_function != NULL);
}

// Update sleep states for all objects
void UpdateSleepStates(PhysicsWorld* world) {
    if (!world) return;
    
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
    
    // Update sphere sleep states
    for (int i = 0; i < world->spheres->count; i++) {
        // Skip static objects - they don't need sleep
        if (world->spheres->is_static[i]) continue;
        
        // Calculate kinetic energy
        Vec3 velocity = world->spheres->velocities[i];
        Vec3 angular_velocity = world->spheres->angular_velocities[i];
        float linear_speed = Vec3Length(velocity);
        float angular_speed = Vec3Length(angular_velocity);
        
        // Check if object is below sleep thresholds
        bool low_energy = (linear_speed < linear_threshold && angular_speed < angular_threshold);
        
        if (low_energy) {
            // Accumulate sleep time
            world->spheres->sleep_timers[i] += dt;
            
            // Put to sleep if timer exceeds threshold
            if (world->spheres->sleep_timers[i] >= time_required && !world->spheres->is_sleeping[i]) {
                world->spheres->is_sleeping[i] = true;
                world->spheres->velocities[i] = (Vec3){ 0.0f, 0.0f, 0.0f };
                world->spheres->angular_velocities[i] = (Vec3){ 0.0f, 0.0f, 0.0f };
            }
        } else {
            // Reset sleep timer and wake up if sleeping
            world->spheres->sleep_timers[i] = 0.0f;
            world->spheres->is_sleeping[i] = false;
        }
    }
    
    // Update box sleep states
    for (int i = 0; i < world->boxes->count; i++) {
        // Skip static objects - they don't need sleep
        if (world->boxes->is_static[i]) continue;
        
        // Calculate kinetic energy
        Vec3 velocity = world->boxes->velocities[i];
        Vec3 angular_velocity = world->boxes->angular_velocities[i];
        float linear_speed = Vec3Length(velocity);
        float angular_speed = Vec3Length(angular_velocity);
        
        // Check if object is below sleep thresholds
        bool low_energy = (linear_speed < linear_threshold && angular_speed < angular_threshold);
        
        // Debug output for first box
        //if (i == 0) {
        //    printf("Box[0]: linear_speed=%.6f, angular_speed=%.6f, low_energy=%d, sleep_timer=%.3f, is_sleeping=%d\n", 
        //           linear_speed, angular_speed, low_energy, world->boxes->sleep_timers[i], world->boxes->is_sleeping[i]);
        //}
        
        if (low_energy) {
            // Accumulate sleep time
            world->boxes->sleep_timers[i] += dt;
            
            // Put to sleep if timer exceeds threshold
            if (world->boxes->sleep_timers[i] >= time_required && !world->boxes->is_sleeping[i]) {
                //printf("BOX[%d] GOING TO SLEEP! timer=%.3f\n", i, world->boxes->sleep_timers[i]);
                world->boxes->is_sleeping[i] = true;
                world->boxes->velocities[i] = (Vec3){ 0.0f, 0.0f, 0.0f };
                world->boxes->angular_velocities[i] = (Vec3){ 0.0f, 0.0f, 0.0f };
            }
        } else {
            // Reset sleep timer and wake up if sleeping
            world->boxes->sleep_timers[i] = 0.0f;
            if (world->boxes->is_sleeping[i]) {
                //printf("BOX[%d] WAKING UP! linear_speed=%.6f, angular_speed=%.6f\n", i, linear_speed, angular_speed);
            }
            world->boxes->is_sleeping[i] = false;
        }
    }
}

// Main physics simulation step - CORRECT Box2D pipeline order
void PhysicsStep(PhysicsWorld* world, float deltaTime) {
    if (!world) return;
    
    // Skip frame if deltaTime is zero or too small (Raylib first frame issue)
    if (deltaTime <= 0.0f || deltaTime < 1e-6f) {
        return;
    }
    
    world->dt = deltaTime;
    
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
    
    // 6. Update sleep states
    UpdateSleepStates(world);
    
    // 7. Cleanup
    CleanupPhysics(world);
}

// Physics sub-systems implementation
void ApplyForces(PhysicsWorld* world) {
    if (!world) return;
    
    // Apply gravity to all dynamic spheres
    for (int i = 0; i < world->spheres->count; i++) {
        if (world->spheres->is_static[i]) continue;  // Skip static objects
        Vec3 gravityForce = Vec3Scale(world->gravity, world->spheres->masses[i]);
        world->spheres->forces[i] = Vec3Add(world->spheres->forces[i], gravityForce);
    }
    
    // Apply gravity to all dynamic boxes
    for (int i = 0; i < world->boxes->count; i++) {
        if (world->boxes->is_static[i]) continue;  // Skip static objects
        Vec3 gravityForce = Vec3Scale(world->gravity, world->boxes->masses[i]);
        world->boxes->forces[i] = Vec3Add(world->boxes->forces[i], gravityForce);
    }
    
    // User forces are already added via AddSphereForce/AddBoxForce calls
}

// BOX2D PIPELINE: Integrate velocities first, positions later after constraint solving
void IntegrateVelocities(PhysicsWorld* world) {
    if (!world) return;
    
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
            
            CollisionContact* contact = &world->collisions->contacts[world->collisions->count];
            if (CheckSphereCollisionWithData(
                world->spheres->positions[i], world->spheres->radii[i],
                world->spheres->positions[j], world->spheres->radii[j], contact)) {
                
                contact->bodyA_type = 0; contact->bodyA_index = i;
                contact->bodyB_type = 0; contact->bodyB_index = j;
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
                    contact->bodyB_type = manifold.bodyB_type;
                    contact->bodyB_index = manifold.bodyB_index;
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
            
            CollisionContact* contact = &world->collisions->contacts[world->collisions->count];
            if (CheckSphereBoxCollisionWithData(
                world->spheres->positions[i], world->spheres->radii[i],
                world->boxes->positions[j], world->boxes->sizes[j], world->boxes->rotations[j], contact)) {
                
                contact->bodyA_type = 0; contact->bodyA_index = i;  // A = sphere
                contact->bodyB_type = 1; contact->bodyB_index = j;  // B = box
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
                
                contact->bodyA_type = 0; contact->bodyA_index = i;  // A = sphere
                contact->bodyB_type = 2; contact->bodyB_index = 0; // B = SDF (type 2, index 0)
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

// Box-SDF collision detection using 3-stage approach
bool CheckBoxSDFCollision(PhysicsWorld* world, Vec3 boxPos, Vec3 boxSize, Quat boxRot, CollisionContact* contact) {
    // Stage 1: Circumscribed sphere test
    // Calculate radius of circumscribed sphere (distance from center to farthest corner)
    float circumscribed_radius = Vec3Length(boxSize);
    
    // Evaluate SDF at box center
    float center_distance = world->SDFCollider(boxPos);
    
    // If SDF distance > circumscribed radius, no collision possible
    if (center_distance > circumscribed_radius) {
        return false;
    }
    
    // Stage 2: Vertex test - check all 8 vertices of the box
    Vec3 vertices[8];
    OBB box_obb = CreateOBB(boxPos, boxSize, boxRot);
    GetOBBVertices(&box_obb, vertices);
    
    float deepest_penetration = -999999.0f;
    Vec3 deepest_point = { 0.0f, 0.0f, 0.0f };
    bool collision_found = false;
    
    for (int i = 0; i < 8; i++) {
        float vertex_distance = world->SDFCollider(vertices[i]);
        
        // If vertex is inside SDF (distance < 0), we have collision
        if (vertex_distance < 0.0f) {
            collision_found = true;
            float penetration = -vertex_distance;
            
            if (penetration > deepest_penetration) {
                deepest_penetration = penetration;
                deepest_point = vertices[i];
            }
        }
    }
    
    // Stage 3: Face center test - check center of each face (check both vertices AND face centers)
    Vec3 face_centers[6];
    
    // Calculate face centers: box center + (half_extent * face_normal)
    Vec3 x_axis = QuatRotateVec3(boxRot, (Vec3){1.0f, 0.0f, 0.0f});
    Vec3 y_axis = QuatRotateVec3(boxRot, (Vec3){0.0f, 1.0f, 0.0f});
    Vec3 z_axis = QuatRotateVec3(boxRot, (Vec3){0.0f, 0.0f, 1.0f});
    
    face_centers[0] = Vec3Add(boxPos, Vec3Scale(x_axis, boxSize.x));   // +X face
    face_centers[1] = Vec3Sub(boxPos, Vec3Scale(x_axis, boxSize.x));   // -X face
    face_centers[2] = Vec3Add(boxPos, Vec3Scale(y_axis, boxSize.y));   // +Y face
    face_centers[3] = Vec3Sub(boxPos, Vec3Scale(y_axis, boxSize.y));   // -Y face
    face_centers[4] = Vec3Add(boxPos, Vec3Scale(z_axis, boxSize.z));   // +Z face
    face_centers[5] = Vec3Sub(boxPos, Vec3Scale(z_axis, boxSize.z));   // -Z face
    
    for (int i = 0; i < 6; i++) {
        float face_distance = world->SDFCollider(face_centers[i]);
        
        // If face center is inside SDF (distance < 0), we have collision
        if (face_distance < 0.0f) {
            float penetration = -face_distance;
            
            if (penetration > deepest_penetration) {
                deepest_penetration = penetration;
                deepest_point = face_centers[i];
            }
            collision_found = true; // Mark as collision found
        }
    }
    
    if (collision_found && contact) {
        contact->penetration = deepest_penetration;
        Vec3 sdf_gradient = SDF_EstimateNormal(world, deepest_point);
        contact->contact_normal = Vec3Scale(sdf_gradient, -1.0f); // Negate to point from box to surface
        contact->contact_point = deepest_point;
        return true;
    }
    
    return false;
}

// Generate multiple contact points for box-SDF collision (for stable stacking)
void GenerateBoxSDFContacts(PhysicsWorld* world, int boxIndex) {
    if (!world || !world->usesSDF || !world->SDFCollider) return;
    if (boxIndex < 0 || boxIndex >= world->boxes->count) return;
    
    Vec3 boxPos = world->boxes->positions[boxIndex];
    Vec3 boxSize = world->boxes->sizes[boxIndex];
    Quat boxRot = world->boxes->rotations[boxIndex];
    
    // Get all 8 vertices of the box
    Vec3 vertices[8];
    OBB box_obb = CreateOBB(boxPos, boxSize, boxRot);
    GetOBBVertices(&box_obb, vertices);
    
    // Check each vertex for SDF collision
    for (int i = 0; i < 8; i++) {
        if (world->collisions->count >= world->collisions->capacity) break;
        
        float vertex_distance = world->SDFCollider(vertices[i]);
        
        // If vertex is inside SDF (distance < 0), generate contact
        if (vertex_distance < -0.001f) { // Small threshold to avoid micro-contacts
            CollisionContact* contact = &world->collisions->contacts[world->collisions->count];
            
            contact->penetration = -vertex_distance;
            Vec3 sdf_gradient = SDF_EstimateNormal(world, vertices[i]);
            contact->contact_normal = Vec3Scale(sdf_gradient, -1.0f); // Normal points from box to SDF
            contact->contact_point = vertices[i];
            
            contact->bodyA_type = 1; contact->bodyA_index = boxIndex;  // A = box
            contact->bodyB_type = 2; contact->bodyB_index = 0;        // B = SDF
            
            world->collisions->count++;
        }
    }
    
    // Also check edge midpoints for better contact coverage
    // Generate contacts at the midpoint of each edge touching the ground
    for (int edge = 0; edge < 12; edge++) {
        if (world->collisions->count >= world->collisions->capacity) break;
        
        // Box edge pairs (vertex indices for each of the 12 edges)
        int edgePairs[12][2] = {
            {0,1}, {1,2}, {2,3}, {3,0}, // Bottom face edges
            {4,5}, {5,6}, {6,7}, {7,4}, // Top face edges  
            {0,4}, {1,5}, {2,6}, {3,7}  // Vertical edges
        };
        
        Vec3 edgeMidpoint = Vec3Scale(Vec3Add(vertices[edgePairs[edge][0]], vertices[edgePairs[edge][1]]), 0.5f);
        float edge_distance = world->SDFCollider(edgeMidpoint);
        
        if (edge_distance < -0.001f) { // Edge is inside SDF
            CollisionContact* contact = &world->collisions->contacts[world->collisions->count];
            
            contact->penetration = -edge_distance;
            Vec3 sdf_gradient = SDF_EstimateNormal(world, edgeMidpoint);
            contact->contact_normal = Vec3Scale(sdf_gradient, -1.0f);
            contact->contact_point = edgeMidpoint;
            
            contact->bodyA_type = 1; contact->bodyA_index = boxIndex;
            contact->bodyB_type = 2; contact->bodyB_index = 0;
            
            world->collisions->count++;
        }
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
    manifold->bodyA_index = boxIndexA;
    manifold->bodyB_type = 1; // box
    manifold->bodyB_index = boxIndexB;
    
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
        constraint->bodyB_type = manifold->bodyB_type;
        constraint->bodyB_index = manifold->bodyB_index;
        
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
