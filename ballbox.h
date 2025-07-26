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
 *   #include "ballbox.h"
 * 
 * License: Public Domain / MIT (choose your preference)
 */

#ifndef BALLBOX_H
#define BALLBOX_H

#include <stdbool.h>
#include <stdlib.h>
#include <math.h>

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
    Vec3* positions;
    float* radii;
    Vec3* velocities;
    Vec3* forces;
    float* masses;
    Vec3* angular_velocities;
    Vec3* torques;
    float* inertias;
    bool* is_static;
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
    int count;
    int capacity;
} BoxSystem;

typedef struct {
    int bodyA_type;     // 0 = sphere, 1 = box
    int bodyA_index;
    int bodyB_type;     // 0 = sphere, 1 = box
    int bodyB_index;
    Vec3 contact_point;
    Vec3 contact_normal; // From A to B
    float penetration;
} CollisionContact;

typedef struct {
    CollisionContact* contacts;
    int count;
    int capacity;
} CollisionCollection;

typedef struct {
    SphereSystem* spheres;
    BoxSystem* boxes;
    CollisionCollection* collisions;
    Vec3 gravity;
    float dt;
} PhysicsWorld;

// Create/destroy physics world
PhysicsWorld* CreatePhysicsWorld(int max_spheres, int max_boxes, int max_collisions);
void DestroyPhysicsWorld(PhysicsWorld* world);

// Main simulation step - call once per frame
void PhysicsStep(PhysicsWorld* world, float deltaTime);

void ApplyForces(PhysicsWorld* world);
void IntegrateMotion(PhysicsWorld* world);
void CollectCollisions(PhysicsWorld* world);
void ResolveCollisions(PhysicsWorld* world);
void CleanupPhysics(PhysicsWorld* world);

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
bool CheckSphereBoxCollision(Vec3 spherePos, float sphereRadius, Vec3 boxPos, Vec3 boxSize, Quat boxRot);
bool CheckSphereBoxCollisionWithData(Vec3 spherePos, float sphereRadius, Vec3 boxPos, Vec3 boxSize, Quat boxRot, CollisionContact* contact);
bool CheckBoxCollision(Vec3 pos1, Vec3 size1, Quat rot1, Vec3 pos2, Vec3 size2, Quat rot2);
bool CheckBoxCollisionWithData(Vec3 pos1, Vec3 size1, Quat rot1, Vec3 pos2, Vec3 size2, Quat rot2, CollisionContact* contact);

// Implementation

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
    system->count = 0;
    system->capacity = capacity;
    
    if (!system->positions || !system->radii || !system->velocities || !system->forces || 
        !system->masses || !system->angular_velocities || !system->torques || !system->inertias || !system->is_static) {
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
    
    system->forces[index] = Vec3Add(system->forces[index], force);
}

void AddSphereTorque(PhysicsWorld* world, int index, Vec3 torque) {
    if (!world || !world->spheres) return;
    SphereSystem* system = world->spheres;
    if (!system || index < 0 || index >= system->count) return;
    
    system->torques[index] = Vec3Add(system->torques[index], torque);
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
    system->count = 0;
    system->capacity = capacity;
    
    if (!system->positions || !system->sizes || !system->rotations || 
        !system->velocities || !system->forces || !system->masses ||
        !system->angular_velocities || !system->torques || !system->inertias || !system->is_static) {
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
    
    system->forces[index] = Vec3Add(system->forces[index], force);
}

void AddBoxTorque(PhysicsWorld* world, int index, Vec3 torque) {
    if (!world || !world->boxes) return;
    BoxSystem* system = world->boxes;
    if (!system || index < 0 || index >= system->count) return;
    
    system->torques[index] = Vec3Add(system->torques[index], torque);
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
    
    // Set default values
    world->gravity = (Vec3){ 0.0f, 0.0f, 0.0f }; // No gravity for now
    world->dt = 1.0f / 60.0f;
    
    // Check for allocation failures
    if (!world->spheres || !world->boxes || !world->collisions || !world->collisions->contacts) {
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
    
    free(world);
}

// Main physics simulation step
void PhysicsStep(PhysicsWorld* world, float deltaTime) {
    if (!world) return;
    
    world->dt = deltaTime;
    
    // 1. Apply Forces (keep your existing)
    ApplyForces(world);
    
    // 2. Integration (keep your existing)
    IntegrateMotion(world);
    
    // 3. Collision Detection (use Randy's approach)
    CollectCollisions(world);
    
    // 4. Collision Resolution (use Randy's approach)
    ResolveCollisions(world);
    
    // 5. Cleanup (keep your existing)
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

void IntegrateMotion(PhysicsWorld* world) {
    if (!world) return;
    
    float dt = world->dt;
    
    // Integrate spheres (linear + angular motion)
    for (int i = 0; i < world->spheres->count; i++) {
        // Skip static objects
        if (world->spheres->is_static[i]) continue;
        
        // Calculate acceleration: a = F / m
        Vec3 acceleration = Vec3Scale(world->spheres->forces[i], 1.0f / world->spheres->masses[i]);
        
        // Integrate velocity: v += a * dt
        world->spheres->velocities[i] = Vec3Add(world->spheres->velocities[i], 
                                               Vec3Scale(acceleration, dt));
        
        // Integrate position: pos += v * dt
        world->spheres->positions[i] = Vec3Add(world->spheres->positions[i], 
                                              Vec3Scale(world->spheres->velocities[i], dt));
        
        // Angular integration for spheres
        // Calculate angular acceleration: α = τ / I
        Vec3 angular_acceleration = Vec3Scale(world->spheres->torques[i], 1.0f / world->spheres->inertias[i]);
        
        // Integrate angular velocity: ω += α * dt
        world->spheres->angular_velocities[i] = Vec3Add(world->spheres->angular_velocities[i], 
                                                       Vec3Scale(angular_acceleration, dt));
        // Note: Spheres don't need orientation integration since they're rotationally symmetric
    }
    
    // Integrate boxes (linear + angular motion)
    for (int i = 0; i < world->boxes->count; i++) {
        // Skip static objects
        if (world->boxes->is_static[i]) continue;
        
        // Linear integration
        // Calculate acceleration: a = F / m
        Vec3 acceleration = Vec3Scale(world->boxes->forces[i], 1.0f / world->boxes->masses[i]);
        
        // Integrate velocity: v += a * dt
        world->boxes->velocities[i] = Vec3Add(world->boxes->velocities[i], 
                                             Vec3Scale(acceleration, dt));
        
        // Integrate position: pos += v * dt
        world->boxes->positions[i] = Vec3Add(world->boxes->positions[i], 
                                            Vec3Scale(world->boxes->velocities[i], dt));
        
        // Angular integration
        // Calculate angular acceleration in world space: α = I⁻¹ * τ 
        // For diagonal inertia tensor: α = τ / I (component-wise)
        Vec3 angular_acceleration;
        angular_acceleration.x = world->boxes->torques[i].x / world->boxes->inertias[i].x;
        angular_acceleration.y = world->boxes->torques[i].y / world->boxes->inertias[i].y; 
        angular_acceleration.z = world->boxes->torques[i].z / world->boxes->inertias[i].z;
        
        // Integrate angular velocity: ω += α * dt
        world->boxes->angular_velocities[i] = Vec3Add(world->boxes->angular_velocities[i], 
                                                     Vec3Scale(angular_acceleration, dt));
        
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
    
    // Apply damping to reduce sliding and spinning
    const float linear_damping = 0.98f;  // Reduces sliding
    const float angular_damping = 0.9f; // Reduces spinning
    
    // Apply damping to spheres
    for (int i = 0; i < world->spheres->count; i++) {
        if (world->spheres->is_static[i]) continue;
        world->spheres->velocities[i] = Vec3Scale(world->spheres->velocities[i], linear_damping);
        world->spheres->angular_velocities[i] = Vec3Scale(world->spheres->angular_velocities[i], angular_damping);
    }
    
    // Apply damping to boxes
    for (int i = 0; i < world->boxes->count; i++) {
        if (world->boxes->is_static[i]) continue;
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
            }
        }
    }
    
    // Box-Box collisions (A=first box, B=second box)
    for (int i = 0; i < world->boxes->count; i++) {
        for (int j = i + 1; j < world->boxes->count; j++) {
            if (world->collisions->count >= world->collisions->capacity) break;
            
            CollisionContact* contact = &world->collisions->contacts[world->collisions->count];
            if (CheckBoxCollisionWithData(
                world->boxes->positions[i], world->boxes->sizes[i], world->boxes->rotations[i],
                world->boxes->positions[j], world->boxes->sizes[j], world->boxes->rotations[j], contact)) {
                
                contact->bodyA_type = 1; contact->bodyA_index = i;
                contact->bodyB_type = 1; contact->bodyB_index = j;
                world->collisions->count++;
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
            }
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
    } else { // Box
        bodyPos = world->boxes->positions[bodyIndex];
        bodyVel = world->boxes->velocities[bodyIndex];
        bodyAngVel = world->boxes->angular_velocities[bodyIndex];
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
                                               : world->boxes->positions[contact->bodyA_index]);
    Vec3 r2 = Vec3Sub(contact->contact_point, 
                     (contact->bodyB_type == 0) ? world->spheres->positions[contact->bodyB_index] 
                                               : world->boxes->positions[contact->bodyB_index]);
    
    // Check for static objects and handle infinite mass
    bool is_static1 = (contact->bodyA_type == 0) ? world->spheres->is_static[contact->bodyA_index] 
                                                 : world->boxes->is_static[contact->bodyA_index];
    bool is_static2 = (contact->bodyB_type == 0) ? world->spheres->is_static[contact->bodyB_index] 
                                                 : world->boxes->is_static[contact->bodyB_index];
    
    // Get masses (use very large mass for static objects)
    //float m1 = is_static1 ? 1e10f : ((contact->bodyA_type == 0) ? world->spheres->masses[contact->bodyA_index] 
    //                                                            : world->boxes->masses[contact->bodyA_index]);
    //float m2 = is_static2 ? 1e10f : ((contact->bodyB_type == 0) ? world->spheres->masses[contact->bodyB_index] 
    //
    float inv_m1 = is_static1 ? 0.0f : 
    1.0f / ((contact->bodyA_type == 0) ? 
            world->spheres->masses[contact->bodyA_index] : 
            world->boxes->masses[contact->bodyA_index]);

    float inv_m2 = is_static2 ? 0.0f : 
        1.0f / ((contact->bodyB_type == 0) ? 
                world->spheres->masses[contact->bodyB_index] : 
                world->boxes->masses[contact->bodyB_index]);
                                                         
    
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
        
    } else { // Box
        // Skip static objects
        if (world->boxes->is_static[bodyIndex]) return;
        
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
    }
}

// Position correction to reduce penetration using Randy's approach
void PositionalCorrection(PhysicsWorld* world, CollisionContact* contact) {
    const float percent = 0.8f;  // Randy uses 20% to 80%
    const float slop = 0.05f;    // Randy uses 0.01 to 0.1
    
    if (contact->penetration <= slop) return;
    
    // Get inverse masses
    bool is_static1 = (contact->bodyA_type == 0) ? world->spheres->is_static[contact->bodyA_index] 
                                                 : world->boxes->is_static[contact->bodyA_index];
    bool is_static2 = (contact->bodyB_type == 0) ? world->spheres->is_static[contact->bodyB_index] 
                                                 : world->boxes->is_static[contact->bodyB_index];
    
    float inv_massA = is_static1 ? 0.0f : 
        1.0f / ((contact->bodyA_type == 0) ? world->spheres->masses[contact->bodyA_index] 
                                           : world->boxes->masses[contact->bodyA_index]);
    float inv_massB = is_static2 ? 0.0f : 
        1.0f / ((contact->bodyB_type == 0) ? world->spheres->masses[contact->bodyB_index] 
                                           : world->boxes->masses[contact->bodyB_index]);
    
    float total_inv_mass = inv_massA + inv_massB;
    if (total_inv_mass < 1e-10f) return;
    
    // RANDY'S FORMULA: correction = (penetration - slop) / total_inv_mass * percent * normal
    Vec3 correction = Vec3Scale(contact->contact_normal, 
                               (contact->penetration - slop) / total_inv_mass * percent);
    
    // RANDY'S APPLICATION: A moves backward, B moves forward
    Vec3 correctionA = Vec3Scale(correction, -inv_massA);
    Vec3 correctionB = Vec3Scale(correction, inv_massB);
    
    // Apply corrections
    if (contact->bodyA_type == 0 && !is_static1) {
        world->spheres->positions[contact->bodyA_index] = 
            Vec3Add(world->spheres->positions[contact->bodyA_index], correctionA);
    } else if (!is_static1) {
        world->boxes->positions[contact->bodyA_index] = 
            Vec3Add(world->boxes->positions[contact->bodyA_index], correctionA);
    }
    
    if (contact->bodyB_type == 0 && !is_static2) {
        world->spheres->positions[contact->bodyB_index] = 
            Vec3Add(world->spheres->positions[contact->bodyB_index], correctionB);
    } else if (!is_static2) {
        world->boxes->positions[contact->bodyB_index] = 
            Vec3Add(world->boxes->positions[contact->bodyB_index], correctionB);
    }
    
    // Cancel out artificial velocity induced by small position corrections
    const float correction_threshold = 0.05f;  // Only for very small corrections
    float correctionMagA = Vec3Length(correctionA);
    float correctionMagB = Vec3Length(correctionB);
    
    if (correctionMagA > 0.0f && correctionMagA < correction_threshold) {
        Vec3 artificial_velocityA = Vec3Scale(correctionA, 1.0f / world->dt);
        if (contact->bodyA_type == 0 && !is_static1) {
            world->spheres->velocities[contact->bodyA_index] = 
                Vec3Add(world->spheres->velocities[contact->bodyA_index], artificial_velocityA);
        } else if (!is_static1) {
            world->boxes->velocities[contact->bodyA_index] = 
                Vec3Add(world->boxes->velocities[contact->bodyA_index], artificial_velocityA);
        }
    }
    
    if (correctionMagB > 0.0f && correctionMagB < correction_threshold) {
        Vec3 artificial_velocityB = Vec3Scale(correctionB, 1.0f / world->dt);
        if (contact->bodyB_type == 0 && !is_static2) {
            world->spheres->velocities[contact->bodyB_index] = 
                Vec3Add(world->spheres->velocities[contact->bodyB_index], artificial_velocityB);
        } else if (!is_static2) {
            world->boxes->velocities[contact->bodyB_index] = 
                Vec3Add(world->boxes->velocities[contact->bodyB_index], artificial_velocityB);
        }
    }
}


void ResolveCollisions(PhysicsWorld* world) {
    if (!world || !world->collisions) return;
    
    const int iterations = 8;        
    const float restitution = 0.0f;  // Much less bouncy  
    
    for (int iter = 0; iter < iterations; iter++) {
        for (int i = 0; i < world->collisions->count; i++) {
            CollisionContact* contact = &world->collisions->contacts[i];
            
            // Get velocities at contact point
            Vec3 velA = GetBodyVelocityAtPoint(world, contact->bodyA_type, contact->bodyA_index, contact->contact_point);
            Vec3 velB = GetBodyVelocityAtPoint(world, contact->bodyB_type, contact->bodyB_index, contact->contact_point);
            
            // RANDY'S WAY: Relative velocity = B - A
            Vec3 relativeVel = Vec3Sub(velB, velA);
            
            // Calculate relative velocity in terms of the normal direction
            float velAlongNormal = Vec3Dot(relativeVel, contact->contact_normal);
            
            // Do not resolve if velocities are separating (Randy's check)
            if (velAlongNormal > 0) continue;
            
            // Calculate effective mass (your existing function should work)
            float effectiveMass = CalculateEffectiveMass(world, contact);
            
            // RANDY'S FORMULA: j = -(1 + e) * velAlongNormal * effectiveMass
            float j = -(1.0f + restitution) * velAlongNormal * effectiveMass;
            
            // Apply impulse (Randy's way)
            Vec3 impulse = Vec3Scale(contact->contact_normal, j);
            
            // RANDY'S APPLICATION: A gets -impulse, B gets +impulse
            ApplyImpulse(world, contact->bodyA_type, contact->bodyA_index, Vec3Scale(impulse, -1.0f), contact->contact_point);
            ApplyImpulse(world, contact->bodyB_type, contact->bodyB_index, impulse, contact->contact_point);
        }
    }
    
    // Position correction
    for (int i = 0; i < world->collisions->count; i++) {
        PositionalCorrection(world, &world->collisions->contacts[i]);
    }
}


void CleanupPhysics(PhysicsWorld* world) {
    if (!world) return;
    
    // Clear collision list for next frame
    if (world->collisions) {
        world->collisions->count = 0;
    }
    
    // Forces are already cleared in IntegrateMotion()
}

#endif // BALLBOX_H
