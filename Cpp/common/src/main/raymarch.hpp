#ifndef COMMON_SRC_MAIN_RAYMARCH_HPP
#define COMMON_SRC_MAIN_RAYMARCH_HPP

#include "datastructures.hpp"

namespace ww
{

/**
 * rm - raymarch
 */
namespace rm
{

/**
 * Default Map for the world of raymarching.
 */
tup MapDefault(tup const &Pos);
tup MapBoxAndSphere(tup const &Pos);

/**
 * RayMarch - move along ray to find hit.
 * @param: R - Ray
 * @param: PtrShape - Shape to compute distance against
 */
float RayMarch(ray const &R, shared_ptr_shape PtrShape);
float RayMarch(ray const &R, world const &World, shared_ptr_shape PtrShape);

/**
 */
tup MainImage(camera const &Camera, world const &World, int X, int Y, shared_ptr_shape PtrShape);

/**
 * Mainimage from https://www.shadertoy.com/view/Xds3zN.
 * Shows a lot of raymarching primitives.
 */
tup MainImage(tup const &FragCoord, mainimage_config const &Cfg);

/**
 */
canvas Render(camera const &Camera, world const &World);

/**
 */
float GetDistance(tup const &P, shared_ptr_shape PtrShape);
float GetDistance(tup const &P, world const &World);

/**
 */
tup GetNormal(tup const &P, shared_ptr_shape PtrShape);
tup CalcNormal(tup const &P, funcPtrMap = nullptr);

float sdBox(tup const &Pos, tup const &Box);
float sdCapsule(tup const &Pos, tup const &A, tup const &B, float Rad);
float sdSphere(tup const &Pos, float S);
float udTriangle(tup const &v1, tup const &v2, tup const &v3, tup const &Pos);
float udQuad(tup const &A, tup const &B, tup const &C, tup const &D, tup const &Pos);

float sdSphereSphereDistance(tup const &PosA, float RadA, tup const &PosB, float RadB);
float sdSphereSpherePenDist(tup const &PosA, float RadA, tup const &PosB, float RadB);
float sdSphereTriangle(tup const &v1, tup const &v2, tup const &v3, tup const &PosSphere, float Radius);
float sdCapsuleCapsule(capsule const &Ca, capsule const &Cb);
bool sdSphereSphereCollision(tup const &PosA, float RadA, tup const &PosB, float RadB);
bool sdSphereTriangleCollision(tup const &v1, tup const &v2, tup const &v3, tup const &PosSphere, float Radius);
tup ClosestPointOnLineSegment(tup const &A, tup const &B, tup const &Pos);

bool CapsuleTriangleIntersect(tup const &CapA, tup const &CapB,  //!<
                              float Radius,                      //!<
                              tup const &V0,                     //!<
                              tup const &V1,                     //!<
                              tup const &V2,                     //!<
                              float *ptrdistance, tup *ptrHitPoint, tup *ptrHitNormal, bool Print = false);

/**
 * Result of Ray towards Triangle hit.
 */
struct rti_result
{
  float t{};
  float u{};
  float v{};
  bool Hit {};
  tup P{};
};

auto RayTriangleIntersect(ray const &R, triangle const &Triangle, bool Print = false) -> rti_result;
};      // end namespace rm
};      // end namespace ww
#endif  // COMMON_SRC_MAIN_RAYMARCH_HPP
