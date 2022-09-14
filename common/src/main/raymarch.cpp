/******************************************************************************
 * * Filename : raymarch.cpp
 * * Date     : 2022 Sep 07
 * * Author   : Willy Clarke (willy@clarke.no)
 * * Version  : Use git you GIT
 * * Copyright: W. Clarke
 * * License  : MIT
 * * Descripti: Implementation of the ray marching algorithm.
 * ******************************************************************************/
#include "raymarch.hpp"

/**
 */
namespace ww
{

namespace rm
{
constexpr int MAX_STEPS = 100;
constexpr float MAX_DIST = 100.f;
constexpr float MIN_DIST = 0.001f;

//------------------------------------------------------------------------------
float SdfSphere(tup const &Center, float Radius, tup const &P)
{
  float const Distance = Mag(P - Center);
  float const DistToSphere = Distance - Radius;
  return DistToSphere;
}

/**
 * GetDistance to objects in the scene.
 * TODO: (Willy Clarke) : this function need to be configured to support all
 *                        objects in a scene.
 *                        For now there is an implicit sphere at origin with
 *                        y equal to sphere diameter.
 */
//------------------------------------------------------------------------------
float GetDistance(tup const &P, shared_ptr_shape PtrShape)
{
  if (PtrShape && PtrShape->isA<sphere>())
  {
    ww::sphere const *pSphere = dynamic_cast<ww::sphere *>(PtrShape.get());
    float const Ds = SdfSphere(pSphere->Center, pSphere->Radius, P);
    // std::cout << __FUNCTION__ << " -> is a sphere with Ds=" << Ds << std::endl;
    return Ds;
  }

  float Distance{};
  tup const Center{0.f, 0.f, 0.f, 0.f};
  float const Radius{0.4f};
  Distance = SdfSphere(Center, Radius, P);
  return Distance;
}

//------------------------------------------------------------------------------
tup GetNormal(tup const &P, shared_ptr_shape PtrShape)
{
  float const E{0.01f};
  tup const EpsilonXYY{E, 0.f, 0.f, 0.f};
  tup const EpsilonYXY{0.f, E, 0.f, 0.f};
  tup const EpsilonYYX{0.f, 0.f, E, 0.f};
  tup const N{GetDistance(P + EpsilonXYY, PtrShape), GetDistance(P + EpsilonYXY, PtrShape),
              GetDistance(P + EpsilonYYX, PtrShape), 0.f};
  return Normalize(N);
}

//------------------------------------------------------------------------------
float GetLight(tup const &P)
{
  tup const LightPos = tup{0.f, 3.f, -2.2f, 0.f};
  tup const LightDir = Normalize(P - LightPos);
  return -Dot(GetNormal(P), LightDir);
}

//------------------------------------------------------------------------------
float RayMarch(ray const &Ray, shared_ptr_shape PtrShape)
{
  float Distance{};
  for (int Idx = 0;          //!<
       Idx < rm::MAX_STEPS;  //!<
       ++Idx                 //!<
  )
  {
    tup iPos = Ray.Origin + Ray.Direction * Distance;

    float const incrDistance = GetDistance(iPos, PtrShape);
    Distance += incrDistance;

    // std::cout << __FUNCTION__ << "-> Step: " << Idx << ". Distance: " << Distance << ". incrDistance: " <<
    // incrDistance
    //           << std::endl;

    if (Distance > ww::rm::MAX_DIST || incrDistance < ww::rm::MIN_DIST) break;
  }
  return Distance;
}

//------------------------------------------------------------------------------
tup MainImage(int X, int Y, int W, int H, shared_ptr_shape PtrShape)
{
  tup UV{float(X) / float(W) - 0.5f, float(Y) / float(H) - 0.5f, 0.f, 0.f};
  UV.X *= float(W) / float(H);

  constexpr float FocalDistance{0.6f};
  ray const R = Ray(Point(0.f, 0.f, -1.6f), Vector(UV.X, UV.Y, FocalDistance));

  tup fragColor{};
  float const Distance = RayMarch(R, PtrShape);
  // std::cout << __FUNCTION__ << "-> Ray: " << R << ". Distance: " << Distance << std::endl;

  if (Distance < MAX_DIST)
  {
    tup const pHit = R.Origin + R.Direction * Distance;
    // std::cout << __FUNCTION__ << "HITHIT HIT ---- > X: " << X << " Y: " << Y << ". Distance: " << Distance
    //           << ". pHit: " << pHit << std::endl;
    fragColor = tup{0.5f, 0.2f, 0.6f, 0.f} * GetLight(pHit) + tup{0.9f, 0.1f, 0.1f, 0.f};
  }
  return fragColor;
}

//------------------------------------------------------------------------------
canvas Render(camera const &Camera, world const &World)
{
  canvas Image(Camera.HSize, Camera.VSize);

  // float const Resolution = Image.W * Image.H;
  for (int X = 0; X < Image.W; ++X)
    for (int Y = 0; Y < Image.H; ++Y)
    {
      for (int ObjIdx = 0; ObjIdx < World.vPtrObjects.size(); ++ObjIdx)
      {
        Image.vXY[X + Image.W * Y] = MainImage(X, Y, Image.W, Image.H, World.vPtrObjects[ObjIdx]);
      }
    }
  return Image;
}

};  // end of namespace rm
};  // end of namespace ww
