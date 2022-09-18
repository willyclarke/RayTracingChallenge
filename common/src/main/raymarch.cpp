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
 * The point in world coordinates is moved into local coordinates by use of
 * the inverse transform of the shape.
 */
//------------------------------------------------------------------------------
float GetDistance(tup const &P, shared_ptr_shape PtrShape)
{
  tup const LocalPoint = ww::Inverse(PtrShape->Transform) * P;
  if (PtrShape->isA<sphere>())
  {
    ww::sphere const *pSphere = dynamic_cast<ww::sphere *>(PtrShape.get());
    float const Ds = SdfSphere(pSphere->Center, pSphere->Radius, LocalPoint);
    return Ds;
  }
  return {};
}

//------------------------------------------------------------------------------
tup GetNormal(tup const &P, shared_ptr_shape PtrShape)
{
  float constexpr E{0.01f};
  tup const EpsilonXYY{E, 0.f, 0.f, 0.f};
  tup const EpsilonYXY{0.f, E, 0.f, 0.f};
  tup const EpsilonYYX{0.f, 0.f, E, 0.f};

  tup const N{GetDistance(P + EpsilonXYY, PtrShape), GetDistance(P + EpsilonYXY, PtrShape),
              GetDistance(P + EpsilonYYX, PtrShape), 0.f};

  tup const WorldNormal = ww::Transpose(ww::Inverse(PtrShape->Transform)) * N;

  return Normalize(WorldNormal);
}

//------------------------------------------------------------------------------
float GetLight(light const &Light, tup const &P, shared_ptr_shape PtrShape = nullptr)
{
  tup const LightDir = Normalize(P - Light.Position);
  return -Dot(GetNormal(P, PtrShape), LightDir);
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
    // ---
    // NOTE: Compute the incremental position.
    // ---
    tup iPos = Ray.Origin + Ray.Direction * Distance;

    // ---
    // NOTE: Get the increment in distance to the shape of interest.
    //      This is called distance aided ray marching.
    //      A useful blogpost is this one:
    //      https://michaelwalczyk.com/blog-ray-marching.html
    // ---
    float const incrDistance = GetDistance(iPos, PtrShape);

    // ---
    // NOTE: Ray march to the new distance.
    // ---
    Distance += incrDistance;

    // std::cout << __FUNCTION__ << "-> Step: " << Idx << ". Distance: " << Distance << ". incrDistance: " <<
    // incrDistance
    //           << std::endl;

    if (Distance > ww::rm::MAX_DIST || incrDistance < ww::rm::MIN_DIST) break;
  }

  return Distance;
}

//------------------------------------------------------------------------------
tup MainImage(camera const &Camera, world const &World, int X, int Y, shared_ptr_shape PtrShape)
{
  ray const R = ww::RayForPixel(Camera, X, Y);
  float const Distance = RayMarch(R, PtrShape);
  tup fragColor{};

  if (Distance < MAX_DIST)  // then march along the ray, yoohoo ... ->->->
  {
    tup const pHit = R.Origin + R.Direction * Distance;
    fragColor = PtrShape->Material.Color * GetLight(*World.vPtrLights[0].get(), pHit, PtrShape);
  }

  return fragColor;
}

//------------------------------------------------------------------------------
canvas Render(camera const &Camera, world const &World)
{
  canvas Image(Camera.HSize, Camera.VSize);

  for (int X = 0; X < Image.W; ++X)
    for (int Y = 0; Y < Image.H; ++Y)
    {
      for (int ObjIdx = 0; ObjIdx < World.vPtrObjects.size(); ++ObjIdx)
      {
        tup const Color = MainImage(Camera, World, X, Y, World.vPtrObjects[ObjIdx]);
        Image.vXY[X + Image.W * Y] = Color + Image.vXY[X + Image.W * Y];
      }
    }
  return Image;
}

};  // end of namespace rm
};  // end of namespace ww
