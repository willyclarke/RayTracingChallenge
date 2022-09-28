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

#include <cmath>
#include <thread>
#include <vector>

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
/**
 * @param: P - Point of interest.
 * @param: B - Box, tuple with dimensions for X, Y, Z.
 * @param: R - Radius for box with rounded corners.
 */
float SdfBox(tup const &P, tup const &B, float const Radius = 0.f)
{
  tup const Q = Abs(P) - B;
  float const Distance = Mag(Max(Q, 0.f)) + std::min(std::max(Q.X, std::max(Q.Y, Q.Z)), 0.f) - Radius;

  // std::cout << "P:\n" << P << "\nBox:\n" << B << "\nDistance:" << Distance << std::endl;

  return Distance;
}

//------------------------------------------------------------------------------
/**
 * @param: P - Point
 * @param: N - Normal (which must be normalized)
 * @param: H - Height to lower the floor. i.e the distance increase by H.
 */
float SdfPlane(tup const &P, tup const &N, float H = 0.f)
{
  float const Result = Dot(P, N) + H;
  return Result;
}

//------------------------------------------------------------------------------
float SdfSphere(tup const &Center, float Radius, tup const &P)
{
  float const Distance = Mag(P - Center);
  float const DistToSphere = Distance - Radius;
  return DistToSphere;
  float const DistToSphere2 = P.Y + Radius;
  return std::min(DistToSphere, DistToSphere2);
}

/**
 * Description: clamp returns the value of x constrained to the range minVal to maxVal.
 * The returned value is computed as min(max(x, minVal), maxVal).
 */
float Clamp(float X, float MinVal, float MaxVal) { return std::min(std::max(X, MinVal), MaxVal); }

/**
 * pow — return the value of the first parameter raised to the power of the second
 * @param: tup X - Specify the value to raise to the power Y.
 * @param: tup Y - Specify the power to which to raise X.
 * @return: tup that maitainse type beeing either Point or Vector. i.e. X.W remains unchanged.
 */
tup Pow(tup const &X, tup const &Y) { return tup{std::powf(X.X, Y.X), std::powf(X.Y, Y.Y), std::powf(X.Z, Y.Z), X.W}; }

/**
 * Description: step generates a step function by comparing x to edge.
 * For element i of the return value, 0.0 is returned if x[i] < edge[i],
 * and 1.0 is returned otherwise.
 */
float Step(float Edge, float X)
{
  if (X < Edge) return 0.f;
  return 1.f;
}

/**
 * GetDistance to one particular object in the scene.
 * The point in world coordinates is moved into local coordinates by use of
 * the inverse transform of the shape.
 * return: Distance to shape or MAX_DIST when shape is not supported.
 */
//------------------------------------------------------------------------------
float GetDistance(tup const &P, shared_ptr_shape PtrShape)
{
  tup const LocalPoint = ww::Inverse(PtrShape->Transform) * P;
  if (PtrShape->isA<sphere>())
  {
    sphere const *pSphere = dynamic_cast<sphere *>(PtrShape.get());
    // ---
    // NOTE: A sphere will have uniform scaling in all directions.
    //       So; for now the scaled X will be used as radius.
    // ---
    float const Radius = pSphere->Transform.R0.X;
    float const Ds = SdfSphere(pSphere->Center, Radius, LocalPoint);
    return Ds;
  }
  else if (PtrShape->isA<cube>())
  {
    cube const *pCube = dynamic_cast<cube *>(PtrShape.get());
    tup const Box = Point(pCube->Transform.R[0].X / 2.f, pCube->Transform.R[1].Y / 2.f, pCube->Transform.R[2].Z / 2.f);
    float const Db = SdfBox(LocalPoint, Box, pCube->R);
    return Db;
  }
  else if (PtrShape->isA<plane>())
  {
    plane const *pPlane = dynamic_cast<plane *>(PtrShape.get());
    float const Dp = SdfPlane(LocalPoint, Point(0.f, 1.f, 0.f), pPlane->H);
    return Dp;
  }
  return MAX_DIST;
}

//------------------------------------------------------------------------------
float GetDistance(tup const &P, world const &World)
{
  float Distance{MAX_DIST};
  // ---
  // NOTE: Set up scene here.
  // ---
  for (int ObjIdx = 0; ObjIdx < World.vPtrObjects.size(); ++ObjIdx)
  {
    shared_ptr_shape PtrShape = World.vPtrObjects[ObjIdx];
    Distance = std::min(GetDistance(P, PtrShape), Distance);
  }

  return Distance;
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
tup GetLight(tup const &P, shared_ptr_shape PtrShape)
{
  // ---
  // NOTE: Test code:
  // ---
  // "normal white" material should be around 0.2 gray
  static tup const Mate = Color(0.2f, 0.2f, 0.2f);
  tup const Nor = GetNormal(P, PtrShape);

  // Lighting
  static tup const SunDir = Normalize(Vector(0.8f, 0.4f, 0.2f));
  float const SunDif = Clamp(Dot(Nor, SunDir), 0.f, 1.f);
  float const SunSha = Step(RayMarch(Ray(P + Nor * 0.001f, SunDir), PtrShape), 0.f);
  float const SkyDif = Clamp(0.5f + 0.5f * Dot(Nor, Vector(0.f, 1.f, 0.f)), 0.f, 1.f);
  float const BouDif = Clamp(0.5f + 0.5f * Dot(Nor, Vector(0.f, -1.f, 0.f)), 0.f, 1.f);

  // ---
  // NOTE: Key light ~- 10.
  //       Fill light ~- 1.
  tup const Color0 = Mate * Point(7.f, 5.f, 3.f) * SunDif * SunSha;
  tup const Color1 = Mate * Point(0.5f, 0.8f, 0.9f) * SkyDif;
  tup const Color2 = Mate * Point(0.7f, 0.3f, 0.2f) * BouDif;
  tup const Color = Color0 + Color1 + Color2;
  return Color;
}

//------------------------------------------------------------------------------
float GetLight(light const &Light, tup const &P, shared_ptr_shape PtrShape)
{
  // ---
  // Original code
  // ---
  tup const LightDir = Normalize(P - Light.Position);
  return -Dot(GetNormal(P, PtrShape), LightDir);
}

//------------------------------------------------------------------------------
tup Lighting(material const &Material, light const &Light, tup const &P, tup const &EyeV, tup const NormalV)
{
  tup const EffectiveColor = Material.Color * Light.Intensity;
  tup const LightV = Normalize(Light.Position - P);
  tup const Ambient = EffectiveColor * Material.Ambient;
  float const LightDotNormal = Dot(LightV, NormalV);

  tup Diffuse{};
  tup Specular{};

  if (LightDotNormal < 0.f)
  {
    return Ambient;
  }
  // std::cout << __FUNCTION__ << ". LightDotNormal: " << LightDotNormal << ". Ambient:" << Ambient << ". Intensity:" <<
  // Light.Intensity << std::endl;

  Diffuse = EffectiveColor * Material.Diffuse * LightDotNormal;

  tup const ReflectV = Reflect(-LightV, NormalV);
  float const ReflectDotEye = Dot(ReflectV, EyeV);
  if (ReflectDotEye > 0.f)  // a negative number means pointing away.
  {
    float const Factor = std::pow(ReflectDotEye, Material.Shininess);
    Specular = Light.Intensity * Material.Specular * Factor;
  }
  return Ambient + Diffuse + Specular;
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
float RayMarch(ray const &Ray, world const &World, shared_ptr_shape PtrShape) { return {}; }

//------------------------------------------------------------------------------
tup MainImage(camera const &Camera, world const &World, int X, int Y, shared_ptr_shape PtrShape)
{
  ray const R = ww::RayForPixel(Camera, X, Y);
  float const Distance = RayMarch(R, PtrShape);

  shared_ptr_light PtrLights = World.vPtrLights[0];
  light const &Light = *PtrLights;

  // tup fragColor = Color(0.65f, 0.75f, 0.9f) * 0.2;
  tup fragColor{};
  if (Distance < MAX_DIST)  // then there is a hit
  {
    tup const pHit = R.Origin + R.Direction * Distance;
    fragColor = PtrShape->Material.Color * GetLight(Light, pHit, PtrShape);
    // fragColor = GetLight(pHit, PtrShape);
    // std::cout << __FUNCTION__ << ": Got hit at X:" << X << " Y:" << Y << ". Color:" << fragColor << std::endl;
  }

  // fragColor = Pow(fragColor, Color(0.4545f, 0.4545f, 0.4545f));
  return fragColor;
}

//------------------------------------------------------------------------------
canvas RenderSingleThread(camera const &Camera, world const &World)
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

/**
 * Render a block of pixels defined in the struct render_block.
 * This function locks on a mutex defined in the render block.
 */
void RenderBlock(render_block const &RB)
{
  for (int Y = RB.VStart;           ///<!
       Y < RB.VStart + RB.VHeigth;  ///<!
       ++Y)
  {
    for (int X = RB.HStart;           ///<!
         X < RB.HStart + RB.HLength;  ///<!
         ++X)
    {
      for (int ObjIdx = 0; ObjIdx < RB.ptrWorld->vPtrObjects.size(); ++ObjIdx)
      {
        world const &World = *RB.ptrWorld;
        canvas &Image = *RB.ptrImage;
        camera const &Camera = *RB.ptrCamera;
        shared_ptr_shape PtrShape = World.vPtrObjects[ObjIdx];

        if (Image.vXY[X + Image.W * Y].W < 1.f)
        {
          tup const Color = MainImage(Camera, World, X, Y, PtrShape);
          if (Mag(Color) > 0.001f)  // then there is an object there
          {
            Image.vXY[X + Image.W * Y] = Color;  //+ Image.vXY[X + Image.W * Y];
            Image.vXY[X + Image.W * Y].W = 1.f;  // store the update
          }
        }
        // Image.vXY[X + Image.W * Y] = MainImage(Camera, World, X, Y, PtrShape);
      }
    }
  }
}

/**
 * Use a number of threads to render the image on to the canvas.
 * NOTE: The function is not thread safe but since the canvas is
 *       split into blocks that are not overlapping there should (?)
 *       not be any undefinded behavior.
 */
void RenderMultiThread(camera const &Camera, world const &World, canvas &Image)
{
  int const NumBlocksH = Camera.NumBlocksH;
  int const NumBlocksV = Camera.NumBlocksV;
  int const HSizeBlock = Camera.HSize / NumBlocksH;
  int const VSizeBlock = Camera.VSize / NumBlocksV;

  std::vector<render_block> vRenderBlocks{};

  for (int HIdx = 0;       //!<
       HIdx < NumBlocksH;  //!<
       ++HIdx)
  {
    for (int VIdx = 0;       //!<
         VIdx < NumBlocksV;  //!<
         ++VIdx)
    {
      render_block RB{};
      RB.HLength = HSizeBlock;
      RB.VHeigth = VSizeBlock;
      RB.HStart = HIdx * HSizeBlock;
      RB.VStart = VIdx * VSizeBlock;
      RB.ptrImage = &Image;
      RB.ptrCamera = &Camera;
      RB.ptrWorld = &World;
      vRenderBlocks.push_back(RB);
    }
  }

  struct WorkerThread
  {
    std::thread T;
    render_block _RB{};
    WorkerThread() { std::cout << __PRETTY_FUNCTION__ << " -> Called XXXXXXXXXXXXXXXXXXXXXXXXXX " << std::endl; }
    WorkerThread(render_block const &RB) { _RB = RB; }

    void Start() { T = std::thread(RenderBlock, _RB); }

    WorkerThread(WorkerThread const &Orig) { _RB = Orig._RB; };
    ~WorkerThread()
    {
      if (T.joinable())
      {
        T.join();
      }
    };
  };

  std::vector<WorkerThread> vWorkerThreads{};

  for (auto const &RB : vRenderBlocks)
  {
    vWorkerThreads.push_back(RB);
  }

  for (auto &T : vWorkerThreads)
  {
    T.Start();
  }
}

//------------------------------------------------------------------------------
canvas Render(camera const &Camera, world const &World)
{
  Assert(World.vPtrLights.size() > 0, __FUNCTION__, __LINE__);

  if (Camera.RenderSingleThread) return RenderSingleThread(Camera, World);

  canvas Image(Camera.HSize, Camera.VSize);
  RenderMultiThread(Camera, World, Image);
  return Image;
}

};  // end of namespace rm
};  // end of namespace ww
