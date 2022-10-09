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
std::mutex gMutexPrint{};
constexpr int MAX_STEPS = 100;
constexpr float MAX_DIST = 100.f;
constexpr float MIN_DIST = 0.001f;
#if HW_PERFORMANCE == 0
#define AA 1
#else
#define AA 2  // make this 2 or 3 for antialiasing
#endif

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

//-
/**
 * sd Functions taken from Inigios Shadertoy primitives.
 */

//------------------------------------------------------------------------------
float sdBox(tup const &Pos, tup const &Box)
{
  tup Distance = Abs(Pos) - Box;
  float const Result = std::min(std::max(Distance.X, std::max(Distance.Y, Distance.Z)),  //!<
                                0.f)                                                     //!<
                       + Mag(Max(Distance, 0.f));
  // std::cout << __FUNCTION__ << ". Distance: " << Distance << ". Lenght: " << Mag(Distance) << ". Result: " << Result
  //           << std::endl;

  return Result;
}

//------------------------------------------------------------------------------
float sdSphere(tup const &Pos, float S)
{
  Assert(IsPoint(Pos), __FUNCTION__, __LINE__);
  return Mag(Pos - Point(0.f, 0.f, 0.f)) - S;
}

//------------------------------------------------------------------------------
// la,lb=semi axis, h=height, ra=corner
float sdRhombus(tup Pos, float la, float lb, float h, float ra)
{
  Pos = Abs(Pos);
  tup B = Vector(la, lb, 0.f);
  float F = Clamp(NDot(B, B - 2.f * tup(Pos.X, Pos.Z)) / Dot(B, B), -1.f, 1.f);
  tup Q = tup(Mag(tup(Pos.X, Pos.Z) - 0.5f * B * tup(1.f - F, 1.f + F))  //!<
                      * Sign(Pos.X * B.Y + Pos.Z * B.X - B.X * B.Y) -
                  ra,  //!<
              Pos.Y - h);
  float Result = std::fmin(std::fmax(Q.X, Q.Y), 0.f) + Mag(Max(Q, 0.f));
  return Result;
}

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
tup OpU(tup const &D1, tup const &D2) { return (D1.X < D2.X) ? D1 : D2; }

//------------------------------------------------------------------------------
/**
 * Map the various Signed Distance Functions of the objects in the scene.
 */
tup Map(tup const &Pos)
{
  Assert(IsPoint(Pos), __FUNCTION__, __LINE__);

  tup Res = Point(Pos.Y, 0.f, 0.f);

  // ---
  // Bounding box.
  // ---
  if (sdBox(Pos - Point(-2.f, 0.3f, 0.25f), Vector(0.3f, 0.3f, 1.f)) < Res.X)
  {
    Res = OpU(Res, Point(sdSphere(Pos - Vector(-2.f, 0.25f, 0.f), 0.25f), 26.9f, 0.f));
    Res = OpU(Res, Point(sdRhombus(Pos - Vector(-2.f, 0.25f, 1.f), 0.15f, 0.25f, 0.04f, 0.08f), 17.f, 0.f));
  }
  return Res;
}

//------------------------------------------------------------------------------
// https://iquilezles.org/articles/rmshadows
float CalcSoftShadow(ray const &R, float tMin, float tMax)
{
  // bounding volume
  float tp = (0.8 - R.Origin.Y) / R.Direction.Y;
  if (tp > 0.f) tMax = std::min(tMax, tp);

  float Res = 1.0;
  float t = tMin;
  for (int Idx = 0; Idx < 24; ++Idx)
  {
    float H = Map(R.Origin + R.Direction * t).X;
    float S = Clamp(8.f * H / t, 0.f, 1.f);
    Res = std::min(Res, S);
    t += Clamp(H, 0.01f, 0.2f);
    if (Res < 0.004f || t > tMax) break;
  }
  Res = Clamp(Res, 0.f, 1.f);
  return Res * Res * (3.f - 2.f * Res);
}

//------------------------------------------------------------------------------
// https://iquilezles.org/articles/normalsSDF
tup CalcNormal(tup const &Pos)
{
#if 1
  tup const e = Vector(1.f, -1.f, 0.f) * 0.5773f * 0.0005f;
  tup const exyy = Vector(e.X, e.Y, e.Y);
  tup const eyyx = Vector(e.Y, e.Y, e.X);
  tup const eyxy = Vector(e.Y, e.X, e.Y);
  tup const exxx = Vector(e.X, e.X, e.X);
  tup const N = Normalize(exyy * Map(Pos + exyy).X +  //!<
                          eyyx * Map(Pos + eyyx).X +  //!<
                          eyxy * Map(Pos + eyxy).X +  //!<
                          exxx * Map(Pos + exxx).X    //!<
  );
  Assert(IsVector(N), __FUNCTION__, __LINE__);
  return N;
#else
  // inspired by tdhooper and klems - a way to prevent the compiler from inlining map() 4 times
  vec3 n = vec3(0.0);
  for (int i = ZERO; i < 4; i++)
  {
    vec3 e = 0.5773 * (2.0 * vec3((((i + 3) >> 1) & 1), ((i >> 1) & 1), (i & 1)) - 1.0);
    n += e * map(pos + 0.0005 * e).x;
    // if( n.x+n.y+n.z>100.0 ) break;
  }
  return normalize(n);
#endif
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

// https://iquilezles.org/articles/boxfunctions
tup iBox(ray const &R, tup const &Rad)
{
  tup const M = 1.f / R.Direction;  // is a Vector.
  tup const N = M * R.Origin;       // becomes a Vector.
  tup const K = Abs(M) * Rad;       // becomes a Vector since M is a Vector.
  tup const t1 = -N - K;
  tup const t2 = -N + K;
  return tup(std::max(std::max(t1.X, t1.Y), t1.Z), std::min(std::min(t2.X, t2.Y), t2.Z));
}

//------------------------------------------------------------------------------
/**
 * Raycast from Inigios primitive's example.
 */
tup RayCast(ray const &R)
{
  tup Res = Vector(-1.f, -1.f, 0.f);
  float Tmin{1.f};
  float Tmax{20.f};

  // ---
  // Raytrace floor plane
  // ---
  float Tp1 = (0.f - R.Origin.Y) / R.Direction.Y;
  if (Tp1 > 0.f)
  {
    Tmax = std::min(Tmax, Tp1);
    Res = Vector(Tp1, 1.f, 0.f);
  }

  //
  // Raymarch primitives
  //
  tup Tb = iBox(Ray(R.Origin - Vector(0.f, 0.4f, -0.5f), R.Direction), Point(2.5f, 0.41f, 3.f));
  if (Tb.X < Tb.Y && Tb.Y > 0.f && Tb.X < Tmax)
  {
    Tmin = std::max(Tb.X, Tmin);
    Tmax = std::min(Tb.Y, Tmax);

    float t = Tmin;
    for (int Idx = 0; Idx < 70 && t < Tmax; ++Idx)
    {
      tup h = Map(R.Origin + R.Direction * t);
      if (std::abs(h.X) < (0.0001f * t))
      {
        Res = Vector(t, h.Y, 0.f);
        break;
      }
      t += h.X;
    }
  }

  return Res;
}

//------------------------------------------------------------------------------
// https://iquilezles.org/articles/rmshadows
float CalcSoftshadow(ray const &R, float mint, float tmax)
{
  // bounding volume
  float tp = (0.8f - R.Origin.Y) / R.Direction.Y;
  if (tp > 0.f) tmax = std::min(tmax, tp);

  float Res = 1.f;
  float t = mint;
  for (int Idx = 0; Idx < 24; ++Idx)
  {
    float h = Map(R.Origin + R.Direction * t).X;
    float s = Clamp(8.f * h / t, 0.f, 1.f);
    Res = std::fmin(Res, s);
    t += Clamp(h, 0.01f, 0.2f);
    if (Res < 0.004f || t > tmax) break;
  }
  Res = Clamp(Res, 0.f, 1.f);
  return Res * Res * (3.f - 2.f * Res);
}

//------------------------------------------------------------------------------
//
// https://iquilezles.org/articles/nvscene2008/rwwtt.pdf
//
//------------------------------------------------------------------------------
/**
 * Calculate the Ambient Occlusion factor.
 * @Pos: Position.
 * @Nor: Normal Vector.
 * @N: Number of samples - optional.
 * @return: float describing the Occlusion after N samples.
 */
float CalcAO(tup const &Pos, tup const &Nor, int N = 5)
{
  float Occ{};
  float Sca{1.f};
  for (int Idx = 0;  //!<
       Idx < N;      //!<
       ++Idx)
  {
    float H = 0.01f + 0.12f * float(Idx) / 4.f;
    float D = Map(Pos + H * Nor).X;
    Occ += (H - D) * Sca;
    Sca *= 0.95f;
    if (Occ > 0.35f) break;
  }
  return Clamp(1.f - 3.f * Occ, 0.f, 1.f) * (0.5f + 0.5 * Nor.Y);
}

//------------------------------------------------------------------------------
// https://iquilezles.org/articles/checkerfiltering
float CheckersGradBox(tup const &Pos, tup const &dPdx, tup const &dPdy)
{
  // filter kernel
  tup W = Abs(dPdx) + Abs(dPdy) + Color(0.001f, 0.001f, 0.f);
  // analytical integral (box filter)
  tup I = 2.f *
          (Abs(Fract((Pos - 0.5f * W) * 0.5f) - Color(0.5f, 0.5f, 0.f)) -
           Abs(Fract((Pos + 0.5f * W) * 0.5f) - Color(0.5f, 0.5f, 0.0f))) /
          W;
  // xor pattern
  return 0.5f - 0.5f * I.X * I.Y;
}

//------------------------------------------------------------------------------
/**
 * Render from Inigios primitive's example.
 */
tup Render(ray const &R, tup const &Rdx, tup const &Rdy)
{
  // Background
  tup Col = Color(0.7f, 0.7f, 0.9f) - 0.3f * Color(std::max(R.Direction.Y, 0.f),  //!<
                                                   std::max(R.Direction.Y, 0.f),  //!<
                                                   std::max(R.Direction.Y, 0.f)   //!<
                                             );

  // ---
  // Raycast scene
  // ---

  tup Res = RayCast(R);
  float t = Res.X;
  float M = Res.Y;

  if (M > -0.5f)
  {
    tup Pos = R.Origin + t * R.Direction;

    tup Nor = (M < 1.5f) ? Vector(0.f, 1.f, 0.f) : CalcNormal(Pos);
    tup Ref = Reflect(R.Direction, Nor);

    // Material
    Col = Color(0.2f, 0.2f, 0.2f) + 0.2f * Sin(2.f * Color(M, M, M) + Color(0.f, 1.f, 2.f));
    float Ks = 1.f;

    if (M < 1.5f)
    {
      // Project pixel footprint into the plane
      tup dPdx = R.Origin.Y * (R.Direction / R.Direction.Y - Rdx / Rdx.Y);
      tup dPdy = R.Origin.Y * (R.Direction / R.Direction.Y - Rdy / Rdy.Y);

      float Floor = CheckersGradBox(3.0 * Vector(Pos.X, Pos.Z, 0.f),    //!<
                                    3.0 * Vector(dPdx.X, dPdx.Z, 0.f),  //!<
                                    3.0 * Vector(dPdy.X, dPdy.Z, 0.f)   //!<
      );
      Col = Color(0.15f, 0.15f, 0.15f) + Floor * Color(0.05f, 0.05f, 0.05f);
      Ks = 0.4;
    }

    // ---
    // Lighting.
    // ---

    // ---
    // Calculate Ambient Occlusion.
    // ---
    float Occ = CalcAO(Pos, Nor);

    tup Lin{};

    // ---
    // Sun
    // ---
    if (1)
    {
      tup Lig = Normalize(Vector(-0.5f, 0.4, -0.6f));
      tup Hal = Normalize(Lig - R.Direction);
      float Dif = Clamp(Dot(Nor, Lig), 0.f, 1.f);

      Dif *= CalcSoftShadow(Ray(Pos, Lig), 0.02f, 2.5f);
      float Spe = std::powf(Clamp(Dot(Nor, Hal), 0.f, 1.f), 16.f);
      Spe *= Dif;
      Spe *= 0.04f + 0.96f * std::powf(Clamp(1.f - Dot(Hal, Lig), 0.f, 1.f), 5.f);
      Lin = Lin + Col * 2.2f * Dif * Vector(1.3f, 1.f, 0.7f);
      Lin = Lin + 5.f * Spe * Vector(1.3f, 1.f, 0.7f) * Ks;
    }

    // ---
    // Sky
    // ---
    if (1)
    {
      float Dif = std::sqrtf(Clamp(0.5f + 0.5f * Nor.Y, 0.f, 1.f));
      Dif *= Occ;
      float Spe = SmoothStep(-0.2f, 0.2f, Ref.Y);
      Spe *= Dif;
      Spe *= 0.04f + 0.96f * std::powf(Clamp(1.f + Dot(Nor, R.Direction), 0.f, 1.f), 5.f);
      // if( spe>0.001 )
      Spe *= CalcSoftshadow(Ray(Pos, Ref), 0.02, 2.5);
      Lin = Lin + Col * 0.60 * Dif * Vector(0.4f, 0.6f, 1.15f);
      Lin = Lin + 2.00 * Spe * Vector(0.4f, 0.6f, 1.3f) * Ks;
    }

    Col = Lin;
  }

  return Clamp(Col, 0.f, 1.f);
}

//------------------------------------------------------------------------------
/**
 * SetCamera
 */
matrix SetCamera(tup const &Origin, tup const &Ta, float const Cr)
{
  tup const Cw = Normalize(Ta - Origin);
  tup const Cp = Vector(std::sinf(Cr), std::cosf(Cr), 0.f);
  tup const Cu = Normalize(Cross(Cw, Cp));
  tup const Cv = Cross(Cu, Cw);
  matrix M{};
  M.R0 = Cu;
  M.R1 = Cv;
  M.R2 = Cw;

  return M;
}

//------------------------------------------------------------------------------
/**
 * MainImage: https://www.shadertoy.com/view/Xds3zN
 * Raymarching - Primitives.
 * A set of raw primitives. All except the ellipsoid are exact euclidean distances.
 * Based on the work by iq in 2013-03-25.
 * @param: FragCoord.X - Horisontal coordinate of the screen.
 * @param: FragCoord.Y - Vertical coordinate of the screen.
 * @param: Resolution.X - Horisontal resolution.
 * @param: Resolution.Y - Vertical resolution.
 * @return: FragColor - The color of the pixel at X,Y.
 */
tup MainImage(tup const &FragCoord, mainimage_config const &Cfg)
{
  ray R{};
  float const Time{};
  tup const Mouse{};

  // ---
  // Camera:
  // ---
  // tup const Ta = Vector(0.25f, -0.75f, -0.75f);
  tup const Ta = Vector(0.f, -1.5f, 0.f);
  R.Origin = Ta + Point(4.5f * std::cosf(0.1f * Time + 7.f * Mouse.X),  //!<
                        2.2f,                                           //!<
                        4.5f * std::sinf(0.1f * Time + 7.f * Mouse.X)   //!<
                  );

  // Camera to World transformation
  // TODO: (Willy Clarke) Move the camera instantiation.
  // static matrix const Ca = SetCamera(R.Origin, Ta, 0.f);
  // static matrix const Ca = TranslateScaleRotate(0.f, 0.f, 0.f, 1.f, 1.f, 1.f, M_PI, -2.f * 0.78f, 0.f);
  matrix const &Ca = Cfg.MCamera;

  tup Tot{};

  // ---
  // Pixel coordinates. NOTE: This becomes a vector.
  // ---
  tup PixelCoord = (2.f * FragCoord - Cfg.Resolution) / Cfg.Resolution.Y;
  PixelCoord.Z = Cfg.FocalLength;

  // ---
  // Ray direction.
  // ---
  R.Direction = Ca * Normalize(PixelCoord);

  // ---
  // Ray differentials
  // ---
  tup Px = (2.f * (FragCoord + Vector(1.f, 0.f, 0.f)) - Cfg.Resolution) / Cfg.Resolution.Y;  // is a vector.
  Px.Z = Cfg.FocalLength;
  Assert(IsVector(Px), __FUNCTION__, __LINE__);

  tup Py = (2.f * (FragCoord + Vector(0.f, 1.f, 0.f)) - Cfg.Resolution) / Cfg.Resolution.Y;  // is a vector.
  Py.Z = Cfg.FocalLength;
  Assert(IsVector(Py), __FUNCTION__, __LINE__);

  tup Rdx = Ca * Normalize(Px);
  tup Rdy = Ca * Normalize(Py);

  // ---
  // Render
  // ---
  tup Color = Render(R, Rdx, Rdy);

  // ---
  // Gain
  // ---
  Color = 3.f * Color / (Vector(2.5f, 2.5f, 2.5f) + Color);

  // ---
  // Gamma.
  // ---
  Color = Pow(Color, Vector(0.4545f, 0.4545f, 0.4545f));

  Tot = Tot + Color;
  return Tot;
}

//------------------------------------------------------------------------------
canvas RenderSingleThread(camera const &Camera, world const &World)
{
  canvas Image(Camera.HSize, Camera.VSize);
  mainimage_config Cfg{};
  Cfg.Resolution = Point(Image.W, Image.H, 0.f);
  Cfg.MCamera = TranslateScaleRotate(0.f, 0.f, 0.f, 1.f, 1.f, 1.f, M_PI, -2.f * 0.78f, 0.f);

  for (int X = 0; X < Image.W; ++X)
    for (int Y = 0; Y < Image.H; ++Y)
    {
      for (int ObjIdx = 0; ObjIdx < World.vPtrObjects.size(); ++ObjIdx)
      {
#if 0
        tup const Color = MainImage(Camera, World, X, Y, World.vPtrObjects[ObjIdx]);
        Image.vXY[X + Image.W * Y] = Color + Image.vXY[X + Image.W * Y];
#endif
        tup const FragCoord = Point(X, Y, 0.f);
        Image.vXY[X + Image.W * Y] = MainImage(FragCoord, Cfg);
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
  // world const &World = *RB.ptrWorld;
  canvas &Image = *RB.ptrImage;
  camera const &Camera = *RB.ptrCamera;
  tup const Resolution = Point(Camera.HSize, Camera.VSize, 0.f);

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
#if 0
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
#endif
        tup const FragCoord = Point(X, Y, 0.f);
        Image.vXY[X + Image.W * Y] = MainImage(FragCoord, RB.Cfg);
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
  mainimage_config Cfg{};
  Cfg.Resolution = Point(Image.W, Image.H, 0.f);
  Cfg.MCamera = TranslateScaleRotate(0.f, 0.f, 0.f, 1.f, 1.f, 1.f, M_PI, -2.f * 0.78f, 0.f);

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
      RB.Cfg = Cfg;
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
