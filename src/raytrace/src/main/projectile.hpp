/******************************************************************************
 * Filename : projectile.hpp
 * Date     : 2018 Sep 01
 * Author   : Willy Clarke (willy@clarke.no)
 * Version  : 0.0.1
 * Copyright: W. Clarke
 * License  : MIT
 * Descripti: Test code for Putting it together part of the first chapter
 *          : in the book 'The Ray Tracer Challenge'
 ******************************************************************************/
#ifndef SRC_RAYTRACE_SRC_MAIN_PROJECTILE_HPP
#define SRC_RAYTRACE_SRC_MAIN_PROJECTILE_HPP

#include "gtest/gtest.h"

#include <datastructures.hpp>

namespace
{
struct projectile
{
  ww::tup Position{ww::Point(0.f, 0.f, 0.f)};   //<! Point.
  ww::tup Velocity{ww::Vector(0.f, 0.f, 0.f)};  //<! Vector.
};

struct world
{
  ww::tup Gravity{ww::Vector(0.f, 0.f, 9.81f)};  //<! Vector.
  ww::tup Wind{ww::Vector(1.0f, 0.f, 0.f)};      //<! Vector.
};

projectile Tick(world const &World, projectile const &P)
{
  projectile Result{};
  Result.Position = P.Position + P.Velocity;
  Result.Velocity = P.Velocity + World.Gravity + World.Wind;
  return (Result);
}

};  // end of anonymous namespace

namespace rtcch1
{
void RunProjectileTest(int argc, char *argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  int const Result = RUN_ALL_TESTS();
  std::cout << "Testing result was " << (Result ? "Failure." : "Success.") << std::endl;

  // ---
  // NOTE: Projectile starts one unit above the origin.
  //     : Velocity is normalized to one unit per tick.
  // ---
  projectile P{ww::Point(0.f, 1.f, 0.f), ww::Normal(ww::Vector(1.f, 1.f, 0.f))};
  // ---
  // NOTE: World gravity of -0.1 unit/tick and wind is 0.01 unit/tick.
  // ---
  world const World{ww::Vector(0.f, -0.01f, 0.f), ww::Vector(-0.01f, 0.f, 0.f)};
  bool Stop{};
  int Count{};
  while ((P.Position.Y > 0.f) && !Stop)
  {
    P = Tick(World, P);
    std::cout << P.Position << std::endl;

    if (Count++ > 10000)
    {
      Stop = true;
    }
  }
}
};  // namespace rtcch1
namespace rtcch2
{
void RunProjectileTest()
{
  ww::tup Start = ww::Point(0.f, 1.f, 0.f);
  ww::tup Velocity = ww::Normal(ww::Vector(1., 1.8f, 0.f)) * 11.25f;

  projectile P{Start, Velocity};

  ww::tup Gravity = ww::Vector(0.f, -0.1f, 0.f);
  ww::tup Wind = ww::Vector(-0.01f, 0.f, 0.f);

  world World{Gravity, Wind};

  ww::canvas Canvas(900, 550);
  ww::tup Color = ww::Color(1.f, 0.f, 0.f);

  bool Stop{};
  int Count{};

  // ---
  // So we need to scale from our world to pixel.
  // Set up pixels per meter.
  // ---
  float const HPixelMeter = Canvas.W / 2000.f;  // use 2000m across.
  auto HPixel = [&](float const Pos) -> int {
    // Pos = 0   -> Pixel = Canvas.W/2
    // Pos = -10 -> Pixel = 0
    // Pos = +10 -> Pixel = Canvas.W
    float const P = HPixelMeter * Pos + Canvas.W / 2.f;
    int const Result = std::min<int>(Canvas.W - 1, std::max<int>(0, int(P)));
    return Result;
  };

  // The vertical axis need to be reversed.
  float const VPixelMeter = -Canvas.H / 2000.f;  // use 2000m across.
  auto VPixel = [&](float const Pos) -> int {
    // Pos y = 0    -> Pixel = Canvas.H
    // Pos y = 2000 -> Pixel = 0
    // Pos y = 1000 -> Pixel = Canvas.H/2
    float const P = VPixelMeter * Pos + Canvas.H;
    int const Result = std::min<int>(Canvas.H - 1, std::max<int>(0, int(P)));
    std::cout << "P:" << P << ". Res:" << Result << std::endl;
    return Result;
  };

  while ((P.Position.Y > 0.f) && !Stop)
  {
    P = Tick(World, P);
    std::cout << "X->" << P.Position << std::endl;
    ww::WritePixel(Canvas, HPixel(P.Position.X), VPixel(P.Position.Y), Color);

    if (Count++ > 10000)
    {
      Stop = true;
    }
  }
  // Finally we write the file
  ww::WriteToPPM(Canvas, "Projectile.ppm");
}

};  // namespace rtcch2
#endif
