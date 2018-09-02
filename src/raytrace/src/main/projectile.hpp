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

namespace rtcch1
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
#endif
