/******************************************************************************
 * Filename : main.cpp
 * Date     : 2018 Aug 25
 * Author   : Willy Clarke (willy@clarke.no)
 * Version  : 0.0.1
 * Copyright: W. Clarke
 * License  : MIT
 * Descripti: Test program for my Raytracing challenge.
 ******************************************************************************/
#include <datastructures.hpp>

#include <iomanip>  // for setw().
#include <iostream>

namespace
{  // anonymous

// ---
// NOTE: Stream operator for this file only: prints a tuple
// ---
std::ostream &operator<<(std::ostream &stream, const ww::tup &T)
{
   stream << ((T.W != 0) ? "Point:" : "Vector:");
   stream << std::setw(4) << T.X         //<!
          << " " << std::setw(4) << T.Y  //<!
          << " " << std::setw(4) << T.Z  //<!
          << " " << std::setw(4) << T.W;
   return stream;
}
};  // namespace

// ---
// NOTE: Test function for the tuples.
//     : Prints out expected results for some simple tests define in the
//     : raytracing book.
// ---
void RunTupleTest()
{
   {
      ww::tup A = ww::Vector(0.f, 0.f, 0.f);
      bool const AIsPoint = ww::IsPoint(A);
      if (!AIsPoint)
      {
         std::cout << "SUCCESS: Tuple a is not a point." << A << "." << std::endl;
      }
      else
      {
         std::cerr << "FAILURE: Tuple A is point. " << A << "." << std::endl;
      }
   }

   {
      ww::tup B = ww::Point(4.3, -4.2, 3.1);
      bool const BIsAPoint = ww::IsPoint(B);
      if (BIsAPoint)
      {
         std::cout << "SUCCESS: B is a Point" << std::endl;
      }
      else
      {
         std::cerr << "FAILURE: B is not a Point" << std::endl;
      }

      bool const BIsAVector = ww::IsVector(B);
      if (!BIsAVector)
      {
         std::cout << "SUCCESS: B is not a Vector" << std::endl;
      }
      else
      {
         std::cerr << "FAILURE: B is a Vector" << std::endl;
      }

      std::cout << B << std::endl;
   }
   {
      ww::tup B = ww::Vector(4, -4, 3);
      if (B.W == 0)
      {
         std::cout << "SUCCESS: " << B << " is a vector." << std::endl;
      }
      else
      {
         std::cout << "FAILURE: " << B << " is not a vector." << std::endl;
      }
   }

   {
      std::cout << "Epsilon:" << ww::EPSILON << std::endl;
   }

   {
      float const A{1.0};
      float const B = A + 0.5f * ww::EPSILON;
      if (ww::Equal(A, B))
      {
         std::cout.precision(27);
         std::cout << "SUCCESS: "      //!<
                   << std::fixed << A  //!<
                   << " is equal to "  //<!
                   << std::fixed << B  //<!
                   << std::endl;
      }
      else
      {
         std::cerr.precision(27);
         std::cerr << "FAILURE: "          //<!
                   << std::fixed << A      //!<
                   << " is not equal to "  //<!
                   << std::fixed << B      //<!
                   << std::endl;
      }
   }
}

// ---
// NOTE: Main function.
// ---
auto main(int argc, char *argv[]) -> int
{
   bool DoTheTest{};

   if (argc > 1)
   {
      std::string const Argv1{argv[1]};
      if ("-t" == Argv1)
      {
         std::cout << "Do the test" << std::endl;
         DoTheTest = true;
      }
   }

   if (DoTheTest)
   {
      RunTupleTest();
   }

   return 0;
}
