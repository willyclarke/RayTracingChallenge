/******************************************************************************
 * Filename : featuretest.hpp
 * Date     : 2018 Aug 29
 * Author   : Willy Clarke (willy@clarke.no)
 * Version  : 0.0.1
 * Copyright: W. Clarke
 * License  : MIT
 * Descripti: Test routines for the Raytracer project.
 ******************************************************************************/
#ifndef SRC_RAYTRACE_SRC_MAIN_FEATURETEST_HPP
#define SRC_RAYTRACE_SRC_MAIN_FEATURETEST_HPP

#include <datastructures.hpp>

#include <iostream>

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
         std::cout << "\nEquality\nSUCCESS: "  //!<
                   << std::fixed << A          //!<
                   << " is equal to "          //<!
                   << std::fixed << B          //<!
                   << std::endl;
      }
      else
      {
         std::cerr.precision(27);
         std::cerr << "\nEquality\nFAILURE: "  //<!
                   << std::fixed << A          //!<
                   << " is not equal to "      //<!
                   << std::fixed << B          //<!
                   << std::endl;
      }
   }

   // ---
   // NOTE: Test for equality
   // ---
   {
      ww::tup const A = ww::Vector(1.f, 2.f, 3.f);
      ww::tup const B = ww::Vector(1.f, 2.f, 3.f);
      if (ww::Equal(A, B))
      {
         std::cout << "\nEquality\nSUCCESS: \n"  //!<
                   << std::fixed << A            //<!
                   << "\n is equal to\n"         //<!
                   << std::fixed << B            //<!
                   << std::endl;
      }
      else
      {
         std::cout << "\nEquality\nFAILURE: \n"  //!<
                   << std::fixed << A            //<!
                   << "\n is not equal to\n"     //<!
                   << std::fixed << B            //<!
                   << std::endl;
      }
   }
   {
      ww::tup const A = ww::Vector(1.f, 2.f, 3.f);
      ww::tup const B = ww::Vector(1.f, 2.000000001f, 3.f);
      if (ww::Equal(A, B))
      {
         std::cout << "\nEquality\nSUCCESS: \n"  //!<
                   << std::fixed << A            //<!
                   << "\n is equal to\n"         //<!
                   << std::fixed << B            //<!
                   << std::endl;
      }
      else
      {
         std::cout << "\nEquality\nFAILURE: \n"  //!<
                   << std::fixed << A            //<!
                   << "\n is not equal to\n"     //<!
                   << std::fixed << B            //<!
                   << std::endl;
      }
   }

   // ---
   // NOTE: Check addition
   // ---
   {
      ww::tup const A = {3.f, -2.f, 5.f, 1.f};
      ww::tup const B = {-2.f, 3.f, 1.f, 0.f};
      ww::tup const Expect = {1.f, 1.f, 6.f, 1.f};
      ww::tup const C = ww::Add(A, B);
      if (ww::Equal(C, Expect))
      {
         std::cout << "\nAddition\nSUCCESS: \n"           //<!
                   << std::fixed << Expect                //<!
                   << "\n is equal to expected result\n"  //<!
                   << std::fixed << C                     //<!
                   << std::endl;
      }
      else
      {
         std::cerr << "\nAddition\nFAILURE: \n"               //<!
                   << std::fixed << Expect                    //<!
                   << "\n is not equal to expected result\n"  //<!
                   << std::fixed << C                         //<!
                   << std::endl;
      }
   }

   // ---
   // NOTE: Check subtraction
   // ---
   {
      ww::tup const P1 = ww::Point(3.f, 2.f, 1.f);
      ww::tup const P2 = ww::Point(5.f, 6.f, 7.f);
      ww::tup const Expect = {-2.f, -4.f, -6.f, 0.f};
      ww::tup const V = ww::Sub(P1, P2);
      if (ww::Equal(V, Expect) && IsVector(V))
      {
         std::cout << "\nSubtration\nSUCCESS: \n"         //<!
                   << std::fixed << Expect                //<!
                   << "\n is equal to expected result\n"  //<!
                   << std::fixed << V                     //<!
                   << std::endl;
      }
      else
      {
         std::cerr << "\nSubtration\nFAILURE: \n"             //<!
                   << std::fixed << Expect                    //<!
                   << "\n is not equal to expected result\n"  //<!
                   << std::fixed << V                         //<!
                   << std::endl;
      }
   }

   // ---
   // NOTE: Check subtraction
   // ---
   {
      ww::tup const P1 = ww::Point(3.f, 2.f, 1.f);
      ww::tup const V1 = ww::Vector(5.f, 6.f, 7.f);
      ww::tup const Expect = {-2.f, -4.f, -6.f, 1.f};
      ww::tup const P = ww::Sub(P1, V1);
      if (ww::Equal(P, Expect) && IsPoint(P))
      {
         std::cout << "\nSubtration Point - Vector\nSUCCESS: \n"  //<!
                   << std::fixed << P1 << "\n"                    //<!
                   << std::fixed << V1 << "\n"                    //<!
                   << "Result is\n"                               //<!
                   << std::fixed << Expect                        //<!
                   << "\n is equal to expected result\n"          //<!
                   << std::fixed << P                             //<!
                   << std::endl;
      }
      else
      {
         std::cerr << "\nSubtration Point - Vector\nFAILURE: \n"  //<!
                   << std::fixed << P1 << "\n"                    //<!
                   << std::fixed << V1 << "\n"                    //<!
                   << "Result is\n"                               //<!
                   << std::fixed << Expect                        //<!
                   << "\n is not equal to expected result\n"      //<!
                   << std::fixed << P                             //<!
                   << std::endl;
      }
   }
}
#endif
