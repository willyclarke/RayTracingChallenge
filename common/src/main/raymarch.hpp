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
 * RayMarch - move along ray to find hit.
 * @param: R - Ray
 * @param: PtrShape - Shape to compute distance against
 */
float RayMarch(ww::ray const &R, shared_ptr_shape PtrShape = nullptr);

/**
 */
tup MainImage(int X, int Y, int W, int H, shared_ptr_shape PtrShape = nullptr);

/**
 */
canvas Render(camera const &Camera, world const &World);

/**
 */
float GetDistance(tup const &P, shared_ptr_shape PtrShape = nullptr);

/**
 */
tup GetNormal(tup const &P, shared_ptr_shape PtrShape = nullptr);

};      // end namespace rm
};      // end namespace ww
#endif  // COMMON_SRC_MAIN_RAYMARCH_HPP
