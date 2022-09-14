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
 */
float RayMarch(ww::ray const &R);

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
