#include "fluid/kernels.hpp"

#include <cmath>
#include <limits>

namespace chemcad::fluid {

double restNumberDensity(double spacing, double support, bool wendland) {
  if (!(spacing > 0.0) || !(support > 0.0) || !std::isfinite(spacing) ||
      !std::isfinite(support)) {
    return 0.0;
  }

  // A fixed integer box, rather than a radius-dependent traversal, makes the
  // summation order part of the numerical contract. The ceil includes the
  // compact-support boundary, where both kernels evaluate to exactly zero.
  const double extentValue = std::ceil(support / spacing);
  if (!std::isfinite(extentValue) ||
      extentValue >= static_cast<double>(std::numeric_limits<int>::max())) {
    return 0.0;
  }
  const int extent = static_cast<int>(extentValue);
  double sum = 0.0;
  for (int iz = -extent; iz <= extent; ++iz) {
    for (int iy = -extent; iy <= extent; ++iy) {
      for (int ix = -extent; ix <= extent; ++ix) {
        const double x = static_cast<double>(ix) * spacing;
        const double y = static_cast<double>(iy) * spacing;
        const double z = static_cast<double>(iz) * spacing;
        const double r = std::sqrt(x * x + y * y + z * z);
        if (r <= support) {
          sum += wendland ? wendlandW(r, support) : cubicW(r, support);
        }
      }
    }
  }
  return sum;
}

}  // namespace chemcad::fluid
