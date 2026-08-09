#pragma once
// SPH smoothing kernels, 3D, parameterised by SUPPORT radius H (not the
// smoothing length h = H/2). Both kernels vanish for r > H.
//
// Cubic spline:
//   W(r,H)      = 8/(pi H^3) * { 1 - 6q^2 + 6q^3      , 0 <= q <= 1/2
//                              { 2 (1-q)^3            , 1/2 < q <= 1
//   grad W      = 48/(pi H^4) * { q(3q-2) e           , 0 <  q <= 1/2
//                              { -(1-q)^2 e           , 1/2 < q <= 1
//
// Wendland C2 (the solver default: strictly non-negative, smooth first
// derivative, cheaper than the quintic spline):
//   W(r,H)      = 21/(2 pi H^3) (1-q)^4 (1+4q)
//   grad W      = -210/(pi H^4) q (1-q)^3 e
//
// with q = r/H and e = (x_i - x_j)/r. Both normalisations are the standard 3D
// ones and match the reference implementation in SPlisHSPlasH's SPHKernels.h.
//
// Cohesion kernel for surface tension (Akinci, Akinci & Teschner 2013,
// "Versatile Surface Tension and Adhesion for SPH Fluids", eq. 2):
//   C(r) = 32/(pi H^9) * { (H-r)^3 r^3                 , H/2 < r <= H
//                        { 2 (H-r)^3 r^3 - H^6/64      , 0 < r <= H/2
// Its coefficient is resolution dependent and MUST be calibrated against the
// Young-Laplace law rather than set equal to a measured interfacial tension;
// see solver.hpp.

#include <cmath>

namespace chemcad::fluid {

constexpr double kPi = 3.14159265358979323846;

inline double wendlandW(double r, double h) {
  if (r >= h || h <= 0.0) return 0.0;
  const double q = r / h;
  const double oneMinusQ = 1.0 - q;
  const double q4 = oneMinusQ * oneMinusQ * oneMinusQ * oneMinusQ;
  return 21.0 / (2.0 * kPi * h * h * h) * q4 * (1.0 + 4.0 * q);
}

// Magnitude of grad W along -e (i.e. dW/dr); multiply by e to get the vector.
inline double wendlandGradMagnitude(double r, double h) {
  if (r >= h || r <= 0.0 || h <= 0.0) return 0.0;
  const double q = r / h;
  const double oneMinusQ = 1.0 - q;
  return -210.0 / (kPi * h * h * h * h) * q * oneMinusQ * oneMinusQ * oneMinusQ;
}

inline double cubicW(double r, double h) {
  if (r >= h || h <= 0.0) return 0.0;
  const double q = r / h;
  const double k = 8.0 / (kPi * h * h * h);
  if (q <= 0.5) return k * (1.0 - 6.0 * q * q + 6.0 * q * q * q);
  const double oneMinusQ = 1.0 - q;
  return k * 2.0 * oneMinusQ * oneMinusQ * oneMinusQ;
}

inline double cubicGradMagnitude(double r, double h) {
  if (r >= h || r <= 0.0 || h <= 0.0) return 0.0;
  const double q = r / h;
  const double k = 48.0 / (kPi * h * h * h * h);
  if (q <= 0.5) return k * q * (3.0 * q - 2.0);
  const double oneMinusQ = 1.0 - q;
  return -k * oneMinusQ * oneMinusQ;
}

// Akinci 2013 cohesion spline.
inline double cohesionC(double r, double h) {
  if (r <= 0.0 || r >= h || h <= 0.0) return 0.0;
  const double k = 32.0 / (kPi * std::pow(h, 9.0));
  const double hr = h - r;
  const double hr3r3 = hr * hr * hr * r * r * r;
  if (r > 0.5 * h) return k * hr3r3;
  const double h6 = std::pow(h, 6.0);
  return k * (2.0 * hr3r3 - h6 / 64.0);
}

// Number density of a perfect infinite rest lattice of spacing dx under this
// kernel: delta0 = sum_j W(|x_i - x_j|, H) over the lattice inside H. The
// solver's density error is (delta_i - delta0)/delta0, so delta0 MUST be
// computed from the same kernel and spacing rather than assumed to be 1/dx^3.
double restNumberDensity(double spacing, double support, bool wendland = true);

}  // namespace chemcad::fluid
