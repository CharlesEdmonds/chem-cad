#include "fluid/frame.hpp"

#include <cmath>

#include "fluid/kernels.hpp"

namespace chemcad::fluid {
namespace {

std::array<double, 3> rotateWorldToBody(const std::array<double, 4>& qInput,
                                        const std::array<double, 3>& v) {
  double w = qInput[0];
  double x = qInput[1];
  double y = qInput[2];
  double z = qInput[3];
  const double norm = std::sqrt(w * w + x * x + y * y + z * z);
  if (norm > 1.0e-15) {
    w /= norm;
    x /= norm;
    y /= norm;
    z /= norm;
  } else {
    w = 1.0;
    x = y = z = 0.0;
  }

  // A pose quaternion maps body coordinates to world coordinates. Applying
  // R^T therefore uses the transpose of the usual quaternion rotation matrix.
  return {
      (1.0 - 2.0 * (y * y + z * z)) * v[0] + 2.0 * (x * y + w * z) * v[1] +
          2.0 * (x * z - w * y) * v[2],
      2.0 * (x * y - w * z) * v[0] + (1.0 - 2.0 * (x * x + z * z)) * v[1] +
          2.0 * (y * z + w * x) * v[2],
      2.0 * (x * z + w * y) * v[0] + 2.0 * (y * z - w * x) * v[1] +
          (1.0 - 2.0 * (x * x + y * y)) * v[2]};
}

}  // namespace

std::array<double, 3> shakeAcceleration(const VesselMotion& motion, double timeS) {
  if (!motion.shaking || motion.shakeRemainingS <= 0.0 || motion.shakeAmplitudeM == 0.0 ||
      motion.shakeFrequencyHz == 0.0) {
    return {0.0, 0.0, 0.0};
  }

  const double axisNorm = std::sqrt(motion.shakeAxis[0] * motion.shakeAxis[0] +
                                    motion.shakeAxis[1] * motion.shakeAxis[1] +
                                    motion.shakeAxis[2] * motion.shakeAxis[2]);
  if (axisNorm <= 1.0e-15) return {0.0, 0.0, 0.0};

  const double omega = 2.0 * kPi * motion.shakeFrequencyHz;
  const double magnitude =
      -motion.shakeAmplitudeM * omega * omega * std::sin(omega * timeS) / axisNorm;
  return {magnitude * motion.shakeAxis[0], magnitude * motion.shakeAxis[1],
          magnitude * motion.shakeAxis[2]};
}

FrameAcceleration frameAcceleration(const VesselMotion& motion, double timeS) {
  constexpr std::array<double, 3> kWorldGravity{0.0, 0.0, -9.80665};
  const std::array<double, 3> gravityBody =
      rotateWorldToBody(motion.pose.orientation, kWorldGravity);
  const std::array<double, 3> shakeWorld = shakeAcceleration(motion, timeS);
  const std::array<double, 3> translationalWorld{
      shakeWorld[0] + motion.manualAcceleration[0],
      shakeWorld[1] + motion.manualAcceleration[1],
      shakeWorld[2] + motion.manualAcceleration[2]};
  const std::array<double, 3> translationalBody =
      rotateWorldToBody(motion.pose.orientation, translationalWorld);

  FrameAcceleration result;
  result.uniform = {gravityBody[0] - translationalBody[0],
                    gravityBody[1] - translationalBody[1],
                    gravityBody[2] - translationalBody[2]};
  result.omega = motion.angularVelocity;
  result.alpha = motion.angularAcceleration;
  constexpr double kRotationEpsilon = 1.0e-12;
  const double omegaSquared = result.omega[0] * result.omega[0] +
                              result.omega[1] * result.omega[1] +
                              result.omega[2] * result.omega[2];
  const double alphaSquared = result.alpha[0] * result.alpha[0] +
                              result.alpha[1] * result.alpha[1] +
                              result.alpha[2] * result.alpha[2];
  result.rotating = omegaSquared > kRotationEpsilon * kRotationEpsilon ||
                    alphaSquared > kRotationEpsilon * kRotationEpsilon;
  return result;
}

double shakePeakVelocity(const VesselMotion& motion) {
  return 2.0 * kPi * std::abs(motion.shakeFrequencyHz) * std::abs(motion.shakeAmplitudeM);
}

double shakeSpecificPower(const VesselMotion& motion) {
  const double velocity = shakePeakVelocity(motion);
  return 0.5 * velocity * velocity * std::abs(motion.shakeFrequencyHz);
}

}  // namespace chemcad::fluid
