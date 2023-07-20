// SPDX-License-Identifier: MIT
pragma solidity >=0.8.0;

/// @title Vesta Elliptic Curve Interface
/// @notice This interface defines the functions for interacting with the Vesta curve
///         using affine and projective coordinates.
interface IVesta {
  struct AffinePoint {
    uint256 x;
    uint256 y;
  }

  struct ProjectivePoint {
    uint256 x;
    uint256 y;
    uint256 z;
  }

  /// @dev Returns generator point of the Vesta curve in affine coordinates.
  function affineGenerator() external pure returns (AffinePoint memory g);

  /// @dev Returns generator point of the Vesta curve in projective coordinates.
  function projectiveGenerator() external pure returns (ProjectivePoint memory g);

  /// @dev Converts a point from projective coordinates to affine coordinates.
  /// @param p The projective point to convert.
  /// @return r The resulting affine point after conversion.
  /// @return success A boolean indicating if the conversion was successful.
  function toAffine(ProjectivePoint memory p) external pure returns (AffinePoint memory r, bool success);

  /// @dev Converts a point from affine coordinates to projective coordinates.
  /// @param p The affine point to convert.
  /// @return r The resulting projective point after conversion.
  /// @return success A boolean indicating if the conversion was successful.
  function toProjective(AffinePoint memory p) external pure returns (ProjectivePoint memory r, bool success);

  /// @dev Checks if a point is in affine coordinates.
  /// @param p The affine point to check.
  /// @return r A boolean indicating if the point is at infinity (false) or not (true).
  function isInfinityAffine(AffinePoint memory p) external pure returns (bool r);

  /// @dev Checks if a point is in projective coordinates.
  /// @param p The projective point to check.
  /// @return r A boolean indicating if the point is at infinity (false) or not (true).
  function isInfinityProjective(ProjectivePoint memory p) external pure returns (bool r);

  /// @dev Negates an affine point.
  /// @param p The affine point to negate.
  /// @return r The negated affine point.
  /// @return success A boolean indicating if the negation was successful.
  function negateAffine(AffinePoint memory p) external pure returns (AffinePoint memory r, bool success);

  /// @dev Negates a projective point.
  /// @param p The projective point to negate.
  /// @return r The negated projective point.
  /// @return success A boolean indicating if the negation was successful.
  function negateProjective(ProjectivePoint memory p) external pure returns (ProjectivePoint memory r, bool success);

  /// @dev Returns the negation of a scalar on the Vesta Curve.
  /// @param scalar The scalar to negate.
  /// @return nScalar The negated scalar value.
  function negateScalar(uint256 scalar) external pure returns (uint256 nScalar);

  /// @dev Doubles an affine point on the Vesta Curve.
  /// @param p The affine point to double.
  /// @return r The resulting affine point after doubling.
  /// @return success A boolean indicating if the doubling was successful.
  function doubleAffine(AffinePoint memory p) external pure returns (AffinePoint memory r, bool success);

  /// @dev Doubles a projective point on the Vesta Curve.
  /// @param p The projective point to double.
  /// @return r The resulting projective point after doubling.
  /// @return success A boolean indicating if the doubling was successful.
  function doubleProjective(ProjectivePoint memory p) external pure returns (ProjectivePoint memory r, bool success);

  /// @dev Adds two affine points on the Vesta Curve.
  /// @param p The first affine point.
  /// @param q The second affine point.
  /// @return r The resulting affine point after addition.
  /// @return success A boolean indicating if the addition was successful.
  function affineAdd(
    AffinePoint memory p,
    AffinePoint memory q
  ) external pure returns (AffinePoint memory r, bool success);

  /// @dev Adds two projective points on the Vesta Curve.
  /// @param p The first projective point.
  /// @param q The second projective point.
  /// @return r The resulting projective point after addition.
  /// @return success A boolean indicating if the addition was successful.
  function projectiveAdd(
    ProjectivePoint memory p,
    ProjectivePoint memory q
  ) external pure returns (ProjectivePoint memory r, bool success);

  /// @dev Multiplies an affine point by a scalar on the Vesta Curve.
  /// @param p The affine point to multiply.
  /// @param scalar The scalar to multiply the point by.
  /// @return r The resulting affine point after multiplication.
  /// @return success A boolean indicating if the multiplication was successful.
  function mulWithScalarAffine(
    AffinePoint memory p,
    uint256 scalar
  ) external pure returns (AffinePoint memory r, bool success);

  /// @dev Multiplies a projective point by a scalar on the Vesta Curve.
  /// @param p The projective point to multiply.
  /// @param scalar The scalar to multiply the point by.
  /// @return r The resulting projective point after multiplication.
  /// @return success A boolean indicating if the multiplication was successful.
  function mulWithScalarProjective(
    ProjectivePoint memory p,
    uint256 scalar
  ) external pure returns (ProjectivePoint memory r, bool success);

  /// @dev Vesta Inverse Mod
  /// @param fr The value to invert modulo m.
  /// @param m The modulus to use for the inversion.
  /// @return iFr The inverse of `fr` modulo `m`.
  function inverseMod(uint256 fr, uint256 m) external pure returns (uint256 iFr);

  /// @dev Checks if a point is on the Vesta Curve.
  /// @param p The affine point to check.
  /// @return r A boolean indicating if the point is on the curve (true) or not (false).
  function isOnCurve(AffinePoint memory p) external pure returns (bool r);

  /// @dev Validates if a scalar is in the valid range [0, R_MOD-1].
  /// @param scalar The scalar to validate.
  /// @return r A boolean indicating if the scalar is in the valid range.
  function validateScalarField(uint256 scalar) external pure returns (bool r);

  /// @dev Converts little-endian bytes to a scalar modulo R_MOD.
  /// @param leBytes The little-endian bytes representing the scalar.
  /// @return r The resulting scalar value.
  function fromLeBytesModOrder(bytes memory leBytes) external pure returns (uint256 r);

  /// @dev Computes `base` raised to the power `exponent` modulo `modulus`.
  /// @param base The base value.
  /// @param exponent The exponent value.
  /// @param modulus The modulus value.
  /// @return r The resulting value of `base` raised to the power `exponent` modulo `modulus`.
  function powSmall(uint256 base, uint256 exponent, uint256 modulus) external pure returns (uint256 r);
}
