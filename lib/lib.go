package lib

import (
	"math/big"
)

const (
	PALLAS_P = "28948022309329048855892746252171976963363056481941560715954676764349967630337"
	PALLAS_R = "28948022309329048855892746252171976963363056481941647379679742748393362948097"
)

type PallasCurve struct {
	P *big.Int
	R *big.Int
}

type VestaCurve struct {
	P *big.Int
	R *big.Int
}

type AffinePoint struct {
	X *big.Int
	Y *big.Int
}

type ProjectivePoint struct {
	X *big.Int
	Y *big.Int
	Z *big.Int
}

// EllipticCurve represents a generic elliptic curve.
type PastaCurve interface {
	AffineGenerator() AffinePoint
	ProjectiveGenerator() ProjectivePoint
	GetP() *big.Int
	GetR() *big.Int
	IsOnCurve(p AffinePoint) bool
	ValidateScalarField(s *big.Int) bool
}

// InitPallas initializes the Pallas curve.
func (c *PallasCurve) InitPallas() {
	c.P = new(big.Int)
	c.R = new(big.Int)
	c.P.SetString(PALLAS_P, 10)
	c.R.SetString(PALLAS_R, 10)
}

// Returns the generator of the Pallas curve in affine coordinates.
func (c PallasCurve) AffineGenerator() AffinePoint {
	generatorX, _ := new(big.Int).SetString("28948022309329048855892746252171976963363056481941560715954676764349967630336", 10)
	return AffinePoint{
		X: generatorX,
		Y: big.NewInt(2),
	}
}

// Returns the generator of the Pallas curve in projective coordinates.
func (c PallasCurve) ProjectiveGenerator() ProjectivePoint {
	generatorX, _ := new(big.Int).SetString("28948022309329048855892746252171976963363056481941560715954676764349967630336", 10)
	return ProjectivePoint{
		X: generatorX,
		Y: big.NewInt(2),
		Z: big.NewInt(1),
	}
}

// GetP returns the prime modulus of the Pallas curve.
func (c PallasCurve) GetP() *big.Int {
	return c.P
}

// GetR returns the order of the Pallas curve.
func (c PallasCurve) GetR() *big.Int {
	return c.R
}

// IsOnCurve checks if the given point is on the Pallas curve.
// The Pallas curve is defined by the equation y^2 = x^3 + 5 mod p, where p is a prime modulus.
// Function checks that:
//  1. x != 0
//  2. y != 0
//  3. x < p
//  4. y < p
//  5. y^2 = x^3 + 5 mod p
func (c PallasCurve) IsOnCurve(p AffinePoint) bool {
	ySquared := new(big.Int).Mul(p.Y, p.Y)
	ySquared.Mod(ySquared, c.P)

	xCubed := new(big.Int).Exp(p.X, big.NewInt(3), c.P)
	xCubedPlusFive := new(big.Int).Add(xCubed, big.NewInt(5))
	xCubedPlusFive.Mod(xCubedPlusFive, c.P)

	isWellFormed := p.X.Cmp(big.NewInt(0)) != 0 &&
		p.Y.Cmp(big.NewInt(0)) != 0 &&
		p.X.Cmp(c.P) < 0 &&
		p.Y.Cmp(c.P) < 0 &&
		ySquared.Cmp(xCubedPlusFive) == 0

	return isWellFormed
}

// Vesta!!!!!

// Returns prime modulus of the Vesta curve.
func (v VestaCurve) GetP() *big.Int {
	return v.P
}

// Returns order of the Vesta curve.
func (v VestaCurve) GetR() *big.Int {
	return v.R
}

// AffineAddition returns the sum of the two given point on the curve in affine form.
func AffineAddition(c PastaCurve, p1, p2 AffinePoint) (AffinePoint, bool) {
	// If p1 is the infinity point, return p2.
	if p1.IsInfinity() {
		return p2, true
	}

	// If p2 is the infinity point, return p1.
	if p2.IsInfinity() {
		return p1, true
	}

	// If both points are the same, return double(p1).
	if p1.X.Cmp(p2.X) == 0 && p1.Y.Cmp(p2.Y) == 0 {
		return DoubleAffine(c, p1)
	}

	// Get the prime modulus of the curve.
	pMod := c.GetP()

	// Calculate lambda
	// lambda = (y2 - y1) / (x2 - x1)
	lambda := new(big.Int).Add(p1.X, pMod)
	lambda.Sub(lambda, p2.X)

	t := new(big.Int).Add(p1.Y, pMod)
	t.Sub(t, p2.Y)

	if lambda.Cmp(pMod) == 0 {
		lambda.Sub(lambda, pMod)
	}

	lambda.ModInverse(lambda, pMod)
	lambda.Mul(lambda, t)
	lambda.Mod(lambda, pMod)

	x3 := new(big.Int).Mul(lambda, lambda)
	x3.Add(x3, pMod)
	x3.Add(x3, pMod)
	x3.Sub(x3, p1.X)
	x3.Sub(x3, p2.X)
	x3.Mod(x3, pMod)

	y3 := new(big.Int).Add(p1.X, pMod)
	y3.Sub(y3, x3)
	y3.Mul(lambda, y3)
	y3.Mod(y3, pMod)
	y3.Add(y3, pMod)
	y3.Sub(y3, p1.Y)
	y3.Mod(y3, pMod)

	ret := AffinePoint{
		X: x3,
		Y: y3,
	}

	// Check if the result is on the curve.
	if onCurve := c.IsOnCurve(ret); onCurve {
		return ret, true
	}

	return AffinePoint{}, false
}

// DoubleAffine returns the double of the given point on the  curve in affine form.
func DoubleAffine(curve PastaCurve, p AffinePoint) (AffinePoint, bool) {
	// If point is infinity, return infinity.
	if p.IsInfinity() {
		return p, false
	}

	pMod := curve.GetP()

	// Calculate the value (3 / 2) mod p.
	threeOverTwo, _ := new(big.Int).SetString("14474011154664524427946373126085988481681528240970780357977338382174983815170", 10)

	// Calculate yInv = y^-1 mod p.
	yInv := new(big.Int).ModInverse(p.Y, pMod)

	// Calculate lambda = (3 * x^2 * yInv / 2) mod p.
	lambda := new(big.Int).Mul(p.X, p.X)
	lambda.Mod(lambda, pMod)
	lambda.Mul(lambda, yInv)
	lambda.Mul(lambda, threeOverTwo)
	lambda.Mod(lambda, pMod)

	// Calculate xPrime = (lambda^2 - 2 * x) mod p.
	xPrime := new(big.Int).Mul(lambda, lambda)
	xPrime.Mod(xPrime, pMod)
	xPrime.Add(xPrime, pMod)
	xPrime.Add(xPrime, pMod)
	xPrime.Sub(xPrime, p.X)
	xPrime.Sub(xPrime, p.X)
	xPrime.Mod(xPrime, pMod)

	// Calculate yPrime = (lambda * (x - xPrime) - y) mod p.
	yPrime := new(big.Int).Add(p.X, pMod)
	yPrime.Sub(yPrime, xPrime)
	yPrime.Mul(lambda, yPrime)
	yPrime.Mod(yPrime, pMod)
	yPrime.Add(yPrime, pMod)
	yPrime.Sub(yPrime, p.Y)
	yPrime.Mod(yPrime, pMod)

	ret := AffinePoint{
		X: xPrime,
		Y: yPrime,
	}

	// Check if the result is on the curve.
	if onCurve := curve.IsOnCurve(ret); onCurve {
		return ret, true
	}

	return AffinePoint{}, false
}

// fromLeBytesModOrder converts a little-endian byte slice to a big integer modulo curve R_MOD.
// The function ensures that the resulting integer is within the valid range [0, R_MOD-1].
// ! Should return 2 values, the second one being a bool indicating if the conversion was successful.
func FromLeBytesModOrder(curve PastaCurve, leBytes []byte) *big.Int {
	res := big.NewInt(0)
	rMod := curve.GetR()

	for i := 0; i < len(leBytes); i++ {
		res.Mul(res, big.NewInt(256))
		res.Mod(res, rMod)
		res.Add(res, big.NewInt(int64(leBytes[len(leBytes)-1-i])))
		res.Mod(res, rMod)
	}

	if !curve.ValidateScalarField(res) {
		return nil
	}

	return res
}

// ValidateScalarField checks if the given scalar is in the field of the Pallas curve.
// The function checks if the scalar value is not in the valid range [0, R_MOD-1].
func (curve PallasCurve) ValidateScalarField(s *big.Int) bool {
	rMod := curve.GetR()
	return s.Cmp(rMod) < 0
}

// IsInfinity checks if the given point represents the infinity point on the Pallas curve.
func (p AffinePoint) IsInfinity() bool {
	return p.X.Sign() == 0 && p.Y.Sign() == 0
}

// IsInfinity checks if the given point represents the infinity point on the Pallas curve.
func (p ProjectivePoint) IsInfinity() bool {
	return p.X.Sign() == 0 && p.Y.Sign() == 0 && p.Z.Sign() == 0
}
