package lib

import (
	"math/big"
)

const (
	PALLAS_P = "28948022309329048855892746252171976963363056481941560715954676764349967630337"
	PALLAS_R = "28948022309329048855892746252171976963363056481941647379679742748393362948097"

	VESTA_P = "28948022309329048855892746252171976963363056481941647379679742748393362948097"
	VESTA_R = "28948022309329048855892746252171976963363056481941560715954676764349967630337"
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

// Vesta!!!!!

// Returns prime modulus of the Vesta curve.
func (v VestaCurve) GetP() *big.Int {
	return v.P
}

// Returns order of the Vesta curve.
func (v VestaCurve) GetR() *big.Int {
	return v.R
}

// Checks if the given scalar is in the field of the curve.
func (v VestaCurve) ValidateScalarField(s *big.Int) bool {
	rMod := v.GetR()
	return s.Cmp(rMod) < 0
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
	if onCurve := IsOnCurve(c, ret); onCurve {
		return ret, true
	}

	return AffinePoint{}, false
}

// IsOnCurve checks if the given point is on the Pallas curve.
// The Pallas curve is defined by the equation y^2 = x^3 + 5 mod p, where p is a prime modulus.
// Function checks that:
//  1. x != 0
//  2. y != 0
//  3. x < p
//  4. y < p
//  5. y^2 = x^3 + 5 mod p
func IsOnCurve(c PastaCurve, p AffinePoint) bool {
	pMod := c.GetP()
	ySquared := new(big.Int).Mul(p.Y, p.Y)
	ySquared.Mod(ySquared, pMod)

	xCubed := new(big.Int).Exp(p.X, big.NewInt(3), pMod)
	xCubedPlusFive := new(big.Int).Add(xCubed, big.NewInt(5))
	xCubedPlusFive.Mod(xCubedPlusFive, pMod)

	isWellFormed := p.X.Cmp(big.NewInt(0)) != 0 &&
		p.Y.Cmp(big.NewInt(0)) != 0 &&
		p.X.Cmp(pMod) < 0 &&
		p.Y.Cmp(pMod) < 0 &&
		ySquared.Cmp(xCubedPlusFive) == 0

	return isWellFormed
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
	if onCurve := IsOnCurve(curve, ret); onCurve {
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

// MulWithScalar returns the product of the given point and scalar on the Pallas curve.
// The function uses the "double-and-add" algorithm to efficiently calculate the result.
func MulWithScalarAffine(curve PastaCurve, p AffinePoint, s *big.Int) (AffinePoint, bool) {
	tmp := p
	var ok bool

	r := AffinePoint{X: big.NewInt(0), Y: big.NewInt(0)}

	for i := 0; i < 256; i++ {
		bit := new(big.Int).And(s, big.NewInt(1))

		s.Div(s, big.NewInt(2))

		if bit.Cmp(big.NewInt(1)) == 0 {
			r, ok = AffineAddition(curve, r, tmp)
			if !ok {
				return AffinePoint{}, false
			}
		}

		// Double the current temporary point (tmp) for the next iteration.
		tmp, ok = DoubleAffine(curve, tmp)
		if !ok {
			return AffinePoint{}, false
		}
	}

	return r, true
}

func MulWithScalarProjective(curve PastaCurve, p ProjectivePoint, s *big.Int) (ProjectivePoint, bool) {
	var ok bool
	tmp := p
	r := ProjectivePoint{big.NewInt(0), big.NewInt(0), big.NewInt(0)}

	for i := 0; i < 256; i++ {
		bit := new(big.Int).And(s, big.NewInt(1))
		s.Div(s, big.NewInt(2))
		if bit.Cmp(big.NewInt(1)) == 0 {
			r, ok = ProjectiveAddition(curve, r, tmp)
			if !ok {
				return ProjectivePoint{}, false
			}
		}
		tmp, ok = DoubleProjective(curve, tmp)
		if !ok {
			return ProjectivePoint{}, false
		}
	}

	return r, true
}

// NegateAffine calculates the negation of a point on the curve.
// Negating a point (x, y) on the curve results in the point (x, -y).
func NegateAffine(curve PastaCurve, p AffinePoint) (AffinePoint, bool) {
	if p.IsInfinity() {
		return p, true
	}

	pMod := curve.GetP()

	// pMod - (p.Y % pMod) = -p.Y mod pMod
	negatedY := new(big.Int).Sub(pMod, new(big.Int).Mod(p.Y, pMod))

	return AffinePoint{X: p.X, Y: negatedY}, true
}

// NegateProjective calculates the negation of a point on the curve in projective coordinates.
// Negating a point (x, y, z) on the curve results in the point (x, -y, z).
func NegateProjective(curve PastaCurve, p ProjectivePoint) (ProjectivePoint, bool) {
	if p.IsInfinity() {
		return p, true
	}

	pMod := curve.GetP()

	return ProjectivePoint{
		X: p.X,
		Y: new(big.Int).Sub(pMod, new(big.Int).Mod(p.Y, pMod)),
		Z: p.Z,
	}, true
}

// NegateScalar returns the negation of the given scalar on the Pallas curve.
func NegateScalar(curve PastaCurve, s *big.Int) *big.Int {
	rMod := curve.GetR()

	// The negation of fr is represented as (R_MOD - fr) mod R_MOD.
	negatedScalar := new(big.Int).Sub(rMod, new(big.Int).Mod(s, rMod))
	return negatedScalar
}

// ToAffine converts the given point from projective to affine coordinates.
func ToAffine(curve PastaCurve, p ProjectivePoint) (AffinePoint, bool) {
	if p.IsInfinity() {
		return AffinePoint{big.NewInt(0), big.NewInt(0)}, false
	}

	pMod := curve.GetP()

	x := new(big.Int).ModInverse(p.Z, pMod)
	y := new(big.Int).Mul(p.Y, x)
	y.Mod(y, pMod)
	x.Mul(p.X, x)
	x.Mod(x, pMod)

	return AffinePoint{x, y}, true
}

// ToProjective converts the given point from affine to projective coordinates.
func ToProjective(p AffinePoint) (ProjectivePoint, bool) {
	if p.IsInfinity() {
		return ProjectivePoint{big.NewInt(0), big.NewInt(0), big.NewInt(0)}, false
	}

	return ProjectivePoint{p.X, p.Y, big.NewInt(1)}, true
}

// DoubleProjective returns the double of the given point on the curve in projective form.
func DoubleProjective(curve PastaCurve, point ProjectivePoint) (ProjectivePoint, bool) {
	// If the point is the infinity point, return the point as it is.
	if point.IsInfinity() {
		return point, false
	}

	// Get the prime modulus of the curve.
	pMod := curve.GetP()

	// Calculate 2P, where P is the given point.
	doubleP := new(big.Int).Lsh(pMod, 1)

	// Calculate intermediate values for the doubling formula.
	a := new(big.Int).Mul(point.X, point.X)
	a.Mod(a, pMod)

	b := new(big.Int).Mul(point.Y, point.Y)
	b.Mod(b, pMod)

	c := new(big.Int).Mul(b, b)
	c.Mod(c, pMod)

	d := new(big.Int).Add(point.X, b)
	d.Mul(d, d)
	d.Mod(d, pMod)
	b.Add(a, c)
	d.Add(d, new(big.Int).Sub(doubleP, b))
	d.Mul(d, big.NewInt(2))
	d.Mod(d, pMod)

	e := new(big.Int).Mul(a, big.NewInt(3))

	f := new(big.Int).Mul(e, e)
	f.Mod(f, pMod)

	z := new(big.Int).Mul(point.Y, big.NewInt(2))
	z.Mul(z, point.Z)
	z.Mod(z, pMod)

	// Calculate the x-coordinate of the doubled point.
	x := new(big.Int).Add(f, new(big.Int).Sub(doubleP, new(big.Int).Mul(d, big.NewInt(2))))
	x.Mod(x, pMod)

	// Calculate the y-coordinate of the doubled point.
	y := new(big.Int).Add(d, new(big.Int).Sub(pMod, x))
	y.Mul(e, y)
	y.Mod(y, pMod)
	y.Add(y, new(big.Int).Sub(pMod, new(big.Int).Mul(c, big.NewInt(8))))
	y.Mod(y, pMod)

	// Return the doubled point in projective coordinates.
	return ProjectivePoint{x, y, z}, true
}

// ProjectiveAddition adds two points in projective coordinates on the curve.
func ProjectiveAddition(curve PastaCurve, p1 ProjectivePoint, p2 ProjectivePoint) (ProjectivePoint, bool) {
	// If p1 is the infinity point, return p2.
	if p1.IsInfinity() {
		return p2, true
	}

	// If p2 is the infinity point, return p1.
	if p2.IsInfinity() {
		return p1, true
	}

	// Get the prime modulus of the curve.
	pMod := curve.GetP()

	// Calculate intermediate values for the addition formula.
	z1z1 := new(big.Int).Mul(p1.Z, p1.Z)
	z1z1.Mod(z1z1, pMod)

	z2z2 := new(big.Int).Mul(p2.Z, p2.Z)
	z2z2.Mod(z2z2, pMod)

	u1 := new(big.Int).Mul(p1.X, z2z2)
	u1.Mod(u1, pMod)

	u2 := new(big.Int).Mul(p2.X, z1z1)
	u2.Mod(u2, pMod)

	s1 := new(big.Int).Mul(p1.Y, p2.Z)
	s1.Mul(s1, z2z2)
	s1.Mod(s1, pMod)

	s2 := new(big.Int).Mul(p2.Y, p1.Z)
	s2.Mul(s2, z1z1)
	s2.Mod(s2, pMod)

	// Check if the points are the same.
	if u1.Cmp(u2) == 0 {
		if s1.Cmp(s2) == 0 {
			// If both points are the same, return double of p1.
			return DoubleProjective(curve, p1)
		}
	}

	// Calculate intermediate values for the addition formula.
	h := new(big.Int).Add(u2, new(big.Int).Sub(pMod, u1))

	i := new(big.Int).Add(h, h)
	i.Mod(i, pMod)
	i.Mul(i, i)
	i.Mod(i, pMod)

	j := new(big.Int).Mul(h, i)
	j.Mod(j, pMod)

	r := new(big.Int).Add(s2, new(big.Int).Sub(pMod, s1))
	r.Add(r, r)
	r.Mod(r, pMod)

	v := new(big.Int).Mul(u1, i)
	v.Mod(v, pMod)

	// Calculate x-coordinate of the result point.
	x3 := new(big.Int).Mul(r, r)
	x3.Mod(x3, pMod)
	tripleP := new(big.Int).Mul(pMod, big.NewInt(3))
	x3.Add(x3, new(big.Int).Sub(tripleP, new(big.Int).Add(j, new(big.Int).Add(v, v))))
	x3.Mod(x3, pMod)

	// Calculate y-coordinate of the result point.
	y3 := new(big.Int).Add(v, new(big.Int).Sub(pMod, x3))
	y3.Mul(r, y3)
	y3.Mod(y3, pMod)
	s1.Mul(s1, big.NewInt(2))
	s1.Mod(s1, j)
	y3.Add(y3, new(big.Int).Sub(pMod, s1))
	y3.Mod(y3, pMod)

	// Calculate z-coordinate of the result point.
	z3 := new(big.Int).Add(p1.Z, p2.Z)
	z3.Mul(z3, z3)
	z3.Mod(z3, pMod)
	doubleP := new(big.Int).Mul(pMod, big.NewInt(2))
	z3.Add(z3, new(big.Int).Sub(doubleP, new(big.Int).Add(z1z1, z2z2)))
	z3.Mul(z3, h)
	z3.Mod(z3, pMod)

	return ProjectivePoint{x3, y3, z3}, true
}

// ValidateScalarField checks if the given scalar is in the field of the Pallas curve.
// The function checks if the scalar value is not in the valid range [0, R_MOD-1].
func (curve PallasCurve) ValidateScalarField(s *big.Int) bool {
	rMod := curve.GetR()
	return s.Cmp(rMod) < 0
}

// IsInfinity checks if the given point represents the infinity point on the curve.
func (p AffinePoint) IsInfinity() bool {
	return p.X.Sign() == 0 && p.Y.Sign() == 0
}

// IsInfinity checks if the given point represents the infinity point on the curve.
func (p ProjectivePoint) IsInfinity() bool {
	return p.X.Sign() == 0 && p.Y.Sign() == 0 && p.Z.Sign() == 0
}

// InverseMod returns the inverse of the given scalar (fr) modulo m.
// ! This function assumes that m is a prime number.
// TODO: Fix this
func InverseMod(fr, m *big.Int) *big.Int {
	output := new(big.Int).Exp(fr, new(big.Int).Sub(m, big.NewInt(2)), m)
	return output
}

// powSmall calculates the modular exponentiation of 'base' raised to the power 'exponent' modulo 'modulus'
// using the binary exponentiation method.
func PowSmall(base *big.Int, exponent *big.Int, modulus *big.Int) *big.Int {
	result := big.NewInt(1)
	input := new(big.Int).Set(base)
	count := big.NewInt(1)

	for count.Cmp(exponent) <= 0 {
		if new(big.Int).And(exponent, count).Cmp(big.NewInt(0)) != 0 {
			result.Mul(result, input)
			result.Mod(result, modulus)
		}
		input.Mul(input, input)
		input.Mod(input, modulus)
		count.Mul(count, big.NewInt(2))
	}

	return result
}
