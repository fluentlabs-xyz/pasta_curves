//! This module provides common utilities, traits and structures for group and
//! field arithmetic.
//!
//! This module is temporary, and the extension traits defined here are expected to be
//! upstreamed into the `ff` and `group` crates after some refactoring.

mod curves;
mod fields;

pub use curves::*;
pub(crate) use fields::*;

use subtle::Choice;

#[repr(transparent)]
pub struct Compat<T>(T);

impl<const N: u8, T: ff::WithSmallOrderMulGroup<N>> FieldExt for T {

    /// Modulus of the field written as a string for display purposes
    const MODULUS: &'static str = <T as ff::PrimeField>::MODULUS;

    /// Inverse of `PrimeField::root_of_unity()`
    const ROOT_OF_UNITY_INV: Self = <T as ff::PrimeField>::ROOT_OF_UNITY_INV;

    /// Generator of the $t-order$ multiplicative subgroup
    const DELTA: Self = <T as ff::PrimeField>::DELTA;

    /// Inverse of $2$ in the field.
    const TWO_INV: Self = <T as ff::PrimeField>::TWO_INV;

    /// Element of multiplicative order $3$.
    const ZETA: Self = <T as ff::WithSmallOrderMulGroup<N>>::ZETA;

    /// Obtains a field element congruent to the integer `v`.
    fn from_u128(v: u128) -> Self { <T as ff::PrimeField>::from_u128(v) }

    /// Obtains a field element that is congruent to the provided little endian
    /// byte representation of an integer.
    fn from_bytes_wide(bytes: &[u8; 64]) -> Self { unimplemented!() }

    /// Exponentiates `self` by `by`, where `by` is a little-endian order
    /// integer exponent.
    fn pow(&self, by: &[u64; 4]) -> Self {
        let mut res = Self::one();
        for e in by.iter().rev() {
            for i in (0..64).rev() {
                res = res.square();
                let mut tmp = res;
                tmp *= self;
                res.conditional_assign(&tmp, (((*e >> i) & 0x1) as u8).into());
            }
        }
        res
    }

    /// Gets the lower 128 bits of this field element when expressed
    /// canonically.
    fn get_lower_128(&self) -> u128 { unimplemented!() }

}

/// This represents an element of a group with basic operations that can be
/// performed. This allows an FFT implementation (for example) to operate
/// generically over either a field or elliptic curve group.
pub trait Group: Copy + Clone + Send + Sync + 'static {
    /// The group is assumed to be of prime order $p$. `Scalar` is the
    /// associated scalar field of size $p$.
    type Scalar: FieldExt;

    /// Returns the additive identity of the group.
    fn group_zero() -> Self;

    /// Adds `rhs` to this group element.
    fn group_add(&mut self, rhs: &Self);

    /// Subtracts `rhs` from this group element.
    fn group_sub(&mut self, rhs: &Self);

    /// Scales this group element by a scalar.
    fn group_scale(&mut self, by: &Self::Scalar);
}

impl<T: ff::PrimeField> Group for T {
}

/// A trait that exposes additional operations related to calculating square roots of
/// prime-order finite fields.
pub trait SqrtRatio: ff::PrimeField {
    /// The value $(T-1)/2$ such that $2^S \cdot T = p - 1$ with $T$ odd.
    const T_MINUS1_OVER2: [u64; 4];

    /// Raise this field element to the power [`Self::T_MINUS1_OVER2`].
    ///
    /// Field implementations may override this to use an efficient addition chain.
    fn pow_by_t_minus1_over2(&self) -> Self {
        ff::Field::pow_vartime(self, &Self::T_MINUS1_OVER2)
    }

    /// Gets the lower 32 bits of this field element when expressed
    /// canonically.
    fn get_lower_32(&self) -> u32 { unimplemented!() }

    fn sqrt_ratio(num: &Self, div: &Self) -> (Choice, Self) { <Self as ff::Field>::sqrt_ratio(num, div) }

    fn sqrt_alt(&self) -> (Choice, Self) { <Self as ff::Field>::sqrt_alt(self) }
}


/// This trait is a common interface for dealing with elements of a finite
/// field.
pub trait FieldExt: SqrtRatio + From<bool> + Ord + Group<Scalar = Self> {

    /// Modulus of the field written as a string for display purposes
    const MODULUS: &'static str;

    /// Inverse of `PrimeField::root_of_unity()`
    const ROOT_OF_UNITY_INV: Self;

    /// Generator of the $t-order$ multiplicative subgroup
    const DELTA: Self;

    /// Inverse of $2$ in the field.
    const TWO_INV: Self;

    /// Element of multiplicative order $3$.
    const ZETA: Self;

    /// Obtains a field element congruent to the integer `v`.
    fn from_u128(v: u128) -> Self;

    /// Obtains a field element that is congruent to the provided little endian
    /// byte representation of an integer.
    fn from_bytes_wide(bytes: &[u8; 64]) -> Self;

    /// Exponentiates `self` by `by`, where `by` is a little-endian order
    /// integer exponent.
    fn pow(&self, by: &[u64; 4]) -> Self {
        let mut res = Self::one();
        for e in by.iter().rev() {
            for i in (0..64).rev() {
                res = res.square();
                let mut tmp = res;
                tmp *= self;
                res.conditional_assign(&tmp, (((*e >> i) & 0x1) as u8).into());
            }
        }
        res
    }

    /// Gets the lower 128 bits of this field element when expressed
    /// canonically.
    fn get_lower_128(&self) -> u128;

}

