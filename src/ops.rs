use crate::{polynomial, Polynomial, Sized};
use auto_impl_ops::auto_ops;
use num_traits::{One, Zero};
use ring_algorithm::RingNormalize;
use std::{
    iter::{Product, Sum},
    ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Rem, RemAssign, Sub, SubAssign},
};

// additive monoid
#[auto_ops]
impl<M> AddAssign<&Polynomial<M>> for Polynomial<M>
where
    M: Sized + Zero + for<'x> AddAssign<&'x M>,
{
    fn add_assign(&mut self, other: &Self) {
        let len = self.len();
        self.extend(other.len());
        self.coeff
            .iter_mut()
            .zip(other.coeff.iter())
            .for_each(|(l, r)| *l += r);
        if len == other.len() {
            self.trim_zero()
        }
    }
}

impl<M> Zero for Polynomial<M>
where
    M: Sized + Zero + for<'x> AddAssign<&'x M>,
{
    fn zero() -> Self {
        Self { coeff: Vec::new() }
    }
    fn is_zero(&self) -> bool {
        self.deg().is_none()
    }
}

impl<M> Sum for Polynomial<M>
where
    M: Sized + Zero + for<'x> AddAssign<&'x M>,
{
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::zero(), Add::add)
    }
}

// additive group
impl<G> Neg for Polynomial<G>
where
    G: Sized + Neg<Output = G>,
{
    type Output = Self;
    fn neg(self) -> Self::Output {
        Polynomial {
            coeff: self.coeff.into_iter().map(|v| -v).collect(),
        }
    }
}

impl<G> Neg for &Polynomial<G>
where
    G: Sized + for<'x> From<<&'x G as Neg>::Output>,
    for<'x> &'x G: Neg,
{
    type Output = Polynomial<G>;
    fn neg(self) -> Self::Output {
        Polynomial {
            coeff: self.coeff.iter().map(|v| G::from(-v)).collect(),
        }
    }
}

#[auto_ops]
impl<G> SubAssign<&Polynomial<G>> for Polynomial<G>
where
    G: Sized + Zero + for<'x> SubAssign<&'x G>,
{
    fn sub_assign(&mut self, other: &Self) {
        let len = self.len();
        self.extend(other.len());
        self.coeff
            .iter_mut()
            .zip(other.coeff.iter())
            .for_each(|(l, r)| *l -= r);
        if len == other.len() {
            self.trim_zero()
        }
    }
}

// unitary ring
fn mul_aux<R>(sum: &mut [R], coeff: &R, vec: &[R])
where
    R: Sized + Clone + for<'x> AddAssign<&'x R> + for<'x> MulAssign<&'x R>,
{
    sum.iter_mut().zip(vec.iter()).for_each(|(l, r)| {
        let mut t = coeff.clone();
        t *= r;
        *l += &t
    });
}
#[auto_ops]
impl<R> Mul for &Polynomial<R>
where
    R: Sized + Clone + Zero + for<'x> AddAssign<&'x R> + for<'x> MulAssign<&'x R>,
{
    type Output = Polynomial<R>;
    fn mul(self, other: Self) -> Self::Output {
        if self.is_zero() || other.is_zero() {
            return Polynomial::<R>::zero();
        }
        let mut coeff = vec![R::zero(); self.len() + other.len() - 1];
        self.coeff
            .iter()
            .enumerate()
            .for_each(|(i, c)| mul_aux::<R>(&mut coeff[i..], c, &other.coeff));
        Polynomial::<R>::new(coeff) // R may not be a domain.
    }
}

impl<R> One for Polynomial<R>
where
    R: Sized + Clone + Zero + for<'x> AddAssign<&'x R> + One + for<'x> MulAssign<&'x R>,
{
    fn one() -> Self {
        polynomial![R::one()]
    }
}

impl<R> Product for Polynomial<R>
where
    R: Sized + Clone + Zero + for<'x> AddAssign<&'x R> + One + for<'x> MulAssign<&'x R>,
{
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::one(), Mul::<Polynomial<R>>::mul)
    }
}

// field
#[auto_ops]
impl<K> Div for &Polynomial<K>
where
    K: Sized
        + Clone
        + Zero
        + for<'x> AddAssign<&'x K>
        + for<'x> SubAssign<&'x K>
        + for<'x> MulAssign<&'x K>
        + for<'x> DivAssign<&'x K>,
{
    type Output = Polynomial<K>;
    fn div(self, other: Self) -> Self::Output {
        let mut f = self.clone();
        f.division(other)
    }
}

#[auto_ops]
impl<K> RemAssign<&Polynomial<K>> for Polynomial<K>
where
    K: Sized
        + Clone
        + Zero
        + for<'x> AddAssign<&'x K>
        + for<'x> SubAssign<&'x K>
        + for<'x> MulAssign<&'x K>
        + for<'x> DivAssign<&'x K>,
{
    fn rem_assign(&mut self, other: &Self) {
        self.division(other);
    }
}

// scalar ops
#[auto_ops]
impl<R> MulAssign<&R> for Polynomial<R>
where
    R: Sized + Zero + for<'x> MulAssign<&'x R>,
{
    fn mul_assign(&mut self, alpha: &R) {
        self.coeff.iter_mut().for_each(|c| *c *= alpha);
        self.trim_zero();
    }
}

#[auto_ops]
impl<R> DivAssign<&R> for Polynomial<R>
where
    R: Sized + for<'x> DivAssign<&'x R>,
{
    fn div_assign(&mut self, alpha: &R) {
        self.coeff.iter_mut().for_each(|c| *c /= alpha);
    }
}

impl<K> RingNormalize for Polynomial<K>
where
    K: Sized
        + Clone
        + Zero
        + One
        + for<'x> AddAssign<&'x K>
        + for<'x> MulAssign<&'x K>
        + for<'x> DivAssign<&'x K>,
{
    fn leading_unit(&self) -> Self {
        if let Some(x) = self.lc() {
            Self::from_monomial(x.clone(), 0)
        } else {
            Self::one()
        }
    }
    fn normalize_mut(&mut self) {
        self.monic();
    }
}
