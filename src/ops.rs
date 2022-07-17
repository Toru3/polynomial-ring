use crate::{Polynomial, Sized};
use auto_impl_ops::auto_ops;
use num_traits::Zero;
use std::ops::{
    Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Rem, RemAssign, Sub, SubAssign,
};

#[auto_ops]
impl<M> AddAssign<&Polynomial<M>> for Polynomial<M>
where
    M: Sized + Zero + for<'x> AddAssign<&'x M>,
{
    fn add_assign(&mut self, other: &Self) {
        let len = self.len();
        self.extend(other.len());
        self.coef
            .iter_mut()
            .zip(other.coef.iter())
            .for_each(|(l, r)| *l += r);
        if len == other.len() {
            self.trim_zero()
        }
    }
}

// Neg
impl<G> Neg for Polynomial<G>
where
    G: Sized + Neg<Output = G>,
{
    type Output = Self;
    fn neg(self) -> Self::Output {
        Polynomial {
            coef: self.coef.into_iter().map(|v| -v).collect(),
        }
    }
}
impl<G> Neg for &Polynomial<G>
where
    G: Sized,
    for<'x> &'x G: Neg<Output = G>,
{
    type Output = Polynomial<G>;
    fn neg(self) -> Self::Output {
        Polynomial {
            coef: self.coef.iter().map(|v| -v).collect(),
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
        self.coef
            .iter_mut()
            .zip(other.coef.iter())
            .for_each(|(l, r)| *l -= r);
        if len == other.len() {
            self.trim_zero()
        }
    }
}

fn mul_aux<R>(sum: &mut [R], coef: &R, vec: &[R])
where
    R: Sized + for<'x> AddAssign<&'x R>,
    for<'x> &'x R: Mul<Output = R>,
{
    sum.iter_mut()
        .zip(vec.iter())
        .for_each(|(l, r)| *l += &(coef * r));
}
#[auto_ops]
impl<R> Mul for &Polynomial<R>
where
    R: Sized + Clone + Zero + for<'x> AddAssign<&'x R>,
    for<'x> &'x R: Mul<Output = R>,
{
    type Output = Polynomial<R>;
    fn mul(self, other: Self) -> Self::Output {
        if self.is_zero() || other.is_zero() {
            return Self::zero();
        }
        let mut coef = vec![R::zero(); self.len() + other.len() - 1];
        self.coef
            .iter()
            .enumerate()
            .for_each(|(i, c)| mul_aux::<R>(&mut coef[i..], c, &other.coef));
        Polynomial::<R>::new(coef) // R may not be a domain.
    }
}

#[auto_ops]
impl<K> Div for &Polynomial<K>
where
    K: Sized + Clone + Zero + for<'x> AddAssign<&'x K> + for<'x> SubAssign<&'x K>,
    for<'x> &'x K: Mul<Output = K> + Div<Output = K>,
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
    K: Sized + Clone + Zero + for<'x> AddAssign<&'x K> + for<'x> SubAssign<&'x K>,
    for<'x> &'x K: Mul<Output = K> + Div<Output = K>,
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
    fn mul_assign(&mut self, other: &R) {
        self.coef.iter_mut().for_each(|c| *c *= alpha);
        self.trim_zero();
    }
}

#[auto_ops]
impl<R> DivAssign<&R> for Polynomial<R>
where
    R: Sized + Zero + for<'x> DivAssign<&'x R>,
{
    fn div_assign(&mut self, other: &R) {
        self.coef.iter_mut().for_each(|c| *c /= alpha);
    }
}
