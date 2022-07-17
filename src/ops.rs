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
        self.add_assign_ref(other)
    }
}

// Neg
impl<G> Neg for Polynomial<G>
where
    G: Sized + Neg<Output = G>,
{
    type Output = Self;
    fn neg(self) -> Self::Output {
        self.neg_impl()
    }
}
impl<G> Neg for &Polynomial<G>
where
    G: Sized,
    for<'x> &'x G: Neg<Output = G>,
{
    type Output = Polynomial<G>;
    fn neg(self) -> Self::Output {
        self.neg_ref()
    }
}

#[auto_ops]
impl<G> SubAssign<&Polynomial<G>> for Polynomial<G>
where
    G: Sized + Zero + for<'x> SubAssign<&'x G>,
{
    fn sub_assign(&mut self, other: &Self) {
        self.sub_assign_ref(other)
    }
}

#[auto_ops]
impl<R> Mul for &Polynomial<R>
where
    R: Sized + Clone + Zero + for<'x> AddAssign<&'x R>,
    for<'x> &'x R: Mul<Output = R>,
{
    type Output = Polynomial<R>;
    fn mul(self, other: Self) -> Self::Output {
        self.mul_impl(other)
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
        self.scalar_mul_assign_impl(other);
    }
}

#[auto_ops]
impl<R> DivAssign<&R> for Polynomial<R>
where
    R: Sized + Zero + for<'x> DivAssign<&'x R>,
{
    fn div_assign(&mut self, other: &R) {
        self.scalar_div_assign_impl(other);
    }
}
