use crate::{Polynomial, Sized};
use num_traits::Zero;
use std::ops::{
    Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Rem, RemAssign, Sub, SubAssign,
};

// AddAssign
impl<'a, M> AddAssign<&'a Polynomial<M>> for Polynomial<M>
where
    M: Sized + Clone + Zero + for<'x> AddAssign<&'x M>,
{
    fn add_assign(&mut self, other: &Self) {
        self.add_assign_ref(other);
    }
}
impl<M> AddAssign for Polynomial<M>
where
    M: Sized + Clone + Zero + for<'x> AddAssign<&'x M>,
{
    fn add_assign(&mut self, other: Self) {
        *self += &other
    }
}

// Add
impl<'a, M> Add for &'a Polynomial<M>
where
    M: Sized + Clone + Zero + for<'x> AddAssign<&'x M>,
{
    type Output = Polynomial<M>;
    fn add(self, other: Self) -> Self::Output {
        let mut f = self.clone();
        f += other;
        f
    }
}
impl<'a, M> Add<Polynomial<M>> for &'a Polynomial<M>
where
    M: Sized + Clone + Zero + for<'x> AddAssign<&'x M>,
{
    type Output = Polynomial<M>;
    fn add(self, other: Polynomial<M>) -> Self::Output {
        let mut f = self.clone();
        f += &other;
        f
    }
}
impl<'a, M> Add<&'a Polynomial<M>> for Polynomial<M>
where
    M: Sized + Clone + Zero + for<'x> AddAssign<&'x M>,
{
    type Output = Self;
    fn add(mut self, other: &Self) -> Self::Output {
        self += other;
        self
    }
}
impl<M> Add for Polynomial<M>
where
    M: Sized + Clone + Zero + for<'x> AddAssign<&'x M>,
{
    type Output = Self;
    fn add(mut self, other: Polynomial<M>) -> Self::Output {
        self += &other;
        self
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
impl<'a, G> Neg for &'a Polynomial<G>
where
    G: Sized,
    for<'x> &'x G: Neg<Output = G>,
{
    type Output = Polynomial<G>;
    fn neg(self) -> Self::Output {
        self.neg_ref()
    }
}

// SubAssign
impl<'a, G> SubAssign<&'a Polynomial<G>> for Polynomial<G>
where
    G: Sized + Clone + Zero + for<'x> SubAssign<&'x G>,
{
    fn sub_assign(&mut self, other: &Self) {
        self.sub_assign_ref(other)
    }
}
impl<G> SubAssign for Polynomial<G>
where
    G: Sized + Clone + Zero + for<'x> SubAssign<&'x G>,
{
    fn sub_assign(&mut self, other: Self) {
        *self -= &other
    }
}

// Sub
impl<'a, G> Sub for &'a Polynomial<G>
where
    G: Sized + Clone + Zero + for<'x> SubAssign<&'x G>,
{
    type Output = Polynomial<G>;
    fn sub(self, other: Self) -> Self::Output {
        let mut f = self.clone();
        f -= other;
        f
    }
}
impl<'a, G> Sub<Polynomial<G>> for &'a Polynomial<G>
where
    G: Sized + Clone + Zero + for<'x> SubAssign<&'x G>,
{
    type Output = Polynomial<G>;
    fn sub(self, other: Polynomial<G>) -> Self::Output {
        let mut f = self.clone();
        f -= &other;
        f
    }
}
impl<'a, G> Sub<&'a Polynomial<G>> for Polynomial<G>
where
    G: Sized + Clone + Zero + for<'x> SubAssign<&'x G>,
{
    type Output = Self;
    fn sub(mut self, other: &Self) -> Self::Output {
        self -= other;
        self
    }
}
impl<G> Sub for Polynomial<G>
where
    G: Sized + Clone + Zero + for<'x> SubAssign<&'x G>,
{
    type Output = Self;
    fn sub(mut self, other: Polynomial<G>) -> Self::Output {
        self -= &other;
        self
    }
}

// Mul
impl<'a, R> Mul for &'a Polynomial<R>
where
    R: Sized + Clone + Zero + for<'x> AddAssign<&'x R>,
    for<'x> &'x R: Mul<Output = R>,
{
    type Output = Polynomial<R>;
    fn mul(self, other: Self) -> Self::Output {
        self.mul_impl(other)
    }
}
impl<'a, R> Mul<Polynomial<R>> for &'a Polynomial<R>
where
    R: Sized + Clone + Zero + for<'x> AddAssign<&'x R>,
    for<'x> &'x R: Mul<Output = R>,
{
    type Output = Polynomial<R>;
    fn mul(self, other: Polynomial<R>) -> Self::Output {
        self * &other
    }
}
impl<R> Mul for Polynomial<R>
where
    R: Sized + Clone + Zero + for<'x> AddAssign<&'x R>,
    for<'x> &'x R: Mul<Output = R>,
{
    type Output = Self;
    fn mul(self, other: Self) -> Self::Output {
        &self * &other
    }
}
impl<'a, R> Mul<&'a Polynomial<R>> for Polynomial<R>
where
    R: Sized + Clone + Zero + for<'x> AddAssign<&'x R>,
    for<'x> &'x R: Mul<Output = R>,
{
    type Output = Self;
    fn mul(self, other: &Self) -> Self::Output {
        &self * other
    }
}

// MulAssign
impl<'a, R> MulAssign<&'a Polynomial<R>> for Polynomial<R>
where
    R: Sized + Clone + Zero + for<'x> AddAssign<&'x R>,
    for<'x> &'x R: Mul<Output = R>,
{
    fn mul_assign(&mut self, other: &Self) {
        *self = &*self * other;
    }
}
impl<R> MulAssign for Polynomial<R>
where
    R: Sized + Clone + Zero + for<'x> AddAssign<&'x R>,
    for<'x> &'x R: Mul<Output = R>,
{
    fn mul_assign(&mut self, other: Self) {
        *self = &*self * &other;
    }
}

// Div
impl<'a, K> Div for &'a Polynomial<K>
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
impl<'a, K> Div<Polynomial<K>> for &'a Polynomial<K>
where
    K: Sized + Clone + Zero + for<'x> AddAssign<&'x K> + for<'x> SubAssign<&'x K>,
    for<'x> &'x K: Mul<Output = K> + Div<Output = K>,
{
    type Output = Polynomial<K>;
    fn div(self, other: Polynomial<K>) -> Self::Output {
        let mut f = self.clone();
        f.division(&other)
    }
}
impl<K> Div for Polynomial<K>
where
    K: Sized + Clone + Zero + for<'x> AddAssign<&'x K> + for<'x> SubAssign<&'x K>,
    for<'x> &'x K: Mul<Output = K> + Div<Output = K>,
{
    type Output = Self;
    fn div(mut self, other: Self) -> Self::Output {
        self.division(&other)
    }
}
impl<'a, K> Div<&'a Polynomial<K>> for Polynomial<K>
where
    K: Sized + Clone + Zero + for<'x> AddAssign<&'x K> + for<'x> SubAssign<&'x K>,
    for<'x> &'x K: Mul<Output = K> + Div<Output = K>,
{
    type Output = Self;
    fn div(mut self, other: &Self) -> Self::Output {
        self.division(other)
    }
}

// DivAssign
impl<'a, K> DivAssign<&'a Polynomial<K>> for Polynomial<K>
where
    K: Sized + Clone + Zero + for<'x> AddAssign<&'x K> + for<'x> SubAssign<&'x K>,
    for<'x> &'x K: Mul<Output = K> + Div<Output = K>,
{
    fn div_assign(&mut self, other: &Self) {
        *self = &*self / other;
    }
}
impl<K> DivAssign for Polynomial<K>
where
    K: Sized + Clone + Zero + for<'x> AddAssign<&'x K> + for<'x> SubAssign<&'x K>,
    for<'x> &'x K: Mul<Output = K> + Div<Output = K>,
{
    fn div_assign(&mut self, other: Self) {
        *self = &*self / &other;
    }
}

// RemAssign
impl<'a, K> RemAssign<&'a Polynomial<K>> for Polynomial<K>
where
    K: Sized + Clone + Zero + for<'x> AddAssign<&'x K> + for<'x> SubAssign<&'x K>,
    for<'x> &'x K: Mul<Output = K> + Div<Output = K>,
{
    fn rem_assign(&mut self, other: &Self) {
        self.division(other);
    }
}
impl<'a, K> RemAssign for Polynomial<K>
where
    K: Sized + Clone + Zero + for<'x> AddAssign<&'x K> + for<'x> SubAssign<&'x K>,
    for<'x> &'x K: Mul<Output = K> + Div<Output = K>,
{
    fn rem_assign(&mut self, other: Self) {
        self.division(&other);
    }
}

// Rem
impl<'a, K> Rem for &'a Polynomial<K>
where
    K: Sized + Clone + Zero + for<'x> AddAssign<&'x K> + for<'x> SubAssign<&'x K>,
    for<'x> &'x K: Mul<Output = K> + Div<Output = K>,
{
    type Output = Polynomial<K>;
    fn rem(self, other: Self) -> Self::Output {
        let mut t = self.clone();
        t %= other;
        t
    }
}
impl<'a, K> Rem<Polynomial<K>> for &'a Polynomial<K>
where
    K: Sized + Clone + Zero + for<'x> AddAssign<&'x K> + for<'x> SubAssign<&'x K>,
    for<'x> &'x K: Mul<Output = K> + Div<Output = K>,
{
    type Output = Polynomial<K>;
    fn rem(self, other: Polynomial<K>) -> Self::Output {
        let mut t = self.clone();
        t %= other;
        t
    }
}
impl<'a, K> Rem<&'a Polynomial<K>> for Polynomial<K>
where
    K: Sized + Clone + Zero + for<'x> AddAssign<&'x K> + for<'x> SubAssign<&'x K>,
    for<'x> &'x K: Mul<Output = K> + Div<Output = K>,
{
    type Output = Self;
    fn rem(mut self, other: &Self) -> Self::Output {
        self %= other;
        self
    }
}
impl<K> Rem for Polynomial<K>
where
    K: Sized + Clone + Zero + for<'x> AddAssign<&'x K> + for<'x> SubAssign<&'x K>,
    for<'x> &'x K: Mul<Output = K> + Div<Output = K>,
{
    type Output = Self;
    fn rem(mut self, other: Self) -> Self::Output {
        self %= &other;
        self
    }
}
