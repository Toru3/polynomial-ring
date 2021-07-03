/*! Polynomial ring $`R[x]`$

```
use num::Rational;
use polynomial_ring::Polynomial;
let p = Polynomial::new(vec![3, 1, 4, 1, 5].into_iter().map(|x| Rational::from_integer(x)).collect());
let q = Polynomial::new(vec![2, 7, 1].into_iter().map(|x| Rational::from_integer(x)).collect());
let mut r = p.clone();
let d = r.division(&q);
assert_eq!(p, d * q + r);
```
*/
#![cfg_attr(feature = "__internal_inject_debug", recursion_limit = "8")]
mod sealed {
    pub trait SizedExt: std::marker::Sized + std::fmt::Debug + std::fmt::Display {}
    impl<T> SizedExt for T where T: std::marker::Sized + std::fmt::Debug + std::fmt::Display {}
    #[cfg(not(feature = "__internal_inject_debug"))]
    pub use std::marker::Sized;
    #[cfg(feature = "__internal_inject_debug")]
    pub use SizedExt as Sized;
}
use num_traits::{One, Zero};
use ring_algorithm::RingNormalize;
use sealed::Sized;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, SubAssign};
mod ops;

/** Polynomial ring $`R[x]`$

```
use num::Rational;
use polynomial_ring::Polynomial;
let p = Polynomial::new(vec![3, 1, 4, 1, 5].into_iter().map(|x| Rational::from_integer(x)).collect());
let q = Polynomial::new(vec![2, 7, 1].into_iter().map(|x| Rational::from_integer(x)).collect());
let mut r = p.clone();
let d = r.division(&q);
assert_eq!(p, d * q + r);
```
*/
#[derive(Clone, Debug, PartialEq, Eq, Default)]
pub struct Polynomial<T> {
    coef: Vec<T>,
}
impl<T: crate::Sized> Polynomial<T> {
    fn len(&self) -> usize {
        self.coef.len()
    }
    /** degree of polynomial

    ```
    use polynomial_ring::Polynomial;
    let p = Polynomial::new(vec![3, 2, 1]); // 3+2x+x^2
    assert_eq!(p.deg(), Some(2));
    let q = Polynomial::new(vec![0]); // 0
    assert_eq!(q.deg(), None);
    ```
    */
    pub fn deg(&self) -> Option<usize> {
        if self.coef.is_empty() {
            None
        } else {
            Some(self.len() - 1)
        }
    }
    /** leading coefficent

    ```
    use polynomial_ring::Polynomial;
    let p = Polynomial::new(vec![3, 2, 1]); // 3+2x+x^2
    assert_eq!(p.lc(), Some(&1));
    let q = Polynomial::new(vec![0]); // 0
    assert_eq!(q.lc(), None);
    ```
    */
    pub fn lc(&self) -> Option<&T> {
        self.deg().map(|d| &self.coef[d])
    }
    /** get coefficents

    ```
    use polynomial_ring::Polynomial;
    let p = Polynomial::new(vec![3, 2, 1]); // 3+2x+x^2
    assert_eq!(p.coefs(), vec![3, 2, 1]);
    let q = Polynomial::new(vec![0]); // 0
    assert_eq!(q.coefs(), vec![]);
    ```
    */
    pub fn coefs(self) -> Vec<T> {
        self.coef
    }
}

// additive monoid
impl<M> Polynomial<M> {
    fn add_assign_ref(&mut self, other: &Self)
    where
        M: Sized + Clone + Zero + for<'x> AddAssign<&'x M>,
    {
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
impl<M> Zero for Polynomial<M>
where
    M: Sized + Clone + Zero + for<'x> AddAssign<&'x M>,
{
    fn zero() -> Self {
        Self { coef: Vec::new() }
    }
    fn is_zero(&self) -> bool {
        self.deg().is_none()
    }
}
impl<M> std::iter::Sum for Polynomial<M>
where
    M: Sized + Clone + Zero + for<'x> AddAssign<&'x M>,
{
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::zero(), Add::add)
    }
}
impl<M: crate::Sized + Zero> Polynomial<M> {
    fn trim_zero(&mut self) {
        let len = self
            .coef
            .iter()
            .rposition(|x| !x.is_zero())
            .map(|pos| pos + 1)
            .unwrap_or(0);
        self.coef.truncate(len);
    }
    /** construct polynomial

    ```
    use polynomial_ring::Polynomial;
    let p = Polynomial::new(vec![3, 2, 1]);
    assert_eq!(p.to_string(), "x^2+2*x+3");
    ```
    */
    pub fn new(coef: Vec<M>) -> Self {
        let mut poly = Self { coef };
        poly.trim_zero();
        poly
    }
    fn extend(&mut self, len: usize) {
        if self.len() < len {
            self.coef.resize_with(len, M::zero);
        }
    }
    /** construct polynomial from monomial $`cx^d`$ ($`c`$=coefficent, $`d`$=degree)

    ```
    use polynomial_ring::Polynomial;
    let p = Polynomial::from_monomial(3, 2);
    let q = Polynomial::new(vec![0, 0, 3]);
    assert_eq!(p, q);
    ```
    */
    pub fn from_monomial(coefficent: M, degree: usize) -> Self {
        let coef = if coefficent.is_zero() {
            Vec::new()
        } else {
            let mut coef = Vec::with_capacity(degree + 1);
            for _ in 0..degree {
                coef.push(M::zero());
            }
            coef.push(coefficent);
            coef
        };
        Self { coef }
    }
}
#[macro_export]
/** make Polynomial (like `vec!`)

```
use polynomial_ring::{Polynomial, polynomial};
let p = Polynomial::new(vec![3, 2, 1]);
let q = polynomial![3, 2, 1];
assert_eq!(p, q);
```
*/
macro_rules! polynomial {
    ($($x:expr),*) => {
        Polynomial::new(vec![$($x), *])
    }
}

// additive group
impl<G> Polynomial<G> {
    fn neg_impl(self) -> Self
    where
        G: Sized + Neg<Output = G>,
    {
        Polynomial {
            coef: self.coef.into_iter().map(|v| -v).collect(),
        }
    }
    fn neg_ref(&self) -> Self
    where
        G: Sized,
        for<'x> &'x G: Neg<Output = G>,
    {
        Polynomial {
            coef: self.coef.iter().map(|v| -v).collect(),
        }
    }
    fn sub_assign_ref(&mut self, other: &Self)
    where
        G: Sized + Clone + Zero + for<'x> SubAssign<&'x G>,
    {
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

// unitary ring
fn mul_aux<R>(sum: &mut [R], coef: &R, vec: &[R])
where
    R: Sized + Clone + Zero + for<'x> AddAssign<&'x R>,
    for<'x> &'x R: Mul<Output = R>,
{
    sum.iter_mut()
        .zip(vec.iter())
        .for_each(|(l, r)| *l += &(coef * r));
}
impl<R> Polynomial<R> {
    fn mul_impl(&self, other: &Self) -> Self
    where
        R: Sized + Clone + Zero + for<'x> AddAssign<&'x R>,
        for<'x> &'x R: Mul<Output = R>,
    {
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
impl<R> One for Polynomial<R>
where
    R: Sized + Clone + Zero + for<'x> AddAssign<&'x R> + One,
    for<'x> &'x R: Mul<Output = R>,
{
    fn one() -> Self {
        polynomial![R::one()]
    }
}
impl<R> std::iter::Product for Polynomial<R>
where
    R: Sized + Clone + Zero + for<'x> AddAssign<&'x R> + One,
    for<'x> &'x R: Mul<Output = R>,
{
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::one(), Mul::mul)
    }
}
impl<R> std::fmt::Display for Polynomial<R>
where
    R: std::cmp::Eq + std::fmt::Display + Zero + One,
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let vec = &self.coef;
        if vec.is_empty() {
            return write!(f, "0");
        }
        let mut is_first = true;
        for (i, c) in vec.iter().enumerate().rev() {
            if c.is_zero() {
                continue;
            }
            if is_first {
                is_first = false;
            } else {
                write!(f, "+")?
            }
            if c.is_one() {
                match i {
                    0 => write!(f, "1")?,
                    1 => write!(f, "x")?,
                    _ => write!(f, "x^{}", i)?,
                }
            } else {
                match i {
                    0 => write!(f, "{}", c)?,
                    1 => write!(f, "{}*x", c)?,
                    _ => write!(f, "{}*x^{}", c, i)?,
                }
            }
        }
        Ok(())
    }
}
impl<R: Sized> Polynomial<R> {
    /** evaluate polynomial by Horner's method

    ```
    use polynomial_ring::Polynomial;
    let p = Polynomial::new(vec![3, 2, 1]); // 3+2x+x^2
    assert_eq!(p.eval(&1), 6);
    assert_eq!(p.eval(&2), 11);
    ```
    */
    pub fn eval<'a>(&self, x: &'a R) -> R
    where
        R: Sized + Clone + Zero + for<'x> AddAssign<&'x R> + MulAssign<&'a R>,
    {
        if self.coef.is_empty() {
            return R::zero();
        }
        let mut sum = self.lc().unwrap().clone();
        for i in (0..self.len() - 1).rev() {
            sum *= x;
            sum += &self.coef[i];
        }
        sum
    }
    /** derivative

    ```
    use polynomial_ring::{Polynomial, polynomial};
    let p = polynomial![1, 2, 3, 2, 1]; // 1+2x+3x^2+2x^3+x^4
    assert_eq!(p.derivative(), polynomial![2, 6, 6, 4]);
    ```
    */
    pub fn derivative(self) -> Self
    where
        R: Sized + Clone + Zero + One + for<'x> AddAssign<&'x R> + Mul<Output = R>,
        for<'x> &'x R: Mul<Output = R>,
    {
        let n = self.coef.len();
        let n = if n > 0 { n - 1 } else { 0 };
        let mut coef = Vec::with_capacity(n);
        let mut i = R::one();
        for c in self.coef.into_iter().skip(1) {
            coef.push(&i * &c);
            i += &R::one();
        }
        Polynomial::new(coef)
    }
}

// field
impl<K> RingNormalize for Polynomial<K>
where
    K: Sized + Clone + Zero + One + for<'x> AddAssign<&'x K> + for<'x> DivAssign<&'x K>,
    for<'x> &'x K: Mul<Output = K>,
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
impl<K: Sized> Polynomial<K> {
    /** make polynomial monic

    ```
    use num::Rational;
    use polynomial_ring::Polynomial;
    let mut p = Polynomial::new(vec![1, 2, 3].into_iter().map(|x| Rational::from_integer(x)).collect());
    p.monic();
    let q = Polynomial::new(vec![(1, 3), (2, 3), (1, 1)].into_iter().map(|(n, d)| Rational::new(n, d)).collect());
    assert_eq!(p, q);
    ```
    */
    pub fn monic(&mut self)
    where
        K: Clone + for<'x> DivAssign<&'x K>,
    {
        if let Some(lc) = self.lc() {
            let lc = lc.clone();
            self.coef.iter_mut().for_each(|v| *v /= &lc);
        }
    }
    /** polynomial division

    ```
    use num::Rational;
    use polynomial_ring::Polynomial;
    let p = Polynomial::new(vec![3, 1, 4, 1, 5].into_iter().map(|x| Rational::from_integer(x)).collect());
    let q = Polynomial::new(vec![2, 7, 1].into_iter().map(|x| Rational::from_integer(x)).collect());
    let mut r = p.clone();
    let d = r.division(&q);
    assert_eq!(p, d * q + r);
    ```
    */
    pub fn division(&mut self, other: &Self) -> Self
    where
        K: Sized + Clone + Zero + for<'x> AddAssign<&'x K> + for<'x> SubAssign<&'x K>,
        for<'x> &'x K: Mul<Output = K> + Div<Output = K>,
    {
        let g_deg = other.deg().expect("Division by zero");
        if self.deg() < other.deg() {
            return Self::zero();
        }
        let lc = other.lc().unwrap();
        let mut coef = vec![K::zero(); self.len() - other.len() + 1];
        while self.deg() >= other.deg() {
            let d = self.deg().unwrap() - g_deg;
            let c = self.lc().unwrap() / lc;
            for i in 0..other.len() - 1 {
                self.coef[i + d] -= &(&c * &other.coef[i]);
            }
            self.coef.pop(); // new deg < prev deg
            self.trim_zero();
            coef[d] = c;
        }
        Self { coef }
    }
    /** square free

    ```
    use polynomial_ring::{Polynomial, polynomial};
    use num::Rational;
    let f = polynomial![Rational::from(1), Rational::from(1)];
    let g = polynomial![Rational::from(1), Rational::from(1), Rational::from(1)];
    let p = &f * &f * &f * &g * &g; // (x+1)^3(x^2+x+1)^2
    assert_eq!(p.square_free(), &f * &g); // (x+1)(x^2+x+1)
    ```
    */
    pub fn square_free(&self) -> Self
    where
        K: Sized
            + Clone
            + Zero
            + One
            + Mul<Output = K>
            + for<'x> AddAssign<&'x K>
            + for<'x> SubAssign<&'x K>
            + for<'x> DivAssign<&'x K>,
        for<'x> &'x K: Mul<Output = K> + Div<Output = K>,
    {
        let d = self.clone().derivative().into_normalize();
        let f = ring_algorithm::gcd::<Self>(self.clone(), d).into_normalize();
        (self / &f).into_normalize()
    }
}
