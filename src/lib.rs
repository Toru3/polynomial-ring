#![doc = include_str!("../README.md")]
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
use ring_algorithm::{gcd, RingNormalize};
use sealed::Sized;
use std::fmt;
use std::ops::{AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Rem, Sub, SubAssign};
mod ops;

/** Polynomial ring $`R[x]`$

```
use num::Rational64;
use polynomial_ring::Polynomial;
let f = Polynomial::new(vec![3, 1, 4, 1, 5].into_iter().map(|x| Rational64::from_integer(x)).collect());
let g = Polynomial::new(vec![2, 7, 1].into_iter().map(|x| Rational64::from_integer(x)).collect());
let mut r = f.clone();
let q = r.division(&g);
assert_eq!(f, q * g + r);
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

    /** The degree of the polynomial.
    This is the highest power, or [None] for the zero polynomial.

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

    /** The leading coefficent.
    This is the coeffient of the highest power, or [None] for the zero polynomial.

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

    /** Get the coefficents.

    ```
    use polynomial_ring::Polynomial;
    let p = Polynomial::new(vec![3, 2, 1]); // 3+2x+x^2
    assert_eq!(p.coeffs(), &[3, 2, 1]);
    let q = Polynomial::new(vec![0]); // 0
    assert_eq!(q.coeffs(), &[]);
    ```
    */
    pub fn coeffs(&self) -> &[T] {
        &self.coef
    }
}

// additive monoid
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

    /** Construct a polynomial from a [Vec] of coefficients.

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

    /** Construct a monomial $`cx^d`$ from a coefficient and a degree ($`c`$=coefficent, $`d`$=degree).

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
/** Make a [Polynomial] from a list of coefficients (like `vec!`).

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

// unitary ring
impl<R> fmt::Display for Polynomial<R>
where
    R: PartialEq + fmt::Display + Zero + One,
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
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
                write!(f, "+")?;
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
    /** Evaluate a polynomial at some point, using Horner's method.

    ```
    use polynomial_ring::Polynomial;
    let p = Polynomial::new(vec![3, 2, 1]); // 3+2x+x^2
    assert_eq!(p.eval(&1), 6);
    assert_eq!(p.eval(&2), 11);
    ```
    */
    pub fn eval<'a>(&self, x: &'a R) -> R
    where
        R: Sized + Zero + for<'x> AddAssign<&'x R> + MulAssign<&'a R>,
    {
        let mut sum = R::zero();
        for coeff in self.coef.iter().rev() {
            sum *= x;
            sum += coeff;
        }
        sum
    }

    /** Calculate the derivative.

    ```
    use polynomial_ring::{Polynomial, polynomial};
    let p = polynomial![1, 2, 3, 2, 1]; // 1+2x+3x^2+2x^3+x^4
    assert_eq!(p.derivative(), polynomial![2, 6, 6, 4]); // 2+6x+6x^2+4x^3
    ```
    */
    #[must_use]
    pub fn derivative(&self) -> Self
    where
        R: Sized + Zero + One + for<'x> AddAssign<&'x R>,
        for<'x> &'x R: Mul<Output = R>,
    {
        let n = self.coef.len();
        let n = if n > 0 { n - 1 } else { 0 };
        let mut coef = Vec::with_capacity(n);
        let mut i = R::one();
        for c in self.coef.iter().skip(1) {
            coef.push(&i * c);
            i += &R::one();
        }
        Polynomial::new(coef)
    }

    /** Pseudo division.

    Let $`R`$ be an [integral domain](https://en.wikipedia.org/wiki/Integral_domain).
    Let $`f, g \in R[x]`$, where $`g \neq 0`$.
    This function calculates $`s \in R`$, $`q, r \in R[x]`$ s.t. $`sf=qg+r`$,
    where $`r=0`$ or $`\deg(r)<\deg(g)`$.
    ```
    use polynomial_ring::{polynomial, Polynomial};

    let f = polynomial![1, 3, 1]; // 1+3x+x^2 ∈ Z[x]
    let g = polynomial![5, 2]; // 5+2x ∈ Z[x]
    let mut r = f.clone();
    let (s, q) = r.pseudo_division(&g);
    assert_eq!(f * s, q * g + r);
    ```
    */
    pub fn pseudo_division(&mut self, other: &Self) -> (R, Self)
    where
        R: Sized + Clone + Zero + One + for<'x> AddAssign<&'x R> + for<'x> MulAssign<&'x R>,
        for<'x> &'x R: Sub<Output = R> + Mul<Output = R>,
    {
        let g_deg = other.deg().expect("Division by zero");
        let f_deg = self.deg();
        if f_deg < other.deg() {
            return (R::one(), Self::zero());
        }
        let f_deg = f_deg.unwrap();
        debug_assert!(f_deg >= g_deg);
        let k = f_deg - g_deg + 1;
        let lc = other.lc().unwrap();
        let mut coef = vec![R::zero(); k];
        let mut scale = R::one();
        while self.deg() >= other.deg() {
            let d = self.deg().unwrap() - g_deg;
            let c = self.lc().unwrap().clone();
            for i in 0..other.len() - 1 {
                self.coef[i + d] = &(lc * &self.coef[i + d]) - &(&c * &other.coef[i]);
            }
            for i in 0..d {
                self.coef[i] *= lc;
            }
            self.coef.pop(); // new deg < prev deg
            self.trim_zero();
            for c_i in coef.iter_mut().skip(d + 1) {
                *c_i *= lc;
            }
            coef[d] = c;
            scale *= lc;
        }
        (scale, Self { coef })
    }

    /** Calculate the [resultant](https://en.wikipedia.org/wiki/Resultant)

    ```
    use polynomial_ring::{polynomial, Polynomial};

    let f = polynomial![0, 2, 0, 1]; // 2x+x^3 ∈ Z[x]
    let g = polynomial![2, 3, 5]; // 2+3x+5x^2 ∈ Z[x]
    let r = f.resultant(g);
    assert_eq!(r, 164);
    ```
    */
    pub fn resultant(mut self, other: Self) -> R
    where
        R: Sized
            + Clone
            + Zero
            + One
            + for<'x> AddAssign<&'x R>
            + for<'x> MulAssign<&'x R>
            + Neg<Output = R>,
        for<'x> &'x R: Sub<Output = R> + Mul<Output = R> + Div<Output = R>,
    {
        let f_deg = self.deg();
        let g_deg = other.deg();
        match (f_deg, g_deg) {
            (Some(0), Some(0)) => R::one(),
            (Some(0), None) => R::one(),
            (None, Some(0)) => R::one(),
            (None, None) => R::zero(),
            (None, Some(_)) => R::zero(),
            (Some(_), None) => R::zero(),
            (Some(0), Some(m)) => ring_algorithm::power::<R>(self.lc().unwrap().clone(), m as u64),
            (Some(n), Some(0)) => ring_algorithm::power::<R>(other.lc().unwrap().clone(), n as u64),
            (Some(n), Some(m)) => {
                debug_assert!(n >= 1);
                debug_assert!(m >= 1);
                let (scale, _) = self.pseudo_division(&other);
                if let Some(l) = self.deg() {
                    let sign = if n * m % 2 == 0 { R::one() } else { -R::one() };
                    let mul =
                        ring_algorithm::power::<R>(other.lc().unwrap().clone(), (n - l) as u64);
                    let div = ring_algorithm::power::<R>(scale, m as u64);
                    &(other.resultant(self) * sign * mul) / &div
                } else {
                    // g | f, gcd(f, g) = g
                    R::zero()
                }
            }
        }
    }

    /** Calculate the primitive part of the input polynomial,
    that is divide the polynomial by the GCD of its coefficents.
    ```
    use polynomial_ring::{polynomial, Polynomial};
    use num_traits::{One, Zero};
    let mut f = polynomial![2, 4, -2, 6]; // 2+4x+2x^2+6x^3 ∈ Z[x]
    f.primitive_part_mut();
    assert_eq!(f, polynomial![1, 2, -1, 3]);// 1+2x+x^2+3x^3 ∈ Z[x]
    let mut g = polynomial![polynomial![1, 1], polynomial![1, 2, 1], polynomial![3, 4, 1], polynomial![-1, -1]]; // (1+x)+(1+2x+x^2)y+(3+4x+x^2)y^2+(-1-x)y^3 ∈ Z[x][y]
    g.primitive_part_mut();
    assert_eq!(g, polynomial![polynomial![1], polynomial![1, 1], polynomial![3, 1], polynomial![-1]]); // 1+(1+x)y+(3+x)y^2-y^3 ∈ Z[x][y]
    ```
    */
    pub fn primitive_part_mut(&mut self)
    where
        R: Sized + Clone + Zero + for<'x> DivAssign<&'x R> + RingNormalize,
        for<'x> &'x R: Rem<Output = R>,
    {
        if self.deg().is_none() {
            return;
        }
        let mut g = self.coef[0].clone();
        for c in &self.coef[1..] {
            g = gcd::<R>(g, c.clone());
        }
        g.normalize_mut();
        *self /= &g;
    }
}

// field
impl<K: Sized> Polynomial<K> {
    /** Make the polynomial monic, that is divide it by its leading coefficient.

    ```
    use num::Rational64;
    use polynomial_ring::Polynomial;
    let mut p = Polynomial::new(vec![1, 2, 3].into_iter().map(|x| Rational64::from_integer(x)).collect());
    p.monic();
    let q = Polynomial::new(vec![(1, 3), (2, 3), (1, 1)].into_iter().map(|(n, d)| Rational64::new(n, d)).collect());
    assert_eq!(p, q);
    ```
    */
    pub fn monic(&mut self)
    where
        K: Clone + for<'x> DivAssign<&'x K>,
    {
        if let Some(lc) = self.lc() {
            let lc = lc.clone();
            *self /= &lc;
        }
    }

    /** Polynomial division.

    ```
    use num::Rational64;
    use polynomial_ring::Polynomial;
    let f = Polynomial::new(vec![3, 1, 4, 1, 5].into_iter().map(|x| Rational64::from_integer(x)).collect());
    let g = Polynomial::new(vec![2, 7, 1].into_iter().map(|x| Rational64::from_integer(x)).collect());
    let mut r = f.clone();
    let q = r.division(&g);
    assert_eq!(f, q * g + r);
    ```
    */
    #[allow(unknown_lints, clippy::return_self_not_must_use)]
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

    /** Calculate the [square-free decomposition](https://en.wikipedia.org/wiki/Square-free_polynomial).

    ```
    use polynomial_ring::{Polynomial, polynomial};
    use num::Rational64;
    let f = polynomial![Rational64::from(1), Rational64::from(1)];
    let g = polynomial![Rational64::from(1), Rational64::from(1), Rational64::from(1)];
    let p = &f * &f * &f * &g * &g; // (x+1)^3(x^2+x+1)^2
    assert_eq!(p.square_free(), &f * &g); // (x+1)(x^2+x+1)
    ```
    */
    #[must_use]
    pub fn square_free(&self) -> Self
    where
        K: Sized
            + Clone
            + Zero
            + One
            + for<'x> AddAssign<&'x K>
            + for<'x> SubAssign<&'x K>
            + for<'x> DivAssign<&'x K>,
        for<'x> &'x K: Mul<Output = K> + Div<Output = K>,
    {
        let d = self.derivative().into_normalize();
        let f = ring_algorithm::gcd::<Self>(self.clone(), d).into_normalize();
        (self / &f).into_normalize()
    }
}
