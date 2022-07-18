use num_traits::{One, Zero};
use polynomial_ring::*;
use ring_algorithm::{extended_euclidian_algorithm, gcd, EuclideanRingOperation, RingNormalize};
use std::fmt::{Debug, Display};

fn make_q_x<T>(v: Vec<(T, T)>) -> Polynomial<num::rational::Ratio<T>>
where
    T: Display + Debug + Clone + num::Integer,
{
    let v = v
        .into_iter()
        .map(|(n, d)| num::rational::Ratio::<T>::new(n, d))
        .collect();
    Polynomial::new(v)
}

#[test]
fn add() {
    let p = polynomial![1, 4, 1, 4, 1, 3, 5, 6];
    let q = polynomial![2, 2, 3, 6, 0, 6, 7, 9];
    assert_eq!(p + q, polynomial![3, 6, 4, 10, 1, 9, 12, 15]);
    let p = polynomial![3, 1, 4];
    let q = polynomial![-3, -1, -4];
    assert_eq!(p + q, polynomial![]);
    let p = polynomial![3, 1, 4];
    let q = polynomial![1, 2, -4];
    assert_eq!(p + q, polynomial![4, 3]);
    let p = polynomial![-3, 1, 4];
    let q = polynomial![3, -1, 1];
    assert_eq!(p + q, polynomial![0, 0, 5]);
}

#[test]
fn add2() {
    let p = make_q_x(vec![(3, 1), (4, 1), (5, 9), (2, 6)]);
    let q = make_q_x(vec![(2, 7), (1, 8), (2, 8)]);
    assert_eq!(p + q, make_q_x(vec![(23, 7), (33, 8), (58, 72), (2, 6)]));
}

#[test]
fn sum() {
    type R = Polynomial<num::Rational64>;
    let v = vec![
        make_q_x(vec![(1, 2), (1, 3)]),
        make_q_x(vec![(1, 4), (1, 5)]),
        make_q_x(vec![(1, 8)]),
        make_q_x(vec![(1, 16), (1, 15), (1, 9)]),
        make_q_x(vec![(1, 32)]),
    ];
    assert_eq!(
        v.into_iter().sum::<R>(),
        make_q_x(vec![(31, 32), (9, 15), (1, 9)])
    );
}

#[test]
fn neg() {
    let p = polynomial![1, 4, 1, 4, 1, 3, 5, 6];
    assert_eq!(-p, polynomial![-1, -4, -1, -4, -1, -3, -5, -6]);
}

#[test]
fn sub() {
    let p = polynomial![1, 4, 1, 4, 1, 3, 5, 6];
    let q = polynomial![2, 2, 3, 6, 0, 6, 7, 9];
    assert_eq!(p - q, polynomial![-1, 2, -2, -2, 1, -3, -2, -3]);
    let p = polynomial![3, 1, 4];
    let q = polynomial![3, 1, 4];
    assert_eq!(p - q, polynomial![]);
    let p = polynomial![3, 1, 4];
    let q = polynomial![1, 2, 4];
    assert_eq!(p - q, polynomial![2, -1]);
    let p = polynomial![3, 1, 4];
    let q = polynomial![3, 1, 1];
    assert_eq!(p - q, polynomial![0, 0, 3]);
}

#[test]
fn sub2() {
    let p = make_q_x(vec![(3, 1), (4, 1), (5, 9)]);
    let q = make_q_x(vec![(2, 7), (1, 8), (2, 8), (1, 1)]);
    assert_eq!(p - q, make_q_x(vec![(19, 7), (31, 8), (22, 72), (-1, 1)]));
}

#[test]
fn mul() {
    let p = polynomial![1, 2, 3];
    let q = polynomial![5, 4, 2];
    assert_eq!(p * q, polynomial![5, 14, 25, 16, 6]);
}

mod residue_class {
    use std::ops::*;
    #[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
    pub struct ResidueClass<const M: i64> {
        n: i64,
    }
    impl<const M: i64> ResidueClass<M> {
        pub fn new(n: i64) -> Self {
            ResidueClass { n }
        }
    }
    impl<const M: i64> std::fmt::Display for ResidueClass<M> {
        fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
            write!(f, "{}", self.n)
        }
    }
    macro_rules! impl_ops {
        ($fn:ident, $trait:ident) => {
            #[auto_impl_ops::auto_ops]
            impl<const M: i64> $trait for &ResidueClass<M> {
                type Output = ResidueClass<M>;
                fn $fn(self, rhs: Self) -> Self::Output {
                    ResidueClass::<M>::new(self.n.$fn(rhs.n).rem(M))
                }
            }
        };
    }
    impl_ops!(add, Add);
    impl_ops!(sub, Sub);
    impl_ops!(mul, Mul);
    impl<const M: i64> num_traits::Zero for ResidueClass<M> {
        fn zero() -> Self {
            Self::new(i64::zero())
        }
        fn is_zero(&self) -> bool {
            self.n.is_zero()
        }
    }
}

#[test]
fn mul2() {
    type R = residue_class::ResidueClass<6>;
    let p = polynomial![R::new(2), R::new(3)];
    let q = polynomial![R::new(1), R::new(2)];
    let r = polynomial![R::new(2), R::new(1)];
    assert_eq!(p * q, r);
}

#[test]
fn product() {
    type R = Polynomial<i64>;
    let v = vec![
        polynomial![-2, 1],
        polynomial![-1, 1],
        polynomial![0, 1],
        polynomial![1, 1],
        polynomial![2, 1],
    ];
    let p = v.into_iter().product::<R>();
    assert_eq!(p, polynomial![0, 4, 0, -5, 0, 1]);
}

#[test]
fn display() {
    assert_eq!(polynomial![0].to_string(), "0");
    assert_eq!(polynomial![3, 2, 1].to_string(), "x^2+2*x+3");
    assert_eq!(polynomial![0, -2, -1, 3].to_string(), "3*x^3+-1*x^2+-2*x");
}

#[test]
fn scalar_mul() {
    let p = polynomial![1, 4, 1, 4, 1, 3, 5, 6];
    assert_eq!(p * 2, polynomial![2, 8, 2, 8, 2, 6, 10, 12]);
}

#[test]
fn scalar_div() {
    let p = make_q_x(vec![(3, 1), (4, 1), (5, 9)]);
    let q = make_q_x(vec![(3, 2), (2, 1), (5, 18)]);
    let two = num::rational::Ratio::new(2, 1);
    assert_eq!(p / two, q);
}

macro_rules! poly {
    ($($x:expr),*) => {
        Polynomial::new(vec![$(num::Rational64::from_integer($x)),*])
    }
}

macro_rules! expand_poly {
    ($([$($x:expr),*]),*) => {
        vec![$(poly![$($x),*]),*].into_iter().product::<Polynomial<num::Rational64>>()
    }
}

fn check_eea<T>(a: T, b: T) -> bool
where
    T: Display + Debug + Zero + One + Clone + Eq + RingNormalize,
    for<'x> &'x T: EuclideanRingOperation<T>,
{
    let g = gcd::<T>(a.clone(), b.clone());
    let (d, x, y) = extended_euclidian_algorithm::<T>(a.clone(), b.clone());
    g.is_similar(&d) && &(&x * &a) + &(&y * &b) == d
}

#[test]
fn test_eea2() {
    type R = Polynomial<num::Rational64>;
    let z = R::zero();
    check_eea::<R>(z.clone(), z.clone());
    let a = expand_poly![[2], [1, 1], [2, 1], [3, 1]];
    let b = expand_poly![[3], [1, 1], [4, 1]];
    let d = expand_poly![[4, 1], [5, 1]];
    assert!(check_eea::<R>(a.clone(), z.clone()));
    assert!(check_eea::<R>(z, a.clone()));
    assert!(check_eea::<R>(a.clone(), b));
    assert!(check_eea::<R>(a, d));
}

#[test]
fn pseudo_division() {
    let f = polynomial![1, -1, -1, 1]; // 1-x-x^2+x^3 ∈ Z[x]
    let g = polynomial![1, 2]; // 1+2x ∈ Z[x]
    let mut r = f.clone();
    let (s, q) = r.pseudo_division(&g);
    assert_eq!(polynomial![s] * f, q * g + r);

    // 1-yx-x^2+yx^3 ∈ Z[y][x]
    let f = polynomial![
        polynomial![1],
        polynomial![0, -1],
        polynomial![-1],
        polynomial![0, 1]
    ];
    // -1+y^2x ∈ Z[y][x]
    let g = polynomial![polynomial![-1], polynomial![0, 0, 1]];
    let mut r = f.clone();
    let (s, q) = r.pseudo_division(&g);
    assert_eq!(f * s, q * g + r);

    // x^3 ∈ Z[y][x]
    let f = polynomial![polynomial![], polynomial![], polynomial![], polynomial![1]];
    // yx ∈ Z[y][x]
    let g = polynomial![polynomial![], polynomial![0, 1]];
    let mut r = f.clone();
    let (s, q) = r.pseudo_division(&g);
    assert_eq!(f * s, q * g + r);
}

#[test]
fn resultant() {
    let f = polynomial![-4, 0, 0, 0, 1]; // -4+x^4 ∈ Z[x]
    let g = polynomial![0, 2, 0, 1]; // 2x+x^3 ∈ Z[x]
    let r = f.resultant(g); // deg(gcd(f, g)) = deg(x^2-2) = 2 ≠ 0
    assert_eq!(r, 0);

    let f = polynomial![polynomial![1], polynomial![0], polynomial![1]]; // 1+x^2 ∈ Z[y][x]
    let g = polynomial![polynomial![1], polynomial![1, 2]]; // 1+(1+2y)x ∈ Z[y][x]
    let r = f.resultant(g);
    assert_eq!(r, polynomial![2, 4, 4]); // 2+4y+4y^2

    let y3 = Polynomial::from_monomial(polynomial![-1], 3); // -y^3 ∈ Z[x][y]
    let y2xy = Polynomial::from_monomial(polynomial![0, -1], 2) + &y3; // -xy^2-y^3 ∈ Z[x][y]
    let x = polynomial![polynomial![0, 1]]; // x ∈ Z[x][y]
    let f = polynomial![y3, polynomial![], x.clone()]; // -y^3+xz^2 ∈ Z[x][y][z]
    let g = polynomial![y2xy, polynomial![], x]; // -xy^2-y^3+xz^2 ∈ Z[x][y][z]
    let r = f.resultant(g);
    assert_eq!(
        r,
        Polynomial::from_monomial(Polynomial::from_monomial(1, 4), 4)
    ); // x^4y^4
}
