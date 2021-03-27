use polynomial_ring::*;
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
    type R = Polynomial<num::Rational32>;
    let mut v = Vec::new();
    v.push(make_q_x(vec![(1, 2), (1, 3)]));
    v.push(make_q_x(vec![(1, 4), (1, 5)]));
    v.push(make_q_x(vec![(1, 8)]));
    v.push(make_q_x(vec![(1, 16), (1, 15), (1, 9)]));
    v.push(make_q_x(vec![(1, 32)]));
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
    use auto_ops::impl_op_ex;
    #[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
    pub struct ResidueClass6 {
        n: i64,
    }
    impl ResidueClass6 {
        pub fn new(n: i64) -> Self {
            ResidueClass6 { n }
        }
    }
    impl std::fmt::Display for ResidueClass6 {
        fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
            write!(f, "{}", self.n)
        }
    }
    impl_op_ex!(+= |a: &mut ResidueClass6, b: &ResidueClass6| { a.n = (a.n + b.n) % 6; });
    impl_op_ex!(*= |a: &mut ResidueClass6, b: &ResidueClass6| { a.n = (a.n * b.n) % 6; });
    impl_op_ex!(+ |a: &ResidueClass6, b: &ResidueClass6| -> ResidueClass6 {
        let mut c = a.clone();
        c += b;
        c
    });
    impl_op_ex!(*|a: &ResidueClass6, b: &ResidueClass6| -> ResidueClass6 {
        let mut c = a.clone();
        c *= b;
        c
    });
    impl num_traits::Zero for ResidueClass6 {
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
    use residue_class::ResidueClass6;
    let p = polynomial![ResidueClass6::new(2), ResidueClass6::new(3)];
    let q = polynomial![ResidueClass6::new(1), ResidueClass6::new(2)];
    let r = polynomial![ResidueClass6::new(2), ResidueClass6::new(1)];
    assert_eq!(p * q, r);
}
#[test]
fn product() {
    type R = Polynomial<i64>;
    let mut v = Vec::new();
    v.push(polynomial![-2, 1]);
    v.push(polynomial![-1, 1]);
    v.push(polynomial![0, 1]);
    v.push(polynomial![1, 1]);
    v.push(polynomial![2, 1]);
    let p = v.into_iter().product::<R>();
    assert_eq!(p, polynomial![0, 4, 0, -5, 0, 1]);
}
#[test]
fn display() {
    assert_eq!(polynomial![0].to_string(), "0");
    assert_eq!(polynomial![3, 2, 1].to_string(), "x^2+2*x+3");
    assert_eq!(polynomial![0, -2, -1, 3].to_string(), "3*x^3+-1*x^2+-2*x");
}
