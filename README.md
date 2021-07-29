# Polynomial Ring
Polynomial ring on Field (or Ring).

```rust
use num::Rational64;
use polynomial_ring::Polynomial;
let f = Polynomial::new(vec![3, 1, 4, 1, 5].into_iter().map(|x| Rational64::from_integer(x)).collect());
let g = Polynomial::new(vec![2, 7, 1].into_iter().map(|x| Rational64::from_integer(x)).collect());
let mut r = f.clone();
let q = r.division(&g);
assert_eq!(f, q * g + r);
```

`Add`, `Sub`, `Mul`, `Div`, and `Rem` of polynomial is implemented.
Derivative, square free, pseudo division, and resultant is implemented.

# Licence
MIT OR Apache-2.0
