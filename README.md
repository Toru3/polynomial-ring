# Polynomial Ring

A polynomial implementation.

```rust
use num::Rational64;
use polynomial_ring::Polynomial;

let f = Polynomial::new(vec![3, 1, 4, 1, 5].into_iter().map(|x| Rational64::from_integer(x)).collect());
let g = Polynomial::new(vec![2, 7, 1].into_iter().map(|x| Rational64::from_integer(x)).collect());
let mut r = f.clone();
let q = r.division(&g);
assert_eq!(f, q * g + r);
let f = Polynomial::new(vec![3, 1, 4, 1, 5].into_iter().map(|x| rug::Rational::from(x)).collect());
let g = Polynomial::new(vec![2, 7, 1].into_iter().map(|x| rug::Rational::from(x)).collect());
let mut r = f.clone();
let q = r.division(&g);
assert_eq!(f, q * g + r);
```

The `Add`, `Sub`, `Mul`, `Div`, and `Rem` traits are implemented for polynomials.
Polynomials also support computing derivative, square free, pseudo division, and resultant.

# Licence

AGPL-3.0-or-later
