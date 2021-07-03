# Polynomial Ring
Polynomial ring on Field (or Ring).

```rust
use num::Rational64;
use polynomial_ring::Polynomial;
let p = Polynomial::new(vec![3, 1, 4, 1, 5].into_iter().map(|x| Rational64::from_integer(x)).collect());
let q = Polynomial::new(vec![2, 7, 1].into_iter().map(|x| Rational64::from_integer(x)).collect());
let mut r = p.clone();
let d = r.division(&q);
assert_eq!(p, d * q + r);
```

# Licence
MIT OR Apache-2.0
