use lambdaworks_math::field::element::FieldElement;
use lambdaworks_math::field::fields::mersenne31::field::{
    Mersenne31Field as M31, MERSENNE_31_PRIME_FIELD_ORDER,
};
use lambdaworks_math::field::traits::IsField;
use std::ops::Add;

pub const MODULUS: u32 = MERSENNE_31_PRIME_FIELD_ORDER;
#[derive(Clone, PartialEq, Copy)]
pub struct CirclePoint {
    x: FieldElement<M31>,
    y: FieldElement<M31>,
}

impl Add for CirclePoint {
    type Output = CirclePoint;
    //(x1,y1)+ (x2,y2) ---> ( x1*x2-y1*y2 , x1*y2+y1*x2 )
    fn add(self, rhs: Self) -> Self::Output {
        let x = self.x.clone() * rhs.x.clone() - self.y.clone() * rhs.y.clone();
        let y = self.x * rhs.y + self.y * rhs.x;
        Self { x, y }
    }
}
pub trait CircleImpl {
    fn new_with_field_elements(x: FieldElement<M31>, y: FieldElement<M31>) -> CirclePoint;
    fn new(x: u32, y: u32) -> CirclePoint;
    fn get_x(&self) -> FieldElement<M31>;
    fn get_y(&self) -> FieldElement<M31>;
    fn zero() -> CirclePoint;
    fn double(&self) -> CirclePoint;
    fn zeroes(shape: usize) -> Vec<CirclePoint>;
    fn inverse_x(&self) -> FieldElement<M31>;
    fn inverse_y(&self) -> FieldElement<M31>;
}

impl CircleImpl for CirclePoint {
    fn new(x: u32, y: u32) -> CirclePoint {
        return CirclePoint {
            x: FieldElement::new(x),
            y: FieldElement::new(y),
        };
    }

    fn new_with_field_elements(x: FieldElement<M31>, y: FieldElement<M31>) -> CirclePoint {
        return CirclePoint { x, y };
    }

    fn get_x(&self) -> FieldElement<M31> {
        return self.x;
    }

    fn get_y(&self) -> FieldElement<M31> {
        return self.y;
    }

    fn zero() -> Self {
        CirclePoint {
            x: FieldElement::new(0),
            y: FieldElement::new(0),
        }
    }
    // (x,y) ->  (2x^2-1 , 2*x*y)
    fn double(&self) -> CirclePoint {
        return CirclePoint {
            x: self.x.square().double() - FieldElement::one(),
            y: (self.y * self.x).double(),
        };
    }

    fn zeroes(shape: usize) -> Vec<CirclePoint> {
        vec![
            CirclePoint {
                x: FieldElement::new(1), // ensure this
                y: FieldElement::new(0)
            };
            shape
        ]
    }

    fn inverse_x(&self) -> FieldElement<M31> {
        return FieldElement::inv(&self.x).unwrap();
    }

    fn inverse_y(&self) -> FieldElement<M31> {
        return FieldElement::inv(&self.y).unwrap();
    }
}

#[allow(non_snake_case)]
pub fn G() -> CirclePoint {
    CirclePoint {
        x: FieldElement::new(1268011823),
        y: FieldElement::new(2),
    }
}

#[allow(non_snake_case)]
pub fn Z() -> CirclePoint {
    CirclePoint {
        x: FieldElement::new(1),
        y: FieldElement::new(0),
    }
}

pub fn inverse<F>(x: FieldElement<F>) -> FieldElement<F>
where
    F: IsField,
{
    return FieldElement::inv(&x).unwrap();
}

pub fn scalar_multiply(c: CirclePoint, n: u32) -> CirclePoint {
    match n {
        0 => CirclePoint::new(1, 0), // confirm if (0,0)
        1 => c,
        _ => {
            let half_result = scalar_multiply(c.clone(), n / 2);
            let doubled = half_result.double();
            if n % 2 == 0 {
                doubled
            } else {
                doubled.add(c)
            }
        }
    }
}

/// n==0 not handled
pub fn scalar_division(c: CirclePoint, n: u32) -> CirclePoint {
    let field_inverse_n: FieldElement<M31> = (FieldElement::new(n)).inv().unwrap();
    return scalar_multiply(c, field_inverse_n.to_raw());
}

// naive implementation for Div , Mul standard operations, rather  implement like `Add` as above
pub fn div(c1: CirclePoint, c2: CirclePoint) -> CirclePoint {
    let new_x = (c1.get_x().inv().unwrap().to_raw()).wrapping_mul(c2.get_x().to_raw());
    let new_y = (c1.get_y().inv().unwrap().to_raw()).wrapping_mul(c2.get_y().to_raw());
    return CirclePoint::new(new_x, new_y);
}

pub fn subtract(c1: CirclePoint, c2: CirclePoint) -> CirclePoint {
    let new_x = c1.get_x().to_raw() - c2.get_x().to_raw() + MODULUS;
    let new_y = c1.get_y().to_raw() - c2.get_y().to_raw() + MODULUS;
    return <CirclePoint as CircleImpl>::new(new_x, new_y);
}
pub fn multiply(c1: CirclePoint, c2: CirclePoint) -> CirclePoint {
    let new_x = (c1.get_x().to_raw()).wrapping_mul(c2.get_x().to_raw());
    let new_y = (c1.get_y().to_raw()).wrapping_mul(c2.get_y().to_raw());
    return CirclePoint::new(new_x, new_y);
}
