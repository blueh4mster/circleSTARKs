use std::ops::Add;

use lambdaworks_math::elliptic_curve::edwards::curves::bandersnatch::field;
use lambdaworks_math::field::element::FieldElement;
use lambdaworks_math::field::fields::mersenne31::field::Mersenne31Field as M31;

#[derive(Clone)]
pub struct CirclePoint {
    x: FieldElement<M31>,
    y: FieldElement<M31>,
}

impl Add for CirclePoint {
    type Output = CirclePoint;
    fn add(self, rhs: Self) -> Self::Output {
        let x = self.x.clone() * rhs.x.clone() - self.y.clone() * rhs.y.clone();
        let y = self.x * rhs.y + self.y * rhs.x;
        Self { x, y }
    }
}
pub trait CircleImpl {
    fn new(x: u32, Y: u32) -> CirclePoint;
    fn zero() -> CirclePoint;
    fn double(&self) -> CirclePoint;
    fn zeroes(shape: usize) -> Vec<CirclePoint>;
}

impl CircleImpl for CirclePoint {
    fn new(x: u32, y: u32) -> CirclePoint {
        return CirclePoint {
            x: FieldElement::new(x),
            y: FieldElement::new(y),
        };
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
}

pub fn G() -> CirclePoint {
    CirclePoint {
        x: FieldElement::new(1268011823),
        y: FieldElement::new(2),
    }
}

pub fn Z() -> CirclePoint {
    CirclePoint {
        x: FieldElement::new(1),
        y: FieldElement::new(0),
    }
}
