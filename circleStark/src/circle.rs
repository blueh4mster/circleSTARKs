use lambdaworks_math::elliptic_curve::edwards::curves::bandersnatch::field;
use lambdaworks_math::field::element::FieldElement;
use lambdaworks_math::field::fields::mersenne31::field::Mersenne31Field as M31;

#[derive(Clone)]
pub struct CirclePoint {
    x: FieldElement<M31>,
    y: FieldElement<M31>,
}
pub trait CircleImpl {
    fn double(&self) -> CirclePoint;
    fn zeroes(&self, shape: usize) -> Vec<CirclePoint>;
}

impl CircleImpl for CirclePoint {
    // (x,y) ->  (2x^2-1 , 2*x*y)
    fn double(&self) -> CirclePoint {
        return CirclePoint {
            x: self.x.square().double() - FieldElement::one(),
            y: (self.y * self.x).double(),
        };
    }
    fn zeroes(&self, shape: usize) -> Vec<CirclePoint> {
        vec![
            CirclePoint {
                x: FieldElement::new(0),
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
