from decimal import Decimal, getcontext
from vector import Vector

# Set the precision
getcontext().prec = 30


class Plane(object):

    NO_NONZERO_ELTS_FOUND_MSG = 'No nonzero elements found'
    NO_INTERSECTION_MSG = 'Lines are parallel and do not intersect'
    INFINITE_INTERSECTION_MSG = 'Lines are equal and have infinite intersections'

    def __init__(self, normal_vector=None, constant_term=None):
        self.dimension = 3

        if not normal_vector:
            all_zeros = ['0'] * self.dimension
            normal_vector = Vector(all_zeros)
        self.normal_vector = normal_vector

        if not constant_term:
            constant_term = Decimal('0')
        self.constant_term = Decimal(constant_term)

        # Initializing here for readability
        self.basepoint = self.set_basepoint()

    def set_basepoint(self):
        try:
            n = self.normal_vector
            c = self.constant_term
            basepoint_coords = ['0'] * self.dimension

            initial_index = Plane.first_nonzero_index(n)
            initial_coefficient = n[initial_index]

            basepoint_coords[initial_index] = c/initial_coefficient
            return Vector(basepoint_coords)

        except Exception as e:
            if str(e) == Plane.NO_NONZERO_ELTS_FOUND_MSG:
                return None
            else:
                raise e

    def __str__(self):

        num_decimal_places = 3

        def write_coefficient(coefficient, is_initial_term=False):
            coefficient = round(coefficient, num_decimal_places)
            if coefficient % 1 == 0:
                coefficient = int(coefficient)

            output = ''

            if coefficient < 0:
                output += '-'
            if coefficient > 0 and not is_initial_term:
                output += '+'

            if not is_initial_term:
                output += ' '

            if abs(coefficient) != 1:
                output += '{}'.format(abs(coefficient))

            return output

        n = self.normal_vector

        try:
            initial_index = Plane.first_nonzero_index(n)
            terms = [write_coefficient(n[i], is_initial_term=(i==initial_index)) + 'x_{}'.format(i+1)
                     for i in range(self.dimension) if round(n[i], num_decimal_places) != 0]
            output = ' '.join(terms)

        except Exception as e:
            if str(e) == self.NO_NONZERO_ELTS_FOUND_MSG:
                output = '0'
            else:
                raise e

        constant = round(self.constant_term, num_decimal_places)
        if constant % 1 == 0:
            constant = int(constant)
        output += ' = {}'.format(constant)

        return output

    def is_parallel_to(self, other):
        return self.normal_vector.is_parallel_to(other.normal_vector)

    def __eq__(self, other):

        if not self.is_parallel_to(other):
            return False

        x0 = self.basepoint
        y0 = other.basepoint

        if x0 is None and y0 is None:
            return self.constant_term == other.constant_term
        elif x0 is None or y0 is None:
            return False

        basepoint_difference = x0 - y0
        n = self.normal_vector

        return basepoint_difference.is_parallel_to(n)

    def intersection_with(self, other):
        try:
            A, B = self.normal_vector.coordinates
            C, D = other.normal_vector.coordinates
            k1, k2 = self.constant_term, other.constant_term

            x_numerator = D * k1 - B * k2
            y_numerator = -C * k1 + A * k2
            one_over_denominator = Decimal('1') / (A * D - B * C)

            return Vector([x_numerator, y_numerator]).scalar_mult(one_over_denominator)

        except ZeroDivisionError:
            if self == other:
                return self
            else:
                return None

    @staticmethod
    def first_nonzero_index(iterable):
        for k, item in enumerate(iterable):
            if not MyDecimal(item).is_near_zero():
                return k
        raise Exception(Plane.NO_NONZERO_ELTS_FOUND_MSG)


class MyDecimal(Decimal):
    def is_near_zero(self, eps=1e-10):
        return abs(self) < eps


if __name__ == '__main__':

    def check_answer(plane1, plane2):
        parallel = plane1.is_parallel_to(plane2)
        print('Plane 1 and Plane 2 are parallel? {}', parallel)
        if parallel:
            equal = plane1 == plane2
            print('Plane 1 and Plane 2 are equal? {}', equal)
        else:
            pass
            # print('Plane 1 and Plane 2 intersect at: {}', plane1.intersection_with(plane2))

    plane_x1 = Plane(normal_vector=Vector([-0.412, 3.806, 0.728]), constant_term=-3.46)
    plane_y1 = Plane(normal_vector=Vector([1.03, -9.515, -1.82]), constant_term=8.65)
    check_answer(plane_x1, plane_y1)

    plane_x1 = Plane(normal_vector=Vector([2.611, 5.528, 0.283]), constant_term=4.6)
    plane_y1 = Plane(normal_vector=Vector([7.715, 8.306, 5.342]), constant_term=3.76)
    check_answer(plane_x1, plane_y1)

    plane_x1 = Plane(normal_vector=Vector([-7.926, 8.625, -7.212]), constant_term=-7.952)
    plane_y1 = Plane(normal_vector=Vector([-2.642, 2.875, -2.404]), constant_term=-2.443)
    check_answer(plane_x1, plane_y1)
