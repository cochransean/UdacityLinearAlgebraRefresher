from decimal import Decimal, getcontext
from copy import deepcopy

from vector import Vector
from plane import Plane

getcontext().prec = 30


class LinearSystem(object):

    ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG = 'All planes in the system should live in the same dimension'
    NO_SOLUTIONS_MSG = 'No solutions'
    INF_SOLUTIONS_MSG = 'Infinitely many solutions'

    def __init__(self, planes):
        try:
            d = planes[0].dimension
            for p in planes:
                assert p.dimension == d

            self.planes = planes
            self.dimension = d

        except AssertionError:
            raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)

    def swap_rows(self, row1, row2):
        self[row2], self[row1] = self[row1], self[row2]

    def multiply_coefficient_and_row(self, coefficient, row):
        n = self[row].normal_vector
        k = self[row].constant_term

        new_normal_vector = n.scalar_mult(coefficient)
        new_constant_term = k * coefficient

        self[row] = Plane(normal_vector=new_normal_vector, constant_term=new_constant_term)

    def add_multiple_times_row_to_row(self, coefficient, row_to_add, row_to_be_added_to):
        n1 = self[row_to_add].normal_vector
        n2 = self[row_to_be_added_to].normal_vector
        k1 = self[row_to_add].constant_term
        k2 = self[row_to_be_added_to].constant_term

        new_normal_vector = n1.scalar_mult(coefficient) + n2
        new_constant_term = k1 * coefficient + k2

        self[row_to_be_added_to] = Plane(normal_vector=new_normal_vector, constant_term=new_constant_term)

    def indices_of_first_nonzero_terms_in_each_row(self):
        num_equations = len(self)
        indices = [-1] * num_equations

        for i, p in enumerate(self.planes):
            try:
                indices[i] = p.first_nonzero_index(p.normal_vector)
            except Exception as e:
                if str(e) == Plane.NO_NONZERO_ELTS_FOUND_MSG:
                    continue
                else:
                    raise e

        return indices

    def compute_triangular_form(self):
        system = deepcopy(self)
        num_equations = len(system)
        num_variables = self.dimension
        current_variable = 0

        for row in range(0, num_equations):
            while current_variable < num_variables:
                d = MyDecimal(system[row].normal_vector[current_variable])
                if d.is_near_zero():
                    swap_succeeded = system.swap_with_row_below_for_nonzero_coefficient(row, current_variable)
                    if not swap_succeeded:
                        current_variable += 1
                        continue

                system.clear_coefficients_above_or_below(row, current_variable)

                current_variable += 1
                break

        return system

    def swap_with_row_below_for_nonzero_coefficient(self, start_row, variable):
        for current_row in range(start_row + 1, len(self)):
            d = MyDecimal(self[current_row].normal_vector[variable])
            if not d.is_near_zero():
                self.swap_rows(start_row, current_row)
                return True

        # If non-zero is not found
        return False

    def clear_coefficients_above_or_below(self, start_row, variable, below=True):
        n1 = self[start_row].normal_vector
        num_equations = len(self)

        if below:
            other_rows = range(start_row + 1, num_equations)
        else:
            other_rows = reversed(range(0, start_row))

        for other_row in other_rows:
            n2 = self[other_row].normal_vector
            coefficient = -(n2[variable] / n1[variable])
            self.add_multiple_times_row_to_row(coefficient, start_row, other_row)

    def scale_row_to_make_coefficient_equal_one(self, row, col):
        n = self[row].normal_vector
        current_coefficient = n[col]
        beta = Decimal('1') / current_coefficient
        self.multiply_coefficient_and_row(beta, row)

    def compute_rref(self):
        tf = self.compute_triangular_form()
        num_equations = len(tf)
        pivot_indicies = tf.indices_of_first_nonzero_terms_in_each_row()

        for row in reversed(range(0, num_equations)):
            current_variable = pivot_indicies[row]
            if pivot_indicies[row] < 0:
                continue
            tf.scale_row_to_make_coefficient_equal_one(row, current_variable)
            tf.clear_coefficients_above_or_below(row, pivot_indicies[row], below=False)

        return tf

    def solve(self):
        tf = self.compute_rref()
        tf.raise_exception_if_contradictory_equation()
        tf.raise_exception_if_too_few_pivots()

        # Otherwise, there must be a unique solution
        num_variables = self.dimension
        return Vector([tf[x].constant_term for x in range(num_variables)])

    def raise_exception_if_contradictory_equation(self):
        for p in self.planes:
            try:
                p.first_nonzero_index(p.normal_vector)

            except Exception as e:
                if str(e) == 'No nonzero elements found':
                    constant_term = MyDecimal(p.constant_term)
                    if not constant_term.is_near_zero():
                        raise Exception(self.NO_SOLUTIONS_MSG)

                else:
                    raise e

    def raise_exception_if_too_few_pivots(self):
        pivot_indicies = self.indices_of_first_nonzero_terms_in_each_row()
        num_pivots = len([x for x in pivot_indicies if x >= 0])

        if num_pivots < self.dimension:
            raise Exception(self.INF_SOLUTIONS_MSG)

    def __len__(self):
        return len(self.planes)

    def __getitem__(self, i):
        return self.planes[i]

    def __setitem__(self, i, x):
        try:
            assert x.dimension == self.dimension
            self.planes[i] = x

        except AssertionError:
            raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)

    def __str__(self):
        ret = 'Linear System:\n'
        temp = ['Equation {}: {}'.format(i+1, p) for i, p in enumerate(self.planes)]
        ret += '\n'.join(temp)
        return ret


class MyDecimal(Decimal):
    def is_near_zero(self, eps=1e-10):
        return abs(self) < eps

if __name__ == '__main__':

    p1 = Plane(normal_vector=Vector(['5.862', '1.178', '-10.366']), constant_term='-8.15')
    p2 = Plane(normal_vector=Vector(['-2.931', '-0.589', '5.183']), constant_term='-4.075')
    s = LinearSystem([p1, p2])
    try:
        print(s.solve())
    except Exception as e:
        print(str(e))

    p1 = Plane(normal_vector=Vector(['8.631', '5.112', '-1.816']), constant_term='-5.113')
    p2 = Plane(normal_vector=Vector(['4.315', '11.132', '-5.27']), constant_term='-6.775')
    p3 = Plane(normal_vector=Vector(['-2.158', '3.01', '-1.727']), constant_term='-0.831')
    s = LinearSystem([p1, p2, p3])
    try:
        print(s.solve())
    except Exception as e:
        print(str(e))

    p1 = Plane(normal_vector=Vector(['5.262', '2.739', '-9.878']), constant_term='-3.441')
    p2 = Plane(normal_vector=Vector(['5.111', '6.358', '7.638']), constant_term='-2.152')
    p3 = Plane(normal_vector=Vector(['2.016', '-9.924', '-1.367']), constant_term='-9.278')
    p4 = Plane(normal_vector=Vector(['2.167', '-13.543', '-18.883']), constant_term='-10.567')
    s = LinearSystem([p1, p2, p3, p4])
    try:
        print(s.solve())
    except Exception as e:
        print(str(e))
