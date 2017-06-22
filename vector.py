from decimal import Decimal, getcontext
from math import acos, pi


# Set the precision
getcontext().prec = 30


class Vector:

    CANNOT_NORMALIZE_ZERO_VECTOR_MSG = 'Cannot normalize zero vector'
    NO_UNIQUE_PARALLEL_COMPONENT_MSG = 'There is no unique parallel component.'
    NO_UNIQUE_ORTHOGONAL_COMPONENT_MSG = 'There is no unique orthogonal component.'
    MUST_BE_2D_OR_3D_MSG = 'Cross products are only defined in two or three dimensions.'

    def __init__(self, coordinates):
        try:
            if not coordinates:
                raise ValueError
            self.coordinates = tuple(Decimal(str(x)) for x in coordinates)
            self.dimension = len(coordinates)

        except ValueError:
            raise ValueError('The coordinates must be non-empty')

        except TypeError:
            raise TypeError('The coordinates must be an iterable')

    def __str__(self):
        return 'Vector: {}'.format(self.coordinates)

    def __eq__(self, other):
        return self.coordinates == other.coordinates

    def __add__(self, other):
        return Vector(tuple(x + other.coordinates[i] for i, x in enumerate(self.coordinates)))

    def __sub__(self, other):
        return Vector(tuple(x - other.coordinates[i] for i, x in enumerate(self.coordinates)))

    def __round__(self, n=None):
        return Vector(tuple(round(x, 3) for x in self.coordinates))

    def __iter__(self):
        return iter(self.coordinates)

    def __getitem__(self, key):
        return self.coordinates[key]

    def __setitem__(self, i, value):
        new_coordinates = list(self.coordinates)
        new_coordinates[i] = value
        self.coordinates = new_coordinates

    def __len__(self):
        return self.dimension

    def scalar_mult(self, scalar):
        return Vector(tuple(x * scalar for x in self.coordinates))

    def magnitude(self):
        return sum(x ** 2 for x in self.coordinates) ** Decimal('0.5')

    def normalized(self):
        try:
            magnitude = self.magnitude()
            return self.scalar_mult(1 / magnitude)

        except ZeroDivisionError:
            raise Exception(Vector.CANNOT_NORMALIZE_ZERO_VECTOR_MSG)

    def dot(self, other):
        try:
            pairs = zip(self.coordinates, other.coordinates)
            return sum([x * y for x, y in pairs])
        except AttributeError:
            raise AttributeError('Other must be another vector.')

    def angle_with(self, other, in_degrees=False):
        try:
            u1 = self.normalized()
            u2 = other.normalized()
            angle_in_radians = acos(round(u1.dot(u2), 12))

            if in_degrees:
                degrees_per_radian = 180 / pi
                return angle_in_radians * degrees_per_radian
            else:
                return angle_in_radians

        except Exception as e:
            if str(e) == Vector.CANNOT_NORMALIZE_ZERO_VECTOR_MSG:
                raise Exception('Cannot compute an angle with the zero vector.')

    def is_parallel_to(self, other, tolerance=1e-10):
        return (self.is_zero() or
                other.is_zero() or
                self.angle_with(other) <= tolerance or
                abs(self.angle_with(other)) - pi <= tolerance)

    def is_zero(self, tolerance=1e-10):
        return self.magnitude() < tolerance

    def is_orthogonal_to(self, other, tolerance=1e-10):
        return abs(self.dot(other)) < tolerance

    def component_parallel_to(self, basis):
        try:
            u = basis.normalized()
            weight = self.dot(u)
            return u.scalar_mult(weight)
        except Exception as e:
            if str(e) == Vector.CANNOT_NORMALIZE_ZERO_VECTOR_MSG:
                raise Exception(Vector.NO_UNIQUE_PARALLEL_COMPONENT_MSG)

    def component_orthogonal_to(self, basis):
        try:
            projection = self.component_parallel_to(basis)
            return self - projection
        except Exception as e:
            if str(e) == Vector.NO_UNIQUE_PARALLEL_COMPONENT_MSG:
                raise Exception(Vector.NO_UNIQUE_ORTHOGONAL_COMPONENT_MSG)
            else:
                raise e

    def cross_product(self, other):
        try:
            x1, y1, z1 = self.coordinates
            x2, y2, z2 = other.coordinates
            new_coordinates = [
                y1 * z2 - y2 * z1,
                -(x1 * z2 - x2 * z1),
                x1 * y2 - x2 * y1
            ]
            return Vector(new_coordinates)
        except ValueError as e:
            msg = str(e)
            if msg == 'not enough values to unpack (expected 3, got 2)':
                self_3dim = Vector(self.coordinates + (0,))
                other_3dim = Vector(other.coordinates + (0,))
                return self_3dim.cross_product(other_3dim)
            else:
                raise Exception(Vector.MUST_BE_2D_OR_3D_MSG)

    def area_parallelogram(self, other):
        return self.cross_product(other).magnitude()

    def area_triangle(self, other):
        return self.area_parallelogram(other) * Decimal('0.5')

if __name__ == '__main__':
    pass
