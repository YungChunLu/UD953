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

    def compute_triangular_form(self):
        system = deepcopy(self)
        # number of equations
        m = len(system)
        # number of variables
        n = system.planes[0].dimension
        for i in range(m):
            j = i
            while j < n:
                plane = system.planes[i]
                # coefficient of variable j in equation i
                c = plane.normal_vector[j]
                if c == 0:
                    indices = system.indices_of_first_nonzero_terms_in_each_row()
                    row = [index for index, coefficient in enumerate(indices) if index > i and coefficient is j][:1]
                    if row:
                        system.swap_rows(i, row[0])
                    else:
                        j += 1
                        continue
                for index in range(i+1, m):
                    p = system.planes[index]
                    c = system.planes[i].normal_vector[j]
                    coefficient = p.normal_vector[j] / (c * -1.0)
                    system.add_multiple_times_row_to_row(coefficient, i, index)
                break
        return system

    def compute_rref(self):
        tf = self.compute_triangular_form()
        # number of equations
        m = len(tf)
        # number of variables
        n = tf.planes[0].dimension
        for index in range(m-1, -1, -1):
            indices = tf.indices_of_first_nonzero_terms_in_each_row()
            if indices[index] is not -1:
                # 1st nonzero term
                i = indices[index]
                if tf.planes[index].normal_vector[i] is not 1:
                    coefficient = 1.0 / tf.planes[index].normal_vector[i]
                    tf.multiply_coefficient_and_row(coefficient, index)
                if index is not 0: 
                    offset = 0
                    for j in range(i, min([n, m])):
                        denominator = tf.planes[index+offset].normal_vector[j]
                        if not (denominator < 1e-10):
                            coefficient = tf.planes[index-1].normal_vector[j] / (denominator * -1)
                            tf.add_multiple_times_row_to_row(coefficient, index+offset, index-1)
                        offset += 1
        return tf

    def compute_solution(self):
        rref = self.compute_rref()
        indices = rref.indices_of_first_nonzero_terms_in_each_row()
        print rref
        if (-1 in indices):
            for i, j in enumerate(indices):
                if j == -1 and abs(rref.planes[i].constant_term) > 0:
                    print "No solution"
                    break
                elif i == len(indices) -1:
                    print "Infinte solutions"
        else:
            print "The solution is : "
            print rref

    def swap_rows(self, row1, row2):
        temp = self[row2]
        self[row2] = self[row1]
        self[row1] = temp


    def multiply_coefficient_and_row(self, coefficient, row):
        plane = self[row]
        plane.constant_term *= coefficient
        plane.normal_vector = plane.normal_vector * coefficient

    def add_multiple_times_row_to_row(self, coefficient, row_to_add, row_to_be_added_to):
        plane_add = self[row_to_add]
        plane_added = self[row_to_be_added_to]
        plane_added.constant_term += (plane_add.constant_term * coefficient)
        plane_added.normal_vector += (plane_add.normal_vector * coefficient)

    def indices_of_first_nonzero_terms_in_each_row(self):
        num_equations = len(self)
        num_variables = self.dimension

        indices = [-1] * num_equations

        for i,p in enumerate(self.planes):
            try:
                indices[i] = p.first_nonzero_index(p.normal_vector)
            except Exception as e:
                if str(e) == Plane.NO_NONZERO_ELTS_FOUND_MSG:
                    continue
                else:
                    raise e

        return indices


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
        temp = ['Equation {}: {}'.format(i+1,p) for i,p in enumerate(self.planes)]
        ret += '\n'.join(temp)
        return ret


class MyDecimal(Decimal):
    def is_near_zero(self, eps=1e-10):
        return abs(self) < eps