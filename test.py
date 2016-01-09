from vector import Vector as v
from line import Line as l
from plane import Plane as p
from linsys import LinearSystem as linSys

# Test for Lines
def checkLine(line1, line2):
	if line1.is_equal(line2):
		print "They are equal"
	elif line1.is_parallel(line2):
		print "They are parallel"
	else:
		print line1.intersection(line2)

a = l(v([4.046, 2.836]), 1.21)
b = l(v([10.115, 7.09]), 3.025)
checkLine(a, b)

c = l(v([7.204, 3.182]), 8.68)
d = l(v([8.172, 4.114]), 9.883)
checkLine(c, d)

e = l(v([1.182, 5.562]), 6.744)
f = l(v([1.773, 8.343]), 9.525)
checkLine(e, f)

# Test for ps
def checkp(p1, p2):
	if p1.is_equal(p2):
		print "They are equal"
	elif p1.is_parallel(p2):
		print "They are parallel"
	else:
		print "They are not parallel"

g = p(v([-0.412, 3.806, 0.728]), -3.46)
h = p(v([1.03, -9.515, -1.82]), 8.65)
checkp(g, h)

i = p(v([2.611, 5.528, 0.283]), 4.6)
j = p(v([7.715, 8.306, 5.342]), 3.76)
checkp(i, j)

k = p(v([-7.926, 8.625, -7.212]), -7.952)
m = p(v([-2.642, 2.875, -2.404]), -2.443)
checkp(k, m)

# Test for Linear System
p0 = p(v([1, 1, 1]), 1)
p1 = p(v([0, 1, 0]), 2)
p2 = p(v([1, 1, -1]), 3)
p3 = p(v([1, 0, -2]), 2)

s = linSys([p0,p1,p2,p3])
s.swap_rows(0,1)
if not (s[0] == p1 and s[1] == p0 and s[2] == p2 and s[3] == p3):
    print 'test case 1 failed'

s.swap_rows(1,3)
if not (s[0] == p1 and s[1] == p3 and s[2] == p2 and s[3] == p0):
    print 'test case 2 failed'

s.swap_rows(3,1)
if not (s[0] == p1 and s[1] == p0 and s[2] == p2 and s[3] == p3):
    print 'test case 3 failed'

s.multiply_coefficient_and_row(1,0)
if not (s[0] == p1 and s[1] == p0 and s[2] == p2 and s[3] == p3):
    print 'test case 4 failed'

s.multiply_coefficient_and_row(-1,2)
if not (s[0] == p1 and
        s[1] == p0 and
        s[2] == p(v([-1, -1, 1]), -3) and
        s[3] == p3):
    print 'test case 5 failed'

s.multiply_coefficient_and_row(10,1)
if not (s[0] == p1 and
        s[1] == p(v([10, 10, 10]), 10) and
        s[2] == p(v([-1, -1, 1]), -3) and
        s[3] == p3):
    print 'test case 6 failed'

s.add_multiple_times_row_to_row(0,0,1)
if not (s[0] == p1 and
        s[1] == p(v([10, 10, 10]), 10) and
        s[2] == p(v([-1, -1, 1]), -3) and
        s[3] == p3):
    print 'test case 7 failed'

s.add_multiple_times_row_to_row(1,0,1)
if not (s[0] == p1 and
        s[1] == p(v([10, 11, 10]), 12) and
        s[2] == p(v([-1, -1, 1]), -3) and
        s[3] == p3):
    print 'test case 8 failed'

s.add_multiple_times_row_to_row(-1,1,0)
if not (s[0] == p(v([-10, -10, -10]), -10) and
        s[1] == p(v([10, 11, 10]), 12) and
        s[2] == p(v([-1, -1, 1]), -3) and
        s[3] == p3):
    print 'test case 9 failed'

# Test for Triangular form
p1 = p(v([1, 1, 1]), 1)
p2 = p(v([0, 1, 1]), 2)
s = linSys([p1,p2])
t = s.compute_triangular_form()
if not (t[0] == p1 and
        t[1] == p2):
    print 'test case 1 failed'

p1 = p(v([1, 1, 1]), 1)
p2 = p(v([1, 1, 1]), 2)
s = linSys([p1,p2])
t = s.compute_triangular_form()
if not (t[0] == p1 and
        t[1] == p(constant_term=1)):
    print 'test case 2 failed'

p1 = p(v([1, 1, 1]), 1)
p2 = p(v([0, 1, 0]), 2)
p3 = p(v([1, 1, -1]), 3)
p4 = p(v([1, 0, -2]), 2)
s = linSys([p1,p2,p3,p4])
t = s.compute_triangular_form()
if not (t[0] == p1 and
        t[1] == p2 and
        t[2] == p(v([0, 0, -2]), 2) and
        t[3] == p()):
    print 'test case 3 failed'

p1 = p(v([0, 1, 1]), 1)
p2 = p(v([1, -1, 1]), 2)
p3 = p(v([1, 2, -5]), 3)
s = linSys([p1,p2,p3])
t = s.compute_triangular_form()
if not (t[0] == p(v([1, -1, 1]), 2) and
        t[1] == p(v([0, 1, 1]), 1) and
        t[2] == p(v([0, 0, -9]), -2)):
    print 'test case 4 failed'

# Test for RREF
p1 = p(v([1, 1, 1]), 1)
p2 = p(v([0, 1, 1]), 2)
s = linSys([p1,p2])
r = s.compute_rref()
if not (r[0] == p(v([1, 0, 0]), -1) and
        r[1] == p2):
    print 'test case 1 failed'

p1 = p(v([1, 1, 1]), 1)
p2 = p(v([1, 1, 1]), 2)
s = linSys([p1,p2])
r = s.compute_rref()
if not (r[0] == p1 and
        r[1] == p(constant_term=1)):
    print 'test case 2 failed'

p1 = p(v([1, 1, 1]), 1)
p2 = p(v([0, 1, 0]), 2)
p3 = p(v([1, 1, -1]), 3)
p4 = p(v([1, 0, -2]), 2)
s = linSys([p1,p2,p3,p4])
r = s.compute_rref()
if not (r[0] == p(v([1, 0, 0]), 0) and
        r[1] == p2 and
        r[2] == p(v([0, 0, 1]), -1) and
        r[3] == p()):
    print 'test case 3 failed'

p1 = p(v([0, 1, 1]), 1)
p2 = p(v([1, -1, 1]), 2)
p3 = p(v([1, 2, -5]), 3)
s = linSys([p1,p2,p3])
r = s.compute_rref()
if not (r[0] == p(v([1, 0, 0]), 23.0 / 9) and
        r[1] == p(v([0, 1, 0]), 7.0 / 9) and
        r[2] == p(v([0, 0, 1]), 2.0 / 9)):
    print 'test case 4 failed'

p1 = p(v([5.862, 1.178, -10.366]), -8.15)
p2 = p(v([-2.931, -0.589, 5.183]), -4.075)
s = linSys([p1,p2])
s.compute_solution()

p1 = p(v([8.631, 5.112, -1.816]), -5.113)
p2 = p(v([4.315, 11.132, -5.27]), -6.775)
p3 = p(v([-2.158, 3.01, -1.727]), -0.831)
s = linSys([p1,p2,p3])
s.compute_solution()

p1 = p(v([5.262, 2.739, -9.878]), -3.441)
p2 = p(v([5.111, 6.358, 7.638]), -2.152)
p3 = p(v([2.016, -9.924, -1.367]), -9.278)
p4 = p(v([2.167, -13.543, -18.833]), -10.567)
s = linSys([p1,p2,p3,p4])
s.compute_solution()

p1 = p(v([0.786, 0.786, 0.588]), -0.714)
p2 = p(v([-0.138, -0.138, 0.244]), 0.319)
s = linSys([p1,p2])
print s.compute_rref()

p1 = p(v([8.631, 5.112, -1.816]), -5.113)
p2 = p(v([4.315, 11.132, -5.27]), -6.775)
p3 = p(v([-2.158, 3.01, -1.727]), -0.831)
s = linSys([p1,p2,p3])
print s.compute_rref()

p1 = p(v([0.935, 1.76, -9.365]), -9.955)
p2 = p(v([0.187, 0.352, -1.873]), -1.991)
p3 = p(v([0.374, 0.704, -3.746]), -3.982)
p4 = p(v([-0.561, -1.056, 5.619]), 5.973)
s = linSys([p1,p2,p3,p4])
print s.compute_rref()