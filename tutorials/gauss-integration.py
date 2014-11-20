from sympy import *

x = Symbol('x')
y = Symbol('y')

exactResult = integrate(x * sin(y), (x, 0, 10), (y, 0, 10))

print "exact: {} = {}".format(exactResult, N(exactResult))
