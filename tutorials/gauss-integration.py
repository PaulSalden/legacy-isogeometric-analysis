from sympy import *

x = Symbol('x')
y = Symbol('y')

exactResult = integrate(x * sin(y), (x, 0, 10), (y, 0, 10))

class BSpline(object):
    def __init__(self, knots, cpoints):
        self.knots = knots
        self.cpoints = cpoints
        # assume an open knot vector
        self.p = knots.count(knots[0]) - 1

    def basisfunc(self, i, p, eta):
        if p == 0:
            if eta >= self.knots[i - 1] and eta < self.knots[i]:
                return 1
            return 0

        An = eta - self.knots[i - 1]
        Ad = self.knots[i + self.p - 1] - self.knots[i - 1]
        # account for repeated knots (I think?)
        A = An / Ad if Ad != 0 else 0

        Bn = self.knots[i + self.p] - eta
        Bd = self.knots[i + self.p] - self.knots[i]
        B = Bn / Bd if Bd != 0 else 0

        return A * self.basisfunc(i, self.p - 1, eta) \
            + B * self.basisfunc(i + 1, self.p - 1, eta)

    def eval(self, eta):
        return sum([self.cpoints[i - 1]
            * self.basisfunc(i, self.p, eta)
            for i in range(len(self.cpoints))])

class NURBS(BSpline):
    def __init__(self, knots, cpoints, weights):
        super(NURBS, self).__init__(knots, cpoints)
        self.weights = weights

    def weightfunc(self, eta):
        return sum([self.weights[i - 1]
            * super(NURBS, self).basisfunc(i, self.p, eta)
            for i in range(len(self.weights))])

    def basisfunct(self, i, p, eta):
        return self.weights[i - 1] \
            * super(NURBS, self).basisfunc(i, self.p, eta) \
            / self.weightfunc(eta)

# integrate a function defined as a linear combination of a
# NURBS basis
def integrate(cvars, NURBS):
    pass

print "exact: {} = {}".format(exactResult, N(exactResult))