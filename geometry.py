import numpy as np

# Cox-De Boor recursion formula, result evaluated at eta
# i = 1, 2, ..., n
def cox_de_boor(knots, i, p, eta):
    if p == 0:
        if eta >= knots[i - 1] and eta < knots[i]:
            return 1
        return 0
     
    Ad = knots[i + p - 1] - knots[i - 1]
    # account for p + 1 repeated knots
    if Ad == 0:
        A = 0
    else:
        An = eta - knots[i - 1]
        A = An / Ad

    Bd = knots[i + p] - knots[i]
    if Bd == 0:
        B = 0
    else:
        Bn = knots[i + p] - eta
        B = Bn / Bd

    return A * cox_de_boor(knots, i, p - 1, eta) \
        + B * cox_de_boor(knots, i + 1, p - 1, eta)

class BSplineSurface(object):
    # kvectors is a list of 2 knot vectors
    # cpoints is a 2D list of control points np.array(x, y)
    def __init__(self, kvectors, cpoints):
        self.kvectors = kvectors
        self.cpoints = cpoints
        # assume an open knot vector
        self.polorders = [v.count(v[0]) - 1 for v in kvectors]

    # basis functions in x direction
    def basis_x(self, i, eta):
        return cox_de_boor(self.kvectors[0], i,
            self.polorders[0], eta)
    # basis function in y direction
    def basis_y(self, j, eta):
        return cox_de_boor(self.kvectors[1], j,
            self.polorders[1], eta)

    # 2D basis function
    def basis(self, i, j, eta):
        return self.basis_x(i, eta) * self.basis_y(j, eta)

    def eval(self, eta):
        result = 0

        for j in range(1, len(self.cpoints) + 1):
            for i in range(1, len(self.cpoints[0]) + 1):
                result += self.cpoints[j - 1][i - 1] * self.basis(i, j, eta)

        return result

class NURBSSurface(BSplineSurface):
    # kvectors is a list of 2 knot vectors
    # cpoints is a 2D list of control points np.array(x, y)
    # weights is a 2D list of weights
    def __init__(self, kvectors, cpoints, weights):
        super(NURBSSurface, self).__init__(kvectors, cpoints)
        self.weights = weights

    # weight function evaluated at eta
    def W(self, eta):
        result = 0

        for j in range(1, len(self.weights) + 1):
            for i in range(1, len(self.weights[0]) + 1):
                result += super(NURBSSurface, self).basis(i, j, eta) * \
                    self.weights[j - 1][i - 1]

        return result

    # 2D basis function
    def basis(self, i, j, eta):
        return super(NURBSSurface, self).basis(i, j, eta) * \
            self.weights[j - 1][i - 1] / self.W(eta)

# allow a test run
if __name__ == "__main__":
    # with n x m basis functions
    xknots = [0., 0., 0., 1., 2., 3., 3., 3.]
    n = 5
    yknots = [0., 0., 1., 2., 3., 3.]
    m = 4

    cpoints = []
    for j in range(m):
        cpoints.append([])
        for i in range(n):
            cpoints[j].append(np.array((i, j)))

    surf = BSplineSurface((xknots, yknots), cpoints)
    print "BSpline at 1.9: {}".format(surf.eval(1.9))

    weights = []
    for j in range(m):
        weights.append([])
        for i in range(n):
            weights[j].append((i*j)**0.5)

    surf2 = NURBSSurface((xknots, yknots), cpoints, weights)
    print "NURBS at 1.9: {}".format(surf2.eval(1.9))
