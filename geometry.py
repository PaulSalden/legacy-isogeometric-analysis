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

    def eval(self, eta):
        result = 0

        for j in range(1, len(self.cpoints) + 1):
            for i in range(1, len(self.cpoints[0]) + 1):
                result += self.cpoints[j - 1][i - 1] * \
                    self.basis_x(i, eta) * self.basis_y(j, eta)

        return result

# allow a test run
if __name__ == "__main__":
    # with n x m basis functions
    xknots = [0., 0., 0., 1., 2., 3., 3., 3.] # n = 5
    yknots = [0., 0., 1., 2., 3., 3.] # m = 4

    cpoints = []
    for j in range(4):
        cpoints.append([])
        for i in range(5):
            cpoints[j].append(np.array((i, j)))

    surf = BSplineSurface((xknots, yknots), cpoints)
    print surf.eval(1.9)
