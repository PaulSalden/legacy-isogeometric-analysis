import numpy as np

# bezier extraction algorithm according to Borden et al.
# input: knot vector U, curve degree p
# output: extraction operators
def extract(U, p):
	a = p + 1
	b = a + 1
	C = np.identity(p + 1)
	while b < len(U):
		Cnext = np.identity(p + 1) # initialize the next extraction operator
		i = b

		# count multiplicity of the knot at location b
		while b < len(U) and U[b] == U[b - 1]: b += 1
		mult = b - i + 1

		if mult < p:
			# use (10) to compute the alphas
			numer = U[b - 1] - U[a - 1]
			alphas = [0] * (p - mult)
			for j in reversed(range(mult + 1, p + 1)):
				alphas[j - mult - 1] = numer / (U[a + j - 1] - U[a - 1])
			r = p - mult
			# update the matrix coefficients for r new knots
			for j in range(1, r + 1):
				save = r - j + 1
				s = mult + j
				for k in reversed(range(s + 1, p + 2)):
					alpha = alphas[k - s - 1]
					# the following line corresponds to (9)
					C[:, k - 1] = alpha * C[:, k - 1] + (1.0 - alpha) * C[:, k - 2]
				if b < len(U):
					# update overlapping coefficients of the next operator
					Cnext[save - 1:j + save, save - 1] = C[p - j:p + 1, p]
			# finished with the current operator
			if b < len(U):
				# update indices for the next operator
				a = b
				b += 1
		yield C
		C = Cnext

# small test coinciding with the example
if __name__ == "__main__":
	U = [0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 4.0, 4.0, 4.0]
	p = 3
	for C in extract(U, p): print C