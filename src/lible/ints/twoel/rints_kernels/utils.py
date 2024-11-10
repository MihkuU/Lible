
def idxCart(i, j, k):
	jk = j + k
	return int(jk * (jk + 1) / 2 + k)

def idxRR(offset, n, t, u, v):
	tuv = t + u + v
	if tuv == 0:
		return int(n)
	else:
		return int(offset + idxCart(t, u, v))

def numHermites(l):
	return int((l + 1) * (l + 2) * (l + 3) / 6)

def hermiteIdxs(l):
	n_hermites = numHermites(l)

	hermite_idxs = []
	for n in range(0, l + 1):
		for i in range(n, -1, -1):
			for j in range(n - i, -1, -1):
				k = n - i - j
				hermite_idxs.append((i, j, k))

	return hermite_idxs

def idxRRollout(lab, t, u, v):
	tuv = t + u + v
	if tuv == 0:
		return 0

	idx_cart = idxCart(t, u, v)
	offset = lab + numHermites(tuv - 1)

	return offset + idx_cart
