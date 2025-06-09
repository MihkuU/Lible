import itertools, math
import numpy as np

def numHermites(l):
	return int((l + 1) * (l + 2) * (l + 3) / 6)

def hermiteIdxPoss(l):

	hermite_idxs_poss = dict()

	pos = 0
	for n in range(0, l + 1):
		for i in range(n, -1, -1):
			for j in range(n - i, -1, -1):
				k = n - i - j
				hermite_idxs_poss[(i, j, k)] = pos
				pos += 1				

	return hermite_idxs_poss

def vm(m):
	if m >= 0:
		return 0.0
	else:
		return 0.5

def non0CartIdxs(l, m):

	t_bound = math.floor((l - abs(m)) / 2)
	v_max = math.floor(abs(m) / 2 - vm(m)) + vm(m)

	v_list = list(np.arange(vm(m), v_max + 1, 1))

	cart_idxs = set()
	for t in range(0, t_bound + 1):
		for u in range(0, t + 1):			
			for v in v_list:
				lx = int(2 * t + abs(m) - 2 * (u + v))
				ly = int(2 * (u + v))
				lz = int(l - 2 * t - abs(m))

				cart_idxs.add((lx, ly, lz))

	return cart_idxs	

def non0HermiteIdxs(cart_idxs):

	hermite_idxs = set()
	for (i, j, k) in cart_idxs:
		for t in range(0, i + 1):
			for u in range(0, j + 1):
				for v in range(0, k + 1):
					hermite_idxs.add((t, u, v))

	return hermite_idxs	

def rolloutSHARKKet1(lbra, lket):

	n_sph_ket = 2 * lket + 1
	n_hermites_bra = numHermites(lbra)
	n_hermites_ket = numHermites(lket)

	m_list_ket = [0] + list(itertools.chain.from_iterable([[m, -m] for m in range(1, lket + 1)]))
	hermite_idxposs_ket = hermiteIdxPoss(lket)
	
	file_str = 'template<> void lible::ints::shark_mm_ket1<{}, {}>(const double *R, const double *ET, double *R_x_ET)\n'.format(lbra, lket)
	file_str += '{\n'

	for tuv_bra in range(0, n_hermites_bra):
		for nu in range(0, len(m_list_ket)):

			mket = m_list_ket[nu]
			cart_idxs_ket = non0CartIdxs(lket, mket)
			hermite_idxs_ket = non0HermiteIdxs(cart_idxs_ket)

			for tuv_ket_triplet in hermite_idxs_ket:
				tuv_ket = hermite_idxposs_ket[tuv_ket_triplet]

				idx_R_x_E = tuv_bra * n_sph_ket + nu
				idx_R = tuv_bra * n_hermites_ket + tuv_ket
				idx_E = tuv_ket * n_sph_ket + nu

				file_str += '    R_x_ET[{}] += R[{}] * ET[{}];\n'.format(idx_R_x_E, idx_R, idx_E)

	file_str += '}\n\n'

	return file_str

def rolloutSHARKKet2(lbra, lc, ld):
	
	n_sph_c = 2 * lc + 1
	n_sph_d = 2 * ld + 1
	n_sph_cd = n_sph_c * n_sph_d
	n_hermites_bra = numHermites(lbra)
	n_hermites_cd = numHermites(lc + ld)
	
	m_list_c = [0] + list(itertools.chain.from_iterable([[m, -m] for m in range(1, lc + 1)]))
	m_list_d = [0] + list(itertools.chain.from_iterable([[m, -m] for m in range(1, ld + 1)]))

	hermite_idxposs_cd = hermiteIdxPoss(lc + ld)

	file_str = 'template<> void lible::ints::shark_mm_ket2<{}, {}, {}>(const double *R, const double *ET, double *R_x_ET)\n'.format(lbra, lc, ld)
	file_str += '{\n'

	for tuv_bra in range(0, n_hermites_bra):
		for ka in range(0, len(m_list_c)):
			for ta in range(0, len(m_list_d)):

				kata = ka * n_sph_d + ta
				mc = m_list_c[ka]
				md = m_list_d[ta]

				cart_idxs_c = non0CartIdxs(lc, mc)
				cart_idxs_d = non0CartIdxs(ld, md)

				cart_idxs_cd = set()

				for [i, j, k] in cart_idxs_c:
					for [i_, j_, k_] in cart_idxs_d:
						ii_ = i + i_
						jj_ = j + j_
						kk_ = k + k_

						cart_idxs_cd.add((ii_, jj_, kk_))

				hermite_idxs_cd = non0HermiteIdxs(cart_idxs_cd)

				for tuv_cd_triplet in hermite_idxs_cd:
					tuv_cd = hermite_idxposs_cd[tuv_cd_triplet]

					idx_R_x_E = tuv_bra * n_sph_cd + kata
					idx_R = tuv_bra * n_hermites_cd + tuv_cd
					idx_E = tuv_cd * n_sph_cd + kata

					file_str += '    R_x_ET[{}] += R[{}] * ET[{}];\n'.format(idx_R_x_E, idx_R, idx_E)


	file_str += '}\n\n'

	return file_str

def rolloutSHARKBra1(lbra, lket):

	n_sph_bra = 2 * lbra + 1
	n_sph_ket = 2 * lket + 1

	n_hermites_bra = numHermites(lbra)
	m_list_bra = [0] + list(itertools.chain.from_iterable([[m, -m] for m in range(1, lbra + 1)]))
	hermite_idxposs_bra = hermiteIdxPoss(lbra)

	file_str = 'template<> void lible::ints::shark_mm_bra1<{}, {}>(const double *E, const double *R_x_ET, double *E_x_R_x_ET)\n'.format(lbra, lket)
	file_str += '{\n'	

	for mu in range(0, n_sph_bra):
		for nu in range(0, n_sph_ket):

			munu = mu * n_sph_ket + nu
			mbra = m_list_bra[mu]
			cart_idxs_bra = non0CartIdxs(lbra, mbra)
			hermite_idxs_bra = non0HermiteIdxs(cart_idxs_bra)

			for tuv_bra_triplet in hermite_idxs_bra:
				tuv_bra = hermite_idxposs_bra[tuv_bra_triplet]

				idx_R_x_ET = tuv_bra * n_sph_ket + nu
				idx_E = mu * n_hermites_bra + tuv_bra

				file_str += '    E_x_R_x_ET[{}] += E[{}] * R_x_ET[{}];\n'.format(munu, idx_E, idx_R_x_ET)

	file_str += '}\n\n'
	
	return file_str

def rolloutSHARKBra2_ERI3(la, lb, lket):

	n_sph_a = 2 * la + 1
	n_sph_b = 2 * lb + 1
	n_sph_ket = 2 * lket + 1

	n_hermites_ab = numHermites(la + lb)
	m_list_a = [0] + list(itertools.chain.from_iterable([[m, -m] for m in range(1, la + 1)]))
	m_list_b = [0] + list(itertools.chain.from_iterable([[m, -m] for m in range(1, lb + 1)]))
	hermite_idxposs_ab = hermiteIdxPoss(la + lb)

	file_str = 'template<> void lible::ints::shark_mm_bra2<{}, {}, {}>(const double *E, const double *R_x_ET, double *E_x_R_x_ET)\n'.format(la, lb, lket)
	file_str += '{\n'

	for mu in range(0, n_sph_a):
		for nu in range(0, n_sph_b):
			munu = mu * n_sph_b + nu
			for ka in range(0, n_sph_ket):

				munuka = munu * n_sph_ket + ka
				ma = m_list_a[mu]
				mb = m_list_b[nu]

				cart_idxs_a = non0CartIdxs(la, ma)
				cart_idxs_b = non0CartIdxs(lb, mb)

				cart_idxs_ab = set()

				for [i, j, k] in cart_idxs_a:
					for [i_, j_, k_] in cart_idxs_b:
						ii_ = i + i_
						jj_ = j + j_
						kk_ = k + k_

						cart_idxs_ab.add((ii_, jj_, kk_))

				hermite_idxs_ab = non0HermiteIdxs(cart_idxs_ab)

				for tuv_ab_triplet in hermite_idxs_ab:
					tuv_ab = hermite_idxposs_ab[tuv_ab_triplet]

					idx_R_x_ET = tuv_ab * n_sph_ket + ka
					idx_E = munu * n_hermites_ab + tuv_ab

					file_str += '    E_x_R_x_ET[{}] += E[{}] * R_x_ET[{}];\n'.format(munuka, idx_E, idx_R_x_ET)

	file_str += '}\n\n'

	return file_str

def rolloutSHARKBra2_ERI4(la, lb, lc, ld):

	n_sph_a = 2 * la + 1
	n_sph_b = 2 * lb + 1
	n_sph_c = 2 * lc + 1
	n_sph_d = 2 * ld + 1
	n_sph_cd = n_sph_c * n_sph_d

	n_hermites_ab = numHermites(la + lb)
	m_list_a = [0] + list(itertools.chain.from_iterable([[m, -m] for m in range(1, la + 1)]))
	m_list_b = [0] + list(itertools.chain.from_iterable([[m, -m] for m in range(1, lb + 1)]))
	hermite_idxposs_ab = hermiteIdxPoss(la + lb)

	file_str = 'template<> void lible::ints::shark_mm_bra2<{}, {}, {}, {}>(const double *E, const double *R_x_ET, double *E_x_R_x_ET)\n'.format(la, lb, lc, ld)
	file_str += '{\n'

	for mu in range(0, n_sph_a):
		for nu in range(0, n_sph_b):
			munu = mu * n_sph_b + nu
			for ka in range(0, n_sph_c):
				for ta in range(0, n_sph_d):
					kata = ka * n_sph_d + ta
					munukata = munu * (n_sph_c * n_sph_d) + kata

					ma = m_list_a[mu]
					mb = m_list_b[nu]

					cart_idxs_a = non0CartIdxs(la, ma)
					cart_idxs_b = non0CartIdxs(lb, mb)

					cart_idxs_ab = set()

					for [i, j, k] in cart_idxs_a:
						for [i_, j_, k_] in cart_idxs_b:
							ii_ = i + i_
							jj_ = j + j_
							kk_ = k + k_

							cart_idxs_ab.add((ii_, jj_, kk_))

					hermite_idxs_ab = non0HermiteIdxs(cart_idxs_ab)

					for tuv_ab_triplet in hermite_idxs_ab:
						tuv_ab = hermite_idxposs_ab[tuv_ab_triplet]

						idx_R_x_ET = tuv_ab * n_sph_cd + kata
						idx_E = munu * n_hermites_ab + tuv_ab

						file_str += '    E_x_R_x_ET[{}] += E[{}] * R_x_ET[{}];\n'.format(munukata, idx_E, idx_R_x_ET)

	file_str += '}\n\n'

	return file_str
