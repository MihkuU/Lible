import itertools, math
import numpy as np

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
	
	hermite_idxs = []
	for n in range(0, l + 1):
		for i in range(n, -1, -1):
			for j in range(n - i, -1, -1):
				k = n - i - j
				hermite_idxs.append((i, j, k))				

	return hermite_idxs

def hermiteIdxsPoss(l):

	hermite_idxs_poss = dict()

	pos = 0
	for n in range(0, l + 1):
		for i in range(n, -1, -1):
			for j in range(n - i, -1, -1):
				k = n - i - j
				hermite_idxs_poss[(i, j, k)] = pos
				pos += 1				

	return hermite_idxs_poss

def idxRRollout(lab, t, u, v):
	tuv = t + u + v
	if tuv == 0:
		return 0

	idx_cart = idxCart(t, u, v)
	offset = lab + numHermites(tuv - 1)

	return offset + idx_cart

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

def non0HermiteIdxs(cart_idxs_ab):

	hermite_idxs = set()
	for (i, j, k) in cart_idxs_ab:
		for t in range(0, i + 1):
			for u in range(0, j + 1):
				for v in range(0, k + 1):
					hermite_idxs.add((t, u, v))

	return hermite_idxs

def rolloutERI4First(la, lb, lc, ld):
	
	n_sph_c = 2 * lc + 1
	n_sph_d = 2 * ld + 1
	n_sph_cd = n_sph_c * n_sph_d

	n_hermites_ab = numHermites(la + lb)
	n_hermites_cd = numHermites(lc + ld)

	m_list_c = [0] + list(itertools.chain.from_iterable([[m, -m] for m in range(1, lc + 1)]))
	m_list_d = [0] + list(itertools.chain.from_iterable([[m, -m] for m in range(1, ld + 1)]))

	hermite_idxs_poss_cd = hermiteIdxsPoss(lc + ld)

	rollout_str = ''
	rollout_str += '                    const double* p_ecoeffs_cd_tsp = &ecoeffs_cd_tsp[icd * n_ecoeffs_cd];\n'
	rollout_str += '                    double* p_rints_x_ecoeffs = &rints_x_ecoeffs[pos_rints_x_ecoeffs];\n\n'

	for tuv_ab in range(0, n_hermites_ab):		
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

				for tuv_cd_3 in hermite_idxs_cd:
					tuv_cd = hermite_idxs_poss_cd[tuv_cd_3]

					idx_x = tuv_ab * n_sph_cd + kata
					idx_r = tuv_ab * n_hermites_cd + tuv_cd					
					idx_e = tuv_cd * n_sph_cd + kata

					rollout_str += '                    p_rints_x_ecoeffs[{}] += rints[{}] * p_ecoeffs_cd_tsp[{}];\n'.format(idx_x, idx_r, idx_e)

	return rollout_str

def rolloutERI4Second(la, lb, lc, ld):

	n_sph_a = 2 * la + 1
	n_sph_b = 2 * lb + 1
	n_sph_c = 2 * lc + 1
	n_sph_d = 2 * ld + 1
	n_sph_cd = n_sph_c * n_sph_d

	n_hermites_ab = numHermites(la + lb)

	m_list_a = [0] + list(itertools.chain.from_iterable([[m, -m] for m in range(1, la + 1)]))
	m_list_b = [0] + list(itertools.chain.from_iterable([[m, -m] for m in range(1, lb + 1)]))

	hermite_idxs_poss_ab = hermiteIdxsPoss(la + lb)

	rollout_str = ''
	rollout_str += '            const double* p_ecoeffs_ab = &ecoeffs_ab[iab * n_ecoeffs_ab];\n'
	rollout_str += '            const double* p_rints_x_ecoeffs = &rints_x_ecoeffs[iab * n_rints_x_ecoeffs];\n\n'

	for mu in range(0, n_sph_a):
		for nu in range(0, n_sph_b):
			munu = mu * n_sph_b + nu
			for ka in range(0, n_sph_c):
				for ta in range(0, n_sph_d):
					kata = ka * n_sph_d + ta

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

					for tuv_ab_3 in hermite_idxs_ab:
						tuv_ab = hermite_idxs_poss_ab[tuv_ab_3]

						idx_ints = munu * n_sph_cd + kata 
						idx_x = tuv_ab * n_sph_cd + kata
						idx_e = munu * n_hermites_ab + tuv_ab

						rollout_str += '            eri4_batch[{}] += p_ecoeffs_ab[{}] * p_rints_x_ecoeffs[{}];\n'.format(idx_ints, idx_e, idx_x)

	return rollout_str

def rolloutERI3First(la, lb, lc):

	n_sph_c = 2 * lc + 1

	n_hermites_ab = numHermites(la + lb)
	n_hermites_c = numHermites(lc)

	m_list_c = [0] + list(itertools.chain.from_iterable([[m, -m] for m in range(1, lc + 1)]))	

	hermite_idxs_poss_c = hermiteIdxsPoss(lc)

	rollout_str = ''
	rollout_str += '                    const double* p_ecoeffs_c = &ecoeffs_c[ic * n_ecoeffs_c];\n'
	rollout_str += '                    double* p_rints_x_ecoeffs = &rints_x_ecoeffs[pos_rints_x_ecoeffs];\n\n'

	for tuv_ab in range(0, n_hermites_ab):		
		for ka in range(0, len(m_list_c)):

			mc = m_list_c[ka]

			cart_idxs_c = non0CartIdxs(lc, mc)

			hermite_idxs_c = non0HermiteIdxs(cart_idxs_c)

			for tuv_c_3 in hermite_idxs_c:
				tuv_c = hermite_idxs_poss_c[tuv_c_3]

				idx_x = tuv_ab * n_sph_c + ka
				idx_r = tuv_ab * n_hermites_c + tuv_c
				idx_e = ka * n_hermites_c + tuv_c

				rollout_str += '                    p_rints_x_ecoeffs[{}] += rints[{}] * p_ecoeffs_c[{}];\n'.format(idx_x, idx_r, idx_e)

	return rollout_str


def rolloutERI3Second(la, lb, lc):

	n_sph_a = 2 * la + 1
	n_sph_b = 2 * lb + 1
	n_sph_c = 2 * lc + 1

	n_hermites_ab = numHermites(la + lb)

	m_list_a = [0] + list(itertools.chain.from_iterable([[m, -m] for m in range(1, la + 1)]))
	m_list_b = [0] + list(itertools.chain.from_iterable([[m, -m] for m in range(1, lb + 1)]))

	hermite_idxs_poss_ab = hermiteIdxsPoss(la + lb)

	rollout_str = ''
	rollout_str += '            const double* p_ecoeffs_ab = &ecoeffs_ab[iab * n_ecoeffs_ab];\n'
	rollout_str += '            const double* p_rints_x_ecoeffs = &rints_x_ecoeffs[iab * n_rints_x_ecoeffs];\n\n'

	for mu in range(0, n_sph_a):
		for nu in range(0, n_sph_b):
			munu = mu * n_sph_b + nu
			for ka in range(0, n_sph_c):

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

				for tuv_ab_3 in hermite_idxs_ab:
					tuv_ab = hermite_idxs_poss_ab[tuv_ab_3]

					idx_ints = munu * n_sph_c + ka
					idx_x = tuv_ab * n_sph_c + ka
					idx_e = munu * n_hermites_ab + tuv_ab

					rollout_str += '            eri3_batch[{}] += p_ecoeffs_ab[{}] * p_rints_x_ecoeffs[{}];\n'.format(idx_ints, idx_e, idx_x)

	return rollout_str

def rolloutERI2First(la, lb):
	
	n_sph_b = 2 * lb + 1

	n_hermites_a = numHermites(la)
	n_hermites_b = numHermites(lb)

	m_list_b = [0] + list(itertools.chain.from_iterable([[m, -m] for m in range(1, lb + 1)]))	

	hermite_idxs_poss_b = hermiteIdxsPoss(lb)

	rollout_str = ''
	rollout_str += '                    const double* p_ecoeffs_b = &ecoeffs_b_tsp[ib * n_ecoeffs_b];\n'
	rollout_str += '                    double* p_rints_x_ecoeffs = &rints_x_ecoeffs[pos_rints_x_ecoeffs];\n\n'

	for tuv_a in range(0, n_hermites_a):		
		for nu in range(0, len(m_list_b)):

			mb = m_list_b[nu]

			cart_idxs_b = non0CartIdxs(lb, mb)

			hermite_idxs_b = non0HermiteIdxs(cart_idxs_b)

			for tuv_b_3 in hermite_idxs_b:
				tuv_b = hermite_idxs_poss_b[tuv_b_3]

				idx_x = tuv_a * n_sph_b + nu
				idx_r = tuv_a * n_hermites_b + tuv_b
				idx_e = tuv_b * n_sph_b + nu

				rollout_str += '                    p_rints_x_ecoeffs[{}] += rints[{}] * p_ecoeffs_b[{}];\n'.format(idx_x, idx_r, idx_e)

	return rollout_str

def rolloutERI2Second(la, lb):
	
	n_sph_a = 2 * la + 1
	n_sph_b = 2 * lb + 1	

	n_hermites_a = numHermites(la)

	m_list_a = [0] + list(itertools.chain.from_iterable([[m, -m] for m in range(1, la + 1)]))	

	hermite_idxs_poss_a = hermiteIdxsPoss(la)

	rollout_str = ''
	rollout_str += '            const double* p_ecoeffs_a = &ecoeffs_a[ia * n_ecoeffs_a];\n'
	rollout_str += '            const double* p_rints_x_ecoeffs = &rints_x_ecoeffs[ia * n_rints_x_ecoeffs];\n\n'

	for mu in range(0, n_sph_a):
		for nu in range(0, n_sph_b):

				ma = m_list_a[mu]

				cart_idxs_a = non0CartIdxs(la, ma)						

				hermite_idxs_a = non0HermiteIdxs(cart_idxs_a)

				for tuv_a_3 in hermite_idxs_a:
					tuv_a = hermite_idxs_poss_a[tuv_a_3]

					idx_ints = mu * n_sph_b + nu
					idx_x = tuv_a * n_sph_b + nu
					idx_e = mu * n_hermites_a + tuv_a

					rollout_str += '            eri2_batch[{}] += p_ecoeffs_a[{}] * p_rints_x_ecoeffs[{}];\n'.format(idx_ints, idx_e, idx_x)

	return rollout_str	