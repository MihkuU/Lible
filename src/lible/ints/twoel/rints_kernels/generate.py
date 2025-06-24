import utils

def writeKernelInstantiate(la, lb):

	file_str = ''

	file_str += '#include <lible/ints/rints_meta.hpp>\n\n'

	file_str += instantiateR_generic(la, lb)

	file_str += instantiateR_ERI2D1(la, lb)

	file_str += instantiateR_ERI2D2(la, lb)

	file_str += instantiateR_ERI3D1(la, lb)

	file_str += instantiateR_ERI4D1(la, lb)

	file_str += 'template void lible::ints::calcRInts_ERI<{}, {}>(const double, const double, const double*, const double*, double*);\n\n'.format(la, lb)

	# ERI2-deriv R-ints kernel instantiationo
	file_str += 'template void lible::ints::calcRInts_ERI2_deriv1<{}, {}>(const double, const double, const double*, const double*, double*);\n'.format(la, lb)	

	with open('rints_kernel_{}_{}.cpp'.format(la, lb), 'w') as file:
		 file.write(file_str)

def instantiateR_generic(la, lb):

	file_str = 'template void lible::ints::calcRInts_ERI_new<{}, {}>(const double alpha, const double fac, const double *fnx,\n'.format(la, lb)
	file_str += '                                                   const double *xyz_pq, const int n_cols, const int ofs_row, \n'
	file_str += '                                                   const int ofs_col, double *rints_out);\n\n'
	
	return file_str

def instantiateR_ERI2D1(la, lb):

	file_str = 'template void lible::ints::calcRInts_ERI2D1<{}, {}>(const double alpha, const double fac, const double *fnx,\n'.format(la, lb)
	file_str += '                                                  const double *xyz_ab, double *rints);\n\n'

	return file_str

def instantiateR_ERI2D2(la, lb):

	file_str = 'template void lible::ints::calcRInts_ERI2D2<{}, {}>(const double alpha, const double fac, const double *fnx,\n'.format(la, lb)
	file_str += '                                                  const double *xyz_ab, double *rints);\n\n'

	return file_str

def instantiateR_ERI3D1(la, lb):

	file_str = 'template void lible::ints::calcRInts_ERI3D1<{}, {}>(const double alpha, const double fac, const double *fnx,\n'.format(la, lb)
	file_str += '                                                  const double *xyz_pc, double *rints);\n\n'

	return file_str

def instantiateR_ERI4D1(la, lb):
	
	file_str = 'template void lible::ints::calcRInts_ERI4D1<{}, {}>(const double alpha, const double fac, const double *fnx,\n'.format(la, lb)
	file_str += '                                                  const double *xyz_pq, const int n_rints, const int ofs_row,\n'
	file_str += '                                                  const int ofs_col, const int n_cols, const int n_rows,\n'
	file_str += '                                                  double *rints_out);\n\n'

	return file_str

def specializeR0(la, lb):

	lab = la + lb

	file_str = 'template<>\n'
	file_str += 'void lible::ints::calcRInts_ERI_new<{}, {}>(const double alpha, const double fac, const double *fnx,\n'.format(la, lb)
	file_str += '                                          const double *xyz_ab, const int n_cols, const int ofs_row, \n'
	file_str += '                                          const int ofs_col, double *rints_out)\n'
	file_str += '{\n'

	file_str += '    constexpr int lab = {};\n'.format(la + lb)
	file_str += '    constexpr int buff_size = {};\n'.format(utils.numHermites(la + lb) + la + lb)
	file_str += '    std::array<double, buff_size> rints_buff{};\n\n'

	file_str += '    rints_buff[0] = fnx[0];\n'
	file_str += '    double x = -2 * alpha;\n'
	file_str += '    double y = x;\n'
	file_str += '    for (int n = 1; n <= lab; n++)\n'
	file_str += '    {\n'
	file_str += '        rints_buff[n] = fnx[n] * y;\n'
	file_str += '        y *= x;\n'
	file_str += '    }\n\n'

	# R-ints recursion
	file_str += '    // R-ints recursion\n'
	for n in range(lab - 1, -1, -1):
		n_ = lab - n
		for m in range(n_, 0, -1):
			offset_lhs = int(lab + utils.numHermites(m - 1))
			offset_rhs1 = int(lab + utils.numHermites(m - 2))
			offset_rhs2 = int(lab + utils.numHermites(m - 3))
			for t in range(m, -1, -1):		
				for u in range(m - t, -1, -1):
					v = m - t - u
					idx_lhs = int(offset_lhs + utils.idxCart(t, u, v))
					if t > 0:
						idx_rhs1 = utils.idxRR(offset_rhs1, n + 1, t - 1, u, v)
						file_str += '    rints_buff[{}] = xyz_ab[0] * rints_buff[{}]'.format(idx_lhs, idx_rhs1)
						if t > 1:
							idx_rhs2 = utils.idxRR(offset_rhs2, n + 1, t - 2, u, v)
							file_str += ' + {} * rints_buff[{}]'.format(t - 1, idx_rhs2)
						file_str += ';\n'
					else:
						if u > 0:
							idx_rhs1 = utils.idxRR(offset_rhs1, n + 1, t, u - 1, v)
							file_str += '    rints_buff[{}] = xyz_ab[1] * rints_buff[{}]'.format(idx_lhs, idx_rhs1)
							if u > 1:
								idx_rhs2 = utils.idxRR(offset_rhs2, n + 1, t, u - 2, v)
								file_str += ' + {} * rints_buff[{}]'.format(u - 1, idx_rhs2)
							file_str += ';\n'
						elif v > 0:
							idx_rhs1 = utils.idxRR(offset_rhs1, n + 1, t, u, v - 1)
							file_str += '    rints_buff[{}] = xyz_ab[2] * rints_buff[{}]'.format(idx_lhs, idx_rhs1)
							if v > 1:
								idx_rhs2 = utils.idxRR(offset_rhs2, n + 1, t, u, v - 2)
								file_str += ' + {} * rints_buff[{}]'.format(v - 1, idx_rhs2)
							file_str += ';\n'

	# R-ints rollout
	file_str += '\n'
	file_str += '    // R-ints rollout\n' 

	idxs_a = utils.hermiteIdxs(la)
	idxs_b = utils.hermiteIdxs(lb)
	for j in range(0, len(idxs_b)):
		t_, u_, v_ = idxs_b[j]

		sign = 1.0
		if (t_ + u_ + v_) % 2 != 0:
			sign = -1.0

		for i in range(0, len(idxs_a)):
			t, u, v = idxs_a[i]

			tt_ = t + t_
			uu_ = u + u_
			vv_ = v + v_

			idx_lhs = '(ofs_row + {}) * n_cols + (ofs_col + {})'.format(i, j)
			idx_rhs = utils.idxRRollout(lab, tt_, uu_, vv_)

			file_str += '    rints_out[{}] = {} * fac * rints_buff[{}];\n'.format(idx_lhs, sign, idx_rhs)
	file_str += '}\n\n'

	return file_str

def writeKernelGenerate(la, lb):
	lab = la + lb

	file_str = ''

	#file_str += '#include <lible/ints/utils.hpp>\n\n'
	file_str += '#include <lible/ints/rints_meta.hpp>\n\n'

	file_str += specializeR0(la, lb)

	file_str += 'namespace lible::ints\n'
	file_str += '{\n'

	file_str += 'template <int la, int lb>\n'
	file_str += 'void calcRInts_ERI(const double alpha, const double fac, const double *fnx, const double *xyz_ab,\n'
	file_str += '                   double *rints_out);\n\n'
	
	file_str += 'template<>\n'
	file_str += 'void calcRInts_ERI<{}, {}>(const double alpha, const double fac, const double *fnx,\n'.format(la, lb)
	file_str += '                         const double *xyz_ab, double *rints_out)\n'
	file_str += '{\n'

	file_str += '    constexpr int lab = {};\n'.format(la + lb)
	file_str += '    constexpr int buff_size = {};\n'.format(utils.numHermites(la + lb) + la + lb)
	file_str += '    std::array<double, buff_size> rints_buff{};\n\n'

	file_str += '    rints_buff[0] = fnx[0];\n'
	file_str += '    double x = -2 * alpha;\n'
	file_str += '    double y = x;\n'
	file_str += '    for (int n = 1; n <= lab; n++)\n'
	file_str += '    {\n'
	file_str += '        rints_buff[n] = fnx[n] * y;\n'
	file_str += '        y *= x;\n'
	file_str += '    }\n'

	file_str += '\n'
 
	recursion_str = '    // R-ints recursion\n'
	
	for n in range(lab - 1, -1, -1):
		n_ = lab - n
		for m in range(n_, 0, -1):
			offset_lhs = int(lab + utils.numHermites(m - 1))
			offset_rhs1 = int(lab + utils.numHermites(m - 2))
			offset_rhs2 = int(lab + utils.numHermites(m - 3))
			for t in range(m, -1, -1):		
				for u in range(m - t, -1, -1):
					v = m - t - u
					idx_lhs = int(offset_lhs + utils.idxCart(t, u, v))
					if t > 0:
						idx_rhs1 = utils.idxRR(offset_rhs1, n + 1, t - 1, u, v)
						recursion_str += '    rints_buff[{}] = xyz_ab[0] * rints_buff[{}]'.format(idx_lhs, idx_rhs1)
						if t > 1:
							idx_rhs2 = utils.idxRR(offset_rhs2, n + 1, t - 2, u, v)
							recursion_str += ' + {} * rints_buff[{}]'.format(t - 1, idx_rhs2)
						recursion_str += ';\n'
					else:
						if u > 0:
							idx_rhs1 = utils.idxRR(offset_rhs1, n + 1, t, u - 1, v)
							recursion_str += '    rints_buff[{}] = xyz_ab[1] * rints_buff[{}]'.format(idx_lhs, idx_rhs1)
							if u > 1:
								idx_rhs2 = utils.idxRR(offset_rhs2, n + 1, t, u - 2, v)
								recursion_str += ' + {} * rints_buff[{}]'.format(u - 1, idx_rhs2)
							recursion_str += ';\n'
						elif v > 0:
							idx_rhs1 = utils.idxRR(offset_rhs1, n + 1, t, u, v - 1)
							recursion_str += '    rints_buff[{}] = xyz_ab[2] * rints_buff[{}]'.format(idx_lhs, idx_rhs1)
							if v > 1:
								idx_rhs2 = utils.idxRR(offset_rhs2, n + 1, t, u, v - 2)
								recursion_str += ' + {} * rints_buff[{}]'.format(v - 1, idx_rhs2)
							recursion_str += ';\n'

	file_str += recursion_str

	file_str += '\n'

	rollout_str = '    // R-ints rollout\n' 

	idxs_a = utils.hermiteIdxs(la)
	idxs_b = utils.hermiteIdxs(lb)
	for j in range(0, len(idxs_b)):
		t_, u_, v_ = idxs_b[j]

		sign = 1.0
		if (t_ + u_ + v_) % 2 != 0:
			sign = -1.0

		for i in range(0, len(idxs_a)):
			t, u, v = idxs_a[i]

			tt_ = t + t_
			uu_ = u + u_
			vv_ = v + v_

			idx_lhs = i * len(idxs_b) + j
			idx_rhs = utils.idxRRollout(lab, tt_, uu_, vv_)

			rollout_str += '    rints_out[{}] = {} * fac * rints_buff[{}];\n'.format(idx_lhs, sign, idx_rhs)

	file_str += rollout_str

	
	file_str += '}\n'	
	file_str += '}\n\n'		

	# ERI2-deriv R-ints kernel instantiation
	file_str += 'template void lible::ints::calcRInts_ERI2_deriv1<{}, {}>(const double, const double, const double*, const double*, double*);\n'.format(la, lb)

	file_str += instantiateR_ERI2D1(la, lb)

	file_str += instantiateR_ERI2D2(la, lb)

	file_str += instantiateR_ERI3D1(la, lb)

	file_str += instantiateR_ERI4D1(la, lb)

	with open('rints_kernel_{}_{}.cpp'.format(la, lb), 'w') as file:
		file.write(file_str)
						
#l_max = 12
l_max = 6
l_max_generate = 6
for la in range(0, l_max + 1):
	for lb in range(0, l_max + 1):
		l_sum = la + lb

		if l_sum > l_max:
			continue

		if l_sum > l_max_generate:	
			writeKernelInstantiate(la, lb)
		else:
			writeKernelGenerate(la, lb)
