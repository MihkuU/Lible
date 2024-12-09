from utils import *

def writeKernelInstantiate(lbra, lket):

	lalb_list = []
	for la in range(0, lbra + 1):
		for lb in range(0, la + 1):	
			if la + lb == lbra:
				lalb_list.append((la, lb))

	lcld_list = []
	for lc in range(0, lket + 1):
		for ld in range(0, lc + 1):
			if lc + ld == lket:
				lcld_list.append((lc, ld))

	with open('eri_kernels_{}_{}.cpp'.format(lbra, lket), 'w') as file:
		file_str = ''
		file_str += '#include <lible/ints/twoel/eri_kernels.hpp>\n\n'

		# Write ERI4 kernels
		for la, lb in lalb_list:
			for lc, ld in lcld_list:
				file_str += 'template void lible::ints::two::eri4Kernel<{}, {}, {}, {}>(const int, const int, const int, const int,\n'.format(la, lb, lc, ld)
				file_str += '                                                       const double*, const double*,\n'
				file_str += '                                                       const double*, const double*,\n'
				file_str += '                                                       const double*, const double*,\n'
				file_str += '                                                       const double*, const double*,\n'
				file_str += '                                                       const double*, const double*,\n'
				file_str += '                                                       double*);\n\n'
		
		# Write ERI3 kernels
		for la, lb in lalb_list:
			file_str += 'template void lible::ints::two::eri3Kernel<{}, {}, {}>(const int, const int, const int,\n'.format(la, lb, lket)
			file_str += '                                                    const double*, const double*, const double*,\n'
			file_str += '                                                    const double*, const double*, const double*,\n'
			file_str += '                                                    const double*, const double*, double*);\n\n'

		# Write ERI2 kernels
		file_str += 'template void lible::ints::two::eri2Kernel<{}, {}>(const int, const int,\n'.format(lbra, lket)
		file_str += '                                                 const double*, const double*,\n'
		file_str += '                                                 const double*, const double*,\n'
		file_str += '                                                 const double*, const double*,\n'
		file_str += '                                                 double*);\n\n'

		file.write(file_str)

def writeKernelGenerate(lbra, lket):

	lalb_list = []
	for la in range(0, lbra + 1):
		for lb in range(0, la + 1):
			if la + lb == lbra:
				lalb_list.append((la, lb))

	lcld_list = []
	for lc in range(0, lket + 1):
		for ld in range(0, lc + 1):
			if lc + ld == lket:
				lcld_list.append((lc, ld))

	file_str = ''
	file_str += '#include <lible/ints/twoel/eri_kernels.hpp>\n\n'

	# Write ERI4 kernels
	for la, lb in lalb_list:
		for lc, ld in lcld_list:
			file_str += 'template<> void lible::ints::two::eri4Kernel<{}, {}, {}, {}>(const int cdepth_a, const int cdepth_b,\n'.format(la, lb, lc, ld)
			file_str += '                                                         const int cdepth_c, const int cdepth_d,\n'
			file_str += '                                                         const double *exps_a, const double *exps_b,\n'
			file_str += '                                                         const double *exps_c, const double *exps_d,\n'
			file_str += '                                                         const double *coords_a, const double *coords_b,\n'
			file_str += '                                                         const double *coords_c, const double *coords_d,\n'
			file_str += '                                                         const double *ecoeffs_ab, const double *ecoeffs_cd_tsp,\n'
			file_str += '                                                         double *eri4_batch)\n'

			file_str += '{\n'
			file_str += '    constexpr int la = {}, lb = {}, lc = {}, ld = {};\n'.format(la, lb, lc, ld)
			file_str += '    constexpr int lab = la + lb;\n'
			file_str += '    constexpr int lcd = lc + ld;\n'
			file_str += '    constexpr int labcd = lab + lcd;\n\n'

			file_str += '    constexpr int n_sph_a = numSphericalsC(la);\n'
			file_str += '    constexpr int n_sph_b = numSphericalsC(lb);\n'
			file_str += '    constexpr int n_sph_c = numSphericalsC(lc);\n'
			file_str += '    constexpr int n_sph_d = numSphericalsC(ld);\n'
			file_str += '    constexpr int n_hermite_ab = numHermitesC(lab);\n'
			file_str += '    constexpr int n_hermite_cd = numHermitesC(lcd);\n'
			file_str += '    constexpr int n_sph_ab = n_sph_a * n_sph_b;\n'
			file_str += '    constexpr int n_sph_cd = n_sph_c * n_sph_d;\n'
			file_str += '    constexpr int n_ecoeffs_ab = n_sph_ab * n_hermite_ab;\n\n'
			file_str += '    constexpr int n_ecoeffs_cd = n_sph_cd * n_hermite_cd;\n\n'

			file_str += '    std::fill(eri4_batch, eri4_batch + n_sph_ab * n_sph_cd, 0);\n\n'

			file_str += '    std::array<double, labcd + 1> fnx;\n'
			file_str += '    BoysF2<labcd> boys_f;\n\n'

			file_str += '    constexpr int n_hermites_abcd = numHermitesC(lab) * numHermitesC(lcd);\n'
			file_str += '    std::array<double, n_hermites_abcd> rints;\n\n'

			file_str += '    constexpr int n_rints_x_ecoeffs = n_sph_cd * n_hermite_ab;\n'
			file_str += '    std::vector<double> rints_x_ecoeffs(cdepth_a * cdepth_b * n_rints_x_ecoeffs, 0);\n\n'

			file_str += '    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)\n'
			file_str += '        for (int ib = 0; ib < cdepth_b; ib++, iab++)\n'
			file_str += '        {\n'
			file_str += '            int pos_rints_x_ecoeffs = iab * n_rints_x_ecoeffs;\n'
			file_str += '            for (int ic = 0, icd = 0; ic < cdepth_c; ic++)\n'
			file_str += '                for (int id = 0; id < cdepth_d; id++, icd++)\n'
			file_str += '                {\n'
			file_str += '                    double a = exps_a[ia];\n'
			file_str += '                    double b = exps_b[ib];\n'
			file_str += '                    double c = exps_c[ic];\n'
			file_str += '                    double d = exps_d[id];\n'

			file_str += '                    double p = a + b;\n'
			file_str += '                    double q = c + d;\n'

			file_str += '                    std::array<double, 3> xyz_p{(a * coords_a[0] + b * coords_b[0]) / p,\n'
			file_str += '                                                (a * coords_a[1] + b * coords_b[1]) / p,\n'
			file_str += '                                                (a * coords_a[2] + b * coords_b[2]) / p};\n\n'

			file_str += '                    std::array<double, 3> xyz_q{(c * coords_c[0] + d * coords_d[0]) / q,\n'
			file_str += '                                                (c * coords_c[1] + d * coords_d[1]) / q,\n'
			file_str += '                                                (c * coords_c[2] + d * coords_d[2]) / q};\n\n'

			file_str += '                    std::array<double, 3> xyz_pq{xyz_p[0] - xyz_q[0], xyz_p[1] - xyz_q[1],\n'
			file_str += '                                                 xyz_p[2] - xyz_q[2]};\n\n'

			file_str += '                    double alpha = p * q / (p + q);\n'
			file_str += '                    double dx{xyz_pq[0]}, dy{xyz_pq[1]}, dz{xyz_pq[2]};\n'
			file_str += '                    double x = alpha * (dx * dx + dy * dy + dz * dz);\n'
			file_str += '                    boys_f.calcFnx(x, &fnx[0]);\n\n'

			file_str += '                    double fac = (2.0 * std::pow(M_PI, 2.5) / (p * q * std::sqrt(p + q)));\n'
			file_str += '                    calcRInts<lab, lcd>(alpha, fac, &fnx[0], &xyz_pq[0], &rints[0]);\n\n'

			# First rollout
			file_str += rolloutERI4First(la, lb, lc, ld)
			file_str += '                }\n'
			file_str += '        }\n\n'

			file_str += '    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)\n'
			file_str += '        for (int ib = 0; ib < cdepth_b; ib++, iab++)\n'
			file_str += '        {\n'

			# Second rollout
			file_str += rolloutERI4Second(la, lb, lc, ld)
			file_str += '        }\n'
			file_str += '}\n\n'
			 
	# print(file_str)

		# Write ERI3 kernels
	for la, lb in lalb_list:
		file_str += 'template void lible::ints::two::eri3Kernel<{}, {}, {}>(const int, const int, const int,\n'.format(la, lb, lket)
		file_str += '                                                    const double*, const double*, const double*,\n'
		file_str += '                                                    const double*, const double*, const double*,\n'
		file_str += '                                                    const double*, const double*, double*);\n\n'

	# Write ERI2 kernels
	file_str += 'template void lible::ints::two::eri2Kernel<{}, {}>(const int, const int,\n'.format(lbra, lket)
	file_str += '                                                 const double*, const double*,\n'
	file_str += '                                                 const double*, const double*,\n'
	file_str += '                                                 const double*, const double*,\n'
	file_str += '                                                 double*);\n\n'	

	with open('eri_kernels_{}_{}.cpp'.format(lbra, lket), 'w') as file:
		file.write(file_str)
	# Write ERI3 kernels

	# Write ERI2 kernels


hermiteIdxs(4)


l_max = 12
l_max_generate = 4 
for lbra in range(0, l_max + 1):
	for lket in range(0, l_max + 1): 
		l_sum = lbra + lket
		writeKernelInstantiate(lbra, lket)
		if l_sum > l_max_generate:
			writeKernelInstantiate(lbra, lket)
		else:
			writeKernelGenerate(lbra, lket)
