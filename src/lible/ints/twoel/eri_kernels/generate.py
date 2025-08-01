from utils import *

def instantiateERI4(la, lb, lc, ld):

	file_str = 'template lible::vec4d lible::ints::eri4KernelFun<{}, {}, {}, {}>(const int ipair_ab, const int ipair_cd,\n'.format(la, lb, lc, ld)
	file_str += '                                                             const ShellPairData &sp_data_ab,\n'
	file_str += '                                                             const ShellPairData &sp_data_cd,\n'
	file_str += '                                                             const ERI4Kernel *eri4_kernel);\n\n'

	return file_str

def instantiateERI3(la, lb, lc):

	file_str = 'template lible::vec3d lible::ints::eri3KernelFun<{}, {}, {}>(const int ipair_ab, const int ishell_c,\n'.format(la, lb, lc)
	file_str += '                                                          const ShellPairData &sp_data_ab,\n'
	file_str += '                                                          const ShellData &sh_data_c,\n'
	file_str += '                                                          const ERI3Kernel *eri3_kernel);\n\n'	

	return file_str

def instantiateERI2(la, lb):

	file_str = 'template lible::vec2d lible::ints::eri2KernelFun<{}, {}>(const int ishell_a, const int ishell_b,\n'.format(la, lb)
	file_str += '                                                       const ShellData &sh_data_a,\n'
	file_str += '                                                       const ShellData &sh_data_b,\n'
	file_str += '                                                       const ERI2Kernel *eri2_kernel);\n\n'
	
	return file_str 

def instantiateERI2D1(la, lb):

	file_str = 'template std::array<lible::vec2d, 6> lible::ints::eri2d1KernelFun<{}, {}>(const int ishell_a, const int ishell_b,\n'.format(la, lb)
	file_str += '                                                                        const ShellData &sh_data_a,\n'
	file_str += '                                                                        const ShellData &sh_data_b,\n'
	file_str += '                                                                        const ERI2D1Kernel *eri2d1_kernel);\n\n'

	return file_str

def instantiateERI2D2(la, lb):

	file_str  = 'template std::array<std::array<lible::vec2d, 6>, 6> \n'
	file_str += 'lible::ints::eri2d2KernelFun<{}, {}>(const int ishell_a, const int ishell_b,\n'.format(la, lb)
	file_str += '                                   const ShellData &sh_data_a,\n'
	file_str += '                                   const ShellData &sh_data_b,\n'
	file_str += '                                   const ERI2D2Kernel *eri2d2_kernel);\n\n'

	return file_str

def instantiateERI3D1(la, lb, lc):

	file_str = 'template std::array<lible::vec3d, 9> lible::ints::eri3d1KernelFun<{}, {}, {}>(const int ipair_ab, const int ishell_c,\n'.format(la, lb, lc)
	file_str += '                                                                           const ShellPairData &sh_data_ab,\n'
	file_str += '                                                                           const ShellData &sh_data_c,\n'
	file_str += '                                                                           const ERI3D1Kernel *eri3d1_kernel);\n\n'

	return file_str

def instantiateERI3D2(la, lb, lc):

	file_str = 'template lible::arr2d<lible::vec3d, 9, 9>\n'
	file_str += 'lible::ints::eri3d2KernelFun<{}, {}, {}>(const int ipair_ab, const int ishell_c,\n'.format(la, lb, lc)
	file_str += '                                      const ShellPairData &sh_data_ab,\n'
	file_str += '                                      const ShellData &sh_data_c,\n'
	file_str += '                                      const ERI3D2Kernel *eri3d2_kernel);\n\n'

	return file_str

def instantiateERI3SOC(la, lb, lc):

	file_str = 'template std::array<lible::vec3d, 3> lible::ints::eri3socKernelFun<{}, {}, {}>(const int ipair_ab, const int ishell_c,\n'.format(la, lb, lc)
	file_str += '                                                                            const ShellPairData &sh_data_ab,\n'
	file_str += '                                                                            const ShellData &sh_data_c,\n'
	file_str += '                                                                            const ERI3SOCKernel *eri3soc_kernel);\n\n'

	return file_str

def instantiateERI4D1(la, lb, lc, ld):

	file_str = 'template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<{}, {}, {}, {}>(const int ipair_ab, const int ipair_cd,\n'.format(la, lb, lc, ld)
	file_str += '                                                                               const ShellPairData &sh_data_ab,\n'
	file_str += '                                                                               const ShellPairData &sp_data_cd,\n'
	file_str += '                                                                               const ERI4D1Kernel *eri4d1_kernel);\n\n'

	return file_str

def instantiateERI4SOC(la, lb, lc, ld):

	file_str = 'template std::array<lible::vec4d, 3> lible::ints::eri4socKernelFun<{}, {}, {}, {}>(const int ipair_ab, const int ipair_cd,\n'.format(la, lb, lc, ld)
	file_str += '                                                                                const ShellPairData &sh_data_ab,\n'
	file_str += '                                                                                const ShellPairData &sp_data_cd,\n'
	file_str += '                                                                                const ERI4SOCKernel *eri4soc_kernel);\n\n'

	return file_str

def specializeERI4(la, lb, lc, ld):

	file_str = 'template<> lible::vec4d lible::ints::eri4KernelFun<{}, {}, {}, {}>(const int ipair_ab, const int ipair_cd,\n'.format(la, lb, lc, ld)
	file_str += '                                                               const ShellPairData &sp_data_ab,\n'
	file_str += '                                                               const ShellPairData &sp_data_cd,\n'
	file_str += '                                                               const ERI4Kernel *eri4_kernel)\n'
	file_str += '{\n'

	file_str += '    // Compile-time data\n'
	file_str += '    constexpr int la = {};\n'.format(la)
	file_str += '    constexpr int lb = {};\n'.format(lb)
	file_str += '    constexpr int lc = {};\n'.format(lc)
	file_str += '    constexpr int ld = {};\n'.format(ld)
	file_str += '    constexpr int lab = la + lb;\n'
	file_str += '    constexpr int lcd = ld + ld;\n'
	file_str += '    constexpr int labcd = lab + lcd;\n\n'

	file_str += '    constexpr int n_sph_a = numSphericalsC(la);\n'
	file_str += '    constexpr int n_sph_b = numSphericalsC(lb);\n'
	file_str += '    constexpr int n_sph_c = numSphericalsC(lc);\n'
	file_str += '    constexpr int n_sph_d = numSphericalsC(ld);\n'
	file_str += '    constexpr int n_sph_ab = n_sph_a * n_sph_b;\n'
	file_str += '    constexpr int n_sph_cd = n_sph_c * n_sph_d;\n'
	file_str += '    constexpr int n_hermite_ab = numHermitesC(lab);\n'
	file_str += '    constexpr int n_hermite_cd = numHermitesC(lcd);\n\n'

	file_str += '    // Read-in data\n'
	file_str += '    const int cdepth_a = sp_data_ab.cdepths[2 * ipair_ab];\n'
	file_str += '    const int cdepth_b = sp_data_ab.cdepths[2 * ipair_ab + 1];\n'
	file_str += '    const int cdepth_c = sp_data_cd.cdepths[2 * ipair_cd];\n'
	file_str += '    const int cdepth_d = sp_data_cd.cdepths[2 * ipair_cd + 1];\n'
	file_str += '    const int cofs_a = sp_data_ab.coffsets[2 * ipair_ab];\n'
	file_str += '    const int cofs_b = sp_data_ab.coffsets[2 * ipair_ab + 1];\n'
	file_str += '    const int cofs_c = sp_data_cd.coffsets[2 * ipair_cd];\n'
	file_str += '    const int cofs_d = sp_data_cd.coffsets[2 * ipair_cd + 1];\n'
	file_str += '    const double *xyz_a = &sp_data_ab.coords[6 * ipair_ab];\n'
	file_str += '    const double *xyz_b = &sp_data_ab.coords[6 * ipair_ab + 3];\n'
	file_str += '    const double *xyz_c = &sp_data_cd.coords[6 * ipair_cd];\n'
	file_str += '    const double *xyz_d = &sp_data_cd.coords[6 * ipair_cd + 3];\n\n'

	file_str += '    const int ofs_E_bra = sp_data_ab.offsets_ecoeffs[ipair_ab];\n'
	file_str += '    const int ofs_E_ket = sp_data_cd.offsets_ecoeffs[ipair_cd];\n'	
	file_str += '    const double *ecoeffs_bra = &eri4_kernel->ecoeffs_bra[ofs_E_bra];\n'
	file_str += '    const double *ecoeffs_ket = &eri4_kernel->ecoeffs_ket[ofs_E_ket];\n\n'

	file_str += '    // R-integrals\n'
	file_str += '    std::array<double, labcd + 1> fnx;\n'
	file_str += '    BoysF2<labcd> boys_f;\n\n'

	file_str += '    int n_r_rows = (cdepth_a * cdepth_b * n_hermite_ab);\n'
	file_str += '    int n_r_cols = (cdepth_c * cdepth_d * n_hermite_cd);\n'
	file_str += '    int n_r_ints = n_r_rows * n_r_cols;\n'
	file_str += '    int n_R_x_E = cdepth_a * cdepth_b * (n_hermite_ab * n_sph_cd);\n'
	file_str += '    std::vector<double> rints(n_r_ints, 0);\n'	
	file_str += '    std::vector<double> R_x_E(n_R_x_E, 0);\n'

	# First loop
	file_str += '    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)\n'
	file_str += '        for (int ib = 0; ib < cdepth_b; ib++, iab++)\n'
	file_str += '            for (int ic = 0, icd = 0; ic < cdepth_c; ic++)\n'
	file_str += '                for (int id = 0; id < cdepth_d; id++, icd++)\n'
	file_str += '                {\n'
	file_str += '                    double a = sp_data_ab.exps[cofs_a + ia];\n'
	file_str += '                    double b = sp_data_ab.exps[cofs_b + ib];\n'
	file_str += '                    double c = sp_data_cd.exps[cofs_c + ic];\n'
	file_str += '                    double d = sp_data_cd.exps[cofs_d + id];\n\n'

	file_str += '                    double p = a + b;\n'
	file_str += '                    double q = c + d;\n'
	file_str += '                    double alpha = p * q / (p + q);\n\n'

	file_str += '                    std::array<double, 3> xyz_p{(a * xyz_a[0] + b * xyz_b[0]) / p,\n'
	file_str += '                                                (a * xyz_a[1] + b * xyz_b[1]) / p,\n'
	file_str += '                                                (a * xyz_a[2] + b * xyz_b[2]) / p};\n\n'

	file_str += '                    std::array<double, 3> xyz_q{(c * xyz_c[0] + d * xyz_d[0]) / q,\n'
	file_str += '                                                (c * xyz_c[1] + d * xyz_d[1]) / q,\n'
	file_str += '                                                (c * xyz_c[2] + d * xyz_d[2]) / q};\n\n'

	file_str += '                    std::array<double, 3> xyz_pq{xyz_p[0] - xyz_q[0],\n'
	file_str += '                                                 xyz_p[1] - xyz_q[1],\n'
	file_str += '                                                 xyz_p[2] - xyz_q[2]};\n\n'

	file_str += '                    double dx{xyz_pq[0]}, dy{xyz_pq[1]}, dz{xyz_pq[2]};\n'
	file_str += '                    double x = alpha * (dx * dx + dy * dy + dz * dz);\n'
	file_str += '                    boys_f.calcFnx(x, &fnx[0]);\n\n'

	file_str += '                    double fac = (2.0 * std::pow(M_PI, 2.5) / (p * q * std::sqrt(p + q)));\n'
	file_str += '                    int ofs_row = iab * n_hermite_ab;\n'
	file_str += '                    int ofs_col = icd * n_hermite_cd;\n\n'

	file_str += '                    calcRInts_ERI_new<lab, lcd>(alpha, fac, fnx.data(), xyz_pq.data(), \n'
	file_str += '                                                n_r_cols, ofs_row, ofs_col, &rints[0]);\n\n'

	file_str += rolloutERI4First(la, lb, lc, ld)
	file_str += '                }\n\n'

	# Second loop	
	file_str += '    int m = n_sph_ab;\n'
	file_str += '    int n = n_sph_cd;\n'
	file_str += '    int k = cdepth_a * cdepth_b * n_hermite_ab;\n\n'

	# file_str += 'vec4d eri4_batch(Fill(0), n_sph_a, n_sph_b, n_sph_c, n_sph_d);'
	file_str += '    vec4d eri4_batch(Fill(0), n_sph_a, n_sph_b, n_sph_c, n_sph_d);\n'
	file_str += '    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0,\n'
	file_str += '                ecoeffs_bra, k, &R_x_E[0], n, 1.0, &eri4_batch[0], n);\n\n'
	# file_str += '    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)\n'
	# file_str += '        for (int ib = 0; ib < cdepth_b; ib++, iab++)\n'
	# file_str += '        {\n'	
	# # file_str += rolloutERI4Second(la, lb, lc, ld)			
	# file_str += '        }\n\n'		

	# Norms
	file_str += '    int ofs_norm_a = sp_data_ab.offsets_norms[2 * ipair_ab];\n'
	file_str += '    int ofs_norm_b = sp_data_ab.offsets_norms[2 * ipair_ab + 1];\n'
	file_str += '    int ofs_norm_c = sp_data_cd.offsets_norms[2 * ipair_cd];\n'
	file_str += '    int ofs_norm_d = sp_data_cd.offsets_norms[2 * ipair_cd + 1];\n'
	file_str += '    for (int mu = 0; mu < n_sph_a; mu++)\n'
	file_str += '        for (int nu = 0; nu < n_sph_b; nu++)\n'
	file_str += '            for (int ka = 0; ka < n_sph_c; ka++)\n'
	file_str += '                for (int ta = 0; ta < n_sph_d; ta++)\n'
	file_str += '                {\n'
	file_str += '                    double norm_a = sp_data_ab.norms[ofs_norm_a + mu];\n'
	file_str += '                    double norm_b = sp_data_ab.norms[ofs_norm_b + nu];\n'
	file_str += '                    double norm_c = sp_data_cd.norms[ofs_norm_c + ka];\n'
	file_str += '                    double norm_d = sp_data_cd.norms[ofs_norm_d + ta];\n'
	file_str += '                    eri4_batch(mu, nu, ka, ta) *= norm_a * norm_b * norm_c * norm_d;\n'
	file_str += '                }\n\n'
		
	file_str += '    return eri4_batch;\n'
	file_str += '}\n\n'	

	return file_str

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
		file_str += '#include <lible/ints/twoel/eri_kernel_funs.hpp>\n\n'

		# Write ERI4 kernels
		for la, lb in lalb_list:
			for lc, ld in lcld_list:
				
				if lbra + lket <= 6:
					file_str += instantiateERI4(la, lb, lc, ld)
					file_str += instantiateERI4D1(la, lb, lc, ld)
					file_str += instantiateERI4SOC(la, lb, lc, ld)
					if la != lb and lc == ld:
						file_str += instantiateERI4(lb, la, lc, ld)
						file_str += instantiateERI4D1(lb, la, lc, ld)
						file_str += instantiateERI4SOC(lb, la, lc, ld)
					elif lc != ld and la == lb:
						file_str += instantiateERI4(la, lb, ld, lc)
						file_str += instantiateERI4D1(la, lb, ld, lc)
						file_str += instantiateERI4SOC(la, lb, ld, lc)
					elif la != lb and lc != ld:
						file_str += instantiateERI4(lb, la, ld, lc)
						file_str += instantiateERI4(lb, la, lc, ld)
						file_str += instantiateERI4(la, lb, ld, lc)
						file_str += instantiateERI4D1(lb, la, ld, lc)
						file_str += instantiateERI4D1(lb, la, lc, ld)
						file_str += instantiateERI4D1(la, lb, ld, lc)
						file_str += instantiateERI4SOC(lb, la, ld, lc)
						file_str += instantiateERI4SOC(lb, la, lc, ld)
						file_str += instantiateERI4SOC(la, lb, ld, lc)
			

				# file_str += 'template lible::vec4d lible::ints::two::eri4Kernel<{}, {}, {}, {}>(const int ipair_ab, const int ipair_cd,\n'.format(la, lb, lc, ld)
				# file_str += '                                                               const std::vector<double> &ecoeffs_ab,\n'
				# file_str += '                                                               const std::vector<double> &ecoeffs_cd_tsp,\n'
				# file_str += '                                                               const ShellPairData &sp_data_ab,\n'
				# file_str += '                                                               const ShellPairData &sp_data_cd);\n\n'
		
		# Write ERI3 kernels
		for la, lb in lalb_list:
			
			if (lbra + lket) <= 6:
				file_str += instantiateERI3(la, lb, lket)
				file_str += instantiateERI3D1(la, lb, lket)
				file_str += instantiateERI3D2(la, lb, lket)
				file_str += instantiateERI3SOC(la, lb, lket)
				if la != lb:
					file_str += instantiateERI3(lb, la, lket)
					file_str += instantiateERI3D1(lb, la, lket)
					file_str += instantiateERI3D2(lb, la, lket)
					file_str += instantiateERI3SOC(lb, la, lket)

			# file_str += 'template lible::vec3d lible::ints::two::eri3Kernel<{}, {}, {}>(const int ipair_ab, const int ishell_c,\n'.format(la, lb, lket)
			# file_str += '                                                             const std::vector<double> &ecoeffs_ab,\n'
			# file_str += '                                                             const std::vector<double> &ecoeffs_c,\n'
			# file_str += '                                                             const ShellPairData &sp_data_ab,\n'
			# file_str += '                                                             const ShellData &sh_data_c);\n\n'			

		# Write ERI2 kernels
		if (lbra + lket) <= 6:
			file_str += instantiateERI2(lbra, lket)
			file_str += instantiateERI2D1(lbra, lket)
			file_str += instantiateERI2D2(lbra, lket)

		# file_str += 'template lible::vec2d lible::ints::two::eri2Kernel<{}, {}>(const int ishell_a, const int ishell_b,\n'.format(lbra, lket)
		# file_str += '                                                         const std::vector<double> &ecoeffs_a,\n'
		# file_str += '                                                         const std::vector<double> &ecoeffs_b_tsp,\n'
		# file_str += '                                                         const ShellData &sh_data_a,\n'
		# file_str += '                                                         const ShellData &sh_data_b);\n\n'

		# file_str += 'template std::array<lible::vec2d, 6> lible::ints::two::eri2d1Kernel<{}, {}>(const int ishell_a, const int ishell_b,\n'.format(lbra, lket)
		# file_str += '                                                                        const std::vector<double> &ecoeffs_a,\n'
		# file_str += '                                                                        const std::vector<double> &ecoeffs_b_tsp,\n'
		# file_str += '                                                                        const ShellData &sh_data_a,\n'
		# file_str += '                                                                        const ShellData &sh_data_b);\n\n'		

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
	file_str += '#include <lible/ints/twoel/eri_kernel_funs.hpp>\n\n'

	# Write ERI4 kernels
	for la, lb in lalb_list:
		for lc, ld in lcld_list:
			# file_str += specializeERI4(la, lb, lc, ld)
			file_str += instantiateERI4(la, lb, lc, ld)			
			file_str += instantiateERI4D1(la, lb, lc, ld)
			file_str += instantiateERI4SOC(la, lb, lc, ld)
			if la != lb and lc == ld:
				file_str += instantiateERI4(lb, la, lc, ld)			
				file_str += instantiateERI4D1(lb, la, lc, ld)
				file_str += instantiateERI4SOC(lb, la, lc, ld)
			elif lc != ld and la == lb:
				file_str += instantiateERI4(la, lb, ld, lc)			
				file_str += instantiateERI4D1(la, lb, ld, lc)
				file_str += instantiateERI4SOC(la, lb, ld, lc)
			elif la != lb and lc != ld:
				file_str += instantiateERI4(lb, la, ld, lc)
				file_str += instantiateERI4(lb, la, lc, ld)
				file_str += instantiateERI4(la, lb, ld, lc)
				file_str += instantiateERI4D1(lb, la, ld, lc)
				file_str += instantiateERI4D1(lb, la, lc, ld)
				file_str += instantiateERI4D1(la, lb, ld, lc)
				file_str += instantiateERI4SOC(lb, la, ld, lc)
				file_str += instantiateERI4SOC(lb, la, lc, ld)
				file_str += instantiateERI4SOC(la, lb, ld, lc)
			
			# file_str += 'template<> lible::vec4d\n'
			# file_str += 'lible::ints::two::eri4Kernel<{}, {}, {}, {}>(const int ipair_ab, const int ipair_cd,\n'.format(la, lb, lc, ld)
			# file_str += '                                         const std::vector<double> &ecoeffs_ab,\n'
			# file_str += '                                         const std::vector<double> &ecoeffs_cd_tsp,\n'
			# file_str += '                                         const ShellPairData &sp_data_ab,\n'
			# file_str += '                                         const ShellPairData &sp_data_cd)\n'
			# file_str += '{\n'
			
			# file_str += '    const int cdepth_a = sp_data_ab.cdepths[2 * ipair_ab];\n'
			# file_str += '    const int cdepth_b = sp_data_ab.cdepths[2 * ipair_ab + 1];\n'
			# file_str += '    const int cdepth_c = sp_data_cd.cdepths[2 * ipair_cd];\n'
			# file_str += '    const int cdepth_d = sp_data_cd.cdepths[2 * ipair_cd + 1];\n'
			# file_str += '    const int cofs_a = sp_data_ab.coffsets[2 * ipair_ab];\n'
			# file_str += '    const int cofs_b = sp_data_ab.coffsets[2 * ipair_ab + 1];\n'
			# file_str += '    const int cofs_c = sp_data_cd.coffsets[2 * ipair_cd];\n'
			# file_str += '    const int cofs_d = sp_data_cd.coffsets[2 * ipair_cd + 1];\n\n'
			
			# file_str += '    const double *exps_a = &sp_data_ab.exps[cofs_a];\n'
			# file_str += '    const double *exps_b = &sp_data_ab.exps[cofs_b];\n'
			# file_str += '    const double *exps_c = &sp_data_cd.exps[cofs_c];\n'
			# file_str += '    const double *exps_d = &sp_data_cd.exps[cofs_d];\n'
			# file_str += '    const double *coords_a = &sp_data_ab.coords[6 * ipair_ab];\n'
			# file_str += '    const double *coords_b = &sp_data_ab.coords[6 * ipair_ab + 3];\n'
			# file_str += '    const double *coords_c = &sp_data_cd.coords[6 * ipair_cd];\n'
			# file_str += '    const double *coords_d = &sp_data_cd.coords[6 * ipair_cd + 3];\n'		
			# file_str += '    const double *pecoeffs_ab = &ecoeffs_ab[sp_data_ab.offsets_ecoeffs[ipair_ab]];\n'
			# file_str += '    const double *pecoeffs_cd_tsp = &ecoeffs_cd_tsp[sp_data_cd.offsets_ecoeffs[ipair_cd]];\n\n'

			# file_str += '    constexpr int la = {}, lb = {}, lc = {}, ld = {};\n'.format(la, lb, lc, ld)
			# file_str += '    constexpr int lab = la + lb;\n'
			# file_str += '    constexpr int lcd = lc + ld;\n'
			# file_str += '    constexpr int labcd = lab + lcd;\n\n'

			# file_str += '    constexpr int n_sph_a = numSphericalsC(la);\n'
			# file_str += '    constexpr int n_sph_b = numSphericalsC(lb);\n'
			# file_str += '    constexpr int n_sph_c = numSphericalsC(lc);\n'
			# file_str += '    constexpr int n_sph_d = numSphericalsC(ld);\n'
			# file_str += '    constexpr int n_hermite_ab = numHermitesC(lab);\n'
			# file_str += '    constexpr int n_hermite_cd = numHermitesC(lcd);\n'
			# file_str += '    constexpr int n_sph_ab = n_sph_a * n_sph_b;\n'
			# file_str += '    constexpr int n_sph_cd = n_sph_c * n_sph_d;\n'
			# file_str += '    constexpr int n_ecoeffs_ab = n_sph_ab * n_hermite_ab;\n'
			# file_str += '    constexpr int n_ecoeffs_cd = n_sph_cd * n_hermite_cd;\n\n'

			# file_str += '    std::array<double, labcd + 1> fnx;\n'
			# file_str += '    BoysF2<labcd> boys_f;\n\n'

			# file_str += '    constexpr int n_hermites_abcd = numHermitesC(lab) * numHermitesC(lcd);\n'
			# file_str += '    std::array<double, n_hermites_abcd> rints;\n\n'

			# file_str += '    constexpr int n_rints_x_ecoeffs = n_sph_cd * n_hermite_ab;\n'
			# file_str += '    std::vector<double> rints_x_ecoeffs(cdepth_a * cdepth_b * n_rints_x_ecoeffs, 0);\n\n'

			# file_str += '    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)\n'
			# file_str += '        for (int ib = 0; ib < cdepth_b; ib++, iab++)\n'
			# file_str += '        {\n'
			# file_str += '            int pos_rints_x_ecoeffs = iab * n_rints_x_ecoeffs;\n'
			# file_str += '            for (int ic = 0, icd = 0; ic < cdepth_c; ic++)\n'
			# file_str += '                for (int id = 0; id < cdepth_d; id++, icd++)\n'
			# file_str += '                {\n'
			# file_str += '                    double a = exps_a[ia];\n'
			# file_str += '                    double b = exps_b[ib];\n'
			# file_str += '                    double c = exps_c[ic];\n'
			# file_str += '                    double d = exps_d[id];\n'

			# file_str += '                    double p = a + b;\n'
			# file_str += '                    double q = c + d;\n'

			# file_str += '                    std::array<double, 3> xyz_p{(a * coords_a[0] + b * coords_b[0]) / p,\n'
			# file_str += '                                                (a * coords_a[1] + b * coords_b[1]) / p,\n'
			# file_str += '                                                (a * coords_a[2] + b * coords_b[2]) / p};\n\n'

			# file_str += '                    std::array<double, 3> xyz_q{(c * coords_c[0] + d * coords_d[0]) / q,\n'
			# file_str += '                                                (c * coords_c[1] + d * coords_d[1]) / q,\n'
			# file_str += '                                                (c * coords_c[2] + d * coords_d[2]) / q};\n\n'

			# file_str += '                    std::array<double, 3> xyz_pq{xyz_p[0] - xyz_q[0], xyz_p[1] - xyz_q[1],\n'
			# file_str += '                                                 xyz_p[2] - xyz_q[2]};\n\n'

			# file_str += '                    double alpha = p * q / (p + q);\n'
			# file_str += '                    double dx{xyz_pq[0]}, dy{xyz_pq[1]}, dz{xyz_pq[2]};\n'
			# file_str += '                    double x = alpha * (dx * dx + dy * dy + dz * dz);\n'
			# file_str += '                    boys_f.calcFnx(x, &fnx[0]);\n\n'

			# file_str += '                    double fac = (2.0 * std::pow(M_PI, 2.5) / (p * q * std::sqrt(p + q)));\n'
			# file_str += '                    calcRInts_ERI<lab, lcd>(alpha, fac, &fnx[0], &xyz_pq[0], &rints[0]);\n\n'

			# # First rollout
			# file_str += rolloutERI4FirstDepr(la, lb, lc, ld)
			# file_str += '                }\n'
			# file_str += '        }\n\n'

			# file_str += '    vec4d eri4_batch(Fill(0), n_sph_a, n_sph_b, n_sph_c, n_sph_d);\n'
			# file_str += '    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)\n'
			# file_str += '        for (int ib = 0; ib < cdepth_b; ib++, iab++)\n'
			# file_str += '        {\n'

			# # Second rollout
			# file_str += rolloutERI4Second(la, lb, lc, ld)			
			# file_str += '        }\n\n'	

			# file_str += '    int ofs_norm_a = sp_data_ab.offsets_norms[2 * ipair_ab];\n'
			# file_str += '    int ofs_norm_b = sp_data_ab.offsets_norms[2 * ipair_ab + 1];\n'
			# file_str += '    int ofs_norm_c = sp_data_cd.offsets_norms[2 * ipair_cd];\n'
			# file_str += '    int ofs_norm_d = sp_data_cd.offsets_norms[2 * ipair_cd + 1];\n'
			# file_str += '    for (int mu = 0; mu < n_sph_a; mu++)\n'
			# file_str += '        for (int nu = 0; nu < n_sph_b; nu++)\n'
			# file_str += '            for (int ka = 0; ka < n_sph_c; ka++)\n'
			# file_str += '                for (int ta = 0; ta < n_sph_d; ta++)\n'
			# file_str += '                {\n'
			# file_str += '                    double norm_a = sp_data_ab.norms[ofs_norm_a + mu];\n'
			# file_str += '                    double norm_b = sp_data_ab.norms[ofs_norm_b + nu];\n'
			# file_str += '                    double norm_c = sp_data_cd.norms[ofs_norm_c + ka];\n'
			# file_str += '                    double norm_d = sp_data_cd.norms[ofs_norm_d + ta];\n'
			# file_str += '                    eri4_batch(mu, nu, ka, ta) *= norm_a * norm_b * norm_c * norm_d;\n'
			# file_str += '                }\n\n'
		
			# file_str += '    return eri4_batch;\n'
			# file_str += '}\n\n'
			 

	# Write ERI3 kernels
	for la, lb in lalb_list:

		if (lbra + lket) <= 6:
			file_str += instantiateERI3(la, lb, lket)
			file_str += instantiateERI3D1(la, lb, lket)
			file_str += instantiateERI3D2(la, lb, lket)
			file_str += instantiateERI3SOC(la, lb, lket)
			if la != lb:
				file_str += instantiateERI3(lb, la, lket)
				file_str += instantiateERI3D1(lb, la, lket)
				file_str += instantiateERI3D2(lb, la, lket)
				file_str += instantiateERI3SOC(lb, la, lket)

		# file_str += 'template<> lible::vec3d\n'
		# file_str += 'lible::ints::two::eri3Kernel<{}, {}, {}>(const int ipair_ab, const int ishell_c,\n'.format(la, lb, lket)
		# file_str += '                                      const std::vector<double> &ecoeffs_ab,\n' 
		# file_str += '                                      const std::vector<double> &ecoeffs_c,\n'
		# file_str += '                                      const ShellPairData &sp_data_ab,\n'
		# file_str += '                                      const ShellData &sh_data_c)\n'		
		# file_str += '{\n'
		# file_str += '    const int cdepth_a = sp_data_ab.cdepths[2 * ipair_ab];\n'
		# file_str += '    const int cdepth_b = sp_data_ab.cdepths[2 * ipair_ab + 1];\n'
		# file_str += '    const int cdepth_c = sh_data_c.cdepths[ishell_c];\n'
		# file_str += '    const int cofs_a = sp_data_ab.coffsets[2 * ipair_ab];\n'
		# file_str += '    const int cofs_b = sp_data_ab.coffsets[2 * ipair_ab + 1];\n'
		# file_str += '    const int cofs_c = sh_data_c.coffsets[ishell_c];\n\n'

		# file_str += '    const double *exps_a = &sp_data_ab.exps[cofs_a];\n'
		# file_str += '    const double *exps_b = &sp_data_ab.exps[cofs_b];\n'
		# file_str += '    const double *exps_c = &sh_data_c.exps[cofs_c];\n'
		# file_str += '    const double *coords_a = &sp_data_ab.coords[6 * ipair_ab];\n'
		# file_str += '    const double *coords_b = &sp_data_ab.coords[6 * ipair_ab + 3];\n'
		# file_str += '    const double *coords_c = &sh_data_c.coords[3 * ishell_c];\n'
		# file_str += '    const double *pecoeffs_ab = &ecoeffs_ab[sp_data_ab.offsets_ecoeffs[ipair_ab]];\n'
		# file_str += '    const double *pecoeffs_c = &ecoeffs_c[sh_data_c.offsets_ecoeffs[ishell_c]];\n\n'
				
		# file_str += '    constexpr int la = {}, lb = {}, lc = {};\n'.format(la, lb, lc)
		# file_str += '    constexpr int lab = la + lb;\n'
		# file_str += '    constexpr int labc = lab + lc;\n'

		# file_str += '    constexpr int n_sph_a = numSphericalsC(la);\n'
		# file_str += '    constexpr int n_sph_b = numSphericalsC(lb);\n'
		# file_str += '    constexpr int n_sph_c = numSphericalsC(lc);\n'
		# file_str += '    constexpr int n_hermite_ab = numHermitesC(lab);\n'
		# file_str += '    constexpr int n_hermite_c = numHermitesC(lc);\n'
		# file_str += '    constexpr int n_sph_ab = n_sph_a * n_sph_b;\n'
		# file_str += '    constexpr int n_ecoeffs_ab = n_sph_ab * n_hermite_ab;\n'
		# file_str += '    constexpr int n_ecoeffs_c = n_sph_c * n_hermite_c;\n\n'

		# file_str += '    std::array<double, labc + 1> fnx;\n'
		# file_str += '    BoysF2<labc> boys_f;\n\n'

		# file_str += '    constexpr int n_hermites_abc = numHermitesC(lab) * numHermitesC(lc);\n'
		# file_str += '    std::array<double, n_hermites_abc> rints;\n\n'

		# file_str += '    constexpr int n_rints_x_ecoeffs = n_sph_c * n_hermite_ab;\n'
		# file_str += '    std::vector<double> rints_x_ecoeffs(cdepth_a * cdepth_b * n_rints_x_ecoeffs, 0);\n\n'

		# file_str += '    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)\n'
		# file_str += '        for (int ib = 0; ib < cdepth_b; ib++, iab++)\n'
		# file_str += '        {\n'
		# file_str += '            int pos_rints_x_ecoeffs = iab * n_rints_x_ecoeffs;\n'
		# file_str += '            for (int ic = 0; ic < cdepth_c; ic++)\n'
		# file_str += '            {\n'
		# file_str += '                double a = exps_a[ia];\n'
		# file_str += '                double b = exps_b[ib];\n'
		# file_str += '                double c = exps_c[ic];\n\n'

		# file_str += '                double p = a + b;\n'
		# file_str += '                double alpha = p * c / (p + c);\n\n'

		# file_str += '                std::array<double, 3> xyz_p{(a * coords_a[0] + b * coords_b[0]) / p,\n'
		# file_str += '                                            (a * coords_a[1] + b * coords_b[1]) / p,\n'
		# file_str += '                                            (a * coords_a[2] + b * coords_b[2]) / p};\n\n'

		# file_str += '                std::array<double, 3> xyz_pc{xyz_p[0] - coords_c[0],\n'
		# file_str += '                                             xyz_p[1] - coords_c[1],\n'
		# file_str += '                                             xyz_p[2] - coords_c[2]};\n\n'

		# file_str += '                double dx{xyz_pc[0]}, dy{xyz_pc[1]}, dz{xyz_pc[2]};\n'
		# file_str += '                double x = alpha * (dx * dx + dy * dy + dz * dz);\n'
		# file_str += '                boys_f.calcFnx(x, &fnx[0]);\n\n'

		# file_str += '                double fac = (2.0 * std::pow(M_PI, 2.5) / (p * c * std::sqrt(p + c)));\n'
		# file_str += '                calcRInts_ERI<lab, lc>(alpha, fac, &fnx[0], &xyz_pc[0], &rints[0]);\n\n'

		# # First rollout
		# file_str += rolloutERI3First(la, lb, lket)
		# file_str += '            }\n'
		# file_str +=         '}\n\n'

		# file_str += '    vec3d eri3_batch(Fill(0), n_sph_a, n_sph_b, n_sph_c);\n'
		# file_str += '    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)\n'
		# file_str += '        for (int ib = 0; ib < cdepth_b; ib++, iab++)\n'
		# file_str += '        {\n'

		# # Second rollout
		# file_str += rolloutERI3Second(la, lb, lket)						
		# file_str += '        }\n\n'

		# file_str += '    int ofs_norm_a = sp_data_ab.offsets_norms[2 * ipair_ab];\n'
		# file_str += '    int ofs_norm_b = sp_data_ab.offsets_norms[2 * ipair_ab + 1];\n'
		# file_str += '    int ofs_norm_c = sh_data_c.offsets_norms[ishell_c];\n'
		# file_str += '    for (int mu = 0; mu < n_sph_a; mu++)\n'
		# file_str += '        for (int nu = 0; nu < n_sph_b; nu++)\n'
		# file_str += '            for (int ka = 0; ka < n_sph_c; ka++)\n'
		# file_str += '            {\n'
		# file_str += '                double norm_a = sp_data_ab.norms[ofs_norm_a + mu];\n'
		# file_str += '                double norm_b = sp_data_ab.norms[ofs_norm_b + nu];\n'
		# file_str += '                double norm_c = sh_data_c.norms[ofs_norm_c + ka];\n'
		# file_str += '                eri3_batch(mu, nu, ka) *= norm_a * norm_b * norm_c;\n'
		# file_str += '            }\n\n'
		
		# file_str += '    return eri3_batch;\n'
		# file_str += '}\n\n'


	# Write ERI2 kernels
	if (lbra + lket) <= 6:
		file_str += instantiateERI2(lbra, lket)
		file_str += instantiateERI2D1(lbra, lket)
		file_str += instantiateERI2D2(lbra, lket)

	# file_str += 'template<> lible::vec2d\n'
	# file_str += 'lible::ints::two::eri2Kernel<{}, {}>(const int ishell_a, const int ishell_b,\n'.format(lbra, lket)
	# file_str += '                                   const std::vector<double> &ecoeffs_a,\n'
	# file_str += '                                   const std::vector<double> &ecoeffs_b_tsp,\n'
	# file_str += '                                   const ShellData &sh_data_a, const ShellData &sh_data_b)\n'	
	# file_str += '{\n'
	# file_str += '    const int cdepth_a = sh_data_a.cdepths[ishell_a];\n'
	# file_str += '    const int cdepth_b = sh_data_b.cdepths[ishell_b];\n'
	# file_str += '    const int cofs_a = sh_data_a.coffsets[ishell_a];\n'
	# file_str += '    const int cofs_b = sh_data_b.coffsets[ishell_b];\n\n'

	# file_str += '    const double *exps_a = &sh_data_a.exps[cofs_a];\n'
	# file_str += '    const double *exps_b = &sh_data_b.exps[cofs_b];\n'
	# file_str += '    const double *coords_a = &sh_data_a.coords[3 * ishell_a];\n'
	# file_str += '    const double *coords_b = &sh_data_b.coords[3 * ishell_b];\n'
	# file_str += '    const double *pecoeffs_a = &ecoeffs_a[sh_data_a.offsets_ecoeffs[ishell_a]];\n'
	# file_str += '    const double *pecoeffs_b_tsp = &ecoeffs_b_tsp[sh_data_b.offsets_ecoeffs[ishell_b]];\n\n'
		
	# file_str += '    constexpr int la = {}, lb = {};\n'.format(lbra, lket)
	# file_str += '    constexpr int lab = la + lb;\n'
	# file_str += '    constexpr int n_sph_a = numSphericalsC(la);\n'
	# file_str += '    constexpr int n_sph_b = numSphericalsC(lb);\n'
	# file_str += '    constexpr int n_hermite_a = numHermitesC(la);\n'
	# file_str += '    constexpr int n_hermite_b = numHermitesC(lb);\n'
	# file_str += '    constexpr int n_ecoeffs_a = n_sph_a * n_hermite_a;\n'
	# file_str += '    constexpr int n_ecoeffs_b = n_sph_b * n_hermite_b;\n\n'

	# file_str += '    std::array<double, lab + 1> fnx;\n'
	# file_str += '    BoysF2<lab> boys_f;\n\n'

	# file_str += '    constexpr int n_hermites_ab = numHermitesC(la) * numHermitesC(lb);\n'
	# file_str += '    std::array<double, n_hermites_ab> rints;\n\n'

	# file_str += '    constexpr int n_rints_x_ecoeffs = n_hermite_a * n_sph_b;\n'
	# file_str += '    std::vector<double> rints_x_ecoeffs(cdepth_a * n_rints_x_ecoeffs, 0);\n\n'

	# file_str += '    std::array<double, 3> xyz_ab{coords_a[0] - coords_b[0],\n'
	# file_str += '    coords_a[1] - coords_b[1],\n'
	# file_str += '    coords_a[2] - coords_b[2]};\n\n'

	# file_str += '    double dx{xyz_ab[0]}, dy{xyz_ab[1]}, dz{xyz_ab[2]};\n'
	# file_str += '    double xyz_ab_dot = dx * dx + dy * dy + dz * dz;	\n\n'

	# file_str += '    for (int ia = 0; ia < cdepth_a; ia++)\n'
	# file_str += '    {\n'
	# file_str += '        int pos_rints_x_ecoeffs = ia * n_rints_x_ecoeffs;\n'
	# file_str += '        for (int ib = 0; ib < cdepth_b; ib++)\n'
	# file_str += '        {\n'
	# file_str += '            double a = exps_a[ia];\n'
	# file_str += '            double b = exps_b[ib];\n\n'

	# file_str += '            double alpha = a * b / (a + b);\n'
	# file_str += '            double x = alpha * xyz_ab_dot;\n'
	# file_str += '            boys_f.calcFnx(x, &fnx[0]);\n\n'

	# file_str += '            double fac = (2.0 * std::pow(M_PI, 2.5) / (a * b * std::sqrt(a + b)));\n'
	# file_str += '            calcRInts_ERI<la, lb>(alpha, fac, &fnx[0], &xyz_ab[0], &rints[0]);\n\n'

	# # First rollout
	# file_str += rolloutERI2First(lbra, lket)
	# file_str += '        }\n'
	# file_str += '    }\n\n'

	# file_str += '    vec2d eri2_batch(Fill(0), n_sph_a, n_sph_b);\n'
	# file_str += '    for (int ia = 0; ia < cdepth_a; ia++)\n'
	# file_str += '    {\n'

	# # Second rollout
	# file_str += rolloutERI2Second(lbra, lket)
	# file_str += '    }\n\n'

	# file_str += '    int ofs_norm_a = sh_data_a.offsets_norms[ishell_a];\n'
	# file_str += '    int ofs_norm_b = sh_data_b.offsets_norms[ishell_b];\n'
	# file_str += '    for (int mu = 0; mu < n_sph_a; mu++)\n'
	# file_str += '        for (int nu = 0; nu < n_sph_b; nu++)\n'
	# file_str += '        {\n'
	# file_str += '            double norm_a = sh_data_a.norms[ofs_norm_a + mu];\n'
	# file_str += '            double norm_b = sh_data_b.norms[ofs_norm_b + nu];\n'
	# file_str += '            eri2_batch(mu, nu) *= norm_a * norm_b;\n'
	# file_str += '        }\n\n'

	# file_str += '    return eri2_batch;\n'
	# file_str += '}\n\n'

	# file_str += 'template std::array<lible::vec2d, 6> lible::ints::two::eri2d1Kernel<{}, {}>(const int ishell_a, const int ishell_b,\n'.format(lbra, lket)
	# file_str += '                                                           const std::vector<double> &ecoeffs_a,\n'
	# file_str += '                                                           const std::vector<double> &ecoeffs_b_tsp,\n'
	# file_str += '                                                           const ShellData &sh_data_a,\n'
	# file_str += '                                                           const ShellData &sh_data_b);\n\n'			
	
	with open('eri_kernels_{}_{}.cpp'.format(lbra, lket), 'w') as file:
		file.write(file_str)

l_max = 6
l_max_generate = 4 # TODO: set to 6 later
for lbra in range(0, l_max + 1):
	for lket in range(0, l_max + 1): 
		l_sum = lbra + lket
		if l_sum > l_max:
			continue
		# writeKernelInstantiate(lbra, lket)
		if l_sum > l_max_generate:
			writeKernelInstantiate(lbra, lket)
		else:
			writeKernelGenerate(lbra, lket)
