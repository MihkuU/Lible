from utils import *

def instantiateERI4(la, lb, lc, ld):

	file_str = 'template lible::vec4d lible::ints::eri4KernelFun<{}, {}, {}, {}>(const size_t ipair_ab, const size_t ipair_cd,\n'.format(la, lb, lc, ld)
	file_str += '                                                             const ShellPairData &sp_data_ab,\n'
	file_str += '                                                             const ShellPairData &sp_data_cd,\n'
	file_str += '                                                             const ERI4Kernel *eri4_kernel);\n\n'

	return file_str

def instantiateERI3(la, lb, lc):

	file_str = 'template lible::vec3d lible::ints::eri3KernelFun<{}, {}, {}>(const size_t ipair_ab, const size_t ishell_c,\n'.format(la, lb, lc)
	file_str += '                                                          const ShellPairData &sp_data_ab,\n'
	file_str += '                                                          const ShellData &sh_data_c,\n'
	file_str += '                                                          const ERI3Kernel *eri3_kernel);\n\n'	

	return file_str

def instantiateERI2(la, lb):

	file_str = 'template lible::vec2d lible::ints::eri2KernelFun<{}, {}>(const size_t ishell_a, const size_t ishell_b,\n'.format(la, lb)
	file_str += '                                                       const ShellData &sh_data_a,\n'
	file_str += '                                                       const ShellData &sh_data_b,\n'
	file_str += '                                                       const ERI2Kernel *eri2_kernel);\n\n'
	
	return file_str 

def instantiateERI2D1(la, lb):

	file_str = 'template std::array<lible::vec2d, 6> lible::ints::eri2d1KernelFun<{}, {}>(const size_t ishell_a, const size_t ishell_b,\n'.format(la, lb)
	file_str += '                                                                        const ShellData &sh_data_a,\n'
	file_str += '                                                                        const ShellData &sh_data_b,\n'
	file_str += '                                                                        const ERI2D1Kernel *eri2d1_kernel);\n\n'

	return file_str

def instantiateERI2D2(la, lb):

	file_str  = 'template std::array<std::array<lible::vec2d, 6>, 6> \n'
	file_str += 'lible::ints::eri2d2KernelFun<{}, {}>(const size_t ishell_a, const size_t ishell_b,\n'.format(la, lb)
	file_str += '                                   const ShellData &sh_data_a,\n'
	file_str += '                                   const ShellData &sh_data_b,\n'
	file_str += '                                   const ERI2D2Kernel *eri2d2_kernel);\n\n'

	return file_str

def instantiateERI3D1(la, lb, lc):

	file_str = 'template std::array<lible::vec3d, 9> lible::ints::eri3d1KernelFun<{}, {}, {}>(const size_t ipair_ab, const size_t ishell_c,\n'.format(la, lb, lc)
	file_str += '                                                                           const ShellPairData &sh_data_ab,\n'
	file_str += '                                                                           const ShellData &sh_data_c,\n'
	file_str += '                                                                           const ERI3D1Kernel *eri3d1_kernel);\n\n'

	return file_str

def instantiateERI3D2(la, lb, lc):

	file_str = 'template lible::arr2d<lible::vec3d, 9, 9>\n'
	file_str += 'lible::ints::eri3d2KernelFun<{}, {}, {}>(const size_t ipair_ab, const size_t ishell_c,\n'.format(la, lb, lc)
	file_str += '                                      const ShellPairData &sh_data_ab,\n'
	file_str += '                                      const ShellData &sh_data_c,\n'
	file_str += '                                      const ERI3D2Kernel *eri3d2_kernel);\n\n'

	return file_str

def instantiateERI3SOC(la, lb, lc):

	file_str = 'template std::array<lible::vec3d, 3> lible::ints::eri3socKernelFun<{}, {}, {}>(const size_t ipair_ab, const size_t ishell_c,\n'.format(la, lb, lc)
	file_str += '                                                                            const ShellPairData &sh_data_ab,\n'
	file_str += '                                                                            const ShellData &sh_data_c,\n'
	file_str += '                                                                            const ERI3SOCKernel *eri3soc_kernel);\n\n'

	return file_str

def instantiateERI4D1(la, lb, lc, ld):

	file_str = 'template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<{}, {}, {}, {}>(const size_t ipair_ab, const size_t ipair_cd,\n'.format(la, lb, lc, ld)
	file_str += '                                                                               const ShellPairData &sh_data_ab,\n'
	file_str += '                                                                               const ShellPairData &sp_data_cd,\n'
	file_str += '                                                                               const ERI4D1Kernel *eri4d1_kernel);\n\n'

	return file_str

def instantiateERI4SOC(la, lb, lc, ld):

	file_str = 'template std::array<lible::vec4d, 3> lible::ints::eri4socKernelFun<{}, {}, {}, {}>(const size_t ipair_ab, const size_t ipair_cd,\n'.format(la, lb, lc, ld)
	file_str += '                                                                                const ShellPairData &sh_data_ab,\n'
	file_str += '                                                                                const ShellPairData &sp_data_cd,\n'
	file_str += '                                                                                const ERI4SOCKernel *eri4soc_kernel);\n\n'

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
		
		# Write ERI3 kernels
		for la, lb in lalb_list:
			
			if (lbra + lket) <= 6:
				file_str += instantiateERI3(la, lb, lket)
				file_str += instantiateERI3D1(la, lb, lket)
				file_str += instantiateERI3SOC(la, lb, lket)
				if la != lb:
					file_str += instantiateERI3(lb, la, lket)
					file_str += instantiateERI3D1(lb, la, lket)
					file_str += instantiateERI3SOC(lb, la, lket)

		# Write ERI2 kernels
		if (lbra + lket) <= 6:
			file_str += instantiateERI2(lbra, lket)
			file_str += instantiateERI2D1(lbra, lket)
			file_str += instantiateERI2D2(lbra, lket)

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

	# Write ERI3 kernels
	for la, lb in lalb_list:

		if (lbra + lket) <= 6:
			file_str += instantiateERI3(la, lb, lket)
			file_str += instantiateERI3D1(la, lb, lket)
			file_str += instantiateERI3SOC(la, lb, lket)
			if la != lb:
				file_str += instantiateERI3(lb, la, lket)
				file_str += instantiateERI3D1(lb, la, lket)
				file_str += instantiateERI3SOC(lb, la, lket)

	# Write ERI2 kernels
	if (lbra + lket) <= 6:
		file_str += instantiateERI2(lbra, lket)
		file_str += instantiateERI2D1(lbra, lket)
		file_str += instantiateERI2D2(lbra, lket)

	with open('eri_kernels_{}_{}.cpp'.format(lbra, lket), 'w') as file:
		file.write(file_str)

l_max = 6
for lbra in range(0, l_max + 1):
	for lket in range(0, l_max + 1): 
		l_sum = lbra + lket
		if l_sum > l_max:
			continue

		writeKernelInstantiate(lbra, lket)
