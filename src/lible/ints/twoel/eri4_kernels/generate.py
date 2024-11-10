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

	with open('eri4_kernel_{}_{}.cpp'.format(lbra, lket), 'w') as file:
		file_str = ''
		file_str += '#include <lible/ints/twoel/eri4_kernel.hpp>\n\n'

		for la, lb in lalb_list:
			for lc, ld in lcld_list:
				file_str += 'template void lible::ints::two::eri4Kernel<{}, {}, {}, {}>(const int, const int, const int, const int,\n'.format(la, lb, lc, ld)
				file_str += '                                                       const double*, const double*,\n'
				file_str += '                                                       const double*, const double*,\n'
				file_str += '                                                       const double*, const double*,\n'
				file_str += '                                                       const double*, const double*,\n'
				file_str += '                                                       const double*, const double*,\n'
				file_str += '                                                       double*);\n\n'

		file.write(file_str)


l_max = 12
l_max_generate = 6
for lbra in range(0, l_max + 1):
	for lket in range(0, l_max + 1): 
		l_sum = lbra + lket
		writeKernelInstantiate(lbra, lket)
		#if l_sum > l_max_generate:
		#	writeKernelInstantiate(lbra, lket)
		#else:
		#	continue
