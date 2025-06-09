import utils

l_max = 6

def writeKernelsBra(lbra, lket):

	file_str = '#include <lible/ints/twoel/shark_mm_kernels.hpp>\n\n'

	file_str +=	utils.rolloutSHARKBra1(lbra, lket)

	for la in range(0, lbra + 1):
		for lb in range(0, lbra + 1):
			if ((la + lb) == lbra):
				file_str += utils.rolloutSHARKBra2_ERI3(la, lb, lket)

	with open('shark_mm_kernels_bra_{}_{}.cpp'.format(lbra, lket), 'w') as file:
		file.write(file_str)

def writeKernelsKet(lbra, lket):

	file_str = '#include <lible/ints/twoel/shark_mm_kernels.hpp>\n\n'

	file_str +=	utils.rolloutSHARKKet1(lbra, lket)

	for lc in range(0, lket + 1):
		for ld in range(0, lket + 1):
			if ((lc + ld) == lket):
				file_str += utils.rolloutSHARKKet2(lbra, lc, ld)

	with open('shark_mm_kernels_ket_{}_{}.cpp'.format(lbra, lket), 'w') as file:
		file.write(file_str)

def writeKernelsBraKet(lbra, lket):

	file_str = '#include <lible/ints/twoel/shark_mm_kernels.hpp>\n\n'

	for la in range(0, lbra + 1):
		for lb in range(0, lbra + 1):
			for lc in range(0, lket + 1):
				for ld in range(0, lket + 1):
					if (la + lb) == lbra and (lc + ld) == lket:
						file_str += utils.rolloutSHARKBra2_ERI4(la, lb, lc, ld)


	with open('shark_mm_kernels_braket_{}_{}.cpp'.format(lbra, lket), 'w') as file:
		file.write(file_str)

for lbra in range(0, l_max + 1):
	for lket in range(0, l_max + 1):
		if (lbra + lket) <= l_max:
			print('lbra = ', lbra, ' lket = ', lket)
			writeKernelsBra(lbra, lket)
			writeKernelsKet(lbra, lket)
			writeKernelsBraKet(lbra, lket)		
