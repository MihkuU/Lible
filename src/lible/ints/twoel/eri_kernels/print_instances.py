
l_max = 6
l_max_bra = 12
l_max_ket = 12

# ERI4 kernels
for la in range(0, l_max + 1):
	for lb in range(0, la + 1):
		for lc in range(0, l_max + 1):
			for ld in range(0, lc + 1):
				#print('{{{{{}, {}, {}, {}}}, eri4KernelFun<{}, {}, {}, {}>}},'.format(la, lb, lc, ld, la, lb, lc, ld))
				continue

# ERI3 kernels
for la in range(0, l_max + 1):
	for lb in range(0, la + 1):
		for lc in range(0, l_max_ket + 1):
			#print('{{{{{}, {}, {}}}, eri3KernelFun<{}, {}, {}>}},'.format(la, lb, lc, la, lb, lc))
			continue
			#print('{{{{{}, {}, {}}}, eri3Kernel<{}, {}, {}>}},'.format(la, lb, lc, la, lb, lc))
		

# ERI2 kernels
for la in range(0, l_max_ket + 1):
	for lb in range(0, l_max_ket + 1):
		print('{{{{{}, {}}}, eri2KernelFun<{}, {}>}},'.format(la, lb, la, lb))
		continue
		#print('{{{{{}, {}}}, eri2Kernel<{}, {}>}},'.format(la, lb, la, lb))
		#print('{{{{{}, {}}}, eri2d1Kernel<{}, {}>}},'.format(la, lb, la, lb))
