
l_max = 6
l_max_bra = 12
l_max_ket = 12

# ERI4 kernels
for la in range(0, l_max + 1):
	for lb in range(0, la + 1):
		for lc in range(0, l_max + 1):
			for ld in range(0, lc + 1):
				continue
				#print('{{{{{}, {}, {}, {}}}, eri4Kernel<{}, {}, {}, {}>}},'.format(la, lb, lc, ld, la, lb, lc, ld))

# ERI3 kernels
for la in range(0, l_max + 1):
	for lb in range(0, la + 1):
		for lc in range(0, l_max_ket + 1):
			continue
			#print('{{{{{}, {}, {}}}, eri3Kernel<{}, {}, {}>}},'.format(la, lb, lc, la, lb, lc))
		

# ERI2 kernels
for la in range(0, l_max_ket + 1):
	for lb in range(0, l_max_ket + 1):
		#continue
		#print('{{{{{}, {}}}, eri2Kernel<{}, {}>}},'.format(la, lb, la, lb))
		print('{{{{{}, {}}}, eri2d1Kernel<{}, {}>}},'.format(la, lb, la, lb))
