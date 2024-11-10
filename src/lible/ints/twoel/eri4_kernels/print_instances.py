
l_max = 6
l_max_bra = 12
l_max_ket = 12

for la in range(0, l_max + 1):
	for lb in range(0, la + 1):
		lbra = la + lb
		for lc in range(0, l_max + 1):
			for ld in range(0, lc + 1):
				print('{{{{{}, {}, {}, {}}}, eri4Kernel<{}, {}, {}, {}>}},'.format(la, lb, lc, ld, la, lb, lc, ld))

