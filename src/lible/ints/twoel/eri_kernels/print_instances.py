
l_max = 6
l_max_bra = 12
l_max_ket = 12

# ERI4 kernels
for la in range(0, l_max + 1):
	for lb in range(0, la + 1):
		for lc in range(0, l_max + 1):
			for ld in range(0, lc + 1):

				if (la + lb + lc + ld) <= 6:
					continue

				"""
					print('{{{{{}, {}, {}, {}}}, eri4socKernelFun<{}, {}, {}, {}>}},'.format(la, lb, lc, ld, la, lb, lc, ld))
					if (la != lb) and (lc == ld):
						print('{{{{{}, {}, {}, {}}}, eri4socKernelFun<{}, {}, {}, {}>}},'.format(lb, la, lc, ld, lb, la, lc, ld))
					elif (lc != ld) and (la == lb):
						print('{{{{{}, {}, {}, {}}}, eri4socKernelFun<{}, {}, {}, {}>}},'.format(la, lb, ld, lc, la, lb, ld, lc))					
					elif (la != lb) and (lc != ld):
						print('{{{{{}, {}, {}, {}}}, eri4socKernelFun<{}, {}, {}, {}>}},'.format(lb, la, ld, lc, lb, la, ld, lc))
						print('{{{{{}, {}, {}, {}}}, eri4socKernelFun<{}, {}, {}, {}>}},'.format(la, lb, ld, lc, la, lb, ld, lc))
						print('{{{{{}, {}, {}, {}}}, eri4socKernelFun<{}, {}, {}, {}>}},'.format(lb, la, lc, ld, lb, la, lc, ld))
				"""

				"""
					print('{{{{{}, {}, {}, {}}}, eri4KernelFun<{}, {}, {}, {}>}},'.format(la, lb, lc, ld, la, lb, lc, ld))
					if (la != lb) and (lc == ld):
						print('{{{{{}, {}, {}, {}}}, eri4KernelFun<{}, {}, {}, {}>}},'.format(lb, la, lc, ld, lb, la, lc, ld))
					elif (lc != ld) and (la == lb):
						print('{{{{{}, {}, {}, {}}}, eri4KernelFun<{}, {}, {}, {}>}},'.format(la, lb, ld, lc, la, lb, ld, lc))					
					elif (la != lb) and (lc != ld):
						print('{{{{{}, {}, {}, {}}}, eri4KernelFun<{}, {}, {}, {}>}},'.format(lb, la, ld, lc, lb, la, ld, lc))
						print('{{{{{}, {}, {}, {}}}, eri4KernelFun<{}, {}, {}, {}>}},'.format(la, lb, ld, lc, la, lb, ld, lc))
						print('{{{{{}, {}, {}, {}}}, eri4KernelFun<{}, {}, {}, {}>}},'.format(lb, la, lc, ld, lb, la, lc, ld))
				"""
				"""
					print('{{{{{}, {}, {}, {}}}, eri4d1KernelFun<{}, {}, {}, {}>}},'.format(la, lb, lc, ld, la, lb, lc, ld))
					if (la != lb) and (lc == ld):
						print('{{{{{}, {}, {}, {}}}, eri4d1KernelFun<{}, {}, {}, {}>}},'.format(lb, la, lc, ld, lb, la, lc, ld))
					elif (lc != ld) and (la == lb):
						print('{{{{{}, {}, {}, {}}}, eri4d1KernelFun<{}, {}, {}, {}>}},'.format(la, lb, ld, lc, la, lb, ld, lc))					
					elif (la != lb) and (lc != ld):
						print('{{{{{}, {}, {}, {}}}, eri4d1KernelFun<{}, {}, {}, {}>}},'.format(lb, la, ld, lc, lb, la, ld, lc))
						print('{{{{{}, {}, {}, {}}}, eri4d1KernelFun<{}, {}, {}, {}>}},'.format(la, lb, ld, lc, la, lb, ld, lc))
						print('{{{{{}, {}, {}, {}}}, eri4d1KernelFun<{}, {}, {}, {}>}},'.format(lb, la, lc, ld, lb, la, lc, ld))
				"""

# ERI3 kernels
for la in range(0, l_max + 1):
	for lb in range(0, la + 1):
		for lc in range(0, l_max_ket + 1):			
			if (la + lb + lc) > 6:
				continue

			#print('{{{{{}, {}, {}}}, eri3KernelFun<{}, {}, {}>}},'.format(la, lb, lc, la, lb, lc))
			#if (la != lb):
			#	print('{{{{{}, {}, {}}}, eri3KernelFun<{}, {}, {}>}},'.format(lb, la, lc, lb, la, lc))
			if (la + lb + lc) <= 6:
				#print('{{{{{}, {}, {}}}, eri3socKernelFun<{}, {}, {}>}},'.format(la, lb, lc, la, lb, lc))
				#print('{{{{{}, {}, {}}}, eri3d1KernelFun<{}, {}, {}>}},'.format(la, lb, lc, la, lb, lc))
				print('{{{{{}, {}, {}}}, eri3d2KernelFun<{}, {}, {}>}},'.format(la, lb, lc, la, lb, lc))
				if (la != lb):
					#print('{{{{{}, {}, {}}}, eri3socKernelFun<{}, {}, {}>}},'.format(lb, la, lc, lb, la, lc))
					#print('{{{{{}, {}, {}}}, eri3d1KernelFun<{}, {}, {}>}},'.format(lb, la, lc, lb, la, lc))
					print('{{{{{}, {}, {}}}, eri3d2KernelFun<{}, {}, {}>}},'.format(lb, la, lc, lb, la, lc))
			#	print('{{{{{}, {}, {}}}, eri3d1KernelFun<{}, {}, {}>}},'.format(la, lb, lc, la, lb, lc))
			#	print('{{{{{}, {}, {}}}, eri3d1KernelFun<{}, {}, {}>}},'.format(lb, la, lc, lb, la, lc))
			#	print('{{{{{}, {}, {}}}, eri3KernelFun<{}, {}, {}>}},'.format(la, lb, lc, la, lb, lc))
			#print('{{{{{}, {}, {}}}, eri3d1KernelFun<{}, {}, {}>}},'.format(lb, la, lc, lb, la, lc))
			continue
			#print('{{{{{}, {}, {}}}, eri3Kernel<{}, {}, {}>}},'.format(la, lb, lc, la, lb, lc))
		

# ERI2 kernels
for la in range(0, l_max_ket + 1):
	for lb in range(0, l_max_ket + 1):
		#print('{{{{{}, {}}}, eri2KernelFun<{}, {}>}},'.format(la, lb, la, lb))

		#if (la + lb) <= 6:			
			#print('{{{{{}, {}}}, eri2KernelFun<{}, {}>}},'.format(la, lb, la, lb))
			#print('{{{{{}, {}}}, eri2d1KernelFun<{}, {}>}},'.format(la, lb, la, lb))
			#print('{{{{{}, {}}}, eri2d2KernelFun<{}, {}>}},'.format(la, lb, la, lb))

		continue
		#print('{{{{{}, {}}}, eri2Kernel<{}, {}>}},'.format(la, lb, la, lb))
		#print('{{{{{}, {}}}, eri2d1Kernel<{}, {}>}},'.format(la, lb, la, lb))
