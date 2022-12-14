extends Node



var matrix_math
### Neo-Hookean...
var f_coefficient : Array
var pk_mu_coefficient : Array
var lambda_log_cofficient : float 
var pk_lambda_coefficient : Array
var P : Array
var j_coefficient : float
var stress_coefficient : Array
var f_transposed_coefficient : Array 
###


func _on_constitutive_models_ready():
	### contains the models to simulate various materials/substances...
	var matrix_math = get_tree().get_root().get_node("Simulation/Matrix Math")
	pass


func Model_of_Water(speck,nub):
	var resulted_stress
	matrix_math = get_tree().get_root().get_node("Simulation/Matrix Math")
	## stress = -p * I ,
	#:: p = k * ( (1 / J ^ y) - 1 )
	#
	# I == Identity Matrix
	# k == bulk modulus
	# J == Determinant(F)
	# y == stiffly penalizes large deviations from incompressibility
	
	# water properties
	# bulk modulus = 2.1
	var y = 3.3 * pow(10,-6)
	#print(nub.particle_mechanics[speck]['J'],' testing J')
	#var deteminant_coefficient = 1.0 / pow(nub.particle_mechanics[speck]['J'],3.0)
	var deteminant_coefficient = 1.0 / pow(nub.particle_mechanics[speck]['J'],y)
	#print(deteminant_coefficient,' deteminant coefficient')
	var p_coefficient = 2.1 * (deteminant_coefficient - 1.0 )
	#print(p_coefficient,' p')
	resulted_stress = matrix_math.Multiply_Matrix_by_Scalar(nub.particle_mechanics[speck]['I'],-p_coefficient,true)
	#print(resulted_stress)
	return resulted_stress
	

func Neo_Hookean(speck,nub):
	### Neo-Hookean model nonlinear hyperelastic models for predicting large deformations of elastic materials. mpm pg. 19
	### example "rubber
	var resulted_stress
	var mu
	var lambda
	matrix_math = get_tree().get_root().get_node("Simulation/Matrix Math")
	#var pk_stress
	
	#var F_transposed = get_tree().get_root().get_node("Simulation/Matrix Math").Transposed_Matrix(nub.particle_mechanics[speck]['F'])
	
	mu = snapped((nub.youngs_modulus / snapped((2.0 * snapped((1.0+nub.poisson_ratio),.01)),.01) ),.01)
	lambda = snapped(snapped((nub.youngs_modulus * nub.poisson_ratio),.01) / ( (1.0 + nub.poisson_ratio) * clamp(snapped((1.0-snapped((2.0 * nub.poisson_ratio),.01)),.01),0.1,snapped((1.0-snapped((2.0 * nub.poisson_ratio),.01)),.01)+1) ),.01)
	#'''
	#  Piola-Kirchoff stress = u*F + (n * In(J) - u )*F^-T
	# note ln() = log()
	# or
	# P = u*(F - F^-T) + n*log (J)*F^-T ; mpm course pg. 19,equation 48
	#f_coefficient =  matrix_math.Subtract_Matrix(nub.particle_mechanics[speck]['F'],matrix_math.Inverse_Matrix(matrix_math.Transposed_Matrix(nub.particle_mechanics[speck]['F'])))
	
	var transposed_f = matrix_math.Transposed_Matrix(nub.particle_mechanics[speck]['F'])
	var inversed_trasposed_f = matrix_math.Inverse_Matrix(transposed_f)
	var f_coefficient = matrix_math.Subtract_Matrix(nub.particle_mechanics[speck]['F'],inversed_trasposed_f)
	
	pk_mu_coefficient =  matrix_math.Multiply_Matrix_by_Scalar(f_coefficient,mu,true)
	
	lambda_log_cofficient = lambda * log(nub.particle_mechanics[speck]['J']) # snapped(n * snapped(log(nub.particle_mechanics[speck]['J']),.001),.001)

	#pk_lambda_coefficient = matrix_math.Multiply_Matrix_by_Scalar(matrix_math.Inverse_Matrix(matrix_math.Transposed_Matrix(nub.particle_mechanics[speck]['F'])),lambda_log_cofficient,true)
	pk_lambda_coefficient =  matrix_math.Multiply_Matrix_by_Scalar(inversed_trasposed_f,lambda_log_cofficient,true)

	P = matrix_math.Add_Matrix(pk_mu_coefficient,pk_lambda_coefficient)

	#stress = 1 / J * P * F^T
	j_coefficient =  clamp(snapped((1.0 / nub.particle_mechanics[speck]['J']),.01),.001,(snapped((1.0 / nub.particle_mechanics[speck]['J']),.01))+1)
	stress_coefficient = matrix_math.Multiply_Matrix_by_Scalar(P,j_coefficient,true)
	#f_transposed_coefficient = matrix_math.Transposed_Matrix(nub.particle_mechanics[speck]['F'])
	resulted_stress = matrix_math.Multiply_Matrix(stress_coefficient,transposed_f)

	return resulted_stress


func Fixed_Corotated(speck,nub):
	### example "snow
	var resulted_stress
	var mu
	var lambda
	var mu_hardened
	var lambda_hardened
	var R
	var harden_coefficient = 10.0 # usually (3 - 10) : The hardening coefficient determines how fast the material breaks once it is plastic (larger = brittle, smaller = ductile)
	matrix_math = get_tree().get_root().get_node("Simulation/Matrix Math")
	
	mu = snapped((nub.youngs_modulus / snapped((2.0 * snapped((1.0+nub.poisson_ratio),.01)),.01) ),.01)
	lambda = snapped(snapped((nub.youngs_modulus * nub.poisson_ratio),.01) / ( (1.0 + nub.poisson_ratio) * clamp(snapped((1.0-snapped((2.0 * nub.poisson_ratio),.01)),.01),0.1,snapped((1.0-snapped((2.0 * nub.poisson_ratio),.01)),.01)+1) ),.01)
	
	var F_inverse = matrix_math.Inverse_Matrix(nub.particle_mechanics[speck]['F'])
	R = nub.particle_mechanics[speck]['stress'].duplicate(true)
	var e = harden_coefficient * (snapped((1.0 - nub.particle_mechanics[speck]['J']),.01))
	mu_hardened = mu * e
	lambda_hardened = lambda * e
	# P = 2 * mu_hardened * ( F - R ) + lambda_hardened * ( J - 1 ) * J * F ^-1 ; mpm course pg 20 eq 52
	
	var f_coefficient = matrix_math.Subtract_Matrix(nub.particle_mechanics[speck]['F'],R)
	var mu_harded_coefficient = snapped((2.0 * mu_hardened),.01)
	var mu_coefficient = matrix_math.Multiply_Matrix_by_Scalar(f_coefficient,mu_harded_coefficient,true)
	
	var lambda_harded_coefficient = snapped(snapped(lambda_hardened * (snapped((nub.particle_mechanics[speck]['J'] - 1.0),.01)),.01) * nub.particle_mechanics[speck]['J'],.01)
	var lambda_coefficient = matrix_math.Multiply_Matrix_by_Scalar(F_inverse,lambda_harded_coefficient,true)
	
	var P = matrix_math.Add_Matrix(mu_coefficient,lambda_coefficient)
	#stress = 1 / J * P * F^T
	var j_coefficient =  clamp(snapped((1.0 / nub.particle_mechanics[speck]['J']),.01),.001,(snapped((1.0 / nub.particle_mechanics[speck]['J']),.01))+1)
	var stress_coefficient = matrix_math.Multiply_Matrix_by_Scalar(P,j_coefficient,true)
	var f_transposed_coefficient = matrix_math.Transposed_Matrix(nub.particle_mechanics[speck]['F'])
	resulted_stress = matrix_math.Multiply_Matrix(stress_coefficient,f_transposed_coefficient)
	
	return resulted_stress


func Drucker_Prager_Elasticity(speck,nub):
	### examples : sand
	var resulted_stress
	var mu
	var lambda
	#var mu_hardened
	#var lambda_hardened
	#var harden_coefficient
	#var diagonalize_helper
	matrix_math = get_tree().get_root().get_node("Simulation/Matrix Math")
	var mu_cofficient
	var lambda_cofficient
	var logged_sigma
	var mu_lambda_sigma
	var P
	
	mu = snapped((nub.youngs_modulus / snapped((2.0 * snapped((1.0+nub.poisson_ratio),.01)),.01) ),.01)
	lambda = snapped(snapped((nub.youngs_modulus * nub.poisson_ratio),.01) / ( (1.0 + nub.poisson_ratio) * clamp(snapped((1.0-snapped((2.0 * nub.poisson_ratio),.01)),.01),0.1,snapped((1.0-snapped((2.0 * nub.poisson_ratio),.01)),.01)+1) ),.01)
	
	var FFtransposed = matrix_math.Multiply_Matrix(nub.particle_mechanics[speck]['F'],matrix_math.Transposed_Matrix(nub.particle_mechanics[speck]['F']))
	nub.particle_mechanics[speck]['U'] = matrix_math.Find_Eigenvectors(FFtransposed)
	var FtransposedF = matrix_math.Multiply_Matrix(nub.particle_mechanics[speck]['F'],matrix_math.Transposed_Matrix(nub.particle_mechanics[speck]['F']))
	nub.particle_mechanics[speck]['V'] = matrix_math.Find_Eigenvectors(FtransposedF)
	#Diagonalize SIGMA = PAP^-1... P = diagonalize_helper ,A = nub.particle_mechanics[speck]['F']
	var diagonalize_helper = matrix_math.Find_Eigenvectors(nub.particle_mechanics[speck]['F'])
	# 
	nub.particle_mechanics[speck]['Sigma'] = matrix_math.Multiply_Matrix(matrix_math.Multiply_Matrix(diagonalize_helper,nub.particle_mechanics[speck]['F']),matrix_math.Inverse_Matrix(diagonalize_helper))
			
	### if any is less than or equal to 0...
	if nub.particle_mechanics[speck]['Sigma'][0][0] <= 0.0 or nub.particle_mechanics[speck]['Sigma'][0][1] <= 0.0 or nub.particle_mechanics[speck]['Sigma'][1][0] <= 0.0 or nub.particle_mechanics[speck]['Sigma'][1][1] <= 0.0:
		if nub.particle_mechanics[speck]['Sigma'][0][0] <= 0.0 and nub.particle_mechanics[speck]['Sigma'][0][1] <= 0.0 and nub.particle_mechanics[speck]['Sigma'][1][0] <= 0.0 and nub.particle_mechanics[speck]['Sigma'][1][1] <= 0.0:
			logged_sigma =[[ snapped(log(.001),.001), snapped(log(.001),.001) ],[snapped(log(.001),.001), snapped(log(.001),.001) ]]
		
		elif nub.particle_mechanics[speck]['Sigma'][0][0] <= 0.0 and nub.particle_mechanics[speck]['Sigma'][0][1] <= 0.0 and nub.particle_mechanics[speck]['Sigma'][1][0] <= 0.0:
			logged_sigma =  [[ snapped(log(.001),.001), snapped(log(.001),.001) ],[ snapped(log(.001),.001), snapped(log(nub.particle_mechanics[speck]['Sigma'][1][1]),.001) ]]
		elif nub.particle_mechanics[speck]['Sigma'][0][0] <= 0.0 and nub.particle_mechanics[speck]['Sigma'][0][1] <= 0.0 and nub.particle_mechanics[speck]['Sigma'][1][1] <= 0.0:
			logged_sigma =  [[ snapped(log(.001),.001), snapped(log(.001),.001) ],[  snapped(log(nub.particle_mechanics[speck]['Sigma'][1][0]),.001), snapped(log(.001),.001) ]]
		elif nub.particle_mechanics[speck]['Sigma'][0][0] <= 0.0 and nub.particle_mechanics[speck]['Sigma'][1][0] <= 0.0 and nub.particle_mechanics[speck]['Sigma'][1][1] <= 0.0:
			logged_sigma =  [[ snapped(log(.001),.001),  snapped(log(nub.particle_mechanics[speck]['Sigma'][0][1]),.001) ],[snapped(log(.001),.001), snapped(log(.001),.001) ]]
		elif nub.particle_mechanics[speck]['Sigma'][0][1] <= 0.0 and nub.particle_mechanics[speck]['Sigma'][1][0] <= 0.0 and nub.particle_mechanics[speck]['Sigma'][1][1] <= 0.0:
			logged_sigma =  [[  snapped(log(nub.particle_mechanics[speck]['Sigma'][0][0]),.001), snapped(log(.001),.001) ],[ snapped(log(.001),.001), snapped(log(.001),.001) ]]
		
		elif nub.particle_mechanics[speck]['Sigma'][0][0] <= 0.0 and nub.particle_mechanics[speck]['Sigma'][0][1] <= 0.0:
			logged_sigma =  [[snapped(log(.001),.001), snapped(log(.001),.001) ],[ log(nub.particle_mechanics[speck]['Sigma'][1][0]),log(nub.particle_mechanics[speck]['Sigma'][1][1])]]
		elif nub.particle_mechanics[speck]['Sigma'][0][0] <= 0.0 and nub.particle_mechanics[speck]['Sigma'][1][0] <= 0.0:
			logged_sigma =  [[snapped(log(.001),.001), log(nub.particle_mechanics[speck]['Sigma'][0][1]) ],[ snapped(log(.001),.001),log(nub.particle_mechanics[speck]['Sigma'][1][1])]]
		elif nub.particle_mechanics[speck]['Sigma'][0][0] <= 0.0 and nub.particle_mechanics[speck]['Sigma'][1][1] <= 0.0:
			logged_sigma =  [[snapped(log(.001),.001),log(nub.particle_mechanics[speck]['Sigma'][0][1]) ],[ log(nub.particle_mechanics[speck]['Sigma'][1][0]), snapped(log(.001),.001)]]
		elif nub.particle_mechanics[speck]['Sigma'][0][1] <= 0.0 and nub.particle_mechanics[speck]['Sigma'][1][0] <= 0.0:
			logged_sigma =  [[log(nub.particle_mechanics[speck]['Sigma'][0][0]),snapped(log(.001),.001) ],[ snapped(log(.001),.001),log(nub.particle_mechanics[speck]['Sigma'][1][1])]]
		elif nub.particle_mechanics[speck]['Sigma'][0][1] <= 0.0 and nub.particle_mechanics[speck]['Sigma'][1][1] <= 0.0:
			logged_sigma =  [[log(nub.particle_mechanics[speck]['Sigma'][0][0]),snapped(log(.001),.001) ],[log(nub.particle_mechanics[speck]['Sigma'][1][0]),snapped(log(.001),.001)]]
		elif nub.particle_mechanics[speck]['Sigma'][1][0] <= 0.0 and nub.particle_mechanics[speck]['Sigma'][1][1] <= 0.0:
			logged_sigma =  [[log(nub.particle_mechanics[speck]['Sigma'][0][0]),log(nub.particle_mechanics[speck]['Sigma'][0][1]) ],[snapped(log(.001),.001),snapped(log(.001),.001)]]
	
	else:
		logged_sigma = [[ snapped(log(nub.particle_mechanics[speck]['Sigma'][0][0]),.001), snapped(log(nub.particle_mechanics[speck]['Sigma'][0][1]),.001) ],[ snapped(log(nub.particle_mechanics[speck]['Sigma'][1][0]),.001), snapped(log(nub.particle_mechanics[speck]['Sigma'][1][1]),.001) ]]
	var V_transposed = matrix_math.Transposed_Matrix(nub.particle_mechanics[speck]['V'])
	var SIGMA_inverse = matrix_math.Inverse_Matrix(nub.particle_mechanics[speck]['Sigma'])
	##print( logged_sigma,' logged sigma')
	#print(V_transposed, ' V^T')
	#print(SIGMA_inverse,' sigma^-1')
	
	### Plasticity Defining...
	### e_sigma = logged_sigma - (trace(logged_sigma)/d*I
	#d : is dimensions, I : Identity Matrix
	var trace_log_sigma = matrix_math.Trace(logged_sigma)
	var trace_identity = matrix_math.Multiply_Matrix_by_Scalar(nub.particle_mechanics[speck]['I'],(trace_log_sigma / 2.0),true)
	var e_sigma =  matrix_math.Subtract_Matrix(logged_sigma,trace_identity)
	
	### amount_of_plastic_deformation = deteminant(e_sigma) + ( d*lambda +2*mu / 2 *mu * trace(log_sigma) )*yield_surface_size
	
	var amount_of_plastic_deformation = matrix_math.Find_Determinant(e_sigma) + (snapped(((( snapped(2.0 * lambda,.01) )  + ( snapped(2.0 * mu,.01) )) / ( snapped((snapped(2.0 * mu,.01) ) * trace_log_sigma,.001) )) * nub.yield_surface,.001))
	
	### Hardening...
	# updated_hardening_state = hardening_state + alter_hardening_state (amount_of_plastic_deformation)
	# internal_coefficient_of_fiction = h_0 + ( h_1 * updated_hardening_state - h_3) * exp( h2 * updated_hardening_state)
	# nub.yield surface = square_root(2/3)*2*sin(internal_coefficient_of_fiction)/3-sin(internal_coefficient_of_fiction)
	
	var hardening_state = 0.0
	var updated_hardening_state = hardening_state + amount_of_plastic_deformation
	
	var h_0 = 35.0
	var h_1 = 0.0
	var h_2 = 0.2
	var h_3 = 10.0
	
	var internal_coefficient_of_fiction = h_0 + ((snapped(h_1 * updated_hardening_state,.001) - h_3) * snapped(exp((snapped(h_2 * updated_hardening_state,.001) )),.01) )
	nub.yield_surface = snapped(snapped(sqrt(snapped(2.0/3.0,.01)),.01) * ( snapped((2.0*snapped(sin(internal_coefficient_of_fiction),.001)),.01) / 3.0 - snapped(sin(internal_coefficient_of_fiction),.01) ),.01)
	
	
	
	# P = U * (2*mu*sigma^-1*log(sigma) + lambda * trace( log(sigma))*sigma^-1)*V^T
	mu_cofficient = matrix_math.Multiply_Matrix_by_Scalar(SIGMA_inverse,(2.0*mu),true)
	#var trace_log_sigma = get_tree().get_root().get_node("Simulation/Matrix Math").Trace(logged_sigma)
	lambda_cofficient = matrix_math.Multiply_Matrix_by_Scalar(SIGMA_inverse,(lambda*trace_log_sigma),true)
	
	mu_lambda_sigma = matrix_math.Add_Matrix(mu_cofficient,lambda_cofficient)
	#print(mu_lambda_sigma,' completed')
	var U_mu_lambda_sigma = matrix_math.Multiply_Matrix(nub.particle_mechanics[speck]['U'],mu_lambda_sigma)
	P = matrix_math.Multiply_Matrix(U_mu_lambda_sigma,V_transposed)
	#print(P,' P ')
	#stress = 1 / J * P * F^T
	var j_coefficient =  clamp(snapped((1.0 / nub.particle_mechanics[speck]['J']),.01),.001,(snapped((1.0 / nub.particle_mechanics[speck]['J']),.01))+1)
	var stress_coefficient = matrix_math.Multiply_Matrix_by_Scalar(P,j_coefficient,true)
	var f_transposed_coefficient = matrix_math.Transposed_Matrix(nub.particle_mechanics[speck]['F'])
	resulted_stress = matrix_math.Multiply_Matrix(stress_coefficient,f_transposed_coefficient)
	#print(resulted_stress,' stress')
	return resulted_stress
	
func Update_Plasticity(case,resurface,nub):
	### .
	matrix_math = get_tree().get_root().get_node("Simulation/Matrix Math")
	if case == 'snow':
		### snow plasicity...
		# ?? and ??s determine when the material starts breaking (larger = chunky, smaller = powdery)
			#critical_compression = 2.5 * 10^-2 : reference base...
		#critical_stretch = 7.5 * 10^-3 : reference base
		# Dry and powdery snow has smaller critical compression and stretch constants, while the opposite is true for 
		# wet and chunky snow. Icy snow has a higher hardening coefficient and Young???s modulus, with the opposite producing muddy snow.
		# material point method for snow simulation pg. 102.5...
		var critical_compression = 2.5 * pow(10.0,-2.0)
		var critical_stretch = 7.5 * pow(10.0,-3.0)
			
		#for row in range(0,len(nub)):
		#	for column in range(0,len(nub[row])):
					
		#		nub[row][column] = clampf(nub[row][column],1.0-critical_compression,1+critical_stretch)
		nub[0] = clampf(nub[0],1.0-critical_compression,1+critical_stretch)
		nub[1] = clampf(nub[1],1.0-critical_compression,1+critical_stretch)
		nub[2] = clampf(nub[2],1.0-critical_compression,1+critical_stretch)
		nub[3] = clampf(nub[3],1.0-critical_compression,1+critical_stretch)
	
	if case == 'sand':
		# if the the sand is wet.
		#if nub.particle_mechanics[speck]['I']s_wet == true:
			#pass
		### if any is less than or equal to 0...
		nub = matrix_math.Multiply_Matrix_by_Scalar(nub,resurface,true)
	
	
	return nub
