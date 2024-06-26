extends Node

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


func Neo_Hookean(speck,nub):
	### Neo-Hookean model nonlinear hyperelastic models for predicting large deformations of elastic materials. mpm pg. 19
	### example "rubber
	var resulted_stress
	var mu
	var lambda
	
	#var pk_stress
	
	#mu = snapped((nub.youngs_modulus / snapped((2.0 * snapped((1.0+nub.poisson_ratio),.01)),.01) ),.01)
	mu = nub.youngs_modulus / (2 * (1.0+nub.poisson_ratio))
	#lambda = snapped(snapped((nub.youngs_modulus * nub.poisson_ratio),.01) / ( (1.0 + nub.poisson_ratio) * clamp(snapped((1.0-snapped((2.0 * nub.poisson_ratio),.01)),.01),0.1,snapped((1.0-snapped((2.0 * nub.poisson_ratio),.01)),.01)+1) ),.01)
	lambda = (nub.youngs_modulus * nub.poisson_ratio) / ((1.0 + nub.poisson_ratio) * (1.0-(2.0*nub.poisson_ratio)))
	#'''
	#  Piola-Kirchoff stress = u*F + (n * In(J) - u )*F^-T
	# note ln() = log()
	# or
	# P = u*(F - F^-T) + lambda*log (J)*F^-T ; mpm course pg. 19,equation 48
	#var transposed_f = get_tree().get_root().get_node('Simulation/Matrix Math').Transposed_Matrix(nub.mechanics[speck]['F'])
	var transposed_f = get_tree().get_root().get_node('Simulation/Matrix Math').Transposed_Matrix(nub.F)
	var inversed_trasposed_f = get_tree().get_root().get_node('Simulation/Matrix Math').Inverse_Matrix(transposed_f)
	
	#f_coefficient = get_tree().get_root().get_node('Simulation/Matrix Math').Subtract_Matrix(nub.mechanics[speck]['F'],inversed_trasposed_f)
	f_coefficient = get_tree().get_root().get_node('Simulation/Matrix Math').Subtract_Matrix(nub.F,inversed_trasposed_f)
	
	pk_mu_coefficient =  get_tree().get_root().get_node('Simulation/Matrix Math').Multiply_Matrix_by_Scalar(f_coefficient,mu,true)
	
	lambda_log_cofficient = lambda * log(nub.J) # snapped(n * snapped(log(nub.mechanics[speck]['J']),.001),.001)
	
	#pk_lambda_coefficient = get_tree().get_root().get_node("Simulation/Matrix Math").Multiply_Matrix_by_Scalar(get_tree().get_root().get_node("Simulation/Matrix Math").Inverse_Matrix(get_tree().get_root().get_node("Simulation/Matrix Math").Transposed_Matrix(nub.mechanics[speck]['F'])),lambda_log_cofficient,true)
	pk_lambda_coefficient =  get_tree().get_root().get_node('Simulation/Matrix Math').Multiply_Matrix_by_Scalar(inversed_trasposed_f,lambda_log_cofficient,true)

	P = get_tree().get_root().get_node('Simulation/Matrix Math').Add_Matrix(pk_mu_coefficient,pk_lambda_coefficient)
	#stress = (1 / J) * P * F^T ; mpm course pg 18 eq 38
	#P = elastic energy density function
	#j_coefficient =  clamp(snapped((1.0 / nub.mechanics[speck]['J']),.01),1.0,(snapped((1.0 / nub.mechanics[speck]['J']),.01))+1)
	j_coefficient = 1.0 / nub.J
	stress_coefficient = get_tree().get_root().get_node('Simulation/Matrix Math').Multiply_Matrix_by_Scalar(P,j_coefficient,true)
	#f_transposed_coefficient = get_tree().get_root().get_node("Simulation/Matrix Math").Transposed_Matrix(nub.mechanics[speck]['F'])
	f_transposed_coefficient = get_tree().get_root().get_node("Simulation/Matrix Math").Transposed_Matrix(nub.F)
	resulted_stress = get_tree().get_root().get_node('Simulation/Matrix Math').Multiply_Matrix(stress_coefficient,transposed_f)
	# inf , nan check...
	for location in range(0,(len(resulted_stress))):
		if is_inf(resulted_stress[location]) == true or is_nan(resulted_stress[location]) == true:
			resulted_stress[location] = 0.0
			#if location == 1 or location == 3:
				#resulted_stress[location] = 0.0
		
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
	
	
	var f_t = get_tree().get_root().get_node('Simulation/Matrix Math').Transposed_Matrix(nub.F)
	#print(f_t,' f_t check')
	var u = get_tree().get_root().get_node('Simulation/Matrix Math').Multiply_Matrix(f_t,nub.F)
	var dia_u = [sqrt(u[0]),0,0,sqrt(u[3])]
	#print(dia_u,' dia_u check')
	var inv_u =  get_tree().get_root().get_node('Simulation/Matrix Math').Inverse_Matrix(dia_u)
	#print(inv_u,' inv_u check')
	#var R = get_tree().get_root().get_node('Simulation/Matrix Math').Multiply_Matrix(material.particle_mechanics[particle]['F'],inv_u)
		
	#print(r,' r check')
	
	mu = snapped((nub.youngs_modulus / snapped((2.0 * snapped((1.0+nub.poisson_ratio),.01)),.01) ),.01)
	lambda = snapped(snapped((nub.youngs_modulus * nub.poisson_ratio),.01) / ( (1.0 + nub.poisson_ratio) * clamp(snapped((1.0-snapped((2.0 * nub.poisson_ratio),.01)),.01),0.1,snapped((1.0-snapped((2.0 * nub.poisson_ratio),.01)),.01)+1) ),.01)
	
	mu_hardened = mu * exp( harden_coefficient * (snapped((1.0 - nub.J),.01)))
	lambda_hardened = lambda * exp( harden_coefficient * (snapped((1.0 - nub.J),.01)))
	
	#var F_inverse = get_tree().get_root().get_node('Simulation/Matrix Math').Inverse_Matrix(nub.particle_mechanics[speck]['F'])
	R = get_tree().get_root().get_node('Simulation/Matrix Math').Multiply_Matrix(nub.F,inv_u)
	
	#print(mu_hardened,' check mu_hardened')
	#print(lambda_hardened,' check lambda_hardened')
	
	#"""
	# energy density function
	# P = 2 * mu_hardened * ( F - R ) + lambda_hardened * ( J - 1 ) * J * F ^-T ; mpm course pg 20 eq 52
	
	f_coefficient = get_tree().get_root().get_node('Simulation/Matrix Math').Subtract_Matrix(nub.F,R)
	
	var mu_harded_coefficient = snapped((2.0 * mu_hardened),.01)
	var mu_coefficient = get_tree().get_root().get_node('Simulation/Matrix Math').Multiply_Matrix_by_Scalar(f_coefficient,mu_harded_coefficient,true)
	
	var lambda_harded_coefficient = snapped(snapped(lambda_hardened * (snapped((nub.J - 1.0),.01)),.01) * nub.J,.01)
	var lambda_coefficient = get_tree().get_root().get_node('Simulation/Matrix Math').Multiply_Matrix_by_Scalar(f_t,lambda_harded_coefficient,true)
	
	P = get_tree().get_root().get_node('Simulation/Matrix Math').Add_Matrix(mu_coefficient,lambda_coefficient)
	#"""
	
	#stress = (1 / J) * P * F^T ; mpm course pg 18 eq 38
	#P = elastic energy density function
	j_coefficient =  clamp(snapped((1.0 / nub.J),.01),.001,(snapped((1.0 / nub.J),.01))+1)
	stress_coefficient = get_tree().get_root().get_node('Simulation/Matrix Math').Multiply_Matrix_by_Scalar(P,j_coefficient,true)
	f_transposed_coefficient = get_tree().get_root().get_node('Simulation/Matrix Math').Transposed_Matrix(nub.F)
	resulted_stress = get_tree().get_root().get_node('Simulation/Matrix Math').Multiply_Matrix(stress_coefficient,f_transposed_coefficient)
	
	if is_inf(resulted_stress[0]) == true or is_nan(resulted_stress[0]) == true:
		resulted_stress[0] = 1.0
	if is_inf(resulted_stress[1]) == true or is_nan(resulted_stress[1]) == true:
		resulted_stress[1] = 1.0
	if is_inf(resulted_stress[2]) == true or is_nan(resulted_stress[2]) == true:
		resulted_stress[2] = 1.0
	if is_inf(resulted_stress[3]) == true or is_nan(resulted_stress[3]) == true:
		resulted_stress[3] = 1.0
		
		
	#print(resulted_stress,' check fixed-corotated')
	
	return resulted_stress





func Model_of_Water(speck,nub):
	var resulted_stress
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
	#var deteminant_coefficient = 1.0 / pow(nub.mechanics[speck]['J'],y)
	var deteminant_coefficient = 1.0 / pow(nub.J,y)
	#print(deteminant_coefficient,' deteminant coefficient')
	#print(nub.J,' nub.J')
	var p_coefficient = 2.1 * (deteminant_coefficient - 1.0 )
	#print(p_coefficient,' p')
	resulted_stress = get_tree().get_root().get_node('Simulation/Matrix Math').Multiply_Matrix_by_Scalar(nub.I,-p_coefficient,true)
	#print(resulted_stress)
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
	
	var mu_cofficient
	var lambda_cofficient
	#var logged_sigma
	#var mu_lambda_sigma
	#var P
	
	mu = snapped((nub.youngs_modulus / snapped((2.0 * snapped((1.0+nub.poisson_ratio),.01)),.01) ),.01)
	lambda = snapped(snapped((nub.youngs_modulus * nub.poisson_ratio),.01) / ( (1.0 + nub.poisson_ratio) * clamp(snapped((1.0-snapped((2.0 * nub.poisson_ratio),.01)),.01),0.1,snapped((1.0-snapped((2.0 * nub.poisson_ratio),.01)),.01)+1) ),.01)
	
	# polar decomposition...
	#print(nub.particle_mechanics[speck]['F'],' nub.particle_mechanics[speck][F] check')
	var f_transposed = get_tree().get_root().get_node('Simulation/Matrix Math').Transposed_Matrix(nub.F)
	#print(f_transposed,' f_transposed check')
	nub.U = get_tree().get_root().get_node('Simulation/Matrix Math').Multiply_Matrix(f_transposed,nub.F)
	#print(nub.particle_mechanics[speck]['U'],' U check')
	nub.sigma = [sqrt(nub.U[0]),0,0,sqrt(nub.U[3])]
	#print(nub.particle_mechanics[speck]['sigma'],' sigma check')
	nub.V = get_tree().get_root().get_node('Simulation/Matrix Math').Multiply_Matrix(nub.F,f_transposed)
	#print(nub.particle_mechanics[speck]['V'],' V check')
	
	
	
	# solving P = U * (2*mu*sigma^-1*log(sigma) + lambda * trace( log(sigma))*sigma^-1)*V^T
	var inverse_sigma =  get_tree().get_root().get_node('Simulation/Matrix Math').Inverse_Matrix(nub.sigma)
	#print(inverse_sigma,' inverse_sigma check')
	var logarithm_sigma = [log(nub.sigma[0]),0,0,log(nub.sigma[3])]
	#print(logarithm_sigma,' logarithm_sigma check')
	var V_transposed = get_tree().get_root().get_node('Simulation/Matrix Math').Transposed_Matrix(nub.V)
	#print(V_transposed,' V_transposed check')
	
	var mu_2x = 2 * mu
	#print(mu_2x,' mu_2x check')
	var mu_inversed = get_tree().get_root().get_node('Simulation/Matrix Math').Multiply_Matrix_by_Scalar(inverse_sigma,mu_2x,true)
	#print(mu_inversed,' mu_inversed check')
	mu_cofficient = get_tree().get_root().get_node('Simulation/Matrix Math').Multiply_Matrix(mu_inversed,logarithm_sigma)
	#print(mu_cofficient,' mu_cofficient check')
	var traced_lambda = lambda * get_tree().get_root().get_node('Simulation/Matrix Math').Trace(logarithm_sigma)
	lambda_cofficient = get_tree().get_root().get_node('Simulation/Matrix Math').Multiply_Matrix_by_Scalar(inverse_sigma,traced_lambda,true)
	#print(lambda_cofficient,' lambda_cofficient check')
	
	var mu_lambda_coefficient = get_tree().get_root().get_node('Simulation/Matrix Math').Add_Matrix(mu_cofficient,lambda_cofficient)
	#print(mu_lambda_coefficient,' mu_lambda_coefficient check')
	var p_component = get_tree().get_root().get_node('Simulation/Matrix Math').Multiply_Matrix(nub.U,mu_lambda_coefficient)
	#print(p_component,' p_component check')
	P = get_tree().get_root().get_node('Simulation/Matrix Math').Multiply_Matrix(p_component,V_transposed)
	#print(P,' Drucker_Prager_Elasticity check')
	
	#stress = (1 / J) * P * F^T ; mpm course pg 18 eq 38
	#P = elastic energy density function
	j_coefficient =  clamp(snapped((1.0 / nub.J),.01),.001,(snapped((1.0 / nub.J),.01))+1)
	stress_coefficient = get_tree().get_root().get_node('Simulation/Matrix Math').Multiply_Matrix_by_Scalar(P,j_coefficient,true)
	f_transposed_coefficient = get_tree().get_root().get_node('Simulation/Matrix Math').Transposed_Matrix(nub.F)
	resulted_stress = get_tree().get_root().get_node('Simulation/Matrix Math').Multiply_Matrix(stress_coefficient,f_transposed_coefficient)
	
	if is_inf(resulted_stress[0]) == true or is_nan(resulted_stress[0]) == true:
		resulted_stress[0] = 1.0
	if is_inf(resulted_stress[1]) == true or is_nan(resulted_stress[1]) == true:
		resulted_stress[1] = 0.0
	if is_inf(resulted_stress[2]) == true or is_nan(resulted_stress[2]) == true:
		resulted_stress[2] = 0.0
	if is_inf(resulted_stress[3]) == true or is_nan(resulted_stress[3]) == true:
		resulted_stress[3] = 1.0
	
	return resulted_stress
	




func Update_Plasticity(nub,name,matter):
	### .
	
	if matter == 'snow':
		### snow plasicity...
		# θ and θs determine when the material starts breaking (larger = chunky, smaller = powdery)
			#critical_compression = 2.5 * 10^-2 : reference base...
		#critical_stretch = 7.5 * 10^-3 : reference base
		# Dry and powdery snow has smaller critical compression and stretch constants, while the opposite is true for 
		# wet and chunky snow. Icy snow has a higher hardening coefficient and Young’s modulus, with the opposite producing muddy snow.
		# material point method for snow simulation pg. 102.5...
		var critical_compression = 2.5 * pow(10.0,-2.0)
		var critical_stretch = 7.5 * pow(10.0,-3.0)
			
		#for row in range(0,len(nub)):
		#	for column in range(0,len(nub[row])):
					
		#		nub[row][column] = clampf(nub[row][column],1.0-critical_compression,1+critical_stretch)
		nub.F[0] = clampf(nub.F[0],1.0-critical_compression,1+critical_stretch)
		nub.F[1] = clampf(nub.F[1],1.0-critical_compression,1+critical_stretch)
		nub.F[2] = clampf(nub.F[2],1.0-critical_compression,1+critical_stretch)
		nub.F[3] = clampf(nub.F[3],1.0-critical_compression,1+critical_stretch)
	
	if matter == 'sand':
		# if the the sand is wet.
		#if nub.mechanics[speck]['I']s_wet == true:
			#pass
		#var case
		var mu
		var lambda
		
		var alpha 
		var hardening_state = 0
		var update_harden
		var project
		#friction angle
		var friction_angle 
		### # hardening parameters
		# Feasible hardening parameters satisfy h0 > h3 ≥ 0 and h1 , h2 ≥ 0::
		var h0 = 40
		var h1 = 10
		var h2 = 0
		var h3 = 20
		
		mu = snapped((nub.youngs_modulus / snapped((2.0 * snapped((1.0+nub.poisson_ratio),.01)),.01) ),.01)
		lambda = snapped(snapped((nub.youngs_modulus * nub.poisson_ratio),.01) / ( (1.0 + nub.poisson_ratio) * clamp(snapped((1.0-snapped((2.0 * nub.poisson_ratio),.01)),.01),0.1,snapped((1.0-snapped((2.0 * nub.poisson_ratio),.01)),.01)+1) ),.01)

		var e_sigma = [log(nub.sigma[0]),0,0,log(nub.sigma[3])]
		#print(e_sigma,' e sigma check')
		var dimension = 2 # 2 or 3...
		update_harden = 0 
		
		var trace_coefficient = get_tree().get_root().get_node("Simulation/Matrix Math").Trace(e_sigma) / dimension
		var identity_trace_coefficient = get_tree().get_root().get_node("Simulation/Matrix Math").Multiply_Matrix_by_Scalar(nub.I,trace_coefficient,true)
		var updated_e_sigma =  get_tree().get_root().get_node("Simulation/Matrix Math").Subtract_Matrix(e_sigma,identity_trace_coefficient)
		
		var determinant_e_sigma = get_tree().get_root().get_node("Simulation/Matrix Math").Find_Determinant(updated_e_sigma)

		if determinant_e_sigma == 0 or get_tree().get_root().get_node("Simulation/Matrix Math").Trace(e_sigma) == 0:
			### case II...
			#print('case II')
			update_harden = get_tree().get_root().get_node("Simulation/Matrix Math").Find_Determinant(e_sigma)
			project = nub.I
		else:
			
			friction_angle = h0 + (h1*update_harden - h3 ) * exp(-h2*update_harden)
			#print(friction_angle,' friction_angle check')
			#alpha = sqrt(2/3) * 2 * rad_to_deg(sin(friction_angle)) / 3 - rad_to_deg(sin(friction_angle))
			alpha = sqrt(2.0/3.0) * ((2.0 * sin(friction_angle)) / (3.0 - sin(friction_angle)))
			#print(alpha,' alpha check')
			
			var mu_lambda_coefficient = ((dimension * lambda) + (2 *  mu)) / (2 *  mu)
			var e_sigma_coefficient = mu_lambda_coefficient * get_tree().get_root().get_node("Simulation/Matrix Math").Trace(e_sigma)
			#var determinant_e_sigma = get_tree().get_root().get_node("Simulation/Matrix Math").Find_Determinant(updated_e_sigma)
			var amount_of_plastic_deformation = determinant_e_sigma + e_sigma_coefficient * alpha
			
			if amount_of_plastic_deformation <= 0:
				### case 1
				#print('case I')
				update_harden = 0
				project = nub.sigma
			else:
				#case 3
				var updated_e_sigma_cofficient = get_tree().get_root().get_node("Simulation/Matrix Math").Divide_Matrix_by_Scalar(updated_e_sigma,determinant_e_sigma,true)
				var H = get_tree().get_root().get_node("Simulation/Matrix Math").Multiply_Matrix_by_Scalar(updated_e_sigma_cofficient,amount_of_plastic_deformation,true)
				update_harden = amount_of_plastic_deformation
				project = H#exp(H)
				
		hardening_state = hardening_state + update_harden
		
		var U_project_coefficient =  get_tree().get_root().get_node("Simulation/Matrix Math").Multiply_Matrix(nub.U,project)
		var V_transposed = get_tree().get_root().get_node("Simulation/Matrix Math").Transposed_Matrix(nub.V)
		nub.F = get_tree().get_root().get_node("Simulation/Matrix Math").Multiply_Matrix(U_project_coefficient,V_transposed)
		### 
		
	# inf,nan check
	if is_inf(nub.F[0]) == true or is_nan(nub.F[0]) == true:
		nub.F[0] = 1.0
	if is_inf(nub.F[1]) == true or is_nan(nub.F[1]) == true:
		nub.F[1] = 0.0
	if is_inf(nub.F[2]) == true or is_nan(nub.F[2]) == true:
		nub.F[2] = 0.0
	if is_inf(nub.F[3]) == true or is_nan(nub.F[3]) == true:
		nub.F[3] = 1.0
	
		
	return nub



