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
	
	mu = snapped((nub.youngs_modulus / snapped((2.0 * snapped((1.0+nub.poisson_ratio),.01)),.01) ),.01)
	#lambda = snapped(snapped((nub.youngs_modulus * nub.poisson_ratio),.01) / ( (1.0 + nub.poisson_ratio) * clamp(snapped((1.0-snapped((2.0 * nub.poisson_ratio),.01)),.01),0.1,snapped((1.0-snapped((2.0 * nub.poisson_ratio),.01)),.01)+1) ),.01)
	lambda = (nub.youngs_modulus * nub.poisson_ratio) / (1 + nub.poisson_ratio) * (1-2*nub.poisson_ratio)
	
	#'''
	#  Piola-Kirchoff stress = u*F + (n * In(J) - u )*F^-T
	# note ln() = log()
	# or
	# P = u*(F - F^-T) + lambda*log (J)*F^-T ; mpm course pg. 19,equation 48
	#f_coefficient =  get_tree().get_root().get_node("Simulation/Matrix Math").Subtract_Matrix(nub.mechanics[speck]['F'],get_tree().get_root().get_node("Simulation/Matrix Math").Inverse_Matrix(get_tree().get_root().get_node("Simulation/Matrix Math").Transposed_Matrix(nub.mechanics[speck]['F'])))
	#print(mu,' mu check')
	#print(lambda,' lambda check')
	var transposed_f = get_tree().get_root().get_node('Simulation/Matrix Math').Transposed_Matrix(nub.mechanics[speck]['F'])
	var inversed_trasposed_f = get_tree().get_root().get_node('Simulation/Matrix Math').Inverse_Matrix(transposed_f)
	#print(transposed_f,' transposed_f check')
	#print(inversed_trasposed_f,' inversed_trasposed_f check')
	f_coefficient = get_tree().get_root().get_node('Simulation/Matrix Math').Subtract_Matrix(nub.mechanics[speck]['F'],inversed_trasposed_f)
	#print(f_coefficient,' f_coefficient check')
	pk_mu_coefficient =  get_tree().get_root().get_node('Simulation/Matrix Math').Multiply_Matrix_by_Scalar(f_coefficient,mu,true)
	
	#print(nub.mechanics[speck]['J'],' J check')
	
	lambda_log_cofficient = lambda * log(nub.mechanics[speck]['J']) # snapped(n * snapped(log(nub.mechanics[speck]['J']),.001),.001)
	#print(inversed_trasposed_f,' inversed_trasposed_f check')
	#print(lambda_log_cofficient,' lambda_log_cofficient check')

	#pk_lambda_coefficient = get_tree().get_root().get_node("Simulation/Matrix Math").Multiply_Matrix_by_Scalar(get_tree().get_root().get_node("Simulation/Matrix Math").Inverse_Matrix(get_tree().get_root().get_node("Simulation/Matrix Math").Transposed_Matrix(nub.mechanics[speck]['F'])),lambda_log_cofficient,true)
	pk_lambda_coefficient =  get_tree().get_root().get_node('Simulation/Matrix Math').Multiply_Matrix_by_Scalar(inversed_trasposed_f,lambda_log_cofficient,true)
	#print(pk_mu_coefficient,' pk_mu_coefficient check')
	#print(pk_lambda_coefficient,' pk_lambda_coefficient check')

	P = get_tree().get_root().get_node('Simulation/Matrix Math').Add_Matrix(pk_mu_coefficient,pk_lambda_coefficient)
	#print(P,' P check')
	#print(nub.mechanics[speck]['J'],' J check')
	#stress = (1 / J) * P * F^T ; mpm course pg 18 eq 38
	#P = elastic energy density function
	#j_coefficient =  clamp(snapped((1.0 / nub.mechanics[speck]['J']),.01),1.0,(snapped((1.0 / nub.mechanics[speck]['J']),.01))+1)
	j_coefficient = 1.0 / nub.mechanics[speck]['J']
	#print(j_coefficient,' j_coefficient check')
	stress_coefficient = get_tree().get_root().get_node('Simulation/Matrix Math').Multiply_Matrix_by_Scalar(P,j_coefficient,true)
	#f_transposed_coefficient = get_tree().get_root().get_node("Simulation/Matrix Math").Transposed_Matrix(nub.mechanics[speck]['F'])
	resulted_stress = get_tree().get_root().get_node('Simulation/Matrix Math').Multiply_Matrix(stress_coefficient,transposed_f)
	#print(resulted_stress,' results stress check')
	# inf , nan check...
	for location in range(0,(len(resulted_stress))):
		if is_inf(resulted_stress[location]) == true or is_nan(resulted_stress[location]) == true:
			resulted_stress[location] = 0.0
			#if location == 1 or location == 3:
				#resulted_stress[location] = 0.0
		
	return resulted_stress



func Update_Plasticity(nub,_name,matter):
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
		nub.mechanics[name]['F'][0] = clampf(nub.mechanics[name]['F'][0],1.0-critical_compression,1+critical_stretch)
		nub.mechanics[name]['F'][1] = clampf(nub.mechanics[name]['F'][1],1.0-critical_compression,1+critical_stretch)
		nub.mechanics[name]['F'][2] = clampf(nub.mechanics[name]['F'][2],1.0-critical_compression,1+critical_stretch)
		nub.mechanics[name]['F'][3] = clampf(nub.mechanics[name]['F'][3],1.0-critical_compression,1+critical_stretch)
	
	if matter == 'sand':
		# if the the sand is wet.
		#if nub.mechanics[speck]['I']s_wet == true:
			#pass
		#var case
		
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
		
		var e_sigma = [log(nub.mechanics[name]['sigma'][0]),0,0,log(nub.mechanics[name]['sigma'][3])]
		#print(e_sigma,' e sigma check')
		var dimension = 2 # 2 or 3...
		
		var trace_coefficient = get_tree().get_root().get_node("Simulation/Matrix Math").Trace(e_sigma) / dimension
		var identity_trace_coefficient = get_tree().get_root().get_node("Simulation/Matrix Math").Multiply_Matrix_by_Scalar(nub.mechanics[name]['I'],trace_coefficient,true)
		var updated_e_sigma =  get_tree().get_root().get_node("Simulation/Matrix Math").Subtract_Matrix(e_sigma,identity_trace_coefficient)
		
		var determinant_e_sigma = get_tree().get_root().get_node("Simulation/Matrix Math").Find_Determinant(updated_e_sigma)

		if determinant_e_sigma == 0 or get_tree().get_root().get_node("Simulation/Matrix Math").Trace(e_sigma) == 0:
			### case II...
			#print('case II')
			update_harden = get_tree().get_root().get_node("Simulation/Matrix Math").Find_Determinant(e_sigma)
			project = nub.mechanics[name]['I']
		else:
			friction_angle = h0 + (h1*update_harden - h3 ) * exp(-h2*update_harden)
			#print(friction_angle,' friction_angle check')
			#alpha = sqrt(2/3) * 2 * rad_to_deg(sin(friction_angle)) / 3 - rad_to_deg(sin(friction_angle))
			alpha = sqrt(2.0/3.0) * ((2.0 * sin(friction_angle)) / (3.0 - sin(friction_angle)))
			#print(alpha,' alpha check')
			
			var mu_lambda_coefficient = ((dimension * nub.mechanics[name]['lambda']) + (2 *  nub.mechanics[name]['mu'])) / (2 *  nub.mechanics[name]['mu'])
			var e_sigma_coefficient = mu_lambda_coefficient * get_tree().get_root().get_node("Simulation/Matrix Math").Trace(e_sigma)
			#var determinant_e_sigma = get_tree().get_root().get_node("Simulation/Matrix Math").Find_Determinant(updated_e_sigma)
			var amount_of_plastic_deformation = determinant_e_sigma + e_sigma_coefficient * alpha
			
			if amount_of_plastic_deformation <= 0:
				### case 1
				#print('case I')
				update_harden = 0
				project = nub.mechanics[name]['sigma']
			else:
				var updated_e_sigma_cofficient = get_tree().get_root().get_node("Simulation/Matrix Math").Divide_Matrix_by_Scalar(updated_e_sigma,determinant_e_sigma,true)
				var H = get_tree().get_root().get_node("Simulation/Matrix Math").Multiply_Matrix_by_Scalar(updated_e_sigma_cofficient,amount_of_plastic_deformation,true)
				update_harden = amount_of_plastic_deformation
				project = exp(H)
				
		hardening_state = hardening_state + update_harden
		
		var U_project_coefficient =  get_tree().get_root().get_node("Simulation/Matrix Math").Multiply_Matrix(nub.mechanics[name]['U'],project)
		var V_transposed = get_tree().get_root().get_node("Simulation/Matrix Math").Transposed_Matrix(nub.mechanics[name]['V'])
		nub.mechanics[name]['F'] = get_tree().get_root().get_node("Simulation/Matrix Math").Multiply_Matrix(U_project_coefficient,V_transposed)
		### 
		
	# inf,nan check
	if is_inf(nub.mechanics[name]['F'][0]) == true or is_nan(nub.mechanics[name]['F'][0]) == true:
		nub.mechanics[name]['F'][0] = 1.0
	if is_inf(nub.mechanics[name]['F'][1]) == true or is_nan(nub.mechanics[name]['F'][1]) == true:
		nub.mechanics[name]['F'][1] = 0.0
	if is_inf(nub.mechanics[name]['F'][2]) == true or is_nan(nub.mechanics[name]['F'][2]) == true:
		nub.mechanics[name]['F'][2] = 0.0
	if is_inf(nub.mechanics[name]['F'][3]) == true or is_nan(nub.mechanics[name]['F'][3]) == true:
		nub.mechanics[name]['F'][3] = 1.0
	
		
	return nub



