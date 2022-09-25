extends Node

var grid_nodes : Dictionary = {}
var gravity : Vector2 = Vector2(0.00,9.80)
#var gravity : Vector2 = Vector2(0.0,9.80)
#var gravity : Vector2 = Vector2(randf_range(-9.80,9.80),randf_range(-9.80,9.80))
var area_multiplier : float = 1.0 
var kernel_distance : Vector2
var kernel_x : float
var kernel_y : float
var flipped_kernel_distance : Vector2
var flipped_kernel_x : float
var flipped_kernel_y : float
var basis : String
var basis_coefficient : float
var basis_function_version : int
# handling particle merging...
var percentage_covered : float
var cycle_of_mergin_particles : Array = []
var merge_particles : bool
var list_of_particles_merged_particles : Array = []
# wall mechanics
var barriers : Dictionary
var wall_data : Dictionary
var window_center : Vector2
var wall_mass : float
var wall_velocity : Vector2 
var wall_momentum : Vector2
var velocity_from_wall_bounce : Vector2
# During ParticletoGrid...
var gradient_weight_term_1 
var gradient_weight_term_2 
var q_timed_volume_coefficient
var q_stressed_coefficient 
var q_gradient_weighted_stressed 
var q_gravity_weighted_stressed 
var mass_c_coefficient 
var flipped_mass_c_coefficient 
var affine_momentum_fuse : Vector2
var transfer_momentum 
var building_weights : Vector2
var weight_interpolation : float
var gradient_weight_interpolation : Vector2
var building_momentum
var forces : Vector2 = Vector2(0.0,0.0)
var Q_x : float 
var Q_y : float 
var Q_term_1
# During GridUpdate...
var standby_momentum_x : float
var standby_momentum_y : float
var placeholder_velocity_x : float
var placeholder_velocity_y : float
# During GridUpdate-Collisions...
var contact_with
var distance_between : Vector2
var particles_set_to_merge
# During GridtoParticles...
var other_particle_inverse_I_coefficient
var c_flipped_velocity_coefficient 
var c_weighted_velocity_matrix
var c_cell_coefficient
var c_inverse_I_coefficient
var c_cell_inverse_coefficient
var sum_of_c_weighted_velocity_matrix
var f_term_1
var f_term_2



func _on_Program_ready():
	### basis fuction setup...
	basis = "quadratic"
	basis_coefficient = 4.0
	basis_function_version = 2
	### setup up the barrier... 
	barriers['window outline'] = get_tree().get_root().get_node("Test Area").simulation_outline
	# allow each wall of the window to be a Constitution...
	#barriers['window outline']['top']['']
	#barriers['window outline']['right']
	#barriers['window outline']['bottom']
	#barriers['window outline']['left']
	### how the simulation window handles the borders...
	# 'None' the particles continue thru and keep going...
	# 'disappear' the particles that breach the window then disappear.
	# '2.0' the particles will bounce at max + 1(explode)
	# '1.0' the particles will bounce at max(perfectly elastic)
	# '0.5' the particles will bounce at partially(partially inelastic)
	# '0.0' the particles will sticky to the border(perfect inelastic)
	#barriers['coefficient of restitution'] = 'disappear'
	window_center = Vector2(ProjectSettings.get_setting('display/window/size/width')/2.0,ProjectSettings.get_setting('display/window/size/height')/2)
	for wall in barriers['window outline'].keys():
		barriers['window outline'][wall]['coefficient of restitution'] = 1.00
		#barriers['window outline'][wall]['coefficient of restitution'] = 0.50
		#barriers['window outline'][wall]['coefficient of restitution'] = 0.0
		barriers['window outline'][wall]['coefficient of static friction'] = 1.0
		barriers['window outline'][wall]['coefficient of kinetic friction'] = 1.0
		#barriers['window outline'][wall]['mass'] = 1000.00
		barriers['window outline'][wall]['mass'] = 1.00
		
		barriers['window outline'][wall]['velocity'] = Vector2(1.0,1.0)
		
		
func _notification(event):
	#print(event)
	if event == Node.NOTIFICATION_WM_CLOSE_REQUEST:
		#print("Exit")
		pass
	if event == Node.NOTIFICATION_WM_SIZE_CHANGED:
		#print('Window size change')
		pass
	
func Substance_Alignment():
	pass


func Weight_Interpolation(splines:String,version:int,kernel:Vector2,cell_size:float):
	#reset parameters
	building_weights = Vector2(0.0,0.0)
	### Weight Interpolation...
	if splines == "quadratic":
		#print('weight interpolation')
		if version == 1:
			pass
		elif version == 2:
			if kernel.x > -((3.0 * cell_size) / 2.0) and kernel.x < -( (1.0*cell_size) / 2.0 ) and  kernel.y > -((3.0 * cell_size) / 2.0) and kernel.y < - ( (1.0*cell_size) / 2.0 ):
			
				#print( -((3.0 * cell_size) / 2.0),' limit of weights')
				#print( -( (1.0*cell_size) / 2.0 ),' limit of weights' )
				
				#print(kernel,' is the kernel')
				#print(snapped(((1.0 * pow(kernel.x,2.0)) /  (2.0 * pow(cell_size,2.0))),.01 ) ,' search of component A')
				#print(snapped(((3.0 * kernel.x) / (2.0 * cell_size)),.01),' searching component b')
				#print(snapped((9.0/8.0),.001),' searching component c')
				
				building_weights.x =  snapped( ((1.0 * pow(kernel.x,2.0)) /  (2.0 * pow(cell_size,2.0))),.01 ) + snapped(((3.0 * kernel.x) / (2.0 * cell_size)),.01) + snapped((9.0/8.0),.001)
				building_weights.y =  snapped( ((1.0 * pow(kernel.y,2.0)) /  (2.0 * pow(cell_size,2.0))),.01 ) + snapped(((3.0 * kernel.y) / (2.0 * cell_size)),.01) + snapped((9.0/8.0),.001)
				#print("1")
				#print(building_weights.x,' finding weights of x')
				#print(building_weights.y,' finding weights of y')
				weight_interpolation = snapped(building_weights.x * building_weights.y,.01)
				
			elif kernel.x > -((1.0 * cell_size) / 2.0) and kernel.x <  ( (1.0*cell_size) / 2.0 ) and kernel.y > -((1.0 * cell_size) / 2.0) and kernel.y < ( (1.0*cell_size) / 2.0 ):
			
				#print( -((1.0 * cell_size) / 2.0),' limit of weights')
				#print( ( (1.0*cell_size) / 2.0 ),' limit of weights' )
				
				#print(kernel,' is the kernel')
				building_weights.x = snapped(((-1.0 * pow(kernel.x,2.0)) / (pow(cell_size,2.0))),.01) + snapped((3.0/4.0),.01)
				building_weights.y = snapped(((-1.0 * pow(kernel.y,2.0)) / (pow(cell_size,2.0))),.01) + snapped((3.0/4.0),.01)
				#print("2")
				#print(building_weights.x,' finding weights of x')
				#print(building_weights.y,' finding weights of y')
				weight_interpolation = snapped(building_weights.x * building_weights.y,.01)
				
			elif kernel.x > ((1.0 * cell_size) / 2.0) and kernel.x <  ( (3.0*cell_size) / 2.0 ) and  kernel.y > ((1 * cell_size) / 2) and kernel.y <  ( (3.0*cell_size) / 2.0 ):
			
				#print( ((1.0 * cell_size) / 2.0),' mid of weights')
				#print( ( (3.0*cell_size) / 2.0 ),' mid of weights' )
				
				#print(kernel,' is the kernel')
				building_weights.x =  snapped( ((1.0 * pow(kernel.x,2.0)) /  (2.0 * pow(cell_size,2.0))),.01 ) - snapped(((3.0 * kernel.x) / (2.0 * cell_size)),.01) + snapped((9.0/8.0),.001)
				building_weights.y =  snapped( ((1.0 * pow(kernel.y,2.0)) /  (2.0 * pow(cell_size,2.0))),.01 ) - snapped(((3.0 * kernel.y) / (2.0 * cell_size)),.01) + snapped((9.0/8.0),.001)
				#print("3")
				#print(kernel,' kernel ')
				
				#print(snapped( ((1.0 * pow(kernel.x,2.0)) /  (2.0 * pow(cell_size,2.0))),.01 ),' first coefficent')
				#print(snapped(((3.0 * kernel.x) / (2.0 * cell_size)),.01),' second coefficent')
				#print(snapped((9.0/8.0),.001),' third coefficent')
				#print(building_weights.x,' finding weights of x')
				#print(building_weights.y,' finding weights of y')
				weight_interpolation = snapped(building_weights.x * building_weights.y,.01)
			else:
				building_weights = Vector2(0.0,0.0)
				#print("4")
				#print(building_weights.x,' finding weights of x')
				#print(building_weights.y,' finding weights of y')
				weight_interpolation = snapped(building_weights.x * building_weights.y,.01)
				
	elif splines == "cubic":
		pass
	
	return weight_interpolation
	
func Simulate(time_passed,grid_domain,material,the_grid):
	### particle simulation...
	
	### reseting the grid data...
	for particle in the_grid.keys():
		
		the_grid[particle] = {'mass': get_tree().get_root().get_node("Test Area/Simulation/Substance").default_mass_of_particle,'velocity':get_tree().get_root().get_node("Test Area/Simulation/Substance").maintain_velocity,'momentum':Vector2(0.0,0.0),'force':Vector2(0.0,0.0)}
		
		#the_grid[particle] = {'mass': get_tree().get_root().get_node("Test Area/Simulation/Substance").default_mass_of_particle,'momentum':Vector2(0.0,0.0),'force':Vector2(0.0,0.0)}
		
		#the_grid[particle] = {'mass': get_tree().get_root().get_node("Test Area/Simulation/Substance").default_mass_of_particle,'velocity':Vector2(0.0,0.0),'momentum':Vector2(0.0,0.0),'force':Vector2(0.0,0.0)}
		
		#the_grid[particle] = {'mass': get_tree().get_root().get_node("Test Area/Simulation/Substance").default_mass_of_particle,'velocity':Vector2(1.0,1.0),'momentum':Vector2(0.0,0.0),'force':Vector2(0.0,0.0)}
		
		
		# establish the particle grid domain...
		#particle.surrounding_area = Rect2(Vector2(particle.position.x - ((grid_domain/2.0)*area_multiplier),particle.position.y - ((grid_domain/2.0)*area_multiplier)),Vector2(grid_domain*area_multiplier,grid_domain*area_multiplier))
		#area_multiplier = 1.0
		#particle.surrounding_area = Rect2(Vector2(particle.position.x - ((particle.get_node("shape").get_size().x/2.0)*area_multiplier),particle.position.y - ((particle.get_node("shape").get_size().y/2.0)*area_multiplier)),Vector2(particle.get_node("shape").get_size().x*area_multiplier,particle.get_node("shape").get_size().y*area_multiplier))
		
		#particle.surrounding_area = Rect2(Vector2(particle.position.x,particle.position.y),Vector2(particle.get_node("shape").get_size().x*area_multiplier,particle.get_node("shape").get_size().y*area_multiplier))
		particle.surrounding_area = Rect2(Vector2((particle.position.x+particle.get_node("shape").get_position().x),(particle.position.y+particle.get_node("shape").get_position().y)),Vector2(particle.get_node("shape").get_size().x*area_multiplier,particle.get_node("shape").get_size().y*area_multiplier))
		
		#print(particle.surrounding_area,' particle domain')
		particle.F = particle.I.duplicate(true)
		particle.J = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Find_Determinant(particle.F)
		particle.volume = particle.volume * particle.J
	#print(' ')
	#print('Sim Start')
	#print('Particles to Grid')
	# Particles to Grid:
	for particle in material:
		#print(particle,' ')
		
		#### the sum of aspects that is resetted...
		standby_momentum_x = 0
		standby_momentum_y = 0
		affine_momentum_fuse = Vector2(0.0,0.0)
		
		
		for other_particle in particle.within_range:
			
			#print(kernel_distance,' testing for floats')
			#print(flipped_kernel_distance,' the flipped test')
			### Weight Interpolation...
			Weight_Interpolation(basis,basis_function_version,particle.relation_to_domain[other_particle],grid_domain)
			
			#"""
			if other_particle.physical_state == 'solid':
				if other_particle.constitutive_model == 'hyperelastic':
					other_particle.stress = get_tree().get_root().get_node("Test Area/Simulation/Constitutive Models").Neo_Hookean(other_particle)
				if other_particle.constitutive_model == 'fixed_corated':
					other_particle.stress = get_tree().get_root().get_node("Test Area/Simulation/Constitutive Models").Fixed_Corotated(other_particle)
				if other_particle.constitutive_model == 'drucker_prager_elasticity':
					other_particle.stress = get_tree().get_root().get_node("Test Area/Simulation/Constitutive Models").Drucker_Prager_Elasticity(other_particle)
			elif other_particle.physical_state == 'liquid':
				if other_particle.constitutive_model == 'water':
					other_particle.stress = get_tree().get_root().get_node("Test Area/Simulation/Constitutive Models").Model_of_Water(other_particle)
					#print(other_particle.stress,' checking stress')
			elif other_particle.physical_state == 'gas':
				pass
			#"""
			
			
			### gradient-weight : weight_interpolation * D ^-1 * flipped_kernel
			#D **-1 = (4/pow(cell_size,2)) * other_particle.I^-1 
			gradient_weight_term_1 = snapped((snapped((weight_interpolation * basis_coefficient),.01) / pow(grid_domain,2.0)),.01)
			gradient_weight_term_2 =  get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Matrix_by_Scalar(particle.I,gradient_weight_term_1,true)
			gradient_weight_interpolation = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Matrix_by_Vector2_to_Vector2(gradient_weight_term_2,particle.domain_relation_to_particle[other_particle])
			
			
			#the_grid[particle]['momentum'] = Affine_Momentum_+_Particle_Force = sum of weight_interpolation * (other_particle_mass * other_particle_velocity + Q * flipped_kernel_distance)
			# Q = - sum of (time_passed * other_particle_volume * other_particle.stress * gradient_weight_interpolation + gravity + mass*other_particle.C)
			# gradient_weight_interpolation = D^-1 * flipped_kernel_distance 
			
			q_timed_volume_coefficient = time_passed * other_particle.volume
			q_stressed_coefficient = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Matrix_by_Scalar(other_particle.stress,q_timed_volume_coefficient,true)
			q_gradient_weighted_stressed = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Matrix_by_Vector2_to_Vector2(q_stressed_coefficient,gradient_weight_interpolation)
			q_gravity_weighted_stressed = (q_gradient_weighted_stressed + gravity)
			#q_gravity_weighted_stressed =  -(q_gradient_weighted_stressed + gravity)
			
			#print(q_gravity_weighted_stressed,' test q gravity_weighted_stress')
			mass_c_coefficient = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Matrix_by_Scalar(other_particle.C,other_particle.mass,true)
			flipped_mass_c_coefficient = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Matrix_by_Vector2_to_Vector2(mass_c_coefficient,particle.domain_relation_to_particle[other_particle])
			
			affine_momentum_fuse = -affine_momentum_fuse + (q_gravity_weighted_stressed + flipped_mass_c_coefficient)
			#affine_momentum_fuse = affine_momentum_fuse + (-1 * (q_gravity_weighted_stressed + flipped_mass_c_coefficient))
			
			#print(affine_momentum_fuse,' affine_momentum_fused ')
			transfer_momentum =  (other_particle.mass * other_particle.velocity)
			#print(transfer_momentum,' transfer_momentum')
			#print(flipped_mass_c_coefficient,' flipped_mass_c_coefficient')
			#standby_momentum_x = standby_momentum_x + weight_interpolation * ((other_particle.mass * other_particle.velocity) + (time_passed * other_particle.volume * other_particle.stress * gradient_weight_interpolation + gravity.x + (the_grid[other_particle]['mass']*other_particle.C * flipped_kernel_x)))
			#standby_momentum_x = standby_momentum_x + (snapped((weight_interpolation  *  (snapped((transfer_momentum.x + q_gravity_weighted_stressed.x + flipped_mass_c_coefficient.x),.01)) ),.01))
			#standby_momentum_y = standby_momentum_y + (snapped((weight_interpolation  *  (snapped((transfer_momentum.y + q_gravity_weighted_stressed.y + flipped_mass_c_coefficient.y),.01)) ),.01))
			standby_momentum_x = standby_momentum_x + (snapped((weight_interpolation  *  (snapped((transfer_momentum.x + affine_momentum_fuse.x),.01)) ),.01))
			standby_momentum_y = standby_momentum_y + (snapped((weight_interpolation  *  (snapped((transfer_momentum.y + affine_momentum_fuse.y),.01)) ),.01))
		
		the_grid[particle]['momentum'] = Vector2(standby_momentum_x,standby_momentum_y)
			
		#print(particle,' ',the_grid[particle]['momentum'],' momentum.')
		

	#Grid Update:
	#print(' ')
	#print('Grid Update')
	for node in the_grid.keys():
		if the_grid[node]['mass'] > 0.0:
			#print(the_grid[node]['mass'],' mass check')
			#the_grid[node]['momentum'] = the_grid[node]['momentum'].normalized() + (time_passed * ((the_grid[node]['mass'] * gravity) + the_grid[node]['force']))
			#the_grid[node]['momentum'] = the_grid[node]['momentum'] + (time_passed * ((the_grid[node]['mass'] * gravity) + the_grid[node]['force']))
			#the_grid[node]['momentum'] = the_grid[node]['momentum'] + (time_passed * ((the_grid[node]['mass'] * gravity) ))
			#standby_momentum_x = snapped((the_grid[node]['momentum'].x + snapped((time_passed * snapped((snapped((the_grid[node]['mass'] * gravity.x),.01) + the_grid[node]['force'].x),.01)),.01)),.01)
			#standby_momentum_y = snapped((the_grid[node]['momentum'].y + snapped((time_passed * snapped((snapped((the_grid[node]['mass'] * gravity.y),.01) + the_grid[node]['force'].y),.01)),.01)),.01)
			#the_grid[node]['momentum'] = Vector2(standby_momentum_x,standby_momentum_y)
			
			#the_grid[node]['velocity'] = (the_grid[node]['momentum'] / the_grid[node]['mass'])
			placeholder_velocity_x = snapped(the_grid[node]['momentum'].x / the_grid[node]['mass'],.01)
			placeholder_velocity_y = snapped(the_grid[node]['momentum'].y / the_grid[node]['mass'],.01)
			the_grid[node]['velocity'] = Vector2(placeholder_velocity_x,placeholder_velocity_y)
			#print(node,' ',the_grid[node]['velocity'],' velocity ')
		else:
			###...
			##print('the particle has no mass.')
			pass
	#"""
	### Collision Detection...
	### how the particles interacts....
	# with other particles,
	# window outline...
	#
	#print(' ')
	#print('collision detect')
	#print(' ')
	for particle in material:
		
		### if the particle has contacted the window outline....
		if particle.surrounding_area.position.x < barriers['window outline']['left']['outline']:
			### the particle left the window out line at the sides...
			#print('breached left')
			the_grid[particle]['velocity'] = get_tree().get_root().get_node("Test Area/Simulation/Particle Interaction").Collision_with_Walls(material,barriers,the_grid,grid_domain,gravity)
		
		elif particle.surrounding_area.end.x >= barriers['window outline']['right']['outline']:
			### the particle left the window out line at the sides...
			#print('breached right')
			the_grid[particle]['velocity'] = get_tree().get_root().get_node("Test Area/Simulation/Particle Interaction").Collision_with_Walls(material,barriers,the_grid,grid_domain,gravity)
		
		elif particle.surrounding_area.position.y < barriers['window outline']['top']['outline']:
			### the particle left the window out line at the sides...
			#print('breached top')
			the_grid[particle]['velocity'] = get_tree().get_root().get_node("Test Area/Simulation/Particle Interaction").Collision_with_Walls(material,barriers,the_grid,grid_domain,gravity)
		
		elif particle.surrounding_area.end.y >= barriers['window outline']['bottom']['outline']:
			### the particle left the window out line at the sides...
			#print('breached bottom')
			the_grid[particle]['velocity'] = get_tree().get_root().get_node("Test Area/Simulation/Particle Interaction").Collision_with_Walls(material,barriers,the_grid,grid_domain,gravity)
		else:
			### the particle is within the window outline...
			pass
		
		#print(the_grid[particle]['velocity'],' checking updated grid velocity after collision')
		#print(' ')
		#print(particle,' checking particle...')
		### if the particle has contacted another particle....
		for other_particle in material:
			
			# particle contact and/or relation to every other particle...
			#contact_with = (particle.get_position() - other_particle.surrounding_area.get_center()) / grid_domain
			contact_with = (particle.surrounding_area.get_center() - other_particle.surrounding_area.get_center()) / grid_domain
			#print("checking if colliding with ",other_particle)
			#print(' ')
			#print(particle,' distance from other particle ',other_particle,' is ',contact_with)#*grid_domain)
			#print(grid_domain/grid_domain,' and ',-(grid_domain/grid_domain),' are the particle barrier')
			#print(particle,' relation to ',other_particle,' is ',particle.relation_to_domain[other_particle],' check before update')
				
			#kernel_distance = (particle.get_position() - other_particle.surrounding_area.get_center())# / cell_size
			kernel_x = snapped(particle.get_position().x - other_particle.surrounding_area.get_center().x,.01)
			kernel_y = snapped(particle.get_position().y - other_particle.surrounding_area.get_center().y,.01)
			kernel_distance = Vector2(kernel_x,kernel_y)
			#flipped_kernel_distance = (other_particle.surrounding_area.get_center() - particle.get_position())# / cell_size
			flipped_kernel_x = snapped(other_particle.surrounding_area.get_center().x - particle.get_position().x,.01)
			flipped_kernel_y = snapped(other_particle.surrounding_area.get_center().y - particle.get_position().y,.01)
			flipped_kernel_distance = Vector2(flipped_kernel_x,flipped_kernel_y)
				
			#print(kernel_distance,' is the kernel between the particles ',particle,' ',other_particle)
			particle.relation_to_domain[other_particle] = kernel_distance
			particle.domain_relation_to_particle[other_particle] = flipped_kernel_distance
			#particle.relation_to_domain[other_particle] = flipped_kernel_distance
			#particle.domain_relation_to_particle[other_particle] = kernel_distance
			#print(particle,' relation to ',other_particle,' is ',particle.relation_to_domain[other_particle],' check after update')
			"""
			if particle != other_particle:
				if particle.surrounding_area.intersects(other_particle.surrounding_area):
					### come into some contact...
				
					#print(particle,' at position ',particle.position,' in contact with ',other_particle,' at position ',other_particle.surrounding_area,' quadrant')
					if particle.within_range.has(other_particle):
						### the particle is already in contact with the other...
						pass
					else:
						### recognized to be within the domain range...
						particle.within_range.append(other_particle)
						
					### particles interaction with other particles...
					get_tree().get_root().get_node("Test Area/Simulation/Particle Interaction").Collision_between_Other_Particles(particle,other_particle,the_grid,grid_domain)
							
				else:
					### particle is not in contact with each other...
					particle.within_range.erase(other_particle)
					other_particle.within_range.erase(particle)
			else:
				### can't collide with itself...
				pass
			#"""
			#"""
			if  contact_with.x == 0 and contact_with.y == 0:
				### always in contact with itself...
				pass
				
			elif abs(contact_with.x) < grid_domain/grid_domain and abs(contact_with.y) < grid_domain/grid_domain:
				### the particles come into contact with each other...
			
			#elif contact_with.x > -(grid_domain/grid_domain) or contact_with.x < grid_domain/grid_domain and contact_with.y > -(grid_domain/grid_domain) or contact_with.y < grid_domain/grid_domain:
				### come into some contact...
				
				#print(particle,' at position ',particle.position,' in contact with ',other_particle,' at position ',other_particle.surrounding_area,' quadrant')
				if particle.within_range.has(other_particle):
					### the particle is already in contact with the other...
					pass
				else:
					### recognized to be within the domain range...
					particle.within_range.append(other_particle)
				
				
				#var test_intersectional = (particle.surrounding_area.intersection(other_particle.surrounding_area)).get_center()
				##print(" ")
				#print(test_intersectional,' formed when the particles collide')
				#print(test_intersectional.get_center(),' is the center of the instersectional')
				### particles interaction with other particles...
				the_grid[particle]['velocity'] = get_tree().get_root().get_node("Test Area/Simulation/Particle Interaction").Collision_between_Other_Particles(particle,other_particle,the_grid,grid_domain)
				#print(the_grid_field[particle])
			else:
				### doesn't come into contact with another... 
				#print(particle,' not in contact with this particle ',other_particle)
				#reset...
				particle.within_range = [particle]
				other_particle.within_range = [other_particle]
			#"""
		
	# Grid to Particle...
	#print(' ')
	#print('Grid to Particle')
	# resetting of the particles...
	for particle in material:
		particle.mass = get_tree().get_root().get_node("Test Area/Simulation/Substance").default_mass_of_particle
		#particle.velocity = get_tree().get_root().get_node("Test Area/Simulation/Substance").maintain_velocity
		#particle.velocity = Vector2(0.0,0.0)
		particle.velocity = Vector2(1.0,1.0)
	
	for particle in material:
		#print(particle,' ')
		#### sum of aspects reset...
		sum_of_c_weighted_velocity_matrix = [[0.0,0.0],[0.0,0.0]]
		
		for other_particle in particle.within_range:
			
			### Weight Interpolation...
			Weight_Interpolation(basis,basis_function_version,particle.relation_to_domain[other_particle],grid_domain)
			#print(weight_interpolation,' wip')
			#MLS MPM 
			#particle.velocity = sum of (the_grid[velocity] * weight_interpolation)
			#particle.C = D^1 * sum of the_grid['velocity'] * flipped_kernel_distance^T * weight_interpolation
			# D^1 = (4/pow(cell_size,2)) * particle.I^-1 
			
			particle.velocity.x = snapped((particle.velocity.x + snapped((the_grid[other_particle]['velocity'].x * weight_interpolation),.01)),.01)
			particle.velocity.y = snapped((particle.velocity.y + snapped((the_grid[other_particle]['velocity'].y * weight_interpolation),.01)),.01)
			
			c_flipped_velocity_coefficient = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Vector2_by_Vector2_to_Matrix(the_grid[other_particle]['velocity'],false,particle.domain_relation_to_particle[other_particle],true)
			c_weighted_velocity_matrix = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Matrix_by_Scalar(c_flipped_velocity_coefficient,weight_interpolation,true)
			sum_of_c_weighted_velocity_matrix = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Add_Matrix(sum_of_c_weighted_velocity_matrix,c_weighted_velocity_matrix)
		
		c_cell_coefficient =  snapped((basis_coefficient / pow(grid_domain,2.0)),.01)
		c_inverse_I_coefficient = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Inverse_Matrix(particle.I)
		c_cell_inverse_coefficient = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Matrix_by_Scalar(c_inverse_I_coefficient,c_cell_coefficient,true)
				
		particle.C = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Matrix(c_cell_inverse_coefficient,sum_of_c_weighted_velocity_matrix)
			
		#print(particle.C,' C updated')
		###particle position update...
		#print(particle.position,' position before')
		#print(time_passed,' time passing')
		#print(particle.velocity,' velocity after')
		
		#particle.position = particle.get_position() + (time_passed * particle.velocity)
		
		particle.position.x = snapped((particle.position.x + snapped((time_passed * particle.velocity.x),.01)  ),.01)
		#if particle.position.x > ProjectSettings.get_setting('display/window/size/width'):
		#	particle.position.x = clampf(particle.position.x,(0.0+(grid_domain/2.0)),(ProjectSettings.get_setting('display/window/size/width')-(grid_domain/2.0) ))
		#particle.position.x = clampf(particle.position.x,(0.0+(grid_domain/2.0)),(ProjectSettings.get_setting('display/window/size/width')-(grid_domain/2.0) ))
		##particle.position.x = clampf(particle.position.x,(0.0),(ProjectSettings.get_setting('display/window/size/width')-(grid_domain*2) ))
		particle.position.y = snapped((particle.position.y + (time_passed * particle.velocity.y)),.01)
		#particle.position.y = clampf(particle.position.y,(0.0+(grid_domain/2.0)),(ProjectSettings.get_setting('display/window/size/height')-(grid_domain/2.0) ))
		#particle.position.y = clampf(particle.position.y,(0.0),(ProjectSettings.get_setting('display/window/size/height')-(grid_domain*2) ))
		
		
		#print(particle.position,' position after')
		###deformation update...
		#particle.F = (particle.I + time_passed * particle.C ) * particle.F
		#print(particle.F,' deformation gradient before')
		#print(particle.C,' checking C')
		#print(time_passed,' checking time passed')
		f_term_1 = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Matrix_by_Scalar(particle.C,time_passed,true)
		#print(f_term_1,' time_passed * particle.C')
		f_term_2 =  get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Add_Matrix(particle.I,f_term_1)
		#print(f_term_2,' particle.I + f term 1')
		particle.F = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Matrix(f_term_2,particle.F)
		### Updating Plasticity...
		particle = get_tree().get_root().get_node("Test Area/Simulation/Constitutive Models").Update_Plasticity(particle)
