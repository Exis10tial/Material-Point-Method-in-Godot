extends Node

var grid_nodes : Dictionary = {}
var gravity : Vector2 = Vector2(0.00,9.80)
#var gravity : Vector2 = Vector2(0.00,-9.80)
#var gravity : Vector2 = Vector2(9.80,0.0)
#var gravity : Vector2 = Vector2(-9.80,0.0)
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
	window_center = Vector2(ProjectSettings.get_setting('display/window/size/viewport_width')/2.0,ProjectSettings.get_setting('display/window/size/viewport_height')/2)
	for wall in barriers['window outline'].keys():
		barriers['window outline'][wall]['coefficient of restitution'] = 1.00
		#barriers['window outline'][wall]['coefficient of restitution'] = 0.50
		#barriers['window outline'][wall]['coefficient of restitution'] = 0.0
		barriers['window outline'][wall]['coefficient of static friction'] = 1.0
		barriers['window outline'][wall]['coefficient of kinetic friction'] = 1.0
		barriers['window outline'][wall]['mass'] = 1000.00
		#barriers['window outline'][wall]['mass'] = 1.00
		barriers['window outline']['top']['outline'] = 0.0
		barriers['window outline']['right']['outline'] = ProjectSettings.get_setting('display/window/size/viewport_width')
		barriers['window outline']['bottom']['outline'] = ProjectSettings.get_setting('display/window/size/viewport_height')
		barriers['window outline']['left']['outline'] = 0.0
		barriers['window outline'][wall]['velocity'] = Vector2(1.0,1.0)
		
		
func _notification(event):
	#print(event)
	if event == Node.NOTIFICATION_WM_CLOSE_REQUEST:
		#print("Exit")
		pass
	if event == Node.NOTIFICATION_WM_MOUSE_EXIT:
		pass
	#if event == Node.NOTIFICATION_WM_MOUSE_ENTER:
		#
		#ProjectSettings.set_setting('display/window/size/width',get_tree().get_root().get_size().x)
		#ProjectSettings.set_setting('display/window/size/width',get_tree().get_root().get_size().y)


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
	
func Simulate(time_passed:float,material:Object,the_grid:Dictionary):
	### particle simulation...
	print(" ")
	print('Cycle Start')
	### reseting the grid data...
	for particle in the_grid.keys():
		
		the_grid[particle] = {'mass': material.default_mass_of_substance,'velocity':Vector2(0.0,0.0),'momentum':Vector2(0.0,0.0),'force':Vector2(0.0,0.0)}
		#print(the_grid[particle], ' grid')
		material.particle_mechanics[particle]['F'] = material.particle_mechanics[particle]['I'].duplicate(true)
		material.particle_mechanics[particle]['J']= get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Find_Determinant(material.particle_mechanics[particle]['F'])
		material.particle_mechanics[particle]['volume'] = material.particle_mechanics[particle].volume * material.particle_mechanics[particle]['J']
		#print(material.particle_mechanics[particle]['F'], ' material.particle_mechanics[particle][F]')
		#print(material.particle_mechanics[particle]['J'], ' material.particle_mechanics[particle][J]')
		#print(material.particle_mechanics[particle]['volume'], ' material.particle_mechanics[particle][volume]')
	
	# Particles to Grid:
	for particle in material.particle_mechanics.keys():

		#### the sum of aspects that is resetted...
		standby_momentum_x = 0
		standby_momentum_y = 0
		affine_momentum_fuse = Vector2(0.0,0.0)
		
		
		for other_particle in material.particle_mechanics[particle]['within_range']:
			### Weight Interpolation...
			Weight_Interpolation(basis,basis_function_version,material.particle_mechanics[particle]['relation_to_domain'][other_particle],material.cell_size)
			print(weight_interpolation,' weight_interpolation')
			#"""
			if material.physical_state == 'solid':
				if material.constitutive_model == 'hyperelastic':
					material.particle_mechanics[other_particle].stress = get_tree().get_root().get_node("Test Area/Simulation/Constitutive Models").Neo_Hookean(particle,material)
				if material.constitutive_model == 'fixed_corated':
					material.particle_mechanics[other_particle].stress = get_tree().get_root().get_node("Test Area/Simulation/Constitutive Models").Fixed_Corotated(particle,material)
				if material.constitutive_model == 'drucker_prager_elasticity':
					material.particle_mechanics[other_particle].stress = get_tree().get_root().get_node("Test Area/Simulation/Constitutive Models").Drucker_Prager_Elasticity(particle,material)
			elif material.physical_state == 'liquid':
				if material.constitutive_model == 'water':
					material.particle_mechanics[other_particle].stress = get_tree().get_root().get_node("Test Area/Simulation/Constitutive Models").Model_of_Water(particle,material)
					#print(other_particle.stress,' checking stress')
			elif material.physical_state == 'gas':
				pass
			print(material.particle_mechanics[other_particle].stress,' checking stress')
			#"""
			### gradient-weight : weight_interpolation * D ^-1 * flipped_kernel
			#D **-1 = (4/pow(cell_size,2)) * other_particle.I^-1 
			gradient_weight_term_1 = snapped((snapped((weight_interpolation * basis_coefficient),.01) / pow(material.cell_size,2.0)),.01)
			gradient_weight_term_2 =  get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Matrix_by_Scalar(material.particle_mechanics[particle]['I'],gradient_weight_term_1,true)
			gradient_weight_interpolation = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Matrix_by_Vector2_to_Vector2(gradient_weight_term_2,material.particle_mechanics[particle]['domain_relation_to_substance'][other_particle])
			
			q_timed_volume_coefficient = time_passed * material.particle_mechanics[other_particle]['volume']
			q_stressed_coefficient = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Matrix_by_Scalar(material.particle_mechanics[other_particle]['stress'],q_timed_volume_coefficient,true)
			q_gradient_weighted_stressed = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Matrix_by_Vector2_to_Vector2(q_stressed_coefficient,gradient_weight_interpolation)
			q_gravity_weighted_stressed = (q_gradient_weighted_stressed + gravity)
			
			mass_c_coefficient = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Matrix_by_Scalar(material.particle_mechanics[other_particle]['C'],material.particle_mechanics[other_particle]['mass'],true)
			flipped_mass_c_coefficient = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Matrix_by_Vector2_to_Vector2(mass_c_coefficient,material.particle_mechanics[particle]['domain_relation_to_substance'][other_particle])
			
			affine_momentum_fuse = -affine_momentum_fuse + (q_gravity_weighted_stressed + flipped_mass_c_coefficient)
			#affine_momentum_fuse = affine_momentum_fuse + (-1 * (q_gravity_weighted_stressed + flipped_mass_c_coefficient))
			
			print(affine_momentum_fuse,' affine_momentum_fused ')
			transfer_momentum =  (material.particle_mechanics[other_particle]['mass'] * material.particle_mechanics[other_particle]['velocity'])
			print(transfer_momentum,' transfer_momentum')
			#print(flipped_mass_c_coefficient,' flipped_mass_c_coefficient')
			#standby_momentum_x = standby_momentum_x + weight_interpolation * ((other_particle.mass * other_particle.velocity) + (time_passed * other_particle.volume * other_particle.stress * gradient_weight_interpolation + gravity.x + (the_grid[other_particle]['mass']*other_particle.C * flipped_kernel_x)))
			#standby_momentum_x = standby_momentum_x + (snapped((weight_interpolation  *  (snapped((transfer_momentum.x + q_gravity_weighted_stressed.x + flipped_mass_c_coefficient.x),.01)) ),.01))
			#standby_momentum_y = standby_momentum_y + (snapped((weight_interpolation  *  (snapped((transfer_momentum.y + q_gravity_weighted_stressed.y + flipped_mass_c_coefficient.y),.01)) ),.01))
			standby_momentum_x = standby_momentum_x + (snapped((weight_interpolation  *  (snapped((transfer_momentum.x + affine_momentum_fuse.x),.01)) ),.01))
			standby_momentum_y = standby_momentum_y + (snapped((weight_interpolation  *  (snapped((transfer_momentum.y + affine_momentum_fuse.y),.01)) ),.01))
		
		the_grid[particle]['momentum'] = Vector2(standby_momentum_x,standby_momentum_y)
			
		print(particle,' ',the_grid[particle]['momentum'],' momentum.')
		
	#Grid Update:
	#print(' ')
	print('Grid Update')
	for particle in the_grid.keys():
		if the_grid[particle]['mass'] > 0.0:
			print(the_grid[particle]['mass'],' mass check')
			#the_grid[node]['momentum'] = the_grid[node]['momentum'].normalized() + (time_passed * ((the_grid[node]['mass'] * gravity) + the_grid[node]['force']))
			#the_grid[node]['momentum'] = the_grid[node]['momentum'] + (time_passed * ((the_grid[node]['mass'] * gravity) + the_grid[node]['force']))
			#the_grid[node]['momentum'] = the_grid[node]['momentum'] + (time_passed * ((the_grid[node]['mass'] * gravity) ))
			#standby_momentum_x = snapped((the_grid[node]['momentum'].x + snapped((time_passed * snapped((snapped((the_grid[node]['mass'] * gravity.x),.01) + the_grid[node]['force'].x),.01)),.01)),.01)
			#standby_momentum_y = snapped((the_grid[node]['momentum'].y + snapped((time_passed * snapped((snapped((the_grid[node]['mass'] * gravity.y),.01) + the_grid[node]['force'].y),.01)),.01)),.01)
			
			#standby_momentum_x = snapped((the_grid[node]['momentum'].x + snapped((time_passed * snapped((snapped((the_grid[node]['mass'] * gravity.x),.01)),.01)),.01)),.01)
			#standby_momentum_y = snapped((the_grid[node]['momentum'].y + snapped((time_passed * snapped((snapped((the_grid[node]['mass'] * gravity.y),.01)),.01)),.01)),.01)
			#the_grid[node]['momentum'] = Vector2(standby_momentum_x,standby_momentum_y)
			
			#the_grid[node]['velocity'] = (the_grid[node]['momentum'] / the_grid[node]['mass'])
			print(the_grid[particle]['momentum'],' the_grid[particle][momentum]')
			print(the_grid[particle]['mass'],' the_grid[particle][mass]')
			placeholder_velocity_x = snapped(the_grid[particle]['momentum'].x / the_grid[particle]['mass'],.01)
			placeholder_velocity_y = snapped(the_grid[particle]['momentum'].y / the_grid[particle]['mass'],.01)
			the_grid[particle]['velocity'] = Vector2(placeholder_velocity_x,placeholder_velocity_y)
			
			### Forces from other objects.
			#the_grid[particle]['velocity']  = the_grid[particle]['velocity']  + Vector2(randf_range(-10.0,10.0),randf_range(-10.0,10.0))
			#the_grid[particle]['velocity']  = the_grid[particle]['velocity'] + Vector2(0.0,0.0)
			
			print(particle,' ',the_grid[particle]['velocity'],' grid velocity ')
		else:
			###...
			##print('the particle has no mass.')
			pass
	#"""
	### Collision Detection...
	### how the particles interacts....
	for particle in material.particle_lineation.keys():
		#print(particle,' particle')
		#print(particle.surrounding_area,' particle.surrounding_area')
		#print(particle.surrounding_area.end.x,' particle.surrounding_area.end.x')
		### if the particle has contacted the window outline....
		#if particle.surrounding_area.end.x >= barriers['window outline']['right']['outline'] and  particle.surrounding_area.end.y >= barriers['window outline']['bottom']['outline']:
		#	the_grid[particle]['velocity'] = Vector2(0.0,0.0)
		if material.particle_lineation[particle].position.y < barriers['window outline']['top']['outline']:
			### the particle left the window out line at the sides...
			print('breached top')
			the_grid[particle]['velocity'] = get_tree().get_root().get_node("Test Area/Simulation/Particle Interaction").Collision_with_Walls('top',material,particle,material.particle_lineation[particle],barriers,the_grid,material.cell_size)
		
		if material.particle_lineation[particle].end.x >= barriers['window outline']['right']['outline']:
			### the particle left the window out line at the sides...
			#print(particle.surrounding_area.end.x,' particle.surrounding_area.end.x')
			#print(barriers['window outline']['right']['outline']," barriers['window outline']['right']['outline']")
			print('breached right')
			the_grid[particle]['velocity'] = get_tree().get_root().get_node("Test Area/Simulation/Particle Interaction").Collision_with_Walls('right',material,particle,material.particle_lineation[particle],barriers,the_grid,material.cell_size)
			
		if material.particle_lineation[particle].end.y >= barriers['window outline']['bottom']['outline']:
			### the particle left the window out line at the sides...
			print('breached bottom')
			the_grid[particle]['velocity'] = get_tree().get_root().get_node("Test Area/Simulation/Particle Interaction").Collision_with_Walls('bottom',material,particle,material.particle_lineation[particle],barriers,the_grid,material.cell_size)
			#print(particle,' ',the_grid[particle]['velocity'],' velocity ')
			
		if material.particle_lineation[particle].position.x < barriers['window outline']['left']['outline']:
			### the particle left the window out line at the sides...
			print('breached left')
			the_grid[particle]['velocity'] = get_tree().get_root().get_node("Test Area/Simulation/Particle Interaction").Collision_with_Walls('left',material,material.particle_lineation[particle],barriers,the_grid,material.cell_size)
		else:
			### the particle is within the window outline...
			pass

		### if the particle has contacted another particle....
		for other_particle in material.particle_lineation.keys():
			
			kernel_x = snapped(material.particle_lineation[particle].get_center().x - material.particle_lineation[other_particle].get_center().x,.01)
			kernel_y = snapped(material.particle_lineation[particle].get_center().y - material.particle_lineation[other_particle].get_center().y,.01)
			kernel_distance = Vector2(kernel_x,kernel_y)
			#flipped_kernel_distance = (other_particle.surrounding_area.get_center() - particle.get_position())# / cell_size
			flipped_kernel_x = snapped(material.particle_lineation[other_particle].get_center().x - material.particle_lineation[particle].get_center().x,.01)
			flipped_kernel_y = snapped(material.particle_lineation[other_particle].get_center().y - material.particle_lineation[particle].get_center().y,.01)
			flipped_kernel_distance = Vector2(flipped_kernel_x,flipped_kernel_y)
				
			#print(kernel_distance,' is the kernel between the particles ',particle,' ',other_particle)
			material.particle_mechanics[particle]['relation_to_domain'][other_particle] = kernel_distance
			material.particle_mechanics[particle]['domain_relation_to_substance'][other_particle] = flipped_kernel_distance
			#particle.relation_to_domain[other_particle] = flipped_kernel_distance
			#particle.domain_relation_to_particle[other_particle] = kernel_distance
			#print(particle,' relation to ',other_particle,' is ',particle.relation_to_domain[other_particle],' check after update')
			contact_with = kernel_distance / sqrt( pow(float(material.particle_lineation[particle].size.x),2.0) + pow(float(material.particle_lineation[particle].size.y),2.0) )
			#print(contact_with, 'contact with')
			#"""
			#print(range(-1.0,2.0))
			if contact_with.x == 0 and contact_with.y == 0:
				### always in contact with itself...
				pass
				
			#elif abs(contact_with.x) < 1.0 abs(contact_with.x) < 1.0and abs(contact_with.y) < 1.0:
			elif contact_with.x < -1.0 or  contact_with.x > 1.0 or contact_with.y < -1.0 or  contact_with.y > 1.0:
				### the particles come into contact with each other...
				#print('touch other particle')
			#elif contact_with.x > -(grid_domain/grid_domain) or contact_with.x < grid_domain/grid_domain and contact_with.y > -(grid_domain/grid_domain) or contact_with.y < grid_domain/grid_domain:
				### come into some contact...
				
				#print(particle,' at position ',particle.position,' in contact with ',other_particle,' at position ',other_particle.surrounding_area,' quadrant')
				if material.particle_mechanics[particle]['within_range'].has(other_particle):
					### the particle is already in contact with the other...
					pass
				else:
					### recognized to be within the domain range...
					material.particle_mechanics[particle]['within_range'].append(other_particle)
				
				
				#var test_intersectional = (particle.surrounding_area.intersection(other_particle.surrounding_area)).get_center()
				##print(" ")
				#print('formed when the particles collide')
				#print(test_intersectional.get_center(),' is the center of the instersectional')
				### particles interaction with other particles...
				the_grid[particle]['velocity'] = get_tree().get_root().get_node("Test Area/Simulation/Particle Interaction").Collision_between_Other_Particles(material[particle],material[other_particle],the_grid,material[particle].cell_size)
				#print(the_grid_field[particle])
			else:
				### doesn't come into contact with another... 
				#print(particle,' not in contact with this particle ',other_particle)
				#reset...
				material.particle_lineation[particle]['within_range'] = [particle]
				material.particle_lineation[other_particle]['within_range'] = [other_particle]
				#print('not touch other particle')
			#"""
		
	# Grid to Particle...
	for particle in material.particle_mechanics.keys():
		material.particle_mechanics[particle]['mass'] = material.default_mass_of_substance
		#particle.velocity = get_tree().get_root().get_node("Test Area/Simulation/Substance").maintain_velocity
		material.particle_mechanics[particle]['velocity'] = Vector2(0.0,0.0)
		#particle.velocity = Vector2(1.0,1.0)
	
	for particle in material.particle_mechanics.keys():
		#print(particle,' ')
		#### sum of aspects reset...
		sum_of_c_weighted_velocity_matrix = [[0.0,0.0],[0.0,0.0]]
		#print(particle,' particle')
		#print(particle.within_range,' particle.within_range')
		for other_particle in material.particle_mechanics[particle]['within_range']:
			
			### Weight Interpolation...
			Weight_Interpolation(basis,basis_function_version,material.particle_mechanics[particle]['relation_to_domain'][other_particle],material.cell_size)
			#print(weight_interpolation,' wip')
			#MLS MPM 
			#particle.velocity = sum of (the_grid[velocity] * weight_interpolation)
			#particle.C = D^1 * sum of the_grid['velocity'] * flipped_kernel_distance^T * weight_interpolation
			# D^1 = (4/pow(cell_size,2)) * particle.I^-1
			print(the_grid[particle],' the_grid[particle] before')
			#print(weight_interpolation,' weight_interpolation ')
			material.particle_mechanics[particle]['velocity'].x = snapped(material.particle_mechanics[particle]['velocity'].x + (snapped((the_grid[other_particle]['velocity'].x * weight_interpolation),.01)),.01)
			material.particle_mechanics[particle]['velocity'].y = snapped(material.particle_mechanics[particle]['velocity'].y + (snapped((the_grid[other_particle]['velocity'].y * weight_interpolation),.01)),.01)
			
			#print()
			c_flipped_velocity_coefficient = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Vector2_by_Vector2_to_Matrix(the_grid[other_particle]['velocity'],false,material.particle_mechanics[particle]['domain_relation_to_substance'][other_particle],true)
			c_weighted_velocity_matrix = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Matrix_by_Scalar(c_flipped_velocity_coefficient,weight_interpolation,true)
			sum_of_c_weighted_velocity_matrix = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Add_Matrix(sum_of_c_weighted_velocity_matrix,c_weighted_velocity_matrix)
		
		c_cell_coefficient =  snapped((basis_coefficient / pow(material.cell_size,2.0)),.01)
		c_inverse_I_coefficient = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Inverse_Matrix(material.particle_mechanics[particle]['I'])
		c_cell_inverse_coefficient = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Matrix_by_Scalar(c_inverse_I_coefficient,c_cell_coefficient,true)
				
		material.particle_mechanics[particle]['C'] = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Matrix(c_cell_inverse_coefficient,sum_of_c_weighted_velocity_matrix)
		
		print(material.particle_mechanics[particle]['velocity'],' material.particle_mechanics[particle][velocity]')
		print(material.particle_lineation[particle],' particle positon before')
		print(time_passed,' time_passed')
		print(material.particle_mechanics[particle]['velocity'],' material.particle_mechanics[particle][velocity]')
		###particle position update...
		#particle.position = particle.get_position() + (time_passed * particle.velocity)
		material.particle_lineation[particle].position.x = snapped((material.particle_lineation[particle].position.x + snapped((time_passed * material.particle_mechanics[particle]['velocity'].x),.01)  ),.01)
		material.particle_lineation[particle].position.y = snapped((material.particle_lineation[particle].position.y + snapped((time_passed * material.particle_mechanics[particle]['velocity'].y),.01)  ),.01)
		
		#particle.position.x = clampf(particle.position.x,1.0,ProjectSettings.get_setting('display/window/size/viewport_width')-1)
		#particle.position.y = clampf(particle.position.y,1.0,ProjectSettings.get_setting('display/window/size/viewport_height')-1)
		print(material.particle_lineation[particle],' particle positon')
		
		#print(particle.position,' position after')
		###deformation update...
		#particle.F = (particle.I + time_passed * particle.C ) * particle.F
		#print(particle.F,' deformation gradient before')
		#print(particle.C,' checking C')
		#print(time_passed,' checking time passed')
		f_term_1 = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Matrix_by_Scalar(material.particle_mechanics[particle]['C'],time_passed,true)
		#print(f_term_1,' time_passed * particle.C')
		f_term_2 =  get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Add_Matrix(material.particle_mechanics[particle]['I'],f_term_1)
		#print(f_term_2,' particle.I + f term 1')
		material.particle_mechanics[particle]['F'] = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Matrix(f_term_2,material.particle_mechanics[particle]['F'])
		### Updating Plasticity...
		material.particle_mechanics[particle]['F'] = get_tree().get_root().get_node("Test Area/Simulation/Constitutive Models").Update_Plasticity(particle,material.type_of_substance,material.yield_surface,material.particle_mechanics[particle]['F'])
	
	#return material
