extends Node


var barriers : Dictionary
var grid_nodes : Dictionary = {}
var gravity : Vector2 = Vector2(0.0,9.8)
# During ParticletoGrid...
var kernel_distance : Vector2
var flipped_kernel_distance : Vector2
var building_weights : Vector2
var weight_interpolation : float
var building_momentum
var forces : Vector2 = Vector2(0.0,0.0)
var force_term_1
# During GridUpdate-Collisions...
var area_multiplier : int
var distance_between : Vector2
var length_covered_area
var width_covered_area
var covered_area
var percentage_covered
var particles_set_to_merge
# During GridtoParticles...
var construct_B

func _on_Program_ready():
	### setup up the barrier... 
	barriers['window outline'] = get_tree().get_root().get_node("Test Area").simulation_outline
	# allow each wall of the window to be a Constitution...
	#barriers['window outline']['top']['']
	#barriers['window outline']['right']
	#barriers['window outline']['bottom']
	#barriers['window outline']['left']
	### how the simulation window handles the borders...
	# 'thru' the particles continue thru and keep going...
	# 'disappear' the particles will breach the window then disappear.
	# '2.0' the particles will bounce at max + 1(explode)
	# '1.0' the particles will bounce at max(perfectly elastic)
	# '0.5' the particles will bounce at partially(partially inelastic)
	# '0.0' the particles will sticky to the border(perfect inelastic)
	barriers['coefficient of reinstition'] = 'disappear'
	

#func _notification(event):
	##print(event)
	
	#if event == Node.NOTIFICATION_WM_CLOSE_REQUEST:
	#	#print("Exit")
	#if event == Node.NOTIFICATION_WM_SIZE_CHANGED:
	#	#print('Window size change')
	#	pass
	
func Substance_Alignment():
	pass

#func B_formation(component,weight,other_particle_grid_velocity):
	### B is the sum of particle.velocity * weight_interpolation.

#	for other_particle in component.within_range:
			
#		kernel_distance = component.get_position() - other_particle.surrounding_area.position.get_center()
#		flipped_kernel_distance = other_particle.surrounding_area.position.get_center() - component.get_position()
			
#		component.velocity += (weight * other_particle_grid_velocity)
	
	#return B_source.velocity

func Weight_Interpolation(kernel:Vector2,cell_size:float):
	#reset parameters
	building_weights = Vector2(0.0,0.0)
	### Weight Interpolation...
	
	if -((3.0 * cell_size) / 2.0) <= kernel.x and kernel.x <= - ( (1.0*cell_size) / 2.0 ) and -((3.0 * cell_size) / 2.0) <= kernel.y and kernel.y <= - ( (1.0*cell_size) / 2.0 ):
		building_weights.x =  ((1.0 * pow(kernel.x,2.0)) / ( 2.0 * pow(cell_size,2.0))) + ((3.0 * kernel.x) / (2.0 * cell_size)) + (9.0/8.0)
		#building_weights.x =  ((1.0 * pow(kernel.x,2.0)) / ( 2.0 * pow(cell_size,2.0))) + ((3.0 * kernel.x) / (2.0 * cell_size)) + (9.0/8.0)
		
		building_weights.y =  ((1.0 * pow(kernel.y,2.0)) / ( 2.0 * pow(cell_size,2.0))) + ((3.0* kernel.y) / (2.0 * cell_size)) + (9.0/8.0)
		weight_interpolation = building_weights.x * building_weights.y
	elif -((1.0 * cell_size) / 2.0) <= kernel.x and kernel.x <=  ( (1.0*cell_size) / 2.0 ) and -((1.0 * cell_size) / 2.0) <= kernel.y and kernel.y <= ( (1.0*cell_size) / 2.0 ):
		building_weights.x = -((1.0 * pow(kernel.x,2.0)) / (pow(cell_size,2.0))) + (3.0/4.0)
		building_weights.y = -((1.0 * pow(kernel.y,2.0)) / (pow(cell_size,2.0))) + (3.0/4.0)
		weight_interpolation = building_weights.x * building_weights.y
	elif ((1.0 * cell_size) / 2.0) <= kernel.x and kernel.x <=  ( (3.0*cell_size) / 2.0 ) and ((1 * cell_size) / 2) <= kernel.y and kernel.y <=  ( (3.0*cell_size) / 2.0 ):
		building_weights.x = ((1.0 * pow(kernel.x,2.0)) / ( 2.0 * pow(cell_size,2))) - ((3.0 * kernel.x) / (2.0 * cell_size)) + (9.0/8.0)
		building_weights.y = ((1.0 * pow(kernel.y,2.0)) / ( 2.0 * pow(cell_size,2))) - ((3.0 * kernel.y) / (2.0 * cell_size)) + (9.0/8.0)
		weight_interpolation = building_weights.x * building_weights.y
	else:
		building_weights = Vector2(0.0,0.0)
		weight_interpolation = building_weights.x * building_weights.y
	
	return weight_interpolation
	
func Simulate(time_passed,grid_domain,material,the_grid):
	### misc...
	
	#var shape_function 
	
	grid_domain = float(grid_domain)
	
	### reseting the grid data...
	for particle in the_grid.keys():
		the_grid[particle] = {'mass':1.0,'velocity':Vector2(0.0,0.0),'momentum':Vector2(0.0,0.0)}
		# the particle comes in contact with each other...
		
		# establish the particle grid domain...
		area_multiplier = 1
		particle.surrounding_area = Rect2i(Vector2(particle.position.x - ((grid_domain/2)*area_multiplier),particle.position.y - ((grid_domain/2)*area_multiplier)),Vector2(grid_domain*area_multiplier,grid_domain*area_multiplier))
		
	#print(' ')
	#print(' Particles to Grid')
	# Particles to Grid:
	for particle in material:
		
		for other_particle in particle.within_range:
			
			### Weight Interpolation...
			Weight_Interpolation(particle.relation_to_domain[other_particle],grid_domain)
			#print(weight_interpolation,' w_ip')
			the_grid[particle]['mass'] = the_grid[particle]['mass'] + (other_particle.mass * weight_interpolation)
			#print(the_grid[particle]['mass'],' the lump sum')
			#the_grid[particle]['momentum'] += weight_interpolation * other_particle.mass * (other_particle.C * flipped_particle.relation_to_domain[other_particle])
			building_momentum = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Matrix_by_Vector2_to_Vector2(other_particle.C,particle.domain_relation_to_particle[other_particle])
			#print(building_momentum,' building momentum.')
			the_grid[particle]['momentum'] = the_grid[particle]['momentum'] + (weight_interpolation * other_particle.mass * (other_particle.velocity + building_momentum))
			#the_grid[particle]['momentum'] = the_grid[particle]['momentum'] + (weight_interpolation * other_particle.mass * building_momentum)
			#print(the_grid[particle]['momentum'],' momentum. during')
			
			# Forces : f =  - sum of particle.volume * particle.stress * flipped_kernel
			force_term_1 = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Matrix_by_Vector2_to_Vector2(other_particle.stress,particle.domain_relation_to_particle[other_particle])
			#print(force_term_1,' build forces')
			forces = -(forces + (other_particle.volume * force_term_1))
		
		#print(the_grid[particle]['momentum'],' momentum. after')
	#Grid Update:
	#print(' ')
	#print('Grid Update')
	for node in the_grid.keys():
		if the_grid[node]['mass'] > 0.0:
			#the_grid[node]['momentum'] = the_grid[node]['momentum'].normalized() + (time_passed * ( the_grid[node]['mass'] * gravity ))
			#the_grid[node]['momentum'] = the_grid[node]['momentum'] + (time_passed * ( the_grid[node]['mass'] * gravity ))
			#the_grid[node]['velocity'] = the_grid[node]['momentum'] / the_grid[node]['mass']
			#print(forces+gravity,' forces check')
			the_grid[node]['velocity'] = the_grid[node]['momentum'] + ((time_passed * (forces+gravity)) / the_grid[node]['mass'])
			#print(the_grid[node]['velocity'],' velocity at the end')
		else:
			###...
			##print('the particle has no mass.')
			pass
	
	### Collision Detection...
	### how the particles interact with the window outline...
	for particle in material:
		
		for other_particle in material:
			
			### if the particles come into contact with each other...
			#var contact_with_x = abs((particle.position.x - other_particle.surrounding_area.get_center().x) / grid_domain)
			#var contact_with_y = abs((particle.position.y - other_particle.surrounding_area.get_center().y) / grid_domain)
			var contact_with = (particle.get_position() - other_particle.surrounding_area.get_center())  / grid_domain
			#print(' ')
			#print(particle,' distance from other particle ',other_particle,' is ',contact_with)
			#print(grid_domain/grid_domain,' and ',-(grid_domain/grid_domain),' are the particle barrier')
			#print(particle,' relation to ',other_particle,' is ',particle.relation_to_domain[other_particle],' check before update')
			
			kernel_distance = (particle.get_position() - other_particle.surrounding_area.position) / grid_domain
			flipped_kernel_distance = (other_particle.surrounding_area.position - particle.get_position()) / grid_domain
			
			#print(kernel_distance,' is the kernel between the particles',particle,' ',other_particle)
			particle.relation_to_domain[other_particle] = kernel_distance
			particle.domain_relation_to_particle[other_particle] = flipped_kernel_distance
			#particle.relation_to_domain[other_particle] = flipped_kernel_distance
			#particle.domain_relation_to_particle[other_particle] = kernel_distance
			#print(particle,' relation to ',other_particle,' is ',particle.relation_to_domain[other_particle],' check after update')
			
			if  contact_with.x == 0 and contact_with.y == 0:
				### always in contact with itself...
				#print(particle,' it recognize itself')
				pass
			elif contact_with.x > -(grid_domain/grid_domain) and contact_with.x < grid_domain/grid_domain and contact_with.y > -(grid_domain/grid_domain) and contact_with.y < grid_domain/grid_domain:
				### come into some contact...
				#print(particle,' at position ',particle.position,' in contact with ',other_particle,' at position ',other_particle.surrounding_area,' quadrant')
				if particle.within_range.has(other_particle):
					### the particle is already in contact with the other...
					pass
				else:
					### recognized to be within the domain range...
					particle.within_range.append(other_particle)
			else:
				###
				#print(particle,' not in contact with this particle ',other_particle)
				particle.within_range.erase(other_particle)
				other_particle.within_range.erase(particle)


	# Grid to Particle...
	#print(' ')
	#print('Grid to Particle')
	# resetting of the particles...
	for particle in material:
		particle.mass = get_tree().get_root().get_node("Test Area/Simulation/Substance").default_mass_of_particle
		#particle.velocity = Vector2(0.0,0.0)
		#particle.B = [[0.0,0.0],[0.0,0.0]]
	
	for particle in material:
		
		for other_particle in particle.within_range:
			
			### Weight Interpolation...
			Weight_Interpolation(particle.relation_to_domain[other_particle],grid_domain)
			#print(weight_interpolation,' during grid to particle.')
			particle.velocity = particle.velocity + (weight_interpolation * the_grid[other_particle]['velocity'])
			#await B_formation(particle,weight_interpolation,the_grid[other_particle]['velocity'])
			#print(particle.velocity,' the particle velocity')
			#particle.B = particle.velocity * flipped_kernerl_distance^Transposed
			construct_B = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Vector2_by_Vector2_to_Matrix(particle.velocity,false,particle.domain_relation_to_particle[other_particle],true)
			#print(construct_B,' build B')
			particle.B = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Add_Matrix(particle.B,construct_B)
			#print(particle.B,' B')
		#C Updated..
		#particle.C = particle.B * D^-1  :: D^-1 = 4/grid_domain^2 * particle.I^-1
		var c_term = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Matrix_by_Scalar(particle.B,4,true)
		var c_term_2 = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Divide_Matrix_by_Scalar(c_term,pow(grid_domain,2),true)
		var c_term_3 = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Inverse_Matrix(particle.I)
		particle.C = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Matrix(c_term_3,c_term_2)
		#print(particle.C,' C')
		#particle.C = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Matrix_by_Vector2_to_Matrix(c_term_4,flipped_kernel_distance,)
		
		
		###particle position update...
		#print(particle.position,' position before')
		#print(time_passed,' time passing')
		#print(particle.velocity,' velocities')
		particle.set_position(particle.get_position() + (time_passed * particle.velocity))
		#print(particle.position,' position after')
		###deformation update...
		#particle.F 
