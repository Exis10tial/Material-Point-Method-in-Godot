extends Node


var particle_interaction
var matrix_math
var constitutive_models

var grid_nodes : Dictionary = {}
#var gravity : Vector2
var gravity : Vector2 = Vector2(0.00,9.80) * 100
#var gravity : Vector2 = Vector2(0.00,-9.80)
#var gravity : Vector2 = Vector2(9.80,0.0)
#var gravity : Vector2 = Vector2(-9.80,0.0)
#var gravity : Vector2 = Vector2(randf_range(-900.80,900.80),randf_range(-900.80,900.80))
#var apply_outside_forces
var identify_number : int
var particle : String
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
var collection_of_velocities : Array
var relation_between_particles
# handling particle merging...
var percentage_covered : float
var cycle_of_mergin_particles : Array = []
var merge_particles : bool
var list_of_particles_merged_particles : Array = []
# wall mechanics
var barriers : Dictionary = {'window outline': {}}
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
	barriers['window outline']['top'] = {'coefficient of restitution' : 0.0,'coefficient of static friction': 0.0,'coefficient of kinetic friction': 0.0,'outline':0,'mass':0.0,'velocity':Vector2(0.0,0.0),'momentum':Vector2(0.0,0.0)}
	barriers['window outline']['right'] = {'coefficient of restitution' : 0.0,'coefficient of static friction': 0.0,'coefficient of kinetic friction': 0.0,'outline':ProjectSettings.get_setting('display/window/size/width'),'mass':0.0,'velocity':Vector2(0.0,0.0),'momentum':Vector2(0.0,0.0)}
	barriers['window outline']['bottom'] = {'coefficient of restitution' : 0.0,'coefficient of static friction': 0.0,'coefficient of kinetic friction': 0.0,'outline':ProjectSettings.get_setting('display/window/size/height'),'mass':0.0,'velocity':Vector2(0.0,0.0),'momentum':Vector2(0.0,0.0)}
	barriers['window outline']['left'] ={'coefficient of restitution' : 0.0,'coefficient of static friction': 0.0,'coefficient of kinetic friction': 0.0,'outline':0,'mass':0.0,'velocity':Vector2(0.0,0.0),'momentum':Vector2(0.0,0.0)}
		
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
	window_center = Vector2(ProjectSettings.get_setting('display/window/size/width')/2.0,ProjectSettings.get_setting('display/window/size/height')/2.0)
	for wall in barriers['window outline'].keys():
		barriers['window outline'][wall]['coefficient of restitution'] = 1.00
		#barriers['window outline'][wall]['coefficient of restitution'] = 0.50
		#barriers['window outline'][wall]['coefficient of restitution'] = 0.0
		barriers['window outline'][wall]['coefficient of static friction'] = 1.0
		barriers['window outline'][wall]['coefficient of kinetic friction'] = 1.0
		barriers['window outline'][wall]['mass'] = 10000.00
		#barriers['window outline'][wall]['mass'] = 1.00
		barriers['window outline']['top']['outline'] = 0.0
		barriers['window outline']['right']['outline'] = ProjectSettings.get_setting('display/window/size/width')
		barriers['window outline']['bottom']['outline'] = ProjectSettings.get_setting('display/window/size/height')
		barriers['window outline']['left']['outline'] = 0.0
		#barriers['window outline'][wall]['velocity'] = Vector2(1.0,1.0)
		
	


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
			
				building_weights.x =  snapped( ((1.0 * pow(kernel.x,2.0)) /  (2.0 * pow(cell_size,2.0))),.01 ) + snapped(((3.0 * kernel.x) / (2.0 * cell_size)),.01) + snapped((9.0/8.0),.001)
				building_weights.y =  snapped( ((1.0 * pow(kernel.y,2.0)) /  (2.0 * pow(cell_size,2.0))),.01 ) + snapped(((3.0 * kernel.y) / (2.0 * cell_size)),.01) + snapped((9.0/8.0),.001)

				weight_interpolation = snapped(building_weights.x * building_weights.y,.01)
				
			elif kernel.x > -((1.0 * cell_size) / 2.0) and kernel.x <  ( (1.0*cell_size) / 2.0 ) and kernel.y > -((1.0 * cell_size) / 2.0) and kernel.y < ( (1.0*cell_size) / 2.0 ):
			
				building_weights.x = snapped(((-1.0 * pow(kernel.x,2.0)) / (pow(cell_size,2.0))),.01) + snapped((3.0/4.0),.01)
				building_weights.y = snapped(((-1.0 * pow(kernel.y,2.0)) / (pow(cell_size,2.0))),.01) + snapped((3.0/4.0),.01)

				weight_interpolation = snapped(building_weights.x * building_weights.y,.01)
				
			elif kernel.x > ((1.0 * cell_size) / 2.0) and kernel.x <  ( (3.0*cell_size) / 2.0 ) and  kernel.y > ((1 * cell_size) / 2) and kernel.y <  ( (3.0*cell_size) / 2.0 ):

				building_weights.x =  snapped( ((1.0 * pow(kernel.x,2.0)) /  (2.0 * pow(cell_size,2.0))),.01 ) - snapped(((3.0 * kernel.x) / (2.0 * cell_size)),.01) + snapped((9.0/8.0),.001)
				building_weights.y =  snapped( ((1.0 * pow(kernel.y,2.0)) /  (2.0 * pow(cell_size,2.0))),.01 ) - snapped(((3.0 * kernel.y) / (2.0 * cell_size)),.01) + snapped((9.0/8.0),.001)

				weight_interpolation = snapped(building_weights.x * building_weights.y,.01)
			else:
				building_weights = Vector2(0.0,0.0)
				weight_interpolation = snapped(building_weights.x * building_weights.y,.01)
				
	elif splines == "cubic":
		pass
	
	return weight_interpolation
	
func Simulate(time_passed:float,material:Object,the_grid:Dictionary):
	#particle_interaction = get_tree().get_root().get_node("Simulation/Particle Interaction")
	#matrix_math = get_tree().get_root().get_node("Simulation/Matrix Math")
	#constitutive_models = get_tree().get_root().get_node("Simulation/Constitutive Models")
	pass
	#"""

func Grid_Reset(material:Object,the_grid:Dictionary):
	#particle_interaction = get_tree().get_root().get_node("Simulation/Particle Interaction")
	matrix_math = get_tree().get_root().get_node("Simulation/Matrix Math")
	constitutive_models = get_tree().get_root().get_node("Simulation/Constitutive Models")
	### particle simulation...
	### reseting the grid data...
	#for particle in the_grid.keys():
	identify_number = 0
	particle = 'null'
	#print('Grid_Reset check')
	while true:
		if identify_number >= len(material.particle_mechanics.keys()):
			break
		
		particle = material.particle_mechanics.keys()[identify_number]
		
		#the_grid[particle] = {'mass': material.mass_in_pieces,'velocity':Vector2(0.0,0.0),'momentum':Vector2(0.0,0.0)}
		the_grid[particle] = {'mass': material.mass_in_pieces,'velocity':Vector2(0.0,0.0),'momentum':Vector2(0.0,0.0)}
		
		#material.particle_mechanics[particle]['F'] = material.particle_mechanics[particle]['I'].duplicate(true)
		material.particle_mechanics[particle]['J'] = matrix_math.Find_Determinant(material.particle_mechanics[particle]['F'])
		material.particle_mechanics[particle]['volume'] = material.volume_in_pieces * material.particle_mechanics[particle]['J']
		
		
		identify_number = wrapi(identify_number+1,0,len(material.particle_mechanics.keys())+1)

func Particles_to_Grid(time_passed:float,material:Object,the_grid:Dictionary):
	#particle_interaction = get_tree().get_root().get_node("Simulation/Particle Interaction")
	matrix_math = get_tree().get_root().get_node("Simulation/Matrix Math")
	constitutive_models = get_tree().get_root().get_node("Simulation/Constitutive Models")
	# Particles to Grid:
	#for particle in material.particle_mechanics.keys():
	identify_number = 0
	particle = 'null'
	#print('Particles_to_Grid check')
	while true:
		if identify_number >= len(material.particle_mechanics.keys()):
			break
		
		particle = material.particle_mechanics.keys()[identify_number]
		#### the sum of aspects that is resetted...
		affine_momentum_fuse = Vector2(0.0,0.0)
		var sum_of_f_forces = Vector2(0.0,0.0)
		var transfer_momentum = Vector2(0.0,0.0)
		var weighted_mass = 0.0
		#print(' ')
		#for other_particle in material.particle_mechanics.keys():
		for other_particle in material.particle_mechanics[particle]['within_range']:
			
			### Weight Interpolation...
				
				#print(other_particle, ' other part')
			relation_between_particles = Vector2(snapped(material.particle_lineation[particle].get_center().x - material.particle_lineation[other_particle].get_center().x,.01),snapped(material.particle_lineation[particle].get_center().y - material.particle_lineation[other_particle].get_center().y,.01))
				
				#Weight_Interpolation(basis,basis_function_version,material.particle_mechanics[particle]['relation_to_domain'][other_particle],material.cell_size)
			Weight_Interpolation(basis,basis_function_version,relation_between_particles,material.cell_size)
				#"""
			if material.physical_state == 'solid':
				if material.constitutive_model == 'hyperelastic':
					material.particle_mechanics[other_particle].stress = constitutive_models.Neo_Hookean(particle,material)
				if material.constitutive_model == 'fixed_corated':
					material.particle_mechanics[other_particle].stress = constitutive_models.Fixed_Corotated(particle,material)
				if material.constitutive_model == 'drucker_prager_elasticity':
					material.particle_mechanics[other_particle].stress = constitutive_models.Drucker_Prager_Elasticity(particle,material)
			elif material.physical_state == 'liquid':
				if material.constitutive_model == 'water':
					material.particle_mechanics[other_particle].stress = constitutive_models.Model_of_Water(particle,material)
					#print(other_particle.stress,' checking stress')
			elif material.physical_state == 'gas':
				pass
			
				#
			### gradient-weight : weight_interpolation * D ^-1 * flipped_kernel
			#D **-1 = (4/pow(cell_size,2)) * other_particle.I^-1
			var q_weighted_time_coefficient = snapped(time_passed * weight_interpolation,.01)
			var q_timed_weighted_intial_volume_coefficient = snapped(q_weighted_time_coefficient * material.volume,.01)
			gradient_weight_term_1 = snapped((snapped((q_timed_weighted_intial_volume_coefficient * basis_coefficient),.01) / pow(material.cell_size,2.0)),.01)
			var q_timed_weighted_volume_d = matrix_math.Multiply_Matrix_by_Scalar(material.particle_mechanics[other_particle]['I'],gradient_weight_term_1,true)
			var particle_volume_q_d = matrix_math.Multiply_Matrix_by_Scalar(q_timed_weighted_volume_d,material.particle_mechanics[particle]['volume'],true)
			var q_timed_volume_weighted_d = matrix_math.Multiply_Matrix_by_Vector2_to_Vector2(particle_volume_q_d,relation_between_particles)
	
			var volumed_stress = matrix_math.Multiply_Matrix_by_Scalar(material.particle_mechanics[other_particle]['stress'],material.particle_mechanics[particle]['volume'],true)
			var q_stressed_timed_volume_weighted_d = matrix_math.Multiply_Matrix_by_Vector2_to_Vector2(volumed_stress,q_timed_volume_weighted_d)
	
			mass_c_coefficient =  matrix_math.Multiply_Matrix_by_Scalar(material.particle_mechanics[other_particle]['C'],material.particle_mechanics[other_particle]['mass'],true)
			var related_mass_c_coeffiecient = matrix_math.Multiply_Matrix_by_Vector2_to_Vector2(mass_c_coefficient,relation_between_particles)
			
			affine_momentum_fuse.x = snapped((affine_momentum_fuse.x + (snapped((-1*snapped((q_stressed_timed_volume_weighted_d.x + related_mass_c_coeffiecient.x),.01)),.01))),.01)
			affine_momentum_fuse.y = snapped((affine_momentum_fuse.y + (snapped((-1*snapped((q_stressed_timed_volume_weighted_d.y + related_mass_c_coeffiecient.y),.01)),.01))),.01)
			
			transfer_momentum.x = snapped((transfer_momentum.x + (snapped((material.particle_mechanics[other_particle]['mass'] * material.particle_mechanics[other_particle]['velocity'].x),.01))),.01)
			transfer_momentum.y = snapped((transfer_momentum.y + (snapped((material.particle_mechanics[other_particle]['mass'] * material.particle_mechanics[other_particle]['velocity'].y),.01))),.01)
			
			weighted_mass = snapped((weighted_mass + (snapped((material.particle_mechanics[particle]['mass'] * weight_interpolation),.01))),.01)
		
		#print(transfer_momentum, ' transfer_momentum check')
		#print(time_passed, ' time_passed check')
		#print(weighted_mass, ' weighted_mass check')
		#print(gravity, ' gravity check')
		#print(affine_momentum_fuse, ' affine_momentum_fuse check')
		
		the_grid[particle]['momentum'] = transfer_momentum +( time_passed * (weighted_mass*gravity + affine_momentum_fuse ))
		#print(the_grid[particle]['momentum'], 'the_grid[particle][momentum] check')
		identify_number = wrapi(identify_number+1,0,len(material.particle_mechanics.keys())+1)

func Grid_Update(material:Object,the_grid:Dictionary):#,outside_forces:Vector2):
	#print('Grid_Update check')
	#Grid Update:
	for particle in the_grid.keys():
		#print(the_grid[particle]['mass'],' the_grid[particle][mass] check')
		if the_grid[particle]['mass'] > 0.0:
			
			the_grid[particle]['velocity'].x = snapped(the_grid[particle]['momentum'].x / the_grid[particle]['mass'],.01)
			the_grid[particle]['velocity'].y = snapped(the_grid[particle]['momentum'].y / the_grid[particle]['mass'],.01)
			
			#print(the_grid[particle]['velocity'],' after update')
		else:
			###...
			##print('the particle has no mass.')
			pass
		

#func Collision_Detection(material:Object,the_grid:Dictionary):
func Collision_with_Wall(material:Object,the_grid:Dictionary):
	#print('Collision_with_Wall check')
	particle_interaction = get_tree().get_root().get_node("Simulation/Particle Interaction")
	#matrix_math = get_tree().get_root().get_node("Simulation/Matrix Math")
	#constitutive_models = get_tree().get_root().get_node("Simulation/Constitutive Models")
	#"""
	### Collision Detection...
	### how the particles interacts....
	collection_of_velocities = []
	for particle in material.particle_lineation.keys():
		### if the particle has contacted the window outline...
		#print(the_grid[particle]['velocity'],' incoming')
		if material.particle_lineation[particle].position.y < barriers['window outline']['top']['outline']:
			### the particle left the window out line at the sides...
			
			the_grid[particle]['velocity'] = particle_interaction.Collision_with_Walls('top',material,particle,material.particle_lineation[particle],barriers,the_grid,material.cell_size)
			#collection_of_velocities.append(particle_interaction.Collision_with_Walls('top',material,particle,material.particle_lineation[particle],barriers,the_grid,material.cell_size))
		if material.particle_lineation[particle].end.x >= barriers['window outline']['right']['outline']:
			### the particle left the window out line at the sides...
			
			the_grid[particle]['velocity'] = particle_interaction.Collision_with_Walls('right',material,particle,material.particle_lineation[particle],barriers,the_grid,material.cell_size)
			#collection_of_velocities.append(particle_interaction.Collision_with_Walls('right',material,particle,material.particle_lineation[particle],barriers,the_grid,material.cell_size))
		if material.particle_lineation[particle].end.y >= barriers['window outline']['bottom']['outline']:
			### the particle left the window out line at the sides...

			the_grid[particle]['velocity'] = particle_interaction.Collision_with_Walls('bottom',material,particle,material.particle_lineation[particle],barriers,the_grid,material.cell_size)
			#collection_of_velocities.append(particle_interaction.Collision_with_Walls('bottom',material,particle,material.particle_lineation[particle],barriers,the_grid,material.cell_size))
		if material.particle_lineation[particle].position.x < barriers['window outline']['left']['outline']:
			### the particle left the window out line at the sides...
			
			the_grid[particle]['velocity'] =  particle_interaction.Collision_with_Walls('left',material,particle,material.particle_lineation[particle],barriers,the_grid,material.cell_size)
			#collection_of_velocities.append(particle_interaction.Collision_with_Walls('left',material,particle,material.particle_lineation[particle],barriers,the_grid,material.cell_size))
		else:
			### the particle is within the window outline...
			pass
		#for velocity in collection_of_velocities:
			#print(particle,' in coming velocity ',velocity)
			#the_grid[particle]['velocity'] = the_grid[particle]['velocity'] + velocity
	
func Collision_with_Other_Particles(material:Object,the_grid:Dictionary):
	### Collision between Particles...
	identify_number = 0
	particle = 'null'
	#print('Collision_with_Other_Particles check')
	while true:
		if identify_number >= len(material.particle_mechanics.keys()):
			break
	
		particle = material.particle_mechanics.keys()[identify_number]
	
		#collection_of_velocities = []
		for other_particle in material.particle_mechanics[particle]['within_range']:
			if particle == other_particle:
				#### it is itself..
				pass
			else:
				#the_grid[particle]['velocity'] = the_grid[particle]['velocity'] + particle_interaction.Collision_between_Other_Particles(particle,material,material.particle_lineation[particle],other_particle,material,material.particle_lineation[other_particle],the_grid)
				the_grid[particle]['velocity'] = particle_interaction.Collision_between_Other_Particles(particle,material,material.particle_lineation[particle],other_particle,material,material.particle_lineation[other_particle],the_grid)
				#collection_of_velocities.append(particle_interaction.Collision_between_Other_Particles(particle,material,material.particle_lineation[particle],other_particle,material,material.particle_lineation[other_particle],the_grid))
		
		### loop counter...
		identify_number = wrapi(identify_number+1,0,len(material.particle_mechanics.keys())+1)
	

func Particle_Reset(material:Object):
	### particle is reset...
	identify_number = 0
	particle = 'null'
	#print('Particle_Reset check')
	while true:
		if identify_number >= len(material.particle_mechanics.keys()):
			break
		
		particle = material.particle_mechanics.keys()[identify_number]
		
		#material.particle_mechanics[particle]['mass'] = material.mass_in_pieces
		#material.particle_mechanics[particle]['velocity'] = material.maintain_velocity
		material.particle_mechanics[particle]['mass'] = material.mass_in_pieces
		material.particle_mechanics[particle]['velocity'] = Vector2(0.0,0.0)
		### loop counter...
		identify_number = wrapi(identify_number+1,0,len(material.particle_mechanics.keys())+1)


func Grid_to_Particle(time_passed:float,material:Object,the_grid:Dictionary):
	#particle_interaction = get_tree().get_root().get_node("Simulation/Particle Interaction")
	matrix_math = get_tree().get_root().get_node("Simulation/Matrix Math")
	#constitutive_models = get_tree().get_root().get_node("Simulation/Constitutive Models")
	#print('Grid_to_Particle check')
	identify_number = 0
	particle = 'null'
	while true:
		if identify_number >= len(material.particle_mechanics.keys()):
			break
		
		particle = material.particle_mechanics.keys()[identify_number]
	
	#for particle in material.particle_mechanics.keys():
		#### sum of aspects reset...
		
		sum_of_c_weighted_velocity_matrix = [0,0,0,0]
		#for other_particle in material.particle_mechanics.keys():
		var setup_grid_velocity : Vector2 = Vector2(0,0)
		for other_particle in material.particle_mechanics[particle]['within_range']:
			
			### Weight Interpolation...
			relation_between_particles = Vector2(snapped(material.particle_lineation[particle].get_center().x - material.particle_lineation[other_particle].get_center().x,.01),snapped(material.particle_lineation[particle].get_center().y - material.particle_lineation[other_particle].get_center().y,.01))
				
			#Weight_Interpolation(basis,basis_function_version,material.particle_mechanics[particle]['relation_to_domain'][other_particle],material.cell_size)
			Weight_Interpolation(basis,basis_function_version,relation_between_particles,material.cell_size)
				
			#print(weight_interpolation,' wip')
			#MLS MPM 
			#particle.velocity = sum of (the_grid[velocity] * weight_interpolation)
			#particle.C = D^-1 * sum of the_grid( weight_interpolation * ['velocity'] * flipped_kernel_distance^T )
			# D^-1 = (4/pow(cell_size,2)) * particle.I^-1
				
			material.particle_mechanics[particle]['velocity'].x = snapped(material.particle_mechanics[particle]['velocity'].x + (snapped((the_grid[other_particle]['velocity'].x * weight_interpolation),.01)),.01)
			material.particle_mechanics[particle]['velocity'].y = snapped(material.particle_mechanics[particle]['velocity'].y + (snapped((the_grid[other_particle]['velocity'].y * weight_interpolation),.01)),.01)
				
			
			setup_grid_velocity.x = snapped((the_grid[other_particle]['velocity'].x * weight_interpolation),.01)
			setup_grid_velocity.y = snapped((the_grid[other_particle]['velocity'].y * weight_interpolation),.01)
					
			c_flipped_velocity_coefficient = matrix_math.Multiply_Vector2_by_Vector2_to_Matrix(setup_grid_velocity,false,relation_between_particles,true)
			#c_weighted_velocity_matrix = matrix_math.Multiply_Matrix_by_Scalar(c_flipped_velocity_coefficient,weight_interpolation,true)
			sum_of_c_weighted_velocity_matrix = matrix_math.Add_Matrix(sum_of_c_weighted_velocity_matrix,c_flipped_velocity_coefficient)
					
		### Compute C...
		c_cell_coefficient =  snapped((basis_coefficient / pow(material.cell_size,2.0)),.01)
		c_inverse_I_coefficient = matrix_math.Inverse_Matrix(material.particle_mechanics[particle]['I'])
		c_cell_inverse_coefficient = matrix_math.Multiply_Matrix_by_Scalar(c_inverse_I_coefficient,c_cell_coefficient,true)
				
		material.particle_mechanics[particle]['C'] = matrix_math.Multiply_Matrix(c_cell_inverse_coefficient,sum_of_c_weighted_velocity_matrix)
		
		###particle position update...
		#print(material.particle_mechanics[particle]['velocity'],' velocity check')
		material.particle_lineation[particle].position.x = snapped((material.particle_lineation[particle].position.x + snapped((time_passed * material.particle_mechanics[particle]['velocity'].x),.01)  ),.01)
		material.particle_lineation[particle].position.y = snapped((material.particle_lineation[particle].position.y + snapped((time_passed * material.particle_mechanics[particle]['velocity'].y),.01)  ),.01)
		#print(material.particle_lineation[particle].position,' material.particle_lineation[particle].position')
		###deformation update...
		#particle.F = (particle.I + time_passed * particle.C ) * particle.F
		
		f_term_1 = matrix_math.Multiply_Matrix_by_Scalar(material.particle_mechanics[particle]['C'],time_passed,true)
		#print(f_term_1,' time_passed * particle.C')
		f_term_2 =  matrix_math.Add_Matrix(material.particle_mechanics[particle]['I'],f_term_1)
		#print(f_term_2,' particle.I + f term 1')
		material.particle_mechanics[particle]['F'] = matrix_math.Multiply_Matrix(f_term_2,material.particle_mechanics[particle]['F'])
		### Updating Plasticity...
		material.particle_mechanics[particle]['F'] = constitutive_models.Update_Plasticity(material.type_of_substance,material.yield_surface,material.particle_mechanics[particle]['F'])
		
		identify_number = wrapi(identify_number+1,0,len(material.particle_mechanics.keys())+1)
	#"""
