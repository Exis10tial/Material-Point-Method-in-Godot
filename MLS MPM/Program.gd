extends Node


var particle_interaction
var matrix_math
var constitutive_models


### Eulerian/Cartesian Grid...
#var Grid_Nodes : Dictionary = {}
var grid_cells : int
var grid_spacing : int
var grid_nodes : Dictionary = {}
var station : Vector2
#var identify_grid_reach : Array = []
#var gravity : Vector2
var gravity : Vector2 = Vector2(0.00,9.80) * 100
#var gravity : Vector2 = Vector2(0.00,-9.80)* 100
#var gravity : Vector2 = Vector2(9.80,0.0) * 100
#var gravity : Vector2 = Vector2(-9.80,0.0) * 100
#var apply_outside_forces
var identify_number : int
var particle : String
var node = null
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
var updating_plasticity : bool = false
var relation_between_particle_and_grid_node : Vector2
var relation_between_grid_node_and_particle : Vector2
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
var weight_interpolation : Dictionary = {'weight':0.0,'gradient':0.0}
var forces = [0,0,0,0]
#var weight_interpolation : float
#var gradient_weight_interpolation : float
var building_momentum
#var forces : Vector2 = Vector2(0.0,0.0)
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
	
	##basis = "cubic"
	##basis_coefficient = 3.0
	
	basis = "quadratic"
	basis_coefficient = 4.0
	basis_function_version = 1
	
	#grid_cells = ProjectSettings.get_setting('display/window/size/viewport_width') * ProjectSettings.get_setting('display/window/size/viewport_height')
	#grid_spacing = 1
	
	grid_cells = 16*9
	grid_spacing = 120
	
	
	### setup up the barrier... 
	barriers['window outline']['top'] = {'coefficient of restitution' : 0.0,'coefficient of static friction': 0.0,'coefficient of kinetic friction': 0.0,'outline':0,'mass':0.0,'velocity':Vector2(0.0,0.0),'momentum':Vector2(0.0,0.0)}
	barriers['window outline']['right'] = {'coefficient of restitution' : 0.0,'coefficient of static friction': 0.0,'coefficient of kinetic friction': 0.0,'outline':ProjectSettings.get_setting('display/window/size/viewport_width'),'mass':0.0,'velocity':Vector2(0.0,0.0),'momentum':Vector2(0.0,0.0)}
	barriers['window outline']['bottom'] = {'coefficient of restitution' : 0.0,'coefficient of static friction': 0.0,'coefficient of kinetic friction': 0.0,'outline':ProjectSettings.get_setting('display/window/size/viewport_height'),'mass':0.0,'velocity':Vector2(0.0,0.0),'momentum':Vector2(0.0,0.0)}
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
	window_center = Vector2(ProjectSettings.get_setting('display/window/size/viewport_width')/2.0,ProjectSettings.get_setting('display/window/size/viewport_height')/2.0)
	for wall in barriers['window outline'].keys():
		barriers['window outline'][wall]['coefficient of restitution'] = 1.00
		#barriers['window outline'][wall]['coefficient of restitution'] = 0.50
		#barriers['window outline'][wall]['coefficient of restitution'] = 0.0
		#barriers['window outline'][wall]['coefficient of static friction'] = 1.0
		barriers['window outline'][wall]['coefficient of static friction'] = 0.80
		#barriers['window outline'][wall]['coefficient of static friction'] = 0.0
		#barriers['window outline'][wall]['coefficient of kinetic friction'] = 1.0
		barriers['window outline'][wall]['coefficient of kinetic friction'] = .650
		#barriers['window outline'][wall]['coefficient of kinetic friction'] = 0.0
		barriers['window outline'][wall]['mass'] = 1.00
		#barriers['window outline'][wall]['mass'] = 1.00
		barriers['window outline']['top']['outline'] = 0.0
		barriers['window outline']['right']['outline'] = ProjectSettings.get_setting('display/window/size/viewport_width')
		barriers['window outline']['bottom']['outline'] = ProjectSettings.get_setting('display/window/size/viewport_height')
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


func Weight_Interpolation(splines:String,version:int,updated_relation:Vector2):#,cell_size:float):
	#reset parameters
	building_weights = Vector2(0.0,0.0)
	### Weight Interpolation...
	
	if splines == "quadratic":
		#print('weight interpolation')
		if version == 1:
			### The Material Point Method for Simulating Continuum Materials pg 33 equation 123
			
			if 1.5 <= abs(updated_relation.x):
				### particles is not in range of the particle to influence.
				#print('3')
				building_weights.x =  0.0
			
			elif abs(updated_relation.x) <= .50 and abs(updated_relation.x) < 1.5:
				#print('2')
				building_weights.x = snappedf((.50 *  pow((1.50 - abs(updated_relation.x)),2.0)),.01)
				
			#if int(relative.x) >= (0*cell_size) and int(relative.x) <= (1*cell_size):
			#if abs(relative.x) == (0*cell_size) and abs(relative.x) <= ((cell_size/4.0)):
			elif 0 <= abs(updated_relation.x) and abs(updated_relation.x) < .50:
				#print('1')
				building_weights.x = snappedf((.75 - snappedf(pow(abs(updated_relation.x),2.0),.01)),.01)
				#print(building_weights,' a building_weights')
			#elif int(relative.x) >= (-.50)*(cell_size) and int(relative.x) <= (.50)*(cell_size):
			#elif abs(relative.x) > ((cell_size/4.0)) and abs(relative.x) <= ((cell_size/2.0)):
			
			#elif int(relative.x) >= (.50)*(cell_size) and int(relative.x) <= (1.50)*(cell_size):
			#elif abs(relative.x) > ((cell_size/2.0)):# and int(relative.x) <= (cell_size/2.0):
			
				
			if 1.5 <= abs(updated_relation.y):
				#print('8')
				building_weights.y = 0.0
			
			elif abs(updated_relation.y) <= .50 and abs(updated_relation.y) < 1.5:
				#print('7')
				building_weights.y = snappedf((.50 * pow((1.50 - abs(updated_relation.y) ),2.0)),.01)
				
			#if int(relative.y) >= (-1.50)*(cell_size) and int(relative.y) <= (-.50)*(cell_size):
			#if int(relative.y) >= (0*cell_size) and int(relative.x) <= (1*cell_size):
			#if abs(relative.y) == (0*cell_size) and abs(relative.y) <= ((cell_size/4.0)):
			elif 0 <= abs(updated_relation.y) and abs(updated_relation.y) < .50:
				#print('6')
				building_weights.y = snappedf((.75 - snappedf(pow(abs(updated_relation.y),2.0),.01)),.01)
				#print(building_weights,' a building_weights')
			#elif int(relative.y) >= (-.50)*(cell_size) and int(relative.y) <= (.50)*(cell_size):
			#elif int(relative.y) >= (1*cell_size) and int(relative.x) <= (2*cell_size):
			#elif abs(relative.y) > ((cell_size/4.0)) and abs(relative.y) <= ((cell_size/2.0)):
			
			#elif int(relative.y) >= (.50)*(cell_size) and int(relative.y) <= (1.50)*(cell_size):
			#elif int(relative.y) >= (2*cell_size):#-(cell_size/4.0) and int(relative.y) <= (cell_size/2.0):
			#elif abs(relative.y) > ((cell_size/2.0)):
			
		elif version == 2:
			### thru Analysis and Reduction of Quadrature Errors in the Material Point Method. steffan .W.pdf-pg(8.equation 16)
			pass
			
		###...
		#print(building_weights.x,' building weights x')
		#print(building_weights.y,' building weights y')
		#weight_interpolation['weight'] = ((1.0/cell_size)*building_weights.x ) * ((1.0/cell_size)*building_weights.y)
		#weight_interpolation['weight'] =(building_weights.x/cell_size) * (building_weights.y/cell_size)
		weight_interpolation['weight'] = snappedf((building_weights.x * building_weights.y),.01)
		
		### unneeded in mls mpm cause of the affine and force combination...
		#weight_interpolation['gradient'] =((1.0/cell_size) * ((1.0/cell_size)*building_weights.y)) * (((1.0/cell_size)*building_weights.x ) * (1.0/cell_size))
		#weight_interpolation['gradient'] = Vector2((( (1.0/cell_size) * ((1.0/cell_size)*building_weights.x ) ) * ((1.0/cell_size)*building_weights.y)),   ((1.0/cell_size)*building_weights.x ) *  ((1.0/cell_size) * ((1.0/cell_size)*building_weights.y)))
		
	elif splines == "cubic":
		pass
	
	return weight_interpolation
	
func Simulate(time_passed:float,material:Object,the_grid:Dictionary):
	#particle_interaction = get_tree().get_root().get_node("Simulation/Particle Interaction")
	#matrix_math = get_tree().get_root().get_node("Simulation/Matrix Math")-2499*
	#constitutive_models = get_tree().get_root().get_node("Simulation/Constitutive Models")
	pass
	#"""

func Grid_Reset(material:Object):#,matrix_math:Object):#,the_grid:Dictionary):
	#particle_interaction = get_tree().get_root().get_node("Simulation/Particle Interaction")
	matrix_math = get_tree().get_root().get_node("Simulation/Matrix Math")
	
	
	#constitutive_models = get_tree().get_root().get_node("Simulation/Constitutive Models")
	### particle simulation...
	### reseting the grid data...
	#for particle in the_grid.keys():
	identify_number = 0
	particle = 'null'
	grid_nodes = {}
	#print('Grid_Reset check')

	while true:
		if identify_number >= len(material.particle_mechanics.keys()):
			break
		
		particle = material.particle_mechanics.keys()[identify_number]
		
		### relation between the particle wth the grid...
		
		# find the nearest node (non-float number)
		var euler_base = Vector2(round(material.particle_lineation[particle].origin.x),round(material.particle_lineation[particle].origin.y))#+Vector2(0.5,0.5)
		
		# from the base node the 9 surrounding nodes...
		#material.particle_mechanics[particle]['eulerian'] = [Vector2(round(material.particle_lineation[particle].origin.x),round(material.particle_lineation[particle].origin.y))]
		#"""
		material.particle_mechanics[particle]['eulerian'] = [Vector2(euler_base.x-1,euler_base.y+1),
		Vector2(euler_base.x,euler_base.y+1),
		Vector2(euler_base.x+1,euler_base.y+1),
		Vector2(euler_base.x-1,euler_base.y),
		Vector2(euler_base.x,euler_base.y),
		Vector2(euler_base.x+1,euler_base.y),
		Vector2(euler_base.x-1,euler_base.y-1),
		Vector2(euler_base.x,euler_base.y-1),
		Vector2(euler_base.x+1,euler_base.y-1),
		]
		#"""
		
		
		#print(material.particle_lineation[particle].origin,' origin check')
		#print(euler_base,' euler_base check')
		#print()
		for surrounding_node in material.particle_mechanics[particle]['eulerian']:
			if grid_nodes.has(surrounding_node) == false:
				### this grid node is to be affected by the particle.
				grid_nodes[surrounding_node] = {'mass': 0.0,'velocity':Vector2(0.0,0.0),'momentum':Vector2(0.0,0.0),'forces':Vector2(0.0,0.0)}.duplicate(true)
				#print(surrounding_node,' surrounding_node check')
				#print(material.particle_lineation[particle].origin-surrounding_node,' relation check')
			else:
				### already recoqnized...
				pass
		
		
		#material.particle_mechanics[particle]['euler data'] = {'mass': material.mass_in_pieces,'velocity':Vector2(0.0,0.0),'momentum':Vector2(0.0,0.0),'forces':Vector2(0.0,0.0)}.duplicate(true)
		#print(material.particle_mechanics[particle]['euler data'],' euler data check')
		
		# inf , nan check...
		for location in range(0,(len(material.particle_mechanics[particle]['F']))):
			if is_inf(material.particle_mechanics[particle]['F'][location]) == true or is_nan(material.particle_mechanics[particle]['F'][location]) == true:
				if location == 0 or location == 2:
					material.particle_mechanics[particle]['F'][location] = 1.0
				if location == 1 or location == 3:
					material.particle_mechanics[particle]['F'][location] = 0.0
		
		#material.particle_mechanics[particle]['F'] = material.particle_mechanics[particle]['I'].duplicate(true)
		material.particle_mechanics[particle]['J'] = matrix_math.Find_Determinant(material.particle_mechanics[particle]['F'])
		#print(material.particle_mechanics[particle]['J'],' reset J check')
		material.particle_mechanics[particle]['volume'] = material.volume * material.particle_mechanics[particle]['J']
		
		identify_number = wrapi(identify_number+1,0,len(material.particle_mechanics.keys())+1)



#func Particles_to_Grid(time_passed:float,material:Object,the_grid:Dictionary,cell_size:int):
func Particles_to_Grid(time_passed:float,material:Object):
	#particle_interaction = get_tree().get_root().get_node("Simulation/Particle Interaction")
	matrix_math = get_tree().get_root().get_node("Simulation/Matrix Math")
	constitutive_models = get_tree().get_root().get_node("Simulation/Constitutive Models")
	
	# Particles to Grid:
	#for particle in material.particle_mechanics.keys():
	identify_number = 0
	particle = 'null'
	#print()
	#print('Particles_to_Grid check')
	while true:
		if identify_number >= len(material.particle_mechanics.keys()):
			break
		
		particle = material.particle_mechanics.keys()[identify_number]
		
		#print(particle,' particle check')
		#### the sum of aspects that is resetted...
		var energy_density_is_known = false
		var energy_denisty_coeficient = [0,0,0,0]
		affine_momentum_fuse = Vector2(0.0,0.0)
		var sum_of_f_forces = Vector2(0.0,0.0)
		var transfer_momentum = Vector2(0.0,0.0)
		var weighted_mass = 0.0
		var stored_relation = []
		var stored_gradient = []
		var stored_node = []
		var stored_affine = []
		var number = 0 
		var different_number = 0
		
		var differentiate_transposed_energy  = [0,0,0,0]
		var q_differentiated_stressed_summed = [0,0,0,0]
		
		var affline_force_contribution = Vector2(0,0)
		#print(' ')
		#for other_particle in material.particle_mechanics.keys():
		
		#print(material.particle_mechanics[particle]['within_range'],' testing')
		
		#constitutive models...
		if material.physical_state == 'solid':
			if material.constitutive_model == 'hyperelastic':
				energy_denisty_coeficient = constitutive_models.Neo_Hookean(particle,material)
				energy_density_is_known = true
			if material.constitutive_model == 'fixed_corated':
				energy_denisty_coeficient = constitutive_models.Fixed_Corotated(particle,material)
				energy_density_is_known = true
				updating_plasticity = true
			if material.constitutive_model == 'drucker_prager_elasticity':
				energy_denisty_coeficient = constitutive_models.Drucker_Prager_Elasticity(particle,material)
				energy_density_is_known = true
				updating_plasticity = true
		elif material.physical_state == 'liquid':
			if material.constitutive_model == 'water':
				energy_denisty_coeficient = constitutive_models.Model_of_Water(particle,material)
				energy_density_is_known = true
					
		elif material.physical_state == 'gas':
			energy_density_is_known = true
			pass
			
			
		if energy_density_is_known == true:
			var volumed_energy = matrix_math.Multiply_Matrix_by_Scalar(energy_denisty_coeficient,material.volume,true)
			var F_transposed = matrix_math.Transposed_Matrix(material.particle_mechanics[particle]['F'])
			
			var transposed_volumed_energy = matrix_math.Multiply_Matrix(volumed_energy,F_transposed)
			
			var basis_coefficient = 4.0 / pow(grid_spacing,2)
			#print(basis_coefficient,' basis_coefficient check')
			var inverse_I_coefficient = matrix_math.Inverse_Matrix(material.particle_mechanics[particle]['I'])
			var differentiate_weight_interpolation = matrix_math.Multiply_Matrix_by_Scalar(inverse_I_coefficient,basis_coefficient,true)
			
			differentiate_transposed_energy = matrix_math.Multiply_Matrix(transposed_volumed_energy,differentiate_weight_interpolation)
			differentiate_transposed_energy = matrix_math.Multiply_Matrix_by_Scalar(differentiate_transposed_energy,-1.0,true)
			
			var q_mass_C = matrix_math.Multiply_Matrix_by_Scalar(material.particle_mechanics[particle]['C'],material.mass,true)
			
			affline_force_contribution = matrix_math.Add_Matrix(differentiate_transposed_energy,q_mass_C) 
			#print(affline_force_contribution,' affine check')
		else:
			### MLS MPM:  fuses the scattering of the affine momentum and particle force contribution
			# weight_interpolation * (Q)* (relation_between_particle_grid_node)
			#Q = time * D^-1 * particle initial volume * differentiate(energy density)/differentiate(F) * F * F_transposed + mass * C (original form)
			#  or ; if energy density is unknown.
			#Q = - sum of ( time * volume of material occupied at p * stess * (differentiate)weight_interpolation) + mass * C
			#(differentiate)weight_interpolation = weight_interpolation * D^-1 * (relation_between_particle_grid_node)
			# D^-1 = (4/pow(cell_size,2)) * particle.I^-1
		
			var q_timed_volumed = material.particle_mechanics[particle]['volume'] * time_passed
			var q_volumed_stress =  matrix_math.Multiply_Matrix_by_Scalar( material.particle_mechanics[particle]['stress'],q_timed_volumed,true)
			
			var basis_coefficient = 4.0 / pow(grid_spacing,2)
			#print(basis_coefficient,' basis_coefficient check')
			var inverse_I_coefficient = matrix_math.Inverse_Matrix(material.particle_mechanics[particle]['I'])
			var differentiate_weight_interpolation = matrix_math.Multiply_Matrix_by_Scalar(inverse_I_coefficient,basis_coefficient,true)
			
			var q_differentiated_stresse_coefficient = matrix_math.Multiply_Matrix(q_volumed_stress,differentiate_weight_interpolation)
			
			q_differentiated_stressed_summed = matrix_math.Add_Matrix(q_differentiated_stressed_summed,q_differentiated_stresse_coefficient)
			q_differentiated_stressed_summed = matrix_math.Multiply_Matrix_by_Scalar(q_differentiated_stressed_summed,-1.0,true)
			
			var q_mass_C = matrix_math.Multiply_Matrix_by_Scalar(material.particle_mechanics[particle]['C'],material.mass,true)
			
			affline_force_contribution = matrix_math.Add_Matrix(q_differentiated_stressed_summed,q_mass_C)
		#print(affline_force_contribution,' affine check')
			
		
		for node in material.particle_mechanics[particle]['eulerian']:
			#print()
			#print(different_number,' number check')
			
			relation_between_particle_and_grid_node = Vector2(snapped(material.particle_lineation[particle].origin.x - node.x,.01),snapped(material.particle_lineation[particle].origin.y - node.y,.01))
			relation_between_grid_node_and_particle = Vector2(snapped(node.x - material.particle_lineation[particle].origin.x,.01),snapped(node.y - material.particle_lineation[particle].origin.y,.01))
			
			var updated_relation = (1.0/grid_spacing) * relation_between_particle_and_grid_node
					
			#print(relation_between_particle_and_grid_node,' relation_between_particle_and_grid_node check')
			stored_relation.append(relation_between_particle_and_grid_node)
		
			Weight_Interpolation(basis,basis_function_version,updated_relation)#,grid_spacing)
			#print(weight_interpolation['weight'],' weight_interpolation check' )
			
			### MPM lumped mass
			#grid_nodes[node]['mass'] = grid_nodes[node]['mass'] + (material.mass * weight_interpolation['weight'])
				
			### MLS MPM lumped mass
			#print(material.mass, ' material.mass check')
			#print(material.particle_mechanics[particle]['euler data']['mass'], ' [euler data][mass] check')
			#print(weight_interpolation['weight'], ' weight_interpolation[weight] check')
			grid_nodes[node]['mass'] = grid_nodes[node]['mass'] + (material.mass * weight_interpolation['weight'])
			
			stored_node.append(node)
			var affine = matrix_math.Multiply_Matrix_by_Vector2_to_Vector2(affline_force_contribution,relation_between_grid_node_and_particle)
			
			#print(material.particle_mechanics[particle]['initial velocity'], ' initial velocity check')
			#print(grid_nodes[node]['mass'], ' mass check')
			#print(affine, ' affine check')
			
			grid_nodes[node]['momentum'] = grid_nodes[node]['momentum'] + (grid_nodes[node]['mass'] * (material.particle_mechanics[particle]['initial velocity'] + affine))
			
			#different_number = different_number + 1
			
		#while number < len(material.particle_mechanics[particle]['eulerian']):
			### MPM momentum.
			#material.particle_mechanics[particle]['euler data']['momentum'] = material.particle_mechanics[particle]['euler data']['momentum'] + (material.particle_mechanics[particle]['euler data']['mass'] * ( material.initial_velocity + ( matrix_math.Multiply_Matrix_by_Vector2_to_Vector2(material.particle_mechanics[particle]['C'],stored_relation[number]))))
				
			### MLS MPM momentum
			#var C_velocity = material.particle_mechanics[particle]['velocity'] + matrix_math.Multiply_Matrix_by_Vector2_to_Vector2(material.particle_mechanics[particle]['C'],stored_relation[number])
			#var C_velocity = material.initial_velocity + matrix_math.Multiply_Matrix_by_Vector2_to_Vector2(material.particle_mechanics[particle]['C'],stored_relation[number])
				
			#material.particle_mechanics[particle]['euler data']['momentum'] = material.particle_mechanics[particle]['euler data']['momentum'] + (C_velocity * material.particle_mechanics[particle]['euler data']['mass'])
			#print(material.particle_mechanics[particle]['euler data']['momentum'],' [euler data][momentum] check')
			
			### MLS MPM thru affine momentum and particle force contribution
			#grid_nodes[stored_node[number]]['momentum'] = grid_nodes[stored_node[number]]['momentum'] + (grid_nodes[stored_node[number]]['mass'] * (material.initial_velocity + (stored_affine[number]*stored_relation[number])))
			
			#number = number + 1
		
		identify_number = wrapi(identify_number+1,0,len(material.particle_mechanics.keys())+1)

#func Grid_Update(material:Object,the_grid:Dictionary):#,outside_forces:Vector2):
func Grid_Update(time_passed:float,material:Object):
	#print()
	#print('Grid_Update check')
	#Grid Update:
	#print()
	for node in grid_nodes.keys():
	
		if grid_nodes[node]['mass'] > 0.0:
			#print(grid_nodes[node], ' [momentum] check')
			#print(grid_nodes[node]['mass'], ' [mass] check')
			#print(time_passed, ' time_passed check')
			#print(gravity, ' gravity check')
			#print((time_passed *( grid_nodes[node]['mass'] * gravity)), ' (time_passed * grid_nodes[node][mass] * gravity) check')
			grid_nodes[node]['momentum'] = grid_nodes[node]['momentum'] + (time_passed * (grid_nodes[node]['mass'] * gravity))
			grid_nodes[node]['velocity'] = grid_nodes[node]['momentum'] / grid_nodes[node]['mass']
			#print(node,' ',grid_nodes[node], ' node update check')
	
#func Collision_Detection(material:Object,the_grid:Dictionary):
#func Collision_with_Wall(material:Object,the_grid:Dictionary):
func Collision_with_Wall(material:Object):
	#print('Collision_with_Wall check')
	particle_interaction = get_tree().get_root().get_node("Simulation/Particle Interaction")
	#matrix_math = get_tree().get_root().get_node("Simulation/Matrix Math")
	#constitutive_models = get_tree().get_root().get_node("Simulation/Constitutive Models")
	#"""
	### Collision Detection...
	### how the particles interacts....
	
	identify_number = 0
	node = null
	#print()
	#print('Particles_to_Grid check')
	while true:
		#print(identify_number,' identify_number')
		if identify_number >= len(grid_nodes.keys()):
			break
		
		node = grid_nodes.keys()[identify_number]
	
		if grid_nodes[node]['mass'] > 0.0:
			
			for particle in material.particle_mechanics.keys():
				#identify the particle/s
				if node == material.particle_mechanics[particle]['eulerian'][4]:
					
					### Non Rotating Particles
					#var topleft = Vector2(snapped(material.particle_lineation[particle].origin.x,.01),snapped(material.particle_lineation[particle].origin.y,.01)) + Vector2(-(material.appearance.x/2.0),-(material.appearance.y/2.0))
					#var topright = Vector2(snapped(material.particle_lineation[particle].origin.x,.01),snapped(material.particle_lineation[particle].origin.y,.01)) + Vector2((material.appearance.x/2.0),-(material.appearance.y/2.0))
					#var bottomright = Vector2(snapped(material.particle_lineation[particle].origin.x,.01),snapped(material.particle_lineation[particle].origin.y,.01)) + Vector2((material.appearance.x/2.0),(material.appearance.y/2.0))
					#var bottomleft = Vector2(snapped(material.particle_lineation[particle].origin.x,.01),snapped(material.particle_lineation[particle].origin.y,.01)) + Vector2(-(material.appearance.x/2.0),(material.appearance.y/2.0))
					
					### Rotating Particles...
					var topright = Vector2(0,0)
					var bottomright = Vector2(0,0)
					var bottomleft = Vector2(0,0)
					var topleft = Vector2(0,0)
					
					topright.x = material.particle_lineation[particle].origin.x + ((material.appearance.x /2.0) * cos(rad_to_deg(material.particle_lineation[particle].get_rotation()))) - ((material.appearance.y /2.0) * sin(rad_to_deg(material.particle_lineation[particle].get_rotation())))
					topright.y = material.particle_lineation[particle].origin.y + ((material.appearance.x /2.0) * sin(rad_to_deg(material.particle_lineation[particle].get_rotation()))) + ((material.appearance.y /2.0) * cos(rad_to_deg(material.particle_lineation[particle].get_rotation())))
					
					bottomright.x = material.particle_lineation[particle].origin.x + ((material.appearance.x /2.0) * cos(rad_to_deg(material.particle_lineation[particle].get_rotation()))) + ((material.appearance.y /2.0) * sin(rad_to_deg(material.particle_lineation[particle].get_rotation())))
					bottomright.y =  material.particle_lineation[particle].origin.y + ((material.appearance.x /2.0) * sin(rad_to_deg(material.particle_lineation[particle].get_rotation()))) - ((material.appearance.y /2.0) * cos(rad_to_deg(material.particle_lineation[particle].get_rotation())))
					
					bottomleft.x = material.particle_lineation[particle].origin.x - ((material.appearance.x /2.0) * cos(rad_to_deg(material.particle_lineation[particle].get_rotation()))) + ((material.appearance.y /2.0) * sin(rad_to_deg(material.particle_lineation[particle].get_rotation())))
					bottomleft.y = material.particle_lineation[particle].origin.y - ((material.appearance.x /2.0) * sin(rad_to_deg(material.particle_lineation[particle].get_rotation()))) - ((material.appearance.y /2.0) * cos(rad_to_deg(material.particle_lineation[particle].get_rotation())))
					
					topleft.x = material.particle_lineation[particle].origin.x - ((material.appearance.x /2.0) * cos(rad_to_deg(material.particle_lineation[particle].get_rotation()))) - ((material.appearance.y /2.0) * sin(rad_to_deg(material.particle_lineation[particle].get_rotation())))
					topleft.y = material.particle_lineation[particle].origin.y -((material.appearance.x /2.0) * sin(rad_to_deg(material.particle_lineation[particle].get_rotation()))) + ((material.appearance.y /2.0) * cos(rad_to_deg(material.particle_lineation[particle].get_rotation())))
					
					
					if topright.x > barriers['window outline']['right']['outline'] or bottomright.x > barriers['window outline']['right']['outline'] or bottomleft.x > barriers['window outline']['right']['outline'] or topleft.x > barriers['window outline']['right']['outline'] :
						### the particle braches the right of the window...
						### collision with the wall...
						var collision_results =  particle_interaction.Collision_with_Walls('right',material,particle,material.particle_lineation[particle],barriers,grid_nodes[node],material.cell_size)
						#print(collision_results,' testing collision')
						if barriers['window outline']['right']['coefficient of restitution'] >= 1.0 and material.coefficient_of_restitution >= 1.0:
							material.particle_mechanics[particle]['initial velocity'] = material.particle_mechanics[particle]['initial velocity'] + collision_results
						grid_nodes[node]['velocity'] =  grid_nodes[node]['velocity'] + collision_results
					elif topright.y > barriers['window outline']['bottom']['outline'] or bottomright.y > barriers['window outline']['bottom']['outline'] or bottomleft.y > barriers['window outline']['bottom']['outline'] or topleft.y > barriers['window outline']['bottom']['outline'] :
						### the particle braches the bottom of the window...
						### collision with the wall...
						var collision_results =  particle_interaction.Collision_with_Walls('bottom',material,particle,material.particle_lineation[particle],barriers,grid_nodes[node],material.cell_size)
						#print(collision_results,' bottom testing collision')
						if barriers['window outline']['bottom']['coefficient of restitution'] >= 1.0 and material.coefficient_of_restitution >= 1.0:
							material.particle_mechanics[particle]['initial velocity'] = material.particle_mechanics[particle]['initial velocity'] + collision_results
						grid_nodes[node]['velocity'] =  grid_nodes[node]['velocity'] + collision_results
					elif topright.x < barriers['window outline']['left']['outline'] or bottomright.x < barriers['window outline']['left']['outline'] or bottomleft.x < barriers['window outline']['left']['outline'] or topleft.x < barriers['window outline']['left']['outline'] :
						### the particle braches the left of the window...
						### collision with the wall...
						var collision_results =  particle_interaction.Collision_with_Walls('left',material,particle,material.particle_lineation[particle],barriers,grid_nodes[node],material.cell_size)
						#print(collision_results,' testing collision')
						if barriers['window outline']['bottom']['coefficient of restitution'] >= 1.0 and material.coefficient_of_restitution >= 1.0:
							material.particle_mechanics[particle]['initial velocity'] = material.particle_mechanics[particle]['initial velocity'] + collision_results
						grid_nodes[node]['velocity'] =  grid_nodes[node]['velocity'] + collision_results
					elif topright.y < barriers['window outline']['top']['outline'] or bottomright.y < barriers['window outline']['top']['outline'] or bottomleft.y < barriers['window outline']['top']['outline'] or topleft.y < barriers['window outline']['top']['outline'] :
						### the particle braches the top of the window...
						### collision with the wall...
						var collision_results =  particle_interaction.Collision_with_Walls('top',material,particle,material.particle_lineation[particle],barriers,grid_nodes[node],material.cell_size)
						#print(collision_results,' testing collision')
						if barriers['window outline']['bottom']['coefficient of restitution'] >= 1.0 and material.coefficient_of_restitution >= 1.0:
							material.particle_mechanics[particle]['initial velocity'] = material.particle_mechanics[particle]['initial velocity'] + collision_results
						grid_nodes[node]['velocity'] =  grid_nodes[node]['velocity'] + collision_results
					
		identify_number = wrapi(identify_number+1,0,len(grid_nodes.keys())+1)
	

#func Collision_with_Other_Particles(material:Object,the_grid:Dictionary):
func Collision_with_Other_Particles(material:Object):
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
				#print('itself')
			else:
				### Collides with another particle...
				#print('collision')
				
				### Enforce Particle Density...
				#relation_between_particle_and_grid_node = Vector2(snapped(material.particle_lineation[particle].get_center().x - material.particle_lineation[other_particle].get_center().x,.01),snapped(material.particle_lineation[particle].get_center().y - material.particle_lineation[other_particle].get_center().y,.01))
				#relation_between_particle_and_grid_node = Vector2(snapped(material.particle_lineation[particle].origin.x - material.particle_lineation[other_particle].origin.x,.01),snapped(material.particle_lineation[particle].origin.y - material.particle_lineation[other_particle].origin.y,.01))
				
				#print(' ')
				#print(relation_between_particles,' relation_between_particles')
				material.particle_lineation[particle].origin.x = material.particle_lineation[particle].origin.x - (relation_between_particle_and_grid_node.x/2.0)
				material.particle_lineation[particle].origin.y = material.particle_lineation[particle].origin.y - (relation_between_particle_and_grid_node.y/2.0)
				
				material.particle_lineation[other_particle].origin.x = material.particle_lineation[other_particle].origin.x - (relation_between_particle_and_grid_node.x/2.0)
				material.particle_lineation[other_particle].origin.y = material.particle_lineation[other_particle].origin.y - (relation_between_particle_and_grid_node.y/2.0)
				
				### Collision between particles...
				#the_grid[particle]['velocity'] = the_grid[particle]['velocity'] + particle_interaction.Collision_between_Other_Particles(particle,material,material.particle_lineation[particle],other_particle,material,material.particle_lineation[other_particle],the_grid)
				#the_grid[particle]['velocity'] = particle_interaction.Collision_between_Other_Particles(particle,material,material.particle_lineation[particle],other_particle,material,material.particle_lineation[other_particle],the_grid)
				
				material.particle_mechanics[particle]['euler data']['velocity'] = particle_interaction.Collision_between_Other_Particles(particle,material,material.particle_lineation[particle],other_particle,material,material.particle_lineation[other_particle],material.particle_mechanics[particle]['euler data'],material.particle_mechanics[other_particle]['euler data'])
				
				#collection_of_velocities.append(particle_interaction.Collision_between_Other_Particles(particle,material,material.particle_lineation[particle],other_particle,material,material.particle_lineation[other_particle],the_grid))
		#print('after particle collision')
		#print(particle,' ',material.particle_mechanics[particle]['euler data']['velocity'],' material.particle_mechanics[euler data][velocity] check')
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
		
		material.particle_mechanics[particle]['mass'] = material.mass
		material.particle_mechanics[particle]['velocity'] = Vector2(0.0,0.0)
		material.particle_mechanics[particle]['B'] = [0,0,0,0]
		
		### loop counter...
		identify_number = wrapi(identify_number+1,0,len(material.particle_mechanics.keys())+1)


#func Grid_to_Particle(time_passed:float,material:Object,the_grid:Dictionary):
func Grid_to_Particle(time_passed:float,material:Object):
	#particle_interaction = get_tree().get_root().get_node("Simulation/Particle Interaction")
	matrix_math = get_tree().get_root().get_node("Simulation/Matrix Math")
	#constitutive_models = get_tree().get_root().get_node("Simulation/Constitutive Models")
	#print('Grid_to_Particle check')
	#print()
	#print(len(the_grid),' len(the_grid) g2p')
	
	identify_number = 0
	particle = 'null'
	while true:
		if identify_number >= len(material.particle_mechanics.keys()):
			break
		
		particle = material.particle_mechanics.keys()[identify_number]
	
		sum_of_c_weighted_velocity_matrix = [0,0,0,0]
		var stored_relation = []
		var number = 0
		#print('')
		#print( 'start g2p')
		#print('')
		#print(material.particle_mechanics[particle]['velocity'],' velocity')
		
		for node in material.particle_mechanics[particle]['eulerian']:
			if grid_nodes[node]['mass'] > 0.0:
				relation_between_particle_and_grid_node = Vector2(snapped(material.particle_lineation[particle].origin.x - node.x,.01),snapped(material.particle_lineation[particle].origin.y - node.y,.01))
				relation_between_grid_node_and_particle = Vector2(snapped(node.x - material.particle_lineation[particle].origin.x,.01),snapped(node.y - material.particle_lineation[particle].origin.y,.01))
				
				stored_relation.append(relation_between_particle_and_grid_node)
				
				var updated_relation = (1.0/grid_spacing) * relation_between_particle_and_grid_node
						
				Weight_Interpolation(basis,basis_function_version,updated_relation)#,grid_spacing)
					
				### MPM,MLS MPM paraticle velocity update
				# particle Velocity = sum of ( grid velocity * weight_interpolation )
				#B = sum of (weight interpolation * grid velocity * relation between particles transposed )
					
				material.particle_mechanics[particle]['velocity'] = material.particle_mechanics[particle]['velocity'] + ( grid_nodes[node]['velocity'] * weight_interpolation['weight'])
				#print(relation_between_particle_and_grid_node,' relation_between_particle_and_grid_node check')
				#print(material.particle_mechanics[particle]['velocity'],' material.particle_mechanics[particle][velocity] check')
				var velocity_relation = matrix_math.Multiply_Vector2_by_Vector2_to_Matrix(material.particle_mechanics[particle]['velocity'],false,relation_between_grid_node_and_particle,true)
				#print(material.particle_mechanics[particle]['B'],' B check')
				#print(velocity_relation,' velocity_relation check')
				material.particle_mechanics[particle]['B'] = matrix_math.Add_Matrix(material.particle_mechanics[particle]['B'],velocity_relation)
					
		
		#while number < len(material.particle_mechanics[particle]['eulerian']):
			# MPM B, MLS MPM
			#material.particle_mechanics[particle]['B'] = material.particle_mechanics[particle]['velocity']
			#print(stored_relation[number],' stored_relation[number] check')
			#var velocity_relation = matrix_math.Multiply_Vector2_by_Vector2_to_Matrix(material.particle_mechanics[particle]['velocity'],false,stored_relation[number],true)
			#print(material.particle_mechanics[particle]['B'],' B check')
			#print(velocity_relation,' velocity_relation check')
			#material.particle_mechanics[particle]['B'] = matrix_math.Add_Matrix(material.particle_mechanics[particle]['B'],velocity_relation)
				
			#number = number + 1
				
		#MLS MPM 
		#particle.C = D^-1 * B
		# B =  sum of the_grid( particle.velocity * flipped_kernel_distance^T )
		# D^-1 = (4/pow(cell_size,2)) * particle.I^-1
		### Compute C
		c_cell_coefficient =  snapped((basis_coefficient / pow(grid_spacing,2.0)),.01)
		c_inverse_I_coefficient = matrix_math.Inverse_Matrix(material.particle_mechanics[particle]['I'])
		c_cell_inverse_coefficient = matrix_math.Multiply_Matrix_by_Scalar(c_inverse_I_coefficient,c_cell_coefficient,true)
		
		#print(material.particle_mechanics[particle]['B'],' B check')
		#print(c_cell_inverse_coefficient,' c_cell_inverse_coefficient check')
		material.particle_mechanics[particle]['C'] = matrix_math.Multiply_Matrix(c_cell_inverse_coefficient,material.particle_mechanics[particle]['B'])
		#print(material.particle_mechanics[particle]['C'],' C check')
		
		
		###particle position update...
		#print(time_passed,' time_passed check')
		#print(material.particle_mechanics[particle]['velocity'],' velocity as particle check')
		material.particle_lineation[particle].origin.x = snapped((material.particle_lineation[particle].origin.x + snapped((time_passed * material.particle_mechanics[particle]['velocity'].x),.01)  ),.01)
		material.particle_lineation[particle].origin.y = snapped((material.particle_lineation[particle].origin.y + snapped((time_passed * material.particle_mechanics[particle]['velocity'].y),.01)  ),.01)
		
		#print(material.particle_lineation[particle].origin,' position check')
		
		
		### MLS-deformation update...
		#particle.F = (particle.I + time_passed * particle.C ) * particle.F
		
		f_term_1 = matrix_math.Multiply_Matrix_by_Scalar(material.particle_mechanics[particle]['C'],time_passed,true)
		#print(f_term_1,' time_passed * particle.C')
		f_term_2 =  matrix_math.Add_Matrix(material.particle_mechanics[particle]['I'],f_term_1)
		#print(f_term_2,' particle.I + f term 1')
		#print(material.particle_mechanics[particle]['F'],' F before final multi check')
		material.particle_mechanics[particle]['F'] = matrix_math.Multiply_Matrix(f_term_2,material.particle_mechanics[particle]['F'])
		#print(material.particle_mechanics[particle]['F'],' F Update check')
		
		"""
		var test_FE = matrix_math.Multiply_Matrix( material.particle_mechanics[particle]['F'],material.particle_mechanics[particle]['I'])
		
		var inverse_FE = matrix_math.Inverse_Matrix(test_FE)
		var test_FP = matrix_math.Multiply_Matrix(test_FE,material.particle_mechanics[particle]['F'])
		
		var test_F = matrix_math.Multiply_Matrix(test_FE,test_FP)
		print(test_F,' F from E and P check')
		#"""
		
		
		if updating_plasticity == true:
			### Updating Plasticity...
			constitutive_models.Update_Plasticity(material,particle,material.type_of_substance)
			#print(material.particle_mechanics[particle]['F'],' F after P check')
			updating_plasticity = false
		
		
		
		
		
		#print(material.particle_mechanics[particle]['F'],' F check')
		var f_t = matrix_math.Transposed_Matrix(material.particle_mechanics[particle]['F'])
		#print(f_t,' f_t check')
		#var u = matrix_math.Multiply_Matrix(f_t,material.particle_mechanics[particle]['F'])
		#var dia_u = [sqrt(u[0]),0,0,sqrt(u[3])]
		#print(dia_u,' dia_u check')
		#var inv_u =  matrix_math.Inverse_Matrix(dia_u)
		#print(inv_u,' inv_u check')
		#var r = matrix_math.Multiply_Matrix(material.particle_mechanics[particle]['F'],inv_u)
		
		#print(r,' r check')
		
		#"""
		#particle rotation...
		var v = matrix_math.Multiply_Matrix(material.particle_mechanics[particle]['F'],f_t)
		#print(v,' v check')
		var principal_orientation = snapped(((2*v[1]) / (v[0]- v[3])),.01)
		if is_inf(principal_orientation) == true or is_nan(principal_orientation) == true:
			principal_orientation = 0
		#print(principal_orientation,' principal_orientation check')
		
		material.particle_lineation[particle].x.x = cos(deg_to_rad(principal_orientation))
		material.particle_lineation[particle].y.x = -sin(deg_to_rad(principal_orientation))
		material.particle_lineation[particle].x.y = sin(deg_to_rad(principal_orientation))
		material.particle_lineation[particle].y.y = cos(deg_to_rad(principal_orientation))
		#"""
		
		
		
		
		identify_number = wrapi(identify_number+1,0,len(material.particle_mechanics.keys())+1)
	#"""
