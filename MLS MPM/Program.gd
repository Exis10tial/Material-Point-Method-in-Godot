extends Node


var particle_interaction
var matrix_math
var constitutive_models


### Eulerian/Cartesian Grid...
#var Grid_Nodes : Dictionary = {}
var grid_nodes : Dictionary = {}
var station : Vector2
#var identify_grid_reach : Array = []
#var gravity : Vector2
var gravity : Vector2 = Vector2(0.00,9.80)# * 100
#var gravity : Vector2 = Vector2(0.00,-9.80#)* 100
#var gravity : Vector2 = Vector2(9.80,0.0)# * 100
#var gravity : Vector2 = Vector2(-9.80,0.0)#* 100
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
		#barriers['window outline'][wall]['coefficient of static friction'] = 1.0
		barriers['window outline'][wall]['coefficient of static friction'] = 0.80
		#barriers['window outline'][wall]['coefficient of static friction'] = 0.0
		#barriers['window outline'][wall]['coefficient of kinetic friction'] = 1.0
		barriers['window outline'][wall]['coefficient of kinetic friction'] = .650
		#barriers['window outline'][wall]['coefficient of kinetic friction'] = 0.0
		barriers['window outline'][wall]['mass'] = 1.00
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


func Weight_Interpolation(splines:String,version:int,relative:Vector2,cell_size:float):
	#reset parameters
	building_weights = Vector2(0.0,0.0)
	### Weight Interpolation...
	
	if splines == "quadratic":
		#print('weight interpolation')
		if version == 1:
			### The Material Point Method for Simulating Continuum Materials pg 33 equation 123
			#if int(relative.x) >= (0*cell_size) and int(relative.x) <= (1*cell_size):
			if abs(relative.x) == (0*cell_size) and abs(relative.x) <= ((cell_size/4.0)):
				#print('1')
				building_weights.x = snappedf((.75 - snappedf(pow(abs(relative.x),2.0),.01)),.01)
				#print(building_weights,' a building_weights')
			#elif int(relative.x) >= (-.50)*(cell_size) and int(relative.x) <= (.50)*(cell_size):
			elif abs(relative.x) > ((cell_size/4.0)) and abs(relative.x) <= ((cell_size/2.0)):
				#print('2')
				building_weights.x = .50 *  pow((1.50 - abs(relative.x)),2.0)
				
			#elif int(relative.x) >= (.50)*(cell_size) and int(relative.x) <= (1.50)*(cell_size):
			elif abs(relative.x) > ((cell_size/2.0)):# and int(relative.x) <= (cell_size/2.0):
				### particles is not in range of the particle to influence.
				#print('3')
				building_weights.x =  0.0
				
			#	pass
			#else:
				
				#building_weights.x = 0.0
			#if int(relative.y) >= (-1.50)*(cell_size) and int(relative.y) <= (-.50)*(cell_size):
			#if int(relative.y) >= (0*cell_size) and int(relative.x) <= (1*cell_size):
			if abs(relative.y) == (0*cell_size) and abs(relative.y) <= ((cell_size/4.0)):
				#print('6')
				building_weights.y = snappedf((.75 - snappedf(pow(abs(relative.y),2.0),.01)),.01)
				#print(building_weights,' a building_weights')
			#elif int(relative.y) >= (-.50)*(cell_size) and int(relative.y) <= (.50)*(cell_size):
			#elif int(relative.y) >= (1*cell_size) and int(relative.x) <= (2*cell_size):
			elif abs(relative.y) > ((cell_size/4.0)) and abs(relative.y) <= ((cell_size/2.0)):
				#print('7')
				building_weights.y =.50 * pow((1.50 - abs(relative.y) ),2.0)
				
			#elif int(relative.y) >= (.50)*(cell_size) and int(relative.y) <= (1.50)*(cell_size):
			#elif int(relative.y) >= (2*cell_size):#-(cell_size/4.0) and int(relative.y) <= (cell_size/2.0):
			elif abs(relative.y) > ((cell_size/2.0)):
				#print('8')
				building_weights.y = 0.0
				
			#	pass
			#else:
				
				#building_weights.y = 0.0
		elif version == 2:
			### thru Analysis and Reduction of Quadrature Errors in the Material Point Method. steffan .W.pdf-pg(8.equation 16)
			
			#if int(relative.x) >= (-1.50)*(cell_size) and int(relative.x) <= (-.50)*(cell_size):
			if relative.x >= -(cell_size/2.0) and relative.x <= -(cell_size/4.0):
				
				building_weights.x = snapped( ((1.0 * pow(relative.x,2.0)) /  (2.0 * pow(cell_size,2.0))),.01 ) + snapped(((3.0 * relative.x) / (2.0 * cell_size)),.01) + snapped((9.0/8.0),.001)
				
			#elif int(relative.x) >= (-.50)*(cell_size) and int(relative.x) <= (.50)*(cell_size):
			elif relative.x >= -(cell_size/4.0) and relative.x <= (cell_size/4.0):
				
				building_weights.x = snapped(((-1.0 * pow(relative.x,2.0)) / (pow(cell_size,2.0))),.01) + snapped((3.0/4.0),.01)
				
			#elif int(relative.x) >= (.50)*(cell_size) and int(relative.x) <= (1.50)*(cell_size):
			elif relative.x >= -(cell_size/4.0) and relative.x <= (cell_size/2.0):
				
				building_weights.x =  snapped( ((1.0 * pow(relative.x,2.0)) /  (2.0 * pow(cell_size,2.0))),.01 ) - snapped(((3.0 * relative.x) / (2.0 * cell_size)),.01) + snapped((9.0/8.0),.001)
				
			else:
				
				building_weights.x = 0.0
			
			#if int(relative.y) >= (-1.50)*(cell_size) and int(relative.y) <= (-.50)*(cell_size):
			if relative.y >= -(cell_size/2.0) and relative.y <= -(cell_size/4.0):
				
				building_weights.y = snapped( ((1.0 * pow(relative.x,2.0)) /  (2.0 * pow(cell_size,2.0))),.01 ) + snapped(((3.0 * relative.y) / (2.0 * cell_size)),.01) + snapped((9.0/8.0),.001)
				
			#elif int(relative.y) >= (-.50)*(cell_size) and int(relative.y) <= (.50)*(cell_size):
			elif relative.y >= -(cell_size/4.0) and relative.y <= (cell_size/4.0):
				
				building_weights.y = snapped(((-1.0 * pow(relative.y,2.0)) / (pow(cell_size,2.0))),.01) + snapped((3.0/4.0),.01)
				
			#elif int(relative.y) >= (.50)*(cell_size) and int(relative.y) <= (1.50)*(cell_size):
			elif relative.y >= -(cell_size/4.0) and relative.y <= (cell_size/2.0):
				
				building_weights.y =  snapped( ((1.0 * pow(relative.y,2.0)) /  (2.0 * pow(cell_size,2.0))),.01 ) - snapped(((3.0 * relative.y) / (2.0 * cell_size)),.01) + snapped((9.0/8.0),.001)
				
			else:
				
				building_weights.y = 0.0
			
		###...
		#print(building_weights.x,' building weights x')
		#print(building_weights.y,' building weights y')
		#weight_interpolation['weight'] = ((1.0/cell_size)*building_weights.x ) * ((1.0/cell_size)*building_weights.y)
		weight_interpolation['weight'] =(building_weights.x/cell_size) * (building_weights.y/cell_size)
	
		### unneeded in mls mpm cause of the affine and force combination...
		#weight_interpolation['gradient'] =((1.0/cell_size) * ((1.0/cell_size)*building_weights.y)) * (((1.0/cell_size)*building_weights.x ) * (1.0/cell_size))
		weight_interpolation['gradient'] = Vector2((( (1.0/cell_size) * ((1.0/cell_size)*building_weights.x ) ) * ((1.0/cell_size)*building_weights.y)),   ((1.0/cell_size)*building_weights.x ) *  ((1.0/cell_size) * ((1.0/cell_size)*building_weights.y)))
		
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
	#print('Grid_Reset check')
	
	### the grid is renewed...
	#grid_nodes = the_grid.duplicate(true)
	#identify_grid_reach = []
	#print(len(grid_nodes),' len(grid_nodes)')
	
	#forces = [0,0,0,0]
	
	while true:
		if identify_number >= len(material.particle_mechanics.keys()):
			break
		
		particle = material.particle_mechanics.keys()[identify_number]
		
		#material.particle_mechanics[particle]['eulerian'] =[Vector2(snapped(material.particle_lineation[particle].position.x,.01),snapped(material.particle_lineation[particle].position.y,.01)),
		#	Vector2(snapped(material.particle_lineation[particle].end.x,.01),snapped(material.particle_lineation[particle].position.y,.01)),
		#	Vector2(snapped(material.particle_lineation[particle].end.x,.01),snapped(material.particle_lineation[particle].end.y,.01)),
		#	Vector2(snapped(material.particle_lineation[particle].position.x,.01),snapped(material.particle_lineation[particle].end.y,.01))
		#	].duplicate(true)
		#print(material.particle_lineation[particle].origin.x,' transform origin check')
		#material.particle_mechanics[particle]['eulerian'] =[ Vector2(snapped(material.particle_lineation[particle].origin.x,.01),snapped(material.particle_lineation[particle].origin.y,.01))
		
		### the relation between the particle with itself and other particles...
		#var topleft = Vector2(snapped(material.particle_lineation[particle].origin.x,.01),snapped(material.particle_lineation[particle].origin.y,.01)) + Vector2(-(material.appearance.x/2.0),-(material.appearance.y/2.0))
		#var topright = Vector2(snapped(material.particle_lineation[particle].origin.x,.01),snapped(material.particle_lineation[particle].origin.y,.01)) + Vector2((material.appearance.x/2.0),-(material.appearance.y/2.0))
		#var bottomleft = Vector2(snapped(material.particle_lineation[particle].origin.x,.01),snapped(material.particle_lineation[particle].origin.y,.01)) + Vector2(-(material.appearance.x/2.0),(material.appearance.y/2.0))
		#var bottomright = Vector2(snapped(material.particle_lineation[particle].origin.x,.01),snapped(material.particle_lineation[particle].origin.y,.01)) + Vector2((material.appearance.x/2.0),(material.appearance.y/2.0))
		#material.particle_mechanics[particle]['eulerian'] = [topleft,topright,bottomleft,bottomright]
		
		### the relation between the particle with itself and other particles...
		var center = Vector2(material.particle_lineation[particle].origin.x,material.particle_lineation[particle].origin.y)
		material.particle_mechanics[particle]['eulerian'] = [center]
		
		
		#print(material.particle_mechanics[particle]['eulerian'],' eulerian check')
		
		material.particle_mechanics[particle]['euler data'] = {'mass': 0.0,'velocity':Vector2(0.0,0.0),'momentum':Vector2(0.0,0.0),'forces':Vector2(0.0,0.0)}.duplicate(true)
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
		material.particle_mechanics[particle]['volume'] = material.volume_in_pieces * material.particle_mechanics[particle]['J']
		
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
		var stored_relation = []
		var stored_gradient = []
		var stored_affine = []
		var number = 0 
		var different_number = 0
		
		
		var q_stressed_volumed_mass_C
		var affline_force_contribution
		#print(' ')
		#for other_particle in material.particle_mechanics.keys():
		
		#print(material.particle_mechanics[particle]['within_range'],' testing')
		
		for other_particle in material.particle_mechanics[particle]['within_range']:
		#print(len(grid_nodes),' len(grid_nodes) p2g')
		#for relation_between_particles in  material.particle_mechanics[particle]['ambit']:
			### Weight Interpolation...
			#print(material.particle_mechanics[other_particle]['eulerian'],' material.particle_mechanics[other_particle][eulerian]')
			stored_relation = []
			material.particle_mechanics[particle]['gradient forces'][other_particle] = []
			stored_gradient = []
			#var number = 0 
			stored_affine = []
			affline_force_contribution = Vector2(0,0)
			
			if material.physical_state == 'solid':
				if material.constitutive_model == 'hyperelastic':
					#material.particle_mechanics[other_particle].stress = constitutive_models.Neo_Hookean(particle,material)
					material.particle_mechanics[particle]['stress'] = constitutive_models.Neo_Hookean(particle,material)
				if material.constitutive_model == 'fixed_corated':
					#material.particle_mechanics[other_particle].stress = constitutive_models.Fixed_Corotated(particle,material)
					material.particle_mechanics[particle]['stress'] = constitutive_models.Fixed_Corotated(particle,material)
				if material.constitutive_model == 'drucker_prager_elasticity':
					#material.particle_mechanics[other_particle].stress = constitutive_models.Drucker_Prager_Elasticity(particle,material)
					material.particle_mechanics[particle]['stress'] = constitutive_models.Drucker_Prager_Elasticity(particle,material)
			elif material.physical_state == 'liquid':
				if material.constitutive_model == 'water':
					material.particle_mechanics[particle]['stress'] = constitutive_models.Model_of_Water(particle,material)
					
			elif material.physical_state == 'gas':
				pass
			
			#print(material.particle_mechanics[particle]['stress'],' stress check')
			for node in material.particle_mechanics[other_particle]['eulerian']:
				
				#relation_between_particles = Vector2(snapped(node.x - material.particle_lineation[particle].get_center().x,.01),snapped(node.y - material.particle_lineation[particle].get_center().y,.01))
				
				relation_between_particles = Vector2(snapped(node.x - material.particle_lineation[particle].origin.x,.01),snapped(node.y - material.particle_lineation[particle].origin.y,.01))
				stored_relation.append(relation_between_particles)
				#relation_between_particles = Vector2(snapped(material.particle_lineation[particle].get_center().x - node.x,.01),snapped(material.particle_lineation[particle].get_center().y -node.y,.01))
				#print('')
				#print(relation_between_particles,' relation_between_particles')
				Weight_Interpolation(basis,basis_function_version,relation_between_particles,material.appearance.x)
				#stored_gradient.append(weight_interpolation['gradient'])
				#"""
				
				#print(weight_interpolation['weight'],' weight_interpolation[weight] P2G')
				#print(material.particle_mechanics[particle]['euler data']['mass'],' material.particle_mechanics[particle][euler data][mass] check')
					
					
				### MPM lumped mass
				#material.particle_mechanics[particle]['euler data']['mass'] = material.particle_mechanics[particle]['euler data']['mass'] + (material.mass * weight_interpolation['weight'])
				
				
				### MLS MPM lumped mass
				#print(material.particle_mechanics[particle]['euler data']['mass'], ' [euler data][mass] check')
				#print(weight_interpolation['weight'], ' weight_interpolation[weight] check')
				material.particle_mechanics[particle]['euler data']['mass'] = material.particle_mechanics[particle]['euler data']['mass'] + (material.mass * weight_interpolation['weight'])
				
				### MLS MPM:  fuses the scattering of the affine momentum and particle force contribution
				# weight_interpolation * (Q)* (relation_between_particle_grid_node)
				#Q = time * D^-1 * initial_volume * differentiate(energy density)/differentiate(F) * F * F_transposed + mass * C (original form)
				#Q = time * volume  D^-1 * stess + mass * C
				# D = (4/pow(cell_size,2)) * particle.I^-1
				var q_timed_volumed = time_passed *  material.volume * basis_coefficient / pow(material.appearance.x,2)
				#print(q_time,' q_time check')
				var q_I = matrix_math.Multiply_Matrix_by_Scalar( material.particle_mechanics[particle]['I'],q_timed_volumed,true)
				#print(q_I,' q_I check')
				#print(q_volumed_I,' q_volumed_I check')
				var q_stressed_volumed_I = matrix_math.Multiply_Matrix(q_I,material.particle_mechanics[particle]['stress'])
				#print(q_stressed_volumed_I,' q_stressed_volumed_I check')
				var q_mass_C = matrix_math.Multiply_Matrix_by_Scalar(material.particle_mechanics[particle]['C'],material.mass,true)
				#print(q_mass_C,' q_mass_C check')
				q_stressed_volumed_mass_C = matrix_math.Add_Matrix(q_stressed_volumed_I,q_mass_C)
				#print(q_stressed_volumed_mass_C,' check')
				var weighted_q_stressed_volumed_mass_C =  matrix_math.Multiply_Matrix_by_Scalar(q_stressed_volumed_mass_C,weight_interpolation['weight'],true)
				affline_force_contribution = matrix_math.Multiply_Matrix_by_Vector2_to_Vector2(weighted_q_stressed_volumed_mass_C,relation_between_particles)
				#print(affline_force_contribution,' affline_force_contribution check')
				
			# used for MPM forces.
			#material.particle_mechanics[particle]['gradient forces'][other_particle] = stored_gradient
			
			#print(material.particle_mechanics[particle]['euler data']['momentum'],' [euler data][momentum] check before')
			### Summation addition...grid nodes
			while number < len(material.particle_mechanics[other_particle]['eulerian']):
				### MPM momentum.
				#material.particle_mechanics[particle]['euler data']['momentum'] = material.particle_mechanics[particle]['euler data']['momentum'] + (material.particle_mechanics[particle]['euler data']['mass'] * ( material.initial_velocity + ( matrix_math.Multiply_Matrix_by_Vector2_to_Vector2(material.particle_mechanics[particle]['C'],stored_relation[number]))))
				
				### MLS MPM momentum
				#var C_velocity = material.particle_mechanics[particle]['velocity'] + matrix_math.Multiply_Matrix_by_Vector2_to_Vector2(material.particle_mechanics[particle]['C'],stored_relation[number])
				#var C_velocity = material.initial_velocity + matrix_math.Multiply_Matrix_by_Vector2_to_Vector2(material.particle_mechanics[particle]['C'],stored_relation[number])
				
				#material.particle_mechanics[particle]['euler data']['momentum'] = material.particle_mechanics[particle]['euler data']['momentum'] + (C_velocity * material.particle_mechanics[particle]['euler data']['mass'])
				#print(material.particle_mechanics[particle]['euler data']['momentum'],' [euler data][momentum] check')
				
				### MLS MPM thru affine momentum and particle force contribution
				material.particle_mechanics[particle]['euler data']['momentum'] = material.particle_mechanics[particle]['euler data']['momentum'] +  (material.initial_velocity + affline_force_contribution)
				#material.particle_mechanics[particle]['euler data']['momentum'] = material.particle_mechanics[particle]['euler data']['momentum'] +  (material.particle_mechanics[particle]['velocity'] + affline_force_contribution)
				
				
				
			### Summation addition......grid nodes
			#while different_number < len(material.particle_mechanics[other_particle]['eulerian']):
				### MPM forces.
				#var stressed_gradient = matrix_math.Multiply_Matrix_by_Vector2_to_Matrix(material.particle_mechanics[particle]['stress'],material.particle_mechanics[particle]['gradient forces'][other_particle][different_number],false,false)
				#var stressed_gradient = matrix_math.Multiply_Matrix_by_Vector2_to_Matrix(material.particle_mechanics[particle]['stress'],stored_gradient[number],false,false)
				
				#print(stressed_gradient,' stressed_gradient check')
				#var volumed_stressed_gradient = matrix_math.Multiply_Matrix_by_Scalar(stressed_gradient,material.particle_mechanics[particle]['volume'],true)
				#var volumed_stressed_gradient_converted =  matrix_math.Convert_to_Vector2(volumed_stressed_gradient)
				
				#material.particle_mechanics[particle]['euler data']['forces'] = material.particle_mechanics[particle]['euler data']['forces'] + volumed_stressed_gradient_converted
				#material.particle_mechanics[particle]['euler data']['forces'] = material.particle_mechanics[particle]['euler data']['forces'] * -1.0
				
				### MLS MPM forces...
				# used if not using Q...
				#material.particle_mechanics[particle]['euler data']['forces'] = material.particle_mechanics[particle]['euler data']['forces'] + matrix_math.Multiply_Matrix_by_Vector2_to_Vector2(stored_affine[number],stored_gradient[number])
				
				number = number + 1
				
				#different_number = different_number + 1
			
			### MPM forces.
			#material.particle_mechanics[particle]['euler data']['forces'] = material.particle_mechanics[particle]['euler data']['forces'] * -1.0
			

		### inf , nan check
		if is_inf(material.particle_mechanics[particle]['euler data']['forces'].x) == true or is_nan(material.particle_mechanics[particle]['euler data']['forces'].y) == true:
			if is_inf(material.particle_mechanics[particle]['euler data']['forces'].x) == true:
				material.particle_mechanics[particle]['euler data']['forces'].x = 1.0
			if is_inf(material.particle_mechanics[particle]['euler data']['forces'].y) == true:
				material.particle_mechanics[particle]['euler data']['forces'].y = 1.0
		if is_nan(material.particle_mechanics[particle]['euler data']['forces'].x) == true or is_nan(material.particle_mechanics[particle]['euler data']['forces'].y) == true:
			if is_nan(material.particle_mechanics[particle]['euler data']['forces'].x) == true:
				material.particle_mechanics[particle]['euler data']['forces'].x = 1.0
			if is_nan(material.particle_mechanics[particle]['euler data']['forces'].y) == true:
				material.particle_mechanics[particle]['euler data']['forces'].y = 1.0
		#print(particle,' ',material.particle_mechanics[particle]['euler data']['forces'],' [euler data][forces] check')
		identify_number = wrapi(identify_number+1,0,len(material.particle_mechanics.keys())+1)
	
#func Grid_Update(material:Object,the_grid:Dictionary):#,outside_forces:Vector2):
func Grid_Update(time_passed:float,material:Object):
	#print('Grid_Update check')
	#Grid Update:
	for particle in material.particle_lineation.keys():
		if material.particle_mechanics[particle]['euler data']['mass'] > 0.0:
			
			### MPM grid update...
			#material.particle_mechanics[particle]['euler data']['velocity'] = material.particle_mechanics[particle]['euler data']['velocity'] + ((time_passed * material.particle_mechanics[particle]['euler data']['forces'])/snapped(material.particle_mechanics[particle]['euler data']['mass'],.01)) #original
			#material.particle_mechanics[particle]['euler data']['velocity'] = material.particle_mechanics[particle]['euler data']['velocity'] + ((time_passed * (material.particle_mechanics[particle]['euler data']['forces'] * gravity) )/snapped(material.particle_mechanics[particle]['euler data']['mass'],.01))
			
			### MLS MPM grid update...
			#print(material.particle_mechanics[particle]['euler data']['momentum'], ' [euler data][momentum] check')
			#print(material.particle_mechanics[particle]['euler data']['mass'], ' [euler data][mass] check')
			#print(material.particle_mechanics[particle]['euler data']['forces'], ' [euler data][forces] check')
			#print(time_passed, ' time_passed check')
			#print(gravity, ' gravity check')
			material.particle_mechanics[particle]['euler data']['momentum'] = material.particle_mechanics[particle]['euler data']['momentum'] + (time_passed * (material.particle_mechanics[particle]['euler data']['mass'] * gravity + material.particle_mechanics[particle]['euler data']['forces']) )
			material.particle_mechanics[particle]['euler data']['velocity'] = material.particle_mechanics[particle]['euler data']['momentum'] / material.particle_mechanics[particle]['euler data']['mass']
			
			#print(material.particle_mechanics[particle]['euler data']['velocity'], ' grid velocity update check')
		else:
			pass
			

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
	collection_of_velocities = []
	for particle in material.particle_lineation.keys():
		### if the particle has contacted the window outline...
		#print(the_grid[particle]['velocity'],' incoming')
		barriers['window outline']['top']['outline'] = 0.0
		barriers['window outline']['right']['outline'] = ProjectSettings.get_setting('display/window/size/width')
		barriers['window outline']['bottom']['outline'] = ProjectSettings.get_setting('display/window/size/height')
		barriers['window outline']['left']['outline'] = 0.0
		
		
		if rad_to_deg(material.particle_lineation[particle].get_rotation()) > 0:
			### if the particle has the slightest rotation...

			#find the position of the mesh corners(square)...
			var topleft_rotated = Vector2(material.particle_lineation[particle].origin.x - ( (material.appearance.x/2.0) * rad_to_deg(cos(material.particle_lineation[particle].get_rotation())) - (material.appearance.y/2.0) * rad_to_deg(sin(material.particle_lineation[particle].get_rotation())) ),
			material.particle_lineation[particle].origin.y - ( (material.appearance.x/2.0) * rad_to_deg(sin(material.particle_lineation[particle].get_rotation())) + (material.appearance.y/2.0) * rad_to_deg(cos(material.particle_lineation[particle].get_rotation())) ) )
			
			var topright_rotated = Vector2(material.particle_lineation[particle].origin.x + ( (material.appearance.x/2.0) * rad_to_deg(cos(material.particle_lineation[particle].get_rotation())) - (material.appearance.y/2.0) * rad_to_deg(sin(material.particle_lineation[particle].get_rotation())) ),
			material.particle_lineation[particle].origin.y + ( (material.appearance.x/2.0) * rad_to_deg(sin(material.particle_lineation[particle].get_rotation())) + (material.appearance.y/2.0) * rad_to_deg(cos(material.particle_lineation[particle].get_rotation())) ) )
			
			var bottomright_rotated = Vector2(material.particle_lineation[particle].origin.x + ( (material.appearance.x/2.0) * rad_to_deg(cos(material.particle_lineation[particle].get_rotation())) + (material.appearance.y/2.0) * rad_to_deg(sin(material.particle_lineation[particle].get_rotation())) ),
			material.particle_lineation[particle].origin.y + ( (material.appearance.x/2.0) * rad_to_deg(sin(material.particle_lineation[particle].get_rotation())) - (material.appearance.y/2.0) * rad_to_deg(cos(material.particle_lineation[particle].get_rotation())) ) )
			
			var bottomleft_rotated = Vector2(material.particle_lineation[particle].origin.x - ( (material.appearance.x/2.0) * rad_to_deg(cos(material.particle_lineation[particle].get_rotation())) + (material.appearance.y/2.0) * rad_to_deg(sin(material.particle_lineation[particle].get_rotation())) ),
			material.particle_lineation[particle].origin.y - ( (material.appearance.x/2.0) * rad_to_deg(sin(material.particle_lineation[particle].get_rotation())) - (material.appearance.y/2.0) * rad_to_deg(cos(material.particle_lineation[particle].get_rotation())) ) )
			
			# finding which is has past a barrier
			if topleft_rotated.y < barriers['window outline']['top']['outline'] or topright_rotated.y < barriers['window outline']['top']['outline'] or bottomright_rotated.y < barriers['window outline']['top']['outline'] or bottomleft_rotated.y < barriers['window outline']['top']['outline']:
				if topleft_rotated.y < barriers['window outline']['top']['outline']:
					### collision with the wall...
					material.particle_mechanics[particle]['euler data']['velocity'] = particle_interaction.Collision_with_Walls('top',material,particle,material.particle_lineation[particle],barriers,material.particle_mechanics[particle]['euler data'],material.cell_size)
				elif topright_rotated.y < barriers['window outline']['right']['outline']:
					### collision with the wall...
					material.particle_mechanics[particle]['euler data']['velocity'] = particle_interaction.Collision_with_Walls('right',material,particle,material.particle_lineation[particle],barriers,material.particle_mechanics[particle]['euler data'],material.cell_size)
				elif bottomright_rotated.y < barriers['window outline']['bottom']['outline']:
					### collision with the wall...
					material.particle_mechanics[particle]['euler data']['velocity'] = particle_interaction.Collision_with_Walls('bottom',material,particle,material.particle_lineation[particle],barriers,material.particle_mechanics[particle]['euler data'],material.cell_size)
				elif bottomleft_rotated.y < barriers['window outline']['left']['outline']:
					### collision with the wall...
					material.particle_mechanics[particle]['euler data']['velocity'] = particle_interaction.Collision_with_Walls('left',material,particle,material.particle_lineation[particle],barriers,material.particle_mechanics[particle]['euler data'],material.cell_size)
			elif topleft_rotated.x > barriers['window outline']['right']['outline'] or topright_rotated.x > barriers['window outline']['right']['outline'] or bottomright_rotated.x > barriers['window outline']['right']['outline'] or bottomleft_rotated.x > barriers['window outline']['right']['outline']:
				if topleft_rotated.x > barriers['window outline']['top']['outline']:
					### collision with the wall...
					material.particle_mechanics[particle]['euler data']['velocity'] = particle_interaction.Collision_with_Walls('top',material,particle,material.particle_lineation[particle],barriers,material.particle_mechanics[particle]['euler data'],material.cell_size)
				elif topright_rotated.x > barriers['window outline']['right']['outline']:
					### collision with the wall...
					material.particle_mechanics[particle]['euler data']['velocity'] = particle_interaction.Collision_with_Walls('right',material,particle,material.particle_lineation[particle],barriers,material.particle_mechanics[particle]['euler data'],material.cell_size)
				elif bottomright_rotated.x > barriers['window outline']['bottom']['outline']:
					### collision with the wall...
					material.particle_mechanics[particle]['euler data']['velocity'] = particle_interaction.Collision_with_Walls('bottom',material,particle,material.particle_lineation[particle],barriers,material.particle_mechanics[particle]['euler data'],material.cell_size)
				elif bottomleft_rotated.x > barriers['window outline']['left']['outline']:
					### collision with the wall...
					material.particle_mechanics[particle]['euler data']['velocity'] = particle_interaction.Collision_with_Walls('left',material,particle,material.particle_lineation[particle],barriers,material.particle_mechanics[particle]['euler data'],material.cell_size)
			elif topleft_rotated.y > barriers['window outline']['bottom']['outline'] or topright_rotated.y > barriers['window outline']['bottom']['outline'] or bottomright_rotated.y > barriers['window outline']['bottom']['outline'] or bottomleft_rotated.y > barriers['window outline']['bottom']['outline']:
				if topleft_rotated.y > barriers['window outline']['top']['outline']:
					### collision with the wall...
					material.particle_mechanics[particle]['euler data']['velocity'] = particle_interaction.Collision_with_Walls('top',material,particle,material.particle_lineation[particle],barriers,material.particle_mechanics[particle]['euler data'],material.cell_size)
				elif topright_rotated.y > barriers['window outline']['right']['outline']:
					### collision with the wall...
					material.particle_mechanics[particle]['euler data']['velocity'] = particle_interaction.Collision_with_Walls('right',material,particle,material.particle_lineation[particle],barriers,material.particle_mechanics[particle]['euler data'],material.cell_size)
				elif bottomright_rotated.y > barriers['window outline']['bottom']['outline']:
					### collision with the wall...
					material.particle_mechanics[particle]['euler data']['velocity'] = particle_interaction.Collision_with_Walls('bottom',material,particle,material.particle_lineation[particle],barriers,material.particle_mechanics[particle]['euler data'],material.cell_size)
				elif bottomleft_rotated.y > barriers['window outline']['left']['outline']:
					### collision with the wall...
					material.particle_mechanics[particle]['euler data']['velocity'] = particle_interaction.Collision_with_Walls('left',material,particle,material.particle_lineation[particle],barriers,material.particle_mechanics[particle]['euler data'],material.cell_size)
			elif topleft_rotated.x < barriers['window outline']['left']['outline'] or topright_rotated.x < barriers['window outline']['left']['outline'] or bottomright_rotated.x < barriers['window outline']['left']['outline'] or bottomleft_rotated.x < barriers['window outline']['left']['outline']:
				if topleft_rotated.x < barriers['window outline']['top']['outline']:
					### collision with the wall...
					material.particle_mechanics[particle]['euler data']['velocity'] = particle_interaction.Collision_with_Walls('top',material,particle,material.particle_lineation[particle],barriers,material.particle_mechanics[particle]['euler data'],material.cell_size)
				elif topright_rotated.x < barriers['window outline']['right']['outline']:
					### collision with the wall...
					material.particle_mechanics[particle]['euler data']['velocity'] = particle_interaction.Collision_with_Walls('right',material,particle,material.particle_lineation[particle],barriers,material.particle_mechanics[particle]['euler data'],material.cell_size)
				elif bottomright_rotated.x < barriers['window outline']['bottom']['outline']:
					### collision with the wall...
					material.particle_mechanics[particle]['euler data']['velocity'] = particle_interaction.Collision_with_Walls('bottom',material,particle,material.particle_lineation[particle],barriers,material.particle_mechanics[particle]['euler data'],material.cell_size)
				elif bottomleft_rotated.x < barriers['window outline']['left']['outline']:
					### collision with the wall...
					material.particle_mechanics[particle]['euler data']['velocity'] = particle_interaction.Collision_with_Walls('left',material,particle,material.particle_lineation[particle],barriers,material.particle_mechanics[particle]['euler data'],material.cell_size)
				
			
		elif material.particle_lineation[particle].get_rotation() == 0:
			### if the particle has no rotation...
			
			var topleft = Vector2(snapped(material.particle_lineation[particle].origin.x,.01),snapped(material.particle_lineation[particle].origin.y,.01)) + Vector2(-(material.appearance.x/2.0),-(material.appearance.y/2.0))
			var topright = Vector2(snapped(material.particle_lineation[particle].origin.x,.01),snapped(material.particle_lineation[particle].origin.y,.01)) + Vector2((material.appearance.x/2.0),-(material.appearance.y/2.0))
			var bottomright = Vector2(snapped(material.particle_lineation[particle].origin.x,.01),snapped(material.particle_lineation[particle].origin.y,.01)) + Vector2((material.appearance.x/2.0),(material.appearance.y/2.0))
			var bottomleft = Vector2(snapped(material.particle_lineation[particle].origin.x,.01),snapped(material.particle_lineation[particle].origin.y,.01)) + Vector2(-(material.appearance.x/2.0),(material.appearance.y/2.0))
			
			#print('before collision')
			#print(particle,' ',material.particle_mechanics[particle]['euler data']['velocity'],' material.particle_mechanics[euler data][velocity] check')
			
			if topleft.y < barriers['window outline']['top']['outline']:
				### the particle breaches the top the window...
				### collision with the wall...
				#material.particle_mechanics[particle]['euler data']['velocity'] = particle_interaction.Collision_with_Walls('top',material,particle,material.particle_lineation[particle],barriers,material.particle_mechanics[particle]['euler data'],material.cell_size)
				material.particle_mechanics[particle]['euler data']['velocity'] = particle_interaction.Collision_with_Walls('top',material,particle,material.particle_lineation[particle],barriers,material.particle_mechanics[particle]['euler data'],material.cell_size)
			elif topright.x > barriers['window outline']['right']['outline']:
				### the particle braches the right of the window...
				### collision with the wall...
				material.particle_mechanics[particle]['euler data']['velocity'] = particle_interaction.Collision_with_Walls('right',material,particle,material.particle_lineation[particle],barriers,material.particle_mechanics[particle]['euler data'],material.cell_size)
			elif bottomright.y > barriers['window outline']['bottom']['outline']:
				### the particle pasts the bottom of the window...
				### collision with the wall...
				material.particle_mechanics[particle]['euler data']['velocity'] = particle_interaction.Collision_with_Walls('bottom',material,particle,material.particle_lineation[particle],barriers,material.particle_mechanics[particle]['euler data'],material.cell_size)
			elif bottomleft.x < barriers['window outline']['left']['outline']:
				### the particle pasts the bottom of the window...
				### collision with the wall...
				material.particle_mechanics[particle]['euler data']['velocity'] = particle_interaction.Collision_with_Walls('left',material,particle,material.particle_lineation[particle],barriers,material.particle_mechanics[particle]['euler data'],material.cell_size)
				
		#print(particle,' ',material.particle_mechanics[particle]['euler data']['velocity'],' material.particle_mechanics[euler data][velocity] check')
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
				#relation_between_particles = Vector2(snapped(material.particle_lineation[particle].get_center().x - material.particle_lineation[other_particle].get_center().x,.01),snapped(material.particle_lineation[particle].get_center().y - material.particle_lineation[other_particle].get_center().y,.01))
				#relation_between_particles = Vector2(snapped(material.particle_lineation[particle].origin.x - material.particle_lineation[other_particle].origin.x,.01),snapped(material.particle_lineation[particle].origin.y - material.particle_lineation[other_particle].origin.y,.01))
				
				#print(' ')
				#print(relation_between_particles,' relation_between_particles')
				material.particle_lineation[particle].origin.x = material.particle_lineation[particle].origin.x - (relation_between_particles.x/2.0)
				material.particle_lineation[particle].origin.y = material.particle_lineation[particle].origin.y - (relation_between_particles.y/2.0)
				
				material.particle_lineation[other_particle].origin.x = material.particle_lineation[other_particle].origin.x - (relation_between_particles.x/2.0)
				material.particle_lineation[other_particle].origin.y = material.particle_lineation[other_particle].origin.y - (relation_between_particles.y/2.0)
				
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
		
		material.particle_mechanics[particle]['mass'] = material.mass_in_pieces
		#material.particle_mechanics[particle]['velocity'] = Vector2(0.0,0.0)
		material.particle_mechanics[particle]['velocity'] = material.initial_velocity
		#material.particle_mechanics[particle].stress = [1,0,0,1]
		#
		material.particle_mechanics[particle]['B'] = [0,0,0,0]
		### loop counter...
		identify_number = wrapi(identify_number+1,0,len(material.particle_mechanics.keys())+1)


#func Grid_to_Particle(time_passed:float,material:Object,the_grid:Dictionary):
func Grid_to_Particle(time_passed:float,material:Object):
	#particle_interaction = get_tree().get_root().get_node("Simulation/Particle Interaction")
	matrix_math = get_tree().get_root().get_node("Simulation/Matrix Math")
	#constitutive_models = get_tree().get_root().get_node("Simulation/Constitutive Models")
	#print('Grid_to_Particle check')
	
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
		for other_particle in material.particle_mechanics[particle]['within_range']:
			
			stored_relation = []
			number = 0
			
			for node in material.particle_mechanics[other_particle]['eulerian']:
				
				relation_between_particles = Vector2(snapped(node.x - material.particle_lineation[particle].origin.x,.01),snapped(node.y - material.particle_lineation[particle].origin.y,.01))
				#print(relation_between_particles,' relation_between_particles check')
				stored_relation.append(relation_between_particles)
				
				Weight_Interpolation(basis,basis_function_version,relation_between_particles,material.appearance.x)
				
				### MPM,MLS MPM paraticle velocity update
				# particle Velocity = sum of ( grid velocity * weight_interpolation )
				# B = sum of (weight interpolation * grid velocity * relation between particles transposed )
				
				material.particle_mechanics[particle]['velocity'] = material.particle_mechanics[particle]['velocity'] + (material.particle_mechanics[other_particle]['euler data']['velocity'] * weight_interpolation['weight'])
				#print(material.particle_mechanics[particle]['velocity'],' at each node')
			
			#print(material.particle_mechanics[particle]['velocity'],' during G2P')
			while number < len(material.particle_mechanics[other_particle]['eulerian']):
				# MPM B, MLS MPM
				#material.particle_mechanics[particle]['B'] = material.particle_mechanics[particle]['velocity']
				
				var velocity_relation = matrix_math.Multiply_Vector2_by_Vector2_to_Matrix(material.particle_mechanics[particle]['velocity'],false,stored_relation[number],true)
				#print(material.particle_mechanics[particle]['B'],' B check')
				#print(velocity_relation,' velocity_relation check')
				material.particle_mechanics[particle]['B'] = matrix_math.Add_Matrix(material.particle_mechanics[particle]['B'],velocity_relation)
				
				number = number + 1
				
				
		#MLS MPM 
		#particle.C = D^-1 * B
		# B =  sum of the_grid( particle.velocity * flipped_kernel_distance^T )
		# D^-1 = (4/pow(cell_size,2)) * particle.I^-1
		### Compute C
		#print(basis_coefficient,' c_cell_inverse_coefficient check')
		#print(material.appearance.x,' material.appearance.x')
		c_cell_coefficient =  snapped((basis_coefficient / pow(material.appearance.x,2.0)),.01)
		c_inverse_I_coefficient = matrix_math.Inverse_Matrix(material.particle_mechanics[particle]['I'])
		c_cell_inverse_coefficient = matrix_math.Multiply_Matrix_by_Scalar(c_inverse_I_coefficient,c_cell_coefficient,true)
		
		#print(material.particle_mechanics[particle]['B'],' B check')
		#print(c_cell_inverse_coefficient,' c_cell_inverse_coefficient check')
		#print(material.particle_mechanics[particle]['C'],' C check')
		material.particle_mechanics[particle]['C'] = matrix_math.Multiply_Matrix(c_cell_inverse_coefficient,material.particle_mechanics[particle]['B'])
		#print(material.particle_mechanics[particle]['C'],' C check')
		
		###particle position update...
		print(material.particle_mechanics[particle]['velocity'],' velocity as particle check')
		#material.particle_lineation[particle].position.x = snapped((material.particle_lineation[particle].position.x + snapped((time_passed * material.particle_mechanics[particle]['velocity'].x),.01)  ),.01)
		#material.particle_lineation[particle].position.y = snapped((material.particle_lineation[particle].position.y + snapped((time_passed * material.particle_mechanics[particle]['velocity'].y),.01)  ),.01)
		
		material.particle_lineation[particle].origin.x = snapped((material.particle_lineation[particle].origin.x + snapped((time_passed * material.particle_mechanics[particle]['velocity'].x),.01)  ),.01)
		material.particle_lineation[particle].origin.y = snapped((material.particle_lineation[particle].origin.y + snapped((time_passed * material.particle_mechanics[particle]['velocity'].y),.01)  ),.01)
		
		### MLS-deformation update...
		#particle.F = (particle.I + time_passed * particle.C ) * particle.F
		
		f_term_1 = matrix_math.Multiply_Matrix_by_Scalar(material.particle_mechanics[particle]['C'],time_passed,true)
		#print(f_term_1,' time_passed * particle.C')
		f_term_2 =  matrix_math.Add_Matrix(material.particle_mechanics[particle]['I'],f_term_1)
		#print(f_term_2,' particle.I + f term 1')
		#print(material.particle_mechanics[particle]['F'],' F before final multi check')
		material.particle_mechanics[particle]['F'] = matrix_math.Multiply_Matrix(f_term_2,material.particle_mechanics[particle]['F'])
		#print(material.particle_mechanics[particle]['F'],' F Update check')
		### Updating Plasticity...
		material.particle_mechanics[particle]['F'] = constitutive_models.Update_Plasticity(material.type_of_substance,material.yield_surface,material.particle_mechanics[particle]['F'])
		#print(material.particle_mechanics[particle]['F'],' F Update check')
		
		
		identify_number = wrapi(identify_number+1,0,len(material.particle_mechanics.keys())+1)
	#"""
