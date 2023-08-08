extends Node


var particle_interaction
var matrix_math
var constitutive_models


### Eulerian/Cartesian Grid...
#var Grid_Nodes : Dictionary = {}
var grid_cells : int
var grid_spacing : float
var grid_nodes : Dictionary = {}
var station : Vector2
#var identify_grid_reach : Array = []
#var gravity : Vector2
var gravity : Vector2 = Vector2(0.00,9.80)# * 10
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
var affine : Vector2
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
	
	grid_cells = (ProjectSettings.get_setting('display/window/size/viewport_width')*2.0) * ((ProjectSettings.get_setting('display/window/size/viewport_height')*2.0))
	grid_spacing = grid_cells #+ (grid_cells/2.0)
	
	
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
		
		#""" null-void
		#barriers['window outline'][wall]['coefficient of restitution'] = 1.0
		barriers['window outline'][wall]['coefficient of restitution'] = 0.50
		#barriers['window outline'][wall]['coefficient of restitution'] = 0.0
		#barriers['window outline'][wall]['coefficient of static friction'] = 1.0
		barriers['window outline'][wall]['coefficient of static friction'] = 0.25
		#barriers['window outline'][wall]['coefficient of static friction'] = 0.0
		#barriers['window outline'][wall]['coefficient of kinetic friction'] = 1.0
		barriers['window outline'][wall]['coefficient of kinetic friction'] = .13
		#barriers['window outline'][wall]['coefficient of kinetic friction'] = 0.0
		barriers['window outline'][wall]['mass'] = 10.00
		#barriers['window outline'][wall]['mass'] = 1.00
		#"""
		""" hyperlaster - natural rubber
		barriers['window outline'][wall]['coefficient of restitution'] = 0.9
		#barriers['window outline'][wall]['coefficient of restitution'] = 0.50
		#barriers['window outline'][wall]['coefficient of restitution'] = 0.0
		#barriers['window outline'][wall]['coefficient of static friction'] = 1.0
		barriers['window outline'][wall]['coefficient of static friction'] = 0.8
		#barriers['window outline'][wall]['coefficient of static friction'] = 0.0
		#barriers['window outline'][wall]['coefficient of kinetic friction'] = 1.0
		barriers['window outline'][wall]['coefficient of kinetic friction'] = 0.6
		#barriers['window outline'][wall]['coefficient of kinetic friction'] = 0.0
		barriers['window outline'][wall]['mass'] = 1000.00
		#barriers['window outline'][wall]['mass'] = 1.00
		#"""
		
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
			
			#if 1.5*grid_spacing <= absf(updated_relation.x):
			if 1.5 <= absf(updated_relation.x):
				### particles is not in range of the particle to influence.
				#print('3')
				building_weights.x =  0.0
			
			#elif .50*grid_spacing <= absf(updated_relation.x) and absf(updated_relation.x) < 1.5*grid_spacing:
			elif .50 <= absf(updated_relation.x) and absf(updated_relation.x) < 1.5:
				#print('2')
				building_weights.x = snappedf((.50 *  pow((1.50 - absf(updated_relation.x)),2.0)),.0001)
				#building_weights.x = .50 *  pow( (1.50 - absf(updated_relation.x)),2.0)
			
			#elif 0 <= absf(updated_relation.x) and absf(updated_relation.x) < .50*grid_spacing:
			elif 0 <= absf(updated_relation.x) and absf(updated_relation.x) < .50:
				#print('1')
				building_weights.x = snappedf((.75 - snappedf(pow(absf(updated_relation.x),2.0),.01)),.0001)
				#print(building_weights,' a building_weights')
				#building_weights.x = .75 - pow(absf(updated_relation.x),2.0)
				
			#if 1.5*grid_spacing <= abs(updated_relation.y):
			if 1.5 <= abs(updated_relation.y):
				#print('8')
				building_weights.y = 0.0
			
			#elif .50*grid_spacing <= absf(updated_relation.y) and absf(updated_relation.y) < 1.5*grid_spacing:
			elif .50 <= absf(updated_relation.y) and absf(updated_relation.y) < 1.5:
				#print('7')
				building_weights.y = snappedf((.50 * pow((1.50 - absf(updated_relation.y) ),2.0)),.0001)
				#building_weights.y = .50 * pow( (1.50 - absf(updated_relation.y) ),2.0)
			
			#elif 0 <= absf(updated_relation.y) and absf(updated_relation.y) < .50*grid_spacing:
			elif 0 <= absf(updated_relation.y) and absf(updated_relation.y) < .50:
				#print('6')
				building_weights.y = snappedf((.75 - snappedf(pow(absf(updated_relation.y),2.0),.01)),.0001)
				#print(building_weights,' a building_weights')
				#building_weights.y = .75 - pow(absf(updated_relation.y),2.0)
		elif version == 2:
			### thru Analysis and Reduction of Quadrature Errors in the Material Point Method. steffan .W.pdf-pg(8.equation 16)
			pass
			
		###...
		#print(building_weights.x,' building weights x')
		#print(building_weights.y,' building weights y')
		
		weight_interpolation['weight'] = snappedf((building_weights.x * building_weights.y),.0001)
		#weight_interpolation['weight'] = (building_weights.x * building_weights.y)
		#print(building_weights.y,' building weights y')
		

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
		#print(euler_base,' euler_base')
		
		"""
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
		"""
		material.particle_mechanics[particle]['eulerian'] = [Vector2(euler_base.x,euler_base.y),
		Vector2(euler_base.x+grid_spacing,euler_base.y),
		Vector2(euler_base.x+(grid_spacing*2),euler_base.y),
		Vector2(euler_base.x,euler_base.y-grid_spacing),
		Vector2(euler_base.x+grid_spacing,euler_base.y-grid_spacing),
		Vector2(euler_base.x+(grid_spacing*2),euler_base.y-grid_spacing),
		Vector2(euler_base.x,euler_base.y-(grid_spacing*(grid_spacing*2))),
		Vector2(euler_base.x+grid_spacing,euler_base.y-(grid_spacing*2)),
		Vector2(euler_base.x+(grid_spacing*2),euler_base.y-(grid_spacing*2)),
		]
		#"""
		#"""
		material.particle_mechanics[particle]['eulerian'] = [Vector2(euler_base.x-grid_spacing,euler_base.y-grid_spacing),
		Vector2(euler_base.x,euler_base.y-grid_spacing),
		Vector2(euler_base.x+grid_spacing,euler_base.y-grid_spacing),
		Vector2(euler_base.x-grid_spacing,euler_base.y),
		Vector2(euler_base.x,euler_base.y),
		Vector2(euler_base.x+grid_spacing,euler_base.y),
		Vector2(euler_base.x-grid_spacing,euler_base.y+grid_spacing),
		Vector2(euler_base.x,euler_base.y+grid_spacing),
		Vector2(euler_base.x+grid_spacing,euler_base.y+grid_spacing),
		]
		#"""
		#print(material.particle_lineation[particle].origin,' origin check')
		#print(euler_base,' euler_base check')
		#print()
		for surrounding_node in material.particle_mechanics[particle]['eulerian']:
			if grid_nodes.has(surrounding_node) == false:
				### this grid node is to be affected by the particle.
				grid_nodes[surrounding_node] = {'mass': 0.0,'velocity':Vector2(0.0,0.0),'momentum':Vector2(0.0,0.0),'forces':Vector2(0.0,0.0),'angular momentum':Vector2(0,0),'normal force':0.0}.duplicate(true)
				#print(surrounding_node,' surrounding_node check')
				#print(material.particle_lineation[particle].origin-surrounding_node,' relation check')
			else:
				### already recoqnized...
				pass
		
		material.particle_mechanics[particle]['w'] = Vector2(0,0)
		
		#print(material.particle_mechanics[particle]['mass'],' mass check')
		#print(gravity,' gravity check')
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
		
		var affline_force_contribution = [0,0,0,0]
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
		#print(affline_force_contribution,' affline_force_contribution check')
		
		if energy_density_is_known == true:
			### MLS MPM:  fuses the scattering of the affine momentum and particle force contribution
			# weight_interpolation * (Q)* (relation_between_grid_node_and_particle)
			#Q = time * particle intial volume * D^-1 * differentiate(energy density)/differentiate(F) * F * F_transposed + partical mass * C
			# D^-1 = (4/pow(grid_spacing,2)) * particle.I^-1
			
			var timed_volume_coefficient = time_passed * material.volume * 4.0 / pow(grid_spacing,2.0)
			var timed_volume_I = matrix_math.Multiply_Matrix_by_Scalar(material.particle_mechanics[particle]['I'],timed_volume_coefficient,true)
			var energy_density_I = matrix_math.Multiply_Matrix(energy_denisty_coeficient,timed_volume_I)
			
			var  F_transposed = matrix_math.Transposed_Matrix(material.particle_mechanics[particle]['F'])
			
			var energy_density_I_transposed = matrix_math.Multiply_Matrix(energy_density_I,F_transposed)
			
			var q_mass_C = matrix_math.Multiply_Matrix_by_Scalar(material.particle_mechanics[particle]['C'],material.mass,true)
				
			affline_force_contribution = matrix_math.Add_Matrix(energy_density_I_transposed,q_mass_C)
			
		else:
			### MLS MPM:  fuses the scattering of the affine momentum and particle force contribution
			# weight_interpolation * (Q)* (relation_between_grid_node_and_particle)
			#  if energy density is unknown.
			#Q = - sum of (volume of material occupied at p * stess * (differentiate)weight_interpolation + mass * C)
			#(differentiate)weight_interpolation = weight_interpolation * D^-1 * (relation_between_particle_grid_node)
			# D^-1 = (4/pow(cell_size,2)) * particle.I^-1
				
			#var q_timed_volumed =  time_passed * material.particle_mechanics[particle]['volume']
			var q_volumed_stress =  matrix_math.Multiply_Matrix_by_Scalar( material.particle_mechanics[particle]['stress'],material.particle_mechanics[particle]['volume'],true)
			var q_volumed_stress_coefficient = matrix_math.Multiply_Matrix_by_Scalar(q_volumed_stress,4.0,true)
			var q_stressed_coefficient = matrix_math.Divide_Matrix_by_Scalar(q_volumed_stress_coefficient,pow(grid_spacing,2),true)
				
			var inverse_I_coefficient = matrix_math.Inverse_Matrix(material.particle_mechanics[particle]['I'])
			var weighted_stressed_I = matrix_math.Multiply_Matrix(inverse_I_coefficient,q_stressed_coefficient)
			var negative_weighted_stressed_I = matrix_math.Divide_Matrix_by_Scalar(weighted_stressed_I,-1.0,true)
			#weighted_stressed_I = weighted_stressed_I +
				
			var q_mass_C = matrix_math.Multiply_Matrix_by_Scalar(material.particle_mechanics[particle]['C'],material.mass,true)
				
			affline_force_contribution = matrix_math.Add_Matrix(negative_weighted_stressed_I,q_mass_C)
				

		#print(material.particle_mechanics[particle]['eulerian'],' material.particle_mechanics[particle][eulerian]')
		for node in material.particle_mechanics[particle]['eulerian']:
			#print()
			#print(different_number,' number check')
			#print(material.particle_lineation[particle].origin,' material.particle_lineation[particle].origin')
			#print(node,' node')
			
			#relation_between_particle_and_grid_node = Vector2(snapped(material.particle_lineation[particle].origin.x - (node.x+(grid_spacing/2.0)),.01),snapped(material.particle_lineation[particle].origin.y - (node.y+(grid_spacing/2.0)),.01))
			#relation_between_grid_node_and_particle = Vector2(snapped((node.x+(grid_spacing/2.0)) - material.particle_lineation[particle].origin.x,.01),snapped((node.y+(grid_spacing/2.0)) - material.particle_lineation[particle].origin.y,.01))
			
			relation_between_particle_and_grid_node = Vector2(snapped(material.particle_lineation[particle].origin.x - node.x,.01),snapped(material.particle_lineation[particle].origin.y - node.y,.01))
			relation_between_grid_node_and_particle = Vector2(snapped(node.x - material.particle_lineation[particle].origin.x,.01),snapped(node.y - material.particle_lineation[particle].origin.y,.01))
			
			#print(relation_between_particle_and_grid_node,' relation_between_particle_and_grid_node check')
			var updated_relation = (1.0/grid_spacing) * relation_between_particle_and_grid_node
			#print((1.0/grid_spacing),' (1.0/grid_spacing) check')
			#print(updated_relation,' updated_relation check')
			Weight_Interpolation(basis,basis_function_version,updated_relation)#,grid_spacing)
			#print(weight_interpolation['weight'],' weight_interpolation[weight] check')

			#affline_force_contribution = matrix_math.Add_Matrix(affline_force_contribution,q_weighted_mass_C)
				
			var weighted_affine_force = matrix_math.Multiply_Matrix_by_Scalar(affline_force_contribution,weight_interpolation['weight'],true)
			affine = matrix_math.Multiply_Matrix_by_Vector2_to_Vector2(weighted_affine_force,relation_between_grid_node_and_particle)
				
			### MPM lumped mass
			#grid_nodes[node]['mass'] = grid_nodes[node]['mass'] + (material.mass * weight_interpolation['weight'])
			
			### MLS MPM lumped mass
			#print(material.mass, ' material.mass check')
			#print(material.particle_mechanics[particle]['euler data']['mass'], ' [euler data][mass] check')
			#print(weight_interpolation['weight'], ' weight_interpolation[weight] check')
			grid_nodes[node]['mass'] = snappedf((grid_nodes[node]['mass'] + (material.mass * weight_interpolation['weight'])),.001)
			
			#print(material.particle_mechanics[particle]['initial velocity'], ' initial velocity check')
			#print(grid_nodes[node]['mass'], ' mass check')
			#print(affine, ' affine check')
			
			grid_nodes[node]['momentum'] = grid_nodes[node]['momentum'] + (grid_nodes[node]['mass'] * (material.particle_mechanics[particle]['initial velocity'] + affine))
			
			
			#grid_nodes[node]['angular momentum'] = grid_nodes[node]['angular momentum'] + (node * grid_nodes[node]['momentum'])
			material.particle_mechanics[particle]['w'] = material.particle_mechanics[particle]['w'] + (node * grid_nodes[node]['momentum'])
			
			#print(material.particle_mechanics[particle]['w'],' material.particle_mechanics[particle][w] check')
			#print(material.particle_mechanics[particle]['normal force'],' material.particle_mechanics[particle][normal force] check')
			
		identify_number = wrapi(identify_number+1,0,len(material.particle_mechanics.keys())+1)
	#print(len(grid_nodes.keys()),' number of gride nodes',grid_nodes )
	

#func Grid_Update(material:Object,the_grid:Dictionary):#,outside_forces:Vector2):
func Grid_Update(time_passed:float,material:Object):
	#print()
	#print('Grid_Update check')
	#Grid Update:
	#print()
	#print(len(grid_nodes.keys()),' number of gride nodes')
	for node in grid_nodes.keys():
	
		if grid_nodes[node]['mass'] > 0.0:
		#	print()
			#print(grid_nodes[node]['momentum'], ' [momentum] check')
			#print(grid_nodes[node]['mass'], ' [mass] check')
			#print(time_passed, ' time_passed check')
			#print(gravity, ' gravity check')
			#print(node,' ',grid_nodes[node]['momentum'], ' momentum before check')
			#print((time_passed *( grid_nodes[node]['mass'] * gravity)), ' (time_passed * grid_nodes[node][mass] * gravity) check')
			grid_nodes[node]['momentum'] = grid_nodes[node]['momentum'] + (time_passed * (grid_nodes[node]['mass'] * gravity))
			#print(grid_nodes[node]['momentum'], ' [momentum] before mass check')
			
			if grid_nodes[node]['mass'] > 1.0:
			#	print('a')
				grid_nodes[node]['velocity'] = grid_nodes[node]['momentum'] / grid_nodes[node]['mass']
			else:
			#	print('b')
			#	print(( 1.0/ grid_nodes[node]['mass']),' check')
				grid_nodes[node]['velocity'] = grid_nodes[node]['momentum'] * ( 1.0/ grid_nodes[node]['mass'])
			#print(node,' ',grid_nodes[node]['momentum'], ' momentum update check')
			#print(node,' ',grid_nodes[node]['velocity'], ' velocity update check')
	
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
					
					### psition of Rotating Particles...
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
					
					var particle_boundary = [topleft,topright,bottomright,bottomleft]
					
					for outline in particle_boundary:
						if outline.y < barriers['window outline']['top']['outline']:
							#print('top)')
							#print(outline,' outline')
							grid_nodes[node]['velocity'] = particle_interaction.Collision_with_Walls('top',material,material.particle_lineation[particle],barriers,grid_nodes[node],material.cell_size)
							
							material.particle_mechanics[particle]['normal force'] = material.particle_mechanics[particle]['mass'] * gravity
							
							break
						if outline.x > barriers['window outline']['right']['outline']:
							#print('right)')
							#print(outline,' outline')
							grid_nodes[node]['velocity'] = particle_interaction.Collision_with_Walls('right',material,material.particle_lineation[particle],barriers,grid_nodes[node],material.cell_size)
							
							material.particle_mechanics[particle]['normal force'] = material.particle_mechanics[particle]['mass'] * gravity
							
							break
						if outline.y > barriers['window outline']['bottom']['outline']:
							#print('bottom)')
							#print(outline,' outline')
							grid_nodes[node]['velocity'] =  particle_interaction.Collision_with_Walls('bottom',material,material.particle_lineation[particle],barriers,grid_nodes[node],material.cell_size)
							
							break
						if outline.x < barriers['window outline']['left']['outline']:
							#print('left)')
							#print(outline,' outline')
							grid_nodes[node]['velocity'] = particle_interaction.Collision_with_Walls('left',material,material.particle_lineation[particle],barriers,grid_nodes[node],material.cell_size)
							
							material.particle_mechanics[particle]['normal force'] = material.particle_mechanics[particle]['mass'] * gravity
							
							break
					
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
		
		material.particle_mechanics[particle]['mass'] = material.mass / len(material.particle_mechanics.keys())
		
		material.particle_mechanics[particle]['velocity'] = Vector2(0.0,0.0)
		
		material.particle_mechanics[particle]['B'] = [0,0,0,0]
		#material.particle_mechanics[particle]['C'] = [0,0,0,0]
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
		#var stored_relation = []
		var number = 0
		#print('')
		#print( 'start g2p')
		#print('')
		#print(material.particle_mechanics[particle]['velocity'],' velocity')
		#print(material.particle_mechanics[particle]['velocity'],' material.particle_mechanics[particle][velocity] check start')
		#print(material.particle_mechanics[particle]['B'],' B before check')
		for node in material.particle_mechanics[particle]['eulerian']:
			#if grid_nodes[node]['mass'] > 0.0:
			#print()
			#print(node,' node check')
			#relation_between_particle_and_grid_node = Vector2(snapped(material.particle_lineation[particle].origin.x - (node.x+(grid_spacing/2.0)),.01),snapped(material.particle_lineation[particle].origin.y - (node.y+(grid_spacing/2.0)),.01))
			#relation_between_grid_node_and_particle = Vector2(snapped((node.x+(grid_spacing/2.0)) - material.particle_lineation[particle].origin.x,.01),snapped((node.y+(grid_spacing/2.0)) - material.particle_lineation[particle].origin.y,.01))
			
			relation_between_particle_and_grid_node = Vector2(snapped(material.particle_lineation[particle].origin.x - node.x,.01),snapped(material.particle_lineation[particle].origin.y - node.y,.01))
			relation_between_grid_node_and_particle = Vector2(snapped(node.x - material.particle_lineation[particle].origin.x,.01),snapped(node.y - material.particle_lineation[particle].origin.y,.01))
			
			
			#print(relation_between_grid_node_and_particle,' relation_between_grid_node_and_particle check')
			
			var updated_relation = (1.0/grid_spacing) * relation_between_grid_node_and_particle
			#print(updated_relation,' updated_relation check')
			Weight_Interpolation(basis,basis_function_version,updated_relation)#,grid_spacing)
			#print(weight_interpolation['weight'],' weight_interpolation[weight] check')
			
			### MPM,MLS MPM paraticle velocity update
			# particle Velocity = sum of ( grid velocity * weight_interpolation )
			#B = sum of (weight interpolation * grid velocity * relation between particles transposed )
			#print(grid_nodes[node]['velocity'],' grid_nodes[node][velocity] check')
			material.particle_mechanics[particle]['velocity'] = material.particle_mechanics[particle]['velocity'] + ( grid_nodes[node]['velocity'] * weight_interpolation['weight'])
			#print(material.particle_mechanics[particle]['velocity'],' material.particle_mechanics[particle][velocity] check just added')
			#print(relation_between_particle_and_grid_node,' relation_between_particle_and_grid_node check')
			#print(material.particle_mechanics[particle]['velocity'],' material.particle_mechanics[particle][velocity] check')
			#var velocity_relation = matrix_math.Multiply_Vector2_by_Vector2_to_Matrix(material.particle_mechanics[particle]['velocity'],false,relation_between_grid_node_and_particle,true)
			var velocity_relation = matrix_math.Multiply_Vector2_by_Vector2_to_Matrix(grid_nodes[node]['velocity'],false,relation_between_grid_node_and_particle,true)
			
			#var angular_I = grid_nodes[node]['velocity'] * relation_between_grid_node_and_particle
			#print(angular_I,' angular_I check')
			#print(material.particle_mechanics[particle]['B'],' B check')
			#print(velocity_relation,' velocity_relation check')
			material.particle_mechanics[particle]['B'] = matrix_math.Add_Matrix(material.particle_mechanics[particle]['B'],velocity_relation)
			#print(material.particle_mechanics[particle]['B'],' B during check')
			
			
		material.particle_mechanics[particle]['initial velocity'] = material.particle_mechanics[particle]['velocity']
		#print()
		#print(material.particle_mechanics[particle]['velocity'],' material.particle_mechanics[particle][velocity] check')
		#print()
		#print(material.particle_mechanics[particle]['w'],' material.particle_mechanics[particle][w] check')
		#MLS MPM 
		#particle.C = B * D^-1
		# B =  sum of the_grid( (particle.velocity*weight_interpolation) * flipped_kernel_distance^T )
		# D^-1 = (4/pow(grid_spacing,2)) * particle.I^-1
		### Compute C
		
		var b_cell_coefficient = matrix_math.Multiply_Matrix_by_Scalar(material.particle_mechanics[particle]['B'],4.0,true)
		var b_coefficient = matrix_math.Divide_Matrix_by_Scalar(b_cell_coefficient,pow(grid_spacing,2.0),true)
		c_inverse_I_coefficient = matrix_math.Inverse_Matrix(material.particle_mechanics[particle]['I'])
		#print(c_inverse_I_coefficient,' c_cell_inverse_coefficient check')
		material.particle_mechanics[particle]['C'] = matrix_math.Multiply_Matrix(b_coefficient,c_inverse_I_coefficient)
		#c_cell_coefficient =  snapped((4.0 / pow(material.appearance.x,2.0)),.01)
		#c_cell_coefficient =  snapped((4.0 / pow(grid_spacing,2.0)),.01)
		#c_inverse_I_coefficient = matrix_math.Inverse_Matrix(material.particle_mechanics[particle]['I'])
		#c_cell_inverse_coefficient = matrix_math.Multiply_Matrix_by_Scalar(c_inverse_I_coefficient,c_cell_coefficient,true)
		
		#print(material.particle_mechanics[particle]['B'],' B after check')
		#print(c_cell_inverse_coefficient,' c_cell_inverse_coefficient check')
		#material.particle_mechanics[particle]['C'] = matrix_math.Multiply_Matrix(c_cell_inverse_coefficient,material.particle_mechanics[particle]['B'])
		#print(material.particle_mechanics[particle]['C'],' C check')
		
		
		###particle position update...
		#print(time_passed,' time_passed check')
		#print(material.particle_mechanics[particle]['velocity'],' velocity as particle check')
		material.particle_lineation[particle].origin.x = snapped((material.particle_lineation[particle].origin.x + snapped((time_passed * material.particle_mechanics[particle]['velocity'].x),.01)  ),.01)
		material.particle_lineation[particle].origin.y = snapped((material.particle_lineation[particle].origin.y + snapped((time_passed * material.particle_mechanics[particle]['velocity'].y),.01)  ),.01)
		
		#print(material.particle_lineation[particle].origin,' position check')
		
		
		#material.particle_lineation[particle]['w'] = material.particle_lineation[particle].origin * (material.mass*material.particle_mechanics[particle]['velocity'])
		
		
		
		
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
