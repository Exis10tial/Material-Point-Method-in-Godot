extends Node



### gravity of the simulation...
var gravity : Vector2 = Vector2(0,9.8)
### establishing grid...
var basis : String 
var basis_coefficient :float
var basis_function_version : int
var grid : Dictionary = {}
var grid_spacing : int
var _by_greatest_common_factor : int
### grid reseet mechanics/components...
var identify_number = 0
var particle = 'null'
### Usage caused by models...]
var updating_plasticity : bool
### for weight interpolation...
var weight_interpolation : Dictionary
### Functions of Ponts to Grid...
var relation_between_particle_and_grid_node : Vector2
var relation_between_grid_node_and_particle : Vector2
### resetting of the particle...
var updated_polygon : PackedVector2Array



func By_Greatest_Common_Factor(multinomial_x,multinomial_y):
	### the cell size will be the greatest common factor between the lenght and width of the substance...
	var x_factors : Array = []
	var y_factors : Array = []
	var size
	
	#find the x factors...
	for number in range(1,multinomial_x+1):
		if int(multinomial_x) % number == 0:
			# that number is a factor...
			x_factors.append(number)
	
	#find the y factors...
	for number in range(1,multinomial_y+1):
		if int(multinomial_y) % number == 0:
			# that number is a factor...
			y_factors.append(number)
	
	# remove itself and 1 from list..
	x_factors.pop_back()
	x_factors.pop_front()
	
	y_factors.pop_back()
	y_factors.pop_front()
	
	#finding/selecting the largest factor between the number..
	if multinomial_y > multinomial_x:
		for number in range(len(y_factors)-1,-1,-1):
			if y_factors[number] in x_factors:
				size = y_factors[number]
				break
			else:
				size = false
	else:
		for number in range(len(x_factors)-1,-1,-1):
			if x_factors[number] in y_factors:
				size = x_factors[number]
				break
			else:
				size = false
	return size


func Weight_Interpolation(splines:String,version:int,updated_relation:Vector2):#,cell_size:float):
	#reset parameters
	var building_weights = Vector2(0.0,0.0)
	### Weight Interpolation...
	
	if splines == "quadratic":
		if version == 1:
			### The Material Point Method for Simulating Continuum Materials pg 33 equation 123 
			### unchanged
			
			#if 1.5*grid_spacing <= absf(updated_relation.x):
			if 1.5*grid_spacing <= absf(updated_relation.x):
				### particles is not in range of the particle to influence.
				building_weights.x =  0.0
			
			#elif .50*grid_spacing <= absf(updated_relation.x) and absf(updated_relation.x) < 1.5*grid_spacing:
			elif .50*grid_spacing <= absf(updated_relation.x) and absf(updated_relation.x) < 1.5*grid_spacing:
				building_weights.x = snappedf((.50 *  pow((1.50 - absf(updated_relation.x)),2.0)),.0001)
				#building_weights.x = .50 *  pow( (1.50 - absf(updated_relation.x)),2.0)
			
			#elif 0 <= absf(updated_relation.x) and absf(updated_relation.x) < .50*grid_spacing:
			elif 0*grid_spacing <= absf(updated_relation.x) and absf(updated_relation.x) < .50*grid_spacing:
				building_weights.x = snappedf((.75 - snappedf(pow(absf(updated_relation.x),2.0),.01)),.0001)
				#building_weights.x = .75 - pow(absf(updated_relation.x),2.0)
				
			if 1.5*grid_spacing <= abs(updated_relation.y):
				building_weights.y = 0.0
			
			elif .50*grid_spacing <= absf(updated_relation.y) and absf(updated_relation.y) < 1.5*grid_spacing:
				building_weights.y = snappedf((.50 * pow((1.50 - absf(updated_relation.y) ),2.0)),.0001)
				#building_weights.y = .50 * pow( (1.50 - absf(updated_relation.y) ),2.0)
			
			elif 0*grid_spacing <= absf(updated_relation.y) and absf(updated_relation.y) < .50*grid_spacing:
				building_weights.y = snappedf((.75 - snappedf(pow(absf(updated_relation.y),2.0),.01)),.0001)
				#building_weights.y = .75 - pow(absf(updated_relation.y),2.0)
		
		elif version == 2:
			### thru Analysis and Reduction of Quadrature Errors in the Material Point Method. steffan .W.pdf-pg(8.equation 16)
			# slight change
			if updated_relation.x > (-0.50*grid_spacing) and updated_relation.x < (0*grid_spacing):
				### location :: inbetween node center edge to halfway inbetween nodes a
				building_weights.x = ((1.0/ (2.0 * (pow(grid_spacing,2)))) * (pow(updated_relation.x,2))) + ((3.0/(2.0*grid_spacing))*updated_relation.x) + (9.0/8.0)
			elif updated_relation.x == (0*grid_spacing):# and updated_relation.x <= (.50*grid_spacing):
				### location :: area around node center... 
				building_weights.x =   -((1 / (pow(grid_spacing,2))) * pow(updated_relation.x,2)) + (3.0/4.0)
			elif updated_relation.x > (0*grid_spacing) and updated_relation.x < (.50*grid_spacing):
				### location :: inbetween node center edge to halfway inbetween nodes b
				building_weights.x = ((1.0/ (2.0 * (pow(grid_spacing,2)))) * (pow(updated_relation.x,2))) - ((3.0/(2.0*grid_spacing))*updated_relation.x) + (9.0/8.0)
			else:
				### location :: space/area between 4 nodes.
				building_weights.x = 0
			
			if updated_relation.y > (-0.50*grid_spacing) and updated_relation.y < (0*grid_spacing):
				### location :: inbetween node center edge to halfway inbetween nodes a
				building_weights.y = ((1.0/ (2.0 * (pow(grid_spacing,2)))) * (pow(updated_relation.y,2))) + ((3.0/(2.0*grid_spacing))*updated_relation.y) + (9.0/8.0)
			elif updated_relation.y == (0*grid_spacing):# and updated_relation.y <= (.50*grid_spacing):
				### location :: area around node center... 
				building_weights.y = -((1 / (pow(grid_spacing,2))) * pow(updated_relation.x,2)) + (3.0/4.0)
			elif updated_relation.y > (0*grid_spacing) and updated_relation.y < (.50*grid_spacing):
				### location :: inbetween node center edge to halfway inbetween nodes b
				building_weights.y = ((1.0/ (2.0 * (pow(grid_spacing,2)))) * (pow(updated_relation.y,2))) - ((3.0/(2.0*grid_spacing))*updated_relation.y) + (9.0/8.0)
			else:
				### location :: space/area between 4 nodes.
				building_weights.y = 0
			
		weight_interpolation['weight'] = snappedf((building_weights.x * building_weights.y),.001)
		
	elif splines == "cubic":
		pass
	
	return weight_interpolation
	



func _ready():
	### grid fuction setup...
	var gathered_into_chunks : bool
	var exact_cells : bool 
	#var _by_greatest_common_factor
	
	
	##basis = "cubic"
	##basis_coefficient = 3.0
	
	basis = "quadratic" #: [piecewise linear,quadratic B-spline, cubic B-spline]
	basis_coefficient = 4.0 #:: [ uadratic B-spline - 4,  cubic B-spline - 3
	basis_function_version = 2
	
	
	### determine how to cut the gridd into even pieces...
	gathered_into_chunks = true
	exact_cells = false
	
	## maybe deifferent shapes/outline for cells...
	
	###determine the size of the grid cells...
	# the size of the cells is parsed into evenly cut squares/cubes = (Length == Width)
	# the siize of window outline is the base...
	# gathered_into_chunks the cells is cut in evenly equal pieces...
	# exact_cells the cellls is 1...
	
	if gathered_into_chunks == true and exact_cells == false:
		#the cell is cut in evenly equal pieces...
		## defualt/lowest number of cells is 4...if possible
		# check to see  a power of 2...
		# check to see if it has a square root...
		# check to see if divisible by 2....
		# then the greatest common factor between the particle length and width.
		
		
		if get_tree().get_root().get_window().size.x == get_tree().get_root().get_window().size.y:
			
			pass
			
		elif get_tree().get_root().get_window().size.x != get_tree().get_root().get_window().size.y:
			
			### using the greastest common factor...
			_by_greatest_common_factor =  By_Greatest_Common_Factor(float(get_tree().get_root().get_window().size.x),float(get_tree().get_root().get_window().size.y))
			if typeof(_by_greatest_common_factor) != TYPE_BOOL:
				grid_spacing = _by_greatest_common_factor 
			else:
				### doesn't have _by_greatest_common_factor
				grid_spacing = 1
			#grid_spacing = (get_tree().get_root().get_window().size.x/_by_greatest_common_factor) * (get_tree().get_root().get_window().size.y/_by_greatest_common_factor)
			
	elif gathered_into_chunks == false and exact_cells == true:
		#the cell size is 1....
		pass




func Grid_Reset(material:Object):#,matrix_math:Object):
	### reseting the grid data...
	#
	identify_number = 0
	particle = 'null'
	grid = {}
	var surrounding_nodes
	while true:
		if identify_number >= len(material.mechanics.keys()):
			break
		
		particle = material.mechanics.keys()[identify_number]
		
		### relation between the particle wth the grid...
		if grid_spacing > 1:
			###...
			#"""
			# find the nearest node (non-float number)
			var euler_base = material.effigy.get_polygon()[particle] 
			#var test_euler = material.effigy.get_polygon()[particle] / grid_spacing
			var results_euler = Vector2(roundi(material.effigy.get_polygon()[particle].x / grid_spacing),roundi(material.effigy.get_polygon()[particle].y / grid_spacing)   )
			var end_euler = results_euler * grid_spacing
			
			surrounding_nodes = [Vector2( (round(euler_base.x/grid_spacing)*grid_spacing) -grid_spacing,(round(euler_base.y/grid_spacing)*grid_spacing)-grid_spacing),
				Vector2( (round(euler_base.x/grid_spacing)*grid_spacing),(round(euler_base.y/grid_spacing)*grid_spacing)-grid_spacing),
				Vector2( (round(euler_base.x/grid_spacing)*grid_spacing)+grid_spacing,(round(euler_base.y/grid_spacing)*grid_spacing)-grid_spacing),
				Vector2((round(euler_base.x/grid_spacing)*grid_spacing)-grid_spacing,(round(euler_base.y/grid_spacing)*grid_spacing)),
				Vector2((round(euler_base.x/grid_spacing)*grid_spacing),(round(euler_base.y/grid_spacing)*grid_spacing)),
				Vector2((round(euler_base.x/grid_spacing)*grid_spacing)+grid_spacing,(round(euler_base.y/grid_spacing)*grid_spacing)),
				Vector2((round(euler_base.x/grid_spacing)*grid_spacing)-grid_spacing,(round(euler_base.y/grid_spacing)*grid_spacing)+grid_spacing),
				Vector2((round(euler_base.x/grid_spacing)*grid_spacing),(round(euler_base.y/grid_spacing)*grid_spacing)+grid_spacing),
				Vector2((round(euler_base.x/grid_spacing)*grid_spacing)+grid_spacing,(round(euler_base.y/grid_spacing)*grid_spacing)+grid_spacing),
			]
		else:
			pass
		
		material.mechanics[particle]['eulerian'] = surrounding_nodes.duplicate(true)
		
		for node in surrounding_nodes:
			if grid.has(node) == false:
				### this grid node is to be affected by the particle.
				grid[node] = {'mass': 0.0,'velocity':Vector2(0.0,0.0),'momentum':Vector2(0.0,0.0),'forces':Vector2(0.0,0.0),'angular momentum':Vector2(0,0),'normal force':0.0}.duplicate(true)
			else:
				### already recoqnized...
				pass
		
		# inf , nan check...
		for location in range(0,(len(material.mechanics[particle]['F']))):
			if is_inf(material.mechanics[particle]['F'][location]) == true or is_nan(material.mechanics[particle]['F'][location]) == true:
				if location == 0 or location == 2:
					material.mechanics[particle]['F'][location] = 1.0
				if location == 1 or location == 3:
					material.mechanics[particle]['F'][location] = 0.0
		
		#material.particle_mechanics[particle]['F'] = material.particle_mechanics[particle]['I'].duplicate(true)
		material.mechanics[particle]['J'] = get_tree().get_root().get_node('Simulation').matrix_math.Find_Determinant(material.mechanics[particle]['F'])
		material.mechanics[particle]['volume'] = material.volume * material.mechanics[particle]['J']
		
		identify_number = wrapi(identify_number+1,0,len(material.mechanics.keys())+1)
		


func Particles_to_Grid(time_passed:float,material:Object):
	#particle_interaction = get_tree().get_root().get_node("Simulation/Particle Interaction")
	# Particles to Grid:
	#for particle in material.particle_mechanics.keys():
	identify_number = 0
	particle = 'null'

	while true:
		if identify_number >= len(material.mechanics.keys()):
			break
		
		particle = material.mechanics.keys()[identify_number]
		
		#### the sum of aspects that is resetted...
		var energy_density_is_known = false
		var energy_denisty_coeficient = [0,0,0,0]
		#var affine_momentum_fuse = Vector2(0.0,0.0)
		
		var affline_force_contribution = [0,0,0,0]
		
		#constitutive models...
		if material.physical_state == 'solid':
			if material.constitutive_model == 'hyperelastic':
				energy_denisty_coeficient = get_tree().get_root().get_node('Simulation').constitutive_models.Neo_Hookean(particle,material)
				energy_density_is_known = true
			if material.constitutive_model == 'fixed_corated':
				energy_denisty_coeficient = get_tree().get_root().get_node('Simulation').constitutive_models.Fixed_Corotated(particle,material)
				energy_density_is_known = true
				updating_plasticity = true
			if material.constitutive_model == 'drucker_prager_elasticity':
				energy_denisty_coeficient = get_tree().get_root().get_node('Simulation').constitutive_models.Drucker_Prager_Elasticity(particle,material)
				energy_density_is_known = true
				updating_plasticity = true
		elif material.physical_state == 'liquid':
			if material.constitutive_model == 'water':
				energy_denisty_coeficient = get_tree().get_root().get_node('Simulation').constitutive_models.Model_of_Water(particle,material)
				energy_density_is_known = true
					
		elif material.physical_state == 'gas':
			energy_density_is_known = true
			pass
		
		if energy_density_is_known == true:
			### MLS MPM:  fuses the scattering of the affine momentum and particle force contribution
			# weight_interpolation * (Q)* (relation_between_grid_node_and_particle)
			#Q = time * particle intial volume * D^-1 * differentiate(energy density)/differentiate(F) * F * F_transposed + partical mass * C
			# D^-1 = (4/pow(grid_spacing,2)) * particle.I^-1
			
			var timed_volume_coefficient = time_passed * material.volume * 4.0 / pow(grid_spacing,2.0)
			var timed_volume_I = get_tree().get_root().get_node('Simulation/Matrix Math').Multiply_Matrix_by_Scalar(material.mechanics[particle]['I'],timed_volume_coefficient,true)
			var energy_density_I = get_tree().get_root().get_node('Simulation/Matrix Math').Multiply_Matrix(energy_denisty_coeficient,timed_volume_I)
			
			var  F_transposed = get_tree().get_root().get_node('Simulation/Matrix Math').Transposed_Matrix(material.mechanics[particle]['F'])
			
			var energy_density_I_transposed = get_tree().get_root().get_node('Simulation/Matrix Math').Multiply_Matrix(energy_density_I,F_transposed)
			
			var q_mass_C = get_tree().get_root().get_node('Simulation/Matrix Math').Multiply_Matrix_by_Scalar(material.mechanics[particle]['C'],material.mass,true)
				
			affline_force_contribution = get_tree().get_root().get_node('Simulation/Matrix Math').Add_Matrix(energy_density_I_transposed,q_mass_C)
			
		else:
			### MLS MPM:  fuses the scattering of the affine momentum and particle force contribution
			# weight_interpolation * (Q)* (relation_between_grid_node_and_particle)
			#  if energy density is unknown.
			#Q = - sum of (volume of material occupied at p * stess * (differentiate)weight_interpolation + mass * C)
			#(differentiate)weight_interpolation = weight_interpolation * D^-1 * (relation_between_particle_grid_node)
			# D^-1 = (4/pow(cell_size,2)) * particle.I^-1
				
			#var q_timed_volumed =  time_passed * material.particle_mechanics[particle]['volume']
			var q_volumed_stress =  get_tree().get_root().get_node('Simulation/Matrix Math').Multiply_Matrix_by_Scalar(material.mechanics[particle]['stress'],material.mechanics[particle]['volume'],true)
			var q_volumed_stress_coefficient = get_tree().get_root().get_node('Simulation/Matrix Math').Multiply_Matrix_by_Scalar(q_volumed_stress,4.0,true)
			var q_stressed_coefficient = get_tree().get_root().get_node('Simulation/Matrix Math').Divide_Matrix_by_Scalar(q_volumed_stress_coefficient,pow(grid_spacing,2),true)
				
			var inverse_I_coefficient = get_tree().get_root().get_node('Simulation/Matrix Math').Inverse_Matrix(material.mechanics[particle]['I'])
			var weighted_stressed_I = get_tree().get_root().get_node('Simulation/Matrix Math').Multiply_Matrix(inverse_I_coefficient,q_stressed_coefficient)
			var negative_weighted_stressed_I = get_tree().get_root().get_node('Simulation/Matrix Math').Divide_Matrix_by_Scalar(weighted_stressed_I,-1.0,true)
				
			var q_mass_C = get_tree().get_root().get_node('Simulation/Matrix Math').Multiply_Matrix_by_Scalar(material.mechanics[particle]['C'],material.mass,true)
				
			affline_force_contribution = get_tree().get_root().get_node('Simulation/Matrix Math').Add_Matrix(negative_weighted_stressed_I,q_mass_C)
		
		for grid_axis in material.mechanics[particle]['eulerian']:
			
			relation_between_particle_and_grid_node = Vector2((snapped(material.effigy.get_polygon()[particle].x - grid_axis.x,.01)),snapped(material.effigy.get_polygon()[particle].y - grid_axis.y,.01)) / grid_spacing
			relation_between_grid_node_and_particle = Vector2(snapped(grid_axis.x - material.effigy.get_polygon()[particle].x,.01),snapped(grid_axis.y - material.effigy.get_polygon()[particle].y,.01)) / grid_spacing
			
			Weight_Interpolation(basis,basis_function_version,relation_between_particle_and_grid_node)#,grid_spacing)
			
			var weighted_affine_force = get_tree().get_root().get_node('Simulation/Matrix Math').Multiply_Matrix_by_Scalar(affline_force_contribution,weight_interpolation['weight'],true)
			var affine = get_tree().get_root().get_node('Simulation/Matrix Math').Multiply_Matrix_by_Vector2_to_Vector2(weighted_affine_force,relation_between_grid_node_and_particle)
				
			### MPM lumped mass
			#grid[node]['mass'] = grid[node]['mass'] + (material.mass * weight_interpolation['weight'])
			
			### MLS MPM lumped mass
			grid[grid_axis]['mass'] = snappedf((grid[grid_axis]['mass'] + (material.mass * weight_interpolation['weight'])),.001)
			
			grid[grid_axis]['momentum'] = grid[grid_axis]['momentum'] + (grid[grid_axis]['mass'] * (material.mechanics[particle]['initial velocity'] + affine))
			
		
		identify_number = wrapi(identify_number+1,0,len(material.mechanics.keys())+1)
	


func Grid_Update(time_passed:float,_material:Object):
	#Grid Update:
	for node_point in grid.keys():
		if grid[node_point]['mass'] > 0.0:
			grid[node_point]['momentum'] = grid[node_point]['momentum'] + (time_passed * (grid[node_point]['mass'] * gravity))
			
			grid[node_point]['velocity'] = grid[node_point]['momentum'] / grid[node_point]['mass']
			

func Collision_with_Wall(material:Object):
	#"""
	### Collision Detection...
	### how the particles interacts....
	
	identify_number = 0
	var node = null
	
	while true:
		if identify_number >= len(grid.keys()):
			break
		
		node = grid.keys()[identify_number]
	
		if grid[node]['mass'] > 0.0:
			
			for identified_particle in material.mechanics.keys():
				#identify the particle/s
				if node == material.mechanics[identified_particle]['eulerian'][4]:
					
					pass
					
		identify_number = wrapi(identify_number+1,0,len(grid.keys())+1)
	






func Particle_Reset(material:Object):
	### particle is reset...
	identify_number = 0
	particle = 'null'
	
	while true:
		if identify_number >= len(material.mechanics.keys()):
			break
		
		particle = material.mechanics.keys()[identify_number]
		
		material.mechanics[particle]['mass'] = material.mass / len(material.mechanics.keys())
		
		material.mechanics[particle]['velocity'] = Vector2(0.0,0.0)
		material.mechanics[particle]['angular momentum'] = [0,0,0,0]
		material.mechanics[particle]['B'] = [0,0,0,0]
		#material.particle_mechanics[particle]['C'] = [0,0,0,0]
		### loop counter...
		identify_number = wrapi(identify_number+1,0,len(material.mechanics.keys())+1)
	
	updated_polygon = PackedVector2Array([])


func Grid_to_Particle(time_passed:float,material:Object):
	
	identify_number = 0
	particle = 'null'
	while true:
		if identify_number >= len(material.mechanics.keys()):
			break
		
		particle = material.mechanics.keys()[identify_number]
		
		for node_location in material.mechanics[particle]['eulerian']:
			
			relation_between_particle_and_grid_node = Vector2(snapped(material.effigy.get_polygon()[particle].x - node_location.x,.01),snapped(material.effigy.get_polygon()[particle].y - node_location.y,.01)) / grid_spacing
			relation_between_grid_node_and_particle = Vector2(snapped(node_location.x - material.effigy.get_polygon()[particle].x,.01),snapped(node_location.y - material.effigy.get_polygon()[particle].y,.01)) / grid_spacing
			
			Weight_Interpolation(basis,basis_function_version,relation_between_particle_and_grid_node)#,grid_spacing)
			
			### MPM,MLS MPM paraticle velocity update
			# particle Velocity = sum of ( grid velocity * weight_interpolation )
			#B = sum of (weight interpolation * grid velocity * relation between particles transposed )
			material.mechanics[particle]['velocity'] = material.mechanics[particle]['velocity'] + ( grid[node_location]['velocity'] * weight_interpolation['weight'])
			
			var velocity_relation = get_tree().get_root().get_node('Simulation/Matrix Math').Multiply_Vector2_by_Vector2_to_Matrix(grid[node_location]['velocity'],false,relation_between_grid_node_and_particle,true)
			
			material.mechanics[particle]['B'] = get_tree().get_root().get_node('Simulation/Matrix Math').Add_Matrix(material.mechanics[particle]['B'],velocity_relation)
			
		#MLS MPM 
		#particle.C = B * D^-1
		# B =  sum of the_grid( (particle.velocity*weight_interpolation) * flipped_kernel_distance^T )
		# D^-1 = (4/pow(grid_spacing,2)) * particle.I^-1
		### Compute C
		
		var b_cell_coefficient = get_tree().get_root().get_node('Simulation/Matrix Math').Multiply_Matrix_by_Scalar(material.mechanics[particle]['B'],4.0,true)
		var b_coefficient = get_tree().get_root().get_node('Simulation/Matrix Math').Divide_Matrix_by_Scalar(b_cell_coefficient,pow(grid_spacing,2.0),true)
		var c_inverse_I_coefficient = get_tree().get_root().get_node('Simulation/Matrix Math').Inverse_Matrix(material.mechanics[particle]['I'])
		
		material.mechanics[particle]['C'] = get_tree().get_root().get_node('Simulation/Matrix Math').Multiply_Matrix(b_coefficient,c_inverse_I_coefficient)
		
		
		###particle position update...
		var new_particle_position = Vector2(0,0)
		
		new_particle_position.x = snapped((material.effigy.get_polygon()[particle].x + snapped((time_passed * material.mechanics[particle]['velocity'].x),.01)  ),.01)
		new_particle_position.y = snapped((material.effigy.get_polygon()[particle].y + snapped((time_passed * material.mechanics[particle]['velocity'].y),.01)  ),.01)
		### set polygon...
		updated_polygon.append(new_particle_position)
		
		### MLS-deformation update...
		#particle.F = (particle.I + time_passed * particle.C ) * particle.F
		
		var f_term_1 = get_tree().get_root().get_node('Simulation/Matrix Math').Multiply_Matrix_by_Scalar(material.mechanics[particle]['C'],time_passed,true)
		var f_term_2 = get_tree().get_root().get_node('Simulation/Matrix Math').Add_Matrix(material.mechanics[particle]['I'],f_term_1)
		material.mechanics[particle]['F'] = get_tree().get_root().get_node('Simulation/Matrix Math').Multiply_Matrix(f_term_2,material.mechanics[particle]['F'])
		
		if updating_plasticity == true:
			### Updating Plasticity...
			get_tree().get_root().get_node('Simulation').constitutive_models.Update_Plasticity(material,particle,material.type_of_substance)
			updating_plasticity = false
		
		identify_number = wrapi(identify_number+1,0,len(material.mechanics.keys())+1)
	
	### the particle is updated...
	material.effigy.set_polygon([])
	material.effigy.set_polygon(updated_polygon)
	

