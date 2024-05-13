extends Node


### Wall Mechanics...
var wall : Dictionary = {}
### gravity of the simulation...
var gravity : Vector2 = Vector2(0.0,9.8)
### establishing grid...
var inverted_grid : float
var basis : String 
var basis_coefficient :float
var basis_function_version : int
var grid : Dictionary = {}
var unusable_grid : Array = []
var grid_spacing : int
var _by_greatest_common_factor : int
### indexing/searching mechanics/components...
var identify_number = 0
var particle = 'null'
var identify_particle = 0
### Usage caused by models...]
var updating_plasticity : bool
### for weight interpolation...
var weight_interpolation : float
### Functions of Ponts to Grid...
#var relation_between_particle_and_grid_node : float
#var relation_between_grid_node_and_particle : float
### resetting of the particle...
var updated_polygon : PackedVector2Array
### maintence of kernel mechanics...
var negative_one_point_five_by_marker#(-1.5*grid_spacing)/grid_spacing
var negative_point_five_by_marker#(-.50*grid_spacing)/grid_spacing
var positive_point_five_by_marker#(.50*grid_spacing)/grid_spacing
var positive_one_point_five_by_marker#(1.5*grid_spacing)/grid_spacing
var one_point_one_two_fiive	#(9.0/8.0)



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


func Kernel(splines:String,version:int,relation:Vector2):#,cell_size:float):
	#reset parameters
	var piecewise : Vector2 = Vector2(0,0)
	
	if splines == "quadratic":
		#if version == :
			#if relation.x >= positive_point_five_by_marker and relation.x < 
		if version == 2:
			### thru Analysis and Reduction of Quadrature Errors in the Material Point Method. steffan .W.pdf-pg(8.equation 16)
			if relation.x >= negative_one_point_five_by_marker and relation.x < negative_point_five_by_marker:
				piecewise.x = ((1.0/ (2.0 * (pow(grid_spacing,2)))) * (pow(relation.x,2))) + ((3.0/(2.0*grid_spacing))*relation.x) + one_point_one_two_fiive
			elif relation.x >= negative_point_five_by_marker and relation.x <= positive_point_five_by_marker:
				piecewise.x =  (1.0/pow(grid_spacing,2) * pow(relation.x,2)) + (3.0/4.0)
			elif relation.x > positive_point_five_by_marker and relation.x <= positive_point_five_by_marker:
				piecewise.x = ((1.0/ (2.0 * (pow(grid_spacing,2)))) * (pow(relation.x,2))) - ((3.0/(2.0*grid_spacing))*relation.x) + one_point_one_two_fiive
			else:
				piecewise.x = 0
				
			if relation.y >= negative_one_point_five_by_marker and relation.y < negative_point_five_by_marker:
				piecewise.y = ((1.0/ (2.0 * (pow(grid_spacing,2)))) * (pow(relation.y,2))) + ((3.0/(2.0*grid_spacing))*relation.y) + one_point_one_two_fiive
			#elif relation.y == (0*grid_spacing)/grid_spacing:# and relation.y <= (.50*grid_spacing):
			elif relation.y >= negative_point_five_by_marker and relation.y <= positive_point_five_by_marker:
				piecewise.y = (1.0/pow(grid_spacing,2) * pow(relation.y,2)) + (3.0/4.0)
			elif relation.y > positive_point_five_by_marker and relation.y <= positive_point_five_by_marker:
				piecewise.y = ((1.0/ (2.0 * (pow(grid_spacing,2)))) * (pow(relation.y,2))) - ((3.0/(2.0*grid_spacing))*relation.y) + one_point_one_two_fiive
			else:
				### location :: 
				piecewise.y = 0
	
	return piecewise
	



func _ready():
	### grid fuction setup...
	var gathered_into_chunks : bool
	var exact_cells : bool 
	#var _by_greatest_common_factor
	
	### Wall Mechanics...defaults
	wall['top'] = {'barrier': 0,'coefficient of restitution' : 0.0,'coefficient of static friction': 0.0,'coefficient of kinetic friction': 0.0,'mass':0.0,'velocity':Vector2(0.0,0.0)}
	wall['right'] = {'barrier':get_tree().get_root().get_window().size.x,'coefficient of restitution' : 0.0,'coefficient of static friction': 0.0,'coefficient of kinetic friction': 0.0,'mass':0.0,'velocity':Vector2(0.0,0.0)}
	wall['bottom'] = {'barrier':get_tree().get_root().get_window().size.y,'coefficient of restitution' : 0.0,'coefficient of static friction': 0.0,'coefficient of kinetic friction': 0.0,'mass':0.0,'velocity':Vector2(0.0,0.0)}
	wall['left'] = {'barrier': 0,'coefficient of restitution' : 0.0,'coefficient of static friction': 0.0,'coefficient of kinetic friction': 0.0,'mass':0.0,'velocity':Vector2(0.0,0.0)}
	
	for outline in wall:
		if outline == 'top':
			#""" null-void
			wall[outline]['coefficient of restitution'] = 1.00
			wall[outline]['coefficient of static friction'] = .43
			wall[outline]['coefficient of kinetic friction'] = .26
			wall[outline]['mass'] = 100.00
			wall[outline]['velocity'] = Vector2(0.0,0.0)
			#"""
			""" hyperelastic  - natural rubber
			wal[outline]['coefficient of restitution'] = 0.9
			wall[outline]['coefficient of static friction'] = 0.8
			wall[outline]['coefficient of kinetic friction'] = 0.6
			wall[outline]['mass'] = 1.00
			wall[outline]['velocity'] = Vector2(0.0,0.0)
			#"""
		if outline == 'right':
			#""" null-void
			wall[outline]['coefficient of restitution'] = 1.00
			wall[outline]['coefficient of static friction'] = .43
			wall[outline]['coefficient of kinetic friction'] = .26
			wall[outline]['mass'] = 100.00
			wall[outline]['velocity'] = Vector2(0.0,0.0)
			#"""
			""" hyperelastic  - natural rubber
			wal[outline]['coefficient of restitution'] = 0.9
			wall[outline]['coefficient of static friction'] = 0.8
			wall[outline]['coefficient of kinetic friction'] = 0.6
			wall[outline]['mass'] = 1.00
			wall[outline]['velocity'] = Vector2(0.0,0.0)
			#"""
		if outline == 'bottom':
			#""" null-void
			wall[outline]['coefficient of restitution'] = 1.00
			wall[outline]['coefficient of static friction'] = .43
			wall[outline]['coefficient of kinetic friction'] = .26
			wall[outline]['mass'] = 100.00
			wall[outline]['velocity'] = Vector2(0.0,0.0)
			#"""
			""" hyperelastic - natural rubber
			wal[outline]['coefficient of restitution'] = 0.9
			wall[outline]['coefficient of static friction'] = 0.8
			wall[outline]['coefficient of kinetic friction'] = 0.6
			wall[outline]['mass'] = 1.00
			wall[outline]['velocity'] = Vector2(0.0,0.0)
			#"""
		if outline == 'left':
			#""" null-void
			wall[outline]['coefficient of restitution'] = 1.00
			wall[outline]['coefficient of static friction'] = .43
			wall[outline]['coefficient of kinetic friction'] = .26
			wall[outline]['mass'] = 100.00
			wall[outline]['velocity'] = Vector2(0.0,0.0)
			#"""
			""" hyperelastic  - natural rubber
			wal[outline]['coefficient of restitution'] = 0.9
			wall[outline]['coefficient of static friction'] = 0.8
			wall[outline]['coefficient of kinetic friction'] = 0.6
			wall[outline]['mass'] = 1.00
			wall[outline]['velocity'] = Vector2(0.0,0.0)
			#"""

	### Grid Mechanics...
	basis = 'quadratic' #"cubic" #: [piecewise linear,quadratic B-spline, cubic B-spline]
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
		grid_spacing = 1
	
	inverted_grid = 1.0/grid_spacing
	
	###mechanic of kernel...
	negative_one_point_five_by_marker = (-1.5*grid_spacing)/grid_spacing
	negative_point_five_by_marker = (-.50*grid_spacing)/grid_spacing
	positive_point_five_by_marker = (.50*grid_spacing)/grid_spacing
	positive_point_five_by_marker = (1.5*grid_spacing)/grid_spacing
	one_point_one_two_fiive = (9.0/8.0)










func Grid_Reset(material:Object):#,matrix_math:Object):
	### reseting the grid data...
	#
	identify_number = 0
	grid.clear()
	var brand_new : Dictionary = {}
	grid = brand_new.duplicate()
	
	var surrounding_nodes
	
	
	# inf , nan check...
	#"""
	for location in range(0,(len(material.F))):
		if is_inf(material.F[location]) == true or is_nan(material.F[location]) == true:
			if location == 0 or location == 3:
				material.F[location] = 1.0
			if location == 1 or location == 2:
				material.F[location] = 0.0
	#"""
		
	
	#"""
	### cycle thru each particle...
	while true:
		if identify_number >= len(material.mechanics.keys()):
			break
		
		var particle = material.mechanics.keys()[identify_number]
	#"""
	#for particle in material.mechanics.keys():
		### relation between the particle wth the grid...
		#var euler_base = material.entity_container[identify_number].position
		
		var euler_base = material.entity_container[particle].position
		
		var whole_x = round(euler_base.x/grid_spacing)*grid_spacing
		var whole_y = round(euler_base.y/grid_spacing)*grid_spacing
		
		### mechanics to find the particle in relation to the whole...
		surrounding_nodes = [Vector2( (whole_x) -grid_spacing,(whole_y)-grid_spacing),
			Vector2( (whole_x),(whole_y)-grid_spacing),
			Vector2( (whole_x)+grid_spacing,(whole_y)-grid_spacing),
			Vector2((whole_x)-grid_spacing,(whole_y)),
			Vector2((whole_x),(whole_y)),
			Vector2((whole_x)+grid_spacing,(whole_y)),
			Vector2((whole_x)-grid_spacing,(whole_y)+grid_spacing),
			Vector2((whole_x),(whole_y)+grid_spacing),
			Vector2((whole_x)+grid_spacing,(whole_y)+grid_spacing),
			]
			
		material.mechanics[particle]['eulerian'] = PackedVector2Array(surrounding_nodes.duplicate(true))
		
		var grid_index
		#"""
		for node in surrounding_nodes:
			grid_index = material.mechanics[particle]['eulerian'].find(node)
			#print(node,' node surrounding check')
			### Grid Mechanics
			if grid.has(node) == false:
				
				### this grid node is to be affected by the particle.
				grid[node] = {'mass': 0.0,'velocity':Vector2(0.0,0.0),'momentum':Vector2(0.0,0.0),'forces':Vector2(0.0,0.0),'angular momentum':Vector2(0,0),'normal force':0.0}.duplicate(true)
			else:
				### already recoqnized...
				pass
			
		
			material.mechanics[particle]['particle to grid'][grid_index] = (euler_base - node) / grid_spacing
			material.mechanics[particle]['grid to particle'][grid_index] = (node - euler_base) / grid_spacing
			
			var relation_inverted = (inverted_grid*material.mechanics[particle]['particle to grid'][grid_index])
			
			var kerneled =  Kernel(basis,basis_function_version,relation_inverted)
			material.mechanics[particle]['weight interpolation'][grid_index] = kerneled.x * kerneled.y
			
		
		#material.particle_mechanics[particle]['F'] = material.particle_mechanics[particle]['I'].duplicate(true)
		#material.mechanics[particle]['J'] = get_tree().get_root().get_node('Simulation').matrix_math.Find_Determinant(material.mechanics[particle]['F'])
		
		#material.mechanics[particle]['J'] = get_tree().get_root().get_node('Simulation').matrix_math.Find_Determinant(material.F)
		material.J = get_tree().get_root().get_node('Simulation').matrix_math.Find_Determinant(material.F)
		material.mechanics[particle]['volume'] = material.volume * material.J
		
		
		### cycle thru each polygon...
		identify_number = wrapi(identify_number+1,0,len(material.mechanics.keys())+1)
		


func Particles_to_Grid(time_passed:float,material:Object):
	#particle_interaction = get_tree().get_root().get_node("Simulation/Particle Interaction")
	# Particles to Grid:
	#for particle in material.particle_mechanics.keys():
	identify_number = 0
	
	#"""
	while true:
		if identify_number >= len(material.mechanics.keys()):
			break
		
		particle = material.mechanics.keys()[identify_number]
	#"""
	#for particle in material.mechanics.keys():
		#### the sum of aspects that is resetted...
		var energy_density_is_known = false
		var energy_denisty_coeficient = [0,0,0,0]
		var affline_force_contribution = [0,0,0,0]
			
		#constitutive models...
		if material.physical_state == 'solid':
			if material.constitutive_model == 'hyperelastic':
				energy_denisty_coeficient =get_tree().get_root().get_node('Simulation/Constitutive Models').Neo_Hookean(particle,material)
				energy_density_is_known = true
			if material.constitutive_model == 'fixed_corated':
				energy_denisty_coeficient = get_tree().get_root().get_node('Simulation/Constitutive Models').Fixed_Corotated(particle,material)
				energy_density_is_known = true
				updating_plasticity = true
			if material.constitutive_model == 'drucker_prager_elasticity':
				energy_denisty_coeficient = get_tree().get_root().get_node('Simulation/Constitutive Models').Drucker_Prager_Elasticity(particle,material)
				energy_density_is_known = true
				updating_plasticity = true
		elif material.physical_state == 'liquid':
			if material.constitutive_model == 'water':
				energy_denisty_coeficient = get_tree().get_root().get_node('Simulation/Constitutive Models').Model_of_Water(particle,material)
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
			#var timed_volume_I = get_tree().get_root().get_node('Simulation/Matrix Math').Multiply_Matrix_by_Scalar(material.mechanics[particle]['I'],timed_volume_coefficient,true)
			var timed_volume_I = get_tree().get_root().get_node('Simulation/Matrix Math').Multiply_Matrix_by_Scalar(material.I,timed_volume_coefficient,true)
			
			#var timed_volume_I = [material.mechanics[particle]['I'][0]*timed_volume_coefficient,
			#material.mechanics[particle]['I'][1]*timed_volume_coefficient,
			#material.mechanics[particle]['I'][2]*timed_volume_coefficient,
			#material.mechanics[particle]['I'][3]*timed_volume_coefficient]
			
			var energy_density_I = get_tree().get_root().get_node('Simulation/Matrix Math').Multiply_Matrix(energy_denisty_coeficient,timed_volume_I)
				
			#var  F_transposed = get_tree().get_root().get_node('Simulation/Matrix Math').Transposed_Matrix(material.mechanics[particle]['F'])
			var  F_transposed = get_tree().get_root().get_node('Simulation/Matrix Math').Transposed_Matrix(material.F)
			
			var energy_density_I_transposed = get_tree().get_root().get_node('Simulation/Matrix Math').Multiply_Matrix(energy_density_I,F_transposed)
				
			var q_mass_C = get_tree().get_root().get_node('Simulation/Matrix Math').Multiply_Matrix_by_Scalar(material.mechanics[particle]['C'],material.mass,true)
			#var q_mass_C = [material.mechanics[particle]['C'][0]*material.mass,
			#	material.mechanics[particle]['C'][1]*material.mass,
			#	material.mechanics[particle]['C'][2]*material.mass,
			#	material.mechanics[particle]['C'][3]*material.mass
			#	]
			
			affline_force_contribution = get_tree().get_root().get_node('Simulation/Matrix Math').Add_Matrix(energy_density_I_transposed,q_mass_C)
			#material.mechanics[particle]['affline force contribution'] = get_tree().get_root().get_node('Simulation/Matrix Math').Add_Matrix(energy_density_I_transposed,q_mass_C)
			
			
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
					
			#var inverse_I_coefficient = get_tree().get_root().get_node('Simulation/Matrix Math').Inverse_Matrix(material.mechanics[particle]['I'])
			var inverse_I_coefficient = get_tree().get_root().get_node('Simulation/Matrix Math').Inverse_Matrix(material.I)
			var weighted_stressed_I = get_tree().get_root().get_node('Simulation/Matrix Math').Multiply_Matrix(inverse_I_coefficient,q_stressed_coefficient)
			var negative_weighted_stressed_I = get_tree().get_root().get_node('Simulation/Matrix Math').Divide_Matrix_by_Scalar(weighted_stressed_I,-1.0,true)
					
			var q_mass_C = get_tree().get_root().get_node('Simulation/Matrix Math').Multiply_Matrix_by_Scalar(material.mechanics[particle]['C'],material.mass,true)
					
			affline_force_contribution = get_tree().get_root().get_node('Simulation/Matrix Math').Add_Matrix(negative_weighted_stressed_I,q_mass_C)
			#material.mechanics[particle]['affline force contribution'] = get_tree().get_root().get_node('Simulation/Matrix Math').Add_Matrix(negative_weighted_stressed_I,q_mass_C)
			
			
		#"""
		#print(len(grid),' grid')
		
		for grid_axis in material.mechanics[particle]['eulerian']:
			var relation_index = material.mechanics[particle]['eulerian'].find(grid_axis)
			weight_interpolation = material.mechanics[particle]['weight interpolation'][relation_index] 
			var relation_between_grid_node_and_particle = material.mechanics[particle]['grid to particle'][relation_index]
			
			var weighted_affine_force = get_tree().get_root().get_node('Simulation/Matrix Math').Multiply_Matrix_by_Scalar(affline_force_contribution,weight_interpolation,true)
			
			#var weighted_affine_force = [affline_force_contribution[0]*weight_interpolation,
			#affline_force_contribution[1] *weight_interpolation,
			#affline_force_contribution[2] *weight_interpolation,
			#affline_force_contribution[3] *weight_interpolation,]
			
			var affine = get_tree().get_root().get_node('Simulation/Matrix Math').Multiply_Matrix_by_Vector2_to_Vector2(weighted_affine_force,relation_between_grid_node_and_particle)
					
			### MPM lumped mass
			
			### MLS MPM lumped mass
			#grid[grid_axis]['mass'] = snappedf((grid[grid_axis]['mass'] + (material.mass * weight_interpolation)),.001)
			
			#print((material.mass * weight_interpolation))
			
			grid[grid_axis]['mass'] = grid[grid_axis]['mass'] + (material.mechanics[particle]['mass'] * weight_interpolation)
			grid[grid_axis]['momentum'] = grid[grid_axis]['momentum'] + (material.mechanics[particle]['mass']*affine)
		
		#"""
		### cycle thru each polygon...
		identify_number = wrapi(identify_number+1,0,len(material.mechanics.keys())+1)
		



func Grid_Update(time_passed:float,_material:Object):
	#Grid Update:
	#print()
	identify_number = 0
	#"""
	### cycle thru each polygon...
	while true:
		if identify_number >= len(grid.keys()):
			break
		
		var node_point = grid.keys()[identify_number]
	#"""
	#for node_point in grid.keys():
		#print(node_point,' node ')
		if grid[node_point]['mass'] > 0.0:
			if node_point.y >= wall['bottom']['barrier'] or node_point.y <= wall['top']['barrier']:
				gravity = Vector2(0,0)
			if node_point.x >= wall['right']['barrier'] or node_point.x <= wall['left']['barrier']:
				gravity = Vector2(0,0)
			grid[node_point]['momentum'] = grid[node_point]['momentum'] + (time_passed * (grid[node_point]['mass'] * gravity))
			grid[node_point]['velocity'] = grid[node_point]['momentum'] / grid[node_point]['mass']
			
			#print(grid[node_point]['velocity'],' node velo')
		### loop counter...
		identify_number = wrapi(identify_number+1,0,len(grid.keys())+1)
	
func Adapt_Constraint():
	### constraint of the soft body applied to the grid nodes...
	
	pass



func Collision_with_Wall(material:Object):

	### Collision Detection...
	### how the particles interacts....
	identify_number = 0
	
	#"""
	### cycle thru each polygon...
	while true:
		if identify_number >= len(material.mechanics.keys()):
			break
		
		particle = material.mechanics.keys()[identify_number]
	#"""
	#for particle in material.mechanics.keys():
		### checks iof the particle comes in contract with the wall...
		if material.entity_container[identify_number].position.y <= wall['top']['barrier']:
			
			material.mechanics[particle]['velocity'] = get_tree().get_root().get_node('Simulation/Interaction').Collision_with_Walls('top',material,particle,wall['top'])
			### Apply friction along the the top...
			if gravity.y < 0:
				material.mechanics[particle]['velocity'] = get_tree().get_root().get_node('Simulation/Interaction').Apply_Friction(wall,'top','left','right',material.mechanics[particle]['mass'],material.mechanics[particle]['velocity'],gravity,material.coefficient_of_static_friction,material.coefficient_of_kinetic_friction)
			
		if material.entity_container[identify_number].position.x >= wall['right']['barrier']:
			material.mechanics[particle]['velocity'] = get_tree().get_root().get_node('Simulation/Interaction').Collision_with_Walls('right',material,particle,wall['right'])
			
			### Apply friction along the the right...
			if gravity.x > 0:
				material.mechanics[particle]['velocity'] = get_tree().get_root().get_node('Simulation/Interaction').Apply_Friction(wall,'right','top','bottom',material.mechanics[particle]['mass'],material.mechanics[particle]['velocity'],gravity,material.coefficient_of_static_friction,material.coefficient_of_kinetic_friction)
			
		if material.entity_container[identify_number].position.y >= wall['bottom']['barrier']:
			
			material.mechanics[particle]['velocity'] = get_tree().get_root().get_node('Simulation/Interaction').Collision_with_Walls('bottom',material,particle,wall['bottom'])
			### Apply friction along the the bottom...
			if gravity.y > 0:
				material.mechanics[particle]['velocity'] = get_tree().get_root().get_node('Simulation/Interaction').Apply_Friction(wall,'bottom','left','right',material.mechanics[particle]['mass'],material.mechanics[particle]['velocity'],gravity,material.coefficient_of_static_friction,material.coefficient_of_kinetic_friction)
			
		if material.entity_container[identify_number].position.x <= wall['left']['barrier']:
			
			material.mechanics[particle]['velocity'] = get_tree().get_root().get_node('Simulation/Interaction').Collision_with_Walls('left',material,particle, wall['left'])
			### Apply friction along the the left...
			if gravity.x < 0:
				material.mechanics[particle]['velocity'] = get_tree().get_root().get_node('Simulation/Interaction').Apply_Friction(wall,'left','top','bottom',material.mechanics[particle]['mass'],material.mechanics[particle]['velocity'],gravity,material.coefficient_of_static_friction,material.coefficient_of_kinetic_friction)
		
		### loop counter...
		identify_number = wrapi(identify_number+1,0,len(material.mechanics.keys())+1)
		


func Particle_Reset(material:Object):
	### particle is reset...
	identify_number = 0
	particle = 'null'
	
	#"""
	while true:
		if identify_number >= len(material.mechanics.keys()):
			break
		
		particle = material.mechanics.keys()[identify_number]
	#"""
	#for particle in material.mechanics.keys():
		material.mechanics[particle]['mass'] = material.mass / len(material.entity_container)
		#material.mechanics[particle]['velocity'] = Vector2(0.0,0.0)
		material.mechanics[particle]['B'] = [0,0,0,0]
		#material.particle_mechanics[particle]['C'] = [0,0,0,0]
		
		### loop counter...
		identify_number = wrapi(identify_number+1,0,len(material.mechanics.keys())+1)
	
	


func Grid_to_Particle(time_passed:float,material:Object):
	identify_number = 0
	particle = 'null'
	
	#"""
	### cycle thru each polygon...
	while true:
		if identify_number >= len(material.mechanics.keys()):
			break
		
		particle = material.mechanics.keys()[identify_number]
	#"""
	#for particle in material.mechanics.keys():
		#"""
		for node_location in material.mechanics[particle]['eulerian']:
				
			var relation_index = material.mechanics[particle]['eulerian'].find(node_location)
			weight_interpolation = material.mechanics[particle]['weight interpolation'][relation_index] 
			var relation_between_grid_node_and_particle = material.mechanics[particle]['grid to particle'][relation_index]
			
			### MPM,MLS MPM paraticle velocity update
			
			# particle Velocity = sum of ( grid velocity * weight_interpolation )
			#B = sum of (weight interpolation * grid velocity * relation between particles transposed )
			material.mechanics[particle]['velocity'] = material.mechanics[particle]['velocity'] + ( grid[node_location]['velocity'] * weight_interpolation)
			
			var velocity_relation = get_tree().get_root().get_node('Simulation/Matrix Math').Multiply_Vector2_by_Vector2_to_Matrix(grid[node_location]['velocity'],false,relation_between_grid_node_and_particle,true)
				
			material.mechanics[particle]['B'] = get_tree().get_root().get_node('Simulation/Matrix Math').Add_Matrix(material.mechanics[particle]['B'],velocity_relation)
		
		#print()
		#print(material.mechanics[particle]['velocity'],' ',particle)
		
		#"""
		#MLS MPM 
		#particle.C = B * D^-1
		# B =  sum of the_grid( (particle.velocity*weight_interpolation) * flipped_kernel_distance^T )
		# D^-1 = (4/pow(grid_spacing,2)) * particle.I^-1
		### Compute C
		var b_cell_coefficient = get_tree().get_root().get_node('Simulation/Matrix Math').Multiply_Matrix_by_Scalar(material.mechanics[particle]['B'],4.0,true)
		var b_coefficient = get_tree().get_root().get_node('Simulation/Matrix Math').Divide_Matrix_by_Scalar(b_cell_coefficient,pow(grid_spacing,2.0),true)
		#var c_inverse_I_coefficient = get_tree().get_root().get_node('Simulation/Matrix Math').Inverse_Matrix(material.mechanics[particle]['I'])
		var c_inverse_I_coefficient = get_tree().get_root().get_node('Simulation/Matrix Math').Inverse_Matrix(material.I)
		#print()
		#print(pow(grid_spacing,2.0),' pow(grid_spacing,2.0)')
		#print(b_cell_coefficient,' b_cell_coefficient')
		#print(b_coefficient,' b_coefficient2')
		#print(c_inverse_I_coefficient,' c_inverse_I_coefficient')
		material.mechanics[particle]['C'] = get_tree().get_root().get_node('Simulation/Matrix Math').Multiply_Matrix(b_coefficient,c_inverse_I_coefficient)
		
		###particle position update...
		var new_particle_position = Vector2(0,0)
		
		#print( material.mechanics[particle]['velocity'],' ')
		new_particle_position.x = snapped((material.entity_container[particle].position.x + snapped((time_passed * material.mechanics[particle]['velocity'].x),.01)  ),.01)
		new_particle_position.y = snapped((material.entity_container[particle].position.y + snapped((time_passed * material.mechanics[particle]['velocity'].y),.01)  ),.01)
		#
		material.entity_container[particle] = Rect2(new_particle_position,Vector2(1,1))
		#print(material.entity_container[particle],' material.entity_container[particle]')
		### MLS-deformation update...
		#particle.F = (particle.I + time_passed * particle.C ) * particle.F
		var f_term_1 = get_tree().get_root().get_node('Simulation/Matrix Math').Multiply_Matrix_by_Scalar(material.mechanics[particle]['C'],time_passed,true)
		#var f_term_2 = get_tree().get_root().get_node('Simulation/Matrix Math').Add_Matrix(material.mechanics[particle]['I'],f_term_1)
		var f_term_2 = get_tree().get_root().get_node('Simulation/Matrix Math').Add_Matrix(material.I,f_term_1)
		#material.mechanics[particle]['F'] = get_tree().get_root().get_node('Simulation/Matrix Math').Multiply_Matrix(f_term_2,material.mechanics[particle]['F'])
		material.F = get_tree().get_root().get_node('Simulation/Matrix Math').Multiply_Matrix(f_term_2,material.F)
		
		if updating_plasticity == true:
			### Updating Plasticity...
			get_tree().get_root().get_node('Simulation/Constitutive Models').Update_Plasticity(material,particle,material.type_of_substance)
			updating_plasticity = false
		
		### cycle thru each polygon...
		identify_number = wrapi(identify_number+1,0,len(material.mechanics.keys())+1)
	

