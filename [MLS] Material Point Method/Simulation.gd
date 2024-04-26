extends Node

''' Base Anchor of the Simulation '''
### the simulation happens here...

#parts of the simulation
# core mathematics- [ MLS Material Point Method]
# supplementary math - Matricies
# particles sent thru the simualtion - in the form of Meshes...
var different_substances : int = 1
var forge : Object
var laboratory : Object
var matrix_math : Object
var constitutive_models : Object
var colliding : Object
var material_method : Object
var substance : Object
### ...
var identify_substance = 0
### establish grid / grid mechanics...
var grid : Dictionary
var complete_grid : Dictionary = {}
var basis
var basis_coefficient
var basis_function_version
### basic math grid mechanics...
var negative_one_point_five_by_marker#(-1.5*grid_spacing)/grid_spacing
var negative_point_five_by_marker#(-.50*grid_spacing)/grid_spacing
var positive_point_five_by_marker#(.50*grid_spacing)/grid_spacing
var positive_one_point_five_by_marker#(1.5*grid_spacing)/grid_spacing
var one_point_one_two_fiive	#(9.0/8.0)
var grid_spacing
var inverted_grid


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





func Quality_of_Grid():
	### determine the size/cells of the grid...
	### Grid Mechanics...
	var basis = 'quadratic' #"cubic" #: [piecewise linear,quadratic B-spline, cubic B-spline]
	var basis_coefficient = 4.0 #:: [ uadratic B-spline - 4,  cubic B-spline - 3
	var basis_function_version = 2
	
	### determine how to cut the gridd into even pieces...
	var gathered_into_chunks = true
	var exact_cells = false
	
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
			var _by_greatest_common_factor =  By_Greatest_Common_Factor(float(get_tree().get_root().get_window().size.x),float(get_tree().get_root().get_window().size.y))
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



func _ready():
	if FileAccess.file_exists("res://MeshForge.tscn"):
		forge = load("res://MeshForge.tscn").instantiate()
	else:
		### file does not exists...
		print(' MeshForge.tscn does not exists...')
	if FileAccess.file_exists("res://PhysicalSubstance.tscn"):
		#laboratory = load("res://PhysicalSubstance.tscn").instantiate()
		pass
	else:
		### file does not exists...
		print(' PhysicalSubstance.tscn does not exists...')
	if FileAccess.file_exists("res://MatrixMath.tscn"):
		#laboratory = load("res://PhysicalSubstance.tscn").instantiate()
		matrix_math = load("res://MatrixMath.tscn").instantiate()
	else:
		### file does not exists...
		print(' MatrixMath.tscn does not exists...')
	if FileAccess.file_exists("res://ConstitutiveModels.tscn"):
		#laboratory = load("res://PhysicalSubstance.tscn").instantiate()
		constitutive_models = load("res://ConstitutiveModels.tscn").instantiate()
	else:
		### file does not exists...
		print(' ConstitutiveModels.tscn does not exists...')
	if FileAccess.file_exists("res://FundamentalColliding.tscn"):
		#laboratory = load("res://PhysicalSubstance.tscn").instantiate()
		colliding = load("res://FundamentalColliding.tscn").instantiate()
	else:
		### file does not exists...
		print(' FundamentalColliding.tscn does not exists...')
	if FileAccess.file_exists("res://MaterialMethod.tscn"):
		#laboratory = load("res://PhysicalSubstance.tscn").instantiate()
		material_method = load("res://MaterialMethod.tscn").instantiate()
	else:
		### file does not exists...
		print(' MaterialMethod.tscn does not exists...')
	
	
	### determine how many different substances...
	different_substances = 1
	
	while len(get_tree().get_root().get_node('Simulation/container').get_children()) < different_substances:
		
		add_child(forge)
		
		forge.free()
		forge = load("res://MeshForge.tscn").instantiate()
	
	#add_child(forge)
	add_child(matrix_math)
	add_child(constitutive_models)
	### used for process...
	#substance = get_tree().get_root().get_node('Simulation/container').get_children()[identify_substance]
	
	add_child(material_method)
	add_child(colliding)
	
	#forge.queue_free()
	
	### determine the size/cells of the grid...
	#Quality_of_Grid()
	
	#print()
	#print(Vector2(0,0) - Vector2((grid_spacing*5),(grid_spacing*5)),' grid start' )
	#print(Vector2(get_tree().get_root().get_window().size.x,get_tree().get_root().get_window().size.y) + Vector2((grid_spacing*5),(grid_spacing*5)),' grid end' )

	#print()
	#print(len(complete_grid),' complete_grid')
	
	set_process(true)
	set_physics_process(false)





func _notification(event):
	###
	#print(event,' event')
	if event == 1003:
		### Window Size Change...
		#print('window change')
		pass

#"""
func _process(delta):
	
	while true:
		if identify_substance >= len(get_tree().get_root().get_node('Simulation/container').get_children()):
			identify_substance = 0
			break
		
		substance = get_tree().get_root().get_node('Simulation/container').get_children()[identify_substance]
		
		$"Material Method".Grid_Reset(substance)
		$"Material Method".Particles_to_Grid(delta,substance)
		$"Material Method".Grid_Update(delta,substance)
		$"Material Method".Collision_with_Wall(substance)
		$"Material Method".Particle_Reset(substance)
		$"Material Method".Grid_to_Particle(delta,substance)
		
		substance.queue_redraw()
		
		### cycle thru each polygon...
		identify_substance = wrapi(identify_substance+1,0,len(get_tree().get_root().get_node('Simulation/container').get_children())+1)
	
#"""
