extends Node

''' The particle sent thru the simulation created here 
		the particle used come from verticies a created Meshes  '''


### the resources/file for Substances
var physical_matter : Object
### used for construction of a mesh
var starting_position : Vector2
var designing_mesh : SurfaceTool
### used for creating a random name...
var letters_list : Array = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
var digits_list : Array = ['0','1','2','3','4','5','6','7','8','9']

func _ready():
	
	physical_matter = load("res://PhysicalSubstance.tscn").instantiate()
	
	### Construction of a partices...*
	#define overall outline shape of polygon.
		#shapes contain a mimumin of points...
		#line = 3 points
		#triangle = 6 points
		#square = 8 points
	# decide pylugon apperance
		#has texture?
		#color verex?
	#define physical properties of the entire polygon...
	#associate material point properties to each polygon vertex...
	
	
	### Creation of Mesh
	physical_matter.effigy = Polygon2D.new()
	
	# create a square mesh...
	var number_of_points = 99 
	var points_per_side
	var collection_of_points = []
	var coloring_points = []
	#make the number of particles a multiple of 4...
	if number_of_points % 4 != 0:
		##just the number of points to equal a multiple of 4.
		
		
		var n = 0
		
		while 4 * n <= number_of_points:
			n = n + 1
		
		number_of_points = 4 * n
		
	points_per_side = number_of_points / 4
	
	var starting_test_position = Vector2(get_tree().get_root().size.x/2.0,get_tree().get_root().size.y/2.0)
	
	var x_point = 0
	var y_point = 0
	
	while len(collection_of_points) < number_of_points:
		
		if len(collection_of_points) >= points_per_side * 1 and len(collection_of_points) < points_per_side * 2:
			### down the right: top to bottom..
			collection_of_points.append(Vector2(starting_test_position.x+(4*x_point),starting_test_position.y+(4*y_point)))
			y_point = y_point + 1
		elif len(collection_of_points) >= points_per_side * 2 and len(collection_of_points) < points_per_side * 3:
			### at the bottom : right to left
			collection_of_points.append(Vector2(starting_test_position.x+(4*x_point),starting_test_position.y+(4*y_point)))
			x_point = x_point - 1
		elif len(collection_of_points) >= points_per_side * 3 and len(collection_of_points) < number_of_points:
			### heading up : bottom to top
			collection_of_points.append(Vector2(starting_test_position.x+(4*x_point),starting_test_position.y+(4*y_point)))
			y_point = y_point - 1
		else:
			### mesh construction starts
			### along the top : left to right...
			collection_of_points.append(Vector2(starting_test_position.x+(4*x_point),starting_test_position.y+(4*y_point)))
			x_point = x_point + 1
		
	physical_matter.effigy.set_polygon(PackedVector2Array(collection_of_points))
	### vertex colors added...
	for point in physical_matter.effigy.get_polygon():
		coloring_points.append(Color(1.0,1.0,1.0))
	
	physical_matter.effigy.set_vertex_colors(PackedColorArray(coloring_points))
	
	"""
	var test_image = Image.new()
	test_image.create(100,100,false,Image.FORMAT_RGBA8)
	test_image.fill(Color(0,0,1.0))
	
	physical_matter.effigy_visual = ImageTexture.new()
	physical_matter.effigy_visual.create_from_image(test_image)
	
	physical_matter.effigy.set_texture(physical_matter.effigy_visual)
	
	physical_matter.effigy.set_uv(PackedVector2Array([Vector2(0,0),
	Vector2(50,0),
	Vector2(100,0),
	Vector2(100,500),
	Vector2(100,100),
	Vector2(0,0),
	Vector2(0,100)]))
	"""
	### the entire properties of the substance...
	### constitution/composition of the substance...
	#constitutive models to simulate...
		
	""" null-void
	physical_matter.coefficient_of_restitution = 1.0
	physical_matter.coefficient_of_static_friction = 0.5
	physical_matter.coefficient_of_kinetic_friction = 0.75
	physical_matter.physical_state = 'none'
	physical_matter.type_of_substance = 'void'
	physical_matter.constitutive_model = 'none'
	physical_matter.poisson_ratio = 0.0
	physical_matter.youngs_modulus = 0.0
	#"""
	#""" hyperelastic model
	physical_matter.coefficient_of_restitution = .9 #rubber
	physical_matter.coefficient_of_static_friction = 0.8 #rubber
	physical_matter.coefficient_of_kinetic_friction = 0.6 #rubber
	physical_matter.physical_state = 'solid'
	physical_matter.constitutive_model = 'hyperelastic'
	physical_matter.poisson_ratio = 0.5 #rubber
	physical_matter.youngs_modulus = .50 #natural rubber
	#"""
	"""
	# testing water model
	substance.coefficient_of_restitution = 0.2
	substance.coefficient_of_static_friction = 0.4
	substance.coefficient_of_kinetic_friction = 0.2
	substance.physical_state = 'liquid'
	substance.constitutive_model = 'water'
	#var flow = randf_range(-10.0,10.0)
	var flow = 100.0
	#maintain_velocity = Vector2(flow,flow)
	#substance.maintain_velocity = Vector2(0.0,flow)
	#substance.volume = 1.0
	#"""
	"""
	# testing fixed-corated model - snow...
	substance.coefficient_of_restitution = randf_range(0.53,1.76) # dry snow:(0.53-1.76) , wet snow:.30-.60 , 
	substance.coefficient_of_static_friction = 0.03 
	substance.coefficient_of_kinetic_friction = 0.015
	substance.physical_state = 'solid'
	substance.constitutive_model = 'fixed_corated'
	substance.type_of_substance = 'snow'
	substance.poisson_ratio = 0.2#0.5
	substance.youngs_modulus = snapped((1.4 * pow(10.0,5.0)),.1)
	substance.volume = snapped((4.0 * pow(10.0,2.0)),.1)
	#
	#"""
	"""
	# testing drucker_prager_elasticity - sand...
	substance.coefficient_of_restitution = randf_range(0.88,0.98) # wet sand (0.05,0.70)
	substance.coefficient_of_static_friction = randf_range(0.3,0.5) # wet sand 0.40 
	substance.coefficient_of_kinetic_friction = 0.15 #random
	substance.physical_state = 'solid'
	substance.constitutive_model = 'drucker_prager_elasticity'
	substance.type_of_substance = 'sand'
	substance.poisson_ratio = 0.29
	substance.youngs_modulus = snapped((3.537 * pow(10.0,7.0)),.1)
	#substance.volume = 1.0#snapped((4.0 * pow(10.0,2.0)),.1)
	#"""
		
	physical_matter.mass = 1.0
	physical_matter.volume = 1
	
	
	### Material Method Properties of each node/particle of substance...
	for node_count in range(0,len(physical_matter.effigy.get_polygon())):
		
		physical_matter.inner_workings['mass'] = physical_matter.mass / len(physical_matter.effigy.get_polygon())
		physical_matter.inner_workings['velocity'] = Vector2(0,0)
		physical_matter.inner_workings['initial velocity'] = Vector2(0,0)
		physical_matter.inner_workings['volume'] = physical_matter.volume / len(physical_matter.effigy.get_polygon())
		physical_matter.inner_workings['stress'] = [1.0,1.0,1.0,1.0]
		physical_matter.inner_workings['B'] = [0,0,0,0]
		physical_matter.inner_workings['C'] = [0,0,0,0]
		physical_matter.inner_workings['I'] = [1,0,0,1]
		physical_matter.inner_workings['F'] = [1,0,0,1]
		physical_matter.inner_workings['J'] = 0
		physical_matter.inner_workings['eulerian'] = {}
		physical_matter.inner_workings['euler data'] = {'mass': physical_matter.inner_workings['mass'],'velocity':Vector2(0.0,0.0),'momentum':[0,0]}.duplicate(true)
		
		physical_matter.mechanics[node_count] = physical_matter.inner_workings.duplicate(true)
	
	
	get_tree().get_root().get_node('Simulation/container').add_child(physical_matter)
	
