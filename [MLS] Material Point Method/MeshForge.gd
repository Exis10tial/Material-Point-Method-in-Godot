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




func polygon_form():
	physical_matter.effigy = Polygon2D.new()
	
	# create a square mesh...
	var number_of_points = 48 
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
	



func mesh_form():
	### Creation of Mesh
	
	var starting_test_position = Vector2(get_tree().get_root().size.x/2.0,get_tree().get_root().size.y/2.0)
	
	var intial_position = [Vector2(starting_test_position.x,starting_test_position.y),
	Vector2(starting_test_position.x+50,starting_test_position.y),
	Vector2(starting_test_position.x,starting_test_position.y+50)]

	
	### determines how many times the inital triangle is subdivide 
	### thru the centroid of the triangle
	
	# the base triangle is 0...
	var centroid__rank = 0
	### number of vertices is determined by the ranks of centroids...
	### for triangle base limit is 3...
	var number_of_vertices = 3
	
	if centroid__rank == 0:
		pass
	elif centroid__rank == 1:
		### number of vertices
		### +3 === median points of the triangle
		### +1 == centroid points of the triangle
		number_of_vertices = number_of_vertices + 3 + 1
	
	var construct_mesh = SurfaceTool.new()
	construct_mesh.begin(Mesh.PRIMITIVE_TRIANGLES)
	for vertex in intial_position:
		construct_mesh.set_color(Color(1,1,1))
		construct_mesh.add_vertex(Vector3(vertex.x,vertex.y,0))
	
	physical_matter.effigy = construct_mesh.commit()
	physical_matter.effigy_data.create_from_surface(physical_matter.effigy,0)
	



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
	physical_matter.effigy = Polygon2D.new()
	var collection_of_points = []
	var coloring_points = []
	var poly_triangle_parts = {'vertices':null,'centroid':null,'medians':null}
	var poly_triangles = {}
	
	var previous_poly_set = {}
	
	var starting_test_position = Vector2(get_tree().get_root().size.x/2.0,get_tree().get_root().size.y/2.0)
	
	###form a initial poly-triangle
	collection_of_points = [Vector2(starting_test_position.x,starting_test_position.y),
	Vector2(starting_test_position.x+100,starting_test_position.y),
	Vector2(starting_test_position.x,starting_test_position.y+100)
	]
	
	### different/seperate ways to break down an original polygon...
	var subdivide_only : bool = false
	var subdivide_into_pieces : bool = true
	
	### if both are false use the default given points...
	#both can be true:: subdivide into pieces first then add subdivide points...
	
	###if sudivide is true:: points are added inbetween points in the exact middle.
	### essentially doubling the amount of points: 3-6-12-24-48-96
	
	### if sudivide into pieces is true::the original polygon is parsed into their own smaller pieces.
	
	
	
	if subdivide_only == true and subdivide_into_pieces == false: 
		###...
		var subdivide_rank = 1
		var median_of_ = Vector2(0,0)  
		var number_of_vertices = 3
		
		poly_triangle_parts['vertices'] = collection_of_points
		
		poly_triangles[0] = poly_triangle_parts
		
		### establish how many total vertices will be created...
		var cycles = 0
		
		while true:
			if cycles >= subdivide_rank:
				break
			
			number_of_vertices = number_of_vertices * 2
			
			cycles = cycles + 1
		
		while len(poly_triangles[0]['vertices']) < number_of_vertices :
			previous_poly_set = poly_triangles.duplicate(true)
			poly_triangles = {}

			### process of establishing median points..
			#var poly_index = 0
			#while true:
				#if poly_index >= len(previous_poly_set.keys()):
					#break
				
			var polygon = previous_poly_set[0]
				
				#median_of_a = Vector2( (previous_poly_set[polygon]['vertices'][0].x + previous_poly_set[polygon]['vertices'][1].x) / 2.0, (previous_poly_set[polygon]['vertices'][0].y + previous_poly_set[polygon]['vertices'][1].y) / 2.0)
				#median_of_b = Vector2( (previous_poly_set[polygon]['vertices'][1].x + previous_poly_set[polygon]['vertices'][2].x) / 2.0, (previous_poly_set[polygon]['vertices'][1].y + previous_poly_set[polygon]['vertices'][2].y) / 2.0)
				#median_of_c = Vector2( (previous_poly_set[polygon]['vertices'][2].x + previous_poly_set[polygon]['vertices'][0].x) / 2.0, (previous_poly_set[polygon]['vertices'][2].y + previous_poly_set[polygon]['vertices'][0].y) / 2.0)
				
				#previous_poly_set[polygon]['medians'] = [median_of_a,median_of_b,median_of_c]
				
				#poly_index = wrapi(poly_index+1,0,len(previous_poly_set.keys())+1)
			var part_c = 0
			var part_d = 1
			var build_medians = []
			while true:
				
				if part_c >= len(polygon['vertices']) :
					break
					
				median_of_ =  Vector2( (polygon['vertices'][part_c].x + polygon['vertices'][part_d].x) / 2.0, (polygon['vertices'][part_c].y + polygon['vertices'][part_d].y) / 2.0)
					
				build_medians.append(median_of_)
				part_c = part_c + 1
				part_d = wrapi(part_d + 1,0,len(polygon['vertices']))
			
			previous_poly_set[0]['medians'] = build_medians.duplicate(true)
			
			### establish a create new vertices for polygons...
			polygon = previous_poly_set[0]
				
			var created_count = 0
			var build_vertices = []
			var part_a = 'vertices'
			var part_b = 'medians'
			var placeholder_a
			var placeholder_b
			part_c = 0
			part_d = 0
				
			while true:
				if len(build_vertices) >= number_of_vertices:
					break
					
				if created_count >= 2:
					part_c = part_c + 1
					created_count = 0
					
				build_vertices.append(polygon[part_a][part_c])
					
				placeholder_a = part_a
				placeholder_b = part_b
				part_a = placeholder_b
				part_b = placeholder_a
					
				created_count = created_count + 1
					
			poly_triangle_parts['vertices'] = build_vertices.duplicate(true)
			poly_triangles[len(poly_triangles.keys())] = poly_triangle_parts.duplicate(true)
					
	
			### the establish polygon from new vertices...
			var color
			for polygon_shape_index in range(0,len(poly_triangles)):
				physical_matter.effigy = Polygon2D.new()
				coloring_points = []
				
				physical_matter.effigy.set_polygon(PackedVector2Array(poly_triangles[polygon_shape_index]['vertices']))
				
				for vertex in poly_triangles[polygon_shape_index]['vertices']:
					coloring_points.append(Color(1.0,1.0,1.0))
		
				physical_matter.effigy.set_vertex_colors(PackedColorArray(coloring_points))
		
				
				physical_matter.effigy_basket.append(physical_matter.effigy)


	elif subdivide_only == false and subdivide_into_pieces == true:
		###...
		###...
		var centroid_rank = 1
		var centroid_point = Vector2(0,0)
		var median_of_a = Vector2(0,0) 
		var median_of_b = Vector2(0,0) 
		var median_of_c = Vector2(0,0) 
		var number_of_polygons = 1
		
		poly_triangle_parts['vertices'] = collection_of_points
		
		poly_triangles[0] = poly_triangle_parts
		
		###.
		
		### establish how many total polygons will be created...
		var cycles = 0
		
		while true:
			if cycles >= centroid_rank:
				break
			
			number_of_polygons = number_of_polygons * 6
			
			cycles = cycles + 1
		
		### mechanics/process of creating new polygons
		### by subdividing existing polygon(triangles)...
		while true:
			if len(poly_triangles) >= number_of_polygons:
				break
		
			previous_poly_set = poly_triangles.duplicate(true)
			poly_triangles = {}

			### process of establishing centroid and median points..
			var poly_index = 0
			while true:
				if poly_index >= len(previous_poly_set.keys()):
					break
				
				var polygon = previous_poly_set.keys()[poly_index]
				
				centroid_point = Vector2( (previous_poly_set[polygon]['vertices'][0].x + previous_poly_set[polygon]['vertices'][1].x + previous_poly_set[polygon]['vertices'][2].x) / 3.0 ,
				(previous_poly_set[polygon]['vertices'][0].y + previous_poly_set[polygon]['vertices'][1].y + previous_poly_set[polygon]['vertices'][2].y)  / 3.0 )
				
				median_of_a = Vector2( (previous_poly_set[polygon]['vertices'][0].x + previous_poly_set[polygon]['vertices'][1].x) / 2.0, (previous_poly_set[polygon]['vertices'][0].y + previous_poly_set[polygon]['vertices'][1].y) / 2.0)
				median_of_b = Vector2( (previous_poly_set[polygon]['vertices'][1].x + previous_poly_set[polygon]['vertices'][2].x) / 2.0, (previous_poly_set[polygon]['vertices'][1].y + previous_poly_set[polygon]['vertices'][2].y) / 2.0)
				median_of_c = Vector2( (previous_poly_set[polygon]['vertices'][2].x + previous_poly_set[polygon]['vertices'][0].x) / 2.0, (previous_poly_set[polygon]['vertices'][2].y + previous_poly_set[polygon]['vertices'][0].y) / 2.0)
				
				previous_poly_set[polygon]['centroid'] = Vector2(int(centroid_point.x),int(centroid_point.y))
				previous_poly_set[polygon]['medians'] = [median_of_a,median_of_b,median_of_c]
				
				poly_index = wrapi(poly_index+1,0,len(previous_poly_set.keys())+1)
			
			### establish a create new vertices for polygons...
			poly_index = 0
			while true:
				if poly_index >= len(previous_poly_set.keys()):
					#print('leaving create new vertices ')
					break
				
				var polygon = previous_poly_set.keys()[poly_index]
				
				var created_count = 0
				var part_a = 'vertices'
				var part_b = 'medians'
				var placeholder_a
				var placeholder_b
				var part_c = 0
				var part_d = 0
				while true:
					if created_count >= 6:
						break
					
					if created_count > 0:
						if created_count in range(1,7,2):
							part_d = wrapi(part_d+1,0,3)
						else:
							part_c = wrapi(part_c+1,0,3)
						
					poly_triangle_parts['vertices'] = [ previous_poly_set[polygon][part_a][part_c],previous_poly_set[polygon]['centroid'],previous_poly_set[polygon][part_b][part_d] ].duplicate(true)
					poly_triangle_parts['centroid'] = null
					poly_triangle_parts['medians'] = []
					
					poly_triangles[len(poly_triangles.keys())] = poly_triangle_parts.duplicate(true)
					
					placeholder_a = part_a
					placeholder_b = part_b
					part_a = placeholder_b
					part_b = placeholder_a
					
					created_count = created_count + 1
					
				poly_index = wrapi(poly_index+1,0,len(previous_poly_set.keys())+1)
			
			### the establish polygon from new vertices...
			var color
			for polygon_shape_index in range(0,len(poly_triangles)):
				physical_matter.effigy = Polygon2D.new()
				coloring_points = []
				
				physical_matter.effigy.set_polygon(PackedVector2Array(poly_triangles[polygon_shape_index]['vertices']))
				
				for vertex in poly_triangles[polygon_shape_index]['vertices']:
					coloring_points.append(Color(1.0,1.0,1.0))
		
				physical_matter.effigy.set_vertex_colors(PackedColorArray(coloring_points))
		
				
				physical_matter.effigy_basket.append(physical_matter.effigy)
	
	elif subdivide_only == true and subdivide_into_pieces == true:
		pass
		
	elif subdivide_only == false and subdivide_into_pieces == false:
		
		physical_matter.effigy.set_polygon(PackedVector2Array(collection_of_points))
		
		### vertex colors added...
		for point in physical_matter.effigy.get_polygon():
			coloring_points.append(Color(1.0,1.0,1.0))
		
		physical_matter.effigy.set_vertex_colors(PackedColorArray(coloring_points))
		
		physical_matter.effigy_basket.append(physical_matter.effigy)
		

	### polygon/particle mechanics...
	var polygon_index = 0
	var particle_index = 0
	var particle_carrier = []
	while true:
		### assoicate polygons to particles.
		if polygon_index == len(physical_matter.effigy_basket):
			break
		
		var polygon = physical_matter.effigy_basket[polygon_index]
		particle_carrier = []
		
		while true:
			if len(particle_carrier) == len(polygon.get_polygon()):
				
				physical_matter.associate_polygon_to_particle[polygon_index] = particle_carrier

				break
			
			particle_carrier.append(particle_index)
			
			particle_index = wrapi(particle_index+1,0,particle_index+2)
		
		###loop 
		polygon_index = wrapi(polygon_index+1,0,len(physical_matter.effigy_basket)+1)
		
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
		
	#""" null-void
	physical_matter.coefficient_of_restitution = 1.0
	physical_matter.coefficient_of_static_friction = 0.5
	physical_matter.coefficient_of_kinetic_friction = 0.75
	physical_matter.physical_state = 'none'
	physical_matter.type_of_substance = 'void'
	physical_matter.constitutive_model = 'none'
	physical_matter.poisson_ratio = 0.0
	physical_matter.youngs_modulus = 0.0
	#"""
	""" hyperelastic model
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
	
	var identify_effigy = 0
	
	while true:
	
		if identify_effigy == len(physical_matter.effigy_basket):
			break
		
		### this is the polygon...
		var effigy_parts = physical_matter.associate_polygon_to_particle[identify_effigy]
		var effigy = physical_matter.effigy_basket[identify_effigy]
		var identify_particle = 0
		
		### cycle thru the particles(vertex) of the polygon...
		while true:
			if identify_particle == len(effigy_parts):
				break
			
			#var particle = effigy.get_polygon()[identify_particle]
			var particle_indicator = effigy_parts[identify_particle]
			var particle = effigy.get_polygon()[identify_particle]
			
			physical_matter.inner_workings['mass'] = physical_matter.mass / len(physical_matter.associate_polygon_to_particle.keys())
			physical_matter.inner_workings['velocity'] = Vector2(0.0,0.0)
			physical_matter.inner_workings['initial velocity'] = Vector2(0,0)
			physical_matter.inner_workings['volume'] = physical_matter.volume / len(physical_matter.associate_polygon_to_particle.keys())
			physical_matter.inner_workings['stress'] = [1.0,1.0,1.0,1.0]
			physical_matter.inner_workings['B'] = [0,0,0,0]
			physical_matter.inner_workings['C'] = [0,0,0,0]
			physical_matter.inner_workings['I'] = [1,0,0,1]
			physical_matter.inner_workings['F'] = [1,0,0,1]
			physical_matter.inner_workings['J'] = 0
			physical_matter.inner_workings['eulerian'] = {}
			physical_matter.inner_workings['euler data'] = {'mass': physical_matter.inner_workings['mass'],'velocity':Vector2(0.0,0.0),'momentum':[0,0]}.duplicate(true)
			physical_matter.inner_workings['correspond'] = identify_effigy
			
			physical_matter.inner_workings['relation within'] = effigy.get_polygon().find(particle)
			
			physical_matter.mechanics[particle_indicator] = physical_matter.inner_workings.duplicate(true)
			
			identify_particle = wrapi(identify_particle+1,0,len(effigy_parts)+1)
		
		identify_effigy = wrapi(identify_effigy+1,0,len(physical_matter.effigy_basket)+1)
		
	
	get_tree().get_root().get_node('Simulation/container').add_child(physical_matter)
	
