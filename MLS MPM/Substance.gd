extends Node


var gather_the_particles : File
var particle : Object
var default_mass_of_particle : float
var maintain_velocity : Vector2
var particle_name : String = ''
var letters_list : Array = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
var digits_list : Array = ['0','1','2','3','4','5','6','7','8','9']
var kernel_distance : Vector2
var kernel_x : float
var kernel_y : float
var flipped_kernel_distance : Vector2
var flipped_kernel_x : float
var flipped_kernel_y : float
var dx : float
var domain : Rect2
var domain_size : Vector2
var gathered_into_chunks : bool = false
var one_particle : bool = false
var x_ : int
var y_ : int
var prime_set_of_numbers : Array = []
var powers_of_2_set_of_numbers : Array = []
var find_x_factors : Array = []
var find_y_factors : Array = []
var defined_rows : int
var defined_columns : int
var particle_limit : int
var cell_size : float 
var substance_starting_point : Vector2
var grid_nodes : Dictionary



func Initial_Collection_Of_Substance():
	### the initial collection of the substance before the simulation...
	
	
	### note ::
	
	#x:100 y:100 = 100 particles
	#x:35 y:30 = 42 particles
	#x:25 y:30 = 30 particles
	#x:140 y:80 = 28 particles
	#x:60 y:80 = 12 particles
	#x:4 y:4 = 2 particles
	#x:1 y:1 = 1 particle
	
	#domain_size = Vector2(100.0,100.0)
	#domain_size = Vector2(48.0,42.0)
	#domain_size = Vector2(35.0,30.0)
	#domain_size = Vector2(25.0,25.0)
	domain_size = Vector2(16.0,16.0)
	#domain_size = Vector2(10.0,10.0)
	#domain_size = Vector2(9.0,9.0)
	#domain_size = Vector2(6.0,6.0)
	#domain_size = Vector2(5.0,5.0)
	#domain_size = Vector2(4.0,4.0)
	#domain_size = Vector2(1.0,2.0)
	#domain_size = Vector2(1.0,1.0)
	#---------------------------------
	# if particles size is always 1...
	#domain_size = Vector2(2,10.0)
	#domain_size = Vector2(10.0,10.0)
	#domain_size = Vector2(100.0,100.0)
	return domain_size


func Cross_Section_Of_Substance(into_pieces:bool):
	### the size of the substance is parsed into evenly cut squares/cubes = (Length == Width)
	
	###if code require into_pieces == false you will get particle at domain size.
	if domain_size.x == domain_size.y and domain_size.x != 1 and domain_size.y != 1 and into_pieces == true :
		
		### Powers of 2 Check...
		var n = 0.0
		
		while pow(2,n) <= domain_size.x:
			
			powers_of_2_set_of_numbers.append(pow(2,n))
			
			n +=1
		
		if powers_of_2_set_of_numbers.has(domain_size.y) == true and into_pieces == true:
			
			cell_size = int(domain_size.x / pow(2,1))
		else:
			### the number is not a power of 2...
			
			#Square Root Check
			if fmod(domain_size.x,sqrt(domain_size.x)) == 0 and into_pieces == true:
				cell_size = sqrt(domain_size.x)
			else:
				### Not a power of 2 or a has a square root...
				
				# find the greatest common factor between the domain sizes::
				for number in range(1,domain_size.x+1):
					if int(domain_size.x) % number == 0:
						find_x_factors.append(number)
				for number in range(1,domain_size.y+1):
					if int(domain_size.y) % number == 0:
						find_y_factors.append(number)
				
				### using the greatest common factor to acquire a cell size...
				if domain_size.x == domain_size.y:
					for number in range(len(find_x_factors)-1,-1,-1):
						if find_x_factors[number] in find_y_factors:
							cell_size = find_x_factors[number] 
							break
				elif domain_size.x < domain_size.y:
					for number in range(len(find_x_factors)-1,-1,-1):
						if find_x_factors[number] in find_y_factors:
							cell_size =  find_x_factors[number] 
							break
				elif domain_size.x > domain_size.y:
					for number in range(len(find_y_factors)-1,-1,-1):
						if find_y_factors[number] in find_x_factors:
							cell_size = find_y_factors[number]
							break
											
	else:
		if into_pieces == false :
			cell_size = 1
		else:
			#"""
			# find the greatest common factor between the domain sizes::
			for number in range(1,domain_size.x+1):
				if int(domain_size.x) % number == 0:
					find_x_factors.append(number)
			for number in range(1,domain_size.y+1):
				if int(domain_size.y) % number == 0:
					find_y_factors.append(number)
			
			### using the greatest common factor to acquire a cell size...
			if domain_size.x == domain_size.y:
				for number in range(len(find_x_factors)-1,-1,-1):
					if find_x_factors[number] in find_y_factors:
						cell_size = find_x_factors[number] 
						break
			elif domain_size.x < domain_size.y:
				for number in range(len(find_x_factors)-1,-1,-1):
					if find_x_factors[number] in find_y_factors:
						cell_size =  find_x_factors[number] 
						break
			elif domain_size.x > domain_size.y:
				for number in range(len(find_y_factors)-1,-1,-1):
					if find_y_factors[number] in find_x_factors:
						cell_size = find_y_factors[number]
						break
			#"""
	return float(cell_size)
	


func Establish_Rate():
	### ...
	#dx = snapped(1.0/float(particle_limit)*1.50,.001)
	#dx = snapped(1.0/float(particle_limit)*5.00,.001)
	dx = snapped(1.0/float(particle_limit),.001)
	#dx = snapped(1.0/1.50,.001)
	#dx = 1.50
	#dx = 1.0
	return dx

func _on_Substance_ready():
	gather_the_particles = File.new()
	
	#if gather_the_particles.file_exists("res://2D Particle.tscn") or gather_the_particles.file_exists("res://2D Particle(thru_rigid_body).tscn"):
	if gather_the_particles.file_exists("res://2D Particle.tscn"):
	#if gather_the_particles.file_exists("res://2D Particle(thru_rigid_body).tscn"):
	#if gather_the_particles.file_exists("res://2D Particle(thru sprite).tscn"):
		Initial_Collection_Of_Substance()
		
		### alters how many particles ...
		# one_particle is 1 particle no matter the domain size...
		# default: Is using the greatest common factor between domain_size x and domain_size.y to acquire a number of particles/cell size..
		# gathered_into_chunks: First a if both domain_size are a power of 2 , then square root is check , then to the default...
		# if one particle and gathered into chunks is false, number of particles is domain_size.x * domain_size.y
		#one_particle = true
		gathered_into_chunks = true
		
		
		if one_particle == true:
			
			defined_columns = int(domain_size.x / domain_size.x)
			defined_rows = int(domain_size.y / domain_size.y)
			
			particle_limit = int(defined_columns * defined_rows)
			cell_size = 1
			
		else:
			cell_size = Cross_Section_Of_Substance(gathered_into_chunks)
			#print(cell_size,' cell size check')
			
			#25x25 : 1 fps
			#10x10 : 13 fps
			#7x8: 33 fps
			#7x5 : 60 fps
			#5x5 : 60 fps
			
			###
			defined_columns = int(domain_size.x / cell_size)
			defined_rows = int(domain_size.y / cell_size)
				
			##print(defined_columns,' checking columns')
			##print(defined_rows,' checking rows')
				
			### determines the size/shape of particles...
			if gathered_into_chunks == false:
				### use if particles to be size 1 always...
				particle_limit = int(defined_columns * defined_rows)
					
			else:
				### the particle will end up as power of two or square root cut...
				particle_limit = int(defined_columns * defined_rows)# / 2
		
		if one_particle == true:
			substance_starting_point = Vector2(
			((ProjectSettings.get_setting('display/window/size/width') / 2.0)),
			((ProjectSettings.get_setting('display/window/size/height') / 2.0))
			)
		else:
			substance_starting_point = Vector2(
			((ProjectSettings.get_setting('display/window/size/width') / 2.0) - ((cell_size * defined_columns) / 2.0) ) + (cell_size/2.0),
			((ProjectSettings.get_setting('display/window/size/height') / 2.0) - ((cell_size * defined_rows) / 2.0)) + (cell_size/2.0)
			)
		
		while get_child_count() < particle_limit:
			###
			
			particle = load("res://2D Particle.tscn").instantiate()
			#particle = load("res://2D Particle(thru_rigid_body).tscn").instantiate()
			#particle = load("res://2D Particle(thru sprite).tscn").instantiate()
			
			
			
			particle_name = '{a}{b}{c}{d}{1}{2}{3}{4}'.format({
				'a':letters_list[int(randi_range(0,len(letters_list)-1))],
				'b':letters_list[int(randi_range(0,len(letters_list)-1))],
				'c':letters_list[int(randi_range(0,len(letters_list)-1))],
				'd':letters_list[int(randi_range(0,len(letters_list)-1))],
				'1':digits_list[int(randi_range(0,len(digits_list)-1))],
				'2':digits_list[int(randi_range(0,len(digits_list)-1))],
				'3':digits_list[int(randi_range(0,len(digits_list)-1))],
				'4':digits_list[int(randi_range(0,len(digits_list)-1))]
				})
			
			particle.set_name(particle_name)
			
			if one_particle == true:
				x_ = domain_size.x
				y_ = domain_size.y
			else:
				x_ = 0
				y_ = 0
			
			particle.adjust_size(cell_size,one_particle,x_,y_)
			
			
			#default_mass_of_particle = 1.00 #/ particle_limit
			default_mass_of_particle = 1.00
			maintain_velocity = Vector2(0.0,0.0)
			#maintain_velocity = Vector2(randf_range(-10.0,10.0),randf_range(-10.0,10.0))
			#maintain_velocity = Vector2(5.0,0.0)
			#maintain_velocity = Vector2(-0.0,0.0)
			#maintain_velocity = Vector2(0.0,-10.0)
			#maintain_velocity = Vector2(0.0,10.0)
			particle.mass = default_mass_of_particle
			particle.stress = [[1.0,1.0],[1.0,1.0]]
			particle.velocity = maintain_velocity
			#"""
			# testing hyperelastic - neohookean
			particle.coefficient_of_restitution = 1.0 #rubber
			particle.coefficient_of_static_friction = 0.9 #rubber
			particle.coefficient_of_kinetic_friction = 0.250 #rubber
			particle.physical_state = 'solid'
			particle.constitutive_model = 'hyperelastic'
			particle.poisson_ratio = 0.5 #rubber
			particle.youngs_modulus = 0.1 #rubber
			particle.volume = 1.0
			#"""
			"""
			# testing water model
			particle.coefficient_of_restitution = 0.2
			particle.coefficient_of_static_friction = 0.4
			particle.coefficient_of_kinetic_friction = 0.2
			particle.physical_state = 'liquid'
			particle.constitutive_model = 'water'
			#var flow = randf_range(-10.0,10.0)
			var flow = 100.0
			#maintain_velocity = Vector2(flow,flow)
			maintain_velocity = Vector2(0.0,flow)
			particle.volume = 1.0
			"""
			"""
			# testing fixed-corated model - snow...
			particle.coefficient_of_restitution = randf_range(0.53,1.76) # dry snow:(0.53-1.76) , wet snow:.30-.60 , 
			particle.coefficient_of_static_friction = 0.03 
			particle.coefficient_of_kinetic_friction = 0.015
			particle.physical_state = 'solid'
			particle.constitutive_model = 'fixed_corated'
			particle.type_of_substance = 'snow'
			particle.poisson_ratio = 0.5
			particle.youngs_modulus = snapped((1.4 * pow(10.0,5.0)),.1)
			particle.volume = 1.0#snapped((4.0 * pow(10.0,2.0)),.1)
			#
			#print(particle.F,' checking F')
			"""
			"""
			# testing fixed-corated model - sand...
			particle.coefficient_of_restitution = randf_range(0.88,0.98) # wet sand (0.05,0.70)
			particle.coefficient_of_static_friction = randf_range(0.3,0.5) # wet sand 0.40 
			particle.coefficient_of_kinetic_friction = 0.15 #random
			particle.physical_state = 'solid'
			particle.constitutive_model = 'drucker_prager_elasticity'
			particle.type_of_substance = 'sand'
			particle.poisson_ratio = 0.29
			particle.youngs_modulus = snapped((3.537 * pow(10.0,7.0)),.1)
			particle.volume = 1.0#snapped((4.0 * pow(10.0,2.0)),.1)
			# Polar SVD :   U from AA^T, V from A^TA...
			var FFtransposed = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Matrix(particle.F,get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Transposed_Matrix(particle.F))
			particle.U = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Find_Eigenvectors(FFtransposed)
			var FtransposedF = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Matrix(particle.F,get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Transposed_Matrix(particle.F))
			particle.V = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Find_Eigenvectors(FtransposedF)
			#Diagonalize SIGMA = PAP^-1... P = diagonalize_helper ,A = artifact.F
			var diagonalize_helper = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Find_Eigenvectors(particle.F)
			# 
			particle.Sigma = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Matrix(get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Matrix(diagonalize_helper,particle.F),get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Inverse_Matrix(diagonalize_helper))
			
			#
			#print(particle.F,' checking F')
			#"""
			### identify other particles with its domain ... itself always...
			particle.within_range.append(particle)
		
			add_child(particle)
		
		var n = 0
		for x in range(1,defined_rows+1):
			for y in range(1,defined_columns+1):
				if one_particle == true:
					get_child(n).set_position(Vector2(substance_starting_point.x + (domain_size.x * y),substance_starting_point.y + (domain_size.y * x)) )
				else:
					get_child(n).set_position(Vector2(substance_starting_point.x + (cell_size * y),substance_starting_point.y + (cell_size * x)) )
				n += 1
		
		### establish the erase board grid
		for created_particle in get_children():
			grid_nodes[created_particle] = {'mass':0.0,'velocity':Vector2(0,0),'momentum':Vector2(0.0,0.0)}
		
		### establish the domain...
		for particle in get_children():
			# establish the particle grid domain...
			var area_multiplier = 1.0
			#particle.surrounding_area = Rect2(Vector2(particle.position.x - ((particle.get_node("shape").get_size().x/2.0)*area_multiplier),particle.position.y - ((particle.get_node("shape").get_size().y/2.0)*area_multiplier)),Vector2(particle.get_node("shape").get_size().x*area_multiplier,particle.get_node("shape").get_size().y*area_multiplier))
			particle.surrounding_area = Rect2(Vector2((particle.position.x+particle.get_node("shape").get_position().x),(particle.position.y+particle.get_node("shape").get_position().y)),Vector2(particle.get_node("shape").get_size().x*area_multiplier,particle.get_node("shape").get_size().y*area_multiplier))
			
			#print(particle.position,' particle positon')
			#print(particle.get_node("shape").get_position())
			#print(particle.surrounding_area.position,' surrounding area')
			
		### establish particle relation to others...
		for particle in get_children():
			for other_particle in get_children():
				
				#kernel_distance = (particle.get_position() - other_particle.surrounding_area.get_center())# / cell_size
				kernel_x = snapped(particle.get_position().x - other_particle.surrounding_area.get_center().x,.01)
				kernel_y = snapped(particle.get_position().y - other_particle.surrounding_area.get_center().y,.01)
				kernel_distance = Vector2(kernel_x,kernel_y)
				#flipped_kernel_distance = (other_particle.surrounding_area.get_center() - particle.get_position())# / cell_size
				flipped_kernel_x = snapped(other_particle.surrounding_area.get_center().x - particle.get_position().x,.01)
				flipped_kernel_y = snapped(other_particle.surrounding_area.get_center().y - particle.get_position().y,.01)
				flipped_kernel_distance = Vector2(flipped_kernel_x,flipped_kernel_y)
				
				particle.relation_to_domain[other_particle] = kernel_distance
				particle.domain_relation_to_particle[other_particle] = flipped_kernel_distance
				#particle.relation_to_domain[other_particle] =flipped_kernel_distance
				#particle.domain_relation_to_particle[other_particle] =kernel_distance
		
		### particles volume...
		
		
		Establish_Rate()

	else:
		print('particles Does not exists')
		pass

