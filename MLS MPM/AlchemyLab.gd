extends Node


#var pile_of_substances : File
var substance : Object
var default_mass_of_substance : float
var maintain_velocity : Vector2
var substance_particle_name : String = ''
var letters_list : Array = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
var digits_list : Array = ['0','1','2','3','4','5','6','7','8','9']
var location_x :int = 0
var location_y :int = 0
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
var one_substance : bool = false
var x_ : int
var y_ : int
var prime_set_of_numbers : Array = []
var powers_of_2_set_of_numbers : Array = []
var find_x_factors : Array = []
var find_y_factors : Array = []
var defined_rows : int
var defined_columns : int
var substance_limit : int
var cell_size : float 
var substance_starting_point : Vector2
var grid_nodes : Dictionary



func Initial_Collection_Of_Substance():
	### the initial collection of the substance before the simulation...
	# Initial Shape of the substance also determine the number of total number of substances of the substance...
	#
	
	#x:100 y:100 = 100 substances
	#x:35 y:30 = 42 substances
	#x:25 y:30 = 30 substances
	#x:140 y:80 = 28 substances
	#x:60 y:80 = 12 substances
	#x:4 y:4 = 2 substances
	#x:1 y:1 = 1 substance
	
	#domain_size = Vector2(1000.0,1000.0)
	#domain_size = Vector2(100.0,100.0)
	#domain_size = Vector2(81.0,81.0)
	#domain_size = Vector2(50.0,50.0)
	#domain_size = Vector2(35.0,30.0)
	#domain_size = Vector2(25.0,25.0)
	#domain_size = Vector2(16.0,16.0)
	domain_size = Vector2(10.0,10.0)
	#domain_size = Vector2(9.0,9.0)
	#domain_size = Vector2(6.0,6.0)
	#domain_size = Vector2(5.0,5.0)
	#domain_size = Vector2(4.0,4.0)
	#domain_size = Vector2(3.0,3.0)
	#domain_size = Vector2(2.0,2.0)
	#domain_size = Vector2(1.0,1.0)
	#---------------------------------
	# if substances size is always 1...
	#domain_size = Vector2(2,10.0)
	#domain_size = Vector2(10.0,10.0)
	#domain_size = Vector2(100.0,100.0)
	return domain_size


func Cross_Section_Of_Substance(into_pieces:bool):
	### the size of the substance is parsed into evenly cut squares/cubes = (Length == Width)
	
	###if code require into_pieces == false you will get substance at domain size.
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
	



func _on_alchemy_lab_ready():
	### the building/construction of substance...
	#pile_of_substances = File.new()
	
	#if pile_of_substances.file_exists("res://Substance.tscn") or pile_of_substances.file_exists("res://Substance(thru_rigid_body).tscn"):
	if FileAccess.file_exists("res://Substance.tscn"):
	#if pile_of_substances.file_exists("res://Substance(thru_rigid_body).tscn"):
	#if pile_of_substances.file_exists("res://Substance(thru sprite).tscn"):
		Initial_Collection_Of_Substance()
		
		### alters how many substances ...
		# one_substance is 1 substance no matter the domain size...
		# default ::one_substance = false , gathered_into_chunks = true:: Is using the greatest common factor between domain_size x and domain_size.y to acquire a number of substances/cell size..
		# gathered_into_chunks: First a if both domain_size are the same first check power of 2 , then square root is check , then to the default...
		# if one substance and gathered into chunks is false, number of substances is domain_size.x * domain_size.y
		#one_substance = true
		#gathered_into_chunks = true
		
		if one_substance == true:
			gathered_into_chunks = false
			defined_columns = int(domain_size.x / domain_size.x)
			defined_rows = int(domain_size.y / domain_size.y)
			
			substance_limit = int(defined_columns * defined_rows)
			cell_size = 1
			
		else:
			one_substance = false
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
				
			### determines the size/shape of substances...
			#if gathered_into_chunks == false:
				### use if substances to be size 1 always...
			#	substance_limit = int(defined_columns * defined_rows)
					
			#else:
			### the substance will end up as power of two or square root cut...
			substance_limit = int(defined_columns * defined_rows)# / 2
		
		#var n = 0
		#for x in range(1,defined_rows+1):
		#	for y in range(1,defined_columns+1):
		#		if one_substance == true:
					#get_child(n).set_position(Vector2(substance_starting_point.x + (domain_size.x * y),substance_starting_point.y + (domain_size.y * x)) )
		#			substance.surrounding_area = Rect2(Vector2(substance_starting_point.x + (domain_size.x * y),substance_starting_point.y + (domain_size.y * x)),Vector2(domain_size.x,domain_size.y) )
		#		else:
					#get_child(n).set_position(Vector2(substance_starting_point.x + (cell_size * y),substance_starting_point.y + (cell_size * x)) )
		#			substance.surrounding_area = Rect2(Vector2(substance_starting_point.x + (cell_size * y),substance_starting_point.y + (cell_size * x)),Vector2(cell_size,cell_size) )
		#		n += 1
			
		if one_substance == true or substance_limit == 1:
			substance_starting_point = Vector2(
			((ProjectSettings.get_setting('display/window/size/width') / 2.0)),
			((ProjectSettings.get_setting('display/window/size/height') / 2.0))
			)
			#print(substance_starting_point,' substance_starting_point a')
		else:
			#sbstance_starting_point = Vector2(
			#((ProjectSettings.get_setting('display/window/size/width') / 2.0) - ((cell_size * defined_columns) / 2.0) ) + (cell_size/2.0),
			#((ProjectSettings.get_setting('display/window/size/height') / 2.0) - ((cell_size * defined_rows) / 2.0)) + (cell_size/2.0)
			#)
			substance_starting_point = Vector2(
			((ProjectSettings.get_setting('display/window/size/width') / 2.0) - ((cell_size * defined_columns) / 2.0) ),
			((ProjectSettings.get_setting('display/window/size/height') / 2.0) - ((cell_size * defined_rows) / 2.0))
			)
		
		substance = load("res://Substance.tscn").instantiate()
		
		substance.default_mass_of_substance = 1.0
		substance.maintain_velocity = Vector2(0.0,0.0)
		
		#substance.mass = default_mass_of_substance
		#substance.velocity = maintain_velocity
		#substance.stress = [[1.0,1.0],[1.0,1.0]]
		
		if one_substance == true:
			substance.appearance = domain_size.x
		else:
			substance.appearance = cell_size
		
		substance.coefficient_of_restitution = 1.0 #rubber
		substance.coefficient_of_static_friction = 0.9 #rubber
		substance.coefficient_of_kinetic_friction = 0.250 #rubber
		substance.physical_state = 'solid'
		substance.constitutive_model = 'hyperelastic'
		substance.poisson_ratio = 0.5 #rubber
		substance.youngs_modulus = 0.1 #rubber
		substance.volume = 1.0
		
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
		maintain_velocity = Vector2(0.0,flow)
		substance.volume = 1.0
		"""
		"""
		# testing fixed-corated model - snow...
		substance.coefficient_of_restitution = randf_range(0.53,1.76) # dry snow:(0.53-1.76) , wet snow:.30-.60 , 
		substance.coefficient_of_static_friction = 0.03 
		substance.coefficient_of_kinetic_friction = 0.015
		substance.physical_state = 'solid'
		substance.constitutive_model = 'fixed_corated'
		substance.type_of_substance = 'snow'
		substance.poisson_ratio = 0.5
		substance.youngs_modulus = snapped((1.4 * pow(10.0,5.0)),.1)
		substance.volume = 1.0#snapped((4.0 * pow(10.0,2.0)),.1)
		#
		"""
		"""
		# testing fixed-corated model - sand...
		substance.coefficient_of_restitution = randf_range(0.88,0.98) # wet sand (0.05,0.70)
		substance.coefficient_of_static_friction = randf_range(0.3,0.5) # wet sand 0.40 
		substance.coefficient_of_kinetic_friction = 0.15 #random
		substance.physical_state = 'solid'
		substance.constitutive_model = 'drucker_prager_elasticity'
		substance.type_of_substance = 'sand'
		substance.poisson_ratio = 0.29
		substance.youngs_modulus = snapped((3.537 * pow(10.0,7.0)),.1)
		substance.volume = 1.0#snapped((4.0 * pow(10.0,2.0)),.1)
		# Polar SVD :   U from AA^T, V from A^TA...
		var FFtransposed = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Matrix(substance.F,get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Transposed_Matrix(substance.F))
		substance.U = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Find_Eigenvectors(FFtransposed)
		var FtransposedF = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Matrix(substance.F,get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Transposed_Matrix(substance.F))
		substance.V = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Find_Eigenvectors(FtransposedF)
		#Diagonalize SIGMA = PAP^-1... P = diagonalize_helper ,A = artifact.F
		var diagonalize_helper = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Find_Eigenvectors(substance.F)
		# 
		substance.Sigma = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Matrix(get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Matrix(diagonalize_helper,substance.F),get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Inverse_Matrix(diagonalize_helper))
		"""
		#while len(substance.particle_mechanics) < substance_limit:
		for x in range(substance_limit):
			### for the position of the starting position of the substance drawn...
			if location_x == domain_size.x/cell_size:
				location_y =  location_y + 1
				location_x = 0
			
			substance_particle_name = '{a}{b}{c}{d}{1}{2}{3}{4}'.format({
				'a':letters_list[int(randi_range(0,len(letters_list)-1))],
				'b':letters_list[int(randi_range(0,len(letters_list)-1))],
				'c':letters_list[int(randi_range(0,len(letters_list)-1))],
				'd':letters_list[int(randi_range(0,len(letters_list)-1))],
				'1':digits_list[int(randi_range(0,len(digits_list)-1))],
				'2':digits_list[int(randi_range(0,len(digits_list)-1))],
				'3':digits_list[int(randi_range(0,len(digits_list)-1))],
				'4':digits_list[int(randi_range(0,len(digits_list)-1))]
				})
			#substance.set_name(substance_name)
			substance.surrounding_area = Rect2(Vector2(substance_starting_point.x-(substance.appearance/2.0) + (substance.appearance * location_x),substance_starting_point.y-(substance.appearance/2.0) + (substance.appearance * location_y)),Vector2(substance.appearance,substance.appearance))
			substance.particle_workings['mass'] = 1.0
			substance.particle_workings['velocity'] = Vector2(0.0,0.0)
			substance.particle_workings['stress'] = [[1.0,1.0],[1.0,1.0]].duplicate(true)
			#substance.particle_workings['stress'] = [1.0,1.0,1.0,1.0]
			substance.particle_workings['volume'] = 1.0
			substance.particle_workings['B'] = substance.B.duplicate(true)
			substance.particle_workings['C'] = substance.C.duplicate(true)
			substance.particle_workings['I'] = substance.I.duplicate(true)
			substance.particle_workings['F'] = substance.F.duplicate(true)
			substance.particle_workings['J'] = substance.J
			
			substance.particle_workings['within_range'] = [substance_particle_name].duplicate(true)
			### 
			substance.particle_mechanics[substance_particle_name] = substance.particle_workings.duplicate(true)
			###
			substance.grid[substance_particle_name] = {'mass':0.0,'velocity':Vector2(0,0),'momentum':Vector2(0.0,0.0)}
			substance.particle_lineation[substance_particle_name] = substance.surrounding_area
			#substance.particle_mechanics[substance_particle_name]['position'] = substance.particle_lineation[substance_particle_name].position
			### 
			location_x =  location_x + 1
			
		
		get_tree().get_root().get_node("Simulation").add_child(substance)
	else:
		print('substances Does not exists')
		pass
