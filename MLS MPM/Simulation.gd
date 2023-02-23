extends Node


#var fund_chemistry : File
var alchemy : Object
#var find_math_book : File
var math_book : Object
var math_book_chapter
var matrix_math_page
#var build_program : File
var program : Object
var collisions : Object
#var contact_with_another : File
#var model_of_substance : File
var models : Object
var check_time =0
var adjust_particles : bool = false
###...
#var Particles : Dictionary = {}
var rate : float
var mechanics_processing : bool
var lore_processing : bool
###..
var Substances : Array = []
### Eulerian/Cartesian Grid...
var Grid : Dictionary = {}
var grid_cell_size : int



func Establish_Rate():
	### ...
	#rate = snapped(1.0/float(len($"Substance".particle_mechanics)),.001)
	#rate = snapped(1.0/1.50,.001)
	rate = 1.0
	return rate


func Switch_Protocol():
	
	if lore_processing == true and mechanics_processing == false:
		lore_processing = false
		mechanics_processing = true
		set_process(lore_processing)
		set_physics_process(mechanics_processing)
		
	elif lore_processing == false and mechanics_processing == true:
		lore_processing = true
		mechanics_processing = false
		set_process(lore_processing)
		set_physics_process(mechanics_processing)
		
func Setup_Outline():
	### the outline of the simulation is set..
	### based of the window size...
	#print($"Simulation".get_child_count(),' nodes check ')
	get_tree().get_root().get_node(".").set_position(Vector2((ProjectSettings.get_setting('display/window/size/width')/2.0),0.0))
	get_tree().get_root().get_node(".").set_position(Vector2(ProjectSettings.get_setting('display/window/size/width'),(ProjectSettings.get_setting('display/window/size/height')/2.0)))
	get_tree().get_root().get_node(".").set_position(Vector2((ProjectSettings.get_setting('display/window/size/width')/2.0),ProjectSettings.get_setting('display/window/size/height')))
	get_tree().get_root().get_node(".").set_position(Vector2(0.0,(ProjectSettings.get_setting('display/window/size/height')/2.0)))




func Determine_Eulerian_Grid(x:int,y:int):
	### the size of the substance is parsed into evenly cut squares/cubes = (Length == Width)
	var size 
	var by_subdivision = false
	var find_x_factors = []
	var find_y_factors = []
	
	if x == y and x != 1 and y != 1:# and into_pieces == true :
		### this is if the length/width of the particle are the same...
		### the substance can be in reduced form...
		var n = 0.0
		
		### Check b y division of 2 first...
		while pow(2,n) <= x:
			#by_subdivision = true
			if pow(2,n) == x:
				by_subdivision = true
				#print('by_subdivision check')
				break
			#subdivision_of_2_set_of_numbers.append(pow(2,n))
			#print('find power of 2')
			n +=1
			
		#if subdivision_of_2_set_of_numbers.has(substance_size.y) == true and into_pieces == true:
		#if fmod(substance_size.x,2.0) == 0.0:
		if by_subdivision == true:
			###.division by 2....
			#print(substance_size.x, ' substance_size.x')
			size = int(x / 2.0)
			#by_subdivision = true
		else:
			### the number is not a subdivision of 2...
			
			#Square Root Check
			if fmod(x,sqrt(x)) == 0:# and into_pieces == true:
				size = sqrt(x)
			else:
				### Not a power of 2 or a has a square root...
				
				# find the greatest common factor between the domain sizes::
				for number in range(1,x+1):
					if int(x) % number == 0:
						find_x_factors.append(number)
				for number in range(1,y+1):
					if int(y) % number == 0:
						find_y_factors.append(number)
				
				### using the greatest common factor to acquire a cell size...
				if x == y:
					for number in range(len(find_x_factors)-1,-1,-1):
						if find_x_factors[number] in find_y_factors:
							size = find_x_factors[number] 
							break
				elif x < y:
					for number in range(len(find_x_factors)-1,-1,-1):
						if find_x_factors[number] in find_y_factors:
							size =  find_x_factors[number] 
							break
				elif x > y:
					for number in range(len(find_y_factors)-1,-1,-1):
						if find_y_factors[number] in find_x_factors:
							size = find_y_factors[number]
							break
											
				#by_square_root= true
	else:
		# find the greatest common factor between the domain sizes::
		for number in range(1,x+1):
			if int(x) % number == 0:
				find_x_factors.append(number)
		for number in range(1,y+1):
			if int(y) % number == 0:
				find_y_factors.append(number)
			
		### using the greatest common factor to acquire a cell size...
		if x == y:
			for number in range(len(find_x_factors)-1,-1,-1):
				if find_x_factors[number] in find_y_factors:
					size = find_x_factors[number] 
					break
		elif x < y:
			for number in range(len(find_x_factors)-1,-1,-1):
				if find_x_factors[number] in find_y_factors:
					size =  find_x_factors[number] 
					break
		elif x > y:
			for number in range(len(find_y_factors)-1,-1,-1):
				if find_y_factors[number] in find_x_factors:
					size = find_y_factors[number]
					break
			#"""
		#by_greastest_common_factor = true
		
	return int(size)
	


func Establish_Grid(grid_cell_size:int):
	###...
	
	#print(ProjectSettings.get_setting('display/window/size/width')," ProjectSettings.get_setting('display/window/size/width')")
	#print(ProjectSettings.get_setting('display/window/size/height')," ProjectSettings.get_setting('display/window/size/height')")
	
	#for width in range(0,(ProjectSettings.get_setting('display/window/size/width'))):
	#	for height in range(0,(ProjectSettings.get_setting('display/window/size/height'))):
	#		Grid[Vector2(width,-height)] = {'mass': $"Substance".mass_in_pieces,'velocity':Vector2(0.0,0.0),'momentum':Vector2(0.0,0.0)}
	print(ProjectSettings.get_setting('display/window/size/height') / grid_cell_size,' height')
	print(ProjectSettings.get_setting('display/window/size/width') / grid_cell_size,' width')
	var width = 0
	var height= 0
	while true:
		if height == (ProjectSettings.get_setting('display/window/size/height') / grid_cell_size):
			break
		if width == (ProjectSettings.get_setting('display/window/size/width') / grid_cell_size):
			height = height + 1
			width = 0
		Grid[Vector2(width,-height)] = {'mass': $"Substance".mass_in_pieces,'velocity':Vector2(0.0,0.0),'momentum':Vector2(0.0,0.0)}
		#
		width = width + 1
	
	$"Program".grid_nodes = Grid.duplicate(true)
	#print(len(Grid),' len(Grid)')
	#print(Grid,' Grid')
	#print(Grid[Vector2(0,0)],' Grid[Vector2(0,0)]')
	
	

#func _notification(event):
	#print(event)
	#if event == CanvasItem.NOTIFICATION_DRAW:
	#	print('test')
		#for particle in Particles:
			#Particles[particle].draw_particle()
	
func _on_Simulation_ready():
	#fund_chemistry = File.new()
	if FileAccess.file_exists("res://AlchemyLab.tscn"):
		alchemy = load("res://AlchemyLab.tscn").instantiate()
	else:
		### file does not exists...
		print(' substance.tscn does not exists...')
	#build_program = File.new()
	if FileAccess.file_exists("res://Program.tscn"):
		program = load("res://Program.tscn").instantiate()
	else:
		### file does not exists...
		print(' program.tscn does not exists...')
	#find_math_book = File.new()
	if FileAccess.file_exists("res://Matrix Math.tscn"):
		math_book = load("res://Matrix Math.tscn").instantiate()
	else:
		### file does not exists...
		print(' matrix_math.tscn does not exists...')
	#contact_with_another = File.new()
	if FileAccess.file_exists("res://Particle Interaction.tscn"):
		collisions = load("res://Particle Interaction.tscn").instantiate()
	else:
		### file does not exists...
		print('collisions are not there...')
	#model_of_substance = File.new()
	if FileAccess.file_exists("res://Constitutive Models.tscn"):
		models = load("res://Constitutive Models.tscn").instantiate()
	else:
		### file does not exists...
		print('models are not there...')
	
	add_child(program)
	add_child(math_book)
	add_child(collisions)
	add_child(models)
	add_child(alchemy)

	
	Setup_Outline()
	
	#grid_cell_size = Determine_Eulerian_Grid(ProjectSettings.get_setting('display/window/size/width'),ProjectSettings.get_setting('display/window/size/height'))
	#print(grid_cell_size,' size of grid cells')
	#Establish_Grid(grid_cell_size)
	
	#Establish_Grid()
	$"Substance".establish_boundary()
	set_process(true)
	
	#set_physics_process(false)
	

func _process(delta):
	#$"Program".Simulate(delta,$"Substance",$"Substance".grid)
	#$"Substance".establish_boundary()
	
	#print(delta,' delta check')
	#print(check_time,' check_time check')
	#if check_time > 0.0:
		#print('start')
		
	
	$"Program".Grid_Reset($"Substance")
	$"Program".Particles_to_Grid(snapped(delta,.0001),$"Substance")
	$"Program".Grid_Update(delta,$"Substance")
	$"Program".Collision_with_Wall($"Substance")
	$"Program".Collision_with_Other_Particles($"Substance")
	$"Program".Particle_Reset($"Substance")
	$"Program".Grid_to_Particle(delta,$"Substance")
	$"Substance".establish_boundary()
	$"Substance".queue_redraw()
	
	#$"Program".Grid_Reset($"Substance",$"Substance".grid)
	#$"Program".Particles_to_Grid(snapped(delta,.0001),$"Substance",$"Substance".grid,$"Substance".cell_size)
	#$"Program".Particles_to_Grid(snapped(delta,.0001),$"Substance",$"Substance".grid_nodes,grid_cell_size)
	#$"Program".Grid_Update($"Substance",$"Substance".grid)#,test_outside_forces )
#	$"Program".Grid_Update($"Substance")
	#$"Program".Collision_with_Wall($"Substance",$"Substance".grid)
	#$"Program".Collision_with_Other_Particles($"Substance",$"Substance".grid)
	#$"Substance".establish_boundary()
		
	#$"Program".Particle_Reset($"Substance")
		
	#$"Program".Grid_to_Particle(delta,$"Substance",$"Substance".grid)
	#$"Program".Grid_to_Particle(delta,$"Substance",$"Program".grid_nodes)
		
	#$"Substance".queue_redraw()
		
		#check_time = 0
	#else:
		### about .0001+ being added normally...
	#	check_time=+ delta #/ 60.0
		
	
func physics_process(_delta):
	#if check_time >= rate:
		#Particle Simulation....
	#$"Program".Grid_Reset($"Substance",$"Substance".grid)
	#$"Program".Particles_to_Grid(check_time,$"Substance",$"Substance".grid)
	
	#$"Program".Grid_Update($"Substance",$"Substance".grid)
	#$"Program".Collision_Detection($"Substance",$"Substance".grid)
	#$"Substance".establish_boundary()
	#$"Program".Grid_to_Particle(check_time,$"Substance",$"Substance".grid)
	
	$"Substance".queue_redraw()
#	else:
	#	collect_time(delta)
	

#func build_time(delta):
#	check_time=+ delta
#	if check_time > 1.0:
#		check_time = 0
