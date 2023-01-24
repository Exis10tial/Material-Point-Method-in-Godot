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
	$"Substance".establish_boundary()
	set_process(true)
	
	#set_physics_process(true)
	

func _process(delta):
	#$"Program".Simulate(delta,$"Substance",$"Substance".grid)
	#$"Substance".establish_boundary()
	#print()
	#print(delta,' delta check')
	#print(check_time,' check_time check')
	#if check_time > 0.0:
		#print('start')
		
	$"Program".Grid_Reset($"Substance",$"Substance".grid)
		
	$"Program".Particles_to_Grid(snapped(delta,.0001),$"Substance",$"Substance".grid)
		
	$"Program".Grid_Update($"Substance",$"Substance".grid)#,test_outside_forces )
		
	$"Program".Collision_with_Wall($"Substance",$"Substance".grid)
	$"Program".Collision_with_Other_Particles($"Substance",$"Substance".grid)
	$"Substance".establish_boundary()
		
	$"Program".Particle_Reset($"Substance")
		
	$"Program".Grid_to_Particle(delta,$"Substance",$"Substance".grid)
		
	$"Substance".queue_redraw()
		
		#check_time = 0
	#else:
		### about .0001+ being added normally...
	#	check_time=+ delta #/ 60.0
		
	
func _physics_process(_delta):
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
