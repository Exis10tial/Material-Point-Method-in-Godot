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
	rate = 0.0
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
	get_tree().get_root().get_node(".").set_position(Vector2((ProjectSettings.get_setting('display/window/size/viewport_width')/2.0),0.0))
	get_tree().get_root().get_node(".").set_position(Vector2(ProjectSettings.get_setting('display/window/size/viewport_width'),(ProjectSettings.get_setting('display/window/size/viewport_height')/2.0)))
	get_tree().get_root().get_node(".").set_position(Vector2((ProjectSettings.get_setting('display/window/size/viewport_width')/2.0),ProjectSettings.get_setting('display/window/size/viewport_height')))
	get_tree().get_root().get_node(".").set_position(Vector2(0.0,(ProjectSettings.get_setting('display/window/size/viewport_height')/2.0)))




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
	
	
	#$"Substance".establish_boundary()
#	print("ready")
	set_process(false)
	set_physics_process(true)
	

func _process(delta):
	
	#print('cycle')
	
	if rate >= 1/60.0:
		$"Program".Grid_Reset($"Substance")
		$"Program".Particles_to_Grid(snapped(rate,.001),$"Substance")
		$"Program".Grid_Update(rate,$"Substance")
		#$"Program".Collision_with_Wall($"Substance")
		#$"Program".Collision_with_Other_Particles($"Substance")
		$"Program".Particle_Reset($"Substance")
		$"Program".Grid_to_Particle(snapped(rate,.001),$"Substance")
		$"Substance".establish_boundary()
		$"Substance".queue_redraw()
		rate = 0.0
	else:
		rate = rate + delta
	
func _physics_process(delta):
	#print()
	#print('turn')
	$"Program".Grid_Reset($"Substance")
	$"Program".Particles_to_Grid(snapped(delta,.001),$"Substance")
	$"Program".Grid_Update(delta,$"Substance")
	#$"Program".Collision_with_Wall($"Substance")
	#$"Program".Collision_with_Other_Particles($"Substance")
	$"Program".Particle_Reset($"Substance")
	$"Program".Grid_to_Particle(snapped(delta,.001),$"Substance")
	#$"Substance".establish_boundary()
	$"Substance".queue_redraw()
#	
#func build_time(delta):
#	check_time=+ delta
#	if check_time > 1.0:
#		check_time = 0
