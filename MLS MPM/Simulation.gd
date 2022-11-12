extends Node


#var fund_chemistry : File
var alchemy : Object
#var find_math_book : File
var math_book : Object
#var build_program : File
var program : Object
var collisions : Object
#var contact_with_another : File
#var model_of_substance : File
var models : Object
var check_time : float
var adjust_particles : bool = false
###...
var Particles : Dictionary = {}
var rate : float
var mechanics_processing : bool
var lore_processing : bool


func Establish_Rate():
	### ...
	rate = snapped(1.0/float(len($"Substance".particle_mechanics)),.001)
	#rate = snapped(1.0/1.50,.001)
	#rate = 1.0
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


	Establish_Rate()
	#$"Substance".establish_boundary()
	
	set_process(true)
	set_physics_process(false)
	

func _process(delta):
	$"Program".Simulate(delta,$"Substance",$"Substance".grid)

	$"Substance".queue_redraw()

	$"Substance".establish_boundary()
	#$"Substance".identify_collisions()

func _physics_process(delta):
	#if check_time >= rate:
	#	#Particle Simulation....
	$"Program".Simulate(delta,$"Substance",$"Substance".grid)

	$"Substance".queue_redraw()
	#else:
	#	collect_time(delta)
	#"""
	pass
func collect_time(delta):
	check_time += snapped(delta,.0001)
	return check_time
