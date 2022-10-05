extends Node


var identify_substance : File
var substance : Object
var find_math_book : File
var math_book : Object
var build_program : File
var program : Object
var collisions : Object
var contact_with_another : File
var model_of_substance : File
var models : Object
var check_time : float
var adjust_particles : bool = false
###...


# Called every frame. 'delta' is the elapsed time since the previous frame.
#func _process(delta):
#	pass

func _on_Simulation_ready():
	###...
	
	identify_substance = File.new()
	if identify_substance.file_exists("res://Substance.tscn"):
		substance = load("res://Substance.tscn").instantiate()
	else:
		### file does not exists...
		print(' substance.tscn does not exists...')
	build_program = File.new()
	if build_program.file_exists("res://Program.tscn"):
		program = load("res://Program.tscn").instantiate()
	else:
		### file does not exists...
		print(' program.tscn does not exists...')
	find_math_book = File.new()
	if find_math_book.file_exists("res://Matrix Math.tscn"):
		math_book = load("res://Matrix Math.tscn").instantiate()
	else:
		### file does not exists...
		print(' matrix_math.tscn does not exists...')
	contact_with_another = File.new()
	if contact_with_another.file_exists("res://Particle Interaction.tscn"):
		collisions = load("res://Particle Interaction.tscn").instantiate()
	else:
		### file does not exists...
		print('collisions are not there...')
	model_of_substance = File.new()
	if model_of_substance.file_exists("res://Constitutive Models.tscn"):
		models = load("res://Constitutive Models.tscn").instantiate()
	else:
		### file does not exists...
		print('models are not there...')
	
	add_child(program)
	add_child(math_book)
	add_child(collisions)
	add_child(models)
	add_child(substance)
	
	### initial grid node data....
	$"Program".grid_nodes = $"Substance".grid_nodes.duplicate(true)
	#for particle in $"Substance".get_children():
	#	print(particle.velocity,' before simulation')
	
	set_process(true)
	#set_physics_process(true)

func _process(delta):
	pass
	#"""
	if check_time >= $'Substance'.dx:
		
		
		### if particles ares set to merge...
		#if $"Program".merge_particles == true:
			
		#	get_tree().get_root().get_node("Test Area/Simulation/Particle Interaction").Merge_Particles($"Program".cycle_of_mergin_particles,$"Program".grid_nodes,$"Substance".cell_size)
			### reset....
		#	$"Program".merge_particles = false
		
		#Particle Simulation....
		$"Program".Simulate($"Substance".dx,$"Substance".get_children(),$"Program".grid_nodes)
		
		### used if the particle travels past the window...
		### and the window is set to 'disappear'...
		#'''
		if adjust_particles == true:
			### if particles are to be removed..
			for particle in $"Substance".get_children():
				if particle.to_remove == true:
					particle.to_remove = false
					particle.free()
					$"Program".grid_nodes.erase(particle)
			
			### reset...
			adjust_particles = false
		#'''
	else:
		collect_time(delta)
	#"""

func collect_time(delta):
	check_time += snapped(delta,.01)
	return check_time
