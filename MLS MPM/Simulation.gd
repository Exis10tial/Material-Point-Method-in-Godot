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

# Called every frame. 'delta' is the elapsed time since the previous frame.
#func _process(delta):
#	pass
func Establish_Rate():
	### ...
	#dx = snapped(1.0/float(particle_limit)*1.50,.001)
	#dx = snapped(1.0/float(particle_limit)*5.00,.001)
	rate = snapped(1.0/float(len($"Substance".particle_mechanics)),.001)
	#dx = snapped(1.0/1.50,.001)
	#dx = 1.50
	#dx = 1.0
	return rate

#func _notification(event):
	#print(event)
	#if event == CanvasItem.NOTIFICATION_DRAW:
	#	print('test')
		#for particle in Particles:
			#Particles[particle].draw_particle()
	
func _on_Simulation_ready():
	###...
	
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
	
	### initial grid node data....
	#$"Program".grid_nodes = $"Substance".grid_nodes.duplicate(true)
	#for particle in $"Substance".get_children():
	#	print(particle.velocity,' before simulation')
	Establish_Rate()
	
	set_process(true)
	#set_physics_process(true)
	
func _process(delta):
	#"""
	if check_time >= rate:

		### if particles ares set to merge...
		#if $"Program".merge_particles == true:
			
		#	get_tree().get_root().get_node("Test Area/Simulation/Particle Interaction").Merge_Particles($"Program".cycle_of_mergin_particles,$"Program".grid_nodes,$"Substance".cell_size)
			### reset....
		#	$"Program".merge_particles = false
		
		#Particle Simulation....
		#$"Program".Simulate($"Substance".dx,$"Substance".get_children(),$"Program".grid_nodes)
		$"Program".Simulate(rate,$"Substance",$"Substance".grid)
		### used if the particle travels past the window...
		### and the window is set to 'disappear'...
		#'''
		if adjust_particles == true:
			### if particles are to be removed..
			for particle in Particles.keys():
				if particle.to_remove == true:
					particle.to_remove = false
					particle.free()
					$"Substance".grid.erase(particle)
				
			### reset...
			adjust_particles = false
		#'''
		
		for particle in $"Substance".particle_lineation.keys():
			print($"Substance".particle_lineation[particle],' end of cycle')
		
		### the particle/s position is updated...
		$"Substance".queue_redraw()
	else:
		collect_time(delta)
	#"""

func collect_time(delta):
	check_time += snapped(delta,.01)
	return check_time
