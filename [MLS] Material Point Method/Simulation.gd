extends Node

''' Base Anchor of the Simulation '''
### the simulation happens here...

#parts of the simulation
# core mathematics- [ MLS Material Point Method]
# supplementary math - Matricies
# particles sent thru the simualtion - in the form of Meshes...

var forge : Object
var laboratory : Object
var matrix_math : Object
var constitutive_models : Object
var colliding : Object
var material_method : Object


func _ready():
	if FileAccess.file_exists("res://MeshForge.tscn"):
		forge = load("res://MeshForge.tscn").instantiate()
	else:
		### file does not exists...
		print(' MeshForge.tscn does not exists...')
	if FileAccess.file_exists("res://PhysicalSubstance.tscn"):
		#laboratory = load("res://PhysicalSubstance.tscn").instantiate()
		pass
	else:
		### file does not exists...
		print(' PhysicalSubstance.tscn does not exists...')
	if FileAccess.file_exists("res://MatrixMath.tscn"):
		#laboratory = load("res://PhysicalSubstance.tscn").instantiate()
		matrix_math = load("res://MatrixMath.tscn").instantiate()
	else:
		### file does not exists...
		print(' MatrixMath.tscn does not exists...')
	if FileAccess.file_exists("res://ConstitutiveModels.tscn"):
		#laboratory = load("res://PhysicalSubstance.tscn").instantiate()
		constitutive_models = load("res://ConstitutiveModels.tscn").instantiate()
	else:
		### file does not exists...
		print(' ConstitutiveModels.tscn does not exists...')
	if FileAccess.file_exists("res://FundamentalColliding.tscn"):
		#laboratory = load("res://PhysicalSubstance.tscn").instantiate()
		colliding = load("res://FundamentalColliding.tscn").instantiate()
	else:
		### file does not exists...
		print(' FundamentalColliding.tscn does not exists...')
	if FileAccess.file_exists("res://MaterialMethod.tscn"):
		#laboratory = load("res://PhysicalSubstance.tscn").instantiate()
		material_method = load("res://MaterialMethod.tscn").instantiate()
	else:
		### file does not exists...
		print(' MaterialMethod.tscn does not exists...')
	
	
	add_child(forge)
	add_child(matrix_math)
	add_child(constitutive_models)
	add_child(material_method)
	add_child(colliding)
	
	
	forge.queue_free()
	
	
	set_process(true)
	set_physics_process(false)



func _notification(event):
	###
	#print(event,' event')
	if event == 1003:
		### Window Size Change...
		#print('window change')
		pass



func _process(delta):
	#"""
	for substance in get_tree().get_root().get_node('Simulation/container').get_children():
		$"Material Method".Grid_Reset(substance)
		$"Material Method".Particles_to_Grid(delta,substance)
		$"Material Method".Grid_Update(delta,substance)
		$"Material Method".Collision_with_Wall(substance)
		$"Material Method".Particle_Reset(substance)
		$"Material Method".Grid_to_Particle(delta,substance)
		
		substance.queue_redraw()
	#"""
	pass
