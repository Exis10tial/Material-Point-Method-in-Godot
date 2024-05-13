extends Node2D


### mechanics of particles with geometrics of 1......
var box_shaped : Array = []
var entity_container : Dictionary
var entity_color : Dictionary
### Physical properties of the entire substance...
var coefficient_of_restitution : float
var coefficient_of_static_friction : float
var coefficient_of_kinetic_friction : float
var physical_state : String
var type_of_substance : String
var constitutive_model : String
var poisson_ratio : float
var youngs_modulus : float
var mass : float
var volume : float
### Material Method Properties Mechanics...
var inner_workings : Dictionary = {'mass':null,'velocity':null,'initial velocity':Vector2(0,0), 'volume':null,'stress':null,'B':null,'C':null,'I':null,'F':null,'J':null,'U':null,'V':null,'sigma':null,'mu':null,'lambda':null, 'eulerian':{},'euler data':null,'particle to grid':null,'grid to particle':null}
var mechanics : Dictionary = {}
### Material Deformation Components...
var I : Array = [1,0,0,1]
var F : Array = I.duplicate(true)
var J : float = 0
var U : Array
var sigma : Array
var V : Array
### partlice relation to grid mechanics...
#var relation_of_particle_to_grid : Dictionary = {}
#var relation_of_grid_to_particle : Dictionary = {}
var weight_interpolation : Dictionary = {}
### allow multiple polygons...
#var centroid_rank : int
#var effigy_basket : Array
#var associate_polygon_to_particle : Dictionary
### Mechanics of Self-Collisions
#var default_distance = 0
#var maximum_distance = default_distance * 2
#var minimum_distance = default_distance * (1.0/256.0)




#"""
func _draw():
	var identified_particle = 0
	
	### drawing a a particle...
	while true:
		if identified_particle >= len(entity_container):
			break
		
		draw_rect(entity_container[identified_particle],entity_color[identified_particle])
		
		### loop
		identified_particle = wrapi(identified_particle+1,0,len(entity_container)+1)
		
