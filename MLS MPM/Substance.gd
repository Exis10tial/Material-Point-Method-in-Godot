extends Node2D

### the particles of a substance...
var particle_workings : Dictionary = {'position':null,'mass':null,'velocity':null,'volume':null,'stress':null,'B':null,'C':null,'I':null,'F':null,'J':null,'eulerian':{},'contact_with_wall':false,'within_range':[],'relation_to_domain':{},'domain_relation_to_substance':{},'U':null,'V':null,'Sigma':null,'contact wall':false,'contact particle':false}
var particle_lineation : Dictionary = {}
var particle_mechanics : Dictionary = {}
var grid_ : Dictionary = {}
#var grid_scope : Array = [{}]
var substance_limit : int = 0
var mass_in_pieces : float
var volume_in_pieces : float
### mechanics of the particle for the servers...
#var particle_body : Dictionary = {}
#var touch_modality : PhysicsMaterial
#var domain : Dictionary
#var effigy : Dictionary
#var physical_structure : Dictionary
#var form : Dictionary
#var building_effigy_material
### particle state,type of state...
var physical_state : String
var type_of_substance : String
var name_of_substance : String
var constitutive_model : String
var youngs_modulus : float
var poisson_ratio : float
### properties...
var mass : float
var initial_velocity : Vector2
var maintain_velocity : Vector2
var volume: float
var density : float = mass / volume
###...
var cell_size : int = 1
var appearance : Vector2
var one_substance_x :int
var one_substance_y : int
###.
var coefficient_of_restitution : float = 0.0
var coefficient_of_static_friction : float = 0.0
var coefficient_of_kinetic_friction : float = 0.0
### used for other mechanical parts...
var identify_number
### scope of how the particles will search the grid...
#center

#var scope = [Vector2(-1,1),Vector2(1,1),Vector2(0,0),Vector2(-1,-1),Vector2(1,-1)]
#var scope = [Vector2(-1,1),Vector2(0,1),Vector2(1,1),Vector2(-1,0),Vector2(1,0),Vector2(-1,-1),Vector2(0,-1),Vector2(1,-1)]


### cauchy stress
var stress : Array = [1.0,1.0,1.0,1.0]
###.
var B : Array = [0.0,0.0,0.0,0.0]
### Affine
var C : Array = [0.0,0.0,0.0,0.0]
### Deformation Identity
var I : Array = [1.0,0.0,0.0,1.0]
#const I : = Transform2D.IDENTITY
### Polar SVD: mainly used for drucker_prager model...
var Sigma : Array
var U : Array
var V : Array
var yield_surface : float = 0.0
### Deformation Gradient...
var F : Array = I.duplicate(true)
var F_elastic : Array
var F_plastic : Array
### determinant of F... 
var J : float
var J_elastic : float
var J_plastic: float


### the surrounding area around the particle/node
var surrounding_area : Rect2
#var surrounding_area : Transform2D
### particle is to be removed from the sim...
var to_remove = false
### parameters used to determine to merge particels...
var magnitude : float 
var direction : float
var covers : float = 0.0
var covered_by_at : Dictionary = {}
### the partilce will merge
var merging : bool = false
var merge_with : Object
### if the particle is a merged particle...
var made_of : Array = []


func _on_substance_ready():
	pass # Replace with function body.

func _on_substance_draw():
	###
	for particle in particle_lineation.keys():
		draw_rect(particle_lineation[particle],Color(1.0,1.0,1.0),true)
		
		#draw_circle(particle_lineation[particle].origin,appearance,Color(1.0,1.0,1.0))
		
		pass

func establish_boundary():
	var copy_lineation 
	var total
	var contact = 0
	var results
	var rotations 
	var switch
	var template = []
	
	copy_lineation = particle_lineation.keys().duplicate(true)
	rotations = 0
	identify_number = 0
	#print('establish_boundary check')
	while true:
		if identify_number >= len(copy_lineation):
		#if identify_number <= len(copy_lineation):
		#	results = []
			identify_number = 0
			copy_lineation.pop_at(0)
			#print('algorithm check')
			if len(copy_lineation) == 0:
				#print(len(copy_lineation),' copy len')
				break

		#contact = snapped(particle_lineation[particle_lineation.keys()[len(particle_lineation)-len(copy_lineation)]].get_center().distance_squared_to(particle_lineation[copy_lineation[identify_number]].get_center()),.0001)
		var relation_between_particles = Vector2(snapped(particle_lineation[copy_lineation[identify_number]].get_center().x - particle_lineation[particle_lineation.keys()[len(particle_lineation)-len(copy_lineation)]].get_center().x,.01),snapped(particle_lineation[copy_lineation[identify_number]].get_center().y - particle_lineation[particle_lineation.keys()[len(particle_lineation)-len(copy_lineation)]].get_center().y,.01))

		#"""
		#if contact <= 1:
		if abs(relation_between_particles/appearance).x >= 0.0 and abs(relation_between_particles/appearance).x < 1.0 and abs(relation_between_particles/appearance).y >= 0.0 and abs(relation_between_particles/appearance).y < 1.0:
			###
			#if contact == 0:
			if relation_between_particles.x == 0.0 and relation_between_particles.y == 0:
				### this is itself...
				if particle_mechanics[particle_lineation.keys()[len(particle_lineation)-len(copy_lineation)]]['within_range'].has(copy_lineation[identify_number]) == false:
					#particle_mechanics[particle_lineation.keys()[len(particle_lineation)-len(copy_lineation)]]['within_range'].append(particle_lineation.keys()[len(particle_lineation)-len(copy_lineation)])
					#print(particle_mechanics[particle_lineation.keys()[len(particle_lineation)-len(copy_lineation)]]['within_range'],' within check a')
					pass
				else:
					### already recognized...
					pass
			else:
					
				#print(particle_lineation.keys()[len(particle_lineation)-len(copy_lineation)],' in contact with ',copy_lineation[identify_number])
				if particle_mechanics[particle_lineation.keys()[len(particle_lineation)-len(copy_lineation)]]['within_range'].has(copy_lineation[identify_number]) == false:
					particle_mechanics[particle_lineation.keys()[len(particle_lineation)-len(copy_lineation)]]['within_range'].append(copy_lineation[identify_number])
					#print(particle_mechanics[particle_lineation.keys()[len(particle_lineation)-len(copy_lineation)]]['within_range'],' within check b')
				else:
					### already recognized...
					pass
		else:
			### The particles is out of range..
			if particle_mechanics[particle_lineation.keys()[len(particle_lineation)-len(copy_lineation)]]['within_range'].has(copy_lineation[identify_number]) == true:
					particle_mechanics[particle_lineation.keys()[len(particle_lineation)-len(copy_lineation)]]['within_range'].erase(copy_lineation[identify_number])
			
					#print(particle_mechanics[particle_lineation.keys()[len(particle_lineation)-len(copy_lineation)]]['within_range'],' within check c')
			else:
				### ...
				pass
			
		#print(particle_mechanics[particle_lineation.keys()[len(particle_lineation)-len(copy_lineation)]]['within_range'],' within check')
		#"""
		
		#identify_number = identify_number + 1
		identify_number = wrapi(identify_number+1,0,len(copy_lineation)+1)
		
	return
