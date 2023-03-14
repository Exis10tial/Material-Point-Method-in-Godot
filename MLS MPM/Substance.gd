extends Node2D

### the particles of a substance...
var particle_workings : Dictionary = {'position':null,'mass':null,'velocity':null,'volume':null,'stress':null,'B':null,'C':null,'I':null,'F':null,'J':null,'eulerian':{},'body':null,'effigy':null,'gradient forces':{},'contact_with_wall':false,'within_range':[],'relation_to_domain':{},'domain_relation_to_substance':{},'U':null,'V':null,'Sigma':null,'contact wall':false,'contact particle':false}
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
#var I : Array = [1.0,0.0,0.0,1.0]
const I : = Transform2D.IDENTITY
### Polar SVD: mainly used for drucker_prager model...
var Sigma : Array
var U : Array
var V : Array
var yield_surface : float = 0.0
### Deformation Gradient...
#var F : Array = I.duplicate(true)
var F : Array 
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
		#draw_rect(particle_lineation[particle],Color(1.0,1.0,1.0),true)
		draw_mesh(particle_mechanics[particle]['body'],particle_mechanics[particle]['effigy'],particle_lineation[particle],Color(1,1,1,1))
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
	
		# copy of particles...
		#particle_lineation[copy_lineation[identify_number]]
		# cycle thru every other particle...
		#particle_lineation.keys()[len(particle_lineation)-len(copy_lineation)]
		
		if rad_to_deg(particle_lineation[copy_lineation[identify_number]].get_rotation()) > 0:
			### if the particle has the slightest rotation...
			
			var topleft_rotated = Vector2(particle_lineation[copy_lineation[identify_number]].origin.x - ( (appearance.x/2.0) * rad_to_deg(cos(particle_lineation[copy_lineation[identify_number]].get_rotation())) - (appearance.y/2.0) * rad_to_deg(sin(particle_lineation[copy_lineation[identify_number]].get_rotation())) ),
			particle_lineation[copy_lineation[identify_number]].origin.y - ( (appearance.x/2.0) * rad_to_deg(sin(particle_lineation[copy_lineation[identify_number]].get_rotation())) + (appearance.y/2.0) * rad_to_deg(cos(particle_lineation[copy_lineation[identify_number]].get_rotation())) ) )
			
			var topright_rotated = Vector2(material.particle_lineation[copy_lineation[identify_number]].origin.x + ( (appearance.x/2.0) * rad_to_deg(cos(particle_lineation[copy_lineation[identify_number]].get_rotation())) - (appearance.y/2.0) * rad_to_deg(sin(particle_lineation[copy_lineation[identify_number]].get_rotation())) ),
			particle_lineation[copy_lineation[identify_number]].origin.y + ( (appearance.x/2.0) * rad_to_deg(sin(particle_lineation[copy_lineation[identify_number]].get_rotation())) + (appearance.y/2.0) * rad_to_deg(cos(particle_lineation[copy_lineation[identify_number]].get_rotation())) ) )
			
			var bottomright_rotated = Vector2(particle_lineation[copy_lineation[identify_number]].origin.x + ( (appearance.x/2.0) * rad_to_deg(cos(particle_lineation[copy_lineation[identify_number]].get_rotation())) + (appearance.y/2.0) * rad_to_deg(sin(particle_lineation[copy_lineation[identify_number]].get_rotation())) ),
			particle_lineation[copy_lineation[identify_number]].origin.y + ( (appearance.x/2.0) * rad_to_deg(sin(particle_lineation[copy_lineation[identify_number]].get_rotation())) - (appearance.y/2.0) * rad_to_deg(cos(particle_lineation[copy_lineation[identify_number]].get_rotation())) ) )
			
			var bottomleft_rotated = Vector2(particle_lineation[copy_lineation[identify_number]].origin.x - ( (appearance.x/2.0) * rad_to_deg(cos(particle_lineation[copy_lineation[identify_number]].get_rotation())) + (appearance.y/2.0) * rad_to_deg(sin(particle_lineation[copy_lineation[identify_number]].get_rotation())) ),
			particle_lineation[copy_lineation[identify_number]].origin.y - ( (appearance.x/2.0) * rad_to_deg(sin(particle_lineation[copy_lineation[identify_number]].get_rotation())) - (appearance.y/2.0) * rad_to_deg(cos(particle_lineation[copy_lineation[identify_number]].get_rotation())) ) )
			
			
			
			
			
			
		elif particle_lineation[copy_lineation[identify_number]].get_rotation() == 0:
			### if the particle has no rotation...
			
			var topleft = Vector2(snapped(particle_lineation[copy_lineation[identify_number]].origin.x,.01),snapped(particle_lineation[copy_lineation[identify_number]].origin.y,.01)) + Vector2(-(appearance.x/2.0),-(appearance.y/2.0))
			var topright = Vector2(snapped(particle_lineation[copy_lineation[identify_number]].origin.x,.01),snapped(particle_lineation[copy_lineation[identify_number]].origin.y,.01)) + Vector2((appearance.x/2.0),-(appearance.y/2.0))
			var bottomright = Vector2(snapped(particle_lineation[copy_lineation[identify_number]].origin.x,.01),snapped(particle_lineation[copy_lineation[identify_number]].origin.y,.01)) + Vector2((appearance.x/2.0),(appearance.y/2.0))
			var bottomleft = Vector2(snapped(particle_lineation[copy_lineation[identify_number]].origin.x,.01),snapped(particle_lineation[copy_lineation[identify_number]].origin.y,.01)) + Vector2(-(appearance.x/2.0),(appearance.y/2.0))
			






		#identify_number = identify_number + 1
		identify_number = wrapi(identify_number+1,0,len(copy_lineation)+1)
		
	return
