extends Node2D


### particle state,type of state...
var physical_state : String
var type_of_substance : String
var name_of_substance : String
var constitutive_model : String
var youngs_modulus : float
var poisson_ratio : float
### properties...
var mass : float
var velocity : Vector2
var volume: float
###...
var coefficient_of_restitution : float = 0.0
var coefficient_of_static_friction : float = 0.0
var coefficient_of_kinetic_friction : float = 0.0
### cauchy stress
var stress : Array = [[1.0,1.0],[1.0,1.0]]


###.
var B : Array = [[0.0,0.0],[0.0,0.0]]
### Affine
var C : Array = [[0.0,0.0],[0.0,0.0]]
### Deformation Identity
var I : Array = [[1.0,0.0],[0.0,1.0]]
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
### particle relation with other particles...
var relation_to_domain : Dictionary = {}
var domain_relation_to_particle : Dictionary = {}
### used for finding other particles around it...
var within_range : Array = []
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

func _on_2D_Particle_ready():
	###...
	#$"shape".set_position($"shape".rect_size/2.0)
	pass

func adjust_size(new_size):
	### adjust the size of the particle...
	$"shape".set_size(Vector2(new_size,new_size))
	$"shape".set_position(Vector2(new_size/2.0,new_size/2.0))
