extends Node2D


#var skeleton : Skeleton2D
### setup thru polygon...
var effigy : Polygon2D
#var effigy_visual : ImageTexture
#var effigy_particles : PackedVector2Array
#var effigy_particle_color : PackedColorArray
#var effify_sensor_uvs : PackedVector2Array
### setup thru Mesh...
#var effigy : Mesh
#var effigy_data = MeshDataTool.new()
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
var inner_workings : Dictionary = {'mass':null,'velocity':null,'initial velocity':Vector2(0,0), 'volume':null,'stress':null,'B':null,'C':null,'I':null,'F':null,'J':null,'U':null,'V':null,'sigma':null,'mu':null,'lambda':null, 'eulerian':{},'euler data':null}
var mechanics : Dictionary = {}
### allow multiple polygons...
var centroid_rank : int
var effigy_basket : Array
var associate_polygon_to_particle : Dictionary


func _draw():
	
	#draw_polygon(effigy.get_polygon(),effigy.get_vertex_colors(),effigy.get_uv())
	for wax in effigy_basket:
		draw_polygon(wax.get_polygon(),wax.get_vertex_colors(),wax.get_uv())
