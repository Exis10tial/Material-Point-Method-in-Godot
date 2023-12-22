extends Node


### determine the type of collision...
var collision_restitution : float = 0.0
var collision_static_friction : float = 0.0
var collision_kinetic_friction : float = 0.0
var wall_center : Vector2
var impact_center : Vector2
var decimal_reminder : float
### particle merging...
var length_covered_area : float 
var width_covered_area : float
var covered_area : float
var percentage_covered : float
### used for elastic collisions >=1.0
var normal_vector : Vector2
var unit_vector : Vector2
var unit_tangent : Vector2
var dotted_unit_mote_velocity : float
var dotted_tangent_mote_velocity : float
var dotted_unit_wall_velocity: float
var dotted_tangent_wall_velocity: float
var final_tangential_mote_velocity : float
var final_tangential_wall_velocity : float
var final_normal_mote_velocity : float
var final_normal_wall_velocity : float
### used for impartial collisions < 1.0 and > 0.0 ...
var line_of_impact_x : float
var line_of_impact_y : float
var line_of_impact : Vector2
var x_component : float
var y_component : float
var updated_artifact_x_component : float
var updated_other_artifact_x_component : float
var static_friction : float
var kinetic_friction : float
var other_static_friction : float
var other_kinetic_friction : float
var incoming_angle 
var outgoing_angle
var y_leg_of_particle_velocity : float
var x_leg_of_particle_velocity : float
### used for perfect inelastic collisions > 0.0 ...
var final_velocity : Vector2
### particles colliding with other particle at elastic collisions...
var dotted_unit_artifact_velocity# : float
var dotted_tangent_artifact_velocity : float
var dotted_unit_other_artifact_velocity: float
var dotted_tangent_other_artifact_velocity: float
var final_tangential_artifact_velocity : float
var final_tangential_other_artifact_velocity : float
var final_normal_artifact_velocity : float
var final_normal_other_artifact_velocity : float
### particles colliding with other particle at impartial collisions...
var incoming_angle_of_artifact
var outgoing_angle_of_artifact
var y_leg_of_artifact_velocity
var incoming_angle_of_other_artifact
var outgoing_angle_of_other_artifact
var x_leg_of_updated_other_artifact_velocity
var x_leg_of_updated_artifact_velocity
### used to create a new-particle when merging particles....
var empirical_particle : Object
var empirical_name : String = ''
var letters_list : Array = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
var digits_list : Array = ['0','1','2','3','4','5','6','7','8','9']
var sum_of_parts_mass : float
var sum_of_parts_kernel : Vector2
var sum_of_parts_velocity : Vector2
var area_multiplier : float
var kernel_distance : Vector2
var kernel_x : float
var kernel_y : float
var flipped_kernel_distance : Vector2
var flipped_kernel_x : float
var flipped_kernel_y : float

###

###
func _on_handle_collisions_ready():
	pass # Replace with function body.


func Collision_with_Walls(breach,mote,particle_boundary,baluster,structure,_cell_size):
	###...
	if breach == 'top':
		
		collision_restitution = baluster['window outline']['top']['coefficient of restitution'] * mote.coefficient_of_restitution
		collision_static_friction = baluster['window outline']['top']['coefficient of static friction'] * mote.coefficient_of_static_friction
		collision_kinetic_friction = baluster['window outline']['top']['coefficient of kinetic friction'] * mote.coefficient_of_kinetic_friction
		
		#if structure['velocity'].y > 0:
		
		#baluster['window outline']['top']['velocity'] = Vector2(0,0)
		#baluster['window outline']['top']['velocity'] = Vector2(0,-structure['velocity'].y)
		
		### for when collision_restitution == 1.0
		baluster['window outline']['top']['velocity'] = Vector2(0,-structure['velocity'].y) * 1.0
		
		wall_center = Vector2(particle_boundary.origin.x,0.0)
		
		if collision_restitution >= 1.0 :
			### the collision is perfect elastic...
			
			#normal_vector = Vector2(wall_center - particle_boundary.get_center())
			normal_vector = Vector2(wall_center - particle_boundary.origin)
			unit_vector = normal_vector / snapped(sqrt((snapped(pow(normal_vector.x,2.0),.01) + snapped(pow(normal_vector.y,2.0),.01))),.01)
			unit_tangent = Vector2(-unit_vector.y,unit_vector.x)
			
			#dotted_unit_mote_velocity = unit_vector.dot(structure['velocity'])
			dotted_unit_mote_velocity = snapped(unit_vector.dot(structure['velocity']),.001)
			#dotted_tangent_mote_velocity = unit_tangent.dot(structure['velocity'])
			dotted_tangent_mote_velocity = snapped(unit_tangent.dot(structure['velocity']),.001)
			dotted_unit_wall_velocity = unit_vector.dot(baluster['window outline']['top']['velocity'])
			dotted_tangent_wall_velocity = unit_tangent.dot(baluster['window outline']['top']['velocity'])
			
			final_tangential_mote_velocity = dotted_tangent_mote_velocity
			final_tangential_wall_velocity = dotted_tangent_wall_velocity
			final_normal_mote_velocity = (dotted_unit_mote_velocity * (structure['mass'] - baluster['window outline']['top']['mass']) + 2.0 * baluster['window outline']['top']['mass'] * dotted_unit_wall_velocity ) / (structure['mass'] + baluster['window outline']['top']['mass'])
			final_normal_wall_velocity = (dotted_unit_wall_velocity * (baluster['window outline']['top']['mass'] - structure['mass']) + 2.0 * baluster['window outline']['top']['mass'] * dotted_unit_mote_velocity ) / (structure['mass'] + baluster['window outline']['top']['mass'])
			
			structure['velocity'] = (final_normal_mote_velocity * unit_vector) + (final_tangential_mote_velocity * unit_tangent)
			#final_tangential_mote_velocity * unit_tangent
			#(final_normal_wall_velocity * unit_vector) + (final_tangential_wall_velocity * unit_tangent)
			
			#print(structure['velocity']," structure['velocity']")
		
		elif collision_restitution < 1.0 and collision_restitution > 0.0:
			
			
			x_component = structure['velocity'].x
			y_component = -baluster['window outline']['top']['coefficient of restitution'] *  structure['velocity'].y 
			
			structure['velocity'] = Vector2(x_component,y_component)
			
		elif collision_restitution == 0.0:
			#perfect inelastic collisions : [m1 / (m1 + m2)] * vi = vf
			
			x_component = structure['mass'] / (structure['mass'] + baluster['window outline']['top']['mass']) * structure['velocity'].x
			y_component = structure['mass'] / (structure['mass'] + baluster['window outline']['top']['mass']) * structure['velocity'].y
			
			structure['velocity'] = Vector2(x_component,y_component)
			
	if breach == 'right':
		
		collision_restitution = baluster['window outline']['right']['coefficient of restitution'] * mote.coefficient_of_restitution
		collision_static_friction = baluster['window outline']['right']['coefficient of static friction'] * mote.coefficient_of_static_friction
		collision_kinetic_friction = baluster['window outline']['right']['coefficient of kinetic friction'] * mote.coefficient_of_kinetic_friction
		
		#baluster['window outline']['right']['velocity'] = Vector2(0,0.0)
		#baluster['window outline']['right']['velocity'] = Vector2(-structure['velocity'].x,0)
		baluster['window outline']['right']['velocity'] = Vector2(-structure['velocity'].x,0) * 1.00
		
		wall_center = Vector2(baluster['window outline']['right']['outline'],particle_boundary.origin.y)
		
		if collision_restitution >= 1.0 :
			### the collision is perfect elastic...
			
			#normal_vector = Vector2(wall_center - particle_boundary.get_center())
			normal_vector = Vector2(wall_center - particle_boundary.origin)
			
			unit_vector = normal_vector / snapped(sqrt((snapped(pow(normal_vector.x,2.0),.01) + snapped(pow(normal_vector.y,2.0),.01))),.01)
			unit_tangent = Vector2(-unit_vector.y,unit_vector.x)
			
			#dotted_unit_mote_velocity = unit_vector.dot(structure['velocity'])
			dotted_unit_mote_velocity = snapped(unit_vector.dot(structure['velocity']),.001)
			#dotted_tangent_mote_velocity = unit_tangent.dot(structure['velocity'])
			dotted_tangent_mote_velocity = snapped(unit_tangent.dot(structure['velocity']),.001)
			dotted_unit_wall_velocity = unit_vector.dot(baluster['window outline']['right']['velocity'])
			dotted_tangent_wall_velocity = unit_tangent.dot(baluster['window outline']['right']['velocity'])
			
			final_tangential_mote_velocity = dotted_tangent_mote_velocity
			final_tangential_wall_velocity = dotted_tangent_wall_velocity
			final_normal_mote_velocity = (dotted_unit_mote_velocity * (structure['mass'] - baluster['window outline']['right']['mass']) + 2.0 * baluster['window outline']['right']['mass'] * dotted_unit_wall_velocity ) / (structure['mass'] + baluster['window outline']['right']['mass'])
			final_normal_wall_velocity = (dotted_unit_wall_velocity * (baluster['window outline']['right']['mass'] - structure['mass']) + 2.0 * baluster['window outline']['right']['mass'] * dotted_unit_mote_velocity ) / (structure['mass'] + baluster['window outline']['right']['mass'])
			
			structure['velocity'] = (final_normal_mote_velocity * unit_vector) + (final_tangential_mote_velocity * unit_tangent)
			
			#final_tangential_mote_velocity * unit_tangent
			#(final_normal_wall_velocity * unit_vector) + (final_tangential_wall_velocity * unit_tangent)
			
			#print(structure['velocity']," structure['velocity']")
		
		elif collision_restitution < 1.0 and collision_restitution > 0.0:
			
			"""
			incoming_angle = rad_to_deg(asin(structure['velocity'].y/structure['velocity'].length()))
			outgoing_angle = snapped(rad_to_deg(atan(baluster['window outline']['right']['coefficient of restitution'] * tan(deg_to_rad(incoming_angle)))),.1)
			"""
			x_component =  -baluster['window outline']['right']['coefficient of restitution'] * structure['velocity'].x
			y_component = structure['velocity'].y
			structure['velocity'] = Vector2(x_component,y_component)
			
		elif collision_restitution == 0.0:
			#perfect inelastic collisions : [m1 / (m1 + m2)] * vi = vf
			
			x_component = structure['mass'] / (structure['mass'] + baluster['window outline']['right']['mass']) * structure['velocity'].x
			y_component = structure['mass'] / (structure['mass'] + baluster['window outline']['right']['mass']) * structure['velocity'].y
			
			structure['velocity'] = Vector2(x_component,y_component)
			
	if breach == 'bottom':
		
		collision_restitution = baluster['window outline']['bottom']['coefficient of restitution'] * mote.coefficient_of_restitution
		collision_static_friction = baluster['window outline']['bottom']['coefficient of static friction'] * mote.coefficient_of_static_friction
		collision_kinetic_friction = baluster['window outline']['bottom']['coefficient of kinetic friction'] * mote.coefficient_of_kinetic_friction
		
		#baluster['window outline']['bottom']['velocity'] = Vector2(0,0)
		#baluster['window outline']['bottom']['velocity'] = Vector2(0,-structure['velocity'].y)
		
		baluster['window outline']['bottom']['velocity'] = Vector2(0,-structure['velocity'].y) * 1.00
		
		#print(collision_restitution,' check collision_restitution')
		
		wall_center = Vector2(particle_boundary.origin.x,baluster['window outline']['bottom']['outline'])
		#wall_center = Vector2(baluster['window outline']['right']['outline']/2,baluster['window outline']['bottom']['outline'])
		
		if collision_restitution >= 1.0 :
			### the collision is perfect elastic...
			#print('bounce back')
			#normal_vector = Vector2(wall_center - particle_boundary.get_center())
			normal_vector = Vector2(wall_center - particle_boundary.origin)
			#print(unit_vector,' unit_vector check')
			unit_vector = normal_vector / snapped(sqrt((snapped(pow(normal_vector.x,2.0),.01) + snapped(pow(normal_vector.y,2.0),.01))),.01)
			unit_tangent = Vector2(-unit_vector.y,unit_vector.x)
			
			#dotted_unit_mote_velocity = unit_vector.dot(structure['velocity'])
			dotted_unit_mote_velocity = snapped(unit_vector.dot(structure['velocity']),.001)
			#dotted_tangent_mote_velocity = unit_tangent.dot(structure['velocity'])
			dotted_tangent_mote_velocity = snapped(unit_tangent.dot(structure['velocity']),.001)
			dotted_unit_wall_velocity = unit_vector.dot(baluster['window outline']['bottom']['velocity'])
			dotted_tangent_wall_velocity = unit_tangent.dot(baluster['window outline']['bottom']['velocity'])
			
			final_tangential_mote_velocity = dotted_tangent_mote_velocity
			final_tangential_wall_velocity = dotted_tangent_wall_velocity
			
			final_normal_mote_velocity = (dotted_unit_mote_velocity * (structure['mass'] - baluster['window outline']['bottom']['mass']) + 2.0 * baluster['window outline']['bottom']['mass'] * dotted_unit_wall_velocity ) / (structure['mass'] + baluster['window outline']['bottom']['mass'])
			final_normal_wall_velocity = (dotted_unit_wall_velocity * (baluster['window outline']['bottom']['mass'] - structure['mass']) + 2.0 * baluster['window outline']['bottom']['mass'] * dotted_unit_mote_velocity ) / (structure['mass'] + baluster['window outline']['bottom']['mass'])
			
			structure['velocity'] = (final_normal_mote_velocity * unit_vector) + (final_tangential_mote_velocity * unit_tangent)
			
		elif collision_restitution < 1.0 and collision_restitution > 0.0:
			#line_of_impact = wall_center
			"""
			print(' ')
			print(structure['velocity'],' velocity')
			#print(wall_center,' wall_center')
			
			print('base ',Vector2(structure['velocity'].x,0))
			print('altitude ',Vector2(0,structure['velocity'].y))
			print('hypotenuse ',structure['velocity'].length() )
			
			print( rad_to_deg(asin(structure['velocity'].y/structure['velocity'].length())),' sin-1')
			print( rad_to_deg(acos(structure['velocity'].x/structure['velocity'].length())),' cos-1')
			print( rad_to_deg(atan(structure['velocity'].y/structure['velocity'].x)),' tan-1')
			
			incoming_angle = rad_to_deg(acos(structure['velocity'].x/structure['velocity'].length()))
			outgoing_angle = snapped(rad_to_deg(atan(baluster['window outline']['bottom']['coefficient of restitution'] * tan(deg_to_rad(incoming_angle)))),.1) 
			#incoming_angle = abs(snapped(rad_to_deg(line_of_impact.angle_to_point(particle_boundary.origin)),1))
			#outgoing_angle = snapped(rad_to_deg(atan(baluster['window outline']['bottom']['coefficient of restitution'] * tan(deg_to_rad(incoming_angle)))),.1)
			#print(incoming_angle,' incoming_angle')
			"""
			
			x_component = structure['velocity'].x
			y_component =  -baluster['window outline']['bottom']['coefficient of restitution'] * structure['velocity'].y
			
			structure['velocity'] = Vector2(x_component,y_component)
			
		elif collision_restitution == 0.0:
			#perfect inelastic collisions : [m1 / (m1 + m2)] * vi = vf
			
			x_component = structure['mass'] / (structure['mass'] + baluster['window outline']['bottom']['mass']) * structure['velocity'].x
			y_component = structure['mass'] / (structure['mass'] + baluster['window outline']['bottom']['mass']) * structure['velocity'].y
			
			structure['velocity'] = Vector2(x_component,y_component)
			#print(structure['velocity']," structure['velocity']")
	
	
	if breach == 'left':
		
		collision_restitution = baluster['window outline']['left']['coefficient of restitution'] * mote.coefficient_of_restitution
		collision_static_friction = baluster['window outline']['left']['coefficient of static friction'] * mote.coefficient_of_static_friction
		collision_kinetic_friction = baluster['window outline']['left']['coefficient of kinetic friction'] * mote.coefficient_of_kinetic_friction
		
		#baluster['window outline']['left']['velocity'] = Vector2(0,0)
		#baluster['window outline']['left']['velocity'] = Vector2(-structure['velocity'].x,0)
		baluster['window outline']['left']['velocity'] = Vector2(-structure['velocity'].x,0) * 1.00
		
		
		wall_center = Vector2(0.0,particle_boundary.origin.y)
		
		if collision_restitution >= 1.0 :
			### the collision is perfect elastic...
			
			#normal_vector = Vector2(wall_center - particle_boundary.get_center())
			normal_vector = Vector2(wall_center - particle_boundary.origin)
			
			unit_vector = normal_vector / snapped(sqrt((snapped(pow(normal_vector.x,2.0),.01) + snapped(pow(normal_vector.y,2.0),.01))),.01)
			unit_tangent = Vector2(-unit_vector.y,unit_vector.x)
			
			#dotted_unit_mote_velocity = unit_vector.dot(structure['velocity'])
			dotted_unit_mote_velocity = snapped(unit_vector.dot(structure['velocity']),.001)
			#dotted_tangent_mote_velocity = unit_tangent.dot(structure['velocity'])
			dotted_tangent_mote_velocity = snapped(unit_tangent.dot(structure['velocity']),.001)
			dotted_unit_wall_velocity = unit_vector.dot(baluster['window outline']['left']['velocity'])
			dotted_tangent_wall_velocity = unit_tangent.dot(baluster['window outline']['left']['velocity'])
			
			final_tangential_mote_velocity = dotted_tangent_mote_velocity
			final_tangential_wall_velocity = dotted_tangent_wall_velocity
			
			final_normal_mote_velocity = (dotted_unit_mote_velocity * (structure['mass'] - baluster['window outline']['left']['mass']) + 2.0 * baluster['window outline']['left']['mass'] * dotted_unit_wall_velocity ) / (structure['mass'] + baluster['window outline']['left']['mass'])
			final_normal_wall_velocity = (dotted_unit_wall_velocity * (baluster['window outline']['left']['mass'] - structure['mass']) + 2.0 * baluster['window outline']['left']['mass'] * dotted_unit_mote_velocity ) / (structure['mass'] + baluster['window outline']['left']['mass'])
			
			structure['velocity'] = (final_normal_mote_velocity * unit_vector) + (final_tangential_mote_velocity * unit_tangent)
			
		elif collision_restitution < 1.0 and collision_restitution > 0.0:
			
			x_component = -baluster['window outline']['left']['coefficient of restitution'] * structure['velocity'].x
			y_component = structure['velocity'].y
			
			structure['velocity'] = Vector2(x_component,y_component)
			
		elif collision_restitution == 0.0:
			#perfect inelastic collisions : [m1 / (m1 + m2)] * vi = vf
			
			x_component = structure['mass'] / (structure['mass'] + baluster['window outline']['left']['mass']) * structure['velocity'].x
			y_component = structure['mass'] / (structure['mass'] + baluster['window outline']['left']['mass']) * structure['velocity'].y
			
			structure['velocity'] = Vector2(x_component,y_component)
			#structure['velocity'] = final_velocity
	
	return structure['velocity']
	
	
