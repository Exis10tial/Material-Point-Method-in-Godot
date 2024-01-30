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
var dotted_unit_particle_composition_velocity : float
var dotted_tangent_particle_composition_velocity : float
var dotted_unit_wall_velocity: float
var dotted_tangent_wall_velocity: float
var final_tangential_particle_composition_velocity : float
var final_tangential_wall_velocity : float
var final_normal_particle_composition_velocity : float
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


func Collision_with_Walls(breach,particle_composition,designation,baluster,):
	###...
	if breach == 'top':
		
		collision_restitution = baluster['coefficient of restitution'] * particle_composition.coefficient_of_restitution
		collision_static_friction = baluster['coefficient of static friction'] * particle_composition.coefficient_of_static_friction
		collision_kinetic_friction = baluster['coefficient of kinetic friction'] * particle_composition.coefficient_of_kinetic_friction
		
		#if particle_composition.mechanics[designation]['velocity'].y > 0:
		
		#baluster['velocity'] = Vector2(0,0)
		#baluster['velocity'] = Vector2(0,-particle_composition.mechanics[designation]['velocity'].y)
		
		### for when collision_restitution == 1.0
		baluster['velocity'] = Vector2(0,-particle_composition.mechanics[designation]['velocity'].y) * 0.50
		
		wall_center = Vector2(particle_composition.effigy_basket[particle_composition.mechanics[designation]['correspond']].get_polygon()[particle_composition.mechanics[designation]['relation within']].x,0.0)
		
		if collision_restitution >= 1.0 :
			### the collision is perfect elastic...
			
			#normal_vector = Vector2(wall_center - particle_boundary.get_center())
			normal_vector = Vector2(wall_center - particle_composition.effigy_basket[particle_composition.mechanics[designation]['correspond']].get_polygon()[particle_composition.mechanics[designation]['relation within']])
			unit_vector = normal_vector / snapped(sqrt((snapped(pow(normal_vector.x,2.0),.01) + snapped(pow(normal_vector.y,2.0),.01))),.01)
			unit_tangent = Vector2(-unit_vector.y,unit_vector.x)
			
			#dotted_unit_particle_composition_velocity = unit_vector.dot(particle_composition.mechanics[designation]['velocity'])
			dotted_unit_particle_composition_velocity = snapped(unit_vector.dot(particle_composition.mechanics[designation]['velocity']),.001)
			#dotted_tangent_particle_composition_velocity = unit_tangent.dot(particle_composition.mechanics[designation]['velocity'])
			dotted_tangent_particle_composition_velocity = snapped(unit_tangent.dot(particle_composition.mechanics[designation]['velocity']),.001)
			dotted_unit_wall_velocity = unit_vector.dot(baluster['velocity'])
			dotted_tangent_wall_velocity = unit_tangent.dot(baluster['velocity'])
			
			final_tangential_particle_composition_velocity = dotted_tangent_particle_composition_velocity
			final_tangential_wall_velocity = dotted_tangent_wall_velocity
			final_normal_particle_composition_velocity = (dotted_unit_particle_composition_velocity * (particle_composition.mechanics[designation]['mass'] - baluster['mass']) + 2.0 * baluster['mass'] * dotted_unit_wall_velocity ) / (particle_composition.mechanics[designation]['mass'] + baluster['mass'])
			final_normal_wall_velocity = (dotted_unit_wall_velocity * (baluster['mass'] - particle_composition.mechanics[designation]['mass']) + 2.0 * baluster['mass'] * dotted_unit_particle_composition_velocity ) / (particle_composition.mechanics[designation]['mass'] + baluster['mass'])
			
			particle_composition.mechanics[designation]['velocity'] = (final_normal_particle_composition_velocity * unit_vector) + (final_tangential_particle_composition_velocity * unit_tangent)
			#final_tangential_particle_composition_velocity * unit_tangent
			#(final_normal_wall_velocity * unit_vector) + (final_tangential_wall_velocity * unit_tangent)
			
			#print(particle_composition.mechanics[designation]['velocity']," particle_composition.mechanics[designation]['velocity']")
		
		elif collision_restitution < 1.0 and collision_restitution > 0.0:
			
			
			x_component = particle_composition.mechanics[designation]['velocity'].x
			y_component = -baluster['coefficient of restitution'] *  particle_composition.mechanics[designation]['velocity'].y 
			
			particle_composition.mechanics[designation]['velocity'] = Vector2(x_component,y_component)
			
		elif collision_restitution == 0.0:
			#perfect inelastic collisions : [m1 / (m1 + m2)] * vi = vf
			
			x_component = particle_composition.mechanics[designation]['mass'] / (particle_composition.mechanics[designation]['mass'] + baluster['mass']) * particle_composition.mechanics[designation]['velocity'].x
			y_component = particle_composition.mechanics[designation]['mass'] / (particle_composition.mechanics[designation]['mass'] + baluster['mass']) * particle_composition.mechanics[designation]['velocity'].y
			
			particle_composition.mechanics[designation]['velocity'] = Vector2(x_component,y_component)
			
	if breach == 'right':
		
		collision_restitution = baluster['coefficient of restitution'] * particle_composition.coefficient_of_restitution
		collision_static_friction = baluster['coefficient of static friction'] * particle_composition.coefficient_of_static_friction
		collision_kinetic_friction = baluster['coefficient of kinetic friction'] * particle_composition.coefficient_of_kinetic_friction
		
		#baluster['velocity'] = Vector2(0,0.0)
		#baluster['velocity'] = Vector2(-particle_composition.mechanics[designation]['velocity'].x,0)
		baluster['velocity'] = Vector2(-particle_composition.mechanics[designation]['velocity'].x,0) * 0.50
		
		wall_center = Vector2(baluster['outline'],particle_composition.effigy_basket[particle_composition.mechanics[designation]['correspond']].get_polygon()[particle_composition.mechanics[designation]['relation within']].y)
		
		if collision_restitution >= 1.0 :
			### the collision is perfect elastic...
			
			#normal_vector = Vector2(wall_center - particle_boundary.get_center())
			normal_vector = Vector2(wall_center - particle_composition.effigy_basket[particle_composition.mechanics[designation]['correspond']].get_polygon()[particle_composition.mechanics[designation]['relation within']])
			
			unit_vector = normal_vector / snapped(sqrt((snapped(pow(normal_vector.x,2.0),.01) + snapped(pow(normal_vector.y,2.0),.01))),.01)
			unit_tangent = Vector2(-unit_vector.y,unit_vector.x)
			
			#dotted_unit_particle_composition_velocity = unit_vector.dot(particle_composition.mechanics[designation]['velocity'])
			dotted_unit_particle_composition_velocity = snapped(unit_vector.dot(particle_composition.mechanics[designation]['velocity']),.001)
			#dotted_tangent_particle_composition_velocity = unit_tangent.dot(particle_composition.mechanics[designation]['velocity'])
			dotted_tangent_particle_composition_velocity = snapped(unit_tangent.dot(particle_composition.mechanics[designation]['velocity']),.001)
			dotted_unit_wall_velocity = unit_vector.dot(baluster['velocity'])
			dotted_tangent_wall_velocity = unit_tangent.dot(baluster['velocity'])
			
			final_tangential_particle_composition_velocity = dotted_tangent_particle_composition_velocity
			final_tangential_wall_velocity = dotted_tangent_wall_velocity
			final_normal_particle_composition_velocity = (dotted_unit_particle_composition_velocity * (particle_composition.mechanics[designation]['mass'] - baluster['mass']) + 2.0 * baluster['mass'] * dotted_unit_wall_velocity ) / (particle_composition.mechanics[designation]['mass'] + baluster['mass'])
			final_normal_wall_velocity = (dotted_unit_wall_velocity * (baluster['mass'] - particle_composition.mechanics[designation]['mass']) + 2.0 * baluster['mass'] * dotted_unit_particle_composition_velocity ) / (particle_composition.mechanics[designation]['mass'] + baluster['mass'])
			
			particle_composition.mechanics[designation]['velocity'] = (final_normal_particle_composition_velocity * unit_vector) + (final_tangential_particle_composition_velocity * unit_tangent)
			
			#final_tangential_particle_composition_velocity * unit_tangent
			#(final_normal_wall_velocity * unit_vector) + (final_tangential_wall_velocity * unit_tangent)
			
			#print(particle_composition.mechanics[designation]['velocity']," particle_composition.mechanics[designation]['velocity']")
		
		elif collision_restitution < 1.0 and collision_restitution > 0.0:
			
			"""
			incoming_angle = rad_to_deg(asin(particle_composition.mechanics[designation]['velocity'].y/particle_composition.mechanics[designation]['velocity'].length()))
			outgoing_angle = snapped(rad_to_deg(atan(baluster['coefficient of restitution'] * tan(deg_to_rad(incoming_angle)))),.1)
			"""
			x_component =  -baluster['coefficient of restitution'] * particle_composition.mechanics[designation]['velocity'].x
			y_component = particle_composition.mechanics[designation]['velocity'].y
			particle_composition.mechanics[designation]['velocity'] = Vector2(x_component,y_component)
			
		elif collision_restitution == 0.0:
			#perfect inelastic collisions : [m1 / (m1 + m2)] * vi = vf
			
			x_component = particle_composition.mechanics[designation]['mass'] / (particle_composition.mechanics[designation]['mass'] + baluster['mass']) * particle_composition.mechanics[designation]['velocity'].x
			y_component = particle_composition.mechanics[designation]['mass'] / (particle_composition.mechanics[designation]['mass'] + baluster['mass']) * particle_composition.mechanics[designation]['velocity'].y
			
			particle_composition.mechanics[designation]['velocity'] = Vector2(x_component,y_component)
			
	if breach == 'bottom':
		
		collision_restitution = baluster['coefficient of restitution'] * particle_composition.coefficient_of_restitution
		collision_static_friction = baluster['coefficient of static friction'] * particle_composition.coefficient_of_static_friction
		collision_kinetic_friction = baluster['coefficient of kinetic friction'] * particle_composition.coefficient_of_kinetic_friction
		
		#baluster['velocity'] = Vector2(0,0)
		#baluster['velocity'] = Vector2(0,-particle_composition.mechanics[designation]['velocity'].y)
		
		baluster['velocity'] = Vector2(0,-particle_composition.mechanics[designation]['velocity'].y) * 0.50
		
		#print(collision_restitution,' check collision_restitution')
		
		wall_center = Vector2(particle_composition.effigy_basket[particle_composition.mechanics[designation]['correspond']].get_polygon()[particle_composition.mechanics[designation]['relation within']].x,baluster['outline'])
		#wall_center = Vector2(baluster['outline']/2,baluster['outline'])
		
		if collision_restitution >= 1.0 :
			### the collision is perfect elastic...
			#print('bounce back')
			#normal_vector = Vector2(wall_center - particle_boundary.get_center())
			normal_vector = Vector2(wall_center - particle_composition.effigy_basket[particle_composition.mechanics[designation]['correspond']].get_polygon()[particle_composition.mechanics[designation]['relation within']])
			#print(unit_vector,' unit_vector check')
			unit_vector = normal_vector / snapped(sqrt((snapped(pow(normal_vector.x,2.0),.01) + snapped(pow(normal_vector.y,2.0),.01))),.01)
			unit_tangent = Vector2(-unit_vector.y,unit_vector.x)
			
			#dotted_unit_particle_composition_velocity = unit_vector.dot(particle_composition.mechanics[designation]['velocity'])
			dotted_unit_particle_composition_velocity = snapped(unit_vector.dot(particle_composition.mechanics[designation]['velocity']),.001)
			#dotted_tangent_particle_composition_velocity = unit_tangent.dot(particle_composition.mechanics[designation]['velocity'])
			dotted_tangent_particle_composition_velocity = snapped(unit_tangent.dot(particle_composition.mechanics[designation]['velocity']),.001)
			dotted_unit_wall_velocity = unit_vector.dot(baluster['velocity'])
			dotted_tangent_wall_velocity = unit_tangent.dot(baluster['velocity'])
			
			final_tangential_particle_composition_velocity = dotted_tangent_particle_composition_velocity
			final_tangential_wall_velocity = dotted_tangent_wall_velocity
			
			final_normal_particle_composition_velocity = (dotted_unit_particle_composition_velocity * (particle_composition.mechanics[designation]['mass'] - baluster['mass']) + 2.0 * baluster['mass'] * dotted_unit_wall_velocity ) / (particle_composition.mechanics[designation]['mass'] + baluster['mass'])
			final_normal_wall_velocity = (dotted_unit_wall_velocity * (baluster['mass'] - particle_composition.mechanics[designation]['mass']) + 2.0 * baluster['mass'] * dotted_unit_particle_composition_velocity ) / (particle_composition.mechanics[designation]['mass'] + baluster['mass'])
			
			particle_composition.mechanics[designation]['velocity'] = (final_normal_particle_composition_velocity * unit_vector) + (final_tangential_particle_composition_velocity * unit_tangent)
			
		elif collision_restitution < 1.0 and collision_restitution > 0.0:
			#line_of_impact = wall_center
			
			x_component = particle_composition.mechanics[designation]['velocity'].x
			y_component =  -baluster['coefficient of restitution'] * particle_composition.mechanics[designation]['velocity'].y
			
			particle_composition.mechanics[designation]['velocity'] = Vector2(x_component,y_component)
			
		elif collision_restitution == 0.0:
			#perfect inelastic collisions : [m1 / (m1 + m2)] * vi = vf
			
			x_component = particle_composition.mechanics[designation]['mass'] / (particle_composition.mechanics[designation]['mass'] + baluster['mass']) * particle_composition.mechanics[designation]['velocity'].x
			y_component = particle_composition.mechanics[designation]['mass'] / (particle_composition.mechanics[designation]['mass'] + baluster['mass']) * particle_composition.mechanics[designation]['velocity'].y
			
			particle_composition.mechanics[designation]['velocity'] = Vector2(x_component,y_component)
			#print(particle_composition.mechanics[designation]['velocity']," particle_composition.mechanics[designation]['velocity']")
	
	
	if breach == 'left':
		
		collision_restitution = baluster['coefficient of restitution'] * particle_composition.coefficient_of_restitution
		collision_static_friction = baluster['coefficient of static friction'] * particle_composition.coefficient_of_static_friction
		collision_kinetic_friction = baluster['coefficient of kinetic friction'] * particle_composition.coefficient_of_kinetic_friction
		
		#baluster['velocity'] = Vector2(0,0)
		#baluster['velocity'] = Vector2(-particle_composition.mechanics[designation]['velocity'].x,0)
		baluster['velocity'] = Vector2(-particle_composition.mechanics[designation]['velocity'].x,0) * 0.50
		
		wall_center = Vector2(0.0,particle_composition.effigy_basket[particle_composition.mechanics[designation]['correspond']].get_polygon()[particle_composition.mechanics[designation]['relation within']].y)
		
		if collision_restitution >= 1.0 :
			### the collision is perfect elastic...
			
			#normal_vector = Vector2(wall_center - particle_boundary.get_center())
			normal_vector = Vector2(wall_center - particle_composition.effigy_basket[particle_composition.mechanics[designation]['correspond']].get_polygon()[particle_composition.mechanics[designation]['relation within']])
			
			unit_vector = normal_vector / snapped(sqrt((snapped(pow(normal_vector.x,2.0),.01) + snapped(pow(normal_vector.y,2.0),.01))),.01)
			unit_tangent = Vector2(-unit_vector.y,unit_vector.x)
			
			#dotted_unit_particle_composition_velocity = unit_vector.dot(particle_composition.mechanics[designation]['velocity'])
			dotted_unit_particle_composition_velocity = snapped(unit_vector.dot(particle_composition.mechanics[designation]['velocity']),.001)
			#dotted_tangent_particle_composition_velocity = unit_tangent.dot(particle_composition.mechanics[designation]['velocity'])
			dotted_tangent_particle_composition_velocity = snapped(unit_tangent.dot(particle_composition.mechanics[designation]['velocity']),.001)
			dotted_unit_wall_velocity = unit_vector.dot(baluster['velocity'])
			dotted_tangent_wall_velocity = unit_tangent.dot(baluster['velocity'])
			
			final_tangential_particle_composition_velocity = dotted_tangent_particle_composition_velocity
			final_tangential_wall_velocity = dotted_tangent_wall_velocity
			
			final_normal_particle_composition_velocity = (dotted_unit_particle_composition_velocity * (particle_composition.mechanics[designation]['mass'] - baluster['mass']) + 2.0 * baluster['mass'] * dotted_unit_wall_velocity ) / (particle_composition.mechanics[designation]['mass'] + baluster['mass'])
			final_normal_wall_velocity = (dotted_unit_wall_velocity * (baluster['mass'] - particle_composition.mechanics[designation]['mass']) + 2.0 * baluster['mass'] * dotted_unit_particle_composition_velocity ) / (particle_composition.mechanics[designation]['mass'] + baluster['mass'])
			
			particle_composition.mechanics[designation]['velocity'] = (final_normal_particle_composition_velocity * unit_vector) + (final_tangential_particle_composition_velocity * unit_tangent)
			
		elif collision_restitution < 1.0 and collision_restitution > 0.0:
			
			x_component = -baluster['coefficient of restitution'] * particle_composition.mechanics[designation]['velocity'].x
			y_component = particle_composition.mechanics[designation]['velocity'].y
			
			particle_composition.mechanics[designation]['velocity'] = Vector2(x_component,y_component)
			
		elif collision_restitution == 0.0:
			#perfect inelastic collisions : [m1 / (m1 + m2)] * vi = vf
			
			x_component = particle_composition.mechanics[designation]['mass'] / (particle_composition.mechanics[designation]['mass'] + baluster['mass']) * particle_composition.mechanics[designation]['velocity'].x
			y_component = particle_composition.mechanics[designation]['mass'] / (particle_composition.mechanics[designation]['mass'] + baluster['mass']) * particle_composition.mechanics[designation]['velocity'].y
			
			particle_composition.mechanics[designation]['velocity'] = Vector2(x_component,y_component)
			#particle_composition.mechanics[designation]['velocity'] = final_velocity
	
	return particle_composition.mechanics[designation]['velocity']
	
	
