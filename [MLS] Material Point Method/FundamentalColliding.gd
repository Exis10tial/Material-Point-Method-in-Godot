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
		
		### used for collision_restitution ( < 1 and > 0 )...
		baluster['velocity'] = Vector2(0,-particle_composition.mechanics[designation]['velocity'].y) * 2.00
		
		wall_center = Vector2(particle_composition.entity_container[designation].position.x,0.0)
		
		if collision_restitution >= 1.0 :
			### the collision is perfect elastic...
			
			#normal_vector = Vector2(wall_center - particle_boundary.get_center())
			normal_vector = Vector2(wall_center - particle_composition.entity_container[designation].position)
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
			# imperflect inelastic collisions :
			### zero_momentum_frame... 
			var y_component
			### acquire the zero momentum velocity...
			var zero_particle_coefficient = particle_composition.mechanics[designation]['mass'] * particle_composition.mechanics[designation]['velocity'].length()
			var zero_wall_coefficient = baluster['mass'] * baluster['velocity'].length()
			var zero_mass_coefficient =  particle_composition.mechanics[designation]['mass'] + baluster['mass'] 
			var zero_momentum_frame_velocity =  (zero_particle_coefficient + zero_wall_coefficient) / zero_mass_coefficient
			### find the zero momentum velocity before the collision...
			var zero_particle_velocity = particle_composition.mechanics[designation]['velocity'].length() - zero_momentum_frame_velocity
			var zero_wall_velocity = baluster['velocity'].length() - zero_momentum_frame_velocity
			
			var incoming_particle_angle = rad_to_deg(atan2(particle_composition.mechanics[designation]['velocity'].y,particle_composition.mechanics[designation]['velocity'].x))
			
			var angle_coefficient = collision_restitution * tan(deg_to_rad(incoming_particle_angle))
			var outgoing_angle = rad_to_deg(atan(angle_coefficient))
			
			var particle_reformed_coefficient = zero_momentum_frame_velocity + zero_particle_velocity
			
			var x_component = ((wall_center.x + particle_reformed_coefficient * cos(deg_to_rad(outgoing_angle))) - wall_center.x )* sign(particle_composition.mechanics[designation]['velocity'].x)
			
			if incoming_particle_angle > -90:
				y_component = (wall_center.y + particle_reformed_coefficient * -sin(deg_to_rad(outgoing_angle))) - wall_center.y
			else:
				y_component = (wall_center.y + particle_reformed_coefficient * sin(deg_to_rad(outgoing_angle))) - wall_center.y
			
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
		
		### used for collision_restitution ( < 1 and > 0 )...
		baluster['velocity'] = Vector2(-particle_composition.mechanics[designation]['velocity'].x,0) * 2.00
		
		wall_center = Vector2(baluster['outline'],particle_composition.entity_container[designation].position.y)
		
		if collision_restitution >= 1.0 :
			### the collision is perfect elastic...
			
			#normal_vector = Vector2(wall_center - particle_boundary.get_center())
			normal_vector = Vector2(wall_center - particle_composition.entity_container[designation].position)
			
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
			#### zero_momentum_frame... 
			var x_component
			### acquire the zero momentum velocity...
			var zero_particle_coefficient = particle_composition.mechanics[designation]['mass'] * particle_composition.mechanics[designation]['velocity'].length()
			var zero_wall_coefficient = baluster['mass'] * baluster['velocity'].length()
			var zero_mass_coefficient =  particle_composition.mechanics[designation]['mass'] + baluster['mass'] 
			var zero_momentum_frame_velocity =  (zero_particle_coefficient + zero_wall_coefficient) / zero_mass_coefficient
			### find the zero momentum velocity before the collision...
			var zero_particle_velocity = particle_composition.mechanics[designation]['velocity'].length() - zero_momentum_frame_velocity
			var zero_wall_velocity = baluster['velocity'].length() - zero_momentum_frame_velocity
			
			var incoming_particle_angle = rad_to_deg(atan2(particle_composition.mechanics[designation]['velocity'].y,particle_composition.mechanics[designation]['velocity'].y))
			
			var angle_coefficient = collision_restitution * tan(deg_to_rad(incoming_particle_angle))
			var outgoing_angle = rad_to_deg(atan(angle_coefficient))
			
			var particle_reformed_coefficient = zero_momentum_frame_velocity + zero_particle_velocity
			
			var y_component = ((wall_center.x + particle_reformed_coefficient * sin(deg_to_rad(outgoing_angle))) - wall_center.x )
			
			if incoming_particle_angle > 0:
				x_component = (wall_center.x + particle_reformed_coefficient * -cos(deg_to_rad(outgoing_angle))) - wall_center.x
			else:
				x_component = (wall_center.x + particle_reformed_coefficient * cos(deg_to_rad(outgoing_angle))) - wall_center.x
			
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
		
		### used for collision_restitution ( < 1 and > 0 )...
		baluster['velocity'] = Vector2(0,-particle_composition.mechanics[designation]['velocity'].y) * 2.0
		### acts as the impact center...
		wall_center = Vector2(particle_composition.entity_container[designation].position.x,baluster['barrier'])
		#wall_center = Vector2(baluster['outline']/2,baluster['outline'])
		
		if collision_restitution >= 1.0 :
			### the collision is perfect elastic...
			#print('bounce back')
			#normal_vector = Vector2(wall_center - particle_boundary.get_center())
			normal_vector = Vector2(wall_center - particle_composition.entity_container[designation].position)
			#print(unit_vector,' unit_vector check')
			unit_vector = normal_vector / snapped(sqrt((snapped(pow(normal_vector.x,2.0),.01) + snapped(pow(normal_vector.y,2.0),.01))),.01)
			#print(unit_vector,' unit_vector check')
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
			### zero_momentum_frame... 
			var y_component
			### acquire the zero momentum velocity...
			var zero_particle_coefficient = particle_composition.mechanics[designation]['mass'] * particle_composition.mechanics[designation]['velocity'].length()
			var zero_wall_coefficient = baluster['mass'] * baluster['velocity'].length()
			var zero_mass_coefficient =  particle_composition.mechanics[designation]['mass'] + baluster['mass'] 
			var zero_momentum_frame_velocity =  (zero_particle_coefficient + zero_wall_coefficient) / zero_mass_coefficient
			### find the zero momentum velocity before the collision...
			var zero_particle_velocity = particle_composition.mechanics[designation]['velocity'].length() - zero_momentum_frame_velocity
			var zero_wall_velocity = baluster['velocity'].length() - zero_momentum_frame_velocity
			
			var incoming_particle_angle = rad_to_deg(atan2(particle_composition.mechanics[designation]['velocity'].y,particle_composition.mechanics[designation]['velocity'].x))
			
			var angle_coefficient = collision_restitution * tan(deg_to_rad(incoming_particle_angle))
			var outgoing_angle = rad_to_deg(atan(angle_coefficient))
			
			var particle_reformed_coefficient = zero_momentum_frame_velocity + zero_particle_velocity
			
			var x_component = ((wall_center.x + particle_reformed_coefficient * cos(deg_to_rad(outgoing_angle))) - wall_center.x )* sign(particle_composition.mechanics[designation]['velocity'].x)
			if incoming_particle_angle < 90:
				y_component = (wall_center.y + particle_reformed_coefficient * -sin(deg_to_rad(outgoing_angle))) - wall_center.y
			else:
				y_component = (wall_center.y + particle_reformed_coefficient * sin(deg_to_rad(outgoing_angle))) - wall_center.y
			
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
		
		### used for collision_restitution ( < 1 and > 0 )...
		baluster['velocity'] = Vector2(particle_composition.mechanics[designation]['velocity'].x,0) * 0
		
		wall_center = Vector2(0.0,particle_composition.entity_container[designation].position.y)
		
		if collision_restitution >= 1.0 :
			### the collision is perfect elastic...
			
			#normal_vector = Vector2(wall_center - particle_boundary.get_center())
			normal_vector = Vector2(wall_center - particle_composition.entity_container[designation].position)
			
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
			#### zero_momentum_frame... 
			
			### acquire the zero momentum velocity...
			var zero_particle_coefficient = particle_composition.mechanics[designation]['mass'] * particle_composition.mechanics[designation]['velocity'].length()
			var zero_wall_coefficient = baluster['mass'] * baluster['velocity'].length()
			var zero_mass_coefficient =  particle_composition.mechanics[designation]['mass'] + baluster['mass'] 
			var zero_momentum_frame_velocity =  (zero_particle_coefficient + zero_wall_coefficient) / zero_mass_coefficient
			### find the zero momentum velocity before the collision...
			var zero_particle_velocity = particle_composition.mechanics[designation]['velocity'].length() - zero_momentum_frame_velocity
			var zero_wall_velocity = baluster['velocity'].length() - zero_momentum_frame_velocity
			
			var incoming_particle_angle = rad_to_deg(atan2(particle_composition.mechanics[designation]['velocity'].y,particle_composition.mechanics[designation]['velocity'].x))
			
			var angle_coefficient = collision_restitution * tan(deg_to_rad(incoming_particle_angle))
			var outgoing_angle = rad_to_deg(atan(angle_coefficient))
			
			var particle_reformed_coefficient = zero_momentum_frame_velocity + zero_particle_velocity
			
			var y_component = ((wall_center.y + particle_reformed_coefficient * -sin(deg_to_rad(outgoing_angle))) - wall_center.y )# * sign(particle_composition.mechanics[designation]['velocity'].y)
			var x_component = (wall_center.x + particle_reformed_coefficient * cos(deg_to_rad(outgoing_angle))) - wall_center.x
			
			particle_composition.mechanics[designation]['velocity'] = Vector2(x_component,y_component)
			
		elif collision_restitution == 0.0:
			#perfect inelastic collisions : [m1 / (m1 + m2)] * vi = vf
			
			x_component = particle_composition.mechanics[designation]['mass'] / (particle_composition.mechanics[designation]['mass'] + baluster['mass']) * particle_composition.mechanics[designation]['velocity'].x
			y_component = particle_composition.mechanics[designation]['mass'] / (particle_composition.mechanics[designation]['mass'] + baluster['mass']) * particle_composition.mechanics[designation]['velocity'].y
			
			particle_composition.mechanics[designation]['velocity'] = Vector2(x_component,y_component)
			#particle_composition.mechanics[designation]['velocity'] = final_velocity
	
	return particle_composition.mechanics[designation]['velocity']
	
	


func Apply_Friction(surrounding_wall,in_contact_with_wall,wall_side_2,wall_side_3,particle_mass,particle_velocity,gravitation,substance_static_coefficient,substance_kinetic_coefficient):
	### find the angle of the wall...
	var angle_of_wall
	var wall_start
	var wall_end
	var particle_force 
	
	if in_contact_with_wall == 'top' or in_contact_with_wall == 'bottom':
		wall_start = Vector2(surrounding_wall[wall_side_2]['barrier'],surrounding_wall[in_contact_with_wall]['barrier'],)
		wall_end = Vector2(surrounding_wall[wall_side_3]['barrier'],surrounding_wall[in_contact_with_wall]['barrier'])
	elif in_contact_with_wall == 'left' or in_contact_with_wall == 'right':
		wall_start = Vector2(surrounding_wall[in_contact_with_wall]['barrier'],surrounding_wall[wall_side_2]['barrier'])
		wall_end = Vector2(surrounding_wall[in_contact_with_wall]['barrier'],surrounding_wall[wall_side_3]['barrier'])
	
	if in_contact_with_wall == 'top' or in_contact_with_wall == 'bottom':
		var midpoint_of_wall = Vector2((wall_start.x + wall_end.x) / 2.0,(wall_start.y + wall_end.y) / 2.0)
		var wall_half_start = wall_start - midpoint_of_wall
		var wall_half_end = wall_end - midpoint_of_wall
					
		angle_of_wall = atan2( (wall_half_end.y*wall_half_start.x) - (wall_end.x*wall_half_start.y) , (wall_half_end.x*wall_half_start.x) + (wall_half_end.y*wall_half_start.y))
	elif in_contact_with_wall == 'left' or in_contact_with_wall == 'right':
		#angle_of_wall = wall_start.angle_to_point(wall_end)
		#print(angle_of_wall)
		var midpoint_of_wall = Vector2((wall_start.x + wall_end.x) / 2.0,(wall_start.y + wall_end.y) / 2.0)
		var wall_half_start = wall_start - midpoint_of_wall
		var wall_half_end = wall_end - midpoint_of_wall
					
		angle_of_wall = atan2( (wall_half_end.y*wall_half_start.x) - (wall_end.x*wall_half_start.y) , (wall_half_end.x*wall_half_start.x) + (wall_half_end.y*wall_half_start.y))
	
	if in_contact_with_wall == 'top' or in_contact_with_wall == 'bottom': 
		
		### find the angle of the particle velocity
		var angle_of_particle = atan2(particle_velocity.y,particle_velocity.x)
		
		#print(int(rad_to_deg(angle_of_particle)),' particle')
		#print(int(rad_to_deg(angle_of_wall)),' angle_of_wall')
		### if the particle is moving parallel to the wall in either direction...
		if int(rad_to_deg(angle_of_particle)) in range(rad_to_deg(angle_of_wall)-2,rad_to_deg(angle_of_wall)+3,1):
			
			### if the particle is moving parallel to the wall and is in contact with the wall...
			var weight = particle_mass * gravitation.length()
			var normal_force = weight
			var wall_friction
			var particle_friction
			### The wall friction depends on the direction of the wall...
			### friction is in the opposite direction of force...
			if surrounding_wall[in_contact_with_wall]['velocity'].x == 0:
				### if the wall isn't moving..
				wall_friction = normal_force * surrounding_wall[in_contact_with_wall]['coefficient of static friction']
			elif surrounding_wall[in_contact_with_wall]['velocity'].x > 0:
				### if the wall is moving positive...
				wall_friction = -1.0 * normal_force * surrounding_wall[in_contact_with_wall]['coefficient of kinetic friction']
			elif surrounding_wall[in_contact_with_wall]['velocity'].x < 0:
				### if the wall is moving negative...
				wall_friction = normal_force * surrounding_wall[in_contact_with_wall]['coefficient of kinetic friction']
			### The particle friction depends on the direction of the wall...
			if particle_velocity.x == 0:
				particle_friction = normal_force * substance_static_coefficient
			elif particle_velocity.x > 0:
				particle_friction = -1.0 * normal_force * substance_kinetic_coefficient
			elif particle_velocity.x < 0:
				particle_friction = normal_force * substance_kinetic_coefficient
			
			var total_friction_force = wall_friction + particle_friction
						
			particle_force = particle_velocity.length() * particle_mass
			
			if particle_force >= total_friction_force:
				particle_force = particle_force - total_friction_force
			elif particle_force < total_friction_force:
				particle_force = 0.0
			###
			return  Vector2( cos(int(rad_to_deg(angle_of_particle))) * particle_force,
			sin(int(rad_to_deg(angle_of_particle))) * particle_force)
		else:
			return particle_velocity

	if in_contact_with_wall == 'left' or in_contact_with_wall == 'right': 
		### if the particle is moving parallel to the wall and is in contact with the wall...
		var weight = particle_mass * gravitation.length()
		var normal_force = weight
		var wall_friction
		var particle_friction
		
		### find the angle of the particle velocity
		var angle_of_particle = atan2(particle_velocity.y,particle_velocity.x)
		
		#print(int(rad_to_deg(angle_of_particle)),' particle')
		#print(int(rad_to_deg(angle_of_wall)),' angle_of_wall')
		### if the particle is moving parallel to the wall in either direction...
		if int(rad_to_deg(angle_of_particle)) in range(rad_to_deg(angle_of_wall)-2,rad_to_deg(angle_of_wall)+3,1) in range(0-2,0+3,1):
				
			### The wall friction depends on the direction of the wall...
			### friction is in the opposite direction of force...
			if surrounding_wall[in_contact_with_wall]['velocity'].y == 0:
				### if the wall isn't moving..
				wall_friction = normal_force * surrounding_wall[in_contact_with_wall]['coefficient of static friction']
			elif surrounding_wall[in_contact_with_wall]['velocity'].y > 0:
				### if the wall is moving positive...
				wall_friction = -1.0 * normal_force * surrounding_wall[in_contact_with_wall]['coefficient of kinetic friction']
			elif surrounding_wall[in_contact_with_wall]['velocity'].y < 0:
				### if the wall is moving negative...
				wall_friction = normal_force * surrounding_wall[in_contact_with_wall]['coefficient of kinetic friction']
			### The particle friction depends on the direction of the wall...
			if particle_velocity.y == 0:
				particle_friction = normal_force * substance_static_coefficient
			elif particle_velocity.y > 0:
				particle_friction = -1.0 * normal_force * substance_kinetic_coefficient
			elif particle_velocity.y < 0:
				particle_friction = normal_force * substance_kinetic_coefficient
				
			var total_friction_force = wall_friction + particle_friction
							
			particle_force = particle_velocity.length() * particle_mass
				
			if particle_force >= total_friction_force:
				particle_force = particle_force - total_friction_force
			elif particle_force < total_friction_force:
				particle_force = 0.0
			
			###
			return  Vector2( cos(int(rad_to_deg(angle_of_particle))) * particle_force,
			sin(int(rad_to_deg(angle_of_particle))) * particle_force)
		else:
			return particle_velocity
