extends Node


### determine the type of collision...
var collision_restitution
var wall_center : Vector2
### particle merging...
var length_covered_area : float 
var width_covered_area : float
var covered_area : float
var percentage_covered : float
### used for elastic collisions <= 1.0 ...
var elastic_term_1_of_vector_1 : float
var elastic_term_1_of_vector_2 : float
var elastic_mass_addition : float
var elastic_dot_vector_1
var elastic_dot_vector_2
var elastic_position_kernel_of_vector_1 : Vector2
var elastic_position_kernel_of_vector_2 : Vector2
var elastic_velocity_kernel_of_vector_1 : Vector2
var elastic_velocity_kernel_of_vector_2 : Vector2
var elastic_magnitude_kernel_of_vector_1
var elastic_magnitude_kernel_of_vector_2
### used for impartial collisions < 1.0 and > 0.0 ...
var line_of_impact_x : float
var line_of_impact_y : float
var line_of_impact : Vector2
var x_component : float
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








func _on_handle_collisions_ready():
	pass # Replace with function body.




func Collision_with_Walls(breach,mote,baluster,structure,cell_size,gravity_force):
	###...
	if breach == 'top':
		##print(' check top')
		collision_restitution = baluster['window outline']['top']['coefficient of restitution'] * mote.coefficient_of_restitution
		if collision_restitution >= 1.0 :
			wall_center = Vector2(mote.position.x,0.0)
			baluster['window outline']['top']['velocity'] = mote.velocity * 0.0
			baluster['window outline']['top']['outline'] = 0 + (cell_size)# / 2.0)
			
			var normal_vector = Vector2(wall_center - mote.surrounding_area.get_center())
			##print(normal_vector,' normal')
			var unit_vector = normal_vector / (snapped(sqrt(pow(normal_vector.x,2.0)),.01) + snapped(sqrt(pow(normal_vector.y,2.0)),.01))
			##print(unit_vector,' unit')
			##print(structure[mote]['velocity'],' right after wall contact 1.0')
			var unit_tangent = Vector2(-unit_vector.y,unit_vector.x)
			##print(unit_tangent,' unit_tangent')
							
			var dotted_unit_mote_velocity = unit_vector.dot(structure[mote]['velocity'])
			var dotted_tangent_mote_velocity = unit_tangent.dot(structure[mote]['velocity'])
			##print(dotted_unit_mote_velocity,' dot unit mote velocity')
			##print(dotted_tangent_mote_velocity,' dotted_tangent_mote_velocity')
			var dotted_unit_wall_velocity = unit_vector.dot(baluster['window outline']['top']['velocity'])
			var dotted_tangent_wall_velocity = unit_tangent.dot(baluster['window outline']['top']['velocity'])
			##print(dotted_unit_wall_velocity,' dotted_unit_wall_velocity')
			##print(dotted_tangent_wall_velocity,' dotted_tangent_wall_velocity')
							
			var final_tangential_mote_velocity = dotted_tangent_mote_velocity
			var final_tangential_wall_velocity = dotted_tangent_wall_velocity
							
							
			var final_normal_mote_velocity = (dotted_unit_mote_velocity * (structure[mote]['mass'] - baluster['window outline']['top']['mass']) + 2.0 * baluster['window outline']['top']['mass'] * dotted_unit_wall_velocity ) / (structure[mote]['mass'] + baluster['window outline']['right']['mass'])
			var final_normal_wall_velocity = (dotted_unit_wall_velocity * (baluster['window outline']['top']['mass'] - structure[mote]['mass']) + 2.0 * baluster['window outline']['top']['mass'] * dotted_unit_mote_velocity ) / (structure[mote]['mass'] + baluster['window outline']['right']['mass'])
							
			##print(final_normal_mote_velocity,' final_normal_mote_velocity')
			##print(final_normal_wall_velocity,' final_normal_wall_velocity')
							
							
			structure[mote]['velocity'] = (final_normal_mote_velocity * unit_vector) + (final_tangential_mote_velocity * unit_tangent)
			#final_tangential_mote_velocity * unit_tangent
							
			#(final_normal_wall_velocity * unit_vector) + (final_tangential_wall_velocity * unit_tangent)
							
				#
			##print(structure[mote]['velocity']," structure[mote]['velocity']")
		
		
		elif collision_restitution < 1.0 and collision_restitution > 0.0:
			line_of_impact = wall_center
			incoming_angle = rad_to_deg(mote.get_position().angle_to(line_of_impact))
			outgoing_angle = rad_to_deg(atan(baluster['window outline']['top']['coefficient of restitution'] * incoming_angle ))
			x_leg_of_particle_velocity = structure[mote]['velocity'].y * tan(outgoing_angle)
			# determine which friction static or kinetic...
			#if structure[particle]['velocity'].x in range(-1,2) and structure[particle]['velocity'].y in range(-1,2):
			if structure[mote]['velocity'].x in range(-cell_size,cell_size+1) and structure[mote]['velocity'].y in range(-cell_size,cell_size+1):
				### using of static friction...
				x_component = x_leg_of_particle_velocity * static_friction
			else:
				### use of kinetic friction...
				x_component = x_leg_of_particle_velocity * kinetic_friction
							
			structure[mote]['velocity'] = structure[mote]['velocity'] - Vector2(x_component,structure[mote]['velocity'].y)
						
			
		elif collision_restitution == 0.0:
			final_velocity = ( (structure[mote]['mass'] * structure[mote]['velocity']) + (baluster['window outline']['top']['mass'] * baluster['window outline']['top']['velocity']) ) / (structure[mote]['mass'] * baluster['window outline']['top']['mass']) 
						
			structure[mote]['velocity'] = structure[mote]['velocity'] - final_velocity
						
			
	if breach == 'right':
		##print(' check right')
		collision_restitution = baluster['window outline']['right']['coefficient of restitution'] * mote.coefficient_of_restitution
		wall_center = Vector2(baluster['window outline']['right']['outline'],mote.position.y)
		if collision_restitution >= 1.0 :
			
			baluster['window outline']['right']['velocity'] = mote.velocity * 0.0
			baluster['window outline']['right']['outline'] = ProjectSettings.get_setting('display/window/size/width') - (cell_size)# / 2.0)
			
			var normal_vector = Vector2(wall_center - mote.surrounding_area.get_center())
			###print(normal_vector,' normal')
			var unit_vector = normal_vector / (snapped(sqrt(pow(normal_vector.x,2.0)),.01) + snapped(sqrt(pow(normal_vector.y,2.0)),.01))
			##print(unit_vector,' unit')
			##print(structure[mote]['velocity'],' right after wall contact 1.0')
			var unit_tangent = Vector2(-unit_vector.y,unit_vector.x)
			##print(unit_tangent,' unit_tangent')
							
			var dotted_unit_mote_velocity = unit_vector.dot(structure[mote]['velocity'])
			var dotted_tangent_mote_velocity = unit_tangent.dot(structure[mote]['velocity'])
			##print(dotted_unit_mote_velocity,' dot unit mote velocity')
			#print(dotted_tangent_mote_velocity,' dotted_tangent_mote_velocity')
			var dotted_unit_wall_velocity = unit_vector.dot(baluster['window outline']['right']['velocity'])
			var dotted_tangent_wall_velocity = unit_tangent.dot(baluster['window outline']['right']['velocity'])
			#print(dotted_unit_wall_velocity,' dotted_unit_wall_velocity')
			#print(dotted_tangent_wall_velocity,' dotted_tangent_wall_velocity')
							
			var final_tangential_mote_velocity = dotted_tangent_mote_velocity
			var final_tangential_wall_velocity = dotted_tangent_wall_velocity
							
							
			var final_normal_mote_velocity = (dotted_unit_mote_velocity * (structure[mote]['mass'] - baluster['window outline']['right']['mass']) + 2.0 * baluster['window outline']['right']['mass'] * dotted_unit_wall_velocity ) / (structure[mote]['mass'] + baluster['window outline']['right']['mass'])
			var final_normal_wall_velocity = (dotted_unit_wall_velocity * (baluster['window outline']['right']['mass'] - structure[mote]['mass']) + 2.0 * baluster['window outline']['right']['mass'] * dotted_unit_mote_velocity ) / (structure[mote]['mass'] + baluster['window outline']['right']['mass'])
							
			#print(final_normal_mote_velocity,' final_normal_mote_velocity')
			#print(final_normal_wall_velocity,' final_normal_wall_velocity')
							
							
			structure[mote]['velocity'] = (final_normal_mote_velocity * unit_vector) + (final_tangential_mote_velocity * unit_tangent)
			#final_tangential_mote_velocity * unit_tangent
							
			#(final_normal_wall_velocity * unit_vector) + (final_tangential_wall_velocity * unit_tangent)
							
				#
			#print(structure[mote]['velocity']," structure[mote]['velocity']")
		
		elif collision_restitution < 1.0 and collision_restitution > 0.0:
			line_of_impact = wall_center
			incoming_angle = rad_to_deg(mote.get_position().angle_to(line_of_impact))
			outgoing_angle = rad_to_deg(atan(baluster['window outline']['right']['coefficient of restitution'] * incoming_angle ))
			x_leg_of_particle_velocity = structure[mote]['velocity'].y * tan(outgoing_angle)
			# determine which friction static or kinetic...
			#if structure[particle]['velocity'].x in range(-1,2) and structure[particle]['velocity'].y in range(-1,2):
			if structure[mote]['velocity'].x in range(-cell_size,cell_size+1) and structure[mote]['velocity'].y in range(-cell_size,cell_size+1):
				### using of static friction...
				x_component = x_leg_of_particle_velocity * static_friction
			else:
				### use of kinetic friction...
				x_component = x_leg_of_particle_velocity * kinetic_friction
							
			structure[mote]['velocity'] = structure[mote]['velocity'] - Vector2(x_component,structure[mote]['velocity'].y)
						
			
		elif collision_restitution == 0.0:
			final_velocity = ( (structure[mote]['mass'] * structure[mote]['velocity']) + (baluster['window outline']['right']['mass'] * baluster['window outline']['right']['velocity']) ) / (structure[mote]['mass'] * baluster['window outline']['right']['mass']) 
						
			structure[mote]['velocity'] = structure[mote]['velocity'] - final_velocity
						
			
	if breach == 'bottom':
		#print(' check bottom')
		collision_restitution = baluster['window outline']['bottom']['coefficient of restitution'] * mote.coefficient_of_restitution
		wall_center = Vector2(mote.position.x,baluster['window outline']['bottom']['outline'])
		if collision_restitution >= 1.0 :
			
			baluster['window outline']['bottom']['velocity'] = mote.velocity * 0.0
			baluster['window outline']['bottom']['outline'] = ProjectSettings.get_setting('display/window/size/height') - (cell_size)# / 2.0)
			
			var normal_vector = Vector2(wall_center - mote.surrounding_area.get_center())
			#print(normal_vector,' normal')
			var unit_vector = normal_vector / (snapped(sqrt(pow(normal_vector.x,2.0)),.01) + snapped(sqrt(pow(normal_vector.y,2.0)),.01))
			#print(unit_vector,' unit')
			##print(structure[mote]['velocity'],' right after wall contact 1.0')
			var unit_tangent = Vector2(-unit_vector.y,unit_vector.x)
			#print(unit_tangent,' unit_tangent')
							
			var dotted_unit_mote_velocity = unit_vector.dot(structure[mote]['velocity'])
			var dotted_tangent_mote_velocity = unit_tangent.dot(structure[mote]['velocity'])
			#print(dotted_unit_mote_velocity,' dot unit mote velocity')
			#print(dotted_tangent_mote_velocity,' dotted_tangent_mote_velocity')
			var dotted_unit_wall_velocity = unit_vector.dot(baluster['window outline']['bottom']['velocity'])
			var dotted_tangent_wall_velocity = unit_tangent.dot(baluster['window outline']['bottom']['velocity'])
			#print(dotted_unit_wall_velocity,' dotted_unit_wall_velocity')
			#print(dotted_tangent_wall_velocity,' dotted_tangent_wall_velocity')
							
			var final_tangential_mote_velocity = dotted_tangent_mote_velocity
			var final_tangential_wall_velocity = dotted_tangent_wall_velocity
							
							
			var final_normal_mote_velocity = (dotted_unit_mote_velocity * (structure[mote]['mass'] - baluster['window outline']['bottom']['mass']) + 2.0 * baluster['window outline']['bottom']['mass'] * dotted_unit_wall_velocity ) / (structure[mote]['mass'] + baluster['window outline']['bottom']['mass'])
			var final_normal_wall_velocity = (dotted_unit_wall_velocity * (baluster['window outline']['bottom']['mass'] - structure[mote]['mass']) + 2.0 * baluster['window outline']['bottom']['mass'] * dotted_unit_mote_velocity ) / (structure[mote]['mass'] + baluster['window outline']['bottom']['mass'])
							
			#print(final_normal_mote_velocity,' final_normal_mote_velocity')
			#print(final_normal_wall_velocity,' final_normal_wall_velocity')
							
							
			structure[mote]['velocity'] = (final_normal_mote_velocity * unit_vector) + (final_tangential_mote_velocity * unit_tangent)
			#final_tangential_mote_velocity * unit_tangent
							
			#(final_normal_wall_velocity * unit_vector) + (final_tangential_wall_velocity * unit_tangent)
							
				#
			#print(structure[mote]['velocity']," structure[mote]['velocity']")
		
		elif collision_restitution < 1.0 and collision_restitution > 0.0:
			line_of_impact = wall_center
			incoming_angle = rad_to_deg(mote.get_position().angle_to(line_of_impact))
			outgoing_angle = rad_to_deg(atan(baluster['window outline']['bottom']['coefficient of restitution'] * incoming_angle ))
			x_leg_of_particle_velocity = structure[mote]['velocity'].y * tan(outgoing_angle)
			# determine which friction static or kinetic...
			#if structure[particle]['velocity'].x in range(-1,2) and structure[particle]['velocity'].y in range(-1,2):
			if structure[mote]['velocity'].x in range(-cell_size,cell_size+1) and structure[mote]['velocity'].y in range(-cell_size,cell_size+1):
				### using of static friction...
				x_component = x_leg_of_particle_velocity * static_friction
			else:
				### use of kinetic friction...
				x_component = x_leg_of_particle_velocity * kinetic_friction
							
			structure[mote]['velocity'] = structure[mote]['velocity'] - Vector2(x_component,structure[mote]['velocity'].y)
						
		elif collision_restitution == 0.0:
			final_velocity = ( (structure[mote]['mass'] * structure[mote]['velocity']) + (baluster['window outline']['bottom']['mass'] * baluster['window outline']['bottom']['velocity']) ) / (structure[mote]['mass'] * baluster['window outline']['bottom']['mass']) 
						
			structure[mote]['velocity'] = structure[mote]['velocity'] - final_velocity
						
	
	if breach == 'left':
		#print(' check left')
		collision_restitution = baluster['window outline']['left']['coefficient of restitution'] * mote.coefficient_of_restitution
		wall_center = Vector2(0.0,mote.position.y)
		if collision_restitution >= 1.0 :
			
			baluster['window outline']['left']['velocity'] = mote.velocity * 0.0
			baluster['window outline']['left']['outline'] = 0 + (cell_size)# / 2.0)
			
			var normal_vector = Vector2(wall_center - mote.surrounding_area.get_center())
			#print(normal_vector,' normal')
			var unit_vector = normal_vector / (snapped(sqrt(pow(normal_vector.x,2.0)),.01) + snapped(sqrt(pow(normal_vector.y,2.0)),.01))
			#print(unit_vector,' unit')
			##print(structure[mote]['velocity'],' right after wall contact 1.0')
			var unit_tangent = Vector2(-unit_vector.y,unit_vector.x)
			#print(unit_tangent,' unit_tangent')
							
			var dotted_unit_mote_velocity = unit_vector.dot(structure[mote]['velocity'])
			var dotted_tangent_mote_velocity = unit_tangent.dot(structure[mote]['velocity'])
			#print(dotted_unit_mote_velocity,' dot unit mote velocity')
			#print(dotted_tangent_mote_velocity,' dotted_tangent_mote_velocity')
			var dotted_unit_wall_velocity = unit_vector.dot(baluster['window outline']['left']['velocity'])
			var dotted_tangent_wall_velocity = unit_tangent.dot(baluster['window outline']['left']['velocity'])
			#print(dotted_unit_wall_velocity,' dotted_unit_wall_velocity')
			#print(dotted_tangent_wall_velocity,' dotted_tangent_wall_velocity')
							
			var final_tangential_mote_velocity = dotted_tangent_mote_velocity
			var final_tangential_wall_velocity = dotted_tangent_wall_velocity
							
							
			var final_normal_mote_velocity = (dotted_unit_mote_velocity * (structure[mote]['mass'] - baluster['window outline']['left']['mass']) + 2.0 * baluster['window outline']['left']['mass'] * dotted_unit_wall_velocity ) / (structure[mote]['mass'] + baluster['window outline']['left']['mass'])
			var final_normal_wall_velocity = (dotted_unit_wall_velocity * (baluster['window outline']['left']['mass'] - structure[mote]['mass']) + 2.0 * baluster['window outline']['left']['mass'] * dotted_unit_mote_velocity ) / (structure[mote]['mass'] + baluster['window outline']['left']['mass'])
							
			#print(final_normal_mote_velocity,' final_normal_mote_velocity')
			#print(final_normal_wall_velocity,' final_normal_wall_velocity')
							
							
			structure[mote]['velocity'] = (final_normal_mote_velocity * unit_vector) + (final_tangential_mote_velocity * unit_tangent)
			#final_tangential_mote_velocity * unit_tangent
							
			#(final_normal_wall_velocity * unit_vector) + (final_tangential_wall_velocity * unit_tangent)
							
				#
			#print(structure[mote]['velocity']," structure[mote]['velocity']")
		
		
		elif collision_restitution < 1.0 and collision_restitution > 0.0:
			line_of_impact = wall_center
			incoming_angle = rad_to_deg(mote.get_position().angle_to(line_of_impact))
			outgoing_angle = rad_to_deg(atan(baluster['window outline']['left']['coefficient of restitution'] * incoming_angle ))
			x_leg_of_particle_velocity = structure[mote]['velocity'].y * tan(outgoing_angle)
			# determine which friction static or kinetic...
			#if structure[particle]['velocity'].x in range(-1,2) and structure[particle]['velocity'].y in range(-1,2):
			if structure[mote]['velocity'].x in range(-cell_size,cell_size+1) and structure[mote]['velocity'].y in range(-cell_size,cell_size+1):
				### using of static friction...
				x_component = x_leg_of_particle_velocity * static_friction
			else:
				### use of kinetic friction...
				x_component = x_leg_of_particle_velocity * kinetic_friction
							
			structure[mote]['velocity'] = structure[mote]['velocity'] - Vector2(x_component,structure[mote]['velocity'].y)
						
			
		elif collision_restitution == 0.0:
			final_velocity = ( (structure[mote]['mass'] * structure[mote]['velocity']) + (baluster['window outline']['left']['mass'] * baluster['window outline']['left']['velocity']) ) / (structure[mote]['mass'] * baluster['window outline']['left']['mass']) 
						
			structure[mote]['velocity'] = structure[mote]['velocity'] - final_velocity
						
	
	
	
	return structure[mote]['velocity']
	
	
	
	
	
	
	
	
	
func Collision_between_Other_Particles(artifact,other_artifact,grid_field,cell_size):
	### Collision between the particles....
	
	collision_restitution = artifact.coefficient_of_restitution * other_artifact.coefficient_of_restitution

	if collision_restitution >= 1.0 :
		### the collision is perfect elastic...
		elastic_term_1_of_vector_1 = 2.0 * grid_field[artifact]['mass']
		elastic_term_1_of_vector_2 = 2.0 * grid_field[other_artifact]['mass']
						
		elastic_mass_addition = grid_field[other_artifact]['mass'] + grid_field[artifact]['mass']
						
		elastic_position_kernel_of_vector_1 = (artifact.get_position() - other_artifact.get_position()) / cell_size
		elastic_position_kernel_of_vector_2 = (other_artifact.get_position() - artifact.get_position()) / cell_size
		### find final-velocity 1::
		elastic_velocity_kernel_of_vector_1 = grid_field[other_artifact]['velocity'] - grid_field[artifact]['velocity']
		elastic_dot_vector_1 = elastic_velocity_kernel_of_vector_1.dot(elastic_position_kernel_of_vector_1)
		### find final-velocit 2::
		elastic_velocity_kernel_of_vector_2  = grid_field[artifact]['velocity'] - grid_field[other_artifact]['velocity']
		elastic_dot_vector_2 = elastic_velocity_kernel_of_vector_2.dot(elastic_position_kernel_of_vector_2)
					
		elastic_magnitude_kernel_of_vector_1 = sqrt( (pow(elastic_position_kernel_of_vector_1.x,2.0)+pow(elastic_position_kernel_of_vector_1.y,2.0)) )
		elastic_magnitude_kernel_of_vector_2 = sqrt( (pow(elastic_position_kernel_of_vector_2.x,2.0)+pow(elastic_position_kernel_of_vector_2.y,2.0)) )
			
		### the velocity of the artifact...
		grid_field[artifact]['velocity'] = grid_field[artifact]['velocity'] - ( elastic_term_1_of_vector_2 / elastic_mass_addition ) * ( elastic_dot_vector_2 / elastic_magnitude_kernel_of_vector_2 ) * ( elastic_position_kernel_of_vector_2 )
		### the velocitiy of the other artifact...
		#grid_field[other_artifact]['velocity'] = grid_field[other_artifact]['velocity'] - ( elastic_term_1_of_vector_1 / elastic_mass_addition ) * ( elastic_dot_vector_1 / elastic_magnitude_kernel_of_vector_1 ) * ( elastic_position_kernel_of_vector_1 )
						
	elif collision_restitution < 1.0 and collision_restitution > 0.0:
		### collision is impartial inelastic...
		#"""
		# find the point of contact between particles...
		##print(artifact.relation_to_domain[other_artifact],' distance between particles.')
		### artifact....
		#line_of_impact = artifact.get_position() + (artifact.relation_to_domain[other_artifact] / 2.0)
		line_of_impact = (artifact.surrounding_area.intersection(other_artifact.surrounding_area)).get_center()
		##print(line_of_impact,' artifact to other artifact')
		incoming_angle_of_artifact = rad_to_deg(artifact.get_position().angle_to(line_of_impact))
		##print(incoming_angle_of_artifact,' incoming angle artifact to other artifact')
		outgoing_angle_of_artifact = rad_to_deg(atan(other_artifact.coefficient_of_restitution * incoming_angle_of_artifact ))
		x_leg_of_updated_artifact_velocity = snapped((grid_field[artifact]['velocity'].y * tan(outgoing_angle_of_artifact)),.01)
		# determine which friction static or kinetic...
		#if grid_field[artifact]['velocity'].x in range(-1,2) and grid_field[artifact]['velocity'].y in range(-1,2):
		if grid_field[artifact]['velocity'].x in range(-cell_size,cell_size+1) and grid_field[artifact]['velocity'].y in range(-cell_size,cell_size+1):
			### using of static friction...
			updated_artifact_x_component = x_leg_of_updated_artifact_velocity * other_artifact.coefficient_of_static_friction
		else:
			### use of kinetic friction...
			updated_artifact_x_component = x_leg_of_updated_artifact_velocity * other_artifact.coefficient_of_kinetic_friction
							
		grid_field[artifact]['velocity'] = grid_field[artifact]['velocity'] - Vector2(updated_artifact_x_component,grid_field[artifact]['velocity'].y)
		#"""
		
		#"""
		### the other artifact...
		##print(other_artifact,' ',grid_field[other_artifact]['velocity'],' checking y value')
		#line_of_impact = other_artifact.get_position() - (other_artifact.relation_to_domain[artifact] / 2.0)
		##print(line_of_impact,' other artifact to artifact')
		incoming_angle_of_other_artifact = rad_to_deg(other_artifact.get_position().angle_to(line_of_impact))
		outgoing_angle_of_other_artifact = rad_to_deg(atan(artifact.coefficient_of_restitution * incoming_angle_of_other_artifact ))
		x_leg_of_updated_other_artifact_velocity = snapped((grid_field[other_artifact]['velocity'].y * tan(outgoing_angle_of_other_artifact)),.01)
		# determine which friction static or kinetic...
		#if grid_field[other_artifact]['velocity'].x in range(-1,2) and grid_field[other_artifact]['velocity'].y in range(-1,2):
		if grid_field[other_artifact]['velocity'].x in range(-cell_size,cell_size+1) and grid_field[other_artifact]['velocity'].y in range(-cell_size,cell_size+1):
			### using of static friction...
			updated_other_artifact_x_component = x_leg_of_updated_other_artifact_velocity * artifact.coefficient_of_static_friction
		else:
			### use of kinetic friction...
			updated_other_artifact_x_component = x_leg_of_updated_other_artifact_velocity * artifact.coefficient_of_kinetic_friction
			
		grid_field[other_artifact]['velocity'] = grid_field[other_artifact]['velocity'] - Vector2(updated_other_artifact_x_component,grid_field[other_artifact]['velocity'].y)
		#"""
	elif collision_restitution == 0.0:
		### Perfect Inelastic Collisions:::
		# objects will stick together...
		
		final_velocity = ( (grid_field[artifact]['mass'] * grid_field[artifact]['velocity']) + (grid_field[other_artifact]['mass'] *grid_field[other_artifact]['velocity']) ) / (grid_field[artifact]['mass'] * grid_field[other_artifact]['mass']) 
						
		grid_field[artifact]['velocity'] = grid_field[artifact]['velocity'] - final_velocity
		
		#grid_field[other_artifact]['velocity'] = grid_field[other_artifact]['velocity'] - final_velocity
	
	return grid_field[artifact]['velocity']
