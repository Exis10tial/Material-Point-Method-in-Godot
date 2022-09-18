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




func Collision_with_Walls(substance,baluster,structure,cell_size):
	###...
	
	#"""
	for particle in substance:
		#print(particle,' particle')
		
		#print(particle,' ',structure[particle]['velocity'],' before collision')
		#print(" ")
		for wall in baluster['window outline'].keys():
			# Determine if the particle breaches the a wall...
			
			
			### suppose to be a wall where the particle hits...
			
			if wall == 'top':
				baluster['window outline'][wall]['velocity'] = particle.velocity * -100.0
				baluster['window outline'][wall]['outline'] = 0 + (cell_size)# / 2.0)
			elif wall == 'right':
				baluster['window outline'][wall]['velocity'] = particle.velocity * -100.0
				baluster['window outline'][wall]['outline'] = ProjectSettings.get_setting('display/window/size/width') - (cell_size)# / 2.0)
			elif wall == 'bottom':
				baluster['window outline'][wall]['velocity'] = particle.velocity * -100.0
				baluster['window outline'][wall]['outline'] = ProjectSettings.get_setting('display/window/size/height') - (cell_size)# / 2.0)
			elif wall == 'left':
				baluster['window outline'][wall]['velocity'] = particle.velocity * -100.0
				baluster['window outline'][wall]['outline'] = 0 + (cell_size)# / 2.0)
				
			if wall == 'top' and particle.surrounding_area.position.y <= baluster['window outline']["top"]['outline']:
				if typeof(baluster['window outline'][wall]['coefficient of restitution']) == TYPE_FLOAT:
					###the wall impacts the particles (elastic explosion,perfect elastic,partial inelastic, total inelastic...
					
					collision_restitution = baluster['window outline'][wall]['coefficient of restitution'] * particle.coefficient_of_restitution
					
					wall_center = Vector2(int(particle.position.x),0.0)
					#print('contact top')
					if collision_restitution >= 1.0 :
						### the collision is perfect elastic...
						
						elastic_term_1_of_vector_1 = 2.0 * baluster['window outline'][wall]['mass']
						elastic_term_1_of_vector_2 = 2.0 * structure[particle]['mass']
						
						elastic_mass_addition = structure[particle]['mass'] + baluster['window outline'][wall]['mass']
						
						elastic_position_kernel_of_vector_1 = (particle.get_position() - wall_center) / cell_size
						elastic_position_kernel_of_vector_2 = (wall_center - particle.get_position()) / cell_size
						### find final-velocity 1::
						elastic_velocity_kernel_of_vector_1 = structure[particle]['velocity'] - baluster['window outline'][wall]['velocity']
						elastic_dot_vector_1 = elastic_velocity_kernel_of_vector_1.dot(elastic_position_kernel_of_vector_1)
						
						### find final-velocit 2::
						elastic_velocity_kernel_of_vector_2  = baluster['window outline'][wall]['velocity'] - structure[particle]['velocity']
						elastic_dot_vector_2 = elastic_velocity_kernel_of_vector_2.dot(elastic_position_kernel_of_vector_2)
					
						elastic_magnitude_kernel_of_vector_1 = sqrt( (pow(elastic_position_kernel_of_vector_1.x,2.0)+pow(elastic_position_kernel_of_vector_1.y,2.0)) )
						elastic_magnitude_kernel_of_vector_2 = sqrt( (pow(elastic_position_kernel_of_vector_2.x,2.0)+pow(elastic_position_kernel_of_vector_2.y,2.0)) )
						
						structure[particle]['velocity'] = structure[particle]['velocity'] - ( elastic_term_1_of_vector_1 / elastic_mass_addition ) * ( elastic_dot_vector_1 / elastic_magnitude_kernel_of_vector_1 ) * ( elastic_position_kernel_of_vector_1 )
						
						return structure[particle]['velocity']
						
					elif collision_restitution < 1.0 and collision_restitution > 0.0:
						### ...
						
						line_of_impact = wall_center
						incoming_angle = rad_to_deg(particle.get_position().angle_to(line_of_impact))
						outgoing_angle = rad_to_deg(atan(baluster['window outline'][wall]['coefficient of restitution'] * incoming_angle ))
						x_leg_of_particle_velocity = structure[particle]['velocity'].y * tan(outgoing_angle)
						# determine which friction static or kinetic...
						#if structure[particle]['velocity'].x in range(-1,2) and structure[particle]['velocity'].y in range(-1,2):
						if structure[particle]['velocity'].x in range(-cell_size,cell_size+1) and structure[particle]['velocity'].y in range(-cell_size,cell_size+1):
							### using of static friction...
							x_component = x_leg_of_particle_velocity * static_friction
						else:
							### use of kinetic friction...
							x_component = x_leg_of_particle_velocity * kinetic_friction
							
						structure[particle]['velocity'] = structure[particle]['velocity'] - Vector2(x_component,structure[particle]['velocity'].y)
						
						return structure[particle]['velocity']
						
					elif collision_restitution == 0.0:
						### Perfect Inelastic Collisions:::
						# objects will stick together...
						
						### final_velocity = (Mass_1_i * Velocity_1_i) + (Mass_2_i * Velocity_2_i)) / (Mass_1 + Mass_2)
						
						final_velocity = ( (structure[particle]['mass'] * structure[particle]['velocity']) + (baluster['window outline'][wall]['mass'] * baluster['window outline'][wall]['velocity']) ) / (structure[particle]['mass'] * baluster['window outline'][wall]['mass']) 
						
						structure[particle]['velocity'] = structure[particle]['velocity'] - final_velocity
						
						return structure[particle]['velocity']
						
				else:
					### establish that there is no wall...### So two different aspect...
					### None or Disappear...
					if baluster['window outline'][wall]['coefficient of restitution'] == 'None':
						### nothing happens no adjustments are made to the particle...
						pass
					elif baluster['window outline'][wall]['coefficient of restitution'] == 'Disappear':
						### once the particle breaches the wall it will disappear from the simulation...
						get_tree().get_root().get_node("Test Area/Simulation").adjust_particles = true
						particle.to_remove = true
				
			elif wall == 'right' and particle.surrounding_area.end.x >= baluster['window outline']['right']['outline']:
				if typeof(baluster['window outline'][wall]['coefficient of restitution']) == TYPE_FLOAT:
					###the wall impacts the particles (pefect elastic,elastic,partial inelastic, perfect inelastic...
					
					collision_restitution = baluster['window outline'][wall]['coefficient of restitution'] * particle.coefficient_of_restitution
					
					wall_center = Vector2(baluster['window outline'][wall]['outline'],int(particle.position.y))
					#print('contact right')
					if collision_restitution >= 1.0 :
						### the collision is perfect elastic...
						
						elastic_term_1_of_vector_1 = 2.0 * baluster['window outline'][wall]['mass']
						elastic_term_1_of_vector_2 = 2.0 * structure[particle]['mass']
						
						elastic_mass_addition = structure[particle]['mass'] + baluster['window outline'][wall]['mass']
						
						elastic_position_kernel_of_vector_1 = (particle.get_position() - wall_center) / cell_size
						elastic_position_kernel_of_vector_2 = (wall_center - particle.get_position()) / cell_size
						### find final-velocity 1::
						elastic_velocity_kernel_of_vector_1 = structure[particle]['velocity'] - baluster['window outline'][wall]['velocity']
						elastic_dot_vector_1 = elastic_velocity_kernel_of_vector_1.dot(elastic_position_kernel_of_vector_1)
						### find final-velocit 2::
						elastic_velocity_kernel_of_vector_2  = baluster['window outline'][wall]['velocity'] - structure[particle]['velocity']
						elastic_dot_vector_2 = elastic_velocity_kernel_of_vector_2.dot(elastic_position_kernel_of_vector_2)
					
						elastic_magnitude_kernel_of_vector_1 = sqrt( (pow(elastic_position_kernel_of_vector_1.x,2.0)+pow(elastic_position_kernel_of_vector_1.y,2.0)) )
						elastic_magnitude_kernel_of_vector_2 = sqrt( (pow(elastic_position_kernel_of_vector_2.x,2.0)+pow(elastic_position_kernel_of_vector_2.y,2.0)) )
						
						
						structure[particle]['velocity'] = structure[particle]['velocity'] - ( ( elastic_term_1_of_vector_1 / elastic_mass_addition ) * ( elastic_dot_vector_1 / elastic_magnitude_kernel_of_vector_1 ) * ( elastic_position_kernel_of_vector_1 ) )
						
						return structure[particle]['velocity']
						
					elif collision_restitution < 1.0 and collision_restitution > 0.0:
						### ...
						line_of_impact = wall_center
						incoming_angle = rad_to_deg(particle.get_position().angle_to(line_of_impact))
						outgoing_angle = rad_to_deg(atan(baluster['window outline'][wall]['coefficient of restitution'] * incoming_angle ))
						x_leg_of_particle_velocity = snapped((structure[particle]['velocity'].y * tan(outgoing_angle)),.01)
						# determine which friction static or kinetic...
						#if structure[particle]['velocity'].x in range(-1,2) and structure[particle]['velocity'].y in range(-1,2):
						#if structure[particle]['velocity'].x in range(-cell_size,cell_size+1) and structure[particle]['velocity'].y in range(-cell_size,cell_size+1):
						if structure[particle]['velocity'].x in range(-cell_size,cell_size+1) and structure[particle]['velocity'].y in range(-cell_size,cell_size+1):
							### using of static friction...
							x_component = x_leg_of_particle_velocity * baluster['window outline'][wall]['coefficient of static friction']
						else:
							### use of kinetic friction...
							x_component = x_leg_of_particle_velocity * baluster['window outline'][wall]['coefficient of kinetic friction']
							
						structure[particle]['velocity'] = structure[particle]['velocity'] - Vector2(x_component,structure[particle]['velocity'].y)
						
						return structure[particle]['velocity']
						
					elif collision_restitution == 0.0:
						### Perfect Inelastic Collisions:::
						# objects will stick together...
						
						### final_velocity = (Mass_1_i * Velocity_1_i) + (Mass_2_i * Velocity_2_i)) / (Mass_1 + Mass_2)
						
						final_velocity = ( (structure[particle]['mass'] * structure[particle]['velocity']) + (baluster['window outline'][wall]['mass'] * baluster['window outline'][wall]['velocity']) ) / (structure[particle]['mass'] * baluster['window outline'][wall]['mass']) 
						
						structure[particle]['velocity'] = structure[particle]['velocity'] - final_velocity
						
						return structure[particle]['velocity']
						
				else:
					### establish that there is no wall...### So two different aspect...
					### None or Disappear...
					if baluster['window outline'][wall]['coefficient of restitution'] == 'None':
						### nothing happens no adjustments are made to the particle...
						pass
					elif baluster['window outline'][wall]['coefficient of restitution'] == 'Disappear':
						### once the particle breaches the wall it will disappear from the simulation...
						get_tree().get_root().get_node("Test Area/Simulation").adjust_particles = true
						particle.to_remove = true
				
				
			elif wall == 'bottom' and particle.surrounding_area.end.y >= baluster['window outline']['bottom']['outline']:
				if typeof(baluster['window outline'][wall]['coefficient of restitution']) == TYPE_FLOAT:
					###the wall impacts the particles (pefect elastic,elastic,partial inelastic, perfect inelastic...
					
					collision_restitution = baluster['window outline'][wall]['coefficient of restitution'] * particle.coefficient_of_restitution
					
					wall_center = Vector2(int(particle.position.x),baluster['window outline'][wall]['outline'])
					#print(wall_center,' wall center check')
					#print('contact bottom')
					if collision_restitution >= 1.0 :
						### the collision is perfect elastic...
						
						elastic_term_1_of_vector_1 = 2.0 * baluster['window outline'][wall]['mass']
						elastic_term_1_of_vector_2 = 2.0 * structure[particle]['mass']
						
						elastic_mass_addition = structure[particle]['mass'] + baluster['window outline'][wall]['mass']
						
						elastic_position_kernel_of_vector_1 = (particle.get_position() - wall_center) / cell_size
						elastic_position_kernel_of_vector_2 = (wall_center - particle.get_position()) / cell_size
						### find final-velocity 1::
						elastic_velocity_kernel_of_vector_1 = structure[particle]['velocity'] - baluster['window outline'][wall]['velocity']
						elastic_dot_vector_1 = elastic_velocity_kernel_of_vector_1.dot(elastic_position_kernel_of_vector_1)
						### find final-velocit 2::
						elastic_velocity_kernel_of_vector_2  = baluster['window outline'][wall]['velocity'] - structure[particle]['velocity']
						elastic_dot_vector_2 = elastic_velocity_kernel_of_vector_2.dot(elastic_position_kernel_of_vector_2)
					
						elastic_magnitude_kernel_of_vector_1 = sqrt( (pow(elastic_position_kernel_of_vector_1.x,2.0)+pow(elastic_position_kernel_of_vector_1.y,2.0)) )
						elastic_magnitude_kernel_of_vector_2 = sqrt( (pow(elastic_position_kernel_of_vector_2.x,2.0)+pow(elastic_position_kernel_of_vector_2.y,2.0)) )
						
						structure[particle]['velocity'] = structure[particle]['velocity'] - ( elastic_term_1_of_vector_1 / elastic_mass_addition ) * ( elastic_dot_vector_1 / elastic_magnitude_kernel_of_vector_1 ) * ( elastic_position_kernel_of_vector_1 )
						
						#print(structure[particle]['velocity'],' right after wall contact 1.0')
						
						return structure[particle]['velocity']
						
					elif collision_restitution < 1.0 and collision_restitution > 0.0:
						### impartial collisions...
						
						# determine which friction static or kinetic...
						#print(structure[particle]['velocity'],' velocity before collision')
						
						line_of_impact = wall_center
						incoming_angle = rad_to_deg(particle.get_position().angle_to(line_of_impact))
						outgoing_angle = rad_to_deg(atan(baluster['window outline'][wall]['coefficient of restitution'] * incoming_angle ))
						y_leg_of_particle_velocity = structure[particle]['velocity'].y
						x_leg_of_particle_velocity = y_leg_of_particle_velocity * tan(outgoing_angle)
						
						# determine which friction static or kinetic...
						# if the particle is moving within a certain range...
						#if structure[particle]['velocity'].x in range(-1.0,2.0) and structure[particle]['velocity'].y in range(-1.0,2.0):
						if structure[particle]['velocity'].x in range(-cell_size,cell_size+1) and structure[particle]['velocity'].y in range(-cell_size,cell_size+1):
							### using of static friction...
							x_component = x_leg_of_particle_velocity * baluster['window outline'][wall]['coefficient of static friction']
						else:
							### use of kinetic friction...
							x_component = x_leg_of_particle_velocity * baluster['window outline'][wall]['coefficient of kinetic friction']
						
						#print(Vector2(x_component,y_leg_of_particle_velocity),' check created velocity')
						structure[particle]['velocity'] = structure[particle]['velocity'] - Vector2(x_component,y_leg_of_particle_velocity)
						#print(structure[particle]['velocity'],' right after wall contact 0.5')
						
						return structure[particle]['velocity']
						
					elif collision_restitution == 0.0:
						### Perfect Inelastic Collisions:::
						# objects will stick together...
						
						### final_velocity = (Mass_1_i * Velocity_1_i) + (Mass_2_i * Velocity_2_i)) / (Mass_1 + Mass_2)
						
						final_velocity = ( (structure[particle]['mass'] * structure[particle]['velocity']) + (baluster['window outline'][wall]['mass'] * baluster['window outline'][wall]['velocity']) ) / (structure[particle]['mass'] * baluster['window outline'][wall]['mass']) 
						
						structure[particle]['velocity'] = structure[particle]['velocity'] - final_velocity
						
						return structure[particle]['velocity']
						
				else:
					### establish that there is no wall...### So two different aspect...
					### None or Disappear...
					if baluster['window outline'][wall]['coefficient of restitution'] == 'None':
						### nothing happens no adjustments are made to the particle...
						pass
					elif baluster['window outline'][wall]['coefficient of restitution'] == 'Disappear':
						### once the particle breaches the wall it will disappear from the simulation...
						get_tree().get_root().get_node("Test Area/Simulation").adjust_particles = true
						particle.to_remove = true
				
			elif wall == 'left' and particle.surrounding_area.position.x <= baluster['window outline']['left']['outline']:
				if typeof(baluster['window outline'][wall]['coefficient of restitution']) == TYPE_FLOAT:
					###the wall impacts the particles (pefect elastic,elastic,partial inelastic, perfect inelastic...

					collision_restitution = baluster['window outline'][wall]['coefficient of restitution'] * particle.coefficient_of_restitution
					
					wall_center = Vector2(0.0,int(particle.position.y))
					#print('contact left')
					if collision_restitution >= 1.0 :
						### the collision is perfect elastic...
						
						elastic_term_1_of_vector_1 = 2.0 * baluster['window outline'][wall]['mass']
						elastic_term_1_of_vector_2 = 2.0 * structure[particle]['mass']
						
						elastic_mass_addition = structure[particle]['mass'] + baluster['window outline'][wall]['mass']
						
						elastic_position_kernel_of_vector_1 = (particle.get_position() - wall_center) / cell_size
						elastic_position_kernel_of_vector_2 = (wall_center - particle.get_position()) / cell_size
						### find final-velocity 1::
						elastic_velocity_kernel_of_vector_1 = structure[particle]['velocity'] - baluster['window outline'][wall]['velocity']
						elastic_dot_vector_1 = elastic_velocity_kernel_of_vector_1.dot(elastic_position_kernel_of_vector_1)
						### find final-velocit 2::
						elastic_velocity_kernel_of_vector_2  = baluster['window outline'][wall]['velocity'] - structure[particle]['velocity']
						elastic_dot_vector_2 = elastic_velocity_kernel_of_vector_2.dot(elastic_position_kernel_of_vector_2)
					
						elastic_magnitude_kernel_of_vector_1 = sqrt( (pow(elastic_position_kernel_of_vector_1.x,2.0)+pow(elastic_position_kernel_of_vector_1.y,2.0)) )
						elastic_magnitude_kernel_of_vector_2 = sqrt( (pow(elastic_position_kernel_of_vector_2.x,2.0)+pow(elastic_position_kernel_of_vector_2.y,2.0)) )
						
						structure[particle]['velocity'] = structure[particle]['velocity'] - ( elastic_term_1_of_vector_1 / elastic_mass_addition ) * ( elastic_dot_vector_1 / elastic_magnitude_kernel_of_vector_1 ) * ( elastic_position_kernel_of_vector_1 )
						
						return structure[particle]['velocity']
						
					elif collision_restitution < 1.0 and collision_restitution > 0.0:
						### ...
						line_of_impact = wall_center
						incoming_angle = rad_to_deg(particle.get_position().angle_to(line_of_impact))
						outgoing_angle = rad_to_deg(atan(baluster['window outline'][wall]['coefficient of restitution'] * incoming_angle ))
						x_leg_of_particle_velocity = structure[particle]['velocity'].y * tan(outgoing_angle)
						# determine which friction static or kinetic...
						#if structure[particle]['velocity'].x in range(-1,2) and structure[particle]['velocity'].y in range(-1,2):
						if structure[particle]['velocity'].x in range(-cell_size,cell_size+1) and structure[particle]['velocity'].y in range(-cell_size,cell_size+1):
							### using of static friction...
							x_component = x_leg_of_particle_velocity * baluster['window outline'][wall]['coefficient of static friction']
						else:
							### use of kinetic friction...
							x_component = x_leg_of_particle_velocity * baluster['window outline'][wall]['coefficient of kinetic friction']
							
						structure[particle]['velocity'] = structure[particle]['velocity'] - Vector2(x_component,structure[particle]['velocity'].y)
						
						return structure[particle]['velocity']
						
					elif collision_restitution == 0.0:
						### Perfect Inelastic Collisions:::
						# objects will stick together...
						
						### final_velocity = (Mass_1_i * Velocity_1_i) + (Mass_2_i * Velocity_2_i)) / (Mass_1 + Mass_2)
						
						final_velocity = ( (structure[particle]['mass'] * structure[particle]['velocity']) + (baluster['window outline'][wall]['mass'] * baluster['window outline'][wall]['velocity']) ) / (structure[particle]['mass'] * baluster['window outline'][wall]['mass']) 
						
						structure[particle]['velocity'] = structure[particle]['velocity'] - final_velocity
						
						return structure[particle]['velocity']
						
				else:
					### establish that there is no wall...### So two different aspect...
					### None or Disappear...
					if baluster['window outline'][wall]['coefficient of restitution'] == 'None':
						### nothing happens no adjustments are made to the particle...
						pass
					elif baluster['window outline'][wall]['coefficient of restitution'] == 'Disappear':
						### once the particle breaches the wall it will disappear from the simulation...
						get_tree().get_root().get_node("Test Area/Simulation").adjust_particles = true
						particle.to_remove = true
				
		
	#"""
	
	
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
		#print(artifact.relation_to_domain[other_artifact],' distance between particles.')
		### artifact....
		#line_of_impact = artifact.get_position() + (artifact.relation_to_domain[other_artifact] / 2.0)
		line_of_impact = (artifact.surrounding_area.intersection(other_artifact.surrounding_area)).get_center()
		#print(line_of_impact,' artifact to other artifact')
		incoming_angle_of_artifact = rad_to_deg(artifact.get_position().angle_to(line_of_impact))
		#print(incoming_angle_of_artifact,' incoming angle artifact to other artifact')
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
		#print(other_artifact,' ',grid_field[other_artifact]['velocity'],' checking y value')
		#line_of_impact = other_artifact.get_position() - (other_artifact.relation_to_domain[artifact] / 2.0)
		#print(line_of_impact,' other artifact to artifact')
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
