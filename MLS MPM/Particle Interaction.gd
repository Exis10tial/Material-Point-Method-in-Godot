extends Node

#class_name ParticleInteraction

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


func Collision_with_Walls(breach,mote,particle_boundary,baluster,structure,cell_size):
	###...
	if breach == 'top':
		
		collision_restitution = baluster['window outline']['top']['coefficient of restitution'] * mote.coefficient_of_restitution
		collision_static_friction = baluster['window outline']['top']['coefficient of static friction'] * mote.coefficient_of_static_friction
		collision_kinetic_friction = baluster['window outline']['top']['coefficient of kinetic friction'] * mote.coefficient_of_kinetic_friction
		
		#baluster['window outline']['top']['velocity'] = Vector2(0.0,0.0)
		
		baluster['window outline']['top']['velocity'] = Vector2(0.0,10) 
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
			final_normal_mote_velocity = (dotted_unit_mote_velocity * (structure['mass'] - baluster['window outline']['top']['mass']) + 2.0 * baluster['window outline']['top']['mass'] * dotted_unit_wall_velocity ) / (structure['mass'] + baluster['window outline']['right']['mass'])
			final_normal_wall_velocity = (dotted_unit_wall_velocity * (baluster['window outline']['top']['mass'] - structure['mass']) + 2.0 * baluster['window outline']['top']['mass'] * dotted_unit_mote_velocity ) / (structure['mass'] + baluster['window outline']['right']['mass'])
			
			structure['velocity'] = (final_normal_mote_velocity * unit_vector) + (final_tangential_mote_velocity * unit_tangent)
			#final_tangential_mote_velocity * unit_tangent
			#(final_normal_wall_velocity * unit_vector) + (final_tangential_wall_velocity * unit_tangent)
			
			#print(structure['velocity']," structure['velocity']")
		
		elif collision_restitution < 1.0 and collision_restitution > 0.0:
			
			line_of_impact = wall_center
			
			incoming_angle = abs(snapped(rad_to_deg(line_of_impact.angle_to_point(particle_boundary.origin)),1))
			outgoing_angle = snapped(rad_to_deg(atan(collision_restitution * tan(deg_to_rad(incoming_angle)))),.1)                                                                                                                                                                                                                                                     
			#print(incoming_angle," incoming_angle")
			#print(outgoing_angle," outgoing_angle")
			if structure['velocity'].x in range(-cell_size,cell_size+1):
				x_component = structure['velocity'].x * baluster['window outline']['top']['coefficient of static friction']
			else:
				x_component = structure['velocity'].x * baluster['window outline']['top']['coefficient of kinetic friction']
			y_component = -collision_restitution *  structure['velocity'].y 
			#print(x_component," x_component")
			#print(y_component," y_component")
			#var component =sqrt( pow(x_component,2) + pow(y_component,2))
			#print(component," component")
			structure['velocity'] = Vector2(x_component,y_component)
			#print("top")			print(structure['velocity']," structure['velocity']")
			
		elif collision_restitution == 0.0:
			#perfect inelastic collisions : [m1 / (m1 + m2)] * vi = vf
			
			x_component = structure['mass'] / (structure['mass'] + baluster['window outline']['top']['mass']) * structure['velocity'].x
			y_component = structure['mass'] / (structure['mass'] + baluster['window outline']['top']['mass']) * structure['velocity'].y
			
			structure['velocity'] = Vector2(x_component,y_component)
			
	if breach == 'right':
		
		collision_restitution = baluster['window outline']['right']['coefficient of restitution'] * mote.coefficient_of_restitution
		collision_static_friction = baluster['window outline']['right']['coefficient of static friction'] * mote.coefficient_of_static_friction
		collision_kinetic_friction = baluster['window outline']['right']['coefficient of kinetic friction'] * mote.coefficient_of_kinetic_friction
		
		baluster['window outline']['right']['velocity'] = Vector2(-10.0,0.0)
		
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
			line_of_impact = wall_center
			
			incoming_angle = abs(snapped(rad_to_deg(line_of_impact.angle_to_point(particle_boundary.origin)),1))
			outgoing_angle = snapped(rad_to_deg(atan(collision_restitution * tan(deg_to_rad(incoming_angle)))),.1)                                                                                                                                                                                                                                                     
			#print(incoming_angle," incoming_angle")
			#print(outgoing_angle," outgoing_angle")
			if structure['velocity'].y in range(-cell_size,cell_size+1):
				y_component = structure['velocity'].y * baluster['window outline']['top']['coefficient of static friction']
			else:
				y_component = structure['velocity'].y * baluster['window outline']['top']['coefficient of kinetic friction']
			x_component =  -collision_restitution * structure['velocity'].x
				
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
		
		baluster['window outline']['bottom']['velocity'] = Vector2(0,-10)
		
		wall_center = Vector2(particle_boundary.origin.x,baluster['window outline']['bottom']['outline'])
		
		if collision_restitution >= 1.0 :
			### the collision is perfect elastic...
			
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
			line_of_impact = wall_center
			
			incoming_angle = abs(snapped(rad_to_deg(line_of_impact.angle_to_point(particle_boundary.origin)),1))
			outgoing_angle = snapped(rad_to_deg(atan(collision_restitution * tan(deg_to_rad(incoming_angle)))),.1)                                                                                                                                                                                                                                                     
			
			if structure['velocity'].x in range(-cell_size,cell_size+1):
				x_component = structure['velocity'].x * baluster['window outline']['top']['coefficient of static friction']
			else:
				x_component = structure['velocity'].x * baluster['window outline']['top']['coefficient of kinetic friction']
			y_component =  -collision_restitution * structure['velocity'].y
			
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
		
		baluster['window outline']['left']['velocity'] = Vector2(10.0,0.0)
		
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
		
			line_of_impact = wall_center
			
			incoming_angle = abs(snapped(rad_to_deg(line_of_impact.angle_to_point(particle_boundary.origin)),1))
			outgoing_angle = snapped(rad_to_deg(atan(collision_restitution * tan(deg_to_rad(incoming_angle)))),.1)                                                                                                                                                                                                                                                     
			#print(incoming_angle," incoming_angle")
			#print(outgoing_angle," outgoing_angle")
			if structure['velocity'].y in range(-cell_size,cell_size+1):
				y_component = structure['velocity'].y * baluster['window outline']['top']['coefficient of static friction']
			else:
				y_component = structure['velocity'].y * baluster['window outline']['top']['coefficient of kinetic friction']
			x_component = -collision_restitution * structure['velocity'].x
				
			structure['velocity'] = Vector2(x_component,y_component)
			
		elif collision_restitution == 0.0:
			#perfect inelastic collisions : [m1 / (m1 + m2)] * vi = vf
			
			x_component = structure['mass'] / (structure['mass'] + baluster['window outline']['left']['mass']) * structure['velocity'].x
			y_component = structure['mass'] / (structure['mass'] + baluster['window outline']['left']['mass']) * structure['velocity'].y
			
			structure['velocity'] = Vector2(x_component,y_component)
			#structure['velocity'] = final_velocity
	
	return structure['velocity']
	
	
	
	
#func Collision_between_Other_Particles(artifact_name,artifact,artifact_zone,other_artifact_name,other_artifact,other_artifact_zone,grid_field):
func Collision_between_Other_Particles(artifact_name,artifact,artifact_zone,other_artifact_name,other_artifact,other_artifact_zone,artifact_grid,other_artifact_grid):
	### Collision between the particles....
	
	collision_restitution = artifact.coefficient_of_restitution * other_artifact.coefficient_of_restitution
	#print(collision_restitution,' between the particles')
	if collision_restitution >= 1.0 :
		
		#normal_vector = other_artifact_zone.get_center() - artifact_zone.get_center()
		normal_vector = other_artifact_zone.origin - artifact_zone.origin
		
		unit_vector = normal_vector / snapped(sqrt((snapped(pow(normal_vector.x,2.0),.001) + snapped(pow(normal_vector.y,2.0),.001))),.001)
		#unit_vector = normal_vector / sqrt((pow(normal_vector.x,2.0) + pow(normal_vector.y,2.0)))
		unit_tangent = Vector2(-unit_vector.y,unit_vector.x)
		
		dotted_unit_artifact_velocity = unit_vector.dot(artifact_grid['velocity'])
		dotted_tangent_artifact_velocity = unit_tangent.dot(artifact_grid['velocity'])
		
		dotted_unit_other_artifact_velocity = unit_vector.dot(other_artifact_grid['velocity'])
		dotted_tangent_other_artifact_velocity = unit_tangent.dot(other_artifact_grid['velocity'])
		
		final_tangential_artifact_velocity = dotted_tangent_artifact_velocity
		#final_tangential_wall_velocity = dotted_tangent_other_artifact_velocity
		final_normal_artifact_velocity = ( (dotted_unit_artifact_velocity * (artifact_grid['mass'] - other_artifact_grid['mass']) )+( 2.0 * other_artifact_grid['mass'] * dotted_unit_other_artifact_velocity ) )/ (artifact_grid['mass'] + other_artifact_grid['mass'])
		
		#final_normal_other_artifact_velocity = (dotted_unit_other_artifact_velocity * (grid_field[other_artifact_name]['mass'] - grid_field[artifact_name]['mass']) + 2.0 * grid_field[other_artifact_name]['mass'] * dotted_unit_artifact_velocity ) / (grid_field[artifact_name]['mass'] + grid_field[other_artifact_name]['mass'])
							
		artifact_grid['velocity'] = (final_normal_artifact_velocity * unit_vector) + (final_tangential_artifact_velocity * unit_tangent)
		
		#final_tangential_mote_velocity * unit_tangent
		#(final_normal_wall_velocity * unit_vector) + (final_tangential_wall_velocity * unit_tangent)
		
		#print(grid_field[artifact_name]['velocity']," grid_field[artifact]['velocity'] right after")
		
	elif collision_restitution < 1.0 and collision_restitution > 0.0:
		### collision is impartial inelastic...
		#"""
		# find the point of contact between particles...
		#line_of_impact = (artifact_zone.intersection(other_artifact_zone)).get_center()
		#line_of_impact = artifact_zone.get_center() - other_artifact_zone.get_center()
		line_of_impact = artifact_zone.origin - other_artifact_zone.origin
		##print(line_of_impact,' artifact to other artifact')
		#incoming_angle_of_artifact = rad_to_deg(artifact_zone.get_center().angle_to(line_of_impact))
		incoming_angle_of_artifact = rad_to_deg(artifact_zone.origin.angle_to(line_of_impact))
		##print(incoming_angle_of_artifact,' incoming angle artifact to other artifact')
		outgoing_angle_of_artifact = rad_to_deg(atan(other_artifact.coefficient_of_restitution * incoming_angle_of_artifact ))
		x_leg_of_updated_artifact_velocity = snapped((artifact_grid['velocity'].y * tan(outgoing_angle_of_artifact)),.01)
		# determine which friction static or kinetic...
		#if grid_field[artifact]['velocity'].x in range(-1,2) and grid_field[artifact]['velocity'].y in range(-1,2):
		#"""
		if artifact_grid['velocity'].x in range(-artifact_zone.size.x,artifact_zone.size.x+1) and artifact_grid['velocity'].y in range(-artifact_zone.size.y,artifact_zone.size.y+1):
			### using of static friction...
			updated_artifact_x_component = x_leg_of_updated_artifact_velocity * other_artifact.coefficient_of_static_friction
			#print(updated_artifact_x_component,' coefficient_of_static_friction')
		else:
			### use of kinetic friction...
			updated_artifact_x_component = x_leg_of_updated_artifact_velocity * other_artifact.coefficient_of_kinetic_friction
			#print(updated_artifact_x_component,' coefficient_of_kinetic_friction')
		#"""
		#grid_field[artifact]['velocity'] = grid_field[artifact]['velocity'] + Vector2(updated_artifact_x_component,grid_field[artifact]['velocity'].y)
		artifact_grid['velocity'] = artifact_grid['velocity'] - Vector2(x_leg_of_updated_artifact_velocity,artifact_grid['velocity'].y)
		#"""
		#print(grid_field[artifact_name]['velocity']," grid_field[artifact]['velocity'] right after")
	elif collision_restitution == 0.0:
		### Perfect Inelastic Collisions:::
		# objects will stick together...
		
		final_velocity = ( (artifact_grid['mass'] * artifact_grid['velocity']) + (other_artifact_grid['mass'] *other_artifact_grid['velocity']) ) / (artifact_grid['mass'] * other_artifact_grid['mass']) 
						
		artifact_grid['velocity'] = final_velocity
		
		#grid_field[other_artifact]['velocity'] = grid_field[other_artifact]['velocity'] - final_velocity
		#print(grid_field[artifact_name]['velocity']," grid_field[artifact]['velocity'] right after")
	#print(' ')
	return artifact_grid['velocity']
