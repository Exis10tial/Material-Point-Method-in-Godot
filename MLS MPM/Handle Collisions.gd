extends Node



func _on_handle_collisions_ready():
	pass # Replace with function body.




func Collision_with_Walls(substance,baluster,structure,center_of_window,cell_size):
	###...
	#print(" ")
	#print('Collision Check')
	#"""
	for particle in substance:
		#print(particle,' particle')
		
		#print(particle,' ',structure[particle]['velocity'],' before collision')
		#print(" ")
		for wall in baluster['window outline'].keys():
			
			# Determine if the particle breaches the a wall...
			if wall == 'top' and particle.surrounding_area.position.y <= baluster['window outline'][wall]['outline']:
				print('breached top')
				if typeof(baluster['window outline'][wall]['coefficient of reinstition']) == TYPE_FLOAT:
					###the wall impacts the particles (elastic explosion,perfect elastic,partial inelastic, total inelastic...
					if baluster['window outline'][wall]['coefficient of reinstition'] == 1.0:
						### the collision is perfect elastic...
						pass
						
					elif baluster['window outline'][wall]['coefficient of reinstition'] == 0.5:
						pass
					elif baluster['window outline'][wall]['coefficient of reinstition'] == 0.0:
						pass
					
				else:
					### establish that there is no wall...### So two different aspect...
					### None or Disappear...
					if baluster['window outline'][wall]['coefficient of reinstition'] == 'None':
						### nothing happens no adjustments are made to the particle...
						pass
					elif baluster['window outline'][wall]['coefficient of reinstition'] == 'Disappear':
						### once the particle breaches the wall it will disappear from the simulation...
						get_tree().get_root().get_node("Test Area/Simulation").adjust_particles = true
						particle.to_remove = true
				
			elif wall == 'right' and particle.surrounding_area.end.x >= baluster['window outline'][wall]['outline']:
				print('breached right')
				if typeof(baluster['window outline'][wall]['coefficient of reinstition']) == TYPE_FLOAT:
					###the wall impacts the particles (pefect elastic,elastic,partial inelastic, perfect inelastic...
					pass
				else:
					### establish that there is no wall...### So two different aspect...
					### None or Disappear...
					if baluster['window outline'][wall]['coefficient of reinstition'] == 'None':
						### nothing happens no adjustments are made to the particle...
						pass
					elif baluster['window outline'][wall]['coefficient of reinstition'] == 'Disappear':
						### once the particle breaches the wall it will disappear from the simulation...
						get_tree().get_root().get_node("Test Area/Simulation").adjust_particles = true
						particle.to_remove = true
				
				
			elif wall == 'bottom' and particle.surrounding_area.end.y >= baluster['window outline'][wall]['outline']:
				#print(particle,' breached bottom')
				if typeof(baluster['window outline'][wall]['coefficient of reinstition']) == TYPE_FLOAT:
					###the wall impacts the particles (pefect elastic,elastic,partial inelastic, perfect inelastic...
					if baluster['window outline'][wall]['coefficient of reinstition'] == 1.0:
						### the collision is perfect elastic...
						pass
						## angle-free :
						# velocity1final = velocity1initial - ( ((2*mass2) * dotproduct(velocitykernel,positionkernel) * positionkernel ) / ( (mass1+mass2) * (|| posiitonkernel ||^2) )
						# for velocity1final notes :: velocitykernel = v1 - v2,positionkernel = x1-x2,|| || = magnitude = squareroot(x^2+y^2)
						
						# velocity2final = velocity2initial - ( ((2*mass1) * dotproduct(velocitykernel,positionkernel) * positionkernel ) / ( (mass1+mass2) * (|| positionkernel ||^2) ) 
						#for velocity2final notes :: velocitykernel = v2 - v1,positionkernel = x2 - x1,|| = magnitude = squareroot(x^2+y^2)
						
						### final velcotiy 1 = particle
						### final velcoity 2 = wall
						var elastic_term_1_of_vector_1 = 2 * baluster['window outline'][wall]['mass']
						var elastic_term_1_of_vector_2 = 2 * structure[particle]['mass']
						
						var elastic_mass_addition = structure[particle]['mass'] + baluster['window outline'][wall]['mass']
						var elastic_dot_vector_1
						var elastic_dot_vector_2
						var elastic_position_kernel_of_vector_1
						var elastic_position_kernel_of_vector_2
						var elastic_velocity_kernel_of_vector_1
						var elastic_velocity_kernel_of_vector_2
						var elastic_magnitude_kernel_of_vector_1
						var elastic_magnitude_kernel_of_vector_2
						
						### wall center is along "y" to the object 
						#var wall_center = Vector2(center_of_window.x,ProjectSettings.get_setting('display/window/size/height'))
						var wall_center = Vector2(particle.position.x,particle.position.y+cell_size)
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
						
						var manual_dot_vector = (elastic_velocity_kernel_of_vector_1.x * elastic_position_kernel_of_vector_1.x) + (elastic_velocity_kernel_of_vector_1.y * elastic_position_kernel_of_vector_1.y)
						#print(wall_center,' wall center')
						#print(elastic_term_1_of_vector_1,' 2 * wall mass')
						#print(elastic_dot_vector_1,' dot product thru engine')
						#print(manual_dot_vector,' dot product thru manual')
						#print(elastic_position_kernel_of_vector_1,' position kernel')
						#print(elastic_mass_addition,' particle mass + wall mass')
						#print(elastic_magnitude_kernel_of_vector_1,' magnitude of kernel')
						structure[particle]['velocity'] = structure[particle]['velocity'] - ( elastic_term_1_of_vector_1 / elastic_mass_addition ) * ( elastic_dot_vector_1 / elastic_magnitude_kernel_of_vector_1 ) * ( elastic_position_kernel_of_vector_1 )
						
						#print(structure[particle]['velocity'],' after collision')
						
						#break
					elif baluster['window outline'][wall]['coefficient of reinstition'] == 0.5:
						### requires 4 seperate equations...
						### line of impact = line between object's masses
						### plane contact = perpendicular to the line of impact
						### i = initial ; f = final ; line of impact = x ; plane contact = y
						### e = coefficient of reinstition
						#
						#1.(mass_1_i * velocity_1_x_i) cos(angle) + (mass_2_i * velocity_2_x_i) cos(angle) =  (mass_1_f * velocity_1_x_f) cos(angle) + (mass_2_f * velocity_2_x_f) cos(angle)
						#
						#2. velocity_2_f - velocity_1_f = e * ( velocity_2_i - velocity_1_i )
						#
						# 3. velocity_y_2_f = velocity_y_2_i
						# 4. velocity_y_1_f = velocity_y_1_i
						#
						######
						#####
						####
						# 1. find angle_i between ( velocity_1_initial : line of impact )
						# .... cos(angle_i) = vector_a.vector_b / |vector_a|*|vector_b|
						#         ... 1. dot_product between vectors
						#         ... 2. calculate the length of the vectors then multiply them together.
						#         ... 3. dot_product results divided by multiplied lengths results.
						#         ... 4. angle_i = arc_cos( results from dot_product / Multipled Lengths)
						# 
						# 2. tan(angle_f) = e * tan(angle_i)
						#    .... angle_f = arc_tan * e * tan(angle_i)
						#
						#3. other leg : velocity_y_2_f * tan(angle_f)
						#
						#4. pythagorean theorem : a^2 + b^2 = c^2
						#
						pass
						
						
						
						
						
					elif baluster['window outline'][wall]['coefficient of reinstition'] == 0.0:
						### Perfect Inelastic Collisions:::
						# objects will stick together...
						
						### final_velocity = square_root((Mass_1_i * Velocity_1_i) + (Mass_2_i * Velocity_2_i)) / (Mass_1 + Mass_2)
						pass
					
				else:
					### establish that there is no wall...### So two different aspect...
					### None or Disappear...
					if baluster['window outline'][wall]['coefficient of reinstition'] == 'None':
						### nothing happens no adjustments are made to the particle...
						pass
					elif baluster['window outline'][wall]['coefficient of reinstition'] == 'Disappear':
						### once the particle breaches the wall it will disappear from the simulation...
						get_tree().get_root().get_node("Test Area/Simulation").adjust_particles = true
						particle.to_remove = true
				
			elif wall == 'left' and particle.surrounding_area.position.x <= baluster['window outline'][wall]['outline']:
				print('breached left')
				if typeof(baluster['window outline'][wall]['coefficient of reinstition']) == TYPE_FLOAT:
					###the wall impacts the particles (pefect elastic,elastic,partial inelastic, perfect inelastic...
					pass
				else:
					### establish that there is no wall...### So two different aspect...
					### None or Disappear...
					if baluster['window outline'][wall]['coefficient of reinstition'] == 'None':
						### nothing happens no adjustments are made to the particle...
						pass
					elif baluster['window outline'][wall]['coefficient of reinstition'] == 'Disappear':
						### once the particle breaches the wall it will disappear from the simulation...
						get_tree().get_root().get_node("Test Area/Simulation").adjust_particles = true
						particle.to_remove = true
				
			
	#"""
	
	
func Collision_between_Other_Particles():
	###...
	pass
