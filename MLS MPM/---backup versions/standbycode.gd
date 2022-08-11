



### particles merging with other ....
#get_tree().get_root().get_node("Test Area/Simulation/Particle Interaction").Particles_Merging()
						
### particles merging with other ....
				#get_tree().get_root().get_node("Test Area/Simulation/Particle Interaction").Particles_Merging()
					
"""
			if particle != other_particle:
				if particle.surrounding_area.intersects(other_particle.surrounding_area):
					### come into some contact...
				
					#print(particle,' at position ',particle.position,' in contact with ',other_particle,' at position ',other_particle.surrounding_area,' quadrant')
					if particle.within_range.has(other_particle):
						### the particle is already in contact with the other...
						pass
					else:
						### recognized to be within the domain range...
						particle.within_range.append(other_particle)
						
					### particles interaction with other particles...
					get_tree().get_root().get_node("Test Area/Simulation/Particle Interaction").Collision_between_Other_Particles(particle,other_particle,the_grid,grid_domain)
							
				else:
					### particle is not in contact with each other...
					particle.within_range.erase(other_particle)
					other_particle.within_range.erase(particle)
			else:
				### can't collide with itself...
				pass
			"""
if  contact_with.x == 0 and contact_with.y == 0:
				### always in contact with itself...
				pass
			elif contact_with.x > -(grid_domain/grid_domain) and contact_with.x < grid_domain/grid_domain and contact_with.y > -(grid_domain/grid_domain) and contact_with.y < grid_domain/grid_domain:
				### come into some contact...
				
				#print(particle,' at position ',particle.position,' in contact with ',other_particle,' at position ',other_particle.surrounding_area,' quadrant')
				if particle.within_range.has(other_particle):
					### the particle is already in contact with the other...
					pass
				else:
					### recognized to be within the domain range...
					particle.within_range.append(other_particle)
					
				### particles merging with other ....
				#get_tree().get_root().get_node("Test Area/Simulation/Particle Interaction").Particles_Merging()
					
				### particles interaction with other particles...
				get_tree().get_root().get_node("Test Area/Simulation/Particle Interaction").Collision_between_Other_Particles(particle,other_particle,the_grid,grid_domain)
						
			else:
				###
				#print(particle,' not in contact with this particle ',other_particle)
				particle.within_range.erase(other_particle)
				other_particle.within_range.erase(particle)
			#"""



#### from grid2particle
#var c_term_1 =  (the_grid['velocity'] * weight_interpolation * 4) / pow(grid_domain,2.0)
			#var c_cell_coefficent =  snapped((basis_coefficent / pow(grid_domain,2.0)),2.0)
			#var c_inverse_I_coefficent = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Inverse_Matrix(other_particle.I)
			#var c_term_1 = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Matrix_by_Scalar(c_inverse_I_coefficent,c_cell_coefficent,true)
			#var c_term_2 = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Matrix_by_Vector2_to_Matrix(c_term_1,the_grid[other_particle]['velocity'],false,false)
			#var c_term_3 = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Matrix_by_Vector2_to_Matrix(c_term_2,particle.domain_relation_to_particle[other_particle],true,false)
			#var sum_of_c = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Matrix_by_Scalar(c_term_3,weight_interpolation,true)
			#particle.C = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Add_Matrix(particle.C,sum_of_c)
			#c_flipped_velocity_coefficient = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Vector2_by_Vector2_to_Matrix(the_grid[other_particle]['velocity'],false,particle.domain_relation_to_particle[other_particle],true)
			#c_weighted_velocity_matrix
			#c_weighted_velocity_matrix = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Matrix_by_Scalar(c_flipped_velocity_coefficient,weight_interpolation,true)
			#sum_of_c_weighted_velocity_matrix = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Add_Matrix(sum_of_c_weighted_velocity_matrix,c_weighted_velocity_matrix)
			#particle.velocity = particle.velocity + (the_grid[other_particle]['velocity'] * weight_interpolation)
			#print(" ")
			#print('loop to find velocity ...')
			#print(the_grid[other_particle]['velocity'],' other particle velocity')
			#print(weight_interpolation,' wip')
			#particle.velocity.x = snapped((particle.velocity.x + snapped((the_grid[other_particle]['velocity'].x * weight_interpolation),.01)),.01)
			#particle.velocity.y = snapped((particle.velocity.y + snapped((the_grid[other_particle]['velocity'].y * weight_interpolation),.01)),.01)
			
			#print(particle.velocity,' velocity')
			### Once the last item of the list is processed...
			#if particle.within_range.back() == other_particle:
			#	c_cell_coefficient =  snapped((basis_coefficient / pow(grid_domain,2.0)),.01)
			#	c_inverse_I_coefficient = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Inverse_Matrix(other_particle.I)
			#	c_cell_inverse_coefficient = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Matrix_by_Scalar(c_inverse_I_coefficient,c_cell_coefficient,true)
				
			#	particle.C = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Matrix(c_cell_inverse_coefficient,sum_of_c_weighted_velocity_matrix)
				






'''
if artifact.sigma[0][0] <= 0.0 or artifact.sigma[0][1] <= 0.0 or artifact.sigma[1][0] <= 0.0 or artifact.sigma[1][1] <= 0.0:
			if artifact.sigma[0][0] <= 0.0 and artifact.sigma[0][1] <= 0.0 and artifact.sigma[1][0] <= 0.0 and artifact.sigma[1][1] <= 0.0:
				e_sigma =[[ log(.001), log(.001) ],[log(.001), log(.001) ]]
			
			elif artifact.sigma[0][0] <= 0.0 and artifact.sigma[0][1] <= 0.0 and artifact.sigma[1][0] <= 0.0:
				e_sigma =  [[ log(.001), log(.001) ],[ log(.001), log(artifact.sigma[1][1]) ]]
			elif artifact.sigma[0][0] <= 0.0 and artifact.sigma[0][1] <= 0.0 and artifact.sigma[1][1] <= 0.0:
				e_sigma =  [[ log(.001), log(.001) ],[ log(artifact.sigma[1][0]), log(.001) ]]
			elif artifact.sigma[0][0] <= 0.0 and artifact.sigma[1][0] <= 0.0 and artifact.sigma[1][1] <= 0.0:
				e_sigma =  [[ log(.001), log(artifact.sigma[0][1]) ],[log(.001), log(.001) ]]
			elif artifact.sigma[0][1] <= 0.0 and artifact.sigma[1][0] <= 0.0 and artifact.sigma[1][1] <= 0.0:
				e_sigma =  [[ log(artifact.sigma[0][0]), log(.001) ],[ log(.001), log(.001) ]]
			
			elif artifact.sigma[0][0] <= 0.0 and artifact.sigma[0][1] <= 0.0:
				e_sigma =  [[log(.001), log(.001) ],[ log(artifact.sigma[1][0]),log(artifact.sigma[1][1])]]
			elif artifact.sigma[0][0] <= 0.0 and artifact.sigma[1][0] <= 0.0:
				e_sigma =  [[log(.001), log(artifact.sigma[0][1]) ],[ log(.001),log(artifact.sigma[1][1])]]
			elif artifact.sigma[0][0] <= 0.0 and artifact.sigma[1][1] <= 0.0:
				e_sigma =  [[log(.001),log(artifact.sigma[0][1]) ],[ log(artifact.sigma[1][0]), log(.001)]]
			elif artifact.sigma[0][1] <= 0.0 and artifact.sigma[1][0] <= 0.0:
				e_sigma =  [[log(artifact.sigma[0][0]),log(.001) ],[ log(.001),log(artifact.sigma[1][1])]]
			elif artifact.sigma[0][1] <= 0.0 and artifact.sigma[1][1] <= 0.0:
				e_sigma =  [[log(artifact.sigma[0][0]),log(.001) ],[log(artifact.sigma[1][0]),log(.001)]]
			elif artifact.sigma[1][0] <= 0.0 and artifact.sigma[1][1] <= 0.0:
				e_sigma =  [[log(artifact.sigma[0][0]),log(artifact.sigma[0][1]) ],[log(.001),log(.001)]]
		
		else:
			e_sigma = [[ log(artifact.sigma[0][0]), log(artifact.sigma[0][1]) ],[ log(artifact.sigma[1][0]), log(artifact.sigma[1][1]) ]]
	
'''





















#### could be added...

#func Particles_Merging():
	###
	### Determine if particles will merge together...
	###
#	if particle.merging == false and other_particle.merging == false:
#		percentage_covered = get_tree().get_root().get_node("Test Area/Simulation/Particle Interaction").Determine_Coverage(particle,other_particle,contact_with)
						
#		if percentage_covered >= .25:
			### if a certain amount is covered now check to see if both particles is heading in the same direction...
				
#			particle.magnitude = snapped(the_grid[particle]['velocity'].length(),.01)
#			particle.direction = snapped( rad2deg( snapped( atan2(the_grid[particle]['velocity'].y,the_grid[particle]['velocity'].x),.01)),.01)
					
#			other_particle.magnitude = snapped(the_grid[other_particle]['velocity'].length(),.01)
#			other_particle.direction = snapped( rad2deg( snapped( atan2(the_grid[other_particle]['velocity'].y,the_grid[other_particle]['velocity'].x),.01)),.01)
								
								
			### checks if both particles....
			# has the same speed.
			# has the same magnitude.
			# has the same direction.
#			if the_grid[particle]['velocity'] == the_grid[other_particle]['velocity'] and particle.magnitude == other_particle.magnitude and particle.direction == other_particle.direction:
				# the pair of particles is set to merge...
									
#				particle.merging = true
#				other_particle.merging = true
									
#				particle.merge_with = other_particle
#				other_particle.merge_with = particle
									
#				merge_particles = true
#				cycle_of_mergin_particles.append([particle,other_particle])
#				list_of_particles_merged_particles.append([particle,other_particle])
#	else:
		### the other particle is already set to merge with another...
#		pass


func Determine_Coverage(artifact,other_artifact,distance_between):
	### find out how much one particle covers the other particle...
	if distance_between.x == 0 and distance_between.y == 0:
		# the particle relation to the other is directly on top...
		percentage_covered = 100
	elif distance_between.x < 0.0 and distance_between.y < 0.0:
		# the particle relation to the other is to the top left...
		length_covered_area = abs(snapped((other_artifact.surrounding_area.position.x - artifact.surrounding_area.end.x),1))
					
		width_covered_area = abs(snapped((other_artifact.surrounding_area.position.y - artifact.surrounding_area.end.y),1))
					
		covered_area = length_covered_area * width_covered_area
		percentage_covered = snapped((float(covered_area) / artifact.surrounding_area.get_area()),.01)
						
	elif distance_between.x < 0.0 and distance_between.y >= 0.0:
		# the particle relation to the other is to the bottom left...
		length_covered_area = abs(snapped((other_artifact.surrounding_area.position.x - artifact.surrounding_area.end.x),1))
					
		width_covered_area = abs(snapped((other_artifact.surrounding_area.end.y - artifact.surrounding_area.position.y),1))
					
		covered_area = length_covered_area * width_covered_area
		percentage_covered = snapped((float(covered_area) / artifact.surrounding_area.get_area()),.01)
					
	elif distance_between.x >= 0.0 and distance_between.y < 0.0:
		# the particle relation to the other is to the top right...
		length_covered_area = abs(snapped((other_artifact.surrounding_area.end.x - artifact.surrounding_area.position.x),1))
					
		width_covered_area = abs(snapped((other_artifact.surrounding_area.position.y - artifact.surrounding_area.end.y),1))
					
		covered_area = length_covered_area * width_covered_area
		percentage_covered = snapped((float(covered_area) / artifact.surrounding_area.get_area()),.01)
					
	elif distance_between.x >= 0.0 and distance_between.y >= 0.0 :
		# the particle realtion to the other is to the bottom right...
		length_covered_area = abs(snapped((other_artifact.surrounding_area.end.x - artifact.surrounding_area.position.x),1))
					
		width_covered_area = abs(snapped((other_artifact.surrounding_area.end.y - artifact.surrounding_area.position.y),1))
					
		covered_area = length_covered_area * width_covered_area
		percentage_covered = snapped((float(covered_area) / artifact.surrounding_area.get_area()),.01)
					
	
	return percentage_covered


func Merge_Particles(set_of_combining_particles,grid_field,grid_area):
	### particle are merge...
	
	
	### accummulate of the particles info in 1 meta-particle....
	
	# cycle thru each set of combining particles...
	for components in set_of_combining_particles:
		### reset...
		sum_of_parts_mass = 0.00
		sum_of_parts_kernel = Vector2(0.0,0.0)
		sum_of_parts_velocity = Vector2(0.0,0.0)
		
		
		### setup to create a new-particle....
		empirical_particle = load("res://2D Particle.tscn").instantiate()
		empirical_name = '{a}{b}{c}{d}{1}{2}{3}{4}'.format({
		'a':letters_list[int(randi_range(0,len(letters_list)-1))],
		'b':letters_list[int(randi_range(0,len(letters_list)-1))],
		'c':letters_list[int(randi_range(0,len(letters_list)-1))],
		'd':letters_list[int(randi_range(0,len(letters_list)-1))],
		'1':digits_list[int(randi_range(0,len(digits_list)-1))],
		'2':digits_list[int(randi_range(0,len(digits_list)-1))],
		'3':digits_list[int(randi_range(0,len(digits_list)-1))],
		'4':digits_list[int(randi_range(0,len(digits_list)-1))]
		})
	
		### name of the meta-particle...
		empirical_particle.set_name(empirical_name)
		empirical_particle.adjust_size(1)
		empirical_particle.coefficient_of_restitution = 0.75
		### recognize that the particle is made of other particles...
		empirical_particle.made_of = components.duplicate(true)
		
		### accummulate of the particless info in 1 meta-particle....
		for part in empirical_particle.made_of:
			sum_of_parts_mass = snapped((sum_of_parts_mass + grid_field[part]["mass"]),.001)
			sum_of_parts_kernel = Vector2(snapped((sum_of_parts_kernel.x + part.get_position().x),1),
										snapped((sum_of_parts_kernel.y + part.get_position().y),1))
			sum_of_parts_velocity = Vector2(snapped((sum_of_parts_velocity.x + (grid_field[part]["velocity"].x * grid_field[part]["mass"])),.01),
										snapped((sum_of_parts_velocity.y + (grid_field[part]["velocity"].y * grid_field[part]["mass"])),.01))
			
			###particles is removed from the simulation...
			get_tree().get_root().get_node("Test Area/Simulation/Substance").remove_child(part)
			#$"Program".grid_nodes.erase(part)
			get_tree().get_root().get_node("Test Area/Simulation/Program").grid_nodes.erase(part)
		
		set_of_combining_particles.erase(components)
		### setting the parameters of the merged particles...
		empirical_particle.mass = sum_of_parts_mass
		empirical_particle.position = Vector2(snapped((sum_of_parts_kernel.x / len(components) ),1),
											snapped((sum_of_parts_kernel.y / len(components) ),1))
		empirical_particle.velocity = Vector2(snapped((sum_of_parts_velocity.x / empirical_particle.mass),.01),
											snapped((sum_of_parts_velocity.y / empirical_particle.mass) ,.01))
		empirical_particle.stress = [[1.0,1.0],[1.0,1.0]]
		empirical_particle.within_range.append(empirical_particle)
		
		area_multiplier = 1.0
		empirical_particle.surrounding_area = Rect2(Vector2(empirical_particle.position.x - ((grid_area/2.0)*area_multiplier),empirical_particle.position.y - ((grid_area/2.0)*area_multiplier)),Vector2(grid_area*area_multiplier,grid_area*area_multiplier))
		
		
		### the merged particle is added to the simulation...
		get_tree().get_root().get_node("Test Area/Simulation/Substance").add_child(empirical_particle)
		get_tree().get_root().get_node("Test Area/Simulation/Program").grid_nodes[empirical_particle] = {'mass':0.0,'velocity':Vector2(0,0),'momentum':Vector2(0.0,0.0)}
	
	### reset....
	for artifacts in get_tree().get_root().get_node("Test Area/Simulation/Substance").get_children():
		artifacts.relation_to_domain = {}
		artifacts.domain_relation_to_particle = {}
			
	
	for artifacts in get_tree().get_root().get_node("Test Area/Simulation/Substance").get_children():
		for other_artifacts in get_tree().get_root().get_node("Test Area/Simulation/Substance").get_children():
			
			kernel_x = snapped(artifacts.get_position().x - other_artifacts.surrounding_area.get_center().x,.01)
			kernel_y = snapped(artifacts.get_position().y - other_artifacts.surrounding_area.get_center().y,.01)
			kernel_distance = Vector2(kernel_x,kernel_y)
			#flipped_kernel_distance = (other_particle.surrounding_area.get_center() - particle.get_position())# / cell_size
			flipped_kernel_x = snapped(other_artifacts.surrounding_area.get_center().x - artifacts.get_position().x,.01)
			flipped_kernel_y = snapped(other_artifacts.surrounding_area.get_center().y - artifacts.get_position().y,.01)
			flipped_kernel_distance = Vector2(flipped_kernel_x,flipped_kernel_y)
					
			artifacts.relation_to_domain[other_artifacts] = kernel_distance
			artifacts.domain_relation_to_particle[other_artifacts] = flipped_kernel_distance
			




















































































































###### old code from old verison....





func Merge_Particles(combining_particles):
	### Merging Particles together...
	var sum_of_parts_mass
	var sum_of_parts_kernel
	var sum_of_parts_velocity
	
	for components in combining_particles:
		### each grouping of merging particles...
		
		### reset...
		sum_of_parts_mass = 0.00
		sum_of_parts_kernel = Vector2(0.0,0.0)
		sum_of_parts_velocity = Vector2(0.0,0.0)
		
		
		### Creating a new particle...
		randomize()
		#used for name of particle...
		var lock_1 : int = int(rand_range(0,len(letters_list)-1))
		var lock_2 : int = int(rand_range(0,len(letters_list)-1))
		var lock_3 : int = int(rand_range(0,len(letters_list)-1))
		var lock_4 : int = int(rand_range(0,len(letters_list)-1))
		var lock_5 : int = int(rand_range(0,len(digits_list)-1))
		var lock_6 : int = int(rand_range(0,len(digits_list)-1))
		var lock_7 : int = int(rand_range(0,len(digits_list)-1))
		var lock_8 : int = int(rand_range(0,len(digits_list)-1))
		var merged_particle = load("res://2d_(particle).tscn").instance()
		particle_name = '{a}{b}{c}{d}{1}{2}{3}{4}'.format({
			'a':letters_list[lock_1],
			'b':letters_list[lock_2],
			'c':letters_list[lock_3],
			'd':letters_list[lock_4],
			'1':digits_list[lock_5],
			'2':digits_list[lock_6],
			'3':digits_list[lock_7],
			'4':digits_list[lock_8]
			})
			
		### setting other parameters of the particle...
		merged_particle.set_name(particle_name)
		##print(merged_particle, ' Merged Particle')
		merged_particle.made_of = components.duplicate(true)
		##print(components,' should be list of merging particles...')
		### cycle thru each particle that will combine...
		for parts in merged_particle.made_of:
			### parameters of the merging particles...
			sum_of_parts_mass = stepify((sum_of_parts_mass + grid_nodes[parts.get_name()]["mass"]),.001)
			sum_of_parts_kernel = Vector2(stepify((sum_of_parts_kernel.x + parts.get_position().x),1),
										stepify((sum_of_parts_kernel.y + parts.get_position().y),1))
			sum_of_parts_velocity = Vector2(stepify((sum_of_parts_velocity.x + (grid_nodes[parts.get_name()]["velocity"].x * grid_nodes[parts.get_name()]["mass"])),.01),
										stepify((sum_of_parts_velocity.y + (grid_nodes[parts.get_name()]["velocity"].y * grid_nodes[parts.get_name()]["mass"])),.01))
			##print(parts,' particle info grabbed')
			##print(sum_of_parts_velocity,' ')
			### adjusting the substance...
			get_tree().get_root().get_node("Test Area/Simulation/program").substance.erase(parts)
			get_tree().get_root().get_node("Test Area/Simulation/program").grid_nodes.erase(parts)
			
			get_tree().get_root().get_node("Test Area/Simulation/substance").remove_child(parts)
		
		merged_particle.mass = sum_of_parts_mass
		merged_particle.position = Vector2(stepify((sum_of_parts_kernel.x / len(components) ),1),
											stepify((sum_of_parts_kernel.y / len(components) ),1))
		merged_particle.velocity = Vector2(stepify((sum_of_parts_velocity.x / merged_particle.mass),.01),
											stepify((sum_of_parts_velocity.y / merged_particle.mass) ,.01))
		##print(merged_particle.mass,' ')
		##print(merged_particle.position, ' ')
		##print(merged_particle.velocity,' ')
		### particles are merged...
		pair_of_merging_particles.erase(components)
		### Merged Particle Added...
		get_tree().get_root().get_node("Test Area/Simulation/program").substance.append(merged_particle)
		get_tree().get_root().get_node("Test Area/Simulation/program").grid_nodes[merged_particle.get_name()] = {'mass':0.0,'velocity':Vector2(0,0),'momentum':Vector2(0.0,0.0)}
		get_tree().get_root().get_node("Test Area/Simulation/substance").add_child(merged_particle)

### also ...

if particle.merging == false:
			for other_particle in material:
				# ...
				if particle != other_particle and other_particle.merging == false and particle.merging == false:
					### if the particles collided with each other...
					distance_between_x = stepify((other_particle.position.x - particle.get_position().x),1)
					distance_between_y = stepify((other_particle.position.y - particle.get_position().y),1)
					#var distance_between_normal = Vector2(distance_between_x,distance_between_y).normalized()
					distance_between = Vector2(distance_between_x,distance_between_y)
					
					#distance_between = Vector2(stepify(distance_between_normal.x,.1),stepify(distance_between_normal.y,.1))
					#distance_between = Vector2(distance_between_x,distance_between_y).normalized()
					##print(distance_between , ' distance between particles')
					#other_particle.surrounding_area = Rect2(Vector2(stepify(other_particle.position.x - (cell_size/2)*area_multiplier,.01),stepify(other_particle.position.y - (cell_size/2)*area_multiplier,.01)),Vector2(cell_size*area_multiplier,cell_size*area_multiplier))
					#if distance_between.x > (cell_size/2) and distance_between.y > (cell_size/2) and 
					#distance_between.x > (cell_size/2) and distance_between_y < (cell_size/2) and 
					#
					### maintenance to determine how to handle math...
					#if distance_between.x == 0.0 and distance_between.y == 0.0:
						### Completely Covered
					#	percentage_covered = 100
					#	particle.within_range.append(other_particle)
					#for y in range(-(cell_size),(cell_size)+1):
					#	for x in range(-(cell_size),(cell_size)+1):
					#		PoolVector2Array.append(Vector2(x,y))
					#	particle.covered_by_at[other_particle] = percentage_covered
						##print(percentage_covered,' is covered by the 2nd particle ')
					if distance_between.x in PoolRealArray(range(-(cell_size),(cell_size)+1)) and distance_between.y in PoolRealArray(range(-(cell_size),(cell_size)+1)):
						### particles is in any contact with another particle ...
						##print('touching')
						
						if particle.within_range.has(other_particle):
							### already contain/recognized
							pass
						else:
							particle.within_range.append(other_particle)
						
						
						### if a particle covers a certain amount of another particle...
						if distance_between.x in PoolRealArray(range(0.0,(cell_size)+1)) and distance_between.y in PoolRealArray(range(0.0,-(cell_size)+1)):
							# bottomleft starting point...
							# warning-ignore:narrowing_conversion
							length_covered_area = abs(stepify((other_particle.surrounding_area.position.x - particle.surrounding_area.end.x),1))
							# warning-ignore:narrowing_conversion
							width_covered_area = abs(stepify((other_particle.surrounding_area.end.y - particle.surrounding_area.position.y),1))
							##print(length_covered_area, ' should be the length')
							##print(width_covered_area, ' should be the width')
							covered_area = length_covered_area * width_covered_area
							percentage_covered = stepify((float(covered_area) / particle.surrounding_area.get_area()),.01)
							
							# if the particle is covered by another by 50 %...
							#if percentage_covered >= .50:
							#	particle.covered_by_at[other_particle.get_name()]  = percentage_covered
							#	identify_possible_merging.append(particle)
							#else:
								### Not schedule to merge...
							#	pass
							##print(percentage_covered,' is covered by the 2nd particle ')
							
						elif distance_between.x in PoolRealArray(range(0.0,(cell_size)+1)) and distance_between.y in PoolRealArray(range(0.0,(cell_size)+1)):
							# topleft
							# warning-ignore:narrowing_conversion
							length_covered_area = abs(stepify((other_particle.surrounding_area.position.x - particle.surrounding_area.end.x),1))
							# warning-ignore:narrowing_conversion
							width_covered_area = abs(stepify((other_particle.surrounding_area.position.y - particle.surrounding_area.end.y),1))
							##print(length_covered_area, ' should be the length')
							##print(width_covered_area, ' should be the width')
							covered_area = length_covered_area * width_covered_area
							percentage_covered = stepify((float(covered_area) / particle.surrounding_area.get_area()),.01)
							
							
							#particle.covered_by_at[other_particle.get_name()]  = percentage_covered
							##print(percentage_covered,' is covered by the 2nd particle ')
							
						elif distance_between.x in PoolRealArray(range(0.0,-(cell_size)+1)) and distance_between.y in PoolRealArray(range(0.0,(cell_size)+1)):
							# topright
							# warning-ignore:narrowing_conversion
							length_covered_area = abs(stepify((other_particle.surrounding_area.end.x - particle.surrounding_area.position.x),1))
							# warning-ignore:narrowing_conversion
							width_covered_area = abs(stepify((other_particle.surrounding_area.position.y - particle.surrounding_area.end.y),1))
							##print(length_covered_area, ' should be the length')
							##print(width_covered_area, ' should be the width')
							covered_area = length_covered_area * width_covered_area
							percentage_covered = stepify((float(covered_area) / particle.surrounding_area.get_area()),.01)
							
							
							#particle.covered_by_at[other_particle.get_name()]  = percentage_covered
							##print(percentage_covered,' is covered by the 2nd particle ')
							
						elif distance_between.x in PoolRealArray(range(0.0,-(cell_size)+1)) and distance_between.y in PoolRealArray(range(0.0,-(cell_size)+1)):
							# bottom right
							# warning-ignore:narrowing_conversion
							length_covered_area = abs(stepify((other_particle.surrounding_area.end.x - particle.surrounding_area.position.x),1))
							# warning-ignore:narrowing_conversion
							width_covered_area = abs(stepify((other_particle.surrounding_area.end.y - particle.surrounding_area.position.y),1))
							##print(length_covered_area, ' should be the length')
							##print(width_covered_area, ' should be the width')
							covered_area = length_covered_area * width_covered_area
							percentage_covered = stepify((float(covered_area) / particle.surrounding_area.get_area()),.01)
							
							
							#particle.covered_by_at[other_particle.get_name()]  = percentage_covered
							##print(percentage_covered,' is covered by the 2nd particle ')
						
						# if the particle is covered by another by 50 %...
						if percentage_covered >= .50:
							#particle.covered_by_at[other_particle]  = percentage_covered
							
							#if identify_possible_merging.has(particle) == false:
							#	identify_possible_merging.append(particle)
							#covering.append(other_particle)
							#determine if particles in heading in the same direction...
							particle.magnitude = stepify(grid_data[particle.get_name()]['velocity'].length(),.01)
							particle.direction = stepify( rad2deg( stepify( atan2(grid_data[particle.get_name()]['velocity'].y,grid_data[particle.get_name()]['velocity'].x),.01)),.01)
							
							other_particle.magnitude = stepify(grid_data[other_particle.get_name()]['velocity'].length(),.01)
							other_particle.direction = stepify( rad2deg( stepify( atan2(grid_data[other_particle.get_name()]['velocity'].y,grid_data[other_particle.get_name()]['velocity'].x),.01)),.01)
							
							if grid_data[particle.get_name()]['velocity'] == grid_data[other_particle.get_name()]['velocity'] and particle.magnitude == other_particle.magnitude and particle.direction == other_particle.direction:
								##print('Merging')
								##print(particle.get_name(), ' with ',other_particle.get_name())
								### mark each particle so dont merge again...
								particle.merging = true
								other_particle.merging = true
								
								particle.within_range = []
								other_particle.within_range = []
								
								pair_of_merging_particles.append([particle,other_particle])
								#merging_particles = true 
						#else:
							### the particles collide by
							### Not covered by another particle by certain percentage...
						#	pass
						
					else:
						# then the particles do not touch...
						##print('not touching')
						if particle.within_range.has(other_particle):
							particle.within_range.erase(other_particle)
						pass
				else:
					###  itself and or already set to merge with another...
					
					
					pass
		else:
			### particle is already merging with another...
			pass




########################################
########################################
########################################




#### misfit cop::: Keep this...

### final_velocity = (Mass_1_i * Velocity_1_i) + (Mass_2_i * Velocity_2_i)) / (Mass_1 + Mass_2)
						


### the collision is perfect elastic...
						
						## angle-free :
						# velocity1final = velocity1initial - ( ((2*mass2) * dotproduct(velocitykernel,positionkernel) * positionkernel ) / ( (mass1+mass2) * (|| posiitonkernel ||^2) )
						# for velocity1final notes :: velocitykernel = v1 - v2,positionkernel = x1-x2,|| || = magnitude = squareroot(x^2+y^2)
						
						# velocity2final = velocity2initial - ( ((2*mass1) * dotproduct(velocitykernel,positionkernel) * positionkernel ) / ( (mass1+mass2) * (|| positionkernel ||^2) ) 
						#for velocity2final notes :: velocitykernel = v2 - v1,positionkernel = x2 - x1,|| = magnitude = squareroot(x^2+y^2)
						
						### final velcotiy 1 = particle
						### final velcoity 2 = wall
						elastic_term_1_of_vector_1 = 2.0 * baluster['window outline'][wall]['mass']
						elastic_term_1_of_vector_2 = 2.0 * structure[particle]['mass']
						
						elastic_mass_addition = structure[particle]['mass'] + baluster['window outline'][wall]['mass']
						
						### wall center is along "y" to the object 
						#var wall_center = Vector2(center_of_window.x,ProjectSettings.get_setting('display/window/size/height'))
						#var wall_center = Vector2(particle.position.x,particle.position.y+cell_size)
						wall_center = Vector2(particle.position.x,baluster['window outline'][wall]['outline'])
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
						




















### requires 4 seperate equations...
						### line of impact = line between object's masses
						### plane contact = perpendicular to the line of impact
						### i = initial ; f = final ; line of impact = x ; plane contact = y
						### e = coefficient of reinstition
						#
						#1.(mass_1_i * velocity_1_x_i) cos(angle_in) + (mass_2_i * velocity_2_x_i) cos(angle_in) =  (mass_1_f * velocity_1_x_f) cos(angle_out) + (mass_2_f * velocity_2_x_f) cos(angle_out)
						#
						#2. velocity_2_f - velocity_1_f = e * ( velocity_2_i - velocity_1_i )
						#  .... velocity_2_f_x,velocity_2_f_y - velocity_1_f_x,velocity_1_f_y = e * (velocity_1_i_x,velocity_1_i_y - velocity_1_i_x,velocity_1_i_y)
						#
						# 3. velocity_y_2_f = velocity_y_2_i
						# 4. velocity_y_1_f = velocity_y_1_i
						#
						######
						#####
						####
						# 1. find angle_i between ( velocity_1_initial : line of impact )
						# .... cos(angle_i) = vector_a.dot_product(vector_b) / |vector_a|*|vector_b|
						#         ... 1. dot_product between vectors
						#         ... 2. calculate the length of the vectors then multiply them together.
						#         ... 3. dot_product results divided by multiplied lengths results.
						#         ... 4. angle_i = arc_cos( results from dot_product / Multipled Lengths)
						# 
						# 2. tan(angle_f) = e * tan(angle_i)
						#    .... angle_f = arc_tan(e * tan(angle_i))
						#
						#3. other leg : velocity_y_2_f * tan(angle_f)
						#
						#4. pythagorean theorem : a^2 + b^2 = c^2
						#
						# find angle....
						var wall_center = Vector2(particle.position.x,(baluster['window outline'][wall]['outline']))
						var line_of_impact = wall_center
						print(" ")
						#var angle = rad2deg(structure[particle]['velocity'].angle_to(line_of_impact))
						var incoming_angle = rad2deg(particle.get_position().angle_to(line_of_impact))
						print(incoming_angle,' in coming angle')
						var outgoing_angle = rad2deg(atan(baluster['window outline'][wall]['coefficient of restitution'] * incoming_angle ))
						print(outgoing_angle,' out going angle')
						var other_leg = structure[particle]['velocity'].y * tan(outgoing_angle)
						
						print(other_leg,' x leg')
						structure[particle]['velocity'] = structure[particle]['velocity'] - Vector2(other_leg,structure[particle]['velocity'].y)
						
####
####











#var manual_dot_vector = (elastic_velocity_kernel_of_vector_1.x * elastic_position_kernel_of_vector_1.x) + (elastic_velocity_kernel_of_vector_1.y * elastic_position_kernel_of_vector_1.y)
						














##################
################## Down below does not matter .....(maybe...)








#####
##### from mls-mpm

#the_grid[particle]['mass'] = the_grid[particle]['mass'] + (the_grid[other_particle]['mass'] * weight_interpolation)
			
			#print(the_grid[particle]['mass'],' the lump sum')
			#building_momentum = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Matrix_by_Vector2_to_Vector2(other_particle.C,particle.domain_relation_to_particle[other_particle])
			#print(building_momentum,' building momentum.')
			#print(other_particle.velocity,' checking velocity')
			#the_grid[particle]['momentum'] = the_grid[particle]['momentum'] + ((the_grid[other_particle]['mass'] * weight_interpolation) * (other_particle.velocity + building_momentum))
			#standby_momentum_x = snapped((the_grid[particle]['momentum'].x + ( snapped((the_grid[other_particle]['mass'] * weight_interpolation),.01) * (snapped((other_particle.velocity.x + building_momentum.x),.01)))),.01)
			#standby_momentum_y = snapped((the_grid[particle]['momentum'].y + ( snapped((the_grid[other_particle]['mass'] * weight_interpolation),.01) * (snapped((other_particle.velocity.y + building_momentum.y),.01)))),.01)
			#the_grid[particle]['momentum'] = Vector2(standby_momentum_x,standby_momentum_y)
			
			#print(the_grid[particle]['momentum'],' momentum. during')
			


			#print(the_grid[particle]['mass'],' the lump sum')
			#building_momentum = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Matrix_by_Vector2_to_Vector2(other_particle.C,particle.domain_relation_to_particle[other_particle])
			#print(building_momentum,' building momentum.')
			#print(other_particle.velocity,' checking velocity')
			#the_grid[particle]['momentum'] = the_grid[particle]['momentum'] + ((the_grid[other_particle]['mass'] * weight_interpolation) * (other_particle.velocity + building_momentum))
			#standby_momentum_x = snapped((the_grid[particle]['momentum'].x + ( snapped((the_grid[other_particle]['mass'] * weight_interpolation),.01) * (snapped((other_particle.velocity.x + building_momentum.x),.01)))),.01)
			#standby_momentum_y = snapped((the_grid[particle]['momentum'].y + ( snapped((the_grid[other_particle]['mass'] * weight_interpolation),.01) * (snapped((other_particle.velocity.y + building_momentum.y),.01)))),.01)
			#the_grid[particle]['momentum'] = Vector2(standby_momentum_x,standby_momentum_y)
			
			#print(the_grid[particle]['momentum'],' momentum. during')
			
			# Forces : f = - sum of: weight_interpolation * particle.volume * D^1 * (deviation_energy_density / deviation_particle.F ) * particle.F * particle.F^T * flipped_kernel_distance
			# or if energy density is unknown
			# Forces : f =  - sum of: particle.volume * particle.stress * gradient_flipped_kernel + gravity
			# additional : f = weight_interpolation * Q * gradient_weight
			# Q = time_passed * particle.volume * D^1 * particle.volume * D^1 * (deviation_energy_density / deviation_particle.F ) * particle.F * particle.F^T + particle.mass * particle.C
			# or
			### Q = time_passed * particle.volume * particle.stress * gradient_weight + gravity + particle.mass * particle.C
			
			#Q_term_1 = weight_interpolation * time_passed * other_particle.volume
			#var Q_term_2 = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Matrix_by_Scalar(other_particle.stress,Q_term_1,true)
			#var Q_term_3 = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Matrix_by_Vector2_to_Vector2(Q_term_2,gradient_weight_interpolation)
			#var Q_term_4 = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Matrix_by_Scalar(other_particle.C,other_particle.mass,true)
			#var Q_term_5 =  get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Matrix_by_Vector2_to_Vector2(Q_term_4,particle.domain_relation_to_particle[other_particle])
			#the_grid[particle]['Q'] = the_grid[particle]['Q'] + -(Q_term_3 + Q_term_5)
			#Q_x = snapped((the_grid[particle]['force'].x + snapped(-(Q_term_3.x + Q_term_5.x),.01)),.01)
			#Q_y = snapped((the_grid[particle]['force'].y + snapped(-(Q_term_3.y + Q_term_5.y),.01)),.01)
			#the_grid[particle]['force'] = Vector2(Q_x,Q_y)























# establish the particle grid domain...
			#area_multiplier = 1
			#particle.surrounding_area = Rect2i(Vector2(particle.position.x - ((grid_domain/2)*area_multiplier),particle.position.y - ((grid_domain/2)*area_multiplier)),Vector2(grid_domain*area_multiplier,grid_domain*area_multiplier))
			#other_particle.surrounding_area = Rect2i(Vector2(snapped(other_particle.position.x - (grid_domain/2)*area_multiplier,1),snapped(other_particle.position.y - (grid_domain/2)*area_multiplier,1)),Vector2(grid_domain*area_multiplier,grid_domain*area_multiplier))
		
#func B_formation(component,weight,other_particle_grid_velocity):
	### B is the sum of particle.velocity * weight_interpolation.

#	for other_particle in component.within_range:
			
#		kernel_distance = component.get_position() - other_particle.surrounding_area.position.get_center()
#		flipped_kernel_distance = other_particle.surrounding_area.position.get_center() - component.get_position()
			
#		component.velocity += (weight * other_particle_grid_velocity)
	
	#return B_source.velocity


						## angle-free :
						
						# velocity1final = velocity1initial - (2*mass2) * dotproduct(velocitykernel,positionkernel) / ( (mass1+mass2) * (|| posiitonkernel ||^2) ) * positionkernel 
						# for velocity1final notes :: velocitykernel = v1 - v2,positionkernel = x1-x2,|| || = magnitude = squareroot(x^2+y^2)
						
						# velocity2final = velocity2initial - (2*mass2) * dotproduct(velocitykernel,positionkernel) / ( (mass1+mass2) * (|| posiitonkernel ||^2) ) * positionkernel 
						#for velocity2final notes :: velocitykernel = v2 - v1,positionkernel = x2 - x1,|| = magnitude = squareroot(x^2+y^2)
						
						#particle = velocity 1
						#wall = velocity 2
						var wall_center = Vector2(window_center.x,0.0)
						
						#var wall_center = Vector2(particle.get_position().x,(0.0-(grid_domain/2))) # causes the particle to bounce in place.
						#var wall_center = Vector2(window_center.x,-(window_center.y/2.0))
						print(wall_center,' wall center')
						var elastic_term_1_of_vector_1 = 2 * barriers['window outline'][wall]['mass']
						print(elastic_term_1_of_vector_1, ' 2 * wall_mass')
						var elastic_mass_addition = the_grid[particle]['mass'] + barriers['window outline'][wall]['mass']
						print(elastic_mass_addition,' particle mass + wall_mass')
						var elastic_position_kernel_of_vector_1 = (particle.get_position() - wall_center)# / grid_domain
						print(elastic_position_kernel_of_vector_1,' ',particle.get_position(),' particle position ','- ',wall_center,' wall-position')
						var elastic_velocity_kernel_of_vector_1 = the_grid[particle]['velocity'] - barriers['window outline'][wall]['velocity']
						print(elastic_velocity_kernel_of_vector_1,' particle velocity - wall velocity')
						var elastic_dot_vector_1 = elastic_velocity_kernel_of_vector_1.dot(elastic_position_kernel_of_vector_1)
						var test_switch = elastic_position_kernel_of_vector_1.dot(elastic_velocity_kernel_of_vector_1)
						print(elastic_dot_vector_1,' velocity dot position')
						print(test_switch,' position dot velocity')
						#var elastic_magnitude_kernel_of_vector_1 = pow(sqrt( (pow(elastic_position_kernel_of_vector_1.x,2.0)+pow(elastic_position_kernel_of_vector_1.y,2.0))),2.0 )
						var elastic_magnitude_kernel_of_vector_1 = pow(elastic_position_kernel_of_vector_1.length(),2.0)
						print(elastic_magnitude_kernel_of_vector_1,' magnitude')
						print(the_grid[particle]['velocity'],' before collision')
						#the_grid[particle]['velocity'] = the_grid[particle]['velocity']  - ( elastic_term_1_of_vector_1 / elastic_mass_addition ) * ( elastic_dot_vector_1 / elastic_magnitude_kernel_of_vector_1 ) * ( elastic_position_kernel_of_vector_1 )
						the_grid[particle]['velocity'] = the_grid[particle]['velocity']  - ((elastic_term_1_of_vector_1 * elastic_dot_vector_1 * elastic_position_kernel_of_vector_1) /( elastic_mass_addition * elastic_magnitude_kernel_of_vector_1) )
						print(the_grid[particle]['velocity'],' after collision adjusted')
						



















				if barriers['window outline'][wall]['coefficient of reinstition'] == 2.0:
					### the collision is perfect elastic +1... (explosion)
					pass
				elif barriers['window outline'][wall]['coefficient of reinstition'] == 1.0:
					### the collision is perfect elastic...
					## angle-free :
					
					# velocity1final = velocity1initial - ( ((2*mass2) * dotproduct(velocitykernel,positionkernel) * positionkernel ) / ( (mass1+mass2) * (|| posiitonkernel ||^2) )
					# for velocity1final notes :: velocitykernel = v1 - v2,positionkernel = x1-x2,|| || = magnitude = squareroot(x^2+y^2)
					
					# velocity2final = velocity2initial - ( ((2*mass1) * dotproduct(velocitykernel,positionkernel) * positionkernel ) / ( (mass1+mass2) * (|| positionkernel ||^2) ) 
					#for velocity2final notes :: velocitykernel = v2 - v1,positionkernel = x2 - x1,|| = magnitude = squareroot(x^2+y^2)
					
					### final velcotiy 1 = particle
					### final velcoity 2 = wall
					var elastic_term_1_of_vector_1 = 2 * barriers['window outline'][wall]['mass']
					var elastic_term_1_of_vector_2 = 2 * the_grid[particle]['mass']
					
					var elastic_mass_addition = the_grid[particle]['mass'] + barriers['window outline'][wall]['mass']
					var elastic_dot_vector_1
					var elastic_dot_vector_2
					var elastic_position_kernel_of_vector_1
					var elastic_position_kernel_of_vector_2
					var elastic_velocity_kernel_of_vector_1
					var elastic_velocity_kernel_of_vector_2
					var elastic_magnitude_kernel_of_vector_1
					var elastic_magnitude_kernel_of_vector_2
					
					var touch_wall : bool = false
					
					
					print(the_grid[particle]['velocity'],' before collision')
					
					if wall == 'top' and particle.surrounding_area.position.y <= barriers['window outline'][wall]['outline']:
						var wall_center = Vector2(window_center.x,0)
						elastic_position_kernel_of_vector_1 = particle.get_position() - wall_center
						elastic_position_kernel_of_vector_2 = wall_center - particle.get_position()
						### find final-velocity 1::
						elastic_velocity_kernel_of_vector_1 = the_grid[particle]['velocity'] - barriers['window outline'][wall]['velocity']
						elastic_dot_vector_1 = elastic_velocity_kernel_of_vector_1.dot(elastic_position_kernel_of_vector_1)
						### find final-velocit 2::
						elastic_velocity_kernel_of_vector_2  = barriers['window outline'][wall]['velocity'] - the_grid[particle]['velocity']
						elastic_dot_vector_2 = elastic_velocity_kernel_of_vector_2.dot(elastic_position_kernel_of_vector_2)
					
						elastic_magnitude_kernel_of_vector_1 = sqrt( (pow(elastic_position_kernel_of_vector_1.x,2.0)+pow(elastic_position_kernel_of_vector_1.y,2.0)) )
						elastic_magnitude_kernel_of_vector_2 = sqrt( (pow(elastic_position_kernel_of_vector_2.x,2.0)+pow(elastic_position_kernel_of_vector_2.y,2.0)) )
						
						touch_wall = true
						
					elif wall == 'right' and particle.surrounding_area.end.x >= barriers['window outline'][wall]['outline']:
						var wall_center = Vector2(ProjectSettings.get_setting('display/window/size/width'),window_center.y)
						elastic_position_kernel_of_vector_1 = particle.get_position() - wall_center
						elastic_position_kernel_of_vector_2 = wall_center - particle.get_position()
						### find final-velocity 1::
						elastic_velocity_kernel_of_vector_1 = the_grid[particle]['velocity'] - barriers['window outline'][wall]['velocity']
						elastic_dot_vector_1 = elastic_velocity_kernel_of_vector_1.dot(elastic_position_kernel_of_vector_1)
						### find final-velocit 2::
						elastic_velocity_kernel_of_vector_2  = barriers['window outline'][wall]['velocity'] - the_grid[particle]['velocity']
						elastic_dot_vector_2 = elastic_velocity_kernel_of_vector_2.dot(elastic_position_kernel_of_vector_2)
					
						elastic_magnitude_kernel_of_vector_1 = sqrt( (pow(elastic_position_kernel_of_vector_1.x,2.0)+pow(elastic_position_kernel_of_vector_1.y,2.0)) )
						elastic_magnitude_kernel_of_vector_2 = sqrt( (pow(elastic_position_kernel_of_vector_2.x,2.0)+pow(elastic_position_kernel_of_vector_2.y,2.0)) )
					
						touch_wall = true
					
					elif wall == 'bottom' and particle.surrounding_area.end.y >= barriers['window outline'][wall]['outline']:
						var wall_center = Vector2(window_center.x,ProjectSettings.get_setting('display/window/size/height'))
						elastic_position_kernel_of_vector_1 = particle.get_position() - wall_center
						elastic_position_kernel_of_vector_2 = wall_center - particle.get_position()
						### find final-velocity 1::
						elastic_velocity_kernel_of_vector_1 = the_grid[particle]['velocity'] - barriers['window outline'][wall]['velocity']
						elastic_dot_vector_1 = elastic_velocity_kernel_of_vector_1.dot(elastic_position_kernel_of_vector_1)
						### find final-velocit 2::
						elastic_velocity_kernel_of_vector_2  = barriers['window outline'][wall]['velocity'] - the_grid[particle]['velocity']
						elastic_dot_vector_2 = elastic_velocity_kernel_of_vector_2.dot(elastic_position_kernel_of_vector_2)
					
						elastic_magnitude_kernel_of_vector_1 = sqrt( (pow(elastic_position_kernel_of_vector_1.x,2.0)+pow(elastic_position_kernel_of_vector_1.y,2.0)) )
						elastic_magnitude_kernel_of_vector_2 = sqrt( (pow(elastic_position_kernel_of_vector_2.x,2.0)+pow(elastic_position_kernel_of_vector_2.y,2.0)) )
					
						touch_wall = true
					
					elif wall == 'left' and particle.surrounding_area.position.x <= barriers['window outline'][wall]['outline']:
						var wall_center = Vector2(0,window_center.y)
						elastic_position_kernel_of_vector_1 = particle.get_position() - wall_center
						elastic_position_kernel_of_vector_2 = wall_center - particle.get_position()
						### find final-velocity 1::
						elastic_velocity_kernel_of_vector_1 = the_grid[particle]['velocity'] - barriers['window outline'][wall]['velocity']
						elastic_dot_vector_1 = elastic_velocity_kernel_of_vector_1.dot(elastic_position_kernel_of_vector_1)
						### find final-velocit 2::
						elastic_velocity_kernel_of_vector_2  = barriers['window outline'][wall]['velocity'] - the_grid[particle]['velocity']
						elastic_dot_vector_2 = elastic_velocity_kernel_of_vector_2.dot(elastic_position_kernel_of_vector_2)
					
						elastic_magnitude_kernel_of_vector_1 = sqrt( (pow(elastic_position_kernel_of_vector_1.x,2.0)+pow(elastic_position_kernel_of_vector_1.y,2.0)) )
						elastic_magnitude_kernel_of_vector_2 = sqrt( (pow(elastic_position_kernel_of_vector_2.x,2.0)+pow(elastic_position_kernel_of_vector_2.y,2.0)) )
					
						touch_wall = true
					
					### 
					if touch_wall == true:
						
						print(elastic_term_1_of_vector_1,' 2 * wall mass')
						print(elastic_dot_vector_1,' dot product')
						print(elastic_position_kernel_of_vector_1,' position kernel')
						print(elastic_mass_addition,' particle mass + wall mass')
						print(elastic_magnitude_kernel_of_vector_1,' magnitude of kernel')
						the_grid[particle]['velocity'] = the_grid[particle]['velocity']  - ( elastic_term_1_of_vector_1 / elastic_mass_addition ) * ( elastic_dot_vector_1 / elastic_magnitude_kernel_of_vector_1 ) * ( elastic_position_kernel_of_vector_1 )
						print(the_grid[particle]['velocity'],' after collision')
						touch_wall = false
						break
					
				elif barriers['window outline'][wall]['coefficient of reinstition'] == 0.5:
					### the collision is partial inelastic...
					pass
				elif barriers['window outline'][wall]['coefficient of reinstition'] == 0.0:
					### the collision is complete inelastic...
					#:: mass_1*velocity_1+mass2*velocity2 = (mass1+mass2)*velocity
				# velocity = mass1*velocity1+mass2*velocity2 /  (mass1+mass2)
				#	print(the_grid[particle]['velocity'],' velocity at the end before wall')
				#for wall in barriers['window outline'].keys():
					if wall == 'top' and particle.surrounding_area.position.y <= barriers['window outline'][wall]['outline']:
						velocity_from_wall_bounce = (the_grid[particle]['momentum'] + barriers['window outline'][wall]['momentum']) /  (the_grid[particle]['mass'] / barriers['window outline'][wall]['mass'])
				#		print(velocity_from_wall_bounce,' checks')
						the_grid[particle]['velocity'] = the_grid[particle]['velocity' ] -(velocity_from_wall_bounce / time_passed)
						break
					if wall == 'right' and particle.surrounding_area.end.x >= barriers['window outline'][wall]['outline'] :
						velocity_from_wall_bounce = (the_grid[particle]['momentum'] + barriers['window outline'][wall]['momentum']) /  (the_grid[particle]['mass'] / barriers['window outline'][wall]['mass'])
				#		print(velocity_from_wall_bounce,' checks')
						the_grid[particle]['velocity'] = the_grid[particle]['velocity'] - (velocity_from_wall_bounce / time_passed)
						break
					if wall == 'bottom' and particle.surrounding_area.end.y >= barriers['window outline'][wall]['outline'] :
						velocity_from_wall_bounce = (the_grid[particle]['momentum'] + barriers['window outline'][wall]['momentum']) /  (the_grid[particle]['mass'] / barriers['window outline'][wall]['mass'])
				#		print(velocity_from_wall_bounce,' checks')
						the_grid[particle]['velocity'] = the_grid[particle]['velocity'] - (velocity_from_wall_bounce / time_passed)
						break
					if wall == 'left' and particle.surrounding_area.position.x <= barriers['window outline'][wall]['outline']:
						velocity_from_wall_bounce = (the_grid[particle]['momentum'] + barriers['window outline'][wall]['momentum']) /  (the_grid[particle]['mass'] / barriers['window outline'][wall]['mass'])
				#		print(velocity_from_wall_bounce,' checks')
						the_grid[particle]['velocity'] = the_grid[particle]['velocity'] - (velocity_from_wall_bounce / time_passed)
						break
				#	print(the_grid[particle]['velocity'],' velocity at the end after wall')
					#velocity_from_wall_bounce = (the_grid[particle]['momentum'] + wall_momentum) /  (the_grid[particle]['mass']+wall_mass)
					#print(velocity_from_wall_bounce,' checks')
					#the_grid[particle]['velocity'] = the_grid[particle]['velocity'] - velocity_from_wall_bounce
			
			#for wall in barriers['window outline'].keys():
			#	if particle.surrounding_area.intersects(barriers['window outline'][wall]):
					###print(' Touching ', wall,' wall')
					### test...
					#particle.velocity = Vector2(0.0,0.0)
			#		the_grid[particle]['velocity'] = Vector2(0.0,0.0)
			#		particle.velocity = Vector2(0.0,0.0)
			#		particle.B = [[0.0,0.0],[0.0,0.0]]
			#		particle.direction = snapped( rad2deg( snapped( atan2(the_grid[particle]['velocity'].y,the_grid[particle]['velocity'].x),.01)),.01)
			else:
				### the particle will travel thru the border...
				# will either continue or disappear from the simulation...
				if barriers['window outline'][wall]['coefficient of reinstition'] == 'thru':
					###print('exist out of sight')
					pass
				elif barriers['window outline'][wall]['coefficient of reinstition'] == 'disappear':
					###print('particles should be removed from the substance...')
					pass
					#grid_nodes.erase(particle.get_name())
						
					#material.erase(particle)
					#get_tree().get_root().get_node("Test Area/Simulation/substance").remove_child(particle)
				

	for particle in material:
		if particle.merging == false:
			for other_particle in particle.within_range:
				if particle != other_particle and other_particle.merging == false:
					
					distance_between = particle.get_position() - other_particle.get_position()
					print(distance_between, ' checking collisions')
					if distance_between.x in range(-(grid_domain),(grid_domain)+1) and distance_between.y in range(-(grid_domain),(grid_domain)+1):
						
						if particle.within_range.has(other_particle):
							pass
						else:
							### recognized to come into contact.
							particle.within_range.append(other_particle)
							# establish the particle grid domain...
							#area_multiplier = 1
							#particle.surrounding_area = Rect2i(Vector2(particle.position.x - ((grid_domain/2)*area_multiplier),particle.position.y - ((grid_domain/2)*area_multiplier)),Vector2(grid_domain*area_multiplier,grid_domain*area_multiplier))
							#other_particle.surrounding_area = Rect2i(Vector2(snapped(other_particle.position.x - (grid_domain/2)*area_multiplier,1),snapped(other_particle.position.y - (grid_domain/2)*area_multiplier,1)),Vector2(grid_domain*area_multiplier,grid_domain*area_multiplier))
			
						if distance_between.x in range(0.0,(grid_domain)+1) and distance_between.y in range(0.0,-(grid_domain)-1):
							length_covered_area = other_particle.surrounding_area.position.x - particle.surrounding_area.end.x
							width_covered_area = other_particle.surrounding_area.end.y - particle.surrounding_area.position.y
							covered_area = length_covered_area * width_covered_area
							percentage_covered = snapped((float(covered_area) / particle.surrounding_area.get_area()),.01)
							
						elif distance_between.x in range(0.0,(grid_domain)+1) and distance_between.y in range(0.0,(grid_domain)+1):
							length_covered_area = other_particle.surrounding_area.position.x - particle.surrounding_area.end.x
							width_covered_area = other_particle.surrounding_area.position.y - particle.surrounding_area.end.y
							covered_area = length_covered_area * width_covered_area
							percentage_covered = snapped((float(covered_area) / particle.surrounding_area.get_area()),.01)
							
						elif distance_between.x in range(0.0,-(grid_domain)-1) and distance_between.y in range(0.0,(grid_domain)+1):
							length_covered_area = other_particle.surrounding_area.end.x - particle.surrounding_area.position.x
							width_covered_area = other_particle.surrounding_area.position.y - particle.surrounding_area.end.y
							covered_area = length_covered_area * width_covered_area
							percentage_covered = snapped((float(covered_area) / particle.surrounding_area.get_area()),.01)
							
						elif distance_between.x in range(0.0,-(grid_domain)-1) and distance_between.y in range(0.0,-(grid_domain)-1):
							length_covered_area = other_particle.surrounding_area.end.x - particle.surrounding_area.position.x
							width_covered_area = other_particle.surrounding_area.end.y - particle.surrounding_area.position.y
							covered_area = length_covered_area * width_covered_area
							percentage_covered = snapped((float(covered_area) / particle.surrounding_area.get_area()),.01)
							
						if percentage_covered >= .50:
							### Aquire the magnitude and direction of both particles...
							particle.magnitude = snapped(the_grid[particle]['velocity'].length(),.01)
							particle.direction = snapped( rad2deg( snapped( atan2(the_grid[particle]['velocity'].y,the_grid[particle]['velocity'].x),.01)),.01)
							
							other_particle.magnitude = snapped(the_grid[other_particle]['velocity'].length(),.01)
							other_particle.direction = snapped( rad2deg( snapped( atan2(the_grid[other_particle]['velocity'].y,the_grid[other_particle]['velocity'].x),.01)),.01)
							
							### check if both particles is set to be merged...
							if the_grid[particle]['velocity'] == the_grid[other_particle]['velocity'] and particle.magnitude == other_particle.magnitude and particle.direction == other_particle.direction:
								particle.merging = true
								other_particle.merging = true
								
								particle.within_range = []
								other_particle.within_range = []
								
								particles_set_to_merge.append([particle,other_particle])
							else:
								# particle does not have the same magnitude or direction to merge...
								pass
						else:
							### comes in contact but doesnt have amount of coverage to merge...
							pass
					else:
						### no longer in contact with an another particle...
						if particle.within_range.has(other_particle):
							particle.within_range.remove(other_particle)
							#other_particle.within_range.remove(particle)
						else:
							### does not contain the other particle..
							pass
				else:
					#
					##print('particle found itself or the other particle is merging elsewhere.')
					pass
		else:
			##print('the particle is schedule to merge with another particle.')
			pass

		### how the particles interact with the window outline...
		if barriers['coefficient of reinstition'] != 'thru' or barriers['coefficient of reinstition'] != 'disappear':
			for wall in barriers['window outline'].keys():
				if particle.surrounding_area.intersects(barriers['window outline'][wall]):
					###print(' Touching ', wall,' wall')
					### test...
					#particle.velocity = Vector2(0.0,0.0)
					the_grid[particle]['velocity'] = Vector2(0.0,0.0)
					particle.direction = snapped( rad2deg( snapped( atan2(the_grid[particle]['velocity'].y,the_grid[particle]['velocity'].x),.01)),.01)
		else:
			### the particle will travel thru the border...
			# will either continue or disappear from the simulation...
			if barriers['coefficient of reinstition'] == 'thru':
				###print('exist out of sight')
				pass
			elif barriers['coefficient of reinstition'] == 'disappear':
				###print('particles should be removed from the substance...')
				pass
				#grid_nodes.erase(particle.get_name())
					
				#substance.erase(particle)
				#get_tree().get_root().get_node("Test Area/Simulation/substance").remove_child(particle)
			


####...
### particle comes in contact with other particles...
			#if particle.surrounding_area.intersects(other_particle.surrounding_area,true):
			### if the particles come into contact with each other...
			#var contact_with_x = abs((particle.position.x - other_particle.surrounding_area.get_center().x) / grid_domain)
			#var contact_with_y = abs((particle.position.y - other_particle.surrounding_area.get_center().y) / grid_domain)
			var contact_with = (particle.get_position() - other_particle.surrounding_area.get_center())  / grid_domain
			#print(' ')
			#print(particle,' distance from other particle ',other_particle,' is ',contact_with)#*grid_domain)
			#print(grid_domain/grid_domain,' and ',-(grid_domain/grid_domain),' are the particle barrier')
			#print(particle,' relation to ',other_particle,' is ',particle.relation_to_domain[other_particle],' check before update')
				
			#kernel_distance = (particle.get_position() - other_particle.surrounding_area.get_center())# / cell_size
			kernel_x = snapped(particle.get_position().x - other_particle.surrounding_area.get_center().x,.01)
			kernel_y = snapped(particle.get_position().y - other_particle.surrounding_area.get_center().y,.01)
			kernel_distance = Vector2(kernel_x,kernel_y)
			#flipped_kernel_distance = (other_particle.surrounding_area.get_center() - particle.get_position())# / cell_size
			flipped_kernel_x = snapped(other_particle.surrounding_area.get_center().x - particle.get_position().x,.01)
			flipped_kernel_y = snapped(other_particle.surrounding_area.get_center().y - particle.get_position().y,.01)
			flipped_kernel_distance = Vector2(flipped_kernel_x,flipped_kernel_y)
					
				
			#print(kernel_distance,' is the kernel between the particles ',particle,' ',other_particle)
			particle.relation_to_domain[other_particle] = kernel_distance
			particle.domain_relation_to_particle[other_particle] = flipped_kernel_distance
			#particle.relation_to_domain[other_particle] = flipped_kernel_distance
			#particle.domain_relation_to_particle[other_particle] = kernel_distance
			#print(particle,' relation to ',other_particle,' is ',particle.relation_to_domain[other_particle],' check after update')
				
			if  contact_with.x == 0 and contact_with.y == 0:
				### always in contact with itself...
				pass
			elif contact_with.x > -(grid_domain/grid_domain) and contact_with.x < grid_domain/grid_domain and contact_with.y > -(grid_domain/grid_domain) and contact_with.y < grid_domain/grid_domain:
				### come into some contact...
				
				#print(particle,' at position ',particle.position,' in contact with ',other_particle,' at position ',other_particle.surrounding_area,' quadrant')
				if particle.within_range.has(other_particle):
					### the particle is already in contact with the other...
					pass
				else:
					### recognized to be within the domain range...
					particle.within_range.append(other_particle)
			else:
				###
				#print(particle,' not in contact with this particle ',other_particle)
				particle.within_range.erase(other_particle)
				other_particle.within_range.erase(particle)
				
		#print(particle,' has ',len(particle.within_range),' number of particle in contact.')
		
		get_tree().get_root().get_node("Test Area/Simulation/Particle Interaction").Collision_with_Walls(material,barriers,the_grid,window_center,grid_domain)















#################
################# collision with walls

#get_tree().get_root().get_node("Test Area/Simulation/Handle Collisions").Collision_with_Walls(material,barriers,the_grid,window_center,grid_domain)
	




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
						
						
						var wall_center = Vector2(center_of_window.x,0)
						elastic_position_kernel_of_vector_1 = particle.get_position() - wall_center
						elastic_position_kernel_of_vector_2 = wall_center - particle.get_position()
						### find final-velocity 1::
						elastic_velocity_kernel_of_vector_1 = structure[particle]['velocity'] - baluster['window outline'][wall]['velocity']
						elastic_dot_vector_1 = elastic_velocity_kernel_of_vector_1.dot(elastic_position_kernel_of_vector_1)
						### find final-velocit 2::
						elastic_velocity_kernel_of_vector_2  = baluster['window outline'][wall]['velocity'] - structure[particle]['velocity']
						elastic_dot_vector_2 = elastic_velocity_kernel_of_vector_2.dot(elastic_position_kernel_of_vector_2)
					
						elastic_magnitude_kernel_of_vector_1 = sqrt( (pow(elastic_position_kernel_of_vector_1.x,2.0)+pow(elastic_position_kernel_of_vector_1.y,2.0)) )
						elastic_magnitude_kernel_of_vector_2 = sqrt( (pow(elastic_position_kernel_of_vector_2.x,2.0)+pow(elastic_position_kernel_of_vector_2.y,2.0)) )
						
						print(elastic_term_1_of_vector_1,' 2 * wall mass')
						print(elastic_dot_vector_1,' dot product')
						print(elastic_position_kernel_of_vector_1,' position kernel')
						print(elastic_mass_addition,' particle mass + wall mass')
						print(elastic_magnitude_kernel_of_vector_1,' magnitude of kernel')
						structure[particle]['velocity'] = structure[particle]['velocity'] - ( elastic_term_1_of_vector_1 / elastic_mass_addition ) * ( elastic_dot_vector_1 / elastic_magnitude_kernel_of_vector_1 ) * ( elastic_position_kernel_of_vector_1 )
						print(structure[particle]['velocity'],' after collision')
						
						
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
						
						
						var wall_center = Vector2(center_of_window.x,ProjectSettings.get_setting('display/window/size/height'))
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
						
						
						print(elastic_term_1_of_vector_1,' 2 * wall mass')
						print(elastic_dot_vector_1,' dot product ',elastic_velocity_kernel_of_vector_1,' velocity_kernel dot to position_kernel ',elastic_position_kernel_of_vector_1)
						
						print(elastic_position_kernel_of_vector_1,' position kernel',' ',particle.get_position(),' particle position ',' - ',wall_center,' wall center')
						print(elastic_mass_addition,' particle mass + wall mass')
						print(elastic_magnitude_kernel_of_vector_1,' magnitude of kernel')
						structure[particle]['velocity'] = structure[particle]['velocity'] -( ( elastic_term_1_of_vector_1 / elastic_mass_addition ) * ( elastic_dot_vector_1 / elastic_magnitude_kernel_of_vector_1 ) * ( elastic_position_kernel_of_vector_1 ) )
						print(structure[particle]['velocity'],' after collision')
						




					if baluster['window outline'][wall]['coefficient of reinstition'] == 2.0:
						### Explosion ::
						# 
						# mass_1_i*velocity_1_i + mass_2_i*velocity_2_i = mass_1_f*velocity_1_f + mass_2_f*velocity_2_f
						# if mass_1_i*velocity_1_i is at rest: 0 = mass_1_f *velocity_1_f + mass_2_f*velocity_2_f
						#                                    : -(mass_1_f*velocity_1_f) = mass_2_f*velocity_2_f
						#                                    :  -(mass_1_f*velocity_1_f) / mass_2_f = velocity_2_f
						#
						#
						#
						#
						pass
						
					el












#for other_particle in particle.within_range:
			
		#	kernel_distance = (particle.get_position() - other_particle.surrounding_area.position) / grid_domain
		#	flipped_kernel_distance = (other_particle.surrounding_area.position - particle.get_position()) / grid_domain
			
			### Weight Interpolation...
		#	Weight_Interpolation(kernel_distance,grid_domain)
		#	particle.velocity += (weight_interpolation * the_grid[other_particle]['velocity'])
		
#print(weight_interpolation,' during grid to particle.')
			#particle.velocity = particle.velocity + (weight_interpolation * the_grid[other_particle]['velocity'])
			#await B_formation(particle,weight_interpolation,the_grid[other_particle]['velocity'])
			#print(particle.velocity,' the particle velocity')
			#particle.B = particle.velocity * flipped_kernerl_distance^Transposed
			#construct_B = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Vector2_by_Vector2_to_Matrix(particle.velocity,false,flipped_kernel_distance,true)
			#print(construct_B,' build B')
			#particle.B = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Add_Matrix(particle.B,construct_B)
			#print(particle.B,' B')
		
		
		#var c_term = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Matrix_by_Scalar(the_grid['velocity'],4,true)
			#var c_term_2 = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Divide_Matrix_by_Scalar(c_term,pow(grid_domain,2),true)
			#var c_term_3 = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Inverse_Matrix(particle.I)
			#particle.C = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Matrix(c_term_3,c_term_2)
		
		
		
		
		
		
		
		#C Updated..
		#particle.C = particle.B * D^-1  :: D^-1 = 4/grid_domain^2 * particle.I^-1
		#var c_term = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Matrix_by_Scalar(particle.B,4,true)
		#var c_term_2 = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Divide_Matrix_by_Scalar(c_term,pow(grid_domain,2),true)
		#var c_term_3 = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Inverse_Matrix(particle.I)
		#particle.C = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Matrix(c_term_3,c_term_2)
		#print(particle.C,' C')
		#particle.C = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Matrix_by_Vector2_to_Matrix(c_term_4,flipped_kernel_distance,)
		
#########################################################################
########################################################################
######################################   MLSMPM- version -still usable....


######### from p2g
#MLS-MPM : 
			#shape_function = weight_interpolation * D^-1 * flipped_kernel_distance^T * flipped_kernel_distance
			#D^1 = (4/pow(cell_size,2)) * other_particle.I^-1 
			var shape_weight_coefficent = (other_particle.mass * weight_interpolation  * 4.0)/pow(grid_domain,2.0)
			var shape_inverse_I_coefficent = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Inverse_Matrix(other_particle.I)
			var shape_weighted_I_coefficent = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Matrix_by_Scalar(shape_inverse_I_coefficent,shape_weight_coefficent,true)
			var shape_term_1 = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Matrix_by_Vector2_to_Matrix(shape_weighted_I_coefficent,particle.domain_relation_to_particle[other_particle],true,true)
			var shape_term_2 = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Matrix_by_Vector2_to_Matrix(shape_term_1,particle.domain_relation_to_particle[other_particle],false,false)
			var shape_function = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Convert_to_Scalar(shape_term_2)
			#print(shape_function,' shape function..with weights.')
			#print(shape_function,' shape function..with weights.')
			
			#the_grid[particle]['mass'] = the_grid[particle]['mass'] + (other_particle.mass * shape_function)
			#the_grid[particle]['mass'] = the_grid[other_particle]['mass'] + shape_function
			
			
		
#### during g2p
var shape_weight_coefficent = (the_grid[other_particle]['velocity'] * weight_interpolation  * 4.0)/pow(grid_domain,2.0)
			var shape_inverse_I_coefficent = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Inverse_Matrix(other_particle.I)
			var shape_weighted_I_coefficent = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Matrix_by_Vector2_to_Matrix(shape_inverse_I_coefficent,shape_weight_coefficent,false,false)
			var shape_term_1 = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Matrix_by_Vector2_to_Matrix(shape_weighted_I_coefficent,particle.domain_relation_to_particle[other_particle],true,false)
			var shape_function = get_tree().get_root().get_node("Test Area/Simulation/Matrix Math").Multiply_Matrix_by_Vector2_to_Vector2(shape_term_1,particle.domain_relation_to_particle[other_particle])
			
			#particle.velocity = particle.velocity + (the_grid[other_particle]['velocity'] * shape_function)
			#particle.velocity = particle.velocity + shape_function
			








#####################
#####################
######   MPM-verison


#"""
	### Collision Detection...
	### how the particles interact with the window outline...
	for particle in material:
		if barriers['coefficient of reinstition'] != 'thru' or barriers['coefficient of reinstition'] != 'disappear':
			for wall in barriers['window outline'].keys():
				if particle.surrounding_area.intersects(barriers['window outline'][wall]):
					###print(' Touching ', wall,' wall')
					### test...
					#particle.velocity = Vector2(0.0,0.0)
					the_grid[particle]['velocity'] = Vector2(0.0,0.0)
					particle.velocity = Vector2(0.0,0.0)
					particle.B = [[0.0,0.0],[0.0,0.0]]
					particle.direction = snapped( rad2deg( snapped( atan2(the_grid[particle]['velocity'].y,the_grid[particle]['velocity'].x),.01)),.01)
		else:
			### the particle will travel thru the border...
			# will either continue or disappear from the simulation...
			if barriers['coefficient of reinstition'] == 'thru':
				###print('exist out of sight')
				pass
			elif barriers['coefficient of reinstition'] == 'disappear':
				###print('particles should be removed from the substance...')
				pass
				#grid_nodes.erase(particle.get_name())
					
				#substance.erase(particle)
				#get_tree().get_root().get_node("Test Area/Simulation/substance").remove_child(particle)
			
	"""
	for particle in material:
		if particle.merging == false:
			for other_particle in particle.within_range:
				if particle != other_particle and other_particle.merging == false:
					
					distance_between = particle.get_position() - other_particle.get_position()
					print(distance_between, ' checking collisions')
					if distance_between.x in range(-(grid_domain),(grid_domain)+1) and distance_between.y in range(-(grid_domain),(grid_domain)+1):
						
						if particle.within_range.has(other_particle):
							pass
						else:
							### recognized to come into contact.
							particle.within_range.append(other_particle)
							# establish the particle grid domain...
							#area_multiplier = 1
							#particle.surrounding_area = Rect2i(Vector2(particle.position.x - ((grid_domain/2)*area_multiplier),particle.position.y - ((grid_domain/2)*area_multiplier)),Vector2(grid_domain*area_multiplier,grid_domain*area_multiplier))
							#other_particle.surrounding_area = Rect2i(Vector2(snapped(other_particle.position.x - (grid_domain/2)*area_multiplier,1),snapped(other_particle.position.y - (grid_domain/2)*area_multiplier,1)),Vector2(grid_domain*area_multiplier,grid_domain*area_multiplier))
			
						if distance_between.x in range(0.0,(grid_domain)+1) and distance_between.y in range(0.0,-(grid_domain)-1):
							length_covered_area = other_particle.surrounding_area.position.x - particle.surrounding_area.end.x
							width_covered_area = other_particle.surrounding_area.end.y - particle.surrounding_area.position.y
							covered_area = length_covered_area * width_covered_area
							percentage_covered = snapped((float(covered_area) / particle.surrounding_area.get_area()),.01)
							
						elif distance_between.x in range(0.0,(grid_domain)+1) and distance_between.y in range(0.0,(grid_domain)+1):
							length_covered_area = other_particle.surrounding_area.position.x - particle.surrounding_area.end.x
							width_covered_area = other_particle.surrounding_area.position.y - particle.surrounding_area.end.y
							covered_area = length_covered_area * width_covered_area
							percentage_covered = snapped((float(covered_area) / particle.surrounding_area.get_area()),.01)
							
						elif distance_between.x in range(0.0,-(grid_domain)-1) and distance_between.y in range(0.0,(grid_domain)+1):
							length_covered_area = other_particle.surrounding_area.end.x - particle.surrounding_area.position.x
							width_covered_area = other_particle.surrounding_area.position.y - particle.surrounding_area.end.y
							covered_area = length_covered_area * width_covered_area
							percentage_covered = snapped((float(covered_area) / particle.surrounding_area.get_area()),.01)
							
						elif distance_between.x in range(0.0,-(grid_domain)-1) and distance_between.y in range(0.0,-(grid_domain)-1):
							length_covered_area = other_particle.surrounding_area.end.x - particle.surrounding_area.position.x
							width_covered_area = other_particle.surrounding_area.end.y - particle.surrounding_area.position.y
							covered_area = length_covered_area * width_covered_area
							percentage_covered = snapped((float(covered_area) / particle.surrounding_area.get_area()),.01)
							
						if percentage_covered >= .50:
							### Aquire the magnitude and direction of both particles...
							particle.magnitude = snapped(the_grid[particle]['velocity'].length(),.01)
							particle.direction = snapped( rad2deg( snapped( atan2(the_grid[particle]['velocity'].y,the_grid[particle]['velocity'].x),.01)),.01)
							
							other_particle.magnitude = snapped(the_grid[other_particle]['velocity'].length(),.01)
							other_particle.direction = snapped( rad2deg( snapped( atan2(the_grid[other_particle]['velocity'].y,the_grid[other_particle]['velocity'].x),.01)),.01)
							
							### check if both particles is set to be merged...
							if the_grid[particle]['velocity'] == the_grid[other_particle]['velocity'] and particle.magnitude == other_particle.magnitude and particle.direction == other_particle.direction:
								particle.merging = true
								other_particle.merging = true
								
								particle.within_range = []
								other_particle.within_range = []
								
								particles_set_to_merge.append([particle,other_particle])
							else:
								# particle does not have the same magnitude or direction to merge...
								pass
						else:
							### comes in contact but doesnt have amount of coverage to merge...
							pass
					else:
						### no longer in contact with an another particle...
						if particle.within_range.has(other_particle):
							particle.within_range.remove(other_particle)
							#other_particle.within_range.remove(particle)
						else:
							### does not contain the other particle..
							pass
				else:
					#
					##print('particle found itself or the other particle is merging elsewhere.')
					pass
		else:
			##print('the particle is schedule to merge with another particle.')
			pass

		### how the particles interact with the window outline...
		if barriers['coefficient of reinstition'] != 'thru' or barriers['coefficient of reinstition'] != 'disappear':
			for wall in barriers['window outline'].keys():
				if particle.surrounding_area.intersects(barriers['window outline'][wall]):
					###print(' Touching ', wall,' wall')
					### test...
					#particle.velocity = Vector2(0.0,0.0)
					the_grid[particle]['velocity'] = Vector2(0.0,0.0)
					particle.direction = snapped( rad2deg( snapped( atan2(the_grid[particle]['velocity'].y,the_grid[particle]['velocity'].x),.01)),.01)
		else:
			### the particle will travel thru the border...
			# will either continue or disappear from the simulation...
			if barriers['coefficient of reinstition'] == 'thru':
				###print('exist out of sight')
				pass
			elif barriers['coefficient of reinstition'] == 'disappear':
				###print('particles should be removed from the substance...')
				pass
				#grid_nodes.erase(particle.get_name())
					
				#substance.erase(particle)
				#get_tree().get_root().get_node("Test Area/Simulation/substance").remove_child(particle)
			
	"""




