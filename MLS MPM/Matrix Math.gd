extends Node


var new_matrix : Array = []
var check_matrix : String
var matrix_toolbox  : Dictionary
var new_vector : Vector2
var new_scalar : float
var converted_vector : Array
var new_determinant : float

func _on_Matrix_Math_ready():
	### Collection of Math Formulas:(Matrix) used/needed...
	pass


func Double_Distribute_Property(a,b,c,d):
	var property_ac
	var property_ad
	var property_bc
	var property_bd
	var double_distribute_result
	### (a %s b)(c %s d); %s can be +,-,*,/ 
	# (a * c) + (a * d) +(b * c) + (b * d)
	double_distribute_result = (a * c) + (a * d) +(b * c) + (b * d)
	if typeof(a) == TYPE_STRING or typeof(b) == TYPE_STRING or typeof(c) == TYPE_STRING or typeof(d) == TYPE_STRING:
		### contains an variable of an unknown...
		a = String(a)
		b = String(b)
		c = String(c)
		d = String(d)
		
		# (a * c) + (a * d) +(b * c) + (b * d)
		if a.match(c):
			### both contain the same variable
			# so it becomes variable squared...
			
			property_ac = a + c
			property_ad = a + d
			property_bc = b + c
			property_bd = b + d
			
			print(property_ac,' ac')
			print(property_ad,' ad')
			print(property_bc,' bc')
			print(property_bd,' bd')
	else:
		### each is a known number ...
		double_distribute_result = (a * c) + (a * d) + (b * c) + (b * d)
	
	return double_distribute_result


func Sum_Of(cumulate):
	return cumulate
	
func Make_Matrix(m:int,n:int):
	### Construct an Matrix M x N....
	var placeholder : Array = []
	new_matrix = []
	
	for row in range(m):
		for column in range(n):
			placeholder.append(0)
		
		new_matrix.append(placeholder)
		placeholder = []
	
	# a template of the matrix is stored...
	check_matrix = '{m} by {n}'.format({'m':m,'n':n})
	matrix_toolbox[check_matrix] = new_matrix
	
	return new_matrix

func Identify_Rows(m:Array):
	### Reveals the number of rows of the matrix...
	return len(m)

func Identify_Columns(m:Array):
	### Reveals the number of columns of the matrix...
	return len(m[0])

func Convert_to_Vector2(m:Array):
	#Converts a Matrix in to a Vector2 """
	
	if len(m)== 2 and len(m[0]) == 1:
		new_vector = Vector2(m[0][0],m[1][0])
	else:
		#print('converts a 2 by 1 matrix into a vector...')
		pass
	return new_vector

func Convert_to_Scalar(m:Array):
	#Converts a 1x1 matrix into a scalar."
	new_scalar = 0.0
	
	if len(m) == 1 and len(m[0]) == 1:
		
		new_scalar = m[0][0]
	
	return new_scalar

func Convert_Vector_2_to_Matrix(n:Vector2,transpose:bool):
	#Converts a Vector2 into a  2x1-Matrix"
	
	if transpose == false:
		check_matrix = "{m} by {n}".format({'m':2,'n':1})
	
		### storing a template of a created matrix...
		if matrix_toolbox.has(check_matrix) == true:
			### Matrix templated exists...
			new_matrix = matrix_toolbox[check_matrix].duplicate(true)
		else:
			### create a placeholding matrix.
			Make_Matrix(2,1)
		new_matrix = [[n.x],[n.y]]
	else:
		check_matrix = "{m} by {n}".format({'m':1,'n':2})
	
		### storing a template of a created matrix...
		if matrix_toolbox.has(check_matrix) == true:
			### Matrix templated exists...
			new_matrix = matrix_toolbox[check_matrix].duplicate(true)
		else:
			### create a placeholding matrix.
			Make_Matrix(1,2)
		
		new_matrix = [[n.x,n.y]]
	
	return new_matrix
	

func Convert_Vector_3_to_Matrix(n:Vector3,transpose:bool):
	#Converts a Vector2 into a  2x1-Matrix"
	
	if transpose == false:
		check_matrix = "{m} by {n}".format({'m':3,'n':1})
	
		### storing a template of a created matrix...
		if matrix_toolbox.has(check_matrix) == true:
			### Matrix templated exists...
			new_matrix = matrix_toolbox[check_matrix].duplicate(true)
		else:
			### create a placeholding matrix.
			Make_Matrix(3,1)
		new_matrix = [[n.x],[n.y],[n.z]]
	else:
		check_matrix = "{m} by {n}".format({'m':1,'n':3})
	
		### storing a template of a created matrix...
		if matrix_toolbox.has(check_matrix) == true:
			### Matrix templated exists...
			new_matrix = matrix_toolbox[check_matrix].duplicate(true)
		else:
			### create a placeholding matrix.
			Make_Matrix(1,3)
		
		new_matrix = [[n.x,n.y,n.z]]
	
	return new_matrix

func Trace(m:Array):
	#The trace of a square matrix is the sum of its (main)-diagonal elements. "
	
	new_scalar = 0.0
	
	if len(m) == len(m[0]):
		### checks if the matrix is a square matrix...
		
		for square in range(0, len(m[0])):
			
			new_scalar = snapped(new_scalar + m[square][square],.001)
	else:
		#print('not a square matrix... m by n where m equals n')
		pass
	return new_scalar

func Find_Determinant(m:Array):
	### a number calculated from a square matrix( rows == columns ):
	if len(m) == len(m[0]):
		
		if len(m) == 2:
			### if the square-matrix is 2...
			# | Determinant | = (a*d) - (b*c)
			new_determinant = snapped(((m[0][0] * m[1][1]) - (m[0][1] * m[1][0])),.01)
		elif len(m) == 3:
			pass
		elif len(m) == 4:
			pass
		pass
	else:
		###
		print('Not a square matrix...')
	
	return new_determinant
	
	
func Transposed_Matrix(m:Array):
	#a transposed matrix is a flipped verison of the original matrix... "
	#To flip the matrix between its rows and columns..."
	new_matrix = []
	
	check_matrix = "{m} by {n}".format({'m':len(m[0]),'n':len(m)})
	
	### storing a template of a created matrix...
	if matrix_toolbox.has(check_matrix) == true:
		### Matrix templated exists...
		new_matrix = matrix_toolbox[check_matrix].duplicate(true)
	
	else:
		### create a placeholding matrix.
		Make_Matrix(len(m[0]),len(m))
	
	for row in range(0,len(m)):
		for column in range(0,len(m[row])):
			new_matrix[column][row] = m[row][column]
	
	return new_matrix


func Inverse_Matrix(m:Array):
	#An Inverse the matrix must be "square" (same number of rows and columns).
	#and also the determinant cannot be zero (or we end up dividing by zero)'
	# matrix to the -1 power...
	new_scalar = 0.0
	var placeholder2x2_1 = 0.0
	var placeholder2x2_2 = 0.0
	var placeholder2x2_3 = 0.0
	var placeholder2x2_4 = 0.0
	var determinant_scalar = 0.0
	
	
	if len(m[0]) == len(m):
		### checks if the matrix is a square matrix...
		if len(m[0]) == 2:
			### a 2x2 square matrix...
			
			#new_scalar = stepify( float( ( m[0][0]*m[1][1] ) ) - float( ( m[0][1]*m[1][0] ) ),.01)
			check_matrix = "{m} by {n}".format({'m':len(m),'n':len(m[0])})
	
			### storing a template of a created matrix...
			if matrix_toolbox.has(check_matrix) == true:
				### Matrix templated exists...
				new_matrix = matrix_toolbox[check_matrix].duplicate(true)
			else:
				### create a placeholding matrix.
				Make_Matrix(len(m),len(m[0]))
			
			### place holding the matrix info...
			placeholder2x2_1 = m[0][0]
			placeholder2x2_2 = m[0][1]
			placeholder2x2_3 = m[1][0]
			placeholder2x2_4 = m[1][1]
			
			new_matrix[0][0] = float(placeholder2x2_4)
			if placeholder2x2_2 >= 1:
				### the number is positive...
				new_matrix[0][1] = -float(placeholder2x2_2)
			elif placeholder2x2_2 <= -1:
				### the number is negative...
				new_matrix[0][1] = float(placeholder2x2_2)
			else:
				### number is 0
				new_matrix[0][1] = float(placeholder2x2_2)
			if placeholder2x2_3 >= 1:
				### the number is positive...
				new_matrix[0][1] = -float(placeholder2x2_3)
			elif placeholder2x2_3 <= -1:
				### the number is negative...
				new_matrix[0][1] = float(placeholder2x2_3)
			else:
				### number is 0
				new_matrix[0][1] = float(placeholder2x2_3)
			new_matrix[1][1] = float(placeholder2x2_1)
			
			if (new_matrix[0][0] * new_matrix[1][1]) - (new_matrix[0][1] * new_matrix[1][0]) > 0:
				
				determinant_scalar = 1.0 / (new_matrix[0][0] * new_matrix[1][1]) - (new_matrix[0][1] * new_matrix[1][0])
				Multiply_Matrix_by_Scalar(new_matrix,determinant_scalar,true)
			
			else:
				#print('The matrix does not have an Inverse.')
				pass
		elif len(m[0]) == 3:
			### a 3x3 square matrix...
			pass
		elif len(m[0]) == 4:
			### a 4x4 square matrix...
			pass
	else:
		#print('to find the inverse matrix the original matrix has to be square')
		#print('m by n where m = n')
		pass
	return new_matrix

#func Find_Eigenvectors(a:Array):
	### Eigenvector formula : AV = evV
	# A = square matrix, ev = eigenvalues, V = eigenvectors
	### eigenvalues formula : | A - Iev |...
	# A = sqare matrix , I = Identity Matrix of A,ev = eigenvalues
#	var eigenvectors_01 : Vector2
	#var eigenvectors_02 : Vector2
	#var solve_by_determinant : Array
	#var find_eigenvector : Array
	#var eigen_identity : Array
	#var quadratic_form 
	#var eigenvalues_01 : float
	#var eigenvalues_02 : float
	
#	var quadratic_coefficient_a
	#var quadratic_coefficient_b
	#var quadratic_constant_c
	
	#var helper_diagonalize_matrix

	#if Identify_Rows(a) == Identify_Columns(a):
		### It is a square matrix...
		#if Identify_Rows(a) == 2:
			
			### Find Eigenvalues...
			# I time v
			# A minus Iv matrix
			# determinant from A minus Iv 
			#eigen_identity = [['ev',0],[0,'ev']]
			#solve_by_determinant = [['{0}'.format({'0':a[0][0]}) + ' - ' + eigen_identity[0][0],a[0][1] - eigen_identity[0][1]],[a[1][0] - eigen_identity[1][0],'{0}'.format({'0':a[1][1]}) + ' - ' + eigen_identity[1][1] ]]
			
			#print(solve_by_determinant,' solving')
			
		#	quadratic_coefficient_a = - 1.0 * -1.0
			#quadratic_coefficient_b = (a[0][0] * -1.0) + (a[1][1] * -1.0)
			#quadratic_constant_c = (a[0][0] * a[1][1]) + (solve_by_determinant[0][1] * solve_by_determinant[1][0])
			
		
			#solving quadratic equation
			# x = -b +- square root ( power(b,2.0) - 4ac  ) / 2a
			
			#eigenvalues_01 = (-1.0*quadratic_coefficient_b) + ( (sqrt(  pow(quadratic_coefficient_b,2.0) - (4.0 *quadratic_coefficient_a*quadratic_constant_c) )) / (2.0 * quadratic_coefficient_a) )
			
			#eigenvalues_02 = (-1.0*quadratic_coefficient_b) - ( (sqrt(  pow(quadratic_coefficient_b,2.0) - (4.0 *quadratic_coefficient_a*quadratic_constant_c) )) / (2.0 * quadratic_coefficient_a) )
	
			#print(eigenvalues_01,' plus')
			#print(eigenvalues_02,' minus')

			# Find Eigenvectors...
			#eigenvectors_01.x = (a[0][0] - eigenvalues_01) / (-1.0 * (a[0][1] - eigenvalues_01))
			#eigenvectors_01.y = (a[1][0] - eigenvalues_01) / (-1.0 * (a[1][1] - eigenvalues_01))
			
			#eigenvectors_02.x = (a[0][0] - eigenvalues_02) / (-1.0 * (a[0][1] - eigenvalues_02))
			#eigenvectors_02.y = (a[1][0] - eigenvalues_02) / (-1.0 * (a[1][1] - eigenvalues_02))
			
			# Eigenvectors form a mastrix...
			#helper_diagonalize_matrix = [[eigenvectors_01.x,eigenvectors_02.x ],[eigenvectors_01.y,eigenvectors_02.y ]]
			
			#return helper_diagonalize_matrix
#	else:
		
		### not a square matrix...
	#	print('not a square matrix')
	


func Add_Matrix(m:Array,n:Array):
	#Add two matrices if they are the same shape(rows) and size(columns)."""
	new_matrix = []

	if len(m) == len(n):
		### rows check...
		if len(m[0]) == len(n[0]):
			### columns check...
			
			check_matrix = "{m} by {n}".format({'m':len(m),'n':len(m[0])})
	
			### storing a template of a created matrix...
			if matrix_toolbox.has(check_matrix) == true:
				### Matrix templated exists...
				new_matrix = matrix_toolbox[check_matrix].duplicate(true)
			else:
				### create a placeholding matrix.
				Make_Matrix(len(m),len(m[0]))
				
			### The Math happens...
			for shape in range(0,len(m)):
				for size in range(0,len(m[shape])):
					#new_matrix[shape][size] = stepify(m[shape][size] + n[shape][size],.001)
					new_matrix[shape][size] = m[shape][size] + n[shape][size]
			
			return new_matrix
		else:
			#print('not the same size')
			pass
	
	else:
		#print('not the same shape')
		pass

func Subtract_Matrix(m:Array,n:Array):
	#Add two matrices if they are the same shape(rows) and size(columns)."""
	new_matrix = []

	if len(m) == len(n):
		### rows check...
		if len(m[0]) == len(n[0]):
			### columns check...
			
			check_matrix = "{m} by {n}".format({'m':len(m),'n':len(m[0])})
	
			### storing a template of a created matrix...
			if matrix_toolbox.has(check_matrix) == true:
				### Matrix templated exists...
				new_matrix = matrix_toolbox[check_matrix].duplicate(true)
			else:
				### create a placeholding matrix.
				Make_Matrix(len(m),len(m[0]))
				
			### The Math happens...
			for shape in range(0,len(m)):
				for size in range(0,len(m[shape])):
					#new_matrix[shape][size] = stepify(m[shape][size] - n[shape][size],.001)
					new_matrix[shape][size] = m[shape][size] - n[shape][size]
			
			return new_matrix
		else:
			#print('not the same size')
			pass
	
	else:
		#print('not the same shape')
		pass


func Add_Matrix_by_Scalar(m:Array,n,is_float:bool):
	#you add a matrix (ab) by a (n) number, 
	#you add every element in the matrix (ab) by the (n) number """
	if is_float == true:
		n = float(n)
	else:
		n = int(n)
		
	check_matrix = "{m} by {n}".format({'m':len(m),'n':len(m[0])})
	
	### storing a template of a created matrix...
	if matrix_toolbox.has(check_matrix) == true:
		### Matrix templated exists...
		new_matrix = matrix_toolbox[check_matrix].duplicate(true)
	else:
		### create a placeholding matrix.
		
		Make_Matrix(len(m),len(m[0]))
	
	### The Math happens...
	for shape in range(0,len(m)):
		for size in range(0,len(m[shape])):
			#new_matrix[shape][size] = stepify(m[shape][size] + n,.001)
			new_matrix[shape][size] = m[shape][size] + n
			
	return new_matrix

func Multiply_Matrix_by_Scalar(m:Array,n,is_float:bool):
	#you multiply a matrix (ab) by a (n) number, 
	#you multiply every element in the matrix (ab) by the (n) number """
	if is_float == true:
		n = float(n)
	else:
		n = int(n)
		
	check_matrix = "{m} by {n}".format({'m':len(m),'n':len(m[0])})
	
	### storing a template of a created matrix...
	if matrix_toolbox.has(check_matrix) == true:
		### Matrix templated exists...
		new_matrix = matrix_toolbox[check_matrix].duplicate(true)
	else:
		### create a placeholding matrix.
		
		Make_Matrix(len(m),len(m[0]))
				
	
	### The Math happens...
	for shape in range(0,len(m)):
		for size in range(0,len(m[shape])):
			#new_matrix[shape][size] = stepify(m[shape][size] * n,.001)
			new_matrix[shape][size] = snapped((m[shape][size] * n),.001)
			
	return new_matrix


func Divide_Matrix_by_Scalar(m:Array,n,is_float:bool):
	#you divide a matrix (ab) by a (n) number, 
	#you divide every element in the matrix (ab) by the (n) number """
	if is_float == true:
		n = float(n)
	else:
		n = int(n)
		
	check_matrix = "{m} by {n}".format({'m':len(m),'n':len(m[0])})
	
	### storing a template of a created matrix...
	if matrix_toolbox.has(check_matrix) == true:
		### Matrix templated exists...
		new_matrix = matrix_toolbox[check_matrix].duplicate(true)
	else:
		### create a placeholding matrix.
		
		Make_Matrix(len(m),len(m[0]))
				

	### The Math happens...
	for shape in range(0,len(m)):
		for size in range(0,len(m[shape])):
			#new_matrix[shape][size] = stepify(m[shape][size] / n,.001)
			new_matrix[shape][size] = snapped((m[shape][size] / n),.001)
			
	return new_matrix

func Multiply_Matrix(m:Array,n:Array):
	#In order to multiply two matrices A and B to get AB the number of
	#columns of A must equal the number of rows of B"""
	#The Matrix AB size will be the row size of A by column size of B

	#new_matrix = []
	
	if len(m[0]) == len(n):
		
		check_matrix = "{m} by {n}".format({'m':len(m),'n':len(n[0])})
		
		### storing a template of a created matrix...
		if matrix_toolbox.has(check_matrix) == true:
			### Matrix templated exists...
			new_matrix = matrix_toolbox[check_matrix].duplicate(true)
		else:
			### create a placeholding matrix.
			
			Make_Matrix(len(m),len(n[0]))
		
		#### The Math happens
		for product_row in range(0,len(new_matrix)):
			for product_column in range(0,len(new_matrix[0])):
				for add_index in range(0,len(n)):
					
					new_matrix[product_row][product_column] = snapped(new_matrix[product_row][product_column] + snapped((m[product_row][add_index] * n[add_index][product_column]),.001),.001)
		
		# Dot Product...
		# the resulted matrix will be the number of rows by the number of columns...
		# the row will pair off the corresponding column...
		# the index of the row is multiplied to each index of the column
		# each of the resulted index are added together...
		# the result from that added-index will go in the resulted matrix[corresponding_row][corresponding_column] 
		# the row will procede to the next column and repeat the process until no more columns.
		# the next row will follow the same procedure as the previous...
		
		return new_matrix
		
	else:
		#print('Columns of A must equal the number of rows of B
		#	The Matrix AB size will be the row size of A by column size of B')
		pass

func Multiply_Matrix_by_Vector2_to_Vector2(m:Array,v:Vector2):
	# a matrix is multipled by a vector2 and a vector is return.

	### matrices of any sizes...
	converted_vector = [[v.x],[v.y]]
	
	Multiply_Matrix(m,converted_vector)
	
	new_vector = Vector2(new_matrix[0][0],new_matrix[1][0])
	
	return new_vector

func Multiply_Matrix_by_Vector2_to_Matrix(m:Array,v:Vector2,vector_transposed:bool,flip:bool):
	#A Matrix multiplied to a Vector , the vector is transformed into a 2x1-Matrix " 
	#Result into a New Matrix..."'''

	if vector_transposed == true:
		converted_vector = [[v.x,v.y]]
	else:
		### matrices of any sizes...
		converted_vector = [[v.x],[v.y]]
	
	if flip == false:
		Multiply_Matrix(m,converted_vector)
	elif flip == true:
		Multiply_Matrix(converted_vector,m)
		
	return new_matrix

func Multiply_Vector2_by_Vector2_to_Matrix(m:Vector2,m_transposed:bool,n:Vector2,n_transposed:bool):
	# Two vectors are turned to matrices then multiplied...
	# Not possible to multiply to vectors after turning into matrix unless 1 of the matrix is transposed...
	
	var matrix_a 
	var matrix_b
	
	#both vectors is converted to a Matrix.
	if m_transposed == true:
		matrix_a = Convert_Vector_2_to_Matrix(m,true).duplicate(true)
	else:
		matrix_a = Convert_Vector_2_to_Matrix(m,false).duplicate(true)
	if n_transposed == true:
		matrix_b = Convert_Vector_2_to_Matrix(n,true).duplicate(true)
	else:
		matrix_b = Convert_Vector_2_to_Matrix(n,false).duplicate(true)
	
	# find if the columns of A and the rows of B are the same...
	if len(matrix_a[0]) == len(matrix_b):
		Multiply_Matrix(matrix_a,matrix_b)
	elif len(matrix_a) == len(matrix_b[0]):
		Multiply_Matrix(matrix_b,matrix_a)
	
	return new_matrix
