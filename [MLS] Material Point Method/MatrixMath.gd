extends Node


var new_matrix : Array = []
var check_matrix : String
var matrix_toolbox  : Dictionary
var matrix_row_element : int
var matrix_column_element :int
var dot_index :int
var new_vector : Vector2
var new_scalar : float
var converted_vector : Array
var new_determinant : float


func _on_MatrixMath_ready():
	### Collection of Math Formulas:(Matrix) used/needed...
	pass


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
	
	#if len(m)== 2 and len(m[0]) == 1:
		#new_vector = Vector2(m[0][0],m[1][0])
	#else:
		#print('converts a 2 by 1 matrix into a vector...')
		#pass
	if len(m) == 2:
		new_vector = Vector2(m[0],m[1])
	#elif len(m) == 4:
		#pass
	return new_vector

func Convert_to_Scalar(m:Array):
	#Converts a 1x1 matrix into a scalar."
	new_scalar = 0.0
	
	if len(m) == 1 and len(m[0]) == 1:
		
		new_scalar = m[0][0]
	
	return new_scalar

func Convert_Vector_2_to_Matrix(_n:Vector2,_transpose:bool):
	#Converts a Vector2 into a  2x1-Matrix"
	pass
	#new_matrix = [[n.x,n.y]][[n.x],[n.y]]
	
	return new_matrix
	

func Convert_Vector_3_to_Matrix(n:Vector3,transpose:bool):
	#Converts a Vector2 into a  2x1-Matrix"...
	
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
	
	if len(m) == 4:
		new_scalar = m[1] + m[3]
	
	#if len(m) == len(m[0]):
	
		### checks if the matrix is a square matrix...
		
		#for square in range(0, len(m[0])):
			
			#new_scalar = snapped(new_scalar + m[square][square],.001)
	#else:
		#print('not a square matrix... m by n where m equals n')
		#pass
	return new_scalar

func Find_Determinant(m:Array):

	### a number calculated from a square matrix( rows == columns ):
	#if len(m) == len(m[0]):
	if len(m) == 4:

		new_determinant = abs(snapped(((m[0] * m[3]) - (m[1] * m[2])),.01))

	
	return new_determinant
	
	
func Transposed_Matrix(m:Array):
	#a transposed matrix is a flipped verison of the original matrix... "
	#To flip the matrix between its rows and columns..."
	
	#new_matrix = [[0,0],[0,0]]
	#new_matrix[0][0] = m[0][0]
	#new_matrix[0][1] = m[1][0]
	#new_matrix[1][0] = m[0][1]
	#new_matrix[1][1] = m[1][1]
	if len(m) == 2:
		new_matrix = [0,0]
		new_matrix[0] = m[1]
		new_matrix[1] = m[0]
	elif len(m) == 4:
		new_matrix = [0,0,0,0]
		new_matrix[0] = m[0]
		new_matrix[1] = m[2]
		new_matrix[2] = m[1]
		new_matrix[3] = m[3]
	
	
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
	
	var updated_new_matrix
	
	#if len(m[0]) == len(m):
		### checks if the matrix is a square matrix...
		#if len(m[0]) == 2:
	if len(m) == 4:
		### a 2x2 square matrix...
			#new_matrix = [[0,0],[0,0]]
		new_matrix = [0,0,0,0]
			### place holding the matrix info...
			#placeholder2x2_1 = m[0][0]
			#placeholder2x2_2 = m[0][1]
			#placeholder2x2_3 = m[1][0]
			#placeholder2x2_4 = m[1][1]
		placeholder2x2_1 = m[0]
		placeholder2x2_2 = m[1]
		placeholder2x2_3 = m[2]
		placeholder2x2_4 = m[3]
			
		new_matrix[0] = float(placeholder2x2_4)
		if placeholder2x2_2 >= 1:
			### the number is positive...
			new_matrix[1] = -float(placeholder2x2_2)
		elif placeholder2x2_2 <= -1:
			### the number is negative...
			new_matrix[1] = float(placeholder2x2_2)
		else:
			### number is 0
			new_matrix[1] = float(placeholder2x2_2)
		if placeholder2x2_3 >= 1:
			### the number is positive...
			new_matrix[2] = -float(placeholder2x2_3)
		elif placeholder2x2_3 <= -1:
			### the number is negative...
			new_matrix[2] = float(placeholder2x2_3)
		else:
			### number is 0
			new_matrix[2] = float(placeholder2x2_3)
		new_matrix[3] = float(placeholder2x2_1)
			
		#if (new_matrix[0][0] * new_matrix[1][1]) - (new_matrix[0][1] * new_matrix[1][0]) > 0:
		if (new_matrix[0] * new_matrix[3]) - (new_matrix[1] * new_matrix[2]) > 0:
			#determinant_scalar = snapped(1.0 / snapped((new_matrix[0][0] * new_matrix[1][1]) - (new_matrix[0][1] * new_matrix[1][0]),.001),.001)
			determinant_scalar = snapped(1.0 / snapped((new_matrix[0] * new_matrix[3]) - (new_matrix[1] * new_matrix[2]),.001),.001)
				
			#Multiply_Matrix_by_Scalar(new_matrix,determinant_scalar,true)
			
			#updated_new_matrix = [[0,0],[0,0]]
			#updated_new_matrix[0][0] = snapped((new_matrix[0][0] * determinant_scalar),.01)
			#updated_new_matrix[0][1] = snapped((new_matrix[0][1] * determinant_scalar),.01)
			#updated_new_matrix[1][0] = snapped((new_matrix[1][0] * determinant_scalar),.01)
			#updated_new_matrix[1][1] = snapped((new_matrix[1][1] * determinant_scalar),.01)
				
			updated_new_matrix = [0,0,0,0]
			updated_new_matrix[0] = snapped((new_matrix[0] * determinant_scalar),.01)
			updated_new_matrix[1] = snapped((new_matrix[1] * determinant_scalar),.01)
			updated_new_matrix[2] = snapped((new_matrix[2] * determinant_scalar),.01)
			updated_new_matrix[3] = snapped((new_matrix[3] * determinant_scalar),.01)
				

		else:
			#print('The matrix does not have an Inverse.')
			return m
		#elif len(m[0]) == 3:
			### a 3x3 square matrix...
		#	pass
		#elif len(m[0]) == 4:
			### a 4x4 square matrix...
		#	pass
	else:
		#print('to find the inverse matrix the original matrix has to be square')
		#print('m by n where m = n')
		pass
	
	return updated_new_matrix

func Add_Matrix(m:Array,n:Array):
	#Add two matrices if they are the same shape(rows) and size(columns)."""
	new_matrix = []

	#if len(m) == len(n) and len(m[0]) == len(n[0]):
	if len(m) == len(n):
		### rows check...### columns check...
		#new_matrix = [[0,0],[0,0]]
		#new_matrix[0][0] = snapped(m[0][0] + n[0][0],.01)
		#new_matrix[0][1] = snapped(m[0][1] + n[0][1],.01)
		#new_matrix[1][0] = snapped(m[1][0] + n[1][0],.01)
		#new_matrix[1][1] = snapped(m[1][1] + n[1][1],.01)
		
		new_matrix = [0,0,0,0]
		new_matrix[0] = snapped(m[0] + n[0],.01)
		new_matrix[1] = snapped(m[1] + n[1],.01)
		new_matrix[2] = snapped(m[2] + n[2],.01)
		new_matrix[3] = snapped(m[3] + n[3],.01)
		
		
		return new_matrix
	
	else:
		#print('not the same shape or same size')
		pass

func Subtract_Matrix(m:Array,n:Array):
	#Add two matrices if they are the same shape(rows) and size(columns)."""
	new_matrix = []

	#if len(m) == len(n) and len(m[0]) == len(n[0]):
	if len(m) == len(n):
		#new_matrix = [[0,0],[0,0]]
		#new_matrix[0][0] = snapped(m[0][0]- n[0][0],.01)
		#new_matrix[0][1] = snapped(m[0][1] - n[0][1],.01)
		#new_matrix[1][0] = snapped(m[1][0] - n[1][0],.01)
		#new_matrix[1][1] = snapped(m[1][1] - n[1][1],.01)
		new_matrix = [0,0,0,0]
		new_matrix[0] = snapped(m[0] - n[0],.01)
		new_matrix[1] = snapped(m[1] - n[1],.01)
		new_matrix[2] = snapped(m[2] - n[2],.01)
		new_matrix[3] = snapped(m[3] - n[3],.01)
		
		
		return new_matrix

	else:
		#print('not the same shape or same size')
		pass


func Add_Matrix_by_Scalar(m:Array,n,_is_float:bool):
	#you add a matrix (ab) by a (n) number, 
	#you add every element in the matrix (ab) by the (n) number """
	
	#new_matrix = [[0,0],[0,0]]
	#new_matrix[0][0] = snapped((m[0][0] + n),.01)
	#new_matrix[0][1] = snapped((m[0][1] + n),.01)
	#new_matrix[1][0] = snapped((m[1][0] + n),.01)
	#new_matrix[1][1] = snapped((m[1][1] + n),.01)
	new_matrix =[0,0,0,0]
	new_matrix[0] = snapped((m[0] + n),.01)
	new_matrix[1] = snapped((m[1] + n),.01)
	new_matrix[2] = snapped((m[2] + n),.01)
	new_matrix[3] = snapped((m[3] + n),.01)
	
	
	return new_matrix

func Multiply_Matrix_by_Scalar(m:Array,n,_is_float:bool):
	
	#you multiply a matrix (ab) by a (n) number, 
	#you multiply every element in the matrix (ab) by the (n) number """
	#new_matrix = [[0,0],[0,0]]
	#new_matrix[0][0] = snapped((m[0][0] * n),.01)
	#new_matrix[0][1] = snapped((m[0][1] * n),.01)
	#new_matrix[1][0] = snapped((m[1][0] * n),.01)
	#new_matrix[1][1] = snapped((m[1][1] * n),.01)
	if len(m)==2:
		new_matrix =[0,0]
		new_matrix[0] = snapped((m[0] * n),.01)
		new_matrix[1] = snapped((m[1] * n),.01)
	elif len(m)== 4:
		new_matrix =[0,0,0,0]
		new_matrix[0] = snapped((m[0] * n),.01)
		new_matrix[1] = snapped((m[1] * n),.01)
		new_matrix[2] = snapped((m[2] * n),.01)
		new_matrix[3] = snapped((m[3] * n),.01)
		

	#if len(new_matrix) == 2 and len(new_matrix[0]) == 1:
	#if len(new_matrix) == 2:
	#	var new_vector = Vector2(0,0)
	#	#new_vector = Vector2(new_matrix[0][0],new_matrix[1][0])
	#	new_vector = Vector2(new_matrix[0],new_matrix[1])
	#	return new_vector
	#else:
	return new_matrix


func Divide_Matrix_by_Scalar(m:Array,n,_is_float:bool):
	#you divide a matrix (ab) by a (n) number, 
	#you divide every element in the matrix (ab) by the (n) number """
	#new_matrix = [[0,0],[0,0]]
	#new_matrix[0][0] = snapped((m[0][0] / n),.01)
	#new_matrix[0][1] = snapped((m[0][1] / n),.01)
	#new_matrix[1][0] = snapped((m[1][0] / n),.01)
	#new_matrix[1][1] = snapped((m[1][1] / n),.01)
	new_matrix =[0,0,0,0]
	new_matrix[0] = snapped((m[0] / n),.01)
	new_matrix[1] = snapped((m[1] / n),.01)
	new_matrix[2] = snapped((m[2] / n),.01)
	new_matrix[3] = snapped((m[3] / n),.01)
	
	return new_matrix

func Multiply_Matrix(m:Array,n:Array):
	#In order to multiply two matrices A and B to get AB the number of
	#columns of A must equal the number of rows of B"""
	#The Matrix AB size will be the row size of A with column size of B
	
	if len(m) == 2 and len(n) == 1:
		## 2x1* 1x1
		### results in a 2x1 matrix...
		pass
	elif len(m) == 2 and len(n) == 2:
		## 2x1 * 1x2
		### results in a 2x2 matrix...
		new_matrix = [0,0,0,0]
		new_matrix[0] = snapped( (m[0] * n[0]) ,.01)
		new_matrix[1] = snapped( (m[0] * n[1]) ,.01)
		new_matrix[2] = snapped( (m[1] * n[0]) ,.01)
		new_matrix[3] = snapped( (m[1] * n[1]) ,.01)
		
	elif len(m) == 4 and len(n) == 2:
		## 2x2 * 2x1
		### results in a 2x1 matrix...
		new_matrix = [0,0]
		new_matrix[0] = snapped( ( (m[0] * n[0]) + (m[1]*n[1]) ),.01)
		new_matrix[1] = snapped( ( (m[2] * n[0]) + (m[3]*n[1]) ),.01)
	
	elif len(m) == 4 and len(n) == 4:
		## 2x2 * 2x2
		### results in a 2x2 matrix...
		new_matrix = [0,0,0,0]
		new_matrix[0] = snapped( ( (m[0] * n[0]) + (m[1]*n[2]) ),.01)
		new_matrix[1] = snapped( ( (m[0] * n[1]) + (m[1]*n[3]) ),.01)
		new_matrix[2] = snapped( ( (m[2] * n[0]) + (m[3]*n[2]) ),.01)
		new_matrix[3] = snapped( ( (m[2] * n[1]) + (m[3]*n[3]) ),.01)
		
		
	return new_matrix
	
func Multiply_Matrix_by_Vector2_to_Vector2(m:Array,v:Vector2):
	# a matrix is multipled by a vector2 and a vector is return.

	converted_vector = [v.x,v.y]

	if len(m) == 4 and len(converted_vector) == 2:
		## 2x2 * 2x1
		### results in a 2x1 matrix...
		new_matrix = [0,0]
		new_matrix[0] = snapped( ( (m[0] * converted_vector[0]) + (m[1]*converted_vector[1]) ),.01)
		
		new_matrix[1] = snapped( ( (m[2] * converted_vector[0]) + (m[3]*converted_vector[1]) ),.01)
		
	new_vector = Vector2(new_matrix[0],new_matrix[1])
	
	return new_vector

func Multiply_Matrix_by_Vector2_to_Matrix(m:Array,v:Vector2,vector_transposed:bool,flip:bool):
	#A Matrix multiplied to a Vector , the vector is transformed into a 2x1-Matrix " 
	#Result into a New Matrix..."'''

	if vector_transposed == true:
		converted_vector = [v.x,v.y]
	else:
		### matrices of any sizes...
		converted_vector = [v.x,v.y]
		
	if flip == false:
		#Multiply_Matrix(m,converted_vector)
		
		if len(m) == 2 and len(converted_vector) == 2 and vector_transposed == true:
			## 2x1 * 1x2
			### results in a 2x2 matrix...
			new_matrix = [0,0,0,0]
			
			new_matrix[0] = snapped( (m[0] * converted_vector[0]) ,.01)
			new_matrix[1] = snapped( (m[0] * converted_vector[1]) ,.01)
			new_matrix[2] = snapped( (m[1] * converted_vector[0]) ,.01)
			new_matrix[3] = snapped( (m[1] * converted_vector[1]) ,.01)
						
		elif len(m) == 4 and len(converted_vector) == 2 and vector_transposed == false:
			## 2x2 * 2x1
			### results in a 2x1 matrix...
			new_matrix = [0,0]
			new_matrix[0] = snapped( ( (m[0] * converted_vector[0]) + (m[1]*converted_vector[1]) ),.01)
			new_matrix[1] = snapped( ( (m[2] * converted_vector[0]) + (m[3]*converted_vector[1]) ),.01)
						
	elif flip == true:
		#Multiply_Matrix(converted_vector,m)
		if vector_transposed == true:
			if len(m) == 4 and len(converted_vector) == 2:
				## 1x2  2x1
				### results in a 1x1 matrix...
				pass
			elif len(converted_vector) == 2 and len(m) == 4 :
				## 1x2 2x2 
				### results in a 1x2 matrix...
				pass
		if vector_transposed == false:
			if len(converted_vector) == 2 and len(m) == 1:
				### 2x1 1x1
				### results in a 2x1 matrix...
				new_matrix = [0,0]
				new_matrix[0] = snapped( ( (converted_vector[0] * m[0]) + (converted_vector[1]*m[1]) ),.01)
				new_matrix[1] = snapped( ( (converted_vector[2] * m[0]) + (converted_vector[3]*m[1]) ),.01)
				
			elif len(converted_vector) == 2 and len(m) == 2:
				### 2x1 1x2
				### results in a 2x2 matrix...
				new_matrix = [0,0,0,0]
			
				new_matrix[0] = snapped( (converted_vector[0] * m[0]) ,.01)
				new_matrix[1] = snapped( (converted_vector[0] * m[1]) ,.01)
				new_matrix[2] = snapped( (converted_vector[1] * m[0]) ,.01)
				new_matrix[3] = snapped( (converted_vector[1] * m[1]) ,.01) 
			
	return new_matrix

func Multiply_Vector2_by_Vector2_to_Matrix(m:Vector2,m_transposed:bool,n:Vector2,n_transposed:bool):
	# Two vectors are turned to matrices then multiplied...
	# Not possible to multiply to vectors after turning into matrix unless 1 of the matrix is transposed...
	
	var matrix_a 
	var matrix_b
	
	matrix_a = [m.x,m.y]
	matrix_b = [n.x,n.y]
	
	if m_transposed == true and n_transposed == true:
		### 2x1 2x1
		# can't multiply to matrix
		pass
	elif m_transposed == true and n_transposed == false:
		###1x2 2x1
		### results in a 1x1 matrix...
		pass
	elif m_transposed == false and n_transposed == true:
		###2x1 1x2
		### results in a 2x2 matrix...
		new_matrix = [0,0,0,0]
			
		new_matrix[0] = snapped( (matrix_a[0] * matrix_b[0]) ,.01)
		new_matrix[1] = snapped( (matrix_a[0] * matrix_b[1]) ,.01)
		new_matrix[2] = snapped( (matrix_a[1] * matrix_b[0]) ,.01)
		new_matrix[3] = snapped( (matrix_a[1] * matrix_b[1]) ,.01) 
			
	elif m_transposed == false and n_transposed == false:
		#2x1 2x1
		# can't multiply to matrix
		pass
	return new_matrix
