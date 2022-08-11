extends Node


var confirm_simulation : File #to search for the corresponding file.
var simulation : Object #the scene containing the simulation
var simulation_window : Rect2 #window outline
var simulation_outline : Dictionary #collision barriers of the window outline
var outline_data : Dictionary
### wall definition
var wall_define : float
###

func Get_Simulation_Domain():
	### the window is the outline of the simulation...
	#simulation_window = Rect2i(Vector2(0,0),
	#Vector2(ProjectSettings.get_setting('display/window/size/width'),
	#ProjectSettings.get_setting('display/window/size/height')))
	
	simulation_outline['top'] = {'coefficient of restitution' : 0.0,'coefficient of static friction': 0.0,'coefficient of kinetic friction': 0.0,'outline':0,'mass':0.0,'velocity':Vector2(0.0,0.0),'momentum':Vector2(0.0,0.0)}
	simulation_outline['right'] = {'coefficient of restitution' : 0.0,'coefficient of static friction': 0.0,'coefficient of kinetic friction': 0.0,'outline':ProjectSettings.get_setting('display/window/size/width'),'mass':0.0,'velocity':Vector2(0.0,0.0),'momentum':Vector2(0.0,0.0)}
	simulation_outline['bottom'] = {'coefficient of restitution' : 0.0,'coefficient of static friction': 0.0,'coefficient of kinetic friction': 0.0,'outline':ProjectSettings.get_setting('display/window/size/height'),'mass':0.0,'velocity':Vector2(0.0,0.0),'momentum':Vector2(0.0,0.0)}
	simulation_outline['left'] ={'coefficient of restitution' : 0.0,'coefficient of static friction': 0.0,'coefficient of kinetic friction': 0.0,'outline':0,'mass':0.0,'velocity':Vector2(0.0,0.0),'momentum':Vector2(0.0,0.0)}
	
# Called every frame. 'delta' is the elapsed time since the previous frame.
#func _process(delta):
#	pass

func _on_Test_Area_ready():
	set_owner(get_tree().get_root())
	Get_Simulation_Domain()
	
	confirm_simulation = File.new()
	if confirm_simulation.file_exists("res://Simulation.tscn"):
		simulation = load("res://Simulation.tscn").instantiate()
		
	else:
		### file does not exists..
		print('simulation does not exists...')
	
	add_child(simulation)
	
