# Material-Point-Method-in-Godot

Looking to implementing MLS version of Material Point Method in Godot Engine.

up to Version.10.6.7:

A project I wanted to work on and for that project, I ended up wanting/needed to use Material Point Method. 
This is the MLS-version of Material Point Method I'm looking to implement (so far thru nodes}. 
I have been using Godot 4 , v4.2.stable.official [46dc27791] being the latest

some progress clips.

video 1
![1 6 0  random scattering -2024-03-11 15-22-33](https://github.com/Exis10tial/Material-Point-Method-in-Godot/assets/62639345/1258e0d3-94e6-4b7c-a7f9-009a39401e16)


video 2
![1 6 0   drop 2024-03-11 15-30-13](https://github.com/Exis10tial/Material-Point-Method-in-Godot/assets/62639345/a0951725-9d79-4f47-8fc9-9fa777ea1259)


[Video Notes]
Window Size is 1152 by 648

Now the particles is represented as rectangles(Rect2)-size of 1.

No models is being used.

Both videos contains particle count of 100.

video 1- 100 particles is at 20-23 fps.

video 2- 100 particles is at 18-20 fps.

1st video there is no gravity or velocity affecting these particles. The particles just randomly jitter in place.

2nd video does have gravity (0.0,9.8), but still no velocity on the particles. The particles drop cause of gravity then bounce cause of walls.
