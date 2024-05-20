import numpy as np
import plotly.graph_objects as go
import turtle
import warnings

while True:
    try:
        sides = np.zeros(3)
        sides[0] = float(input('Side no 1: '))
        sides[1] = float(input('Side no 2: '))
        sides[2] = float(input('Side no 3: '))
        rotation = float(input('Starting angle from horizontal: '))
        sides = sorted(sides)
        break
    except ValueError:
        print('Sides cannot be strings')

a,b,c = sides[0], sides[1], sides[2]
# If triangle sides are wrong then it will throw up a run time error
# when calculating the angles

with warnings.catch_warnings(record=True) as w:
    
    alpha = np.arccos((b**2+c**2-a**2)/(2*b*c))*(180/np.pi)
    beta = np.arccos((a**2+c**2-b**2)/(2*a*c))*(180/np.pi)
    gamma = np.arccos((a**2+b**2-c**2)/(2*a*b))*(180/np.pi)
   
    if len(w) > 0:
        print('These side lengths cannot form a triangle')
        exit(0)
    else:
        area = 0.5*a*b*np.sin(gamma*np.pi/180)
        height = 2*area/a
        

# Rearrange angles for turtle drawing
# each angle is opposite its corresponding side so needs to be shifted
angles = [gamma, alpha, beta]

# Define the Triangle Type
if a == b and b == c:
    tritype = 'equilateral triangle'
elif (a == b and a != c) or (b == c and b != a) or (c == a and c != b):
    tritype = 'isoceles triangle'
elif a != b and a != c:
    tritype = 'scalene triangle'

# Last check involves angles so should be seperate so it isn't skipped
if np.isclose(90, gamma, atol = 0.001):
    tritype = 'Right Angled triangle'

t = turtle.Turtle()
screen = turtle.Screen()
rootwindow = screen.getcanvas().winfo_toplevel()
rootwindow.call('wm', 'attributes', '.', '-topmost', '1')
rootwindow.call('wm', 'attributes', '.', '-topmost', '0')
screen.setup(800,800)
# Adjust starting position
t.penup()
t.goto(-10*a, -10*height)
t.pendown()
t.left(rotation)
for i in range(3):
    #distance = sides[i]
    t.forward(distance=50*sides[i])
    t.left(180-angles[i])
    t.write(sides[i], align= 'right', font = ('Verdana', 10, 'normal'))

# Print out area and triangle type
t.hideturtle()
t.penup()
t.goto((-10*a)-50,-10*height-50)
t.write(f'Area = {np.round(area,3)} units squared', 
        font = ('Verdana', 10, 'normal'))
t.penup()
t.goto((-10*a)-50,-10*height-65)
t.write(f'Triangle is a {tritype}', 
        font = ('Verdana', 10, 'normal'))

turtle.done()