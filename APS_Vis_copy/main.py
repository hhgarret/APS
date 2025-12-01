# An attempt to visualise the APS demo using VPython
from vpython import *
import numpy as np
import math
import time
import pandas as pd
from pynput import keyboard

# TODO: Show #s and state information in side of canvas
# TODO: Show different color for boardcasting nodes

# pyqt can't be used on windows? Look into work around
# set_browser(type='pyqt')
scene = canvas(width=1200, height=900)
scene.align = 'left'
# scene.caption = """In GlowScript programs:
# Right button drag or Ctrl-drag to rotate "camera" to view scene.
# To zoom, drag with middle button or Alt/Option depressed, or use scroll wheel.
#   On a two-button mouse, middle is left + right.
# Touch screen: pinch/extend to zoom, swipe or two-finger rotate."""
scene.background = color.white
scene.forward = vector(0,-.3,-1)
scene.autoscale = False
# scene.userzoom = False

def square_matrix(arr):
    retM = np.zeros((len(arr), len(arr)))
    for index in range(len(arr)):
        retM[index][index] = np.real(arr[index])
    return retM

def mds(dists):
    numPoints = len(dists)
    C = np.eye(len(dists)) - (np.ones((len(dists), len(dists)))/len(dists))  # centering matrix C
    Dsquared = np.zeros((len(dists), len(dists))) # piece-wise squared distance matrix
    for i in range(len(dists)):
        for j in range(len(dists)):
            Dsquared[i][j] = (dists[i][j]*dists[i][j])
    B = -0.5*np.linalg.multi_dot((C,Dsquared,C))  # B = -0.5*C*D^(2)*C

    ret = np.linalg.eig(B)
    # Get m largest eigenvalues
    theoreticalEigenValues = np.argsort(ret.eigenvalues)[-2:]  # Get largest eigenvalues/vectors
    # print(theoreticalEigenValues)
    eigenValues = ret.eigenvalues[theoreticalEigenValues]
    # And corresponding eigenvectors
    eigenVectors = ret.eigenvectors[:, theoreticalEigenValues]
    eigenVectors = eigenVectors / np.linalg.norm(eigenVectors, axis=0, ord=2)
    d2c = (np.dot(eigenVectors, np.sqrt(square_matrix(eigenValues))))  # Recalculate new points
    # d2c = np.asmatrix(d2c)

    newdists = np.zeros((numPoints, numPoints))
    for i in range(numPoints):
        for j in range(numPoints):
            newdists[i][j] = (np.linalg.norm(d2c[i] - d2c[j]))

    rmse = 0  # Calculate difference between initial data and answer to geometry
    for i in range(numPoints):
        for j in range(numPoints):
            rmse += ((dists[i][j]-newdists[i][j])**2)/(numPoints**2)

    rmse = math.sqrt(rmse)
    d2c = np.real(d2c)
    # print(d2c)
    return d2c, rmse

class Tree:
    def __init__(self, pos):
        self.pos = pos
        self.stump = cylinder(pos = self.pos, size = vector(10, 2.5, 2.5), axis=vector(0,1,0), color= vector(150/255, 75/255, 0))
        self.leaves = []
        self.numleaves = 2
        for i in range(self.numleaves):
            self.leaves.append(cone(pos = self.pos+vector(0, (1+(i/self.numleaves))*5, 0), length=5, radius=0.75*5, axis=vector(0,1,0), color=color.green))
        self.tree = compound([self.stump] + self.leaves)

class Node:
    def __init__(self, pos, scalefactor = 1, alpha=1, opacity=1):
        self.pos = pos
        self.scalefactor = scalefactor
        self.base1 = cylinder(pos = self.pos, size = vector(0.04, 0.13, 0.13)*self.scalefactor, axis = vector(0,1,0), color=color.red, alpha=alpha, opacity=opacity)
        self.base2 = cone(pos=self.pos+vector(0, 0.04*self.scalefactor, 0), length=0.014*self.scalefactor, radius=0.065*self.scalefactor, axis = vector(0,1,0), color=color.black, alpha=alpha, opacity=opacity)
        self.base3 = cylinder(pos=self.pos+vector(0, 0.04*self.scalefactor, 0), size=vector(.014, 0.114, 0.114)*self.scalefactor, axis=vector(0,1,0), color=color.white, alpha=alpha, opacity=opacity)
    def diff(self, node2):
        pos1 = np.asarray((self.pos.x, self.pos.y, self.pos.z))
        pos2 = np.asarray((node2.pos.x, node2.pos.y, node2.pos.z))
        return np.linalg.norm(pos1 - pos2)
    def update_pos(self, newpos):
        self.base1.pos = newpos
        self.base2.pos = newpos+vector(0, 0.04*self.scalefactor, 0)
        self.base3.pos = newpos+vector(0, 0.04*self.scalefactor, 0)
    def set_emissive(self, emissive_bool):
        self.base1.emissive = emissive_bool
        self.base2.emissive = emissive_bool
        self.base3.emissive = emissive_bool



islandwidth = 100
phi = (1 + np.sqrt(5)) / 2


# Add in an environment and populate it with trees
ground = cylinder(pos=vector(0, -1, 0), size = vector(1, islandwidth, islandwidth), axis=vector(0, 1,0), color = color.yellow)
dirt = cylinder(pos=vector(0, -1.1, 0), size = vector(1, islandwidth+1, islandwidth+1), axis=vector(0, 1, 0), color=vector(150/255, 75/255, 0))
trees = []
numtrees = 15
for i in range(numtrees):
    trees.append(Tree(vector(.45*islandwidth*np.cos(i*phi*2*np.pi/numtrees),0, .45*islandwidth*np.sin(i*phi*2*np.pi/numtrees))))

# Instantiate each node on the grass
nodes = []
nodesize = 5

original_points = pd.read_csv("Sample 25 Node Data/1_coordinates_original.txt",sep="\t", header=None).values + np.asarray([-6, -6])
# nxn square shape
for i in range(25):
        # nodes.append(Node(vector(3*i - 6, 0, 3*j - 6), 5))
        nodes.append(Node(vector(original_points[i][0], 0, original_points[i][1]), nodesize))


class Skybox:
    def __init__(self):
        # Sky boxes
        hshape1 = shapes.arc(angle1=0, angle2=pi/2, radius=islandwidth*1.5, thickness=0.01)
        hpath1 = paths.circle(radius=0.5)
        demihemi1 = extrusion( shape=hshape1, path=hpath1, color=color.cyan)

        hshape2 = shapes.arc(angle1=0, angle2=-1*pi/2, radius=islandwidth*1.5, thickness=0.01)
        hpath2 = paths.circle(radius=0.5)
        demihemi2 = extrusion( shape=hshape2, path=hpath2, color=color.black)
Skybox()

rad = sphere(pos = nodes[0].pos, opacity=0.2, radius = 0, color=color.red)


# Logic for animating the demonstration
print(nodes[0].diff(nodes[1]))
dists = [ [0.0 for i in range(len(nodes))] for i in range(len(nodes))]
def mat_to_str(mat):
    retstr = ""
    for row in mat:
        for column in row:
            retstr += f"{column:05.2f}" + ", "
        retstr += "\n"
    return retstr
# scene2 = canvas(width=1200, height=750)
# # scene2.forward = vector(0,-.3,-1)
# scene2.autoscale = False
# scene2.userzoom = False
# scene2.userspin = False
# scene2.userpan = False
#
# l = label(text=mat_to_str(dists), scene=scene2, yoffset=1, pos=vector(0,-7,0))

visibility = 3*sqrt(2)

mds_coords = pd.read_csv("Sample 25 Node Data/2_coordinates_MDS.txt",sep="\t", header=None).values
print(mds_coords, len(mds_coords))
opt_coords = pd.read_csv("Sample 25 Node Data/3_coordinates_optimization.txt", sep="\t", header=None).values
opt_coords_rot = pd.read_csv("Sample 25 Node Data/4_coordinates_optimization_rotated.txt", sep="\t", header=None).values + np.asarray([-6, -6])

lerpmax = 200
lerpcurrent = 0
lerpstep = 0
def lerp_function(p, x, y):
    output = []
    for i in range(len(x)):
        output.append(x[i] * (1-p) +y[i]*p)
    return output
# exit()


running = True

def toggle_pause(b):
    global running
    running = not running
    if running:
        b.text = "Pause"
    else:
        b.text = "Play"
# Create a listener
b = button(bind=toggle_pause, text="Pause")

radpos = 0
nodesmade = False
newnodes = []
time.sleep(5)
# text(text="testing, testing, 123")
while True:
    rate(60)
    if not running:
        continue
    if rad.radius < 10:
        if rad.radius == 0:
            for i, node in enumerate(nodes):
                if nodes[radpos].diff(node) <= visibility:
                    node.set_emissive(True)
                else:
                    node.set_emissive(False)
        rad.radius += (2/6)
    elif radpos+1 < len(nodes):
        rad.radius = 0
        for i, node in enumerate(nodes):
            dists[radpos][i] = nodes[radpos].diff(node)
        # print(dists[radpos])
        # l.text = mat_to_str(dists)
        radpos += 1
        rad.pos = nodes[radpos].pos
    elif rad.visible:
        rad.visible = False
        time.sleep(3)
        for i, node in enumerate(nodes):
            dists[radpos][i] = nodes[radpos].diff(node)
            node.set_emissive(False)
        # print(dists[radpos])
        # l.text = mat_to_str(dists)
        # d2c, rmse = mds(dists)
        # print("RMSE: ", rmse)
        d2c = mds_coords
        for coord in d2c:
            newnodes.append(Node(vector(coord[0], 0, coord[1]), scalefactor=nodesize, opacity=0.5))
        nodesmade = True
        # time.sleep(3)
        toggle_pause(b)
    elif nodesmade and lerpstep == 0:
        # if nodes[0].diff(newnodes[0]) > .1:
        #     for newnode in newnodes:
        #         possave = newnode.pos
        #         rotationr = -0.25 * np.pi / 180;
        #
        #         newnode.pos = vector((newnode.pos.x * cos(rotationr) - newnode.pos.z * sin(rotationr)), newnode.pos.y, (newnode.pos.z * cos(rotationr) + newnode.pos.x * sin(rotationr)))
        #         newnode.update_pos(newnode.pos)
        if lerpcurrent < lerpmax:
            for i, newnode in enumerate(newnodes):
                # possave = newnode.pos
                newpos = lerp_function(lerpcurrent/lerpmax, mds_coords[i], opt_coords[i])
                newnode.pos = vector(newpos[0], newnode.pos.y, newpos[1])
                newnode.update_pos(newnode.pos)
            lerpcurrent += 1
        else:
            lerpstep += 1
            lerpcurrent = 0
            toggle_pause(b)
    elif nodesmade and lerpstep == 1:
        if lerpcurrent < lerpmax:
            for i, newnode in enumerate(newnodes):
                # possave = newnode.pos
                newpos = lerp_function(lerpcurrent/lerpmax, opt_coords[i], opt_coords_rot[i])
                newnode.pos = vector(newpos[0], newnode.pos.y, newpos[1])
                newnode.update_pos(newnode.pos)
            lerpcurrent += 1
        else:
            lerpstep += 1

if __name__ == "__main__":
    print("running")
