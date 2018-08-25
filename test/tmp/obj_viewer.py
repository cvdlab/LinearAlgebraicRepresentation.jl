from pyplasm import *
batches=[]
batches+=Batch.openObj("two_3x3cubes.obj")
octree=Octree(batches)
glcanvas=GLCanvas()
glcanvas.setOctree(octree)
glcanvas.runLoop()

