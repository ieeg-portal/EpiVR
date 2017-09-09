import os
import bpy
import math
import mathutils
from mathutils import Vector, Matrix

################################################## #####
# scene stype:
# NEW New, Add new scene.
# EMPTY Copy Settings, Make a copy without any objects.
# LINK_OBJECTS Link Objects, Link to the objects from the current scene.
# LINK_OBJECT_DATA Link Object Data, Copy objects linked to data from the current scene.
# FULL_COPY Full Copy, Make a full copy of the current scene.

def addScene(name, stype):
    bpy.ops.scene.new(type=stype)
    ob = bpy.context.scene
    ob.name = name
    return ob


def makeMaterial(name, diffuse, specular, alpha):
    mat = bpy.data.materials.new(name)
    mat.diffuse_color = diffuse
    mat.diffuse_intensity = 1.0
    mat.specular_color = specular
    mat.specular_intensity = 0.5
    mat.alpha = alpha
    mat.ambient = 1
    return mat

def createMeshFromPrimitive(name, origin, scale):
    bpy.ops.mesh.primitive_uv_sphere_add(location=origin, segments=64, ring_count=64)

    ob = bpy.context.object
    ob.scale = (scale,scale,scale)
    bpy.ops.object.transform_apply(scale=True)
    ob.name = name
    ob.show_name = True
    bpy.ops.object.shade_smooth()
    me = ob.data
    me.name = name+'Mesh'
    return ob

def setMaterial(ob, mat):
    me = ob.data
    me.materials.append(mat)
    return


def addTrackToConstraint(ob, name, target):
    cns = ob.constraints.new('TRACK_TO')
    print(dir(cns))
    cns.name = name
    cns.target = target
    cns.track_axis = 'TRACK_NEGATIVE_Z'
    cns.up_axis = 'UP_Y'
    cns.owner_space = 'LOCAL'
    cns.target_space = 'LOCAL'
    return


def addCamera(name, place):
    bpy.ops.object.camera_add(location=place)
    ob = bpy.context.object
    ob.name = name
    return ob

################################################## #####
# Lamp ltype:
# POINT Point, Omnidirectional point light source.
# SUN Sun, Constant direction parallel ray light source.
# SPOT Spot, Directional cone light source.
# HEMI Hemi, 180 degree constant light source.
# AREA Area, Directional area light source.

def addLamp(name, position, ltype):
    bpy.ops.object.lamp_add(location=position, type = ltype)
    ob = bpy.context.object
    ob.name = name
    return ob

