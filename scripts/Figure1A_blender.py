from blender_util import *

def run(csv_fn, pial_fn, res_seg_fn, scale, x0, y0, z0, lesion_fn=None):
    scene1=addScene('scene1', 'NEW')

    bpy.context.screen.scene=bpy.data.scenes[scene1.name]

    # lamp1 = addLamp('lamp1', (0,0,200*scale), 'SPOT')
    # lamp1.data.energy = 5.

    cam1 = addCamera('cam1', (scale,-scale,0))

    bpy.ops.import_scene.obj(filepath=res_seg_fn)
    res_seg_obj = bpy.data.objects[res_seg_fn.split('/')[-1].split('.obj')[0]]
    res_seg_obj.scale = (scale,scale,scale)
    # Deselect All
    bpy.ops.object.select_all(action='DESELECT')
    # Select the new object.
    res_seg_obj.select = True
    bpy.context.scene.objects.active = res_seg_obj


    bpy.ops.import_scene.obj(filepath=pial_fn)
    pial_obj = bpy.data.objects[pial_fn.split('/')[-1].split('.obj')[0]]
    pial_obj.scale = (scale,scale,scale)

    if('lh' in pial_fn):
        bpy.ops.import_scene.obj(filepath=pial_fn.replace('lh','rh'))
        pial_obj2 = bpy.data.objects[pial_fn.replace('lh','rh').split('/')[-1].split('.obj')[0]]
    else:
        bpy.ops.import_scene.obj(filepath=pial_fn.replace('rh','lh'))
        pial_obj2 = bpy.data.objects[pial_fn.replace('rh','lh').split('/')[-1].split('.obj')[0]]
    pial_obj2.scale = (scale,scale,scale)

    # Deselect All
    bpy.ops.object.select_all(action='DESELECT')
    # Select the new object.
    pial_obj.select = True
    bpy.context.scene.objects.active = pial_obj

    if lesion_fn is not None:
        bpy.ops.import_scene.obj(filepath=lesion_fn)
        lesion_seg_obj = bpy.data.objects[lesion_fn.split('/')[-1].split('.obj')[0]]
        lesion_seg_obj.scale = (scale,scale,scale)

    # mod = pial_obj.modifiers.new("Boolean", type='BOOLEAN')
    # mod.operation = 'INTERSECT'
    # mod.object = res_seg_obj

    # # large cube has context.
    # bpy.ops.object.modifier_apply(apply_as='DATA', modifier=mod.name)

    # bpy.context.scene.objects.unlink(res_seg_obj)
    # bpy.data.objects.remove(res_seg_obj)

    # bpy.ops.import_scene.obj(filepath=res_seg_fn)
    # res_seg_obj = bpy.data.objects[res_seg_fn.split('/')[-1].split('.obj')[0]+'.001']
    # res_seg_obj.scale = (scale,scale,scale)

    addTrackToConstraint(cam1,"myTrackTo",pial_obj)
    # addTrackToConstraint(lamp1,"myTrackTo1",pial_obj)

    RR = 3
    TT = 9600

    lines = open(csv_fn,'r').readlines()
    for line in lines:
        if('rh' in pial_fn):
            origo = [scale*(float(line.split(',')[0])-x0+3),scale*(float(line.split(',')[1])-y0),scale*(float(line.split(',')[2])-z0)]
        else:
            origo = [scale*(float(line.split(',')[0])-x0-3),scale*(float(line.split(',')[1])-y0),scale*(float(line.split(',')[2])-z0)]
        label = line.split(',')[3].replace('\n','')
        value = float(line.split(',')[4].replace('\n',''))
        chCol = makeMaterial('chCol',(1,1,1), (1,1,1), 1)
        origin = Vector(origo)

        cone3 = createMeshFromPrimitive(label.replace('\n',''), origin+Vector((0,0,0)), scale*1.25)
        cone3.active_material = chCol
        if(value < 0):
            chCol.diffuse_color = (0,0,0)
        elif(value == 0):
            chCol.diffuse_color = (1,1,1)
        else:
            chCol.diffuse_color = (1,0,0)

    if('rh' in pial_fn):
        x = RR * math.cos(1.0/TT * 2 * math.pi)
    else:
        x = -RR * math.cos(1.0/TT * 2 * math.pi)
    y = RR * math.sin(1.0/TT * 2 * math.pi)
    z = 1
    cam1.location = (x,y,z)

    return


if __name__ == "__main__":
    comp_dir = '/Users/lkini/Documents/LittLab/aim1/results'
    patient_id = 'HUP075'
    soz_coord_fn = ''
    scale = 0.01
    x0 = 8.94377
    y0 = -11.531066
    z0 = 3.9223785
    lesion_fn = None
    for fn in os.listdir('%s/%s/aim1'%(comp_dir,patient_id)):
        if '_soz_coord.csv' in fn:
            soz_coord_fn = fn
        if 'LESION' in fn and '.obj' in fn:
            lesion_fn = fn
    if lesion_fn is not None:
        run('%s/%s/aim1/%s'%(comp_dir,patient_id,soz_coord_fn),'%s/%s/aim1/lh.pial.native.obj'%(comp_dir,patient_id),'%s/%s/aim1/%s_T1_19990604_RES_SEG_BINARY_DILATED.obj'%(comp_dir,patient_id,patient_id), scale, x0, y0, z0, lesion_fn='%s/%s/aim1/%s'%(comp_dir,patient_id,lesion_fn))
    else:
        run('%s/%s/aim1/%s'%(comp_dir,patient_id,soz_coord_fn),'%s/%s/aim1/lh.pial.native.obj'%(comp_dir,patient_id),'%s/%s/aim1/%s_T1_19990604_RES_SEG_BINARY_DILATED.obj'%(comp_dir,patient_id,patient_id), scale, x0, y0, z0, lesion_fn=None)