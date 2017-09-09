from blender_util import *

def hub_run(csv_fn, pial_fn, res_seg_fn, scale, x0, y0, z0):
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
    # Deselect All
    bpy.ops.object.select_all(action='DESELECT')
    # Select the new object.
    pial_obj.select = True
    bpy.context.scene.objects.active = pial_obj

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
        
        chCol.diffuse_color = (value,0,1.0-value)

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
    patient_id = 'HUP074'
    hub_coord_fn = ''
    scale = 0.01
    x0 = 3.7572
    y0 = -7.1747
    z0 = -2.9175
    lesion_fn = None
    for fn in os.listdir('%s/%s/aim3'%(comp_dir,patient_id)):
        if '_hub_coord.csv' in fn:
            hub_coord_fn = fn
        if 'LESION' in fn and '.obj' in fn:
            lesion_fn = fn
    clip = 'B'
    pre_post = 'post'
    fconn = 'broadband_CC'
    hub_coord_fn = 'HUP074_T1_19990327_Ictal_%s_%s_%s_hub_coord.csv'%(pre_post,clip,fconn)
    print(hub_coord_fn)

    hub_run('%s/%s/aim3/%s'%(comp_dir,patient_id,hub_coord_fn),'%s/%s/aim1/lh.pial.native.obj'%(comp_dir
    ,patient_id),'%s/%s/aim1/RES_SEG_BINARY_DILATED.obj'%(comp_dir,patient_id), scale, x0, y0, z0)
