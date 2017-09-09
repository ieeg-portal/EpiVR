from blender_util import *

def get_rgb(val, minval=-1, maxval=1):
    if val <= minval:
        val = -0.99
    elif val >= maxval:
        val = 0.99
    normval = (val + 1.0)/2.0*1000.0
    cmap_lines = open('/Users/lkini/Documents/LittLab/aim3/code/scripts/coolwarm.csv','r').readlines()
    return map(float,cmap_lines[int(normval)].split(','))



def run(csv_fn, pial_fn, res_seg_fn, scale, x0, y0, z0, radius, cres_csv_fn, lesion_fn=None):

    scene1=addScene('scene1', 'NEW')

    bpy.context.screen.scene=bpy.data.scenes[scene1.name]
    if not scene1.render.engine == 'CYCLES':
        scene1.render.engine = 'CYCLES'

    bpy.ops.world.new()
    scene1.world = bpy.data.worlds['World.001']
    scene1.world.node_tree.nodes['Background'].inputs['Color'].default_value = (1.0,1.0,1.0,1.0)

    cam1 = addCamera('cam1', (scale,-scale,0))

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

    bpy.ops.import_scene.obj(filepath=res_seg_fn)
    res_seg_obj = bpy.data.objects[res_seg_fn.split('/')[-1].split('.obj')[0]]
    res_seg_obj.scale = (scale,scale,scale)
    # Deselect All
    bpy.ops.object.select_all(action='DESELECT')
    # Select the new object.
    res_seg_obj.select = True
    bpy.context.scene.objects.active = res_seg_obj
    # Add a modifier
    bpy.ops.object.modifier_add(type='SHRINKWRAP')
    mod = res_seg_obj.modifiers
    mod[0].name = "Shrinkwrap_sphere_brain"
    if ('rh' in pial_fn):
        mod[0].target = bpy.data.objects['rh.pial.native']
    else:
        mod[0].target = bpy.data.objects['lh.pial.native']
    mod[0].offset = 1
    mod[0].use_keep_above_surface = True
    # Apply modifier
    bpy.ops.object.modifier_apply(apply_as='DATA', modifier=mod[0].name)
    res_seg_material = makeMaterial('res_seg_material',(1,0,0),(1,1,1),0.1)
    # get the material
    res_seg_material = bpy.data.materials['res_seg_material']
    res_seg_material.use_nodes = True
    nodes = res_seg_material.node_tree.nodes
    for node in nodes:
        nodes.remove(node)
    # create emission node
    node_transparent = nodes.new(type='ShaderNodeBsdfTransparent')
    node_transparent.inputs[0].default_value = (0.5,0.5,0.5,0.1)  # RGBA
    # create output node
    node_output = nodes.new(type='ShaderNodeOutputMaterial')   
    node_output.location = 400,0
    # link nodes
    links = res_seg_material.node_tree.links
    link = links.new(node_transparent.outputs[0], node_output.inputs[0])
    res_seg_obj.active_material = res_seg_material


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


        ###### Make spheres
        cres = makeMaterial('cres_%s'%label,(1,1,1),(1,1,1),0.1)
        # get the material
        cres = bpy.data.materials['cres_%s'%label]
        cres.use_nodes = True
        nodes = cres.node_tree.nodes
        for node in nodes:
            nodes.remove(node)
        # create emission node
        node_transparent = nodes.new(type='ShaderNodeBsdfTransparent')
        node_transparent.inputs[0].default_value = (0.5,0.5,0.5,0.5)  # RGBA
        # create output node
        node_output = nodes.new(type='ShaderNodeOutputMaterial')   
        node_output.location = 400,0
        # link nodes
        links = cres.node_tree.links
        link = links.new(node_transparent.outputs[0], node_output.inputs[0])


        sphere = createMeshFromPrimitive('Sphere_%s'%label, origin+Vector((0,0,0)), scale*radius)
        sphere.active_material = cres

        # Deselect All
        bpy.ops.object.select_all(action='DESELECT')

        # Select the new object.
        target = sphere
        target.select = True
        bpy.context.scene.objects.active = target

        # Add a modifier
        bpy.ops.object.modifier_add(type='SHRINKWRAP')
        mod = target.modifiers
        mod[0].name = "Shrinkwrap_sphere_brain"
        if ('rh' in pial_fn):
            mod[0].target = bpy.data.objects['rh.pial.native']
        else:
            mod[0].target = bpy.data.objects['lh.pial.native']
        mod[0].offset = 1
        mod[0].use_keep_above_surface = True
        # Apply modifier
        bpy.ops.object.modifier_apply(apply_as='DATA', modifier=mod[0].name)

        timeseries_lines = open(cres_csv_fn,'r').readlines()
        label_found = False
        for timeseries_line in timeseries_lines:
            electrode_name = timeseries_line.split(',')[3].replace('"','')
            if electrode_name == label:
                label_found = True
                break
        if label_found:
            print(electrode_name, label)
            T = len(timeseries_line.split(',')[4:])
            print(T)
            # Create color animation frame by frame
            for frame in range(1,T+1):
                timeseries_val = float(timeseries_line.split(',')[4:][frame-1].replace('"',''))
                if timeseries_val < 0:
                    r = 0.0
                    g = 0.0
                    b = -timeseries_val/0.1
                else:
                    r = timeseries_val/0.1
                    g = 0.0
                    b = 0.0
                # r,g,b = get_rgb(timeseries_val,minval=-0.01,maxval=0.01)
                # print(timeseries_val,r,g,b)
                cres.node_tree.nodes['Transparent BSDF'].inputs['Color'].default_value = (r,g,b,1.0)
                cres.node_tree.nodes['Transparent BSDF'].inputs['Color'].keyframe_insert(data_path="default_value", frame=frame)
                # print(timeseries_val,cres.node_tree.nodes['Transparent BSDF'].inputs['Color'].default_value[0])


    if('rh' in pial_fn):
        x = RR * math.cos(1.0/TT * 2 * math.pi)
    else:
        x = -RR * math.cos(1.0/TT * 2 * math.pi)
    y = RR * math.sin(1.0/TT * 2 * math.pi)
    z = 1
    cam1.location = (x,y,z)
    scene1.frame_start = 1
    scene1.frame_end = 1346
    return


if __name__ == "__main__":
    comp_dir = '/Users/lkini/Documents/LittLab/aim1/results'
    patient_id = 'HUP075'
    soz_coord_fn = ''
    cres_csv_fn = 'HUP075_5_node_coordinates_cres.csv'
    scale = 0.01
    x0 = 8.94377
    y0 = -11.531066
    z0 = 3.9223785
    radius = 5
    lesion_fn = None
    for fn in os.listdir('%s/%s/aim1'%(comp_dir,patient_id)):
        if '_soz_coord.csv' in fn:
            soz_coord_fn = fn
        if 'LESION' in fn and '.obj' in fn:
            lesion_fn = fn
    if lesion_fn is not None:
        run('%s/%s/aim1/%s'%(comp_dir,patient_id,soz_coord_fn),'%s/%s/aim1/lh.pial.native.obj'%(comp_dir,patient_id),'%s/%s/aim1/%s_T1_19990604_RES_SEG_BINARY_DILATED.obj'%(comp_dir,patient_id,patient_id), scale, x0, y0, z0, radius, '%s/%s/aim3/%s'%(comp_dir, patient_id, cres_csv_fn), lesion_fn='%s/%s/aim1/%s'%(comp_dir,patient_id,lesion_fn))
    else:
        run('%s/%s/aim1/%s'%(comp_dir,patient_id,soz_coord_fn),'%s/%s/aim1/lh.pial.native.obj'%(comp_dir,patient_id),'%s/%s/aim1/%s_T1_19990604_RES_SEG_BINARY_DILATED.obj'%(comp_dir,patient_id,patient_id), scale, x0, y0, z0, radius, '%s/%s/aim3/%s'%(comp_dir, patient_id, cres_csv_fn), lesion_fn=None)

