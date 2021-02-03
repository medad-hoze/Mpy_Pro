import arcpy

import os

arcpy.env.overwriteOutput = True

def Mxd_format_to_10_3(MXD_to_convert,temp_aprx):
    folder      = os.path.dirname(MXD_to_convert)
    name        = os.path.basename(MXD_to_convert).split('.')[0] + '.aprx'
    New_aprx     = folder +'\\'+name

    aprx = arcpy.mp.ArcGISProject(temp_aprx)
    aprx.importDocument(MXD_to_convert)
    aprx.saveACopy(New_aprx)

def Fix_fodler_raster(Fodler,search_endwith):
    save_mxd = []
    already_exists = []
    for root,dirs,files in os.walk(Fodler):
        for file in files:
            if file.endswith(search_endwith):
                print (root + '\\' + file)
                save_mxd.append(root + '\\' + file)
    return save_mxd


save_mxd  = Fix_fodler_raster(r'C:\Users\Administrator\Desktop\medad\python\Work\Mpy_Pro\Source_files','.mxd')
temp_aprx = r"C:\Users\Administrator\Desktop\medad\python\Work\Mpy_Pro\aprx_temp\Edit89689.aprx"

for i in save_mxd:
    Mxd_format_to_10_3(i)


