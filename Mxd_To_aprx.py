import arcpy

import os

arcpy.env.overwriteOutput = True

def Mxd_format_to_10_3(MXD_to_convert,temp_aprx,temp_gdb):
    folder      = os.path.dirname(MXD_to_convert[0])
    name        = os.path.basename(MXD_to_convert[0]).split('.')[0] + '.aprx'
    New_aprx     = folder +'\\'+name

    aprx = arcpy.mp.ArcGISProject(temp_aprx)
    aprx.importDocument(MXD_to_convert[0])

    m = aprx.listMaps('Layers')[0]
    m.updateConnectionProperties(temp_gdb,MXD_to_convert[1])
    aprx.saveACopy(New_aprx)


def Fix_fodler_raster(Fodler):
    save_mxd,save_gdb,= [],[]

    for root,dirs,files in os.walk(Fodler):
        for file in files:
            if file.endswith('.mxd') and 'Copy' not in file:
                print (root + '\\' + file)
                save_mxd.append(root + '\\' + file)
            if file.endswith('gdb'):
                gdb_name = os.path.dirname(root + '\\' + file)
                save_gdb.append(gdb_name)
                print (gdb_name)


    return list(zip(save_mxd,save_gdb))


folder    = r'C:\Users\Administrator\Desktop\medad\python\Work\Mpy_Pro\Source_files'
temp_aprx = r"C:\Users\Administrator\Desktop\medad\python\Work\Mpy_Pro\aprx_temp\EditTazar89689\Edit89689.aprx"
temp_gdb  = r"C:\Users\Administrator\Desktop\medad\python\Work\Mpy_Pro\aprx_temp\EditTazar89689\CadasterEdit_Tazar.gdb"

save_mxd  = Fix_fodler_raster(folder)

for i in save_mxd:
    Mxd_format_to_10_3(i,temp_aprx,temp_gdb)

