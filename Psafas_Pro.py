# -*- coding: utf-8 -*-



import json,sqlite3
import arcpy
import os,math,sys
import pandas as pd
import datetime,ast
import numpy as np


arcpy.env.overwriteOutput = True

# -*- coding: utf-8 -*-

# -*- coding: utf-8 -*-

# # # # # # # # # # # # # # # P S A F A S    T O O L S # # # # # # # # # # # # # # # 



def NewGushim(parcel_tazar, parcel_all_bankal,layer_f):

    print_arcpy_message("START Func: NewGushim",1)

    layer_finish = layer_f+'Temp'
    
    arcpy.Select_analysis                  (layer_f,layer_finish)
    arcpy.MakeFeatureLayer_management      (parcel_all_bankal, "in_memory\\bankal")
    arcpy.SelectLayerByLocation_management ("in_memory\\bankal", 'INTERSECT', parcel_tazar,"2000 Meters")
        
    arcpy.Dissolve_management              ("in_memory\\bankal", "in_memory\\sub_gush_all","GUSH_NUM", "", "SINGLE_PART")
    arcpy.Dissolve_management              (layer_finish, "in_memory\\layer_finish_gush","GUSH_NUM", "", "SINGLE_PART")     
        
    Gush_stay   = []
    Gush_change = []
        
    parcel_GUSH_bankal = [[i.GUSH_NUM,i.SHAPE] for i in arcpy.SearchCursor("in_memory\\sub_gush_all")]
    with arcpy.da.SearchCursor("in_memory\\layer_finish_gush",["GUSH_NUM","SHAPE@"]) as cursor:
                for row in cursor:
                        if parcel_GUSH_bankal:
                            if str(row[0]) == str(parcel_GUSH_bankal[0][0]):
                                    if row[1].union(parcel_GUSH_bankal[0][1]) == row[1]:
                                            Gush_stay.append(row[0])
                                    else:
                                            Gush_change.append(row[0])
        
    if len(Gush_change) > 0:
                print_arcpy_message ("there is gush that moved, name: {}".format(Gush_change),status = 1)   

        
    gushim_tazar = list(set([row.GUSH_NUM for row in arcpy.SearchCursor(parcel_tazar)]))
    gushim_all = list(set([row.GUSH_NUM for row in arcpy.SearchCursor("in_memory\\sub_gush_all")]))
    new_gushim = [g for g in gushim_tazar if g not in gushim_all]
    if len(new_gushim) > 0:
                print_arcpy_message("there is  New Gush in the tazar name {}".format(new_gushim),status = 1) 
    else:
                print_arcpy_message("there is no New Gush in the tazar",status = 1) 

    arcpy.Delete_management(layer_finish)


def Sub_Processing(bankal,modad_c,points,pnt_modad,Lines,line_modad,border,copy_tazar,gdb):
    
    def Check_If_Sub_Processing(after_del,bankal):
        No_geo_changes = True
        dic_AREA_ID = {i.PARCEL_ID:i.SHAPE_Area for i in arcpy.SearchCursor(after_del)}
        with arcpy.da.SearchCursor(bankal,['PARCEL_ID','SHAPE@AREA']) as cursor:
            for row in cursor:
                if dic_AREA_ID.get(row[0]):
                    changed = abs(row[1] - dic_AREA_ID[row[0]])
                    if changed > 0.1:
                        No_geo_changes = False
                        break
        return No_geo_changes

    after_del  = gdb + '\\' + 'afetr_del'

    if not arcpy.Exists(after_del):
        Delete_polygons (bankal,border,after_del)

    No_geo_changes = Check_If_Sub_Processing(after_del,bankal)
    arcpy.Delete_management (after_del)

    if No_geo_changes:
        print_arcpy_message     ("Sub Processing was acctiveted",status = 1)
        path_before      = gdb + '\\' + 'Copy_bankal'
        arcpy.Select_analysis (bankal,path_before)


        # insert: polygons
        Update_Polygons (bankal,modad_c,'')

        Parcel_data                (bankal,path_before,copy_tazar)
          
        # insert: Points
        Delete_polygons          (points,border)
        arcpy.Append_management  (pnt_modad,points,"NO_TEST")


                # מחיקת כפילויות אם יש
        add_field                 (points,'X_Y','TEXT')
        arcpy.CalculateField_management(points, "X_Y", "calc(!shape.centroid.X!,!shape.centroid.Y!)", "PYTHON_9.3", "def calc (x,y):\\n    return str(round(x,1)) + '-' + str(round(y,1))")

        del_identical   (points,'X_Y')

        # insert Lines

        Layer_Management        (Lines).Select_By_Location('HAVE_THEIR_CENTER_IN',border)

        arcpy.Append_management (line_modad,Lines,"NO_TEST")

        # Delete Duplicate Lines
        Delete_Duplic_Line(Lines)

        print_arcpy_message     ("# # # # # # # F I N I S H # # # # # #",status = 1)
        sys.exit()


def CheckIfSkipProcess(parcel_bankal,PARCELS_inProc_edit,gdb):

    conti = True
    intersect  = gdb + '\\' + 'inter'
    Bankal_cut = gdb + '\\' + 'cut_bankal'
    Layer_Management(parcel_bankal).Select_By_Location("INTERSECT",PARCELS_inProc_edit,"",Bankal_cut)
    name_ID    = 'FID_' + os.path.basename(Bankal_cut)
    arcpy.Intersect_analysis([Bankal_cut,PARCELS_inProc_edit],intersect)
    dic      = {i[0]:[round(i[1],1),i[2],i[3]] for i in arcpy.da.SearchCursor(intersect,[name_ID,'SHAPE_Area','GUSH_NUM','PARCEL'])}
    list_ref = [[i.OBJECTID,round(i.SHAPE_Area,1)] for i in arcpy.SearchCursor(Bankal_cut)]
    for n in list_ref:
        if dic.get(n[0]):
            if n[1] == dic[n[0]][0] and dic[n[0]][1] == dic[n[0]][2]:
                print (n)
                pass
            else:
                conti = False
    arcpy.Delete_management(intersect)
    arcpy.Delete_management(Bankal_cut)
    return conti

def PrePare_Data(parcel_bankal,parcels_copy,points_copy,Point_bankal,GDB,name_bankal,name_tazar):

    '''
    INPUTS
    1) parcel_bankal - שכבת החלקות של הבנק"ל
    2) parcels_copy  - שכבת החלקות של המודד
    3) points_copy   - שכבת נקודות המודד
    4) Point_bankal  - שכבת נקודות הבנק"ל
    5) GDB           - בסיס הנתונים בו ישמרו השכבות
    6) name_bankal   - שם השדה של שם הנקודה בבנק"ל
    7) name_bankal   - שם השדה של שם הנקודה בתצ"ר

    OUTPUTS
    1) AOI               - אזור העבודה החדש
    2) tazar_border      - גבול התצ"ר
    3) Curves            - קשתות של אזור העבודה
    4) parcel_Bankal_cut - חיתוך של הבנק"ל כל חלקה בטווח 10 מטר מהתצ"ר תיכנס
    5) Point_bankal_Cut  - חיתוך של נקודות הבנק"ל, כל נקודה בטווח 10 מטר מהתצ"ר תיכנס
    '''

    # #Prepare data

    parcel_Bankal_cut  = GDB + '\\' + 'Bankal_Cut'
    tazar_border       = GDB + '\\' + 'tazar_border'
    AOI                = GDB + '\\' + 'AOI'
    Point_bankal_Cut   = GDB + '\\' + 'Point_bankal_Cut'
    Holes_data         = GDB + '\\' + 'Holes_Prepare_data'


    # Create Tazar Border, Curves
    arcpy.Dissolve_management                  (parcels_copy,tazar_border)

    # Create Parcel Bankal For AOI
    Layer_Management(parcel_bankal).Select_By_Location("INTERSECT",tazar_border,"10 Meters",parcel_Bankal_cut)

    add_field                       (parcel_Bankal_cut, "AREA_Orig","DOUBLE")
    arcpy.CalculateField_management (parcel_Bankal_cut, "AREA_Orig","!shape.area!", "PYTHON_9.3")

    # Cut Points From Bankal
    Layer_Management(Point_bankal).Select_By_Location("INTERSECT",tazar_border,"10 Meters",Point_bankal_Cut)

    Move_Vertices_By_Name                      (parcel_Bankal_cut,Point_bankal_Cut,name_bankal,points_copy,name_tazar) # לשים לב לשדות שמות הנקודות

    Delete_polygons         (parcel_Bankal_cut,parcels_copy,AOI)
    arcpy.Append_management (parcels_copy,AOI,'NO_TEST')

    # מחיקה של חלקים הקטנים מ-20 אחוז של גודלם המקורי
    Multi_to_single                         (AOI)
    arcpy.AddField_management               (AOI, "OVERLAP_PRCT", "DOUBLE")
    arcpy.CalculateField_management         (AOI,"OVERLAP_PRCT", "((!SHAPE_Area!  / !AREA_Orig!) * 100)", "PYTHON")
    arcpy.MakeFeatureLayer_management       (AOI,'parcel_Bankal_cut_Layer',"\"OVERLAP_PRCT\" < 20") 
    arcpy.Select_analysis                   ('parcel_Bankal_cut_Layer',Holes_data)
    arcpy.DeleteFeatures_management         ('parcel_Bankal_cut_Layer')

    # Update_Polygons                         (AOI,parcels_copy)
    Curves = generateCurves                 (AOI)
    Update_Polygons                         (AOI,parcels_copy)

    Multi = Multi_to_single                         (AOI)
    if Multi:
        print_arcpy_message("You have Multi layer after insert new tazar")

    return AOI,tazar_border,Curves,parcel_Bankal_cut,Point_bankal_Cut


def CheckResultsIsOK(parcel_all,tazar_border,num):

    GDB       = os.path.dirname(tazar_border)
    Holes     = GDB +'\\' + 'Holes_Check_'     + str(num)
    Intersect = GDB +'\\' + 'Intersect_Check_' + str(num)

    holes,intersect = topology         (parcel_all)

    if holes:
        Layer_Management    (holes).Select_By_Location   ('INTERSECT',tazar_border,0,Holes)
        Delete_By_area        (Holes,0.1)
        holes_count         = int(str(arcpy.GetCount_management(Holes)))
        print_arcpy_message  ("holes: {}".format(holes_count,1))
    else:
        holes_count = 0
    
    if intersect:
        Layer_Management      (intersect).Select_By_Location  ('INTERSECT',tazar_border,0,Intersect)
        Delete_By_area        (Intersect,0.1)
        intersect_count       = int(str(arcpy.GetCount_management(Intersect)))
        print_arcpy_message   ("Intersect: {}".format(intersect_count),1)
    else:
        intersect_count = 0

    return True if ((holes_count == 0)  and (intersect_count == 0)) else False



def Calculate_Area_Rashum(PARCEL_ALL_FINAL):
    
    def find_problem(Area_rasum,Shape_area,delta):

        minus = abs(Area_rasum - Shape_area)
        return 'Warning, Delta is to big' if minus > delta else 'Ok'

    def math_delta_rashum(area_rashum):

        area_rashum = float(area_rashum)
        delta1 = (0.3 * (math.sqrt(area_rashum)) + (0.005 * area_rashum))
        delta2 = (0.8 * (math.sqrt(area_rashum)) + (0.002 * area_rashum))

        return delta1 if delta1 > delta2 else delta2

    fields = [add_field(PARCEL_ALL_FINAL,i[0],i[1]) for i in [["GAP", "DOUBLE"],["delta", "DOUBLE"],["Check", "TEXT"]]]
        
    with arcpy.da.UpdateCursor(PARCEL_ALL_FINAL,["LEGAL_AREA","SHAPE_Area","GAP","delta","Check"]) as cursor:
        for row in cursor:
            if row[0]:
                delta  = math_delta_rashum(row[0])
                row[3] = delta
                row[2] = abs(row[1] - row[0])- delta
                row[4] = find_problem(row[0],row[1],delta)
                cursor.updateRow (row)
    del cursor


def Get_Attr_From_parcel(parcel_all_final,tazar_copy):

    Uni_data = len(list(set([i.LOCALITY_ID for i in arcpy.SearchCursor(parcel_all_final) if i.LOCALITY_ID])))
    fields = [["REGION_NAME","TEXT"],["REGION_ID","LONG"],["COUNTY_NAME","TEXT"],["COUNTY_ID","LONG"],["REG_MUN_ID","LONG"],["LOCALITY_NAME","TEXT"],["LOCALITY_ID","LONG"],["REG_MUN_NAME","TEXT"],["WP","LONG"]]
    for i in fields: add_field(tazar_copy,i[0],i[1])
    if Uni_data == 1:
        data = [[i.REGION_NAME,i.REGION_ID,i.COUNTY_NAME,i.COUNTY_ID,i.REG_MUN_ID,i.LOCALITY_NAME,i.LOCALITY_ID,i.REG_MUN_NAME,i.WP] for i in arcpy.SearchCursor(parcel_all_final) if i.LOCALITY_ID][0]
        with arcpy.da.UpdateCursor(tazar_copy,["REGION_NAME","REGION_ID","COUNTY_NAME","COUNTY_ID","REG_MUN_ID","LOCALITY_NAME","LOCALITY_ID","REG_MUN_NAME","WP"]) as Ucursor:
            for row in Ucursor:
                row[:] = data[:]
                Ucursor.updateRow(row)
    else:
        print_arcpy_message ("Get_Attr_From_parcel coudnt give names",2)

def connect_parcel_to_sett(layer,sett,bankal_c,sett_fields = ['MUN_HEB','MACHOZ','NAFA1','SETTEL_NAM']):

    '''
    layer       = שכבה אליה יכנסו מספר מזהה ושם הישוב
    sett        = שכבת הישובים של מפ"י
    sett_fields = 1) CODE field , 2) NAME field
    '''

    list_fields   = [['REG_MUN_NAME','TEXT'],['REGION_NAME','TEXT'],['COUNTY_NAME','TEXT'],['LOCALITY_NAME','TEXT']]

    fields_name   = [i[0] for i in list_fields]

    add_op = [add_field(layer,i[0],i[1]) for i in list_fields]

    if sett != '':
        sett_fields.insert                      (0, 'SHAPE@') 
        fields_name.insert                      (0, 'SHAPE@')
        arcpy.MakeFeatureLayer_management       (sett,'sett_layer')
        arcpy.SelectLayerByLocation_management  ('sett_layer','INTERSECT',layer)
        data = [[i[0],i[1],i[2],i[3],i[4]] for i in arcpy.da.SearchCursor('sett_layer',sett_fields)]
        with arcpy.da.UpdateCursor (layer, fields_name) as cursor:
            for row in cursor:
                geom   = row[0]
                midpnt = geom.labelPoint
                for i in data:
                    if i[0].distanceTo(midpnt) == 0:
                        row[1:] = i[1:]
                        cursor.updateRow(row)

        del cursor

    arcpy.JoinField_management (layer, 'REGION_NAME',bankal_c , 'REGION_NAME',['REG_MUN_ID','REGION_ID','COUNTY_ID','LOCALITY_ID'])


def get_no_node_vertex(AOI,tazar_border,Modad_node,PARCEL_ALL_node):

    '''
    [INFO] -  בודק אם המודד שכח לתת נקודת חיצונית, ולכן אין אנו יכולים לחבר אותה
    Inputs:
    1) AOI             - Result of the tool after runing
    2) tazar_border    - border of the tazar
    3) Modad_node      - points coming from the modad
    4) PARCEL_ALL_node - Point of the PARCAL_ALL before starting to change him
    OUTPUT:
    1) points layer of missing vertxes
    '''

    GDB = os.path.dirname(AOI)

    print_arcpy_message                       ("START Func: get no node vertex",1)

    Point_AOI = Layer_Management(AOI).Get_vertxs_As_Point ()

    arcpy.MakeFeatureLayer_management         (Point_AOI,"AOI_lyr")
    arcpy.SelectLayerByLocation_management    ("AOI_lyr","BOUNDARY_TOUCHES",tazar_border,0.01)
    arcpy.SelectLayerByLocation_management    ("AOI_lyr","INTERSECT",Modad_node,0.01,"REMOVE_FROM_SELECTION")
    arcpy.SelectLayerByLocation_management    ("AOI_lyr","INTERSECT",PARCEL_ALL_node,0.01,"REMOVE_FROM_SELECTION")
    arcpy.Select_analysis                     ("AOI_lyr",GDB + "\\" + "Possible_Error_points")

    return GDB + "\\" + "Possible_Error_points"



def clean_slivers_by_vertex(PARCEL_ALL,SLIVERS_CLEAN,border,Dis_search,PARCEL_ALL_lyr):

    print_arcpy_message ("START Func: clean slivers by vertex")

    '''
    [INFO] -  מוחק לפי ליניאריות ומרחק את הוורטקסים שנמצאים ליד החורים של התצ"ר
    INPUT-
    1) PARCEL_ALL     - שכבת רצף
    2) SLIVERS_CLEAN  - חורים של הרצף
    3) border         - גבול התצ"ר
    4) Dis_search     - מרחק חיפוש הוורטקסים
    5) PARCEL_ALL_lyr - שכבת הפלט
    '''

    gdb = os.path.dirname(border)

    tazar_border = 'in_memory\TazarBorderDiss'
    arcpy.Dissolve_management (border,tazar_border)

    conn = sqlite3.connect(':memory:')
    c    = conn.cursor()
    c.execute('''CREATE TABLE old_vertices(pnt_num real, x real, y real, xy text, part real, oid real)''')
    c.execute('''CREATE TABLE new_vertices(pnt_num real, x real, y real, xy text, part real, oid real)''')


    c.execute('''CREATE TABLE vertices(pnt_num real, x real, y real, xy text, part real, oid real, junction real, linearity real)''')

    c.execute('''CREATE TABLE sliver_vertices(pnt_num real, x real, y real, xy text, part real, oid real, junction real, linearity real)''')

    c.execute('''CREATE TABLE border_vertices(pnt_num real, x real, y real, xy text, part real, oid real, junction real, linearity real)''')

    arcpy.Select_analysis            (PARCEL_ALL, PARCEL_ALL_lyr)
    
    arcpy.CopyFeatures_management        (PARCEL_ALL_lyr,gdb + "\\PARCEL_ALL_lyr_COPY_DEL")
    
    VerticesToTable2(PARCEL_ALL_lyr, "vertices",c)
    VerticesToTable2(SLIVERS_CLEAN, "sliver_vertices",c)
    VerticesToTable2(tazar_border, "border_vertices",c)


    parcel_common_vertices = [row for row in c.execute('''SELECT * FROM vertices
                                                                                                left join sliver_vertices
                                                                                                on vertices.xy = sliver_vertices.xy
                                                                                                where  sliver_vertices.xy is not null''')]
                                                                                            
    border_common_vertices = [row for row in c.execute('''SELECT * FROM border_vertices
                                                                                                left join sliver_vertices
                                                                                            on border_vertices.xy = sliver_vertices.xy
                                                                                                where  sliver_vertices.xy is not null''')]
                                                                                                

    distance_vertices = [[p[:8] + b[:8],math.sqrt(((p[1]-b[1])**2)+((p[2]-b[2])**2))] for p in parcel_common_vertices for b in border_common_vertices if math.sqrt(((p[1]-b[1])**2)+((p[2]-b[2])**2))\
         < Dis_search or (float("{0:.2f}".format(p[1])) == float("{0:.2f}".format(b[1])) and float("{0:.2f}".format(p[2])) == float("{0:.2f}".format(b[2])))]

                                           
    rows = arcpy.UpdateCursor(PARCEL_ALL_lyr)
    for row in rows:
            geometry = row.Shape
            oid = row.OBJECTID
            pts = []
            poly_vertices = [r for r in distance_vertices if r[0][5] == oid]
            for part in geometry:
                    for pt in part:
                            if str(type(pt)) != "<type 'NoneType'>":
                                    num_point = 0
                                    #print str(pt.X) + "--" + str(pt.Y)
                                    this_x = float("{0:.2f}".format(pt.X))
                                    this_y = float("{0:.2f}".format(pt.Y))      
                                    this_vertex = [p for p in poly_vertices if float("{0:.2f}".format(p[0][1])) == this_x and float("{0:.2f}".format(p[0][2])) == this_y]
                                    if this_vertex:
                                            if this_vertex[0][0][8] == None:
                                                    if this_vertex[0][0][7] < 0.7 and this_vertex[0][0][6] == 1:
                                                        print ("pseodo: delete vertex")
                                                    else:
                                                            #print "pseodo, but important: keep the vertex"
                                                            point = pt
                                                            pts.append(point)
                                            # tazar point in buffer
                                            else:
                                                    # check minimum distance
                                                    the_minimum_vertex = [v for v in this_vertex if v[1] == min([i[1] for i in this_vertex])]
                                                    point = arcpy.Point(the_minimum_vertex[0][0][9], the_minimum_vertex[0][0][10])
                                                    pts.append(point)
                                    # point not on sliver: keep the vertex
                                    else:
                                            point = pt
                                            pts.append(point)
                                    if num_point == 0:
                                            first_point = point
                                    num_point = num_point + 1
            polygon = PtsToPolygon(pts)
            if pts[0] != pts[-1] and first_point:
                    #print "ooops.... - polygon not closed"
                    pts.append(first_point)
            row.Shape       = polygon
            rows.updateRow(row)

    arcpy.Delete_management(gdb + "\\PARCEL_ALL_lyr_COPY_DEL")
    return PARCEL_ALL_lyr




def get_default_Snap_border(point_bankal,tazar,Distance_min):

    '''
    [INFO] - בודק את נקודות הבנקל ליד התצ"ר במידה ויש 2 נקודות בנקל קרובות אחת לשניה, הוא נותן את המרחק הזה כברירת מחדל
    INPUT-
    1) point_bankal - שכבת נוקודת בנקל
    2) tazar        - שכבת תצ"ר של המודד
    3) Distance_min - מה המינימום, יופעל במידה שמרחק הנקודות גדול מידי
    OUTPUT-
    1) המרחק בקטן ביותר בין שני נקודות בנק"ל בסמוך לתצר
    '''

    GDB = os.path.dirname(tazar)

    PntTmp = r'in_memory' + '\\' + 'PntTmp'
    buffer = r'in_memory' + '\\' + 'buffer'
    dissol = r'in_memory' + '\\' + 'dissol'
    multiP = r'in_memory' + '\\' + 'multiP'

    def Getmin(list1,Dis_min = 2):
        li = [i[2] for i in list1 if i[2] < 1]
        return min(li) - 0.01 if li else Dis_min

    arcpy.MakeFeatureLayer_management           (point_bankal,'path_lyr')
    arcpy.SelectLayerByLocation_management      ('path_lyr','WITHIN_A_DISTANCE',tazar,'5 Meters')
    arcpy.Select_analysis                       ('path_lyr',PntTmp)
    arcpy.MakeFeatureLayer_management           (PntTmp,'PntTmp_lyr')
    arcpy.SelectLayerByLocation_management      ('PntTmp_lyr',"COMPLETELY_WITHIN",tazar)
    arcpy.SelectLayerByAttribute_management     ('PntTmp_lyr',"SWITCH_SELECTION")

    arcpy.Buffer_analysis                   ('PntTmp_lyr',buffer,0.5)
    arcpy.Dissolve_management               (buffer,dissol)
    arcpy.MultipartToSinglepart_management  (dissol,multiP)


    with arcpy.da.UpdateCursor(multiP,['SHAPE@AREA']) as cursor:
        for row in cursor:
            if row[0] < 0.8:
                cursor.deleteRow()

    arcpy.MakeFeatureLayer_management       (PntTmp,'path2_lyr')
    arcpy.SelectLayerByLocation_management  ('path2_lyr','INTERSECT',multiP)

    dis_point  = [[row[0],row[1]] for row in arcpy.da.SearchCursor('path2_lyr',['OBJECTID','SHAPE@'])]
    list_dis   = [[row[1],n[0],row[0].distanceTo(n[1])] for n in dis_point for row in arcpy.da.SearchCursor('path2_lyr',['SHAPE@','OID@']) if row[0].distanceTo(n[1]) > 0]


    Min_dist = Getmin(list_dis,Distance_min)
    print_arcpy_message(Min_dist, status=1)
    return Min_dist



def ChangeFieldNames(parcel,line,point):
    '''
        Take 3 layers, Changing fields from Source layers to bankal format
    '''

    wrong = {'TALAR_NUM':'TALAR_NUMBER','GushNum':'GUSH_NUM','GushSuffix':'GUSH_SUFFIX','ParcelName':'PARCEL','LegalArea':'LEGAL_AREA','GUSHNUM':'GUSH_NUM','GUSHSUFFIX':'GUSH_SUFFIX','PARCEL_FINAL':'PARCEL','LEGALAREA':'LEGAL_AREA'}

    List_fields = [[str(i.name),wrong[str(i.name)]] for i in arcpy.ListFields(parcel) if str(i.name) in list(wrong.keys())]

    if List_fields:

        print_arcpy_message("Changing Fields",status = 1)
        list_layers = [parcel,line,point]

        for lyr in list_layers:
            for field in List_fields:
                layer   = os.path.basename(lyr)
                parcels = [os.path.basename(parcel)]
                others  = [os.path.basename(line),os.path.basename(point)]
                if layer in others:
                    if field[0] in ['GushNum','GUSHNUM']:
                        add_field(lyr,field[1],'LONG')
                        arcpy.CalculateField_management  (lyr, field[1]  , "!"+field[0]+"!"  , "PYTHON", "")
                    if field[0] in ['GushSuffix','GUSHSUFFIX']:
                        add_field(lyr,field[1],'SHORT')
                        arcpy.CalculateField_management  (lyr, field[1]  , "!"+field[0]+"!"  , "PYTHON", "")
                if layer in parcels:
                    if field[0] in ['GushNum','GUSHNUM']:
                        add_field(lyr,field[1],'LONG')
                        arcpy.CalculateField_management  (lyr, field[1]  , "!"+field[0]+"!"  , "PYTHON", "")

                    if field[0] in ['GushSuffix','ParcelName','PARCEL_FINAL','GUSHSUFFIX']:
                        add_field(lyr,field[1],'SHORT')
                        arcpy.CalculateField_management  (lyr, field[1]  , "!"+field[0]+"!"  , "PYTHON", "")

                    if field[0] in ['LegalArea','LEGALAREA']:
                        add_field(lyr,field[1],'DOUBLE')
                        arcpy.CalculateField_management  (lyr, field[1]  , "!"+field[0]+"! * 1000"  , "PYTHON", "")

                    if field[0] in ['TALAR_NUM']:
                        add_field(lyr,field[1],'LONG')
                        arcpy.CalculateField_management  (lyr, field[1]  , "!"+field[0]+"!"  , "PYTHON", "")

    try:
        arcpy.CalculateField_management  (parcel, 'PARCEL', "int( ''.join ([i for i in !ParcelName! if i.isdigit()]))", "PYTHON" ) 
    except:
        arcpy.CalculateField_management  (parcel, 'PARCEL', "int( ''.join ([i for i in !PARCEL_FINAL! if i.isdigit()]))", "PYTHON" ) 


def add_err_pts_to_mxd(our_gdb, folder, data_source,CURRENT):

    # copy 3 error fcs from data_source (demo.gdb) to our_gdb
    err_fc_names = ["Errors_Line", "Errors_Point", "Errors_Polygon"]
    for err_fc_name in err_fc_names:
        arcpy.DeleteRows_management(data_source + "\\" + err_fc_name)
        arcpy.Copy_management(data_source + "\\" + err_fc_name, our_gdb + "\\" + err_fc_name)
    
    mxd = arcpy.mp.ArcGISProject (CURRENT)
    df  = mxd.listMaps('Layers')[0]
    for root, dir, files in os.walk(folder):
        for file in files:
            file_full_path  = root + "\\" + file
            if file == "Errors_Line.lyr" or file == "Errors_Point.lyr" or file == "Errors_Polygon.lyr" or file == "Possible_Error_points.lyr" or file == "PARCEL_ALL_EDIT_copy.lyr" or file == "PARCEL_NODE_EDIT_copy.lyr" or file == "PARCEL_ARC_EDIT_copy.lyr":
                addLayer   = arcpy.mp.LayerFile(file_full_path)
                df.addLayer(addLayer, "TOP")
                try:
                    mxd.updateConnectionProperties(data_source, our_gdb)
                except:
                    print ("Coudnt replace Data Source")


def Parcel_data(path_after,path_before,copy_tazar):
    
    def Get_Runing_numbers(data1):
        for i in range(len(data1)):
            if data1[i][1] == data1[i-1][1]:
                if data1[i][0]+1 == data1[i-1][0]:
                    print ("its ok ,in {} value is equal with: {}".format(data1[i][2],data1[i-1][2]))
                else:
                    print ("in {} value is not equal with: {}".format(data1[i][2],data1[i-1][2]))
            else:
                pass


    conn = sqlite3.connect(":memory:")
    c = conn.cursor()
    c.execute("""CREATE TABLE Before_Table (
                        PARCEL     INTEGER,
                        GUSH_NUM   INTEGER,
                        KEY        text
                        )""")

    c = conn.cursor()
    c.execute("""CREATE TABLE Table_After (
                        PARCEL     INTEGER,
                        GUSH_NUM   INTEGER,
                        KEY        text
                        )""")

    for i in arcpy.SearchCursor(path_before):
        c.execute ("INSERT INTO Before_Table VALUES (" + str(i.PARCEL) +','+ str(i.GUSH_NUM) + ",'"+str(i.PARCEL)+"-"+str(i.GUSH_NUM)+"-"+ str(i.GUSH_SUFFIX)+"')")

    for i in arcpy.SearchCursor(path_after):
        c.execute ("INSERT INTO Table_After VALUES (" + str(i.PARCEL) +','+ str(i.GUSH_NUM) + ",'"+str(i.PARCEL)+"-"+str(i.GUSH_NUM) +"-"+ str(i.GUSH_SUFFIX)+ "')")


    count_before = [row for row in c.execute ('''SELECT * FROM  (SELECT *, COUNT(*) as count FROM Before_Table group by KEY) t1 WHERE t1.count > 1;''')]
    count_after  = [row for row in c.execute ('''SELECT * FROM  (SELECT *, COUNT(*) as count FROM Table_After group by KEY) t1 WHERE t1.count > 1;''')]

    if count_before:
        msg  =  "Found identical parcels on orig parcels : {}".format(count_before)
        print_arcpy_message(msg, status=2)

    if count_before:
        msg2 = "Found identical parcels on new parcels : {}".format(count_after)
        print_arcpy_message(msg2, status=2)

    #data1 = [row for row in c.execute ('''SELECT * FROM  Before_Table ORDER BY GUSH_NUM DESC, PARCEL DESC;''')]
    #Get_Runing_numbers(data1)

    add_parcels = [str(row[0]) for row in c.execute ('''SELECT A.KEY FROM Table_After A LEFT JOIN Before_Table B ON A.KEY = B.KEY WHERE B.KEY is NULL;''')]
    del_parcels = [str(row[0]) for row in c.execute ('''SELECT A.KEY FROM Before_Table A LEFT JOIN Table_After B ON A.KEY = B.KEY WHERE B.KEY is NULL;''')]

    gdb = os.path.dirname(path_after)
    Insert_to_table(path_before,copy_tazar,gdb)

    msg2 = "added parcels: {}  ".format(add_parcels)
    msg3 = "Deleted parcels: {}".format(del_parcels)

    print_arcpy_message(msg2, status=1)
    print_arcpy_message(msg3, status=1)

    data = {str(i.PARCEL) +'-' +str(i.GUSH_NUM)+'-'+ str(i.GUSH_SUFFIX):[i.LOCALITY_ID,i.LOCALITY_NAME,i.LEGAL_AREA] for i in arcpy.SearchCursor(copy_tazar)}
    with arcpy.da.UpdateCursor(path_after,['PARCEL','GUSH_NUM','GUSH_SUFFIX','LOCALITY_ID','LOCALITY_NAME','LEGAL_AREA']) as ucursor:
        for row in ucursor:
            key = str(row[0]) +'-' +str(row[1])+'-'+ str(row[2])
            if data.get(key):
                row[3] = data[key][0]
                row[4] = data[key][1]
                row[5] = data[key][2]
                ucursor.updateRow(row)
    del ucursor

    # Make Dic:  {GUSH_NUM:[STATUS,STATUS_TEXT]} and copy this values to PARCEL LAYER if with Nulls
    Get_Status_Field(path_after)



def Insert_to_table(bankal,tazar_copy,GDB):

    arcpy.MakeFeatureLayer_management      (bankal,'bankal_lyr')
    arcpy.SelectLayerByLocation_management ('bankal_lyr','INTERSECT',tazar_copy, '1 Meters')

    None_me = [i for i in arcpy.SearchCursor(tazar_copy) if i.PARCEL_FINAL == None]
    if None_me:
        arcpy.CalculateField_management  (tazar_copy, 'PARCEL_FINAL', "int( ''.join ([i for i in !PARCELNAME! if i.isdigit()]))", "PYTHON" ) 

    new_shamce = False
    try:
        data = [[i.shape,str(i.PARCEL_FINAL)+'-'+str(i.GUSHNUM)+'-'+str(i.GUSHSUFFIX)] for i in arcpy.SearchCursor(tazar_copy)]
    except:
        data = [[i.shape,str(i.PARCEL)+'-'+str(i.GUSH_NUM)+'-'+str(i.GUSH_SUFFIX)] for i in arcpy.SearchCursor(tazar_copy)]
        new_shamce = True

    in_tazar_copy = []
    with arcpy.da.SearchCursor ('bankal_lyr', ['SHAPE@','PARCEL','GUSH_NUM','GUSH_SUFFIX']) as cursor:
        for row in cursor:
            geom   = row[0]
            midpnt = geom.labelPoint
            key    = str(row[1])+'-'+str(row[2])+'-'+str(row[3])
            for i in data:
                if i[0].distanceTo(midpnt) == 0:
                    in_tazar_copy.append([i[1],key])
        del cursor


    data          = [[i.shape,str(i.PARCEL)+'-'+str(i.GUSH_NUM)+'-'+str(i.GUSH_SUFFIX)] for i in arcpy.SearchCursor('bankal_lyr')]
    in_bankal     = []
    fields_parcel = ['SHAPE@','PARCEL_FINAL','GUSHNUM','GUSHSUFFIX']
    if new_shamce == True:
        fields_parcel = ['SHAPE@','PARCEL','GUSH_NUM','GUSH_SUFFIX']

    with arcpy.da.SearchCursor (tazar_copy, fields_parcel) as cursor:
        for row in cursor:
            geom   = row[0]
            midpnt = geom.labelPoint
            key    = str(row[1])+'-'+str(row[2])+'-'+str(row[3])
            for i in data:
                if i[0].distanceTo(midpnt) == 0:
                    in_bankal.append([key,i[1]])
        del cursor

    data1 = [ast.literal_eval(i) for i in list(set([str(i) for i in in_tazar_copy + in_bankal]))]


    path      = GDB
    name1     = 'CANCEL_PARCEL_EDIT'
    full_path = path +'\\'+ name1

    if not arcpy.Exists(full_path):
        arcpy.CreateTable_management(path,name1)

    fields = ['F_PARCEL_NUM','F_GUSH_SUFFIX','T_PARCEL_NUM','T_GUSH_SUFFIX']
    add_field(full_path,"F_GUSH_NUM",'LONG')
    add_field(full_path,"T_GUSH_NUM",'LONG')
    for i in fields:
        add_field(full_path,i,'SHORT')

    for row in data1:
            insert = arcpy.InsertCursor (full_path)
            in_row               = insert.newRow()
            in_row.F_GUSH_NUM    = int(row[1].split('-')[1])
            in_row.F_PARCEL_NUM  = int(row[1].split('-')[0])
            in_row.F_GUSH_SUFFIX = int(row[1].split('-')[2])
            in_row.T_GUSH_NUM    = int(row[0].split('-')[1])
            in_row.T_PARCEL_NUM  = int(row[0].split('-')[0])
            in_row.T_GUSH_SUFFIX = int(row[0].split('-')[2])
            insert.insertRow   (in_row)


def Get_Status_Field(layer):
    data_status = {i.GUSH_NUM:[i.STATUS,i.STATUS_TEXT] for i in arcpy.SearchCursor(layer) if i.GUSH_NUM != None and i.STATUS != None}

    with arcpy.da.UpdateCursor(layer,['GUSH_NUM',"STATUS","STATUS_TEXT"]) as cursor:
        for row in cursor:
            if not row[1]:
                if data_status.get(row[0]):
                    row[1] = data_status[row[0]][0]
                    row[2] = data_status[row[0]][1]
                    cursor.updateRow(row)
    del cursor


def Delete_curves_out_AOI(parcel_new,bankal_old):

    diss_old = r'in_memory' + '\\' + 'Diss_old'
    diss_new = r'in_memory' + '\\' + 'Diss_new'
    lyr_old  = r'in_memory' + '\\' + 'lyr_old'
    lyr_new  = r'in_memory' + '\\' + 'lyr_new'


    arcpy.Dissolve_management    (bankal_old,diss_old)
    arcpy.Dissolve_management    (parcel_new,diss_new)

    Feature_to_polygon           (diss_old,lyr_old)
    Feature_to_polygon           (diss_new,lyr_new)

    Delete_polygons (lyr_new,lyr_old)
    Delete_polygons (parcel_new,lyr_new)



def update_curves(fc,curve):
        for row in arcpy.SearchCursor(curve):
                upd_rows = arcpy.UpdateCursor(fc)
                curve_g = row.Shape
                midpnt = curve_g.centroid
                for upd_row in upd_rows:
                        if upd_row.Shape.distanceTo(midpnt)== 0:
                                diff = upd_row.Shape.difference (curve_g)
                                new_geometry = curve_g.union(diff)
                                upd_row.Shape = new_geometry
                                upd_rows.updateRow(upd_row) 

def names_curves(fc,tazar,curve):

        '''
        [INFO] -  חיבור של הקשתות לשמות של חלקות הבנקל שמתחת, צריך להיות במצב בו הקשתות נמצאות מעל הבנקל
        INPUT:
        1) fc     = שכבת הרצף
        2) curves = קשתות
        OUTPUT:
        1) curves = שכבה חדשה שמכוונת לאיזה שכבת בנקל הוא מתאים
        '''

        add_field                              (curve,'PARCEL_ID','LONG')
        arcpy.MakeFeatureLayer_management      (curve,'curvs_lyr')
        arcpy.SelectLayerByLocation_management ('curvs_lyr',"HAVE_THEIR_CENTER_IN",tazar,'',"NEW_SELECTION","INVERT")

        upd_rows = arcpy.UpdateCursor('curvs_lyr')
        for row in upd_rows:
                for Search_row in arcpy.SearchCursor(fc):
                    try: 
                        if Search_row.Shape.distanceTo(row.Shape.centroid)== 0:
                                if Search_row.PARCEL_ID > 0:
                                    if row.PARCEL_ID == None:
                                        row.PARCEL_ID = Search_row.PARCEL_ID
                                        upd_rows.updateRow(row)
                    except:
                        print ("Coudnt make Names for curves")
                        pass


def Update_Layer_Curves_By_ID(fc,tazar,curve):

    '''
    [INFO] - מחבר את הקשתות לחלקות המתאימות ומעדכן את התצ"ר בהתאם, מתשמש בשמות מפונקציה
            names_curves
    INPUT:
    1) fc     = שכבת הרצף
    2) tazar  =  שכבת התצר
    3) curves =  שכבת הקשתות

    INPUT     = שכבת הרצף מעודכנת עם השמות
    '''

    curve_temp = curve + '_temp'
    curve_diss = r'in_memory\\curve_diss'

    add_field      (curve,'PARCEL_ID','LONG')
    Delete_polygons(fc,curve)
    Update_Polygons(fc,tazar)
    Fix_curves     (fc,tazar,curve)


    Search_data = {row.PARCEL_ID:row.Shape for row in arcpy.SearchCursor(curve) if row.PARCEL_ID}

    # חיבור בין הבנקל לשכבת הקשתות

    arcpy.Select_analysis    (curve,curve_temp,"\"PARCEL_ID\" > 0 ")
    arcpy.Dissolve_management(curve_temp,curve_diss,['PARCEL_ID'])
    Search_data = {row.PARCEL_ID:row.Shape for row in arcpy.SearchCursor(curve_diss) if row.PARCEL_ID}

    Delete_polygons(fc,curve_temp)

    upd_rows = arcpy.UpdateCursor(fc)
    for row in upd_rows:
        if Search_data.get(row.PARCEL_ID):
            row.shape = Search_data[row.PARCEL_ID].union(row.Shape)
            upd_rows.updateRow(row) 


def Fix_Multi_part_Bankal(layer,tazar_border,parcel_Bankal_cut):

	'''
	[INFO] -  מתקן חלקות בנקל שמחוברות לתצר אבל נהרסו במהל העבודה, יתקן על מה ש-250 מטר מהתצר

	INPUT:
	1) layer             = שכבת העבודה אחרי עריכה
	2) tazar_border      = גבול התצ"ר
	3) parcel_Bankal_cut = בנקל חתוך לפני עריכה

	OUTPUT:
	1) שכבת עבודה בה כל מה ש-250 מטר מהתצר חזר להיות כמו המקור
	'''

	diss_temp        = r'in_memory'     +'\\'+'diss_temp'
	Temp_inter       = layer            +'Temp'
	after_del        = r'in_memory'     +'\\'+'after_del'
	Multi_part_inter = layer            +'Temp2'
	save_name        = layer

	arcpy.Buffer_analysis     (tazar_border,diss_temp,250)

	Delete_polygons            (parcel_Bankal_cut,diss_temp,after_del)

	arcpy.Clip_analysis       (layer,diss_temp,Temp_inter)

	fields       = Layer_Management(Temp_inter).fields()
	fields_layer = [n for n in fields if n not in ["SHAPE_Area","SHAPE_Length","OBJECTID","SHAPE@","SHAPE","GAP","delta","Check"]]

	arcpy.Dissolve_management (Temp_inter,Multi_part_inter,fields_layer)

	data = {str(i.PARCEL) +'-' +str(i.GUSH_NUM)+'-'+ str(i.GUSH_SUFFIX):i.shape for i in arcpy.SearchCursor(after_del)}

	with arcpy.da.UpdateCursor(Multi_part_inter,['PARCEL','GUSH_NUM','GUSH_SUFFIX','SHAPE@']) as cursor:
		for row in cursor:
			geom = row[-1]
			key  = str(row[0]) +'-' +str(row[1])+'-'+ str(row[2])
			if data.get(key):
				row[-1] = geom.union(data[key])
				cursor.updateRow(row)


	arcpy.Delete_management                (Temp_inter)
	arcpy.Delete_management                (layer)
	arcpy.Rename_management                (Multi_part_inter,save_name)




def Get_Point_AOI(AOI_final,point_bankal,point_modad,AOI_Point):

    '''
    [INFO] - מייצר שכבת נקודות מתאימה ל
    AOI Final
    חיבור של שכבת הבנקל לשכבת המודד, ללא הנקודות שנמחקו
    '''
    Layer_Management(point_bankal).Select_By_Location('INTERSECT',AOI_final,'0.001 Meters',AOI_Point)

    layer1 = Layer_Management    (AOI_Point)
    Fix_Pnt_Tolerance            (AOI_final,AOI_Point)
    layer1.Select_By_Location    ('COMPLETELY_WITHIN',AOI_final)

    pnt_save = [str(round(pt.X,1)) + '-' + str(round(pt.Y,1)) for row in arcpy.SearchCursor(AOI_final) for part in row.shape for pt in part if pt]
    with arcpy.da.UpdateCursor(AOI_Point,['SHAPE@']) as cursor:
        for row in cursor:
            if str(round(row[0].centroid.X,1)) + '-' + str(round(row[0].centroid.Y,1)) not in pnt_save:
                cursor.deleteRow()

    arcpy.Append_management      (point_modad,AOI_Point,'NO_TEST')

    layer1.delete_identical()
    arcpy.DeleteField_management(layer1.layer,'X_Y')

    return AOI_Point

    

def delete_Line_by_polygon(AOI_Line,tazar_border,Dissolve = False,num = 0.001):

    save_name = tazar_border

    if Dissolve:
        tazar_diss = 'in_memory\\' + os.path.basename(tazar_border+"_Temp")
        arcpy.Dissolve_management(tazar_border,tazar_diss)
        save_name = tazar_border
    
    if int(str(arcpy.GetCount_management(save_name))) > 0:
        data = [i.shape for i in arcpy.SearchCursor(save_name) if i.shape][0]
        with arcpy.da.UpdateCursor(AOI_Line,['SHAPE@']) as cursor:
            for row in cursor:
                if row[0]:
                    geom      = row[0]
                    new_geom  = geom.difference(data.buffer(num))
                    row[0]    = new_geom
                    cursor.updateRow(row)

def Create_Line_AOI(aoi,tazar_border,curves,bankal_line,modad_line,New_Line):

    print_arcpy_message('START Func: Create_Line_AOI',1)

    gdb = os.path.dirname(tazar_border)
    bankal_cut   = gdb + '\\' + 'bankal_cut'
    curves_temp  = gdb + '\\' + 'curves_temp'
    Return_line  = gdb + '\\' + 'Return_line'

    Polygon_To_Line_holes    (aoi,New_Line)
    Split_Line_By_Vertex     (New_Line)
    delete_Line_by_polygon   (New_Line,tazar_border)

    Layer_Management         (curves).Select_By_Location('COMPLETELY_WITHIN',tazar_border,'1 Meters',curves_temp,'invert')

    Layer_Management         (New_Line).Select_By_Location('INTERSECT',curves_temp)

    Layer_Management         (bankal_line).Select_By_Location ('INTERSECT',aoi,0,bankal_cut)

    Layer_Management         (bankal_cut).Select_By_Location ('INTERSECT',curves_temp,0,Return_line)  # החזרה של הקווים שנמחקו בעקבות הקשתות

    Layer_Management         (bankal_cut).Select_By_Location ('INTERSECT',aoi,'1 Meters',None,'invert')

    delete_Line_by_polygon   (bankal_cut,tazar_border,False,1)

    arcpy.Append_management  (bankal_cut,New_Line,'NO_TEST')
    Layer_Management         (Return_line).Select_By_Location('INTERSECT',tazar_border,0,None,'invert') 

    arcpy.Append_management  (Return_line,New_Line,'NO_TEST')

    Multi_to_single          (New_Line)

    del_line_Not_on_parcels  (New_Line,aoi)

    Delete_Duplic_Line       (New_Line)

    fix_tolerance_line       (New_Line,tazar_border)

    Delete_layers_after_use([bankal_cut,curves_temp])


    return New_Line


def Find_stubbern_lines(bankal_arc,aoi,tazar_border):

    gdb = os.path.dirname(tazar_border)

    Line_bankal = gdb + '\\' + 'Line_bankal_cut'
    new_line    = gdb + '\\' + 'AOI_Bankal_line'

    Layer_Management         (bankal_arc).Select_By_Location('INTERSECT',aoi,0,Line_bankal)
    Connect_Lines            (Line_bankal,new_line,50)

    # מחיקות קווים שעלולים להיות בעיתיים
    delete_Line_by_polygon   (new_line,tazar_border)
    Layer_Management         (new_line).Select_By_Location('INTERSECT',tazar_border,'0.1 Meters',None,'invert')
    Layer_Management         (new_line).Select_By_Location("WITHIN_CLEMENTINI",aoi)

    arcpy.Append_management  (new_line,bankal_arc,'NO_TEST')

    Delete_layers_after_use([new_line,Line_bankal])




def add_field(fc,field,Type = 'TEXT'):
    try:
        TYPE = [i.name for i in arcpy.ListFields(fc) if i.name == field]
        if not TYPE:
            arcpy.AddField_management (fc, field, Type, "", "", 500)
    except:
        arcpy.AddField_management (fc, field, Type, "", "", 500)

def del_identical(points,field):

    before = int(str(arcpy.GetCount_management(points)))

    data       = [[row[0],row[1]] for row in arcpy.da.SearchCursor(points,["OBJECTID",field])]

    df         = pd.DataFrame(data,columns= ["OBJECTID",field])
    df["RANK"] = df.groupby(field).rank(method='first',ascending=False)
    df         = df[df['RANK'] > 1]

    data_to_gis = []
    for row in df.itertuples(index=True, name='Pandas'):
            data_to_gis.append([getattr(row, "OBJECTID")])

    flat_list = [item for sublist in data_to_gis for item in sublist]
    

    with arcpy.da.UpdateCursor(points,["OBJECTID"]) as cursor:
            for row in cursor:
                    if int(row[0]) in flat_list:
                            cursor.deleteRow()
    del cursor
                            

    after = int(str(arcpy.GetCount_management(points)))
    deleted = before - after

    return deleted



def Split_List_by_value(list1,value,del_value = False):

    list_index = [n for n,val in enumerate(list1) if val == value]

    list_index.append(len(list1))

    list_val = []
    num = 0
    for i in list_index:
        list_val.append(list1[num:i])
        num = + i

    if del_value:
        for i in list_val:
            for n in i:
                if n is None:
                        i.remove(value)
    return list_val


def Feature_to_polygon(path,Out_put):

    path_diss = arcpy.Dissolve_management(path,r'in_memory\path_diss')

    polygon = []
    cursor = arcpy.SearchCursor(path_diss)
    for row in cursor:
        geom = row.shape
        for part in geom:
            num = 0
            for pt in part:
                if pt:
                    polygon.append([pt.X,pt.Y])
                else:
                    polygon.append(None)

    poly    = Split_List_by_value(polygon,None,True)            
    feature = arcpy.CopyFeatures_management(path,Out_put)

    for i in poly[1:]:
        array = arcpy.Array()
        for n in i:
            array.add(arcpy.Point(n[0],n[1]))
        polygon      = arcpy.Polygon(array, arcpy.SpatialReference("Israel TM Grid"))
        in_rows      = arcpy.InsertCursor(feature)
        in_row       = in_rows.newRow()
        in_row.Shape = polygon
        in_rows.insertRow(in_row)
        
    arcpy.RepairGeometry_management(Out_put)
    return Out_put   


# # # # # # # # # # # # # # # Layer Class # # # # # # # # # # # # # # # 

class Layer_Management():

    '''
    1)  Calc_XY
    2)  Select_By_Location
    3)  Multi_to_single
    4)  add_field
    5)  is_Curves
    6)  Erase
    7)  Destroy_layer
    8)  Get_vertxs_As_Point
    9)  None_in_fields
    10) Get_Closest_Distance
    11) Fill_Holes_in_Polygon
    12) Get_Label_Point_As_Point
    13) delete_identical
    '''

    def __init__(self,Layer):
        if arcpy.Exists(Layer):
            self.gdb          = os.path.dirname  (Layer)
            self.name         = os.path.basename (Layer)
            self.layer        = Layer
            self.oid          = str(arcpy.Describe(Layer).OIDFieldName)

            desc = arcpy.Describe(Layer)
            if str(desc.shapeType) == 'Point':
                self.Geom_type = 'Point'
            elif str(desc.shapeType) == 'Polyline':
                self.Geom_type = 'Polyline'
            else:
                self.Geom_type = 'Polygon'
        else:
            print ("Layer is not exist")
            pass

    def Calc_XY(self):
        add_field                       (self.layer,'X_Y','TEXT')
        with arcpy.da.UpdateCursor(self.layer,['SHAPE@','X_Y']) as Ucursor:
            for row in Ucursor:
                row[1] = str(round(row[0].centroid.X,1)) + '-' + str(round(row[0].centroid.Y,1))
                Ucursor.updateRow(row)
        return 'X_Y'

    def is_Curves(self):
        for row in arcpy.da.SearchCursor(self.layer,['SHAPE@']):
            geom = row[0]
            if geom:
                j    = json.loads(geom.JSON)
                if 'curve' in str(j):
                    return True
                else:
                    return False


    def fields(self):
        return [str(f.name) for f in arcpy.ListFields(self.layer)]

    def len(self):
        return int(str(arcpy.GetCount_management(self.layer)))

    def Get_Field_And_Type(self):
        return {str(f.name):str(f.type) for f in arcpy.ListFields(self.layer)}

    def vertxs_Count(self):
        if self.Geom_type != 'Point':
            return sum([row.objectid for row in arcpy.SearchCursor(self.layer) for part in row.shape for pt in part if pt])

    def  Select_By_Location(self,Connection,ref_layer,distance = 0,New_Layer = None,invert = ''):

        '''
        Connection:
        1) ARE_IDENTICAL_TO 
        2) BOUNDARY_TOUCHES
        3) ARE_IDENTICAL_TO
        4) INTERSECT
        5) HAVE_THEIR_CENTER_IN
        6) WITHIN
        7) WITHIN_A_DISTANCE
        8) COMPLETELY_WITHIN
        9) SHARE_A_LINE_SEGMENT_WITH
        DISTANCE:
        exp: "5 Meters"
        '''
        if invert != '':
            invert = "INVERT"
    
        FeatureLyr   = arcpy.MakeFeatureLayer_management(self.layer,self.name +'_Layer')
        arcpy.SelectLayerByLocation_management  (FeatureLyr,Connection,ref_layer,distance,'',invert)
        if not New_Layer:
            arcpy.DeleteFeatures_management (FeatureLyr)
        else:
            arcpy.Select_analysis           (FeatureLyr,New_Layer)

        if New_Layer:
            return New_Layer


    def Multi_to_single(self,temp_lyer = ''):
        
        temp_lyer = self.gdb + '\\' + 'Temp'
        save_name = self.layer
        arcpy.MultipartToSinglepart_management (self.layer,temp_lyer)
        arcpy.Delete_management                (self.layer)
        arcpy.Rename_management                (temp_lyer,save_name)

        return save_name


    def Erase(self,del_geom,Out_put = ''):

        if self.Geom_type == 'Polygon':
            if Out_put == '':
                Out_put = self.gdb + '\\' + self.name + 'Erased'
            Delete_polygons(self.layer,del_geom,Out_put)

            return Out_put

    def Destroy_layer(self):

        arcpy.Delete_management(self.layer)

    def Generate_Name(self,out_put = '' ,add_name= 'temp'):
        if arcpy.Exists(out_put):
            arcpy.Delete_management(out_put)
        return self.gdb + '\\' + self.name + '_'+ add_name if out_put == '' else out_put

    def Get_vertxs_As_Point(self,out_put = ''):

        if self.Geom_type != 'Point':
            out_put  = self.Generate_Name(out_put,'Point')
            arcpy.CopyFeatures_management([arcpy.PointGeometry(arcpy.Point(j.X,j.Y)) for i in arcpy.SearchCursor (self.layer) for n in i.shape for j in n if j],out_put)
            return out_put

    def Get_Label_Point_As_Point(self,out_put = ''):

        out_put  = self.Generate_Name(out_put,'Label_Point')
        arcpy.CopyFeatures_management([arcpy.PointGeometry(i.shape.labelPoint) for i in arcpy.SearchCursor (self.layer) if i.shape],out_put)
        return out_put

    def None_in_fields(self,fields_to_Check=[]):

        if not fields_to_Check:
            fields_to_Check = self.fields()

        fields_None = [i for i in arcpy.da.SearchCursor(self.layer,fields_to_Check)]
        New_list    = set([str(fields_to_Check[n]) for i in range(len(fields_None)) for n in range(len(fields_to_Check)) if fields_None[i][n] == None])
        
        return New_list

    def Get_Closest_Distance(self,layer2 = '',Return_list = False):

        if layer2 =='':
            layer2 = self.layer

        layer2_oid = str(arcpy.Describe(layer2).OIDFieldName)    
        list1      = [[i[0],i[1]] for i in arcpy.da.SearchCursor(self.layer,    ["SHAPE@",self.oid])]
        all_list   = [[row[1],round(n[0].distanceTo(row[0]),2),n[1]] for row in arcpy.da.SearchCursor(layer2,["SHAPE@",layer2_oid])for n in list1 if n[0].distanceTo(row[0]) != 0]

        df         = pd.DataFrame     (all_list,columns= ['KEY','NUM','KEY2'])
        gb         = df.groupby       ('KEY').agg({'NUM':'min'}).reset_index()
        df2        = pd.merge         (gb,df,how = 'inner',on=['NUM','KEY'])

        data_to_gis = {getattr(row, "KEY") : [getattr(row, "KEY2"),getattr(row, "NUM")] for row in df2.itertuples(index=True, name='Pandas')}
            
        add_field(self.layer,'ID_ref','LONG')
        add_field(self.layer,'Dis','DOUBLE')
        with arcpy.da.UpdateCursor(self.layer,[self.oid,'ID_ref','Dis']) as cursor:
            for row in cursor:
                if data_to_gis.get(row[0]):
                    row[1] = data_to_gis[row[0]][0]
                    row[2] = data_to_gis[row[0]][1]
                    cursor.updateRow (row)

        return self.layer if Return_list == False else df2.values.tolist()

    def delete_identical(self,field = ''):
        if field == '':
            field = self.Calc_XY()
        del_identical(self.layer,field)

    def Fill_Holes_in_Polygon(self,Out_put = '' , delete_layer = False ,Return_holes = False):

        if self.Geom_type == 'Polygon':
            if Out_put == '':
                Out_put = self.gdb + '\\' + self.name + '_filled'

        add_field                        (self.layer,'Holes','SHORT')
        arcpy.CalculateField_management  (self.layer, 'Holes', "0" , "PYTHON", "")

        Feature_to_polygon(self.layer,Out_put)
        New_polygons = int(str(arcpy.GetCount_management(Out_put))) -  self.len()
        print ('New polygons : {}'.format(str(New_polygons)))
        with arcpy.da.UpdateCursor(Out_put,['Holes']) as ucursor:
            for row in ucursor:
                if row[0] != 0:
                    row[0] = 1
                    ucursor.updateRow(row)


        if delete_layer:
            arcpy.Delete_management(Out_put)
            return New_polygons

        if Return_holes:
            arcpy.MakeFeatureLayer_management(Out_put,'Out_put_lyr',"\"Holes\" = 0")
            arcpy.DeleteFeatures_management  ('Out_put_lyr')
            return Out_put

        return Out_put


# # # # # # # # # # # # # # # Advanced Func # # # # # # # # # # # # # # # 

def stubborn_parts(path,bankal,tazar,Out_put,curves = ''):

    '''
    [INFO]- חוסם חורים ע"י מילוי של האזור הבעייתי, עובד במידה ואין חורים שמשותפים לכמה חלקות
    INPUT - 
    1) path     =  שכבת רצף
    2) bankal   -  שכבת ההבנקל הרחבה 
    3) tazar    - חלקות המודד
    4) curves   - במידה ויש קשתות, ניתן למחוק חלקות שנראות בדיוק כמו קשת

    OUTPUT - 
    1) Out_put  - שכבת הרצף לאחר התיקון, במידה ויש פגיעה, יחזור אחרוה
    '''

    print_arcpy_message("START Func: stubborn parts",1)

    memory = r'in_memory'
    gdb    = os.path.dirname(path)

    ## Create_layers
    inter            = memory + "\\" + "inter"
    sliver_curves    = memory + '\\'+ 'sliver_curves'
    path2            = memory + '\\'+ 'COPY_TEMP'
    Featur_to_poly   = memory + '\\' + 'Featur_to_poly'
    paracel_around   = memory +"\\"+ "paracel2_around"
    Parcel_deleted   = memory +"\\"+ "Parcel_deleted"
    parcal_all_Final = Out_put

    # מחיקה של קשתות שהן בדיוק כמו חלקות - יכול להיות בעיה אם החלקה היא קשת
    if curves:
        arcpy.MakeFeatureLayer_management(path,'path_lyr')
        arcpy.SelectLayerByLocation_management('path_lyr',"ARE_IDENTICAL_TO",curves)
        if int(str(arcpy.GetCount_management('path_lyr'))) > 0:
                        #print_arcpy_message("found identical rings")
                        arcpy.DeleteFeatures_management('path_lyr')             
            

    # מחיקה של חפיפות בשביל שהיו על תקן חורים
    arcpy.Intersect_analysis          ([path], inter)
    Delete_polygons                   (path,inter)

    ## יצירה של החורים לפיהם נמחק את הפוליגונים של הבנקל 
    Feature_to_polygon             (path,Featur_to_poly)
    Delete_polygons                (Featur_to_poly,path,sliver_curves)

    arcpy.CopyFeatures_management  (path,path2)

    arcpy.Select_analysis (path2,path2+"_path2_born") # DELETE
    # בדיקה שהחורים נוגעים התצ"ר והם רלוונטים
    arcpy.MakeFeatureLayer_management      (sliver_curves, 'sliver_curves_lyr')
    arcpy.SelectLayerByLocation_management ('sliver_curves_lyr', 'BOUNDARY_TOUCHES', tazar)

    num_slivers = int(str(arcpy.GetCount_management('sliver_curves_lyr')))
    if num_slivers > 0:

            print_arcpy_message ("you still have {} slivers, rebuild geometry and attributes".format(num_slivers),1)

            # #  #   #  #  #  #  יצירה של חלקות במעטפת נוספת מסביב לחלקה

            # מחיקה של החלקה מהבנקל שנוגעת בחור
            arcpy.MakeFeatureLayer_management     (path2, "FINAL2_lyr", "\"PARCEL_ID\" > 0") 
            arcpy.SelectLayerByLocation_management("FINAL2_lyr","SHARE_A_LINE_SEGMENT_WITH",sliver_curves)
            arcpy.Select_analysis                 ("FINAL2_lyr",Parcel_deleted) # שמירה על החלקה שאנחנו מוחקים
            arcpy.DeleteFeatures_management       ("FINAL2_lyr")

            # שליפת החלקות החדשות מסביב לחלקה שנעלמה
            arcpy.MakeFeatureLayer_management     (bankal, "paracel2_lyr") 
            arcpy.SelectLayerByLocation_management("paracel2_lyr","INTERSECT",Parcel_deleted)
            arcpy.Select_analysis                 ("paracel2_lyr", paracel_around)
            delete_parts_if_inside                (paracel_around,path2)          # מחיקת חלקות העבודה, השארת רק החלקות החדשות שנוספו מהבנקל 
            delete_parts_if_inside                (paracel_around,Parcel_deleted) # מחיקת החלקה עם החור, השארת רק החלקות החדשות שנוספו מהבנקל

            # חיבור השכבות החדשות לשכבת הרצף
            Update_Polygons                       (path2,paracel_around) 
            arcpy.Select_analysis (path2,path2 +'before_Feat_To_Poly')  # DELETEEEEEEEEEEEEEEEEE

            Feature_to_polygon                    (path2,Out_put)

            # חיבור המידע של השכבות לשכבה החדשה
            Spatial_Connection_To_LabelPoint      (Out_put,path)
            Spatial_Connection_To_LabelPoint      (Out_put,paracel_around)

            delete_parts_if_inside                (Out_put,paracel_around) # מחיקת שכבות הבנק"ל מסביב לאזור העבודה

            # בדיקה אם התיקון פגם בתוצאה, במידה וכן, יחזור אחורה
            before    = int(str(arcpy.GetCount_management(path)))
            after     = int(str(arcpy.GetCount_management(Out_put)))
            None_in_f = len(Layer_Management(Out_put).None_in_fields(['PARCEL']))

            if (before > after) or (None_in_f > 0):
                print_arcpy_message     ("Stubbern seems to delete features, Cancel and return 1 step back",1)
                arcpy.Select_analysis   (Out_put,Out_put + 'AFTER STUBBURN_PARTS') # #################delete
                arcpy.Delete_management (Out_put)
                arcpy.Select_analysis   (path,parcal_all_Final)

    else:
            print_arcpy_message          ("No stubborn parts")
            arcpy.CopyFeatures_management(path,Out_put)
      

def fix_holes_Overlaps_By_Length(path,tazar,path2):

    GDB = os.path.dirname(path)

    in_memory          = r'in_memory' 
    path2              = arcpy.CopyFeatures_management(path,path2)
    FEATURE_TO_POLYGON = in_memory + '\FEATURE_TO_POLYGON'
    slivers            = GDB + '\slivers'
    PARACELS_Only      = in_memory + '\PARACELS_Only'
    inter              = in_memory + '\inter'
    line               = GDB + '\Line'
    slivers_Intersect  = GDB + '\slivers_Intersect'

    arcpy.Intersect_analysis          ([path2], inter)
    if int(str(arcpy.GetCount_management(inter))) > 0:
        Delete_polygons                   (path2,inter)

    arcpy.AddField_management        (path2, "KEY_parcel", "LONG")
    arcpy.CalculateField_management  (path2, "KEY_parcel", "!OBJECTID!", "PYTHON", "")
        
    Feature_to_polygon(path2, FEATURE_TO_POLYGON)
    Delete_polygons             (FEATURE_TO_POLYGON, path2, slivers)
        
    number_of_slivers = int(str(arcpy.GetCount_management(slivers)))
    if number_of_slivers > 0:
            print("there is {} holes, start working to fix them".format(str(number_of_slivers)))

            arcpy.AddField_management        (slivers, "KEY_sliv", "LONG")
            arcpy.CalculateField_management  (slivers, "KEY_sliv", "!OBJECTID!", "PYTHON", "")

            Delete_polygons             (path2, tazar, PARACELS_Only)

            Polygon_To_Line   (PARACELS_Only, line)

            sliver_feature_layer = GDB + '\\' + 'sliver_feature_layer'
            arcpy.MakeFeatureLayer_management      (slivers, sliver_feature_layer)
            arcpy.SelectLayerByLocation_management (sliver_feature_layer, 'BOUNDARY_TOUCHES', tazar)
            intersect_list = [sliver_feature_layer,line]

            arcpy.Intersect_analysis    (intersect_list, slivers_Intersect, "ALL", ".001 Meters", "INPUT")

            data       = [[row[0],row[1],row[2]] for row in arcpy.da.SearchCursor(slivers_Intersect,['KEY_sliv','FID_Line','SHAPE@LENGTH'])]
                
            df         = pd.DataFrame(data,columns= ['KEY_sliv','KEY_parcel_1','SHAPE@LENGTH'])
            df["RANK"] = df.groupby('KEY_sliv')['SHAPE@LENGTH'].rank(method='first',ascending=False)
            df         = df[df['RANK'] == 1]

            data_to_gis = [[getattr(row, "KEY_sliv"), getattr(row, "KEY_parcel_1")]for row in df.itertuples(index=True, name='Pandas')]

            arcpy.AddField_management (slivers, "ID_KEY_par", "LONG")
            for data in data_to_gis:
                with arcpy.da.UpdateCursor(slivers,['KEY_sliv','ID_KEY_par']) as cursor:
                    for row in cursor:
                        if row[0] == data[0]:
                            row[1] = data[1]
                            cursor.updateRow (row)
                            

            x = [[x[0],x[1]] for x in arcpy.da.SearchCursor(slivers,['ID_KEY_par','SHAPE@'])]
            for i in x:
                with arcpy.da.UpdateCursor(path2,['OID@','SHAPE@']) as icursor:
                    for row in icursor:
                            if row[0] == i[0]:
                                    new = row[1].union(i[1])
                                    row[1] = new
                                    icursor.updateRow(row)

                                
            arcpy.Delete_management(line)
            arcpy.Delete_management(slivers_Intersect)
    else:
            print("no holes found".format(str(number_of_slivers)))


def Snap_border_pnts(ws,border,parcel_all,Dis_search = 1):


    print_arcpy_message('START Func: Snap border pnts',1)

    tazar_border = 'in_memory\Tazar_Border_diss'
    arcpy.Dissolve_management(border,tazar_border)

    arcpy.MakeFeatureLayer_management(parcel_all, "parcel_all_lyr")
    arcpy.SelectLayerByLocation_management("parcel_all_lyr", "WITHIN_A_DISTANCE", tazar_border, '5 Meters')

    conn = sqlite3.connect(':memory:')
    c = conn.cursor()

    c.execute('''CREATE TABLE border(pnt_num real, x real, y real, xy text, part real, oid real, junction real, linearity real)''')


    c.execute('''CREATE TABLE parcels(pnt_num real, x real, y real, xy text, part real, oid real, junction real, linearity real)''')


    VerticesToTable2(tazar_border, "border",c)
    border_vertices = [row for row in c.execute('''SELECT * FROM border''')]
    VerticesToTable2("parcel_all_lyr", "parcels",c)

    parcel_non_common_vertices = [row for row in c.execute('''SELECT * FROM parcels
                                                                     left join border
                                                                     on parcels.xy = border.xy
                                                                     where  border.xy is null''')]
    
    
    border_geom = arcpy.CopyFeatures_management(tazar_border, arcpy.Geometry())[0]

    vertices_on_border_outline = [row for row in parcel_non_common_vertices if border_geom.distanceTo (arcpy.Point(row[1], row[2])) < 5]

    distance_vertices = [[p[:8] + b[:8],math.sqrt(((p[1]-b[1])**2)+((p[2]-b[2])**2))] for p in vertices_on_border_outline for b in border_vertices if math.sqrt(((p[1]-b[1])**2)+((p[2]-b[2])**2))\
        < Dis_search or (float("{0:.2f}".format(p[1])) == float("{0:.2f}".format(b[1])) and float("{0:.2f}".format(p[2])) == float("{0:.2f}".format(b[2])))]


    rows = arcpy.UpdateCursor("parcel_all_lyr")
    for row in rows:
        geometry = row.Shape
        oid = row.OBJECTID
        pts = []
        poly_vertices = [r for r in distance_vertices if r[0][5] == oid]
        for part in geometry:
            for pt in part:
                if str(type(pt)) != "<type 'NoneType'>":
                    num_point = 0
                    #print str(pt.X) + "--" + str(pt.Y)
                    this_x = float("{0:.2f}".format(pt.X))
                    this_y = float("{0:.2f}".format(pt.Y))      
                    this_vertex = [p for p in poly_vertices if float("{0:.2f}".format(p[0][1])) == this_x and float("{0:.2f}".format(p[0][2])) == this_y]
                    if this_vertex:
                        if this_vertex[0][0][8] == None:
                            if this_vertex[0][0][7] < 0.5 and this_vertex[0][0][6] == 1:
                                print ("pseodo: delete vertex")
                            else:
                                #print "pseodo, but important: keep the vertex"
                                point = pt
                                pts.append(point)
                        # tazar point in buffer
                        else:
                            # check minimum distance
                            the_minimum_vertex = [v for v in this_vertex if v[1] == min([i[1] for i in this_vertex])]
                            point = arcpy.Point(the_minimum_vertex[0][0][9], the_minimum_vertex[0][0][10])
                            pts.append(point)
                    # point not on sliver: keep the vertex
                    else:
                        point = pt
                        pts.append(point)
                    if num_point == 0:
                        first_point = point
                    num_point = num_point + 1
        polygon = PtsToPolygon(pts)
        if pts[0] != pts[-1] and first_point:
            #print "ooops.... - polygon not closed"
            pts.append(first_point)
        row.Shape       = polygon
        rows.updateRow(row)


def clean_pseudo(parcel_all, border,curves):

    def clean_pseudo_vertices(polygon_before, border_geom, nodes_pts):
        for part in polygon_before:
            pts_final = []
            pts_trio = []
            deleted = []
            for pt in part:
                    if str(type(pt)) != "<type 'NoneType'>":
                            pts_trio.append([pt.X, pt.Y])
                            if len(pts_trio) == 3:
                                x = pts_trio[1][0]
                                y = pts_trio[1][1]
                                if collinearity(pts_trio[0], pts_trio[1], pts_trio[2]) < 0.9 and [float("{0:.2f}".format(x)), float("{0:.2f}".format(y))] not in nodes_pts and border_geom_buffer.contains(arcpy.Point(x,y)):
                                        #print_arcpy_message("delete vertex",1)
                                        deleted.append([x,y])
                                        #print_arcpy_message([x,y],1)
                                else:
                                    pts_final.append(arcpy.Point(x,y))
                                pts_trio = [pts_trio[1], pts_trio[2]]
                            else:
                                pts_final.append(pt)
            if deleted:
                polygon_after = PtsToPolygon(pts_final)
                return polygon_after
            else:
                return polygon_before


    print_arcpy_message("START Func: clean pseudo",1)

    before_vertxs = Layer_Management(parcel_all).vertxs_Count()
    
    arcpy.MakeFeatureLayer_management(parcel_all, "parcel_all_lyr")
    arcpy.SelectLayerByLocation_management("parcel_all_lyr", "SHARE_A_LINE_SEGMENT_WITH", border)
        
    node = "in_memory" +'\\'+"node"
    arcpy.CreateFeatureclass_management("in_memory", "node", "POINT", "", "", "",border)
    border_geom = arcpy.CopyFeatures_management(border, arcpy.Geometry())
    border_geom_buffer =  border_geom[0].buffer(0.05)

    Layer_Management(border).Get_vertxs_As_Point(node)     

    # להכניס גם נקודות של הבנקל לרשימת נקודות שלא צריכות להימחק
    # למה לא להשתמש בנקודות קיימות?            

    nodes_pts = [[float("{0:.2f}".format(row.Shape.centroid.X)), float("{0:.2f}".format(row.Shape.centroid.Y))] for row in arcpy.SearchCursor(node)]
    upd_rows = arcpy.UpdateCursor("parcel_all_lyr")
    for upd_row in upd_rows:
        polygon_before = upd_row.Shape
        polygon_after = clean_pseudo_vertices(polygon_before, border_geom, nodes_pts)
        if polygon_after:
            upd_row.Shape = polygon_after
            upd_rows.updateRow(upd_row)


    Fix_curves              (parcel_all,border,curves)

    after_vertxs    = Layer_Management(parcel_all).vertxs_Count()
    deleted_vertexs = before_vertxs - after_vertxs
    print_arcpy_message('Total Vertexs Deleted: {}'.format(deleted_vertexs))



def Move_Vertices_By_Name(polygon,points,field_name_points,points_to_move,field_name_to_move = 'POINT_NAME',Dis_limit_to_move = 10):

    '''
    [INFO] -  מזיז וורטקסים של פוליגון לפי השם של של הוורטקסים שלו (שכבת נקודות נפרדת), לעומת שם של שכבת נקודות אחרת
    INPUT - 
    1) polygon            - שכבת הפוליגון שתזוז
    2) points             - שכבת נקודות על הוורטקסים של השכבה הפוליגונלית, עם שמות
    3) field_name_points  - שם השדה של שכבת הנקודות בו מופיע שם הנקודה
    4) points_to_move     - שכבת הנקודות אליהם אנחנו רוצים להזיז
    5) field_name_to_move - שם השדה בו נמצת שם הנקודה אליה אנו רוצים להזיז

    OUTPUT - שכבת הפוליגונים תזוז אל הנקודות לה יש שם דומה
    '''

    Save_Source  = polygon + 'Save'
    before_inter = r'in_memory\\before_inter'
    after_inter  = r'in_memory\\after_inter'
    save_name    = polygon

    arcpy.Select_analysis(polygon,Save_Source)

    data_points = {str(round(n.X,2))    +'-' + str(round(n.Y,2)):str(i.getValue(field_name_points)) for i in arcpy.SearchCursor(points) for n in i.shape if i.getValue(field_name_points) != None and i.getValue(field_name_points) != ''}
    move_points = {str(i.getValue(field_name_to_move)):[n.X,n.Y] for i in arcpy.SearchCursor(points_to_move) for n in i.shape if i.getValue(field_name_to_move) != None and i.getValue(field_name_to_move) != ''}

    if (bool(data_points) == True) and (bool(move_points)==True):
        data_polygons = [[str(round(pts.X,2)) +'-' + str(round(pts.Y,2)),round(pts.X,2),round(pts.Y,2),'',0,0] for i in arcpy.SearchCursor(polygon) for n in i.shape for part in i.shape for pts in part if pts]

        for i in data_polygons:
            if data_points.get(i[0]):
                i[-3] = data_points[i[0]]
                
        for i in data_polygons:
            if move_points.get(i[-3]):
                i[-2] = move_points[i[-3]][0]
                i[-1] = move_points[i[-3]][1]

        data_polygons = {i[0]:i[1:] for i in data_polygons if i[-3] != '' and i[-2] != 0}

        num_inter = int(str(arcpy.GetCount_management(arcpy.Intersect_analysis(polygon,before_inter))))

        if data_polygons:
            with arcpy.da.UpdateCursor(polygon,["SHAPE@"]) as cursor:
                for row in cursor:
                    geom = row[0]
                    array = arcpy.Array()
                    count = 0
                    for part in geom:
                        for pt in part:
                            if pt:
                                key = str(round(pt.X,2))+'-' + str(round(pt.Y,2))
                                if data_polygons.get(key):
                                    distance = dis(data_polygons[key][-2],data_polygons[key][-5],data_polygons[key][-1],data_polygons[key][-4])
                                    if distance < Dis_limit_to_move:
                                        Point =  arcpy.Point(data_polygons[key][-2],data_polygons[key][-1])
                                    else:
                                        Point = arcpy.Point(pt.X,pt.Y)
                                else:
                                    Point = arcpy.Point(pt.X,pt.Y)
                                if count == 0:
                                    first = Point
                                array.add(Point)
                                count += 1
                            else:
                                array.add(first)
                                array.add(None)
                    New_poly = arcpy.Polygon(array)
                    row[0] = New_poly
                    cursor.updateRow(row)
            del cursor 

        num_inter_AFTER = int(str(arcpy.GetCount_management(arcpy.Intersect_analysis(polygon,after_inter))))

        if num_inter_AFTER > num_inter:
            print_arcpy_message     ("Coudnt Move By Name, Polygon Building Error",2)
            arcpy.Delete_management (polygon)
            arcpy.Rename_management (Save_Source,save_name)

def Update_Polygons(Layer_To_Update,New_Item,To_New_Layer = ''):

    Delete_polygons                            (Layer_To_Update ,New_Item,To_New_Layer)
    if To_New_Layer == '':
        To_New_Layer = Layer_To_Update
    arcpy.Append_management                    (New_Item        ,To_New_Layer, 'NO_TEST')


def Connect_Lines(layer,layer_new,min_dis):

    '''
    [INFO] - יצירת קווים המחברים בין יישויות במידה ויש ישות שלא מחוברת "עם קצה" יש לבחור במרחק מקסימלי לחיבור
    '''

    new_list = Layer_To_Edge_list(layer)

    Diss = 'in_memory\Diss_layer'

    ws, fc_name = os.path.split (layer_new)
    s_r         = arcpy.Describe (layer).spatialReference

    if arcpy.Exists(layer_new):
        arcpy.Delete_management(layer_new)

    line = arcpy.CreateFeatureclass_management (ws, fc_name, 'POLYLINE', spatial_reference=s_r)

    insert = arcpy.InsertCursor(line)

    for i in range(len(new_list)):
        if new_list[i][2] < min_dis:
            points   = [arcpy.Point(new_list[i][0][1],new_list[i][0][2]),arcpy.Point(new_list[i][1][1],new_list[i][1][2])]
            array    = arcpy.Array(points)
            polyline = arcpy.Polyline(array)
            feat     = insert.newRow ()
            feat.shape    = polyline
            insert.insertRow(feat)

    arcpy.RepairGeometry_management(layer_new)

    Delete_By_length                (layer_new,0.2)
    return layer_new


def Clean_non_exist_pnts(AOI,border,bankal,tazar_copy):

    print_arcpy_message("START FUNC: Clean_non_exist_pnts")

    gdb            = os.path.dirname(border)
    bankal_cut     = gdb + '\\' + 'Bankal_Cut_inter'
    holes_to_keep  = gdb + '\\' + 'holes_to_keep'
    holes_keeping  = gdb + '\\' + 'holes_keeping'
    pts_on_border  = gdb + '\\' + 'pnts_on_border'

    save_source    = gdb + '\\' + 'Save_Source'
    save_name      = AOI
    arcpy.Select_analysis(AOI,save_source)

    Layer_Management(bankal).Select_By_Location('INTERSECT',border,'10 Meters',bankal_cut)
        
    border_xy     = [[round(j.X,1),round(j.Y,1)] for i in arcpy.SearchCursor (border) for n in i.shape for j in n if j]     # לא נוגעים בנקודות המודד
    bankal_cut_xy = [[round(j.X,1),round(j.Y,1)] for i in arcpy.SearchCursor (bankal_cut) for n in i.shape for j in n if j] # לא נוגעים בנקודות הבנקל
    Layer_Management(Layer_Management(AOI).Get_vertxs_As_Point()).Select_By_Location('INTERSECT',border,'0.01 Meters',pts_on_border)
    pts_border_xy = [[round(row.Shape.centroid.X,1),round(row.Shape.centroid.Y,1)] for row in arcpy.SearchCursor(pts_on_border)] # לא נוגעים בנקודות החדשות על התצר
    arcpy.Dissolve_management(AOI,r'in_memory\dissolve_me')
    AOI_border_xy = [[round(j.X,1),round(j.Y,1)] for i in arcpy.SearchCursor (r'in_memory\dissolve_me') for n in i.shape for j in n if j] # לא נוגעים בנקודות גבול של אזור העבודה
  
    saved_vertexs = border_xy + bankal_cut_xy + pts_border_xy + AOI_border_xy

    Feature_to_polygon (AOI,holes_to_keep)
    Delete_polygons    (holes_to_keep,AOI)
    holes_befotre = int(str(arcpy.GetCount_management(holes_to_keep)))

    Layer_Management        (holes_to_keep).Select_By_Location('INTERSECT',border,0,holes_keeping,'invert')

    Delete_polygons(AOI,border)

    Ucursor = arcpy.UpdateCursor(AOI)
    for row in Ucursor:
        geom = row.shape
        array = arcpy.Array()
        for part in geom:
            count = 0
            for pt in part:
                if pt:
                    if [round(pt.X,1),round(pt.Y,1)] in saved_vertexs:
                        array.add(pt)
                    if count == 0:
                        first = pt
                    else:
                        pass
                    count += 1
                else:
                    array.add(first)
                    array.add(None)
                    count = 0

        New_poly = arcpy.Polygon(array)
        row.shape = New_poly
        Ucursor.updateRow(row)

    del Ucursor 

    Delete_polygons         (AOI,holes_keeping)
    arcpy.Append_management (tazar_copy,AOI,'NO_TEST')

    Feature_to_polygon (AOI,holes_to_keep)
    Delete_polygons    (holes_to_keep,AOI)
    after_befotre = int(str(arcpy.GetCount_management(holes_to_keep)))

    if holes_befotre < after_befotre:
        print_arcpy_message("There is modad  point that may be missing, return 1 step back")
        arcpy.Delete_management(AOI)
        arcpy.Rename_management(save_source,save_name)



def print_arcpy_message(msg,status = 1):
    '''
    return a message :
    
    print_arcpy_message('sample ... text',status = 1)
    [info][08:59] sample...text
    '''
    msg = str(msg)
    
    if status == 1:
        prefix = '[info]'
        msg = prefix + str(datetime.datetime.now()) +"  "+ msg
        # print (msg)
        arcpy.AddMessage(msg)
        
    if status == 2 :
        prefix = '[!warning!]'
        msg = prefix + str(datetime.datetime.now()) +"  "+ msg
        print (msg)
        arcpy.AddWarning(msg)
            
    if status == 0 :
        prefix = '[!!!err!!!]'
        
        msg = prefix + str(datetime.datetime.now()) +"  "+ msg
        print (msg)
        arcpy.AddWarning(msg)
        msg = prefix + str(datetime.datetime.now()) +"  "+ msg
        print (msg)
        arcpy.AddWarning(msg)
            
        warning = arcpy.GetMessages(1)
        error   = arcpy.GetMessages(2)
        arcpy.AddWarning(warning)
        arcpy.AddWarning(error)
            
    if status == 3 :
        prefix = '[!FINISH!]'
        msg = prefix + str(datetime.datetime.now()) + " " + msg
        print (msg)
        arcpy.AddWarning(msg) 


def del_line_Not_on_parcels(ARC_bankal,Parcel_makor):

    #  # cuting layer , to work on less data # #

    # # Check Arc points\ID

    dicLine = [[str(round(pt.X,1)) + '-' + str(round(pt.Y,1)),row.objectid] for row in arcpy.SearchCursor(ARC_bankal) for part in row.shape for pt in part if pt]
    data_p  = [str(round(pts.X,1)) +'-' + str(round(pts.Y,1)) for i in arcpy.SearchCursor(Parcel_makor) for n in i.shape for part in i.shape for pts in part if pts]
    del_line = list(set([i[1] for i in dicLine if i[0] not in data_p]))
    
    if del_line:
        arcpy.MakeFeatureLayer_management      (ARC_bankal,'ARC_bankal_lyr')
        arcpy.SelectLayerByAttribute_management('ARC_bankal_lyr',"NEW_SELECTION","\"OBJECTID\" IN ("+str(del_line)[1:-1]+")")
        arcpy.DeleteFeatures_management        ('ARC_bankal_lyr')

  


def Delete_polygons(fc,del_layer,Out_put = ''):

    '''
    fc        = השכבה הראשית- שכבה ממנה רוצים למחוק
    del_layer = שכבה שתמחק את השכבה הראשית
    Out_put   = שכבת הפלט, במידה ולא תוכנס שכבה, ימחק מהשכבה הראשית
    '''
    
    desc = arcpy.Describe(fc)

    if not Out_put == '':
        fc = arcpy.CopyFeatures_management(fc,Out_put)
    else:
        Out_put = fc
    
    if desc.ShapeType == u'Point':
        del_layer_temp = 'in_memory' + '\\' + 'Temp'
        arcpy.Dissolve_management(del_layer,del_layer_temp)

        if desc.ShapeType == u'Point':
            geom_del = [row.shape for row in arcpy.SearchCursor (del_layer_temp)][0]
            Ucursor  = arcpy.UpdateCursor (Out_put)
            for row in Ucursor:
                point_shape = row.shape.centroid
                if geom_del.distanceTo(point_shape)== 0:
                    Ucursor.deleteRow(row)
                else:
                    pass
            del Ucursor
        del del_layer_temp
                        
    else:
        count_me = int(str(arcpy.GetCount_management(del_layer)))
        if count_me > 0:
            temp = 'in_memory' +'\\'+'_temp'
            arcpy.Dissolve_management(del_layer,temp)
            if int(str(arcpy.GetCount_management(temp))) > 0:
                geom_del = [row.shape for row in arcpy.SearchCursor (temp)][0]
                Ucursor  = arcpy.UpdateCursor (Out_put)
                for row in Ucursor:
                    geom_up     = row.shape
                    new_geom    = geom_up.difference(geom_del)
                    try:
                        row.shape = new_geom
                        Ucursor.updateRow (row)
                    except:
                        pass
                del Ucursor
            arcpy.Delete_management(temp)
        else:
            pass

                        
    if desc.ShapeType == u'Point':
        pass
    else:
        up_cursor = arcpy.UpdateCursor(Out_put)
        for row in up_cursor:
            geom = row.shape
            if geom.area == 0:
                up_cursor.deleteRow(row)
        del up_cursor
        
    arcpy.RepairGeometry_management(Out_put)
    return Out_put


def topology(final):
        
    gdb                 = os.path.dirname(final)
    Dissolve_temp       = r'in_memory' + '\\'+ 'dissolve_me'
    Feature_to_poly     = r'in_memory' + '\\'+'Feature_to_poly'
    Topolgy_Check_holes = gdb + '\\'+'Topolgy_Check_holes'
    Topolgy__intersect  = gdb + '\\'+'Topolgy_Check_intersect'

    arcpy.Dissolve_management                 (final,Dissolve_temp)
    Feature_to_polygon                        (Dissolve_temp,Feature_to_poly)
    Delete_polygons                           (Feature_to_poly,Dissolve_temp,Topolgy_Check_holes)
    count = int(str(arcpy.GetCount_management (Topolgy_Check_holes)))
    if count == 0:
        arcpy.Delete_management(Topolgy_Check_holes)
        Topolgy_Check_holes = None

    over_lap       = arcpy.Intersect_analysis([final],Topolgy__intersect)
    over_lap_count = int(str(arcpy.GetCount_management (over_lap)))
    if over_lap_count == 0:
        arcpy.Delete_management(Topolgy__intersect)
        Topolgy__intersect = None

    arcpy.Delete_management(Dissolve_temp)
    del Dissolve_temp
    return Topolgy_Check_holes,Topolgy__intersect




def Fix_Pnt_Tolerance(AOI_final,AOI_Point):

    pnt_save = {str(round(pt.X,1)) + '-' + str(round(pt.Y,1)):[pt.X,pt.Y] for row in arcpy.SearchCursor(AOI_final) for part in row.shape for pt in part if pt}

    with arcpy.da.UpdateCursor(AOI_Point,['SHAPE@']) as cursor:
        for row in cursor:
            X_pt = str(round(row[0].centroid.X,1))
            Y_pt = str(round(row[0].centroid.Y,1))
            key = X_pt + '-' + Y_pt
            if pnt_save.get(key):
                point = arcpy.Point(pnt_save[key][0],pnt_save[key][1])
                row[0] = point
                cursor.updateRow(row)


def getLayerPath(fc,CURRENT = 'CURRENT'):

    '''
    [INFO]
        מקבל שכבה מהמשתמש, מחזיר את בסיס הנתונים בו הוא נמצא
    '''

    MXD  = arcpy.mp.ArcGISProject (CURRENT)
    df   = MXD.listMaps('Layers')[0]
    lyrs = df.listLayers()
    if lyrs:
        for i in lyrs:
            if i.isFeatureLayer:
                if 'demo.gdb' not in os.path.dirname(i.dataSource) and 'Scratch' not in os.path.dirname(i.dataSource):
                    return os.path.dirname(i.dataSource)
    else:
        print_arcpy_message(os.path.dirname(fc))
        return os.path.dirname(fc)


def CreateWorkingGDB(gdb,Folder,copy,fc_name,CURRENT):


    '''
        [INFO]
            Create GDB and copy all the layers from the source GDB to the new GDB, also except layer from mxd
        INPUT - GDB     which all his layer will be copyed
              - Folder  New GDB will be created, with copyed layers
              - Names of the layer that you need to copy, in list
              - fc_name , name of 1 layer from the mxd, in case layer from MXD will enter the tool

        RETURN - NEW gdb with all the "workimg on" layers
    '''

    def get_fc_from_mxd(fc_name,CURRENT):
    #CURRENT
        MXD  = arcpy.mp.ArcGISProject(CURRENT)
        df   = MXD.listMaps('Layers')[0]
        lyrs = lyrs = df.listLayers()
        fc = None
        for lyr in lyrs:
                if lyr.isFeatureLayer:
                        if os.path.basename(lyr.dataSource) == fc_name:
                                fc = lyr.dataSource
        return fc
    
    folder_source = os.path.dirname  (gdb)
    name          = os.path.basename (folder_source)
    tazar_num     = ''.join([i for i in name if i.isdigit()])
    ws = Folder + '\\' + 'Tazar_{}.gdb'.format(tazar_num)
    if arcpy.Exists(ws):
                    arcpy.Delete_management(ws)

    print ('Tazar_{}.gdb'.format(tazar_num))
    arcpy.CreateFileGDB_management(Folder,'Tazar_{}.gdb'.format(tazar_num))
    
    return_list = []
    for fc in copy:
        try:
            arcpy.CopyFeatures_management    (gdb + '\\' + fc ,ws + '\\' + fc + '_copy')
        except:
            copy_me = get_fc_from_mxd        (fc_name,CURRENT)
            gdb     = os.path.dirname(copy_me)
            arcpy.CopyFeatures_management    (gdb + '\\' + fc ,ws + '\\' + fc + '_copy')
        return_list.append(ws + '\\' + fc + '_copy')

    return ws


def generateCurves(fc):
    desc    = arcpy.Describe(fc)
    fc_name = desc.name
    fc_gdb  = desc.path
    Curves  = fc_gdb + "\\" + fc_name + "_curves_polygon"
    #print "generateCurves("+fc_name+")..."
    arcpy.CreateFeatureclass_management(fc_gdb, fc_name + "_curves_polygon", "POLYGON", "", "", "",fc)
    curveFeatureList = []
    for row in arcpy.SearchCursor(fc):
        pts = []
        geom = row.Shape
        j = json.loads(geom.JSON)
        if 'curve' in str(j):
            coords = geom.__geo_interface__['coordinates']

            for i in coords:
                if i:
                    for f in i:
                        if f:
                            pts.append(arcpy.Point(f[0],f[1]))
                        else:
                            pts.append(None)

        poly    = Split_List_by_value(pts,None,True) 

        if pts:
            array        = arcpy.Array(None)
            polygon = arcpy.Polygon(array, arcpy.SpatialReference("Israel TM Grid"))
            if len(poly) > 1:
                for part in poly:
                    poly_part = PtsToPolygon(part)
                    polygon   = polygon.symmetricDifference(poly_part)
            else:
                polygon = PtsToPolygon(poly[0])

            diff    = polygon.symmetricDifference(geom)
            diff_sp = arcpy.MultipartToSinglepart_management(diff, arcpy.Geometry())
            if len(diff_sp) > 0:
                arcpy.Append_management(diff_sp, Curves, "NO_TEST")
    return Curves


def dis(x1,x2,y1,y2):
    dist = math.sqrt(((x1-x2)**2) + ((y1-y2)**2))
    return dist

def PtsToPolygon(pts):
    point = arcpy.Point()
    array = arcpy.Array()
    for point in pts:
        array.add(point)
    array.add(array.getObject(0))

    polygon = arcpy.Polygon(array, arcpy.SpatialReference("Israel TM Grid"))
    return polygon


def Delete_layers_after_use(layers):
    for i in layers:
        try:
                arcpy.Delete_management (i)
        except:
            pass 


def collinearity(p1, p2, p3):
    """return True if 3 points are collinear.
    tolerance value will decide whether lines are collinear; may need
    to adjust it based on the XY tolerance value used for feature class"""
    x1, y1 = p1[0], p1[1]
    x2, y2 = p2[0], p2[1]
    x3, y3 = p3[0], p3[1]
    res = x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)
    return abs(res)  



def Polygon_To_Line_holes(Polygon,New_Line):

    Multi_to_single(Polygon)

    ws, fc_name = os.path.split (New_Line)
    s_r = arcpy.Describe (Polygon).spatialReference
    arcpy.CreateFeatureclass_management (ws, fc_name, 'POLYLINE', spatial_reference=s_r)

    ins_cursor = arcpy.da.InsertCursor (New_Line, ["SHAPE@"])

    New_Lines = []
    with arcpy.da.SearchCursor(Polygon,['SHAPE@','OBJECTID']) as cursor:
        for row in cursor:
            geom = row[0]
            array_temp = arcpy.Array() 
            conti      = True
            for part in geom:
                for pt in part:
                    if pt:
                        if conti:
                            array_temp.append(pt)
                        else:
                            New_Lines.append(arcpy.Point(pt.X,pt.Y))
                    else:
                        New_Lines.append(None)
                        conti = False



            polyline = arcpy.Polyline (array_temp,s_r)
            ins_cursor.insertRow ([polyline])

    del cursor

    # Insert rings in polygon, and make them lines

    InsertCursor = arcpy.InsertCursor(New_Line)
    insert       = InsertCursor.newRow()

    New_Lines = Split_List_by_value(New_Lines,None,True)

    for row in New_Lines:
        if row:
            row            = arcpy.Array(row)
            line           = arcpy.Polyline(row)
            insert.shape   = line
            InsertCursor.insertRow  (insert)

    arcpy.RepairGeometry_management(New_Line)


def Polygon_To_Line(fc,layer_new):
    
    ws, fc_name = os.path.split (layer_new)
    s_r = arcpy.Describe (fc).spatialReference

    if arcpy.Exists(layer_new):
        arcpy.Delete_management(layer_new)
        
    line = arcpy.CreateFeatureclass_management (ws, fc_name, 'POLYLINE', spatial_reference=s_r)

    insert = arcpy.da.InsertCursor(line,"SHAPE@")

    Search = arcpy.da.SearchCursor(fc,"SHAPE@"  )
    Get_Line_list = []
    pid = 0
    for row in Search:
        for part in row[0]:
            for pt in part:
                if pt:
                    Get_Line_list.append([pid,pt.X,pt.Y])
                else:
                    pass
        pid +=1

    for i in range(pid):
        points   = [arcpy.Point(n[1],n[2]) for n in Get_Line_list if n[0] == i]
        array    = arcpy.Array(points)
        polyline = arcpy.Polyline(array)
        insert.insertRow([polyline])

    arcpy.RepairGeometry_management(layer_new)

def Fix_curves(fc,tazar_border,curves):

        print_arcpy_message("START Func: Fix curves",1)

        name       = fc
        gdb        = os.path.dirname(fc)
        curves_cut = gdb + '\\' + 'curves_cut'
        fc2        = gdb + '\\' + 'temp'
        
        Delete_polygons(curves,tazar_border,curves_cut)
        Delete_polygons(fc,curves_cut,fc2)      
        arcpy.MakeFeatureLayer_management(fc2,'ARCEL_ALL_FINAL_lyr')
        
        list_Upd = []
        cursor = arcpy.SearchCursor(curves_cut)
        for i in cursor:
                arcpy.SelectLayerByLocation_management('ARCEL_ALL_FINAL_lyr',"SHARE_A_LINE_SEGMENT_WITH",i.shape)
                layer_ID = [row.OBJECTID for row in arcpy.SearchCursor('ARCEL_ALL_FINAL_lyr',['OBJECTID','PARCEL_ID']) if row.PARCEL_ID is not None]
                if layer_ID:
                    list_Upd.append([layer_ID[0],i.shape])


        for i in list_Upd:
                upd_cursor = arcpy.UpdateCursor(fc2)
                for up_row in upd_cursor:
                        geom = up_row.shape
                        id   = up_row.OBJECTID
                        if str(id) == str(i[0]):
                                new_geom     = geom.union (i[1])
                                up_row.shape = new_geom
                                upd_cursor.updateRow(up_row)  


        arcpy.Delete_management (fc)
        arcpy.Delete_management (curves_cut)
        arcpy.Rename_management (fc2, name)
        

def VerticesToTable2(fc, table,c):
    #read fc vertices
    xys = []
    fc_rows = [[row.OBJECTID, row.SHAPE] for row in arcpy.SearchCursor(fc)]
    for item in fc_rows:
        fc_oid = item[0]
        geometry = item[1]
        fc_part_num = 1
        for fc_part in geometry:
            fc_pnt_num = 1
            xys_for_linerity = []
            for fc_pnt in fc_part:
                if fc_pnt:
                    x = fc_pnt.X
                    y = fc_pnt.Y
                    c.execute("INSERT INTO "+table+" VALUES ("+str(fc_pnt_num) + ","
                              +str(x) + ","
                              +str(y) + ","
                              +"'"+str(float("{0:.2f}".format(x)))+"-"+str(float("{0:.2f}".format(y)))+"',"
                              +str(fc_part_num) + ","
                              +str(fc_oid) +", 0, 0)")
                    xys.append(str(float("{0:.2f}".format(x)))+"-"+str(float("{0:.2f}".format(y))))
                    xys_for_linerity.append([x,y])
                    fc_pnt_num = fc_pnt_num + 1
            counter = 0
            for xy_l in xys_for_linerity:
                if counter > 0 and counter < len(xys_for_linerity) - 1:
                    xy1 = xys_for_linerity[counter  - 1]
                    xy2 = xy_l
                    xy3 = xys_for_linerity[counter  + 1]
                    linearity = collinearity(xy1, xy2, xy3)
                    c.execute("UPDATE "+table+" SET linearity = "+str(linearity) + " WHERE x = "+str(xy_l[0])+" AND y = "+str(xy_l[1]))
                if counter == 0:
                    xy1 = xys_for_linerity[-1]
                    xy2 = xy_l
                    xy3 = xys_for_linerity[1]
                    linearity = collinearity(xy1, xy2, xy3)
                    c.execute("UPDATE "+table+" SET linearity = "+str(linearity) + " WHERE x = "+str(xy_l[0])+" AND y = "+str(xy_l[1]))
                if counter == len(xys_for_linerity) - 1:
                    xy1 = xys_for_linerity[-2]
                    xy2 = xy_l
                    xy3 = xys_for_linerity[0]
                    linearity = collinearity(xy1, xy2, xy3)
                    c.execute("UPDATE "+table+" SET linearity = "+str(linearity) + " WHERE x = "+str(xy_l[0])+" AND y = "+str(xy_l[1]))
                counter = counter + 1
                              
    for xy in xys:
        xy_count = len([i for i in xys if i == xy])
        c.execute("UPDATE "+table+" SET junction = "+str(xy_count) + " WHERE xy = '"+xy+"'")




def fix_tolerance(layer_path,border):
    
    gdb           = os.path.dirname(border)
    border_diss   = gdb + '\\' + 'border_diss'
    holes_to_keep = gdb + '\\' + 'holes_to_keep'
    arcpy.Dissolve_management(border,border_diss)

    dic_point = {str([float('{0:.0f}'.format(pt.X)),float('{0:.0f}'.format(pt.Y))]):[pt.X,pt.Y] for i in arcpy.SearchCursor(border_diss) for part in i.Shape for pt in part if pt}

    arcpy.Delete_management (gdb + '\\' + 'border_diss')

    lyr_management = Layer_Management    (layer_path)
    lyr_management.Multi_to_single       ()
    lyr_management.Fill_Holes_in_Polygon (holes_to_keep,False,True)

    Ucursor = arcpy.UpdateCursor(layer_path)
    for i in Ucursor:
        ring = arcpy.Array()
        geom = i.Shape
        j = json.loads(geom.JSON)
        if 'curve' not in str(j):
            for part in geom:
                counter = 0
                for pt in part:
                    if pt:
                        if counter == 0:
                            first_pt = pt
                        key = str([float('{0:.0f}'.format(pt.X)),float('{0:.0f}'.format(pt.Y))])
                        if dic_point.get(key):
                            ring.add(arcpy.Point(dic_point[key][0],dic_point[key][1]))
                        else:
                            ring.add(pt)
                        counter = counter + 1
                    else:
                        ring.add(first_pt)
                        ring.add(None)
                        counter = 0

            polygon = arcpy.Polygon(ring)
            i.Shape = polygon
            Ucursor.updateRow(i) 
        else:
            pass

    
    Delete_polygons(layer_path,holes_to_keep)


def Spatial_Connection_To_LabelPoint(layer,ref,field_to_pass = []):

    '''
    [INFO] - חיבור מרחבי בין יישויות
    INPUT:
    1) layer         - השכבה אליה יכנסו הערכים לשדות חדשים
    2) ref           - השכבה ממנה ישלפו הערכים והשדות
    3) field_to_pass - שדות אותם רוצים להעביר, במידה ולא ינתן, יעביר את כל השדות
    OUTPUT:
    1) layer         - אותה שכבה עם השדות החדשים
    '''

    if field_to_pass == []:
        field_to_pass  = [i for i in Layer_Management(ref).fields() if i != 'SHAPE']

    # שכבת עזר 
    temp_intersect = 'in_memory\intersectAOI3'

    # שכבה אליה יכנסו השדות החדשות
    layer       = Layer_Management(layer)
    # המספר מזהה של השכבה
    OID         = layer.oid

    # יצירת נקודות במרכז היישות, 
    LabelPoint  = layer.Get_Label_Point_As_Point()
    LabelPoint  = Layer_Management(LabelPoint)

    # שם השכבה הראשית אצל שכבת הנקודות, הקשר בין השכבות
    ID_Field    = 'FID_'+LabelPoint.name # OID after intersect

    # חיבור בין שכבת הנקודות לשכבה חדשה
    arcpy.Intersect_analysis([ref,LabelPoint.layer],temp_intersect)

    # הכנה של השדות לפני הכנסה לחיבור הטבלאי בין הנקודות לשכבה הראשית
    field_to_pass = ''.join([i+';' for i in field_to_pass])[:-1]

    arcpy.JoinField_management(layer.layer, OID, temp_intersect, ID_Field, field_to_pass)

    LabelPoint.Destroy_layer()
    del temp_intersect

def Average(lst): 
    return sum(lst) / len(lst) 


def delete_parts_if_inside(layer,delete):

    '''
    [INFO] - Delete all Polygons 
    '''

    layer          = Layer_Management(layer)
    LabelPoint     = layer.Get_Label_Point_As_Point()
    LabelPoint     = Layer_Management(LabelPoint)

    LabelPoint.Select_By_Location('INTERSECT',delete,invert = "INVERT")

    layer.Select_By_Location('INTERSECT',LabelPoint.layer)


def Split_Line_By_Vertex(aoi_line):

    Multi_to_single(aoi_line)
    New_Line  = aoi_line + '_Temp'
    save_name = aoi_line

    arcpy.Select_analysis(aoi_line, New_Line, "\"OBJECTID\" < 0")
    iCursor = arcpy.da.InsertCursor(New_Line, ["SHAPE@"])
    with arcpy.da.SearchCursor(aoi_line,["SHAPE@"]) as sCursor:
        for row in sCursor:
            for part in row[0]:
                prevX = None
                prevY = None
                for pnt in part:
                    if pnt:
                        if prevX:
                            array = arcpy.Array([arcpy.Point(prevX, prevY),
                                                arcpy.Point(pnt.X, pnt.Y)])
                            polyline = arcpy.Polyline(array)
                            iCursor.insertRow([polyline])
                        prevX = pnt.X
                        prevY = pnt.Y
                    else:
                        pass

    del iCursor

    arcpy.Delete_management                (aoi_line)
    arcpy.Rename_management                (New_Line,save_name)


def Multi_to_single(layer):
    
    multi = False
    len_before = int(str(arcpy.GetCount_management(layer)))
    temp_lyer = layer  + 'Temp'
    save_name = layer
    arcpy.MultipartToSinglepart_management (layer,temp_lyer)
    arcpy.Delete_management                (layer)
    arcpy.Rename_management                (temp_lyer,save_name)
    len_after = int(str(arcpy.GetCount_management(layer)))
    if len_after > len_before:
        multi = True

    return multi


def Delete_By_area(layer,dis=0.2):
    with arcpy.da.UpdateCursor(layer,['SHAPE@AREA']) as ucursor:
        for row in ucursor:
            if row[0] < dis:
                ucursor.deleteRow()
    del ucursor


def Delete_Duplic_Line(fc):

	del_layer    = 'in_memory' + '\\' + 'arc_inter'
	diss_layer   = 'in_memory' + '\\' + 'diss_layer'
	Append_layer = 'in_memory' + '\\' + 'Append_layer'

	arcpy.Intersect_analysis          ([fc],del_layer)

	if int(str(arcpy.GetCount_management(del_layer))) > 0:

                del_layer_temp = 'in_memory' + '\\' + 'Temp'
                arcpy.Dissolve_management(del_layer,del_layer_temp)

                geom_del = [row.shape for row in arcpy.SearchCursor (del_layer_temp)][0]
                Ucursor  = arcpy.UpdateCursor (fc)
                for row in Ucursor:
                        for row in Ucursor:
                            if geom_del:
                                if row.shape:
                                    geom_up     = row.shape
                                    new_geom    = geom_up.difference(geom_del)
                                    row.shape = new_geom
                                    Ucursor.updateRow (row)


                arcpy.Dissolve_management              (del_layer,diss_layer)
                arcpy.MultipartToSinglepart_management (diss_layer,Append_layer)
                arcpy.Append_management                (Append_layer,fc,"NO_TEST")


def Layer_To_Edge_list(layer):

    '''
    [INFO] -  מוצא את זוג הוורטקסים הקרובים ביותר לכל קצה של ישות, ומחזיר רשימה עם הקארדינטות של שניהם והמרחק
    INPUT: 
    1) layer = שכבה קווית
    OUTPUT = [[ID_1,x1,y1],[ID_2,x2,y2],distance] 
    '''

    def dis(x1,y1,x2,y2):
        dist = math.sqrt(((x1-x2)**2) + ((y1-y2)**2))
        return dist

    data = [[row.OBJECTID,pt.X,pt.Y] for row in arcpy.SearchCursor(layer) for part in row.shape for pt in part]

    df            = pd.DataFrame(data,columns=["OBJECTID","X","Y"])
    df['index1']  = df.index


    gb_obj = df.groupby(by = 'OBJECTID')

    df_min = gb_obj.agg({'index1' : np.min})
    df_max = gb_obj.agg({'index1' : np.max})

    df_edge = pd.concat([df_min,df_max])

    df2 = pd.merge(df,df_edge, how='inner', on='index1')

    df_edge = df2.values.tolist()
    df_list = df.values.tolist()

    new_list = []
    for n in range(len(df_edge)):
        min_list = []
        dict1    = {}
        for i in range(len(df_list)):
            if df_edge[n][0] != df_list[i][0]:
                dist = dis(df_edge[n][1],df_edge[n][2],data[i][1],data[i][2])
                if dist != 0:
                    min_list.append(dist)
                    dict1[dist] = df_list[i]

        if min_list:
            min_l = min(min_list)
            new_list.append([df_edge[n][:-1],dict1[min_l][:-1],min_l])
        else:
            print ("part have no match type")

    return new_list

def Delete_By_length(layer,dis=0.2):
    with arcpy.da.UpdateCursor(layer,['SHAPE@LENGTH']) as ucursor:
        for row in ucursor:
            if row[0] < dis:
                ucursor.deleteRow()
    del ucursor


def fix_tolerance_line(layer_path,border):

    Multi_to_single(layer_path)

    dic_point = {str([float('{0:.1f}'.format(pt.X)),float('{0:.1f}'.format(pt.Y))]):[pt.X,pt.Y] for i in arcpy.SearchCursor(border) for part in i.Shape for pt in part if pt}

    Ucursor = arcpy.UpdateCursor(layer_path)
    for i in Ucursor:
        geom = i.Shape
        array = arcpy.Array()
        j = json.loads(geom.JSON)
        if 'curve' not in str(j):
            for part in geom:
                for pt in part:
                    if pt:
                        key = str([float('{0:.1f}'.format(pt.X)),float('{0:.1f}'.format(pt.Y))])
                        if dic_point.get(key):
                            # dis_moved = dis(pt.X, dic_point[key][0],  pt.Y, dic_point[key][1])
                            # print dis_moved # בדיקה כמה זה מתן
                            array.add(arcpy.Point(dic_point[key][0],dic_point[key][1]))
                        else:
                            array.add(arcpy.Point(pt.X, pt.Y))
                    else:
                        array.add(None)
            polyline = arcpy.Polyline(array)
            i.Shape = polyline
            Ucursor.updateRow(i)
            
        else:
            pass
       

def del_Non_Boundery_Line(line_layer,aoi_Final,tazar_border):
    '''
    [INFO] - מוחק קווים שנקודות הלייבל שלהם לא נמצאת על החלקות
    INPUT:
    1) line_layer = שכבת הקווים ממנה נמחק את הקווים שלא יושבים על החלקות
    2) aoi_Final  = שכבת אזור העבודה אליה תהיה השוואה
    3) tazar_border = הכלי יעבוד רק על קווים שיחתכו את שכבה זו
    ''' 
    Line_cut    = line_layer + '_Cut'
    label_point = line_layer + '_labelPoint'

    #get point from line in AOI 
    Layer_Management(Layer_Management(line_layer).Select_By_Location('INTERSECT',tazar_border,'5 Meters',Line_cut)).Get_Label_Point_As_Point(label_point)
    Layer_Management(label_point).Select_By_Location("BOUNDARY_TOUCHES",aoi_Final)

    Layer_Management(line_layer).Select_By_Location("INTERSECT",label_point)

    arcpy.Delete_management(Line_cut)
    arcpy.Delete_management(label_point)

def Get_layer_gdb(gdb):

    parcel_bankal = gdb + '\\' + 'PARCEL_ALL_EDIT'
    arc_bankal    = gdb + '\\' + 'PARCEL_ARC_EDIT'
    point_bankal  = gdb + '\\' + 'PARCEL_NODE_EDIT'

    parcel_modad  = gdb + '\\' + 'PARCELS_inProc_edit'
    arc_modad     = gdb + '\\' + 'LINES_inProc_edit'
    point_modad   = gdb + '\\' + 'POINTS_inProc_edit'

    return parcel_bankal,arc_bankal,point_bankal,parcel_modad,arc_modad,point_modad

def Get_layer_gdb_Copy(gdb):

    parcel_bankal_c = GDB + '\\' + 'PARCEL_ALL_EDIT_copy'
    arc_bankal_c    = GDB + '\\' + 'PARCEL_ARC_EDIT_copy'
    point_bankal_c  = GDB + '\\' + 'PARCEL_NODE_EDIT_copy'

    parcel_modad_c  = GDB + '\\' + 'PARCELS_inProc_edit_copy'
    arc_modad_c     = GDB + '\\' + 'LINES_inProc_edit_copy'
    point_modad_c   = GDB + '\\' + 'POINTS_inProc_edit_copy'

    return parcel_bankal_c,arc_bankal_c,point_bankal_c,parcel_modad_c,arc_modad_c,point_modad_c



scriptPath = os.path.abspath(__file__)
Scripts    = os.path.dirname(scriptPath)
ToolShare  = os.path.dirname(Scripts)
Scratch    = ToolShare + "\\Scratch"
ToolData   = ToolShare + "\\ToolData"

parcels_bankal         = arcpy.GetParameterAsText(0)
Folder                 = Scratch
Dis_limit_border_pnts  = 1
sett                   = ToolData + '\\' + r'Set.gdb\Sett'
CURRENT                = r'CURRENT'


print_arcpy_message     ("# # # # # # # S T A R T # # # # # #",status = 1)


layers_to_Copy  = ['PARCEL_ALL_EDIT','PARCEL_ARC_EDIT','PARCEL_NODE_EDIT','PARCELS_inProc_edit','LINES_inProc_edit','POINTS_inProc_edit']

GDB_Source      = getLayerPath     (parcels_bankal,CURRENT)
GDB             = CreateWorkingGDB (GDB_Source,Folder,layers_to_Copy,'Parcels_inProc_edit',CURRENT)

# קריאה של שכבות העבודה מבסיס הנתונים המקורי והחדש
parcel_bankal    ,arc_bankal    ,point_bankal    ,parcel_modad     , arc_modad    ,  point_modad    = Get_layer_gdb       (GDB_Source)
print_arcpy_message (parcel_bankal)
parcel_bankal_c  ,arc_bankal_c  ,point_bankal_c  ,parcel_modad_c   , arc_modad_c  ,  point_modad_c  = Get_layer_gdb_Copy  (GDB)
print_arcpy_message (parcel_bankal_c)
# שכבות עזר

Continue   = True                         # בדיקת צורך בהמשך פעולות גאומטריות, יהיה שלילי כאשר לא היו בעיות גאומטריות 
holes_2    = GDB + '\\' + 'Holes_Check_2' # return from func: CheckResultsIsOK(AOI,tazar_border,2)
inter2     = GDB + '\\' + 'Intersect_Check_2'
AOI2       = GDB + '\\' + 'AOI2'          # clean_slivers_by_vertex
AOI3       = GDB + '\\' + 'AOI3'          # after using: fix_holes_by_neer_length
AOI_best   = ''                           # return the last AOI Layer
AOI_final  = GDB + '\\' + 'AOI_Final'
AOI_Point  = GDB + '\\' + 'AOI_Point'
AOI_Line   = GDB + '\\' + 'AOI_Line'
AOI_Fix    = GDB + '\\' + 'Fix_holes'

diss_aoi   = r'in_memory' + '\\' +'diss_aoi'

if CheckIfSkipProcess(parcel_bankal_c,parcel_modad_c,GDB):
    print_arcpy_message     ("Exit",status = 1)
    sys.exit(0)

ChangeFieldNames                   (parcel_modad_c,arc_modad_c,point_modad_c)  # שינוי שמות השדות לפורמט בנק"ל
connect_parcel_to_sett             (parcel_modad_c,sett,parcel_bankal_c)       # הזנה של שמות ומספרי מפתח של ישובים לחלקות ניסיון ראשון

AOI,tazar_border,Curves,parcel_Bankal_cut,Point_bankal_Cut = PrePare_Data    (parcel_bankal_c,parcel_modad_c,point_modad_c,point_bankal_c,GDB,'POINT_NAME','POINT_NAME')

Get_Attr_From_parcel               (AOI,parcel_modad_c)                        #במידה והכלי מאתר שכל הישובים מסביב אותו דבר connect_parcel_to_sett גורס את הפעולה של
Fix_curves                         (AOI,tazar_border,Curves)
add_err_pts_to_mxd                 (GDB, ToolData + "\\lyr_files", ToolData + "\\demo.gdb",CURRENT) # parcels_bankal[1] = Current


if CheckResultsIsOK(AOI,tazar_border,1):
    AOI_best  = AOI
    Continue  = False

    # # # # # # # # בודק אם יש שינוי בגבולת אך ללא שינוי בתוך התצר, במידה ואין שינוי בגבולות, מחליף בין היישויות
Sub_Processing(parcel_bankal,parcel_modad_c,point_bankal,point_modad,arc_bankal,arc_modad,tazar_border,AOI,GDB)

if Continue:
    Dis_border_pnts = get_default_Snap_border (Point_bankal_Cut,parcel_modad_c,Dis_limit_border_pnts) 
    Snap_border_pnts        (GDB , tazar_border ,AOI,Dis_border_pnts) # סתימת חורים ע"י הזזת נקודות גבול
    Update_Polygons         (AOI , parcel_modad_c)
    Fix_curves              (AOI,tazar_border,Curves)
    AOI_best  = AOI

    if CheckResultsIsOK(AOI,tazar_border,2):                
        Continue   = False

if arcpy.Exists(holes_2):
    if int(str(arcpy.GetCount_management(holes_2))) == 0 and arcpy.Exists(inter2):
        holes_2 = inter2
else:
    holes_2 = inter2

    # # # # # # # במידה ואין חורים, ויש אינטרסקטים קטנים ממוצע של 1, הכלי ידלג על המשך הפעולות הגאומטריות

if Continue:
    names_curves              (AOI ,parcel_modad_c,Curves)
    clean_slivers_by_vertex   (AOI ,holes_2,tazar_border,2,AOI2)
    Update_Layer_Curves_By_ID (AOI2,parcel_modad_c,Curves)
    Clean_non_exist_pnts      (AOI2,tazar_border ,parcel_bankal_c,parcel_modad_c)
    Update_Layer_Curves_By_ID (AOI2,parcel_modad_c,Curves)
    AOI_best  = AOI2
        # # # # # # Check if all holes are closed after the Geometry tool
    if CheckResultsIsOK(AOI2,tazar_border,3):
        Continue   = False

if Continue:
        # # # # # # # Work Only if there is still Holes
    fix_holes_Overlaps_By_Length  (AOI2,tazar_border   ,AOI3)      # סתימת חורים ע"י שיוך לפי קו גדול משותף גדול ביותר
    clean_pseudo                  (AOI3,tazar_border   ,Curves)
    Update_Layer_Curves_By_ID     (AOI3,parcel_modad_c ,Curves)

    AOI_best = AOI3

CheckResultsIsOK(AOI_best,tazar_border,5)

fix_holes_Overlaps_By_Length  (AOI_best,tazar_border   ,AOI_Fix) 
stubborn_parts                (AOI_Fix,parcel_bankal_c,parcel_modad_c,AOI_final,Curves)
Update_Layer_Curves_By_ID     (AOI_final,parcel_modad_c,Curves)

fix_tolerance                 (AOI_final,tazar_border)
get_no_node_vertex            (AOI_final,tazar_border,point_modad_c,Point_bankal_Cut)
Delete_curves_out_AOI         (AOI_final,parcel_bankal)
Fix_Multi_part_Bankal         (AOI_final,tazar_border,parcel_Bankal_cut) # מתקן חלקות רחוקות שנפגעו בגלל שיש בהן חורים
Update_Polygons               (AOI_final,parcel_modad_c)
CheckResultsIsOK              (AOI_final,tazar_border,6)                # בדיקת סופית, כמה חורים נשארו

#  #  #  #  #  # # Prepare Insert to Razaf  #  #  #  #  #  #  # 
print_arcpy_message ("  #   #   #    # Preper data For Insert  #   #   #   #  ")

Calculate_Area_Rashum   (AOI_final)
NewGushim               (parcel_modad_c, parcel_bankal,AOI_final)
Get_Point_AOI           (AOI_final,point_bankal_c,point_modad_c,AOI_Point)
Create_Line_AOI         (AOI_final,tazar_border,Curves,arc_bankal_c,arc_modad_c,AOI_Line)

#  #  #  #  #  #  # insert To Razaf #  #  #  #  #  #  #  # #

print_arcpy_message ("  #   #   #    # insert To Razaf  #   #   #   #  ")

# # # # # Polygons 

Update_Polygons           (parcel_bankal,AOI_final)

# # # # Lines

arcpy.Dissolve_management (AOI_final,diss_aoi)
Layer_Management          (arc_bankal).Select_By_Location ('COMPLETELY_WITHIN',diss_aoi)
arcpy.Append_management   (AOI_Line,arc_bankal,'NO_TEST')

Multi_to_single           (arc_bankal)

arcpy.Append_management  (arc_modad_c,arc_bankal,'NO_TEST')

del_Non_Boundery_Line    (arc_bankal,AOI_final,tazar_border)

Find_stubbern_lines      (arc_bankal,AOI_final,tazar_border)
Delete_Duplic_Line       (arc_bankal)

# # # # Points

bankal_pnts = Layer_Management (point_bankal).Select_By_Location('INTERSECT',AOI_final)

arcpy.Append_management        (AOI_Point,point_bankal)

#    #   #    #  #

print_arcpy_message ("  #   #   #    #  Last Checks  #   #   #   #  ")

Parcel_data         (parcel_bankal,parcel_bankal_c,parcel_modad_c)


print_arcpy_message     ("# # # # # # # #     F I N I S H     # # # # # # # # #",1)
