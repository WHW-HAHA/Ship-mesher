import load_linesplan
import pro_ship
import panel_size
from interpolated_shipnew import interpolated_ship as Iship
from Patran import Patran
from Pat2VTK_custom import vtk

'''
Hanwei Wangd RHDHV 7_2_2019
'''
import numpy as np


# scaling
scaling_ratio = [1, 1, 1]
water_depth = 20
linesplan = load_linesplan.linesplan()
linesplan.get_info()
linesplan.get_fraction()
linesplan.get_coordinate()

# write coordinate in .txt file
text_name = 'Ship.txt'
with open(text_name, 'w') as ship_text:
    for crossection in linesplan.coordinate:
        for point in crossection:
            point = list(point)
            ship_text.write(str(point[0])+',' + str(point[2])+',' + str(point[1]) + '\n')
    ship_text.close()


# calculat the ship main particulars
model = pro_ship.Pro_ship()
model.main_particulars(linesplan.coordinate, linesplan.B/linesplan.L, linesplan.D/linesplan.L, linesplan.name,linesplan.L, linesplan.B, linesplan.T, linesplan.D)

# calculate the panel size(length of the panel)
Panel_size = panel_size.panel_size(water_depth)
Panel_size.get_max_frequency(model.omega_heave, model.omega_roll, model.omega_pitch)
Panel_size.panel_size = Panel_size.panel_size()

# interpolate the points on the ship frame
Post_ship = Iship(linesplan.L, linesplan.B, linesplan.T, model.num_crossection, model.out_line, model.half_crossection, Panel_size.panel_size )
Post_ship.ship_interpolate()

# Write Patran file
Pat = Patran(Post_ship.X, Post_ship.Y, Post_ship.Z, Post_ship.order, Post_ship.coor)
Pat.Pat_write()

# # convert pat to vtk
VTK_file = vtk()
VTK_file.create_vtk(Post_ship.Lookup_table, Pat.Allpoints)


# model.plot_frame()
scaled_x = np.array(model.x) * scaling_ratio[0]  # column list
scaled_y = np.array(model.y) * scaling_ratio[1]
scaled_z = np.array(model.z) * scaling_ratio[2]
linesplan.scaled_coordinates = list(zip(scaled_x, scaled_y , scaled_z))

# # calculate the panel size(length of the panel)
# Panel_size = panel_size.panel_size(water_depth)
# Panel_size.get_max_frequency(model.omega_heave, model.omega_roll, model.omega_pitch)
# Panel_size.panel_size = Panel_size.panel_size()


# main particulars of the scaled ship
linesplan.scaled_L = linesplan.L * scaling_ratio[0]
linesplan.scaled_B = linesplan.B * scaling_ratio[1]
linesplan.scaled_D = linesplan.D * scaling_ratio[2]
linesplan.scaled_T = linesplan.T * scaling_ratio[2]
linesplan.scaled_kxx = 0.35 * linesplan.scaled_B
linesplan.scaled_kyy = 0.25 * linesplan.scaled_L
linesplan.scaled_kzz = 0.25 * linesplan.scaled_L
linesplan.scaled_KG = 0.5 * linesplan.scaled_D

# calculate the panel size


