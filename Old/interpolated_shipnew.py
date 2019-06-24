from scipy import interpolate
from pro_ship import Pro_ship
from load_linesplan import linesplan
from panel_size import panel_size
from scipy.interpolate import pchip_interpolate
import math
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from copy import copy
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import matplotlib.pyplot as plt

'''
Hanwei Wang RHDHV 19-02-2019
'''


class interpolated_ship():
    '''
    This class will interpolate the points on ship hull before editing draft
    '''

    def __init__(self, L, B, T, num_crossection, outlines, half_crossection, panel_size):
        self.L = L
        self.B = B
        self.T = T
        self.num_crossection = num_crossection
        self.outlines = outlines
        self.half_crossection = half_crossection
        self.panel_size = panel_size
        self.panel_size = 2

    def cosspace(self, startpoint, endpoint, numpoints):
        '''
        COSSPACE(X1, X2) generates a row vector of 100 cosine-spaced points between X1 and X2.
        This method of spacing concentrates samples at the ends while producing fewer sample points in the center
        '''
        start_angle = math.pi / 2
        x_fore = []
        x_aft = []
        x_fore.append(startpoint)
        # x.append(endpoint)
        if endpoint <= startpoint:
            print('End point must greater than start point')
            x = []
            return x
        else:
            midpoint = (endpoint - startpoint) / 2
        angleInc = (2 * start_angle) / (numpoints - 1)
        curAngle = -start_angle + angleInc

        for idx in range(1, math.ceil(numpoints / 2) + 1):
            x_fore.append(startpoint + midpoint + (1 / math.sin(start_angle)) * midpoint * math.sin(curAngle))
            x_aft.append(
                endpoint - (startpoint + midpoint + (1 / math.sin(start_angle)) * midpoint * math.sin(curAngle)))
            curAngle = curAngle - angleInc
        x_aft.sort()
        x = x_fore + x_aft
        return x

    def equal_interpolate(self, x_series, y_series, z_series):
        '''
        interpolate data in 3D curve with equal step
        '''

        M = 1000
        t = np.linspace(0, len(y_series), M)
        y_new1 = np.interp(t, np.arange(len(y_series)), y_series)
        z_new1 = np.interp(t, np.arange(len(z_series)), z_series)
        tol = self.panel_size  # panel_size
        k, idx = 0, [0]
        while k < len(y_new1):
            total_dist = 0
            for j in range(k + 1, len(y_new1)):
                total_dist += math.sqrt((y_new1[j] - y_new1[j - 1]) ** 2 + (z_new1[j] - z_new1[j - 1]) ** 2)
                if total_dist > tol:
                    idx.append(j)
                    break
            k = j + 1
        y_new1 = y_new1[idx]
        z_new1 = z_new1[idx]
        x_new1 = np.ones(len(y_new1)) * x_series[0]
        return x_new1, y_new1, z_new1

    def equal_interpolate_3D(self, x_series, y_series, z_series):
        '''
        interpolate data in 3D curve with equal step
        '''
        M = 1000
        t = np.linspace(0, len(y_series), M)
        x_new1 = np.interp(t, np.arange(len(x_series)), x_series)
        y_new1 = np.interp(t, np.arange(len(y_series)), y_series)
        tol = self.panel_size  # panel_size
        k, idx = 0, [0]
        while k < len(x_new1):
            total_dist = 0
            for j in range(k + 1, len(x_new1)):
                total_dist += math.sqrt((x_new1[j] - x_new1[j - 1]) ** 2 + (y_new1[j] - y_new1[j - 1]) ** 2)
                if total_dist > tol:
                    idx.append(j)
                    break
            k = j + 1
        x_new1 = x_new1[idx]
        y_new1 = y_new1[idx]
        z_new1 = np.ones(len(y_new1)) * max(z_series)
        return x_new1, y_new1, z_new1

    def arclength(self, x, y, z):
        '''
        return the length of curve
        '''
        length = []
        for i in range(1, len(x)):
            dl_i = np.sqrt((x[i] - x[i - 1]) ** 2 + (y[i] - y[i - 1]) ** 2 + (z[i] - z[i - 1]) ** 2)
            length.append(dl_i)
        self.length_i = sum(length)
        return sum(length)

    def return_values(self, series):
        try:
            print(len(series))
        except:
            print(series)
        data = []
        try:
            for item in series:
                if item != None:
                    data.append(item)
        except:
            pass
        return data

    def ship_interpolate(self):
        # define smooth keel line
        cuts = np.linspace(0, self.L, self.L / self.panel_size)
        keelpoints = np.zeros([3, len(cuts)])
        keel = self.outlines[:, 3:6]
        keelpoints[2, :] = pchip_interpolate(keel[:, 0], keel[:, 2], cuts)
        keelpoints[0, :] = cuts
        # new_keel_x = np.linspace(max(keel[:,0]) - min(keel[:,0])
        x_slices = self.half_crossection[0:-2:3]
        y_slices = self.half_crossection[1:-1:3]
        z_slices = self.half_crossection[2::3]
        x_crossection = {}
        y_crossection = {}
        z_crossection = {}
        Name = []
        for i in range(self.num_crossection):
            name = 'crossection' + str(i)
            Name.append(name)
            x_crossection[name], y_crossection[name], z_crossection[name] = self.equal_interpolate(x_slices[i],
                                                                                                   y_slices[i],
                                                                                                   z_slices[i])
        mid_num_crossection = np.ceil(self.num_crossection / 2)
        # new_x_crossection use the same ram as x_crossecton, therefore changes together
        # copy to duplicate with 2 ram
        new_x_crossection = copy(x_crossection)
        new_y_crossection = copy(y_crossection)
        new_z_crossection = copy(z_crossection)

        for o in range(self.num_crossection - 1):
            name_1 = Name[o]
            name_2 = Name[o + 1]
            # first insert more points in between cross-section if possible
            # if num_crossection is odd number , neglect the middle crossection
            # if num_crossection is even number
            extra = len(x_crossection[name_2]) - len(x_crossection[name_1])
            dx = abs((x_crossection[name_2][0] - x_crossection[name_1][0])) / (abs(extra) + 1)
            dz = abs((min(z_crossection[name_1]) - min(z_crossection[name_2]))) / (abs(extra) + 1)

            if extra > 0:  # in first half crossection
                for i in range(1, extra + 1):
                    new_x_crossection[name_1] = np.insert(new_x_crossection[name_1], 0,
                                                          min(x_crossection[name_1]) + dx * i)
                    new_y_crossection[name_1] = np.insert(new_y_crossection[name_1], 0, 0)
                    # interplate Z
                    # new_z_crossection[name_1] =  pchip_interpolate(keelpoints[0,:], keelpoints[2,:], new_x_crossection[name_1])

                    new_z_crossection[name_1] = np.insert(new_z_crossection[name_1], 0,
                                                          pchip_interpolate(keelpoints[0, :], keelpoints[2, :],
                                                                            min(x_crossection[name_1]) + dx * i))
            elif extra < 0:  # in next half crossection
                for i in range(1, abs(extra) + 1):
                    new_x_crossection[name_2] = np.insert(new_x_crossection[name_2], 0,
                                                          min(x_crossection[name_2]) - dx * i)
                    new_y_crossection[name_2] = np.insert(new_y_crossection[name_2], 0, 0)
                    new_z_crossection[name_2] = np.insert(new_z_crossection[name_2], 0,
                                                          pchip_interpolate(keelpoints[0, :], keelpoints[2, :],
                                                                            min(x_crossection[name_2]) - dx * i))
        fig = plt.figure()
        ax = Axes3D(fig)
        K = 0

        for i, j in zip(Name[0:-1], Name[1:]):
            K += 1
            # ax.scatter(new_x_crossection[i], new_y_crossection[i], new_z_crossection[i], c='r')
            # if K < mid_num_crossection:
            #     '''
            #     new_crossection[i] is longer
            #     '''
            #     for k in range(len(new_x_crossection[i])):
            #         ax.plot([new_x_crossection[i][k], x_crossection[j][k]],
            #                 [new_y_crossection[i][k], y_crossection[j][k]],
            #                 [new_z_crossection[i][k], z_crossection[j][k]], c='black')
            # else:
            #     for k in range(len(x_crossection[i])):
            #         ax.plot([x_crossection[i][k], new_x_crossection[j][k]],
            #                 [y_crossection[i][k], new_y_crossection[j][k]],
            #                 [z_crossection[i][k], new_z_crossection[j][k]],
            #                 c='black')

            ax.scatter(new_x_crossection[j], new_y_crossection[j], new_z_crossection[j], c='r')
        ax.get_proj = lambda: np.dot(Axes3D.get_proj(ax), np.diag([1.5, 1.5 * 0.2, 1.5 * 0.1, 1]))
        plt.show()

        # redefine spline
        fig = plt.figure()
        ax = Axes3D(fig)
        num_vertical_layers = len(
            new_x_crossection['crossection{}'.format(int(mid_num_crossection))])  # !!!! need to be changed
        spline_x = []
        spline_y = []
        spline_z = []
        Spline_x = {}
        Spline_y = {}
        Spline_z = {}
        Hor_name = []
        for i in range(num_vertical_layers):
            ver_name = 'ver_layer' + str(i)
            Hor_name.append(ver_name)
            for j in range(self.num_crossection):
                name = Name[j]
                # if j < mid_num_crossection:
                try:
                    spline_x.append(new_x_crossection[name][-1 - i])
                    spline_y.append(new_y_crossection[name][-1 - i])
                    spline_z.append(new_z_crossection[name][-1 - i])

                except:
                    pass

            # calculate the length of spline
            new_spline_x, new_spline_y, new_spline_z = self.equal_interpolate_3D(spline_x, spline_y, spline_z)
            # second half
            temp_x = []
            temp_y = []
            temp_z = []

            for i in cuts:
                if i >= min(new_spline_x) and i <= max(new_spline_x):
                    if abs(i-min(new_spline_x)) <= self.panel_size:
                        temp_x.append(min(new_spline_x))
                        temp_y.append(float(interpolate.pchip_interpolate(new_spline_x, new_spline_y, min(new_spline_x))))
                        temp_z.append(new_spline_z[0])
                    elif abs(i-max(new_spline_x)) <= self.panel_size:
                        temp_x.append(max(new_spline_x))
                        temp_y.append(float(interpolate.pchip_interpolate(new_spline_x, new_spline_y, max(new_spline_x))))
                        temp_z.append(new_spline_z[0])
                    else:
                        temp_x.append(i)
                        temp_y.append(float(interpolate.pchip_interpolate(new_spline_x, new_spline_y, i)))
                        temp_z.append(new_spline_z[0])
                else:
                    temp_x.append(None)
                    temp_y.append(None)
                    temp_z.append(None)

            # ax.scatter(list(filter(lambda v: v != None, temp_x)), list(filter(lambda v: v != None, temp_y)),list(filter(lambda v: v != None, temp_z)), c='red')
            # # # # # # plt.show()
            # lock the start and end point on new keel line
            print(len(cuts))
            print(len(temp_y))
            Spline_x[ver_name] = temp_x
            for i in range(len(temp_y)):
                if temp_y[i] != None:
                    temp_y[i] = 0
                    # temp_z[i] = new_keel[i]
                    break
            for i in range(len(temp_y)):
                if temp_y[len(temp_y) - i - 1] != None:
                    temp_y[len(temp_y) - i - 1] = 0
                    # temp_z[len(temp_y) - i-1] = new_keel[len(temp_y) - i-1]
                    break
            ax.scatter(list(filter(lambda v: v != None, temp_x)), list(filter(lambda v: v != None, temp_y)),list(filter(lambda v: v != None, temp_z)), c='red')
            Spline_y[ver_name] = temp_y
            Spline_z[ver_name] = temp_z
            spline_x = []
            spline_y = []
            spline_z = []
        # # # # # plt.show()
        # reorganize the points in crossections（cuts）
        Hcross_name = ['cross' + str(i) for i in range(len(cuts))]  # strat with 'vc0'
        Hor_coord = {}
        num_cor = 1
        for i in range(len(cuts)):
            temp = []
            for ver_name in Spline_x.keys():
                temp.append(Spline_x[ver_name][i])
                temp.append(Spline_y[ver_name][i])
                temp.append(Spline_z[ver_name][i])
                if Spline_x[ver_name][i] != None:
                    temp.append(num_cor)
                    num_cor += 1
                else:
                    temp.append(None)
            Hor_coord[Hcross_name[i]] = np.array(list(filter(lambda x: x != None, temp))).reshape(-1, 4)
            try:
                Hor_coord[Hcross_name[i]][-1][1] = 0.0
            except:
                print(i)
            ax.scatter(Hor_coord[Hcross_name[i]][:][:, 0], Hor_coord[Hcross_name[i]][:][:, 1],
                       Hor_coord[Hcross_name[i]][:][:, 2], c='red')

            self.coor = []
            # collect the coordinate
        for name in Hor_coord.keys():
            for x, y, z, index in zip(Hor_coord[name][:, 0], Hor_coord[name][:, 1], Hor_coord[name][:, 2],
                                      Hor_coord[name][:, 3]):
                self.coor.append([x, y, z, index])

        # sort the data and patch, stern and bow are included in the ship body.
        # x,y,z for visualization
        x = []
        y = []
        z = []
        # X, Y, Z store data
        X = []
        Y = []
        Z = []
        order = []
        alpha = 0.9
        fc = "C0"
        for j in range(len(cuts) - 1):
            cross_fore = Hor_coord[Hcross_name[j]]
            cross_aft = Hor_coord[Hcross_name[j + 1]]
            len1 = len(cross_fore)
            len2 = len(cross_aft)
            for i in range(len1):
                try:
                    x.append(cross_fore[i][0])
                    y.append(cross_fore[i][1])
                    z.append(cross_fore[i][2])
                    X.append(cross_fore[i][0])
                    Y.append(cross_fore[i][1])
                    Z.append(cross_fore[i][2])
                    order.append(cross_fore[i][3])
                except:
                    X.append(None)
                    Y.append(None)
                    Z.append(None)
                    order.append(None)

                try:
                    x.append(cross_aft[i][0])
                    y.append(cross_aft[i][1])
                    z.append(cross_aft[i][2])
                    X.append(cross_aft[i][0])
                    Y.append(cross_aft[i][1])
                    Z.append(cross_aft[i][2])
                    order.append(cross_aft[i][3])
                except:
                    X.append(None)
                    Y.append(None)
                    Z.append(None)
                    order.append(None)

                try:
                    x.append(cross_aft[i + 1][0])
                    y.append(cross_aft[i + 1][1])
                    z.append(cross_aft[i + 1][2])
                    X.append(cross_aft[i + 1][0])
                    Y.append(cross_aft[i + 1][1])
                    Z.append(cross_aft[i + 1][2])
                    order.append(cross_aft[i + 1][3])
                except:
                    X.append(None)
                    Y.append(None)
                    Z.append(None)
                    order.append(None)

                try:
                    x.append(cross_fore[i + 1][0])
                    y.append(cross_fore[i + 1][1])
                    z.append(cross_fore[i + 1][2])
                    X.append(cross_fore[i + 1][0])
                    Y.append(cross_fore[i + 1][1])
                    Z.append(cross_fore[i + 1][2])
                    order.append(cross_fore[i + 1][3])
                except:
                    X.append(None)
                    Y.append(None)
                    Z.append(None)
                    order.append(None)

                # close the panel, only for visualization
                try:
                    x.append(cross_fore[i][0])
                    y.append(cross_fore[i][1])
                    z.append(cross_fore[i][2])
                except:
                    pass

                verts_1st_half = [list(zip(x, y, z))]
                x = []
                y = []
                z = []

                # for i in range(len1):
                #     try:
                #         x.append(cross_fore[i][0])
                #         y.append(-cross_fore[i][1])
                #         z.append(cross_fore[i][2])
                #         X.append(cross_fore[i][0])
                #         Y.append(-cross_fore[i][1])
                #         Z.append(cross_fore[i][2])
                #         order.append(cross_fore[i][3])
                #     except:
                #         pass
                #     try:
                #         x.append(cross_aft[i][0])
                #         y.append(-cross_aft[i][1])
                #         z.append(cross_aft[i][2])
                #         X.append(cross_aft[i][0])
                #         Y.append(-cross_aft[i][1])
                #         Z.append(cross_aft[i][2])
                #         order.append(cross_aft[i][3])
                #     except:
                #         pass
                #     try:
                #         x.append(cross_aft[i + 1][0])
                #         y.append(-cross_aft[i + 1][1])
                #         z.append(cross_aft[i + 1][2])
                #         X.append(cross_aft[i + 1][0])
                #         Y.append(-cross_aft[i + 1][1])
                #         Z.append(cross_aft[i + 1][2])
                #         order.append(cross_aft[i + 1][3])
                #     except:
                #         pass
                #     try:
                #         x.append(cross_fore[i + 1][0])
                #         y.append(-cross_fore[i + 1][1])
                #         z.append(cross_fore[i + 1][2])
                #         X.append(cross_fore[i + 1][0])
                #         Y.append(-cross_fore[i + 1][1])
                #         Z.append(cross_fore[i + 1][2])
                #         order.append(cross_fore[i + 1][3])
                #     except:
                #         pass
                #
                #
                # verts_sec_half = [list(zip(x,y,z))]
                pc_1st_half_face = Poly3DCollection(verts_1st_half, alpha=alpha, facecolors=fc, linewidths=1)
                pc_1st_half_wire_frame = Line3DCollection(verts_1st_half, colors='k', linewidths=1)
                # pc_2ec_half_face = Poly3DCollection(verts_sec_half, alpha=alpha, facecolors=fc, linewidths=3)
                # pc_2ec_half_wire_frame = Line3DCollection(verts_sec_half, colors='k', linewidths=1)
                # ax.add_collection3d(pc_1st_half_face)
                # ax.add_collection3d(pc_2ec_half_face)
                ax.add_collection3d(pc_1st_half_wire_frame)
                # ax.add_collection3d(pc_2ec_half_wire_frame)
        # ax.get_proj = lambda: np.dot(Axes3D.get_proj(ax), np.diag([1.5, 1.5 * 0.2, 1.5 * 0.1, 1]))
        print('Body points {}'.format(len(X)))
        self.Body_dot = len(X)
        # # # # # plt.show()

        # # Stern + water_lid + bow
        # stern_x = Hor_coord[Hcross_name[0]][0][0]
        # stern_y = sorted(Hor_coord[Hcross_name[0]][:,1])
        # stern_z = Hor_coord[Hcross_name[0]][:,2]
        # stern_cross_name = ['stern_cross' + str(i) for i in range(len(stern_y))]
        # stern_cross_x = {}
        # stern_cross_y = {}
        # stern_cross_z = {}
        # for i in range(len(stern_y)):
        #     name = stern_cross_name[i]
        #     stern_cross_x[name] = []
        #     stern_cross_y[name] = []
        #     stern_cross_z[name] = []
        #     for j in range(len(stern_z)-i):
        #         stern_cross_x[name].append(stern_x)
        #         stern_cross_y[name].append(stern_y[i])
        #         stern_cross_z[name].append(stern_z[j])
        #
        # x = []
        # y = []
        # z = []
        # for i in range(len(stern_y)-1):
        #     cross_fore = stern_cross_name[i]
        #     cross_aft = stern_cross_name[i+1]
        #     len1 = len(stern_cross_x[cross_fore])
        #     for i in range(len1-1):
        #         try:
        #             x.append(stern_cross_x[cross_fore][i])
        #             y.append(stern_cross_y[cross_fore][i])
        #             z.append(stern_cross_z[cross_fore][i])
        #         except:
        #             pass
        #         try:
        #             x.append(stern_cross_x[cross_aft][i])
        #             y.append(stern_cross_y[cross_aft][i])
        #             z.append(stern_cross_z[cross_aft][i])
        #         except:
        #             pass
        #         try:
        #             x.append(stern_cross_x[cross_aft][i+1])
        #             y.append(stern_cross_y[cross_aft][i+1])
        #             z.append(stern_cross_z[cross_aft][i+1])
        #         except:
        #             pass
        #         try:
        #             x.append(stern_cross_x[cross_fore][i+1])
        #             y.append(stern_cross_y[cross_fore][i+1])
        #             z.append(stern_cross_z[cross_fore][i+1])
        #         except:
        #             pass
        #         try:
        #             x.append(stern_cross_x[cross_fore][i])
        #             y.append(stern_cross_y[cross_fore][i])
        #             z.append(stern_cross_z[cross_fore][i])
        #         except:
        #             pass
        #
        #         verts_1st_stern = [list(zip(x, y, z))]
        #
        #         x = []
        #         y = []
        #         z = []
        #         # '''Second stern'''
        #         # for i in range(len1 - 1):
        #         #     try:
        #         #         x.append(stern_cross_x[cross_fore][i])
        #         #         y.append(-stern_cross_y[cross_fore][i])
        #         #         z.append(stern_cross_z[cross_fore][i])
        #         #         # X.append(cross_fore[i][0])
        #         #         # Y.append(cross_fore[i][1])
        #         #         # Z.append(cross_fore[i][2])
        #         #         # order.append(cross_fore[i][3])
        #         #     except:
        #         #         pass
        #         #     try:
        #         #         x.append(stern_cross_x[cross_aft][i])
        #         #         y.append(-stern_cross_y[cross_aft][i])
        #         #         z.append(stern_cross_z[cross_aft][i])
        #         #         # X.append(cross_aft[i][0])
        #         #         # Y.append(cross_aft[i][1])
        #         #         # Z.append(cross_aft[i][2])
        #         #         # order.append(cross_aft[i][3])
        #         #     except:
        #         #         pass
        #         #     try:
        #         #         x.append(stern_cross_x[cross_aft][i + 1])
        #         #         y.append(-stern_cross_y[cross_aft][i + 1])
        #         #         z.append(stern_cross_z[cross_aft][i + 1])
        #         #         # Y.append(cross_aft[i+1][1])
        #         #         # Z.append(cross_aft[i+1][2])
        #         #         # order.append(cross_aft[i+1][3])
        #         #     except:
        #         #         pass
        #         #     try:
        #         #         x.append(stern_cross_x[cross_fore][i + 1])
        #         #         y.append(-stern_cross_y[cross_fore][i + 1])
        #         #         z.append(stern_cross_z[cross_fore][i + 1])
        #         #         # X.append(cross_fore[i + 1][0])
        #         #         # Y.append(cross_fore[i + 1][1])
        #         #         # Z.append(cross_fore[i + 1][2])
        #         #         # order.append(cross_fore[i + 1][3])
        #         #     except:
        #         #         pass
        #         #     verts_2nd_stern = [list(zip(x, y, z))]
        #         pc_1st_stern_wire_frame = Line3DCollection(verts_1st_stern, colors='k', linewidths=1)
        #         pc_1st_stern_face = Poly3DCollection(verts_1st_stern, alpha=alpha, facecolors=fc, linewidths=1)
        #         # pc_2nd_stern_wire_frame = Line3DCollection(verts_2nd_stern, colors='k', linewidths=1)
        #         # pc_2nd_stern_face = Poly3DCollection(verts_2nd_stern, alpha=alpha, facecolors=fc, linewidths=3)
        #         # ax.add_collection3d(pc_1st_stern_wire_frame)
        #         # ax.add_collection3d(pc_2nd_stern_face)
        #         # ax.add_collection3d(pc_1st_stern_face)
        #         # ax.add_collection3d(pc_2nd_stern_wire_frame)

        lid_cross_name = ['lid_cross' + str(i) for i in range(len(cuts))]
        lid_cross_x = {}
        lid_cross_y = {}
        lid_cross_z = {}
        lid_index = {}
        # skip -1, since the last crossection is wrong

        for j in range(len(cuts) - 1):
            name = lid_cross_name[j]
            lid_cross_x[name] = []
            lid_cross_y[name] = []
            lid_cross_z[name] = []
            lid_index[name] = []
            if j == 0:
                lid_y = 0
                lid_cross_x[name].append(Hor_coord[Hcross_name[j]][0][1])
                lid_cross_y[name].append(lid_y)
                lid_cross_z[name].append(Hor_coord[Hcross_name[j]][0][2])
                lid_index[name].append(Hor_coord[Hcross_name[j]][0][3])

            elif j == len(cuts) - 2:
                lid_y = 0
                lid_cross_x[name].append(Hor_coord[Hcross_name[j]][0][0])
                lid_cross_y[name].append(lid_y)
                lid_cross_z[name].append(Hor_coord[Hcross_name[j]][0][2])
                lid_index[name].append(Hor_coord[Hcross_name[j]][0][3])

            else:
                lid_y = np.linspace(0, Hor_coord[Hcross_name[j]][0][1],
                                    math.ceil(Hor_coord[Hcross_name[j]][0][1] / self.panel_size))
                for y in lid_y[::-1]:  # from y = 0 to y = contour
                    if y != lid_y[-1]:
                        lid_cross_x[name].append(Hor_coord[Hcross_name[j]][0][0])
                        lid_cross_y[name].append(y)
                        lid_cross_z[name].append(Hor_coord[Hcross_name[j]][0][2])
                        lid_index[name].append(num_cor)  # skip the last one, since is alredy counted in ship body
                        self.coor.append([Hor_coord[Hcross_name[j]][0][0], y, Hor_coord[Hcross_name[j]][0][2], num_cor])
                        num_cor += 1
                    else:
                        lid_cross_x[name].append(Hor_coord[Hcross_name[j]][0][0])
                        lid_cross_y[name].append(y)
                        lid_cross_z[name].append(Hor_coord[Hcross_name[j]][0][2])
                        lid_index[name].append(Hor_coord[Hcross_name[j]][0][3])

        # The 1st lid is the first ver-crossection

        for o in range(0, len(cuts) - 2):
            # for o in [len(cuts) - 4]:
            # skip the last water lid
            x = []
            y = []
            z = []
            lid_fore = lid_cross_name[o]
            lid_aft = lid_cross_name[o + 1]
            len1 = len(lid_cross_x[lid_fore])
            len2 = len(lid_cross_x[lid_aft])
            # triangle panel at the first and last water lid cross
            # at stern
            if len1 == 1:
                try:
                    x.append(lid_cross_x[lid_fore][0])
                    y.append(lid_cross_y[lid_fore][0])
                    z.append(lid_cross_z[lid_fore][0])
                    X.append(lid_cross_x[lid_fore][0])
                    Y.append(lid_cross_y[lid_fore][0])
                    Z.append(lid_cross_z[lid_fore][0])
                    order.append(lid_index[lid_fore][0])
                except:
                    X.append(None)
                    Y.append(None)
                    Z.append(None)
                    order.append(None)
                try:
                    x.append(lid_cross_x[lid_aft][0])
                    y.append(lid_cross_y[lid_aft][0])
                    z.append(lid_cross_z[lid_aft][0])
                    X.append(lid_cross_x[lid_aft][0])
                    Y.append(lid_cross_y[lid_aft][0])
                    Z.append(lid_cross_z[lid_aft][0])
                    order.append(lid_index[lid_aft][0])
                except:
                    X.append(None)
                    Y.append(None)
                    Z.append(None)
                    order.append(None)
                try:
                    x.append(lid_cross_x[lid_aft][-1])
                    y.append(lid_cross_y[lid_aft][-1])
                    z.append(lid_cross_z[lid_aft][-1])
                    X.append(lid_cross_x[lid_aft][-1])
                    Y.append(lid_cross_y[lid_aft][-1])
                    Z.append(lid_cross_z[lid_aft][-1])
                    order.append(lid_index[lid_aft][-1])
                except:
                    X.append(None)
                    Y.append(None)
                    Z.append(None)
                    order.append(None)
                try:
                    x.append(lid_cross_x[lid_fore][0])
                    y.append(lid_cross_y[lid_fore][0])
                    z.append(lid_cross_z[lid_fore][0])
                    X.append(None)
                    Y.append(None)
                    Z.append(None)
                    order.append(None)
                except:
                    pass

            # at bow
            elif len2 == 1:
                try:
                    x.append(lid_cross_x[lid_fore][0])
                    y.append(lid_cross_y[lid_fore][0])
                    z.append(lid_cross_z[lid_fore][0])
                    X.append(lid_cross_x[lid_fore][0])
                    Y.append(lid_cross_y[lid_fore][0])
                    Z.append(lid_cross_z[lid_fore][0])
                    order.append(lid_index[lid_fore][0])
                except:
                    X.append(None)
                    Y.append(None)
                    Z.append(None)
                    order.append(None)
                try:
                    x.append(lid_cross_x[lid_aft][0])
                    y.append(lid_cross_y[lid_aft][0])
                    z.append(lid_cross_z[lid_aft][0])
                    X.append(lid_cross_x[lid_aft][0])
                    Y.append(lid_cross_y[lid_aft][0])
                    Z.append(lid_cross_z[lid_aft][0])
                    order.append(lid_index[lid_aft][0])
                except:
                    X.append(None)
                    Y.append(None)
                    Z.append(None)
                    order.append(None)
                try:
                    x.append(lid_cross_x[lid_fore][-1])
                    y.append(lid_cross_y[lid_fore][-1])
                    z.append(lid_cross_z[lid_fore][-1])
                    X.append(lid_cross_x[lid_fore][-1])
                    Y.append(lid_cross_y[lid_fore][-1])
                    Z.append(lid_cross_z[lid_fore][-1])
                    order.append(lid_index[lid_fore][-1])
                except:
                    X.append(None)
                    Y.append(None)
                    Z.append(None)
                    order.append(None)

                try:
                    x.append(lid_cross_x[lid_fore][0])
                    y.append(lid_cross_y[lid_fore][0])
                    z.append(lid_cross_z[lid_fore][0])
                    X.append(None)
                    Y.append(None)
                    Z.append(None)
                    order.append(None)
                except:
                    pass

            elif len1 <= len2:
                for i in range(len1):
                    try:
                        x.append(lid_cross_x[lid_fore][i])
                        y.append(lid_cross_y[lid_fore][i])
                        z.append(lid_cross_z[lid_fore][i])
                        X.append(lid_cross_x[lid_fore][i])
                        Y.append(lid_cross_y[lid_fore][i])
                        Z.append(lid_cross_z[lid_fore][i])
                        order.append(lid_index[lid_fore][i])
                    except:
                        X.append(None)
                        Y.append(None)
                        Z.append(None)
                        order.append(None)
                    try:
                        x.append(lid_cross_x[lid_aft][i])
                        y.append(lid_cross_y[lid_aft][i])
                        z.append(lid_cross_z[lid_aft][i])
                        X.append(lid_cross_x[lid_aft][i])
                        Y.append(lid_cross_y[lid_aft][i])
                        Z.append(lid_cross_z[lid_aft][i])
                        order.append(lid_index[lid_aft][i])
                    except:
                        X.append(None)
                        Y.append(None)
                        Z.append(None)
                        order.append(None)
                    try:
                        x.append(lid_cross_x[lid_aft][i + 1])
                        y.append(lid_cross_y[lid_aft][i + 1])
                        z.append(lid_cross_z[lid_aft][i + 1])
                        X.append(lid_cross_x[lid_aft][i + 1])
                        Y.append(lid_cross_y[lid_aft][i + 1])
                        Z.append(lid_cross_z[lid_aft][i + 1])
                        order.append(lid_index[lid_aft][i + 1])
                    except:
                        X.append(None)
                        Y.append(None)
                        Z.append(None)
                        order.append(None)
                    try:
                        x.append(lid_cross_x[lid_fore][i + 1])
                        y.append(lid_cross_y[lid_fore][i + 1])
                        z.append(lid_cross_z[lid_fore][i + 1])
                        X.append(lid_cross_x[lid_fore][i + 1])
                        Y.append(lid_cross_y[lid_fore][i + 1])
                        Z.append(lid_cross_z[lid_fore][i + 1])
                        order.append(lid_index[lid_fore][i + 1])
                    except:
                        X.append(None)
                        Y.append(None)
                        Z.append(None)
                        order.append(None)

                    try:
                        x.append(lid_cross_x[lid_fore][i])
                        y.append(lid_cross_y[lid_fore][i])
                        z.append(lid_cross_z[lid_fore][i])
                    except:
                        pass
            else:
                for i in range(len2):
                    try:
                        x.append(lid_cross_x[lid_fore][i])
                        y.append(lid_cross_y[lid_fore][i])
                        z.append(lid_cross_z[lid_fore][i])
                        X.append(lid_cross_x[lid_fore][i])
                        Y.append(lid_cross_y[lid_fore][i])
                        Z.append(lid_cross_z[lid_fore][i])
                        order.append(lid_index[lid_fore][i])
                    except:
                        X.append(None)
                        Y.append(None)
                        Z.append(None)
                        order.append(None)
                    try:
                        x.append(lid_cross_x[lid_aft][i])
                        y.append(lid_cross_y[lid_aft][i])
                        z.append(lid_cross_z[lid_aft][i])
                        X.append(lid_cross_x[lid_aft][i])
                        Y.append(lid_cross_y[lid_aft][i])
                        Z.append(lid_cross_z[lid_aft][i])
                        order.append(lid_index[lid_aft][i])
                    except:
                        X.append(None)
                        Y.append(None)
                        Z.append(None)
                        order.append(None)
                    try:
                        x.append(lid_cross_x[lid_aft][i + 1])
                        y.append(lid_cross_y[lid_aft][i + 1])
                        z.append(lid_cross_z[lid_aft][i + 1])
                        X.append(lid_cross_x[lid_aft][i + 1])
                        Y.append(lid_cross_y[lid_aft][i + 1])
                        Z.append(lid_cross_z[lid_aft][i + 1])
                        order.append(lid_index[lid_aft][i + 1])
                    except:
                        X.append(None)
                        Y.append(None)
                        Z.append(None)
                        order.append(None)
                    try:
                        x.append(lid_cross_x[lid_fore][i + 1])
                        y.append(lid_cross_y[lid_fore][i + 1])
                        z.append(lid_cross_z[lid_fore][i + 1])
                        X.append(lid_cross_x[lid_fore][i + 1])
                        Y.append(lid_cross_y[lid_fore][i + 1])
                        Z.append(lid_cross_z[lid_fore][i + 1])
                        order.append(lid_index[lid_fore][i + 1])
                    except:
                        X.append(None)
                        Y.append(None)
                        Z.append(None)
                        order.append(None)

                    try:
                        x.append(lid_cross_x[lid_fore][i])
                        y.append(lid_cross_y[lid_fore][i])
                        z.append(lid_cross_z[lid_fore][i])
                    except:
                        pass

            lid_1st_half = [list(zip(x, y, z))]
            pc_1st_lid_wire_frame = Line3DCollection(lid_1st_half, colors='k', linewidths=1)
            # pc_1st_lid_face = Poly3DCollection(lid_1st_half, alpha=alpha, facecolors=fc, linewidths=1)
            # ax.add_collection3d(pc_1st_lid_face)
            ax.add_collection3d(pc_1st_lid_wire_frame)
        ax.get_proj = lambda: np.dot(Axes3D.get_proj(ax), np.diag([1.5, 1.5 * 0.2, 1.5 * 0.1, 1]))

        print('Total points: {}'.format(len(X)))
        self.Lid_dot = len(X) - self.Body_dot
        # # # # # plt.show()

        # collecting the data and write them in a text file.
        # remove the whole panel with more than 2 rows Nones
        # assign look up table
        i_start = 0
        i_end = i_start + 4
        X_new = []
        Y_new = []
        Z_new = []
        order_new = []
        self.Lookup_table = []
        num = 1
        while num <= len(X) / 4:
            num += 1
            count = 0
            for x, y, z, index in zip(X[i_start: i_end], Y[i_start: i_end], Z[i_start: i_end], order[i_start:i_end]):
                if x == None:
                    count += 1
            if count < 2:
                for x, y, z, index in zip(X[i_start: i_end], Y[i_start: i_end], Z[i_start: i_end],
                                          order[i_start:i_end]):
                    X_new.append(x)
                    Y_new.append(y)
                    Z_new.append(z)
                    order_new.append(index)
                if num <= self.Body_dot/4:
                    self.Lookup_table.append(0)
                else:
                    self.Lookup_table.append(2)
            i_start += 4
            i_end = i_start + 4
        self.X = X_new
        self.Y = Y_new
        self.Z = Z_new
        self.order = order_new
