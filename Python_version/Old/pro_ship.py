import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import pchip_interpolate
import math

'''
Hanwei Wang RHDHV 15_2_2019
'''
class Pro_ship():
    def __init__(self):
        self.phase = None
        self.x = []
        self.y = []
        self.z = []

    def plot_all(self):
        fig = plt.figure()
        ax = Axes3D(fig)
        pass

    def simpsons(self, arguements, values):
        ''''
        Use simpsons iteration to calculate the area and volume
        '''
        s1 = 0
        s2 = 0
        h = (arguements[-1] - arguements[0]) / (2 * len(arguements))
        index_1 = np.arange(0, len(arguements)-1)
        index_2 = np.arange(1, len(arguements))
        for index in index_1:
            s1 += values[index]
        for index in index_2:
            s2 = s2 + values[index]
        s = h * (values[0] + values[-1] + 4 * s1 + 2 * s2) / 3
        return s

    def main_particulars(self, data, ratio_y_x, ratio_z_x, name, Lpp, Breadth, Draft, Depth ):
        x = []
        y = []
        z = []
        fig = plt.figure()
        ax = Axes3D(fig)
        crossections = []
        half_crossection = []
        self. num_crossection = 0
        j = 0
        # k = 0
        for cross_section in data:
            self.num_crossection += 1
            for point in list(cross_section):
                # print(point)
                self.x.append(point[0])
                if point[2] in y:
                    j += 1
                    self.y.append(point[2]+ 0.1*j)
                    y.append(point[2]+ 0.1*j)
                else:
                    self.y.append(point[2])
                    y.append(point[2])
                self.z.append(point[1])
                x.append(point[0])
                z.append(point[1])
            j = 0
            Crossection = np.vstack((x, y, z))
            Crossection = Crossection[:, Crossection[1, :].argsort()]  # Crossection is the half ship body
            ##############
            Crossection_next = np.vstack((Crossection[0,1:], -Crossection[1,1:], Crossection[2,1:]))
            Crossection_next = Crossection_next[:, Crossection_next[1,:].argsort()]
            Crossection_whole = np.hstack((Crossection_next[:, :], Crossection))

            # for point in list(cross_section):
            #     # print(point)
            #     if point[2] != 0:
            #         self.x.append(point[0])
            #         if -point[2] in self.y:
            #             k += 1
            #             self.y.append(-point[2])
            #             y.append(-point[2]-0.01*k)
            #         else:
            #             self.y.append(-point[2])
            #             y.append(-point[2])
            #         self.z.append(point[1])
            #         x.append(point[0])
            #         z.append(point[1])
            #
            # crossection = np.vstack((x, y, z)) # crossection is the whole body, but only for visualization
            # crossection = crossection[:, crossection[1, :].argsort()]  # 按照数组的一个维度排序
            # print('j = {}'.format(j))
            # print('k = {}'.format(k))
            x = []
            y = []
            z = []

            # plot the crossections
            # ax.plot(Crossection_whole[0, :], Crossection_whole[1, :], Crossection_whole[2, :], c='r')

            # plot the water lid
            # ax.plot([Crossection_whole[0, 0], Crossection_whole[0, -1]], [Crossection_whole[1, 0], Crossection_whole[1, -1]],
            #         [Crossection_whole[2, 0], Crossection_whole[2, -1]], c='g')
            half_crossection.extend(Crossection)
            crossections.extend(Crossection_whole)

            # plot the vertical water plane

            ax.scatter(self.x[0: :2], self.y[0: :2], self.z[0: :2], c='b', marker="o", alpha=0.5)
            # ax.plot(self.x[0::2], self.y[0::2], self.z[0::2], c='b', alpha=0.5)
            # ax.plot(  )
        self.half_crossection = half_crossection
        self.crossections = crossections

        '''
        Data collecting ends here
        Start data Processing 
        '''
        '''
        # 1, calculate the water plane area
        '''
        # First create a nutural line at surface
        nutural_line_1 = np.ones([self.num_crossection, 3])
        nutural_line_1[:,2] = self.crossections[2][-1]
        x_slices = self.crossections[0:-2:3]
        for i in range(self.num_crossection):
            nutural_line_1[i,0] = x_slices[i][0]
        nutural_line_1[:,1] = 0

        # Then create nutural line at keel
        nutural_line_2 = np.ones([self.num_crossection, 3])
        nutural_line_2[:, 0] = nutural_line_1[:,0]
        z_slices = self.crossections[2: :3]
        for i in range(self.num_crossection):
            nutural_line_2[i, 2] = min(z_slices[i])
        nutural_line_2[:, 1] = 0

        # Then create contour line
        contour_line_1 = np.ones([self.num_crossection, 3])
        contour_line_1[:,2] = self.crossections[2][-1]
        contour_line_1[:,0] = nutural_line_1[:,0]
        y_slices = self.crossections[1:-1:3]
        for i in range(self.num_crossection):
            contour_line_1[i,1] = max(y_slices[i])

        # Next half of the contour line
        contour_line_2 = np.ones([self.num_crossection, 3])
        contour_line_2[:,2] = self.crossections[2][-1]
        contour_line_2[:,0] = nutural_line_1[:,0]
        y_slices = self.crossections[1:-1:3]
        for i in range(self.num_crossection):
            contour_line_2[i,1] = -max(y_slices[i])

        # plot
        # plt.plot(nutural_line_1[:, 0], nutural_line_1[:, 1], nutural_line_1[:, 2], c='g', linewidth='1')
        # plt.plot(nutural_line_2[:, 0], nutural_line_2[:, 1], nutural_line_2[:, 2], c='red')
        # plt.plot(contour_line_1[:, 0], contour_line_1[:, 1], contour_line_1[:, 2], c='g', linewidth='1')
        # plt.plot(contour_line_2[:, 0], contour_line_2[:, 1], contour_line_2[:, 2], c='g', linewidth='1')
        self.out_line = np.hstack((nutural_line_1, nutural_line_2, contour_line_1, contour_line_2))

        # calculate the water plane area by simpson iteration
        self.water_plane = self.simpsons(nutural_line_1[:, 0], contour_line_1[:, 1]) * 2

        print('Water plane area is {} m^2'.format(int(self.water_plane)))

        '''
        #2, calculate the areas of all the cross-sections and then the displacement of the ship
        '''
        area_cross_sections = []
        for i in range(self.num_crossection):
            area_cross_sections.append(nutural_line_1[0,2] * contour_line_1[i,1] - self.simpsons(half_crossection[3*i+1], half_crossection[3*i+2]))
        self.displacement = abs(self.simpsons(nutural_line_1[:, 0], area_cross_sections)*2)
        print('Displacement is {} m^3'.format(self.displacement))

        # calculate the LCF
        # 垂向计算法， 先计算各个水线面面积， 然后将水线面面积按照吃水方向进行积分拉计算浮心位置
        # Yf = 0, Xf
        self.LCF = 2 * self.simpsons(nutural_line_1[:,0],nutural_line_1[:,0]*contour_line_1[:,1])/self.water_plane
        '''
        #3 calculate the ship main particulars
        '''
        # Inertia for all the crossections
        self.I = 2/3 * self.simpsons(nutural_line_1[:,0], contour_line_1[:,1]*contour_line_1[:,1]*contour_line_1[:,1])
        # CB
        self.CB = self.displacement/(Lpp*Breadth*Draft)
        # BM
        self.BM = self.I/ self.displacement
        # KB
        '''
        some doubts
        '''
        self.KB = abs(self.simpsons(contour_line_1[:,1], self.water_plane * contour_line_1[:,1])/self.displacement)
        # LCB
        # self.LCB = self.simpsons(contour_line_1[:,1], self.water_plane * (Xf - Lpp/2))/self.displacement
        ''' May not necessary, change it later. '''
        self.LCB = 0.2
        # GM
        ''' KG is user input, here assume is half of depth'''
        self.GM = self.KB + self.BM - Depth/2
        # k_xx, k_yy
        self.k_yy = 0.3 * Lpp
        self.k_xx = 0.25 * Breadth

        '''
        #4 calculate the natural frequencies of ship motions 
        '''
        self.omega_heave = math.sqrt(self.water_plane * 9.81 / self.displacement)  # pay attention for the unit
        self.omega_roll = math.sqrt(self.GM * 9.81/ self.k_xx ** 2)
        self.I_l = 2*self.simpsons(nutural_line_1[:,0], contour_line_1[:,1] * nutural_line_1[:,0] ** 2) - self.LCF ** 2 * self.water_plane
        self.BM_l = self.I_l/self.displacement
        self.GM_l = self.KB + self.BM_l - Depth/2
        self.omega_pitch = math.sqrt(self.GM_l * 9.81/(self.k_yy ** 2))

        print('Heave $\omega$ = {} rad/s\nRoll $\omega$ = {} rad/s\nPitch $\omega$ = {} rad/s'.format(self.omega_heave, self.omega_roll, self.omega_pitch))

        ax.get_proj = lambda: np.dot(Axes3D.get_proj(ax), np.diag([1.5, 1.5 * ratio_y_x, 1.5 * ratio_z_x, 1]))
        ax.set_xlabel('Length [m]')
        ax.set_ylabel('Breadth [m]')
        ax.set_zlabel('Draft [m]')
        ax.set_title(name)
        plt.show()



