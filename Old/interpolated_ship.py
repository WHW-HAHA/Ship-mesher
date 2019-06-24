from scipy import interpolate
from pro_ship import Pro_ship
from load_linesplan import linesplan
from panel_size import panel_size
from scipy.interpolate import pchip_interpolate
import math
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


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

    def cosspace(self, startpoint, endpoint, numpoints):
        '''
        COSSPACE(X1, X2) generates a row vector of 100 cosine-spaced points between X1 and X2.
        This method of spacing concentrates samples at the ends while producing fewer sample points in the center
        '''
        start_angle = math.pi/2
        x_fore = []
        x_aft = []
        x_fore.append(startpoint)
        # x.append(endpoint)
        if endpoint <= startpoint:
            print('End point must greater than start point')
            x = []
            return x
        else:
            midpoint= (endpoint - startpoint)/2
        angleInc = (2*start_angle)/(numpoints-1)
        curAngle = -start_angle + angleInc
        for idx in range(1, math.ceil(numpoints/2)+1):
            x_fore.append(startpoint + midpoint + (1/math.sin(startpoint))* midpoint*math.sin(curAngle))
            curAngle = curAngle - angleInc
        x_aft = (-1* x_fore + endpoint).sort()
        print(x_fore)
        print(x_aft)
        x = x_fore.extend(x_aft)
        print(x)
        return x


    def xdevide(self):
        pass


    def arclength(self, x, y, z):
        length = []
        for i in range(1, len(x)):
            dl_i = np.sqrt((x[i] - x[i - 1]) ** 2 + (y[i] - y[i - 1]) ** 2 + (z[i] - z[i - 1]) ** 2)
            length.append(dl_i)
        self.length_i = sum(length)
        self.x = x
        self.y = y
        self.z = z
        return sum(length)

    def ship_interpolate(self):
        # define smooth keel line
        keel = self.outlines[:,6:9]

        """
        x_new = np.linspace(keel[0,0], keel[-1,0], 100)
        ynew1 = interpolate.pchip_interpolate(keel[:,0], keel[:,1], x_new)
        tck = interpolate.splrep(keel[:,0],keel[:,1])
        ynew2 = interpolate.splev(x_new,tck)
        fig = plt.figure()
        plt.plot(keel[:,0], keel[:,1])
        plt.plot(x_new, ynew1, 'o')
        plt.plot(x_new, ynew2 , '*')

        plt.legend(['ori', 'pchip','Bspline'])
        # plt.show()
        """
        num_per_crossection = []

        # new_keel_x = np.linspace(max(keel[:,0]) - min(keel[:,0])
        x_slices = self.half_crossection[0:-2:3]
        y_slices = self.half_crossection[1:-1:3]
        z_slices = self.half_crossection[2::3]
        X_slices = x_slices
        Y_slices = y_slices
        Z_slices = z_slices

        fig = plt.figure()
        ax = Axes3D(fig)

        # using dictionary to store coordinates in each crossection, each crossection is 3 * n numpy matrix

        for i in range(self.num_crossection):

            y_new1 = np.linspace(y_slices[i][0],y_slices[i][-1],20)
            z_new1 = interpolate.pchip_interpolate(y_slices[i], z_slices[i], y_new1)
            z_new2 = np.linspace(z_slices[i][0],y_slices[i][-1],5)
            y_new2 = interpolate.pchip_interpolate(z_slices[i], y_slices[i], z_new2)

            y_new = np.append(y_new1, y_new2)
            z_new = np.append(z_new1, z_new2)

            print(y_new)
            print(z_new)
            x_new = np.ones(len(z_new)) * x_slices[i][0]
            ax.scatter(x_new, y_new, z_new)
        ax.get_proj = lambda: np.dot(Axes3D.get_proj(ax), np.diag([1.5, 1.5 * 0.2, 1.5 * 0.1, 1]))
        plt.show()




        for o in range(self.num_crossection):
            num_per_crossection.append(len(x_slices[o]))
        print(num_per_crossection)

        # Extra data are necessary to keep the number of points in
        # all crossections identical, and then connect lines

        # for o in range(self.num_crossection-1):
        #     # first insert more points in between cross-section if possible
        #     # and directly add these point in original crossection
        #     extra = len(x_slices[o+1]) - len(x_slices[o])
        #     if extra>0:
        #         dx = (x_slices[o+1][0] - x_slices[o][0])/(extra+1)
        #         dz = (z_slices[o][0])/(extra+1)
        #         # dz = 0
        #         for i in range(1, extra+1):
        #             if o < self.num_crossection/2:
        #                 X_slices[o] = np.append(X_slices[o], X_slices[o][0] + dx*i)
        #                 Y_slices[o] = np.append(Y_slices[o], [0])
        #                 '''
        #                 Approximate the z values
        #                 '''
        #                 Z_slices[o] = np.append(Z_slices[o], [Z_slices[o][0] - dz*i])
        #             elif o > self.num_crossection/2+1:
        #                 X_slices[o+1] = np.append(X_slices[o+1], [X_slices[o+1][0] -i*dx])
        #                 Y_slices[o+1] = np.append(Y_slices[o+1], 0)
        #                 Z_slices[o+1] = np.append(Z_slices[o+1], [Z_slices[o][0] - dz*i])
        #
        #     elif extra<0:
        #         extra = abs(extra)
        #         dx = (x_slices[o+1][0] - x_slices[o][0])/(extra+1)
        #         dz = (z_slices[o][0])/(extra+1)
        #         # dz = 0
        #         for i in range(1, extra+1):
        #             if o < self.num_crossection/2:
        #                 X_slices[o + 1] = np.append(X_slices[o + 1], X_slices[o][0] - dx*i)
        #                 Y_slices[o + 1] = np.append(Y_slices[o], [0])
        #                 '''
        #                 Approximate the z values
        #                 '''
        #                 Z_slices[o + 1] = np.append(Z_slices[o], [Z_slices[o][0] - dz*i])
        #             elif o > self.num_crossection/2+1:
        #                 X_slices[o] = np.append(X_slices[o], [X_slices[o][0] -i*dx])
        #                 Y_slices[o] = np.append(Y_slices[o], 0)
        #                 Z_slices[o] = np.append(Z_slices[o], [Z_slices[o][0] - dz*i])

        for o in range(self.num_crossection-1):
            # first insert more points in between cross-section if possible
            # and directly add these point in original crossection
            extra = len(x_slices[o+1]) - len(x_slices[o])
            if extra>0:
                extra = extra-1
                dx = (x_slices[o+1][0] - x_slices[o][0])/(extra+1)
                dz = (z_slices[o][0]-z_slices[o+1][0])/(extra+1)
                for i in range(1, extra+1):
                        X_slices[o] = np.append(x_slices[o],  x_slices[o][0] + dx*i)
                        Y_slices[o] = np.append(y_slices[o],  [0])
                        Z_slices[o] = np.append(z_slices[o],  [z_slices[o][0] - dz*i])
            elif extra<0:
                extra = abs(extra)-1
                dx = (x_slices[o+1][0] - x_slices[o][0])/(extra+1)
                dz = (z_slices[o+1][0]-z_slices[o][0])/(extra+1)
                # dz = 0x0
                for i in range(1, extra+1):
                        X_slices[o+1] = np.append(x_slices[o+1],  [x_slices[o+1][0] - i*dx])
                        Y_slices[o+1] = np.append(y_slices[o+1],  0)
                        Z_slices[o+1] = np.append(z_slices[o+1],  [z_slices[o+1][0] - dz*i])


        # Get number of new lengthwise splines needed
        m= []
        for o in range(self.num_crossection):
            # m.append(len(X_slices[o]))
            m.append(len(x_slices[o]))

        print('New_num_points_per_section is \b')
        print(m)
        Figure = plt.figure()
        ax = Axes3D(Figure)
        for o in range(self.num_crossection):
            ax.scatter(X_slices[o], Y_slices[o], Z_slices[o], 'b')
        for o in range(19):
            for j in range(self.num_crossection-1):
                try:
                    ax.plot([x_slices[j][-1-o],x_slices[j+1][-1-o]], [y_slices[j][-1-o],y_slices[j+1][-1-o]], [z_slices[j][-1-o],z_slices[j+1][-1-o]])
                except:
                    pass
        ax.get_proj = lambda: np.dot(Axes3D.get_proj(ax), np.diag([1.5, 1.5 * 0.2, 1.5 * 0.1, 1]))
        plt.show()




        for i in range(max(m)):
            c = 0
            x = []
            y = []
            z = []
            points = []
            strook = []
            for o in range(self.num_crossection):
                if len(X_slices[o]) > i:
                    x.append(X_slices[o])
                    y.append(Y_slices[o])
                    z.append(Z_slices[o])
                    c += 1
                ax.scatter(x_slices[o], y_slices[o], z_slices[o], 'b')
            ax.get_proj = lambda: np.dot(Axes3D.get_proj(ax), np.diag([1.5, 1.5 * 0.2, 1.5 * 0.1, 1]))
            plt.show()

            if len(x[0]) > 1:
                # Create the vertical slices
                points.extend(list(zip(x, y, z)))
                '''
                Strook is the interpolated spline in vertical direction
                不知道是在干嘛 pending
                '''
                strook.append(None)

        # use cosin spacing to add more inter-sections
        # new crossection
        cuts = self.cosspace(min(X_slices[0]), max(X_slices[-1]), 30)
        keelpoints = np.zeros([3, len(cuts)])
        keelpoints[0,:] = cuts
        keelpoints[2,:] = pchip_interpolate(keel[:,0], keel[:,2], cuts)
        keelpoints[1,:] = 0

        #
        # divs_1 = []
        # divs_2 = []
        # length_spline = []
        # s = []
        # x_frame = self.half_crossection[0:-2:3]
        # y_frame = self.half_crossection[1:-1:3]
        # z_frame = self.half_crossection[2::3]
        # for i in range(self.num_crossection):
        #     length = self.arclength(x_frame[i], y_frame[i], z_frame[i])
        #     length_spline.append(length)
        #     num = round(length/self.Panel_size) + 1
        #     divs_1.append(num)
        #     y_inter = np.linspace(min(y_frame[i]), max(y_frame[i]), num)
        #     # evaluate the z values at y_inter
        #
        # pass
        #







