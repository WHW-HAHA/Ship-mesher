import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import PchipInterpolator
from scipy.interpolate import pchip_interpolate
from scipy.interpolate import griddata
from math import *


'''
Hanwei Wang RHDHV 8_2_2019
'''
class Visual():
    def __init__(self):
        self.phase = None
        self.x = []
        self.y = []
        self.z = []

    def plot_all(self):
        fig = plt.figure()
        ax = Axes3D(fig)
        pass

    def plot_points(self, data, ratio_y_x, ratio_z_x, name):
        x = []
        y = []
        z = []
        fig = plt.figure()
        ax = Axes3D(fig)
        crossections = []
        self. num_crossection = 0
        j = 0
        k = 0
        for cross_section in data:
            self.num_crossection += 1
            for point in list(cross_section):
                # print(point)
                self.x.append(point[0])
                if point[2] in self.y:
                    j += 1
                    self.y.append(point[2])
                else:
                    self.y.append(point[2])
                self.z.append(point[1])
                x.append(point[0])
                y.append(point[2])
                z.append(point[1])
            Crossection = np.vstack((x, y, z))
            Crossection = Crossection[:, Crossection[1, :].argsort()]

            for point in list(cross_section):
                # print(point)
                if point[2] != 0:
                    self.x.append(point[0])
                    if point[2] in self.y:
                        k += 1
                        self.y.append(-point[2])
                        y.append(-point[2]-0.1)
                    else:
                        self.y.append(-point[2])
                        y.append(-point[2])
                    self.z.append(point[1])
                    x.append(point[0])
                    z.append(point[1])

            crossection = np.vstack((x, y, z))
            crossection = crossection[:, crossection[1, :].argsort()]  # 按照数组的一个维度排序
            x = []
            y = []
            z = []
            ax.plot(crossection[0, :], crossection[1, :], crossection[2, :], c='r')
            # plot the water lid
            ax.plot([crossection[0, 0], crossection[0, -1]], [crossection[1, 0], crossection[1, -1]],
                    [crossection[2, 0], crossection[2, -1]], c='g')
            crossections.extend(Crossection)
            ax.scatter(self.x, self.y, self.z, c='b', marker="o", alpha=0.5)

        self.crossections = crossections
        self.data_interpolation()

        for i in range(self.num_crossection-1):
            print(i)
            self.frame[3*i] = np.hstack((self.frame[3 * i], self.frame[3 * i]))
            self.frame[3*i+1] = np.hstack((self.frame[3 * i + 1],-self.frame[3 * i + 1]))
            self.frame[3*i+2] = np.hstack((self.frame[3 * i + 2], self.frame[3 * i + 2]))

        for i in range(self.num_crossection):
            try:
                ax.scatter(self.frame[3 * i], self.frame[3 * i + 1], self.frame[3 * i + 2], c='black', marker=".",
                           alpha=0.5)
            except:
                pass

        ax.get_proj = lambda: np.dot(Axes3D.get_proj(ax), np.diag([1.5, 1.5 * ratio_y_x, 1.5 * ratio_z_x, 1]))
        ax.set_xlabel('Length [m]')
        ax.set_ylabel('Breadth [m]')
        ax.set_zlabel('Draft [m]')
        ax.set_title(name)
        plt.show()



    def data_interpolation(self):
        self.frame = []
        # for i in range(self.num_crossection)[-1]:
        for i in range(self.num_crossection):
            try:
                self.fram_z = pchip_interpolate(self.crossections[i * 3 + 1], self.crossections[i * 3 + 2], np.linspace(self.crossections[i * 3 + 1][0], self.crossections[i * 3 + 1][-1]))
                self.fram_y = np.linspace(self.crossections[i * 3 + 1][0], self.crossections[i * 3 + 1][-1])
                self.fram_x = np.ones((1, self.fram_y.size)) * self.crossections[i * 3][0]
                self.frame.extend((self.fram_x, self.fram_y, self.fram_z))
            except:
                print('Wrong again')

    def simpsons_iteration(self, data, width, n):
        # width
        sum = data[0] + data[n-1]
        for i in range(2, n):
            if i%2 == 0:
                sum =  sum + 4*data[i-1]
            else:
                sum = sum + 2*data[i-1]
            return sum*width/3.0
        pass










'''
    def plot_frame(self):
            self.fram_y_z = PchipInterpolator(np.array(self.y), np.array(self.z))
            self.fram_y = self.fram_y_z.c
            self.fram_z = self.fram_y_z(self.fram_y)
            self.fram_y = np.reshape(self.fram_y, [self.fram_y.size, 1])
            self.fram_z = np.reshape(self.fram_z, [self.fram_z.size, 1])
            self.fram_x = self.x[0] * np.ones(self.fram_y.size)
            # self.fram.append(np.concatenate(self.fram_x, self.fram_y, self.fram_z, axis = 1))
            self.fram.append(list(zip(self.fram_x, self.fram_y, self.fram_z)))




            # print(self.x)
            fig = plt.figure()
            ax = Axes3D(fig)
            ax.scatter(self.x, self.y, self.z)
            ax.get_proj = lambda: np.dot(Axes3D.get_proj(ax), np.diag([1.5, 1.5*ratio_y_x, 1.5*ratio_z_x, 1]))
            ax.set_xlabel('Length [m]')
            ax.set_ylabel('Breadth [m]')
            ax.set_zlabel('Draft [m]')
            ax.set_title(name)
            # ax.set_aspect('equal')
        except:
            pass
        plt.show()
        '''





