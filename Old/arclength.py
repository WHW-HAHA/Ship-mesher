import numpy as np
import math

'''
Hanwei Wang RHDHV 19-02-2019
'''


class arclength():
    def Length_cal(self, data):
        '''
        :param data: 3 * n matrix
        :return: length of curve
        '''
        length = []
        x = data[0, :]
        y = data[1, :]
        z = data[2, :]
        for i in range(1, len(x)):
            dl_i = np.sqrt((x[i] - x[i - 1]) ** 2 + (y[i] - y[i - 1]) ** 2 + (z[i] - z[i - 1]) ** 2)
            length.append(dl_i)
        self.length_i = sum(length)
        self.x = x
        self.y = y
        self.z = z
        return sum(length)

    def new_coor(self, panel_size):
        num_steps = math.ceil(self.length_i/panel_size)
        new_x = np.linspace(min(self.x), max(self.x), num_steps)
        new_y = np.linspace(min(self.y), max(self.y), num_steps)
        new_z = np.linspace(min(self.z), max(self.z), num_steps)
        new_data = np.vstack((new_x, new_y, new_z))
        return new_data
