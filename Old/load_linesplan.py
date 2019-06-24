import numpy as np
import tkinter as tk
import linecache
import math
import re

'''
Hanwei Wangd RHDHV 7_2_2019
'''

class linesplan():
    def __init__(self):
        self.L = None
        self.B = None
        self.T = None  # draft of the ship
        self.D = None  # depth of the ship
        self.name = None
        self.type = None
        self.coordinate = None
        self.num_cross_section = None
        self.fractions = None
        self.scaling_flag = None
        self.kxx = None
        self.kyy = None
        self.KG = None
        self.displacement = None
        self.file_path = r'C:\Users\907932\Desktop\Software\Diffrac_related\Diffrac database\001_linesplan database\linesplans\coaster_v08.linesplan'
        self.file_path = r'C:\Users\907932\Desktop\Software\Diffrac_related\Diffrac database\001_linesplan database\Original linesplan database Journee Versluis\Versluis.004'


    def get_info(self):
        # get the originald ship dimensions(L, B, T, D)
        with open(self.file_path, 'r') as plan:
            # line1 = linecache.getline(self.file_path,2).split(' ')
            line1 = linecache.getline(self.file_path,2)
            # LBTD = re.search('\d{1,3}.\d{2}', line1)
            regex = re.compile(r'\d{1,3}.\d{2}.*\d{2}')
            LBTD = regex.search(line1).group(0)
            LBTD = LBTD.split('x')
            self.L = float(LBTD[0])
            self.B = float(LBTD[1])
            TD = LBTD[2].split('(')
            self.T = float(TD[0])
            self.D = float(TD[1])
            self.name = line1.split(' ')[0]
            self.type = line1.split(' ')[1]
            self.num_cross_section = int(linecache.getline(self.file_path,4).strip('\n').split(' ')[-1])
        plan.close()


    def get_fraction(self):

        fractions = []
        fraction_start = 4
        fraction_end = 4 + 1 + math.ceil((self.num_cross_section + 1)/6)
        for i in np.arange(fraction_start+1,fraction_end):
            line = linecache.getline(self.file_path, i)
            for fraction in line.strip('\n').split(' '):
                # print(fraction)
                if fraction != '':
                    fractions.append(float(fraction))
        self.fractions = fractions


    def get_coordinate(self):
        '''

        '''
        # First find the origin of the ship
        # X_0 is the stern of ship, y_0 is central line, z_0 is keel
        fraction_sum = [0]
        for i in range(1, len(self.fractions)+1):   # 千万注意啊
            fraction_sum.append(sum(self.fractions[0: i]))

        # print(' fraction: {}'.format(fraction_sum))
        x_cors = np.array(fraction_sum) * self.L

        with open(self.file_path, 'r') as plan:
             content = plan.read().split('\n')
             Content = []
             effective_content = content[4 + math.ceil((self.num_cross_section)/6) + 1: -3]
             for items in effective_content:
                 for item in items.split(' '):
                     if item != '':
                        Content.append(float(item))
             count = 0
             Coordinates = []

             while Content:
                 index = 1
                 num_points_per_section = int(Content[index]) + 1
                 y_cor = np.array(Content[index + 3: index + num_points_per_section*2 + 3 : 2]) * self.T
                 y_cor = np.sort(y_cor)
                 z_cor = np.array(Content[index + 2: index + num_points_per_section*2 + 2 : 2]) * self.B # !!!!!!!
                 z_cor = np.sort(z_cor)
                 # if count == 0:
                 #     y_cor = np.insert( y_cor, 0, self.T)
                 #     z_cor = np.insert(z_cor,0, 0)
                 Content = Content[index + num_points_per_section*2 + 2: ] # 不是-1, 因为列表不会包含最后一个值
                 # organize coordinate per section
                 try:
                     x_cor = np.ones(len(y_cor)) * x_cors[count]
                     coordinates = list(zip(x_cor, y_cor, z_cor))
                     Coordinates.append(coordinates)
                     count += 1
                 except:
                     print('The ship model is wrong, you may lose the last cross-section!')

             plan.close()
             # Coordinates[0] = np.insert(Coordinates[0], 0, (0 , 0, self.T))
             # Coordinates[-1] = np.append(Coordinates[-1], (self.L, 0, self.T))
             # self.num_cross_section = self.num_cross_section
             self.coordinate = Coordinates


    











