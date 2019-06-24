from scipy import optimize
import math

'''
Hanwei Wangd RHDHV 18_2_2019
'''

class panel_size():
    def __init__(self, water_depth):
        self.water_depth = water_depth
    '''
    calculate the panelsize, with input of wave_length
    '''
    def get_max_frequency(self, w_heave, w_roll, w_pitch):
        max_w = max(w_heave, w_roll, w_pitch)
        self.max_w = max_w

    def wave_dispersion(self, Lambda):
        return 9.81 * (2 * math.pi/ self.max_w) ** 2 / 2 / math.pi * math.tanh(2 * math.pi * self.water_depth / Lambda) - Lambda

    def panel_size(self):
        # for very shallow water, the simplified wave length follows formula: L = T * sqrt(g*d)
        length_app = 2*math.pi/self.max_w * math.sqrt(9.81 * self.water_depth)
        self.Wave_length = optimize.fsolve(self.wave_dispersion, length_app)
        self.Panel_size = int(math.floor(self.Wave_length/4))
        print('Wave length is {} m, panel size is {} m'.format(self.Wave_length, self.Panel_size))
        return self.Panel_size










