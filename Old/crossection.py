import numpy as np

# separate the ship drat into number of slices
class crossection():
    def __init__(self):
        self.g = 9.81
        self.wls = 10


    def water_line_y(self, draft, coordinates):
        Y = np.linspace(0.01, draft, self.wls)
        pass





