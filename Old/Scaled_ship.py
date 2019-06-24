from .load_linesplan import linesplan

class scaled_ship(linesplan):

    def scaling(self, scaling_ratio):
        linesplan.L = linesplan.L * scaling_ratio[0]
        linesplan.D = linesplan.D * scaling_ratio[2]
        linesplan.T = linesplan.T * scaling_ratio[2]
        linesplan.B = linesplan.B * scaling_ratio[1]


