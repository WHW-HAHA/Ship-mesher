import math
class Patran():
    def __init__(self, X, Y, Z, order, coor):
        self.X = X
        self.Y = Y
        self.Z = Z
        self.order = order
        self.num_panel = len(self.X)/4
        self.num_cor = len(coor)
        self.coor = coor
        self.info = 'info'


    def Pat_write(self):
        path = r'C:\Users\907932\Desktop\out_put.pat'
        with open(path, 'w') as output:
            output.write('25       0       0       1       0       0       0       0       0\n')
            output.write('ISYM=0 nbod= 1%f\n'%math.floor(self.num_panel))
            output.write('26       0       0       1%8.0f%8.0f       0       0       0\n'%(self.num_cor, self.num_panel))
            output.write(self.info + '\n')
            '''
            for item in self.coor:
                output.write('%f       %f       %f       %d\n'%(item[0], item[1], item[2], item[3] ))
            '''
            # first part, write the coordinate of points
            for item in self.coor:
                output.write(' 1       %d       0       2       0       0       0       0       0\n'%(item[3]))
                output.write('  %f   %f    %f\n'%(item[0], item[1], item[2]))
                output.write('1G       6       0       0  000000\n')

            # second part, write the connectivity
            Count = 0 # start with 0
            i = 0
            Num = 0
            while i < len(self.X):
                if None in [self.order[i], self.order[i+1], self.order[i+2], self.order[i+3]]:
                    num = 3
                else:
                    num = 4
                output.write(' 2       %d       %d       2       0       0       0       0       0\n' %(Count,num))
                output.write('       %d       0       0       0 0.000000000E+00 0.000000000E+00 0.000000000E+00\n'%(num))
                for item in [self.order[i], self.order[i+1], self.order[i+2], self.order[i+3]]:
                    if item != None:
                        output.write('%8.0f' % item)
                output.write('\n')
                i  += 4
                Count += 1
                Num += num
            print(Num+self.num_panel)
            self.Allpoints = Num+self.num_panel
            output.write('99       0       0       1       0       0       0       0       0')
            output.close()







