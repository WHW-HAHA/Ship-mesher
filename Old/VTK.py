class VTK():
    def __init__(self, X, Y, Z, order, coor):
        self.X = X
        self.Y = Y
        self.Z = Z
        self.order = order
        self.num_panel = len(self.X) / 4
        self.num_cor = len(coor)
        self.coor = coor
        self.info = 'info'

    def Pat_write(self):
        with open(r'C:\Users\907932\Desktop\out_put.vtk', 'w') as outfile:
            # header
            outfile.write('# vtk DataFile Version 2.0\n')
            outfile.write('An case related object\n')
            outfile.write('ASCII\n')
            outfile.write('\n')
            outfile.write('DATASET UNSTRUCTURED_GRID\n')
            outfile.write('POINTS' + ' ' + str(self.num_cor) + ' ' + 'float\n')

            # vertex
            for vertex in self.coor:
                outfile.write(str(vertex[0]) + '\t' + str(vertex[1]) + '\t' + str(vertex[2]) + '\n')
            outfile.write('\n')

            # Cells
            outfile.write('CELLS' + ' ' + str(self.num_panel) + ' ' + str(5 * int(self.num_panel)) + '\n')

            outfile.write('25       0       0       1       0       0       0       0       0\n')
            outfile.write('ISYM=0 nbod= 1%f\n' % math.floor(self.num_panel))
            outfile.write(
                '26       0       0       1%8.0f%8.0f       0       0       0\n' % (self.num_cor, self.num_panel))
            outfile.write(self.info + '\n')
            '''
            for item in self.coor:
                outfile.write('%f       %f       %f       %d\n'%(item[0], item[1], item[2], item[3] ))
            '''
            # first part, write the coordinate of points
            for item in self.coor:
                outfile.write(' 1       %d       0       2       0       0       0       0       0\n' % (item[3]))
                outfile.write('  %f   %f    %f\n' % (item[0], item[1], item[2]))
                outfile.write('1G       6       0       0  000000\n')

            # second part, write the connectivity
            Count = 0  # start with 0
            i = 0
            Num = 0
            while i < len(self.X):
                if None in [self.order[i], self.order[i + 1], self.order[i + 2], self.order[i + 3]]:
                    num = 3
                else:
                    num = 4
                outfile.write(' 2       %d       %d       2       0       0       0       0       0\n' % (Count, num))
                outfile.write(
                    '       %d       0       0       0 0.000000000E+00 0.000000000E+00 0.000000000E+00\n' % (num))
                for item in [self.order[i], self.order[i + 1], self.order[i + 2], self.order[i + 3]]:
                    if item != None:
                        outfile.write('       %d' % item)
                outfile.write('\n')
                i += 4
                Count += 1
                Num += num
            print(Num)
            outfile.write('99       0       0       1       0       0       0       0       0')
            outfile.close()
