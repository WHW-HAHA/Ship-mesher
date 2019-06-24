# import rhinoscriptsyntax as rs


class ship():

    def ImportPoints(self):
        # filter = "Text file (*.txt)|*.txt|All Files (*.*)|*.*||"
        # filename = rs.OpenFileName("Open Point File", filter)
        # if not filename: return

        file = open(r'C:\Users\907932\Desktop\DMA_Work_Flow\DIFFRAC_work_flow\PAT_related\Python_pat\Ship.txt', "r")
        contents = file.readlines()
        X = []
        Y = []
        Z = []
        X_crossection = {}
        Y_crossection = {}
        Z_crossection = {}
        tempx = None
        i = 1
        name = 'crossection' + '1'

        for line in contents:
            items = line.strip('()\n').split(',')
            x = float(items[0])
            y = float(items[1])
            z = float(items[2])
            if tempx == None:
                tempx = x

            if x == tempx:
                # in the same crossection
                X.append(x)
                Y.append(y)
                Z.append(z)
            else:
                X_crossection[name] = X
                Y_crossection[name] = Y
                Z_crossection[name] = Z
                i += 1
                name = 'crossection' + str(i)
                X = []
                Y = []
                Z = []
                tempx = x
        file.close()
        self.X_crossection = X_crossection
        self.Y_crossection = Y_crossection
        self.Z_crossection = Z_crossection
        self.L = max(self.X_crossection)
        self.D = max(self.Z_crossection)

    def other_lines(self):
        point_fore_deck = [self.L,0, self.D]
        point_aft_deck = [0, 0 ,self.D ]
        point_fore_keel = [self.L, 0, 0]
        point_aft_keel = [0,0,0]
        self.nutural_line1 = [point_fore_deck, point_aft_deck]
        self.nutural_line2 = [point_fore_keel, point_aft_keel]
        self.contour_line = []
        for key in self.X_crossection:
            self.contour_line.append([self.X_crossection[key][-1], self.Y_crossection[key][-1], self.Z_crossection[key][-1]])


if __name__ == '__main__':
    ship = ship()
    contents = ship.ImportPoints()
    list = []
    for key in ship.X_crossection:
        for (x, y, z) in zip(ship.X_crossection[key], ship.Y_crossection[key], ship.Z_crossection[key]):
            # print("{}, {}, {}".format(x, y, z))
            list.append([x,y,z])
        print(list)
        print('***********************')
        list= []

    ship.other_lines()
    print(ship.contour_line)



    # contents = [ship.point_from_string(line) for line in contents]
    # rs.AddPoints(contents)
    # rs.AddPoints([0,0,ship.Z])
    # rs.AddPoints([ship.X, 0, ship.Z])
    # rs.AddLine([0,0,ship.Z], [ship.X, 0, ship.Z])


