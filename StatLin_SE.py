

class Measurement:

    def __init__(self, meas_data):
        self.origin = meas_data[0]
        self.destiny = meas_data[1]
        self.value = float(meas_data[2])
        self.std_dev = float(meas_data[3])


class Line:

    def __init__(self, line_data):
        self.origin = 1
        self.destiny = 1
        self.R = 1
        self.X = 1
        self.B = 1


class Bus:

    def __init__(self, bus_data):
        self.ID = 1
        self.bustype = 1
        self.V = 1
        self.theta = 1


measures = {}

data_measures = open("Measurements.txt", "r").read().split("\n")

for each_line in data_measures:
    if each_line and each_line[0] is not '%':
        measures[each_line.split("\t")[0] + '_' + each_line.split("\t")[1]] = Measurement(each_line.split("\t"))


for key in measures.keys():
    print(measures[key].value)
