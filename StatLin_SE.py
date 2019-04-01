

class Measurement:

    def __init__(self, meas_data):
        self.origin = meas_data[0]
        self.destiny = meas_data[1]
        self.value = float(meas_data[2])
        self.std_dev = float(meas_data[3])


class Line:

    def __init__(self, dataline):

        self.origin = int(dataline[0:4].strip())
        self.destiny = int(dataline[4:12].strip())

        if dataline[16:23].strip():
            self.R = float(dataline[16:23])
        else:
            self.R = 0.

        if dataline[23:29].strip():
            self.X = float(dataline[23:29])
        else:
            self.X = 0.

        if dataline[29:35].strip():
            self.B = float(dataline[29:35])
        else:
            self.B = 0.


class Bus:

    def __init__(self, databus):

        bustypes = {'0': 'PQ', '1': 'PV', '2': 'Vθ'}

        self.ID = int(databus[0:4])
        self.name = databus[8:22]
        self.bustype = bustypes[databus[4:8].strip()]

        self.V = float(databus[22:26]) / 1000

        if databus[26:30].strip():
            self.theta = float(databus[26:30])
        else:
            self.theta = 0.

        p = []
        q = []
        for item in [databus[30:35], databus[56:60]]:
            if not item.strip():
                p.append(0.)
            else:
                p.append(float(item))

        for item in [databus[35:40], databus[60:65]]:
            if not item.strip():
                q.append(0.)
            else:
                q.append(float(item))

        self.P = (p[0] - p[1])
        self.Q = (q[0] - q[1])


# Importing data
data_measures = open("Measurements.txt", "r").read().split("\n")
datasets = open("System.txt", "r").read().split("9999\n")

# Create measurement objects
measures = dict()


for row in data_measures:
    if row.strip() and row[0] is not '%':
        measures[row.split("\t")[0] + '_' + row.split("\t")[1]] = Measurement(row.split("\t"))


# Create bus objects
buses = dict()
bus_set = datasets[0].split('\n')

for row in bus_set:
    if row.strip():
        buses[str(int(row[0:4]))] = Bus(row)

# Create line objects
lines = dict()
line_set = datasets[1].split("\n")

for row in line_set:
    if row.strip():
        lines[row[0:4].strip() + "-" + row[4:12].strip()] = Line(row)
