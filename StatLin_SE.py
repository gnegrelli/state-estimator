

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
        print(row)
        # buses[str(int(row[0:4]))] = Bus(row)

# Create line objects
lines = dict()
line_set = datasets[1].split("\n")

for row in line_set:
    if row.strip():
        print(row)
        # lines[row[0:4].strip() + "-" + row[4:12].strip()] = Line(row)
