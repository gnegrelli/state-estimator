import numpy as np


class Measurement:

    def __init__(self, meas_data):

        self.origin = int(meas_data[0])
        if meas_data[1] == "*":
            self.destiny = meas_data[1]
        else:
            self.destiny = int(meas_data[1])

        self.value = float(meas_data[2])
        self.std_dev = float(meas_data[3])

        self.type = meas_data[4]


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

        self.V = float(databus[22:26])/1000

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


# Flag for WLS or LS state estimator
wls = True

# Importing data
data_measures = open("Measurements.txt", "r").read().split("\n")
datasets = open("System.txt", "r").read().split("9999\n")

# Create measurement objects
measures = dict()

for row in data_measures:
    if row.strip() and row[0] is not '%':
        measures[row.split("\t")[0] + '-' + row.split("\t")[1]] = Measurement(row.split("\t"))

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

# Nodal Admitance Matrix
Ybus = np.zeros((len(buses.keys()), len(buses.keys())), dtype=complex)

# Shunt Elements Vector
Bshunt = np.zeros(len(buses), dtype=complex)

for key in lines.keys():
    Ybus[lines[key].origin - 1][lines[key].destiny - 1] = -1/(lines[key].R + 1j*lines[key].X)
    Bshunt[lines[key].origin - 1] += 1j*lines[key].B/2
    Bshunt[lines[key].destiny - 1] += 1j*lines[key].B/2

Ybus += Ybus.T

np.fill_diagonal(Ybus, Bshunt - np.sum(Ybus, axis=1))

# Create Jacobian Matrix
H = np.zeros((len(measures.keys()), len(buses.keys())))

# Create measurement vector
z = np.zeros((len(measures.keys()), 1))

# Create weight matrix
W = np.zeros((len(measures.keys()), len(measures.keys())))

# Fill Jacobian matrix, measurement vector and weight matrix
aux = 0
for key in measures.keys():

    # Fill Jacobian matrix
    if measures[key].type == 'u':
        if key.split("-")[1] is not "*":
            H[aux, measures[key].origin - 1] = 1/np.imag(1/-Ybus[measures[key].origin - 1, measures[key].destiny - 1])
            H[aux, measures[key].destiny - 1] = -1/np.imag(1/-Ybus[measures[key].origin - 1, measures[key].destiny - 1])
        else:
            for lkey in lines.keys():
                if measures[key].origin == lines[lkey].origin:
                    H[aux, measures[key].origin - 1] += 1/lines[lkey].X
                    H[aux, lines[lkey].destiny - 1] += -1 / lines[lkey].X
                elif measures[key].origin == lines[lkey].destiny:
                    H[aux, measures[key].origin - 1] += 1 / lines[lkey].X
                    H[aux, lines[lkey].origin - 1] += -1 / lines[lkey].X
    else:
        H[aux, measures[key].origin - 1] = 1

    # Fill measurement vector
    z[aux] = measures[key].value

    # Fill weight matrix
    if wls:
        W[aux, aux] = 1/measures[key].std_dev**2

    aux += 1

# Exclude column for Vθ buses
for key in sorted(buses.keys(), reverse=True):
    if buses[key].bustype == 'Vθ':
        H = np.delete(H, int(key)-1, 1)

# Calculate states
if not wls:
    x_hat = np.linalg.solve(np.dot(H.T, H), np.dot(H.T, z))
    print("Least-Square State Estimator")
else:
    H_til = np.dot(np.sqrt(W), H)
    z_til = np.dot(np.sqrt(W), z)

    x_hat = np.linalg.solve(np.dot(H_til.T, H_til), np.dot(H_til.T, z_til))
    print("Weighted Least-Square State Estimator")

print("\nEstimated States:")
for item in x_hat:
    print("%.6f" % item[0])

# Residue calculation
r = z - np.dot(H, x_hat)

# Gain Matrix
G = np.dot(H.T, np.dot(W, H))

# Covariance matrix of residues
Omega = np.linalg.inv(W) - np.dot(H, np.linalg.solve(G, H.T))

# Normalized residue calculation
r_n = np.abs(r.T/np.sqrt(np.diag(Omega))).T

print("\nResidue:")
for item in r:
    print("%.5f" % item[0])

print("\nNormalized Residue:")
for item in r_n:
    if item == max(r_n):
        print('\033[31m' + "%.5f" % item[0] + '\033[0m')
    else:
        print("%.5f" % item[0])

# Treating Gross Error in measurements

# Flags for eliminate or recuperate measures with error
eliminate = not True
recuperate = True

# Counter of number of iterations
counter = 0

while max(r_n) > 3 and counter < len(z):

    counter += 1

    if counter % 10 == 1:
        print("\n%dst Iteration" % counter)
    elif counter % 10 == 2:
        print("\n%dnd Iteration" % counter)
    elif counter % 10 == 3:
        print("\n%drd Iteration" % counter)
    else:
        print("\n%dth Iteration" % counter)
    print(50*"-")

    # Obtain position of maximum normalized residue
    m = r_n.argmax()

    if eliminate:
        # Delete row of z and column of H correspondent to measurement with error
        print("\nEliminating measure with highest normalized residue")
        z = np.delete(z, m, axis=0)
        H = np.delete(H, m, axis=0)
        W = np.delete(np.delete(W, m, axis=1), m, axis=0)

        # Recalculate Gain and Covariance of Residues Matrices
        G = np.dot(H.T, np.dot(W, H))
        Omega = np.linalg.inv(W) - np.dot(H, np.linalg.solve(G, H.T))

    elif recuperate:
        # Recuperate measurement with error
        print("\nRecuperating measure with highest normalized residue")
        delta_z = np.ones_like(z)
        delta_z[m] = - r[m]/(np.diag(W)[m]*np.diag(Omega)[m])
        z += delta_z
        H = H
        W = W

    # Recalculate states
    if not wls:
        x_hat = np.linalg.solve(np.dot(H.T, H), np.dot(H.T, z))
    else:
        H_til = np.dot(np.sqrt(W), H)
        z_til = np.dot(np.sqrt(W), z)
        x_hat = np.linalg.solve(G, np.dot(H_til.T, z_til))

    print("\nNew Measurements:")
    for item in z:
        print("%.6f" % item[0])

    print("\nNew Estimated States:")
    for item in x_hat:
        print("%.6f" % item[0])

    # Residue recalculation
    r = z - np.dot(H, x_hat)

    # Normalized residue recalculation
    r_n = np.abs(r.T/np.sqrt(np.diag(Omega))).T

    print("\nNew Residue:")
    for item in r:
        print("%.5f" % item[0])

    print("\nNew Normalized Residue:")
    for item in r_n:
        if item == max(r_n):
            print('\033[31m' + "%.5f" % item[0] + '\033[0m')
        else:
            print("%.5f" % item[0])
