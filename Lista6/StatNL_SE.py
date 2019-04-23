import numpy as np


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

        # Power flowing from origin to destiny and destiny to origin
        self.S_od = 0
        self.S_do = 0

        # Measurements of Power
        self.P_m = 0
        self.sd_P = 0
        self.Q_m = 0
        self.sd_Q = 0

    def add_measure(self, data):

        if not self.P_m:
            self.P_m = float(data[2])
            self.sd_P = float(data[3])
        else:
            self.Q_m = float(data[2])
            self.sd_Q = float(data[3])

    # Function to save values of power flowing through line
    def save_flow(self, p, q, bus_origin):
        if bus_origin == self.origin:
            self.S_od = p + 1j*q
        elif bus_origin == self.destiny:
            self.S_do = p + 1j*q
        else:
            print("ERROR: The bus must be at one of the ends of the line")


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

        self.P_m = 0
        self.sd_P = 0
        self.Q_m = 0
        self.sd_Q = 0
        self.V_m = 0
        self.sd_V = 0

    def add_measure(self, data):

        if data[4] == 'd':
            self.V_m = float(data[2])
            self.sd_V = float(data[3])
        elif not self.P:
            self.P_m = float(data[2])
            self.sd_P = float(data[3])
        else:
            self.Q_m = float(data[2])
            self.sd_Q = float(data[3])

    # Function to save values of power
    def save_power(self, p, q):

        self.P = p
        self.Q = q


# Flag for WLS or LS state estimator
wls = True

# Importing data
data_measures = open("Measurements2.txt", "r").read().split("\n")
datasets = open("System2.txt", "r").read().split("9999\n")

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

# Add measurements to objects
for row in data_measures:
    if row.strip() and row[0] is not '%':
        if row.split("\t")[1] == '*':
            buses[row.split("\t")[0]].add_measure(row.split("\t"))
        else:
            lines[row.split("\t")[0] + '-' + row.split("\t")[1]].add_measure(row.split("\t"))

# Flat Start
for bus in buses.values():
    bus.V = 1
    bus.theta = 0

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


# Calculation of power injected to each bus
for bus in buses.values():
    p, q = 0, 0

    for otherbus in buses.values():

        # Calculate angle difference
        theta_km = bus.theta - otherbus.theta

        # Calculate active and reactive power reaching bus
        p += bus.V*otherbus.V*(np.real(Ybus[bus.ID - 1, otherbus.ID - 1])*np.cos(theta_km) + np.imag(Ybus[bus.ID - 1, otherbus.ID - 1]*np.sin(theta_km)))
        q += bus.V*otherbus.V*(np.real(Ybus[bus.ID - 1, otherbus.ID - 1])*np.sin(theta_km) - np.imag(Ybus[bus.ID - 1, otherbus.ID - 1]*np.cos(theta_km)))

    # Save calculated values of power
    bus.save_power(p, q)


# Calculation of power flowing through lines
for line in lines.values():

    theta_km = buses[str(line.origin)].theta - buses[str(line.destiny)].theta

    Vk = buses[str(line.origin)].V
    Vm = buses[str(line.destiny)].V

    Y = 1/(line.R + 1j*line.X)

    pkm = (Vk**2)*np.real(Y) - Vk*Vm*(np.real(Y)*np.cos(theta_km) + np.imag(Y)*np.sin(theta_km))
    qkm = -(Vk**2)*(np.imag(Y) + line.B/2) - Vk*Vm*(np.real(Y)*np.sin(theta_km) - np.imag(Y)*np.cos(theta_km))

    line.save_flow(pkm, qkm, line.origin)

    pmk = (Vm**2)*np.real(Y) - Vm*Vk*(np.real(Y)*np.cos(-theta_km) + np.imag(Y)*np.sin(-theta_km))
    qmk = -(Vm**2)*(np.imag(Y) + line.B/2) - Vk*Vm*(np.real(Y)*np.sin(-theta_km) - np.imag(Y)*np.cos(-theta_km))

    line.save_flow(pmk, qmk, line.destiny)


h_p = np.array([])
h_q = np.array([])
h_v = np.array([])

z_p = np.array([])
z_q = np.array([])
z_v = np.array([])

w_p = np.array([])
w_q = np.array([])
w_v = np.array([])

H_p = np.array([])
H_q = np.array([])
H_v = np.array([])

for key in lines.keys():
    if lines[key].P_m is not 0 and lines[key].Q_m is not 0:
        z_p = np.hstack((z_p, np.array([lines[key].P_m])))
        w_p = np.hstack((w_p, np.array([lines[key].sd_P])))
        h_p = np.hstack((h_p, np.array([np.real(lines[key].S_od)])))

        z_q = np.hstack((z_q, np.array([lines[key].Q_m])))
        w_q = np.hstack((w_q, np.array([lines[key].sd_Q])))
        h_q = np.hstack((h_q, np.array([np.imag(lines[key].S_od)])))

for key in buses.keys():
    if buses[key].P_m is not 0 and buses[key].Q_m is not 0:
        z_p = np.hstack((z_p, np.array([buses[key].P_m])))
        w_p = np.hstack((w_p, np.array([buses[key].sd_P])))
        h_p = np.hstack((h_p, np.array([np.real(buses[key].P)])))

        z_q = np.hstack((z_q, np.array([buses[key].Q_m])))
        w_q = np.hstack((w_q, np.array([buses[key].sd_Q])))
        h_q = np.hstack((h_q, np.array([np.real(buses[key].Q)])))

    if buses[key].V_m is not 0:
        z_v = np.hstack((z_v, np.array([buses[key].V_m])))
        w_v = np.hstack((w_v, np.array([buses[key].sd_V])))
        h_v = np.hstack((h_v, np.array([buses[key].V])))

# print(z_p.reshape((len(z_p), 1)))

z = np.hstack((z_p, z_q, z_v))

h = np.hstack((h_p, h_q, h_v))

W = np.zeros((len(z), len(z)))
np.fill_diagonal(W, np.hstack((w_p, w_q, w_v)))

print(z)

print(h)

'''
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
        z[m] += -r[m]/(np.diag(W)[m]*np.diag(Omega)[m])
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
'''