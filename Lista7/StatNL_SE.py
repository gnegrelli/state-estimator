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

    # Function to refresh values of voltage and angle on buses
    def refresh(self, v, ang):

        self.theta = ang
        self.V = v


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

for line in lines.values():
    Ybus[line.origin - 1][line.destiny - 1] = -1/(line.R + 1j*line.X)
    Bshunt[line.origin - 1] += 1j*line.B/2
    Bshunt[line.destiny - 1] += 1j*line.B/2

np.fill_diagonal(Ybus, Bshunt - np.sum(Ybus, axis=1))

# Tolerance
tolerance = 0.01

# Maximum absolute value of delta_x will be compared to tolerance. Delta_x is initiated bigger than tolerance
delta_x = np.array([tolerance + 1, tolerance + 1])

# Iteration counter
counter = 0

while max(abs(delta_x)) > tolerance and counter < 5:

    print("\nIteration #%d" % counter)
    for bus in buses.values():
        print("V%d: %.6f < %.4f" % (bus.ID, bus.V, bus.theta))

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

        Y = Ybus[line.origin - 1, line.destiny - 1]

        pkm = (Vk**2)*np.real(Y) - Vk*Vm*(np.real(Y)*np.cos(theta_km) + np.imag(Y)*np.sin(theta_km))
        qkm = -(Vk**2)*(np.imag(Y) + line.B/2) - Vk*Vm*(np.real(Y)*np.sin(theta_km) - np.imag(Y)*np.cos(theta_km))

        line.save_flow(pkm, qkm, line.origin)

        pmk = (Vm**2)*np.real(Y) - Vm*Vk*(np.real(Y)*np.cos(-theta_km) + np.imag(Y)*np.sin(-theta_km))
        qmk = -(Vm**2)*(np.imag(Y) + line.B/2) - Vk*Vm*(np.real(Y)*np.sin(-theta_km) - np.imag(Y)*np.cos(-theta_km))

        line.save_flow(pmk, qmk, line.destiny)

    # Submatrices of vector h (Calculated values of measurements)
    h_p = np.array([])
    h_q = np.array([])
    h_v = np.array([])

    # Submatrices of vector z (Measured values)
    z_p = np.array([])
    z_q = np.array([])
    z_v = np.array([])

    # Submatrices of vector w (Diagonal of Covariance Matrix)
    w_p = np.array([])
    w_q = np.array([])
    w_v = np.array([])

    # Submatrices of matrix H (Jacobian Matrix)
    H_p = np.array([])
    H_q = np.array([])
    H_v = np.hstack((np.zeros((len(buses), len(buses) - 1)), np.eye(len(buses))))

    for line in lines.values():
        if line.P_m is not 0 and line.Q_m is not 0:
            z_p = np.hstack((z_p, np.array([line.P_m])))
            w_p = np.hstack((w_p, np.array([line.sd_P])))
            h_p = np.hstack((h_p, np.array([np.real(line.S_od)])))

            z_q = np.hstack((z_q, np.array([line.Q_m])))
            w_q = np.hstack((w_q, np.array([line.sd_Q])))
            h_q = np.hstack((h_q, np.array([np.imag(line.S_od)])))

            Vk = buses[str(line.origin)].V
            Vm = buses[str(line.destiny)].V

            # theta_km =

            # Partial derivatives on theta
            for bus in buses.values():
                if bus.bustype != 'Vθ':
                    if bus.ID == line.origin:
                        H_p = np.hstack((H_p, buses[str(line.origin)].V*buses[str(line.destiny)].V*(np.real(Ybus[line.origin - 1, line.destiny - 1])*np.sin(buses[str(line.origin)].theta - buses[str(line.destiny)].theta) - np.imag(Ybus[line.origin - 1, line.destiny - 1])*np.cos(buses[str(line.origin)].theta - buses[str(line.destiny)].theta))))
                        H_q = np.hstack((H_q, buses[str(line.origin)].V*buses[str(line.destiny)].V*(-np.real(Ybus[line.origin - 1, line.destiny - 1])*np.cos(buses[str(line.origin)].theta - buses[str(line.destiny)].theta) - np.imag(Ybus[line.origin - 1, line.destiny - 1])*np.sin(buses[str(line.origin)].theta - buses[str(line.destiny)].theta))))

                    else:
                        H_p = np.hstack((H_p, buses[str(line.origin)].V*buses[str(line.destiny)].V*(-np.real(Ybus[line.origin - 1, line.destiny - 1])*np.sin(buses[str(line.origin)].theta - buses[str(line.destiny)].theta) + np.imag(Ybus[line.origin - 1, line.destiny - 1])*np.cos(buses[str(line.origin)].theta - buses[str(line.destiny)].theta))))
                        H_q = np.hstack((H_q, buses[str(line.origin)].V*buses[str(line.destiny)].V*(np.real(Ybus[line.origin - 1, line.destiny - 1])*np.cos(buses[str(line.origin)].theta - buses[str(line.destiny)].theta) + np.imag(Ybus[line.origin - 1, line.destiny - 1])*np.sin(buses[str(line.origin)].theta - buses[str(line.destiny)].theta))))

            # Partial derivatives on V
            for bus in buses.values():
                if bus.ID == line.origin:
                    H_p = np.hstack((H_p, 2*buses[str(line.origin)].V*np.real(Ybus[line.origin - 1, line.destiny - 1]) - buses[str(line.destiny)].V*(np.real(Ybus[line.origin - 1, line.destiny - 1])*np.cos(buses[str(line.origin)].theta - buses[str(line.destiny)].theta) + np.imag(Ybus[line.origin - 1, line.destiny - 1])*np.sin(buses[str(line.origin)].theta - buses[str(line.destiny)].theta))))
                    H_q = np.hstack((H_q, -2*buses[str(line.origin)].V*(np.imag(Ybus[line.origin - 1, line.destiny - 1]) + line.B/2) + buses[str(line.destiny)].V*(-np.real(Ybus[line.origin - 1, line.destiny - 1])*np.sin(buses[str(line.origin)].theta - buses[str(line.destiny)].theta) + np.imag(Ybus[line.origin - 1, line.destiny - 1])*np.cos(buses[str(line.origin)].theta - buses[str(line.destiny)].theta))))

                else:
                    H_p = np.hstack((H_p, -buses[str(line.origin)].V*(np.real(Ybus[line.origin - 1, line.destiny - 1])*np.cos(buses[str(line.origin)].theta - buses[str(line.destiny)].theta) + np.imag(Ybus[line.origin - 1, line.destiny - 1])*np.sin(buses[str(line.origin)].theta - buses[str(line.destiny)].theta))))
                    H_q = np.hstack((H_q, buses[str(line.origin)].V*(-np.real(Ybus[line.origin - 1, line.destiny - 1])*np.sin(buses[str(line.origin)].theta - buses[str(line.destiny)].theta) + np.imag(Ybus[line.origin - 1, line.destiny - 1])*np.cos(buses[str(line.origin)].theta - buses[str(line.destiny)].theta))))

    for bus in buses.values():
        if bus.P_m is not 0 and bus.Q_m is not 0:
            z_p = np.hstack((z_p, np.array([bus.P_m])))
            w_p = np.hstack((w_p, np.array([bus.sd_P])))
            h_p = np.hstack((h_p, np.array([np.real(bus.P)])))

            z_q = np.hstack((z_q, np.array([bus.Q_m])))
            w_q = np.hstack((w_q, np.array([bus.sd_Q])))
            h_q = np.hstack((h_q, np.array([np.real(bus.Q)])))

            # I should calculate H_p and H_q here. I don't want though

        if bus.V_m is not 0:
            z_v = np.hstack((z_v, np.array([bus.V_m])))
            w_v = np.hstack((w_v, np.array([bus.sd_V])))
            h_v = np.hstack((h_v, np.array([bus.V])))

    # Assembling submatrices

    # Measurements vector
    z = np.hstack((z_p, z_q, z_v))

    # Calculated values of measurements
    h = np.hstack((h_p, h_q, h_v))

    # Covariance Matrix
    W = np.zeros((len(z), len(z)))
    np.fill_diagonal(W, np.hstack((w_p, w_q, w_v)))

    # Jacobian Matrix
    H = np.hstack((H_p, H_q, H_v.flatten())).reshape((len(z), 2*len(buses)-1))

    # Submatrices of states
    x_theta = np.array([])
    x_V = np.array([])

    for bus in buses.values():
        if bus.bustype != 'Vθ':
            x_theta = np.hstack((x_theta, bus.theta))

        x_V = np.hstack((x_V, bus.V))

    # Assembling submatrices to state matrix
    x = np.hstack((x_theta, x_V))

    # Gain Matrix
    G = np.dot(H.T, np.dot(W, H))

    # Update State Values
    delta_x = np.linalg.solve(G, np.dot(H.T, np.dot(W, (z-h).T)))
    x += delta_x

    vtheta = 0
    for bus in buses.values():
        if bus.bustype == 'Vθ':
            bus.refresh(x[len(buses) + bus.ID - 2], 0)
            vtheta += 1
        else:
            bus.refresh(x[len(buses) + bus.ID - 2], x[bus.ID - vtheta - 1])

    counter += 1

print("\n" + 30*"@" + "\n")
print("Final States")

for bus in buses.values():
    print("V%d: %.6f < %.4f"% (bus.ID, bus.V, bus.theta))