import numpy as np


class Line:

    def __init__(self, dataline):

        self.origin = int(dataline[0:4].strip())
        self.destiny = int(dataline[8:12].strip())

        if dataline[17:23].strip():
            self.R = float(dataline[17:23])/100
        else:
            self.R = 0.

        if dataline[23:29].strip():
            self.X = float(dataline[23:29])/100
        else:
            self.X = 0.

        if dataline[29:35].strip():
            self.B = float(dataline[29:35])
        else:
            self.B = 0.

        # Power flowing from origin to destiny and destiny to origin
        self.S_od = 0.
        self.S_do = 0.

        # Measurement of active power flowing from origin to destiny bus
        self.flagPkm = int(dataline[80])
        if self.flagPkm:
            self.Pkm_m = float(dataline[81:88])
            self.sd_Pkm = float(dataline[88:94])

        # Measurement of reactive power flowing from origin to destiny bus
        self.flagQkm = int(dataline[95])
        if self.flagQkm:
            self.Qkm_m = float(dataline[96:103])
            self.sd_Qkm = float(dataline[103:109])

        # Measurement of active power flowing from destiny to origin bus
        self.flagPmk = int(dataline[116])
        if self.flagPmk:
            self.Pmk_m = float(dataline[117:124])
            self.sd_Pmk = float(dataline[124:130])

        # Measurement of reactive power flowing from destiny to origin bus
        self.flagQmk = int(dataline[131])
        if self.flagQmk:
            self.Qmk_m = float(dataline[132:139])
            self.sd_Qmk = float(dataline[139:145])

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

        self.ID = int(databus[:5])
        self.name = databus[10:20]
        self.bustype = bustypes[databus[7]]

        self.V = float(databus[22:26])/1000

        if databus[26:30].strip():
            self.theta = float(databus[26:30])
        else:
            self.theta = 0.

        # Read value of shunt on bus
        if databus[65:70].strip():
            self.shunt = float(databus[65:70])
        else:
            self.shunt = 0.

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

        # Read measurement of active power injected to bus
        self.flagP = int(databus[80])
        if self.flagP:
            self.P_m = float(databus[81:88])
            self.sd_P = float(databus[88:94])

        # Read measurement of reactive power injected to bus
        self.flagQ = int(databus[95])
        if self.flagQ:
            self.Q_m = float(databus[96:103])
            self.sd_Q = float(databus[103:109])

        # Read measurement of bus voltage
        self.flagV = int(databus[110])
        if self.flagV:
            self.V_m = float(databus[111:118])/1000
            self.sd_V = float(databus[118:124])

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


# Flag for WLS (True) or LS (False) state estimator
wls = True

# Flag for Fast Decoupled (True) or Traditional (False) Method
decoupled = True

# Importing data
# data_measures = open("Measurements2.txt", "r").read().split("\n")
datasets = open("2Barras.txt", "r").read().split("9999\n")

# Create bus objects
buses = dict()
bus_set = datasets[0].split('\n')

for row in bus_set[15:]:
    if row.strip():
        buses[str(int(row[0:4]))] = Bus(row)

# Create line objects
lines = dict()
line_set = datasets[1].split("\n")

for row in line_set[6:]:
    if row.strip():
        lines[row[0:4].strip() + "-" + row[8:12].strip()] = Line(row)

# Add measurements to objects
# for row in data_measures:
#     if row.strip() and row[0] is not '%':
#         if row.split("\t")[1] == '*':
#             buses[row.split("\t")[0]].add_measure(row.split("\t"))
#         else:
#             lines[row.split("\t")[0] + '-' + row.split("\t")[1]].add_measure(row.split("\t"))

# Flat Start
for bus in buses.values():
    bus.V = 1
    bus.theta = 0

# Nodal Admittance Matrix
Ybus = np.zeros((len(buses.keys()), len(buses.keys())), dtype=complex)

# Shunt Elements Vector
Bshunt = np.zeros(len(buses), dtype=complex)

for line in lines.values():
    Ybus[line.origin - 1][line.destiny - 1] = -1/(line.R + 1j*line.X)
    Bshunt[line.origin - 1] += 1j*line.B/2
    Bshunt[line.destiny - 1] += 1j*line.B/2

Ybus += Ybus.T
np.fill_diagonal(Ybus, Bshunt - np.sum(Ybus, axis=1))

# Tolerance
tolerance = 0.01

# Maximum absolute value of delta_x will be compared to tolerance. Delta_x is initiated with value higher than tolerance
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
        if line.flagPkm:
            z_p = np.hstack((z_p, np.array([line.Pkm_m])))
            w_p = np.hstack((w_p, np.array([line.sd_Pkm])))
            h_p = np.hstack((h_p, np.array([np.real(line.S_od)])))

            # Partial derivatives of Pkm on theta
            for bus in buses.values():
                if bus.bustype != 'Vθ':
                    if bus.ID == line.origin:
                        H_p = np.hstack((H_p, buses[str(line.origin)].V*buses[str(line.destiny)].V*(np.real(Ybus[line.origin - 1, line.destiny - 1])*np.sin(buses[str(line.origin)].theta - buses[str(line.destiny)].theta) - np.imag(Ybus[line.origin - 1, line.destiny - 1])*np.cos(buses[str(line.origin)].theta - buses[str(line.destiny)].theta))))
                    else:
                        H_p = np.hstack((H_p, buses[str(line.origin)].V*buses[str(line.destiny)].V*(-np.real(Ybus[line.origin - 1, line.destiny - 1])*np.sin(buses[str(line.origin)].theta - buses[str(line.destiny)].theta) + np.imag(Ybus[line.origin - 1, line.destiny - 1])*np.cos(buses[str(line.origin)].theta - buses[str(line.destiny)].theta))))

            # Partial derivatives of Pkm on V
            for bus in buses.values():
                if bus.ID == line.origin:
                    H_p = np.hstack((H_p, 2*buses[str(line.origin)].V*np.real(Ybus[line.origin - 1, line.destiny - 1]) - buses[str(line.destiny)].V*(np.real(Ybus[line.origin - 1, line.destiny - 1])*np.cos(buses[str(line.origin)].theta - buses[str(line.destiny)].theta) + np.imag(Ybus[line.origin - 1, line.destiny - 1])*np.sin(buses[str(line.origin)].theta - buses[str(line.destiny)].theta))))
                else:
                    H_p = np.hstack((H_p, -buses[str(line.origin)].V*(np.real(Ybus[line.origin - 1, line.destiny - 1])*np.cos(buses[str(line.origin)].theta - buses[str(line.destiny)].theta) + np.imag(Ybus[line.origin - 1, line.destiny - 1])*np.sin(buses[str(line.origin)].theta - buses[str(line.destiny)].theta))))

        if line.flagQkm:
            z_q = np.hstack((z_q, np.array([line.Qkm_m])))
            w_q = np.hstack((w_q, np.array([line.sd_Qkm])))
            h_q = np.hstack((h_q, np.array([np.imag(line.S_od)])))

            # Partial derivatives of Qkm on theta
            for bus in buses.values():
                if bus.bustype != 'Vθ':
                    if bus.ID == line.origin:
                        H_q = np.hstack((H_q, buses[str(line.origin)].V*buses[str(line.destiny)].V*(-np.real(Ybus[line.origin - 1, line.destiny - 1])*np.cos(buses[str(line.origin)].theta - buses[str(line.destiny)].theta) - np.imag(Ybus[line.origin - 1, line.destiny - 1])*np.sin(buses[str(line.origin)].theta - buses[str(line.destiny)].theta))))
                    else:
                        H_q = np.hstack((H_q, buses[str(line.origin)].V*buses[str(line.destiny)].V*(np.real(Ybus[line.origin - 1, line.destiny - 1])*np.cos(buses[str(line.origin)].theta - buses[str(line.destiny)].theta) + np.imag(Ybus[line.origin - 1, line.destiny - 1])*np.sin(buses[str(line.origin)].theta - buses[str(line.destiny)].theta))))

            # Partial derivatives of Qkm on V
            for bus in buses.values():
                if bus.ID == line.origin:
                    H_q = np.hstack((H_q, -2*buses[str(line.origin)].V*(np.imag(Ybus[line.origin - 1, line.destiny - 1]) + line.B/2) + buses[str(line.destiny)].V*(-np.real(Ybus[line.origin - 1, line.destiny - 1])*np.sin(buses[str(line.origin)].theta - buses[str(line.destiny)].theta) + np.imag(Ybus[line.origin - 1, line.destiny - 1])*np.cos(buses[str(line.origin)].theta - buses[str(line.destiny)].theta))))

                else:
                    H_q = np.hstack((H_q, buses[str(line.origin)].V*(-np.real(Ybus[line.origin - 1, line.destiny - 1])*np.sin(buses[str(line.origin)].theta - buses[str(line.destiny)].theta) + np.imag(Ybus[line.origin - 1, line.destiny - 1])*np.cos(buses[str(line.origin)].theta - buses[str(line.destiny)].theta))))

        if line.flagPmk:
            z_p = np.hstack((z_p, np.array([line.Pmk_m])))
            w_p = np.hstack((w_p, np.array([line.sd_Pmk])))
            h_p = np.hstack((h_p, np.array([np.real(line.S_do)])))

            # Partial derivatives of Pmk on theta
            for bus in buses.values():
                if bus.bustype != 'Vθ':
                    if bus.ID == line.destiny:
                        H_p = np.hstack((H_p, buses[str(line.destiny)].V*buses[str(line.origin)].V*(np.real(Ybus[line.destiny - 1, line.origin - 1])*np.sin(buses[str(line.destiny)].theta - buses[str(line.origin)].theta) - np.imag(Ybus[line.destiny - 1, line.origin - 1])*np.cos(buses[str(line.destiny)].theta - buses[str(line.origin)].theta))))
                    else:
                        H_p = np.hstack((H_p, buses[str(line.destiny)].V*buses[str(line.origin)].V*(-np.real(Ybus[line.destiny - 1, line.origin - 1])*np.sin(buses[str(line.destiny)].theta - buses[str(line.origin)].theta) + np.imag(Ybus[line.destiny - 1, line.origin - 1])*np.cos(buses[str(line.destiny)].theta - buses[str(line.origin)].theta))))

            # Partial derivatives of Pmk on V
            for bus in buses.values():
                if bus.ID == line.destiny:
                    H_p = np.hstack((H_p, 2*buses[str(line.destiny)].V*np.real(Ybus[line.destiny - 1, line.origin - 1]) - buses[str(line.origin)].V*(np.real(Ybus[line.destiny - 1, line.origin - 1])*np.cos(buses[str(line.destiny)].theta - buses[str(line.origin)].theta) + np.imag(Ybus[line.destiny - 1, line.origin - 1])*np.sin(buses[str(line.destiny)].theta - buses[str(line.origin)].theta))))
                else:
                    H_p = np.hstack((H_p, -buses[str(line.destiny)].V*(np.real(Ybus[line.destiny - 1, line.origin - 1])*np.cos(buses[str(line.destiny)].theta - buses[str(line.origin)].theta) + np.imag(Ybus[line.destiny - 1, line.origin - 1])*np.sin(buses[str(line.destiny)].theta - buses[str(line.origin)].theta))))

        if line.flagQmk:
            z_q = np.hstack((z_q, np.array([line.Qmk_m])))
            w_q = np.hstack((w_q, np.array([line.sd_Qmk])))
            h_q = np.hstack((h_q, np.array([np.imag(line.S_do)])))

            # Partial derivatives of Qmk on theta
            for bus in buses.values():
                if bus.bustype != 'Vθ':
                    if bus.ID == line.destiny:
                        H_q = np.hstack((H_q, buses[str(line.destiny)].V*buses[str(line.origin)].V*(-np.real(Ybus[line.destiny - 1, line.origin - 1])*np.cos(buses[str(line.destiny)].theta - buses[str(line.origin)].theta) - np.imag(Ybus[line.destiny - 1, line.origin - 1])*np.sin(buses[str(line.destiny)].theta - buses[str(line.origin)].theta))))
                    else:
                        H_q = np.hstack((H_q, buses[str(line.destiny)].V*buses[str(line.origin)].V*(np.real(Ybus[line.destiny - 1, line.origin - 1])*np.cos(buses[str(line.destiny)].theta - buses[str(line.origin)].theta) + np.imag(Ybus[line.destiny - 1, line.origin - 1])*np.sin(buses[str(line.destiny)].theta - buses[str(line.origin)].theta))))

            # Partial derivatives of Qmk on V
            for bus in buses.values():
                if bus.ID == line.destiny:
                    H_q = np.hstack((H_q, -2*buses[str(line.destiny)].V*(np.imag(Ybus[line.destiny - 1, line.origin - 1]) + line.B/2) + buses[str(line.origin)].V*(-np.real(Ybus[line.destiny - 1, line.origin - 1])*np.sin(buses[str(line.destiny)].theta - buses[str(line.origin)].theta) + np.imag(Ybus[line.destiny - 1, line.origin - 1])*np.cos(buses[str(line.destiny)].theta - buses[str(line.origin)].theta))))
                else:
                    H_q = np.hstack((H_q, buses[str(line.destiny)].V*(-np.real(Ybus[line.destiny - 1, line.origin - 1])*np.sin(buses[str(line.destiny)].theta - buses[str(line.origin)].theta) + np.imag(Ybus[line.destiny - 1, line.origin - 1])*np.cos(buses[str(line.destiny)].theta - buses[str(line.origin)].theta))))

    for bus in buses.values():
        if bus.flagP:
            z_p = np.hstack((z_p, np.array([bus.P_m])))
            w_p = np.hstack((w_p, np.array([bus.sd_P])))
            h_p = np.hstack((h_p, np.array([np.real(bus.P)])))

        if bus.flagQ:
            z_q = np.hstack((z_q, np.array([bus.Q_m])))
            w_q = np.hstack((w_q, np.array([bus.sd_Q])))
            h_q = np.hstack((h_q, np.array([np.real(bus.Q)])))

            # I should calculate H_p and H_q here. I don't want though

        if bus.flagV:
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
