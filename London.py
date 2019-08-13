import numpy as np

# H = np.array([[1, 0, -1, 0, 0],
#               [-1, 0, 1, 0, 0],
#               [2, -1, -1, 0, 0],
#               [0, -1, -1, 3, -1]], dtype=float)

# H = np.array([[1, -1, 0, 0, 0, 0],
#               [1, 0, -1, 0, 0, 0],
#               [-1, 1, 0, 0, 0, 0],
#               [0, -1, 1, 0, 0, 0],
#               [0, 0, 0, -1, 1, 0],
#               [0, 0, 0, -1, 0, 1],
#               [2, -1, -1, 0, 0, 0],
#               [0, 0, 0, -1, 1, 0],
#               [0, 0, 0, -1, 0, 1]], dtype=float)

# H = np.array([[1, -1, 0, 0, 0, 0],
#               [0, 1, -1, 0, 0, 0],
#               [0, 0, 0, 1, -1, 0],
#               [0, 0, 0, -1, 1, 0],
#               [0, 0, 0, 1, 0, -1],
#               [2, -1, -1, 0, 0, 0],
#               [-1, -1, 3, -1, 0, 0],
#               [0, 0, 0, -1, 1, 0],
#               [0, 0, 0, -1, 0, 1]], dtype=float)

# H = np.array([[5, -3, -2, 0, 0, 0],
#               [-3, 2, 1, 0, 0, 0],
#               [-2, 1, 2, -3, 1, 1],
#               [0, 0, -3, 10, -4, -3],
#               [0, 0, 1, -4, 2, 1],
#               [0, 0, 1, -3, 1, 1]], dtype=float)

H = np.array([[1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
              [-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
              [1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
              [-1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
              [0, 1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
              [0, -1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
              [0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
              [0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
              [0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0],
              [0, 0, 0, -1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0],
              [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -1, 0],
              [0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 1, 0],
              [0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1],
              [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0],
              [0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0],
              [0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 1, 0, 0],
              [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0],
              [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0],
              [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1],
              [2, -1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
              [0, -1, -1, 5, -1, 0, -1, 0, -1, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 2, -1, 0],
              [0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, -1, 3, -1]], dtype=float)

H = H.T

# Nametag list
measures = ["P12", "P21", "P15", "P51", "P25", "P52", "P34", "P43", "P47", "P74", "P611", "P613", "P136", "P78", "P194", "P1011", "P1110", "P126", "P1213", "P1312", "P1413", "P1", "P4", "P8", "P12", "P13"]

factors = np.zeros((H.shape[0], H.shape[0] - 1))

# List of pseudomeasurements
pseudomeas = [np.array([[0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0]]).T,
              np.array([[0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1]]).T,
              np.array([[0, 0, 0, 0, 0, -1, 0, 0, 0, -1, 2, 0, 0, 0]]).T,
              np.array([[0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]]).T,
              np.array([[0, 0, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0, 0]]).T]
added_meas = 0

# pseudomeas = []

# Gaussian decomposition adding pseudo-measurements
for i in range(H.shape[0] - 1):

    # Check if diagonal is null
    if H[i, i] == 0:

        # While the entire row is null add a pseudo-measurement from list
        while not H[i].any():
            try:
                add_col = pseudomeas.pop(0)
                added_meas += 1
                print("Adding pseudo-measurement #%d" % added_meas)
            except IndexError:
                print("All pseudo-measurements available were added. System is not observable.")
                break
            for j in range(i):
                add_col[j] = np.dot(np.hstack((factors, add_col))[j], add_col)
            H = np.hstack((H, add_col))

        # Swap columns and continue decomposition
        try:
            col = np.where(H[i] != 0)[0][0]
        except IndexError:
            break
        aux = np.array(H[:, i])
        H[:, i] = H[:, col]
        H[:, col] = aux

        # Swap columns on nametag list
        aux = measures[i]
        measures[i] = measures[col]
        measures[col] = aux

    factors[i, i] = 1/H[i, i]
    H[i] = H[i]/H[i, i]

    for j in range(i + 1, H.shape[0]):
        factors[j, i] = -H[j, i]/H[i, i]
        H[j] = -H[j, i]/H[i, i]*H[i] + H[j]

# Backward processing
for i in range(H.shape[0] - 2, 0, -1):
    for j in range(i - 1, -1, -1):
        factors[j, i] = -H[j, i]/H[i, i]
        H[j] = -H[j, i]/H[i, i]*H[i] + H[j]

print(measures)
print(H)
print(factors)
