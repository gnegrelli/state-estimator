import numpy as np

H = np.array([[-1, 0, 1, 0, 0],
              [1, 0, -1, 0, 0],
              [-2, 1, 1, 0, 0],
              [0, 1, 1, -3, 1]])

H = np.array([[1, -1, 0, 0, 0, 0],
              [0, 1, -1, 0, 0, 0],
              [0, 0, 0, 1, -1, 0],
              [0, 0, 0, -1, 1, 0],
              [0, 0, 0, 1, 0, -1],
              [2, -1, -1, 0, 0, 0],
              [-1, -1, 3, -1, 0, 0],
              [0, 0, 0, -1, 1, 0],
              [0, 0, 0, -1, 0, 1]])

H = H.T

factors = np.zeros((H.shape[0], H.shape[0] - 1))

# Gaussian decomposition adding pseudo-measurements
for i in range(H.shape[0] - 1):

    # Check if diagonal is null
    if H[i, i] == 0:

        # If the entire row is null add a pseudo-measurement
        if not H[i].any():
            zeros = np.zeros((H.shape[0], 1))
            zeros[i] = 1
            H = np.hstack((H, zeros))

        # Swap columns and continue decomposition
        col = np.where(H[i] != 0)[0][0]
        aux = H[:, i]
        H[:, i] = H[:, col]
        H[:, col] = aux

    else:
        H[i] = H[i]/H[i, i]

    for j in range(i + 1, H.shape[0]):
        H[j] = -H[j, i]/H[i, i]*H[i] + H[j]

print(H)
print(factors)
