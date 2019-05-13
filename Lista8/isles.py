import numpy as np
import scipy.linalg

z = np.array([0, 0, 0, 0])

H = np.array([[1., 0., -1., 0., 0.],
              [-1., 0., 1., 0., 0.],
              [2., -1., -1., 0., 0.],
              [0., -1., -1., 3., -1.]])

G = np.dot(H.T, H)

G_f = np.around(scipy.linalg.lu(G)[2], 2)

print(G_f)

counter = 0
for row in range(5):
    if np.diag(G_f)[row] != 0:
        G_f[row] = G_f[row]/np.diag(G_f)[row]
    elif row != 4:
        G_f[row, row] = 1
        aux = np.zeros(5)
        aux[row] = 1
        H = np.vstack((H, aux))
        z = np.hstack((z, counter*np.ones(1)))
        counter += 1

print(G_f)

print(H)
print(z)

print(np.dot(H.T, H))
