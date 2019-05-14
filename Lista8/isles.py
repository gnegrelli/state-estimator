import numpy as np
import scipy.linalg

z = np.array([[0, 0, 0, 0]])

H = np.array([[1., 0., -1., 0., 0.],
              [-1., 0., 1., 0., 0.],
              [2., -1., -1., 0., 0.],
              [0., -1., -1., 3., -1.]])

G = np.dot(H.T, H)
G_f = np.zeros_like(G)

P, L, U = scipy.linalg.lu(G)

null_pivot = 0
for row in range(len(U)):
    if np.diag(U)[row] != 0:
        G_f[row] = U[row]/np.diag(U)[row]
    else:
        G_f[row, row] = 1

        aux = np.zeros(5)
        aux[row] = 1
        H = np.vstack((H, aux))

        z = np.hstack((z, np.array([[null_pivot]])))
        null_pivot += 1

print(G_f)

print(H)
print(z)

print(np.dot(H.T, H))
