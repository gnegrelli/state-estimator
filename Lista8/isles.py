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

        aux = np.zeros(5)
        aux[row] = 1
        H = np.vstack((H, aux))

        z = np.hstack((z, np.array([[null_pivot]])))
        null_pivot += 1

print(G_f)

print(H)
print(z)

print(np.dot(H.T, H))


G = np.dot(H.T, H)
G_f = np.zeros_like(G)

P, L, U = scipy.linalg.lu(G)

null_pivot = 0
for row in range(len(U)):

    if np.diag(U)[row] != 0:

        G_f[row] = U[row]/np.diag(U)[row]

    else:

        aux = np.zeros(5)
        aux[row] = 1
        H = np.vstack((H, aux))

        z = np.hstack((z, np.array([[null_pivot]])))
        null_pivot += 1

print(np.around(G_f, 2))

theta = np.around(np.linalg.solve(G, np.dot(H.T, z.T)), 2)

for b1 in range(5):
    for b2 in range(b1 + 1, 5):
        print("P%d%d = %.4f" % (b1 + 1, b2 + 1, theta[b1] - theta[b2]))
