import scipy.stats
import numpy as np

scipy.stats.norm(0, 1)

std_dev = np.sqrt(25)

mean = 20

x1 = 10
x2 = -1

z1 = (x1 - mean)/std_dev

z2 = (x2 - mean)/std_dev

print("p(z < %.2f) = %.4f%%" % (z1, (scipy.stats.norm(mean, std_dev).cdf(x1))*100))
print("p(%.2f < z < %.2f) = %.4f%%" % (z2, z1, (scipy.stats.norm(mean, std_dev).cdf(x1) - scipy.stats.norm(mean, std_dev).cdf(x2))*100))

# print((scipy.stats.norm(0, 1).cdf(z1))*100)
# print((scipy.stats.norm(mean, std_dev).cdf(x1))*100)

# print((scipy.stats.norm(0, 1).cdf(z1) - scipy.stats.norm(0, 1).cdf(z2))*100)
# print((scipy.stats.norm(mean, std_dev).cdf(x1) - scipy.stats.norm(mean, std_dev).cdf(x2))*100)
