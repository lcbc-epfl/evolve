import numpy as np


# f1 = np.asarray([-2.312E-02, -2.401E-03, 2.386E-03]) * (627.5 / 0.5291772083)
# f2 = np.asarray([-2.358E-02, -2.702E-03, 2.227E-03]) * (627.5 / 0.5291772083)

f1 = np.asarray([-2.465E-02, -1.058E-02, 9.966E-05]) * (627.5 / 0.5291772083)
f2 = np.asarray([-2.089E-02, -3.969E-03, 1.346E-03]) * (627.5 / 0.5291772083)


print np.abs(np.linalg.norm(f1) - np.linalg.norm(f2))


print np.arccos(np.dot(f1 / np.linalg.norm(f1), f2 / np.linalg.norm(f2))) / np.pi
