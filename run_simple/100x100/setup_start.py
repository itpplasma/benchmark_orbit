import numpy as np

nth = 100
nph = 100

s = 0.5
th = np.linspace(0, 2 * np.pi, 101)[:-1]
ph = np.linspace(0, 2 * np.pi, 101)[:-1]

vnorm = 1.0
vthnorm = 0.0

with open("start.dat", "w") as f:
    for i in range(nth):
        for j in range(nph):
            f.write(
                f"{s: .17f} {th[i]: .17f} {ph[j]: .17f} {vnorm: .17f} {vthnorm: .17f}\n"
            )
