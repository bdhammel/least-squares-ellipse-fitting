

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2578663.svg)](https://doi.org/10.5281/zenodo.2578663)


# Least Squares fitting of ellipses, python routine 

based on the  publication 
[Halir, R., Flusser, J.: 'Numerically Stable Direct Least Squares 
            Fitting of Ellipses'](./media/WSCG98.pdf)

## Example execution

```python
import ellipses as el
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

data = el.make_test_ellipse()

lsqe = el.LSqEllipse()
lsqe.fit(data)
center, width, height, phi = lsqe.parameters()

plt.close('all')
fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)
ax.axis('equal')
ax.plot(data[0], data[1], 'ro', label='test data', zorder=1)

ellipse = Ellipse(xy=center, width=2*width, height=2*height, angle=np.rad2deg(phi),
               edgecolor='b', fc='None', lw=2, label='Fit', zorder = 2)
ax.add_patch(ellipse)

plt.legend()
plt.show()
```

![ellipse fit](./media/ellipse_fit.png)


**Cite this work**
```
@misc{https://doi.org/10.5281/zenodo.2578663,
  doi = {10.5281/zenodo.2578663},
  url = {https://zenodo.org/record/2578663},
  author = {Hammel,  Ben and Sullivan-Molina,  Nick},
  title = {bdhammel/least-squares-ellipse-fitting: Initial release},
  publisher = {Zenodo},
  year = {2019}
}
```
B. Hammel and N. Sullivan-Molina. bdhammel/least-squares-ellipse-fitting: Initial release, 2019.
