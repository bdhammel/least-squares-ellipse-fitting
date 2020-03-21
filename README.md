[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3723294.svg)](https://doi.org/10.5281/zenodo.3723294)
[![bdhammel](https://circleci.com/gh/bdhammel/least-squares-ellipse-fitting.svg?style=shield)](https://app.circleci.com/pipelines/github/bdhammel/least-squares-ellipse-fitting)


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

X1, X2 = el.make_test_ellipse()

X = np.array(list(zip(X1, X2)))
reg = LsqEllipse().fit(X)
center, width, height, phi = reg.as_parameters()

plt.close('all')
fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)
ax.axis('equal')
ax.plot(X1, X2, 'ro', label='test data', zorder=1)

ellipse = Ellipse(xy=center, width=2*width, height=2*height, angle=np.rad2deg(phi),
               edgecolor='b', fc='None', lw=2, label='Fit', zorder = 2)
ax.add_patch(ellipse)

plt.legend()
plt.show()
```

![ellipse fit](./media/ellipse_fit.png)


**Cite this work**
```
@software{ben_hammel_2020_3723294,
  author       = {Ben Hammel and Nick Sullivan-Molina},
  title        = {bdhammel/least-squares-ellipse-fitting: v2.0.0},
  month        = mar,
  year         = 2020,
  publisher    = {Zenodo},
  version      = {v2.0.0},
  doi          = {10.5281/zenodo.3723294},
  url          = {https://doi.org/10.5281/zenodo.3723294}
}
```
Ben Hammel, & Nick Sullivan-Molina. (2020, March 21). bdhammel/least-squares-ellipse-fitting: v2.0.0 (Version v2.0.0). Zenodo. http://doi.org/10.5281/zenodo.3723294
