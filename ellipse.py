import numpy as np
import numpy.linalg as la

__version__ = '3.0.1'


class LsqEllipse:
    """Lest Squares fitting of Elliptical data

    Attributes
    ----------
    coef_ : array
        Estimated coefficients for the Least squares fit to the elliptical data
        containing the values [a,b,c,d,f,g].T corresponding to
        ax**2 + 2bxy + cy**2 + 2dx + 2fy + g

    References
    ----------
    (*) Halir R., Flusser J. 'Numerically Stable Direct Least Squares
    Fitting of Ellipses'
    (**) Weisstein, Eric W. "Ellipse." From MathWorld--A Wolfram Web Resource.
    http://mathworld.wolfram.com/Ellipse.html

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.linear_model import LsqEllipse
    >>> x = np.array([ 1.,  0., -1., -0.,  1.])
    >>> y = np.array([ 0. ,  0.5,  0. , -0.5, -0. ])
    >>> X = np.c_[x, y]
    >>> reg = LsqEllipse().fit(X)
    >>> reg.as_parameters()
    """
    ALLOWED_FEATURES = 2

    def __init__(self):
        self.coef_ = None

    def _check_data(self, X):

        n_samples, n_features = X.shape
        if not n_features == self.ALLOWED_FEATURES:
            raise ValueError("Incorrect number of features. "
                             f"Got {n_features} features, expected 2. ")

        if n_samples < 5:
            raise ValueError("Received too few samples"
                             f"Got {n_samples} features, 5 or more required. ")

        return X

    def fit(self, X):
        """Fit the data

        Parameters
        ----------
        X : array, shape (n_points, 2)
            Data values for the x-y data pairs to fit

        Returns
        -------
        self : returns an instance of self.
        """
        X = self._check_data(X)

        # extract x-y pairs
        x, y = X.T

        # Quadratic part of design matrix [eqn. 15] from (*)
        
        D1 = np.vstack([x*x, x * y, y*y]).T     
        
        # Linear part of design matrix [eqn. 16] from (*)
        D2 = np.vstack([x, y, np.ones_like(x)]).T

        # Forming scatter matrix [eqn. 17] from (*)
        S1 = D1.T @ D1
        S2 = D1.T @ D2
        S3 = D2.T @ D2

        # Constraint matrix [eqn. 18]
        C1 = np.array([[0., 0., 2.], [0., -1., 0.], [2., 0., 0.]])

        # Reduced scatter matrix [eqn. 29]
        M = la.inv(C1) @ (S1 - S2 @ la.inv(S3) @ S2.T)

        # M*|a b c >=l|a b c >. Find eigenvalues and eigenvectors from this
        # equation [eqn. 28]
        eigval, eigvec = np.linalg.eig(M)

        # Eigenvector must meet constraint 4ac - b^2 to be valid.
        cond = (
            4*np.multiply(eigvec[0, :], eigvec[2, :])
            - np.power(eigvec[1, :], 2)
        )
        a1 = eigvec[:, np.nonzero(cond > 0)[0]]

        # |d f g> = -S3^(-1) * S2^(T)*|a b c> [eqn. 24]
        a2 = la.inv(-S3) @ S2.T @ a1

        # Eigenvectors |a b c d f g>
        # list of the coefficients describing an ellipse [a,b,c,d,f,g]
        # corresponding to ax**2 + 2bxy + cy**2 + 2dx + 2fy + g
        if (a1[0] < 0):                                                     #change the signum in case < 0
            a1 = -a1
            a2 = -a2
        self.coef_ = np.vstack([a1, a2])

        return self

    @property
    def coefficients(self):
        """
        List of the coefficients describing the fitted ellipse

        Returns
        -------
        [a,b,c,d,f,g] corresponding to ax**2 + 2bxy + cy**2 + 2dx + 2fy + g
        """
        return np.asarray(self.coef_).ravel()

    def as_parameters(self):
        """Returns the definition of the fitted ellipse as localized parameters

        Returns
        _______
        center : list
            [x0, y0]
        width : float
            Semimajor axis
        height : float
            Semiminor axis
        phi : float
            The counterclockwise angle of rotation from the x-axis to the major
            axis of the ellipse
        """

        # Eigenvectors are the coefficients of an ellipse in general form
        # a*x^2 + 2*b*x*y + c*y^2 + 2*d*x + 2*f*y + g = 0
        # [eqn. 15) from (**) or (***)
        a = self.coefficients[0]
        b = self.coefficients[1] / 2.0
        c = self.coefficients[2]
        d = self.coefficients[3] / 2.0
        f = self.coefficients[4] / 2.0
        g = self.coefficients[5]
        
        # Finding center of ellipse [eqn.19 and 20] from (**)
        x0 = (c*d - b*f) / (b*b - a*c)          
        y0 = (a*f - b*d) / (b*b - a*c)          
        center = [x0, y0]

        # Find the semi-axes lengths [eqn. 21 and 22] from (**)
        numerator = 2.0 * (a*f*f + c*d*d + g*b*b - 2.0*b*d*f - a*c*g)   
        
        denominator1 = (b*b - a*c) * (np.sqrt((a - c)*(a - c) + 4.0*b*b) - (a + c))                 
        
        denominator2 = (b*b - a*c) * (-np.sqrt((a - c)*(a - c) + 4.0*b*b) - (a + c))              
        
        width  = np.sqrt(numerator / denominator1)
        height = np.sqrt(numerator / denominator2)

        # Angle of counterclockwise rotation of major-axis of ellipse to x-axis
        # [eqn. 23] from (**) or [eqn. 26] from (***).

        if (b == 0 and a < c):                  #calculate with the b == 0 condition
            phi = 0.0

        if (b == 0 and a > c):
            phi = 0.5 * np.pi

        if (b != 0 and a != c):                 #calculate with the a > c condition
            if (a < c):
                phi = 0.5 * np.arctan((2.0 * b) / (a - c))
            if (a > c):
                phi = 0.5 * np.pi + 0.5 * np.arctan((2.0 * b) / (a - c))

        if (a == c):                            # add the case of a circle, "undefined" is the right value
            phi = 0.0
      

        return center, width, height, phi
        
    # Add a function that shows the ellipse coefficients
    def ellipse_coeff(self):
        """
        List of the coefficients describing the fitted ellipse

        Returns
        -------
        [a,b,c,d,f,g] corresponding to ax**2 + 2bxy + cy**2 + 2dx + 2fy + g
        """
        a = self.coefficients[0]
        b = self.coefficients[1] / 2.0
        c = self.coefficients[2]
        d = self.coefficients[3] / 2.0
        f = self.coefficients[4] / 2.0
        g = self.coefficients[5]
        
        return a, b, c, d, f, g
        
        
    def return_fit(self, n_points=None, t=None):
        """Return the X, Y values of the predicted ellipse

        Points are returned along the parametric curve of the ellipse as evenly
        spaced points starting at t=0 to t=2pi

        Parameters
        ---------
        n_points : int
            Number of points to return
        t : array
            Parametric points used to generate x-y pairs, If provided,
            `n_points` will be ignored

        Returns
        -------
        X : array, shape (n_points, 2)
            data values for the x-y data pairs
        """
        if self.coef_ is None:
            raise ValueError("Must call .fit() before using .return_fit()")

        if n_points is None and t is None:
            raise AttributeError("A value for `n_points` or `t` must be ",
                                 "provided")

        if t is None:
            t = np.linspace(0, 2 * np.pi, n_points)

        center, width, height, phi = self.as_parameters()

        x = (center[0] 
             + width * np.cos(t) * np.cos(phi) 
             - height * np.sin(t) * np.sin(phi))
        y = (center[1] 
             + width * np.cos(t) * np.sin(phi) 
             + height * np.sin(t) * np.cos(phi))

        return np.c_[x, y]
