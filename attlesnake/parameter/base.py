"""Base classes for attitude parameterizations."""

import numpy as np


class BaseAttitudeParameter:
    """Base class for attitude parameterizations."""

    def __init__(self, array: np.ndarray) -> None:
        """
        Abstract initialization based on numpy array of values.
        Derived classes should define how to initialize based on
        these values.
        """
