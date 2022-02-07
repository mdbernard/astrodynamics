"""Implementation of Euler Angles."""

import numpy as np

import attitude.parameter.base as base
import attitude.parameter.constants as constants
import attitude.parameter.direction_cosine_matrix as direction_cosine_matrix


class BaseEulerAngle(base.BaseAttitudeParameter):
    """
    Base class for Euler Angle parameterizations.
    All angles are handled in radians by default except in
    methods explicitly denoted as working with degrees.
    """

    def __init__(
        self, a1: float = None, a2: float = None, a3: float = None
    ) -> None:
        """
        Initialize the Euler Angles based on the three angles (rad)
        about the sequence of axes (e.g. 3-2-1, 3-1-3, or others).
        """
        self.initialized = False
        self.a1 = self.a2 = self.a3 = None
        if all([a is not None for a in (a1, a2, a3)]):
            self.a1 = a1
            self.a2 = a2
            self.a3 = a3
            self.initialized = True

    def __repr__(self):
        return f"{self.a1:.4f}, {self.a2:.4f}, {self.a3:.4f}"

    @classmethod
    def from_dcm(
        cls, dcm: "direction_cosine_matrix.DCM"
    ) -> "BaseEulerAngle":
        """Initialize the Euler Angles instance from a DCM."""
        raise NotImplementedError("Derived classes should define this method")

    @property
    def rate_matrix(self) -> np.ndarray:
        """
        The rate matrix for the differential kinematic equation
        d(EA_321)/dt = rate_matrix @ angular_velocity_vector.
        """
        raise NotImplementedError("Derived classes should define this property")

    @property
    def array(self):
        return np.array([self.a1, self.a2, self.a3])

    @property
    def is_singular(self):
        """Whether the set of values is singular."""
        raise NotImplementedError("Derived classes should define this property")

    def print_rad(self):
        print(self)

    def print_deg(self):
        _a1, _a2, _a3 = self.array
        self.a1 = np.rad2deg(self.a1)
        self.a2 = np.rad2deg(self.a2)
        self.a3 = np.rad2deg(self.a3)
        print(self)
        self.a1 = _a1
        self.a2 = _a2
        self.a3 = _a3


class BaseSymmetricEulerAngle(BaseEulerAngle):
    @property
    def is_singular(self):
        """
        Symmetric Euler Angles are singular when the second angle
        is 0 or pi radians.
        """
        return (
            abs(self.a2) < constants.SINGULARITY_EPS
            or abs(self.a2 - constants.PI) < constants.SINGULARITY_EPS
        )


class BaseAsymmetricEulerAngle(BaseEulerAngle):
    @property
    def is_singular(self):
        """
        Asymmetric Euler Angles are singular when the second angle
        is +/- pi/2 radians.
        """
        return abs(abs(self.a2) - constants.PI/2) < constants.SINGULARITY_EPS


class EulerAngle321(BaseAsymmetricEulerAngle):
    """3-2-1 Euler Angle (yaw, pitch, roll)."""

    @property
    def rate_matrix(self) -> np.ndarray:
        return (1/np.cos(self.a2))*np.array([
            [
                0,
                np.sin(self.a3),
                np.cos(self.a3)
            ],
            [
                0,
                np.cos(self.a2)*np.cos(self.a3),
                -np.cos(self.a2)*np.sin(self.a3)
            ],
            [
                np.cos(self.a2),
                np.sin(self.a2)*np.sin(self.a3),
                np.sin(self.a2)*np.cos(self.a3)
            ]
        ])

    @classmethod
    def from_dcm(
        cls, dcm: "direction_cosine_matrix.DCM"
    ) -> "EulerAngle321":
        C = dcm.array
        a1 = np.arctan2(C[0,1], C[0,0])
        a2 = -np.arcsin(C[0,2])
        a3 = np.arctan2(C[1,2], C[2,2])
        return cls(a1, a2, a3)


class EulerAngle313(BaseSymmetricEulerAngle):
    """3-1-3 Euler Angle (yaw, roll, yaw)."""

    def from_dcm(
        cls, dcm: "direction_cosine_matrix.DCM"
    ) -> "EulerAngle313":
        C = dcm.array
        a1 = np.arctan2(C[2,0], -C[2,1])
        a2 = np.arcnp.cos(C[2,2])
        a3 = np.arctan2(C[0,2], C[1,2])
        return cls(a1, a2, a3)
