"""Implements Direction Cosine Matrix (DCM) attitude parameter."""

from typing import Optional

import numpy as np

import attitude.parameter.base as base
import attitude.parameter.euler_angle as euler_angle


class DCM(base.BaseAttitudeParameter):
    """Direction Cosine Matrix (DCM)."""
    def __init__(
        self, array_3x3: Optional[float] = None) -> None:
        if array_3x3 is None:
            self.array = np.eye(3)
        else:
            if not np.shape(array_3x3) == (3, 3):
                raise ValueError("DCM must be initialized with a 3x3 array of floats")
            # TODO: add validation that given array represents a valid DCM
            # (i.e. that columns are orthonormal with determinant +1)
            self.array = np.array(array_3x3)

    @classmethod
    def from_321_euler_angles(
        cls, ea_321: "euler_angle.EulerAngle321"
    ) -> None:
        """
        Initialize the DCM from a 3-2-1 set of Euler
        angles, where the angles are expressed in radians.
        3-2-1: a1*yaw, a2*pitch, a3*roll.
        """
        a1, a2, a3 = ea_321.array
        array = np.array([
            [
                np.cos(a2)*np.cos(a1),
                np.cos(a2)*np.sin(a1),
                -np.sin(a2)
            ],
            [
                np.sin(a3)*np.sin(a2)*np.cos(a1) - np.cos(a3)*np.sin(a1),
                np.sin(a3)*np.sin(a2)*np.sin(a1) + np.cos(a3)*np.cos(a1),
                np.sin(a3)*np.cos(a2)
            ],
            [
                np.cos(a3)*np.sin(a2)*np.cos(a1) + np.sin(a3)*np.sin(a1),
                np.cos(a3)*np.sin(a2)*np.sin(a1) - np.sin(a3)*np.cos(a1),
                np.cos(a3)*np.cos(a2)
            ]
        ])
        return cls(array)

    @classmethod
    def from_313_euler_angles(
        cls, ea_313: "euler_angle.EulerAngle313"
    ) -> "DCM":
        """
        Initialize the DCM from a 3-1-3 set of Euler
        angles, where the angles are expressed in radians.
        3-1-3: a1*yaw, a2*roll, a3*yaw.
        """
        a1, a2, a3 = ea_313.array
        array = np.array([
            [
                np.cos(a3)*np.cos(a1) - np.sin(a3)*np.cos(a2)*np.sin(a1),
                np.cos(a3)*np.sin(a1) + np.sin(a3)*np.cos(a2)*np.cos(a1),
                np.sin(a3)*np.sin(a2)
            ],
            [
                -np.sin(a3)*np.cos(a1) - np.cos(a3)*np.cos(a2)*np.sin(a1),
                -np.sin(a3)*np.sin(a1) + np.cos(a3)*np.cos(a2)*np.cos(a1),
                np.cos(a3)*np.sin(a2)
            ],
            [
                np.sin(a2)*np.sin(a1),
                -np.sin(a2)*np.cos(a1),
                np.cos(a2)
            ]
        ])
        return cls(array)
