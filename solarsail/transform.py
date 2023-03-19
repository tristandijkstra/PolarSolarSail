import math
import numpy as np

def simplifiedSailToInertial(inclination, argPeriapsis, trueAnomaly, RAAN):
    cosArgTrue = math.cos(argPeriapsis + trueAnomaly)
    sinArgTrue = math.sin(argPeriapsis + trueAnomaly)
    cosRAAN = math.cos(RAAN)
    sinRAAN = math.sin(RAAN)
    cosIncl = math.cos(inclination)
    sinIncl = math.sin(inclination)

    res = np.array(
        [
            [
                (
                    cosArgTrue * cosRAAN * cosIncl * cosIncl
                    - sinArgTrue * sinRAAN * cosIncl
                    + cosArgTrue * cosRAAN * sinIncl * sinIncl
                )
                / (
                    cosArgTrue * cosArgTrue * cosRAAN * cosRAAN * cosIncl * cosIncl
                    + cosArgTrue * cosArgTrue * cosRAAN * cosRAAN * sinIncl * sinIncl
                    + cosArgTrue * cosArgTrue * sinRAAN * sinRAAN * cosIncl * cosIncl
                    + cosArgTrue * cosArgTrue * sinRAAN * sinRAAN * sinIncl * sinIncl
                    + sinArgTrue * sinArgTrue * cosRAAN * cosRAAN * cosIncl * cosIncl
                    + sinArgTrue * sinArgTrue * cosRAAN * cosRAAN * sinIncl * sinIncl
                    + sinArgTrue * sinArgTrue * sinRAAN * sinRAAN * cosIncl * cosIncl
                    + sinArgTrue * sinArgTrue * sinRAAN * sinRAAN * sinIncl * sinIncl
                ),
                -(
                    sinArgTrue * cosRAAN * cosIncl * cosIncl
                    + cosArgTrue * sinRAAN * cosIncl
                    + sinArgTrue * cosRAAN * sinIncl * sinIncl
                )
                / (
                    cosArgTrue * cosArgTrue * cosRAAN * cosRAAN * cosIncl * cosIncl
                    + cosArgTrue * cosArgTrue * cosRAAN * cosRAAN * sinIncl * sinIncl
                    + cosArgTrue * cosArgTrue * sinRAAN * sinRAAN * cosIncl * cosIncl
                    + cosArgTrue * cosArgTrue * sinRAAN * sinRAAN * sinIncl * sinIncl
                    + sinArgTrue * sinArgTrue * cosRAAN * cosRAAN * cosIncl * cosIncl
                    + sinArgTrue * sinArgTrue * cosRAAN * cosRAAN * sinIncl * sinIncl
                    + sinArgTrue * sinArgTrue * sinRAAN * sinRAAN * cosIncl * cosIncl
                    + sinArgTrue * sinArgTrue * sinRAAN * sinRAAN * sinIncl * sinIncl
                ),
                (sinRAAN * sinIncl)
                / (
                    cosRAAN * cosRAAN * cosIncl * cosIncl
                    + cosRAAN * cosRAAN * sinIncl * sinIncl
                    + sinRAAN * sinRAAN * cosIncl * cosIncl
                    + sinRAAN * sinRAAN * sinIncl * sinIncl
                ),
            ],
            [
                (
                    cosArgTrue * sinRAAN * cosIncl * cosIncl
                    + sinArgTrue * cosRAAN * cosIncl
                    + cosArgTrue * sinRAAN * sinIncl * sinIncl
                )
                / (
                    cosArgTrue * cosArgTrue * cosRAAN * cosRAAN * cosIncl * cosIncl
                    + cosArgTrue * cosArgTrue * cosRAAN * cosRAAN * sinIncl * sinIncl
                    + cosArgTrue * cosArgTrue * sinRAAN * sinRAAN * cosIncl * cosIncl
                    + cosArgTrue * cosArgTrue * sinRAAN * sinRAAN * sinIncl * sinIncl
                    + sinArgTrue * sinArgTrue * cosRAAN * cosRAAN * cosIncl * cosIncl
                    + sinArgTrue * sinArgTrue * cosRAAN * cosRAAN * sinIncl * sinIncl
                    + sinArgTrue * sinArgTrue * sinRAAN * sinRAAN * cosIncl * cosIncl
                    + sinArgTrue * sinArgTrue * sinRAAN * sinRAAN * sinIncl * sinIncl
                ),
                -(
                    sinArgTrue * sinRAAN * cosIncl * cosIncl
                    - cosArgTrue * cosRAAN * cosIncl
                    + sinArgTrue * sinRAAN * sinIncl * sinIncl
                )
                / (
                    cosArgTrue * cosArgTrue * cosRAAN * cosRAAN * cosIncl * cosIncl
                    + cosArgTrue * cosArgTrue * cosRAAN * cosRAAN * sinIncl * sinIncl
                    + cosArgTrue * cosArgTrue * sinRAAN * sinRAAN * cosIncl * cosIncl
                    + cosArgTrue * cosArgTrue * sinRAAN * sinRAAN * sinIncl * sinIncl
                    + sinArgTrue * sinArgTrue * cosRAAN * cosRAAN * cosIncl * cosIncl
                    + sinArgTrue * sinArgTrue * cosRAAN * cosRAAN * sinIncl * sinIncl
                    + sinArgTrue * sinArgTrue * sinRAAN * sinRAAN * cosIncl * cosIncl
                    + sinArgTrue * sinArgTrue * sinRAAN * sinRAAN * sinIncl * sinIncl
                ),
                -(cosRAAN * sinIncl)
                / (
                    cosRAAN * cosRAAN * cosIncl * cosIncl
                    + cosRAAN * cosRAAN * sinIncl * sinIncl
                    + sinRAAN * sinRAAN * cosIncl * cosIncl
                    + sinRAAN * sinRAAN * sinIncl * sinIncl
                ),
            ],
            [
                (sinArgTrue * sinIncl)
                / (
                    cosArgTrue * cosArgTrue * cosIncl * cosIncl
                    + cosArgTrue * cosArgTrue * sinIncl * sinIncl
                    + sinArgTrue * sinArgTrue * cosIncl * cosIncl
                    + sinArgTrue * sinArgTrue * sinIncl * sinIncl
                ),
                (cosArgTrue * sinIncl)
                / (
                    cosArgTrue * cosArgTrue * cosIncl * cosIncl
                    + cosArgTrue * cosArgTrue * sinIncl * sinIncl
                    + sinArgTrue * sinArgTrue * cosIncl * cosIncl
                    + sinArgTrue * sinArgTrue * sinIncl * sinIncl
                ),
                cosIncl / (cosIncl * cosIncl + sinIncl * sinIncl),
            ],
        ]
    )

    return res