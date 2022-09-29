import numpy as np
import enum

_FOA = float | np.ndarray[np.float64]

# ----------------------------------------------------------------------------------------
# Core Two-Body Functions
# ----------------------------------------------------------------------------------------
def semi_minor_axis(sma, ecc):
    """Compute the semi-minor axis of a conic orbit."""
    return sma * np.sqrt(1.0 - ecc**2)


def semi_latus_rectum(sma, ecc):
    """Compute the semi-latus rectum (parameter) of a conic orbit."""
    return sma * (1 - ecc**2)


def periapsis_radius(sma, ecc):
    return sma * (1.0 - ecc)


def apoapsis_radius(sma, ecc):
    return sma * (1.0 + ecc)


def apsis_radii(sma, ecc):
    ae = sma * ecc
    return (sma - ae, sma + ae)


def circular_orbit_speed(gm, sma):
    return np.sqrt(gm / sma)


def two_body_radius(sma, ecc, ta):
    return semi_latus_rectum(sma, ecc) / (1 + ecc * np.cos(ta))


def specific_angular_momentum(gm, sma, ecc):
    return np.sqrt(gm * sma * (1.0 - ecc**2))


def areal_velocity(gm, sma, ecc):
    return specific_angular_momentum(gm, sma, ecc) / 2.0


# ----------------------------------------------------------------------------------------
# Keplerian Element Abstraction
# ----------------------------------------------------------------------------------------
class KeplerianElements:
    def __init__(self, sma, *, ecc=0.0, inc=0.0, aop=0.0, raan=0.0, ta=0.0):
        def _as_np_array(v, *, dtype=np.float64):
            return np.asarray(v, dtype=dtype).reshape(1, -1)[0, :]

        self._sma = _as_np_array(sma)
        self._ecc = _as_np_array(ecc)
        self._inc = _as_np_array(inc)
        self._aop = _as_np_array(aop)
        self._raan = _as_np_array(raan)
        self._ta = _as_np_array(ta)

        self.count = np.max(
            np.array(
                [
                    len(self._sma),
                    len(self._ecc),
                    len(self._inc),
                    len(self._aop),
                    len(self._raan),
                    len(self._ta),
                ]
            )
        )

        if np.any(self._ecc < 0.0):
            raise ValueError("negative eccentricity")

        if ((self._sma > 0.0) & self.is_hyperbolic()).any():
            raise ValueError("positive semi-major axis for hyperbolic eccentricity")

        if ((self._sma < 0.0) & self.is_closed()).any():
            raise ValueError(
                "negative semi-major axis for elliptical/circular eccentricity"
            )

        if ((self._sma != np.inf) & self.is_parabolic()).any():
            raise ValueError("finite semi-major axis for parabolic eccentricity")

    def is_closed(self) -> np.ndarray[bool]:
        return self._ecc < 1.0

    def is_open(self) -> np.ndarray[bool]:
        return not self.is_closed()

    def is_circular(self) -> np.ndarray[bool]:
        return self._ecc == 0.0

    def is_elliptical(self) -> np.ndarray[bool]:
        return self._ecc < 1.0

    def is_parabolic(self) -> np.ndarray[bool]:
        return self._ecc == 1.0

    def is_hyperbolic(self) -> np.ndarray[bool]:
        return self._ecc > 1.0

    def semi_major_axis(self) -> np.ndarray[np.float64]:
        return self._sma

    def eccentricity(self) -> np.ndarray[np.float64]:
        return self._ecc

    def inclination(self) -> np.ndarray[np.float64]:
        return self._inc

    def right_ascension(self) -> np.ndarray[np.float64]:
        return self._raan

    def argument_of_periapsis(self) -> np.ndarray[np.float64]:
        return self._aop

    def true_anomaly(self) -> np.ndarray[np.float64]:
        return self._ta

    def semi_minor_axis(self) -> np.ndarray[np.float64]:
        return semi_minor_axis(self.semi_major_axis(), self.eccentricity())

    def semi_latus_rectum(self) -> np.ndarray[np.float64]:
        return semi_latus_rectum(self.semi_major_axis(), self.eccentricity())

    def periapsis_radius(self) -> np.ndarray[np.float64]:
        return periapsis_radius(self.semi_major_axis(), self.eccentricity())

    def apoapsis_radius(self) -> np.ndarray[np.float64]:
        return apoapsis_radius(self.semi_major_axis(), self.eccentricity())

    def apsis_radii(self) -> np.ndarray[np.float64]:
        return apsis_radii(self.semi_major_axis(), self.eccentricity())

    def radius_at(self, ta) -> np.ndarray[np.float64]:
        return self.semi_latus_rectum() / (1.0 + self.eccentricity() * np.cos(ta))

    def radius(self) -> np.ndarray[np.float64]:
        return self.radius_at(self.true_anomaly())

    def mean_motion(self, gm: float) -> np.ndarray[np.float64]:
        return np.sqrt(gm / self.semi_major_axis())

    def period(self, gm: float) -> np.ndarray[np.float64]:
        return 2 * np.pi / self.mean_motion(gm)

    def specific_energy(self, gm: float) -> np.ndarray[np.float64]:
        return -gm / (2.0 * self.semi_major_axis())

    def specific_angular_momentum(self, gm: float) -> np.ndarray[np.float64]:
        return specific_angular_momentum(
            gm, self.semi_major_axis(), self.eccentricity()
        )

    def speed_at(self, gm: float, ta: float) -> np.ndarray[np.float64]:
        energy = self.specific_energy(gm)
        radius = self.radius_at(ta)
        return np.sqrt(2 * (energy + gm / radius))

    def speed(self, gm: float) -> np.ndarray[np.float64]:
        return self.speed_at(gm, self.true_anomaly())

    def to_perifocal_states(self, gm: float):
        ta = self.true_anomaly()
        h = self.specific_angular_momentum(gm)
        radius = self.radius()

        cta = np.cos(ta)
        sta = np.sin(ta)

        states = np.zeros((6, self.count))

        states[0] = cta * radius  # Radial component
        states[1] = sta * radius  # Tangential component
        states[3] = gm / h * self.eccentricity() * sta  # Radial component
        states[4] = gm / h * (1 + self.eccentricity() * cta)  # Tangential component

        return states

    def perifocal_to_inertial_rotation_matrix(self):
        sAop = np.sin(self.argument_of_periapsis())
        cAop = np.cos(self.argument_of_periapsis())
        sRaan = np.sin(self.right_ascension())
        cRaan = np.cos(self.right_ascension())
        sInc = np.sin(self.inclination())
        cInc = np.cos(self.inclination())
        
        out = np.zeros((self.count, 3, 3))

        out[:, 0, 0] =  cRaan * cAop - sRaan * sAop * cInc
        out[:, 0, 1] = -cRaan * sAop - sRaan * cAop * cInc
        out[:, 0, 2] =  sRaan * sInc

        out[:, 1, 0] =  sRaan * cAop + cRaan * sAop * cInc
        out[:, 1, 1] = -sRaan * sAop + cRaan * cAop * cInc
        out[:, 1, 2] = -cRaan * sInc

        out[:, 2, 0] = sAop * sInc
        out[:, 2, 1] = cAop * sInc
        out[:, 2, 2] = cInc

        return out
        

    def to_inertial_states(self, gm: float):
        pass

    @classmethod
    def from_states(cls, gm: float, states: np.ndarray):
        pos = states[:3, :]
        vel = states[3:, :]

        drv = np.sum(pos * vel, axis=0)

        pos_mag = np.linalg.norm(pos, axis=0)
        pos_hat = pos / pos_mag
        vel_mag = np.linalg.norm(vel, axis=0)

        # Semi-major axis
        energy = 0.5 * vel_mag**2 - gm / pos_mag
        sma = -gm / (2 * energy)

        # Inclination
        ang_mom = np.cross(pos, vel, axis=0)
        inc = np.arccos(ang_mom[2] / np.linalg.norm(ang_mom, axis=0))

        # Eccentricity
        ev = (1 / gm) * (pos * vel_mag**2 - drv * vel) - pos_hat
        ecc = np.linalg.norm(ev, axis=0)
        ev_hat = ev / ecc

        # Nodal Axis
        nodal_axis = np.zeros(ang_mom.shape)
        nodal_axis[0] = -ang_mom[1]
        nodal_axis[1] = ang_mom[0]
        nodal_axis_hat = nodal_axis / np.linalg.norm(nodal_axis, axis=0)
        nodal_axis_hat[0, ang_mom[1] == 0.0] = 1.0
        nodal_axis_hat[1, ang_mom[1] == 0.0] = 0.0

        # Right Ascension
        raan = np.arccos(nodal_axis_hat[0])
        raan[nodal_axis_hat[1] < 0.0] = 2 * np.pi - raan[nodal_axis_hat[1] < 0.0]

        # Argument of periapsis
        aop = np.arccos(np.sum(nodal_axis_hat * ev_hat, axis=0))
        aop[ev[2] < 0.0] = 2 * np.pi - aop[ev[2] < 0.0]

        # True Anomaly
        ta = np.arccos(np.clip(np.sum(ev_hat * pos_hat, axis=0), -1, 1))
        ta[drv < 0] = 2 * np.pi - ta[drv < 0]

        return cls(sma, ecc=ecc, inc=inc, aop=aop, raan=raan, ta=ta)
