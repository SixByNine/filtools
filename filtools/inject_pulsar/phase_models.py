from __future__ import division
import numpy as np

class simple_phase_model:
    """
    A very simple model of the rotational phase of a pulsar

    This model is that the pulsar is in the rest frame of the observer
    and follows a polynomial spin where the phase is given by:

    .. math:: \phi(t) = f_0*t + (f_1*t^2)/2.0 + (f_2*t^3)/6.0$

    phase is therefore defined as being from 0 to 1, rather than as an angle.

    :param epoch: The reference epoch (:math:`t=0`)
    :param f0:    Frequency at t=0
    :param f1:    Frequency derivative at t=0
    :param f2:    Frequency second derivative at t=0
    :param accn:  Doppler acceleration at t=0 (overrides choice of f1).
    """

    def __init__(self, epoch, f0, f1=0.0, f2=0.0,accn=None):
        self.f0=f0
        if accn is None:
            self.f1=f1
        else:
            self.f1 = -accn*f0/3e8   
        self.f2=f2
        self.epoch=epoch

    def get_phase(self,time,mjdref):
        """Get the pulse phase at time in seconds relative to refmjd"""
        t = time - (self.epoch-mjdref)*86400.0
        return self.f0*t + self.f1*t*t/2.0 + self.f2*t*t*t/6.0

    def get_frequency(self,time,mjdref):
        """Get the pulse frequency at time in seconds relative to refmjd"""
        t = time - (self.epoch-mjdref)*86400.0
        return self.f0*t + self.f1*t + self.f2*t*t/2.0

