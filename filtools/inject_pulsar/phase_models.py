from __future__ import division
import numpy as np
from scipy import interpolate as interp

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
        return self.f0 + self.f1*t + self.f2*t*t/2.0


class precision_predictor_phase_model:
    """
    A phase model that uses a phase predictor (or polyco) from
    tempo or tempo2. It uses the t2pred library to convert between
    time and phase, so can model a wide range of pulsar behaviour.
    This model calls t2pred on every request so is slow, but very
    accurate.

    :param predictor: The phase predictor to use
    :type predictor:  A t2pred.phase_predictor object

    .. note::
       This module will only function if you have installed the t2pred python library in tempo2.
       This can be installed using setup.py in the python/t2pred directory in the tempo2 source code.
    """
    
    def __init__(self,predictor):
        self.predictor=predictor


    def get_phase(self, time, mjdref):
        """Get the pulse phase at time in seconds relative to refmjd"""
#        phase=np.zeros(np.product(time.shape))
#        i=0
#        for t in (time.flatten().astype(np.longdouble)/86400.0)+mjdref:
#            phase[i] = self.predictor.getPhase(t,1e9)
#            i+=1
#        return phase.reshape(time.shape)
        return self.predictor.getPhase_array(time.astype(np.longdouble)/86400.0 + mjdref,1400)
#        return self.predictor.getPhase_array(time/86400.0 + mjdref,1e9)


    def get_frequency(self, time, mjdref):
        """Get the pulse frequency at time in seconds relative to refmjd"""
        freqs=np.zeros_like(time)
        i=0
        for t in (time.astype(np.longdouble)/86400.0)+mjdref:
            freqs[i] = self.predictor.getPhase(t,1400)
            i+=1
        return freqs

class predictor_phase_model:
    """
    A phase model that uses a phase predictor (or polyco) from
    tempo or tempo2. It uses the t2pred library to convert between
    time and phase, so can model a wide range of pulsar behaviour.
    This version keeps its own polynomial interpolation of the
    predictor over short spans to reduce the number of calls to t2pred
    and greatly increase the performance for a very minor precision
    reduction.

    :param predictor: The phase predictor to use
    :type predictor:  t2pred.phase_predictor
    :param dt: The interpolation window size in seconds

    .. note::
       This module will only function if you have installed the t2pred python library in tempo2.
       This can be installed using setup.py in the python/t2pred directory in the tempo2 source code.
    """

    def __init__(self,predictor,dt=0.1):
        self.predictor=predictor
        self.dt=dt
        self.default_tspan=10
        self.mjdref=None

    def make_poly(self,t):
        """
        This routine is used internally to make the polynomial approimation of the timing model.
        If a new-ish version of scipy is installed it uses scipy.interp.CubicHermiteSpline, or
        otherwise uses a scipy.interp.BPoly.
        """
        phs = self.predictor.getPhase_array(t.astype(np.longdouble)/86400.0 + self.mjdref,1400)
        freq= self.predictor.getFrequency_array(t.astype(np.longdouble)/86400.0 + self.mjdref,1400)
        try:
            # Newer versions of scipy
            pp = interp.CubicHermiteSpline(t,phs,freq)
        except:
            # fall back if CubicHermiteSpline is not existing
            yi=np.vstack((phs,freq)).T
            pp = interp.BPoly.from_derivatives(t,yi)
        return pp
        
    def get_phase(self, time, mjdref):
        """
        Get the pulse phase at time in seconds relative to refmjd.
        """
        if self.mjdref==None:
            self.mjdref=mjdref
            self.tstart=np.amin(time)
            t = np.linspace(self.tstart,self.tstart+self.default_tspan,self.default_tspan/self.dt,endpoint=True)
            self.tend=t[-1]
            self.pp = self.make_poly(t)

        time += (mjdref-self.mjdref)*86400.0
        while np.amin(time) < self.tstart:
            self.tstart -= self.default_tspan
            t = np.linspace(self.tstart,self.tstart+self.default_tspan,self.default_tspan/self.dt,endpoint=True)
            pp = self.make_poly(t)
            self.pp.extend(pp.c,pp.x[:-1])
        while np.amax(time) > self.tend:
            self.tend += self.default_tspan
            t = np.linspace(self.tend-self.default_tspan,self.tend,self.default_tspan/self.dt,endpoint=True)
            pp = self.make_poly(t)
            self.pp.extend(pp.c,pp.x[1:])
        return self.pp(time)



    def get_frequency(self, time, mjdref):
        """Get the pulse frequency at time in seconds relative to refmjd"""
        freqs=np.zeros_like(time)
        i=0
        for t in (time.astype(np.longdouble)/86400.0)+mjdref:
            freqs[i] = self.predictor.getPhase(t,1400)
            i+=1
        return freqs



