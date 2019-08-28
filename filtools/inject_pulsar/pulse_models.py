from __future__ import division
import numpy as np
from scipy import interpolate as interp

class delta_pulse_model:
    """
    Models a pulsar as a delta function.

    Only has flux at phase zero. Phase average flux density is 1 unit.

    :param phase_model: Object that implements the get_phase method.
    :param freq: An array of observing frequnecies (ignored)
    """

    def __init__(self,phase_model,freq):
        self.phase_model = phase_model
        self.freq=freq



    def __call__(self,t1,t2,mjdref):
        """
        Compute average flux density between t1 and t2 in each frequency channel.

        Output is in flux density units.

        :param t1: Start time in seconds, relative to mjdref
        :param t2: End time in seconds, relative to mjdref
        :returns:  Average flux density between t1 and t2.
        """

        p1 = self.phase_model.get_phase(t1,mjdref)
        p2 = self.phase_model.get_phase(t2,mjdref)

        pr =  np.floor(p1)

        p1 -= pr
        p2 -= pr

        # Integrate over delta function at zero phase
        val=np.zeros_like(p2,dtype=np.float32)
        val += np.floor(p2)

        # Get back to flux units
        val /= (p2-p1)


        # This is a simple model that has no frequency structure.
        return val



class sinusoid_pulse_model:
    r"""
    Models a pulsar as a sinusoid function.

    .. math:: I(\phi) = \cos{(2\pi\phi)+1}

    Where :math:`\phi` is the phase in turns (i.e. between 0 and 1).

    :param phase_model: Object that implements the get_phase method.
    :param freq: An array of observing frequnecies (ignored)
    """

    def __init__(self,phase_model,freq):
        self.phase_model = phase_model
        self.freq = freq


    def __call__(self,t1,t2,mjdref):
        r"""
        Compute average flux density between t1 and t2 in each frequency channel.

        Output is in flux density units.

        :param t1: Start time in seconds, relative to mjdref
        :param t2: End time in seconds, relative to mjdref
        :returns:  Average flux density between t1 and t2.
        """


        p1 = self.phase_model.get_phase(t1,mjdref)
        p2 = self.phase_model.get_phase(t2,mjdref)



        #spin_freq = self.phase_model.get_freq(t2)

        # Integrate cos(2*pi*x).dx between p1 and p2
        # This gives result in flux * phase units
        val = (np.sin(2*np.pi*p2) - np.sin(2*np.pi*p1))/2.0/np.pi + p2 - p1

        # Get back to flux units
        val /= (p2-p1)


        # This is a simple model that has no frequency structure.

        return val


class piecewise_pulse_model:
    r"""
    Models a pulsar as an arbitrary piecewise pulse shape


    :param phase_model: Object that implements the get_phase method.
    :param freq: An array of observing frequnecies (ignored)
    :param profile: An array of pulse intensity as a function of time
    """

    def __init__(self,phase_model,freq,profile):
        self.phase_model = phase_model
        self.freq = freq
        self.orig_profile=profile
        self.orig_profile /= np.mean(self.orig_profile)
        self.nbin = len(profile)
        self.repad(10)

    def repad(self,n):
        r"""
        Pads the internal array with multiple copies of the pulse profile. Used internally to
        ensure that the integration can cope with profiles spanning multiple pulse periods.

        :param n: The number of times to pad the profile.
        """
        self.profile = np.pad(self.orig_profile,pad_width=(0,self.nbin*(n-1)),mode='wrap')
        self.pad=n
        y=self.profile
        x=np.arange(self.nbin*n)/self.nbin
        dx=(x[1]-x[0])

        m = (y[1:] - y[:-1])/dx
        c = y[:-1]

        #self.pp = interp.PPoly([m,c],x)

        cc=np.zeros_like(m)
        cc[1:] = np.cumsum(m*dx*dx/2.0 + c*dx)[:-1]

        self.ipp = interp.PPoly([m/2.0, c, cc],x)


    def __call__(self,t1,t2,mjdref):
        r"""
        Compute average flux density between t1 and t2 in each frequency channel.

        Output is in flux density units.

        :param t1: Start time in seconds, relative to mjdref
        :param t2: End time in seconds, relative to mjdref
        :returns:  Average flux density between t1 and t2.
        """

        p1 = self.phase_model.get_phase(t1,mjdref)
        p2 = self.phase_model.get_phase(t2,mjdref)

        intpart = np.floor(p1)
        p1 -= intpart
        p2 -= intpart

        val = self.ipp(p2)-self.ipp(p1)
        # This is a simple model that has no frequency structure.

        return val/(p2-p1)
