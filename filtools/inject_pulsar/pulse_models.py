from __future__ import division
import numpy as np

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

    .. math:: I(\phi) = \sin{(2\pi\phi)}

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



