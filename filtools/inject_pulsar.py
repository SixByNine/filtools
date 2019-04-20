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

    def get_spectrum(self,t1,t2,mjdref):
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
    
class simple_propagation_model:
    r"""Models simple cold-plasma dispersion delays.

       Optionally includes smearing due to intra-channel delay.
       Disable to model coherent dedispersion.

       The DM delay is modeled as 

       .. math:: t_\mathrm{DM} = \frac{\nu^{-2}\mathrm{DM}}{k_\mathrm{DM}}


       :param freq:        array of centre frequency for each channel in MHz. The algorithm assumes equally spaced channels.
       :param input_model: an object which implements input_model methods
       :param dm:          the simulated DM in cm$^-3$pc.
       :param intra_chan:  if True, will model the smearing in the pulse due to the intra-channel delays
       :param dm_const:    override to change the DM constant (:math:`k_\mathrm{DM}`) used to compute the DM delays.
    """

    def __init__(self,freq,input_model,dm,intra_chan=True,dm_const=2.41e-4):
        self.input_model = input_model
        
        ch_w = freq[1]-freq[0]
        self.freq_start = freq-ch_w/2.0
        self.freq_end   = freq+ch_w/2.0
        print(freq)

        self.nchan = len(freq)
        
        
        if intra_chan:
            self.dm_delays_start = np.power(self.freq_start,-2)*dm/dm_const
            self.dm_delays_end = np.power(self.freq_end,-2)*dm/dm_const
        else:
            self.dm_delays_start = np.power(freq,-2)*dm/dm_const
            self.dm_delays_end   = self.dm_delays_start.copy()
       

    
    
    def get_spectrum(self,t1,t2,mjdref):
        """
        Compute average flux density between t1 and t2 in each frequency channel.
        
        Output is in flux density units.
        
        :param t1: Start time in seconds, relative to mjdref
        :param t2: End time in seconds, relative to mjdref
        :returns:  Average flux density between t1 and t2.
        """

        spec = np.zeros(self.nchan,dtype=np.float32)
        spec += self.input_model.get_spectrum(t1-self.dm_delays_start,\
                                          t2-self.dm_delays_end,\
                                          mjdref)
        return spec

class simple_scattering_model:
    r"""Approximates thin-screen scattering in the ionised interstellar medium.

        Effectively convolves the incoming signal with an exponential of the form

        .. math:: \exp{\left(-\frac{t}{t_\mathrm{scatt}}\left(\frac{\nu_\mathrm{ref}}{\nu}\right)^\alpha\right)}

        .. warning::
            The simple scattering model assumes that each sample to be generated is called in sequence.
            I.e. no calls to get_spectrum occur out of sequence. This will prevent multi-threaded operation
            of codes using this module.

        :param input_model: the input signal generator
        :type input_model: input_model
        :param freq:        array of frequency channels in MHz.
        :param t_scatt:     the scattering timescale.
        :param dx:          the sampling interval in the data file.
        :param scatt_idx:   the spectral index of the scattering timescale (:math:`\alpha`)

    """
    def __init__(self,t_scatt,dx,input_model,freq,scatt_idx=-4.0,ref_freq=1400.0):
        
        self.t_scatt = t_scatt
        
        self.const=np.exp(-dx/(self.t_scatt*np.power(freq/ref_freq,scatt_idx)))
        self.rconst = 1.0 - self.const
        self.input_model=input_model
        nchan = len(freq)
        self.store=np.zeros(nchan)
    
    def get_spectrum(self,t1,t2,mjdref):    
        """
        Compute average flux density between t1 and t2 in each frequency channel.
        
        Output is in flux density units.

        .. warning::
            The simple scattering model assumes that each sample to be generated is called in sequence.
            I.e. no calls to get_spectrum occur out of sequence. This will prevent multi-threaded operation
            of codes using this module.
        
        :param t1: Start time in seconds, relative to mjdref
        :param t2: End time in seconds, relative to mjdref
        :returns:  Average flux density between t1 and t2.
        """

        spec = self.input_model.get_spectrum(t1,t2,mjdref)
        self.store += spec
        spec = self.store * self.rconst
        self.store *= self.const
        return spec
        


class sinusoid_input_model:
    def __init__(self,phase_model,freq):
        self.phase_model = phase_model
        self.freq = freq


    def get_spectrum(self,t1,t2,mjdref):
        """
        Average flux density between t1 and t2
        
        Output is in flux density units.
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
        
    

class simple_integerise:
    def __init__(self,thresh=1e-6):
        self.thresh=thresh
        
    def __call__(self,block):
        mask=np.abs(block)>self.thresh
        block[mask] += np.random.uniform(size=np.sum(mask))
        return np.floor(block)

        
    
class fluxscale:
    def __init__(self, scale=1, offset=0):
        self.scale  = scale
        self.offset = offset
        
    def set_by_gain(self,gain,tsys,tsamp,bw,rms,offset):
        
        sefd = tsys/gain
        sefd /= np.sqrt(tsamp*bw*1e6)
        
        self.scale=rms/sefd
        
        self.offset
        
    
    def set_by_params(self, ndish, D, efficiency, tsys, tsamp, bw, rms,offset):
        effective_area = ndish * np.pi*np.power(D/2.0,2)*efficiency
        gain = effective_area / (2 * 1.38064852e3)
        self.set_by_gain(gain,tsys,tsamp,bw,rms,offset)
    
    def to_flux_units(self,buffer):
        return (buffer-self.offset)/self.scale
 
    def to_file_units(self,buffer):
        return buffer * self.scale + self.offset
        
    
