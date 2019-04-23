from __future__ import division
import numpy as np

    
class simple_propagation_model:
    r"""Models simple cold-plasma dispersion delays.

       Optionally includes smearing due to intra-channel delay.
       Disable to model coherent dedispersion.

       The DM delay is modeled as 

       .. math:: t_\mathrm{DM} = \frac{\nu^{-2}\mathrm{DM}}{k_\mathrm{DM}}


       :param freq:        array of centre frequency for each channel in MHz. The algorithm assumes equally spaced channels.
       :param input:       Callable or function.
       :param dm:          the simulated DM in cm$^-3$pc.
       :param intra_chan:  if True, will model the smearing in the pulse due to the intra-channel delays
       :param dm_const:    override to change the DM constant (:math:`k_\mathrm{DM}`) used to compute the DM delays.
    """

    def __init__(self,freq,input,dm,intra_chan=True,dm_const=2.41e-4):
        self.input = input
        
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
       

    
    
    def __call__(self,t1,t2,mjdref):
        """
        Compute average flux density between t1 and t2 in each frequency channel.
        
        Output is in flux density units.
        
        :param t1: Start time in seconds, relative to mjdref
        :param t2: End time in seconds, relative to mjdref
        :returns:  Average flux density between t1 and t2.
        """

        spec = np.zeros(self.nchan,dtype=np.float32)
        spec += self.input(t1-self.dm_delays_start, t2-self.dm_delays_end, mjdref)
        return spec

class simple_scattering_model:
    r"""Approximates thin-screen scattering in the ionised interstellar medium.

        Effectively convolves the incoming signal with an exponential of the form

        .. math:: \exp{\left(-\frac{t}{t_\mathrm{scatt}}\left(\frac{\nu_\mathrm{ref}}{\nu}\right)^\alpha\right)}

        .. warning::
            The simple scattering model assumes that each sample to be generated is called in sequence.
            I.e. no calls to get_spectrum occur out of sequence. This will prevent multi-threaded operation
            of codes using this module.

        :param input: the input signal generator
        :param freq:        array of frequency channels in MHz.
        :param t_scatt:     the scattering timescale.
        :param dx:          the sampling interval in the data file.
        :param scatt_idx:   the spectral index of the scattering timescale (:math:`\alpha`)

    """
    def __init__(self,freq,input,t_scatt,dx,scatt_idx=-4.0,ref_freq=1400.0):
        
        self.t_scatt = t_scatt
        
        self.const=np.exp(-dx/(self.t_scatt*np.power(freq/ref_freq,scatt_idx)))
        self.rconst = 1.0 - self.const
        self.input=input
        nchan = len(freq)
        self.store=np.zeros(nchan)
    

    def __call__(self,t1,t2,mjdref):    
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

        spec = self.input_model(t1,t2,mjdref)
        self.store += spec
        spec = self.store * self.rconst
        self.store *= self.const
        return spec
        

