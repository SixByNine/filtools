from __future__ import division
import numpy as np
import yaml
import filtools
import os


class SimpleIntegerise:
    """
        Converts floats to integers assuming previously digitised data had a uniform
        underlying distribution.

        :param thresh: Ignore values below this to improve performance (default 1e-8).
    """
    def __init__(self,thresh=1e-8):
        self.thresh=thresh

    def __call__(self,block):
        """
            Convert a block of floating point data to integers.
        """
        mask=np.abs(block)>self.thresh
        block[mask] += np.random.uniform(size=np.sum(mask))
        block[np.logical_not(mask)] = 0 # Stops negative values causing issues
        return np.floor(block)


class FluxScale:
    """
        Helps scaling flux density data to the data format.

        This is used either by directly specifying a scale and offset to apply,
        or by setting a gain and Tsys value, or by specifying physical telescope
        parameters to define the gain and Tsys.

    """

    def __init__(self, scale=1, offset=0):
        self.scale  = scale
        self.offset = offset
        self.rms=0.5
        self.gain=1.0
        self.tsys=20.

    def init(self,tsamp,bw,npol_avg=2):
        """
        Call this once you have set the gain and tsys to
        set up the scaling.

        :param tsamp: sample time in seconds.
        :param bw:    channel bandwidth in MHz.
        :param npol_avg: Number of polarisations that have been averaged already (1 or 2).
        """

        sefd = self.tsys/self.gain
        sefd /= np.sqrt(npol_avg*tsamp*abs(bw)*1e6)
        self.scale=self.rms/sefd

    def set_gain_from_physical(self, ndish, D, efficiency):
        """
        Specify physical parameters to determine gain. 

        :param ndish: Number of antennas
        :param D: Diameter of each antenna (m).
        :param efficiency: Efficiency of the observing system (antenna efficiency).
        """

        effective_area = ndish * np.pi*np.power(D/2.0,2)*efficiency
        gain = effective_area / (2 * 1.38064852e3)
        self.gain=gain


    def to_flux_units(self,buffer):
        """ Convert from file units to flux density units """
        return (buffer-self.offset)/self.scale

    def to_file_units(self,buffer):
        """ Convert from flux density units to file units """
        return buffer * self.scale + self.offset

    def set_gain_from_file(self,file):
        """
        Read a yaml file with a telescope definition. Properties are: tsys, n_antenna, diameter, efficiency and/or gain.
        """
        params={}
        lfile = os.path.join(os.path.dirname(filtools.__file__),"site_data" ,file)
        try:
            with open(file, 'r') as f:
                try:
                    params = yaml.safe_load(f)
                except yaml.YAMLError as exc:
                    print(exc)
        except:
            try:
                with open(lfile, 'r') as f:
                    try:
                        params = yaml.safe_load(f)
                    except yaml.YAMLError as exc:
                        print(exc)
            except:
                pass

        if "tsys" in params:
            self.tsys=float(params['tsys'])
        if "n_antenna" not in params:
            params["n_antenna"] = 1
        if "diameter" in params and "efficiency" in params and "n_antenna" in params:
            self.set_gain_from_physical(ndish=int(params['n_antenna']),D=float(params['diameter']),efficiency=float(params['efficiency']))
        elif "gain" in params:
            self.gain=float(params['gain'])
        else:
            print(params)
            raise Exception("Cannot load telescope info from '{}' or '{}'",file,lfile)

    
