from __future__ import division
import numpy as np
import yaml
import filtools
import os


class SimpleIntegerise:
    def __init__(self,thresh=1e-6):
        self.thresh=thresh
    def __call__(self,block):
        mask=np.abs(block)>self.thresh
        block[mask] += np.random.uniform(size=np.sum(mask))
        return np.floor(block)


class FluxScale:
    def __init__(self, scale=1, offset=0):
        self.scale  = scale
        self.offset = offset
        self.rms=0.5
        self.gain=1.0
        self.tsys=20.

    def init(self,tsamp,bw):
        sefd = self.tsys/self.gain
        sefd /= np.sqrt(tsamp*abs(bw)*1e6)
        self.scale=self.rms/sefd

    def set_gain_from_physical(self, ndish, D, efficiency):
        effective_area = ndish * np.pi*np.power(D/2.0,2)*efficiency
        gain = effective_area / (2 * 1.38064852e3)
        self.gain=gain

    def to_flux_units(self,buffer):
        return (buffer-self.offset)/self.scale

    def to_file_units(self,buffer):
        return buffer * self.scale + self.offset

    def set_gain_from_file(self,file):
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

    
