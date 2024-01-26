import numpy as np
import yaml
class inject_rfi:

    def configure_db(self,filename,nchans,tsamp,fch1,foff):
        r"""
        This function will be used to setup the RFI injection configuration with desired nature of RFI. 
        A YAML file containing all the information about RFI is read, samples that are going to be affected are figured out in this function.

        :param filename:    Path to the YAML configuration file
        :param bsize:   Blocksize of operation (default from ft_inject_pulsar = 256)
        :param nchans:  Number of Frequency Channels
        :param tsamp:   Sampling time
        :param fch1:    Frequency of the first channel
        :param foff:    Frequency resolution or difference in frequency between two channels
        :returns:       True (if code executes successfully) or False (otherwise)
        """


        self.nchans = nchans        
        rfi_file = yaml.safe_load(open(filename))                   #loading the RFI yaml files into list of dictionary
        rfi_file = sorted(rfi_file, key=lambda d: d['tstart'])      #sorting the dictionaries based on 'tstart' of every RFI

        #The elements are sorted because referring to frequency indices would be difficult (some mismatch was observed)
        
        self.rfi_begin, self.rfi_end, self.rfi_freq_start, self.rfi_freq_stop, self.rfi_mag = [],[],[],[],[]
        try:
            for i in rfi_file:
                if i['period'] > 0:
                    #periodic RFI is treated as many discrete RFI chunks. Hence, the entry is appended as many chunks of discrete RFI and would be treated in regular fashion
                    tnow=i['tstart']
                    while tnow<i['tstop']:
                        self.rfi_begin.append(int(np.floor(tnow/tsamp)))  
                        self.rfi_end.append(int(np.floor((tnow+i['period']*i['duty'])/tsamp)))
                        self.rfi_freq_start.append(int((fch1-i['fstart']+foff*nchans)/foff))
                        self.rfi_freq_stop.append(int((fch1-i['fstop']+foff*nchans)/foff))
                        self.rfi_mag.append(i['mag'])
                        tnow+=i['period']
                else:
                    self.rfi_begin.append(int(np.floor(i['tstart']/tsamp)))     #Converts units of time to corresponding sample number before appending
                    self.rfi_end.append(int(np.ceil(i['tstop']/tsamp)))
                    self.rfi_freq_start.append(int((fch1-i['fstart']+foff*nchans)/foff))    #converts units of frequency into corresponding channel number before appending
                    self.rfi_freq_stop.append(int((fch1-i['fstop']+foff*nchans)/foff))
                    self.rfi_mag.append(i['mag'])
        
            #Converting the lists to numpy array for easy operation
            self.rfi_begin = np.array(self.rfi_begin)
            self.rfi_end = np.array(self.rfi_end)
            self.rfi_freq_start = np.array(self.rfi_freq_start)
            self.rfi_freq_stop = np.array(self.rfi_freq_stop)
            self.rfi_mag = np.array(self.rfi_mag)
            
        except:
            return False

        return True
    
    
    def injector(self,bin_start,block_length):
        r"""
        Generates a mask which will be added to the filterbank being operated
        
        :param bin_start:   Time-sample number of the block that is being worked on
        :returns:   Mask of size same as block with affected channels and time filled with respective magnitude (:math:`\times\sigma`)
        """
        mask = np.empty((self.nchans,block_length))      #Creating an empty mask (local: clearing up memory)

        #Creating a list of RFI chunks that needs to be added in the current block
        rfi_list = np.argwhere(bin_start < self.rfi_end[np.argwhere(self.rfi_begin<bin_start+block_length)])[:,0].T
        

        for i in rfi_list:
            #Iterating over the RFI chunks that needs to be injected, adjusting the 'tstart' and 'tstop' in order to fit them within the selected block
            tstart = max(self.rfi_begin[i],bin_start) - bin_start    #If the start of RFI was much before the block, then adjust the RFI start to the start of current block            
            tstop = min(self.rfi_end[i]-bin_start,block_length)    #If the end of RFI chunk is after the end of chosen block, then adjust it to the end of current block
            
            fstart = self.rfi_freq_start[i]
            fstop = self.rfi_freq_stop[i]
            #Fill the entries of affected channel with ones, scaled with respective magnitude
            mask[fstart:fstop,tstart:tstop] += np.ones((fstop-fstart,tstop-tstart)) * self.rfi_mag[i]
        
        return mask.T