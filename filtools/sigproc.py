from __future__ import division, print_function
import struct
import os
import numpy as np

def _check_modifiable(meth):
    def inner(self,*args,**kwargs):
        if self.state==0:
            meth(self,*args,**kwargs)
        else:
            raise Exception("Cannot modify file whilst file is open.")
    return inner

class FilterbankIO:
    """
    This class deals with reading from and writing to sigproc filterbank files.

    To read a filterbank file, first set up a new :class:`FilterbankIO` object, and use the
    :meth:`read_header` method to read the header parameters. The FilterbankIO object is then ready to 
    be used for reading. To read the file, we use the :meth:`read_block` method to read samples from the file.

    Example:
    
    .. code-block:: python

        inputfile = filtools.FilterbankIO()
        inputfile.read_header("example.fil")
        data = inputfile.read_block(1024) # Reads 1024 samples from the file

    .. note::
        Currently :class:`FilterbankIO` supports sigproc data files with the following properties:
         * Single IF (i.e. single polarisation)
         * Fixed channel bandwidth (i.e. each channel has same bandwidth/offset)
        And, data type is one of: 
         * 1-bit unsigned integer
         * 2-bit unsigned integer
         * 4-bit unsigned integer
         * 8-bit unsigned integer
         * 16-bit unsigned integer
         * 32-bit floating point
 

    .. py:attribute:: header

        The header dict can be used to provide direct access to the sigproc header parameters. Parameters changed here
        before calling :meth:`write_header` can change the ouput header parameters. However, it is usually best to access
        information about the data using the other methods provided below, as the structure can be re-used for other data
        types.
    
       
    """
    _inttypes = ['machine_id', 'telescope_id', 'data_type', 'nchans',\
                    'nbits', 'nifs', 'scan_number', 'barycentric','pulsarcentric',\
                    'nbeams', 'ibeam']
    _strtypes = ['source_name','rawdatafile']
    _dbltypes = ['tstart','tsamp','fch1','foff','refdm','az_start','za_start','src_raj','src_dej']
    _chrtypes = ['signed']
    
    def __init__(self):
        self.header={}
        self._debug=False
        self.state=0
        self._read_block = self.read_block
        self._write_block=self.write_block

    
    def read_header(self,filename):
        """
        Read header from sigproc filterbank file

        :param filename: filename to read
        :type filename: str
        :returns: Number of bytes read (i.e. size of header)

        This method will also put the :class:`FilterbankIO` object in a state for reading the 
        file, i.e. enabling the :meth:`read_block` method.

        Header parameters can be accessed by the methods of the :class:`FilterbankIO` object, or
        by direct interface with the sigproc API using the :data:`header` dict member.
        """
        def read_string(file):
            nchar = np.fromfile(file,dtype=np.int32,count=1)[0]
            if nchar == 0:
                return ""
            if nchar < 1 or nchar > 80:
                raise Exception("Cannot parse filterbank header (Nchar was {} when reading string).".format(nchar))
            byte_data = file.read(nchar)
            string_data = byte_data.decode("UTF-8")
            return string_data
        self.file = open(filename, "rb")
        self.header_position=self.file.tell()

        key = read_string(self.file)
        bytes_read = len(key)+4
        if key=="HEADER_START":
            key = read_string(self.file)
            bytes_read += len(key)+4
            while key != "HEADER_END":
                if self._debug:
                    print("FilterbankIO: Got '{}'".format(key))
                if key in FilterbankIO._inttypes:
                    self.header[key] = np.fromfile(self.file,dtype=np.int32,count=1)[0]
                    bytes_read += 4
                elif key in FilterbankIO._strtypes:
                    self.header[key] = read_string(self.file)
                    bytes_read += len(self.header[key])+4
                elif key in FilterbankIO._dbltypes:
                    self.header[key] = np.fromfile(self.file,dtype=np.float64,count=1)[0]
                    bytes_read += 8
                elif key in FilterbankIO._chrtypes:
                    self.header[key] = np.fromfile(self.file,dtype=np.int8,count=1)[0]
                    bytes_read += 1
                else:
                    raise Exception("Cannot parse filterbank header, key '{}' not understood".format(key))
                key = read_string(self.file)
                bytes_read += len(key)+4
        else:
            raise Exception("Cannot parse filterbank header, HEADER_START was not found")
        self.state=1
        self.data_position=self.file.tell()

        
        return bytes_read
        
    def write_header(self,filename):
        """
        Write the current state to a new sigproc filterbank file

        :param filename: filename to write
        :returns: Number of bytes written (i.e. size of header)

        This method will also put the :class:`FilterbankIO` object in a state for writing data to the 
        file, i.e. enabling the :meth:`write_block` method.
        """

        def write_string(file,string):
            l=len(string)
            byte_data= struct.pack("i{}s".format(l),l,string.encode("UTF-8"))
            file.write(byte_data)
        
        
        self.file = open(filename, "wb")
        self.header_position=self.file.tell()

        write_string(self.file,"HEADER_START")
        for key in self.header:
            write_string(self.file,key)
            if key in FilterbankIO._inttypes:
                self.file.write(struct.pack("i",self.header[key]))
            elif key in FilterbankIO._strtypes:
                write_string(self.file,self.header[key])
            elif key in FilterbankIO._dbltypes:
                self.file.write(struct.pack("d",self.header[key]))
            elif key in FilterbankIO._chrtypes:
                self.file.write(struct.pack("b",self.header[key]))
            else:
                raise Exception("Cannot understand filterbank header for writing, key '{}' not understood".format(key))
        write_string(self.file,"HEADER_END")
        
        self.data_position=self.file.tell()
        self.state=2
  
    def clone_header(self):
        """
        Create a new :class:`FilterbankIO` object with the header parameters copied from this object

        :returns: A new :class:`FilterbankIO` object
        """
        ret = FilterbankIO()
        ret.header = self.header.copy()
        return ret

    def seconds_since_start(self,sample):
        """ Get the number of seconds since start of observation at sample number ``sample``."""
        return sample * self.header['tsamp']

    def seconds_since(self,sample, mjdref):
        """
        Get the number of seconds since ``mjdref`` at sample number ``sample``.

        :param sample: Integer sample number
        :param refmjd: Reference mjd
        :returns: seconds since ``mjdref`` at sample ``sample``
        """
        return 86400.0*(self.header['tstart'] - mjdref) + self.seconds_since_start(sample)

    def mjd_at_sample(self,sample):
        """Get the MJD at given sample number"""
        return self.header['tstart'] + (sample * self.header['tsamp'])/86400.0
    
    @_check_modifiable
    def set_mjd_at_start(self,mjd):
        """Set the MJD at the start of the observation."""
        self.header['tstart'] = mjd
    
    def frequency_table(self):
        """Get an array of centre frequency of each channel"""
        return np.linspace(self.header['fch1'],self.header['fch1'] + self.header['nchans']*self.header['foff'],self.header['nchans'],endpoint=False)


    def seek_to_sample(self,sample):
        """
        Position the file pointer at given sample number (index from zero).

        The next call to :meth:`read_block` or :meth:`write_block` will start at the given sample. E.g.

        .. code-block:: python

            fil.seek_to_sample(100)
            block = fil.read_block(128) # block[0] is sample number 100
        """
        bytes_from_data_start = (sample * self.header['nchans'] * 8)//self.header['nbits']
        self.file.seek(self.data_position+bytes_from_data_start,os.SEEK_SET)
    
    def current_position(self):
        """
        Print the current position of the file pointer in samples.

        The next call to :meth:`read_block` or :meth:`write_block` will read or write from this sample.
        """
        return ((self.file.tell()-self.data_position) * self.header['nbits']) // (8*self.header['nchans'])
        
    def bits_per_sample(self):
        """Number of bits per sample in data file"""
        return self.header['nbits']

    @_check_modifiable
    def set_bits_per_sample(self,val):
        """Set the number of bits per sample to write"""
        self.header['nbits'] = int(val)

    def total_bytes(self):
        """Compute the total number of bytes in the file"""
        return os.fstat(self.file.fileno()).st_size

    def total_channels(self):
        """The total number of frequency channels in the file"""
        return self.header['nchans']

    @_check_modifiable
    def set_total_channels(self,val):
        self.header['nchans'] = val

    def total_samples(self):
        """The total number of samples in the file"""
        header_size = self.data_position - self.header_position
        return int(8*(self.total_bytes()-header_size) / self.header['nchans'] / self.header['nbits'])
        
    def sample_interval(self):
        """The time between each sample"""
        return self.header['tsamp']

    @_check_modifiable
    def set_sample_interval(self,val):
        self.header['tsamp'] = val


    def reference_dm(self):
        if 'refdm' in self.header:
            return self.header['refdm']
        else:
            return 0.0

    @_check_modifiable
    def set_reference_dm(self,val):
        self.header['refdm'] = val

    def channel_bandwidth(self):
        """The channel bandwidth."""
        return self.header['foff']

    def centre_frequency(self):
        """The centre frequency"""
        return self.header['fch1'] + self.header['foff'] * (self.header['nchans']-1)/2.0

    @_check_modifiable
    def set_frequency_table(self,new_centre_freq,new_ch_bw,new_nchans):
        self.header['foff'] = new_ch_bw
        self.header['nchans'] = new_nchans
        self.header['fch1'] = new_centre_freq - self.header['foff'] * (self.header['nchans']-1)/2.0


    def read_block(self,nsamples,dtype=np.float32):
        """
        Read a block of data. 

        This will read up to ``nsamples`` from the file, returning a 2-d array indexed first by sample and then by channel.  
        i.e. the shape of the returned array will be (nsamples_read,nchannels).

        If the file has reached the end of file, then less than nsamples will be read, and if past the end of the file, an array of size (0,nchannels) will be returned.

        .. note::
            On first call this method is dynamically replaced using :meth:`get_reader` to define an optimised verison of the reader method.
            This should be invisiable to most users.

        :param nsamples: The number of samples to read.
        :param dtype: A numpy compatible data type. Data will be converted to this type using the astype numpy method, and should be at least large enough to hold the read data types.
        :returns: A 2-d numpy array of requested type.
        """

        self.read_block=self.get_reader()
        return self.read_block(nsamples,dtype)

    def get_reader(self):
        """This method returns a function optimised to read the input data.`"""
        
        def read_block_32bit(nsamples,dtype=np.float32):
            nchans=  self.header['nchans']
            nbytes = nsamples * nchans
            raw = np.fromfile(self.file,dtype=np.float32,count=nbytes)
            return raw.reshape(-1,nchans).astype(dtype)
        
        def read_block_16bit(nsamples,dtype=np.float32):
            nchans=  self.header['nchans']
            nbytes = nsamples * nchans
            raw = np.fromfile(self.file,dtype=np.uint16,count=nbytes)
            return raw.reshape(-1,nchans).astype(dtype)
        
        def read_block_8bit(nsamples,dtype=np.float32):
            nchans=  self.header['nchans']
            nbytes = nsamples * nchans
            raw = np.fromfile(self.file,dtype=np.uint8,count=nbytes)
            return raw.reshape(-1,nchans).astype(dtype)
        
        def read_block_4bit(nsamples,dtype=np.float32):
            nchans=  self.header['nchans']
            nbytes = int((nsamples * nchans)/2)
            raw = np.fromfile(self.file,dtype=np.uint8,count=nbytes)    
            out = np.zeros(len(raw) * 2, dtype=dtype)
            out[::2] += np.bitwise_and(raw,0x0f)
            out[1::2] += np.right_shift(raw,4)
            return out.reshape(-1,nchans).astype(dtype)       
        
        def read_block_2bit(nsamples,dtype=np.float32):
            nchans=  self.header['nchans']
            nbytes = int((nsamples * nchans)/4)
            raw = np.fromfile(self.file,dtype=np.uint8,count=nbytes)
            out = np.zeros(len(raw) * 4, dtype=dtype)
            out[::4] += np.bitwise_and(raw,0x03)
            out[1::4] += np.right_shift(np.bitwise_and(raw,0x0f),2)
            out[2::4] += np.right_shift(np.bitwise_and(raw,0x3f),4)
            out[3::4] += np.right_shift(np.bitwise_and(raw,0xff),6)
            return out.reshape(-1,nchans).astype(dtype)
        
        def read_block_1bit(nsamples,dtype=np.float32):
            nchans=  self.header['nchans']
            nbytes = int((nsamples * nchans)/8)
            raw = np.fromfile(self.file,dtype=np.uint8,count=nbytes)
            
            out = np.zeros(len(raw) * 8, dtype=dtype)
            out[::8] += np.bitwise_and(raw,0x01)
            out[1::8] += np.right_shift(np.bitwise_and(raw,0x03),1)
            out[2::8] += np.right_shift(np.bitwise_and(raw,0x07),2)
            out[3::8] += np.right_shift(np.bitwise_and(raw,0x0f),3)
            out[4::8] += np.right_shift(np.bitwise_and(raw,0x1f),4)
            out[5::8] += np.right_shift(np.bitwise_and(raw,0x3f),5)
            out[6::8] += np.right_shift(np.bitwise_and(raw,0x7f),6)
            out[7::8] += np.right_shift(np.bitwise_and(raw,0xff),7)
            return out.reshape(-1,nchans).astype(dtype)
        
        if self.state==1:
            if self.header['nbits'] == 32:
                return read_block_32bit
            
            elif 'signed' in self.header and self.header['signed']!= 1:
                    raise Exception("FilterbankIO Only works on unsigned data")
            elif self.header['nbits'] == 16:
                return read_block_16bit
            elif self.header['nbits'] == 8:
                return read_block_8bit
            elif self.header['nbits'] == 4:
                return read_block_4bit
            elif self.header['nbits'] == 2:
                return read_block_2bit
            elif self.header['nbits'] == 1:
                return read_block_1bit
            else:
                raise Exception("FilterbankIO does not work with {} bits data".format(self.header['nbits']))
        else:
            raise Exception("File not opened for reading")
        

    def write_block(self,block):
        """
        Write a block of data to a file opened for writing.

        The data block should be a numpy array in the same format as would be obtained from :meth:`read_block`.
        The data will be truncated and cast to the required type to write to the file using the ndarray.astype method from numpy.

        The block should have shape (nsamples,nchannels), where nsamples is the number of samples to write.

        :param block: A 2-d numpy array to write.

        .. note::
            On first call this this method is dynamically replaced using :meth:`get_writer` to define an optimised verison of the writer method.
            This should be invisiable to most users.
        """
        self.write_block = self.get_writer()
        return self.write_block(block)

    def get_writer(self):
        """Returns an optimised write function. Used by :meth:`write_block`."""
        def write_block_32bit(block):
            raw = block.flatten(order='C').astype('float32')
            self.file.write(raw.tobytes())
            
        def write_block_16bit(block):
            block[block>65535] = 65535
            block[block < 0] = 0
            raw = block.flatten(order='C').astype('uint16')

            self.file.write(raw.tobytes())
            
        def write_block_8bit(block):
            block[block>255] = 255
            block[block < 0] = 0
            raw = block.flatten(order='C').astype('uint8')
            
            self.file.write(raw.tobytes())

        def write_block_4bit(block):
            block[block > 15] = 3
            block[block < 0] = 0
            raw = block.flatten(order='C').astype('uint8')
            
            out = np.zeros(int(len(raw)/2),dtype=np.uint8)
            out += raw[::2]
            out += np.left_shift(raw[1::2],4)
            self.file.write(out.tobytes())
            
        def write_block_2bit(block):
            block[block > 3] = 3
            block[block < 0] = 0
            
            raw = block.flatten(order='C').astype('uint8')
            out = np.zeros(int(len(raw)/4),dtype=np.uint8)
            out += raw[::4]
            out += np.left_shift(raw[1::4],2)
            out += np.left_shift(raw[2::4],4)
            out += np.left_shift(raw[3::4],6)
            self.file.write(out.tobytes())
            
        def write_block_1bit(block):
            block[block > 1] = 1
            block[block < 0] = 0
            raw = block.flatten(order='C').astype('uint8')
            out = np.zeros(int(len(raw)/8),dtype=np.uint8)
            out += raw[::8]
            out += np.left_shift(raw[1::8],1)
            out += np.left_shift(raw[2::8],2)
            out += np.left_shift(raw[3::8],3)
            out += np.left_shift(raw[4::8],4)
            out += np.left_shift(raw[5::8],5)
            out += np.left_shift(raw[6::8],6)
            out += np.left_shift(raw[7::8],7)
            self.file.write(out.tobytes())

            
        if self.state == 2:
            if self.header['nbits'] == 32:
                return write_block_32bit
            
            elif 'signed' in self.header and self.header['signed']!= 1:
                    raise Exception("FilterbankIO Only works on unsigned data")      
            elif self.header['nbits'] == 16:
                return write_block_16bit
            elif self.header['nbits'] == 8:
                return write_block_8bit
            elif self.header['nbits'] == 4:
                return write_block_4bit
            elif self.header['nbits'] == 2:
                return write_block_2bit
            elif self.header['nbits'] == 1:
                return write_block_1bit
            else:
                raise Exception("FilterbankIO does not work with {} bits data".format(self.header['nbits']))
        else:
            raise Exception("File not opened for writing")
            
    
    def close(self):
        """Close the input file, reset the read/write state of the :class:`FilterbankIO` object"""
        if self.state>0:
            self.file.close()
            self.state=0
        self.write_block=self._write_block
        self.read_block = self._read_block
