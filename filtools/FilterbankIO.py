from __future__ import division, print_function
import struct
import os
import numpy as np

class FilterbankIO:
    inttypes = ['machine_id', 'telescope_id', 'data_type', 'nchans',\
                    'nbits', 'nifs', 'scan_number', 'barycentric','pulsarcentric',\
                    'nbeams', 'ibeam']
    strtypes = ['source_name','rawdatafile']
    dbltypes = ['tstart','tsamp','fch1','foff','refdm','az_start','za_start','src_raj','src_dej']
    chrtypes = ['signed']
    
    def __init__(self):
        self.header={}
        self._debug=False
        self.state=0
    
    def read_header(self,filename):
        def read_string(file):
            nchar = np.fromfile(file,dtype=np.int32,count=1)[0]
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
                if key in FilterbankIO.inttypes:
                    self.header[key] = np.fromfile(self.file,dtype=np.int32,count=1)[0]
                    bytes_read += 4
                elif key in FilterbankIO.strtypes:
                    self.header[key] = read_string(self.file)
                    bytes_read += len(self.header[key])+4
                elif key in FilterbankIO.dbltypes:
                    self.header[key] = np.fromfile(self.file,dtype=np.float64,count=1)[0]
                    bytes_read += 8
                elif key in FilterbankIO.chrtypes:
                    self.header[key] = np.fromfile(self.file,dtype=np.int8,count=1)[0]
                    bytes_read += 1
                else:
                    raise Exception("Cannot parse filterbank header, key '{}' not understood".format(key))
                key = read_string(self.file)
                bytes_read += len(key)+4
        else:
            raise Exception("Cannot parse filterbank header, HEADER_START was not fount")
        self.state=1
        self.data_position=self.file.tell()

        
        return bytes_read
        
    def write_header(self,filename):
        def write_string(file,string):
            l=len(string)
            byte_data= struct.pack("i{}s".format(l),l,string.encode("UTF-8"))
            file.write(byte_data)
        
        
        self.file = open(filename, "wb")
        self.header_position=self.file.tell()

        write_string(self.file,"HEADER_START")
        for key in self.header:
            write_string(self.file,key)
            if key in FilterbankIO.inttypes:
                self.file.write(struct.pack("i",self.header[key]))
            elif key in FilterbankIO.strtypes:
                write_string(self.file,self.header[key])
            elif key in FilterbankIO.dbltypes:
                self.file.write(struct.pack("d",self.header[key]))
            elif key in FilterbankIO.chrtypes:
                self.file.write(struct.pack("b",self.header[key]))
            else:
                raise Exception("Cannot understand filterbank header for writing, key '{}' not understood".format(key))
        write_string(self.file,"HEADER_END")
        
        self.data_position=self.file.tell()
        self.state=2
  
    def clone_header(self):
        ret = FilterbankIO()
        ret.header = self.header.copy()
        return ret

    def seconds_since(self,sample, mjdref):
        return 86400.0*(self.header['tstart'] - mjdref) + sample * self.header['tsamp']

    def mjd_at_sample(self,sample):
        return self.header['tstart'] + (sample * self.header['tsamp'])/86400.0
    
    def frequency_table(self):
        return np.linspace(self.header['fch1'],self.header['fch1'] + self.header['nchans']*self.header['foff'],self.header['nchans'],endpoint=False)

    def seek_to_sample(self,sample):
        bytes_from_data_start = (sample * self.header['nchans'] * 8)//self.header['nbits']
        self.file.seek(self.data_position+bytes_from_data_start,os.SEEK_SET)
    
    def current_position(self):
        return ((self.file.tell()-self.data_position) * self.header['nbits']) // (8*self.header['nchans'])
        
    def total_bytes(self):
        return os.fstat(self.file.fileno()).st_size
    
    def total_channels(self):
        return self.header['nchans']
    def total_samples(self):
        header_size = self.data_position - self.header_position
        return int(8*(self.total_bytes()-header_size) / self.header['nchans'] / self.header['nbits'])
        
    def sample_interval(self):
        return self.header['tsamp']

    def channel_bandwidth(self):
        return self.header['foff']

    def get_reader(self):
        
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
            out[1::2] += np.right_shift(np.bitwise_and(raw,0x0f),4)
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
        
    def get_writer(self):
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
        if self.state>0:
            self.file.close()
            self.state=0
