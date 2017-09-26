import os
import struct

class Mappability(object):
    def __init__(self,fasta_name,fileAccess="memory"):
        self._lowFN=fasta_name+".bin/map.low.bin"
        self._highFN=fasta_name+".bin/map.high.bin"
        if os.path.getsize(self._lowFN)!=os.path.getsize(self._highFN):
           raise Error("The low and high Mappability files should be the same size")

        self._lowF=open(self._lowFN)
        self._highF=open(self._highFN)        
        if fileAccess=="memory":
            self._lowFMem=bytearray(os.path.getsize(self._lowFN))
            self._lowF.readinto(self._lowFMem)
            self._lowF.close();

            self._highFMem=bytearray(os.path.getsize(self._highFN))
            self._highF.readinto(self._highFMem)
            self._highF.close();

    def low(pos):
        return struct.unpack("B",self._lowFMem[pos])

    def high(pos):
        return struct.unpack("B",self._highFMem[pos])
        
        
        
