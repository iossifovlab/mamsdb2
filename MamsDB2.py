import os
import struct
import array

class BinaryFile(object):
    def __init__(self,fileName,fileAccess="file"):
        self.fileName=fileName
        self.fileAccess=fileAccess
        self._file=open(fileName,'rb')
        if fileAccess=="memory":
            self._data=bytearray(os.path.getsize(fileName))
            self._file.readinto(self._data)
            self._file.close()

    def readIndex(pos):
        if fileAccess=="memory":
            return self_.data[pos]
        elif fileAccess=="file":
            self._file.seek(pos)
            self._file.read(1)

    def close():
        if fileAccess=="memory":
            self._data=None
        elif fileAccess=="file":
            self._file.close()
        
            
        
    

class Mappability(object):
    def __init__(self,fasta_name,fileAccess="file"):
        self._lowFN=fasta_name+".bin/map.low.bin"
        self._highFN=fasta_name+".bin/map.high.bin"
        if os.path.getsize(self._lowFN)!=os.path.getsize(self._highFN):
           raise Error("The low and high Mappability files should be the same size")

        self._lowF=BinaryFile(self._lowFN,fileAccess)
        self._highF=BinaryFile(self._highFN,fileAccess)
        self.fasta=fasta_name

    def low(pos):
        return struct.unpack("B",self._lowF.readIndex(pos))

    def high(pos):
        return struct.unpack("B",self._highF.readIndex(pos))

    @staticmethod
    def createFromMumdexDir(mumdexDir):
        fasta_name=""
        with open(os.path.join(mumdexDir,"ref.txt")) as refF:
            fasta_name=refF.readline().strip()

        return Mappability(fasta_name)
    
        
        
class Reference(object):
    def __init__(self,fasta_name,fileAccess="file"):
        self._seqFN=fasta_name+".bin/ref.seq.bin"
        self._seqF=BinaryFile(self._seqFN,fileAccess)
        self._chrLenFN=fasta_name+".bin/ref.chr_len.bin"
        self._chrNameFN=fasta_name+".bin/ref.chr_name.bin"
        self.chrLen=array.array("I")
        with open(self._chrLenFN) as chrLenF:
            self.chrLen.fromstring(chrLenF.read())

        self.chrName=[]
        with open(self._chrNameFN) as chrNameF:
            for line in chrNameF:
                self.chrName.append(line.strip())

        self.fasta=fasta_name

    @staticmethod
    def createFromMumdexDir(mumdexDir):
        fasta_name=""
        with open(os.path.join(mumdexDir,"ref.txt")) as refF:
            fasta_name=refF.readline().strip()

        return Reference(fasta_name)
        
    def name(self,chromIndex):
        return self.chrName[chromIndex]
        
