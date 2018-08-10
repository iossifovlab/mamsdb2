from BinaryFile import BinaryFile,BinaryCDataFile
import ctypes
import os

class lazy_property(object):
    def __init__(self,fget):
        self.fget = fget
        self.func_name = fget.__name__
    def __get__(self,obj,cls):
        if obj is None:
            return None
        value = self.fget(obj)
        setattr(obj,self.func_name,value)
        return value


class Bridge_CData(ctypes.Structure):
    _fields_=[
        ("pos1",ctypes.c_uint64,31),
        ("pos2",ctypes.c_uint64,31),
        ("high1",ctypes.c_uint64,1),
        ("high2",ctypes.c_uint64,1),
        ("chr1",ctypes.c_uint64,8),
        ("chr2",ctypes.c_uint64,8),
        ("padding",ctypes.c_uint64,1),
        ("offset",ctypes.c_int64,11),
        ("anchor1_length",ctypes.c_uint64,9),
        ("anchor2_length",ctypes.c_uint64,9),
        ("mate_anchor1_length",ctypes.c_uint64,9),
        ("mate_anchor2_length",ctypes.c_uint64,9)
    ]

    @lazy_property
    def invariant(self):
        return (1 if self.high1 else -1)*self.pos1+(1 if self.high2 else -1)*self.pos2 + self.offset

class BridgeDB():
    def __init__(self,bridgesDir,fileAccess="mmap"):
        self.bridgesDir=bridgesDir
        self.fileAccess=fileAccess

    def getBridgesByChrom(self,chrom):
        chrom=chrom.lower()
        bridgeFilePath=os.path.join(self.bridgesDir,"newbridges.%s.bin" % (chrom))
        
        bridgeFile=BinaryCDataFile(bridgeFilePath,Bridge_CData,self.fileAccess)
        
        for i in range(0,bridgeFile.numRecords):
            yield bridgeFile.readIndex(i)
        
