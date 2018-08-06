from BinaryFile import BinaryFile,BinaryCDataFile
import ctypes

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
        ("anchor_length1",ctypes.c_uint64,9),
        ("anchor_length2",ctypes.c_uint64,9),
        ("mate_anchor1_length",ctypes.c_uint64,9),
        ("mate_anchor2_length",ctypes.c_uint64,9)
    ]


class BridgeDataFile(BinaryCDataFile):
    def __init__(self,bridgeFile,fileAccess="mmap"):
        BinaryCDataFile.__init__(self,bridgeFile,Bridge_CData,fileAccess)
