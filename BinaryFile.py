import os,mmap,ctypes

class BinaryFile(object):
    def __init__(self,fileName,fileAccess="file"):
        self.fileName=fileName
        self.fileAccess=fileAccess
        self._file=open(fileName,'rb')
        if fileAccess=="memory":
            # create an anonymous mmap so the data can be shared across multiple processes
            self._data=mmap.mmap(-1,os.path.getsize(self.fileName))
            self._file.readinto(self._data)
            self._file.close()

        elif fileAccess=="mmap":
            self._data=mmap.mmap(self._file.fileno(),0,prot=mmap.PROT_READ)

    def readIndex(self,pos):
        if self.fileAccess=="memory" or self.fileAccess=="mmap":
            return self_.data[pos]
        elif self.fileAccess=="file":
            self._file.seek(pos)
            return self._file.read(1)

    def readRange(self,start,end):        
        if self.fileAccess=="memory" or self.fileAccess=="mmap":
            return self_.data[start:end]
        elif self.fileAccess=="file":
            self._file.seek(start)
            return self._file.read(end-start)

    def close(self):
        if self.fileAccess=="memory":
            self._data=None
        elif self.fileAccess=="file":            
            self._file.close()
        elif self.fileAccess=="mmap":
            self._data.close()
            self._file.close()

class BinaryCDataFile(BinaryFile):
    '''
    A binary data file where each record is defined by a ctypes Structure class
    '''
    def __init__(self,fileName,cDataClass,fileAccess="file"):
        BinaryFile.__init__(self,fileName,fileAccess)
        self.cDataClass=cDataClass
        self.sizeOfRecord=ctypes.sizeof(self.cDataClass)
        self.numRecords=os.path.getsize(self.fileName)/ctypes.sizeof(self.cDataClass)

    def readIndex(self,index):
        pos=self.sizeOfRecord*index
        if self.fileAccess=="memory":
            data=self.cDataClass.from_buffer(self._data,pos)
        elif self.fileAccess=="file" or self.fileAccess=="mmap":
            data=self.cDataClass()       
            self._file.seek(pos)
            self._file.readinto(data)            

        return data

    def readRange(self,startIndex,endIndex):
        records=[]
        for i in xrange(startIndex,endIndex):
            records.append(self.readIndex(i))
        return records
        
