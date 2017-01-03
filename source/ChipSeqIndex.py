import os
import json
import bisect
import hashlib

class ChipSeqIndex:
    '''Class to index and retrieve data from Chip_Seq files
    '''
    def __init__(self, filename):
        self.identity=hashlib.md5(filename).hexdigest()
        self.name=filename
        self.dir='/'.join(filename.split('/')[:-1])
        chrlist=[x.split('\t')[0] for x in os.popen('samtools idxstats %s'%self.name).readlines()]
        if '*' in chrlist:
            chrlist.remove('*')
        chrlist.sort()
        self.CHRMULT=1000000000
        self.chrmap1=dict(zip(chrlist,range(51,51+len(chrlist))))
        self.chrmap2=dict(zip(range(51,51+len(chrlist)),chrlist))
        self.idxname=filename[0:-4]+'_chip.idx' 
        self.chrposdict={}
        if not self.checkindex():
            self.buildindex()
    
    def checkindex(self):
        if os.path.isfile(self.idxname):
            if os.path.getmtime(self.idxname) > os.path.getmtime(self.name):
                return 1
        return 0
    
    def buildindex(self):
        """
        convert the bam to readlist of chrpos to integer
        save the list
        """
        idxdir='/tmp/%s_%s_chip'%(self.name.split('/')[-1].split('.')[0],self.identity)
        if not os.path.exists(idxdir):
            os.system('mkdir %s'%(idxdir))
        chrlist=self.chrmap1.keys()
        for chrnm in chrlist:
            chrposlist=[]
            for lntxt in os.popen('samtools view %s %s'%(self.name,chrnm)):
                ln=lntxt.split('\t')
                chrpos=[ln[2],int(ln[3])]
                chrposlist.append(self._chrpos2int(chrpos))
            json.dump(chrposlist,open('%s/%s_data.json'%(idxdir,chrnm),'w'))
        os.chdir(self.dir)
        os.system('tar -cf %s.tar %s'%(self.name,idxdir))
        os.system('mv %s.tar %s'%(self.name,self.idxname))
        os.system('rm -r %s'%(idxdir))
        return 0
    
    def _chrpos2int(self,chrpos):
        """
        chromosome mapped to integer and attached to position
        chromosome names in sorted order from samtools idxstats
        append to position
        e.g.
            chrpos=[chr1,3201912]
        """
        chrnm=chrpos[0]
        pos=chrpos[1]
        return self.chrmap1[chrnm]*self.CHRMULT+pos
    
    def _int2chrpos(self,chrposint):
        """
        extract position, convert chromosome
        """
        chrnm=self.chrmap2[chrposint/self.CHRMULT]
        pos=chrposint % self.CHRMULT
        return [chrnm,pos]
    
    def _getindexesinrange(self,chrrange):
        """
        Fix: the name of the object and test existance before loading
        """
        idxdir='/tmp/%s_%s_chip'%(self.name.split('/')[-1].split('.')[0],self.identity)
        chrnm=chrrange.split(':')[0]
        pos1=int(chrrange.split(':')[1].split('-')[0])
        pos2=int(chrrange.split(':')[1].split('-')[1])
        if chrnm in self.chrposdict:
            chrposlist=self.chrposdict[chrnm]
        else:
            if not os.path.exists(idxdir):
                os.system('mkdir %s'%idxdir)
                print 'tar xf %s'%(self.idxname)
                os.system('cd / && tar xvf %s'%(self.idxname))
            chrjson='%s/%s_data.json'%(idxdir,chrnm)
            if os.path.exists(chrjson):
                chrposlist=json.load(open(chrjson))
                self.chrposdict[chrnm]=chrposlist
            else:
                return ([],0,0)
        chrpos1int=self._chrpos2int([chrnm,pos1])
        chrpos2int=self._chrpos2int([chrnm,pos2])
        leftidx=bisect.bisect(chrposlist,chrpos1int-1)
        rightidx=bisect.bisect(chrposlist,chrpos2int)
        return (chrposlist,leftidx,rightidx)
    
    def listreadsinrange(self,chrrange):
        chrposlist,leftidx,rightidx=self._getindexesinrange(chrrange)
        readlist=[]
        for chrposint in chrposlist[leftidx:rightidx]:
            readlist.append(self._int2chrpos(chrposint))
        return readlist
    
    def countreadsinrange(self,chrrange):
        _,leftidx,rightidx=self._getindexesinrange(chrrange)
        return rightidx-leftidx
        
        
        
    