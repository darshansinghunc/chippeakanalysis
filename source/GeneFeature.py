"""
Created on Oct 20, 2015
@author: darshan
"""
import sys
import os
import HTSeq
sys.path.insert(0, '..')
import common

class GeneFeature():
    """
    Finding Gene Features for a peak
    """
    def __init__(self, tssgtffile, allgtffile):
        """
        gtffile: contains information about genes in gtf format
        """
        if os.path.exists(tssgtffile):
            self.tssgtf=tssgtffile
            self.genetssdict={}
            self.genetssindex={}
        else:
            message='TSS GTF file not found'
            common.printstatus(message,'F',common.func_name())
        if os.path.exists(allgtffile):
            self.allgtf=allgtffile
            self.geneexondict={}
            self.exonindex={}
            self.genebodydict={}
            self.genebodyindex={}
        else:
            message='ALL GTF file not found'
            common.printstatus(message,'F',common.func_name())
            
    def chrintervaldistance(self,chrinterval1,chrinterval2):
        """finds distance between two intervals (chr,start,end)
        return 0 if overlap
        return 1000000000 if different chr
        else distance
        """ 
        chr1,start1,end1=chrinterval1
        chr2,start2,end2=chrinterval2
        if chr1!=chr2:
            return 1000000000
        else:
            return max(0,max(start1,start2)-min(end1,end2))
        
    def getGeneTSS(self,all_flg=1):
        """
        gets the gene TSS from the gtf file
        if all_flg==0: gets first TSS in the order of transcription 
        if ignore_weird_chr=1: only chr len <=5 are used
        returns dictionary of gene, chrposfeature (chrom,pos,strand) of tss
        """
        ignore_weird_chr=1
        
        gtffile = HTSeq.GFF_Reader(self.tssgtf)
        genetssdict={}
        for feature in gtffile:
            chrposfeature=feature.iv.start_d_as_pos
            if ignore_weird_chr and len(chrposfeature.chrom)>5:
                continue
            gene=feature.attr["gene_id"]
            if gene in genetssdict:
                genetssdict[gene].add(chrposfeature)
            else:
                genetssdict[gene]=set([chrposfeature])

        if all_flg==1:
            self.genetssdict=genetssdict
        else:
            for gene in genetssdict.keys():
                allpos=[]
                tsslist=genetssdict[gene]
                for chrposfeature in tsslist:
                    allpos.append(chrposfeature.pos)
                    pstrand=chrposfeature.strand
                if pstrand=="+":
                    usepos=min(allpos)
                else:
                    usepos=max(allpos)
                genetssdict[gene]=[p for p in tsslist if p.pos==usepos]
            self.genetssdict=genetssdict
            
    def _getGeneTSSindex(self):
        """gets index on genes as a hundred-thousandth position
        
        returns the dictionary of genes per key of 100000th position
        """
        genetssindex={}
        for gene in self.genetssdict:
            for tss in self.genetssdict[gene]:
                key='%s-%04d'%(tss.chrom,int(tss.pos/100000))
                if key in genetssindex:
                    genetssindex[key].add(gene)
                else:
                    genetssindex[key]=set([gene])
        self.genetssindex=genetssindex
        
    def peak2promotor(self,peakinterval,distance):
        """
        returns a list of genes and tss location near the peak
        """
        nearbygenesdict={}
        if len(self.genetssdict)==0:
            self.getGeneTSS()
        if len(self.genetssindex)==0:
            self._getGeneTSSindex()
        chrpeak,peakstart,peakend=peakinterval
        keyrange=[(peakstart-distance)/100000,(peakend+distance)/100000+1]
        candidategenes=[]
        for keyidx in range(keyrange[0],keyrange[1]):
            key='%s-%04d'%(chrpeak,keyidx)
            if key in self.genetssindex:
                candidategenes+=self.genetssindex[key]
        for gene in candidategenes:
            tsslist=self.genetssdict[gene]
            for tss in tsslist:
                chrom=tss.chrom
                pos=tss.pos
                chrinterval1=[chrom,pos,pos]
                tsspeakdistance=self.chrintervaldistance(chrinterval1,peakinterval)
                if tsspeakdistance<distance:
                    if (peakend<tss.pos and tss.strand=='+') or (tss.strand=='-' and peakstart>tss.pos):
                        signmultiplier=-1
                    else:
                        signmultiplier=1
                    if gene in nearbygenesdict:
                        nearbygenesdict[gene].append(['%s:%d/%s'%(tss.chrom,tss.pos,tss.strand),tsspeakdistance*signmultiplier])
                    else:
                        nearbygenesdict[gene]=[['%s:%d/%s'%(tss.chrom,tss.pos,tss.strand),tsspeakdistance*signmultiplier]]
        return nearbygenesdict
    
    def _getgeneexondict(self):
        geneexondict={}
        gtffile = HTSeq.GFF_Reader(self.allgtf)
        for feature in gtffile:
            if feature.type=='exon':
                gene=feature.name
                if len(feature.iv.chrom)>5:
                    continue
                exon=[feature.iv.chrom,feature.iv.start,feature.iv.end]
                if gene in geneexondict:
                    geneexondict[gene]+=[exon]
                else:
                    geneexondict[gene]=[exon]
        genebodydict={}  
        for gene in geneexondict:
            genechrlist=list(set([x[0] for x in geneexondict[gene]]))
            genebodydict[gene]=[]
            for chrnm in genechrlist:
                genebodydict[gene].append([chrnm,min([x[1] for x in geneexondict[gene] if x[0]==chrnm]),max([x[2] for x in geneexondict[gene] if x[0]==chrnm])])     
        self.geneexondict=geneexondict
        self.genebodydict=genebodydict
        
    def _getexonindex(self):
        exonindex={}
        for gene in self.geneexondict:
            for exon in self.geneexondict[gene]:
                for poskey in range(exon[1]/100000,exon[2]/100000+1):
                    key='%s-%04d'%(exon[0],poskey)
                    if key in exonindex:
                        exonindex[key].append([exon,gene])
                    else:
                        exonindex[key]=[[exon,gene]]
        self.exonindex=exonindex
    
    def peak2exon(self,peakinterval):
        if len(self.geneexondict)==0:
            self._getgeneexondict()
        if len(self.exonindex)==0:
            self._getexonindex()
        chrpeak,peakstart,peakend=peakinterval
        keyrange=[peakstart/100000,peakend/100000+1]
        candidateexons=[]
        for keyidx in range(keyrange[0],keyrange[1]):
            key='%s-%04d'%(chrpeak,keyidx)
            if key in self.exonindex:
                candidateexons+=self.exonindex[key]
        unique_candidateexons=[] 
        for x in candidateexons:
            if x not in unique_candidateexons:
                unique_candidateexons.append(x)
        unique_intervals=[x[0] for x in unique_candidateexons]
        overlapintervals=common.interval_join(unique_intervals,[peakinterval],5,overlappct=0.5)
        overlapexonlist=[]
        for interval in overlapintervals:
            for exon in unique_candidateexons:
                if exon[0]==interval:
                    overlapexonlist.append(exon)
                    break
        return overlapexonlist       

    def _getgenebodyindex(self):
        genebodyindex={}
        for gene in self.genebodydict:
            if gene[0:3]=='MIR':
                continue
            for genebody in self.genebodydict[gene]:
                for poskey in range(genebody[1]/100000,genebody[2]/100000+1):
                    key='%s-%04d'%(genebody[0],poskey)
                    if key in genebodyindex:
                        genebodyindex[key].append([genebody,gene])
                    else:
                        genebodyindex[key]=[[genebody,gene]]
        self.genebodyindex=genebodyindex
    
    def peak2genebody(self,peakinterval):
        if len(self.geneexondict)==0:
            self._getgeneexondict()
        if len(self.genebodyindex)==0:
            self._getgenebodyindex()
        chrpeak,peakstart,peakend=peakinterval
        keyrange=[peakstart/100000,peakend/100000+1]
        candidategenebody=[]
        for keyidx in range(keyrange[0],keyrange[1]):
            key='%s-%04d'%(chrpeak,keyidx)
            if key in self.genebodyindex:
                candidategenebody+=self.genebodyindex[key]
        unique_candidategenebody=[]
        for x in candidategenebody:
            if x not in unique_candidategenebody:
                unique_candidategenebody.append(x)
        unique_intervals=[x[0] for x in unique_candidategenebody]
        overlapintervals=common.interval_join(unique_intervals,[peakinterval],4,overlappct=0.5)
        overlapgenebody=[]
        for interval in overlapintervals:
            for genebody in unique_candidategenebody:
                if genebody[0]==interval:
                    overlapgenebody.append(genebody)
                    break        
        return overlapgenebody         
    
    def peakto3prime(self,peakinterval,bodypeakdistancethreshold):
        genesizeconstant=4000000
        if len(self.genetssdict)==0:
            self.getGeneTSS()
        if len(self.genetssindex)==0:
            self._getGeneTSSindex()            
        if len(self.geneexondict)==0:
            self._getgeneexondict()
        if len(self.genebodyindex)==0:
            self._getgenebodyindex()
        chrpeak,peakstart,peakend=peakinterval
        keyrange=[(peakstart-genesizeconstant)/100000,(peakend+genesizeconstant)/100000+1]
        candidategenes=[]
        for keyidx in range(max(keyrange[0],0),keyrange[1]):
            key='%s-%04d'%(chrpeak,keyidx)
            if key in self.genetssindex:
                candidategenes+=self.genetssindex[key]
        candidategenes=list(set(candidategenes))
        threeprimegenedict=[]
        for gene in candidategenes:
            tss=list(self.genetssdict[gene])[0]
            if tss.strand=='-':
                multiplier=-1
            else:
                multiplier=1
            bodyinterval=self.genebodydict[gene][0]
            bodypeakdistance=self.chrintervaldistance(bodyinterval,peakinterval)
            if bodypeakdistance>0 and bodypeakdistance<bodypeakdistancethreshold:
                #print gene, bodyinterval, peakinterval, bodypeakdistance, multiplier
                if (peakinterval[1]-bodyinterval[2])*multiplier>0:
                    if multiplier>0:
                        lastpos=bodyinterval[2]
                    else:
                        lastpos=bodyinterval[1]
                    threeprimegenedict.append([gene,lastpos,bodypeakdistance])
        return threeprimegenedict
                    
    