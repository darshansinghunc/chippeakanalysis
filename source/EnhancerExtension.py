import common
import ChipSeqIndex as csi
import EnhancerPipeline as ep
from ChipSeqIndex import ChipSeqIndex
import json
import math
import scipy.stats
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch
import numpy as np
from matplotlib import rcParams
#from superenhancers import peaksummarydict
rcParams.update({'figure.autolayout': True})

def upanddown(valuelist):
    if valuelist[1]>0.2:
        return valuelist[1]*valuelist[1]/(max(valuelist[2],0.05)*max(valuelist[0],0.05))
    else:
        return 0

def upanddownminusinput(valuelist):
    if valuelist[1]<valuelist[0] or valuelist[1]<valuelist[2]:
        return 0
    if valuelist[1]>0.2:
        return (valuelist[1]-valuelist[3])*(valuelist[1]-valuelist[3])/(max(valuelist[2]-valuelist[3],0.05)*max(valuelist[0]-valuelist[3],0.05))
    else:
        return 0
    
def upandupminusinput(valuelist):
    if valuelist[1]<valuelist[0] or valuelist[1]>valuelist[2]:
        return 0
    if valuelist[1]>0.2:
        return (valuelist[2]-valuelist[3])/max(valuelist[0]-valuelist[3],0.05)
    else:
        return 0    

def upminusinput(valuelist):
    if valuelist[1]<valuelist[0]:
        return 0
    if valuelist[1]>0.2:
        return (valuelist[1]-valuelist[2])/max(valuelist[0]-valuelist[2],0.05)
    else:
        return 0
    
def downminusinput(valuelist):
    if valuelist[1]>valuelist[0]:
        return 0
    if valuelist[1]>0.2:
        return (valuelist[0]-valuelist[2])/max(valuelist[1]-valuelist[2],0.05)
    else:
        return 0    

def topchipsubregions(inreport, columnlist, maxexpressionfunc, numberofregions, sizeofregion, slidingsize, index_file_parameters):
    peakobjectivedict={}
    lncnt=0
    columnidxlist=[-1 for x in columnlist]
    for lntxt in open(inreport):
        ln=lntxt.rstrip('\n').split('\t')
        lncnt+=1
        if lncnt==1:
            for columnidx in range(len(columnlist)):
                for idx in range(len(ln)):
                    if columnlist[columnidx]==ln[idx]:
                        columnidxlist[columnidx]=idx
            if len([x for x in columnidxlist if x>-1])!=len(columnlist):
                message="Some columns are not found in the report"
                funcname="topchipsubregions"
                print common.printstatus(message, 'W',funcname)
        else:
            valuelist=[]
            for idx in columnidxlist:
                if idx==-1:
                    valuelist.append(0.0)
                else:
                    valuelist.append(float(ln[idx]))
            objectivevalue=maxexpressionfunc(valuelist)
            peakobjectivedict[ln[0]]=objectivevalue
    objectivevaluelist=peakobjectivedict.values()
    objectivevaluelist.sort()
    valuethreshold=objectivevaluelist[max(0,len(objectivevaluelist)-5*numberofregions)]
    print 'Screening Cutoff:',len(objectivevaluelist),valuethreshold
    candidatepeaklist=[x for x in peakobjectivedict.keys() if peakobjectivedict[x]>valuethreshold]
    columncsidict={}
    columnreadcntdict={}
    for column in columnlist:
        bamfile=ep.getcodetofilename(index_file_parameters, column)
        columncsidict[column]=ChipSeqIndex(bamfile)
        columnreadcntdict[column]=ep.getmappedreadcount(bamfile)
    topregiondict={}
    for candidatepeak in candidatepeaklist:
        peakinterval=[candidatepeak.split(':')[0],int(candidatepeak.split(':')[1].split('-')[0]),int(candidatepeak.split(':')[1].split('-')[1])]
        numwindows=(max(peakinterval[2]-peakinterval[1]-sizeofregion,slidingsize))/slidingsize
        for i in range(numwindows):
            candidateinterval='%s:%d-%d'%(peakinterval[0],peakinterval[1]+i*slidingsize,min(peakinterval[1]+i*slidingsize+sizeofregion,peakinterval[2]))
            valuelist=[]
            for column in columnlist:
                numreads=columncsidict[column].countreadsinrange(candidateinterval)
                valuelist.append(numreads*1000000.0/columnreadcntdict[column])
            objectivevalue=maxexpressionfunc(valuelist)
            topregiondict[candidateinterval]=objectivevalue
    topvalueintervallist=[(topregiondict[x],x) for x in topregiondict.keys()]
    topvalueintervallist.sort()
    topvalueintervallist.reverse()
    topintervallist=[]
    for topvalueinterval in topvalueintervallist:
        intervalstr=topvalueinterval[1]
        interval=[intervalstr.split(':')[0],int(intervalstr.split(':')[1].split('-')[0]),int(intervalstr.split(':')[1].split('-')[1])]
        if len(common.interval_join(topintervallist,[interval],0))==0:
            topintervallist.append(interval)
        if len(topintervallist)==numberofregions:
            break
    print 'Final Cutoff:',len(topvalueintervallist),topvalueinterval[0]
    return topintervallist
    
def mergeChipRNA(injson,RNAcolumnlist,Chipcolumnlist,index_file_parameters,outfile):
    """
    denormalizers the peak and gene and adds columns for RNA-seq
    """
    rnasamplegeneexpressiondict={}
    for rna_column in RNAcolumnlist:
        rnasamplegeneexpressiondict[rna_column]={}
        rna_filename=ep.getcodetofilename(index_file_parameters,rna_column)
        for lntxt in open(rna_filename):
            ln=lntxt.rstrip('\n').split(',')
            if ln[0]=='gene_id':
                continue
            gene=ln[0].split('|')[0]
            expression=float(ln[-1])
            rnasamplegeneexpressiondict[rna_column][gene]='%6.4f'%expression
    peak_attr_dict=json.load(open(injson))
    fout=open(outfile,'w')
    fout.write("peak\tgene\tnearest gene tss\tpeak gene tss distance")
    for column in Chipcolumnlist:
        fout.write("\tchip_%s"%column)
    for column in RNAcolumnlist:
        fout.write("\trna_%s"%column)    
    fout.write('\n')
    for peak in peak_attr_dict:
        chip_expression=['%6.2f'%peak_attr_dict[peak]['chip_coverage'][column]['RPM'] for column in Chipcolumnlist]
        for gene in peak_attr_dict[peak]['peak_attr']['nearby_genes']:
            min_dist=min([abs(x[1]) for x in peak_attr_dict[peak]['peak_attr']['nearby_genes'][gene]])
            nearest_tss=[x[0] for x in peak_attr_dict[peak]['peak_attr']['nearby_genes'][gene] if abs(x[1])==min_dist][0]
            min_dist=[x[1] for x in peak_attr_dict[peak]['peak_attr']['nearby_genes'][gene] if x[0]==nearest_tss][0]
            rna_expr_str='\t'.join([rnasamplegeneexpressiondict[rna_column].get(gene,'0.0000') for rna_column in RNAcolumnlist])
            fout.write('%s\t%s\t%s\t%d\t%s\t%s\n'%(peak,gene,nearest_tss,min_dist,'\t'.join(chip_expression),rna_expr_str))
    fout.close()
            
            
def testChipReplicate(chip_id1,chip_id2,index_file_parameters,blocksize,chromsize_file,imgfile,peak_file):
    chromsizedict=dict([(x.rstrip().split('\t')[0],int(x.rstrip().split('\t')[1])) for x in open(chromsize_file).readlines()])
    allblocklist=[]
    for chrnm in chromsizedict:
        for i in range(int(chromsizedict[chrnm]/blocksize)+1):
            allblocklist.append([chrnm,i*blocksize,min((i+1)*blocksize,chromsizedict[chrnm])])
    peaklist=[]
    lncnt=0
    for lntxt in open(peak_file):
        lncnt+=1
        if lncnt==1:
            continue
        peak=lntxt.split('\t')[0]
        peaklist.append([peak.split(':')[0],int(peak.split(':')[1].split('-')[0]),int(peak.split(':')[1].split('-')[1])])
        
    candidateintvlist=common.interval_join(allblocklist,peaklist,3,overlappct=0.1)
    chip_bam1=ep.getcodetofilename(index_file_parameters,chip_id1)
    chip_bam2=ep.getcodetofilename(index_file_parameters,chip_id2)
    chip1=ChipSeqIndex(chip_bam1)
    bam1count=ep.getmappedreadcount(chip_bam1)
    chip2=ChipSeqIndex(chip_bam2)
    bam2count=ep.getmappedreadcount(chip_bam2)   
    pointxlist=[] ; pointylist=[] 
    morethan2cnt=0
    for candidateintvtup in candidateintvlist:
            candidateinterval='%s:%d-%d'%(candidateintvtup[0],candidateintvtup[1],candidateintvtup[2])
            S1=math.log10((chip1.countreadsinrange(candidateinterval)+0.01)*1000000.0/bam1count)
            S2=math.log10((chip2.countreadsinrange(candidateinterval)+0.01)*1000000.0/bam2count)
            if S1+S2>-2:
                morethan2cnt+=1
            pointxlist.append((S1+S2)/2.0)
            pointylist.append(S1-S2)
    print len(pointxlist), morethan2cnt
    diffmean=np.mean(pointylist)
    diffstd=np.std(pointylist)
    print diffmean,diffstd
    cntoutside=len([x for x in pointylist if x<diffmean-1.96*diffstd or x>diffmean+1.96*diffstd])
    outpct= cntoutside*1.0/len(pointylist)
    fig = plt.figure(figsize=(8,8))
    plt.plot(pointxlist,pointylist,'k.')
    plt.plot([-4,4],[diffmean,diffmean],'b-')   
    plt.plot([-4,4],[diffmean+1.96*diffstd,diffmean+1.96*diffstd],'r-') 
    plt.plot([-4,4],[diffmean-1.96*diffstd,diffmean-1.96*diffstd],'r-') 
    plt.axis((-4,4,-2,2))
    plt.title('%s vs %s\npercentage outside limits %5.2f'%(chip_id1,chip_id2,outpct*100))
    plt.savefig('%s__%s_%s.png'%(imgfile,chip_id1,chip_id2))

def foldchangecolumnupdate(csvfilename,minvalue,compare_column_list,reference_column,outfile,filtertuple=(0.1,2)):
    """
    filtertup=(x,y): count y more than x
    """
    matrix=[]
    fout=open(outfile,'w')
    lncnt=0
    columnindexdict={}
    for lntxt in open(csvfilename):
        lncnt+=1
        ln=lntxt.rstrip('\n').split('\t')
        if lncnt==1:
            for col in compare_column_list:
                if col in ln:
                    columnindexdict[col]=ln.index(col)
                else:
                    print col
                    print ln
                    message="Some columns are not found in the report"
                    funcname="foldchangecolumnupdate"
                    print common.printstatus(message, 'F',funcname)
                if reference_column in ln:
                    columnindexdict[reference_column]=ln.index(reference_column)
                else:
                    message="Reference column is not found in the report"
                    funcname="foldchangecolumnupdate"
                    print common.printstatus(message, 'F',funcname)                    
            fout.write('%s\t%s\t%s\n'%('\t'.join(ln[0:4]),'\t'.join(compare_column_list),'\t'.join(ln[-5:])))
        else:
            logfoldchangelist=[]
            rpmlist=[]
            refrpm=float(ln[columnindexdict[reference_column]])
            rpmlist.append(refrpm)
            for col in compare_column_list:
                colrpm=float(ln[columnindexdict[col]])
                logfoldchange=math.log(max(colrpm,minvalue)/max(refrpm,minvalue),2)
                logfoldchangelist.append('%6.2f'%logfoldchange)
                rpmlist.append(colrpm)
            rpmthreshold,cntthreshold=filtertuple
            if len([x for x in rpmlist if x>rpmthreshold])>=cntthreshold:
                fout.write('%s\t%s\t%s\n'%('\t'.join(ln[0:4]),'\t'.join(logfoldchangelist),'\t'.join(ln[-5:])))
                matrix.append([float(x) for x in logfoldchangelist])
    return (compare_column_list,matrix)       
            
def matrixprocessing(labels,matrix,dist_name,imgfile):
    npmatrix=np.array(matrix)
    txnpmatrix=np.transpose(npmatrix)
    D=sch.distance.pdist(txnpmatrix,dist_name)
    fig = plt.figure(figsize=(8,8))
    Z= sch.linkage(D,method='complete')
    P =sch.dendrogram(Z, labels=labels,leaf_rotation=90)
    plt.savefig('%s_%s.pdf'%(imgfile,dist_name))
    return D
            
def twosamplepeakcomparison(project1_name, sample1, project2_name, sample2):
    """
    Part of the same project
    """
    project1root='/media/ddrive/project/%s'%project1_name
    project2root='/media/ddrive/project/%s'%project2_name
    peak1dict=json.load(open('%s/output/peaks/%s__%s__sample__peaks_attr.json'%(project1root,project1_name,sample1)))
    peak2dict=json.load(open('%s/output/peaks/%s__%s__peaks.json'%(project2root,project2_name,sample2)))
    peak1list=[[x.split(':')[0],int(x.split(':')[1].split('-')[0]),int(x.split(':')[1].split('-')[1])] for x in peak1dict.keys()]
    peak2list=[[x.split(':')[0],int(x.split(':')[1].split('-')[0]),int(x.split(':')[1].split('-')[1])] for x in peak2dict.keys()]
    
    peak1list.sort()
    peak2list.sort()
 
    overlaplist=common.interval_join2(peak1list,peak2list,3,overlappct=0.01)
    
    extrafields=['compare_peak','overlapsize']
    peak_attr_dict=peak1dict.copy()
    for interval in overlaplist:
        print interval
        peak1=[interval[0],interval[1],interval[2]]
        peak2=[interval[0],interval[3],interval[4]]
        overlapsize=(min(interval[2],interval[4])-max(interval[1],interval[3])+1)*1.0/(interval[2]-interval[1]+1)
        peak_key='%s:%d-%d'%(peak1[0],peak1[1],peak1[2])
        peak_attr_dict[peak_key]['extra']={}
        peak_attr_dict[peak_key]['extra']['compare_peak']='%s:%d-%d'%(peak2[0],peak2[1],peak2[2])
        peak_attr_dict[peak_key]['extra']['overlapsize']='%5.2f'%overlapsize
    
    samplelist=[sample1]
    reportfilterdict={'target':'%s'%sample1,'thresholds':{},'constraints':{}}
    csvfile='%s/output/report/%s_overlap_%s_%s.csv'%(project1root,project1_name, sample1, sample2)
    ep.preparecsvanalysis(peak_attr_dict,samplelist,reportfilterdict,csvfile,extrafields)


            
    
     
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
        
    
            
                
                        