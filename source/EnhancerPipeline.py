#!/usr/bin/python

import os
import time
import ChipSeqIndex
import json
from json import encoder
encoder.FLOAT_REPR = lambda o: format(o, '.2f')
import common
import GeneFeature
import HTSeq
import numpy
from matplotlib import pyplot
import bisect
import random

# PLOTs
def CHIPseq_averageTSSplots(tsslist,bamfiledict,halfwinwidth,imgfile,readfragmentsize): 
    for sampleid in bamfiledict:
        bamfilename=bamfiledict[sampleid]
        flog.write('%s: Processing %s\n'%(time.asctime(),bamfilename))
        multiplier=1000000.0/(getmappedreadcount(bamfilename)*len(tsslist))
        sortedbamfile = HTSeq.BAM_Reader(bamfilename)
        profile = numpy.zeros( 2*halfwinwidth, dtype='i' )   
        for tss in tsslist:
            tss=tss[0]
            window =  HTSeq.GenomicInterval( tss.chrom, tss.pos - halfwinwidth - readfragmentsize, tss.pos + halfwinwidth + readfragmentsize, "." )
            bamstream=sortedbamfile[ window ]
            while True:
                try:
                    almnt = bamstream.next() 
                    almnt.iv.length = readfragmentsize
                    if tss.strand == "+":
                        start_in_window = almnt.iv.start - tss.pos + halfwinwidth 
                        end_in_window   = almnt.iv.end   - tss.pos + halfwinwidth 
                    else:
                        start_in_window = tss.pos + halfwinwidth - almnt.iv.end
                        end_in_window   = tss.pos + halfwinwidth - almnt.iv.start
                    start_in_window = max( start_in_window, 0 )
                    end_in_window = min( end_in_window, 2*halfwinwidth )
                    if start_in_window >= 2*halfwinwidth or end_in_window < 0:
                        continue
                    profile[ start_in_window : end_in_window ] += 1
                except StopIteration:
                    break
                except Exception:
                    pass
        pyplot.plot( numpy.arange( -halfwinwidth, halfwinwidth ), profile*multiplier, label=sampleid)   
    pyplot.legend(prop={'size':9})
    pyplot.savefig(imgfile) 
    pyplot.clf()

def CHIPseq_geneplots(gene,genetss,bamfiledict,halfwinwidth,imgfile,readfragmentsize): 
    """Plots chip-seq coverage data for given bamfiles 
    
    Input:
        gene: Around the gene TSS
    """
    profilelist=[]
    samplelist=bamfiledict.keys()
    samplelist.sort()
    for sampleid in samplelist:
        bamfilename=bamfiledict[sampleid]
        multiplier=1000000.0/getmappedreadcount(bamfilename)
        sortedbamfile = HTSeq.BAM_Reader(bamfilename)
        profile = numpy.zeros( 2*halfwinwidth, dtype='i' )   
        p=genetss[0]
        window =  HTSeq.GenomicInterval( p.chrom, p.pos - halfwinwidth - readfragmentsize, p.pos + halfwinwidth + readfragmentsize, "." )
        bamstream=sortedbamfile[ window ]
        while True:
            try:
                almnt = bamstream.next() 
                almnt.iv.length = readfragmentsize
                if p.strand == "+":
                    start_in_window = almnt.iv.start - p.pos + halfwinwidth 
                    end_in_window   = almnt.iv.end   - p.pos + halfwinwidth 
                else:
                    start_in_window = p.pos + halfwinwidth - almnt.iv.end
                    end_in_window   = p.pos + halfwinwidth - almnt.iv.start
                start_in_window = max( start_in_window, 0 )
                end_in_window = min( end_in_window, 2*halfwinwidth )
                if start_in_window >= 2*halfwinwidth or end_in_window < 0:
                    continue
                profile[ start_in_window : end_in_window ] += 1
            except StopIteration:
                break
            except Exception:
                #print e
                pass
        profilelist.append(profile*multiplier)
    numplots=len(bamfiledict)
    fig = pyplot.figure(figsize=(8,2*len(profilelist)))
    fig.subplots_adjust(hspace=1)
    fileidx=0
    max_coverage=max([max(profile) for profile in profilelist])
    y_limit=int(max_coverage)+1
    
    for profile in profilelist:
        pyplot.subplot(numplots,1,fileidx+1)
        pyplot.ylim([0,y_limit])      
        pyplot.title(samplelist[fileidx])
        pyplot.plot( numpy.arange( -halfwinwidth, halfwinwidth ), profile)   
        fileidx+=1
    
    pyplot.savefig(imgfile) 
    pyplot.close(fig)        


def plot_superenhancer(sampleid,peak_attr_dict,imgfile,annotated_genelist,super_cutoff):
    ignore_promoter_cnt=3
    
    chip_sig_y=[]
    annotated_gene_dict={}
    for peak in sorted(peak_attr_dict.keys()):
        if len(peak_attr_dict[peak]['peak_attr']['promotor'])>=ignore_promoter_cnt:
            continue
        chip_signal=max(0.0,peak_attr_dict[peak]['chip_coverage'][sampleid]['RPM'])
        chip_sig_y.append(chip_signal)
        for gene in peak_attr_dict[peak]['peak_attr']['nearby_genes']:
            if gene in annotated_genelist:
                if gene not in peak_attr_dict[peak]['peak_attr']['promotor']:
                    annotated_gene_dict[gene]=max(chip_signal,annotated_gene_dict.get(gene,0.0))

    chip_sig_y.sort()
                    
    fig = pyplot.figure()
    ax = fig.add_subplot(111)
    ax.plot(range(len(chip_sig_y)), chip_sig_y, lw=2)
    idx=bisect.bisect(chip_sig_y,super_cutoff)
    num_superenhancers=len(chip_sig_y)-idx
    ax.plot(idx-1,chip_sig_y[idx-1],'*')
    max_x=len(chip_sig_y)
    max_y=2000
    #max_y=int(chip_sig_y[-1])
    ax.set_ylim([-max_y/20,21*max_y/20])
    ax.set_xlim([-max_x/20,21*max_x/20])  
    ax.set_xlabel('Enhancers ranked by increasing chip-seq signal')
    ax.set_ylabel('magnitude of chip-seq signal')
    ax.set_title('%s'%sampleid)    
    ax.text(11*max_x/10, max_y/2, "%d Super-Enhancers"%num_superenhancers, ha="right", va="center", rotation=270)
    for gene in annotated_gene_dict:
        y_val=annotated_gene_dict[gene]
        xy_val=(bisect.bisect(chip_sig_y,y_val),y_val)
        ax.annotate(gene, xy=xy_val, xytext=([xy_val[0]-random.randint(1,2)*max_x/10,xy_val[1]+random.randint(1,2)*max_y/15]),
             arrowprops=dict(arrowstyle="->"))
    pyplot.savefig(imgfile)
    pyplot.close("all")

#########

def preparecsvanalysis(peak_attr_dict,samplelist,reportfilterdict,csvfile,extrafields=[]):
    """
    analysisdict keys:
        targettreatment: Treatment around which the analysis moves
        targettreatmentrange: Threshold of coverage per million for target
        dict of constraints:
            name as key and list of ranges for treatments
    e.g.
        analysisdict={'target':'SUM159BRD_GSK','targetthreshold':5,'constraints':{'SUM159BRD_DMSO':[0,0.3],'SUM159BRD_COMB':[0,0.6],'SUM159K27AC_DMSO':[0.0.5]}}
    """ 
    print csvfile
    fout=open(csvfile,'w')
    fout.write('Peak\tPeaksize\tClassification\tClassification gene')
    fout.write('\t%s'%'\t'.join(samplelist))
    fout.write('\tPromotor\tEnhancer\tnearby Kinases\tnearby Transcription Factors\tnearby genes\t%s\n'%'\t'.join(extrafields))    
    for peak in peak_attr_dict:
        unselectflag=0
        targetid=reportfilterdict['target']
        targetRPM=peak_attr_dict[peak]['chip_coverage'][targetid]['RPM']
        thresholddict=reportfilterdict['thresholds']
        for sampleid in thresholddict:
            threshold=thresholddict[sampleid]
            sampleRPM=peak_attr_dict[peak]['chip_coverage'][sampleid]['RPM']
            if sampleRPM<threshold:
                unselectflag=1
                continue
        if unselectflag==1:
            continue
        for sampleid in reportfilterdict['constraints']:
            constraintrange=reportfilterdict['constraints'][sampleid]
            sampleRPM=peak_attr_dict[peak]['chip_coverage'][sampleid]['RPM']
            sampletargetratio=sampleRPM/targetRPM
            if sampletargetratio<constraintrange[0] or sampletargetratio>constraintrange[1]:
                unselectflag=1
                continue
        if unselectflag==1:
            continue
        else:
            peaksize=int(peak.split(':')[1].split('-')[1])-int(peak.split(':')[1].split('-')[0])
            peakclassification=peak_attr_dict[peak]['peak_attr']['classification']
            if 'classification_gene' not in peak_attr_dict[peak]['peak_attr']:
                peak_attr_dict[peak]['peak_attr']['classification_gene']=[]
            fout.write('%s\t%d\t%s\t%s'%(peak,peaksize,peakclassification,'|'.join(peak_attr_dict[peak]['peak_attr']['classification_gene'])))
            for sampleid in samplelist:
                sampleRPM=peak_attr_dict[peak]['chip_coverage'][sampleid]['RPM']
                fout.write('\t%6.2f'%sampleRPM)
            nearbygenedict=peak_attr_dict[peak]['peak_attr']['nearby_genes']
            nearbygenestrlist=[]
            for gene in nearbygenedict:
                genetsslist=nearbygenedict[gene]
                for genetss in genetsslist:
                    nearbygenestrlist.append('%s::%s:%010d'%(gene,genetss[0],genetss[1]))            
            fout.write('\t%s\t%s\t%s\t%s\t%s'%('|'.join(peak_attr_dict[peak]['peak_attr']['promotor']),
                                         '|'.join(peak_attr_dict[peak]['peak_attr']['enhancer']),
                                         '|'.join(peak_attr_dict[peak]['peak_attr']['kinase']),
                                         '|'.join(peak_attr_dict[peak]['peak_attr']['xcription_factor']),
                                         '|'.join(nearbygenestrlist)))  
            if len(extrafields)==0:
                fout.write('\n')
            else:
                if 'extra' in peak_attr_dict[peak]:
                    datafields=[peak_attr_dict[peak]['extra'][x] for x in extrafields]
                    fout.write('\t%s\n'%'\t'.join(datafields))
                else:
                    fout.write('\n')


def getcodetofilename(index_file_parameters,bamfile_id):
    """finds the bamfile name from the index file provided in the config dictionary
    
    The first two columns of the index file are:
        col1: bamfile_id: <tissue>_<enrichment>_<treatment>
        col2: location of the bam file
    """
    index_file=index_file_parameters['index']
    relative_flg=index_file_parameters['relative']
    
    index_dict=dict([(lntxt.rstrip().split(',')[0],lntxt.rstrip().split(',')[1]) for lntxt in open(index_file).readlines()])
    
    if bamfile_id not in index_dict:
        return ''
    
    if relative_flg==0:
        return index_dict[bamfile_id]
    else:
        relative_dir='/'.join(index_file.split('/')[0:-1])
        return '%s/%s'%(relative_dir,index_dict[bamfile_id])

def genPeakToolRunCommands(project_name,treatment_id, treatment_bamfile, control_bamfile, tool_parameters_dict, temp_dir ):
    """generates the run command for all the tools
    """
    cmd_dict={}
    if 'MACS' in tool_parameters_dict:
        macs_out_dir='%s/MACS'%(temp_dir)
        macs_exec=tool_parameters_dict['MACS']['exec']
        macs_options=tool_parameters_dict['MACS']['options']
        if not os.path.isdir(macs_out_dir):
            os.system('mkdir -p %s'%macs_out_dir)
        if treatment_bamfile==control_bamfile:
            macs_cmd='cd %s && %s -t %s -n %s %s'%(macs_out_dir,macs_exec,treatment_bamfile,treatment_id,macs_options)
        else:
            macs_cmd='cd %s && %s -t %s -c %s -n %s %s'%(macs_out_dir,macs_exec,treatment_bamfile,control_bamfile, treatment_id,macs_options)
        cmd_dict['MACS']=macs_cmd

    hmcan_out_dir='%s/HMCan'%(temp_dir)
    if not os.path.isdir(hmcan_out_dir):
        os.system('mkdir -p %s'%hmcan_out_dir)
    HMCan_config=tool_parameters_dict['HMCan']['config_file']
    HMCan_exec=tool_parameters_dict['HMCan']['exec']
    hmcan_cmd='cd %s && %s %s %s %s %s'%(hmcan_out_dir,HMCan_exec,treatment_bamfile,control_bamfile, HMCan_config, treatment_id)
    cmd_dict['HMCan']=hmcan_cmd 
    return cmd_dict

def getmappedreadcount(bamfile):
    return int(os.popen("samtools view -c -F 4  %s"%bamfile).read().strip('\n'))

def getchrsizes(bamfile):
    chrsizedict={}
    cmd='samtools idxstats %s'%bamfile
    stdo=os.popen(cmd)
    for linetxt in stdo:
        line=linetxt.rstrip('\n').split('\t')
        if line[0]=='*':
            continue
        chrsizedict[line[0]]=int(line[1])
    return chrsizedict              

def getmissedoutregions(peakfile,treatment_bamfile, min_size, min_coverage_gain_over_average,window_size):
    """
    The peakfile is assumed to be sorted
    The first three columns of the file are chr, peakstart, peakend
    Need chromosome sizes
    """
    total_read_count=getmappedreadcount(treatment_bamfile)
    chrsizedict=getchrsizes(treatment_bamfile)
    total_chr_size=sum([chrsizedict[chrnm] for chrnm in chrsizedict])
    avg_reads=1.0*total_read_count/total_chr_size
    indexed_bamfile=ChipSeqIndex.ChipSeqIndex(treatment_bamfile)
    tempmissedoutregionslist=[]
    prev_chrinterval=['chr1',0,0]
    lncnt=0
    for lntxt in open(peakfile):
        lncnt+=1
        ln=lntxt.rstrip('\n').split('\t')
        chrinterval=[ln[0],int(ln[1]),int(ln[2])]
        cadidategaplist=[]
        if prev_chrinterval[0]==chrinterval[0]:
            if (chrinterval[1]-prev_chrinterval[2]-1)>min_size:
                cadidategaplist.append([chrinterval[0],prev_chrinterval[2]+1,chrinterval[1]-1])
        else:
            if (chrsizedict[prev_chrinterval[0]]-prev_chrinterval[2])>min_size:
                cadidategaplist.append([prev_chrinterval[0],prev_chrinterval[2]+1,chrsizedict[prev_chrinterval[0]]])
            if (chrinterval[1]-1)>min_size:
                cadidategaplist.append([chrinterval[0],1,chrinterval[1]-1])
        for cadidatelargegap in cadidategaplist:
            for i in range((cadidatelargegap[2]-cadidatelargegap[1])/window_size+1):
                cadidategap=[cadidatelargegap[0],cadidatelargegap[1]+i*window_size,min(cadidatelargegap[1]+(i+1)*window_size,cadidatelargegap[2])]
                chrrange='%s:%d-%d'%tuple(cadidategap)
                readcounts=indexed_bamfile.countreadsinrange(chrrange)
                if readcounts*1.0/(chrinterval[1]-prev_chrinterval[2])> min_coverage_gain_over_average*avg_reads:
                    tempmissedoutregionslist.append(cadidategap)
        prev_chrinterval=chrinterval[0:]

    #join adjacent regions
    missedoutregionslist=[]
    if len(tempmissedoutregionslist)!=0:
        prev_region=tempmissedoutregionslist[0]
        for region in tempmissedoutregionslist[1:]:
            if region[0]==prev_region[0] and (region[1]-prev_region[2])<=1:
                prev_region[2]=region[2]
            else:
                missedoutregionslist.append(prev_region)
                prev_region=region[0:]
        missedoutregionslist.append(prev_region)
    
    return missedoutregionslist
    

def findpeaks(project_name, treatment_id, control_id, index_file_parameters, tool_parameters_dict, temp_dir, macs_cnv_region_identifiers, output_dir):
    """
    Merge the peaks from both MACS and HMCan and output the results as json file
    Input:
        treatment_bam_file
        control_bam_file
        tool_parameters_dict=options for the tools
        macs_cnv_region_identifiers=(min_region_size, min_coverage_gain_over_average,window_size)
        temp_directory=output files of MACS and HMCan
    method:
        Run MACS and HMCan
        Find the regions in MACS output that could be copy number regions
        In the macs gaps where size>min_size and coverage>min_coverage_ratio*average average
        use HMCan peaks
    Output:
        json file of peaks: dictionary 
        {
        chr1-0000713858-0000714238": {
        "RPM": 0.80, 
        "called_by":'MACS' }, 
    """
    treatment_bamfile=getcodetofilename(index_file_parameters,treatment_id)
    control_bamfile=getcodetofilename(index_file_parameters,control_id)
    
    cmd_dict=genPeakToolRunCommands(project_name,treatment_id,treatment_bamfile,control_bamfile, tool_parameters_dict, temp_dir )
    MACSpeakfile='%s/MACS/%s_peaks.bed'%(temp_dir,treatment_id)
    HMCanpeakfile='%s/HMCan/%s_regions.bed'%(temp_dir,treatment_id)
    
    if not os.path.exists(MACSpeakfile): 
        flog.write('%s: Running %s\n'%(time.asctime(),cmd_dict['MACS']))
        os.system(cmd_dict['MACS'])
    else:
        flog.write('%s: No need to run %s\nMACS peaks already there\n'%(time.asctime(),cmd_dict['MACS']))
    
    if not os.path.exists(HMCanpeakfile): 
        flog.write('%s: Running %s\n'%(time.asctime(),cmd_dict['HMCan']))    
        os.system(cmd_dict['HMCan'])
    else:
        flog.write('%s: No need to run %s\nHMCan peaks already there'%(time.asctime(),cmd_dict['HMCan']))    
     
    min_size,min_coverage_gain_over_average,window_size=macs_cnv_region_identifiers
    
    MACSpeaklist=[]
    for lntxt in open(MACSpeakfile):
        ln=lntxt.rstrip('\n').split('\t')
        MACSpeaklist.append([ln[0],int(ln[1]),int(ln[2])])   
    flog.write('%s: Info: number of MACS peaks %d\n'%(time.asctime(),len(MACSpeaklist)))
    missedoutregionslist=getmissedoutregions(MACSpeakfile,treatment_bamfile, min_size, min_coverage_gain_over_average,window_size)
       
    
    HMCanpeaklist=[]
    for lntxt in open(HMCanpeakfile):
        ln=lntxt.rstrip('\n').split('\t')
        HMCanpeaklist.append([ln[0],int(ln[1]),int(ln[2])])
    flog.write('%s: Info: number of HMCan peaks %d\n'%(time.asctime(),len(HMCanpeaklist)))
       
    HMCanadditions=common.interval_join(HMCanpeaklist, missedoutregionslist,3)
    flog.write('%s: Info: number of HMCan added peaks %d\n'%(time.asctime(),len(HMCanadditions)))
       
    all_peaklist=[]
    for peak in MACSpeaklist:
        all_peaklist.append(peak+['MACS'])
    for peak in HMCanadditions:
        all_peaklist.append(peak+['HMCan'])    
    all_peaklist.sort()
       
    outcsv='%s/peaks/%s__%s__peaks.bed'%(output_dir,project_name,treatment_id)
    outjson='%s/peaks/%s__%s__peaks.json'%(output_dir,project_name,treatment_id)
     
    fout=open(outcsv,'w')
    jsondict={}
       
    for peak in all_peaklist:
        fout.write('%s\t%d\t%d\t%s\n'%tuple(peak))
        jsondict['%s:%d-%d'%tuple(peak[0:3])]={}
        jsondict['%s:%d-%d'%tuple(peak[0:3])]['called_by']=peak[3]
       
    fout.close()
    json.dump(jsondict, open(outjson,'w'),indent=4,sort_keys=True)

def stitchpeaklist(inpeak_list,mergethreshold):
    """
    Stitch the peaks of the peaklist within given distance
    Used for stitching peaks for a sample or for union
    """
    peak_list=[]
    prev_peak=['chr0',0,1]
    inpeak_list.sort()
    for curr_peak in inpeak_list:
        if curr_peak[0]==prev_peak[0] and prev_peak[2]+mergethreshold>=curr_peak[1]:
            curr_peak[1]=min(prev_peak[1],curr_peak[1])
            curr_peak[2]=max(prev_peak[2],curr_peak[2])
        else:
            if prev_peak!=['chr0',0,1]:
                peak_list.append(prev_peak)
        prev_peak=curr_peak[:]
    peak_list.append(prev_peak)
    return peak_list   

def getpeakcoverages(treatment_id,control_id,peaklist,coveragepeakdict,readcountsdict):
    treatment_bamfile=getcodetofilename(index_file_parameters,treatment_id)
    control_bamfile=getcodetofilename(index_file_parameters,control_id)
    if treatment_id not in readcountsdict:
        readcountsdict[treatment_id]=getmappedreadcount(treatment_bamfile)
    if control_id not in readcountsdict:
        readcountsdict[control_id]=getmappedreadcount(control_bamfile)
    indexed_treatmentbamfile=ChipSeqIndex.ChipSeqIndex(treatment_bamfile)
    indexed_controlbamfile=ChipSeqIndex.ChipSeqIndex(control_bamfile)
    flog.write('%s: chip_indexes created : %s\n'%(time.asctime(),treatment_bamfile))
    for peak in peaklist:
        chrrange='%s:%d-%d'%tuple(peak)
        treatmentreadcounts=indexed_treatmentbamfile.countreadsinrange(chrrange)
        controlreadcounts=indexed_controlbamfile.countreadsinrange(chrrange)
        coveragepeakdict[chrrange]['chip_coverage'][treatment_id]={}
        treatment_RPM=treatmentreadcounts*1000000.0/readcountsdict[treatment_id]
        control_RPM=controlreadcounts*1000000.0/readcountsdict[control_id]
        coveragepeakdict[chrrange]['chip_coverage'][treatment_id]['treatment_RPM']=treatment_RPM
        coveragepeakdict[chrrange]['chip_coverage'][treatment_id]['Control_RPM']=control_RPM
        coveragepeakdict[chrrange]['chip_coverage'][treatment_id]['RPM']=treatment_RPM-control_RPM
    return coveragepeakdict

def getpeakattributes(peakdict,distancethreshold,promotorthreshold,genefeature,tflist,kinaselist):
    for peak in peakdict:
        peakenhancerlist=[]
        peakpromotorlist=[]
        peakkinaselist=[]
        peaktflist=[]
        peakinterval=[peak.split(':')[0],int(peak.split(':')[1].split('-')[0]),int(peak.split(':')[1].split('-')[1])]
        peakoverlapexon=genefeature.peak2exon(peakinterval)
        peakoverlapgenebody=genefeature.peak2genebody(peakinterval)
        nearbygenes=genefeature.peak2promotor(peakinterval, distancethreshold)
        nearby3primegenes=genefeature.peakto3prime(peakinterval,distancethreshold)
        for gene in nearbygenes:
            min_dist=min([x[1] for x in nearbygenes[gene]])
            min_abs_dist=min([abs(x[1]) for x in nearbygenes[gene]])
            if min_dist<-1*promotorthreshold and min_dist>-1*distancethreshold:
                peakenhancerlist.append(gene)
            if min_abs_dist<promotorthreshold:
                peakpromotorlist.append(gene)
            if gene in kinaselist:
                peakkinaselist.append(gene)
            if gene in tflist:
                peaktflist.append(gene)
        peakdict[peak]['peak_attr']['nearby_genes']=nearbygenes
        peakdict[peak]['peak_attr']['enhancer']=peakenhancerlist
        peakdict[peak]['peak_attr']['promotor']=peakpromotorlist
        peakdict[peak]['peak_attr']['kinase']=peakkinaselist
        peakdict[peak]['peak_attr']['xcription_factor']=peaktflist
        if len(peakpromotorlist)>0:
            peakdict[peak]['peak_attr']['classification']='promotor'
            peakdict[peak]['peak_attr']['classification_gene']=peakpromotorlist
        elif len(peakenhancerlist)>0:
            peakdict[peak]['peak_attr']['classification']='enhancer'
            peakdict[peak]['peak_attr']['classification_gene']=peakenhancerlist
        elif len(peakoverlapexon)>0: 
            peakdict[peak]['peak_attr']['classification']='genebody_exon'
            peakdict[peak]['peak_attr']['classification_gene']=[x[1] for x in peakoverlapexon]
        elif len(peakoverlapgenebody)>0: 
            peakdict[peak]['peak_attr']['classification']='genebody_intron'    
            peakdict[peak]['peak_attr']['classification_gene']=[x[1] for x in peakoverlapgenebody]   
        elif len(nearby3primegenes)>0:
            peakdict[peak]['peak_attr']['classification']='3 prime'    
            peakdict[peak]['peak_attr']['classification_gene']=[x[0] for x in nearby3primegenes]   
        else:
            peakdict[peak]['peak_attr']['classification']='other'
            peakdict[peak]['peak_attr']['classification_gene']=[]
        
    return peakdict
    
if __name__=='__main__':
    #config_filename='/media/ddrive/project/config/config_jon_proj01.json'
    import argparse

    parser = argparse.ArgumentParser(description='Process Chip-seq data')
    parser.add_argument('-c', '--config', dest='config_file',default='config.json',
                        help='config file name')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true',
                        default=0,help='write status on stdout')
    args = parser.parse_args()
    
    config=eval(open(args.config_file).read())
    
    project_name=config['project']['name']
    index_file_parameters=config['datasets']['index_file_parameters']
    tool_parameters_dict=config['peak_calling']['tools']
    mcri_dict=config['peak_calling']['macs_cnv_region_identifiers']
    macs_cnv_region_identifiers=(mcri_dict['min_size'],
                                 mcri_dict['min_coverage_gain_over_average'],
                                 mcri_dict['window_size'])
    comparison_list=config['datasets']['comparison']
    readfragmentsize=config['analysis']['parameters']['readfragmentsize']
    
    output_dir='%s/%s/output'%(config['output']['folder'],project_name)
    temp_dir='%s/%s/temp'%(config['output']['folder'],project_name)
    log_dir='%s/%s/log'%(config['output']['folder'],project_name)
    log_file='%s/%s.log'%(log_dir,project_name)
    if not os.path.isdir('%s/peaks'%output_dir):
        os.system('mkdir -p %s/peaks'%output_dir)
    if not os.path.isdir(temp_dir):
        os.system('mkdir -p %s'%temp_dir)
    if not os.path.isdir(log_dir):
        os.system('mkdir -p %s'%log_dir)    
    os.system('cp %s %s'%(args.config_file,log_dir))
    flog=open(log_file,'a')
    flog.write('\n\n\n\n######################\n%s: %s started\n'%(time.asctime(),project_name))

    redo_peak_calling_flag=config['analysis']['redo']['peak_calling']
    driverpeaksamplelist=config['analysis']['union']
    for comparison in comparison_list:
        treatment_id=comparison['treatment']
        control_id=comparison['control']
        peakfile='%s/peaks/%s__%s__peaks.bed'%(output_dir,project_name,treatment_id)
        if not (os.path.exists(peakfile) or redo_peak_calling_flag) and (treatment_id in driverpeaksamplelist):
            findpeaks(project_name, treatment_id, control_id, index_file_parameters, tool_parameters_dict, temp_dir, macs_cnv_region_identifiers, output_dir)
    flog.write('%s: Peak Finding done\n'%(time.asctime()))
         
    redo_project_analysis_flag=config['analysis']['redo']['projectjson']
    union_json='%s/peaks/%s__union__peaks.json'%(output_dir,project_name)
    union_bed='%s/peaks/%s__union__peaks.bed'%(output_dir,project_name)
    if not os.path.exists(union_json) or redo_project_analysis_flag:
        peakdistancethreshold=config['analysis']['parameters']['stitchpeakdistance']
        driverpeaksamplelist=config['analysis']['union']
         
        #Stitch Sample Peaks
        allstitchedpeaklist=[]
        for sample in driverpeaksamplelist:
            peakfile='%s/peaks/%s__%s__peaks.bed'%(output_dir,project_name,sample)
            inpeaklist=[[lntxt.split('\t')[0],int(lntxt.split('\t')[1]),int(lntxt.split('\t')[2])] for lntxt in open(peakfile)]
            stitched_peaklist=stitchpeaklist(inpeaklist,peakdistancethreshold)
            allstitchedpeaklist+=stitched_peaklist
            stitchedpeakfile='%s/peaks/%s__%s__stitchedpeaks.bed'%(output_dir,project_name,sample)
            fstitchout=open(stitchedpeakfile,'w')
            for peak in stitched_peaklist:
                fstitchout.write('%s\t%d\t%d\n'%tuple(peak))
            fstitchout.close()
         
        allstitchedpeaklist.sort()
        flog.write('%s: Peak Stitching done\n'%(time.asctime()))
         
        #Make Union Peaks
        union_peaklist=stitchpeaklist(allstitchedpeaklist,0)
        union_peakdict={}
        funionout=open(union_bed,'w')
        for peak in union_peaklist:
            union_peakdict['%s:%d-%d'%tuple(peak)]={'chip_coverage':{},'peak_attr':{}}
            funionout.write('%s\t%d\t%d\n'%tuple(peak))
        funionout.close()
        readcountsdict={}
        reportsamplelist=config['analysis']['report']
        for comparison in comparison_list:
            flog.write('%s: Computing RPM %s\n'%(time.asctime(), comparison))
            treatment_id=comparison['treatment']
            control_id=comparison['control']
            if treatment_id in reportsamplelist:
                union_peakdict=getpeakcoverages(treatment_id,control_id,union_peaklist,union_peakdict,readcountsdict)
        json.dump(union_peakdict, open(union_json,'w'),indent=4,sort_keys=True)
        flog.write('%s: json RPM done\n'%(time.asctime()))
     
    #get peak attributes
    union_attr_json='%s/peaks/%s__union__peaks_attr.json'%(output_dir,project_name)
    redo_project_peakattr_flag=config['analysis']['redo']['projectpeakattrjson']   
    if not os.path.exists(union_attr_json) or redo_project_peakattr_flag:
        flog.write('%s: Started peak attr json\n'%(time.asctime()))
        union_json='%s/peaks/%s__union__peaks.json'%(output_dir,project_name)
        union_peakdict=json.load(open(union_json))
        genetssfile=config['genefiles']['genetssgtf']
        genegtffile=config['genefiles']['genegtf']
        genefeature=GeneFeature.GeneFeature(genetssfile,genegtffile)
        distancethreshold=config['analysis']['parameters']['enhancerdistance']
        promotorthreshold=config['analysis']['parameters']['promotordistance']
        kinase_file=config['genefiles']['kinasefile']
        kinaselist=[lntxt.rstrip().split(',')[0] for lntxt in open(kinase_file).readlines()[1:]]
        transcriptionfactor_file=config['genefiles']['transcriptionfactorfile']
        tflist=[lntxt.rstrip().split(',')[0] for lntxt in open(transcriptionfactor_file).readlines()[1:]]
        union_peakdict=getpeakattributes(union_peakdict,distancethreshold,promotorthreshold,genefeature,tflist,kinaselist)
        json.dump(union_peakdict, open(union_attr_json,'w'),indent=4,sort_keys=True)
        
    peak_attr_dict=json.load(open('%s/peaks/%s__union__peaks_attr.json'%(output_dir,project_name,)))
    analysisreportdict=config['output']['report']
    reportfolder='%s/report'%output_dir
    reportsamplelist=config['analysis']['report']
    if not os.path.isdir(reportfolder):
        os.system('mkdir -p %s'%reportfolder)    
    for reportid in analysisreportdict:
        reportfilterdict=analysisreportdict[reportid]
        csvfile='%s/%s_report_%s.csv'%(reportfolder,project_name,reportid)
        redo_report_flag=config['analysis']['redo']['csv']
        if redo_report_flag==1 or not os.path.isfile(csvfile):
            flog.write('%s: Preparing %s\n'%(time.asctime(),csvfile))
            preparecsvanalysis(peak_attr_dict,reportsamplelist,reportfilterdict,csvfile)
    flog.write('%s: Done: Report processing\n'%(time.asctime()))
    
    redo_plot_avg_TSS=config['analysis']['redo']['plot_avg_TSS'] 
    redo_plot_superenhancer=config['analysis']['redo']['plot_superenhancer'] 
    redo_plot_chip_coverage=config['analysis']['redo']['plot_chip_coverage'] 
    plot_avg_halfwindowwidth=config['analysis']['parameters']['plot_avg_halfwindowwidth']
    plot_gene_halfwindowwidth=config['analysis']['parameters']['plot_gene_halfwindowwidth']
    avg_plot_samplelists=config['output']['plot']['avg_plot_samplelists']
    plotcoveragegenelist=config['output']['plot']['plotcoveragegenelist']
    plotsuperenhancerdict=config['output']['plot']['superenhancer']
    plotfolder='%s/plot'%output_dir
    if not os.path.isdir(plotfolder):
        os.system('mkdir -p %s'%plotfolder)

    #Plot
    try:
        genefeature
    except:
        genetssfile=config['genefiles']['genetssgtf']
        genegtffile=config['genefiles']['genegtf']
        genefeature=GeneFeature.GeneFeature(genetssfile,genegtffile)
    genefeature.getGeneTSS(all_flg=0)    
    genetssdict=genefeature.genetssdict
    
    #Plot Super Enhancers
    samplelist=plotsuperenhancerdict.keys()
    peak_attr_dict=json.load(open('%s/peaks/%s__union__peaks_attr.json'%(output_dir,project_name,)))
    superenhancercutoffdict=config['analysis']['parameters']['superenhancercutoff']
    for sampleid in samplelist:
        annotated_genelist=plotsuperenhancerdict[sampleid]
        if sampleid in superenhancercutoffdict:
            super_cutoff=superenhancercutoffdict[sampleid]
        elif 'default' in superenhancercutoffdict:
            super_cutoff=superenhancercutoffdict['default']
        else:
            super_cutoff=35
        imgfile='%s/%s_union_superenhancer.pdf'%(plotfolder,sampleid)
        if redo_plot_superenhancer==1 or not os.path.isfile(imgfile):
            flog.write('%s: Preparing %s\n'%(time.asctime(),imgfile))
            plot_superenhancer(sampleid,peak_attr_dict,imgfile,annotated_genelist,super_cutoff)
        imgfile='%s/%s_sample_superenhancer.pdf'%(plotfolder,sampleid)
        if redo_plot_superenhancer==1 or not os.path.isfile(imgfile):
            flog.write('%s: Preparing %s\n'%(time.asctime(),imgfile))
            sample_peak_attr_json='%s/peaks/%s__%s__sample__peaks_attr.json'%(output_dir,project_name,sampleid)
            if os.path.isfile(sample_peak_attr_json):
                sample_peak_attr_dict=json.load(open(sample_peak_attr_json))
            else:
                comparison_list=config['datasets']['comparison']
                control_id=[x['control'] for x in comparison_list if x['treatment']==sampleid][0]
                stitchedpeakfile='%s/peaks/%s__%s__stitchedpeaks.bed'%(output_dir,project_name,sampleid)
                stitchedpeaklist=[[lntxt.rstrip().split('\t')[0],int(lntxt.rstrip().split('\t')[1]),int(lntxt.rstrip().split('\t')[2])] for lntxt in open(stitchedpeakfile)]
                sample_peak_cov_dict={}
                for peak in stitchedpeaklist:
                    sample_peak_cov_dict['%s:%d-%d'%tuple(peak)]={'chip_coverage':{},'peak_attr':{}}
                try:
                    readcountsdict
                except:
                    readcountsdict={}
                sample_peak_cov_dict=getpeakcoverages(sampleid,control_id,stitchedpeaklist,sample_peak_cov_dict,readcountsdict)
                distancethreshold=config['analysis']['parameters']['enhancerdistance']
                promotorthreshold=config['analysis']['parameters']['promotordistance']
                kinase_file=config['genefiles']['kinasefile']
                kinaselist=[lntxt.rstrip().split(',')[0] for lntxt in open(kinase_file).readlines()[1:]]
                transcriptionfactor_file=config['genefiles']['transcriptionfactorfile']
                tflist=[lntxt.rstrip().split(',')[0] for lntxt in open(transcriptionfactor_file).readlines()[1:]]
                sample_peak_attr_dict=getpeakattributes(sample_peak_cov_dict,distancethreshold,promotorthreshold,genefeature,tflist,kinaselist)
                json.dump(sample_peak_attr_dict, open(sample_peak_attr_json,'w'),indent=4,sort_keys=True)
            plot_superenhancer(sampleid,sample_peak_attr_dict,imgfile,annotated_genelist,super_cutoff)            
            
    flog.write('%s: Done: Super Enhancer Plot\n'%(time.asctime()))
    
    #Plot Genes
    for gene in plotcoveragegenelist:
        genetss=genetssdict[gene]
        for samplelist in avg_plot_samplelists:
            bamfiledict={}
            for sampleid in samplelist:
                bamfiledict[sampleid]=getcodetofilename(index_file_parameters,sampleid)
            imgfile='%s/%s_%s_%02d.pdf'%(plotfolder,project_name,gene,avg_plot_samplelists.index(samplelist)+1)
            if redo_plot_chip_coverage==1 or not os.path.isfile(imgfile):
                flog.write('%s: Preparing %s\n'%(time.asctime(),imgfile))
                CHIPseq_geneplots(gene,genetss,bamfiledict,plot_gene_halfwindowwidth,imgfile,readfragmentsize)
    flog.write('%s: Done: Gene Coverage Plot\n'%(time.asctime()))
         
     
    #Plot Average TSS
    tsslist=genetssdict.values()
    for samplelist in avg_plot_samplelists:
        bamfiledict={}
        for sampleid in samplelist:
            bamfiledict[sampleid]=getcodetofilename(index_file_parameters,sampleid)
        imgfile='%s/%s_average_TSS_%02d.pdf'%(plotfolder,project_name,avg_plot_samplelists.index(samplelist)+1)
        if redo_plot_avg_TSS==1 or not os.path.isfile(imgfile):
            flog.write('%s: Preparing %s\n'%(time.asctime(),imgfile))
            CHIPseq_averageTSSplots(tsslist,bamfiledict,plot_avg_halfwindowwidth,imgfile,readfragmentsize)
    flog.write('%s: Done: Average TSS Plot\n'%(time.asctime()))
    flog.close()
    

    