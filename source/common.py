import os
import math
import time
import sys
import random
import re
import csv
        
def func_name():
    """
    Returns the name of function its called from
    """
    return sys._getframe(1).f_code.co_name

def printstatus(message,message_type,func_nm,verbose_flg=0):
    """
    Prints the Status/Warning/Error message along with the function and time of the message 
    """
    if message_type=='E':
        print '%s : Error   : %30s : %s'%(time.asctime(),func_nm,message)
    if message_type=='W' and verbose_flg:
        print '%s : Warning : %30s : %s'%(time.asctime(),func_nm,message)
    if message_type=='S':
        print '%s : Status  : %30s : %s'%(time.asctime(),func_nm,message)     
    if message_type=='F':
        print '%s : Fatal Error  : %30s : %s'%(time.asctime(),func_nm,message) 
        sys.exit()
    
def fl2str(self,inlist,fmt='%8.4f'):
    """
    Prints a list of floats as comma separated string with given float format
    """
    return self._tfl2str(inlist,fmt)[:-1]

def _tfl2str(self,inlist,fmt='%8.4f'):
    """
    A recursive function used by function fl2str
    """
    outstr='['
    for el in inlist:
        if isinstance(el,list):
            outstr+=self._tfl2str(el,fmt)
        else:
            outstr+=fmt%el+','
    outstr=outstr[:-1]+'],'
    return outstr

def str2fl(self,istr,stype='float'):
    """
    Converts a list string (int or float) to a list
    """
    istr=istr.replace(' ','')
    strlist=istr[1:-1].split(',')
    if stype=='float':
        l1=[float(elem) for elem in strlist]
    if stype=='int':
        l1=[int(elem) for elem in strlist]        
    return l1

def median(self,l1):
    l1.sort()
    l=len(l1)
    if l%2==0:
        return 1.0/2*(l1[l/2]+l1[l/2-1])
    else:
        return l1[l/2]

def toss(self,pct=0.5):
    x=random.random()
    if x<pct:
        return 1
    else:
        return 0

def getlistidx(self,l1,element): 
    for i in range(len(l1)):
        if l1[i] > element:
            return i-1
    return len(l1)-1

def mergedict(self,a, b, path=None):
    "merges b into a"
    if path is None: path = []
    for key in b:
        if key in a:
            if isinstance(a[key], dict) and isinstance(b[key], dict):
                self.mergedict(a[key], b[key], path + [str(key)])
            elif a[key] == b[key]:
                pass # same leaf value
            else:
                raise Exception('Conflict at %s' % '.'.join(path + [str(key)]))
        else:
            a[key] = b[key]
    return a

def mergedictlist(self,dictlist):
    out={}
    for dicti in dictlist:
        self.mergedict(out,dicti)
    return out

def natural_sort(self,l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

def concatfiles(self,outfile,filelist,headerflg=1):
    fout=open(outfile,'w')
    filecnt=0
    for infile in filelist:
        filecnt+=1
        fin=open(infile)
        if filecnt!=1:
            fin.readline()
        for lntxt in fin:
            fout.write(lntxt)
        fin.close()

def csv2html(self,csvfile,htmlfile,delim='\t'):
    reader = csv.reader(csvfile,delimiter=delim)
    rownum = 0
    htmlfile.write('<table border="1">\n')
    for row in reader:
        if rownum == 0:
            htmlfile.write('<tr>') 
            for column in row:
                htmlfile.write('<th>' + column + '</th>')
            htmlfile.write('</tr>\n')
        else:
            htmlfile.write('<tr>')    
            for column in row:
                htmlfile.write('<td>' + column + '</td>')
            htmlfile.write('</tr>\n')
        rownum += 1
    htmlfile.write('</table>\n')
    print "Created " + str(rownum) + " row table." 

def lndistance(self,point1,point2,n=1):
    if isinstance(point1,list):
        if n==1:
            distance=1.0*sum([abs(point1[i]-point2[i]) for i in range(len(point1))])
        elif n==2:
            distance=math.sqrt(sum([(point1[i]-point2[i])**2 for i in range(len(point1))]))
    else:
        distance=abs(point1-point2)
    return max(distance,0.000001)
             
def interval_join(interval_list1,interval_list2,join_type,overlappct=0.7):
    """
    Input:
        list: list of (chr,pos1,pos2)
        join_type=0,1,2,3,4,5:intersection,list1 in list2, union, item1 70% covered in item2, item1 that covers 70% item2,item1 if either item1 or item2 is covered 70%
    Output:
        list
    Todo: Make efficient
    """
    out_interval_list=[]
    if join_type==0:
        for interval1 in interval_list1:
            for interval2 in interval_list2:
                if interval1[0]==interval2[0]:
                    if max(interval1[1],interval2[1]) < min(interval1[2],interval2[2]):
                        out_interval_list.append([interval1[0],max(interval1[1],interval2[1]),min(interval1[2],interval2[2])])        
    if join_type==1:
        for interval1 in interval_list1:
            for interval2 in interval_list2:
                if interval1[0]==interval2[0]:
                    if interval1[1]>=interval2[1] and interval1[2]<=interval2[2]:
                        out_interval_list.append(interval1)
    if join_type==3:
        for interval1 in interval_list1:
            for interval2 in interval_list2:
                if interval1[0]==interval2[0]:
                    if max(interval1[1],interval2[1]) <  min(interval1[2],interval2[2]):
                        if (min(interval1[2],interval2[2])-max(interval1[1],interval2[1]))>overlappct*(interval1[2]-interval1[1]):
                            out_interval_list.append(interval1)      
    if join_type==4:
        for interval1 in interval_list1:
            for interval2 in interval_list2:
                if interval1[0]==interval2[0]:
                    if max(interval1[1],interval2[1]) <  min(interval1[2],interval2[2]):
                        if (min(interval1[2],interval2[2])-max(interval1[1],interval2[1]))>overlappct*(interval2[2]-interval2[1]):
                            out_interval_list.append(interval1)                
    if join_type==5:
        for interval1 in interval_list1:
            for interval2 in interval_list2:
                if interval1[0]==interval2[0]:
                    if max(interval1[1],interval2[1]) <  min(interval1[2],interval2[2]):
                        overlapsize=min(interval1[2],interval2[2])-max(interval1[1],interval2[1])
                        if overlapsize>overlappct*min((interval2[2]-interval2[1]),(interval1[2]-interval1[1])):
                            out_interval_list.append(interval1)                                                                         
    return out_interval_list


def interval_join2(interval_list1,interval_list2,join_type,overlappct=0.7):
    """
    Input:
        list: list of (chr,pos1,pos2)
        join_type=
            0:intersection
            1:list1 in list2
            2:Union (not coded)
            3:item1 70% covered in item2
            4:item1 that covers 70% item2
            5:item1 if either item1 or item2 is covered 70%
    Output:
        list
    Todo: Make efficient
    """
    out_interval_list=[]
    if join_type==0:
        for interval1 in interval_list1:
            for interval2 in interval_list2:
                if interval1[0]==interval2[0]:
                    if max(interval1[1],interval2[1]) < min(interval1[2],interval2[2]):
                        out_interval_list.append(interval1+[max(interval1[1],interval2[1]),min(interval1[2],interval2[2])])        
    if join_type==1:
        for interval1 in interval_list1:
            for interval2 in interval_list2:
                if interval1[0]==interval2[0]:
                    if interval1[1]>=interval2[1] and interval1[2]<=interval2[2]:
                        out_interval_list.append(interval1+[interval2[1],interval2[2]])
    if join_type==3:
        for interval1 in interval_list1:
            for interval2 in interval_list2:
                if interval1[0]==interval2[0]:
                    if max(interval1[1],interval2[1]) <  min(interval1[2],interval2[2]):
                        if (min(interval1[2],interval2[2])-max(interval1[1],interval2[1]))>overlappct*(interval1[2]-interval1[1]):
                            out_interval_list.append(interval1+[interval2[1],interval2[2]])      
    if join_type==4:
        for interval1 in interval_list1:
            for interval2 in interval_list2:
                if interval1[0]==interval2[0]:
                    if max(interval1[1],interval2[1]) <  min(interval1[2],interval2[2]):
                        if (min(interval1[2],interval2[2])-max(interval1[1],interval2[1]))>overlappct*(interval2[2]-interval2[1]):
                            out_interval_list.append(interval1+[interval2[1],interval2[2]])                
    if join_type==5:
        for interval1 in interval_list1:
            for interval2 in interval_list2:
                if interval1[0]==interval2[0]:
                    if max(interval1[1],interval2[1]) <  min(interval1[2],interval2[2]):
                        overlapsize=min(interval1[2],interval2[2])-max(interval1[1],interval2[1])
                        if overlapsize>overlappct*min((interval2[2]-interval2[1]),(interval1[2]-interval1[1])):
                            out_interval_list.append(interval1+[interval2[1],interval2[2]])                                                                         
    return out_interval_list

