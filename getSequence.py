#!/usr/bin/python
# -*- coding=utf8 -*-

import os
import re
import sys
#from baseFunctions import getSeq,raiseError
import math

###get start position in fasta
def getBtyes(start,length=60):
    startByte=start+int(start/length)
    return startByte

###get reference index info
def getIndex(fai):
    faiInfo={}
    lenInfo={}
    for line in open(fai,"r"):
        lineinfo=line.strip("\n").split("\t")
        faiInfo[lineinfo[0]]=int(lineinfo[2])
        lenInfo[lineinfo[0]]=int(lineinfo[1])
    return faiInfo,lenInfo	

def raiseError(string):
    sys.exit(string)

###get sequence
def getSeq(ref,chr,start,end,faiInfo,nearbySize=2):###1-based
    startBytes=getBtyes(start-nearbySize-1) ###byte positon for variant start
    endBytes=getBtyes(end+nearbySize-1) ###byte position for variant end
    if chr not in dict.keys(faiInfo):
        chr=re.sub("chr","",chr) if re.match("chr",chr) else "chr"+chr
    posBytes=faiInfo[chr] ###byte position for chr
    length=end-start+1  ###variant length
    ###get sequence
    fa=open(ref,"r")
    fa.seek(posBytes+startBytes,0) ###left 2bp
    sequence=fa.read(endBytes-startBytes+1) ###right 2bp
    sequence=re.sub("\n","",sequence)
    fa.close()
    return sequence[nearbySize:-(nearbySize)],sequence[:nearbySize],sequence[-(nearbySize):]

###iteration to get left dup region list,需要调整
def findleftDup(reference,chr,start,length,ref,alt,faiInfo,regionList,seqList,left2bpList,right2bpList,nearbySize):
    #print(chr,start,ref,alt)
    stemp=start-length
    etemp=start-1
    if stemp<1:
        return regionList,seqList,left2bpList,right2bpList
    varType,seq="",""
    if ref=="-":
        seq=alt
        varType="INS"
    elif alt=="-":
        seq=ref
        varType="DEL"
    else:
        seq=ref
        print("don't need to find dup")
        return regionList,seqList,left2bpList,right2bpList
    #print(stemp,etemp,seq)
    rt,st,et,seqt,left2bp,right2bp=finddupForleft(reference,chr,stemp,etemp,length,seq,faiInfo,nearbySize,varType)
    #print(rt,st,et)
    status=False
    if et!=start and ref=="-":
        regionList.append(str(st)+"-"+str(et))
        seqList.append(ref+":"+rt)
        left2bpList.append(left2bp)
        right2bpList.append(right2bp)
        status=True
    elif st!=start and alt=="-":
        regionList.append(str(st)+"-"+str(et))
        seqList.append(rt+":"+alt)
        left2bpList.append(left2bp)
        right2bpList.append(right2bp)
        status=True
    if not status:
        return regionList,seqList,left2bpList,right2bpList
    else:
        ref=seqList[-1].split(":")[0]
        alt=seqList[-1].split(":")[1]
        stemp=int(regionList[-1].split("-")[1]) if ref=="-" else int(regionList[-1].split("-")[0])
        #print(ref,alt,etemp)
        return findleftDup(reference,chr,stemp,length,ref,alt,faiInfo,regionList,seqList,left2bpList,right2bpList,nearbySize)

###findrightDup for INS
def finddupForright(reference,chr,s,e,length,seq,faiInfo,nearbySize,varType):
    seqtemp,left2bp,right2bp=getSeq(reference,chr,s,e,faiInfo,nearbySize)
    i=0
    while i < length:
        st=seq[i]
        stt=seqtemp[i]
        if st!=stt:
            break
        i+=1
    if varType=="INS":
        seqt,left2bp,right2bp=getSeq(reference,chr,s+i-1,s+i,faiInfo,nearbySize)
        return seq[i:]+seqtemp[:i],s+i-1,s+i,seqt,left2bp,right2bp
    elif varType=="DEL":
        seqt,left2bp,right2bp=getSeq(reference,chr,s+i-length,e+i-length,faiInfo,nearbySize) 
        return seq[i:]+seqtemp[:i],s+i-length,e+i-length,seqt,left2bp,right2bp
  
##findleftdup for INS
def finddupForleft(reference,chr,s,e,length,seq,faiInfo,nearbySize,varType):
    seqtemp,left2bp,right2bp=getSeq(reference,chr,s,e,faiInfo,nearbySize)
    #print(seqtemp,length,seq)
    i=length
    while i>0:
        st=seq[i-1]
        stt=seqtemp[i-1]
        if st!=stt:
            break
        i-=1
    if varType=="INS":
        seqt,left2bp,right2bp=getSeq(reference,chr,e-length+i,e-length+i+1,faiInfo,nearbySize)
        return seqtemp[i:]+seq[:i],e-length+i,e-length+i+1,seqt,left2bp,right2bp
    elif varType=="DEL":
        seqt,left2bp,right2bp=getSeq(reference,chr,s+i,e+i,faiInfo,nearbySize)
        return seqtemp[i:]+seq[:i],s+i,e+i,seqt,left2bp,right2bp

###iterration to get rigreference,chr,etemp,length,ref,alt,faiInfo,lenInfo,regionList,seqList,left2bpList,right2bpList,nearbySizeht dup region list
def findrightDup(reference,chr,end,length,ref,alt,faiInfo,lenInfo,regionList,seqList,left2bpList,right2bpList,nearbySize):
    stemp=end+1
    etemp=end+length
    if etemp > lenInfo[chr]:
        return regionList,seqList,left2bpList,right2bpList
    varType,seq="",""
    if ref=="-":
        seq=alt
        varType="INS"
    elif alt=="-":
        seq=ref
        varType="DEL"
    else:
        seq=ref
        print("don't need to find dup")
        return regionList,seqList,left2bpList,right2bpList
    #print(seq)
    rt,st,et,seqt,left2bp,right2bp=finddupForright(reference,chr,stemp,etemp,length,seq,faiInfo,nearbySize,varType)
    ##print(rt,st,et)
    status=False
    if et!=end+1 and ref=="-":
        regionList.append(str(st)+"-"+str(et))
        seqList.append(ref+":"+rt)
        left2bpList.append(left2bp)
        right2bpList.append(right2bp)
        status=True
    elif et!=end and alt=="-":
        regionList.append(str(st)+"-"+str(et))
        seqList.append(rt+":"+alt)
        left2bpList.append(left2bp)
        right2bpList.append(right2bp)
        status=True
    if not status:
        return regionList,seqList,left2bpList,right2bpList
    else:
        ref=seqList[-1].split(":")[0]
        alt=seqList[-1].split(":")[1]
        etemp=int(regionList[-1].split("-")[0]) if ref=="-" else int(regionList[-1].split("-")[1])
        #print(ref,alt,etemp)
        return findrightDup(reference,chr,etemp,length,ref,alt,faiInfo,lenInfo,regionList,seqList,left2bpList,right2bpList,nearbySize)

###get ref sequence for variant
def getVariantInfo(chr,start,end,reference,ref,alt,strand,direction):
    print(start,end,ref,alt)
    if not os.path.isfile(reference):
        raiseError("ReferenceError: no ref!")
    refIndex=reference+".fai"
    if not os.path.isfile(refIndex):
        raiseError("ReferenceError: no index file!")
    ###get nearby dup sequence
    faiInfo,lenInfo=getIndex(refIndex) ###bytes position for each chromsome
    ###get bytes position for variant in reference
    length=end-start+1
    dupList=[] 
    left2bpList=[]
    right2bpList=[]
    seqList=[]
    nearbySize=max(len(ref),len(alt))+4
    ###fix chr
    if chr not in lenInfo:
        chr=re.sub("chr","",chr) if "chr" in chr else "chr"+str(chr)
    ###get variant sequence
    seq,left2bp,right2bp=getSeq(reference,chr,start,end,faiInfo,nearbySize)
    #print(seq,left2bp,right2bp)
    if seq!=ref and ref!="-":
        raiseError("VariantError: "+chr+":"+str(start)+"-"+str(end)+":"+ref)
    ###对于SNP和等长度的ref及alt，不寻找dup
    if len(ref)==len(alt) and ref!="-" and alt!="-":
        return start,end,ref,alt,seq,left2bp,right2bp,nearbySize,faiInfo
    ###ins
    st,et=start,end
    if ref=="-":
        length=len(alt)
        st=end
        et=start
#        start=st
#        end=et
    ###优化点：根据参数进行查找，正链只需要找右边，负链向左找，位置
    if direction=="gd5":
        dupList,seqList,left2bpList,right2bpList=findleftDup(reference,chr,st,length,ref,alt,faiInfo,dupList,seqList,left2bpList,right2bpList,nearbySize)
    elif direction=="gd3":
        dupList,seqList,left2bpList,right2bpList=findrightDup(reference,chr,et,length,ref,alt,faiInfo,lenInfo,dupList,seqList,left2bpList,right2bpList,nearbySize)
    elif direction=="td5":
        if strand=="-":
            dupList,seqList,left2bpList,right2bpList=findrightDup(reference,chr,et,length,ref,alt,faiInfo,lenInfo,dupList,seqList,left2bpList,right2bpList,nearbySize)
        else: 
            dupList,seqList,left2bpList,right2bpList=findleftDup(reference,chr,st,length,ref,alt,faiInfo,dupList,seqList,left2bpList,right2bpList,nearbySize)
    elif direction=="td3":
        if strand=="+":
            dupList,seqList,left2bpList,right2bpList=findrightDup(reference,chr,et,length,ref,alt,faiInfo,lenInfo,dupList,seqList,left2bpList,right2bpList,nearbySize)
        else:
            dupList,seqList,left2bpList,right2bpList=findleftDup(reference,chr,st,length,ref,alt,faiInfo,dupList,seqList,left2bpList,right2bpList,nearbySize)
    if len(dupList)>=1:
        start=int(dupList[-1].split("-")[0])
        end=int(dupList[-1].split("-")[1])
        ref=seqList[-1].split(":")[0]
        alt=seqList[-1].split(":")[1]
        seq,left2bp,right2bp=getSeq(reference,chr,start,end,faiInfo,nearbySize)
 #   return seq,left2bp,right2bp,dupList,seqList,left2bpList,right2bpList,faiInfo
    #print(start,end,ref,alt,seq)
    return start,end,ref,alt,seq,left2bp,right2bp,nearbySize,faiInfo
