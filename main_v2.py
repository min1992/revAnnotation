 #!/usr/bin/python
# -*- coding=utf8 -*-
 ###function: reverse annotation
import re
import os
import sys
from argparse import ArgumentParser
import getSequence

STARTCODON="ATG"
ENDCONDON=["TAA","TAG","TGA"]
CONDON2AA={"TT[T|C]":"Phe","TT[A|G]":"Leu","CT.":"Leu","AT[T|C|A]":"Ile","ATG":"Met","GT.":"Val","TC.":"Ser","CC.":"Pro","AC.":"Thr","GC.":"Ala",
"TA[T|C]":"Tyr","CA[T|C]":"His","CA[A|G]":"Gln","AA[T|C]":"Asn","AA[A|G]":"Lys","GA[T|C]":"Asp","GA[A|G]":"Glu","TG[T|C]":"Cys",
"TGG":"Trp","CG.":"Arg","AG[T|C]":"Ser","AG[A|G]":"Arg","GG.":"Gly","T[AA|AG|GA]":"Ter"}
abbCONDON2AA={"Phe":"F","Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C","Gln":"Q","Glu":"E","Gly":"G","His":"H","Ile":"I",
"Leu":"L","Lys":"K","Met":"M","Pro":"P","Ser":"S","Thr":"T","Trp":"W","Tyr":"Y","Val":"V","Ter":"*"}
reverseCONDON2AA=dict(zip(CONDON2AA.values(),CONDON2AA.keys()))
reverseabbCONDON2AA=dict(zip(abbCONDON2AA.values(),abbCONDON2AA.keys()))

def obtainPars():
    parser=ArgumentParser()
    parser.add_argument("--refseq",dest="refseq",action="store",help="refseq file from UCSC",required=True)
    parser.add_argument("--variant",dest="variant",action="store",help="variant format transcript:cchange or transcript:pchange",required=True)
    parser.add_argument("--reference",dest="reference",action="store",help="reference",required=True)
    parser.add_argument("--direction",dest="direction",action="store",choices=['gd5','gd3','td5','td3'],help="move indel to different directionss",default="gd5")
    pars=parser.parse_args()
    return pars
###judge if two sections are intersectant
def judgeIntersect(start1,start2,end1,end2):
    maxStart=max(int(start1),int(start2))
    minEnd=min(int(end1),int(end2)) 
    if maxStart<=minEnd:
        return True,"middle"
    elif int(start2) > int(start1):
        return False,"right"
    elif int(start2) < int(start1):
        return False,"left"

##locate CDS region
def locateCDS(CDSstart,CDSend,exonstart,exonend):
    #print(exonstart)
    i,j=0,len(exonstart)-1
    while i < len(exonstart):
        status,ori=judgeIntersect(CDSstart,exonstart[i],CDSend,exonend[i])
        if status:
            break
        i+=1
    while j >=0:
        status,ori=judgeIntersect(CDSstart,exonstart[j],CDSend,exonend[j])
        if status:
            break
        j-=1
    leftUTRexonstart=exonstart[:i+1]##0-based
    leftUTRexonend=exonend[:i]+[CDSstart-1]###1-based
    rightUTRexonstart=[CDSend+1-1]+exonstart[j+1:]###0-based
    rightUTRexonend=exonend[j:]###1-based
    return exonstart[i:j+1],exonend[i:j+1],i,len(exonstart)-1-j,leftUTRexonstart,leftUTRexonend,rightUTRexonstart,rightUTRexonend

###judge the distance between the base in intron and leftCDNA or rightCDNA
def getDis(pos,leftborder,rightborder,length):
    leftdis=abs(pos-leftborder)
    rightdis=abs(pos-rightborder)
    leftCDNA=abs(leftborder-length)
    rightCDNA=leftCDNA+1
    if leftdis<=rightdis:
        return str(leftCDNA)+"+"+str(leftdis)
    else:
        return str(rightCDNA)+"-"+str(rightdis)

###judge distance for intron
def judgeD(pos,leftborder,rightborder):
    leftdis=abs(pos-leftborder)+1
    rightdis=abs(pos-rightborder)+1
    if leftdis<rightdis:
        return "left",leftdis
    elif leftdis>rightdis:
        return "right",rightdis
    else:
        return "both",leftdis
###get intron length
def getintronlength(exonstart,exonend,strand,label):
    length=[]
    if (strand=="+" and label=="5'UTR") or (strand=="-" and label=="3'UTR"):
        end,lengthtmp=0,0
        for i in list(range(len(exonstart)))[::-1]:
            if i==len(exonstart)-1:
                length.append(0)
                lengthtmp+=0
                end=int(exonstart[i])
            else:
                ###intron in UTR
                leftborder=int(exonend[i])+1
                rightborder=end
                lengthtmp+=(rightborder-leftborder+1)
                length.append(lengthtmp)
                end=int(exonstart[i])
        return length[::-1]
    elif (strand=="+" and label=="3'UTR") or (strand=="-" and label=="5'UTR"):
        start,lengthtmp=0,0
        for i in range(len(exonstart)):
            if i==0:
                lengthtmp+=0
                length.append(lengthtmp)
                start=int(exonend[i])+1
            else:
                leftborder=start
                rightborder=exonstart[i]
                lengthtmp+=(rightborder-leftborder+1)
                length.append(lengthtmp)
                start=int(exonend[i])+1
        return length

###parse UTR regions
def getUTRPosition(UTRexonstart,UTRexonend,boundary,strand,label,regionInfo,totalExons): 
    temp,left,right="","",""
    if label=="5'UTR":
        temp="-"
    elif label=="3'UTR":
        temp="*"
    else:
        os._exit("cann't recognise label:"+label)
    if strand=="+":
        left="+"
        right="-"
    else:
        left="-"
        right="+"
    if len(UTRexonstart)!=len(UTRexonend):
        os._exit("please check UTR start positions and UTR end positions!")
    intronlength=getintronlength(UTRexonstart,UTRexonend,strand,label) 
    start,end=float('inf'),float('inf')
    for i in range(len(UTRexonstart)):
        if i==0:
            if UTRexonstart[i]==UTRexonend[i]:
                continue
            for p in range(int(UTRexonstart[i])+1,int(UTRexonend[i])+1):
                regionInfo[temp+str(abs(p-boundary+intronlength[i]))]=str(p)+":"+label+"|"+str(totalExons)
                start=p+1
        else:
            if UTRexonstart[i]==UTRexonend[i]:
                end=float('inf')
            else:
                end=UTRexonstart[i]
            ###intron in UTR 
            for p in range(start,int(UTRexonstart[i])+1):
                direction,dis=judgeD(p,start,end)
                if direction=='left' or direction=="both":
                    regionInfo[temp+str(abs(start-1-boundary+intronlength[i-1]))+left+str(dis)]=str(p)+":"+label+"|"+str(totalExons)
                else:
                    regionInfo[temp+str(abs(end+1-boundary+intronlength[i]))+right+str(dis)]=str(p)+":"+label+"|"+str(totalExons)
            if UTRexonstart[i]==UTRexonend[i]:
                continue
            ###exon in UTR
            for p in range(int(UTRexonstart[i])+1,int(UTRexonend[i]+1)):
                regionInfo[temp+str(abs(p-boundary+intronlength[i]))]=str(p)+":"+label+"|"+str(totalExons)
                start=p+1
    return regionInfo

                

###get transcript and record CDNA position
def parseTranscript(lineinfo,reference,faiInfo):
    regionInfo={}
    chrom=lineinfo[2]
    strand=lineinfo[3]
    transtart=int(lineinfo[4])+1
    tranend=int(lineinfo[5])
    cdsstart=int(lineinfo[6])+1
    cdsend=int(lineinfo[7])
    exonstart=map(lambda x:int(x),re.split(",",lineinfo[9])[:-1])
    exonend=map(lambda x: int(x),re.split(",",lineinfo[10])[:-1]) 
    rawtotalExons=int(lineinfo[8])
    exonstart,exonend,leftexons,rightexons,leftUTRexonstart,leftUTRexonend,rightUTRexonstart,rightUTRexonend=locateCDS(cdsstart,cdsend,exonstart,exonend)
    #print(lineinfo[1],leftexons,rightexons,leftUTRexonstart,leftUTRexonend,rightUTRexonstart,rightUTRexonend)
    totalExons=len(exonstart)
    leftborder=0 ###left exon border for CDNA position
    rightborder=0 ###right exon border for CDNA position
    length=0 ###length for calculating CDNA change
    cdnaseq=""
    if strand=="+":
        if transtart!=cdsstart: 
            leftborder=transtart
            rightborder=cdsstart-1
            length=0
            if leftexons>=1:
                regionInfo=getUTRPosition(leftUTRexonstart,leftUTRexonend,cdsstart,strand,"5'UTR",regionInfo,totalExons)
            else:
                for p in range(leftborder,cdsstart):
                    regionInfo[str(p-cdsstart)]=str(p)+":5'UTR|"+str(totalExons)
        if tranend!=cdsend:
            leftborder=cdsend+1
            rightborder=tranend
            length=0
            if rightexons>=1:
                regionInfo=getUTRPosition(leftUTRexonstart,leftUTRexonend,cdsend,strand,"3'UTR",regionInfo,totalExons)
            else:
                for p in range(leftborder,tranend+1):
                    regionInfo["*"+str(p-cdsend)]=str(p)+":3'UTR|"+str(totalExons)
        start=0
        for i in range(totalExons):
            if i==0:
                leftborder=cdsstart ###left border
                rightborder=exonend[i]
                if totalExons==1:
                    rightborder=cdsend
                start=exonend[i]+1
                length+=cdsstart-1
                exonNO="exon"+str(i+1+leftexons)
                for p in range(leftborder,rightborder+1):
                    regionInfo[str(p-length)]=str(p)+":"+exonNO+"|"+str(totalExons)
                cdnaseq+=getSeq(reference,chrom,cdsstart,exonend[i],faiInfo)[0]
            else:
                leftborder=start
                rightborder=exonstart[i]
                exonNO="intron"+str(i+leftexons)
                for p in range(leftborder,rightborder+1):
                    cdna=getDis(p,leftborder-1,rightborder+1,length)
                    regionInfo[cdna]=str(p)+":"+exonNO+"|"+str(totalExons)
                length+=exonstart[i]-start+1
                if i<totalExons-1:
                    leftborder=exonstart[i]+1
                    rightborder=exonend[i]
                    exonNO='exon'+str(i+1+leftexons)
                    for p in range(leftborder,rightborder+1):
                        regionInfo[str(p-length)]=str(p)+":"+exonNO+"|"+str(totalExons)
                    cdnaseq+=getSeq(reference,chrom,exonstart[i]+1,exonend[i],faiInfo)[0]
                else:
                    leftborder=exonstart[i]+1
                    rightborder=cdsend
                    cdnaseq+=getSeq(reference,chrom,exonstart[i]+1,cdsend,faiInfo)[0]
                    exonNO='exon'+str(i+1+leftexons)
                    for p in range(leftborder,rightborder+1):
                        regionInfo[str(p-length)]=str(p)+":"+exonNO+"|"+str(totalExons)
                start=exonend[i]+1
    else:
        if transtart!=cdsstart:
            leftborder=transtart
            rightborder=cdsstart-1
            length=0
            exonNO="3'UTR"
            if leftexons>=1:
                regionInfo=getUTRPosition(leftUTRexonstart,leftUTRexonend,cdsstart,strand,"3'UTR",regionInfo,totalExons)
            else:
                for p in range(leftborder,rightborder+1):
                    regionInfo["*"+str(cdsstart-p)]=str(p)+":"+exonNO+"|"+str(totalExons)
        if tranend!=cdsend:
            leftborder=cdsend+1
            rightborder=tranend
            length=0
            exonNO="5'UTR"
            if rightexons>=1:
                regionInfo=getUTRPosition(rightUTRexonstart,rightUTRexonend,cdsend,strand,"5'UTR",regionInfo,totalExons)
            else:
                for p in range(leftborder,rightborder+1):
                    regionInfo[str(cdsend-p)]=str(p)+":"+exonNO+"|"+str(totalExons)
        end=0
        for i in range(totalExons-1,-1,-1): 
            if i==totalExons-1:
                leftborder=exonstart[i]+1
                rightborder=cdsend
                if totalExons==1:
                    leftborder=cdsstart
                end=exonstart[i]
                seqt=getSeq(reference,chrom,leftborder,rightborder,faiInfo)[0] 
                cdnaseq+=transRevAndComple(seqt)
                length=cdsend+1 ###leftCDNA=length-e 
                exonNO="exon"+str(totalExons-i+rightexons)
                for p in range(leftborder,rightborder+1):
                    regionInfo[str(length-p)]=str(p)+":"+exonNO+"|"+str(totalExons)
            else: 
                leftborder=exonend[i]+1
                rightborder=end
                exonNO="intron"+str(totalExons-i-1+rightexons)
                for p in range(leftborder,rightborder+1):
                    cdna=getDis(p,rightborder+1,leftborder-1,length)
                #    print(p,cdna)
                    regionInfo[cdna]=str(p)+":"+exonNO+"|"+str(totalExons)
                length-=(end-exonend[i])
                if i>0:
                    seqt=getSeq(reference,chrom,exonstart[i]+1,exonend[i],faiInfo)[0]
                    cdnaseq+=transRevAndComple(seqt)
                    leftborder=exonstart[i]+1
                    rightborder=exonend[i]
                    exonNO="exon"+str(totalExons-i+rightexons)
                    for p in range(leftborder,rightborder+1):
                        regionInfo[str(length-p)]=str(p)+":"+exonNO+"|"+str(totalExons)
                else:
                    seqt=getSeq(reference,chrom,cdsstart,exonend[i],faiInfo)[0]
                    cdnaseq+=transRevAndComple(seqt)
                    leftborder=cdsstart
                    rightborder=exonend[i]
                    exonNO="exon"+str(totalExons-i+rightexons)
                    for p in range(leftborder,rightborder+1):
                        regionInfo[str(length-p)]=str(p)+":"+exonNO+"|"+str(totalExons)
                end=exonstart[i]
    return regionInfo,chrom,strand,cdnaseq

###search transcript
def transcriptSearch(refseq,transcript):
    lines=[x for x in open(refseq).readlines() if transcript in x.strip("\n").split("\t")]
    lineinfo=lines[0].strip("\n").split("\t")
    return lineinfo

###trans reverse and complementary
def transRevAndComple(seq): 
    if seq=="none":
        return "none"
    basecomple={"A":"T","G":"C","C":"G","T":"A","-":"-",".":".","[":"]","]":"[","|":"|"}
    revcompleSeq=[basecomple[base] for base in seq[::-1]]
    return "".join(revcompleSeq)
###get start position in fasta
def getBtyes(start,length=60):
    startByte=start+int(start/length)
    return startByte

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

###get reference index info
def getIndex(fai):
    faiInfo={}
    lenInfo={}
    for line in open(fai,"r"):
        lineinfo=line.strip("\n").split("\t")
        faiInfo[lineinfo[0]]=int(lineinfo[2])
        lenInfo[lineinfo[0]]=int(lineinfo[1])
    return faiInfo,lenInfo

###parse transcript:c.
def parsecchange(transcript,change,refseq,reference,direction):
    #print(variant)
    fai=reference+".fai"
    if not os.path.isfile(fai):
        os._exit("no fai for reference!")
    faiInfo,lenInfo=getIndex(fai)
    ###get transcript line
    lineinfo=transcriptSearch(refseq,transcript)
    gene=lineinfo[-4]
    ##parse transcript info
    regionInfo,chrom,strand,cdnaseq=parseTranscript(lineinfo,reference,faiInfo)
    leftCDNA,rightCDNA,leftbase,rightbase,seq,varType="","","","","",""
    snppattern="c.(?P<leftCDNA>[\*\+\-]*[0-9]*[\*\+\-]*[0-9]*)(?P<leftbase>[ATGCNagctn]*)(?P<varType>\>)(?P<rightbase>[ATGCNatgcn]*)"
    indelpattern="c.(?P<leftCDNA>[\*\+\-]*[0-9]*[\*\+\-]*[0-9]*)(|_(?P<rightCDNA>[\*\+\-]*[0-9]*[\*\+\-]*[0-9]*))(?P<varType>(delins|del|ins|dup|inv))(?P<seq>[ATGCNatgcn]*)"
    #print(change)
    if ">" in change:###snv 
        try:
            res=re.search(snppattern,change)
            groups=res.groupdict()
        except:
            os._exit("variant format is wrong: "+transcript+":"+change)
        leftCDNA=groups['leftCDNA']
        leftbase=groups['leftbase']
        rightbase=groups['rightbase']
        if strand=="-":
            leftbase=transRevAndComple(leftbase)
            rightbase=transRevAndComple(rightbase)
        try:
            position=regionInfo[leftCDNA].split(":")[0]
        except:
            return None,gene
        s,e,ref,alt,seq,left5bp,right5bp,nearbySize,faiInfo=getSequence.getVariantInfo(chrom,int(position),int(position),reference,leftbase,rightbase,strand,direction)
        return chrom+":"+str(s)+"-"+str(e)+":"+ref+":"+alt,gene
    else:
        #print(change)
        try:
            res=re.search(indelpattern,change)
            groups=res.groupdict()
        except:
            os._exit("variant format id wrong: "+variant)
        leftCDNA=groups['leftCDNA']
        try:
            leftposition=regionInfo[leftCDNA].split(":")[0]
        except:
            return None,gene
        #print(groups)
        if ("seq" in groups.keys()) and (groups["varType"]!="del"):
            insseq=groups['seq'] if strand=="+" else transRevAndComple(groups['seq'])
        else:
            insseq="-"
        if groups['varType']=="del" or groups['varType']=="delins":##long del and one-based del
            if 'rightCDNA' in groups.keys() and groups['rightCDNA']!="" and groups['rightCDNA']:#long del
                rightCDNA=groups['rightCDNA']
                try:
                    rightposition=regionInfo[rightCDNA].split(":")[0]
                except:
                    return None,gene
                seq=getSeq(reference,chrom,int(leftposition),int(rightposition),faiInfo)[0] if int(leftposition)<=int(rightposition) else getSeq(reference,chrom,int(rightposition),int(leftposition),faiInfo)[0]
                #print(seq)
                if int(leftposition)<=int(rightposition):
                    s,e,ref,alt,seq,left5bp,right5bp,nearbySize,faiInfo=getSequence.getVariantInfo(chrom,int(leftposition),int(rightposition),reference,seq,insseq,strand,direction)
                else:
                    s,e,ref,alt,seq,left5bp,right5bp,nearbySize,faiInfo=getSequence.getVariantInfo(chrom,int(rightposition),int(leftposition),reference,seq,insseq,strand,direction)
                #print(s,e,ref,alt)
#                posinfo=chrom+":"+leftposition+"-"+rightposition+":"+seq+":"+insseq if int(leftposition)<=int(rightposition) else chrom+":"+rightposition+"-"+leftposition+":"+seq+":"+insseq
                posinfo=chrom+":"+str(s)+"-"+str(e)+":"+ref+":"+alt
                return posinfo,gene
            else:#one-based del
                seq=getSeq(reference,chrom,int(leftposition),int(leftposition),faiInfo)[0]
                s,e,ref,alt,seq,left5bp,right5bp,nearbySize,faiInfo=getSequence.getVariantInfo(chrom,int(leftposition),int(leftposition),reference,seq,insseq,strand,direction)
                return chrom+":"+str(s)+"-"+str(e)+":"+ref+":"+alt,gene
        elif groups['varType']=="ins":###
            rightCDNA=groups['rightCDNA']
            try:
                rightposition=regionInfo[rightCDNA].split(":")[0]
            except:
                return None,gene
            if int(leftposition)<=int(rightposition):
                s,e,ref,alt,seq,left5bp,right5bp,nearbySize,faiInfo=getSequence.getVariantInfo(chrom,int(leftposition),int(rightposition),reference,"-",insseq,strand,direction)
            else:
                s,e,ref,alt,seq,left5bp,right5bp,nearbySize,faiInfo=getSequence.getVariantInfo(chrom,int(rightposition),int(leftposition),reference,"-",insseq,strand,direction)
#            posinfo=chrom+":"+leftposition+"-"+rightposition+":-:"+insseq if int(leftposition)<=int(rightposition) else chrom+":"+rightposition+"-"+leftposition+":-:"+insseq
            posinfo=chrom+":"+str(s)+"-"+str(e)+":"+ref+":"+alt
            return posinfo,gene
        elif groups['varType']=="dup":###dup 
            if 'rightCDNA' in groups.keys() and groups['rightCDNA']!="" and groups['rightCDNA']:
                rightCDNA=groups['rightCDNA']
                try:
                    leftposition=regionInfo[rightCDNA].split(":")[0] if strand=="+" else int(regionInfo[rightCDNA].split(":")[0])-1
                except:
                    return None,gene
                rightposition=str(int(leftposition)+1)
                if insseq=="-" or insseq=="":
                    seq=getSeq(reference,chrom,int(regionInfo[leftCDNA].split(":")[0]),int(regionInfo[rightCDNA].split(":")[0]),faiInfo)[0] if strand=="+" else getSeq(reference,chrom,int(regionInfo[rightCDNA].split(":")[0]),int(regionInfo[leftCDNA].split(":")[0]),faiInfo)[0]
                    insseq=seq
#                    print(insseq)
                s,e,ref,alt,seq,left5bp,right5bp,nearbySize,faiInfo=getSequence.getVariantInfo(chrom,int(leftposition),int(rightposition),reference,"-",insseq,strand,direction)
                return chrom+":"+str(s)+"-"+str(e)+":"+ref+":"+alt,gene
            else:
                leftposition=leftposition if strand=="+" else int(leftposition)-1
                rightposition=str(int(leftposition)+1)
                if insseq=="-" or insseq=="":
                    seq=getSeq(reference,chrom,int(leftposition),int(leftposition),faiInfo)[0] if strand=="+" else getSeq(reference,chrom,int(rightposition),int(rightposition),faiInfo)[0]
                    insseq=seq
                s,e,ref,alt,seq,left5bp,right5bp,nearbySize,faiInfo=getSequence.getVariantInfo(chrom,int(leftposition),int(rightposition),reference,"-",insseq,strand,direction)
                return chrom+":"+str(s)+"-"+str(e)+":"+ref+":"+alt,gene
        elif groups['varType']=="inv":###inversion
            rightCDNA=groups['rightCDNA']
            try:
                rightposition=regionInfo[rightCDNA].split(":")[0] if strand=="+" else leftposition
            except:
                return None,gene
            leftposition=leftposition if starnd=="+" else regionInfo[rightCDNA].split(":")[0]
            seq=getSeq(reference,chrom,int(leftposition),int(rightposition),faiInfo)[0]
            insseq=transRevAndComple(insseq)
            if insseq=="-" or insseq=="":
                insseq=transRevAndComple(seq)
            s,e,ref,alt,seq,left5bp,right5bp,nearbySize,faiInfo=getSequence.getVariantInfo(chrom,int(leftposition),int(rightposition),reference,seq,insseq,strand,direction)
            return chrom+":"+str(s)+"-"+str(e)+":"+ref+":"+alt,gene
        else:
            os._exit("unknown variant type: "+variant)

##trans abberivation into fullname
def transAbbintoFullname(aas):
    newaas=""
    firstaa=aas[0:3]
    if firstaa not in reverseabbCONDON2AA.values():
        while i < len(aas):
            newaas+=reverseabbCONDON2AA[aas[i]]
    else:
        return aas
    return newaas

##trans aas into bases
def transAAintobases(aas):
    bases=""
    i=0
    while i<len(aas):
        bases+=reverseCONDON2AA[aas[i:i+3]]
        i+=3
    return bases

###compare refbases and altbases, exclude same bases in header or tail
def comparebases(refbases,altbases):
    i,leftexcnum,rightexcnum=0,0,0
    length=min(len(refbases),len(altbases))
    if altbases=="":
        return refbases,altbases,leftexcnum,rightexcnum
    while i <length:
        reft=refbases[i:i+1]
        altt=altbases[i:i+1]
        if reft!=altt:
            break
        i+=1
    refbases=refbases[i:]
    altbases=altbases[i:]
    leftexcnum=i
    length=min(len(refbases),len(altbases))
    i=0
    while i < length:
        reft=refbases[-i-1:-i] if i>0 else refbases[-i-1:]
        altt=altbases[-i-1:-i] if i>0 else altbases[-i-1:]
        if reft!=altt:
            break
        i+=1
    rightexcnum=i
    refbases=refbases[0:-i] if i>0 else refbases[0:]
    altbases=altbases[0:-i] if i>0 else altbases[0:] 
    return refbases,altbases,leftexcnum,rightexcnum
###prolong variant
#def extendVariant(leftposition,rightposition,strand,chrom,ref,alt,reference):
#    fai=reference+".fai"
#    if not os.path.isfile(fai):
#        os._exit("no fai for reference!")
#    faiInfo,lenInfo=getIndex(fai)
#    for 
#    seqt=getSeq(reference,chrom,leftposition,leftposition,faiInfo,0)

###find nearest endcondon
def findnearerstENDCONDON(cdnaseq):
    i=0
    while i<len(cdnaseq)-3:
        temp=cdnaseq[i:i+3]
        if temp in ENDCONDON:
            return i
        i+=1
    return False

###parse trancript:p.
def parsepchange(transcript,change,refseq,reference):
    fai=reference+".fai"
    if not os.path.isfile(fai):
        os._exit("no fai for reference!")
    faiInfo,lenInfo=getIndex(fai)
    ###get transcript line
    lineinfo=transcriptSearch(refseq,transcript)
    gene=lineinfo[-4]
    ##parse transcript info
    regionInfo,chrom,strand,cdnaseq=parseTranscript(lineinfo,reference,faiInfo)
    ##single aa
    singleAApattern="p.(?P<leftaa>[a-zA-Z]*)(?P<leftaaNO>[0-9]*)(?P<rightaa>[a-zA-Z]*)"
    ##frameshift: del or ins
    singleAAfspattern="p.(?P<leftaa>[a-zA-Z]*)(?P<leftaaNO>[0-9]*)(?P<rightaa>[a-zA-Z]*)(fs\*(?P<fsNO>[0-9]*)|fsTer(?P<fsNO>[0-9]*))"
    ##nonframeshift: del or ins aas
    delinsAApattern="p.(?P<leftaa>[a-zA-Z]*)(?P<leftaaNO>[0-9]*)(|_(?P<rightaa>[a-zA-Z]*)(?P<rightaaNO>[0-9]*))(?P<varType>(delins|del|ins|dup))(?P<insaa>[a-zA-Z]*)"
    if "del" in change or "ins" in change or "delins" in change:
        res=re.search(delinsAApattern,change)
        groups=res.groupdict()
        varType=groups['varType']
        leftaa=transAbbintoFullname(groups['leftaa'])
        leftaaNO=groups['leftaaNO']
        if ("insaa" in groups.keys()) and (varType!="del"):
            insaa=transAbbintoFullname(groups['insaa'])
            insseq=transAAintobases(insaa) if strand=="+" else transRevAndComple(transAAintobases(insaa))
        else:
            insseq="-"
        if 'rightaa' in groups.keys():###del >=2 aas or ins aa, need to exclude ref from condon[] of aa
            rightaa=transAbbintoFullname(groups['rightaa'])
            rightaaNO=groups['rightaaNO']
            if varType=="ins":
                leftCDNA=int(leftaaNO)*3
                rightaaNO=int(rightaaNO)*3-2
                if strand=="+":
                    try:
                        leftposition=regionInfo[str(leftCDNA)].split(":")[0]
                        rightposition=regionInfo[str(rightCDNA)].split(":")[0]
                    except:
                        return None,gene
                else:
                    try:
                        leftposition=regionInfo[str(rightCDNA)].split(":")[0]
                        rightposition=regionInfo[str(leftCDNA)].split(":")[0]
                    except:
                        return None,gene
#                ref=""
#                alt=insseq
                return chrom+":"+str(leftposition)+"-"+str(rightposition)+":-:"+insseq,gene
            elif varType=="dup":
                rightCDNA=int(rightaaNO)*3
                if strand=="+":
                    try:
                        leftposition=int(regionInfo[str(rightCDNA)].split(":")[0])
                        rightposition=int(leftposition)+1
                    except:
                        return None,gene
                else:
                    try:
                        rightposition=int(regionInfo[str(rightCDNA)].split(":")[0])
                        leftposition=rightposition-1
                    except:
                        return None,gene
                return chrom+":"+str(leftposition)+"-"+str(rightposition)+":-:"+insseq,gene
            else:
                leftCDNA=int(leftaaNO)*3-2
                rightCDNA=int(rightaaNO)*3
                if strand=="+":
                    try:
                        leftposition=int(regionInfo[str(leftCDNA)].split(":")[0])
                        rightposition=int(regionInfo[str(rightCDNA)].split(":")[0])
                    except:
                        return None,gene
                else:
                    try:
                        leftposition=int(regionInfo[str(rightCDNA)].split(":")[0])
                        rightposition=int(regionInfo[str(leftCDNA)].split(":")[0])
                    except:
                        return None,gene
                seq=getSeq(reference,chrom,leftposition,rightposition,faiInfo)[0]
                seqbases,insbases,leftexcnum,rightexcnum=comparebases(seq,insseq)
                leftposition+=leftexcnum
                rightposition-=rightexcnum
                return chrom+":"+str(leftposition)+"-"+str(rightposition)+":"+seqbases+":"+insbases,gene
        else:######del one aa
            if varType=="dup":
                rightCDNA=int(leftaaNO)*3
                if strand=="+":
                    try:
                        leftposition=int(regionInfo[str(rightCDNA)].split(":")[0])
                        rightposition=leftposition+1
                    except:
                        return None,gene
                else:
                    try:
                        rightposition=int(regionInfo[str(rightCDNA)].split(":")[0])
                        leftposition=rightposition-1
                    except:
                        return None,gene
                return chrom+":"+str(leftposition)+"-"+str(rightposition)+":-:"+insseq,gene
            else:
                leftCDNA=int(leftaaNO)*3-2
                rightCDNA=int(leftaaNO)*3
                if strand=="+":
                    try:
                        leftposition=int(regionInfo[str(leftCDNA)].split(":")[0])
                        rightposition=int(regionInfo[str(rightCDNA)].split(":")[0])
                    except:
                        return None,gene
                else:
                    try:
                        leftposition=int(regionInfo[str(rightCDNA)].split(":")[0])
                        rightposition=int(regionInfo[str(leftCDNA)].split(":")[0])
                    except:
                        return None,gene
                seq=getSeq(reference,chrom,leftposition,rightposition,faiInfo)[0]
                seqbases,insbases,leftexcnum,rightexcnum=comparebases(seq,insseq)
                leftposition+=leftexcnum
                rightposition-=rightexcnum
                return chrom+":"+str(leftposition)+"-"+str(rightposition)+":"+seqbases+":"+insbases,gene
    elif "fs" in change:##framshift variant, which is complicated and put away
        res=re.search(singleAAfspattern,change)
        groups=res.groupdict()
        leftaa=transAbbintoFullname(groups['leftaa'])
        leftaaNO=groups['leftaaNO']
        rightaa=transAbbintoFullname(groups['rightaa'])
        fsNO=int(groups['fsNO'])*3
        leftCDNA=int(leftaaNO)*3-2
        rightCDNA=int(leftaaNO)*3
        refcdnaseq=cdnaseq[leftCDNA:]
        if strand=="+":
            try:
                leftposition=int(regionInfo[str(leftCDNA)].split(":")[0])
                rightposition=int(regionInfo[str(rightCDNA)].split(":")[0])
            except:
                return None,gene
        else:
            try:
                leftposition=int(regionInfo[str(rightCDNA)].split(":")[0])
                rightposition=int(regionInfo[str(leftCDNA)].split(":")[0])
            except:
                return None,gene
        seq=getSeq(reference,chrom,leftposition,rightposition,faiInfo)[0]
        altseq=reverseCONDON2AA[rightaa]
        endCondonIndex=findnearerstENDCONDON(refcdnaseq)
       # if not endCondonIndex:###can't find endcondon, insertion endcondon
            
        #elif fsNO < endCondonIndex:###ins leading to fs, greed strategy

        refseq,altseq,leftexcnum,rightexcnum=comparebases(seq,altseq)
        leftposition+=leftexcnum
        rightposition-=rightexcnum
        return chrom+":"+str(leftposition)+"-"+str(rightposition)+":"+refseq+":"+altseq,gene
    else: ###single aa 
        res=re.search(singleAApattern,change)
        groups=res.groupdict() 
        leftaa=transAbbintoFullname(groups['leftaa'])
        leftaaNO=groups['leftaaNO']
        rightaa=transAbbintoFullname(groups['rightaa'])
        leftCDNA=int(leftaaNO)*3-2
        rightCDNA=int(leftaaNO)*3
        if strand=="+":
            try:
                leftposition=int(regionInfo[str(leftCDNA)].split(":")[0])
                rightposition=int(regionInfo[str(rightCDNA)].split(":")[0])
            except:
                return None,gene
        else:
            try:
                leftposition=int(regionInfo[str(rightCDNA)].split(":")[0])
                rightposition=int(regionInfo[str(leftCDNA)].split(":")[0])
            except:
                return None,gene
        seq=getSeq(reference,chrom,leftposition,rightposition,faiInfo)[0]
        altseq=reverseCONDON2AA[rightaa]
        refseq,altseq,leftexcnum,rightexcnum=comparebases(seq,altseq)
        leftposition+=leftexcnum
        rightposition-=rightexcnum
        return chrom+":"+str(leftposition)+"-"+str(rightposition)+":"+refseq+":"+altseq,gene

###get same variant in 3' direction
#s,e,ref,alt,seq,left5bp,right5bp,nearbySize,faiInfo=getSequence.getVariantInfo(c,s,e,fasta,ref,alt,strand)


###find other transcripts when can't search correct position in provided transcript
def searchOT(gene,refseq):
    transinfo=[x.strip("\n").split("\t")[1] for x in open(refseq).readlines() if gene in x.strip("\n").split("\t")]
    return transinfo

def main():
    pars=obtainPars()
    refseq,variant,reference,direction = pars.refseq,pars.variant,pars.reference,pars.direction
    if "c." in variant:
        transcript,cchange=variant.split(":")[0],variant.split(":")[1]
#        print(parsecchange(transcript,cchange,refseq,reference))
        info,gene=parsecchange(transcript,cchange,refseq,reference,direction)
        transinfo=searchOT(gene,refseq)
        if not info:
            temp=""
            for trans in transinfo:
                if re.match("NR",trans):
                    continue
                info,gene=parsecchange(trans,cchange,refseq,reference,direction)
                if info:
                    temp+=trans+":"+info+";"
            print(temp)
        else:
            print(transcript+":"+info)
    elif "p." in variant:
        transcript,pchange=variant.split(":")[0],variant.split(":")[1]
        info,gene=parsepchange(transcript,pchange,refseq,reference,direction)
        transinfo=searchOT(gene,refseq)
        if not info:
            temp=""
            for trans in transinfo:
                if re.match("NR",trans):
                    continue
                info,gene=parsecchange(trans,pchange,refseq,reference,direction)
                if info:
                    temp+=trans+":"+info+";"
            print(temp)
        else:
            print(transcript+":"+info)
    else:
        os._exit("the variant format is wrong: "+variant)

if __name__ == "__main__":
    main()
