#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import pandas as pd
from Bio import SeqIO
import re
import numpy as np
import time

anno_file = sys.argv[1]
exonskip_file = sys.argv[2]
ref_file = sys.argv[3]
geneRef_file = sys.argv[4]
intron_ref_file = sys.argv[5]
output_file = sys.argv[6]

'''Part 0: import files'''
anno_df = pd.read_csv(anno_file, index_col=0, sep=' ')

#join on the info for exon skipping to find the intron number
exonskip_df = pd.read_csv(exonskip_file, index_col=0, sep=' ')

anno_df = anno_df.join(exonskip_df[['tags_skipping', 'curTag_skipping', 'relStartExon', 'relEndExon']], how='left', on='intron_junction')

records = list(SeqIO.parse(ref_file, 'fasta'))

geneRef_df = pd.read_csv(geneRef_file, index_col=False)
geneRef_df = geneRef_df[geneRef_df['type'] == 'transcript']

intron_df = pd.read_csv(intron_ref_file, sep='\t', header=None)


'''Part 1: look for stop codons in the reference file'''
def nucleotideSeqFinder(start, end, strand, sequence, gene, intronDataFrame, tag, startExon, endExon, exonDataFrame):
    '''builds the nucleotide sequence for a gene with a given AS junction'''
    try:
        if (np.isnan(startExon) == False) & (np.isnan(endExon) == False):
            intronStarts = list(intronDataFrame[1][(intronDataFrame[3] == gene) & (intronDataFrame[9] == tag) & (intronDataFrame[7] != startExon)])
            intronEnds = list(intronDataFrame[2][(intronDataFrame[3] == gene) & (intronDataFrame[9] == tag) & (intronDataFrame[7] != endExon)])

            intronEnds = [x-1 for x in intronEnds]

            intronStarts.insert(int(startExon), int(start))
            intronEnds.insert(int(endExon), int(end))
            transcripts = exonDataFrame[(exonDataFrame['gene'] == gene)]
            transcripts['len'] = exonDataFrame.apply(lambda x: abs(x['end'] - x['start']), axis=1)
            transcripts = transcripts.sort_values('len', axis=0, ascending=False)
            geneStart = [transcripts.iloc[0,list(transcripts.columns.values).index('start')]]
            geneEnd = [transcripts.iloc[0,list(transcripts.columns.values).index('end')]]
    

            exonStarts = geneStart
            exonStarts.extend(intronEnds)
            exonEnds = intronStarts
            exonStarts.sort()
            exonEnds.sort()
            exonEnds.extend(geneEnd)
            exons = [(exonStarts[i],exonEnds[i]) for i in range(len(exonStarts))]

            outputSeq = ''
            for exon in exons:
                outputSeq = outputSeq + sequence[exon[0]:exon[1]]
                outputSeq = str(outputSeq)
            if strand == '-':
                outputSeq = outputSeq[::-1]
            if strand == '+': 
                startToJuncDist = sum([sum(abs(exon[0] - exon[1]) for exon in exons if exon[0] < start)])
            if strand == '-':
                startToJuncDist = sum([sum(abs(exon[0] - exon[1]) for exon in exons if exon[1] > end)])
            if strand == '+':
                seqStarts = [m.start(0) for m in re.finditer('ATG', outputSeq[:startToJuncDist])]
            
            elif strand == '-':
                seqStarts = [m.start(0) for m in re.finditer('TAC', outputSeq[:startToJuncDist])]
            if len(seqStarts) == 0:
                if strand == '+':
                    seqStarts = [m.start(0) for m in re.finditer('ATG', outputSeq)]
                elif strand == '-':
                    seqStarts = [m.start(0) for m in re.finditer('TAC', outputSeq)]
            if strand == '-':
                last_exon_len = abs(exons[0][1] - exons[0][0])
            elif strand == '+':
                last_exon_len = abs(exons[-1][1] - exons[-1][0])
    
            if strand == '-':
                outputSeq = outputSeq[::-1]
            
            return(pd.Series([str(outputSeq), last_exon_len, seqStarts, startToJuncDist]))
            
        else:
            return(pd.Series(['', '', '', '']))
    except:
        return(pd.Series(['', '', '', '']))

def correctStartCodonFinder(nucleotideSeq, lastExonLen, potStartCodonStarts, distToJunc, strand, codonsPos, codonsNeg):
    '''finds the correct start codon for a nucleotide sequence (if a start codon causes a PTC between the start codon and the AS junction, it is not correct)'''
    if strand == '-':
        nucleotideSeq = nucleotideSeq[::-1]
        codons = [x[::-1] for x in codonsNeg]
    else:
        codons = codonsPos
    for i in potStartCodonStarts:
        potSeq = nucleotideSeq[i:distToJunc]
        codonList = [potSeq[0+i:3+i] for i in range(0, len(potSeq), 3)]
        if len(set(codonList).intersection(codons)) > 0:       
            pass
        else:
            finalSeq = nucleotideSeq[i:len(nucleotideSeq)-lastExonLen-50]
            if strand == '-':
                finalSeq = finalSeq[::-1]
            return(finalSeq)

def allInOneNMDFinder(strand, nucleotideSeq, codonsPos, codonsNeg):
    '''from the nucleotide sequence with the correct start, if there is a stop codon prior to the last exon-exon junction + 50 nt, it likely goes through NMD'''
    try:
        if strand == '-':
            nucleotideSeq = nucleotideSeq[::-1]
            codons = [x[::-1] for x in codonsNeg]
        else:
            codons = codonsPos
        codonList = [nucleotideSeq[0+i:3+i] for i in range(0, len(nucleotideSeq), 3)]
        if strand == '-':
            codonList = codonList[::-1]
        codonList = set(codonList)
        if len(set(codonList).intersection(codons)) > 0:
            return('NMD')
        else:
            return('not NMD')
    except:
        return('')
        
def postFixer(gene, finalNucleotideSeq, seqStarts, verdict):
    '''There are some AS junctions where the correctStartCodonFinder fails because the junction occurs before the healthy start codon. These are deemed to go through NMD'''
    if (verdict == '') & (finalNucleotideSeq is None) & (len(seqStarts) > 0):
        return('NMD')
    else:
        return(verdict)
    
def stopCodonIdentifier(strand, nucleotideSeq, codonsPos, codonsNeg):
    '''identifies codons of interest in the extended sequences and returns their index'''
    if strand == '-':
        nucleotideSeq = nucleotideSeq[::-1]
    codonList = [] 
    while len(nucleotideSeq) > 0:
        try:
            codon = nucleotideSeq[:3]
            nucleotideSeq = nucleotideSeq[3:]
        except:
            codon = nucleotideSeq[0:]
        codonList.append(codon)
    outputList = []
    if strand == '+':
        codons = codonsPos
    elif strand == '-':
        codonsNeg = [x[::-1] for x in codonsNeg]
        codons = codonsNeg
    for i in range(len(codonList)):
        if codonList[i] in codons:
            outputList.append(i)
    if len(outputList) == 0:
        return(False)
    else:
        return(outputList)

stopCodonsPos = ['TAG', 'TAA', 'TGA']  
stopCodonsNeg = ['CTA', 'TTA', 'TCA']  

'''call the previous functions to find PTCs in our junction file'''
dfList = []
chrs = ['chr' + str(x) for x in range(1,23)]
start_time = time.time()
for chrZ in chrs:
    curr_geneRef_Df = geneRef_df[geneRef_df['chromosome'] == chrZ]
    chrZ_df = anno_df[anno_df['chr'] == chrZ]

    '''get the relevant reference file sequence'''
    recordZ = records[chrs.index(chrZ)]
    '''this is in case that the reference file is not in order for whatever
    reason'''
    if recordZ.id != chrZ:
        for record in records:
            if record.id == chrZ:
                recordZ = record
                break
    seqZ = recordZ.seq

    '''now find PTCs in every junction'''
    chrZ_df[['bpString', 'lastExonLen', 'potStarts', 'distToJunc']] = chrZ_df.apply(lambda x: nucleotideSeqFinder(x['start'], x['end'], x['strand'], seqZ, x['gene'], intron_df, x['curTag_skipping'], x['relStartExon'], x['relEndExon'], geneRef_df), axis=1)
    chrZ_df['finalSeq'] = chrZ_df.apply(lambda x: correctStartCodonFinder(x['bpString'], x['lastExonLen'], x['potStarts'], x['distToJunc'], x['strand'], stopCodonsPos, stopCodonsNeg), axis=1)
    chrZ_df['NMDVerdict'] = chrZ_df.apply(lambda x: allInOneNMDFinder(x['strand'], x['finalSeq'], stopCodonsPos, stopCodonsNeg), axis=1)
    chrZ_df['NMDVerdict'] = chrZ_df.apply(lambda x: postFixer(x['gene'], x['finalSeq'], x['potStarts'], x['NMDVerdict']), axis=1)
    dfList.append(chrZ_df)
    print('chromosome: %s' %chrZ)
    print('Time elapsed: %ss' %round((time.time() - start_time), 0))
output_df = pd.concat(dfList)

del output_df['bpString']
del output_df['lastExonLen']
del output_df['potStarts']
del output_df['distToJunc']
del output_df['NMDVerdict']

output_df.to_csv(output_file)
