#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import re
import numpy as np
import math
from Bio import SeqIO
import time
import sys


'''code for finding frameshifts'''
def frameShiftFinder(fivep_dist, threep_dist):
    '''checks if there is a frameshift. That is, if the sum of three prime distance and tive prime distance is divisible by 3, it is not a frameshift'''
    if (abs(fivep_dist) + abs(threep_dist))%3 == 0:
        return('not_frame_shift')
    else:
        return('frame_shift')
        
'''find non NMD causing PTCs'''
def nucleotideSeqFinderNonAs(strand, sequence, gene, intronDataFrame, tag, exonDataFrame):
    '''get the healthy nucleotide sequence for a gene'''
    try:
        intronStarts = list(intronDataFrame[1][(intronDataFrame[3] == gene) & (intronDataFrame[9] == tag)])
        intronEnds = list(intronDataFrame[2][(intronDataFrame[3] == gene) & (intronDataFrame[9] == tag)])

        intronEnds = [x-1 for x in intronEnds]
        
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
        if strand == '-':
            last_exon_len = abs(exons[0][1] - exons[0][0])
        elif strand == '+':
            last_exon_len = abs(exons[-1][1] - exons[-1][0])
        if strand == '+':
            seqStarts = [m.start(0) for m in re.finditer('ATG', outputSeq[:-last_exon_len])]
        elif strand == '-':
            seqStarts = [m.start(0) for m in re.finditer('TAC', outputSeq[:-last_exon_len])]
        if len(seqStarts) == 0:
            if strand == '+':
                seqStarts = [m.start(0) for m in re.finditer('ATG', outputSeq)]
            elif strand == '-':
                seqStarts = [m.start(0) for m in re.finditer('TAC', outputSeq)]

        if strand == '-':
            outputSeq = outputSeq[::-1]
        
        return(pd.Series([str(outputSeq), last_exon_len, seqStarts]))
            
    except:
        return(pd.Series(['', '', '']))

def correctStartCodonFinderNonAs(nucleotideSeq, lastExonLen, potStartCodonStarts, strand, codonsPos, codonsNeg):
    '''get the correct start codon for a healthy gene'''
    try:
        if strand == '-':
            nucleotideSeq = nucleotideSeq[::-1]
            codons = [x[::-1] for x in codonsNeg]
        else:
            codons = codonsPos
        for i in potStartCodonStarts:
            potSeq = nucleotideSeq[i:-(lastExonLen)]
            codonList = [potSeq[0+i:3+i] for i in range(0, len(potSeq), 3)]
            if len(set(codonList).intersection(codons)) > 0:       
                pass
            else:
                finalSeq = nucleotideSeq[i:]
                if strand == '-':
                    finalSeq = finalSeq[::-1]
                return(finalSeq)
    except:
        return('broke at correctStartCodonFinderNonAs')

def codonFinder(nucleotideSeq):
    '''take a nucleotide sequence and returns codons'''
    codonList = [nucleotideSeq[0+i:3+i] for i in range(0, len(nucleotideSeq), 3)]
    return(codonList)

def correctStopCodonFinderNonAs(finalSeq, strand, codonsPos, codonsNeg):
    '''finds the correct stop codon for a healthy sequence'''
    try:
        if strand == '-':
            finalSeq = finalSeq[::-1]
            codonsNeg = [x[::-1] for x in codonsNeg]
        codons = codonFinder(finalSeq)
        if strand == '+':
            correctStopCodonIndex = min([codons.index(x) for x in codonsPos if x in codons])
        elif strand == '-':
            correctStopCodonIndex = min([codons.index(x) for x in codonsNeg if x in codons])
        nucleotideSeqToMatch = ''.join(codons[correctStopCodonIndex:correctStopCodonIndex+4])
        if strand == '-':
            nucleotideSeqToMatch = nucleotideSeqToMatch[::-1]
        return(nucleotideSeqToMatch)
    except:
        return('broke at correctStopCodonFinderNonAs')

def nucleotideSeqFinderAs(start, end, strand, sequence, gene, intronDataFrame, tag, startExon, endExon, exonDataFrame):
    '''generates the nucleotide sequence of a gene with a given AS junction'''
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
    '''find the correct start codon'''
    try:
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
                finalSeq = nucleotideSeq[i:]
                if strand == '-':
                    finalSeq = finalSeq[::-1]
                return(finalSeq)
    except:
        return('broke at correctStartCodonFinder')

def checkCodonMatch(bigCodonList, smallCodonList):
    '''checks if a list of codons contains a smaller list of 4 codons'''
    try:
        for i in range(len(bigCodonList) - len(smallCodonList) + 1):
            if bigCodonList[i] == smallCodonList[0]:
                try:
                    if bigCodonList[i+1] == smallCodonList[1]:
                        try:
                            if bigCodonList[i+2] == smallCodonList[2]:
                                try:
                                    if bigCodonList[i+3] == smallCodonList[3]:
                                        return(i)
                                except:
                                    return(i)
                        except:
                            return(i)
                except:
                    return(i)
        return('No_match')
    except:
        return('broke at checkCodonMatch')
        
def nonNmdTriggeringPtcsFinder(NmdVerdict, frameShiftVerdict, finalSeqAs, ntSeqToMatch, lastExonLenAs, strand, stopCodonsPos, stopCodonsNeg):
    '''checks the region of the nt rule for premature stop codons'''
    try:
        if strand == '-':
            finalSeqAs = finalSeqAs[::-1]
            ntSeqToMatch = ntSeqToMatch[::-1]
        if (NmdVerdict == 'not NMD') & (frameShiftVerdict == 'not_frame_shift'):
            codonsAs = codonFinder(finalSeqAs)
            codonsOfInterest = codonsAs[-math.floor((lastExonLenAs + 50)/3):]
            ntSeq = ''.join(codonsOfInterest)
            if re.search(ntSeqToMatch, ntSeq):
                if strand == '+':
                    stopCodon = stopCodonsPos
                elif strand == '-':
                    stopCodon = [x[::-1] for x in stopCodonsNeg]
                firstStopCodon = min([codonsOfInterest.index(x) for x in stopCodon if x in codonsOfInterest])
                correctCodons = codonFinder(ntSeqToMatch)
                matchIndex = checkCodonMatch(codonsOfInterest, correctCodons)
                if firstStopCodon < matchIndex:
                    return('PTC_not_causing_NMD')
                else:
                    return('No_PTC')
            else:
                return('edge_case')
        else:
            return(None)
    except:
        return('edge_case')
    
'''import files'''

junction_file = sys[1]
fasta_file = sys[2]
gene_ref_file = sys[3]
intron_ref_file = sys[4]
output_file = sys[5]

#file = '/Users/lkluegel/Documents/Splicing/Fede_Ally_Paulina/NMD/logOR_within_cell_type_ONLY_CRYPTIC_3P_Junctions_with_threshold_info_w_NMD_faster_with_start_codon_correction.csv'
df = pd.read_csv(file, index_col=0)

#ref_file = '/Users/lkluegel/Documents/Splicing/Fede_Ally_Paulina/Annotator/Data/GRCh38.p12.genome.fa'
records = list(SeqIO.parse(ref_file, 'fasta'))

#geneRef_file = '/Users/lkluegel/Documents/Splicing/Fede_Ally_Paulina/Annotator/Data/gencode.v31.basic.annotation.csv'
geneRef_df = pd.read_csv(geneRef_file, index_col=False)
geneRef_df = geneRef_df[geneRef_df['type'] == 'transcript']

#intron_ref_file = '/Users/lkluegel/Documents/Splicing/Fede_Ally_Paulina/Annotator/Data/leafviz_all_introns_cleaned.bed'
intron_df = pd.read_csv(intron_ref_file, sep='\t', header=None)     

stopCodonsPos = ['TAG', 'TAA', 'TGA']  
stopCodonsNeg = ['CTA', 'TTA', 'TCA']  

'''call the previous functions to find PTCs in our junction file'''
dfList = []
chrs = ['chr' + str(x) for x in range(1,23)]

start_time = time.time()
for chrZ in chrs:
    curr_geneRef_Df = geneRef_df[geneRef_df['chromosome'] == chrZ]
    chrZ_df = df[df['chr'] == chrZ]

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
    chrZ_df['frameShift'] = chrZ_df.apply(lambda x: frameShiftFinder(x['fivep_distance'], x['threep_distance']), axis=1)
    chrZ_df[['bpStringNonAs', 'lastExonLenNonAs', 'potStartsNonAs']] = chrZ_df.apply(lambda x: nucleotideSeqFinderNonAs(x['strand'], seqZ, x['gene'], intron_df, x['curTag_skipping'], geneRef_df), axis=1)
    chrZ_df['finalSeqNonAs'] = chrZ_df.apply(lambda x: correctStartCodonFinderNonAs(x['bpStringNonAs'], x['lastExonLenNonAs'], x['potStartsNonAs'], x['strand'], stopCodonsPos, stopCodonsNeg), axis=1)
    chrZ_df['nucleotideSeqToMatch'] = chrZ_df.apply(lambda x: correctStopCodonFinderNonAs(x['finalSeqNonAs'], x['strand'], stopCodonsPos, stopCodonsNeg), axis=1)
    chrZ_df[['bpStringAs', 'lastExonLenAs', 'potStartsAs', 'distToJuncAs']] = chrZ_df.apply(lambda x: nucleotideSeqFinderAs(x['start'], x['end'], x['strand'], seqZ, x['gene'], intron_df, x['curTag_skipping'], x['relStartExon'], x['relEndExon'], geneRef_df), axis=1)
    chrZ_df['finalSeqAs'] = chrZ_df.apply(lambda x: correctStartCodonFinder(x['bpStringAs'], x['lastExonLenAs'], x['potStartsAs'], x['distToJuncAs'], x['strand'], stopCodonsPos, stopCodonsNeg), axis=1)
    chrZ_df['nonNmdPtc'] = chrZ_df.apply(lambda x: nonNmdTriggeringPtcsFinder(x['NMDVerdict'], x['frameShift'], x['finalSeqAs'], x['nucleotideSeqToMatch'], x['lastExonLenAs'], x['strand'], stopCodonsPos, stopCodonsNeg), axis=1)
    dfList.append(chrZ_df)
    print('chromosome: %s' %chrZ)
    print('Time elapsed: %ss' %round((time.time() - start_time), 0))
    
output_df = pd.concat(dfList)

output_df.drop(['bpStringNonAs','lastExonLenNonAs','potStartsNonAs','finalSeqNonAs','nucleotideSeqToMatch','bpStringAs','lastExonLenAs','potStartsAs','distToJuncAs','finalSeqAs'], inplace=True, axis=1)

#output_df.to_csv('/Users/lkluegel/Documents/Splicing/Fede_Ally_Paulina/NMD/logOR_within_cell_type_ONLY_CRYPTIC_3P_Junctions_with_threshold_info_w_NMD_faster_with_start_codon_correction_frame_shift_and_nonNmdPtc.csv')

output_df.to_csv(output_file)
