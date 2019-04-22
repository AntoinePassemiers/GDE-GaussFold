# -*- coding: utf-8 -*-
# parsers.py
# author : Antoine Passemiers

import os
import warnings
import numpy as np


class Parser:

    def parse(self, filepath):
        _, file_ext = os.path.splitext(filepath)
        return self.__parse__(filepath)


class SS3Parser(Parser):

    def __init__(self, target_indices=[3, 4, 5]):
        self.target_indices = target_indices

    def __parse__(self, filepath):
        data = list()
        with open(filepath, 'r') as f:
            lines = f.readlines()
        for line in lines:
            if line.startswith('#'):
                pass
            elif len(line) > 2:
                els = line.split()
                data.append([float(els[i]) for i in self.target_indices])
        return np.asarray(data)


class ContactParser(Parser):

    def __init__(self, sequence_length, delimiter=' ', target_cols=[2]):
        self.sequence_length = sequence_length
        self.delimiter = delimiter
        self.target_cols = target_cols

    def is_comment(self, elements):
        is_comment = False
        for element in elements:
            try:
                float(element)
            except ValueError:
                is_comment = True
                break
        return is_comment

    def get_features(self, elements):
        features = list()
        if self.target_cols:
            for i in self.target_cols:
                features.append(float(elements[i]))
        else:
            target_col = len(elements) - 1
            features.append(float(elements[target_col]))
        return features

    def __parse__(self, filepath):
        with open(filepath, 'r') as f:
            lines = f.readlines()
            n_features = 1 if not self.target_cols else len(self.target_cols)
            data = np.full((n_features, self.sequence_length, self.sequence_length), np.nan, dtype=np.double)
            data[:, np.arange(self.sequence_length), np.arange(self.sequence_length)] = 0
            for line in lines:
                if self.delimiter:
                    elements = line.rstrip('\r\n').split(self.delimiter)
                else:
                    elements = line.rstrip('\r\n').split()
                if not self.is_comment(elements):
                    i, j = int(elements[0]) - 1, int(elements[1]) - 1
                    data[:, i, j] = data[:, j, i] = self.get_features(elements)
            if np.isnan(data).any():
                # sequence_length is wrong or the input file has missing pairs
                warnings.warn('Warning: Pairs of residues are missing from the contacts text file')
                warnings.warn('Number of missing pairs: %i ' % np.isnan(data).sum())
            return np.squeeze(data)


class FastaParser(Parser):

    def __parse__(self, filepath):
        with open(filepath, 'r') as f:
            data = { 'sequences' : list(), 'comments' : list() }
            lines = f.readlines()
            sequence = ""
            for line in lines:
                if line[0] == '>':
                    data['comments'].append(line[1:].rstrip("\r\n"))
                    if len(sequence) > 0:
                        sequence = sequence.rstrip("\r\n")
                        data['sequences'].append(sequence)
                        sequence = ''
                else:
                    line = line.replace('\n', '').replace('\r', '').strip()
                    if len(line) > 1:
                        sequence += line
            sequence = sequence.rstrip("\r\n")
            data['sequences'].append(sequence)
            return data




RES_NAME_TO_SYM = {
    'ALA': 'A',
    'ARG': 'R',
    'ASN': 'N',
    'ASP': 'D',
    'CYS': 'C',
    'GLU': 'E',
    'GLN': 'Q',
    'GLY': 'G',
    'HIS': 'H',
    'ILE': 'I',
    'LEU': 'L',
    'LYS': 'K',
    'MET': 'M',
    'MSE': 'M',
    'PHE': 'F',
    'PRO': 'P',
    'SER': 'S',
    'THR': 'T',
    'TRP': 'W',
    'TYR': 'Y',
    'VAL': 'V'
}


class ResidueId:

    def __init__(self, idf):
        self.idf = idf.strip()
        self.number = ''
        self.suffix = ''
        for c in self.idf:
            if c.isdigit() or c in ['+', '-']:
                self.number += c
            else:
                self.suffix += c
        self.number = int(self.number)
        assert(len(self.suffix) <= 1)

    def __float__(self):
        f = self.number * 1000
        if self.suffix:
            f += ord(self.suffix)
        return f

    def __le__(self, other):
        return self.__float__() <= other.__float__()

    def __ge__(self, other):
        return self.__float__() >= other.__float__()

    def __eq__(self, other):
        return self.__float__() == other.__float__()

    def __ne__(self, other):
        return self.__float__() != other.__float__()

    def __lt__(self, other):
        return self.__float__() < other.__float__()

    def __gt__(self, other):
        return self.__float__() > other.__float__()
    
    def __hash__(self):
        return hash(self.__float__())


class PDBParser(Parser):

    def __init__(self, sequence, prot_name, method='CASP'):
        self.sequence = sequence
        self.L = len(sequence)
        self.prot_name = prot_name
        self.method = method

    def isAppropriateAtom(self, res_name, atom):
        valid = False
        if self.method == 'CA' and atom == 'CA':
            valid = True
        elif self.method == 'CB' and atom == 'CB':
            valid = True
        elif self.method == 'CASP':
            if res_name == 'GLY' and atom == 'CA':
                valid = True
            elif res_name != 'GLY' and atom == 'CB':
                valid = True
        return valid

    def align_query_sequence(self, query_seq, whole_seq):
        import minineedle # https://github.com/scastlara/minineedle
        import miniseq # https://github.com/scastlara/miniseq
        if len(query_seq) != len(whole_seq):
            query_seq = miniseq.Protein('', query_seq)
            whole_seq = miniseq.Protein('', whole_seq)
            alignment = minineedle.Needleman(query_seq, whole_seq)
            alignment.align()
            assert(len(whole_seq) == len(alignment.alseq2))
            indices = [i for i, res_name in enumerate(alignment.alseq1) if res_name != '-']
        else:
            indices = [i for i, res_name in enumerate(query_seq) if res_name != '-']
        return indices

    def res_name_to_sym(self, res_name):
        try:
            sym = RES_NAME_TO_SYM[res_name]
        except KeyError:
            sym = 'X' # TODO
        return sym

    def __parse__(self, filepath):
        if self.prot_name[0] == 'T' and '-' in self.prot_name:
            het_chain = None
        else:
            het_chain = self.prot_name[-1]
        distances = np.full((self.L, self.L), np.nan, dtype=np.float)
        
        residues = dict()
        with open(filepath, 'r') as f:
            lines = f.readlines()
        for line in lines:
            if line[:10] == 'REMARK 465' and line[19] == het_chain:
                res_name = line[15:19].strip()
                res_id = line[21:27].strip()
                if any(c.isdigit() for c in res_id):
                    res_id = ResidueId(res_id)
                    residues[res_id] = (res_name, None)

        for line in lines:
            if (line[:4] == 'ATOM' or line[:6] == 'HETATM') and (line[21] == het_chain or het_chain is None):
                res_id = ResidueId(line[23:27].strip())
                res_name = line[17:20].strip()
                atom = line[13:16].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                if self.isAppropriateAtom(res_name, atom):
                    residues[res_id] = (res_name, np.asarray([x, y, z]))
                elif atom == 'CA' and not res_id in residues.keys():
                    residues[res_id] = (res_name, np.asarray([x, y, z]))

        residues = [x[1] for x in list(sorted(residues.items(), key=lambda kv: kv[0]))]
        whole_seq = ''.join([self.res_name_to_sym(pos[0]) for pos in residues])
        query_seq = str(self.sequence)
       
        indices = self.align_query_sequence(query_seq, whole_seq)
        coordinates = [residues[i][1] for i in indices]
        assert(len(coordinates) == self.L)

        for i, coords_i in enumerate(coordinates):
            for j, coords_j in enumerate(coordinates):
                if coords_i is not None and coords_j is not None:
                    distances[i, j] = np.sqrt(np.sum((coords_i - coords_j) ** 2.))
        return distances, np.asarray(coordinates)
