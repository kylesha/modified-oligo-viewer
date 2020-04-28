class ModifiedOligo:
    '''
    A class to facilitate viewing of modified bases present in an oligonucleotide. It does this by
        (1) Replacing the modified base position(s) with a user-defined symbol.
        (2) Colorizing each DNA base (and its modified base symbol) with the same color.

    Attributes
    ----------
    SeqRecord: Bio.SeqRecord
        A Biopython SeqRecord object
    modifications: dict
        A dictionary of the form {(symbol, name):(positions)}. It maps a user-defined symbol and the modified base name to the base position(s) in which the modified base occurs.
        The base position is zero-based (so as to be consistent with SeqRecord numbering convention in Biopython). For example, suppose we have an oligo (say ATGTCAGTC) in which the second and eighth T are deoxyUracil (dU) bases,
        and we wish to represent it using the pound '#' symbol. We would specify that as {('#', 'dU'):(1, 7)}
    '''

    #############| class constructor |#############
    def __init__(self, SeqRecord, modifications=None):
        self.__seq_rec = SeqRecord
        self.__modifications = modifications
        self.__oligo = [base for base in self.__seq_rec]
        self.__info_table = {}

        '''
        Map each modified base's symbol and its position(s) to the canonical base.
        The new dict has structure {(symbol, name, canonical_base):(positions)}
        '''
        for symbol, positions in self.__modifications.items():
            if len(positions) == 1:
                canonical_base = self.__seq_rec[positions[0]]
                self.__info_table.update({(symbol[0], symbol[1], canonical_base):positions})
            else:   # ensure symbol does NOT map to mulitple (different) canonical bases
                canonical_base = set([self.__seq_rec[i] for i in positions])
                if len(canonical_base) == 1:
                    self.__info_table.update({(symbol[0], symbol[1], canonical_base.pop()):positions})
                else:
                    print("ERROR!: \'{}\' is assigned to multiple canonical bases {}".format(symbol[0], canonical_base))

        # reconstruct oligo object: replacing modified base positions with use-defined symbols
        for symbol, positions in self.__info_table.items():
            if len(positions) == 1:
                self.__oligo[positions[0]] = symbol[0]
            else:
                for i in positions:
                   self.__oligo[i] = symbol[0]


    #############| private class methods |#############
    def __colorize(self, base):
        for key in self.__info_table.keys():
            if base.upper() == 'A' or ('A' in key and base in key):
                return ''.join(['\033[0;31m', base, '\033[0m'])  # red
            elif base.upper() == 'G' or ('G' in key and base in key):
                return ''.join(['\033[0;32m', base, '\033[0m'])  # green
            elif base.upper() == 'T' or ('T' in key and base in key):
                return ''.join(['\033[0;34m', base, '\033[0m'])  # blue
            elif base.upper() == 'C' or ('C' in key and base in key):
                return ''.join(['\033[0;33m', base, '\033[0m'])  # yellow

    def __legend(self):
        legend = []
        for key in self.__info_table.keys():
            if key[2] == 'A':
                legend.append(''.join(['\033[0;31m', ' = '.join([key[0], key[1]]), '\033[0m']))  # red
            elif key[2] == 'G':
                legend.append(''.join(['\033[0;32m', ' = '.join([key[0], key[1]]), '\033[0m']))  # green
            elif key[2] == 'T':
                legend.append(''.join(['\033[0;34m', ' = '.join([key[0], key[1]]), '\033[0m']))  # blue
            elif key[2] == 'C':
                legend.append(''.join(['\033[0;33m', ' = '.join([key[0], key[1]]), '\033[0m']))  # yellow
        return ''.join([i+'\n' for i in legend]).rstrip()


   #############| instance methods |#############
    def view53(self, modified=True, showLegend = False):
        # whether to view the unmodified (self.__seq_rec) or modified (self.__oligo) oligo
        oligo_colored = [self.__colorize(base) for base in self.__oligo] if modified else [self.__colorize(base) for base in self.__seq_rec]
        if showLegend:
            return ''.join([self.__legend(), '\n',
                            "5' ", ''.join(oligo_colored), " 3'"])
        else:
            return ''.join(["5' ", ''.join(oligo_colored), " 3'"])

    def view35(self, modified=True, showLegend = False):
        # oligo_colored = [self.__colorize(base) for base in self.__oligo]
        oligo_colored = [self.__colorize(base) for base in self.__oligo] if modified else [self.__colorize(base) for base in self.__seq_rec]
        if showLegend:
            return ''.join([self.__legend(), '\n',
                            "3' ", ''.join(oligo_colored[::-1]), " 5'"])
        else:
            return ''.join(["3' ", ''.join(oligo_colored[::-1]), " 5'"])
