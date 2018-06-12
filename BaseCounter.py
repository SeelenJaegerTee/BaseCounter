import linecache
import os
import datetime
from io import *

# ---------------------------------------------------------------------------------------------------------------------
# Global Variables
# ---------------------------------------------------------------------------------------------------------------------
# path = 'C:/Users/NeugebauerC/Desktop/20170630_run187/Analysen_JenB/'
path = 'C:/Users/NeugebauerC/Desktop/python kurs/BaseCounter/'
out_path = path + "output.txt"
log_path = path + 'log.txt'
fileEndings = ['.fasta']
bases = 'ACGTUWSMKRYBDHVN'
bases += bases.lower()


# ---------------------------------------------------------------------------------------------------------------------
# Classes
# ---------------------------------------------------------------------------------------------------------------------
class Config:
    DEBUG = True
    CHUNK_SIZE = 65535


# can provide lines from the file by buffering a large chunk of the file
# should prevent filling memory in case input file is incredibly large
# warning when requesting a line and if the file has no line-breaks (\n)
# the entire file will be successively stored into a string
# at this point it would probably have been faster to just readlines the entire file
class ChunkProvider:
    def __init__(self, file_path):
        self.file_path = file_path
        self.pos = 0
        self.chunk_buffer = ''
        self.EOF_reached = False
        self.file_exhausted = False

        self.line = 0
        self.chunk_pos = 0
        self.read_chunk()

    def read_chunk(self, chunk_size=Config.CHUNK_SIZE):
        # if Config.DEBUG: print('before call self.chunk_buffer', repr(self.chunk_buffer))
        if not self.EOF_reached:
            with open(self.file_path) as fp:
                fp.seek(self.pos)
                self.chunk_buffer = fp.read(chunk_size)
                # if Config.DEBUG is True: print('chunk loaded: ' + repr(self.chunk_buffer))
                self.pos = fp.tell()
                # check whether EOF is reached
                if self.pos == fp.seek(0, SEEK_END):
                    self.EOF_reached = True
                    if Config.DEBUG is True:
                        print('EOF reached')
                # we need to reset the filepointer to where it was
                fp.seek(self.pos)
        else:
            self.chunk_buffer = ''
            if Config.DEBUG is True:
                print('requested to read chunk despite EOF reached, returning empty string')
        # if Config.DEBUG: print('after call self.chunk_buffer', repr(self.chunk_buffer))

    def next_line(self):
        new_line_pos = self.chunk_buffer.find('\n', self.chunk_pos)
        if Config.DEBUG: print('new_line_pos', new_line_pos)
        if new_line_pos != -1:
            ret = self.chunk_buffer[self.chunk_pos: new_line_pos+1]
            self.line += 1
            if Config.DEBUG: print('found line')
            self.chunk_pos = new_line_pos + 1
        elif not self.EOF_reached:
            if Config.DEBUG: print('still more file to buffer ...')
            ret = self.chunk_buffer[self.chunk_pos: ]
            if Config.DEBUG: print('reloading')
            self.read_chunk()
            self.chunk_pos = 0
            ret += self.next_line()
        elif not self.file_exhausted:
            if Config.DEBUG: print('last bit of file already loaded in buffer')
            ret = self.chunk_buffer[self.chunk_pos:]
            self.line += 1
            self.file_exhausted = True
        else:
            if Config.DEBUG: print('finished')
            return ''
        return ret


# class LineProvider:
#     def __init__(self, file_path):
#         self.line = 0
#         self.file_path = file_path
#         self.pos = 0
#
#     def next_line(self):
#         with open(self.file_path) as fp:
#             fp.seek(self.pos)
#             rtrn = fp.readline()
#             self.line += 1
#             self.pos = fp.tell()
#             return rtrn


class File:
    def __init__(self, file_path, name):
        self.name = name
#        self.lp = LineProvider(file_path)
        self.lp = ChunkProvider(file_path)
        self.sequences = []

    def analyze(self):
        counter = 0
        nuc_types = {'-': 0}
        for c in bases:
            nuc_types[c] = 0
        position = 0
        name = ''
        buff = self.lp.next_line()
        while buff != '':
            if buff[0] is '>':
                if counter > 0:
                    self.sequences.append(Sequence(name, counter, nuc_types))
                counter = 0
                nuc_types = {'-': 0}
                for c in bases:
                    nuc_types[c] = 0
                position = 0
                name = buff
            if buff[0] is not '>':
                for c in buff:
                    position += 1
                    if c in bases:
                        counter += 1
                        nuc_types[c] += 1
                    elif c is '-':
                        nuc_types['-'] += 1
                    elif c is '\n':
                        position = 0
                    else:
                        if c in nuc_types:
                            nuc_types[c] += 1
                        else:
                            nuc_types[c] = 1
                        warn_string = 'encountered non-nucleotide: ' + repr(c) + ' in ' + self.name \
                                      + ' \tline: ' + str(self.lp.line) + ' \tposition: ' + str(position) \
                                      + ' \t...' + buff[(position - 5):(position + 4)] + '...\n'
                        with open(log_path, "a") as log:
                            log.write(warn_string)
            buff = self.lp.next_line()
        self.sequences.append(Sequence(name, counter, nuc_types))

    def print_size(self, verbose=True):
        print(self.name)
        total_size = 0
        for seq in self.sequences:
            total_size = total_size + seq.bases
            if verbose:
                seq.print_size()
        print('total size: \t', total_size / 1000000, 'Mbp')
        print('---------------------------------------------')

    def write_size(self, out_file, verbose=True):
        out_file.write(self.name + '\n')
        total_size = 0
        for seq in self.sequences:
            total_size = total_size + seq.bases
            if verbose:
                seq.write_size(out_file)
        out_string = '\n' + 'total size: \t' + str(total_size / 1000000) \
                     + ' Mbp' + '\t\t' + str(total_size) + ' bp' \
                     + '\n----------------------------------------------------------\n\n'
        out_file.write(out_string)


class Sequence:
    def __init__(self, name, base_count, nuc_types):
        self.name = name.rstrip()
        self.bases = base_count
        self.nuc_types = nuc_types

    def print_size(self):
        print(self.name, ':\t', self.bases / 1000000, ' Mbp')

    def write_size(self, out_file):
        out_string = '\n' + self.name + ':\t' + str(self.bases / 1000000) + ' Mbp' + '\t\t' + str(self.bases) + ' bp\n'
        for key in self.nuc_types:
            if key is '-':
                out_string += key + ': ' + str(self.nuc_types[key]) + ' (not counted in bp)\n'
            elif key is 'G' or key is 'C' or key is 'A' or key is 'T':
                out_string += key+': '+str(self.nuc_types[key])+'\n'
            elif self.nuc_types[key] != 0:
                out_string += key + ': ' + str(self.nuc_types[key]) + '\n'
        out_file.write(out_string)


# ---------------------------------------------------------------------------------------------------------------------
# Script
# ---------------------------------------------------------------------------------------------------------------------
files_index = os.listdir(path)
t1 = datetime.datetime.now()
with open(out_path, 'w') as fp:
    fp.write('%s\n\n' % datetime.datetime.now())

with open(log_path, 'w') as fp:
    fp.write('%s\n\n' % datetime.datetime.now())

tmp = []
for filename in files_index:
    for ending in fileEndings:
        if ending in filename:
            tmp.append(filename)
files_index = tmp
with open(out_path, 'a') as fp:
    fp.write('evaluated files: ' + str(files_index) + '\n\n\n')
inFiles = []
for filename in files_index:
    filePath = path + filename
    inFiles.append(File(filePath, filename))

for file in inFiles:
    file.analyze()
    file.print_size(verbose=False)
print('\ndetailed:\n')
for file in inFiles:
    file.print_size(verbose=True)

for file in inFiles:
    print('writing to file')
    with open(out_path, "a") as fp:
        file.write_size(fp, verbose=False)
with open(out_path, 'a') as fp:
    fp.write('\n#############\n# detailed: #\n#############\n\n')
for file in inFiles:
    with open(out_path, "a") as fp:
        file.write_size(fp, verbose=True)
t2 = datetime.datetime.now()
with open(out_path, "a") as fp:
    fp.write('\n\n time to finish:'+str(t2 - t1))
print('job is finished, have a nice day!')
