"""read reference genome fasta file with seek
"""
import os
import logging

logger = logging.getLogger()


class Chromosome:
    """Chromesome
    """
    def __init__(self, f, index, linelen):
        self.f = f
        self.index = index
        self.linelen = linelen

    def __len__(self):
        return self.index[1] - self.index[0]

    def __getitem__(self, key):
        if isinstance(key, slice):
            slice_start = max(key.start, 0)
            slice_stop = min(key.stop, self.__len__())
        elif isinstance(key, int):
            if key < 0 or key >= self.__len__():
                raise IndexError('chromsome index out of range')
            slice_start, slice_stop = key, key + 1
        else:
            raise TypeError('Chromsome index must be integer or slices')

        length = slice_stop - slice_start
        if length <= 0:
            return ''

        first = slice_start % self.linelen
        start_offset = slice_start // self.linelen * (self.linelen + 1) + first

        lines = (length + first) // self.linelen
        last = (length + first) % self.linelen if lines else length

        self.f.seek(self.index[0] + start_offset)
        s = ''.join(self.f.readline().strip('\n') for i in range(lines)) + self.f.read(last)
        return s


class SequenceFileDB:
    """read reference genome fasta file with seek or load all into memory
    """
    def __init__(self, filename, load_all=False):
        if load_all:
            self._read_all(filename)
            self.file = None
        else:
            offset_filename = f'{filename}.offset'
            if not os.path.exists(offset_filename):
                self._create_offset(filename, offset_filename)

            self._read_offset(offset_filename)
            self.file = open(filename)
            self.chrom = None

    def _read_all(self, filename):
        self.chrom = dict()
        with open(filename) as f:
            for line in f:
                if line.startswith('>chr'):
                    chrom_name = line[1:].strip()
                    self.chrom[chrom_name] = []
                else:
                    self.chrom[chrom_name].append(line.strip())
        for chrom_name in self.chrom:
            self.chrom[chrom_name] = ''.join(self.chrom[chrom_name])

    def _read_offset(self, filename):
        self.chrom_idx = dict()
        with open(filename) as f:
            self.linelen = int(f.readline())
            for line in f:
                key, start, end = line.split('\t')
                self.chrom_idx[key] = (int(start), int(end))

    def _create_offset(self, infile, outfile):
        logger.info('creating ucsc.hg19.fasta.offset ...')
        linelen = None
        lines = []
        with open(infile) as f:
            line = f.readline()
            chrom, start_offset = None, None
            while line:
                if line.startswith('>chr'):
                    if chrom:
                        end_offset = f.tell() - len(line)
                        lines.append('{}\t{}\t{}'.format(chrom, start_offset, end_offset))

                    start_offset = f.tell()
                    chrom = line[1:].strip()

                elif not linelen:
                    linelen = len(line.strip())

                line = f.readline()

            if chrom:
                end_offset = f.tell()
                lines.append('{}\t{}\t{}'.format(chrom, start_offset, end_offset))

        lines = [str(linelen)] + lines

        try:
            with open(outfile, 'w') as fout:
                fout.write('\n'.join(lines))
        except Exception:
            os.unlink(outfile)
            raise

    def __getitem__(self, chrom):
        if self.chrom:
            return self.chrom[chrom]
        return Chromosome(self.file, self.chrom_idx[chrom], self.linelen)

    def __del__(self):
        if self.file:
            self.file.close()
