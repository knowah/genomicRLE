import gzip

class ChromRLE:
    def __init__(self, chrom):
        self.chromosome = chrom
        self.starts = []
        self.lengths = []
    
    def _within_run(self, idx, pos):
        return pos >= self.starts[idx] and pos <= (self.starts[idx] + self.lengths[idx] - 1)

    def contains(self, pos):
        n_runs = len(self.starts)
        if n_runs == 0:
            return False
        if n_runs == 1:
            return self._within_run(0, pos)

        # binary search for where this position would be 
        s = 0
        e = n_runs-1
        while s < e:
            mid = int((s+e)/2)
            if pos < self.starts[mid]:
                e = mid-1
                continue
            if pos == self.starts[mid] or self._within_run(mid, pos):
                return True
            else:
                s = mid+1
        if s == e:
            return self._within_run(s, pos)

        # not sure this will ever happen, but check both in this instance?
        # TODO: check if this ever e > s and figure what to do efficiently if so
        return self._within_run(s, pos) or self._within_run(e, pos)

    def append_run(self, start, length, check=True):
        if check and len(self.starts) > 0:
            if start < self.starts[-1]:
                raise ValueError("Attempting to append run out of order.")
            if self._within_run(len(self.starts)-1, start):
                raise ValueError("Attempting to append overlapping run.")
        self.starts.append(start)
        self.lengths.append(length)

    def extend_runs(self, start_lst, len_lst):
        self.starts.extend(start_lst)
        self.lengths.extend(len_lst)
        
        
class GenomicRLE:
    def __init__(self, rlefile, strict=False):
        self.strict = strict
        self.chrom_dict = {}
        
        buffer_size = 1000
        curr_buff = 0
        entry_buffer = [None] * buffer_size
        
        def eval_buffer(b, sz, chrom):
            st = [int(x[0]) for x in b[0:sz]]
            ln = [int(x[1]) for x in b[0:sz]]
            self.chrom_dict[chrom].extend_runs(st, ln)
        
        fopen = open
        if rlefile.endswith(".gz"):
            fopen = gzip.open
        
        with fopen(rlefile, 'rt') as rlef:
            curr_chrom = None
            for line in rlef:
                line = line.rstrip()
                if line.startswith("chrom="):
                    if curr_buff > 0:
                        eval_buffer(entry_buffer, curr_buff, curr_chrom)
                        entry_buffer = [None] * buffer_size
                        curr_buff = 0
                    curr_chrom = line.split('=')[1]
                    if curr_chrom in self.chrom_dict:
                        raise ValueError("RLE file {}: Chromsome {} seen multiple times.".format(rlefile, curr_chrom))
                    self.chrom_dict[curr_chrom] = ChromRLE(curr_chrom)
                else:
                    if curr_buff == buffer_size:
                        eval_buffer(entry_buffer, buffer_size, curr_chrom)
                        entry_buffer = [None] * buffer_size
                        curr_buff = 0
                    entry_buffer[curr_buff] = line.split('\t')
                    curr_buff += 1
                    #pos_start, run_len = line.split('\t')
                    #self.chrom_dict[curr_chrom].append_run(int(pos_start), int(run_len), False)


    def __getitem__(self, key):
        if key not in self.chrom_dict:
            if self.strict:
                raise KeyError("Chromosome {} not present in this GenomicRLE".format(key))
            else:
                return ChromRLE(key)
        return self.chrom_dict[key]

