import h5py
import numpy as np
from itertools import product


quadrants = np.array([
  # rotx=0, roty=0-3
  [[[0, 7], [1, 6]], [[3, 4], [2, 5]]],
  [[[7, 4], [6, 5]], [[0, 3], [1, 2]]],
  [[[4, 3], [5, 2]], [[7, 0], [6, 1]]],
  [[[3, 0], [2, 1]], [[4, 7], [5, 6]]],
  # rotx=1, roty=0-3
  [[[1, 0], [6, 7]], [[2, 3], [5, 4]]],
  [[[0, 3], [7, 4]], [[1, 2], [6, 5]]],
  [[[3, 2], [4, 5]], [[0, 1], [7, 6]]],
  [[[2, 1], [5, 6]], [[3, 0], [4, 7]]],
  # rotx=2, roty=0-3
  [[[6, 1], [7, 0]], [[5, 2], [4, 3]]],
  [[[1, 2], [0, 3]], [[6, 5], [7, 4]]],
  [[[2, 5], [3, 4]], [[1, 6], [0, 7]]],
  [[[5, 6], [4, 7]], [[2, 1], [3, 0]]],
  # rotx=3, roty=0-3
  [[[7, 6], [0, 1]], [[4, 5], [3, 2]]],
  [[[6, 5], [1, 2]], [[7, 4], [0, 3]]],
  [[[5, 4], [2, 3]], [[6, 7], [1, 0]]],
  [[[4, 7], [3, 0]], [[5, 6], [2, 1]]],
  # rotx=4, roty=0-3
  [[[6, 7], [5, 4]], [[1, 0], [2, 3]]],
  [[[7, 0], [4, 3]], [[6, 1], [5, 2]]],
  [[[0, 1], [3, 2]], [[7, 6], [4, 5]]],
  [[[1, 6], [2, 5]], [[0, 7], [3, 4]]],
  # rotx=5, roty=0-3
  [[[2, 3], [1, 0]], [[5, 4], [6, 7]]],
  [[[3, 4], [0, 7]], [[2, 5], [1, 6]]],
  [[[4, 5], [7, 6]], [[3, 2], [0, 1]]],
  [[[5, 2], [6, 1]], [[4, 3], [7, 0]]]
])

rotxmap_table = np.array([4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 17, 18, 19, 16, 23, 20, 21, 22])
rotymap_table = np.array([1, 2, 3, 0, 16, 17, 18, 19, 11, 8, 9, 10, 22, 23, 20, 21, 14, 15, 12, 13, 4, 5, 6, 7])

rotx_table = np.array([3, 0, 0, 2, 2, 0, 0, 1])
roty_table = np.array([0, 1, 1, 2, 2, 3, 3, 0])

sense_table = [-1, -1, -1, +1, +1, -1, -1, -1]


class EagleSnapshotClosedException(Exception):
    pass


class EagleSnapshot(object):
    """Class to represent an open Eagle snapshot"""

    def __init__(self, fname, verbose=False):
        """Open a new snapshot"""
        self.fname = fname
        self.verbose = verbose
        if self.verbose:
            print("__init__() called")
        self.sampling_rate = 1.0
        with h5py.File(self.fname, 'r') as f:
            if self.verbose:
                print("  - Opened file: {:s}".format(self.fname))
            self.boxsize = f['/Header'].attrs['BoxSize']
            self.numfiles = f['/Header'].attrs['NumFilesPerSnapshot']
            nptot = f['/Header'].attrs['NumPart_Total']
            nptot_hw = f['/Header'].attrs['NumPart_Total_HighWord']
            self.hashbits = f['/HashTable'].attrs['HashBits']
            self.ncell = 1 << self.hashbits
            self.nhash = 1 << 3 * self.hashbits
            self.numpart_total = [nptot[i] + nptot_hw[i] << 32 for i in range(6)]
            if self.verbose:
                print("  - Read in file header")
            # initialize hashmap all false
            self.hashmap = np.zeros(self.nhash)
            name_parts = self.fname.split('.')
            if name_parts[-1] != 'hdf5' or not name_parts[-2].isdigit():
                raise RuntimeError("Don't understand snapshot file name!")
            else:
                self.basename = ".".join(name_parts[:-2])
                if self.verbose:
                    print("  - Base name is {:s}".format(self.basename))
            if self.verbose:
                for itype in range(6):
                    if self.numpart_total[itype] > 0:
                        print("  - Have particles of type {:d}".format(itype))
            # These three datasets repeated in each file
            self.first_key_in_file = [
                f['/HashTable/PartType{:d}/FirstKeyInFile'.format(itype)][...]
                if self.numpart_total[itype] > 0 else None
                for itype in range(6)
            ]
            self.last_key_in_file = [
                f['/HashTable/PartType{:d}/LastKeyInFile'.format(itype)][...]
                if self.numpart_total[itype] > 0 else None
                for itype in range(6)
            ]
            self.num_keys_in_file = [
                f['/HashTable/PartType{:d}/NumKeysInFile'.format(itype)][...]
                if self.numpart_total[itype] > 0 else None
                for itype in range(6)
            ]
        self.num_part_in_file = [[None for ifile in range(self.numfiles)] for itype in range(6)]
        for ifile in range(self.numfiles):
            fname = '{:s}.{:d}.hdf5'.format(self.basename, ifile)
            with h5py.File(fname, 'r') as f:
                numpart_thisfile = f['/Header'].attrs['NumPart_ThisFile']
                for itype in range(6):
                    self.num_part_in_file[itype][ifile] = numpart_thisfile[itype]

        # These two datasets different in each file
        self.part_per_cell = [[None for ifile in range(self.numfiles)] for itype in range(6)]
        self.first_in_cell = [[None for ifile in range(self.numfiles)] for itype in range(6)]
        # don't force read here, read on demand instead

        self._collect_dataset_names()
        
        self.split_rank = -1
        self.split_size = -1

        return


    def __del__(self):
        self.close()
        return

    def select_region(self, xmin, xmax, ymin, ymax, zmin, zmax):
        """Select a region to read in"""
        if self.verbose:
            print('select_region() called')
        ixmin = int(np.floor(xmin / self.boxsize * self.ncell))
        ixmax = int(np.floor(xmax / self.boxsize * self.ncell))
        iymin = int(np.floor(ymin / self.boxsize * self.ncell))
        iymax = int(np.floor(ymax / self.boxsize * self.ncell))
        izmin = int(np.floor(zmin / self.boxsize * self.ncell))
        izmax = int(np.floor(zmax / self.boxsize * self.ncell))
        return self.select_grid_cells(ixmin, ixmax, iymin, iymax, izmin, izmax)


    def select_grid_cells(self, ixmin, ixmax, iymin, iymax, izmin, izmax):
        """Select hash grid cells to read in"""
        if self.verbose:
            print('select_grid_cells() called')
        n = 0
        for ixyz in product(
                range(ixmin, ixmax + 1),
                range(iymin, iymax + 1),
                range(izmin, izmax + 1)
        ):
            ixyz = np.array(ixyz)
            while (ixyz < 0).any():
                ixyz[ixyz < 0] += self.ncell
            while (ixyz >= self.ncell).any():
                ixyz[ixyz >= self.ncell] -= self.ncell
            self.hashmap[self._peano_hilbert_key(*tuple(ixyz))] = 1
            n += 1
        if self.verbose:
            print("  - Selected {:d} cells of {:d}".format(n, self.nhash))
        return

    def select_rotated_region(self, centre, xvec, yvec, zvec, length):
        """Select a non axis aligned region to read in"""
        diagonal = np.sqrt(3) * self.boxsize / self.ncell
        for ixyz in product(range(self.ncell), range(self.ncell), range(self.ncell)):
            ixyz = np.array(ixyz)
            cell_centre = self.boxsize / self.ncell * (ixyz + .5) - np.array(centre)
            while (cell_centre > .5 * self.boxsize).any():
                cell_centre[cell_centre > .5 * self.boxsize] -= self.boxsize
            while (cell_centre < .5 * self.boxsize).any():
                cell_centre[cell_centre < .5 * self.boxsize] += self.boxsize
            pos = np.array([
                np.sum(cell_centre * xvec),
                np.sum(cell_centre * yvec),
                np.sum(cell_centre * zvec)
            ])
            if (pos > -.5 * length - .5 * diagonal).all() \
               and (pos < .5 * length + .5 * diagonal).all():
                self.hashmap[self._peano_hilbert_key(*tuple(ixyz))] = 1
        return

    def set_sampling_rate(self, rate):
        """Set the sampling rate for subsequent reads"""
        self.sampling_rate = rate
        return

    def clear_selection(self):
        """Clear the current selection"""
        if self.verbose:
            print('clear_selection() called')
        self.hashmap = np.zeros(self.nhash)
        # it's ok to call split selection() again after clearing the selection
        self.split_size = -1
        self.split_rank = -1
        return

    def count_particles(self, itype):
        """Return the number of particles in the selected region"""
        return self.get_particle_locations(itype, _count=True)

    def get_particle_locations(self, itype, _count=False):
        """Return the locations of particles in the selected region"""
        nmax = 0
        file_index = []
        file_offset = []
        if itype < 0 or itype > 5:
            raise ValueError("Particle type index is out of range")
        if self.verbose:
            print("count_particles() called")
        if self.numpart_total[itype] == 0:
            return 0
        np.random.seed(1)
        n_to_read = 0
        end_key = -1
        for ifile in range(self.numfiles):
            if self.num_keys_in_file[itype][ifile] > 0:
                for key in range(self.first_key_in_file[itype][ifile],
                                 self.last_key_in_file[itype][ifile] + 1):
                    if key <= end_key:
                        continue
                    if self.hashmap[key]:
                        if self.part_per_cell[itype][ifile] is None:
                            self._load_hash_table(itype, ifile)
                        end_key = key
                        while end_key <= self.last_key_in_file[itype][ifile] and self.hashmap[end_key]:
                            end_key += 1
                        end_key -= 1

                        for this_key in range(key, end_key + 1):
                            if self.sampling_rate >= 1.0:
                                offset = this_key - self.first_key_in_file[itype][ifile]
                                n = self.part_per_cell[itype][ifile][offset]
                                n_to_read += n
                                if not _count:
                                    file_index.append(np.ones(n, dtype=np.int) * ifile)
                                    start = self.first_in_cell[itype][ifile][offset]
                                    end = start + n
                                    file_offset.append(np.arange(start, end, dtype=np.int))
                            else:
                                pass
        if _count:
            return n_to_read
        else:
            return np.concatenate(file_index), np.concatenate(file_offset)

    def read_dataset(self, itype, name):
        """Read a dataset and return it as a Numpy array"""
        return self.read_extra_dataset(itype, name, None)

    def read_extra_dataset(self, itype, name, basename):
        """Read a dataset from an auxiliary file and return it as a np.array"""
        return self._read_extra_dataset(itype, name, basename)

    def datasets(self, itype):
        """Return a list of datasets for the specified particle type"""
        if itype < 0 or itype > 5:
            raise ValueError("Particle type index is out of range")
        return self.dataset_names[itype]

    def split_selection(self, ThisTask, NTask):
        """Split the selected region(s) between processors"""
        self._split_selection(ThisTask, NTask)

    def close(self):
        """Close the snapshot and deallocate associated memory"""
        # should probably delete the arrays here?
        return

    def _collect_dataset_names(self):
        self.num_datasets = [0 for i in range(6)]
        self.dataset_names = [None for i in range(6)]
        #locate a file which actually has data for given particle type
        for itype in range(6):
            if self.numpart_total[itype] > 0:
                ifile = 0
                while self.num_keys_in_file[itype][ifile] == 0:
                    ifile += 1
                filename = '{:s}.{:d}.hdf5'.format(self.basename, ifile)
                # open this file and record dataset names
                with h5py.File(filename, 'r') as f:
                    try:
                        g = f['/PartType{:d}'.format(itype)]
                    except KeyError:
                        self.num_datasets[itype] = 0
                    else:
                        self.dataset_names[itype] = _get_dataset_list(g)
                        self.num_datasets[itype] = len(self.dataset_names[itype])
        return  

    def _read_extra_dataset(self, itype, name, basename):
        if (itype < 0) or (itype > 5):
            raise ValueError('Particle type itype is outside range 0-5!')
        if self.verbose:
            print('read_dataset() called')
        if self.numpart_total[itype] == 0:
            return
        name = "PartType{:d}/{:s}".format(itype, name)
        retval = []
        np.random.seed(1)
        n_to_read = 0
        end_key = -1
        for ifile in range(self.numfiles):
            if self.num_keys_in_file[itype][ifile] > 0:
                fname = '{:s}.{:d}.hdf5'.format(basename if basename else self.basename, ifile)
                with h5py.File(fname, 'r') as f:
                    if self.verbose:
                        print('  - Opened file {:d}'.format(ifile))
                    # boolean masks apparently may be slow in h5py? need to use slices?
                    s = np.zeros(self.num_part_in_file[itype][ifile], dtype=np.bool)
                    for key in range(self.first_key_in_file[itype][ifile],
                                     self.last_key_in_file[itype][ifile]):
                        if key <= end_key:
                            continue
                        if self.hashmap[key]:
                            # this key is in the current file, so make sure the hash
                            # table for the file is in memory
                            if self.part_per_cell[itype][ifile] is None:
                                self._load_hash_table(itype, ifile)
                            # check if there's a series of consecutive selected cells we can read
                            end_key = key
                            while end_key <= self.last_key_in_file[itype][ifile] and self.hashmap[end_key]:
                                end_key += 1
                            end_key -= 1
                            try:
                                d = f[name]
                            except KeyError:
                                raise KeyError('Unable to open dataset: {:s}'.format(name))
                            nread_file = 0
                            count = 0
                            for this_key in range(key, end_key + 1):
                                count += self.part_per_cell[itype][ifile][
                                    this_key - self.first_key_in_file[itype][ifile]]
                            start = int(self.first_in_cell[itype][ifile][
                                    key - self.first_key_in_file[itype][ifile]])
                            rank = d.ndim
                            nread_file += count
                            s[start:start + count] = True
                            key = end_key
                    if nread_file > 0:
                        if not rank in (1, 2):
                            raise RuntimeError('Can only read 1D or 2D datasets!')
                        if self.sampling_rate >= 1.0:
                            if rank == 2:
                                s = np.s_[s, :]
                            retval.append(f[name][s])
                        else:
                            pass
        return np.concatenate(retval)

    def _split_selection(self, ThisTask, NTask):
        if (ThisTask < 0) or (Ntask < 1) or (ThisTask >= Ntask):
            raise ValueError('Invalid paramters')
        # don't allow repeat splitting
        if (self.split_rank >= 0) or self.split_size >= 0:
            raise RuntimeError('Selection has already been split!')
        # count how many hash cels have been selected
        selected_keys = np.sum(self.hashmap)
        # decide how any keys to assign to this and previous processors
        Task = 0
        nkey_prev = 0
        nkey_this = 0
        for key in range(selected_keys):
            if Task < ThisTask:
                nkey_prev += 1
            if Task == ThisTask:
                nkey_this += 1
            Task = (Task + 1) % NTask
        # un-select keys outside range to read on this processor
        selected_keys = 0
        for key in range(self.nhash):
            if self.hashmap[key] != 0:
                selected_keys += 1
                if selected_keys <= nkey_prev:
                    self.hashmap[key] = 0
                if selected_keys > nkey_prev + nkey_this:
                    self.hashmap[key] = 0
        # store split parameters
        self.split_size = NTask
        self.split_rank = ThisTask
        return

    def _load_hash_table(self, itype, ifile):
        if self.verbose:
            print('  - Reading hash table for file {:d} particle type {:d}'.format(ifile, itype))
        if (self.numpart_total[itype] > 0) and (self.num_keys_in_file[itype][ifile] > 0):
            with h5py.File('{:s}.{:d}.hdf5'.format(self.basename, ifile), 'r') as f:
                self.part_per_cell[itype][ifile] = \
                    f['/HashTable/PartType{:d}/NumParticleInCell'.format(itype)][...]
                self.first_in_cell[itype][ifile] = np.r_[
                    0,
                    np.cumsum(self.part_per_cell[itype][ifile])[:-1]
                ]
        return

    def _peano_hilbert_key(self, x, y, z):
        mask = 1 << (self.hashbits - 1)
        key = 0
        rotation = 0
        sense = 1

        for i in range(self.hashbits):
            bitx = np.where(x & mask, 1, 0)
            bity = np.where(y & mask, 1, 0)
            bitz = np.where(z & mask, 1, 0)

            quad = quadrants[rotation][bitx][bity][bitz]

            key <<= 3
            key += quad if sense == 1 else 7 - quad

            rotx = rotx_table[quad]
            roty = roty_table[quad]
            sense *= sense_table[quad]

            while rotx > 0:
                rotation = rotxmap_table[rotation]
                rotx -= 1
            while roty > 0:
                rotation = rotymap_table[rotation]
                roty -= 1

            mask >>= 1

        return key


def _get_dataset_list(grp, prefix=""):
    all_objs = list()
    if prefix:
        grp[prefix].visit(all_objs.append)
    else:
        grp.visit(all_objs.append)
    all_dsets = ['/' + obj for obj in all_objs if isinstance(grp[obj], h5py.Dataset)]
    return all_dsets