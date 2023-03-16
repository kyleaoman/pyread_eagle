import h5py
import numpy as np
from functools import wraps


_random_seed = 1  # doesn't matter which, but ensure consistency

_quadrants = np.array(
    [
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
        [[[5, 2], [6, 1]], [[4, 3], [7, 0]]],
    ]
)

_rotxmap_table = np.array(
    [
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        0,
        1,
        2,
        3,
        17,
        18,
        19,
        16,
        23,
        20,
        21,
        22,
    ]
)
_rotymap_table = np.array(
    [
        1,
        2,
        3,
        0,
        16,
        17,
        18,
        19,
        11,
        8,
        9,
        10,
        22,
        23,
        20,
        21,
        14,
        15,
        12,
        13,
        4,
        5,
        6,
        7,
    ]
)

_rotx_table = np.array([3, 0, 0, 2, 2, 0, 0, 1])
_roty_table = np.array([0, 1, 1, 2, 2, 3, 3, 0])

_sense_table = np.array([-1, -1, -1, +1, +1, -1, -1, -1])


def _get_dataset_list(grp, prefix=""):
    all_dsets = list()
    if prefix:
        all_dsets = [
            "/{:s}{:s}".format(prefix, dset_name)
            for dset_name in _get_dataset_list(grp[prefix])
        ]
    else:
        for key in grp.keys():
            if isinstance(grp[key], h5py._hl.group.Group):
                all_dsets.extend(_get_dataset_list(grp, prefix=key))
            else:
                all_dsets.append("/{:s}".format(key))
    return all_dsets


class EagleSnapshotClosedException(Exception):
    pass


def check_open(method):
    @wraps(method)
    def _check_open(self, *method_args, **method_kwargs):
        if self.isclosed:
            raise EagleSnapshotClosedException
        else:
            return method(self, *method_args, **method_kwargs)

    return _check_open


class EagleSnapshot(object):
    """Class to represent an open Eagle snapshot"""

    def __init__(self, fname, verbose=False):
        """Open a new snapshot"""
        self.isclosed = False
        self.fname = fname
        self.verbose = verbose
        if self.verbose:
            print("__init__() called")
        self.sampling_rate = 1.0
        self.split_rank = -1
        self.split_size = -1
        with h5py.File(self.fname, "r") as f:
            if self.verbose:
                print("  - Opened file: {:s}".format(self.fname))
            self.boxsize = f["/Header"].attrs["BoxSize"]
            self.numfiles = f["/Header"].attrs["NumFilesPerSnapshot"]
            nptot = f["/Header"].attrs["NumPart_Total"]
            nptot_hw = (f["/Header"].attrs["NumPart_Total_HighWord"]).astype("int64")
            self.hashbits = f["/HashTable"].attrs["HashBits"]
            self.ncell = 1 << self.hashbits
            self.nhash = 1 << 3 * self.hashbits
            self.numpart_total = [nptot[i] + (nptot_hw[i] << 32) for i in range(6)]
            if (np.array(self.numpart_total) < 0).any():
                print(
                    "EagleSnapshot WARNING: " "negative total counts (overflow?).",
                    self.numpart_total,
                )
            if self.verbose:
                print("  - Read in file header")
            # initialize hashmap all false
            self.hashmap = np.zeros(self.nhash, dtype=bool)
            name_parts = self.fname.split(".")
            if name_parts[-1] != "hdf5" or not name_parts[-2].isdigit():
                raise RuntimeError("Don't understand snapshot file name!")
            else:
                self.basename = ".".join(name_parts[:-2])
                if self.verbose:
                    print("  - Base name is {:s}".format(self.basename))
            if self.verbose:
                for itype in range(6):
                    if self.numpart_total[itype] != 0:
                        print("  - Have particles of type {:d}".format(itype))
            # These three datasets repeated in each file
            self.first_key_in_file = [
                f["/HashTable/PartType{:d}/FirstKeyInFile".format(itype)][...]
                if self.numpart_total[itype] != 0
                else None
                for itype in range(6)
            ]
            self.last_key_in_file = [
                f["/HashTable/PartType{:d}/LastKeyInFile".format(itype)][...]
                if self.numpart_total[itype] != 0
                else None
                for itype in range(6)
            ]
            self.num_keys_in_file = [
                f["/HashTable/PartType{:d}/NumKeysInFile".format(itype)][...]
                if self.numpart_total[itype] != 0
                else None
                for itype in range(6)
            ]
        self.num_part_in_file = [
            [None for ifile in range(self.numfiles)] for itype in range(6)
        ]
        for ifile in range(self.numfiles):
            fname = "{:s}.{:d}.hdf5".format(self.basename, ifile)
            with h5py.File(fname, "r") as f:
                numpart_thisfile = f["/Header"].attrs["NumPart_ThisFile"]
                for itype in range(6):
                    self.num_part_in_file[itype][ifile] = numpart_thisfile[itype]

        # These two datasets different in each file
        self.part_per_cell = [
            [None for ifile in range(self.numfiles)] for itype in range(6)
        ]
        self.first_in_cell = [
            [None for ifile in range(self.numfiles)] for itype in range(6)
        ]
        # don't force read here, read on demand instead

        self._collect_dataset_names()

        return

    def __del__(self):
        self.close()
        return

    @check_open
    def select_region(self, xmin, xmax, ymin, ymax, zmin, zmax):
        """Select a region to read in"""
        if self.verbose:
            print("select_region() called")
        ixmin = int(np.floor(xmin / self.boxsize * self.ncell))
        ixmax = int(np.floor(xmax / self.boxsize * self.ncell))
        iymin = int(np.floor(ymin / self.boxsize * self.ncell))
        iymax = int(np.floor(ymax / self.boxsize * self.ncell))
        izmin = int(np.floor(zmin / self.boxsize * self.ncell))
        izmax = int(np.floor(zmax / self.boxsize * self.ncell))
        return self.select_grid_cells(ixmin, ixmax, iymin, iymax, izmin, izmax)

    @check_open
    def select_grid_cells(self, ixmin, ixmax, iymin, iymax, izmin, izmax):
        """Select hash grid cells to read in"""
        if self.verbose:
            print("select_grid_cells() called")
        ixyz = np.vstack(
            [
                _i.ravel()
                for _i in np.mgrid[
                    ixmin : ixmax + 1, iymin : iymax + 1, izmin : izmax + 1
                ]
            ]
        )
        while (ixyz < 0).any():
            ixyz[ixyz < 0] += self.ncell
        while (ixyz >= self.ncell).any():
            ixyz[ixyz >= self.ncell] -= self.ncell
        self.hashmap[self._peano_hilbert_key(*tuple(ixyz))] = True
        if self.verbose:
            print(
                "  - Selected {:d} cells of {:d}".format(
                    np.sum(self.hashmap), self.nhash
                )
            )
        return

    @check_open
    def select_rotated_region(self, centre, xvec, yvec, zvec, length):
        """Select a non axis aligned region to read in"""
        centre = np.array(centre)
        xvec = np.array(xvec)
        yvec = np.array(yvec)
        zvec = np.array(zvec)
        length = np.array(length)
        if any([param.shape != (3,) for param in (centre, xvec, yvec, zvec, length)]):
            raise ValueError(
                "select_rotated_region parameters must all have " "shape (3, )"
            )
        if any(
            [
                np.dot(xvec, yvec) != 0.0,
                np.dot(xvec, zvec) != 0,
                np.linalg.norm(xvec) != 1.0,
                np.linalg.norm(yvec) != 1.0,
                np.linalg.norm(zvec) != 1.0,
            ]
        ):
            raise ValueError(
                "select_rotated_region parameters xvec, yvec &"
                " zvec must be mutually orthogonal unit vectors"
            )
        diagonal = np.sqrt(3) * self.boxsize / self.ncell
        ixyz = np.vstack(
            [
                _i.ravel()
                for _i in np.mgrid[0 : self.ncell, 0 : self.ncell, 0 : self.ncell]
            ]
        ).T
        cell_centre = self.boxsize / self.ncell * (ixyz + 0.5) - centre
        while (cell_centre > 0.5 * self.boxsize).any():
            cell_centre[cell_centre > 0.5 * self.boxsize] -= self.boxsize
        while (cell_centre < -0.5 * self.boxsize).any():
            cell_centre[cell_centre < -0.5 * self.boxsize] += self.boxsize
        pos = np.vstack(
            (
                np.sum(cell_centre * xvec, axis=1),
                np.sum(cell_centre * yvec, axis=1),
                np.sum(cell_centre * zvec, axis=1),
            )
        ).T
        accept_cells = np.logical_and(
            (pos > -0.5 * length - 0.5 * diagonal).all(axis=1),
            (pos < 0.5 * length + 0.5 * diagonal).all(axis=1),
        )
        self.hashmap[self._peano_hilbert_key(*tuple(ixyz[accept_cells].T))] = True
        return

    @check_open
    def set_sampling_rate(self, rate):
        """Set the sampling rate for subsequent reads"""
        self.sampling_rate = rate
        return

    @check_open
    def clear_selection(self):
        """Clear the current selection"""
        if self.verbose:
            print("clear_selection() called")
        self.hashmap = np.zeros(self.nhash, dtype=bool)
        # it's ok to call split selection() again after clearing the selection
        self.split_size = -1
        self.split_rank = -1
        return

    @check_open
    def count_particles(self, itype):
        """Return the number of particles in the selected region"""
        return self.get_particle_locations(itype, _count=True)

    @check_open
    def get_particle_locations(self, itype, _count=False):
        """Return the locations of particles in the selected region"""
        # revise this function to be similar to read_extra_dataset to speed up
        file_offsets = list()
        if itype < 0 or itype > 5:
            raise ValueError("Particle type index is out of range")
        if self.verbose:
            print("count_particles() called")
        if self.numpart_total[itype] == 0:
            if _count:
                return 0
            else:
                return np.empty((0,)), np.empty((0,))
        counts = list()
        for ifile in range(self.numfiles):
            if self.num_keys_in_file[itype][ifile] == 0:
                counts.append(0)
                continue
            cell_mask = self.hashmap[
                self.first_key_in_file[itype][ifile] : self.last_key_in_file[itype][
                    ifile
                ]
                + 1
            ]
            if not cell_mask.any():
                counts.append(0)
                continue
            if self.part_per_cell[itype][ifile] is None:
                self._load_hash_table(itype, ifile)
            lengths = self.part_per_cell[itype][ifile][cell_mask]
            if self.sampling_rate < 1.0:
                counts.append(int(np.floor(np.sum(lengths) * self.sampling_rate)))
            else:
                counts.append(np.sum(lengths))
            if self.verbose:
                print("  ", counts[-1])
            if not _count:
                starts = self.first_in_cell[itype][ifile][cell_mask]
                ends = starts + lengths
                fos = np.array(
                    np.repeat(ends - np.cumsum(lengths), lengths)
                    + np.arange(np.sum(lengths)),
                    dtype=int,
                )
                if self.sampling_rate < 1.0:
                    sub = np.arange(fos.shape[0])
                    np.random.seed(_random_seed)
                    np.random.shuffle(sub)
                    sub = sub[: int(np.floor(counts[-1]))]
                    fos = fos[sub]
                file_offsets.append(fos)
        if _count:
            # count
            return np.sum(counts)
        else:
            # file index, file_offset
            return np.repeat(np.arange(self.numfiles), counts), np.concatenate(
                file_offsets
            )

    @check_open
    def read_dataset(self, itype, name):
        """Read a dataset and return it as a Numpy array"""
        return self.read_extra_dataset(itype, name, None)

    @check_open
    def read_extra_dataset(self, itype, name, basename):
        """Read a dataset from an auxiliary file and return it as a np.array"""
        if (itype < 0) or (itype > 5):
            raise ValueError("Particle type itype is outside range 0-5!")
        if self.verbose:
            print("read_dataset() called")
        if self.numpart_total[itype] == 0:
            return
        name = "PartType{:d}/{:s}".format(itype, name)
        all_retval = list()
        for ifile in range(self.numfiles):
            if self.num_keys_in_file[itype][ifile] == 0:
                continue
            fname = "{:s}.{:d}.hdf5".format(
                basename if basename else self.basename, ifile
            )
            with h5py.File(fname, "r") as f:
                if self.verbose:
                    print("  - Opened file {:d}".format(ifile))
                try:
                    d = f[name]
                except KeyError:
                    raise KeyError("Unable to open dataset: {:s}".format(name))
                ndim = d.ndim
                ncol = 1 if ndim == 1 else d.shape[1]
                if ndim not in (1, 2):
                    raise RuntimeError("Can only read 1D or 2D datasets!")
                cell_mask = self.hashmap[
                    self.first_key_in_file[itype][ifile] : self.last_key_in_file[itype][
                        ifile
                    ]
                    + 1
                ]
                if not cell_mask.any():
                    continue
                if self.part_per_cell[itype][ifile] is None:
                    self._load_hash_table(itype, ifile)
                lowers = np.argwhere(np.diff(cell_mask.astype(int)) > 0) + 1
                uppers = np.argwhere(np.diff(cell_mask.astype(int)) < 0) + 1
                if cell_mask[0]:
                    lowers = np.r_[[[0]], lowers]
                if cell_mask[-1]:
                    uppers = np.r_[uppers, [[len(cell_mask)]]]
                cell_intervals = np.hstack((lowers, uppers))
                counts = list()
                starts = list()
                retval = list()
                for interval in cell_intervals:
                    counts.append(
                        int(
                            np.sum(
                                self.part_per_cell[itype][ifile][
                                    interval[0] : interval[1]
                                ]
                            )
                        )
                    )
                    starts.append(int(self.first_in_cell[itype][ifile][interval[0]]))
                if np.sum(counts) > 0:
                    for start, count in zip(starts, counts):
                        # the reading here is a current bottleneck
                        dat = f[name][start : start + count]
                        if self.sampling_rate >= 1.0:
                            all_retval.append(dat)
                        else:
                            retval.append(dat)
                if self.sampling_rate < 1.0:
                    # cut file-by-file to conserve memory
                    retval = np.concatenate(retval)
                    sub = np.arange(retval.shape[0])
                    np.random.seed(_random_seed)
                    np.random.shuffle(sub)
                    sub = sub[: int(np.floor(retval.shape[0] * self.sampling_rate))]
                    all_retval.append(retval[sub])
        if len(all_retval):  # check that data not empty
            return np.concatenate(all_retval)
        elif ndim == 2:  # guaranteed to be 1 or 2 by check above
            return np.empty((0, ncol))  # empty is ok with axis length 0
        else:  # ndim == 1
            return np.empty((0,))

    @check_open
    def datasets(self, itype):
        """Return a list of datasets for the specified particle type"""
        if itype < 0 or itype > 5:
            raise ValueError("Particle type index is out of range")
        return self.dataset_names[itype]

    @check_open
    def split_selection(self, ThisTask, NTask, new=False):
        """Split the selected region(s) between processors"""
        if (ThisTask < 0) or (NTask < 1) or (ThisTask >= NTask):
            raise ValueError("Invalid paramters")
        # don't allow repeat splitting
        if (self.split_rank >= 0) or self.split_size >= 0:
            raise RuntimeError("Selection has already been split!")
        # count how many hash cels have been selected
        selected_keys = np.sum(self.hashmap)
        # decide how any keys to assign to this and previous processors
        nkey_split = selected_keys // NTask + np.asarray(
            selected_keys % NTask > np.arange(NTask), dtype=int
        )
        nkey_this = nkey_split[ThisTask]
        nkey_prev = np.sum(nkey_split[:ThisTask])
        # un-select keys outside range to read on this processor
        self.hashmap[self.hashmap] = np.r_[
            np.zeros(nkey_prev, dtype=bool),
            np.ones(nkey_this, dtype=bool),
            np.zeros(selected_keys - nkey_prev - nkey_this, dtype=bool),
        ]
        # store split parameters
        self.split_size = NTask
        self.split_rank = ThisTask
        return

    @check_open
    def _load_hash_table(self, itype, ifile):
        if self.verbose:
            print(
                "  - Reading hash table for file "
                "{:d} particle type {:d}".format(ifile, itype)
            )
        if (self.numpart_total[itype] != 0) and (
            self.num_keys_in_file[itype][ifile] > 0
        ):
            fname = "{:s}.{:d}.hdf5".format(self.basename, ifile)
            with h5py.File(fname, "r") as f:
                hpath = "/HashTable/PartType{:d}/" "NumParticleInCell".format(itype)
                self.part_per_cell[itype][ifile] = f[hpath][...]
                self.first_in_cell[itype][ifile] = np.r_[
                    0, np.cumsum(self.part_per_cell[itype][ifile])[:-1]
                ]
        return

    @check_open
    def close(self):
        """Close the snapshot and deallocate associated memory"""
        del self.fname
        del self.verbose
        del self.sampling_rate
        del self.boxsize
        del self.numfiles
        del self.hashbits
        del self.ncell
        del self.nhash
        del self.numpart_total
        del self.hashmap
        del self.basename
        del self.first_key_in_file
        del self.last_key_in_file
        del self.num_keys_in_file
        del self.num_part_in_file
        del self.part_per_cell
        del self.first_in_cell
        del self.num_datasets
        del self.dataset_names
        del self.split_rank
        del self.split_size
        self.isclosed = True

        return

    @check_open
    def _collect_dataset_names(self):
        self.num_datasets = [0 for i in range(6)]
        self.dataset_names = [None for i in range(6)]
        # locate a file which actually has data for given particle type
        for itype in range(6):
            if self.numpart_total[itype] == 0:
                continue
            ifile = 0
            while self.num_part_in_file[itype][ifile] == 0:
                ifile += 1
            filename = "{:s}.{:d}.hdf5".format(self.basename, ifile)
            # open this file and record dataset names
            with h5py.File(filename, "r") as f:
                try:
                    g = f["/PartType{:d}".format(itype)]
                except KeyError:
                    self.num_datasets[itype] = 0
                else:
                    self.dataset_names[itype] = _get_dataset_list(g)
                    self.num_datasets[itype] = len(self.dataset_names[itype])
        return

    @check_open
    def _peano_hilbert_key(self, x, y, z):
        x = np.atleast_1d(x)
        y = np.atleast_1d(y)
        z = np.atleast_1d(z)
        mask = 1 << (self.hashbits - 1)  # leave scalar
        key = np.zeros(x.shape, dtype=int)
        rotation = np.zeros(x.shape, dtype=int)
        sense = np.ones(x.shape, dtype=int)

        for i in range(self.hashbits):
            bitx = np.where(x & mask, 1, 0)
            bity = np.where(y & mask, 1, 0)
            bitz = np.where(z & mask, 1, 0)

            quad = _quadrants[rotation, bitx, bity, bitz]

            key = key << 3
            key += np.where(sense > 0, quad, 7 - quad)

            rotx = _rotx_table[quad]
            roty = _roty_table[quad]
            sense = sense * _sense_table[quad]
            while (rotx > 0).any():
                rotation = np.where(rotx > 0, _rotxmap_table[rotation], rotation)
                rotx[rotx > 0] -= 1
            while (roty > 0).any():
                rotation = np.where(roty > 0, _rotymap_table[rotation], rotation)
                roty[roty > 0] -= 1

            mask >>= 1

        return key
