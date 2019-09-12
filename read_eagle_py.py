import h5py
import numpy as np
from itertools import product


class EagleSnapshotClosedException(Exception):
    pass


class EagleSnapshot(object):
    """Class to represent an open Eagle snapshot"""

    def __init__(self, fname, verbose=False):
        self.verbose = verbose
        self.F = None
        """Open a new snapshot"""
        (self.snap, n0, n1, n2, n3, n4, n5, self.boxsize, self.numfiles,
         self.hashbits) = self._open_snapshot(fname)
        # Store particle number
        self.numpart_total = (n0, n1, n2, n3, n4, n5)
        # Get names of datasets
        self._dataset = []
        for itype in range(6):
            self._dataset.append([])
            for iset in range(self._get_dataset_count(itype)):
                self._dataset[-1].append(
                    self._get_dataset_name(itype, iset)
                )

    def __del__(self):
        """Close the snapshot and deallocate, if we haven't already"""
        if self.F is not None:
            self._close_snapshot()

    def select_region(self, xmin, xmax, ymin, ymax, zmin, zmax):
        """Select a region to read in"""
        if not self.F:
            raise EagleSnapshotClosedException(
                "Cannot select region in closed snapshot!"
            )
        self._select_region(self.snap, xmin, xmax, ymin, ymax, zmin, zmax)

    def select_rotated_region(self, centre, xvec, yvec, zvec, length):
        """Select a non axis aligned region to read in"""
        if not self.F:
            raise EagleSnapshotClosedException(
                "Cannot select region in closed snapshot!"
            )
        cx, cy, cz = [float(c) for c in centre]
        xx, xy, xz = [float(c) for c in xvec]
        yx, yy, yz = [float(c) for c in yvec]
        zx, zy, zz = [float(c) for c in zvec]
        lx, ly, lz = [float(c) for c in length]
        self._select_rotated_region(
            cx, cy, cz,
            xx, xy, xz,
            yx, yy, yz,
            zx, zy, zz,
            lx, ly, lz
        )

    def select_grid_cells(self, ixmin, ixmax, iymin, iymax, izmin, izmax):
        """Select hash grid cells to read in"""
        if not self.F:
            raise EagleSnapshotClosedException(
                "Cannot select region in closed snapshot!"
            )
        self._select_grid_cells(
            ixmin, ixmax,
            iymin, iymax,
            izmin, izmax
        )

    def set_sampling_rate(self, rate):
        """Set the sampling rate for subsequent reads"""
        if not self.F:
            raise EagleSnapshotClosedException(
                "Cannot set sample rate in closed snapshot!"
            )
        self._set_sampling_rate(self.snap, rate)

    def clear_selection(self):
        """Clear the current selection"""
        if not self.F:
            raise EagleSnapshotClosedException(
                "Cannot clear selection in closed snapshot!"
            )
        self._clear_selection()

    def count_particles(self, itype):
        """Return the number of particles in the selected region"""
        if not self.F:
            raise EagleSnapshotClosedException(
                "Cannot count particles in closed snapshot!"
            )
        return self._count_particles(itype)

    def get_particle_locations(self, itype):
        """Return the locations of particles in the selected region"""
        if not self.F:
            raise EagleSnapshotClosedException(
                "Cannot count particles in closed snapshot!"
            )
        return self._get_particle_locations(itype)

    def read_dataset(self, itype, name):
        """Read a dataset and return it as a Numpy array"""
        if not self.F:
            raise EagleSnapshotClosedException(
                "Cannot read dataset from closed snapshot!"
            )
        return self._read_dataset(itype, name)

    def read_extra_dataset(self, itype, name, basename):
        """Read a dataset from an auxiliary file and return it as a np.array"""
        if not self.F:
            raise EagleSnapshotClosedException(
                "Cannot read dataset from closed snapshot!"
            )
        return self._read_extra_dataset(itype, name, basename)

    def datasets(self, itype):
        """Return a list of datasets for the specified particle type"""
        if itype < 0 or itype > 5:
            raise ValueError("Particle type index is out of range")
        return self._dataset[itype]

    def split_selection(self, ThisTask, NTask):
        """Split the selected region(s) between processors"""
        if not self.F:
            raise EagleSnapshotClosedException(
                "Cannot split selection in closed snapshot!"
            )
        self._split_selection(ThisTask, NTask)

    def close(self):
        """Close the snapshot and deallocate associated memory"""
        if not self.F:
            raise EagleSnapshotClosedException("Snapshot is already closed!")
        self._close_snapshot()
        self.F = None

    def _open_snapshot(self, fname):
        if self.verbose:
            print("open_snapshot() called")
        self.hashmap = None
        self.first_key_in_file = [None for i in range(6)]
        self.last_key_in_file = [None for i in range(6)]
        self.num_keys_in_file = [None for i in range(6)]
        self.part_per_cell = [None for i in range(6)]
        self.first_in_cell = [None for i in range(6)]
        self.dataset_name = [None for i in range(6)]
        try:
            self.F = h5py.File(fname, 'r')
        except OSError:
            raise IOError("Unable to open file: {:s}".format(fname))
        else:
            if self.verbose:
                print("  - Opened file: {:s}".format(fname))
        try:
            self.boxsize = self.F['/Header'].attrs['BoxSize']
            self.numfiles = self.F['/Header'].attrs['NumFilesPerSnapshot']
            nptot = self.F['/Header'].attrs['NumPart_Total']
            nptot_hw = self.F['/Header'].attrs['NumPart_Total_HighWord']
            self.hashbits = self.F['/HashTable'].attrs['HashBits']
        except KeyError:
            raise KeyError("Unable to read hash table information from file: "
                           "{:s}".format(fname))
        self.ncell = 1 << self.hashbits
        self.nhash = 1 << 3 * self.hashbits
        self.numpart_total = [nptot[i] + nptot_hw[i] << 32 for i in range(6)]
        if self.verbose:
            print("  - Read in file header")
        self.hashmap = np.zeros(self.nhash)
        name_parts = fname.split('.')
        if name_parts[-1] != 'hdf5' or not name_parts[-2].isdigit():
            raise RuntimeError("Don't understand snapshot file name!")
        else:
            self.basename = ".".join(name_parts[:-2])
            if self.verbose:
                print("  - Base name is {:s}".format(self.basename))
        for itype in range(6):
            if self.numpart_total[itype] > 0:
                if self.verbose:
                    print("  - Have particles of type {:d}".format(itype))
                self.first_key_in_file[itype] = self.F[
                    '/HashTable/PartType{:d}/FirstKeyInFile'.format(itype)
                ]
                self.last_key_in_file[itype] = self.F[
                    '/HashTable/PartType{:d}/LastKeyInFile'.format(itype)
                ]
                self.num_keys_in_file[itype] = self.F[
                    '/HashTable/PartType{:d}/NumKeysInFile'.format(itype)
                ]
        self.F.close()  # this will mess up open checks, I guess
        # need to distinguish between open *file* and open *snapshot*
        self.part_per_cell = [[None for j in range(self.numfiles)]
                              for i in range(6)]
        self.first_in_cell = [[None for j in range(self.numfiles)]
                              for i in range(6)]

        self.num_datasets = [0 for i in range(6)]
        self.dataset_name = [None for i in range(6)]

        for itype in range(6):
            if self.numpart_total[itype] > 0:
                ifile = 0
                while self.num_keys_in_file[itype][ifile] == 0:
                    ifile += 1
                filename = '{:s}.{:d}.hdf5'.format(self.basename, ifile)
                with h5py.File(filename, 'r') as f:
                    try:
                        g = f['/PartType{:d}'.format(itype)]
                    except KeyError:
                        self.num_datasets[itype] = 0
                    else:
                        self._get_dataset_list(g, itype, "")
                        self.dataset_name[itype] = [None for i in range(self.num_datasets[itype])]
                        self.num_datasets[itype] = 0
                        self._get_dataset_list(g, itype, "")
                        # two calls to _get_dataset_list here, I think the first counts and
                        # the second actually stores the names, can condense this in python

        self.split_rank = -1
        self.split_size = -1

        return

    def _get_dataset_list(self, grp, itype, prefix):
        all_objs = list()
        if prefix:
            grp[prefix].visit(all_objs.append)
        else:
            grp.visit(all_objs.append)
        all_dsets = [obj for obj in all_objs if isinstance(obj, h5py.Dataset)]
        itype_dsets = [ds for ds in all_dsets if 'PartType{:d}'.format(itype) in ds]
        self.num_datasets[itype] = len(itype_dsets)
        self.dataset_name[itype] = itype_dsets
        return

    def _get_dataset_count(self, itype):
        if (itype > 5) or (itype < 0):
            raise ValueError("Particle type itype is outside range 0-5!")
        if self.numpart_total[itype] == 0:
            return 0
        return self.num_datasets[itype]

    def _get_dataset_name(self, itype, iset):
        if (itype < 0) or (itype > 5):
            raise ValueError("Particle type itype is outside range 0-5!")
        if (iset < 0) or iset >= self.num_datasets[itype]:
            raise ValueError("Dataset index is out of range!")
        # this function turns out to be pointless in python, refactor aggressively
        # to just get the names and store them during init (or maybe not even that)

    def _select_region(self, xmin, xmax, ymin, ymax, zmin, zmax):
        ixmin = int(np.floor(xmin / self.boxsize * self.ncell))
        ixmax = int(np.floor(xmax / self.boxsize * self.ncell))
        iymin = int(np.floor(ymin / self.boxsize * self.ncell))
        iymax = int(np.floor(ymax / self.boxsize * self.ncell))
        izmin = int(np.floor(zmin / self.boxsize * self.ncell))
        izmax = int(np.floor(zmax / self.boxsize * self.ncell))
        if self.verbose:
            print('select_region() called')
        n = 0
        for ix in range(ixmin, ixmax + 1):
            iix = ix
            while iix < 0:
                iix += self.ncell
            while iix >= self.ncell:
                iix -= self.ncell
            for iy in range(iymin, iymax + 1):
                iiy = iy
                while iiy < 0:
                    iiy += self.ncell
                while iiy >= self.ncell:
                    iiy -= self.ncell
                for iz in range(izmin, izmax + 1):
                    iiz = iz
                    while iiz < 0:
                        iiz += self.ncell
                    while iiz >= self.ncell:
                        iiz -= self.ncell
                    self.hashmap[peano_hilbert_key(iix, iiy, iiz, self.hashbits)] = 1
                    n += 1
        if self.verbose:
            print("  - Selected {:d} cells of {:d}".format(n, self.nhash))
        return

    def _close_snapshot(self):
        # revise once init is finalized
        pass

    def _select_rotated_region(self, cx, cy, cz, xx, xy, xz, yx, yy, yz,
                               zx, zy, zz, lx, ly, lz):
        diagonal = np.sqrt(3) * self.boxsize / self.ncell
        for ix, iy, iz in product(range(self.ncell), range(self.ncell), range(self.ncell)):
            cell_centre = [
                self.boxsize / self.ncell * (ix + .5) - cx,
                self.boxsize / self.ncell * (iy + .5) - cy,
                self.boxsize / self.ncell * (iz + .5) - cz
            ]
            for i in range(3):
                while cell_centre[i] > .5 * self.boxsize:
                    cell_centre[i] -= self.boxsize
                while cell_centre[i] < .5 * self.boxsize:
                    cell_centre[i] += self.boxsize
            pos = [
                cell_centre[0] * xx + cell_centre[1] * xy + cell_centre[2] * xz,
                cell_centre[0] * yx + cell_centre[1] * yy + cell_centre[2] * yz,
                cell_centre[0] * zx + cell_centre[1] * zy + cell_centre[2] * zz
            ]
            if (pos[0] > -.5 * lx - .5 * diagonal) and (pos[0] < .5 * lx + .5 * diagonal) \
               and (pos[1] > -.5 * ly - .5 * diagonal) and (pos[1] < .5 * ly + .5 * diagonal) \
               and (pos[2] > -.5 * lz - .5 * diagonal) and (pos[2] < .5 * lz + .5 * diagonal):
                self.hashmap[peano_hilbert_key(ix, iy, iz, self.hashbits)] = 1
        return
            

    def _select_grid_cells(self, ixmin, ixmax, iymin, iymax, izmin, izmax):
        # first arg was self.snap, no need to pass
        pass

    def _set_sampling_rate(self, rate):
        # first arg was self.snap, no need to pass
        pass

    def _clear_selection():
        # first arg was self.snap, no need to pass
        pass

    def _count_particles(self, itype):
        # first arg was self.snap, no need to pass
        pass

    def _get_particle_locations(self, itype):
        # first arg was self.snap, no need to pass
        pass

    def _read_dataset(self, itype, name):
        # first arg was self.snap, no need to pass
        pass

    def _read_extra_dataset(self, itype, name, basename):
        # first arg was self.snap, no need to pass
        pass

    def _split_selection(self, ThisTask, NTask):
        # first arg was self.snap, no need to pass
        pass

    def _close(self):
        # first arg was self.snap, no need to pass
        pass
