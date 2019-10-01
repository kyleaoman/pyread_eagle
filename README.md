## Read routines for Peano-Hilbert key sorted Eagle snapshots
Pure-python port of J. Helly's [read_eagle](https://github.com/jchelly/read_eagle).

### Installation
 - Download via web UI, or 'git clone https://github.com/kyleaoman/pyread_eagle.git'
 - Install dependencies if necessary (see 'setup.py'), some may be found in other repositories by kyleaoman.
 - Global install (Linux): 
   - cd to directory with 'setup.py'
   - run 'sudo pip install -e .' (-e installs via symlink, so pulling repository will do a 'live' update of the installation)
 - User install (Linux):
   - cd to directory with 'setup.py'
   - ensure '~/lib/python3.7/site-packages' or similar is on your PYTHONPATH (e.g. 'echo $PYTHONPATH'), if not, add it (perhaps in .bash_profile or similar)
   - run 'pip install --prefix ~ -e .' (-e installs via symlink, so pulling repository will do a 'live' update of the installation)
 - cd to a directory outside the module and launch python; you should be able to do 'from pyread_eagle import *'

### Documentation

The documentation below is adapted from the original [read_eagle readme file](https://github.com/jchelly/read_eagle/blob/master/README).

These routines can be used if Gadget was run with the `-DEAGLE_SORT_OUTPUT`
flag enabled. They provide a quick way to read in spatial regions
without having to read all the data or even open all of the files.

It works by splitting the simulation volume into a grid and ensuring
that particles in the same grid cell are stored consecutively in the 
snapshot files. Some extra datasets are added to the snapshot which
specify which grid cells are stored in which files and the location
of the start of each cell in the file.

The procedure to read in data is as follows:

- Open the snapshot by specifying the name of one file
- Flag the grid cells to be read in by calling `select_region` one or
  more times.
- Call `read_dataset` once for each quantity to be read, specifying
  which particle type to read and the hdf5 dataset name.

You can then either close/delete the snapshot object or call `clear_selection` to
read in a different region.

Any of the datasets in the `PartTypeX` groups can be read in. The code
makes no assumptions about what datasets are present in the snapshot.

Unlike the original `read_eagle`, this version is written directly in python
and therefore does not depend on having a properly configured C compiler
or hdf5 libraries, all that is needed is python, numpy and h5py.

Python2 is not supported. I *think* `pyread_eagle` will run in python2.7,
but since python2 is deprecated as of Jan 2020, I have not and have no intention
to test this.

### Simple example usage

```
from pyread_eagle import EagleSnapshot

snap = EagleSnapshot("./snap_020.0.hdf5")
snap.select_region(4.0, 5.0, 2.0, 3.0, 3.0, 4.0)
xyz = snap.read_dataset(0, "Coordinates")
ids = snap.read_dataset(0, "ParticleIDs")
del snap
```

### Speed

Currently pyread_eagle is around 20x slower than read_eagle. I think this can
still be improved by further refactoring to vectorize some loops and especially
making better use of the h5py API. Currently a boolean mask is used to apply
the particle selection to the hdf5 dataset, which in h5py is a known-issue
slowdown. For now, though, the hit in speed is the price to pay for ease of
installation.

### Error handling

If a routine fails it will raise an exception. Because of differences in 
the structure of the code relative to `read_eagle`, not all the same sanity
checks are made, not all the same errors can arise, and therefore not all
the same error messages may appear. I have also chosen to be more lax and
not try to catch every edge case by hand, rather letting failures happen
and be reported organically by the underlying python libraries. Please
report any errors whose messages are obscure, and of course any bugs.


### Description of the routines and parameters

#### Open snapshot

`snap = pyread_eagle.EagleSnapshot(fname, verbose=False)`

This opens the snapshot which contains the specified file. Closing
explicitly is not necessary: in python the snapshot will be automatically 
deallocated once there are no more references to it.

Parameters
  - `fname`: name of any one file in the snapshot
  - `verbose`: if `True`, print extra messages during function calls (default: False)

Return value
  An instance of the `EagleSnapshot` class

#### Close snapshot

`del snap`
`snap.close()`

Deallocates memory associated with the snap object.
Not strictly necessary - will happen automatically when the snap variable goes
out of scope. Attempting to call methods of a closed snapshot raises an
`EagleSnapshotClosedException`.

#### Region selection

`snap.select_region(xmin, xmax, ymin, ymax, zmin, zmax)`

All grid cells overlapping the specified region are flagged to be read
in by subsequent `read_dataset calls`. You can call `select_region` multiple
times to make oddly shaped or disjoint selections.

If selected regions overlap or the same region is selected multiple times
particles in these regions will still only be read in once.

Parameters
  - `xmin`: the minimum x coordinate of the region to read
  - `xmax`: the maximum x coordinate of the region to read
  - `ymin`: the minimum y coordinate of the region to read
  - `ymax`: the maximum y coordinate of the region to read
  - `zmin`: the minimum z coordinate of the region to read
  - `zmax`: the maximum z coordinate of the region to read

`snap.select_grid_cells(ixmin, ixmax, iymin, iymax, izmin, izmax)`

All grid cells in the specified range of grid coordinates are flagged to
be read in by subsequent `read_dataset` calls. You can call `select_grid_cells`
multiple times to make oddly shaped or disjoint selections.

The coordinates ixmin, ixmax etc are integer coordinates in the hash grid,
starting from zero. The maximum coordinate is `(2**hashbits)-1`. The value of
hashbits is stored in `snap.hashbits`.

If selected regions overlap or the same region is selected multiple times
particles in these regions will still only be read in once.

Parameters
  - `ixmin`: the minimum x coordinate of the region to read
  - `ixmax`: the maximum x coordinate of the region to read
  - `iymin`: the minimum y coordinate of the region to read
  - `iymax`: the maximum y coordinate of the region to read
  - `izmin`: the minimum z coordinate of the region to read
  - `izmax`: the maximum z coordinate of the region to read
  
`snap.select_rotated_region(centre, xvec, yvec, zvec, length)`

???

Parameters
  - `centre`: ?
  - `xvec`: ?
  - `yvec`: ?
  - `zvec`: ?
  - `length`: ?

#### Count particles

`n = snap.count_particles(itype)`

This returns the number of particles of the specified type which will be
read by the next `read_dataset` call. Note that only whole grid cells can
be read so some particles outside the selected region may be read in.
These are included in the count.

In python this routine is not usually needed, but it is included for
compatibility with codes formerly using the read_eagle python wrapper.

Parameters
  - `itype`: which particle type to count (integer, `0-5`)

Return value
  The number of particles to be read in

#### Get particle locations

`file_index, file_offset = snap.get_particle_locations(itype)`

This returns two arrays which each have one element for each selected
particle. `file_index` contains the index of the file each particle is in.
`file_offset` contains the position in the file, numbering from zero.

Parameters
  - `itype`: which particle type to count (integer, `0-5`)

Return value
  `file_index`  - integer array with index of the file containing each particle
  `file_offset` - integer array with position of each particle in its file


#### Read datasets

`data = snap.read_dataset(itype, name)`

This reads in the specified dataset for all particles of type itype in the
selected region(s). Use repeated calls to read in multiple datasets.

The type of array you get back reflects the type of the dataset in the file.

Parameters
  - `itype`: which particle type to read (integer, `0-5`)
  - `name`: the HDF5 name of the dataset, relative to the `PartTypeX` group

Return value
  The contents of the dataset for the selected particles
  
`data = snap.read_extra_dataset(itype, name, basename)`

This is the same as `read_dataset`, except that it reads from a set of
"auxiliary files" comprised of the same number of parts as the snapshot and 
which contain additional datasets with the same sorting, distribution across
files, etc. as the datasets in the main snapshot files.

Parameters
  - `itype`: which particle type to read (integer, `0-5`)
  - `name`: the HDF5 name of the dataset, relative to the `PartTypeX` group
  - `basename`: name (including path, if necessary) of the auxiliary files, omitting the `.X.hdf5` portion

#### Clear selection

`snap.clear_selection()`

This clears the flags which specify which grid cells should be
read in on the next read_dataset() call. If you've already read
in a region and you want to read a different region you should call 
clear_selection() before calling select_region() again.


#### Split selection

`snap.split_selection(ThisTask, NTask)`

This is for use in MPI programs to allow parallel I/O.
When running under MPI all read_eagle routines are collective -
they must be called on all processes.

The procedure to read a region in parallel is as follows:

1. All processes open the snapshot
2. All processes select the SAME region
3. Call `split_selection`. This causes each processor to select
   a subset of the originally selected particles.
4. Call `read_dataset()` on all processors to read the required 
   quantities

This results in all of the particles in the specified region 
being read in exactly once with the data spread across the MPI 
processes.

I have not made any attempt to use pyread_eagle with e.g. mpi4py
or other MPI-ish implementations in python. If you succeed in this,
please get in touch!

Parameters
  - `ThisTask`: rank of this process in an MPI program
  - `NTask`: number of processes in an MPI program

#### Set sampling rate

`snap.set_sampling_rate(rate)`

If rate is set < 1 before reading datasets, a random subsample
of particles will be read, with the fraction of particles read
set by the rate.

Parameters
  - `rate`: fraction of particles to be read (rate >= 1 for all particles)
  
#### List datasets

`dsets = snap.datasets(itype)`

Returns a list of the available datasets for a given particle type.
This can be iterated to read all datasets, for instance.

Parameters
  - `itype`: which particle type to read (integer, `0-5`)
  
 
