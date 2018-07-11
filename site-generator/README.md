## Requirements

The following python packages must be installed to run `data_preprocessor.py`:

- `h5py` - for data storing;
- `biopandas` - I use this to read mol2 files.
- `numpy`, `scipy` - for working with arrays
- `tqdm` - to show fancy status bar for long-running tasks (mostly for data preprocessing)
- `pyvdwsurface` - this is [my fork](https://github.com/latticetower/pyvdwsurface) 
of small, but good package used to get points on protein's Van der Waals surface.
