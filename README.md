# coralscape

Custom scripts for measuring coral reef microhabitat structural complexity on 3D point cloud data (structure-from-motion photogrammetry) using open3D and numpy packages in python.

## Running analyses

1. Rotate point cloud and annotations according to the 'up vector'

```bash
$ rotate_plot.py [ply filname] [annotations filename] [json file]
```
The annotations filename is the coordinates and names of your samples (or points of interest) and the json file is the metadata from Viscore.

2. Get the environment points surrounding your points of interest

```bash
$ export_environment_points.py [ply filname] [annotations filename] [json file] [env radius]
```

3. Get measures

```bash
$ python get_measures.py [environment points ply filename] [annotations filename]
```

## HPC workflow

I prefer to use the software `parallel` for parallel computation.