## BackTrackBB
Multi-band array detection and location of seismic sources

(c) 2015-2018  Natalia Poiata <poiata@ipgp.fr>, Claudio Satriano <satriano@ipgp.fr>;

(c) 2013-2014  Natalia Poiata <poiata@ipgp.fr>, Claudio Satriano <satriano@ipgp.fr>, Pierre Romanet <romanet@ipgp.fr>

BackTrackBB is a program for detection and space-time location of seismic sources
based on multi-scale, frequency-selective statistical coherence of the wave field
recorded by dense large-scale seismic networks and local antennas.
The method is designed to enhance coherence of the signal statistical features
across the array of sensors and consists of three steps:
  * signal processing;
  * space-time imaging;
  * detection and location.



## Getting Started

Clone or download the project from GitHub, if needed uncompress the archive.

### Installation:

Run following command from within the main directory:

    pip install .

or to install developer mode use:

    pip install -e .


### Running examples:
First, download the file [examples.zip](https://www.dropbox.com/s/emlz4lbd6dpu9a9/examples.zip?dl=1) containing additional data (seismograms and theoretical travel-time grids).


Run the main detection and location code on an example dataset:

    btbb  examples/BT_ChileExample.conf

Run an example illustrating the procedure of Multi-Band Filter Characteristic Function calculation:

    mbf_plot  examples/MBF_ChileExample.conf


## Documentation

A detailed documentation is available here: [backtrackbb.readthedocs.io](http://backtrackbb.readthedocs.io/en/latest/)


### Contact Information:
  * [Natalia Poiata](mailto:poiata@ipgp.fr)
  * [Claudio Satriano](mailto:satriano@ipgp.fr)


### References
Poiata, N., C. Satriano, J.-P. Vilotte, P. Bernard, and K. Obara (2016). Multi-band array detection and location of seismic sources recorded by dense seismic networks, Geophys. J. Int., 205(3), 1548-1573, doi:[10.1093/gji/ggw071](https://doi.org/10.1093/gji/ggw071).
