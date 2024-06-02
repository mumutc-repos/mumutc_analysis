# mumutc_analysis

This a program for analyzing tc final state at muon collider. 
It can read events from a HepMC file and output some features of the events. 
The ouput data were used in this publication: [https://arxiv.org/pdf/2302.01143](https://arxiv.org/pdf/2302.01143).

## Prerequisites
The following packages are necessary to run this program:
* ROOT ([https://root.cern/](https://root.cern/))
* HepMC2 ([http://hepmc.web.cern.ch/hepmc/](http://hepmc.web.cern.ch/hepmc/))
* FastJet ([https://fastjet.fr/](https://fastjet.fr/))

The installation and configuration of these codes please refer their homepages. 

## Usage
Please edit your pathes to these tools in the Makefile to compile the codes by running
```
make
```
You should get an executable **main**, and then you can run 
```
main [input_file] [output_file] [R] [E]
```
where input file must be HepMC2 format and the output file will be ROOT format.
**R** and **E** are two float numbers. 
**R** is the R parameter for kt jet algorithm, and **E** is the energy cut for jets. 
