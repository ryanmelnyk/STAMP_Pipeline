# Processing reads for STAMP

Collaboration with Finlay Lab

### Pipeline

ProcessReads.py takes two arguments: the directory of raw reads and an output directory to store files.

### Naming

The folder containing the reads should be renamed to have no spaces.  Also the name of each sample should be unique and contain no periods. The sample name will be everything up to the final dash.

i.e 1-3_Feces_Day3-57133183 will become 1-3_Feces_Day3

### Final output

The main file is "barcode_matrix.txt".  This contains a line for each unique barcode cluster (with >97% identity) with each column corresponding to a sample.
