# maf_matrix

### Create matrix plots from TCGA MAF files

Simply provide a path to a maf file and a list of genes, and the script can either:

<ol type="a">
  <li>save the matrix figure (to an optional path). (The "-a" flag shows all patients.), or</li>
  <li>look up the amino acid changes from Broad's Oncotator and print out text that can be used in the CBioPortal oncoprinter tool.</li>
</ol>

## Using the executable
```
plot_maf_matrix [-h] -i MAF_PATH [-t | -o OUT_PATH] [-a]
                       [-g GENES [GENES ...]]

optional arguments:
  -h, --help            show this help message and exit
  -i MAF_PATH, --maf_path MAF_PATH
                        Path to maf file.
  -t, --oncoprinter
  -o OUT_PATH, --out_path OUT_PATH
                        Output figure path including extension.
  -a, --all_patients    Show all patients in matrix plot.
  -g GENES [GENES ...], --genes GENES [GENES ...]
                        Gene hugo symbols, case sensitive.
```
e.g. for first matrix plot shown above:
```
plot_maf_matrix -i hgsc.bcm.edu_CHOL.IlluminaGA_DNASeq.1.somatic.maf -o test_matrix.pdf -g MUC4 TP53 KRAS BRAF NF1
```

## Using the python module
https://dl.dropboxusercontent.com/u/5141228/matrix_demo.html


## Installation
- unzip the archive and run
```
python setup.py install
```

## License

    Copyright (C) 2015 Stephen Gaffney

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.