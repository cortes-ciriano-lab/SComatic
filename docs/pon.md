## Panel of normals (PoN)
To build your Panel of normals (PoN), the user needs to provide (--in_tsv) a tsv (or txt) file listing the full path to all the desired files/samples obtained by BaseCellCalling.step1.py. The file should look like next:

```
/path/to/sample1.basecalling.step1.tsv
/path/to/sample2.basecalling.step1.tsv
/path/to/sample3.basecalling.step1.tsv
/path/to/sample4.basecalling.step1.tsv
/path/to/sample5.basecalling.step1.tsv
```

The python script takes these parameters:

```
python scripts/PoN/PoN.py --help
usage: PoN.py [-h] --in_tsv IN_TSV --out_file OUT_FILE
              [--min_samples MIN_SAMPLES] [--rm_prefix {Yes,No}]

Script to obtain a Panel Of Normals (PoN) file for scRNA somatic variant
caller

optional arguments:
  -h, --help            show this help message and exit
  --in_tsv IN_TSV       File with tsv files to be used for final PoN
                        construction (ideally files obtained in
                        BaseCellCalling.step1.py)
  --out_file OUT_FILE   PoN output file name
  --min_samples MIN_SAMPLES
                        Minimum number of significant samples to consider a
                        site in the PoN. [Default: 2]
  --rm_prefix {Yes,No}  Remove chr prefix from input files (Yes) or no (No)
                        [Default: Yes]
```
