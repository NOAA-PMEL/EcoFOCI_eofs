# EcoFOCI_eofs

This is a wrapper to smoothly analyze EcoFOCI timeseries data via EOF analysis.

Input files are EPIC standard netCDF files.  Output files are EPIC standard netCDF files and summary text file.

Tested on python 3.6, 2.7

## Based on python EOF analysis package

*** eofs V1.3.0***

[http://ajdawson.github.io/eofs/index.html](http://ajdawson.github.io/eofs/index.html)  
[https://github.com/ajdawson/eofs](https://github.com/ajdawson/eofs)  
[http://doi.org/10.5334/jors.122](http://doi.org/10.5334/jors.122) - journal article  

Many results are not relevant here for these purposes from this package and are available as analysis extensions.

### Science Notes

We expect to have the same number of EOF's / PCS's as we have timeseries

### Usage

```text
usage: EcoFOCI_eof.py [-h] [-o OUTFILE] [-s START_DATE] [-e END_DATE]
                      [--eof_num EOF_NUM] [--epic]
                      pfile varname config_file_name

Analyze timeseries EOF

positional arguments:
  pfile                 pointer file with full paths on each line
  varname               name of variable, may be EPIC name
  config_file_name      full path to config file - eof_config.yaml

optional arguments:
  -h, --help            show this help message and exit
  -o OUTFILE, --outfile OUTFILE
                        name of output file (run-name)
  -s START_DATE, --start_date START_DATE
                        yyyymmddhhmmss
  -e END_DATE, --end_date END_DATE
                        yyyymmddhhmmss
  --eof_num EOF_NUM     number of eofs. default is all (arbitrarily large
                        number)
  --epic                assume EPIC time format
  --plots               output some basic plots - TODO
  --summary             output summary only
```

#### Example

```bash
EcoFOCI_eof.py test/14ck9a_curr_brg_subset_f35.point U_320 config/eof_config.yaml -s 20141005000000 -e 20150501000000 --epic
```

Will run the analysis files listed in: `test/14ck9a_curr_brg_subset_f35.point`

for the variable: U_320 (which will also select u_1205)  
with the EPIC .nc file having properties specified by `config/eof_config.yaml`

starting at 2014-10-05 00:00:00  
ending at 2015-05-01 00:00:00  
and expecting  EPIC format

```bash
EcoFOCI_eof.py test/14ck9a_curr_brg_subset_f35.point U_320 config/eof_config.yaml -s 20141005000000 -e 20150501000000 -eof_num=1 -o=test --epic
```

is the same as the above but only returning the first EOF and naming the run 'test'

PJS
---

```bash
progdir=/full/path/to/dir/  
pointer=14ck9a_curr_brg_subset_f35.pavlof.point

python ${progdir}EcoFOCI_eof.py ${pointer} U_320 ${progdir}config/eof_config.yaml -s 20141005000000 -e 20150501000000 --epic -o eof_results
```

## Requirements

Found in the ci/requirements.txt

## Legal Disclaimer

This repository is a scientific product and is not official communication of the National Oceanic and Atmospheric Administration (NOAA), or the United States Department of Commerce (DOC).
All NOAA GitHub project code is provided on an 'as is' basis and the user assumes responsibility for its use.
Any claims against the DOC or DOC bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law.
Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation, or favoring by the DOC.
The DOC seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by the DOC or the United States Government.
