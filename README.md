# nmd-escape
Code to determine NMD escaping regions of transcripts and calculate NMD annotation for frameshift variants

![alt text](https://github.com/ritorene/nmd-escape/blob/main/images/nmd.png?raw=true)

Nonsense-mediated decay (NMD) escaping protein truncating variants can cause Mendelian disease.

The code in this repo contains utility functions for:
* defining NMD escaping regions based on coding sequence bed files
* determining total NMD escape size per transcript
* annotating truncating variants for whether the shifted stop codon would be NMD escaping


# Getting Started
```
git clone git@github.com:rebeccaito/nmd-escape.git
cd nmd-escape
python install -e .

# run unit tests to see whether install worked
 python -m unittest ./tests/test_annotating_nmd.py
```

# Examples
## 1. Determine NMD regions from a 6-column bed file of coding sequence regions
*Note:* the `cds_id` column must follow the conventions for RefSeq transcript CDS downloaded from the UCSC Genome Browser:
e.g., `NM_138576.4_cds_0_0_chr14_99640488_r` is a valid `cds_id` 
```
from annotating_nmd import * 

cds_bed_df = pd.read_table('/path/to/nmd-escape/tests/data/test_cds.bed', names=['chrom', 'start', 'end', 'cds_id', 'score', 'strand'])
print(cds_bed_df)
nmd_bed_df = make_boundaries_df(cds_bed_df)
print(nmd_bed_df)
```

## 2. Determine NMD size per transcript
```
from annotating_nmd import * 

cds_bed_df = pd.read_table('/path/to/nmd-escape/tests/data/test_cds.bed', names=['chrom', 'start', 'end', 'cds_id', 'score', 'strand'])
print(cds_bed_df)

sizes_df = make_cds_size_df(cds_bed_df)
print(sizes_df)
```

## 3. Determine whether frameshifts would be NMD escaping
```
from annotating_nmd import * 

cds_bed_df = pd.read_table('/path/to/nmd-escape/tests/data/test_cds.bed', names=['chrom', 'start', 'end', 'cds_id', 'score', 'strand'])
sizes_df = make_cds_size_df(cds_bed_df)

# variant annotations (from VEP for example)
data = [
    {
        'transcript_name': 'NM_003620.4',
        'HGVSp': 'NP_003611.1:p.Cys300LeufsTer129',
    }, {
        'transcript_name': 'NM_003620.4',
        'HGVSp':'NP_003611.1:p.Cys81LeufsTer129',
    }
]
variant_df = pd.DataFrame(data)
result = get_upstream_frameshift(variant_df, sizes_df)
print(result)
```

