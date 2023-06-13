# nmd-escape
Code to determine NMD escaping regions of transcripts and calculate NMD annotation for frameshift variants

![alt text](https://github.com/ritorene/nmd-escape/blob/main/images/nmd.png?raw=true)

Nonsense-mediated decay (NMD) escaping protein truncating variants can cause Mendelian disease.

The code in this repo contains utility functions for:
* defining NMD escaping regions based on coding sequence bed files
* determining total NMD escape size per transcript
* annotating truncating variants for whether the shifted stop codon would be NMD escaping
