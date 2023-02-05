# NOLAN
TENET(Transfer Entropy-based causal gene NETwork) refined

## Original
Please visit the orignal [TENET Github](https://github.com/neocaleb/TENET).
Also take a look at the original [TENET Paper](https://doi.org/10.1093/nar/gkaa1014).

## Improvements
- refined and packaged as python3 module for improved usability
- explicit parallelization by python multiprocessing, does not depend on openmpi
- utilizes shared memory to avoid overhead from collecting results from threads
- core java code remains unchanged to guarentee same results

## Dependencies
- Python >= 3.10
- [JPype](https://github.com/jpype-project/jpype) >= 1.4.0
- pandas >= 1.4
- numpy >= 1.23
- tqdm >= 4.0

## Installation
This module is not on pypi.
```bash
git clone --recurse-submodules https://github.com/Stfort52/NOLAN.git
```

## Reference
Please refer to the repository [wiki](https://github.com/Stfort52/NOLAN/wiki).