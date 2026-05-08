# Computational supplement for cycle-ratio obstructions

This repository contains optional computational materials for the manuscript
``Cycle-ratio obstructions to low-period cycles of
psi_{d,a,b}(z)=z/(a z^d+b)''.

The computations are reproducibility and sanity checks only. The manuscript's
proofs are self-contained and do not rely on these scripts as black boxes.

## Contents

- `code/verify_theorems.py`: Python/SymPy checks for the displayed polynomial
  identities and finite illustrative orbit searches.
- `code/magma_verify_and_explore.m`: Magma checks for the corresponding
  polynomial identities and exploratory computations.
- `outputs/verify_log.txt`: recorded output from the Python/SymPy checks.
- `outputs/magma_verify_log.txt`: recorded output from the Magma checks.
- `requirements.txt`: minimal optional Python dependency list.

## Reproduce the Python checks

```sh
python3 -m pip install -r requirements.txt
python3 code/verify_theorems.py
```

The included Python log was generated with Python 3.9.6 and SymPy 1.14.0.
The script is version-light and should also run in current SymPy installations.

## Reproduce the Magma checks

```sh
magma -b code/magma_verify_and_explore.m
```

The included Magma log was generated with Magma V2.29-6.

## Notes

The symbolic portions check the displayed polynomial identities used for
comparison in the manuscript. The finite orbit enumeration is illustrative and
is not used as evidence for the theorems.
