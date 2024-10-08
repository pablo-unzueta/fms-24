# Home
A reimplementation of FMS90 in PyTorch

## Installation
Recommended installation:

`pip install fms24`

One can also install from github:

`pip install git+https://github.com/pablo-unzueta/fms-24/fms-24.git`

## Example Usage
Similar to FMS90's `Control.dat`, the `config.yaml` file is used to specify the parameters of the AIMS run.

A tutorial and a full example of a `config.yaml` file can be found in the examples directory.


## Electronic Structure
Currently, only TeraChem and Q-Chem are supported. If you would like a specific program to be supported, please open an issue.

## ML Potentials
We currently support [Nequip-style potentials](https://github.com/mir-group/nequip). Create potentials for the adiabatic states of interest add the path(s) to your `config.yaml`.

```
ml_potentials:
  - type: "nequip"
    paths:
      - {path: "path/to/first/ml/potential", order: 0}
      - {path: "path/to/second/ml/potential", order: 1}
      - {path: "path/to/third/ml/potential", order: 2}
```

## Acknowledgements
Written by Pablo Unzueta and Todd Martinez. We'd like to thank the NSF MPS-ASCEND (Grant No. X) for financial support.


