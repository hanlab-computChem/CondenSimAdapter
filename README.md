# CondenSimAdapter

CondenSimAdapter is an automated workflow for protein condensate simulations, covering the main stages from CG to AA. This README only keeps the installation guide.

## Installation

### 1. Create a conda environment

```bash
conda create -n conden python=3.11 -y
conda activate conden
```

### 2. Run the install script

```bash
bash install.sh
```

### 3. Install the local package

```bash
pip install .
```

### 4. Only for COCOMO multidomain protein simulations

```bash
pip install git+https://github.com/feiglab/mdsim.git
pip install numpy==1.26.4 mdtraj==1.11.0
```

Note: `mdsim` claims it needs higher dependency versions, but downgrading `numpy` and `mdtraj` to the versions above does not affect COCOMO simulations.

## Requirements

- Python = 3.11
- CUDA >= 12.1
- GROMACS >= 2023 (install it yourself)