# CondenSimAdapter

CondenSimAdapter is an automated workflow for protein condensate simulations, covering the main stages from CG to AA. This README only keeps the installation guide.

## Installation

### 1. Create a conda environment

```bash
conda create -n conden python=3.11 -y
conda activate conden
```

### 2. Download and install dependency

Download the zip from GitHub webpage or run:

```bash
wget https://github.com/pearlDingzhen/CondenSimAdapter/archive/refs/heads/master.zip
```

Unzip and enter the directory:

```bash
unzip CondenSimAdapter-master.zip
cd CondenSimAdapter-master
```

Run the installation script:

```bash
bash install.sh
```

### 3. Install CondenSimAdapter

```bash
pip install .
```

### 4. For COCOMO multidomain protein simulations (Optional)

```bash
pip install git+https://github.com/feiglab/mdsim.git
pip install numpy==1.26.4 mdtraj==1.11.0
```

Note: `mdsim` claims it needs higher dependency versions, but downgrading `numpy` and `mdtraj` to the versions above does not affect COCOMO simulations.

### 5. Run test 

```bash
python -m pytest
```

### 6. Run with CondenSimAdapter command adapter

```bash
adapter -h 
```

## Requirements

- Python = 3.11
- CUDA >= 12.1
- GROMACS >= 2023 (install it yourself)