# CondenSimAdapter API Documentation

## Overview

This document describes the API for the CondenSimAdapter package, including all public classes, methods, and their usage.

## Core Classes

### CalvadosGenerator

The main class for generating and running CALVADOS simulations.

#### Constructor

```python
CalvadosGenerator(config_path)
```

**Parameters:**
- `config_path` (str): Path to the YAML configuration file

**Example:**
```python
from CondenSimAdapter.calvados_generator import CalvadosGenerator

generator = CalvadosGenerator("config.yaml")
```

#### Methods

##### `setup_calvados_environment(output_dir)`

Sets up the CALVADOS simulation environment by creating input directories and copying necessary files.

**Parameters:**
- `output_dir` (str): Directory where the simulation will be run

**Returns:**
- `str`: Path to the input directory

**Example:**
```python
input_dir = generator.setup_calvados_environment("output_simulation")
```

##### `generate_idp_prepare_script(output_dir)`

Generates the prepare.py script for IDP (Intrinsically Disordered Proteins) simulations.

**Parameters:**
- `output_dir` (str): Directory where the script will be written

**Returns:**
- `str`: Path to the generated prepare.py script

**Example:**
```python
prepare_script = generator.generate_idp_prepare_script("output_idp")
```

##### `generate_mdp_prepare_script(output_dir)`

Generates the prepare.py script for MDP (Membrane-Disrupting Peptides) simulations.

**Parameters:**
- `output_dir` (str): Directory where the script will be written

**Returns:**
- `str`: Path to the generated prepare.py script

**Example:**
```python
prepare_script = generator.generate_mdp_prepare_script("output_mdp")
```

##### `generate_and_run(output_dir, protein_name=None, gpu_id=0, replica=1)`

Generates all necessary scripts and returns the command to run the CALVADOS simulation.

**Parameters:**
- `output_dir` (str): Directory where the simulation will be run
- `protein_name` (str, optional): Name of the protein. If None, uses the name from config
- `gpu_id` (int, optional): GPU device ID to use (default: 0)
- `replica` (int, optional): Replica number (default: 1)

**Returns:**
- `str`: Command string to run the simulation

**Example:**
```python
cmd = generator.generate_and_run("output_sim", protein_name="FUS_LC", gpu_id=0, replica=1)
```

#### Properties

##### `config`

The loaded configuration dictionary.

**Type:** `dict`

**Example:**
```python
print(generator.config['cg_calvados']['task_type'])
```

##### `calvados_config`

The CALVADOS-specific configuration section.

**Type:** `dict`

**Example:**
```python
print(generator.calvados_config['steps'])
```

##### `protein_name`

The name of the protein from the configuration.

**Type:** `str`

**Example:**
```python
print(f"Simulating protein: {generator.protein_name}")
```

##### `protein_nmol`

The number of protein molecules from the configuration.

**Type:** `int`

**Example:**
```python
print(f"Number of molecules: {generator.protein_nmol}")
```

## Configuration Structure

### Input Files Configuration

```python
input_files = {
    'sequence_fasta': 'protein.fasta',  # For IDP tasks
    'structure_pdb': 'protein.pdb',     # For MDP tasks
}
```

### Protein Configuration

```python
protein = {
    'name': 'protein_name',
    'nmol': 20,
}
```

### CALVADOS Configuration

```python
cg_calvados = {
    'task_type': 'IDP',  # or 'MDP'
    'platform': 'CUDA',
    'temp': 310.0,
    'box': [25.0, 25.0, 50.0],
    'steps': 200000,
    'wfreq': 5000,
    'topol': 'slab',  # 'random', 'slab', or 'grid'
    'slab_width': 20,
    'ionic': 0.15,
    'pH': 7.0,
    'restraint': False,  # True for MDP, False for IDP
    'fdomains': 'domains.yaml',  # Required for MDP tasks
    # ... other parameters
}
```

## Usage Examples

### Basic IDP Simulation

```python
from CondenSimAdapter.calvados_generator import CalvadosGenerator
import subprocess

# Initialize generator
generator = CalvadosGenerator("config_idp.yaml")

# Generate and get command
cmd = generator.generate_and_run("output_idp", protein_name="FUS_LC", gpu_id=0, replica=1)

# Run the simulation
result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

if result.returncode == 0:
    print("Simulation completed successfully!")
else:
    print("Simulation failed:", result.stderr)
```

### Basic MDP Simulation

```python
from CondenSimAdapter.calvados_generator import CalvadosGenerator
import subprocess

# Initialize generator
generator = CalvadosGenerator("config_mdp.yaml")

# Generate and get command
cmd = generator.generate_and_run("output_mdp", protein_name="TDP43", gpu_id=0, replica=1)

# Run the simulation
result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

if result.returncode == 0:
    print("Simulation completed successfully!")
else:
    print("Simulation failed:", result.stderr)
```

### Custom Configuration

```python
from CondenSimAdapter.calvados_generator import CalvadosGenerator
import yaml

# Create custom configuration
config = {
    'input_files': {
        'sequence_fasta': 'custom_protein.fasta'
    },
    'protein': {
        'name': 'CustomProtein',
        'nmol': 25
    },
    'cg_calvados': {
        'task_type': 'IDP',
        'platform': 'CUDA',
        'temp': 300.0,
        'box': [30.0, 30.0, 60.0],
        'steps': 500000,
        'wfreq': 10000,
        'topol': 'slab',
        'slab_width': 25,
        'ionic': 0.1,
        'pH': 6.5
    }
}

# Write configuration to file
with open('custom_config.yaml', 'w') as f:
    yaml.dump(config, f)

# Use the custom configuration
generator = CalvadosGenerator("custom_config.yaml")
cmd = generator.generate_and_run("custom_output", protein_name="CustomProtein", gpu_id=0, replica=1)
```

### Multiple Replicas

```python
from CondenSimAdapter.calvados_generator import CalvadosGenerator
import subprocess

generator = CalvadosGenerator("config.yaml")

# Run multiple replicas
for replica in range(1, 6):  # 5 replicas
    output_dir = f"output_replica_{replica}"
    cmd = generator.generate_and_run(output_dir, protein_name="FUS_LC", gpu_id=0, replica=replica)
    
    print(f"Running replica {replica}...")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    if result.returncode == 0:
        print(f"Replica {replica} completed successfully!")
    else:
        print(f"Replica {replica} failed:", result.stderr)
```

## Error Handling

### Common Exceptions

#### `FileNotFoundError`

Raised when required input files are not found.

```python
try:
    generator = CalvadosGenerator("config.yaml")
except FileNotFoundError as e:
    print(f"Input file not found: {e}")
```

#### `ValueError`

Raised when configuration parameters are invalid.

```python
try:
    generator = CalvadosGenerator("config.yaml")
    cmd = generator.generate_and_run("output")
except ValueError as e:
    print(f"Configuration error: {e}")
```

#### `ImportError`

Raised when CALVADOS is not available.

```python
try:
    from CondenSimAdapter.calvados_generator import CalvadosGenerator
except ImportError as e:
    print(f"CALVADOS not available: {e}")
```

## Best Practices

### 1. Configuration Validation

Always validate your configuration before running simulations:

```python
generator = CalvadosGenerator("config.yaml")

# Check essential parameters
print(f"Task type: {generator.calvados_config['task_type']}")
print(f"Protein: {generator.protein_name}")
print(f"Molecules: {generator.protein_nmol}")
print(f"Steps: {generator.calvados_config['steps']}")
```

### 2. Output Directory Management

Create unique output directories for different simulations:

```python
import os
from datetime import datetime

# Create timestamped output directory
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
output_dir = f"simulation_{timestamp}"

cmd = generator.generate_and_run(output_dir, protein_name="FUS_LC", gpu_id=0, replica=1)
```

### 3. Error Handling

Implement proper error handling for long-running simulations:

```python
import subprocess
import time

cmd = generator.generate_and_run("output_sim", protein_name="FUS_LC", gpu_id=0, replica=1)

try:
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=3600)  # 1 hour timeout
    
    if result.returncode == 0:
        print("Simulation completed successfully!")
    else:
        print("Simulation failed:", result.stderr)
        
except subprocess.TimeoutExpired:
    print("Simulation timed out after 1 hour")
except Exception as e:
    print(f"Unexpected error: {e}")
```

### 4. Resource Management

Monitor GPU usage and system resources:

```python
import subprocess

# Check GPU availability
try:
    result = subprocess.run(['nvidia-smi', '--query-gpu=index,name,memory.used,memory.total', '--format=csv'], 
                          capture_output=True, text=True)
    print("GPU Status:")
    print(result.stdout)
except FileNotFoundError:
    print("nvidia-smi not found. GPU monitoring not available.")
```
