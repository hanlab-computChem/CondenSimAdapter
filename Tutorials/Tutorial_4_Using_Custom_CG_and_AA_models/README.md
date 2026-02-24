# Tutorial 4: Using Custom CG and AA Models

## Tutorial 4.1: Using a Custom All-Atom Force Field

This tutorial demonstrates how to register and use a custom GROMACS all-atom force field with CondenSimAdapter.

## Overview

`adapter minimize` relies on GROMACS (`gmx pdb2gmx`, `gmx solvate`, etc.) to generate an all-atom topology and (optionally) build an explicit-solvent system. By default, CondenSimAdapter provides several built-in all-atom force fields. When you need a force field that is not built in, you can register a local `*.ff` directory (AMBER- or CHARMM-family) and then refer to it by an assigned ID (`a1`, `a2`, ...).

## Prerequisites

- `CondenSimAdapter` installed and `adapter` available in your environment
- A working GROMACS installation (`gmx` command)
- A GROMACS force field directory `*.ff` — this tutorial provides an example in `Tutorials/Tutorial_4_Using_Custom_CG_and_AA_models/AA/amber03.ff`

## Step 0: Enter the tutorial directory

```bash
cd Tutorials/Tutorial_4_Using_Custom_CG_and_AA_models/AA
```

## Step 1: Inspect the `forcefield` CLI

The `adapter forcefield` subcommand manages the registry of all-atom force fields. Run the help to see available subcommands:

```bash
adapter forcefield -h
```

```text
Usage: adapter forcefield [OPTIONS] COMMAND [ARGS]...

  Manage all-atom force fields (add/remove/list custom entries).

Options:
  -h, --help  Show this message and exit.

Commands:
  add     Register a custom all-atom force field and assign an aN id.
  list    List built-in and custom all-atom force fields.
  remove  Remove a custom all-atom force field by aN id.
```

The `add` command registers a new force field. Its key options are:

```bash
adapter forcefield add -h
```

```text
Usage: adapter forcefield add [OPTIONS]

  Register a custom all-atom force field and assign an aN id.

Options:
  --ff-dir DIRECTORY       Path to custom force field directory (*.ff).
                           [required]
  --pdb2gmx-name TEXT      Name used for gmx pdb2gmx -ff (default: inferred
                           from ff-dir).
  --family [amber|charmm]  Force field family.  [default: AMBER]
  --water-model TEXT       Water model for pdb2gmx.  [default: tip3p]
  --solvate-cs TEXT        Water model used by gmx solvate -cs.  [default:
                           spc216]
  --description TEXT       Optional human-readable description.
  -h, --help               Show this message and exit.
```

## Step 2: Register a custom force field

Register `amber03.ff` as a custom force field. A few notes on the key options:

- **`--family`**: Categorizes the force field as AMBER-like or CHARMM-like, enabling compatible defaults and checks.
- **`--water-model`**: The water model name used by `gmx pdb2gmx` (what `pdb2gmx -water` expects).
- **`--solvate-cs`**: The coordinate set used by `gmx solvate -cs` (a pre-built water box structure).

These do not need to be identical, but they must be consistent in resolution (how many points represent one water molecule). Typically, use `spc216` for 3-point water models and `tip4p` for 4-point water models.

Now, run the registration:

```bash
adapter forcefield add --ff-dir amber03.ff/ --family amber --water-model tip3p --solvate-cs spc216
```

```text
Registered custom force field: a1
  pdb2gmx name: amber03
  family: AMBER
  water model: tip3p
  gbsa mapping: AMBER99SB-ILDN
```

(If you run the same command again, you may see a slightly different summary depending on your version.)

You can verify the force field was registered by listing all available force fields:

```bash
adapter forcefield list
```

```text
All-atom force fields:
------------------------------------------------------------
1-a99SBdisp        | builtin | a99SBdisp
2-amber03wsc       | builtin | amber03wsc
3-amber99sbws-stqp | builtin | amber99sbws-STQp
4-amber99sbws-stq  | builtin | amber99sbws-stq
5-des-amber        | builtin | des-amber
6-des-amber-sf1.0  | builtin | des-amber-SF1.0
7-amber99sb-ildn   | builtin | amber99sb-ildn
8-amber14sb        | builtin | amber14sb_parmbsc1
9-charmm36m        | builtin | charmm36-jul2021
a1                 | custom  | amber03
------------------------------------------------------------
Total: 10
```

The newly registered `amber03` now appears as a custom force field with ID `a1`.

## Step 3: Use the custom force field in `adapter minimize`

Once registered, you can use the ID `a1` anywhere the `-ff/--forcefield` option is accepted:

```bash
adapter minimize -f FUS_LC.yaml -l medium -ff a1 --salt-conc 0.15 --solvate
```

From here on, the workflow is identical to using built-in all-atom force fields.

## (Optional) Clean up: Remove the custom force field

When you no longer need this custom force field, remove it from the registry:

```bash
adapter forcefield remove a1
```

```text
Removed custom force field: a1 (amber03)
```


## Tutorial 4.2: Using a Custom CG Model

This section shows how to use a CG structure generated outside CondenSimAdapter and convert it to an all-atom structure with `adapter backmap`.

## Overview

For custom CG workflows, you only need to:

1. Run your own CG simulation with your preferred package.
2. Prepare a CG snapshot (`.pdb`) as the backmapping input.
3. Choose the correct CG model type in `adapter backmap`.

Supported CG model types are:

- `com` (COM one-bead model)
- `ca` (CA one-bead model)
- `martini2`
- `martini3`

## Step 0: Enter the CG tutorial directory

```bash
cd Tutorials/Tutorial_4_Using_Custom_CG_and_AA_models/CG
```

## Step 1: Run backmapping from CG to all-atom

Use the provided Martini3 example:

```bash
adapter backmap -i m3.pdb -m martini3 -f FUS_LC.yaml
```

Example output:

```text
============================================================
Backmap CG to All-Atom
============================================================

  Input: m3.pdb
  Config: FUS_LC.yaml
  Output: FUS_LC_backmap
  Device: cpu (fixed)
  Model: martini3

Preparing input...
Loading checkpoint from: .../model/Martini3.ckpt
  SLAB topology detected: centering condensate in z direction...
    Protein COM z: 22.49 nm
    Box z center: 21.64 nm
    Offset (z only): -0.85 nm
    ✓ Condensate centered in z direction
  ✓ Backmap completed
  Input PDB: m3.pdb
  Output PDB: FUS_LC_backmap/final.aa.pdb
  Model type: Martini3
```

The generated all-atom structure is:

- `FUS_LC_backmap/final.aa.pdb`

## Step 2: Generate a production-ready all-atom system

After backmapping, run `adapter minimize` to build topology and generate an MD-ready all-atom system:

```bash
adapter minimize -f FUS_LC.yaml -l medium
```

At this stage, the flow is the same as a normal all-atom workflow.