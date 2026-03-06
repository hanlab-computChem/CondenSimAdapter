#!/usr/bin/env python

import os
import sys
import mdtraj
import pathlib
import numpy as np
import functools

import dgl
import torch

from .libconfig import MODEL_HOME
from .model_downloader import ensure_model_available
from .libdata import (
    PredictionData,
    create_trajectory_from_batch,
    create_topology_from_data,
)
from .residue_constants import read_coarse_grained_topology
from . import libcg
from .libpdb import write_SSBOND
from .libter import patch_termini
from . import libmodel

import warnings

warnings.filterwarnings("ignore")

# Models that support fix_atom parameter (FIX version)
MODELS_WITH_FIX = {"CalphaBasedModel"}

# Supported models
SUPPORTED_MODELS = {
    "CalphaBasedModel": libcg.CalphaBasedModel,
    "ResidueBasedModel": libcg.ResidueBasedModel,
    "Martini": libcg.Martini,
    "Martini3": libcg.Martini3,
}


def convert_cg2all(
    in_pdb_fn,
    out_fn,
    model_type="CalphaBasedModel",
    in_dcd_fn=None,
    ckpt_fn=None,
    fix_atom=False,
    device=None,
    n_proc=int(os.getenv("OMP_NUM_THREADS", 1)),
    write_ssbond=False,
):
    """
    Convert coarse-grained protein structure to all-atom model.
    
    Args:
        in_pdb_fn: Input CG PDB file path
        out_fn: Output AA PDB file path
        model_type: Model type (CalphaBasedModel, ResidueBasedModel, Martini, Martini3)
        in_dcd_fn: Optional DCD trajectory file
        ckpt_fn: Optional custom checkpoint path
        fix_atom: For CalphaBasedModel only - preserve CA coordinates from input
        device: Device to use ("cuda" or "cpu")
        n_proc: Number of processes for data loading
        write_ssbond: Whether to write SSBOND records for disulfide bonds (default: False)
    """
    # Validate model type
    if model_type not in SUPPORTED_MODELS:
        raise ValueError(
            f"Unsupported model type: {model_type}\n"
            f"Supported models: {list(SUPPORTED_MODELS.keys())}"
        )
    
    # set device
    if device is None:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    else:
        device = torch.device(device)

    # load model ckpt file
    if ckpt_fn is None:
        # Only use FIX version for CalphaBasedModel when fix_atom=True
        if model_type in MODELS_WITH_FIX and fix_atom:
            model_name = f"{model_type}-FIX"
            ckpt_fn = MODEL_HOME / f"{model_type}-FIX.ckpt"
        else:
            model_name = model_type
            ckpt_fn = MODEL_HOME / f"{model_type}.ckpt"
        
        # Try to download if not exists
        if not ckpt_fn.exists() or ckpt_fn.stat().st_size < 1024*1024:
            try:
                ckpt_fn = ensure_model_available(model_name, MODEL_HOME)
            except Exception as e:
                # Download failed, show helpful error
                error_msg = f"Checkpoint file not found: {ckpt_fn}\n"
                error_msg += f"\nDownload error: {e}\n"
                error_msg += "\nTo use backmapping, you need to either:\n"
                error_msg += "1. Install huggingface_hub: pip install huggingface_hub\n"
                error_msg += f"2. Manually download models from: https://huggingface.co/hanlab/condensimadapter-cg2all-models\n"
                error_msg += f"3. Place models in: {MODEL_HOME}\n"
                raise FileNotFoundError(error_msg)
    
    # Check if checkpoint file exists
    if not ckpt_fn.exists():
        error_msg = f"Checkpoint file not found: {ckpt_fn}\n"
        
        # Show available models
        error_msg += "\nAvailable local models:\n"
        for model in SUPPORTED_MODELS:
            normal_ckpt = MODEL_HOME / f"{model}.ckpt"
            fix_ckpt = MODEL_HOME / f"{model}-FIX.ckpt"
            
            normal_status = "✓" if normal_ckpt.exists() and normal_ckpt.stat().st_size > 1024*1024 else "✗"
            if model in MODELS_WITH_FIX:
                fix_status = "✓" if fix_ckpt.exists() and fix_ckpt.stat().st_size > 1024*1024 else "✗"
                error_msg += f"  {model}: normal={normal_status}, FIX={fix_status}\n"
            else:
                error_msg += f"  {model}: normal={normal_status}\n"
        
        raise FileNotFoundError(error_msg)
    
    # Check file size
    file_size = ckpt_fn.stat().st_size
    if file_size < 1000:
        raise ValueError(
            f"Checkpoint file is too small or empty: {ckpt_fn} ({file_size} bytes)"
        )
    
    print(f"Loading checkpoint from: {ckpt_fn}")
    ckpt = torch.load(ckpt_fn, map_location=device)
    config = ckpt["hyper_parameters"]

    # configure model
    cg_model = SUPPORTED_MODELS.get(config["cg_model"])
    if cg_model is None:
        raise ValueError(f"Unknown CG model in checkpoint: {config['cg_model']}")
    
    config = libmodel.set_model_config(config, cg_model, flattened=False)
    model = libmodel.Model(config, cg_model, compute_loss=False)

    # update state_dict
    state_dict = ckpt["state_dict"]
    for key in list(state_dict):
        state_dict[".".join(key.split(".")[1:])] = state_dict.pop(key)
    model.load_state_dict(state_dict)
    model = model.to(device)
    model.set_constant_tensors(device)
    model.eval()

    # Read unitcell information from input PDB (if available)
    input_unitcell_lengths = None
    input_unitcell_angles = None
    if in_dcd_fn is None:
        # Try to read unitcell info from input PDB
        try:
            input_traj = mdtraj.load(in_pdb_fn, standard_names=False)
            if input_traj.unitcell_lengths is not None and input_traj.unitcell_angles is not None:
                # mdtraj stores unitcell_lengths/angles as (n_frames, 3) arrays
                # For single frame, take the first frame
                if input_traj.n_frames > 0:
                    input_unitcell_lengths = input_traj.unitcell_lengths[0]
                    input_unitcell_angles = input_traj.unitcell_angles[0]
        except Exception:
            # If reading fails, continue without unitcell info
            pass
    
    # prepare input data
    input_s = PredictionData(
        in_pdb_fn,
        cg_model,
        dcd_fn=in_dcd_fn,
        radius=config.globals.radius,
        fix_atom=config.globals.fix_atom,
    )
    if in_dcd_fn is not None:
        unitcell_lengths = input_s.cg.unitcell_lengths
        unitcell_angles = input_s.cg.unitcell_angles
    if len(input_s) > 1 and n_proc > 1:
        input_s = dgl.dataloading.GraphDataLoader(
            input_s, batch_size=1, num_workers=n_proc, shuffle=False
        )

    if in_dcd_fn is None:  # PDB file
        batch = input_s[0].to(device)
        #
        with torch.no_grad():
            R = model.forward(batch)[0]["R"]
        #
        traj_s, ssbond_s = create_trajectory_from_batch(batch, R)
        output = patch_termini(traj_s[0])
        
        # Set unitcell information if available from input
        if input_unitcell_lengths is not None and input_unitcell_angles is not None:
            # Convert to numpy arrays with shape (1, 3) for single frame
            output.unitcell_lengths = np.array([input_unitcell_lengths])
            output.unitcell_angles = np.array([input_unitcell_angles])
        
        output.save(out_fn)
        
        # Write SSBOND records only if write_ssbond=True
        if write_ssbond and len(ssbond_s[0]) > 0:
            write_SSBOND(out_fn, output.top, ssbond_s[0])

    else:  # DCD file
        xyz = []
        for batch in input_s:
            batch = batch.to(device)
            #
            with torch.no_grad():
                R = model.forward(batch)[0]["R"].cpu().detach().numpy()
                mask = batch.ndata["output_atom_mask"].cpu().detach().numpy()
                xyz.append(R[mask > 0.0])
        #
        top, atom_index = create_topology_from_data(batch)
        xyz = np.array(xyz)[:, atom_index]
        traj = mdtraj.Trajectory(
            xyz=xyz,
            topology=top,
            unitcell_lengths=unitcell_lengths,
            unitcell_angles=unitcell_angles,
        )
        output = patch_termini(traj)
        output.save(out_fn)

    return output
