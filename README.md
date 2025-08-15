# Minimal and Customizable Relative Free Energy with OpenFE

This repository provides a lightweight and customizable wrapper around the `openfe` toolkit, designed for performing Relative Free Energy (RFE) calculations. It streamlines the core functionalities of `openfe` while allowing users to easily customize key aspects, from ligand and system setup to detailed simulation parameters.

This project is inspired by and built upon the concepts shown in the [OpenFE showcase Colab notebook](https://colab.research.google.com/github/OpenFreeEnergy/ExampleNotebooks/blob/main/showcase/openfe_showcase.ipynb).

## Features

- **Customizable:** Easily modify ligand preparation, system setup, and simulation parameters.
- **Minimalist:** A streamlined approach focusing on the core RFE calculation pipeline.
- **Lightweight:** Designed to be a simple wrapper for common use cases.

## Installation

This project requires a `mamba` or `conda` environment. The recommended approach is to create a dedicated environment to manage dependencies.

1.  **Create and activate the `mamba` environment:**

    ```bash
    mamba create -c conda-forge -n min_fe_env openfe=1.6
    mamba activate min_fe_env
    ```

2.  **Install the repository in editable mode:**

    This allows you to make changes to the source code without needing to reinstall the package.

    ```bash
    pip install -e .
    ```

> **Note:** For the latest installation instructions for the `openfe` library itself, please refer to the official [OpenFE documentation](https://docs.openfree.energy/en/stable/installation.html).

## Usage

To perform a default RFE calculation, you only need to provide a single `.sdf` file containing your ligands and a `.pdb` file for the protein.

**Running the calculation with default parameters:**

```bash
python scripts/compute_rbfep.py --ligands_path path/to/ligands.sdf --protein_path path/to/proteins.pdb
```

Sample ligand and protein files are provied at [files/ligands.sdf](files/ligands.sdf) and [files/protein.pdb](files/protein.pdb)