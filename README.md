# MPBuild/Mp_build Tool Documentation

## English Version

### Tool Overview
**MPBuild/Mp_build** is a highly automated membrane protein system construction tool designed for AMBER force field-based molecular dynamics (MD) simulations. It integrates multiple computational biology tools (e.g., Modeller, Schrodinger Suite, AmberTools24, PyMOL) and offers the following core features:

1. **End-to-End Membrane Protein System Construction**  
   - One-click generation from PDB input to complete simulation system
   - Automatic repair of missing residues, mutation sites, and disulfide bonds
   - Compatibility with non-standard residues and small molecule parameterization
   - Transmembrane cavity water generation (based on 3D-RISM theory)

2. **Key Features**  
   - Dual operation modes: Simplified (auto-repair) and Advanced (user-defined templates)
   - Support for AlphaFold-predicted structures as templates
   - Batch system construction capability
   - Built-in binding free energy calculation module (MM-PBSA/GBSA)

3. **Application Scenarios**  
   - Drug target-membrane protein complex simulations
   - Mutant conformation studies
   - Peptide/small molecule interaction analysis with membrane proteins

---

### Installation Guide

#### 1. Unpack Toolkit
```bash
# Navigate to target directory (example)
cd ~/Downloads/software/

# Extract toolkit
tar -xvf mp_build_v9.tar.gz
```

#### 2. Configure Database Files
```bash
# Create database directories
mkdir -p mp_build_v9/databases/{opm_pdbs,nonaa}

# Download and deploy essential databases
## OPM membrane protein orientation database (~500MB)
wget https://opm.phar.umich.edu/opm_pdbs.tar.gz -P mp_build_v9/databases/
tar -xzf mp_build_v9/databases/opm_pdbs.tar.gz -C mp_build_v9/databases/opm_pdbs

## Modeller sequence database (pdball.pir)
wget https://salilab.org/modeller/supplemental.html/pdball.pir -P mp_build_v9/databases/
```

#### 3. Install Dependencies
```bash
| Software               | Version       | Verification Command          | Key Configuration                                                                 |
|------------------------|---------------|--------------------------------|-----------------------------------------------------------------------------------|
| **AmberTools24**       | 24.0+         | `tleap -v`                    | Set environment variables:<br>`export AMBERHOME=/path/to/amber24`<br>`source $AMBERHOME/amber.sh` |
| **Modeller**           | 10.6          | `mod9.26 -h`                  | Configure license key:<br>`export KEY_MODELLER=XXXX-XXXX-XXXX-XXXX`                |
| **PyMOL**              | 2.5.5+        | `pymol -cq`                   | Ensure executable path is in system PATH                                         |
| **Schrodinger Suite**  | 2018+         | `$SCHRODINGER/run prepwizard` | Installation path must contain `utilities/` subdirectory:<br>`export SCHRODINGER=/path/to/schrodinger` |
| **Gaussian**           | G09+          | Optional                      | Configure node path:<br>`export GAUSS_EXEDIR=/path/to/g16/bsd`                   |
| **R**                  | 4.0+          | `R --version`                 | Install required packages:<br>`install.packages(c("ggplot2", "dplyr", "MASS"))`  |
```

#### 4. Verify Directory Structure
```bash
mp_build_v9/
├── build.sh                # Main execution script
├── databases/
│   ├── nonaa/              # User-defined non-standard residue/molecule parameters
│   ├── opm_pdbs/           # OPM database files (*.pdb)
│   └── pdball.pir          # Modeller homology sequence database
├── scripts/
│   ├── homo_build/         # Homology modeling & structural repair scripts
│   ├── sys_build/          # System assembly scripts
│   ├── make_nonaa_lig/     # Small molecule/residue parameter generation scripts
│   ├── run_mmpbsa_amber/   # Free energy calculation scripts
│   └── build.sh            # Main execution script
└── sample/                 # Example files (with test cases)
```

#### 5. Core Parameters
```bash
| Parameter      | Required/Optional | Description                                                              |
|----------------|-------------------|--------------------------------------------------------------------------|
| `-m_path`      | Required          | MPBuild toolkit directory path (e.g., `~/mp_build_v9`)                 |
| `-s_path`      | Required          | Schrodinger installation path (must contain `utilities/` subdirectory) |
| `-p_com`       | Required          | Target PDB file path (must contain hydrogen atoms and complete structure) |
| `-p_tmm`       | Required          | Transmembrane chain ID (single chain only, e.g., `A`)                  |
| `-p_cid`       | Optional          | Specify modeling chain IDs (repeatable, e.g., `-p_cid A -p_cid B`)     |
| `-p_seq`       | Optional          | FASTA template sequence file (must include chain ID tags, e.g., `>A`)  |
| `-p_tpt`       | Optional          | Template structure PDB file (recommend AlphaFold-predicted full-length structures) |
| `-c_lig`       | Optional          | Ligand PDB file (requires pre-generated parameters via `make_lig_v5.pl`) |
| `-c_pep`       | Optional          | Peptide PDB file (requires non-standard residue parameters via `make_nonaa_v5.pl`) |
| `-w_inh`       | Optional          | Cavity water control (`0`=disable/`1`=enable/`*.pdb`=pre-generated water file) |
```

#### 6. Basic Usage Examples
```bash
# Simplified mode (auto-repair all domains)
bash build.sh -m_path ~/mp_build_v9 \
              -s_path /opt/schrodinger2023 \
              -p_com 7xau.pdb \
              -p_tmm A

# Advanced mode (with templates and ligands)
bash build.sh -m_path ~/mp_build_v9 \
              -s_path /opt/schrodinger2023 \
              -p_com 9bi6.pdb \
              -p_tmm R \
              -p_cid R -p_cid A -p_cid B \
              -p_seq complex_seq.txt \
              -p_tpt template.pdb \
              -c_pep oct.pdb \
              -w_inh 1

# Generate intermediate residue parameters (requires capped structures)
perl make_nonaa_inr_v5.pl -pdb residue_capped.pdb \
                          -mol2 residue_nocap.mol2 \
                          -name RES_NAME \
                          -gau g16 \
                          -np 8 -ac 0

# Generate small molecule parameters (GAFF2 force field)
perl make_lig_v5.pl -pdb ligand.pdb \
                    -mol2 ligand.mol2 \
                    -name LIG \
                    -gau g16 \
                    -np 4
```

# Case Study: G Protein-Coupled Receptor Complex System Construction
## Case Background

### This case demonstrates how to construct a G protein-coupled receptor complex system containing a transmembrane domain (Chain R) and three intracellular subunits (Chains A/B/G). The original PDB file 9bi6.pdb from RCSB requires mutation repair and missing structure completion.

## File Preparation
### 1. Obtain Original Structure
```bash
wget https://files.rcsb.org/download/9bi6.pdb
```
### 2. Prepare Template Sequence File
```bash
Create complex_seq.txt containing Uniprot sequences for four domains:
>R
MELTIV... (Full Q15743 sequence)
>A
MGSKGE... (Full P50148 sequence)
>B
MAAVAG... (Full P62873 sequence)
>G
MGLQDS... (Full P59768 sequence)

​Format Requirements​​:

    Each domain starts with >[Chain ID]
    Sequence length must match actual domains (verify residue numbering via PyMOL)
```
### 3. Generate AlphaFold Template Structures
```bash
    Access Uniprot for predicted structures:
        Q15743: https://www.uniprot.org/uniprotkb/Q15743/entry#structure

    Repeat for other IDs

Merge downloaded AF2 structures into template.pdb:

cat AF2_Q15743.pdb AF2_P50148.pdb AF2_P62873.pdb AF2_P59768.pdb > template.pdb
```

### 4. Pre-generate Cavity Waters
```bash
First-time run to generate water file:

bash build.sh -m_path ~/mp_build_v9 \
              -s_path /opt/schrodinger2023 \
              -p_com 9bi6.pdb \
              -p_tmm R \
              -w_inh 1 \
              -p_cid R

Copy generated wats_inhole_del_2.pdb to working directory.
```

### 5. Execute Modeling Command
```bash
bash build.sh -m_path ~/mp_build_v9 \
              -s_path /opt/schrodinger2023 \
              -p_com 9bi6.pdb \
              -p_tmm R \
              -p_cid R -p_cid A -p_cid B -p_cid G \
              -p_seq complex_seq.txt \
              -p_tpt template.pdb \
              -w_inh wats_inhole_del_2.pdb
```

### 6. Key Parameter Analysis (Case-Specific)
```bash
| Parameter       | Type         | Description                                                             | Input Requirements                                                     |
|-----------------|--------------|-------------------------------------------------------------------------|-------------------------------------------------------------------------|
| `-p_cid`        | Optional     | Specify modeling chain IDs                                             | Requires repeated use (e.g., `-p_cid R -p_cid A`)                      |
| `-p_seq`        | Optional     | Provide native sequence template for mutation/defect repair            | FASTA format with `>[ChainID]` headers                                 |
| `-p_tpt`        | Optional     | Use AlphaFold-predicted structures for missing region modeling         | Requires merged PDB of all domains (via `cat` command)                 |
| `-w_inh`        | Optional     | Control transmembrane cavity water generation                          | Supports three modes:<br>• `0`=disable<br>• `1`=first-time generation<br>• `*.pdb`=reuse pre-generated file |
```

### 7. Parameter Interactions
```bash
| Parameter Combination          | System Behavior                                                                 |
|---------------------------------|---------------------------------------------------------------------------------|
| `-p_seq` + `-p_tpt`            | Prioritizes user template sequences for repair with AlphaFold structural guidance |
| `-p_cid` + `-p_tmm`            | Models specified chains only (transmembrane chain must be included in `-p_cid`) |
| `-w_inh 1` + `-p_tmm`          | Automatically calculates cavity water positions based on specified TM chain (requires OPM database) |
```
