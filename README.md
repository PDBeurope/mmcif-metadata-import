# mmCIF Metadata Importer

This tool imports metadata from mmCIF files into new metadata-only files or into existing models. It uses the gemmi library, with automatic method detection and method-specific CSV specification files.

**Protein Data Bank in Europe (PDBe)** · [pdbe.org](https://www.ebi.ac.uk/pdbe)

## Installation

1. Install the required dependency:
```bash
pip install -r requirements.txt
```

### Jupyter notebook (interactive UI)

A **Jupyter notebook** provides an interactive form (file upload, checkboxes, run button)—no command line or web hosting needed.

**Run in browser (no install):**  
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/PDBeurope/mmcif-metadata-import/HEAD?labpath=metadata_import.ipynb)  
Click the badge to open the notebook on [mybinder.org](https://mybinder.org) ([repo](https://github.com/PDBeurope/mmcif-metadata-import)). The first launch may take a few minutes while the environment builds (gemmi install). **Download any output files** from the notebook links before closing the tab—Binder sessions are temporary.  
**Note:** Binder only works with **public** repos. If [mmcif-metadata-import](https://github.com/PDBeurope/mmcif-metadata-import) is private, the badge link will not work; use **Run locally** instead.

**Run locally:**  
1. Install notebook dependencies:
```bash
pip install -r requirements-notebook.txt
```
2. Start Jupyter (`jupyter notebook` or `jupyter lab`), open `metadata_import.ipynb`, and run all cells. Use the widgets to upload mmCIF files, select specifications, and run import. Outputs are saved in `notebook_output/`.

## Usage

```bash
python import_metadata.py <input_file> [--xray] [--xray_serial] [--em] [--nmr] [--macromolecules] [--citation] [--authors] [--funding] [--keywords] [-o output_file] [--merge_to_file target_file] [--log]
```

### Arguments

- `input_file`: Input mmCIF file (supports `.cif` and `.cif.V[ordinal]` extensions)
- `--xray`: Optional flag to include X-ray specific categories from specs/XRAY.csv
- `--xray_serial`: Optional flag to include X-ray serial specific categories from specs/XRAY_SERIAL.csv
- `--em`: Optional flag to include electron microscopy specific categories from specs/EM.csv
- `--nmr`: Optional flag to include NMR specific categories from specs/NMR.csv
- `--macromolecules`: Optional flag to include macromolecules categories from specs/MACROMOLECULES.csv
- `--citation`: Optional flag to include citation categories from specs/CITATION.csv
- `--authors`: Optional flag to include author categories from specs/AUTHORS.csv
- `--funding`: Optional flag to include funding categories from specs/FUNDING.csv
- `--keywords`: Optional flag to include keyword categories from specs/KEYWORDS.csv
- `-o, --output`: Optional output file name (default: `[input_name]_metadata.cif`)
- `--merge_to_file`: Optional file path to merge imported metadata into (instead of creating a new file). Metadata will be added to the first data block of the target file. The output file will be named `<originalname>_merged_with_<inputfilename>` in the same directory as the target file.
- `--log`: Optional flag to generate a log file with detailed information about the import process. The log file is automatically named based on the output file (same name with `.log` extension) and placed in the same directory as the output file.

**Note**: At least one specification file must be provided.

**Merge Mode**: When `--merge_to_file` is provided, the imported metadata will be merged into the first data block of the specified file. The metadata items will be added at the end of the first data block, before any subsequent data blocks. A new file will be created with the name pattern `<originalname>_merged_with_<inputfilename>` (e.g., if merging `6qvt_nometdata.cif` with metadata from `6qvt.cif`, the output will be `6qvt_nometdata_merged_with_6qvt.cif`). The original target file is not modified. **Important**: Categories and items that already exist in the target file will not be merged to avoid overwriting existing data. These will be reported in the log file as "Categories not imported" and "Items not imported". If `--merge_to_file` is not provided, a new metadata file will be created as specified by `-o/--output`.

**Method Validation**: The script automatically detects the input file's method and validates method-specific flags. If you try to use `--xray` on an EM file, the script will warn you and skip the X-ray specification to prevent importing incompatible metadata.

### Examples

```bash
# Basic usage with method-specific files
python import_metadata.py input.cif --xray
python import_metadata.py input.cif --xray_serial
python import_metadata.py input.cif --em
python import_metadata.py input.cif --nmr

# With custom output name
python import_metadata.py input.cif --xray -o custom_output.cif

# Using only optional specification files
python import_metadata.py input.cif --macromolecules
python import_metadata.py input.cif --citation --authors
python import_metadata.py input.cif --funding --keywords

# Combine method-specific with optional files
python import_metadata.py input.cif --em --macromolecules
python import_metadata.py input.cif --xray --citation --authors
python import_metadata.py input.cif --nmr --funding --keywords

# Multiple method-specific files
python import_metadata.py input.cif --xray --xray_serial --em --nmr

# All optional categories
python import_metadata.py input.cif --macromolecules --citation --authors --funding --keywords

# Everything together
python import_metadata.py input.cif --xray --em --nmr --macromolecules --citation --authors --funding --keywords

# Method validation example (EM file with X-ray flag - X-ray will be skipped)
python import_metadata.py em_file.cif --em --xray --macromolecules
# Output: "Warning: Skipping X-ray specification - input file method (EM_MAP_ONLY) doesn't match X-ray method"

# Merge metadata into an existing file (single data block)
python import_metadata.py 6qvt.cif --xray --merge_to_file test_files/xray/6qvt_nometdata.cif

# Merge metadata into an existing file with multiple data blocks
python import_metadata.py 6qvs.cif --xray --merge_to_file test_files/xray/6qvt_nometdata_extra_datablock.cif
# Metadata will be added to the first data block, before the second data block

# Generate a log file with detailed import information (automatically named input.log)
python import_metadata.py input.cif --xray --log

# Combine merge with log file (log file automatically named based on merge output)
python import_metadata.py 6qvt.cif --xray --merge_to_file test_files/xray/6qvt_nometdata.cif --log
# Log file will be: test_files/xray/6qvt_nometdata_merged_with_6qvt.log
```

## Method Detection

The script automatically detects the source method (FROM) from the input mmCIF file based on:

- **XRAY**: `exptl.method = "X-RAY DIFFRACTION"`
- **NMR**: `exptl.method = "SOLUTION NMR"`
- **EM_MAP_ONLY**: `exptl.method = "ELECTRON MICROSCOPY"` + `database_2.database_id` contains "WWPDB" and "EMDB"
- **EM_MODEL_ONLY**: `exptl.method = "ELECTRON MICROSCOPY"` + `database_2.database_id` contains "WWPDB" and "PDB"
- **EM_MAP_MODEL**: `exptl.method = "ELECTRON MICROSCOPY"` + `database_2.database_id` contains "WWPDB", "PDB", and "EMDB"

## Specification Files

All specification CSV files are located in the `specs/` subdirectory.

The script uses simplified method-specific CSV files:

### Method-Specific Files:
- `specs/XRAY.csv` - X-ray crystallography specific categories
- `specs/XRAY_SERIAL.csv` - X-ray serial specific categories
- `specs/EM.csv` - Electron microscopy specific categories  
- `specs/NMR.csv` - Nuclear magnetic resonance specific categories

## Optional Specification Files

The script supports several optional flags that add additional categories from separate CSV files. These are merged with the method-specific specification file to provide comprehensive metadata import.

### Available Optional Files:

#### `--macromolecules` (specs/MACROMOLECULES.csv)
Contains macromolecule-related categories:
- `_entity`, `_entity_name_com`, `_entity_poly`, `_entity_poly_seq`
- `_entity_src_nat`, `_entity_src_gen`, `_pdbx_entity_src_syn`
- `_struct_ref`, `_struct_ref_seq`, `_struct_ref_seq_dif`

#### `--citation` (specs/CITATION.csv)
Contains citation-related categories:
- `_citation`, `_citation_author`

#### `--authors` (specs/AUTHORS.csv)
Contains author-related categories:
- `_pdbx_contact_author`, `_em_author_list`

#### `--funding` (specs/FUNDING.csv)
Contains funding-related categories:
- `_pdbx_audit_support`

#### `--keywords` (specs/KEYWORDS.csv)
Contains keyword-related items:
- `_struct_keywords.text`, `_struct_keywords.pdbx_keywords`, `_struct_keywords.pdbx_details`

All optional categories are merged with the method-specific specification file to provide comprehensive metadata information in the output.

## CSV Specification File Format

Each CSV specification file should contain the following columns:

- `category`: The mmCIF category name (e.g., `_pdbx_contact_author`)
- `item`: The specific item name within the category (e.g., `id`, `name_first`). Leave empty for category-level specifications.
- `should_import`: Whether to include this category/item (`Y` for yes, `N` for no)
- `type`: Either `category` (for entire category) or `item` (for specific items)

### Example CSV structure:

```csv
category,item,should_import,type
_pdbx_contact_author,,Y,category
_citation,,Y,category
_struct_keywords,text,Y,item
_struct_keywords,pdbx_keywords,Y,item
_database_2,,N,category
_struct_keywords,entry_id,N,item
```

### Annotated Example:

```csv
# Header row
category,item,should_import,type

# Include entire _pdbx_contact_author category (all items)
_pdbx_contact_author,,Y,category

# Include entire _citation category (all items)
_citation,,Y,category

# Include only specific items from _struct_keywords category
_struct_keywords,text,Y,item                    # Include _struct_keywords.text
_struct_keywords,pdbx_keywords,Y,item           # Include _struct_keywords.pdbx_keywords
_struct_keywords,entry_id,N,item                # Exclude _struct_keywords.entry_id

# Exclude entire _database_2 category (no items)
_database_2,,N,category
```

**Key Points:**
- **Empty `item` column** = entire category (use `type=category`)
- **Filled `item` column** = specific item (use `type=item`)
- **`Y`** = include this category/item
- **`N`** = exclude this category/item

## Output

The script creates a new mmCIF file containing only the specified categories and items from the input file. The output filename follows the pattern `[input_name]_metadata.cif`.

**Output Format**: The output file does not include a `data_` block declaration line at the beginning. This allows the metadata content to be easily appended to the first data block of an existing mmCIF file. The file starts directly with the metadata categories and items.

## Log File

When using the `--log` flag, a detailed log file is automatically generated with the same name as the output file but with a `.log` extension, placed in the same directory as the output file. For example:
- If output file is `input_metadata.cif`, the log file will be `input_metadata.log`
- If merge output is `test_files/xray/6qvt_nometdata_merged_with_6qvt.cif`, the log file will be `test_files/xray/6qvt_nometdata_merged_with_6qvt.log`

The log file contains:

- **Requested Categories and Items**: Lists all categories and items that were requested to be imported based on the specification files
- **Skipped Specifications**: Lists any specification files that were skipped (e.g., due to method mismatch) with the reason
- **Imported Categories and Items**: Lists all categories and items that were successfully imported
- **Categories Not Found**: Lists categories that were requested but not found in the input file
- **Items Not Found**: Lists items that were requested but not found in the input file
- **Categories Not Imported** (merge mode only): Lists categories that were not imported because they already exist in the target file
- **Items Not Imported** (merge mode only): Lists items that were not imported because they already exist in the target file
- **Summary**: Provides counts of requested vs imported categories/items, skipped specifications, categories not found, items not found, and (for merge mode) categories/items not imported

This log file is useful for debugging and understanding what metadata was imported and what was skipped.

## Features

- Supports both `.cif` and `.cif.V[ordinal]` input file extensions
- Processes only the first data block in multi-block mmCIF files
- Handles both single items and loop structures in mmCIF files
- Uses CSV format for easy specification management
- Provides detailed error messages for file reading/writing issues
- Optional log file generation for detailed import tracking

## Author & Affiliation

**Deborah Harrus** — [Protein Data Bank in Europe (PDBe)](https://www.ebi.ac.uk/pdbe)

## License

This project is licensed under the **Apache License 2.0**. See [LICENSE](LICENSE) for the full text.