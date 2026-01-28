# Web Interface for mmCIF Metadata Importer

This web interface provides a user-friendly way to use the `import_metadata.py` script (mmCIF metadata importer) without needing to use the command line.

**Protein Data Bank in Europe (PDBe)** · [pdbe.org](https://www.ebi.ac.uk/pdbe)

## Features

- **Easy-to-use web form** with checkboxes for all options
- **File upload** for input files and merge target files
- **Method validation** - automatically validates that your input file matches the selected method
- **Automatic file downloads** - results are packaged in a ZIP file
- **Beautiful, modern UI** - responsive design that works on all devices

## Installation

1. Install the required dependencies:
   ```bash
   cd web_app
   pip install -r requirements.txt
   ```

2. Make sure all specification CSV files are in the `../specs/` directory (relative to web_app):
   - `XRAY.csv`
   - `XRAY_SERIAL.csv`
   - `EM.csv`
   - `NMR.csv`
   - `MACROMOLECULES.csv`
   - `CITATION.csv`
   - `AUTHORS.csv`
   - `FUNDING.csv`
   - `KEYWORDS.csv`

## Running the Web Application

1. Navigate to the web_app directory:
   ```bash
   cd web_app
   ```

2. Start the Flask server:
   ```bash
   python app.py
   ```

3. Open your web browser and navigate to:
   ```
   http://localhost:5000
   ```

4. The web interface will be available at that address.

**Note:** The web app needs to be run from the `web_app` directory, but it will automatically find the `import_metadata.py` script and `specs/` directory in the parent directory.

## Usage

1. **Upload Input File**: Select your input mmCIF file (required)
   - Supports `.cif` and `.cif.V[ordinal]` files

2. **Upload Merge File** (optional): If you want to merge metadata into an existing file

3. **Select Specifications**: 
   - Choose at least one method-specific specification (X-ray, EM, or NMR)
   - Optionally select additional specifications (macromolecules, citation, authors, funding, keywords)

4. **Output Options**:
   - Optionally specify a custom output filename
   - Check "Generate Log File" to get detailed import information

5. **Submit**: Click "Import Metadata" to process your file

6. **Download Results**: The results will be automatically downloaded as a ZIP file containing:
   - The output file (merged or imported metadata)
   - The log file (if requested)

## Notes

- The web application validates that your input file method matches the selected method-specific specification
- If validation fails, you'll see an error message explaining why
- All processing happens server-side, so large files may take a moment to process
- Maximum file size is 100MB (configurable in `app.py`)

## Deployment

For production deployment, consider:

1. Using a production WSGI server like Gunicorn:
   ```bash
   pip install gunicorn
   gunicorn -w 4 -b 0.0.0.0:5000 app:app
   ```

2. Setting up proper security (HTTPS, authentication, etc.)

3. Configuring proper file cleanup and storage

4. Adding rate limiting to prevent abuse

## Troubleshooting

- **"No specification file found"**: Make sure all CSV files are in the `specs/` directory
- **"Method validation failed"**: Your input file's method doesn't match the selected specification
- **"File too large"**: The default limit is 100MB. You can increase this in `app.py` by modifying `MAX_CONTENT_LENGTH`

