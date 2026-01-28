#!/usr/bin/env python3
"""
Flask web application for the import_metadata.py script.
Provides a web interface for users who don't want to use the command line.

Author: Deborah Harrus, Protein Data Bank in Europe (PDBe)
       https://www.ebi.ac.uk/pdbe
"""

import os
import sys
import subprocess
import tempfile
import shutil
from pathlib import Path
from flask import Flask, render_template, request, send_file, flash, redirect, url_for, after_this_request
from werkzeug.utils import secure_filename

# Add parent directory to path to import import_metadata
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

# Import the import_metadata function from the script
from import_metadata import import_metadata

app = Flask(__name__)
app.secret_key = os.urandom(24)  # For flash messages
app.config['MAX_CONTENT_LENGTH'] = 100 * 1024 * 1024  # 100MB max file size
app.config['UPLOAD_FOLDER'] = tempfile.mkdtemp()

# Allowed file extensions
ALLOWED_EXTENSIONS = {'cif'}


def allowed_file(filename):
    """Check if file has an allowed extension."""
    # Handle .cif.V[ordinal] extensions
    if '.cif.V' in filename:
        return True
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS


def cleanup_temp_files(temp_dir):
    """Remove temporary directory and its contents."""
    try:
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)
    except Exception as e:
        print(f"Warning: Could not clean up temp directory {temp_dir}: {e}")


@app.route('/')
def index():
    """Render the main form page."""
    return render_template('index.html')


@app.route('/process', methods=['POST'])
def process():
    """Handle form submission and process the metadata import."""
    temp_dir = None
    try:
        # Create a temporary directory for this request
        temp_dir = tempfile.mkdtemp()
        
        # Check if input file was provided
        if 'input_file' not in request.files:
            flash('Error: No input file provided', 'error')
            return redirect(url_for('index'))
        
        input_file_obj = request.files['input_file']
        if input_file_obj.filename == '':
            flash('Error: No input file selected', 'error')
            return redirect(url_for('index'))
        
        if not allowed_file(input_file_obj.filename):
            flash('Error: Invalid file type. Only .cif files are allowed', 'error')
            return redirect(url_for('index'))
        
        # Save input file
        input_filename = secure_filename(input_file_obj.filename)
        input_path = os.path.join(temp_dir, input_filename)
        input_file_obj.save(input_path)
        
        # Handle merge_to_file if provided
        merge_to_path = None
        merge_output_path = None
        if 'merge_to_file' in request.files:
            merge_file_obj = request.files['merge_to_file']
            if merge_file_obj.filename != '':
                if not allowed_file(merge_file_obj.filename):
                    flash('Error: Invalid merge file type. Only .cif files are allowed', 'error')
                    cleanup_temp_files(temp_dir)
                    return redirect(url_for('index'))
                
                merge_filename = secure_filename(merge_file_obj.filename)
                merge_to_path = os.path.join(temp_dir, merge_filename)
                merge_file_obj.save(merge_to_path)
        
        # Collect specification files based on checkboxes
        spec_files = []
        skipped_specs = []
        
        # Get the base directory (parent of web_app)
        base_dir = Path(__file__).parent.parent
        
        # Method-specific specs (need validation)
        method_flags = {
            'xray': str(base_dir / 'specs' / 'XRAY.csv'),
            'xray_serial': str(base_dir / 'specs' / 'XRAY_SERIAL.csv'),
            'em': str(base_dir / 'specs' / 'EM.csv'),
            'nmr': str(base_dir / 'specs' / 'NMR.csv')
        }
        
        # Additional specs (no validation needed)
        additional_specs = {
            'macromolecules': str(base_dir / 'specs' / 'MACROMOLECULES.csv'),
            'citation': str(base_dir / 'specs' / 'CITATION.csv'),
            'authors': str(base_dir / 'specs' / 'AUTHORS.csv'),
            'funding': str(base_dir / 'specs' / 'FUNDING.csv'),
            'keywords': str(base_dir / 'specs' / 'KEYWORDS.csv')
        }
        
        # Check if at least one spec is selected
        has_method_spec = any(request.form.get(flag) == 'on' for flag in method_flags.keys())
        has_additional_spec = any(request.form.get(flag) == 'on' for flag in additional_specs.keys())
        
        if not has_method_spec and not has_additional_spec:
            flash('Error: At least one specification option must be selected', 'error')
            cleanup_temp_files(temp_dir)
            return redirect(url_for('index'))
        
        # Detect method from input file for validation
        detected_method = None
        if has_method_spec:
            try:
                import gemmi
                doc = gemmi.cif.read(input_path)
                from import_metadata import detect_method_from_input
                detected_method = detect_method_from_input(doc)
            except Exception as e:
                flash(f'Error detecting method from input file: {e}', 'error')
                cleanup_temp_files(temp_dir)
                return redirect(url_for('index'))
        
        # Process method-specific specs with validation
        for flag, spec_file in method_flags.items():
            if request.form.get(flag) == 'on':
                spec_path = Path(spec_file)
                if not spec_path.exists():
                    flash(f'Error: Specification file {spec_file} does not exist', 'error')
                    cleanup_temp_files(temp_dir)
                    return redirect(url_for('index'))
                
                # Validate method compatibility
                if detected_method:
                    expected_methods = {
                        'xray': ['XRAY'],
                        'xray_serial': ['XRAY'],
                        'em': ['EM_MAP_MODEL', 'EM_MAP_ONLY', 'EM_MODEL_ONLY'],
                        'nmr': ['NMR']
                    }
                    if detected_method not in expected_methods.get(flag, []):
                        reason = f"Input file method ({detected_method}) doesn't match {flag} method"
                        skipped_specs.append((spec_file, reason))
                        continue
                
                spec_files.append(spec_file)
        
        # Process additional specs
        for flag, spec_file in additional_specs.items():
            if request.form.get(flag) == 'on':
                spec_path = Path(spec_file)
                if not spec_path.exists():
                    flash(f'Error: Specification file {spec_file} does not exist', 'error')
                    cleanup_temp_files(temp_dir)
                    return redirect(url_for('index'))
                spec_files.append(spec_file)
        
        # Determine output file path
        if merge_to_path:
            # When merging, output is automatically generated
            output_path = None
            # Generate merge output filename
            merge_name = os.path.basename(merge_to_path)
            input_name = os.path.basename(input_path)
            
            # Extract stems
            if '.cif.V' in merge_name:
                merge_stem = merge_name.split('.cif.V')[0]
                merge_suffix = '.cif.V' + merge_name.split('.cif.V')[1]
            elif merge_name.endswith('.cif'):
                merge_stem = merge_name[:-4]
                merge_suffix = '.cif'
            else:
                merge_stem = Path(merge_name).stem
                merge_suffix = Path(merge_name).suffix
            
            if '.cif.V' in input_name:
                input_stem = input_name.split('.cif.V')[0]
            elif input_name.endswith('.cif'):
                input_stem = input_name[:-4]
            else:
                input_stem = Path(input_name).stem
            
            merge_output_path = os.path.join(temp_dir, f"{merge_stem}_merged_with_{input_stem}{merge_suffix}")
        else:
            # Regular output
            if request.form.get('output_filename'):
                output_filename = secure_filename(request.form.get('output_filename'))
                if not output_filename.endswith('.cif'):
                    output_filename += '.cif'
                output_path = os.path.join(temp_dir, output_filename)
            else:
                # Generate default output name
                input_stem = Path(input_filename).stem
                if input_stem.endswith('.cif'):
                    input_stem = input_stem[:-4]
                output_path = os.path.join(temp_dir, f"{input_stem}_metadata.cif")
            merge_output_path = None
        
        # Determine log file path
        log_path = None
        if request.form.get('log') == 'on':
            if merge_output_path:
                log_base = Path(merge_output_path).stem
            elif output_path:
                log_base = Path(output_path).stem
            else:
                log_base = Path(input_filename).stem
            log_path = os.path.join(temp_dir, f"{log_base}.log")
        
        # Call the import_metadata function
        success = import_metadata(
            input_path,
            spec_files,
            output_path,
            merge_to_path,
            merge_output_path,
            log_path,
            skipped_specs if skipped_specs else None
        )
        
        if not success:
            flash('Error: Metadata import failed. Check the log file for details.', 'error')
            cleanup_temp_files(temp_dir)
            return redirect(url_for('index'))
        
        # Prepare files for download
        files_to_download = []
        
        if merge_output_path and os.path.exists(merge_output_path):
            files_to_download.append(('merged_output', merge_output_path, os.path.basename(merge_output_path)))
        elif output_path and os.path.exists(output_path):
            files_to_download.append(('output', output_path, os.path.basename(output_path)))
        
        if log_path and os.path.exists(log_path):
            files_to_download.append(('log', log_path, os.path.basename(log_path)))
        
        if not files_to_download:
            flash('Error: No output files were generated', 'error')
            cleanup_temp_files(temp_dir)
            return redirect(url_for('index'))
        
        # Store temp_dir in session or pass to download handler
        # For simplicity, we'll create a zip file with all results
        import zipfile
        zip_path = os.path.join(temp_dir, 'results.zip')
        with zipfile.ZipFile(zip_path, 'w') as zipf:
            for file_type, file_path, filename in files_to_download:
                zipf.write(file_path, filename)
        
        # Clean up temp directory after response is sent
        @after_this_request
        def cleanup(response):
            cleanup_temp_files(temp_dir)
            return response
        
        # Return the zip file
        return send_file(
            zip_path,
            as_attachment=True,
            download_name='metadata_import_results.zip',
            mimetype='application/zip'
        )
        
    except Exception as e:
        flash(f'Error processing request: {str(e)}', 'error')
        if temp_dir:
            cleanup_temp_files(temp_dir)
        return redirect(url_for('index'))


if __name__ == '__main__':
    # Clean up temp directory on exit
    import atexit
    atexit.register(lambda: cleanup_temp_files(app.config['UPLOAD_FOLDER']))
    
    app.run(debug=True, host='0.0.0.0', port=5000)

