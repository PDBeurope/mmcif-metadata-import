#!/usr/bin/env python3
"""
mmCIF Metadata Importer: import specific categories and items from mmCIF files using gemmi.
Creates metadata-only files or merges metadata into existing models based on CSV specifications.

Author: Deborah Harrus, Protein Data Bank in Europe (PDBe)
       https://www.ebi.ac.uk/pdbe
"""

import sys
import os
import re
import argparse
import csv
from pathlib import Path
import gemmi

PROJECT_ROOT = Path(__file__).resolve().parent


def resolve_spec_path(spec_path):
    """
    Resolve a specification file path from either:
    - current working directory (repo usage), or
    - installed package location (PyPI usage).
    """
    candidate = Path(spec_path)
    if candidate.exists():
        return candidate

    package_candidate = PROJECT_ROOT / candidate
    if package_candidate.exists():
        return package_candidate

    return candidate


def parse_specification_file(spec_file_path):
    """
    Parse the CSV specification file to determine which categories and items to include/exclude.
    
    Returns:
        tuple: (included_categories, included_items, excluded_categories, excluded_items)
    """
    included_categories = set()
    included_items = set()
    excluded_categories = set()
    excluded_items = set()
    
    spec_file_path = resolve_spec_path(spec_file_path)
    with open(spec_file_path, 'r', newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            category = row['category']
            item = row['item']
            should_import = row['should_import'].upper() == 'Y'
            item_type = row['type']
            
            if item_type == 'category':
                if should_import:
                    included_categories.add(category)
                else:
                    excluded_categories.add(category)
            elif item_type == 'item':
                if item:  # Only process if item is not empty
                    full_item_name = f"{category}.{item}"
                    if should_import:
                        included_items.add(full_item_name)
                    else:
                        excluded_items.add(full_item_name)
    
    return included_categories, included_items, excluded_categories, excluded_items


def parse_multiple_specification_files(spec_files):
    """
    Parse multiple CSV specification files and merge their specifications.
    
    Args:
        spec_files (list): List of CSV file paths to parse
    
    Returns:
        tuple: (included_categories, included_items, excluded_categories, excluded_items)
    """
    all_included_categories = set()
    all_included_items = set()
    all_excluded_categories = set()
    all_excluded_items = set()
    
    for spec_file in spec_files:
        included_categories, included_items, excluded_categories, excluded_items = parse_specification_file(spec_file)
        
        all_included_categories.update(included_categories)
        all_included_items.update(included_items)
        all_excluded_categories.update(excluded_categories)
        all_excluded_items.update(excluded_items)
    
    return all_included_categories, all_included_items, all_excluded_categories, all_excluded_items


def should_include_block(block_name, included_categories, excluded_categories):
    """Determine if a block should be included based on category rules."""
    # For mmCIF files, we typically have one data block, so we need to check
    # if any items in the block should be included
    if block_name in excluded_categories:
        return False
    if block_name in included_categories:
        return True
    # If not explicitly included or excluded, we need to check individual items
    return None  # Special value to indicate we need to check items


def get_category_from_item(item_name):
    """Extract category name from item name (e.g., '_pdbx_contact_author.id' -> '_pdbx_contact_author')"""
    if '.' in item_name:
        return item_name.split('.')[0]
    return item_name


def detect_method_from_input(doc):
    """
    Detect the method FROM the input mmCIF file based on exptl.method and database_2.database_id.
    
    Returns:
        str: The detected method (EM_MAP_MODEL, EM_MAP_ONLY, EM_MODEL_ONLY, NMR, XRAY)
    """
    if len(doc) == 0:
        raise ValueError("No data blocks found in input file")
    
    block = doc[0]  # Get the first data block
    
    # Get exptl.method value
    exptl_method = None
    for item in block:
        if item.pair is not None and item.pair[0] == "_exptl.method":
            exptl_method = item.pair[1]
            break
    
    if exptl_method is None:
        raise ValueError("Could not find _exptl.method in input file")
    
    # Remove quotes if present (both single and double quotes)
    exptl_method = exptl_method.strip('"').strip("'")
    
    # Check for XRAY
    if exptl_method == "X-RAY DIFFRACTION":
        return "XRAY"
    
    # Check for NMR
    if exptl_method == "SOLUTION NMR":
        return "NMR"
    
    # Check for EM methods
    if exptl_method == "ELECTRON MICROSCOPY":
        # Get database_2.database_id values
        database_ids = []
        for item in block:
            if item.loop is not None:
                loop = item.loop
                if "_database_2.database_id" in [tag for tag in loop.tags]:
                    # Find the index of database_id in the loop
                    db_id_index = None
                    for i, tag in enumerate(loop.tags):
                        if tag == "_database_2.database_id":
                            db_id_index = i
                            break
                    
                    if db_id_index is not None:
                        # Extract database_id values
                        num_tags = len(loop.tags)
                        for i in range(db_id_index, len(loop.values), num_tags):
                            if i < len(loop.values):
                                database_ids.append(loop.values[i])
        
        # Determine EM method based on database_ids
        if "WWPDB" in database_ids and "EMDB" in database_ids and "PDB" in database_ids:
            return "EM_MAP_MODEL"
        elif "WWPDB" in database_ids and "EMDB" in database_ids:
            return "EM_MAP_ONLY"
        elif "WWPDB" in database_ids and "PDB" in database_ids:
            return "EM_MODEL_ONLY"
        else:
            raise ValueError(f"Unexpected database_2.database_id values for EM method: {database_ids}")
    
    raise ValueError(f"Unknown exptl.method: {exptl_method}")


def get_spec_file_path(from_method, to_method):
    """
    Get the path to the specification CSV file based on FROM and TO methods.
    
    Args:
        from_method (str): The source method
        to_method (str): The target method
    
    Returns:
        str: Path to the specification CSV file
    """
    spec_filename = f"FROM_{from_method}_TO_{to_method}.csv"
    return spec_filename


def should_include_item(item_name, included_items, excluded_items):
    """Determine if an item should be included based on item rules."""
    if item_name in excluded_items:
        return False
    if item_name in included_items:
        return True
    # If not explicitly included or excluded, exclude by default
    return False


def merge_metadata_to_file(metadata_block, merge_to_file, output_file):
    """
    Merge imported metadata into the first data block of a target file.
    Works with text files directly.
    Skips categories/items that already exist in the target file.

    Args:
        metadata_block (gemmi.cif.Block): Block containing imported metadata
        merge_to_file (str): Path to target file to merge metadata into
        output_file (str): Path to output file where merged result will be written
    
    Returns:
        tuple: (success (bool), already_present_categories (set), already_present_items (set))
    """
    try:
        # Read the target file to check what already exists
        target_doc = gemmi.cif.read(merge_to_file)
        if len(target_doc) == 0:
            print("Error: No data blocks found in merge target file")
            return False, set(), set()
        
        target_block = target_doc[0]
        
        # Track what categories and items already exist in target file
        existing_categories = set()
        existing_items = set()
        
        for item in target_block:
            if item.pair is not None:
                item_name = item.pair[0]
                category = get_category_from_item(item_name)
                existing_items.add(item_name)
                existing_categories.add(category)
            elif item.loop is not None:
                for tag in item.loop.tags:
                    category = get_category_from_item(tag)
                    existing_items.add(tag)
                    existing_categories.add(category)
        
        # Filter metadata_block to exclude items/categories that already exist
        filtered_metadata_block = gemmi.cif.Block(metadata_block.name)
        already_present_categories = set()
        already_present_items = set()
        
        for item in metadata_block:
            if item.pair is not None:
                item_name = item.pair[0]
                category = get_category_from_item(item_name)
                
                if item_name in existing_items:
                    already_present_items.add(item_name)
                    already_present_categories.add(category)
                    continue  # Skip this item
                else:
                    filtered_metadata_block.set_pair(item.pair[0], item.pair[1])
                    
            elif item.loop is not None:
                # Check if any tags in the loop already exist
                loop = item.loop
                loop_tags = [tag for tag in loop.tags]
                
                # Check which tags should be included (not already present)
                tags_to_include = []
                for tag in loop_tags:
                    category = get_category_from_item(tag)
                    if tag in existing_items:
                        already_present_items.add(tag)
                        already_present_categories.add(category)
                    else:
                        tags_to_include.append(tag)
                
                # If we have tags to include, create a filtered loop
                if tags_to_include:
                    # For simplicity, if any tag is new, include the whole loop
                    # This is a simplified approach
                    filtered_metadata_block.add_item(item)
                else:
                    # All tags already exist, mark categories
                    for tag in loop_tags:
                        category = get_category_from_item(tag)
                        already_present_categories.add(category)
        
        # Write the filtered metadata block to a temporary string
        metadata_doc = gemmi.cif.Document()
        metadata_doc.add_copied_block(filtered_metadata_block)
        metadata_text = metadata_doc.as_string()
        
        # Read the target file as text
        with open(merge_to_file, 'r', encoding='utf-8') as f:
            target_lines = f.readlines()
    except Exception as e:
        print(f"Error reading files: {e}")
        return False, set(), set()
    
    # Find where the first data block ends (before the second "data_" line)
    first_data_block_end = None
    data_block_count = 0
    
    for i, line in enumerate(target_lines):
        if line.strip().startswith('data_'):
            data_block_count += 1
            if data_block_count == 2:
                # Found the second data block - insert before this line
                first_data_block_end = i
                break
    
    # If no second data block found, insert at the end of the file
    if first_data_block_end is None:
        first_data_block_end = len(target_lines)
    
    # Extract metadata content (skip the first "data_" line and any leading whitespace)
    metadata_lines = metadata_text.split('\n')
    # Find the first "data_" line in metadata
    metadata_start = 0
    for i, line in enumerate(metadata_lines):
        if line.strip().startswith('data_'):
            metadata_start = i + 1
            break
    
    # Get the actual metadata content (everything after the data_ line)
    metadata_content = '\n'.join(metadata_lines[metadata_start:])
    # Remove trailing newlines but keep the content
    metadata_content = metadata_content.rstrip()
    
    # Build the merged content
    # Take everything up to (but not including) the second data block
    merged_lines = target_lines[:first_data_block_end]
    
    # Add a newline if the last line doesn't end with one
    if merged_lines and not merged_lines[-1].endswith('\n'):
        merged_lines[-1] = merged_lines[-1] + '\n'
    
    # Add the metadata content
    if metadata_content:
        merged_lines.append(metadata_content)
        merged_lines.append('\n')
    
    # Add the rest of the file (second data block and beyond)
    merged_lines.extend(target_lines[first_data_block_end:])
    
    # Write the merged content to the output file
    try:
        with open(output_file, 'w', encoding='utf-8') as f:
            f.writelines(merged_lines)
        print(f"Successfully merged metadata into: {output_file}")
        return True, already_present_categories, already_present_items
    except Exception as e:
        print(f"Error writing merged file {output_file}: {e}")
        return False, set(), set()


def write_log_file(log_file, input_file, spec_files, requested_categories, requested_items, 
                   imported_categories, imported_items, skipped_specs, not_found_items, not_found_categories,
                   not_imported_categories=None, not_imported_items=None):
    """
    Write a log file with detailed information about the import process.
    
    Args:
        log_file (str): Path to log file
        input_file (str): Path to input file
        spec_files (list): List of specification files used
        requested_categories (set): Categories requested to be imported
        requested_items (set): Items requested to be imported
        imported_categories (set): Categories actually imported
        imported_items (set): Items actually imported
        skipped_specs (list): List of tuples (spec_file, reason) for skipped specifications
        not_found_items (set): Items requested but not found in input file
        not_found_categories (set): Categories requested but not found in input file
        not_imported_categories (set, optional): Categories not imported because already present in target file
        not_imported_items (set, optional): Items not imported because already present in target file
    """
    try:
        with open(log_file, 'w', encoding='utf-8') as f:
            f.write("=" * 80 + "\n")
            f.write("Metadata Import Log\n")
            f.write("=" * 80 + "\n\n")
            
            f.write(f"Input file: {input_file}\n")
            f.write(f"Specification files used: {', '.join(spec_files)}\n\n")
            
            f.write("-" * 80 + "\n")
            f.write("REQUESTED CATEGORIES AND ITEMS\n")
            f.write("-" * 80 + "\n\n")
            
            if requested_categories:
                f.write("Categories requested to be imported:\n")
                for cat in sorted(requested_categories):
                    f.write(f"  - {cat}\n")
                f.write("\n")
            else:
                f.write("No categories requested.\n\n")
            
            if requested_items:
                f.write("Items requested to be imported:\n")
                for item in sorted(requested_items):
                    f.write(f"  - {item}\n")
                f.write("\n")
            else:
                f.write("No specific items requested.\n\n")
            
            if skipped_specs:
                f.write("-" * 80 + "\n")
                f.write("SKIPPED SPECIFICATIONS\n")
                f.write("-" * 80 + "\n\n")
                for spec_file, reason in skipped_specs:
                    f.write(f"Specification file: {spec_file}\n")
                    f.write(f"Reason: {reason}\n\n")
            
            f.write("-" * 80 + "\n")
            f.write("IMPORTED CATEGORIES AND ITEMS\n")
            f.write("-" * 80 + "\n\n")
            
            if imported_categories:
                f.write("Categories successfully imported:\n")
                for cat in sorted(imported_categories):
                    f.write(f"  - {cat}\n")
                f.write("\n")
            else:
                f.write("No categories imported.\n\n")
            
            if imported_items:
                f.write("Items successfully imported:\n")
                for item in sorted(imported_items):
                    f.write(f"  - {item}\n")
                f.write("\n")
            else:
                f.write("No items imported.\n\n")
            
            if not_found_categories:
                f.write("-" * 80 + "\n")
                f.write("CATEGORIES NOT FOUND IN INPUT FILE\n")
                f.write("-" * 80 + "\n\n")
                f.write("The following categories were requested but not found in the input file:\n")
                for category in sorted(not_found_categories):
                    f.write(f"  - {category}\n")
                f.write("\n")
            
            if not_found_items:
                f.write("-" * 80 + "\n")
                f.write("ITEMS NOT FOUND IN INPUT FILE\n")
                f.write("-" * 80 + "\n\n")
                f.write("The following items were requested but not found in the input file:\n")
                for item in sorted(not_found_items):
                    f.write(f"  - {item}\n")
                f.write("\n")
            
            if not_imported_categories:
                f.write("-" * 80 + "\n")
                f.write("CATEGORIES NOT IMPORTED\n")
                f.write("-" * 80 + "\n\n")
                f.write("The following categories were not imported because they already exist in the target file:\n")
                for category in sorted(not_imported_categories):
                    f.write(f"  - {category}\n")
                f.write("\n")
            
            if not_imported_items:
                f.write("-" * 80 + "\n")
                f.write("ITEMS NOT IMPORTED\n")
                f.write("-" * 80 + "\n\n")
                f.write("The following items were not imported because they already exist in the target file:\n")
                for item in sorted(not_imported_items):
                    f.write(f"  - {item}\n")
                f.write("\n")
            
            # Summary
            f.write("-" * 80 + "\n")
            f.write("SUMMARY\n")
            f.write("-" * 80 + "\n\n")
            f.write(f"Categories requested: {len(requested_categories)}\n")
            f.write(f"Categories imported: {len(imported_categories)}\n")
            f.write(f"Categories not found: {len(not_found_categories)}\n")
            if not_imported_categories:
                f.write(f"Categories not imported (already present): {len(not_imported_categories)}\n")
            f.write(f"Items requested: {len(requested_items)}\n")
            f.write(f"Items imported: {len(imported_items)}\n")
            f.write(f"Items not found: {len(not_found_items)}\n")
            if not_imported_items:
                f.write(f"Items not imported (already present): {len(not_imported_items)}\n")
            f.write(f"Specifications skipped: {len(skipped_specs)}\n")
        
        print(f"Log file written to: {log_file}")
        return True
    except Exception as e:
        print(f"Error writing log file {log_file}: {e}")
        return False


def import_metadata(input_file, spec_files, output_file, merge_to_file=None, merge_output_file=None, 
                   log_file=None, skipped_specs=None):
    """
    Import metadata from input mmCIF file based on specification.
    
    Args:
        input_file (str): Path to input mmCIF file
        spec_files (str or list): Path to specification file(s)
        output_file (str): Path to output metadata file (used if merge_to_file is None)
        merge_to_file (str, optional): Path to target file to merge metadata into
        merge_output_file (str, optional): Path to output file for merged result
        log_file (str, optional): Path to log file to write
        skipped_specs (list, optional): List of tuples (spec_file, reason) for skipped specifications
    """
    # Handle both single file and multiple files
    if isinstance(spec_files, str):
        spec_files = [spec_files]
    
    # Parse specification(s)
    if len(spec_files) == 1:
        included_categories, included_items, excluded_categories, excluded_items = parse_specification_file(spec_files[0])
    else:
        included_categories, included_items, excluded_categories, excluded_items = parse_multiple_specification_files(spec_files)
    
    # Read input file
    try:
        doc = gemmi.cif.read(input_file)
    except Exception as e:
        print(f"Error reading input file {input_file}: {e}")
        return False
    
    # Create new document for output
    output_doc = gemmi.cif.Document()
    
    # Process only the first data block
    if len(doc) == 0:
        print("No data blocks found in input file")
        return False
    
    block = doc[0]  # Get the first data block
    block_name = block.name
    print(f"Processing data block: {block_name}")
    
    # Create new block for output
    new_block = gemmi.cif.Block(block_name)
    has_items = False
    
    # Track what was imported for logging
    imported_categories = set()
    imported_items = set()
    all_items_in_file = set()  # Track all items found in the file
    all_categories_in_file = set()  # Track all categories found in the file
    
    # Process each item in the data block
    for item in block:
        if item.pair is not None:
            # Single item
            item_name = item.pair[0]  # The item name already includes the category prefix
            category = get_category_from_item(item_name)
            all_items_in_file.add(item_name)
            all_categories_in_file.add(category)
            
            # Check if this category should be included/excluded
            if category in excluded_categories:
                continue  # Skip this item
            elif category in included_categories:
                # Include entire category, so include this item
                new_block.set_pair(item.pair[0], item.pair[1])
                imported_categories.add(category)
                imported_items.add(item_name)
                has_items = True
            else:
                # Check if this specific item should be included
                if should_include_item(item_name, included_items, excluded_items):
                    new_block.set_pair(item.pair[0], item.pair[1])
                    imported_items.add(item_name)
                    has_items = True
                    
        elif item.loop is not None:
            # Loop items - check if any tags in this loop should be included
            loop = item.loop
            loop_tags = [tag for tag in loop.tags]
            
            # Track all items in loops
            for tag in loop_tags:
                all_items_in_file.add(tag)
                tag_category = get_category_from_item(tag)
                all_categories_in_file.add(tag_category)
            
            # Check if any tags in this loop should be included
            included_tags = []
            for tag in loop_tags:
                category = get_category_from_item(tag)
                
                # Check if this category should be included/excluded
                if category in excluded_categories:
                    continue  # Skip this tag
                elif category in included_categories:
                    # Include entire category, so include this tag
                    included_tags.append(tag)
                    imported_categories.add(category)
                    imported_items.add(tag)
                else:
                    # Check if this specific item should be included
                    if should_include_item(tag, included_items, excluded_items):
                        included_tags.append(tag)
                        imported_items.add(tag)
            
            # If we have included tags, copy the entire loop (simpler approach)
            if included_tags:
                # For now, include the entire loop if any tag should be included
                # This is a simplified approach - in a more sophisticated version,
                # we would filter the loop data
                new_block.add_item(item)
                has_items = True
    
    if not has_items:
        print("No items to include in output")
        # Still write log if requested
        if log_file:
            requested_categories = included_categories.copy()
            requested_items = included_items.copy()
            not_found_items = requested_items - all_items_in_file
            not_found_categories = requested_categories - all_categories_in_file
            write_log_file(log_file, input_file, spec_files, requested_categories, requested_items,
                          imported_categories, imported_items, skipped_specs or [], not_found_items, not_found_categories,
                          None, None)
        return False
    
    # If merge_to_file is provided, merge into that file
    not_imported_categories = set()
    not_imported_items = set()
    if merge_to_file:
        result, not_imported_categories, not_imported_items = merge_metadata_to_file(new_block, merge_to_file, merge_output_file)
    else:
        # Otherwise, create a new metadata file
        output_doc.add_copied_block(new_block)
        print(f"Added items to output")
        
        # Write output file
        try:
            output_doc.write_file(output_file)
            
            # Remove the first "data_" line from the output file
            with open(output_file, 'r', encoding='utf-8') as f:
                lines = f.readlines()
            
            # Find and remove the first line that starts with "data_"
            filtered_lines = []
            data_line_removed = False
            for line in lines:
                if not data_line_removed and line.strip().startswith('data_'):
                    data_line_removed = True
                    continue  # Skip this line
                filtered_lines.append(line)
            
            # Write the file back without the data_ line
            with open(output_file, 'w', encoding='utf-8') as f:
                f.writelines(filtered_lines)
            
            print(f"Successfully created metadata file: {output_file}")
            result = True
        except Exception as e:
            print(f"Error writing output file {output_file}: {e}")
            result = False
    
    # Write log file if requested
    if log_file and result:
        requested_categories = included_categories.copy()
        requested_items = included_items.copy()
        not_found_items = requested_items - all_items_in_file
        not_found_categories = requested_categories - all_categories_in_file
        
        # Adjust imported counts for merge mode: exclude what was not imported (already present)
        if merge_to_file:
            # Only count what was actually merged (not skipped)
            actually_imported_categories = imported_categories - not_imported_categories
            actually_imported_items = imported_items - not_imported_items
        else:
            # For non-merge mode, use the original counts
            actually_imported_categories = imported_categories
            actually_imported_items = imported_items
        
        write_log_file(log_file, input_file, spec_files, requested_categories, requested_items,
                      actually_imported_categories, actually_imported_items, skipped_specs or [], not_found_items, not_found_categories,
                      not_imported_categories if merge_to_file else None, not_imported_items if merge_to_file else None)
    
    return result


def main():
    parser = argparse.ArgumentParser(description="mmCIF Metadata Importer: import metadata from mmCIF files using gemmi")
    parser.add_argument("input_file", help="Input mmCIF file (.cif or .cif.V[ordinal])")
    parser.add_argument("--xray", action="store_true", 
                       help="Include X-ray specific categories from XRAY.csv")
    parser.add_argument("--xray_serial", action="store_true", 
                       help="Include X-ray serial specific categories from XRAY_SERIAL.csv")
    parser.add_argument("--em", action="store_true", 
                       help="Include electron microscopy specific categories from EM.csv")
    parser.add_argument("--nmr", action="store_true", 
                       help="Include NMR specific categories from NMR.csv")
    parser.add_argument("--macromolecules", action="store_true", 
                       help="Include macromolecules categories from MACROMOLECULES.csv")
    parser.add_argument("--citation", action="store_true", 
                       help="Include citation categories from CITATION.csv")
    parser.add_argument("--authors", action="store_true", 
                       help="Include author categories from AUTHORS.csv")
    parser.add_argument("--funding", action="store_true", 
                       help="Include funding categories from FUNDING.csv")
    parser.add_argument("--keywords", action="store_true", 
                       help="Include keyword categories from KEYWORDS.csv")
    parser.add_argument("-o", "--output", help="Output file name (default: [input_name]_metadata.cif)")
    parser.add_argument("--merge_to_file", help="File to merge imported metadata into (instead of creating a new file)")
    parser.add_argument("--log", action="store_true", help="Generate a log file with detailed import information (automatically named based on output file)")
    
    args = parser.parse_args()
    
    # Validate input file
    input_path = Path(args.input_file)
    if not input_path.exists():
        print(f"Error: Input file {args.input_file} does not exist")
        sys.exit(1)
    
    # Validate merge_to_file if provided
    merge_output_file = None
    if args.merge_to_file:
        merge_path = Path(args.merge_to_file)
        if not merge_path.exists():
            print(f"Error: Merge target file {args.merge_to_file} does not exist")
            sys.exit(1)
        
        # Generate output filename: <originalname>_merged_with_<inputfilename>
        merge_dir = merge_path.parent
        merge_name = merge_path.name
        
        # Extract stem and extension from merge_to_file
        # Handle .cif.V[ordinal] extension
        if '.cif.V' in merge_name:
            parts = merge_name.split('.cif.V')
            merge_stem = parts[0]
            merge_suffix = '.cif.V' + parts[-1]
        elif merge_name.endswith('.cif'):
            merge_stem = merge_name[:-4]  # Remove .cif
            merge_suffix = '.cif'
        else:
            # Fallback: use Path methods
            merge_stem = merge_path.stem
            merge_suffix = merge_path.suffix
        
        # Extract stem from input_file
        input_name = input_path.name
        if '.cif.V' in input_name:
            parts = input_name.split('.cif.V')
            input_stem = parts[0]
        elif input_name.endswith('.cif'):
            input_stem = input_name[:-4]  # Remove .cif
        else:
            input_stem = input_path.stem
        
        merge_output_file = str(merge_dir / f"{merge_stem}_merged_with_{input_stem}{merge_suffix}")
    
    # Determine output file name (only used if merge_to_file is not provided)
    if args.merge_to_file:
        output_file = None  # Not used when merging
    elif args.output:
        output_file = args.output
    else:
        # Remove .cif.V[ordinal] or .cif extension and add _metadata.cif
        stem = input_path.stem
        # If the stem ends with .cif, remove it
        if stem.endswith('.cif'):
            stem = stem[:-4]  # Remove .cif
        output_file = f"{stem}_metadata.cif"
    
    # Generate log file name if requested
    log_file = None
    if args.log:
        # Determine which output file to base the log name on
        if args.merge_to_file and merge_output_file:
            # Use merge output file
            log_base = Path(merge_output_file)
        elif output_file:
            # Use regular output file
            log_base = Path(output_file)
        else:
            # Fallback: use input file name
            log_base = input_path
        
        # Generate log file name: same directory, same stem, .log extension
        log_file = str(log_base.parent / f"{log_base.stem}.log")
    
    # Prepare list of specification files
    spec_files = []
    skipped_specs = []  # Track skipped specifications for logging
    
    # Detect input file method for validation
    detected_method = None
    if any([args.xray, args.xray_serial, args.em, args.nmr]):
        # Only read the file if method-specific flags are used
        try:
            doc = gemmi.cif.read(args.input_file)
            detected_method = detect_method_from_input(doc)
            print(f"Detected input file method: {detected_method}")
        except Exception as e:
            print(f"Error reading input file {args.input_file}: {e}")
            sys.exit(1)
    
    # Handle method-specific specification files with validation
    method_specs = [
        (args.xray, "specs/XRAY.csv", "X-ray", "XRAY"),
        (args.xray_serial, "specs/XRAY_SERIAL.csv", "X-ray serial", "XRAY"),
        (args.em, "specs/EM.csv", "electron microscopy", ["EM_MAP_MODEL", "EM_MAP_ONLY", "EM_MODEL_ONLY"]),
        (args.nmr, "specs/NMR.csv", "NMR", "NMR")
    ]
    
    for flag, filename, description, expected_methods in method_specs:
        if flag:
            # Validate method compatibility
            if detected_method is not None:
                if isinstance(expected_methods, str):
                    expected_methods = [expected_methods]
                
                if detected_method not in expected_methods:
                    reason = f"Input file method ({detected_method}) doesn't match {description} method"
                    print(f"Warning: Skipping {description} specification - {reason}")
                    skipped_specs.append((filename, reason))
                    continue
            
            spec_path = resolve_spec_path(filename)
            if not spec_path.exists():
                print(f"Error: {description} specification file {filename} does not exist")
                sys.exit(1)
            spec_files.append(str(spec_path))
            print(f"Using {description} specification file: {filename}")
    
    # Check if at least one specification file is provided
    if not spec_files and not any([args.macromolecules, args.citation, args.authors, args.funding, args.keywords]):
        print("Error: At least one specification file must be provided")
        print("Use --xray, --xray_serial, --em, --nmr for method-specific files, or --macromolecules, --citation, --authors, --funding, --keywords for optional files")
        sys.exit(1)
    
    # Add additional specification files if requested
    additional_specs = [
        (args.macromolecules, "specs/MACROMOLECULES.csv", "macromolecules"),
        (args.citation, "specs/CITATION.csv", "citation"),
        (args.authors, "specs/AUTHORS.csv", "authors"),
        (args.funding, "specs/FUNDING.csv", "funding"),
        (args.keywords, "specs/KEYWORDS.csv", "keywords")
    ]
    
    for flag, filename, description in additional_specs:
        if flag:
            spec_path = resolve_spec_path(filename)
            if not spec_path.exists():
                print(f"Error: {description.title()} specification file {filename} does not exist")
                sys.exit(1)
            spec_files.append(str(spec_path))
            print(f"Also using {description} specification file: {filename}")
    
    # Run metadata import
    success = import_metadata(args.input_file, spec_files, output_file, args.merge_to_file, merge_output_file,
                             log_file, skipped_specs)
    
    if not success:
        sys.exit(1)


if __name__ == "__main__":
    main() 