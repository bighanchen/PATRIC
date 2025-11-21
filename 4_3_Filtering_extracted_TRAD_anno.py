import argparse
import pandas as pd
import openpyxl
from openpyxl.utils import get_column_letter
import sys

def parse_alleles(allele_str):
    """Parse allele string, return set of repeat numbers"""
    if not allele_str or pd.isna(allele_str) or allele_str == '':
        return set()
    return set(map(int, allele_str.split(',')))

def has_valid_common_allele(row, ref_repeat_idx=5):
    """Check if a row meets the conditions:
    1. All samples have overlapping repeat numbers
    2. The overlap contains numbers other than the reference repeat number
    """
    # Get reference repeat number
    try:
        ref_repeat = int(row[ref_repeat_idx])
    except (ValueError, IndexError):
        return False
    
    alleles_sets = []
    
    # Parse column 8 (index 7)
    if len(row) > 7 and row[7]:
        alleles_sets.append(parse_alleles(row[7]))
    
    # Parse column 10 (index 9)
    if len(row) > 9 and row[9]:
        alleles_sets.append(parse_alleles(row[9]))
    
    # Parse column 12 (index 11), if not empty
    if len(row) > 11 and row[11]:
        alleles_sets.append(parse_alleles(row[11]))
    
    # Skip row if no allele data
    if not alleles_sets:
        return False
    
    # Calculate intersection of all sets
    common_alleles = set.intersection(*alleles_sets)
    
    # Check conditions:
    # 1. Intersection is not empty
    # 2. Intersection contains numbers other than the reference repeat number
    if len(common_alleles) > 0:
        # If intersection has numbers other than reference repeat
        if len(common_alleles - {ref_repeat}) > 0:
            return True
        # Or if intersection has multiple numbers (even if including reference repeat)
        elif len(common_alleles) > 1:
            return True
    
    return False

def filter_file(input_file, output_txt, output_xlsx):
    """Filter file and output in two formats"""
    
    # Read file
    print(f"Reading file: {input_file}")
    try:
        # Use flexible method to read file, handling possible column count inconsistencies
        data = []
        with open(input_file, 'r') as f:
            for line in f:
                # Remove line breaks and split by tabs
                row = line.strip().split('\t')
                data.append(row)
        
        if not data:
            print("Error: Input file is empty")
            return
        
        print(f"Successfully read {len(data)} rows of data")
        
    except Exception as e:
        print(f"Error reading file: {e}")
        return
    
    # Filter data
    print("Filtering data...")
    filtered_data = [row for row in data if has_valid_common_allele(row)]
    print(f"Remaining {len(filtered_data)} rows after filtering")
    
    if not filtered_data:
        print("Warning: No data passed the filtering criteria")
        return
    
    # Output TXT file
    print(f"Outputting TXT file: {output_txt}")
    try:
        with open(output_txt, 'w') as f:
            for row in filtered_data:
                f.write('\t'.join(row) + '\n')
        print("TXT file output completed")
    except Exception as e:
        print(f"Error outputting TXT file: {e}")
        return
    
    # Output Excel file
    print(f"Outputting Excel file: {output_xlsx}")
    try:
        # Determine maximum column count
        max_cols = max(len(row) for row in filtered_data)
        
        # Create DataFrame
        columns = [f'Col_{i+1}' for i in range(max_cols)]
        df = pd.DataFrame(filtered_data, columns=columns)
        
        # Create Excel writer
        with pd.ExcelWriter(output_xlsx, engine='openpyxl') as writer:
            df.to_excel(writer, index=False, header=False)
            
            # Get workbook and worksheet
            workbook = writer.book
            worksheet = writer.sheets['Sheet1']
            
            # Set columns 7, 9, 11 (indices 6,8,10) as text format to prevent conversion to dates
            text_columns = [6, 8, 10]  # 0-based indices
            
            for col_idx in text_columns:
                if col_idx < max_cols:
                    col_letter = get_column_letter(col_idx + 1)  # +1 because Excel is 1-based
                    for row in range(2, len(filtered_data) + 2):  # +2 because Excel first row is header, data starts from row 2
                        cell = worksheet[f'{col_letter}{row}']
                        # Set cell format to text
                        cell.number_format = '@'
        
        print("Excel file output completed")
        print(f"Processing completed successfully! Input: {len(data)} rows, Output: {len(filtered_data)} rows")
        
    except Exception as e:
        print(f"Error outputting Excel file: {e}")
        return

def main():
    parser = argparse.ArgumentParser(description='Filter TRAD data, retain rows where all samples have overlapping repeat numbers and the overlap contains non-reference repeat numbers')
    parser.add_argument('input_file', help='Input file path')
    parser.add_argument('output_txt', help='Output TXT file path')
    parser.add_argument('output_xlsx', help='Output Excel file path')
    
    args = parser.parse_args()
    
    # Check if input file exists
    try:
        with open(args.input_file, 'r') as f:
            pass
    except FileNotFoundError:
        print(f"Error: Input file '{args.input_file}' does not exist")
        sys.exit(1)
    except Exception as e:
        print(f"Error: Cannot read input file '{args.input_file}': {e}")
        sys.exit(1)
    
    filter_file(args.input_file, args.output_txt, args.output_xlsx)

if __name__ == "__main__":
    main()