import re
import math
import argparse

def parse_info_field(info):
    """Parse VCF INFO field into a dictionary"""
    info_dict = {}
    for item in info.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            info_dict[key] = value
    return info_dict

def process_files(file1_path, file2_path, output_path):
    # Read file1 and build dictionary
    file1_dict = {}
    with open(file1_path, 'r') as f1:
        for line in f1:
            fields = line.strip().split(',')
            if len(fields) < 5:
                continue
                
            # Extract key fields
            tr_id = fields[0]
            motif = fields[1]
            chrom = fields[2]
            start = fields[3]
            end = fields[4]
            
            # Build composite key
            key = f"{chrom}-{start}-{end}-{motif}"
            
            # Parse d-prefixed fields into dictionary
            data_dict = {}
            for item in fields[5:]:
                if item.startswith('d') and ':' in item:
                    k, v = item.split(':', 1)
                    try:
                        data_dict[k] = float(v)
                    except ValueError:
                        pass  # Ignore values that cannot be converted to float
            
            file1_dict[key] = (tr_id, data_dict)
    
    # Process file2 and generate output
    with open(file2_path, 'r') as f2, \
         open(output_path, 'w') as out:        
           
        for line in f2:
            raw_line = line  # Save original line
            line = line.strip()
            
            # Process comment lines
            if line.startswith('#'):
                out.write(raw_line)  # Write original line (with newline)
                continue
                
            # Split VCF line
            parts = line.split('\t')
            if len(parts) < 10:
                continue
                
            chrom = parts[0]
            pos = parts[1]
            
            # Parse INFO field
            info_dict = parse_info_field(parts[7])
            
            # Get END and RU values
            end_val = info_dict.get('END', '')
            ru_val = info_dict.get('RU', '').upper()  # Convert to uppercase
            
            # Build composite key
            variant_key = f"{chrom}-{pos}-{end_val}-{ru_val}"
            
            # Parse FORMAT and SAMPLE fields
            format_fields = parts[8].split(':')
            sample_fields = parts[9].split(':')
            
            # Get REPCN index and value
            if 'REPCN' not in format_fields:
                continue
            repcn_idx = format_fields.index('REPCN')
            if repcn_idx >= len(sample_fields):
                continue
                
            repcn_vals = sample_fields[repcn_idx].split(',')
            if len(repcn_vals) < 2:
                continue
                
            try:
                a = int(repcn_vals[0])
                b = int(repcn_vals[1])
            except ValueError:
                continue
                
            # Look for matching entry in file1 dictionary
            if variant_key in file1_dict:
                tr_id, data_dict = file1_dict[variant_key]
                
                # Build d-value key lists
                a_keys = [f"d{i}_{a}" for i in range(1, 7)]
                b_keys = [f"d{i}_{b}" for i in range(1, 7)]
                
                # Get corresponding values
                a_values = [data_dict.get(k, float('nan')) for k in a_keys]
                b_values = [data_dict.get(k, float('nan')) for k in b_keys]
                
                # Check if d4_a and d4_b meet the condition
                d4_a = data_dict.get(f"d4_{a}", float('nan'))
                d4_b = data_dict.get(f"d4_{b}", float('nan'))
                
                # Condition: d4_a and d4_b are not empty and at least one is less than 0.005
                #if not (isnan(d4_a) or isnan(d4_b)) and (d4_a < 0.005 or d4_b < 0.005):
                if math.isnan(d4_a) or math.isnan(d4_b) or d4_a < 0.005 or d4_b < 0.005:
                    # Write to output file (original VCF line)
                    out.write(raw_line)
                    

def isnan(value):
    """Check if value is NaN"""
    try:
        return math.isnan(value)
    except:
        return False

# Execute processing
if __name__ == "__main__":
    # Set up command line argument parsing
    parser = argparse.ArgumentParser(description='Process TR files and VCF files to generate filtered outputs.')
    parser.add_argument('--file1', required=True, help='Path to the first input file (TR data in CSV format)')
    parser.add_argument('--file2', required=True, help='Path to the second input file (VCF format)')
    parser.add_argument('--output', required=True, help='Path for the output file (filtered VCF)')
    
    # Parse command line arguments
    args = parser.parse_args()
    
    # Execute processing
    process_files(args.file1, args.file2, args.output)