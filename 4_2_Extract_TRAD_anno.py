import sys
import re

def extract_info_field(info, field_name):
    """Extract the value of the specified field from the INFO field"""
    pattern = f"{field_name}=([^;]+)"
    match = re.search(pattern, info)
    return match.group(1) if match else ""

def extract_sample_fields(sample_data, format_keys):
    """Extract GT and REPCN fields from sample data"""
    if sample_data == "":
        return "", ""
    
    sample_values = sample_data.split(':')
    format_dict = dict(zip(format_keys, sample_values))
    
    gt = format_dict.get('GT', '')
    repcn = format_dict.get('REPCN', '')
    
    return gt, repcn

def main():
    if len(sys.argv) != 4:
        print("Usage: python script.py <family_id> <input_file> <output_file>")
        sys.exit(1)
    
    family_id = sys.argv[1]
    input_file = sys.argv[2]
    output_file = sys.argv[3]
    
    try:
        with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
            # Write header
            header = [
                "FamilyID", "CHROM", "POS", "END", "RU", "REF",
                "Sample1_GT", "Sample1_REPCN",
                "Sample2_GT", "Sample2_REPCN", 
                "Sample3_GT", "Sample3_REPCN",
                "Col13", "Col14", "Col15", "Col16", "Col17"
            ]
            #f_out.write('\t'.join(header) + '\n')
            
            for line in f_in:
                if line.startswith('#'):  # Skip comment lines
                    continue
                    
                fields = line.strip().split('\t')
                
                # Determine number of samples - key fix section
                # Total columns minus first 9 columns (VCF standard + FORMAT) and last 10 columns (annotation columns)
                total_cols = len(fields)
                n_samples = total_cols - 9 - 10  # 9 VCF columns + 10 annotation columns
                
                # Extract basic fields
                chrom = fields[0]
                pos = fields[1]
                info = fields[7]
                
                # Extract END, RU, REF from INFO field
                end = extract_info_field(info, "END")
                ru = extract_info_field(info, "RU")
                ref = extract_info_field(info, "REF")
                
                # Extract format fields
                format_keys = fields[8].split(':')
                
                # Extract sample data
                sample_data = []
                # Samples start from column 9 (index 9)
                for i in range(min(3, n_samples)):  # Process up to 3 samples
                    sample_idx = 9 + i
                    if sample_idx < len(fields):
                        gt, repcn = extract_sample_fields(fields[sample_idx], format_keys)
                        sample_data.extend([gt, repcn])
                    else:
                        sample_data.extend(["", ""])
                
                # If fewer than 3 samples, fill remaining positions with empty strings
                while len(sample_data) < 6:  # 3 samples * 2 fields = 6
                    sample_data.extend(["", ""])
                
                # Extract the 7th, 5th, 3rd, 2nd, and 1st columns from the end
                col13 = fields[total_cols - 7] if total_cols >= 7 else ""
                col14 = fields[total_cols - 5] if total_cols >= 5 else ""
                col15 = fields[total_cols - 3] if total_cols >= 3 else ""
                col16 = fields[total_cols - 2] if total_cols >= 2 else ""
                col17 = fields[total_cols - 1] if total_cols >= 1 else ""
                
                # Build output row
                output_fields = [
                    family_id, chrom, pos, end, ru, ref,
                    sample_data[0], sample_data[1],  # Sample1 GT and REPCN
                    sample_data[2], sample_data[3],  # Sample2 GT and REPCN
                    sample_data[4], sample_data[5],  # Sample3 GT and REPCN
                    col13, col14, col15, col16, col17
                ]
                
                f_out.write('\t'.join(output_fields) + '\n')
                
        print(f"File processed successfully, output saved to: {output_file}")
        
    except FileNotFoundError:
        print(f"Error: Input file not found {input_file}")
        sys.exit(1)
    except Exception as e:
        print(f"Error occurred while processing file: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()