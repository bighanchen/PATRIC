import gzip
import sys
from collections import defaultdict

def main():
    if len(sys.argv) != 4:
        print("Usage: python script.py <file1.vcf.gz> <file2.txt> <output.txt>")
        sys.exit(1)
    
    file1_path = sys.argv[1]
    file2_path = sys.argv[2]
    output_path = sys.argv[3]
    
    # Step 1: Read file2 and build position index dictionary
    file2_dict = defaultdict(list)
    
    try:
        with open(file2_path, 'r') as f2:
            for line in f2:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    chrom = parts[0]
                    pos = parts[1]
                    # Use tuple (chrom, pos) as key, store entire line
                    file2_dict[(chrom, pos)].append(line.strip())
    
        # Step 2: Read file1 (gzipped vcf file) and write output
        with gzip.open(file1_path, 'rt') as f1, open(output_path, 'w') as out_file:
            for line in f1:
                # Skip comment lines
                if line.startswith('#'):
                    continue
                    
                parts = line.strip().split('\t')
                if len(parts) < 2:
                    continue
                    
                chrom1 = parts[0]
                pos1 = parts[1]  # Second column of file1
                key = (chrom1, pos1)
                
                # Find matches - multiple file2 lines may match the same position
                if key in file2_dict:
                    for match_line in file2_dict[key]:
                        # Output combination of file1 line and file2 matching line
                        out_file.write(f"{line.strip()}\t{match_line}\n")
    
    except FileNotFoundError as e:
        print(f"Error: File not found - {e.filename}")
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()