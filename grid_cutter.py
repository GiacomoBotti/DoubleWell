#############################################
# Code wrote by Chatgpt to parse SOFT wfn   #
# plots and make cuts                       #
#############################################


import sys

def filter_rows(filename):
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            
            # Skip blank lines
            if not line:
                continue
            
            parts = line.split()
            
            # Ensure there are at least 2 columns
            if len(parts) < 2:
                continue
            
            # Compare first two columns
            # x=y cut
            if parts[0] == parts[1]:
            # y=0 cut
            #if parts[1] == 0:
                print(line)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py input.dat")
    else:
        filter_rows(sys.argv[1])
