import csv
import os
import sys
import glob

# call the file using python on the command line, then enter the input and output directories
# Does not use the input {filename}.txt as argument, but rather the directory only containing the .txt files
# Usage looks something like this. 
# python txt_to_csv.py ../../spacegroup/P121/temp_test/ ../../spacegroup/P121/input_data/

def txt_to_csv(txt_file_path, csv_file_path):
    print(f"Processing: {txt_file_path} -> {csv_file_path}")
    with open(txt_file_path, 'r') as txt_file:
        lines = txt_file.readlines()
        with open(csv_file_path, 'w', newline='') as csv_file:
            writer = csv.writer(csv_file)
            for line in lines:
                # Detect delimiter (space or tab)
                delimiter = '\t' if '\t' in line else ' '
                writer.writerow(line.split(delimiter))

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python txt_to_csv.py <input_txt_directory> <output_csv_directory>")
        sys.exit(1)

    input_txt_directory = sys.argv[1]
    output_csv_directory = sys.argv[2]

    # Ensure the output directory exists
    os.makedirs(output_csv_directory, exist_ok=True)

    # Get all txt files in the input directory
    txt_files = glob.glob(os.path.join(input_txt_directory, '*.txt'))

    if not txt_files:
        print(f"No .txt files found in {input_txt_directory}")
        sys.exit(1)

    for txt_file_path in txt_files:
        # Construct the corresponding csv file path
        csv_file_name = os.path.basename(txt_file_path).replace('.txt', '.csv')
        csv_file_path = os.path.join(output_csv_directory, csv_file_name)

        # Convert the txt file to csv
        txt_to_csv(txt_file_path, csv_file_path)

    print("Conversion complete.")
