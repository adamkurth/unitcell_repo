import os
import glob

# file must be within the current working directory of the files you want to rename

def rename_files(file_extension):
    # Get the list of files matching the pattern in the current working directory
    files = glob.glob(f"*{file_extension}")

    for file_path in files:
        # Extract the filename and extension
        filename = os.path.basename(file_path)
        file_extension = os.path.splitext(filename)[1]

        # Extract the protein ID (first 4 characters of the filename)
        pdb_id = filename[:4]

        # Determine the type (input or output) and construct the new filename
        if "input" in filename:
            new_filename = f"{pdb_id}_input_formatted{file_extension}"
        elif "output" in filename:
            new_filename = f"{pdb_id}_output_formatted{file_extension}"
        else:
            print(f"Skipping {file_path} - type (input/output) not found in filename")
            continue

        # Construct the new file path
        new_file_path = os.path.join(os.getcwd(), new_filename)

        # Rename the file
        os.rename(file_path, new_file_path)
        print(f"Renamed '{file_path}' to '{new_file_path}'")

if __name__ == "__main__":
    file_extension = input("Enter the file extension (.csv or .txt): ")
    rename_files(file_extension)
