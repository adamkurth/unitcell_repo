import sys
import os
import pandas as pd

def is_header(line):
    # Define a list of header keywords
    header_keywords = ["num", "h", "k", "l", "FREE", "FP", "SIGFP", "FCalc", "PHICalc", "FC", "PHIC", "FWT", "PHWT", "DELFWT", "PHDELWT", "FOM", "FC_ALL", "PHIC_ALL_LS"]
    # Check if any keyword is present in the line
    return any(keyword in line for keyword in header_keywords)

def get_skiprows(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()
    return [i for i, line in enumerate(lines) if is_header(line)]

def format_output(output_file_path, drop_questionable_rows):
    if not output_file_path.endswith(".txt"):
        return 0

    skiprows = get_skiprows(output_file_path)
    df = pd.read_csv(output_file_path, delim_whitespace=True, header=None, skiprows=skiprows)

    num_question_marks = df.applymap(lambda x: '?' in str(x)).sum().sum()

    if drop_questionable_rows:
        df = df[~df.applymap(lambda x: '?' in str(x)).any(axis=1)]

    df = df.applymap(lambda x: str(x).replace('nan', '').replace(',', ''))

    header = ["num", "h", "k", "l", "FREE", "FP", "SIGFP", "FCalc", "PHICalc"]
    output_file_name = os.path.join(os.path.dirname(output_file_path), os.path.splitext(os.path.basename(output_file_path))[0] + "_formatted.csv")

    while df.shape[1] < len(header):
        df[df.shape[1]] = ''

    df = df.iloc[:, :len(header)]
    df.to_csv(output_file_name, index=False, header=header)

    return num_question_marks

def format_input(input_file_path, drop_questionable_rows):
    skiprows = get_skiprows(input_file_path)
    df = pd.read_csv(input_file_path, delim_whitespace=True, header=None, skiprows=skiprows)

    header = ["num", "h", "k", "l", "FREE", "FP", "SIGFP", "FC", "PHIC", "FC_ALL", "PHIC_ALL", "FWT", "PHWT", "DELFWT", "PHDELWT", "FOM", "FC_ALL_LS", "PHIC_ALL_LS"]

    if drop_questionable_rows:
        df = df[~df.applymap(lambda x: '?' in str(x))]

    while df.shape[1] < len(header):
        df[df.shape[1]] = ''

    df = df.iloc[:, :len(header)]
    output_file_name = os.path.join(os.path.dirname(input_file_path), os.path.splitext(os.path.basename(input_file_path))[0] + "_formatted.csv")
    df.to_csv(output_file_name, index=False, header=header)

    num_question_marks = df.applymap(lambda x: '?' in str(x)).sum().sum()
    return num_question_marks

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python script.py file1_input.txt file2_input.txt ... file1_output.txt file2_output.txt ...")
        sys.exit(1)

    input_files = []
    output_files = []

    for arg in sys.argv[1:]:
        if arg.endswith("_input.txt") or arg.endswith("_input_formatted.txt"):
            input_files.append(arg)
        elif arg.endswith("_output.txt") or arg.endswith("_out.txt") or arg.endswith("_output_formatted.txt"):
            output_files.append(arg)
        else:
            print(f"Skipping {arg} - not an input or output file")

    total_question_marks = 0

    for input_file in input_files:
        num_question_marks = format_input(input_file, False)
        total_question_marks += num_question_marks
        print(f"{input_file}: type INPUT")
        print(f"{input_file}: {num_question_marks} question marks")

    for output_file in output_files:
        num_question_marks = format_output(output_file, False)
        total_question_marks += num_question_marks
        print(f"{output_file}: type OUTPUT")
        print(f"{output_file}: {num_question_marks} question marks")

    print(f"Total question marks in all files: {total_question_marks}")
    drop_rows = input("Would you like to drop rows with question marks (yes/no)? ").lower()
    if drop_rows == "yes" or drop_rows == "y":
        total_question_marks = 0
        for input_file in input_files:
            num_question_marks = format_input(input_file, True)
            total_question_marks += num_question_marks
            print(f"{input_file}: type INPUT (after dropping rows with question marks)")
            print(f"{input_file}: {num_question_marks} question marks removed")
        for output_file in output_files:
            num_question_marks = format_output(output_file, True)
            total_question_marks += num_question_marks
            print(f"{output_file}: type OUTPUT (after dropping rows with question marks)")
            print(f"{output_file}: {num_question_marks} question marks removed")
        print(f"Total question marks in all files: {total_question_marks}")
