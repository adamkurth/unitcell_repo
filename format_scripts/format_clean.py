import sys
import os

def format_output(output_file_path, drop_questionable_rows):
    if not output_file_path.endswith(".txt"):
        return 0  # Return 0 if the file is not a .txt file

    with open(output_file_path, 'r') as f:
        lines = f.readlines()

    num_question_marks = 0
    formatted_lines = []

    for line in lines[1:]:
        # Check if the line contains a question mark
        if "?" in line:
            if drop_questionable_rows:
                continue  # Skip the line if it contains a question mark and drop_questionable_rows is True
            else:
                num_question_marks += line.count("?")  # Count the number of "?" in the line
        else:
            num_question_marks += line.count("?")  # Count the number of "?" in the line

        # Split the line into columns
        columns = line.split()

        # Format the columns
        formatted_columns = []
        for i, column in enumerate(columns):
            formatted_column = '{:>12}'.format(column)
            formatted_columns.append(formatted_column)

        # Join the formatted columns into a single line
        formatted_line = ' '.join(formatted_columns)

        # Add the formatted line to the list of formatted lines
        formatted_lines.append(formatted_line)

    # Write the formatted text to a new file
    output_file_name = os.path.basename(output_file_path)
    output_file_name = os.path.splitext(output_file_name)[0] + "_formatted.txt"
    with open(output_file_name, 'w') as f:
        f.write('\n'.join(formatted_lines))
    return num_question_marks

def format_input(input_file_path, drop_questionable_rows):
    if not input_file_path.endswith("_input.txt"):
        return 0  # Return 0 if the file is not an input file

    with open(input_file_path, 'r') as f:
        lines = f.readlines()

    num_question_marks = 0
    formatted_lines = []
    header_row = "         h          k          l         FP      SIGFP         FC       PHIC     FC_ALL   PHIC_ALL        FWT       PHWT     DELFWT    PHDELWT        FOM  FC_ALL_LS PHIC_ALL_LS"

    for line in lines[1:]:
        # Check if the line contains a question mark
        if "?" in line:
            if drop_questionable_rows:
                continue  # Skip the line if it contains a question mark and drop_questionable_rows is True
            else:
                num_question_marks += line.count("?")  # Count the number of "?" in the line
        else:
            num_question_marks += line.count("?")  # Count the number of "?" in the line

        # Split the line into columns
        columns = line.split()

        # Extract the "h" value
        h_value = columns[0]

        # Format the columns (excluding the first column and "FREE")
        formatted_columns = [h_value]  # Add the "h" value
        for i, column in enumerate(columns[1:]):
            if i < 7:
                formatted_column = '{:>10}'.format(column)
            else:
                formatted_column = '{:>12}'.format(column)
            formatted_columns.append(formatted_column)

        # Join the formatted columns into a single line
        formatted_line = ' '.join(formatted_columns)

        # Add the formatted line to the list of formatted lines
        formatted_lines.append(formatted_line)

    formatted_lines.insert(0, header_row)

    # Write the formatted text to a new file
    output_file_name = os.path.basename(input_file_path)
    output_file_name = os.path.splitext(output_file_name)[0] + "_formatted.txt"
    with open(output_file_name, 'w') as f:
        f.write('\n'.join(formatted_lines))
    return num_question_marks

def insert_header(file_type, file_path):
    if file_type == "input":
        header_row = "         h          k          l       FREE         FP      SIGFP         FC       PHIC     FC_ALL   PHIC_ALL        FWT       PHWT     DELFWT    PHDELWT        FOM  FC_ALL_LS PHIC_ALL_LS"
    elif file_type == "output":
        header_row = "         h          k          l       FREE         FP      SIGFP         FC       PHIC"
    else:
        return "Invalid file type"

    # Read the formatted lines from the input file
    with open(file_path, 'r') as f:
        formatted_lines = f.readlines()

    formatted_lines.insert(0, header_row)

    # Write the formatted text back to the input file
    with open(file_path, 'w') as f:
        f.write(''.join(formatted_lines))

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
        num_question_marks = format_input(input_file, False)  # Set to True if you want to drop rows with question marks
        total_question_marks += num_question_marks
        print(f"{input_file}: type INPUT")
        print(f"{input_file}: {num_question_marks} question marks")
        insert_header("input", input_file)

    for output_file in output_files:
        num_question_marks = format_output(output_file, False)  # Set to True if you want to drop rows with question marks
        total_question_marks += num_question_marks
        print(f"{output_file}: type OUTPUT")
        print(f"{output_file}: {num_question_marks} question marks")
        insert_header("output", output_file)

    print(f"Total question marks in all files: {total_question_marks}")
    drop_rows = input("Would you like to drop rows with question marks (yes/no)? ").lower()
    if drop_rows == "yes" or drop_rows == "y":
        total_question_marks = 0
        for input_file in input_files:
            num_question_marks = format_input(input_file, True)  # Drop rows with question marks
            total_question_marks += num_question_marks
            print(f"{input_file}: type INPUT (after dropping rows with question marks)")
            print(f"{input_file}: {num_question_marks} question marks removed")
        for output_file in output_files:
            num_question_marks = format_output(output_file, True)  # Drop rows with question marks
            total_question_marks += num_question_marks
            print(f"{output_file}: type OUTPUT (after dropping rows with question marks)")
            print(f"{output_file}: {num_question_marks} question marks removed")
        print(f"Total question marks in all files: {total_question_marks}")
