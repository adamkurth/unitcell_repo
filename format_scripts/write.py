import os

def write(final, intensity, phases, spacegroup, output_dir):
    # Create a directory to store the output files
    destination = f"{spacegroup}_output"
    output_path = os.path.join(output_dir, destination)
    os.makedirs(output_path, exist_ok=True)

    # Write the intensity data to a file
    intensity_filename = f"{output_path}/{spacegroup}_intensity_data.txt"
    with open(intensity_filename, "w") as f:
        f.write(intensity.to_string())

    # Write the phases data to a file
    phases_filename = f"{output_path}/{spacegroup}_phases_data.txt"
    with open(phases_filename, "w") as f:
        f.write(phases.to_string())
    
    # Write the final crystal data to a file
    crystal_filename = f"{output_path}/{spacegroup}_crystal_data.txt"
    with open(crystal_filename, "w") as f:
        f.write(final.to_string())

    # Print the output directory
    print(f"Output directory: {os.path.abspath(output_path)}")
