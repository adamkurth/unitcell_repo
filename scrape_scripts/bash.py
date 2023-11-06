from rcsbsearchapi.search import AttributeQuery
import random

# Create a query for space group P1211
query = AttributeQuery("rcsb_struct_symmetry.space_group_name_H-M", "exact_match", "P 1 2 1 1")

# Execute the query to get PDB IDs
pdb_ids = query("entry_id")

# Print the list of PDB IDs before randomization
print("PDB IDs before randomization:")
pdb_id_list = []

# Suppose 'pdb_ids' is an enumerate object
for index, pdb_id in pdb_ids:
    # Add the PDB ID to the list
    pdb_id_list.append(pdb_id)

# Now 'pdb_id_list' contains all the PDB IDs
print("PDB IDs after randomization:")
print(pdb_id_list)


# # Shuffle the list of PDB IDs to randomize the order
# random.shuffle(pdb_ids)

# # Select the first N IDs (e.g., N = 10 for 10 random IDs)
# n = 10  # Adjust this number as needed
# random_pdb_ids = pdb_ids[:n]

# # Print the list of PDB IDs after randomization
# print("PDB IDs after randomization:")
# print(random_pdb_ids)