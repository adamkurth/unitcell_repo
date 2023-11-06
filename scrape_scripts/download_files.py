import requests
import os

class ProteinDownloader:
    def __init__(self, base_dir="/Users/adamkurth/Documents/vscode/CXFEL_Image_Analysis/CXFEL/unitcell_project/run_sfall/pdb"):
        self.base_dir = base_dir
        os.makedirs(base_dir, exist_ok=True)

    def is_url_accessible(self, url):
        """Check if a URL is accessible."""
        response = requests.head(url)
        return response.status_code == 200

    def download_files(self, protein_ids, space_group):
        """Download PDB files for specified proteins."""
        PDB_URL_TEMPLATE = "https://files.rcsb.org/download/{}.pdb"
        
        # Convert space_group to lowercase for consistent processing
        space_group = space_group.lower()

        for protein_id in protein_ids:
            print(f"Processing Protein: {protein_id}, Space Group: {space_group}")

            pdb_url = PDB_URL_TEMPLATE.format(protein_id)

            # Check if PDB file is accessible
            if self.is_url_accessible(pdb_url):
                # Create directory for the space group if it doesn't exist
                space_group_dir = os.path.join(self.base_dir, f"pdb_{space_group}")
                os.makedirs(space_group_dir, exist_ok=True)

                # Download PDB file
                response = requests.get(pdb_url)
                with open(os.path.join(space_group_dir, f"{protein_id}.pdb"), 'wb') as f:
                    f.write(response.content)
                print(f"Downloaded {protein_id}.pdb")
            else:
                print(f".pdb file not available for {protein_id}")
