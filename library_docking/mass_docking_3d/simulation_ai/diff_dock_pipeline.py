import requests
import csv

url = "https://health.api.nvidia.com/v1/biology/mit/diffdock"
header_auth = "Bearer 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"

def _upload_asset(input):
    assets_url = "https://api.nvcf.nvidia.com/v2/nvcf/assets"

    headers = {
        "Authorization": header_auth,
        "Content-Type": "application/json",
        "accept": "application/json",
    }

    s3_headers = {
        "x-amz-meta-nvcf-asset-description": "diffdock-file",
        "content-type": "text/plain",
    }

    payload = {
        "contentType": "text/plain", 
        "description": "diffdock-file"
    }

    response = requests.post(
        assets_url, headers=headers, json=payload, timeout=30
    )

    response.raise_for_status()

    asset_url = response.json()["uploadUrl"]
    asset_id = response.json()["assetId"]

    response = requests.put(
        asset_url,
        data=input,
        headers=s3_headers,
        timeout=300,
    )

    response.raise_for_status()
    return asset_id

def read_file_as_string(file_path):
    with open(file_path, 'r') as file:
        return file.read()

def process_ligand_requests(protein_file, ligand_list_file, output_csv):
    # Upload the protein file once
    protein_id = _upload_asset(read_file_as_string(protein_file))

    # Read ligand filenames from the text file
    with open(ligand_list_file, 'r') as f:
        ligand_files = ["./../../library/3d_sdf/" + line.strip() for line in f if line.strip()]

    # Open the CSV file for writing
    with open(output_csv, mode='w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        
        # Write the initial header with Ligand File
        csv_writer.writerow(["Ligand File", "Confidence Array (Dynamic Columns)"])

        for ligand_file in ligand_files:
            try:
                # Upload the ligand file
                ligand_id = _upload_asset(read_file_as_string(ligand_file))

                # Prepare headers
                headers = {
                    "Content-Type": "application/json",
                    "NVCF-INPUT-ASSET-REFERENCES": ",".join([protein_id, ligand_id]),
                    "Authorization": header_auth
                }

                # Send the API request
                response = requests.post(
                    url,
                    headers=headers,
                    json={
                        "ligand": ligand_id,
                        "ligand_file_type": "sdf",
                        "protein": protein_id,
                        "num_poses": 5,
                        "time_divisions": 20,
                        "steps": 18,
                        "save_trajectory": True,
                        "is_staged": True,
                    },
                )

                # Check the response and write to CSV
                if response.status_code == 200:
                    response_data = response.json()
                    if "position_confidence" in response_data:
                        confidence = response_data["position_confidence"]
                        # Make sure the confidence is an array
                        if isinstance(confidence, list):
                            # Write ligand file name followed by confidence values
                            csv_writer.writerow([ligand_file] + confidence)
                        else:
                            csv_writer.writerow([ligand_file, "Invalid confidence format"])
                    else:
                        csv_writer.writerow([ligand_file, "Position Confidence not found"])
                else:
                    csv_writer.writerow([ligand_file, f"Error: {response.status_code}"])
            except Exception as e:
                csv_writer.writerow([ligand_file, f"Error: {e}"])

# Update these file paths
protein_file_path = "./../../../raw_files/hCAR1.pdb"
ligand_list_file_path = "./ligand_files.txt"
output_csv_path = "./diff_dock_scores_5.csv"

# Process the requests and save results to CSV
process_ligand_requests(protein_file_path, ligand_list_file_path, output_csv_path)
