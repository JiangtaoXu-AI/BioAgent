from langchain.tools import BaseTool
import csv
import paramiko
import requests
import os
from typing import Dict, List, Set, Optional
import json
class GenomeCollectorTool(BaseTool):
    name = "ProteinGenomeCollector"
    description = """Input several seed sequences of a protein along with an optional species name to build a genome database. Example input: {"ncbi_ids": ["OIR22820.1", "RZD35475.1".........], "species": "Archaea"})"""

    def __init__(self):
        super().__init__()

    def _run(self, input_dict: Dict) -> str:
        """
        Executes BLAST searches for the provided NCBI IDs, extracts Hit IDs, and sends them to a Linux server.
        :param input_dict: A dictionary containing a species name and a list of NCBI IDs.
                          Example: {'species': 'Methanobacteriota archaeon', 'ncbi_ids': ['RZD39016.1', 'RZD35475.1']}
        :return: A message indicating the number of queried sequences, retrieved sequences, and the Linux server output.
        """
        try:
            print(type(input_dict))
            print(input_dict)
            if isinstance(input_dict, str):
                try:
                    input_dict = json.loads(input_dict)
                except json.JSONDecodeError:
                    return "Input must be a valid JSON string or dictionary."

            # Extract species and NCBI IDs from the input dictionary
            species_name = input_dict.get("species")  # 获取 species 值
            if not species_name:  # 如果 species 未提供或为空
                species_name = "Archaea"  # 设置为默认值 # Default to 'Archaeon' if species is not provided
            ncbi_ids = input_dict.get("ncbi_ids", [])

            if not ncbi_ids:
                return "No valid NCBI IDs provided."

            # Perform BLAST searches for each NCBI ID and extract Hit IDs
            hit_ids = set()
            for accession in ncbi_ids:
                print(f"\nSearching for similar sequences of {accession} based on BLAST")
                test_data = {
                    "query": accession,
                    "eq_menu": species_name,
                    "max_num_seq": 500
                }
                print(test_data)
                result = self._perform_blast_search(test_data)
                print(len(hit_ids))
                if result:
                    # Extract Hit IDs from the result
                    for row in result:
                        hit_id = row.get("Hit_id", "")
                        if hit_id:
                            # Extract the content between | |
                            parts = hit_id.split("|")
                            if len(parts) >= 2:
                                hit_ids.add(parts[1])

            # Send Hit IDs to the Linux server
            if hit_ids:
                print(hit_ids)
                if "RLI38395.1" in hit_ids:
                    print("11")
                self._send_to_linux_server(hit_ids)
                genome_number=self._receive_genome_ids_file()
                # Receive the genome IDs file from the Linux server
                return f"Queried {len (ncbi_ids)} seed sequences, Retrieved {len (hit_ids)} protein sequences and obtain {genome_number} genome sequence."
            else:
                return "No Hit IDs found for the provided NCBI IDs."

        except Exception as e:
            return f"Error occurred {str(e)}"

    async def _arun(self, input_dict: Dict) -> str:
        """Use the tool asynchronously."""
        raise NotImplementedError("Asynchronous execution is not supported for this tool.")

    def _perform_blast_search(self, test_data: Dict) -> List[Dict]:
        """Execute BLAST search"""
        base_url = "http://183.131.35.27:57100/api/v1"
        sync_search_url = f"{base_url}/blast/sync-search"
        try:
            response = requests.post(sync_search_url, json=test_data)
            if response.status_code == 200:
                result = response.json()
                if "csv_content" in result:
                    return self._parse_csv_content(result["csv_content"])
            return []
        except Exception as e:
            print(f"BLAST search failed: {str(e)}")
            return []

    def _parse_csv_content(self, csv_content: str) -> List[Dict]:
        """Parse CSV content"""
        lines = csv_content.split("\n")
        headers = lines[0].split(",")
        results = []
        for line in lines[1:]:
            if line:
                values = line.split(",")
                results.append(dict(zip(headers, values)))
        return results

    def _send_to_linux_server(self, hit_ids: Set[str]) -> str:
        """
        Send Hit IDs to the Linux server for further processing.
        :param hit_ids: A set of Hit IDs.
        :return: The output from the Linux server.
        """
        hostname = "localhost"  # Replace with your Linux server IP
        username = "xjtfi"  # Replace with your username
        password = "123456"  # Replace with your password
        script_path = "/home/xjtfi/blast/ProteinGenomeCollector.py"  # Replace with your script path

        # Convert Hit IDs to a comma-separated string
        hit_ids_str = ",".join(hit_ids)

        try:
            # Connect to the Linux server and run the script
            ssh = paramiko.SSHClient()
            ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
            ssh.connect(hostname, username=username, password=password)

            # Ensure Conda is initialized and activate the bioinfo environment
            command = f"source ~/miniconda3/etc/profile.d/conda.sh && conda activate bioinfo && python {script_path} {hit_ids_str}"

            stdin, stdout, stderr = ssh.exec_command(command)
            output = stdout.read().decode().strip()
            error = stderr.read().decode().strip()
            print(output)

            ssh.close()

            if error:
                raise Exception(f"Error occurred on Linux server: {error}\nOutput: {output}")
            return f"Successfully executed program in Linux system hostname"
        except Exception as e:
            raise Exception(f"Failed to communicate with Linux server: {str(e)}")

    def _receive_genome_ids_file(self):
        """
        Receive the Plsc_genome_ID.txt file from the Linux server and save it to the local result folder.
        """
        hostname = "localhost"  # Replace with your Linux server IP
        username = "xjtfi"  # Replace with your username
        password = "123456"  # Replace with your password
        remote_file_path = "/home/xjtfi/Plsc_genome_ID.txt"  # Path to the file on the Linux server

        local_file_path = os.path.join(r"D:\Desktop\new_start\knowledge_discovery\code\agent\my_agent2\results", "Plsc_genome_ID.txt")  # Local file path


        try:
            # Connect to the Linux server using SFTP
            transport = paramiko.Transport((hostname, 22))
            transport.connect(username=username, password=password)
            sftp = paramiko.SFTPClient.from_transport(transport)

            # Download the file
            sftp.get(remote_file_path, local_file_path)
            sftp.close()
            transport.close()
            with open(local_file_path, 'r', encoding='utf-8') as file:
                line_count = sum(1 for line in file)
            return line_count
        except Exception as e:
            print(f"Failed to receive file from Linux server: {str(e)}")
class GenomeQueryTool(BaseTool):
    name = "GenomeDatabaseQuery"

    description ="""Manually input or open a file to obtain some seed sequences of a protein. These sequences will be used to search the established genome database for genomes containing the protein.Example input in JSON format: {"ncbi_ids": ["RZD39016.1", "RZD35475.1", ...], "file_input": "Yes/No"}."""
    def __init__(self):
        super().__init__()


    def _run(self, input_dict: Dict) -> str:
        """
        Sends a list of NCBI protein IDs to a Linux server for genome database querying.
        :param input_dict: A dictionary containing a list of NCBI protein IDs.
                          Example: {'ncbi_ids': ['RZD39016.1', 'RZD35475.1']}
        :return: The output from the Linux server after querying the genome database.
        """
        try:
            # If input_dict is a string, try to parse it as JSON
            if isinstance(input_dict, str):
                try:
                    input_dict = json.loads(input_dict)
                except json.JSONDecodeError:
                    return "Input must be a valid JSON string or dictionary."

            # Ensure input_dict is a dictionary
            if not isinstance(input_dict, dict):
                return "Input must be a dictionary or a JSON string."

            # Check if input is from a file or manual input
            file_input = input_dict.get("file_input", "No").strip().lower()
            if file_input == "yes":
                # Read NCBI IDs from file
                file_path = r"D:\Desktop\new_start\knowledge_discovery\code\agent\my_agent2\results\GGGPS_query_ID.txt"

                with open(file_path, "r") as f:
                    ncbi_ids = [line.strip() for line in f if line.strip()]
            else:
                # Get NCBI IDs from manual input
                ncbi_ids = input_dict.get("ncbi_ids", [])

            if not ncbi_ids:
                return "No valid NCBI IDs provided."

            # Send NCBI IDs to the Linux server for genome database querying
            self._send_to_linux_server(ncbi_ids)
            row_count=self._receive_genome_ids_file()
            output_file= os.path.join(r"D:\Desktop\new_start\knowledge_discovery\code\agent\my_agent2\results", "genome_classification.txt")
            return f"The relevant genome has been saved to {output_file} and there are {row_count} genomes"

        except Exception as e:
            return f"Error occurred during genome database query: {str(e)}"

    async def _arun(self, input_dict: Dict) -> str:
        """Use the tool asynchronously."""
        raise NotImplementedError("Asynchronous execution is not supported for this tool.")

    def _send_to_linux_server(self, ncbi_ids: List[str]) -> tuple:
        """
        Sends a list of NCBI protein IDs to the Linux server for genome database querying.
        :param ncbi_ids: A list of NCBI protein IDs.
        :return: A tuple containing the output file path and the number of rows in the xlsx file.
        """
        hostname = "localhost"  # Replace with your Linux server IP
        username = "xjtfi"      # Replace with your username
        password = "123456"    # Replace with your password
        script_path = "/home/xjtfi/blast/GenomeDatabaseQuery.py"  # Replace with your script path

        # Convert NCBI IDs to a comma-separated string
        ncbi_ids_str = ",".join(ncbi_ids)

        try:
            # Connect to the Linux server and run the script
            ssh = paramiko.SSHClient()
            ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
            ssh.connect(hostname, username=username, password=password)

            # Ensure Conda is initialized and activate the bioinfo environment
            command = f"source ~/miniconda3/etc/profile.d/conda.sh && conda activate bioinfo && python {script_path} {ncbi_ids_str}"

            stdin, stdout, stderr = ssh.exec_command(command)
            output = stdout.read().decode().strip()
            error = stderr.read().decode().strip()

            ssh.close()

            if error:
                raise Exception(f"Error occurred on Linux server: {error}\nOutput: {output}")
            return f"Successfully retrieved genomes that meet the requirements"
        except Exception as e:
            raise Exception(f"Failed to communicate with Linux server: {str(e)}")
    def _receive_genome_ids_file(self):
        """
        Receive the Plsc_genome_ID.txt file from the Linux server and save it to the local result folder.
        """
        hostname = "localhost"  # Replace with your Linux server IP
        username = "xjtfi"  # Replace with your username
        password = "123456"  # Replace with your password
        remote_file_path = "/home/xjtfi/genome_classification.txt"  # Path to the file on the Linux server

        local_file_path = os.path.join(r"D:\Desktop\new_start\knowledge_discovery\code\agent\my_agent2\results", "genome_classification.txt")  # Local file path


        try:
            # Connect to the Linux server using SFTP
            transport = paramiko.Transport((hostname, 22))
            transport.connect(username=username, password=password)
            sftp = paramiko.SFTPClient.from_transport(transport)

            # Download the file
            sftp.get(remote_file_path, local_file_path)
            sftp.close()
            transport.close()
            with open(local_file_path, 'r', encoding='utf-8') as file:
                line_count = sum(1 for line in file)
            return line_count-1
        except Exception as e:
            print(f"Failed to receive file from Linux server: {str(e)}")
# Example usage
# if __name__ == "__main__":
#     tool = GenomeCollectorTool()
#     input_dict = {
#         "species": "Archaea",  # 物种名称（可选，默认为 'Archaea'）
#         "ncbi_ids": [
#             "OIR22820.1",
#             "OIR22374.1",
#         ]
#     }
#     input_dict = {"ncbi_ids": [ "OIR22820.1","OIR13649.1","OUX24843.1"], "species": ""}
#     #input_dict = '{"ncbi_ids": ["OIR22820.1", "RZD35475.1"], "species": "Archaea"}'
#     result = tool._run(input_dict)
#     print(result)
#
#     tool = GenomeQueryTool()
#     input_dict = {
#         "ncbi_ids": [], "file_input": "Yes"
#     }
#     result = tool._run(input_dict)
#     print(result)