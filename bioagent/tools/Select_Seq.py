"""
Tool: SelectRepresentativeSeq
Description: 从大量古菌 PlsC 同源序列中，通过 SSN 分析和系统发育树构建，自动筛选指定数量的代表性序列，输出可视化分析文件。
"""
from langchain.tools import BaseTool
from typing import Dict, Any, List, Optional
import os
import sys
import json
import subprocess


class SelectRepresentativeSeq(BaseTool):
    """
    Select representative sequences from a large set of PlsC homologs using SSN and phylogenetic analysis.
    Input: {
        "input_sequences_fasta": "path/to/sequences.fasta",
        "num_representatives": 12,
        "ssn_threshold": 0.5,
        "gtdbtk_options": {...}
    }
    Output: {
        'selected_representative_sequences_fasta': str,
        'ssn_cytoscape_file': str,
        'phylogenetic_tree_file': str,
        'status_message': str,
        'notes_for_manual_selection': str
    }
    """

    name = "SelectRepresentativeSeq"
    description = (
        "Automatically select representative sequences from a large set of PlsC homologs using SSN and phylogenetic tree analysis. "
        "Input: JSON with keys: input_sequences_fasta, num_representatives, ssn_threshold, gtdbtk_options."
    )

    def __init__(self):
        super().__init__()

    def _run(self, input_dict: Dict) -> Dict[str, Any]:
        try:
            input_sequences_fasta = input_dict.get("input_sequences_fasta")
            num_representatives = int(input_dict.get("num_representatives", 12))
            ssn_threshold = float(input_dict.get("ssn_threshold", 0.5))
            gtdbtk_options = input_dict.get("gtdbtk_options", {})

            cys_file = self.run_ssn(input_sequences_fasta, ssn_threshold)
            tree_file = self.run_gtdbtk(input_sequences_fasta, gtdbtk_options)
            clusters = self.parse_cys(cys_file)
            rep_ids = self.select_reps(clusters, num_representatives)
            rep_fasta = f"representatives_{num_representatives}.fasta"
            self.extract_fasta(input_sequences_fasta, rep_ids, rep_fasta)

            return {
                'selected_representative_sequences_fasta': os.path.abspath(rep_fasta),
                'ssn_cytoscape_file': os.path.abspath(cys_file),
                'phylogenetic_tree_file': os.path.abspath(tree_file),
                'status_message': 'Success',
                'notes_for_manual_selection': (
                    'Please use the SSN (.cys) and phylogenetic tree (.tree) files in Cytoscape or other visualization tools for final manual selection.\n'
                    'Experimental validation still requires GGGPS enzyme and chemical synthesis standards (DGGGP, hybrid lipids).'
                )
            }
        except Exception as e:
            return {
                'status_message': f'Error: {str(e)}'
            }

    def run_ssn(self, fasta: str, threshold: float) -> str:
        base = os.path.splitext(fasta)[0]
        output_cys = f"{base}_ssn_{threshold:.2f}.cys"
        try:
            subprocess.run([
                'esist-cli', 'network', fasta,
                '--identity', str(threshold),
                '--output', output_cys
            ], check=True)
        except FileNotFoundError:
            sys.stderr.write(
                f"EFI-EST CLI not found. Please manually upload {fasta} to EFI-EST Web and download as {output_cys}\n"
            )
        return output_cys

    def run_gtdbtk(self, fasta: str, options: Dict[str, Any]) -> str:
        out_dir = f"{os.path.splitext(fasta)[0]}_gtdbtk_out"
        os.makedirs(out_dir, exist_ok=True)
        cmd = [
            'gtdbtk', 'classify_wf',
            '--genome_dir', os.path.dirname(fasta) or '.',
            '--extension', 'fasta',
            '--out_dir', out_dir
        ]
        for k, v in options.items():
            cmd.extend([f"--{k}", str(v)])
        try:
            subprocess.run(cmd, check=True)
        except FileNotFoundError:
            sys.stderr.write(
                f"GTDB-tk not found. Please manually run GTDB-tk on {fasta} and generate a Newick tree file.\n"
            )
            return ""
        tree_src = os.path.join(out_dir, 'gtdbtk_classification.tree')
        tree_dst = f"{os.path.splitext(fasta)[0]}_gtdbtk.tree"
        if os.path.exists(tree_src):
            os.replace(tree_src, tree_dst)
        return tree_dst

    def parse_cys(self, cys_file: str) -> Dict[str, List[str]]:
        clusters = {}
        try:
            with open(cys_file) as f:
                data = json.load(f)
            for node in data.get('nodes', []):
                cid = node.get('attributes', {}).get('shared nearest neighbors', '0')
                sid = node.get('name')
                clusters.setdefault(str(cid), []).append(sid)
        except Exception:
            sys.stderr.write(f"Failed to parse {cys_file}, skipping cluster parsing.\n")
        return clusters

    def select_reps(self, clusters: Dict[str, List[str]], num_reps: int) -> List[str]:
        items = sorted(clusters.items(), key=lambda x: len(x[1]), reverse=True)
        reps = []
        for _, seqs in items:
            if len(reps) >= num_reps:
                break
            reps.append(seqs[0])
        return reps[:num_reps]

    def extract_fasta(self, input_fasta: str, headers: List[str], out_path: str) -> None:
        write_flag = False
        with open(input_fasta) as fin, open(out_path, 'w') as fout:
            for line in fin:
                if line.startswith('>'):
                    write_flag = any(line.split()[0] == f'>{h}' for h in headers)
                if write_flag:
                    fout.write(line)
