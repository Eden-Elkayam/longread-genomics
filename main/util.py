import os
from pathlib import Path
from datetime import datetime
import subprocess

def get_reads() -> str:
    input_path = input("Enter the absolute path to your file:")
    if not os.path.isfile(input_path):
        print("File was not found, please try again.")
        get_reads()
    # continue when found a legitimate path to file
    if not input_path.endswith((".fasta", ".fasta.gz", ".fastq", ".fastq.gz")):
        print("File is not a fasta or fastq file, please try again.")
        get_reads()
    return input_path

# returns the path to the result directory
def get_output_dir(input_path:str) -> str:
    output_path = input("Enter the absolute path to your output folder:")
    if output_path:
        output_path = Path(output_path)
    else:
        cur = os.getcwd()
        output_path = Path(cur + "/..")
    output_path.mkdir(parents=True, exist_ok=True)

    # Now that output_path is all good, we need to create a new directory for the results

    last_slash = input_path.rfind("/")
    file_name = input_path[last_slash + 1:]
    timestamp = datetime.now().strftime("%Y%m%d-%H%M")
    output_dir = str(output_path) + "/assembly_results/" + file_name + "_" + str(timestamp)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    return str(output_dir)


def make_dnadiff(reference_path, sequence_path,output_dir,prefix):
    dnadiff_output_dir = Path(str(output_dir) + "/dnadiffout")
    # create directory
    dnadiff_output_dir.mkdir(parents=True, exist_ok=True)
    # Run dnadiff
    subprocess.run(
        [
            "dnadiff",
            f"-p", str(dnadiff_output_dir)+"/"+prefix,
            reference_path,
            sequence_path
        ],
        check=True
    )
    dnadiff_report = str(dnadiff_output_dir) + f"/{prefix}.report"
    print("dnadiff report generated. path to file:", dnadiff_report)
    return dnadiff_report

def dna_diff_check():
    reference_path = "/Users/edenelkayam/flamholz_lab/assembly/SRR1302084/sequence.fasta"
    sequence_path = "/Users/edenelkayam/flamholz_lab/assembly/assembly_pipeline/assembly_results/SRR31302084.fastq_20251024-0931/medaka_out/consensus.fasta"
    output_dir = "/Users/edenelkayam/flamholz_lab/assembly/assembly_pipeline/assembly_results/SRR31302084.fastq_20251024-0931/medaka_out/"
    prefix = "ecoli"
    make_dnadiff(reference_path, sequence_path, output_dir, prefix)

#dna_diff_check()