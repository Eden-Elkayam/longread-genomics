import os
from pathlib import Path
from datetime import datetime
import subprocess
from Bio import SeqIO

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

def flye_assembly(input_reads, output_dir, genome_size=4800000, threads=8, quality="hq"):
    flye_output_dir = Path(str(output_dir) + "/flye_out")

    # Make sure output directory exists
    flye_output_dir.mkdir(parents=True, exist_ok=True)

    # Run Flye
    subprocess.run(
        [
            "flye",
            f"--nano-{quality}", str(input_reads),
            "--out-dir", str(flye_output_dir),
            "--threads", str(threads),
            "--iterations",  str(2),
            "--genome-size", str(genome_size),
            "--plasmids"
        ],
        check=True
    )
    flye_file = str(flye_output_dir) + "/assembly.fasta"
    print("Flye assembly complete. path to file:", flye_file)
    return flye_file

def medaka_polishing(reads, assembly, out_dir, threads=8, m="r1041_e82_400bps_sup_v4.3.0"):
    medaka_output_dir = Path(str(out_dir) + "/medaka_out")
    # Make sure output directory exists
    medaka_output_dir.mkdir(parents=True, exist_ok=True)

    subprocess.run(
        [
            "medaka_consensus",
            "-i", str(reads),
            "-d", str(assembly),
            "-o", str(medaka_output_dir),
            "-t", str(threads),
            "-m", m
        ],
        check=True
    )
    output_file = str(medaka_output_dir) + "/consensus.fasta"
    print(f"Medaka polishing complete. Output:{output_file}")
    return output_file


def make_dnadiff(sequence_path,output_dir,prefix,reference_path=None):
    while reference_path is None:
        reference_path = input("Enter the absolute path to your reference file:")
        if not os.path.isfile(reference_path):
            reference_path = None
            print("File was not found, please try again.")
        if not reference_path.endswith(".fasta"):
            reference_path = None
            print("File is not a fasta file, please try again.")

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

def reference_gbk_to_fasta(input_path=None):
    if input_path is None:
        input_path = input("Enter the absolute path to your reference file:")
    i = input_path.rfind("/")
    dir = input_path[:i]
    name = input_path[i + 1:(input_path.rfind("."))]
    output_file = f"{dir}/{name}.fasta"
    count = SeqIO.convert(input_path, "genbank", output_file, "fasta")
    print("Converted %i records" % count)
    return output_file


