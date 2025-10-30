# from paper: https://doi.org/10.1016/j.btre.2025.e00931

import subprocess
from pathlib import Path
import util

def chop(input_file:str, output_dir:str):
    chop_dir = Path(output_dir + "/porechop_output")
    chop_dir.mkdir(parents=True, exist_ok=True)
    output_file_path = f"{chop_dir}/trimmed.fastq"
    subprocess.run(["porechop",
     "-i", input_file,
     "-o", output_file_path], check=True)
    return  output_file_path


def canu(reads,output_dir ):
    canu_dir = Path(output_dir + "/canu_output")
    canu_dir.mkdir(parents=True, exist_ok=True)
    canu_path = "/Users/edenelkayam/canu-2.3/bin/canu" #manually installed so..

    subprocess.run([canu_path, "-correct",
                    "-p", "corrected",
                    "-d", canu_dir,
                    "genomeSize=4.8m",
                    "-nanopore", reads])
    print("Canu is done! Output Folder is:", canu_dir)
    return canu_dir

def chopper_filter(input_file, output_dir):
    chopper_dir = Path(output_dir + "/chopper_output")
    chopper_dir.mkdir(parents=True, exist_ok=True)
    output_file_path = f"{chopper_dir}/filtered_reads.fastq"
    with open(output_file_path, "w") as out:
        subprocess.run(
            ["chopper",
             "-i", input_file,
             "-t", "8",
             "--trim-approach", "trim-by-quality",
             "--cutoff", "10"],
            stdout=out,
        )
        print("Chopper is done! Output file:", output_file_path)
    return output_file_path

def main():
    raw = "/Users/edenelkayam/flamholz_lab/assembly/SRR1302084/SRR31302084.fastq"
    output_dir =util.get_output_dir(raw)
    chopped = "/Users/edenelkayam/flamholz_lab/assembly/assembly_pipeline/assembly_results/SRR31302084.fastq_20251030-1045/porechop_output/trimmed.fastq"
    #corrected = canu(chopped, output_dir) # takes so long
    filtered = chopper_filter(raw, output_dir)
    return filtered

main()


