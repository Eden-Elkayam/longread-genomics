import os
import subprocess
from pathlib import Path
from datetime import datetime

# user provides name of input file.

# return the full path to the input file, making sure it exists and in the right format
def get_reads():
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
def get_output_dir(input_path):
    output_path = input("Enter the absolute path to your output folder:")
    if output_path:
        output_path = Path(output_path)
        output_path.mkdir(parents=True, exist_ok=True)
    else:
        output_path = os.getcwd()
    # Now that output_path is all good, we need to create a new directory for the results

    last_slash = input_path.rfind("/")
    file_name = input_path[last_slash + 1:]
    timestamp = datetime.now().strftime("%Y%m%d-%H%M")
    output_dir = str(output_path) + "/assembly_results/" + file_name + "_" + str(timestamp)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    return str(output_dir)

# Filter reads based on x100 coverage, reads that undergo this are ready to be assembled
def filter_reads(input_file, output_dir):
    os.chdir(output_dir)
    # New files to be created:
    best95 = output_dir + "/reads_keep95.fastq"
    best250mb = output_dir + "/reads_250Mb.fastq"
    overlaps = output_dir + "/overlaps.paf.gz"
    mini_assembly = output_dir + "/assembly_sketch.gfa"
    x100coverage = output_dir + "/reads_x100coverage.fastq"

    # Remove worst 5%
    with open(best95, "w") as out:
        subprocess.run(["filtlong", "--keep_percent", "95", input_file], stdout=out)

    # Create map tp assess genome size:

    # keep best 250 Mb for the map
    with open(best250mb, "w") as out:
        subprocess.run(["filtlong", "--target_bases", "250000000", best95], stdout=out)

    # Overlap graph (Node = read, edge = overlap)
    with open(overlaps, "wb") as out:
        p1 = subprocess.Popen(["minimap2", "-x", "ava-ont", "-t", "8", best250mb, best250mb], stdout=subprocess.PIPE)
        p2 = subprocess.Popen(["gzip"], stdin=p1.stdout, stdout=out)
        p1.stdout.close()  # allow p1 to receive a SIGPIPE if p2 exits
        p2.communicate()

# Feed the overlap graph and the 250Mb file to create a sketch assembly for knowing genome size
    with open(mini_assembly, "w") as out:
        subprocess.run(
            ["miniasm", "-f", best250mb, overlaps], stdout=out, check=True)

    # finds genome size by the sketch we just made
    def find_sketch_size(sketch):
        count = 0
        with open(sketch, "r") as file:
            for line in file:
                if line.startswith("S"):
                    l = line.split('LN:i:')
                    count += int(l[1])
        return count

    # according to size, re-filter for x100 coverage
    genome_size = find_sketch_size(mini_assembly)
    target = genome_size * 100
    with open(x100coverage, "w") as out:
        subprocess.run(
            ["filtlong", "--mean_q_weight", "10",
             "--target_bases", str(target),
             str(best95)],
            stdout=out,
            check=True
        )
    print(f"reads were succesfully filtered. path to filtered read is: {x100coverage}, Genome size: {genome_size}")
    return x100coverage, genome_size


def assemble(input_filtered, output_dir, genome_size, threads=8):
    flye_output_dir = Path(str(output_dir) + "/flye_out")

    # Make sure output directory exists
    flye_output_dir.mkdir(parents=True, exist_ok=True)

    # Run Flye
    subprocess.run(
        [
            "flye",
            "--nano-hq", str(input_filtered),
            "--out-dir", str(flye_output_dir),
            "--threads", str(threads),
            "--iterations",  str(2),
            "--genome-size", str(genome_size),
            "--plasmids"
        ],
        check=True
    )
    assembly_flye = str(flye_output_dir) + "/assembly.fasta"
    print("Flye assembly complete. path to file:", assembly_flye)
    return assembly_flye

def polishing(input_filtered, assembly_flye, out_dir, threads=8, m="r1041_e82_400bps_sup_v4.3.0"):
    medaka_output_dir = Path(str(out_dir) + "/medaka_out")
    # Make sure output directory exists
    medaka_output_dir.mkdir(parents=True, exist_ok=True)

    subprocess.run(
        [
            "medaka_consensus",
            "-i", str(input_filtered),
            "-d", str(assembly_flye),
            "-o", str(medaka_output_dir),
            "-t", str(threads),
            "-m", m
        ],
        check=True
    )
    print(f"Medaka polishing complete. Output: {out_dir}/consensus.fasta")

def main():
    raw_reads = get_reads()
    output_dir = get_output_dir(raw_reads)
    filtered_reads, genome_size = filter_reads(raw_reads, output_dir)
    assembly = assemble(filtered_reads, output_dir, genome_size)
    polishing(filtered_reads, assembly, output_dir)

main()
