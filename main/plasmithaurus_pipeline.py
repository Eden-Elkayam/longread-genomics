import os
import subprocess
import util

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

def main():
    raw_reads = util.get_reads()
    output_dir = util.get_output_dir(raw_reads)
    filtered_reads, genome_size = filter_reads(raw_reads, output_dir)
    assembly = util.flye_assembly(filtered_reads, output_dir, genome_size)
    polished = util.medaka_polishing(filtered_reads, assembly, output_dir)
    reference_file = util.reference_gbk_to_fasta()
    comparison = util.make_dnadiff(polished, output_dir, "plasmithaurus_pipeline", reference_file)
    return comparison

main()
