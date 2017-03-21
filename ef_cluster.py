# Common functions for clustering

# Needs to be done in chunks to prevent massive command lines
def run_mash_sketch(file_names, kmer_size, sketch_size):

    mash_sep = " "
    for chunk in range(((len(file_names) - 1) // mash_chunk_size) + 1):
        chunk_start = chunk * mash_chunk_size
        chunk_end = (chunk+1) * mash_chunk_size
        if chunk_end > len(file_names):
            chunk_end = len(file_names)

        mash_command = str(args.mash_exec) + " sketch -k " + kmer_size + " -s " + sketch_size + " -o reference" + str(chunk) + " " + mash_sep.join(file_names[chunk_start:chunk_end])
        retcode = subprocess.call(mash_command, shell=True)
        if retcode < 0:
            sys.stderr.write("Mash sketch failed with signal " + str(retcode) + "\n")
            sys.exit(1)

    if (len(file_names) // mash_chunk_size) > 0:
        paste_join = ".msh reference"
        mash_paste_command = str(args.mash_exec) + " paste reference reference" + paste_join.join([str(chunk) for chunk in range(((len(file_names) - 1) // mash_chunk_size) + 1)]) + ".msh"
        retcode = subprocess.call(mash_paste_command, shell=True)
        if retcode < 0:
            sys.stderr.write("Mash paste failed with signal " + str(retcode) + "\n")
            sys.exit(1)
        else:
            for chunk in range(((len(file_names) - 1) // mash_chunk_size) + 1):
                os.remove("reference" + str(chunk) + ".msh")
    else:
        os.rename("reference0.msh", "reference.msh")

# Read samples from vcf header
def read_header(vcf_file):
    separator = "\t"
    p = subprocess.Popen(['bcftools view -h ' + vcf_file], stdout=subprocess.PIPE, shell=True)
    for line in iter(p.stdout.readline, ''):
        line = line.decode('utf-8')
        line = line.rstrip()
        header_line = re.match("^#CHROM", line)
        if header_line != None:
            samples = line.split(separator)
            del samples[0:9]
            break

    return samples

