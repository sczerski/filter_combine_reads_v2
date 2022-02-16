#Author: Sam Czerski
#Date: 2/16/2022
'''
Description of Program: This script accepts two required input files. 1) a fasta file of an OTU from Pacbio post-ccs/demultiplex, 2) a fasta file of a genus downloaded from NCBI
via entrez edirect utils. GI number is not output from entrez edirect and therefore v2 will work with NR/accessions rather than GI numbers. The NR numbers will
be repeated in two columns because the pacbio IDs are different, and while for ncbi, one column maps to tree nodes and other acts as tree labels.
Output is a combined file of the Pacbio and NCBI sequences filtered to only contain reads where more than one identical read was found (Pacbio only)
the Pacbio file must contain this "size" information in the header. This threshold can be changed depending on your data.
Lastly, an identifier is added, either "pacbio" or "ncbi" so read's origin can be easily known for analysis.
'''
#Imports
from Bio import SeqIO
import glob

#Takes the OTU fasta file as input and returns a fasta file of filtered records based on a chosen "size" parameter
def filter_single_reads(otu_fasta_file):
    handles = open(otu_fasta_file, "r")
    handles_list = []
    pre_filt_count = 0
    for handle in SeqIO.parse(handles, "fasta"):
        handles_list.append(handle.id)
        pre_filt_count+=1
    print("Number of reads before filtering:", pre_filt_count)


    post_filt_count = 0
    filtered_reads = get_read_size(handles_list)
    for i in filtered_reads:
        post_filt_count+=1
    print("Number of reads after filtering:", post_filt_count)

    handles.close()

    return get_pacbio_seqs(filtered_reads)

#Takes a list of handles parsed from BioPython SeqIO.parse and returns a list of filtered reads based on "size"
def get_read_size(handles_list):
    filtered_reads = []
    for element in handles_list:
        #Delimiter can be changed depending on which is used in your file
        handle_id_items = element.split(sep=';')
        #depending on the number of elements in the handle, this index must be changed
        #"size" must be present on pacbio headers post CCS/Demultiplexing
        if handle_id_items[3] not in ["size=1", "size=2", "size=3", "size=4", "size=5"]:
            filtered_reads.append(element)
        

    return filtered_reads

#Takes a list of filtered reads, finds the associated sequence, and writes to a new output file
def get_pacbio_seqs(filtered_reads):
    handles = open(otu_fasta_file[0], "r")
    total_records = 0

    with open("filtered_records.fasta", "w") as fr:
        for handle in SeqIO.parse(handles, "fasta"):
            if handle.id in filtered_reads:
                SeqIO.write(handle, fr, "fasta")
                total_records+=1

    #print(total_records)
    handles.close()

    return "filtered_records.fasta"

#Accepts two fasta files as input, the filtered pacbio reads and the ncbi downloaded reads, and combines them into one file
def combine_pacbio_w_ncbi(new_records, ncbi_records):
    with open(ncbi_records, "r") as ncbi_file, open("ncbi_seqs_no_blanks.fasta", "w") as ncbi_out_file, open(new_records, "r"):
    #first remove blank lines from ncbi fastas
        for line in ncbi_file:
            if line.strip():
                ncbi_out_file.write(line)

        files_to_combine = ["ncbi_seqs_no_blanks.fasta", new_records]
        #now combine the files in a new file containing all of the filtered records
        with open("all_filtered_seqs.fasta", "w") as all_record_file:
            for file in files_to_combine:
                with open(file) as in_file:
                    for line in in_file:
                        all_record_file.write(line)



        return "all_filtered_seqs.fasta"

#Function to print the number of reads after combining pacbio and ncbi reads, may be helpful for keeping track of records
def print_total_num_of_seqs(all_filtered_records_file):
    with open(all_filtered_records_file, "r") as f:
        total_seqs = 0
        for line in f:
            if line.startswith(">"):
                total_seqs+=1

        print("Total number of records in all_filtered_seqs:", total_seqs)

#This function accepts the combined fasta of pacbio and ncbi reads, and adds an "origin identifier" to each header.
def add_seq_origin_to_id(all_filtered_seqs_file):
    pre_seqs_ncbi = []
    pre_seqs_pacbio = []
    post_all_seqs = []
    with open(all_filtered_seqs_file, "r") as f:
        for header in f:
            if header.startswith(">NR_") == True:
                pre_seqs_ncbi.append(header)

            elif header.startswith(">m") == True:
                pre_seqs_pacbio.append(header)

        id_ncbi = "ncbi"
        for i in pre_seqs_ncbi:
            post_seqs_ncbi = i.split(sep=" ")
            x = 1
            post_seqs_ncbi.insert(x, id_ncbi)
            x += 8
            post_seqs_ncbi = [' '.join(i for i in post_seqs_ncbi)]
            for seq in post_seqs_ncbi:
                post_all_seqs.append(seq)


        id_pacbio = "pacbio"
        for j in pre_seqs_pacbio:
            post_seqs_pacbio = j.split(sep=";")
            x = 2
            post_seqs_pacbio.insert(x, id_pacbio)
            x += 3
            post_seqs_pacbio = [';'.join(j for j in post_seqs_pacbio)]
            for seq in post_seqs_pacbio:
                post_all_seqs.append(seq)

        #post_all_seqs now contains all the sequence headers with the correct id
        #note: for ncbi, it appears that the the NR_# and ncbi_id are delimited by only a space, but it is in fact a tab,
        # I was able to upload into R successfully as a tsv, so not sure why visually it only looks like a space.

        #write new headers to a new file with associated sequences
        all_filtered_seqs_w_id = get_all_seqs(all_filtered_seqs_file, post_all_seqs)

        return all_filtered_seqs_w_id

def make_all_reads_info_tsv(all_filtered_reads_w_id):
    pre_reads_info_ncbi = []
    pre_reads_info_pacbio = []
    all_reads_info_tsv = []
    with open(all_filtered_reads_w_id, "r") as f:
        for header in f:
            if header.startswith(">NR_") == True:
                pre_reads_info_ncbi.append(header[1:])
            elif header.startswith(">m") == True:
                pre_reads_info_pacbio.append(header[1:])

        with open("all_reads_info.tsv",  'w') as ari_tsv:
            #ncbi
            for i in pre_reads_info_ncbi:
                post_reads_info_ncbi = i.split(sep=" ")
                extra_info = post_reads_info_ncbi[1:]

                extra_info = [item.strip() for item in extra_info]

                genus_species = extra_info[1] + "_" + extra_info[2]
                nr_num = post_reads_info_ncbi[0]
                no_quotes = ""
                i=3

                while i < len(extra_info):
                    no_quotes += " " + extra_info[i]
                    i += 1
                no_quotes = no_quotes.strip()
                w_quotes = '"{}"\n'.format(no_quotes)


                #removing redundant elements
                post_reads_info_ncbi = post_reads_info_ncbi[0:2]

                #appending new consolidated elements
                post_reads_info_ncbi.append(nr_num)
                post_reads_info_ncbi.append(genus_species)
                post_reads_info_ncbi.append(w_quotes)


                #reordering elements for ease of adding annotation layers
                swap_element_positions(post_reads_info_ncbi, 2, 3)


                #Add blank elements to list so when formatting tsv in R, everything is lined up correctly
                post_reads_info_ncbi.insert(3, "")
                post_reads_info_ncbi.insert(4, "")

                #joining elements into single entries in list
                post_reads_info_ncbi = ['\t'.join(i for i in post_reads_info_ncbi)]

                for handle in post_reads_info_ncbi:
                    all_reads_info_tsv.append(handle)

        #pacbio
            for j in pre_reads_info_pacbio:
                post_reads_info_pacbio = j.split(sep=";")
                #get elements to be edited
                barcode_label = post_reads_info_pacbio[1]
                ccs_passes = post_reads_info_pacbio[3]
                size = post_reads_info_pacbio[4]

                #I think lstrip should work here but it doesn't, this feels like a long-winded way of doing this...
                #keep everything after the "="
                barcode_label = barcode_label[13:]
                ccs_passes = ccs_passes[4:]
                size = size[5:]
                size_w_space = "{}\n".format(size)
                pacbio_id = "pacbio"


                #remove old elements and add new ones
                post_reads_info_pacbio.pop(4)
                post_reads_info_pacbio.pop(3)
                post_reads_info_pacbio.pop(2)
                post_reads_info_pacbio.pop(1)

                post_reads_info_pacbio.append(pacbio_id)
                post_reads_info_pacbio.append(barcode_label)
                post_reads_info_pacbio.append(ccs_passes)
                post_reads_info_pacbio.append(size_w_space)

                post_reads_info_pacbio[0] = post_reads_info_pacbio[0].strip()
                post_reads_info_pacbio[1] = post_reads_info_pacbio[1].strip()

                post_reads_info_pacbio = ['\t'.join(j for j in post_reads_info_pacbio)]

                #somehow a space is getting added in as the second element in the list which is messing up the format
                for j in post_reads_info_pacbio:
                    post_post_reads_info_pacbio = j.split(sep='\t')
                    post_post_reads_info_pacbio.pop(1)

                post_reads_info_pacbio = ['\t'.join(j for j in post_post_reads_info_pacbio)]

                for handle in post_reads_info_pacbio:
                    all_reads_info_tsv.append(handle)


            for new_header in all_reads_info_tsv:
                ari_tsv.write(new_header)



def get_all_seqs(all_filtered_seqs_file, post_all_seqs_list):
    handles = open(all_filtered_seqs_file, "r")
    i = 0
    total_records = 0
    with open("all_filtered_seqs_w_id.fasta", "w") as fr:
        for line in handles:
            if line.startswith(">") == True:
                line_w_id = line.replace(line, post_all_seqs_list[i])
                fr.write(line_w_id)
                i += 1
                total_records += 1
            else:
                fr.write(line)

    handles.close()
    print("Number of records wrote to all_filtered_seqs_w_id.fasta", total_records)

    return "all_filtered_seqs_w_id.fasta"


#Swap function
def swap_element_positions(list, pos1, pos2):
    list[pos1], list[pos2] = list[pos2], list[pos1]
    return list

#main function
if __name__ == '__main__':
    #Get Input
    otu_fasta_file = glob.glob("*_eefilt_subset.fasta")
    ncbi_records = glob.glob("*_sequence.fasta") 

    #Filter out singleton reads where only one identical sequence was found (size=1)
    #NOTE: THIS FUNCTION REMOVES SINGLETON READS, BUT CAN (AND POSSIBLY SHOULD) BE CHANGED TO YOUR DESIRED SIZE
    #ie. size=1 may not be appropriate for your data
    new_records = filter_single_reads(otu_fasta_file[0])

    #Combine filtered reads with the reads downloaded from ncbi
    all_filtered_records = combine_pacbio_w_ncbi(new_records, ncbi_records[0])
    print("\nAll pacbio and ncbi reads have been written to all_filtered_seqs.fasta")
    #Print total number of records for sanity
    print_total_num_of_seqs(all_filtered_records)
    #Add sequence identifier to headers ("pacbio" or "ncbi")
    print("\nAdding identifier to sequence headers")
    all_filtered_reads_w_id = add_seq_origin_to_id(all_filtered_records)
    #Make tsv of read info for tree annotations
    make_all_reads_info_tsv(all_filtered_reads_w_id)
    print("\nAll reads info wrote to all_reads_info.tsv")
    print("Fin")

