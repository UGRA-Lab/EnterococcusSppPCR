#!/usr/bin/python
from os import path
import argparse
import subprocess
from Bio import SeqIO
from Bio import SearchIO

#This must be edited accordingly to your server structure
blastn_path = "/usr/bin/"

def create_directory(outdir):
    """Function to create the output directory if it is not already created."""
    if not path.isdir(outdir):
        mkoutdir = ["mkdir", outdir]
        subprocess.call(mkoutdir)

def check_primer_hit(hit, match_length):
    for hsp in hit:
        if abs(hsp.hit_start - hsp.hit_end) == match_length:
            if "rev" in hit.id:
                return hsp.query_id, 'rev'
            elif "fwd" in hit.id:
                return hsp.query_id, 'fwd'
        else:
            return False, False

def check_primer_hit_Ehirae(hit, match_length):
    for hsp in hit:
        if abs(hsp.hit_start - hsp.hit_end) == match_length:
            if "mur2" in hit.id:
                 if "rev" in hit.id:
                     return hsp.query_id, 'mur2', 'rev'
                 elif "fwd" in hit.id:
                     return hsp.query_id, 'mur2', 'fwd'
            if "copY" in hit.id:
                 if "rev" in hit.id:
                     return hsp.query_id, 'copY', 'rev'
                 elif "fwd" in hit.id:
                     return hsp.query_id, 'copY', 'fwd'
            if "murG" in hit.id:
                 if "rev" in hit.id:
                     return hsp.query_id, 'murG', 'rev'
                 elif "fwd" in hit.id:
                     return hsp.query_id, 'murG', 'fwd'
        else:
            return False, False, False


def define_enterococcus_specie(fasta_file, ddlE_db, samplename, outdir):
    print("##### Determining specie")
    blastn_cmd = ["{0}/blastn".format(blastn_path), "-query", fasta_file, "-db", ddlE_db, "-outfmt", "6", "-out", "{0}/{1}_blastn.tab".format(outdir, samplename), "-word_size", "16", "-dust", "no", "-soft_masking", "false"]
    process = subprocess.Popen(blastn_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    blastfile = "{0}/{1}_blastn.tab".format(outdir, samplename)
    with open(blastfile) as blast_results:
        species = ['Efs', 'Efm', 'Egallinarum', 'Ecasseliflavus', 'Edurans', 'Eavium']
        species_primers = {sp: {'fwd':False, 'rev':False} for sp in species}
        species_primers['Ehirae'] = {'mur2':{'fwd':False, 'rev':False}, 'copY':{'fwd':False, 'rev':False}, 'murG':{'fwd':False, 'rev':False}}
        contig_names = set()
        for query in SearchIO.parse(blast_results, "blast-tab"):
            for hit in query:
                match_length = int(hit.id.split("_")[-1])
                for sp in species:
                    if sp in hit.id and sp != 'Ehirae':
                        contig_name, direction = check_primer_hit(hit, match_length)
                        if contig_name:
                            contig_names.add(contig_name)
                            species_primers[sp][direction] = True
                    elif sp in hit.id and sp == 'Ehirae':
                        contig_name, gene, direction = check_primer_hit_Ehirae(hit, match_length)                            
                        if contig_name:
                            contig_names.add(contig_name)
                            species_primers[sp][gene][direction] = True
    species_name = "Enterococcus "
    sp_names = {'Efm': 'faecium', 'Efs': 'faecalis', 'Edurans': 'durans', 'Egallinarum': 'gallinarum',
                'Ecasseliflavus':'casseliflavus', 'Eavium':'avium', 'Ehirae':'hirae'}
    species_found = {sp:False for sp in species}
    for sp in species_primers:
        if sp == "Ehirae":
            if all([species_primers[sp][gene]['rev'] and species_primers[sp][gene]['fwd'] for gene in species_primers[sp].keys()]):
                species_name += sp_names[sp]
                species_found[sp] = True
        elif species_primers[sp]['rev'] and species_primers[sp]['fwd']:
            species_name += sp_names[sp]
            species_found[sp] = True
    number_of_species_found = sum([int(v) for v in species_found.values()])
    if number_of_species_found == 0 or number_of_species_found > 1:
       species_name = 'not determined'
    labels = ['SpeciesDetected', 'contigs']
    values = [species_name, "::".join(contig_names)]
    for sp in species_primers.keys():
        if sp == 'Ehirae': continue
        labels.append('{0}Fwd'.format(sp))
        labels.append('{0}Rev'.format(sp))
        if species_primers[sp]['fwd']: values.append('pos')
        else: values.append('neg')
        if species_primers[sp]['rev']: values.append('pos')
        else: values.append('neg')
    for gene in species_primers['Ehirae'].keys():
        labels.append('EhiraeFwd{0}'.format(gene))
        labels.append('EhiraeRevFwd{0}'.format(gene))
        if species_primers['Ehirae'][gene]['fwd']: values.append('pos')
        else: values.append('neg')
        if species_primers['Ehirae'][gene]['rev']: values.append('pos')
        else: values.append('neg')
    return species_name, labels, values

def main():
    parser = argparse.ArgumentParser(description="Script to identify Enterococcus species by insilico PCR for Enterococcus sp. assembled genomes.")
    parser.add_argument('--samplename', help='Name of the sample to use in the outfiles as preffix', type=str, default=None, required=True)
    parser.add_argument("--fasta_file", help="Path to the fasta file with the genome assembly.", type=str, default=None, required=True)
    parser.add_argument('--outdir', help='Name of the output directory', type=str, default=None, required=True)
    parser.add_argument("--ddlE_db", type=str, help="Path to the ddlE_primers.fasta file with the primers sequences for specie identification", required=True)
    args = parser.parse_args()
    create_directory(args.outdir)
    print("##### Performing analysis on sample: {0}".format(args.samplename))
    fasta_file = args.fasta_file
    specie, species_labels, species_values = define_enterococcus_specie(fasta_file, args.ddlE_db, args.samplename, args.outdir)
    with open("{0}/{1}_stats.tab".format(args.outdir, args.samplename), "w") as handler:
        header = ["Sample"] + species_labels
        values = [args.samplename] + species_values
        handler.write("\t".join(map(str, header)))
        handler.write("\n")
        handler.write("\t".join(map(str,values)))      
     
if __name__ == "__main__": main()
