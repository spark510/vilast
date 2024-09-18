import argparse
import os
import pandas as pd
from blast_handler import BlastProcessor
from blast_handler import BlastParser
from taxonomy import Acc2TaxidLoader
from taxonomy import TaxonomyProcessor
from datetime import datetime

# from taxon import Taxon
# from taxon_utils import TaxonUtils

def find_pairs(inputdir, findex, rindex=None):
    files = os.listdir(inputdir)
    forward_files = [f for f in files if findex in f]
    file_pairs = []
    unmatched_forward = []
    
    if rindex:
        reverse_files = [f for f in files if rindex in f]
        for f_file in forward_files:
            # Derive the reverse file name
            expected_reverse_file = f_file.replace(findex, rindex)
            if expected_reverse_file in reverse_files:
                file_pairs.append((os.path.abspath(os.path.join(inputdir, f_file)), 
                                   os.path.abspath(os.path.join(inputdir, expected_reverse_file))))
            else:
                unmatched_forward.append(os.path.abspath(os.path.join(inputdir, f_file)))
        return file_pairs, unmatched_forward
    else:
        # If rindex is not provided, return forward files only
        for f_file in forward_files:
            file_pairs.append((os.path.abspath(os.path.join(inputdir, f_file)),))
        return file_pairs, []
    
def get_cpu_cores():
    num_cores = os.cpu_count()
    if num_cores is None:
        print("Unable to determine the number of CPU cores.")
        return 1  # Fallback to 1 core if the number cannot be determined
    return num_cores

def main():
    parser = argparse.ArgumentParser(description="Run MEGABLAST, parse the results, and simplify taxonomy in a pipeline.")
    parser.add_argument("-i","--inputdir", type=str, required=True, help="Directory containing input fasta files.")
    parser.add_argument("-o","--outdir", type=str, required=True, help="Directory to save output files.")
    parser.add_argument("-f","--findex", type=str, required=True, help="Pattern to target forward strand files. ex. _1, _F ")
    parser.add_argument("-r","--rindex", type=str, required=False, help="Pattern to target reverse strand files. ex. _2, _R ")
    parser.add_argument("-d","--db", type=str, required=False, help="Path to the BLAST database file.")
    parser.add_argument("-w","--worker", type=int, default=10, help="Number of worker processes to use.")
    parser.add_argument("-m","--midentity", type=float, default=75, required=False, help="cutoff for valid hit identity")
    parser.add_argument("-c","--mcoverage", type=float, default=90, required=False, help="cutoff for valid hit coverage")
    parser.add_argument("-s","--steps", type=str, required=False, default='1-2', help="Steps to execute. Use format '1-2' to run steps 1 to 2, or '1,2' to run steps 1 and 2.")
    
    args = parser.parse_args()

    steps = []
    if '-' in args.steps:
        start, end = map(int, args.steps.split('-'))
        steps = list(range(start, end + 1))
    else:
        steps = list(map(int, args.steps.split(',')))

    input_dir = args.inputdir
    out_dir = args.outdir
    os.makedirs(out_dir,exist_ok=True)
    
    findex = args.findex
    rindex = args.rindex
    n_worker = args.worker

    acc2taxid_file='/home/share/programs/LPBvirus/necessaries/acc2taxid.tsv'
    vacc2taxid_file='/home/share/programs/LPBvirus/necessaries/virus_acc2taxid.tsv.gz'
    vacc2taxid_pickle='/home/share/programs/LPBvirus/necessaries/virus_acc2taxid.pkl'
    taxonkit_pickle = '/home/share/programs/LPBvirus/necessaries/taxonkit.pkl'
    
    # Step 1: Process FASTA files
    if 0 in steps:
        print(0)
        qcoverage = args.midentity
        pidentity = args.mcoverage
  
        acc2taxid_loader = Acc2TaxidLoader(acc2taxid_file,vacc2taxid_file,vacc2taxid_pickle)  # Assuming this class is implemented correctly
        blast_parser = BlastParser(out_dir, acc2taxid_loader, taxonkit_pickle,qcov=qcoverage,pident=pidentity)
        blast_parser.process_all_files(findex=findex, rindex=rindex, max_workers=n_worker)

        sqlite3_file = "/home/share/programs/LPBvirus/necessaries/virus_taxonomy.db"
        db_name = "taxonomy"
        processor = TaxonomyProcessor(sqlite3_file, db_name)
        full_df = processor.process_dataframe(blast_parser.result_dict)
        file_suffix = datetime.today().strftime("%Y%m%d")+"_"+str(qcoverage)+"_"+str(pidentity)
        processor.process_and_save_aggregations(full_df, out_dir, file_suffix)
        
    if 1 in steps:
        if not (args.inputdir and args.outdir and args.findex):
            parser.error("Step 1 requires --inputdir, --outdir, --findex arguments.")
        
        # Default database if not provided
        db = args.db if args.db else '/home/share/databases/ncbi_nt_viruses_custom/ncbi_virus_host_mammalia_prok_wo_covid.rmdup.fasta'
        # Find file pairs
        file_pairs, unmatched_forward = find_pairs(input_dir, findex, rindex)
        print("File pairs found:", file_pairs)
        if unmatched_forward:
            print("Unmatched forward files:", unmatched_forward)
        
        num_cores = get_cpu_cores()
        blast_worker =  min(n_worker // 10, num_cores // 10)
                
        # Initialize BlastProcessor with the database
        bp = BlastProcessor(db)
        bp.process_fasta_multithread(file_pairs, out_dir, max_workers=blast_worker)
        print("STEP1. DONE: megablast")

    # Step 2: Further processing (Placeholder, implementation required)
    if 2 in steps:
        qcoverage = args.midentity
        pidentity = args.mcoverage
        acc2taxid_loader = Acc2TaxidLoader(acc2taxid_file,vacc2taxid_file,vacc2taxid_pickle)  # Assuming this class is implemented correctly
        blast_parser = BlastParser(out_dir, acc2taxid_loader, taxonkit_pickle,qcov=qcoverage,pident=pidentity)
        print(out_dir)
        blast_parser.process_all_files(findex=findex, rindex=rindex, max_workers=n_worker)

        sqlite3_file = "/home/share/programs/LPBvirus/necessaries/virus_taxonomy.db"
        db_name = "taxonomy"
        processor = TaxonomyProcessor(sqlite3_file, db_name)
        full_df = processor.process_dataframe(blast_parser.result_dict)
        file_suffix = datetime.today().strftime("%Y%m%d")+"_"+str(qcoverage)+"_"+str(pidentity)
        processor.process_and_save_aggregations(full_df, out_dir, file_suffix)

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"An error occurred: {e}")