from concurrent.futures import ThreadPoolExecutor, as_completed
import os
import pickle
import subprocess
import re
import csv
import pandas as pd
import pytaxonkit as ptk
from taxonomy import TaxonUtils

class BlastProcessor:
    def __init__(self, db):
        self.db = db
        
    def run_megablast(self, fasta_file, out_file):
        if not os.path.exists(fasta_file):
            print(f"No such file: {fasta_file}")
            return None
        
        cmd = (
            f'blastn -db {self.db} -num_threads 10 -query {fasta_file} '
            f'-task megablast -out {out_file} '
            f'-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sgi saccver slen qlen staxids sscinames stitle qcovs qcovhsp"'
        )
        
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, encoding='utf-8')
        while True:
            output = process.stdout.readline()
            if output == '' and process.poll() is not None:
                break
            if output:
                print(output.strip())
        return out_file
    
    def is_fasta(self, file_path):
        try:
            with open(file_path, 'r') as file:
                first_line = file.readline().strip()
                return first_line.startswith('>')
        except Exception as e:
            print(f"Error reading {file_path}: {e}")
            return False

    def process_fasta(self, fasta_file, output_dir):
        print(fasta_file)
        if not self.is_fasta(fasta_file):
            print(f"File {fasta_file} is not a valid FASTA file. Skipping.")
            return None
        
        out_file = os.path.join(output_dir, os.path.basename(fasta_file) + ".blast6.tsv")
        self.run_megablast(fasta_file, out_file)
        return out_file
    
    def flatten_fasta_files(self, fasta_files):
        """Flatten a list of tuples into a flat list of files."""
        flat_list = []
        for item in fasta_files:
            if isinstance(item, tuple):
                flat_list.extend(item)  # Add each element in the tuple to the list
            else:
                flat_list.append(item)  # Add single files directly
        return flat_list


    def process_fasta_multithread(self, fasta_files, output_dir, max_workers=4):
        results = []
        
        flat_fasta_files = self.flatten_fasta_files(fasta_files)
        results = []
        
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_fasta = {
                executor.submit(self.process_fasta, fasta, output_dir): fasta
                for fasta in flat_fasta_files
            }
            for future in as_completed(future_to_fasta):
                fasta_file = future_to_fasta[future]
                try:
                    result = future.result()
                    results.append(result)
                    print(f"Completed processing: {fasta_file}")
                except Exception as exc:
                    print(f"FASTA file {fasta_file} generated an exception: {exc}")
        
        return results


class BlastParser:
    EXCLUDE_PATTERN =  re.compile(r'MH590414\.[0-9]|MH590388\.[0-9]|MG934417\.[0-9]|OR951374\.[0-9]')
    EXCLUDE_PATTERN2 =  re.compile(r'Uncultured human fecal virus')

    def __init__(self, bout_dir, acc2taxid_loader, taxon_structure_pickle='taxon_structure.pkl',best_hit_only = False,qcov=75,pident=90, intersect = True):
        self.bout_dir = bout_dir
        self.acc2taxid_loader = acc2taxid_loader
        self.taxon_structure_pickle = taxon_structure_pickle
        self.utils = self.initialize_taxon_utils()
        self.result_dict = {}
        self.best_hit_only = best_hit_only
        self.qcov = qcov
        self.pident = pident
        self.intersect = intersect


    def initialize_taxon_utils(self):
        taxon_structure = self.load_taxon_structure()
        taxon_utils = TaxonUtils()
        taxon_utils.build_taxon_structure(taxon_structure)
        return taxon_utils

    def load_taxon_structure(self):
        if os.path.exists(self.taxon_structure_pickle):

            print(f"Loading taxon structure from {self.taxon_structure_pickle}")
            with open(self.taxon_structure_pickle, 'rb') as f:
                taxon_structure = pickle.load(f)
        else:
            print("Generating taxon structure...")
            taxon_structure = ptk.list([10239], threads=20, raw=True)
            with open(self.taxon_structure_pickle, 'wb') as f:
                pickle.dump(taxon_structure, f)
            print(f"Taxon structure saved to {self.taxon_structure_pickle}")
        
        return taxon_structure

    def load_acc2taxid(self):
        return self.acc2taxid_loader.load()

    def should_skip_row(self, row, qcov, pident):
        return float(row[19]) < qcov or float(row[2]) < pident or re.search(self.EXCLUDE_PATTERN, row[1]) is not None or re.search(self.EXCLUDE_PATTERN2, row[18]) is not None

    def get_match_info(self, row, acc2taxid):
        match_info = {
            'sseqid': row[1],
            'staxids': row[16].split(';')[0],
            'sscinames': row[17].split(';')[0],
            'stitle': row[18].split(';')[0],
            'qcov': row[19],
            'pident': float(row[2])
        }
        if int(match_info["staxids"]) == 0:
            acc = re.sub(r"\.[0-9]+", "", match_info["sseqid"])
            match_info['staxids'] = int(acc2taxid.get(acc, 0))
        return match_info

    def update_results(self, query_id, match_info, results):
        if self.best_hit_only:
            if query_id not in results:
                results[query_id] = [match_info]
            else:
                current_best_pident = results[query_id][0]['pident']
                if match_info['pident'] > current_best_pident:
                    results[query_id] = [match_info]  # Replace with new best hit
                elif match_info['pident'] == current_best_pident:
                    results[query_id].append(match_info)  # Add to the list of best hits
        else:
            results.setdefault(query_id, []).append(match_info)

    def parse_blast_output_add_tax(self, blast_output_file, results=None):
        if results is None:
            results = {}
        acc2taxid = self.load_acc2taxid()

        with open(blast_output_file, 'r') as file:
            reader = csv.reader(file, delimiter='\t')
            for row in reader:
                if self.should_skip_row(row, self.qcov, self.pident):
                    continue
                query_id = row[0]
                match_info = self.get_match_info(row, acc2taxid)
                if int(match_info["staxids"]) == 0:
                    continue
                self.update_results(query_id, match_info, results)
        return results

    def categorize_matches(self, results):
        categories = {}
        paired_queries = {query_id.rsplit('/', 1)[0] if '/' in query_id else query_id for query_id in results}
        for read_uid in paired_queries:
            # for_read, rev_read = f"{read_uid}/1", f"{read_uid}/2"
            for_read = f"{read_uid}/1" if f"{read_uid}/1" in results else read_uid
            rev_read = f"{read_uid}/2" if f"{read_uid}/2" in results else read_uid

            for_hits = results.get(for_read, [])
            rev_hits = results.get(rev_read, [])
            for_hits_subjid = {hit['sseqid'] for hit in for_hits}
            rev_hits_subjid = {hit['sseqid'] for hit in rev_hits}
            inter_sseqids = for_hits_subjid.intersection(rev_hits_subjid)
            if self.intersect:
            # If intersection=True, only add entries with intersecting sseqid values
                if inter_sseqids:
                    common_hits = [hit for hit in for_hits if hit['sseqid'] in inter_sseqids]
                    categories[read_uid] = {'hits': common_hits, 'mtype': 2}
            else:
                # Original behavior: handle both intersection and non-intersection cases
                if inter_sseqids:
                    common_hits = [hit for hit in for_hits if hit['sseqid'] in inter_sseqids]
                    categories[read_uid] = {'hits': common_hits, 'mtype': 2}
                else:
                    if len(for_hits) > 0 and len(rev_hits) > 0:
                        categories[for_read] = {'hits': for_hits, 'mtype': 1}
                        categories[rev_read] = {'hits': rev_hits, 'mtype': 1}
                    else:
                        categories[read_uid] = {'hits': for_hits + rev_hits, 'mtype': 1}
        return categories
    
    def append_lineage(self, tutils, hits_list, level='species'):
        try:
            list_taxa = []
            for hit in hits_list:      
                found_taxon = tutils.find_by_taxid(int(hit['staxids']))
                if found_taxon:
                    list_taxa.append(found_taxon)
                # else:
                #     print(hit)
            majority_taxon = tutils.find_majority_taxon(list_taxa, level)
            return majority_taxon

        except Exception as err:
            print(f"Unexpected {err=}, {type(err)=} {hit}")

    def update_sample_dict(self, sample_dict, taxon):
        if taxon not in sample_dict:
            sample_dict[taxon] = {
                'count': 0,
                'taxid': taxon.taxid,
                'lineage': taxon.get_all_lineage()  # Assuming this method returns the lineage as a string
            }   
        sample_dict[taxon]['count'] += 1

    def process_file_pair(self, idx, for_bout, rev_bout):
        print(f'Processing: {for_bout} and {rev_bout}')
        blast_res_forward = self.parse_blast_output_add_tax(for_bout)
        blast_res_pair = self.parse_blast_output_add_tax(rev_bout, results=blast_res_forward)        
        
        paired_hits = self.categorize_matches(blast_res_pair)
        
        pair_result = {}
        for ruid, hits in paired_hits.items():
            major_taxon = self.append_lineage(self.utils, hits['hits'], 'species')
            self.update_sample_dict(pair_result, major_taxon)
        
        return idx, pair_result

    def process_all_files(self, findex="_kneaddata_paired_R1", rindex="_kneaddata_paired_R2", max_workers=4):
        tasks = []
        print (len([name for name in os.listdir(self.bout_dir) if os.path.isfile(name)]))

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            for parent, _, files in os.walk(self.bout_dir):
                for f in files:
                    if not ".blast6.tsv" in f:
                        continue
                    # if not "vir17_" in f:
                    #     continue
                    if rindex in f:
                        continue  # Skip reverse index files
                    idx = os.path.basename(f).replace(findex, "")
                    for_bout = os.path.join(parent, f)
                    rev_f = os.path.basename(f).replace(findex, rindex)
                    rev_bout = os.path.join(parent, rev_f)
                    tasks.append(executor.submit(self.process_file_pair, idx, for_bout, rev_bout))

            for future in as_completed(tasks):
                result = future.result()
                if result:
                    idx, pair_result = result
                    self.result_dict[idx] = pair_result
        print("All files processed.")