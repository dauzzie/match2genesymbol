import pandas as pd
import subprocess as sub
import gffutils
from datetime import datetime
import argparse
import os

parser = argparse.ArgumentParser(description='Add gene symbols to GFF files.')
parser.add_argument('gff3', help='the gff3 file to be modified')
parser.add_argument('db', help='the gene2accession file.')
parser.add_argument('tax_id', help='the tax id of the species. Integers only.')
parser.add_argument('-o', '--output', help='Name of output file (Default: output_<datetime>.gff3)')
args = parser.parse_args()

id_fname = 'id_{}.txt'.format(datetime.now().strftime('%Y%m%d_%H%M%S'))
sliced_fname = 'extracted_{}.txt'.format(datetime.now().strftime('%Y%m%d_%H%M%S'))
with open(id_fname, 'w') as fout:
    fout.write("#tax_id\n" + args.tax_id)
print('Extracting the NCBI database...')
COMMAND_1 = "grep -w -F -f " + id_fname + " " + args.db + " > " + sliced_fname
sub.call(COMMAND_1, shell=True)
print('Done!')

print('Formatting refined extractions...')
df = pd.read_csv(sliced_fname, sep="\t")
idx = []
for index, row in df.iterrows():
    if row['RNA_nucleotide_accession.version'] == '-' or row['#tax_id'] != int(args.tax_id):
        idx.append(index)
df = df.drop(df.iloc[:, 4:15],axis=1).drop(df.iloc[:, 1:3],axis=1).drop(labels=idx, axis=0)
print('Done!')

print('Creating GFF database...')
gff3_file = args.gff3
db_name = 'gff3_{}.db'.format(datetime.now().strftime('%Y%m%d_%H%M%S'))
db = gffutils.create_db(gff3_file, db_name, keep_order=True)
print('Done!')

print('Processing GFF database and extracted data...')
match_list = pd.DataFrame(columns=['match_id','gene_symbol'])
match_id = []
print('Searching IDs with gene symbols...')
for a in db.features_of_type('match'):
    match_id.append(a.attributes['match_id'][0][:len(a.attributes['match_id'][0])-2])
match_id = sorted(match_id)
df = df.sort_values('RNA_nucleotide_accession.version')
print('Matching the accession numbers with the symbol...')
for row in df.itertuples():
    if(row['RNA_nucleotide_accession.version'][:len(row['RNA_nucleotide_accession.version'])-2] in match_id):
        match_list = match_list.append({'match_id': row['RNA_nucleotide_accession.version'], 'gene_symbol': row['Symbol']}, ignore_index=True)
match_list['match_id'] = match_list['match_id'].str[:-2]
match_list = match_list.set_index('match_id')
print('Done!')
            
output_file = args.output if args.output is not None else 'output_{}.gff3'.format(
    datetime.now().strftime('%Y%m%d_%H%M%S'))
print('Done!')
print('Sending output to file...')
with open(output_file, 'w') as fout:
    for d in db.directives:
        fout.write('##{0}\n'.format(d))

    for feature in db.all_features():
        if feature.featuretype == 'match':

            # sanity checking
            f_id = feature.attributes['match_id'][0][:-2]
            try:
                feature.attributes['symbol'] = match_list.loc[f_id]['gene_symbol']
            except:
                fout.write(str(feature) + '\n')
                continue
        fout.write(str(feature) + '\n')
print('Done!')
os.remove(db_name)
os.remove(sliced_fname)
os.remove(id_fname)
