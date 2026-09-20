#!/usr/bin/env python3
"""
Convert Liftoff GFF3 to featureCounts-compatible format
Adds gene_id and transcript_id attributes to all features
"""
import sys
import re

gff3_file = sys.argv[1] if len(sys.argv) > 1 else \
    "/Volumes/Expansion/NCBI/PRJNA934294/liftoff/Mcoruscus_HiC_liftoff.gff3"
output_file = sys.argv[2] if len(sys.argv) > 2 else \
    "/Volumes/Expansion/NCBI/PRJNA934294/liftoff/Mcoruscus_HiC_liftoff_fc.gff3"

print(f"Converting {gff3_file} -> {output_file}")

# First pass: build mRNA→gene mapping
mRNA_to_gene = {}
gene_ids = set()
with open(gff3_file) as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        if len(parts) < 9:
            continue
        attrs_str = parts[8]
        attrs = {}
        for item in attrs_str.split(';'):
            if '=' in item:
                k, v = item.split('=', 1)
                attrs[k.strip()] = v.strip()

        feat_type = parts[2]
        feat_id = attrs.get('ID', '')
        parent = attrs.get('Parent', '')

        if feat_type == 'gene':
            gene_ids.add(feat_id)
        elif feat_type == 'mRNA' and parent.startswith('gene-'):
            mRNA_to_gene[feat_id] = parent

# Second pass: add gene_id and transcript_id
out_lines = 0
with open(gff3_file) as f, open(output_file, 'w') as out:
    for line in f:
        if line.startswith('#'):
            out.write(line)
            continue

        parts = line.strip().split('\t')
        if len(parts) < 9:
            out.write(line)
            continue

        attrs_str = parts[8]
        attrs = {}
        for item in attrs_str.split(';'):
            if '=' in item:
                k, v = item.split('=', 1)
                attrs[k.strip()] = v.strip()

        feat_type = parts[2]
        feat_id = attrs.get('ID', '')
        parent = attrs.get('Parent', '')

        # Determine gene_id and transcript_id
        gene_id = ''
        transcript_id = ''

        if feat_type == 'gene':
            gene_id = feat_id
        elif feat_type == 'mRNA':
            gene_id = parent if parent.startswith('gene-') else ''
            transcript_id = feat_id
        elif feat_type in ('exon', 'CDS', 'five_prime_UTR', 'three_prime_UTR'):
            if parent.startswith('gene-'):
                gene_id = parent
            elif parent.startswith('rna-') and parent in mRNA_to_gene:
                gene_id = mRNA_to_gene[parent]
                transcript_id = parent
            else:
                gene_id = parent

        # Add gene_id and transcript_id
        new_attrs = []
        has_gene_id = False
        has_transcript_id = False
        for item in attrs_str.split(';'):
            item = item.strip()
            if '=' in item:
                k, v = item.split('=', 1)
                if k.strip() == 'gene_id':
                    has_gene_id = True
                    new_attrs.append(f'gene_id={gene_id}')
                    continue
                elif k.strip() == 'transcript_id':
                    has_transcript_id = True
                    new_attrs.append(f'transcript_id={transcript_id}')
                    continue
            new_attrs.append(item)

        if not has_gene_id and gene_id:
            new_attrs.append(f'gene_id={gene_id}')
        if not has_transcript_id and transcript_id:
            new_attrs.append(f'transcript_id={transcript_id}')

        parts[8] = ';'.join(new_attrs)
        out.write('\t'.join(parts) + '\n')
        out_lines += 1
        if out_lines % 100000 == 0:
            print(f"  Processed {out_lines} lines...")

print(f"Conversion complete: {out_lines} features written")
print(f"Genes: {len(gene_ids)}, mRNA-gene mappings: {len(mRNA_to_gene)}")
