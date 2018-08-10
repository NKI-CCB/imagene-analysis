import click
import numpy as np
import pandas as pd


def only(x):
    assert len(x) == 1, f"Length is {len(x)}"
    return list(x)[0]


click_in_path = click.Path(exists=True, dir_okay=False, resolve_path=True)
click_out_path = click.Path(exists=False, dir_okay=False, resolve_path=True)


@click.command()
@click.argument('xls', type=click_in_path)
@click.argument('refseq', type=click_in_path)
@click.argument('ensembl', type=click_in_path)
@click.argument('out', type=click_out_path)
def map_genes(xls, refseq, ensembl, out):
    genes_df = pd.read_excel(xls)
    genes_df = genes_df.rename(columns={'#name': 'refseq_id',
                                        'name2': 'original_gene_name'})

    # Map Refseq identifier (NM_XXXX) to HGNC identifier
    refseq_annot = pd.read_csv(refseq, sep='\t')
    refseq_to_hgnc = dict()
    for refseq_id, hgnc_id in zip(refseq_annot['refseq_id'],
                                  refseq_annot['hgnc_id']):
        refseq_to_hgnc.setdefault(refseq_id, set()).add(hgnc_id)
    genes_df['hgnc_id'] = [only(refseq_to_hgnc.get(r, set([''])))
                           for r in genes_df['refseq_id']]

    # Map HGNC identifier to Ensembl identifier
    ensembl_annot = pd.read_csv(ensembl, sep='\t')
    hgnc_to_ensembl = dict()
    for ensembl_id, hgnc_id in zip(ensembl_annot['Ensembl Gene ID'],
                                   ensembl_annot['HGNC ID(s)']):
        if hgnc_id:
            hgnc_to_ensembl.setdefault(hgnc_id, set()).add(ensembl_id)
    ensembl_annot = ensembl_annot.set_index('Ensembl Gene ID')

    genes_df['ensembl_id'] = np.full(genes_df.shape[0], None)
    ensembl_col = genes_df.columns.get_loc('ensembl_id')
    for row_idx, hgnc_id in enumerate(genes_df['hgnc_id']):
        ensembl_ids = hgnc_to_ensembl.get(hgnc_id, set())
        if len(ensembl_ids) > 1:
            chroms = ensembl_annot['Chromosome Name'].loc[list(ensembl_ids)]
            ensembl_ids = set([e for e, c in chroms.items() if len(c) < 2])
        if len(ensembl_ids) == 0:
            genes_df.iloc[row_idx, ensembl_col] = ''
        else:
            genes_df.iloc[row_idx, ensembl_col] = only(ensembl_ids)

    # Output
    genes_df.to_csv(out, sep='\t', index=False)


if __name__ == '__main__':
    map_genes()
