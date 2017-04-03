import argparse
import click
from pathlib import Path

import requests

biomart_url = 'http://sep2015.archive.ensembl.org/biomart/martservice'

fields = ["ensembl_gene_id",
          "chromosome_name",
          "start_position",
          "end_position",
          "strand",
          "entrezgene",
          "hgnc_symbol",
          "hgnc_id"]

query_xml = '''<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "TSV" header = "1"
        uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
    <Dataset name = "hsapiens_gene_ensembl" interface = "default" >
''' + '\n'.join(['<Attribute name = "{}" />'.format(f) for f in fields]) + '''
    </Dataset>
</Query>'''


def parse_args():
    parser = argparse.ArgumentParser(
        description="Download ensemble reference"
    )
    parser.add_argument('out', help='Output tsv file.')
    args = parser.parse_args()

    # Check for existence for paths
    args.out = Path(args.out)
    args.out = args.out.parent.resolve() / args.out.name

    return args


if __name__ == '__main__':
    args = parse_args()

    r = requests.get(biomart_url, params={'query': query_xml}, stream=True)
    r.raise_for_status()

    r_length = r.headers.get('content-length')
    if r_length is not None:
        r_length = int(r_length)

    with args.out.open('wb') as file:
        with click.progressbar(r.iter_content(1024), length=r_length) as bar:
            for chunk in bar:
                file.write(chunk)
