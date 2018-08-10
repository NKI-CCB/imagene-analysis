import gzip

import click


def split_loci(lines):
    locus = None
    section_title = None
    section = None
    for l in lines:
        if l == '//\n':
            locus.append((section_title, section))
            yield locus
            locus = None
            section_title = None
            section = None
            continue
        if l.startswith('LOCUS'):
            locus = list()
        if not l.startswith(' '):
            if section is not None:
                locus.append((section_title, section))
            section_title = l[0:12].strip()
            section = list()
        section.append(l)
    if locus is not None:
        yield locus


def split_features(lines):
    header = next(lines)
    assert header == 'FEATURES             Location/Qualifiers\n'
    feature_title = None
    feature = None
    for l in lines:
        assert l.startswith(' ' * 5)
        l = l[5:]
        title = l[0:16].strip()
        if title:
            if feature is not None:
                yield (feature_title, feature)
            feature_title = title
            feature = list()
        feature.append(l[16:].strip())
    yield (feature_title, feature)


click_in_path = click.Path(exists=True, dir_okay=False, resolve_path=True)
click_out_path = click.Path(exists=False, dir_okay=False, resolve_path=True)


@click.command()
@click.argument('gbff', type=click_in_path)
@click.argument('out', type=click_out_path)
def parse_gbff(gbff, out, sep="\t"):
    with open(out, 'w') as out:
        out.write("refseq_id")
        out.write(sep)
        out.write("hgnc_id")
        out.write("\n")
        with gzip.open(gbff, 'rt') as f:
            loci = split_loci(f)
            for l in loci:
                for s, sc in l:
                    if s == 'LOCUS':
                        assert len(sc) == 1
                        sc = sc[0].split()
                        refseq_id = sc[1]
                    if s == 'FEATURES':
                        features = split_features(iter(sc))
                        for f, fc in features:
                            if f == 'gene':
                                for fl in fc:
                                    if fl.startswith('/db_xref="HGNC:'):
                                        hgnc_id = fl[15:-1]
                                        out.write(refseq_id)
                                        out.write(sep)
                                        out.write(hgnc_id)
                                        out.write("\n")


if __name__ == '__main__':
    parse_gbff()
