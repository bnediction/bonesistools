import bonesistools as bt

GENE_IDS = [
    "20375",
    "78688",
    "17869",
    "435318",
    "18033",
    "12394",
    "12393",
    "12305",
    "22059",
]

GENE_NAMES = [
    "PU.1",
    "B430311C09Rik",
    "Myc",
    "Gm5657",
    "NF-kappaB",
    "Aml1",
    "Aml3",
    "Ddr1",
    "Tp53",
]

NCBI_NAMES = [
    "Spi1",
    "Nol3",
    "Myc",
    "Serpina3d-ps",
    "Nfkb1",
    "Runx1",
    "Runx2",
    "Ddr1",
    "Trp53",
]


def test_mouse_gene_synonyms_convert_gene_ids_to_ncbi_names():

    genesyn = bt.dbs.ncbi.GeneSynonyms(organism="mouse")

    converted = genesyn.convert_sequence(
        GENE_IDS,
        input_identifier_type="gene_id",
        output_identifier_type="ncbi_name",
    )

    assert converted == NCBI_NAMES


def test_mouse_gene_synonyms_convert_gene_names_to_ncbi_names():

    genesyn = bt.dbs.ncbi.GeneSynonyms(organism="mouse")

    converted = genesyn.convert_sequence(
        GENE_NAMES,
        input_identifier_type="name",
        output_identifier_type="ncbi_name",
    )

    assert converted == NCBI_NAMES
