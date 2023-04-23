from biomart import BiomartServer
from os import mkdir, path

def main():
    print("Retrieving database information from ensembl.org...")

    servers = [ ("Ensembl server (asia)", BiomartServer("http://asia.ensembl.org/biomart"), 'default'),
                ("Ensembl plants server", BiomartServer("http://plants.ensembl.org/biomart"), 'plants_mart_26') ]

    datasets = []
    
    for servername, server, vsname in servers:
        print("Retrieving dataset information from {0}...".format(servername))
        for key in server.datasets:
            if 'eg_gene' in key or 'gene_ensembl' in key:
                datasets.append((server.datasets[key], vsname, key))

    print('')
    print("Available datasets:")
    for i, dataset in enumerate(datasets, 1):
        print(str(i)+". "+str(dataset[0]))

    dataset_id = int(raw_input("Select a dataset to create RGEN database: "))

    dataset = datasets[dataset_id-1]
    print('')
    print("Selected dataset: "+str(dataset[0]))
    print('')
    print("Retrieving gene/transcript/CDS information...")

    if 'external_gene_name' in dataset[0].attributes:
        attributes = [
                'ensembl_gene_id',
                'ensembl_transcript_id',
                'external_gene_name',
                'chromosome_name',
                'description',
                'genomic_coding_start',
                'genomic_coding_end',
            ]
    else:
        attributes = [
                'ensembl_gene_id',
                'ensembl_transcript_id',
                'external_gene_id',
                'chromosome_name',
                'description',
                'genomic_coding_start',
                'genomic_coding_end',
            ]
    response = dataset[0].search({ 'attributes': attributes }, vsname=dataset[1])
    if not path.isdir(dataset[2]):
        mkdir(dataset[2])
    with open(dataset[2]+"/gene_data", "w") as f:
        f.write("#"+'\t'.join(attributes)+"\n")
        for i, line in enumerate(response.iter_lines()): # Response is TSV
            entries = line.split('\t')
            if not '' in entries[5:]: # force results to have CDS informations
                f.write(line+"\n")
            if i%10000 == 0:
                print(i)

if __name__ == "__main__":
    main()
