"""
Simple snakefile to download and extract part of the GTEx public expression
data.

Usage:
    `snakemake` for the default subset.
    `snakemake gtex_<tissue>_top<N>.ft` for a custom subset.

    Valid tissues:

    ['Whole_Blood',
     'Brain_-_Frontal_Cortex_(BA9)',
     'Adipose_-_Subcutaneous',
     'Muscle_-_Skeletal',
     'Artery_-_Tibial',
     'Artery_-_Coronary',
     'Heart_-_Atrial_Appendage',
     'Adipose_-_Visceral_(Omentum)',
     'Ovary',
     'Uterus',
     'Vagina',
     'Breast_-_Mammary_Tissue',
     'Skin_-_Not_Sun_Exposed_(Suprapubic)',
     'Minor_Salivary_Gland',
     'Brain_-_Cortex',
     'Adrenal_Gland',
     'Thyroid',
     'Lung',
     'Spleen',
     'Pancreas',
     'Esophagus_-_Muscularis',
     'Esophagus_-_Mucosa',
     'Esophagus_-_Gastroesophageal_Junction',
     'Stomach',
     'Colon_-_Sigmoid',
     'Small_Intestine_-_Terminal_Ileum',
     'Colon_-_Transverse',
     'Prostate',
     'Testis',
     'Skin_-_Sun_Exposed_(Lower_leg)',
     'Nerve_-_Tibial',
     'Heart_-_Left_Ventricle',
     'Pituitary',
     'Brain_-_Cerebellum',
     'Cells_-_Cultured_fibroblasts',
     'Artery_-_Aorta',
     'Cells_-_EBV-transformed_lymphocytes',
     'Brain_-_Cerebellar_Hemisphere',
     'Brain_-_Caudate_(basal_ganglia)',
     'Brain_-_Nucleus_accumbens_(basal_ganglia)',
     'Brain_-_Putamen_(basal_ganglia)',
     'Brain_-_Hypothalamus',
     'Brain_-_Spinal_cord_(cervical_c-1)',
     'Liver',
     'Brain_-_Hippocampus',
     'Brain_-_Anterior_cingulate_cortex_(BA24)',
     'Brain_-_Substantia_nigra',
     'Kidney_-_Cortex',
     'Brain_-_Amygdala',
     'Cervix_-_Ectocervix',
     'Fallopian_Tube',
     'Cervix_-_Endocervix',
     'Bladder',
     'Kidney_-_Medulla',
     'Cells_-_Leukemia_cell_line_(CML)']


     To import the resulting dataframe in python:
     ```
     import pandas as pd

     # As a dataframe:
     pd.read_feather('gtex_Whole_Blood_top500.ft')

     # As a numpy array of just the numerical matrix:
     pd.read_feather('gtex_Whole_Blood_top500.ft').drop('Name', 1).to_numpy()
     ```
"""
import io
import gzip
import requests
import pandas as pd
import pyarrow as pa
import pyarrow.csv as pacsv
import pyarrow.compute as pac
from tqdm.auto import tqdm


# Constants

## Url of public RNA-seq counts
URL = 'https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz'
## Url of corresponding meta data, mostly tissue type
URL_META = 'https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt'

## Buffer size for streaming raw RNA-seq
BUF_SIZE = 10 * 2**20

# Default subset
TISSUE = 'Whole Blood'
TOP_N = 500

rule all:
    input: "gtex_" + TISSUE.replace(' ', '_') + "_top" + str(TOP_N) + ".ft"
 
class ResFile:
    """
    Helper class to interface request and pyarrow
    """
    def __init__(self, res, chunk_size=BUF_SIZE):
        self.it = res.iter_content(chunk_size=chunk_size)
        self.buf = b''
    def read(self, size):
        # Feed if starved
        if not self.buf:
            try:
                self.buf = next(self.it)
            except StopIteration:
                return b''
        buf = self.buf
        # cut if more than asked
        if size <= len(buf):
            self.buf = buf[size:]
            return buf[:size]
        # else return current buf
        self.buf = ''
        return buf
    @property
    def closed(self):
        return False

rule download_gex:
    output: 'gtex.ft'
    run:
        # Start HTTP request
        res = requests.get(URL, stream=True)
        # Wrap response in file object with gz decompression
        rf = ResFile(res)
        text = gzip.open(rf)

        # Build Arrow csv reader
        reader = pacsv.open_csv(
            text,
            read_options=pacsv.ReadOptions(skip_rows=2, block_size=BUF_SIZE),
            parse_options=pacsv.ParseOptions(delimiter='\t'),
        )

        # Get batch #1 and infer schema
        batch = reader.read_next_batch()

        # Initialize writer
        with pa.RecordBatchFileWriter(output[0], batch.schema,
                    options=pa.ipc.IpcWriteOptions(compression='lz4')
                    ) as writer:
            writer.write(batch)

            i = 1
            while True:
                if not i % 10:
                    print('Imported:', i, 'chunks')
                try:
                    batch = reader.read_next_batch()
                    writer.write(batch)
                    #chunk = reader.get_chunk(1000)
                    #table = pa.Table.from_pandas(chunk)
                    #writer.write(table)
                except StopIteration:
                    break
                i += 1

        text.close()
        res.close()


rule filter_gex:
    input: "gtex.ft"
    output: "gtex_{tissue}_top{top_n}.ft"
    run:
        # Load meta data
        pmeta = pd.read_csv(URL_META, sep='\t').set_index('SAMPID')

        # Open dataset
        readback = pa.ipc.open_file(input[0])

        # Targeted subset
        tissue = wildcards.tissue.replace('_', ' ')
        top_n = int(wildcards.top_n)

        # Filter samples ids by tissue
        tissue_names = [
            field.name
            for i, field in enumerate(readback.schema)
            if field.name in pmeta.index
            if pmeta.loc[field.name]['SMTSD'] == tissue
            ]
       
        # Scan dataset one, compute median expressions
        print('Scan 1: computing medians')
        medians = pd.concat([
            readback
                .get_batch(i)
                .to_pandas()
                .set_index('Name')
                [tissue_names]
                .median(axis=1)
            for i in tqdm(range(readback.num_record_batches))
        ])


        # Define threshold for inclusion
        thr = medians.sort_values(ascending=False)[top_n]

        # Define second scan
        def filter_batch(batch):
            batch = batch.to_pandas().set_index('Name')[tissue_names]
            return pa.Table.from_pandas(batch.loc[
                        batch
                        .median(axis=1)
                        > thr
                    ].reset_index())

        
        # Filter first batch to get schema
        b0 = filter_batch(readback.get_batch(0))

        # Create sub dataset and etl filtered dataset to it
        print('Scan 2: building subset')
        with pa.RecordBatchFileWriter(output[0], b0.schema, 
                    options=pa.ipc.IpcWriteOptions(compression='lz4')) as writer:
            for i in tqdm(range(readback.num_record_batches)):
                writer.write(filter_batch(readback.get_batch(i)))

