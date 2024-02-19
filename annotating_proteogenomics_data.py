import rdflib as rl
import pandas as pd
from rdflib.namespace import FOAF, RDF

## load ontology
g = rl.Graph()
g.parse("C:/Users/Finder/Downloads/ImPO_new.owl",format="xml")
impo = "https://raw.githubusercontent.com/liseda-lab/ImPO/main/ImPO.owl#"
placeholder = rl.Namespace("http://placeholder.com/")


#######################################################
    ##################################################
    # proteogenomics DATA FILES INTEGRATION          #
    ##################################################
#######################################################

## load proteogenomics data files - Ashwin
epitope_contig = pd.read_csv("C:/Users/Finder/Downloads/24681513/EpitopeContigs.tsv", sep="\t", index_col=None)
immunopeptide = pd.read_csv("C:/Users/Finder/Downloads/24681513/immunopeptides_library_new.tsv", sep="\t", index_col=None)
mappings = pd.read_csv("C:/Users/Finder/Downloads/24681513/mappings_new.tsv", sep="\t", index_col=None)
metadata = pd.read_csv("C:/Users/Finder/Downloads/24681513/dataset_metadata.csv", sep=";", index_col=None)
epitope_pep = pd.read_csv("C:/Users/Finder/Downloads/24681513/contigs_peptides.tsv", sep="\t", index_col=None)

#######################################################
#######################################################
dict_study = {}

for k in range(len(metadata["S.No"])): # add sample individual uris
    study_uri = rl.URIRef(placeholder + "study_" + str(metadata["S.No"][k])) ## study individuals uris
    g.add((study_uri, RDF.type, rl.URIRef(impo + "study"))) # adds as individual of class study
    g.add((study_uri, rl.URIRef(impo + "has_id"), rl.Literal(metadata["S.No"][k]))) # adds study id
    dataset_uri = rl.URIRef(placeholder + "dataset_" + metadata["Dataset_ID"][k]) ## dataset individuals uris
    g.add((dataset_uri, RDF.type, rl.URIRef(impo + "dataset"))) # adds as individual of class dataset
    g.add((dataset_uri, rl.URIRef(impo + "has_id"), rl.Literal(metadata["Dataset_ID"][k]))) # adds dataset id
    g.add((dataset_uri, rl.URIRef(impo + "has_publication"), rl.Literal(metadata["Source publication"][k]))) # adds Source publication to dataset
    g.add((dataset_uri, rl.URIRef(impo + "has_database"), rl.Literal(metadata["Source"][k]))) # adds database to dataset
    g.add((study_uri, rl.URIRef(impo + "has_outcome"), dataset_uri)) # adds dataset as outcome of study
    dict_study[metadata["Dataset_ID"][k]] = []
    dict_study[metadata["Dataset_ID"][k]].append(str(metadata["S.No"][k]))
    dict_study[metadata["Dataset_ID"][k]].append(metadata["MS instrument"][k])
print("------ after study and dataset instantiation ------")

dict_spec = {}
sample_uri = rl.URIRef(placeholder)
for k in range(len(immunopeptide["study ID"])):
    if dict_study.get(immunopeptide["study ID"][k]) != None:
        study_uri = rl.URIRef(placeholder + str(dict_study[immunopeptide["study ID"][k]][0])) ## study individuals uris
        if rl.URIRef(placeholder + immunopeptide["sample ID"][k]) != sample_uri:
            #sample_uri = rl.URIRef(placeholder + immunopeptide["sample ID"][k]) ## sample individuals uris
            sample_uri = rl.URIRef(placeholder + "sample_" + str(k)) ## sample individuals uris
            g.add((sample_uri, RDF.type, rl.URIRef(impo + "sample"))) # adds as individual of class sample
            g.add((sample_uri, rl.URIRef(impo + "has_id"), rl.Literal(immunopeptide["sample ID"][k]))) # adds sample id
            assay_uri = rl.URIRef(placeholder + "assay_" + immunopeptide["sample ID"][k]) ## make assay uri per sample
            g.add((assay_uri, RDF.type, rl.URIRef(impo + "mass_spectrometry"))) # added as instance of mass spec because file has only for spec data
            g.add((assay_uri, rl.URIRef(impo + "has_id"), rl.Literal(immunopeptide["sample ID"][k]))) # adds an internal assay id
            g.add((assay_uri, rl.URIRef(impo + "part_of"), study_uri)) # change to add assay as is_part of study; and assay has_input sample (how to model participation constraint?)
            g.add((assay_uri, rl.URIRef(impo + "has_input"), sample_uri)) # adds sample as input to assay
            g.add((assay_uri, rl.URIRef(impo + "has_instrument"), rl.Literal(dict_study[immunopeptide["study ID"][k]][1]))) # adds an internal assay id
            for j in immunopeptide["Spectrum"][k].split(","): # because 1 sample has more than one assay
                spectrum_uri = rl.URIRef(placeholder + "spectrum_" + j) ## make spectrum uri per MS assay
                g.add((spectrum_uri, RDF.type, rl.URIRef(impo + "spectrum"))) # added as instance of spectrum
                g.add((spectrum_uri, rl.URIRef(impo + "has_name"), rl.Literal(j))) # adds spectrum name
                g.add((assay_uri, rl.URIRef(impo + "has_output"), spectrum_uri))
                if dict_spec.get(immunopeptide["Peptide sequence"][k]) == None: dict_spec[immunopeptide["Peptide sequence"][k]] = []
                dict_spec[immunopeptide["Peptide sequence"][k]].append(j)
del dict_study

print("------ after spectrum instantiation ------")


for k in range(len(immunopeptide["sample ID"])):
    assay_uri = rl.URIRef(placeholder + "assay_" + immunopeptide["sample ID"][k])
    peptide_uri = rl.URIRef(placeholder + "peptide_" + immunopeptide["Peptide sequence"][k])  # + str(count_pep_id))
    g.add((peptide_uri, RDF.type, rl.URIRef(impo + "peptide")))
    g.add((peptide_uri, rl.URIRef(impo + "has_sequence"), rl.Literal(immunopeptide["Peptide sequence"][k])))
    g.add((peptide_uri, rl.URIRef(impo + "has_tag"), rl.Literal(immunopeptide["tag"][k])))
    g.add((assay_uri, rl.URIRef(impo + "identifies"), peptide_uri))


#### PEPTIDE - SPECTRUM IDENTIFICATION  - W/ properties score and analysis strategy
for k in range(len(immunopeptide["Peptide sequence"])):
    for j in dict_spec[immunopeptide["Peptide sequence"][k]]:
        spec_pep_identification = rl.URIRef(placeholder + "spec_pep_id_" + immunopeptide["Peptide sequence"][k] + "_" + j)
        g.add((spec_pep_identification, RDF.type, rl.URIRef(impo + "mass_spectrometry-peptide_identification")))
        g.add((spec_pep_identification, rl.URIRef(impo + "has_source"), rl.URIRef(placeholder + "spectrum_" + j)))
        g.add((spec_pep_identification, rl.URIRef(impo + "has_target"), rl.URIRef(placeholder + "peptide_" + immunopeptide["Peptide sequence"][k])))
        g.add((spec_pep_identification, rl.URIRef(impo + "has_analysis_strategy"), rl.Literal(immunopeptide["mass spectrometry analysis strategy"][k])))
del dict_spec
print("------ after assignment ------")

## in file metadata, link samples to source tissue individual through the study/dataset sample is input to
for k in range(len(immunopeptide["sample ID"])):
    for j in range(len(metadata["Dataset_ID"])):
        if immunopeptide["sample ID"][k] == metadata["Dataset_ID"][j]:
            sample_uri = rl.URIRef(placeholder + immunopeptide["sample ID"][k])
            tissue_uri = rl.URIRef(placeholder + metadata["Tissue"][j]) ## make assay uri per sample  
            if metadata["Disease_type"][j] == "Cancer": 
                g.add((tissue_uri, RDF.type, rl.URIRef(impo + "tumor_tissue")))
                cancer_uri = rl.URIRef(placeholder + metadata["Disease_type"][j]) ## make assay uri per sample
                g.add((tissue_uri, rl.URIRef(impo + "part_of"), cancer_uri))
            else: g.add((tissue_uri, RDF.type, rl.URIRef(impo + "normal_tissue")))
            g.add((sample_uri, rl.URIRef(impo + "has_source"), tissue_uri))
            g.add((sample_uri, rl.URIRef(impo + "has_type"), rl.Literal(metadata["Source_Tissue"][k])))



#### INSTANCIATE EPITOPE CONTIG
for k in range(len(epitope_contig["id"])):
    contig_uri = rl.URIRef(placeholder + "contig_" + epitope_contig["id"][k])
    g.add((contig_uri, RDF.type, rl.URIRef(impo + "epitope_contig")))
    g.add((contig_uri, rl.URIRef(impo + "has_id"), rl.Literal(epitope_contig["id"][k])))
    g.add((contig_uri, rl.URIRef(impo + "has_number_of_unique_peptides"), rl.Literal(epitope_contig["unique peptides"][k])))
    g.add((contig_uri, rl.URIRef(impo + "has_population_coverage"), rl.Literal(epitope_contig["population coverage"][k])))
    g.add((contig_uri, rl.URIRef(impo + "has_conservation_score"), rl.Literal(epitope_contig["conservation score"][k])))
    g.add((contig_uri, rl.URIRef(impo + "has_hotspot_score"), rl.Literal(epitope_contig["hotspot score"][k])))
    g.add((contig_uri, rl.URIRef(impo + "has_max_depth"), rl.Literal(epitope_contig["max depth"][k])))
    g.add((contig_uri, rl.URIRef(impo + "has_mean_depth"), rl.Literal(epitope_contig["mean depth"][k])))
    g.add((contig_uri, rl.URIRef(impo + "has_mutation_frequency"), rl.Literal(epitope_contig["mutation frequency"][k])))
    g.add((contig_uri, rl.URIRef(impo + "has_spectral_counts"), rl.Literal(epitope_contig["spectral counts"][k])))

## EPITOPE CONTIG <belong_to> PEPTIDE
for k in range(len(epitope_pep["peptide"])):
    g.add((rl.URIRef(placeholder + "contig_" + epitope_pep["contig"][k]), rl.URIRef(impo + "belong_to"), rl.URIRef(placeholder + "peptide_" + epitope_pep["peptide"][k])))


print("------ after contig instantiation ------")

#### INSTANCIATE AND LINK PROTEIN <product_of> TRANSCRIPT <product_of> GENE
for k in range(len(mappings["gene name"])):
    protein_uri = rl.URIRef(placeholder +  mappings["protein id"][k])
    g.add((protein_uri, RDF.type, rl.URIRef(impo + "protein")))
    g.add((protein_uri, rl.URIRef(impo + "has_id"), rl.Literal(mappings["protein id"][k])))

    transcript_uri = rl.URIRef(placeholder + mappings["transcript id"][k])
    g.add((transcript_uri, RDF.type, rl.URIRef(impo + "transcript")))
    g.add((transcript_uri, rl.URIRef(impo + "has_id"), rl.Literal(mappings["transcript id"][k])))
    g.add((transcript_uri, rl.URIRef(impo + "has_name"), rl.Literal(mappings["transcript name"][k])))

    gene_uri = rl.URIRef(placeholder + mappings["gene id"][k])
    g.add((gene_uri, RDF.type, rl.URIRef(impo + "gene")))
    g.add((gene_uri, rl.URIRef(impo + "has_id"), rl.Literal(mappings["gene id"][k])))
    g.add((gene_uri, rl.URIRef(impo + "has_name"), rl.Literal(mappings["gene name"][k])))

    g.add((protein_uri, rl.URIRef(impo + "product_of"), transcript_uri))
    g.add((transcript_uri, rl.URIRef(impo + "product_of"), gene_uri))

print("------ after all entities' instantiation ------")


#### PEPTIDE - PROTEIN ASSIGNMENT
dict_peptide_pos = {}
for k in range(len(immunopeptide["Peptide sequence"])):
    dict_peptide_pos[immunopeptide["Peptide sequence"][k]] = [immunopeptide["start"][k],immunopeptide["end"][k]]

for k in range(len(mappings["peptide"])):
    protein_uri = rl.URIRef(placeholder +  mappings["protein id"][k])
    peptide_uri = rl.URIRef(placeholder + "peptide_" + mappings["peptide"][k])
    pep_prot_assignment = rl.URIRef(placeholder + "pep_prot_id_" + mappings["peptide"][k] + "_" + mappings["protein id"][k])
    g.add((pep_prot_assignment, RDF.type, rl.URIRef(impo + "peptide-protein_assignment")))
    g.add((pep_prot_assignment, rl.URIRef(impo + "has_source"), protein_uri))
    g.add((pep_prot_assignment, rl.URIRef(impo + "has_target"), peptide_uri))
    if dict_peptide_pos.get(mappings["peptide"][k]) != None:
        g.add((pep_prot_assignment, rl.URIRef(impo + "has_start_position"), rl.Literal(dict_peptide_pos[mappings["peptide"][k]][0])))
        g.add((pep_prot_assignment, rl.URIRef(impo + "has_end_position"), rl.Literal(dict_peptide_pos[mappings["peptide"][k]][1])))
del dict_peptide_pos


g.serialize("C:/Users/laura/Downloads/proteogenomics_impo.owl") # generates owl file

#### TO CHECK:
## USAGE OF COLUMNS IN FILE METADATA
## START AND END IN IMMUNOPEPTIDE_LIBRARY ARE FOR PEPTIDE - PRODUCT_OF - GENOMIC_REGION OR PEPTIDE - ASSIGNED_TO - PROTEIN ?