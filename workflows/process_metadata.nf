/*

 A nextflow workflow to wrangle, extract, and fetch metadata 
 related to this repository's paper 

 Mark Stenglein 10/20/2026

 */

workflow process_metadata {

  take:
    genbank // [path_to_genbank, value: segment_ID]

  main:

    // extract metadata from input genbank files
    extract_metadata(genbank)
    combine_metadata(extract_metadata.out.metadata.collect())

    // download info about reptarenavirus taxa from NCBI
    fetch_reptarenavirus_taxa()

    // output ch with all generated files
    metadata_output = combine_metadata.out.metadata.concat(fetch_reptarenavirus_taxa.out.map)

  emit:
    metadata  = metadata_output
}

// pull metadata out of one genbank file
process extract_metadata {
  tag "$genbank"
  label 'process_low'

  // a singularity image with biopython
  container 'https://depot.galaxyproject.org/singularity/biopython:1.81'

  input:
    tuple path (genbank), val (segment)

  output:
    path "${segment}.accession_map.txt", emit: metadata

  script:
    """
    extract_taxonomy_from_genbank.py $genbank $segment > ${segment}.accession_map.txt
    """
}

// concatenate individual metadata files
process combine_metadata {
  label 'process_low'

  // a basic linux image for cat and awk
  container 'https://depot.galaxyproject.org/singularity/ubuntu:24.04' 

  input:
    path genbanks

  output:
    path "genbank_metadata.txt", emit: metadata

  script:
    """
    # concatenate individual files into one file (cat), removing all but first instance of header line (awk)
    cat $genbanks | awk '/^accession\t/ { if (seen++) next } 1' - > genbank_metadata.txt
    """
}

// this process uses R API to NCBI E-utils to pull down info
// about reptarenavirus taxa from NCBI 
process fetch_reptarenavirus_taxa {
  label 'process_low'

  // this custom docker/singularity container includes
  // necessary tools (R, R packages, etc)
  container 'ghcr.io/stenglein-lab/reptarenavirus_taxonomy:1.0'

  output:
    path "accession_species_taxid_map.txt", emit: map

  script:
  """
  fetch_reptarenavirus_taxa_accessions.R
  """
}



