
/*
 A nextflow workflow to extract metadata from some genbank files

 Mark Stenglein 11/3/2025
 */

workflow {

  main:
    // make alignments
    genbank_ch = Channel.fromPath(params.genbank)
    segment_ch = Channel.value(params.segment)
    extract_metadata(genbank_ch, segment_ch)

  publish:
    metadata     = extract_metadata.out.metadata
}

output {
    metadata {
        mode 'link'
    }
}

// pull metadata out of one genbank file
process extract_metadata {
  tag "$genbank"
  label 'process_low'

  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/biopython:1.81':
    'quay.io/biocontainers/biopython:1.81' }"

  input:
    path genbank
    val segment

  output:
    path "${segment}.accession_map.txt", emit: metadata

  script:
  """
  extract_taxonomy_from_genbank.py $genbank $segment > ${segment}.accession_map.txt
  """
}
