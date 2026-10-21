/*
 A nextflow workflow to reproduce the paper described in this repository

 Mark Stenglein 10/20/2026 
 */

include { process_metadata } from './workflows/process_metadata.nf'
include { align_and_tree   } from './workflows/align_and_tree.nf'
include { collapse_cds     } from './workflows/collapse_cds.nf'

// the main entry workflow
workflow {

  main:

    // create little channels with input files
    genbank_input = Channel.of(
    [ file("${projectDir}/input/all_reptarenavirus_L_seqs_2025_10_31.gb"), "L" ],
    [ file("${projectDir}/input/all_reptarenavirus_S_seqs_2025_10_31.gb"), "S" ] )
    
    // some Metadata about snakes in Hepojoki et al 2015
    hepojoki_metadata = Channel.fromPath("${projectDir}/input/Hepojoki_2015_table_1.xlsx", checkIfExists: true)

    // The ICTV MSL
    ictv_msl = Channel.fromPath("${projectDir}/input/ICTV_Master_Species_List_2025_MSL41.v1.xlsx", checkIfExists: true)
     
    // run a workflow to collect and wrangle metadata
    process_metadata(genbank_input)

    // create multiple sequence alignments and trees
    align_and_tree_input_ch = Channel.fromPath("${projectDir}/input/*_{CDS,prot}.fasta", checkIfExists: true )

    align_and_tree(align_and_tree_input_ch)

    // collapse NP and L sequences at 80% identity threshold using cd-hit
    cds_input = Channel.fromPath("${projectDir}/input/*_CDS.fasta", checkIfExists: true)
    cd_hit_cutoff = Channel.value("0.8")
    collapse_cds(cds_input, cd_hit_cutoff)

    // run code that generates figures and programmatically-generated values
    generate_figures_and_text(hepojoki_metadata, 
                              ictv_msl,
                              process_metadata.out.metadata.collect(), 
                              align_and_tree.out.treefile.collect(),
                              align_and_tree.out.alignment.collect(),
                              collapse_cds.out.fasta.collect() 
                              )


    // run code to generate the .doc format paper 
    Channel.of( 
     tuple(
        file(params.qmd), 
        file(params.refs), 
        file(params.csl), 
        file(params.custom_doc)
     ) ).set { paper_files_ch }

    edited_figures_ch = Channel.fromPath(params.edited_figures, checkIfExists: true)
    render_paper(paper_files_ch,
                 edited_figures_ch,
                 generate_figures_and_text.out.figures,
                 generate_figures_and_text.out.tables,
                 generate_figures_and_text.out.rds )

  publish:
    paper    = render_paper.out.paper

    figures  = generate_figures_and_text.out.figures
    tables   = generate_figures_and_text.out.tables
    rds      = generate_figures_and_text.out.rds   

    metadata      = process_metadata.out.metadata

    collapsed_cds = collapse_cds.out.fasta

    alignment   = align_and_tree.out.alignment
    contree     = align_and_tree.out.contree
    iqtree      = align_and_tree.out.iqtree
    treefile    = align_and_tree.out.treefile
    treelog     = align_and_tree.out.treelog

}

// workflow main output files 
output {
    // main output of quarto rendering
    paper {
        mode 'link'
    }

    // outputs of figure/table/text generation
    figures {
        path 'figures'
        mode 'link'
    }
    tables {
        path 'tables'
        mode 'link'
    }
    rds {
        path 'RDS'
        mode 'link'
    }

    // outputs of process_metadata
    metadata {
        path 'metadata'
        mode 'link'
    }

    // outputs of collapse_cds
    collapsed_cds {
        path 'collapsed_cds'
        mode 'link'
    }

    // outputs of align_and_tree
    alignment {
        path 'alignments'
        mode 'link'
    }
    contree {
        path 'trees'
        mode 'link'
    }
    iqtree {
        path 'trees'
        mode 'link'
    }
    treefile {
        path 'trees'
        mode 'link'
    }
    treelog {
        path 'trees'
        mode 'link'
    }

}

process generate_figures_and_text {

  // this custom docker/singularity container includes
  // necessary tools (R, R packages, etc)
  container 'ghcr.io/stenglein-lab/reptarenavirus_taxonomy:1.0'

  input:
    path hepojoki_metadata
    path ictv_msl

    // include all of the rest of these files so they will be staged in the pwd of the R script
    path metadata_files
    path tree_files
    path msa_files
    path collapsed_cds_files

  output:
    path "RDS",     emit: rds    
    path "figures", emit: figures
    path "tables",  emit: tables

  script:
  """
    # run R script that does main analysis heavy-lifting 
    # this script is in the bin dir so nextflow will be able to find it
    # and it has executable permissions so can be run 
    # it has an Rscript shebang
    reptarenavirus_taxonomy.R $hepojoki_metadata $ictv_msl
  """
}

process render_paper {

  // this custom docker/singularity container includes
  // necessary tools (R, R packages, etc)
  container 'ghcr.io/stenglein-lab/reptarenavirus_taxonomy:1.0'

  input:
    tuple path (qmd), path (refs), path (csl), path (custom_doc)
    // including the following paths as inputs to this process 
    // will cause the corresponding directories to be available
    // as relative paths to the .qmd document
    path edited_figures
    path figures
    path tables
    path rds

  output:
    path "paper.docx", emit: paper

  script:
  """
    # render quarto markdown file to docx
    quarto render $qmd --to docx -P baseDir:"$projectDir"
  """
}


