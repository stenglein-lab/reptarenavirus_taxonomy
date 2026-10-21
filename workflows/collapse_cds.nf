    
// run cd-hit-est to collapse coding sequences at a particular identity threshold  
process collapse_cds {
    
  container 'https://depot.galaxyproject.org/singularity/cd-hit:4.8.1--h5ca1c30_13'
    
  input:
    path fasta
    val threshold
    
  output:
    path "*.fasta", emit: fasta
    
  script:
    def fasta_base = fasta.baseName
    """
    cd-hit-est -c $threshold -i $fasta -o ${fasta_base}.cdhit_${threshold}.fasta
    """
}   
