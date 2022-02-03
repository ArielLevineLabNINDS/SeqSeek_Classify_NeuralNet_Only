library(tidyverse)
library(Seurat)
library(Matrix)


mapReadGTF <- function(file,hasColumnNames){
  
  cnames = hasColumnNames
  if (!cnames) cnames=c("Chromosome","Source","Feature","Start","End","Score","Strand","Frame","Attributes") 
  
  gtf <- read_tsv(file,col_names = cnames, col_types = "ccccccccc", comment = '#') %>% 
    mutate(
      gene_id    = gsub("gene_id \"(ENSMUSG\\d+)\";.+","\\1",Attributes),
      gene_names = gsub(".+ gene_name \"([[:print:]]+?)\"; .+","\\1",Attributes) 
    ) %>% 
    select(gene_id,gene_names) %>% 
    distinct() %>%
    # Seurat checks that all genes are uniques. Gene ids that map to the same gene name
    # cause the count matrix to have duplicated gene names.  Seurat uses the make unique
    # function to make the gene names unique.  We're going to do the same.
    mutate(gene_names=make.unique(gene_names))
  
  return(gtf)
}

mapGenes <- function(fromGenome,ToGenome,count_matrix){

  
  # before we start make sure the genes in the count matrix match the genes in the fromGenome.
  if ( !all(rownames(count_matrix) %in% fromGenome$gene_names)){
    stop("The genes in the count matrix do not match the genes in the fromGenome.")
  }
    
  
  ## map via gene id's. so create a map from the "fromGenome" to the "toGenome"
  gene_map <- fromGenome %>% left_join(ToGenome,by="gene_id",suffix=c(".from",".to"))

  # if the gene_name.to is NA  then I don;t have a good match.  
  # one way to handle this is to remove these counts...
  ignore_unmapped_genes = TRUE  
  
  if (ignore_unmapped_genes){
    message("dropping ",sum(is.na(gene_map$gene_names.to))," genes")
    gene_map <- gene_map %>% filter(!is.na(gene_names.to))
    
    # remove all genes not in the map...
    row_indx = rownames(count_matrix) %in% gene_map$gene_names.from
    count_matrix = count_matrix[row_indx,]
  }
  
  # other ways to deal with unmapped genes???

    
  # ok ready to map the gene.
  # map the gene names...
  map_vector = pull(gene_map,gene_names.to,gene_names.from)
  rownames(count_matrix)<- recode(rownames(count_matrix),!!!map_vector)
  
  # now lets add the "missing" genes from the toGenome...
  missing_genes = setdiff(ToGenome$gene_names,rownames(count_matrix))
  zeros = Matrix(0,nrow=length(missing_genes),ncol=dim(count_matrix)[2],dimnames = list(missing_genes,colnames(count_matrix)),sparse = T )
  x=rbind(count_matrix,zeros)

  # now lets set the order of the matrix...
  x=x[ToGenome$gene_names,]
  
  return (x)
}

filter_genes <- function (x){
  genes <- read_csv("https://raw.githubusercontent.com/ArielLevineLabNINDS/Seq-Seek-classifyData/master/dirty_neurons_genes.csv.gz",col_names  = "genes") %>%
    pull(genes)
  
  return( x[genes,])
}



