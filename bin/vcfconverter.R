#!/usr/bin/env Rscript

### Transformando dados de VCF em tabela para analise de cluster hierarquico ###

### Prep baseado na frequência de alelos dos libidinosus ou de qualquer subamostragem (ver funcao genotypes)
library(vcfR)
library(dplyr)

args<-commandArgs(trailingOnly=TRUE)

vcfp<-args[1]
ind_to_freq<-args[2]

myvcf<-read.vcfR(vcfp)

twocol<-function(vcfpath, threecols=T){
  phy2<-read.vcfR(
    vcfpath,
    limit = 1e+07,
    nrows = -1,
    skip = 0,
    cols = NULL,
    convertNA = TRUE,
    checkFile = TRUE,
    check_keys = TRUE,
    verbose = TRUE )
  
  
  SNPs_phy2<-extract.haps(phy2, mask = FALSE, unphased_as_NA = FALSE, verbose = TRUE)
  SNPs_phy2<-t(SNPs_phy2)
  
  phy2tidy<-vcfR2tidy(phy2)
  isolates<- unique(phy2tidy$gt$Indiv)
  
  #Extraindo rows impares e pares:
  parouimpar<-seq_len(nrow(SNPs_phy2)) %% 2
  data_odd <- SNPs_phy2[parouimpar == 1, ]
  colnames(data_odd)<-as.matrix(paste(colnames(data_odd), "_0",sep = "")) #renomeando colunas para juntar
  rownames(data_odd)<-as.matrix(isolates)
  
  data_even <- SNPs_phy2[parouimpar == 0, ]
  colnames(data_even)<-as.matrix(paste(colnames(data_even), "_1",sep = ""))
  rownames(data_even)<-as.matrix(isolates)
  
  
  #Construindo matriz com duas colunas por locu:
  
  test<-cbind(data_odd, data_even)
  col.order <- sort(colnames(test))
  test<-test[,col.order]
  
  if (threecols==F){return(test)}
  
  if (threecols==T){
    #Adicionar coluna de genotipos (para representação com uma coluna)
    data_genotype <- SNPs_phy2[parouimpar == 0, ]
    colnames(data_genotype)<-as.matrix(paste(colnames(data_genotype), "_2",sep = ""))
    rownames(data_genotype)<-as.matrix(isolates)
    data_genotype[,]<-NA
    
    #Construindo matriz com uma colunas por locu:
    test<-cbind(data_odd, data_even, data_genotype)
    col.order <- sort(colnames(test))
    test<-test[,col.order]
    
    return(test)
  }
}

matriz<-twocol(vcfp, threecols = T)

inds<-read.table(ind_to_freq)
inds
ind_freqed<-matriz[row.names(matriz)%in%as.vector(inds[,1]),]

genotypes<-function(data=NA, freq=NA){
  #Nomeando com base no genótipo:
  for (i in 1:length(colnames(data))){
    if (endsWith(colnames(data)[i], "2")){
      
      unicos_l1<-unique(freq[,(i-1)])
      unicos_l1<-sort(unicos_l1[!is.na(unicos_l1)])
      unicos_l2<-unique(freq[,(i-2)])
      unicos_l2<-sort(unicos_l2[!is.na(unicos_l2)])
      
      if (length(unicos_l1)==2 & length(unicos_l2)==2){
        l1_al1<-sum(as.numeric(freq[,(i-1)]==unicos_l1[1]),na.rm = T)
        l2_al1<-sum(as.numeric(freq[,(i-2)]==unicos_l2[1]),na.rm = T)
        al1f<-l1_al1+l2_al1
        
        l1_al2<-sum(as.numeric(freq[,(i-1)]==unicos_l1[2]),na.rm = T)
        l2_al2<-sum(as.numeric(freq[,(i-2)]==unicos_l2[2]),na.rm = T)
        al2f<-l1_al2+l2_al2
        
        
        freqs<-c(al1f,al2f)
      }
      
      if (length(unicos_l1)==2 & length(unicos_l2)==1){if (unicos_l2[1]==unicos_l1[2]){
        l1_al2<-sum(as.numeric(freq[,(i-1)]==unicos_l1[2]),na.rm = T)
        l2_al1<-sum(as.numeric(freq[,(i-2)]==unicos_l2[1]),na.rm = T)
        al2f<-l1_al2+l2_al1
        
        l1_al1<-sum(as.numeric(freq[,(i-2)]==unicos_l1[1]),na.rm = T)
        al1f<-l1_al1
        
        freqs<-c(al1f,al2f)
      }
        if (unicos_l2[1]==unicos_l1[1]){
          l1_al1<-sum(as.numeric(freq[,(i-1)]==unicos_l1[1]),na.rm = T)
          l2_al1<-sum(as.numeric(freq[,(i-2)]==unicos_l2[1]),na.rm = T)
          al1f<-l1_al1+l2_al1
          
          l1_al2<-sum(as.numeric(freq[,(i-2)]==unicos_l1[2]),na.rm = T)
          al2f<-l1_al2
          
          freqs<-c(al1f,al2f)
        }
      }
      
      if (length(unicos_l2)==2 & length(unicos_l1)==1){if (unicos_l1[1]==unicos_l2[2]){
        l1_al1<-sum(as.numeric(freq[,(i-1)]==unicos_l1[1]),na.rm = T)
        l2_al2<-sum(as.numeric(freq[,(i-2)]==unicos_l2[2]),na.rm = T)
        al2f<-l1_al1+l2_al2
        
        l2_al1<-sum(as.numeric(freq[,(i-2)]==unicos_l2[1]),na.rm = T)
        al1f<-l2_al1
        
        freqs<-c(al1f,al2f)
      }
        if (unicos_l1[1]==unicos_l2[1]){
          l1_al1<-sum(as.numeric(freq[,(i-1)]==unicos_l1[1]),na.rm = T)
          l2_al1<-sum(as.numeric(freq[,(i-2)]==unicos_l2[1]),na.rm = T)
          al1f<-l1_al1+l2_al1
          
          l2_al2<-sum(as.numeric(freq[,(i-2)]==unicos_l2[2]),na.rm = T)
          al2f<-l2_al2
          
          freqs<-c(al1f,al2f)
        }}
      
      if (length(unicos_l1)==1 & length(unicos_l2)==1){
        l1_al1<-sum(as.numeric(freq[,(i-1)]==unicos_l1[1]),na.rm = T)
        l2_al1<-sum(as.numeric(freq[,(i-2)]==unicos_l2[1]),na.rm = T)
        al1f<-l1_al1+l2_al1
        
        freqs<-al1f
      }
      
      unicos<-sort(unique(c(unicos_l1,unicos_l2)))
      for (j in 1:length(rownames(data))){
        a=data[j, (i-1)]
        b=data[j, (i-2)]
        
        if (is.na(a)){next}
        if (is.na(b)){next}
        if (a==b){
          if (a %in% unicos==F){
            data[j,i]<-0
          }
          
          if (a %in% unicos){
            if (length(unicos)==1){
              data[j,i]<-2
            }
            if (length(unicos)==2){if (a==unicos[which.min(freqs)]){
              data[j,i]<-0
            }
              if (a==unicos[which.max(freqs)]){
                data[j,i]<-2
              }  }
          }}
        
        if (a!=b){
          data[j,i]<-1
        }
      }
      
    }
  }
  data<-data[,endsWith(colnames(data),"2")]
  return(data)
}

testgen<-genotypes(data = matriz, freq = ind_freqed)


remove_invariable<-function(genotyp){
  t<-data.frame(row.names = row.names(genotyp))
  for (i in 1:length(colnames(genotyp))){
    if (length(na.omit(unique(genotyp[,colnames(genotyp)[i]])))==2){
      col<-colnames(genotyp)[i]
      t[col]<-genotyp[,i]
    }
    if (length(na.omit(unique(genotyp[,colnames(genotyp)[i]])))==3){
      col<-colnames(genotyp)[i]
      t[col]<-genotyp[,i]
    }
  }
  return(t)
}

final_matrix<-remove_invariable(testgen)

filename<-basename(vcfp)
outname<-sub(".recode.vcf","",basename(vcfp))
outname<-paste(outname, ".genotype.txt", sep = "")
write.table(final_matrix, outname)
