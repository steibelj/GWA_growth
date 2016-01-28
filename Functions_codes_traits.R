#Sys.setlocale(locale="C")
################################################################################################################################################
# Functions paper:                                                                                                                             #
# "Refining genome wide association for growth and fat deposition traits in a F2 pig population."                                              #
################################################################################################################################################


##############################
## FUNCTIONS FOR PAPER CODE ##
##############################

#1.######## function zstandard ########
#function to compute the matrix Z for create G matrix, this is a generic function is included, 
#which computes allele frequencies from certain genotype matrix and also, its expected frq and respective
#standard deviation. This function allows two options: Homogeneous, which create Z matrix without divide normalized
#counts of B allele by the expected stand deviation Heterogeneous, which create Z matrix dividin normalized counts
#of B allele by the expected stand deviation
#First step is the definition of input genotypes file (genomatrix) and also the procedure to
#compute Z (homogeneous or heterogeneous). In the case of parameter "alfreq" it can be NULL and then, it will be computed by the function. 
#Alternatively, a vector of allele freq can be provided by the user.In our case, we will allow function to compute allele freq


#input: genomatrix : matrix with genotypes (animals in rows and markers in columns)
#       alfreq: frequency of the markers optional, NULLif the frequency for the markers is not given  
#       procedure: specify the standardization method=> "homogeneous"= weights markers by the sum of the expected variance across loci
#                                                       "heteregeneous"= weights markers by reciprocals of their expected variance  
#output:
#        Z: Return a matrix "Z" standardazid by "homogeneous" or "heteregeneous" procedure
#        


zstandard<-function(genomatrix,alfreq=NULL,procedure="homogeneous"){
  # However and before computations, it is convenient to do some controls of input genotype matrix.
  
  # 1. Check that a genotype matrix is provided
  if (!is.matrix(genomatrix))stop("a genotype matrix is needed")
  
  # 2. Check that genotype matrix contains numeric values (for instance, genotypes expressed as 0,1,2)
  if (!is.numeric(genomatrix)) stop ("genomatrix contains non numeric values")
  
  # 3. Check that allele dosages are in the range [0,2]
  if ((max(genomatrix,na.rm=T)>2)|(min(genomatrix,na.rm=T)<0)) stop ("valid allelic dosages should
                                                                     be in the interval [0,2]")  
  # 4. Check that genotype matrix doesn't contain NA
  if (sum(is.na(genomatrix))) warning("genomatrix contains NA values")
  # After quality control checks, and if allele freq are not provided, compute them obtaining
  #means per column divided by 2
  if (is.null(alfreq)) alfreq<-colMeans(genomatrix,na.rm=T)/2
  
  # Define names of each element of vector containing allele freq (same as SNP names)
  names(alfreq)<-colnames(genomatrix)
  
  # Further control steps;  
  # Check for inconsistencies between dimensions
  if (ncol(genomatrix)!=length(alfreq)) stop("length of vector alfreq should equal number of
                                             columns in genomatrix")  
  # Check for inconsistencies between SNP names
  if (sum(names(alfreq)!=colnames(genomatrix))!=0) stop ("SNP names on frequency vector and
                                                         genomatrix don't match")  
  # Check for presence of NAs
  if (is.na(sum(alfreq))) stop ("allelic frequency vector can't have NA values")
  
  # Check for allele freq greater than 1 or less than 0
  if ((min(alfreq)<=0)|(max(alfreq)>=1)) stop ("allelic frequency should be in (0,1)")
  
  # After previous steps, compute expected value of allele dosages
  eg<-2*alfreq  
  # Compute expected standard deviation of allele dosages,required in heterogeneous option
  sdg<-eg*(1-alfreq)
  if (procedure=="homogeneous")  {
    
    # In this case, Z will be the matrix containing normalized allelic dosages (counts of allele
    #?B? minus its expected value)
    
    #     divided by the sum of the standard deviation
    Zn<-as.matrix(sweep(genomatrix,2,FUN="-",STATS=eg))
    Zn<-Zn/sqrt(sum(sdg))
  }
  else {
    
    # Alternatively, and to make G matrix analogous to the numerator relationship matrix A, Z matrix will be the matrix containing
    
    #     but divided by expected standard deviation of allele dosages
    
    Zn<-as.matrix(sweep(genomatrix,2,FUN="-",STATS=eg))
    Zn<-as.matrix(sweep(Zn,2,FUN="/",STATS=sqrt(sdg)))
    Zn<-Zn/sqrt(ncol(Zn))
  }
  return(Zn)
}



#2.######## function gblup ########
#function to compute SNP effects and BLUP
#input: y : is a vector with phenotypes (can have missing values as NAs)
#       mf: is a incidence matrix with fixed effects and covariates
#       zf : is a Matrix with Standarized genotypes for F2 animals
#       g : is a matrix with Genomic relationship for F2 animlas
#       Ginv: is a invert matrix of "g"

#output:
#        llik:  LogLikehood of the model "f"
#        a_hat: Genome breeding values (GEBVs).
#        E_variance: Error variance  
#        A_variance: Additive variance  
#        Heritability: Heritability of the trait 
#        iteratons: Number of iterations to converge. 
#        Ginv: The inverse of G matrix
#        ll0: LogLikehood of the model "f0"


gblup<-function(y,mf,zf,g=NULL,Ginv=NULL){
  if(is.null(g)) g=zf%*%t(zf) #if g is NULL, the matrix g is created from zf
  #Run "regress"
  y[is.na(y)] <- 0
  f<-regress(y~mf,~g,pos=c(T,T))
  f0<-regress(y~mf,pos=c(T))
  #Components variance estimate by REML(Newton-Raphson algorithm to maximise the residual log likelihood)
  ve<-f$sigma["In"] # VarE
  vu<-f$sigma["g"]  # VarU
  h2<-vu/(ve+vu)    #Heredability
  itercon<-f$cycle  #Number of cycles to convergence.
  
  #BLUP and animal effect(ah)
  blupf<-BLUP(f)
  ah<-blupf$Mean
  
  #Inverse of GF2
  if (is.null(Ginv)) Ginv<-solve(g)
  
  return(list(llik=f$llik,a_hat=ah,E_variance=ve,A_variance=vu,Heritability=h2,n_iter=itercon,Ginv=Ginv,ll0=f0$llik)) 
}



#3.######## function abmap ########
#function to compute absolute map
#input: 
#       map: Data frame or matrix
#            first column  -> chromosome name
#            second column -> NUMERIC, position relative to chromosome
#output: 
#       Gives a vector with the Absolute position  

abmap<-function(map){
  dims<-dim(map)[2]
  if (dims!=2){
    stop ("Dimensions incorrect")
  }
  
  if(!is.numeric(map[,2])){
    stop ("column 2 must be numeric")
  }
  
  chrom<-names(table(map[,1]))       
  last_pos<-0                        
  absolute_map<-rep(0,nrow(map))     
  for (j in chrom){                         # j comes from 1 to # total of chromosome
    idx<-map[,1]==j                         # index for id column chromosome
    absolute_map[idx]<-map[idx,2]+last_pos
    last_pos<-absolute_map[idx][sum(idx)]
  }
  return(absolute_map)
}


#4.######## function varsnp ########
#Function to compute elements for Variance SNP effects
#input: 
#       zf : is a Matrix with Standarized genotypes for F2 animals
#       iG: is a invert matrix of "g"
#       iGcG: is the result of iG%*%Cuu%*%iG, Cuu: Extract from MME

#output: Gives two rows (elements to compute the variance of the each marker effect) 
#        First row is the result  of:  zf'%*%iG%*%zf
#        Second row is the result of:  zf'%*%iG%*%Cuu%*%zf%*%iG%*zf

varsnp<-function(zf,iG,iGCG){
  t1<-zf%*%iG%*%zf
  # from Gualdr?n Duarte et al.(2014):Z'(G^-1)Z(VarA^2)-Z'(G^-1)Cuu(G^-1)Z
  t2<-zf%*%iGCG%*%zf 
  return(c(t1,t2))
}


#5.######## function snpe_GWA ######## 
#Function to compute the SNPe or "g_hat", SNP variance and p-values for each SNP
#input: 
#       trout_gblup : the output of function "gblup"
#       x: is a matrix that indicates the sex of the animal in the F2
#       Znf2: is the matrix Z extract for the F2 animals
#output: 
#        beta: column with the SNP effect "SNPe"
#        snp_variance: SNP variance for each SNP  
#        pvalues:p-values for each SNP

snpe_GWA<-function(trout_gblup,x,Znf2){

#########################
## Create Cuu from MME ##
#########################
### Matrix Z
z<-diag(nrow(Znf2))

### Cuu= ve*[z'z - z'x(x'x)^-1 x'z + G^-1*lambda]^-1
#z'z
zzt<-t(z)%*%z
#z'x(x'x)^-1 x'z
z1<-solve(t(x)%*%x)
z2<-t(z)%*%x%*%z1%*%t(x)%*%z
##
ve1<-trout_gblup[[3]]                # error variance "E_variance"
vu1<-trout_gblup[[4]]                # aditive variance "A_variance"
lambda<-as.numeric(ve1/vu1)
G_inv<-trout_gblup[[7]]              # Inverse G matrix "Ginv"
##
Cuu<-(solve(zzt-z2+(G_inv*lambda)))*as.numeric(ve1)

#########################
## SNP effect variance ##
#########################
iGCG<-G_inv%*%Cuu%*%G_inv

vars<-apply(Znf2,2,varsnp,iG=G_inv,iGCG=iGCG)
vsnp<-vars[1,]*vu1-t(vars[2,])                  #Extract the variance effect for each SNP

sdsnp<-sqrt(vsnp)                               # Standar desviation snp effects

uh<-t(Znf2)%*%G_inv%*%trout_gblup[[2]]          # SNPe "g_hat"
SNPe_ad<-uh/t(sdsnp)                            # SNP_effect"SNPe"/Standar_deviaton of snp effects
pvalue<-2*(1-pnorm(abs(SNPe_ad)))               # p-value for each snp 

return(list(beta=uh,snp_variance=vsnp,pvalues=pvalue))

}



#6.######## function delete.NULLS  ########
#Function to delete null/empty entries in a list
#input:
#      list file with null or empty entries
#output:
#      list file without null or empty entries

delete.NULLs  <-  function(x.list){    
  x.list[unlist(lapply(x.list, length) != 0)] 
} 



#7.######## function SNPS_FDR ####### 
#Function to compute the False Discovery Rate (FDR) > "fdr" or thershold
#input: 
#       gwatrait: the output of function "snpe_GWA" 
#       map: marker map with: SNP name as rowname, 1st column: chromosome, 2nd Column: Position Mega-bases 
#       fdr: false discovery rate thershold 
#output: 
#       list:marker map with: SNP name as rowname,
#        1st column:trait
#        2nd column:chromosome
#        3th column:Position of the SNP in Mega -bases
#        4th column:p-value
#        5th column: -log10(p-value)
#        6th column: q-value
#        7th column: the SNP effect

SNPs_FDR<-function(gwatrait,map,fdr=0.05){
 
  ### Read the outpufiles (gwa_"trait") for all growth traits ##
  t_pv<-list.files(pattern="gwa") # put in a list files with pattern gwa
  nt<-length(t_pv)  # [1] number of traits
  namet<-sub("gwa_","",t_pv)
  
  ### LOOP FOR ALL GROWTH TRAITS ###
  
  result_qv <- list() 
  
  for(i in 1:nt){
    
    load(t_pv[i]) 
    qvs<-qvalue(gwa_trait$pvalues) #Apply function qvalue for the p-values in gwa_trait
    qvsidx<-qvs$qvalues<fdr        #Index Filter by False Discovery (FDR) thershold 
    qvst<-qvs$qvalues[qvsidx,]     #Extract the q-values filtered by FDR
    
    #name trait:
    #list the SNP name and the q-value, when multiple traits are in the loop and one SNP  
    #filtered is significant for more than one trait.
    
    nametrait<-rep(namet[i],length(qvst))   
    namet_qvs<-cbind(nametrait,qvst)
    
    #pvalue:
    #Extract the -Log10(p-value) for the SNP filtered by FDR and trait
    pvalue<-gwa_trait$pvalues
    pvalueidx<-rownames(pvalue)%in%rownames(namet_qvs)
    pvalue_names<-pvalue[pvalueidx,]
    logpvalues<--log(pvalue_names,10)
    
    #Beta:
    #Extract the SNP effect for the SNP filtered by FDR and trait
    beta<-gwa_trait$beta
    betaidx<-rownames(beta)%in%rownames(namet_qvs)
    beta_names<-beta[betaidx,]
    
    #map:
    #Extract the SNP map position for the SNP filtered by FDR and trait    
    mapidx<-rownames(map2)%in%rownames(namet_qvs)
    map_names<-map[mapidx,]
    
    # Here joint the information of pvalue, -log10(p-value), Beta(or SNP effect), map position of the SNP selected by FDR and trait
    if (nrow(namet_qvs)==1){
      result<-cbind(rownames(namet_qvs),namet_qvs[,1],map_names[1],map_names[2],pvalue_names,logpvalues,namet_qvs[,2],beta_names)
    }else{  
      
      result<-cbind(rownames(namet_qvs),namet_qvs[,1],map_names[,1],map_names[,2],pvalue_names,logpvalues,namet_qvs[,2],beta_names)
    }
    
    result_qv[[i]]<-result
    
  }
  
  result1<-delete.NULLs(result_qv)        # Apply function and delete blank spaces 
  result2 <- do.call(rbind, result1)      # Joint list (result1)
  rownames(result2)<-result2[,1]          # rownames as SNPname
  snpFDR<-result2[,-1] 
  
  # snpFDR: List per trait the SNP filtered by qvalues<0.05 (FDR<0.05)
  colnames(snpFDR)<-c("trait","Chr","Pos","pvalue","logpvalue","qvalue","beta")
  return(snpFDR)
  
}


#8. ######## FUNCTION pvalue_FDR ########
# Function to select snp by FDR and lowest pvalue by chromosome and trait
#input:
#      output file from function "SNPs_FDR", where the columns are:
#      "trait"-"Chr"-Chromosome-"Pos": Position (Mega-bases) -"pvalue"-"logpvalue":-log10(pvalue) -"qvalue"-"beta": SNP effect    
#output:
#      file with the lowest palue per chromosome/trait, filtering by FDR.
#      Where tthe columns are:
#      "snpID": SNP name - "trait"- "Chr" - "Pos" - "pvalue" - "logpvalue" - "qvalue" - "beta"


pvalue_FDR<-function(snpFDR){
  
  snpID<-rownames(snpFDR)
  snpFDR2<-cbind(snpID,snpFDR) # Add the rownames(SNP names) in the forst column
  
  ##### 
  tr<-unique(snpFDR2[,2])  # Vector with the names trait
  #List for ranges
  range2<-list()
  range3<-list()
  range_trait<-list()
  
  #List for min pvalues
  snp_b<-list()
  snp_c<-list()
  for(i in 1:length(tr)){
 
    snpt<-subset(snpFDR2,snpFDR2[,2]==tr[i])
    chr<-unique(snpt[,3])   # Chromosomes with significant p-value per trait
    
    for(j in 1:length(chr)){ # Select the minimum p-value in each chromosome per trait

      # Select the minimum p-value in each chromosome per trait  
      chr_t<-subset(snpt,snpt[,3]==chr[j])
          
      
      minpvalue<- which.min(chr_t[,5])
      snp_a<-chr_t[minpvalue,]
      
      snp_b[[j]]<-snp_a
      t_minpv<-do.call(rbind,snp_b)
      
      # Positions distance-range and SNP number per trait and chromosome   
      minpos<- min(chr_t[,4])  
      maxpos<- max(chr_t[,4])
      nsnp<-nrow(chr_t)
      range1<-cbind(unique(chr_t[,2]),nsnp,chr[j],minpos,maxpos)  
      range2[[j]]<-range1
      range3<-do.call(rbind,range2)
      
    }
    range_trait[[i]]<-rbind(range3)
    
    #List by trait and chromosome the SNP with lowest p-value and filtered by 0.05<FDR
    snp_c[[i]]<-t_minpv
    
  }
  
  # SNP with minimum p-value in each trait=>
  snp_mpvalueFDR<-do.call(rbind,snp_c)
  snp_mpvalueFDR<-unique(snp_mpvalueFDR[,c(1,2,3,4,5,6,7,8),drop=FALSE])
  
  #List info. for the SNP selected: 
  #snpID, trait, Chromosome(chr), Position (Pos), pvalue,logpvalue, qvalue, beta. 
  snp_mpvalueFDR<-snp_mpvalueFDR[sort.list(snp_mpvalueFDR[,3]), ] # order the output by chromosome
  
  return(snp_mpvalueFDR)
}



#9.######## function summary.ll.fit ########
#Function to order the outputfiles of regress package for models: "y=Xb+a+e" and "y=Xb+a1+a+e" 
#input: 
#      list of output files from regress for model1:y=Xb+a+e and model2:y=Xb+a1+a+e           
#output:
#       Loglikem2:   LogLikehood "model2"  
#       LogLikem1:   LogLikehood "model1"   
#       LRT:         Likehood Ratio Test for "model1" and "model2"
#       LRTseg:      p-value for Likehood Ratio Test for the segment
#       varE1:       Error variance   of "model1"
#       varA1:       Additive variance   of "model1"
#       varE2:       Error variance   of "model2" 
#       VarA2:       Additive variance  of "model2"
#       Varseg:      Additive variance segment of "model2" 
#       Varseg_pr:   Proportion in % of the total variance explained by the segment.

summary.ll.fit<-function(x){
  m1<-x$l1
  m0<-x$l0
  ll1<- as.numeric(m1["llik"])  # llik:Value of the marginal log likelihood at the point of convergence.
  ll0<- as.numeric(m0["llik"])
  print(c(ll1,ll0))
  LRT<- 2*(ll1 -ll0)
  varE1<-m1$sigma["In"]         #Error Variance 
  varE0<-m0$sigma["In"] 
  nms0<- m0$sigma[!(m0$Vnames%in%"In")] # Vnames: Names associated with each variance component, used in print.regress.
  nms1<- m1$sigma[!(m1$Vnames%in%"In")]
  varU0<-nms0[length(nms0)]     #Aditive Variance
  varU1<-nms1[length(nms1)]
  varBR<-nms1[length(nms1)-1]   # Variance of the segment
  q1<-varBR/(varBR+varU1+varE1) # % of variablity explain by QTL or segment
  LRTseg<-1-pchisq(LRT,0.5)     # Likehood ratio test for the segment 
  return(list(Loglikem2=ll0,LogLikem1=ll1,LRT=LRT,LRTseg=LRTseg,varE1=varE0,varA1=varU0,varE2=varE1,VarA2=varU1,Varseg=varBR,Varseg_pr=q1,c<-"#-----------#"))
}



#10.######## function propor_seg ########
# Function to assess the additive variance proportion segment 
#input:
#      trait  : vector with the phenotype (rownames=animal ID)
#      map    : first column: chromosome number, second column: marker Position, (rownames=markerID)
#      genof2 : Genotype matrix, animal(rownames: animal ID) x marker(columnname: marker ID)  
#      frqF0  : Frequency markers for F0 animals   
#      G      : Genomic Relationship matrix (dimention: #rows_genof2 x #rows_genof2)
#      Zf2    : Standardaized matrix of gnoetypes (Zf2 x t(Zf2) = G) 
#      dis_snp: length of the segment, distance in Mb up and down stream from the SNP selected
#      snp_character : matrix(snp_name_selected, trait_name, chromosome_number, Position(Mega-baese)   
#output/per trait:
#       snpname:     SNP name selected (form segment)  
#       trait:       trait name
#       Loglikem2:   LogLikehood "model2"  
#       LogLikem1:   LogLikehood "model1"   
#       LRT:         Likehood Ratio Test for "model1" and "model2"
#       LRTseg:      p-value for Likehood Ratio Test for the segment
#       varE1:       Error variance   of "model1"
#       varA1:       Additive variance   of "model1"
#       varE2:       Error variance   of "model2" 
#       VarA2:       Additive variance  of "model2"
#       Varseg:      Additive variance segment of "model2" 
#       Varseg_pr:   Proportion in % of the total variance explained by the segment.

propor_seg<-function(trait,map,genof2,frqF0,G,Zf2,dis_snp=NULL,snp_character,chr=NULL,nameplow=NULL,namephigh=NULL){  
  
  ######################
  # G1 for the Segment #
  ######################
  
  # Made a particular segment  #
  if(is.null(nameplow)){ 
  
  # distance (dis_snp)  Mb up-down stream from the Reference SNPs   
  ref_snp<-grep(snp_character[1,1],rownames(map))
  chr_snp<-as.numeric(map[ref_snp,1])
  pos_snp<-map[ref_snp,2]
  phigh<-as.numeric(pos_snp)+dis_snp
  plow<-as.numeric(pos_snp)-dis_snp
  #Extract chromosme map from the SNP reference
  chr_map<-subset(map,map[,1]==chr_snp)
  #Extract SNP in this range
  slxmap<-subset(chr_map,(chr_map[,2]<=phigh)&(chr_map[,2]>=plow))
  slxnames<-rownames(slxmap) 
  
  }else{
 
  # Use nameplow and namephigh to demarcate the genome segment   
  idxpl<-rownames(map)==nameplow
  plow<-as.numeric(map[idxpl,2])-dis_snp
  
  idxph<-rownames(map)==namephigh
  phigh<-as.numeric(map[idxph,2])+dis_snp
  
  ### Extract SNP in this range
  mapslx<-subset(map,map[,1]==chr ) # Select the chromose in the map
  slxmap<-subset(mapslx,(mapslx[,2]<=phigh)&(mapslx[,2]>=plow))
  slxnames<-rownames(slxmap) 

  }
  
  
  ### compute matrix Z1 and G1 for snp selected in the segment (normalized separtely) 
  Zidx<-colnames(genof2)%in%slxnames
  genof2fil<-genof2[,Zidx]
  dim(genof2fil) 
  
  idxfrq<-names(frqF0)%in%slxnames # index to Select the frequencies for the SNP segment
  sum(idxfrq==TRUE) 
  fil_frq_f0<-frqF0[idxfrq]        # Extract the F0 frequencies for the SNP segment
  length(fil_frq_f0) 
  
  # Matrices Z1 and G1
  Z1<-zstandard(genof2fil,alfreq=fil_frq_f0,procedure="heterogeneous") #Apply the function for the SNP segment
  G1<-Z1%*%t(Z1)                  # G1 for SNP selected
  
  ### G with out SNP columns belong to the segement
  Zidx<-colnames(Zf2)%in%slxnames
  Zseg<-Zf2[,Zidx]     
  G2<-G-Zseg%*%t(Zseg)
  
  #phenotype  
  names_trait<-snp_character[,2]  
  res<- list() 
  for(i in 1:length(names_trait)){
    pheno<-trait[,as.character(names_trait[i]),1]
    indx<-match(rownames(Zf2),names(pheno))
    pheno<-pheno[indx]
    
    ##############
    ##  Models  ##
    ##############
    
    model1<-regress(pheno~x,~G,pos=c(T,T))     # y=Xb+a+e
    
    model2<-regress(pheno~x,~G1+G2,pos=c(T,T)) # y=Xb+a1+a+e: add segment into the model
    
    ##################################
    ##  Performed models comparison ##
    ##################################
    
    llik<-list(Result=list(l0=model1,l1=model2)) 
    
    #################################################
    ##  Results from Regress for model1 and model2 ##
    #################################################
    model_result<-sapply(llik,summary.ll.fit) #
    
    res[[i]]<-rbind(snpname=snp_character[,1],trait=as.character(names_trait[i]),model_result)
    
  }
  
  sig_seg<-do.call(rbind, res) #merge the results per trait
  
  return(sig_seg)
}



#11.######## function plot_ld ########
# Function to create input files for "snp.plotter" R package, and obtain the Linkage Disequilibrium-plot 
# 
#input:
#      smap: map file for the SNP selected. rownames:SNP ID, 1st column: Chromosome, 2nd column: Position Mega-bases
#      pvalue: p-values for the SNp selected
#      ped: pedigree file with columns:"Family","ID":animla ID,"Sire","Dam","sex": 1(male)-2(female),"pheno": -9  
#      genotype: row-names: animal ID, column-names: Marker-ID, genotype codification as 0, 1 and 2 (counting of the reference allele) 
#      name_trait: name of the trait
#output:
#      input files for "snp.plotter"   
#      snpfile:  marker file
#      genofile: genotype file
#      config:   configuration file to run "snp.plotter"    


plot_ld<-function(smap,pvalue,ped,genotype,name_trait){
  
  ###############  
  #1) SNP FILE ## 
  ###############
  
  #Select the p-values for markers in the segment
  idxp<-rownames(pvalue)%in%rownames(smap)
  pvt1<-pvalue[idxp,] 
  sigplus<-rep("+",length(pvt1)) # Add column with "+" symbol
  snpIDs<-rownames(smap)
  snpfile<-as.data.frame(cbind(sigplus,snpIDs,smap[,2]*1000000,pvt1)) 
  colnames(snpfile)<-c("ASSOC", "SNP.NAME",  "LOC"  , "SS.PVAL") #Add column names 
  
  #2) GENOTYPE FILE ## 
  
  #Extract segment genotype 
  idxplot<-colnames(genotype)%in%rownames(smap)
  geno_seg<-genotype[,idxplot]     
  
  #Round imputed genotypes between [0-2]
  geno_seg<-round(geno_seg)
  
  #Replace values
  geno1<-replace(geno_seg,geno_seg==2,"11")
  geno2<-replace(geno1,geno1==1,"12")
  geno3<-as.data.frame(replace(geno2,geno2==0,"22"))
  
  result <- list() 
  for(i in 1:ncol(geno3)){
    geno4 <- data.frame(do.call('rbind', strsplit(as.character(geno3[,i]),"",fixed=TRUE)))
    result[[i]]<-cbind(geno4)  
    
  }
  
  result2 <- do.call(cbind, result) 
  rownames(result2)<-rownames(geno3)
  
  ## Extract the family-id-idsire-iddam-sex-infec for result2
  id_result2<-subset(ped,ped[,2]%in%rownames(result2))
  #check order
  idr2<-as.numeric(rownames(result2))
  if (sum(as.numeric(id_result2[,2])-idr2)!=0)stop("order is not the same in ped(1st 6 columns) and genotype recouded")# [1] 0 !!!order OK  
  #replace affection status ("-9") column 6 "genofile1" by "1"  
  column6<-replace(id_result2,id_result2==-9,1)
  #GENOFILE=>
  genofile<-cbind(column6,result2)
  
  
  #3) Script to create "config.txt" for  "snp.plotter" ## 
  
  nline1 = paste ("SNP.FILE=snpfile_",name_trait,".txt",sep="") 
  nline2 = paste ("HAP.FILE=NULL")
  nline3 = paste ("PALETTE.FILE=NULL")
  nline4 = paste ("EVEN.SPACED=FALSE")
  nline5 = paste ("USE.GBL.PVAL=FALSE")                            
  nline6 = paste ("DISP.HAP=FALSE" )
  nline7 = paste ("DISP.LDMAP=TRUE")
  nline8 = paste ("LD.COLOR.SCHEME=heat")
  nline9 = paste ("COLOR.LIST=red")
  nline10 = paste ("SYMBOLS=circle-filled")
  nline11 = paste ("PVAL.THRESHOLD=1")
  nline12 = paste ("LAB.Y=log")
  nline13 = paste ("GENOTYPE.FILE=genofile_",name_trait,".dat",sep="") 
  nline14 = paste ("LD.TYPE=rsquare")
  nline15 = paste ("DISP.COLOR.BAR=TRUE")
  nline16 = paste ("DISP.TYPE=symbol")
  nline17 = paste ("DISP.LEGEND=FALSE")
  nline18 = paste ("SAMPLE.LABELS=LD_segment")
  nline19 = paste ("IMAGE.TYPE=pdf")
  nline20 = paste ("DISP.SNP=TRUE")
  nline21 = paste ("DISP.MULT.LAB.X=FALSE")
  nline22 = paste ("IMAGE.TITLE=LD SEGMENT ") 
  nline23 = paste ("IMAGE.NAME=ld_plot_",name_trait,sep="")   
  nline24 = paste ("IMAGE.SIZE=7")
  nline25 = paste ("CONNECTING.LINES.FACTOR=1.5")
  nline26 = paste ("DISP.PHYS.DIST=TRUE")
  nline27 = paste ("CONNECTING.LINES.ADJ=0.02")
  
  config = rbind(nline1,nline2,nline3,nline4,nline5,nline6,nline7,nline8,nline9,nline10,nline11,nline12,nline13,nline14,nline15,nline16,nline17,nline18,nline19,nline20,nline21,nline22,nline23,nline24,nline25,nline26,nline27)
  
  return(list(snpfile=snpfile,genofile=genofile,config=config))
  
}


