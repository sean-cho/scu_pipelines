# This script contains all the modular functions and a pipeline function to
# read, filter, and analyze methylation data per SCU pipeline.

# A custom S4 class is created to handle the object through the pipeline

require(minfi)
require(limma)
require(org.Hs.eg.db)
require(scales)
require(DMRcate)
require(toolkit)
source('code/def_classes.R')
load('annot/il450k.rda')

.makedesigns <- function(scu_meth, phenotype){
  foo <- factor(pData(scu_meth@processed_dat)[,phenotype,drop=TRUE])
  bar <- model.matrix(~ 0 + foo)
  colnames(bar) <- gsub('foo', '', colnames(bar))
  return(bar)
}

.process_beta_table <- function(x, 
                                annot = c('gene_symbol', 'Relation_to_Island',
                                          'UCSC_RefGene_Group'),
                                filter = FALSE){
  if(filter){
    foo <- subset(x, adj.P.Val <= 0.05)
  } else {
    foo <- x
  }
  bar <- data.frame(il450k[,annot], stringsAsFactors = FALSE)
  bar$Relation_to_Island = gsub('S_|N_', '', bar$Relation_to_Island)
  foobar <- merge(foo, bar, by = 0)
  rownames(foobar) <- foobar$Row.names
  foobar$Row.names <- NULL
  foobar <- foobar[rownames(foo),]
  return(foobar)
}

.limma_mule <- function(m, des, comparison){
  res_lm <- lmFit(m, design = des)
  res_contrasts <- contrasts.fit(fit = res_lm, 
                                 contrasts = 
                                   makeContrasts(contrasts = comparison,
                                                 levels = des))
  res_ebayes <- eBayes(res_contrasts)
  res_toptable <- topTable(res_ebayes, number = nrow(m))
  # rownames(res_toptable) <- res_toptable$ID
  res_toptable <- .process_beta_table(res_toptable)
  
  output <- new(Class = 'SCU_Methylation_Limma',
                model.mat = des,
                model = comparison,
                linear_fit = res_lm,
                contrasts_fit = res_contrasts,
                ebayes = res_ebayes,
                toptable = res_toptable)
  return(output)
}

pmm_01_readqc <- function(target_dir,
                          output_dir = 'output',
                          pdat_file = NA,
                          pval_file = NA,
                          th_pval = NULL,
                          th_proberate = 0.95,
                          th_samplerate = 0.8,
                          verbose = 2){
  
  ### Set up
  if(verbose >= 1)
    cat('Checking for data consistency.\n')
  
  # Create output directory if it does not exist
  if(!dir.exists(output_dir)){
    dir.create(output_dir)
  }
  
  # Evaluate if files exist
  foo <- c(pdat_file, pval_file)
  names(foo) <- c('pdat_file', 'pval_file')
  bar <- sapply(foo, function(x){
    if(is.na(x))
      return(TRUE) else
        return(file.exists(x))
  })
  if(!all(as.logical(bar))) {
    missing_files <- paste('The following files are not found:',
                           paste(names(which(!bar)), collapse = '; '), sep = '\n')
    stop(missing_files)
  }
  
  # Define thresholds
  if(missing(pval_file) & missing(th_pval)){
    th_pval <- 1e-8
  } else {
    if(missing(th_pval)) {
      th_pval <- 0.01
    }
  }
  
  ### Read and QC
  
  # Read sample sheet and idat files
  if(verbose >= 1)
    cat('Reading samplesheet.\n')
  ss <- read.metharray.sheet(base = target_dir)
  if(verbose >= 1)
    cat('Reading data.\n')
  rgdat <- read.metharray.exp(targets = ss, verbose = verbose >= 2)
  
  # Get phenodata
  if(missing(pdat_file)){
    pdat <- pData(rgdat)
  } else {
    pdat_rg <- pData(rgdat)
    pdat_ext <- read.delim(pdat_file, stringsAsFactors = FALSE)
    
    # Scan for samples
    foo <- apply(pdat_ext, 2, function(x) mean(pdat_rg$Sample_Name %in% x))
    if(max(foo, na.rm = TRUE) == 0){
      stop('No matching samples identified.')
    } else {
      bar <- names(which.max(foo))
      if(max(foo, na.rm = TRUE) != 1){
        matching_n <- sum(which(pdat_rg$Sample_Name %in% pdat_ext[,bar]))
        warning_pdat <- paste('Only', matching_n, 'of', length(pdat_rg$Sample_Name), 
                              'samples were matched between external and samplesheet samples.\n')
        warning(warning_pdat)
        rm(matching_n)
      }
    }
    pdat <- merge(pdat_rg, pdat_ext, by.x = 'Sample_Name', by.y = bar, all.x = TRUE)
    rm(pdat_rg, pdat_ext)
  }
  
  pData(rgdat) <- pdat
  rownames(pdat) <- pdat$Sample_Name
  sampleNames(rgdat) <- pdat$Sample_Name
  
  ### Run QC
  
  if(verbose >= 1)
    cat('Calculating QC metrics.\n')
  
  if(verbose >= 2)
    cat('Generating QC report.\n')
  # QC report for control probes
  qcReport(rgSet = rgdat, sampNames = pData(rgdat)$Sample_Name,
           pdf = paste0(output_dir, '/qcreport.pdf'))
  
  if(verbose >= 2)
    cat('Estimating beta-value densities.\n')
  # Beta value distribution
  foo <- getBeta(rgdat)
  bar <- apply(foo, 2, density, na.rm = TRUE)
  bar_x <- range(unlist(lapply(bar, function(x) x$x)))
  bar_y <- range(unlist(lapply(bar, function(x) x$y)))
  bar_col <- colorRampPalette(c('orange', 'skyblue3'))(nrow(pData(rgdat)))
  pdf(file = paste0(output_dir, '/beta_distribution_raw.pdf'),
      width = 8, height = 6)
  plot(0, type = 'n', xlim = bar_x, ylim = bar_y, xlab = 'Beta values',
       ylab = 'Density', main = 'Beta value distributions')
  i <- 0
  lapply(bar, function(x){
    i <<- i + 1
    lines(x, lwd = 2, col = alpha(bar_col[i], 0.5))
  })
  # i_group # plan to implement separated legends
  plot.new()
  legend('center', legend = sampleNames(rgdat), 
         lty = 1, col = bar_col, lwd = 2)
  dev.off()
  
  # Detection P
  if(verbose >= 2)
    cat('Estimating detection p-value and call rates for the following thresholds:',
        '\np-value threshold              : ', th_pval,
        '\nprobe-wise callrate threshold  : ', th_proberate,
        '\nsample-wise callrate threshold : ', th_samplerate,
        '\n', sep='')
  detp <- detectionP(rgdat)
  
  # Probe-wise
  probe_callrate <- apply(detp <= th_pval, 1, mean, na.rm = TRUE)
  probe_callrate_bin <- probe_callrate >= th_proberate
  cat(length(which(probe_callrate_bin)), '/485577 probes passed QC.\n', sep = '')
  
  # Sample-wise
  sample_callrate <- apply(detp[probe_callrate_bin,] <= th_pval, 2, 
                           mean, na.rm = TRUE)
  sample_callrate_bin <- sample_callrate >= th_samplerate
  cat(length(which(sample_callrate_bin)), '/', length(sample_callrate_bin),
      ' samples passed QC.\n', sep = '')
  
  output <- new(Class = 'SCU_Methylation', 
                rgdat = rgdat, 
                detp = detp,
                probe_callrate = probe_callrate,
                sample_callrate = sample_callrate,
                pval_threshold = th_pval,
                callrate_thresholds = c(sample_threshold = th_samplerate,
                                        probe_threshold = th_proberate))
  
  return(output)
}

pmm_02_normalization <- function(scu_meth,
                                 output_dir = 'output',
                                 th_pval = NULL,
                                 th_proberate = NULL,
                                 th_samplerate = NULL,
                                 filter_samples = TRUE,
                                 filter_probes = TRUE,
                                 preprocess = c('funnorm', 'raw',
                                                'illumina', 'quantile'),
                                 verbose = 2, 
                                 snp_filter = TRUE,
                                 snp_filter_dist = 2,
                                 snp_filter_maf = 0.05,
                                 getm = TRUE,
                                 ...){
  
  ### Set up
  if(verbose >= 1)
    cat('Checking for data consistency.\n')
  
  if(!inherits(scu_meth, 'SCU_Methylation')){
    stop('Input class is not SCU_Methylation. Please run psmm_01_readqc or wrapper.')
  }
  
  # Create output directory if it does not exist
  if(!dir.exists(output_dir)){
    dir.create(output_dir)
  }
  
  # Replace thresholds and recalculate if there are new user inputs
  if(!missing(th_pval)){
    cat('Recalculating call rates.\n')
    scu_meth@pval_threshold <- th_pval
    scu_meth@sample_callrate <- apply(scu_meth@detp <= scu_meth@pval_threshold,
                                      2, mean, na.rm = TRUE)
    scu_meth@probe_callrate <- apply(scu_meth@detp <= scu_meth@pval_threshold,
                                     1, mean, na.rm = TRUE)
  }
  if(!missing(th_samplerate)){
    scu_meth@callrate_thresholds['sample_threshold'] <- th_samplerate
  }
  if(!missing(th_proberate)){
    scu_meth@callrate_thresholds['probe_threshold'] <- th_proberate
  }
  
  # Apply filters
  if(filter_samples) {
    scu_meth@rgdat <- scu_meth@rgdat[,scu_meth@sample_callrate >= 
                                       scu_meth@callrate_thresholds['sample_threshold']]
  }
  
  ### Preprocessing
  
  if(preprocess[1] == 'funnorm'){
    if(length(sampleNames(scu_meth@rgdat)) < 3){
      preprocess[1] <- preprocess[2]
    } else {
      cat('Running functional normalization.\n')
      processed_dat <- preprocessFunnorm(scu_meth@rgdat, verbose = verbose >= 1, ...)
    }
  }
  if(preprocess[1] == 'raw'){
    cat('Running raw preprocessing.\n')
    processed_dat <- preprocessRaw(scu_meth@rgdat)
  }
  if(preprocess[1] == 'illumina'){
    cat('Running Illumina pre-processing.\n')
    processed_dat <- preprocessIllumina(scu_meth@rgdat, ...)
  }
  if(preprocess[1] == 'quantile'){
    cat('Running quantile normalization.\n')
    processed_dat <- preprocessQuantile(scu_meth@rgdat, ...)
  }
  scu_meth@processed_dat <- processed_dat
  
  ### Get beta and M-values
  b <- getBeta(scu_meth@processed_dat)
  if(filter_probes){
    b <- b[
      scu_meth@probe_callrate >= scu_meth@callrate_thresholds['probe_threshold'],]
    if(snp_filter)
      cat('Removing probes adjacent to SNPs.\n')
    b <- rmSNPandCH(b, dist = snp_filter_dist, mafcut = snp_filter_maf)
  }
  
  if(getm){
    m <- getM(scu_meth@processed_dat)[rownames(b),]
    m <- pmin(m, max(m[m != 'Inf']))
    m <- pmax(m, min(m[m != '-Inf']))
  }
  
  scu_meth@beta_values <- b
  scu_meth@m_values <- m
  
  return(scu_meth)
}

pmm_03_limma <- function(scu_meth,
                         phenotype,
                         comparison = NULL,
                         top_probes = 0.95,
                         limma_mode = c('m', 'b'),
                         output_dir = 'output',
                         verbose = 2){
  ### Set up
  if(verbose >= 1)
    cat('Checking for data consistency.\n')
  
  if(!inherits(scu_meth, 'SCU_Methylation')){
    stop('Input class is not SCU_Methylation. Please run psmm_01_readqc or wrapper.')
  }
  
  if(top_probes < 0){
    stop('Top_probes must be a value between 0 to 1, with 0 = no filtering.')
  }
  
  # Create output directory if it does not exist
  if(!dir.exists(output_dir)){
    dir.create(output_dir)
  }
  
  # Set up phenotype and comparisons
  if(verbose >= 1)
    cat('Setting up designs.\n')
  
  # Check for phenotype(s)
  if(missing(phenotype) | !all(phenotype %in% colnames(pData(scu_meth@processed_dat)))){
    stop('Column name indicating phenotype of interest is required.\n')
  }
  
  # Identify comparison(s)
  if(is.null(comparison)){
    # Populate list for comparison
    comparison <- apply(pData(scu_meth@processed_dat)[,phenotype,drop=FALSE],
                        2,
                        function(x) apply(combn(unique(x), 2), 2,
                                          paste0, collapse = ' - '))
  } else {
    # Check for list
    if(length(comparison)!=length(phenotype))
      stop('Comparison length is not equal to phenotype length.')
    # Check for presence of comparison units in phenotype
    for(i in 1:length(phenotype)){
      unique(as.character(sapply(comparison[[i]], function(x){
        gsub(' ', '', strsplit(x, split = '-')[[1]])}))) -> foo
      if(!all(foo %in% unique(pData(scu_meth@processed_dat[,phenotype[i],drop=TRUE]))))
        stop('Not all members of comparisons found in phenotype columns.')
    }
  }
  
  ### Start limma
  if(verbose >= 1)
    cat('Begin limma analysis.\n')
  
  # Create empty holder 'Res'
  Res <- list()
  
  # Check for presence of mode, if not populate the data
  if(limma_mode[1] == 'm'){
    if(length(scu_meth@m_values) == 0){
      cat('M-values not found. Populating.\n')
      m <- getM(scu_meth@processed_dat)[rownames(scu_meth@beta_values),]
      m <- pmin(m, max(m[m != 'Inf']))
      m <- pmax(m, min(m[m != '-Inf']))
      scu_meth@m_values <- m
    } 
    M <- scu_meth@m_values
  } else {
    M <- scu_meth@beta_values
  }
  # Find most variable probes
  if(top_probes){
    if(verbose >= 1)
      cat('Identifying most variable probes.\n')
    foo <- nrow(M) - round(top_probes*nrow(M))
    M <- dpsd(M, foo)
  }
  
  for(i in 1:length(phenotype)){
    if(verbose >= 1)
      cat('Running limma for phenotype ', i, '.\n', sep = '')
    foo_phenotype <- phenotype[i]
    foo_comparison <- comparison[[i]]
    foo_design <- .makedesigns(scu_meth, foo_phenotype)
    Res[[i]] <- .limma_mule(m = M,
                            des = foo_design,
                            comparison = foo_comparison)
  }
  names(Res) <- phenotype
  return(Res)
}

pipeline_methylation_microarray <- function(target_dir,
                                            phenotype,
                                            project_name = NULL,
                                            comparison = NULL,
                                            output_dir = 'output',
                                            pdat_file = NA,
                                            pval_file = NA,
                                            th_pval = NULL,
                                            th_proberate = 0.95,
                                            th_samplerate = 0.8,
                                            filter_samples = TRUE,
                                            filter_probes = TRUE,
                                            preprocess = c('funnorm', 'raw',
                                                           'illumina', 'quantile'),
                                            snp_filter = TRUE,
                                            snp_filter_dist = 2,
                                            snp_filter_maf = 0.05,
                                            getm = TRUE,
                                            top_probes = 0.95,
                                            limma_mode = c('m', 'b'),
                                            save_data = TRUE,
                                            .cache = TRUE,
                                            verbose = 2, ...){
  ## .cache is not implemented yet
  
  ### Complete checks before running pipeline
  
  if(verbose >= 1)
    cat('Checking for data consistency.\n')
  
  # Create output directory if it does not exist
  if(!dir.exists(output_dir)){
    dir.create(output_dir)
  }
  
  # Make project directory if missing, using date & timestamp
  if(missing(project_name)){
    project_name <- paste0('project_', gsub(' ', '_', gsub('-|:', '', Sys.time())))
  }
  project_dir <- paste0(output_dir, '/', project_name)
  
  if(!dir.exists(project_dir)){
    dir.create(project_dir)
  }
  
  # Evaluate if files exist
  foo <- c(pdat_file, pval_file)
  names(foo) <- c('pdat_file', 'pval_file')
  bar <- sapply(foo, function(x){
    if(is.na(x))
      return(TRUE) else
        return(file.exists(x))
  })
  if(!all(as.logical(bar))) {
    missing_files <- paste('The following files are not found:',
                           paste(names(which(!bar)), collapse = '; '), sep = '\n')
    stop(missing_files)
  }
  
  # Define thresholds
  if(missing(pval_file) & missing(th_pval)){
    th_pval <- 1e-8
  } else {
    if(missing(th_pval)) {
      th_pval <- 0.01
    }
  }
  
  if(top_probes < 0){
    stop('Top_probes must be a value between 0 to 1, with 0 = no filtering.')
  }
  
  if(verbose >= 1)
    cat('\nStep 1: Reading data.\n\n')
  scu_meth <- pmm_01_readqc(target_dir = target_dir, output_dir = project_dir,
                            pdat_file = pdat_file, pval_file = pval_file,
                            th_pval = th_pval, th_proberate = th_proberate,
                            th_samplerate = th_samplerate, verbose = verbose)
  
  if(verbose >= 1)
    cat('\nStep 2: Beginning normalization.\n\n')
  scu_meth <- pmm_02_normalization(scu_meth = scu_meth, output_dir = project_dir,
                                   filter_samples = filter_samples,
                                   filter_probes = filter_probes,
                                   preprocess = preprocess, verbose = verbose,
                                   snp_filter = snp_filter,
                                   snp_filter_dist = snp_filter_dist,
                                   snp_filter_maf = snp_filter_maf,
                                   getm = getm, ...)
  if(save_data){
    if(verbose >= 1)
      cat('Saving rdata, beta, and m-value files.\n')
    save(scu_meth, file = paste0(project_dir, '/scu_methylation.rda'))
    
    foo <- as.data.frame(scu_meth@beta_values)
    foo$probe_id <- rownames(foo)
    foo <- foo[,c(ncol(foo), 1:(ncol(foo) - 1))]
    write.table(foo, file = paste0(project_dir, '/beta_values.txt'),
                sep = '\t', quote = FALSE, row.names = FALSE)
    
    foo <- as.data.frame(scu_meth@beta_values)
    foo$probe_id <- rownames(foo)
    foo <- foo[,c(ncol(foo), 1:(ncol(foo) - 1))]
    write.table(foo, file = paste0(project_dir, '/m_values.txt'),
                sep = '\t', quote = FALSE, row.names = FALSE)
  }
    
  if(verbose >= 1)
    cat('\nStep 3: Running differential methylation analysis by limma.\n\n')
  scu_limma <- pmm_03_limma(scu_meth = scu_meth,
                            phenotype = phenotype,
                            comparison = comparison,
                            top_probes = top_probes,
                            limma_mode = limma_mode,
                            output_dir = project_name,
                            verbose = verbose)
  
  Res <- list()
  Res$scu_methylation <- scu_meth
  Res$scu_limma <- scu_limma
  return(Res)
}


# Use case ----------------------------------------------------------------

td = 'data/' # Contains .csv and folders with .idat
pf = 'data/phenodat.txt' # Contains phenotypic data with 'Sample_Name' identifiers

test = pmm_01_readqc(target_dir = td, pdat_file = pf)
test_norm = pmm_02_normalization(test, preprocess = 'funnorm')
test_limma = pmm_03_limma(scu_meth = scu_meth, limma_mode = 'b',
                          phenotype = 'phenotype')

test = pipeline_methylation_microarray(target_dir = td, pdat_file = pf,
                                       phenotype = 'phenotype', project = 'test')
