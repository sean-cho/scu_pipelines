# Methylation -------------------------------------------------------------

library(minfi)

setClassUnion('SCU_processed_dat', 
              c('GenomicRatioSet', 'MethylSet'))

setClass(Class = 'SCU_Methylation', 
         representation = representation(
           rgdat = 'RGChannelSet',
           processed_dat = 'SCU_processed_dat',
           detp = 'matrix', 
           probe_callrate = 'numeric',
           sample_callrate = 'numeric',
           callrate_thresholds = 'numeric',
           pval_threshold = 'numeric',
           beta_values = 'matrix',
           m_values = 'matrix'))

setClass(Class = 'SCU_Methylation_Limma', 
         representation = representation(
           model.mat = 'matrix',
           model = 'character',
           linear_fit = 'MArrayLM', 
           contrasts_fit = 'MArrayLM',
           ebayes = 'MArrayLM',
           toptable  = 'data.frame'))

setMethod('show', signature = signature(object = 'SCU_Methylation'), 
          function(object){
            show(object@rgdat)
            foo <- object@processed_dat
            if(length(foo) != 0)
              show(foo) else cat('No processed dataset.\n')
            cat('SCU_Extensions\n')
            
            cat('  pval_threshold:', object@pval_threshold, '\n')
            cat('  callrate_thresholds; sample: ', 
                object@callrate_thresholds[1],
                ', probe:', object@callrate_thresholds[2], 
                '\n', sep = '')
            
            cat('  ')
            if(length(object@detp) != 0) 
              cat('detp: TRUE\n') else cat('detp: FALSE\n')
            
            cat('  ')
            if(length(object@probe_callrate) != 0) 
              cat('probe_callrate: TRUE\n') else cat('probe_callrate: FALSE\n')
            
            cat('  ')
            if(length(object@sample_callrate) != 0)
              cat('sample_callrate: TRUE\n') else cat('sample_callrate: FALSE\n')
            
            cat('  ')
            if(length(object@beta_values) != 0)
              cat('beta_values: TRUE\n') else cat('beta_values: FALSE\n')
            
            cat('  ')
            if(length(object@m_values) != 0)
              cat('M_values: TRUE\n') else cat('M_values: FALSE\n')
          })

setMethod('show', signature = signature(object = 'SCU_Methylation_Limma'), 
          function(object){
            cat('Model: ', paste(object@model, collapse = ' ; '), '\n', sep = '')
            cat('Top 10 DMPs:\n')
            print(head(object@toptable))
          })