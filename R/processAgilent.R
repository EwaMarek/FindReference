# na wejscie obiekt klasy EListRaw lub marrayRaw
# na wyjscie dane po korekcji tla, normalizacji i sumaryzacji

processAgilent = function(a){


  if(class(a) == 'marrayRaw'){
    # background correction
    a_back = list(red = a@maRf-a@maRb , green = a@maGf-a@maGb)
    a_back = lapply(a_back, function(x){x[x<0] = 0.0000001; x})

    # vsn normalization
    a_norm = lapply(a_back, vsn2, verbose = FALSE)
    a_norm = lapply(a_norm, exprs)

    #annotation
    controls = a@maGnames@maInfo$ControlType
    annot = annotAgilent(a@maGnames@maInfo, controls)
    a_norm = lapply(a_norm, function(x, controls, annot){x[which(controls==0), ]}, controls=controls, annot=annot)
    a_norm = lapply(a_norm, function(x, annot){row.names(x) = annot; x}, annot=annot)
    a_norm = lapply(a_norm, function(x){x[which(rownames(x)!="NA"), ]})

    # summarization
    a_out = rep(list(list()), 2)
    a_out[[1]] = as.matrix(aggregate.Matrix(a_norm[[1]], groupings=row.names(a_norm[[1]]), fun='mean'))
    a_out[[2]] = as.matrix(aggregate.Matrix(a_norm[[2]], groupings=row.names(a_norm[[2]]), fun='mean'))


  }else if(class(a) == 'EListRaw'){
    # background correction
    a_back = a$E-a$Eb
    a_back[a_back<0] = 0.0000001

    # vsn normalization
    a_norm = vsn2(a_back, verbose = FALSE)
    a_norm = exprs(a_norm)

    #annotation
    controls = a@.Data[[4]]$ControlType
    annot = annotAgilent(a@.Data[[4]], controls)
    a_norm = a_norm[which(controls==0), ]
    row.names(a_norm) = annot
    a_norm = a_norm[which(rownames(a_norm)!="NA"), ]

    # summarization
    a_out = as.matrix(aggregate.Matrix(a_norm, groupings=row.names(a_norm), fun='mean'))

  }else if(class(a_out) == 'try-error'){
    a_out == 'There were no data to normalize'
    warning("There are no data to normalize.")

  }else{
    a_out = "This object class is not supported for chosen platform."
    warning(paste("This object class (", as.character(class(a)), ") could not be normalized with norm_and_annot_function with platform argument Agilent.
                  Supported classes for this function for Agilent platform are ElistRaw, marrayRaw or list of ElistRaw and marrayRaw objects.", sep=' '))
  }

  return(a_out)

}
