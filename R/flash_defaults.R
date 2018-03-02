flash_default_ebnm_param = function(ebnm_fn) {
    if (identical(ebnm_fn, ebnm_ash)) {
        return(list(outputlevel = 5, mixcompdist = "normal", method = "shrink"))
    } else if (identical(ebnm_fn, ebnm_pn)) {
        return(list())
    } else {
        stop("no defaults available for ebnm_param for that ebnm function, please supply them")
    }
}
