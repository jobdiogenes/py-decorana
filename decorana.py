import numpy as np

import numba

@numba.jit
def is_sorted(a):
    for i in range(a.size-1):
         if a[i+1] < a[i] :
               return False
    return True

@numba.jit
def downweigth(veg, fraction=5):
    CONST1 = 1e-10
    if fraction<1:
        fraction= 1/fraction
    if not isinstance(veg,np.ndarray):
        raise TypeError("Matrix must be a numpy array type")
    yeig1 = np.sum(veg,axis=1)
    y2 = np.sum(veg**2,axis=1)
    y2 = yeig1**2/y2
    amax = np.max(y2)/fraction
    v = np.repeat(1,veg.shape[1])
    downers = y2 < amax
    v[downers] = (y2/amax)[downers]
    
@numba.jit          
def decorana(veg, iweigh = 0, iresc = 4, ira = 0, mk = 26, short = 0,
                before = None, after = None):

    CONST1 = 1e-10
    CONST2 = 5
    CONST3 = 1e-11
    if not isinstance(veg,np.ndarray):
       raise TypeError("Values Matrix must be a numpy array type")
    if veg.ndim!=2:
       raise TypeError("Values is not a matrix")
    
    aidot = np.sum(veg,axis=0)
    if np.any(aidot <= 0):
        raise ValueError("all row sums must be >0 in the community matrix: remove empty sites")
    if np.any(veg < 0):
        raise ValueError("'decorana cannot handle negative data entries")
    adotj = np.sum(veg,axis=1)
    if np.any(adotj <= 0):
        print("some species were removed because they were missing in the data")
    if mk < 10:
        mk = 10
    if mk > 46:
        mk = 46
    if ira==1:
        iresc = 0
    if before is not None:
        if not is_sorted(before):
            raise ValueError("'before' must be sorted")
        if len(before) != len(after):
            raise ValueError("'before' and 'after' must have same lengths")
        
        for i in range(np.shape(veg)[0]): 
            tmp = veg[i,veg[i,:] > 0]
            veg[i, tmp] = approx(before, after, veg[i, tmp],rule = 2)
            
    if iweigh:
        veg = downweight(veg, CONST2)
        aidot = np.sum(veg,axis=0)
        adotj = np.sum(veg,axis=1)
    
    v <- attr(veg, "v")
    v.fraction <- attr(veg, "fraction")
    adotj[adotj < Const3] <- Const3
    CA = do_decorana, veg, ira, iresc, short, mk, as.double(aidot), as.double(adotj)
    
    if ira:
        dnames <- paste("RA", 1:4, sep = "")
    else:
         dnames <- paste("DCA", 1:4, sep = "")
    dimnames(CA$rproj) <- list(rownames(veg), dnames)
    dimnames(CA$cproj) <- list(colnames(veg), dnames)
    names(CA$evals) <- dnames
    origin <- apply(CA$rproj, 2, weighted.mean, aidot)
    if (ira) {
        evals.decorana <- NULL
    }
    else {
        evals.decorana <- CA$evals
        var.r <- diag(cov.wt(CA$rproj, aidot, method = "ML")$cov)
        var.c <- diag(cov.wt(CA$cproj, adotj, method = "ML")$cov)
        CA$evals <- var.r/var.c
        if (any(ze <- evals.decorana <= 0))
            CA$evals[ze] <- 0
    }
    additems <- list(evals.decorana = evals.decorana, origin = origin,
                     v = v, fraction = v.fraction, iweigh = iweigh,
                     before = before, after = after,
                     call = match.call())
    CA <- c(CA, additems)
    class(CA) <- "decorana" # c() strips class
    CA

