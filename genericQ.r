## Functions to get Q matrices for the SIR model from known start and end
## states, and for the Moran model over the entire state space.
##
## Chris Sherlock

#' Get Q matrix for SIR model
#'
#' Gets the rate matrix for an SIR model with known start and end states using the reduced statespace of Ho, Crawford and Suchard (2018): births of I and births of R.
#'
#' @param x0 initial vector: (S,I)
#' @param x1 final vector: (S,I)
#' @param thetas rate vector: (infection, recovery)
#' @return Q A rate matrix where the first state corresponds to x0 and the penultimate state is x1
#'
#' @author Chris Sherlock
#'
#' @example
#' Q=SIRgetQHCS(c(100,2),c(90,1),c(.01,.1))
#'
#' @export
SIRgetQHCS<-function(x0,x1,thetas) {
## Think of index 0 to bI and 0 to bR, with inner loop over bI
## But array indices in the R language start at 1.
    nbI=x0[1]-x1[1]
    nbR=nbI+x0[2]-x1[2]
    d=(nbI+1)*(nbR+1)
    maxnentry=(1+length(thetas))*d+1 ## 1 for the coffin
    rows=rep(0,maxnentry)
    cols=rep(0,maxnentry)
    vals=rep(0,maxnentry)
    nentries=0
    for (i in 0:nbI) {
        for (j in 0:nbR) {
            row=i+(nbI+1)*j+1
            S=x0[1]-i
            I=x0[2]+i-j
            sumrats=0
            if (I>0) {
                ratbI=thetas[1]*S*I
                ## S->S-1 and I->I+1 i.e. i->i+1
                rows[nentries+1]=row
                if (i<nbI) {
                    cols[nentries+1]=row+1
                }
                else { ## Coffin
                    cols[nentries+1]=d+1
                }
                vals[nentries+1]=ratbI
                nentries=nentries+1
                sumrats=sumrats+ratbI
            }
            if (I>0) {
                ratbR=thetas[2]*I
                rows[nentries+1]=row
                if (j<nbR) {
                    cols[nentries+1]=row+(nbI+1)
                }
                else {
                    cols[nentries+1]=d+1
                }
                vals[nentries+1]=ratbR
                nentries=nentries+1
                sumrats=sumrats+ratbR
            }
            if (sumrats>0) {
                rows[nentries+1]=row
                cols[nentries+1]=row
                vals[nentries+1]=-sumrats
                nentries=nentries+1
            }
        }
    }
    rows=rows[1:nentries]; cols=cols[1:nentries]; vals=vals[1:nentries]

    Q=Matrix::sparseMatrix(rows,cols,x=vals,dims=c(d+1,d+1))
    return(Q)
}


## Create matrix, the (i+1,j+1) element of which is the state index for
## i I births and j R births
## Removes all states where number of infected < 0
CreateLookup<-function(nbI,nbR,I0) {
    M=matrix(nrow=nbI+1,ncol=nbR+1,data=-1)
    currstart=1
    for (i in 0:nbI) {
        ninrow=min(1+i+I0,nbR+1)
        currend=currstart+ninrow-1
        M[i+1,1:ninrow]=currstart:currend
        currstart=currend+1
    }
    return(M)
}

#' Get Q matrix for SIR model
#'
#' Gets the rate matrix for an SIR model with known start and end states using the reduced statespace of Ho, Crawford and Suchard (2018) [births of I and births of R] but with the additional constraint in Sherlock (2020).
#'
#' @param x0 initial vector: (S,I)
#' @param x1 final vector: (S,I)
#' @param thetas rate vector: (infection, recovery)
#' @return Q A rate matrix where the first state corresponds to x0 and the penultimate state is x1
#'
#' @author Chris Sherlock
#'
#' @example
#' Q=SIRgetQslim(c(100,2),c(90,1),c(.01,.1))
#'
#' @export
SIRgetQslim<-function(x0,x1,thetas) {
    nbI=x0[1]-x1[1]
    nbR=nbI+x0[2]-x1[2]
    M=CreateLookup(nbI,nbR,x0[2])
    d=(nbI+1)*(nbR+1)-sum(M<0)
    maxnentry=(1+length(thetas))*d+1 ## 1 for the coffin
    rows=rep(0,maxnentry)
    cols=rep(0,maxnentry)
    vals=rep(0,maxnentry)
    nentries=0
    for (i in 0:nbI) {
        for (j in 0:nbR) {
            row=M[i+1,j+1]
            if (row>0) {
                S=x0[1]-i
                I=x0[2]+i-j
                sumrats=0
                if (I>0) {
                    ratbI=thetas[1]*S*I
                ## S->S-1 and I->I+1 i.e. i->i+1
                    rows[nentries+1]=row
                    if (i<nbI) {
                        cols[nentries+1]=M[i+2,j+1]
                    }
                    else { ## Coffin
                        cols[nentries+1]=d+1
                    }
                    vals[nentries+1]=ratbI
                    nentries=nentries+1
                    sumrats=sumrats+ratbI
                }
                if (I>0) {
                    ratbR=thetas[2]*I
                    rows[nentries+1]=row
                    if (j<nbR) {
                        cols[nentries+1]=M[i+1,j+2]
                    }
                    else {
                        cols[nentries+1]=d+1
                    }
                    vals[nentries+1]=ratbR
                    nentries=nentries+1
                    sumrats=sumrats+ratbR
                }
                if (sumrats>0) {
                    rows[nentries+1]=row
                    cols[nentries+1]=row
                    vals[nentries+1]=-sumrats
                    nentries=nentries+1
                }
            }
        }
    }
    rows=rows[1:nentries]; cols=cols[1:nentries]; vals=vals[1:nentries]

    Q=Matrix::sparseMatrix(rows,cols,x=vals,dims=c(d+1,d+1))
    return(Q)
}


#' Get Q matrix for Moran model
#'
#' Gets the rate matrix for the full statespace of the Moran model.
#'
#' @param thetas parameter vector: (alpha,beta,u,v,npop)
#' @return Q An (npop+1)*(npop+1) rate matrix
#'
#' @author Chris Sherlock
#'
#' @example
#' Q=MoranGetQ(c(1,.3,.2,.1,1000))
#'
#' @export
MoranGetQ<-function(thetas) {
    alpha=thetas[1]; beta=thetas[2]; u=thetas[3]; v=thetas[4]; npop=thetas[5]
    nentry=3*(npop+1)-2
    rows=rep(0,nentry); cols=rep(0,nentry); vals=rep(0,nentry)
    nentries=0

    for (i in 0:npop) { ## row
        f=1-i/npop
        rowsum=0
        if (i>0) {
            lambda=(1-f)*( alpha*f*(1-u) +beta*(1-f)*v)
            nentries=nentries+1
            rows[nentries]=i+1; cols[nentries]=i; vals[nentries]=lambda
            rowsum=rowsum+lambda
        }
        if (i<npop) {
            mu=f* (alpha*f*u + beta*(1-f)*(1-v))
            nentries=nentries+1
            rows[nentries]=i+1; cols[nentries]=i+2; vals[nentries]=mu
            rowsum=rowsum+mu            
        }
        nentries=nentries+1
        rows[nentries]=i+1; cols[nentries]=i+1; vals[nentries]=-rowsum
    }


    Q=Matrix::sparseMatrix(rows,cols,x=vals)
    return(Q)
}

#' Get Q matrix for BDI model for an interval + coffin
#'
#' Gets the rate matrix for the birth-death-immigration model with a state space of interval union coffin.
#'
#' @param thetas parameter vector: (lambda,mu,gamma)
#' @param interval vector: (nmin,nmax)
#' @return Q An (nmax-nmin+2)*(nmax-nmin+2) rate matrix - includes coffin state
#'
#' @author Chris Sherlock
#'
#' @example
#' Q=BDIGetQ(c(.01,.02,.2),c(45,50))
#'
#' @export
BDIGetQ<-function(thetas,interval) {
    lambda=thetas[1]; mu=thetas[2]; gamma=thetas[3];
    nmin=interval[1]; nmax=interval[2]; 
    nentry=3*(nmax-nmin+1)
    if (nmin==0) { 
        nentry=nentry-1 ## no deaths from x=0
    }    
    rows=rep(0,nentry); cols=rep(0,nentry); vals=rep(0,nentry)

    nentries=0
    nnoncoffin=nmax-nmin+1
    for (i in 1:nnoncoffin) { ## row from 1 ... 
        x=i-1+nmin ## actual state
        rowsum=0
        if (x>0) { ## allow death
            nentries=nentries+1
            rows[nentries]=i; cols[nentries]=i-1; vals[nentries]=mu*x
            if (i==1) { ## but if in lowest state this goes to the coffin
                cols[nentries]=nnoncoffin+1; 
            }
            rowsum=rowsum+vals[nentries]
            nentries=nentries+1
        }
        ## When i=nmax then movement is to the coffin state automatically
        rows[nentries]=i; cols[nentries]=i+1; vals[nentries]=lambda*x+gamma
        rowsum=rowsum+vals[nentries]            

        ## Add the diagonal entry
        nentries=nentries+1
        rows[nentries]=i; cols[nentries]=i; vals[nentries]=-rowsum
    }
    ## all rates are zero for the coffin row.    

    Q=Matrix::sparseMatrix(rows,cols,x=vals,dims=c(nmax-nmin+2,nmax-nmin+2))
    return(Q)
}

