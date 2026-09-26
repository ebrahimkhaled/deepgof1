"""Review 2 check: finite-n null sd of the projection statistic T = e'Ae/n^2 (printed weight, raw covariates,
fitted logistic model, Bernoulli y) against its U-statistic limit nu/(2 pi) from PROJ_LIU_RAW.json 'base'
(nu_F0 = 0.00488, 0.00482 at k = 8; 0.0274, 0.0269 at k = 0, niche design c_irr = 0.15). Four covariate samples
per (n, k), 400 null replicates each."""
import os
os.environ["OMP_NUM_THREADS"]="2"; os.environ["OPENBLAS_NUM_THREADS"]="2"
import numpy as np
from scipy.special import expit
def liu(x):
    n=len(x); A=np.zeros((n,n))
    for l in range(n):
        D=x-x[l]; nm=np.sqrt((D*D).sum(1)); z=nm<=1e-12; U=D/np.where(z,1,nm)[:,None]
        K=(np.pi-np.arccos(np.clip(U@U.T,-1,1)))/(2*np.pi); np.fill_diagonal(K,.5)
        K[z,:]=.5; K[:,z]=.5; K[np.ix_(z,z)]=1; A+=K
    return A
def fit(X,y):
    b=np.zeros(X.shape[1])
    for _ in range(25):
        p=expit(X@b); W=p*(1-p); b=b+np.linalg.solve(X.T@(W[:,None]*X), X.T@(y-p))
    return expit(X@b)
LIM={8:np.mean([0.00488,0.00482]),0:np.mean([0.0274,0.0269])}
for n,k in ((200,8),(500,8),(200,0),(500,0)):
    out=[]
    for s in range(4):
        rng=np.random.default_rng(100*s+k+n); x=rng.normal(size=(n,2+k)); X=np.column_stack([np.ones(n),x])
        pi=expit(-0.5+x[:,0]+0.8*x[:,1]+0.15*x[:,2:].sum(1)); A=liu(x); d=np.diag(A); Ao=A-np.diag(d)
        T=[];To=[]
        for r in range(400):
            y=(rng.random(n)<pi).astype(float); e=y-fit(X,y); T.append(e@A@e/n**2)
        out.append(2*np.pi*np.std(T,ddof=1))
    out=np.array(out)
    print("n=%d k=%d: 2pi*sd(T) = %.5f (samples %s) ; limit nu_F = %.5f ; ratio %.2f"%(n,k,out.mean(),np.round(out,5),LIM[k],out.mean()/LIM[k]))
