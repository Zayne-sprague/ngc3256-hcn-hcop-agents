# Q4c-4: differential N vs S fit. Which single parameter must change between the N and S nuclei to explain the full data set (ladders, R32, R10, 3-2/1-0 per molecule)?
# Family: both thick (tau32>1), common-(T,n) model, kappa=1, physT30 prior; data = 'all' (8 per nucleus; g free per nucleus). Options: nothing changes; kappa(S) free; N(HCN) differs; N(HCO+) differs;
# both columns differ; n differs; Tkin differs; everything independent. Output out/q4c/ns_options.csv
import sys; sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd
import q4c_joint as J
def diagcube(c,mode,thick=True):
    d,use,chi,cls,aux=J.cube(c,mode); sh=aux['sh']; n=np.arange(J.nn)
    msk=cls['A both thick, kappa=1'][0] if thick else None
    out=[]
    for ik in (len(J.KAP)-1,):
        x=chi[ik][:,n,n,:,:]; 
        if thick: x=np.where(msk[:,n,n,:,:],x,1e9)
        out.append(x)
    xk=chi[:,:,n,n][:, :, :, :, :] if False else chi[:,:,n,n,:,:]
    if thick: xk=np.where(msk[None][:,:,n,n,:,:],xk,1e9)
    return out[0],xk.min(0),int(use.sum()),xk.argmin(0)
rows=[]
for mode in ('all','all_no21'):
    for thick in (True,False):
        XN,_,nd,_=diagcube('N',mode,thick); XS,XSk,_,ak=diagcube('S',mode,thick)
        T,n,NH,NO=XN.shape; Tk,Lnn,LNc=J.TK,J.LN,J.LNC
        def info(a,x): i=np.unravel_index(np.argmin(x),x.shape); return i
        res={}
        s=XN+XS; res['nothing changes (4 par)']=(s.min(),4)
        res['kappa(S) differs (5)']=((XN+XSk).min(),5)
        res['N(HCN) differs (5)']=((XN.min(2,keepdims=True)+XS.min(2,keepdims=True)).min(),5)
        res['N(HCO+) differs (5)']=((XN.min(3,keepdims=True)+XS.min(3,keepdims=True)).min(),5)
        res['both columns differ (6)']=((XN.min((2,3),keepdims=True)+XS.min((2,3),keepdims=True)).min(),6)
        res['n differs (5)']=((XN.min(1,keepdims=True)+XS.min(1,keepdims=True)).min(),5)
        res['Tkin differs (5)']=((XN.min(0,keepdims=True)+XS.min(0,keepdims=True)).min(),5)
        res['Tkin and n differ, columns common (6)']=((XN.min((0,1),keepdims=True)+XS.min((0,1),keepdims=True)).min(),6)
        res['independent (8)']=(XN.min()+XS.min(),8)
        base=res['nothing changes (4 par)'][0]
        for k,(c2,npar) in res.items():
            rows.append(dict(mode=mode,family='both thick' if thick else 'any tau (diag)',option=k,chi2=c2,npar=npar,ndata_total=2*nd,dchi2_vs_nothing=base-c2,BIC=c2+npar*np.log(2*nd),dof=2*nd-npar))
        # best parameters for the selected options
        i=np.unravel_index(np.argmin(XN.min(2,keepdims=True)+XS.min(2,keepdims=True)),(T,n,1,NO)); 
        jN=XN[i[0],i[1],:,i[3]].argmin(); jS=XS[i[0],i[1],:,i[3]].argmin()
        rows.append(dict(mode=mode,family='both thick' if thick else 'any tau (diag)',option='  best N(HCN) option: T,logn,logNO,logNH_N,logNH_S = %.0f %.2f %.2f %.2f %.2f'%(Tk[i[0]],Lnn[i[1]],LNc[i[3]],LNc[jN],LNc[jS])))
        i=np.unravel_index(np.argmin(XN+XSk),XN.shape); 
        rows.append(dict(mode=mode,family='both thick' if thick else 'any tau (diag)',option='  best kappa option: T,logn,logNH,logNO,kappaS = %.0f %.2f %.2f %.2f %.3f'%(Tk[i[0]],Lnn[i[1]],LNc[i[2]],LNc[i[3]],J.KAP[ak[i]])))
        i=np.unravel_index(np.argmin(s),s.shape)
        rows.append(dict(mode=mode,family='both thick' if thick else 'any tau (diag)',option='  best common: T,logn,logNH,logNO = %.0f %.2f %.2f %.2f'%(Tk[i[0]],Lnn[i[1]],LNc[i[2]],LNc[i[3]])))
R=pd.DataFrame(rows); R.round(3).to_csv(J.O+'ns_options.csv',index=False)
pd.set_option('display.width',250,'display.max_colwidth',90); print(R.round(2).to_string())
# profile of the HCN column shift between N and S (all other parameters common): chi2(Delta = logN_HCN(S) - logN_HCN(N))
for mode in ('all','all_no21'):
    XN,_,nd,_=diagcube('N',mode); XS,XSk,_,_=diagcube('S',mode); nN=XN.shape[2]; prof={}
    for sh in range(-6,7):
        best=1e30
        for a in range(nN):
            b=a+sh
            if 0<=b<nN: best=min(best,float((XN[:,:,a,:]+XS[:,:,b,:]).min()))
        prof[sh*0.25]=best
    m=min(prof.values()); ok=[k for k,v in prof.items() if v-m<1.0]; ok4=[k for k,v in prof.items() if v-m<4.0]
    print(mode,'profile',{k:round(v-m,1) for k,v in prof.items()},'| 1sigma range of Delta logN_HCN(S-N):',min(ok),max(ok),'| 2sigma',min(ok4),max(ok4))
    pd.DataFrame({'delta_logNHCN_S_minus_N':list(prof),'chi2':list(prof.values())}).to_csv(J.O+f'ns_profile_NHCN_{mode}.csv',index=False)
