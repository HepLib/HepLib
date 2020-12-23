*
*  -------------------------------------------------------
*
*       qgraf-3.1.4
*       a module for Feynman diagram generation
*
*       design and copyright by paulo.nogueira@ist.utl.pt
*
*       reference:
*         [1] J. Comput. Phys. 105 (1993) 279-289
*
*       documentation:
*         [2] files 'qgraf-3.0.pdf' and 'qgraf-3.1.pdf'
*
*  -------------------------------------------------------
*
      program qgraf
      implicit integer(a-z)
      save
      parameter ( maxleg=10, maxrho=6, maxi=10 )
      parameter ( maxdeg=6 )
      parameter ( maxn=maxi+maxi-2, maxli=maxn+maxrho-1 )
      parameter ( ipar1= maxleg*maxleg-maxleg )
      parameter ( ipar2= ipar1/2 )
      parameter ( maxtak=2 )
      parameter ( srec=81 )
      parameter ( lfeed=10, dquote=34, rquote=39 )
      parameter ( sibuff=1048576, scbuff=524288, mimal=2, mcmal=2 )
      common/z1in/leg(1:maxleg),nleg,nloop,incom
      common/z2in/nv(1:maxdeg)
      common/z3in/momep(0:maxleg),momel(0:maxleg)
      common/z4in/mflag(1:21)
      common/z6in/dunit,munit,ounit,sunit,funit
      common/z7in/lsfile,sfilea,sfileb
      common/z8in/lmfile,mfilea,mfileb
      common/z9in/lofile,ofilea,ofileb
      common/z10in/lffile,ffilea,ffileb
      common/z1g/nrho,rho(1:maxdeg),g(1:maxn,1:maxn)
      common/z2g/dsym,dis,ndiag
      common/z3g/vparto(0:0),vdeg(0:0),nvert
      common/z4g/n,nli
      common/z5g/psym(0:0),psyms,nsym
      common/z6g/p1(1:maxleg),invp1(1:maxleg)
      common/z7g/lmap(1:maxn,1:maxdeg),vmap(1:maxn,1:maxdeg),
     :pmap(1:maxn,1:maxdeg),vlis(1:maxn),invlis(1:maxn)
      common/z8g/degree(1:maxn)
      common/z9g/tpc(0:0)
      common/z10g/p1l(1:ipar2),p1r(1:ipar2),ns1
      common/z11g/npart,nblok,nprop
      common/z13g/xn(1:maxn)
      common/z14g/zcho(0:maxli),zbri(0:maxli),zpro(0:maxli),
     :rbri(0:maxli),sbri(0:maxli)
      common/z15g/dpntro(1:maxdeg)
      common/z16g/rdeg(1:maxn),amap(1:maxn,1:maxdeg)
      common/z17g/xtail(1:maxn),xhead(1:maxn),ntadp
      common/z18g/eg(1:maxn,1:maxn),flow(1:maxli,0:maxleg+maxrho)
      common/z19g/vfo(1:maxn)
      common/z20g/tftype(0:0),tfnarg(0:0),tfa(0:0),tfb(0:0),tfo(0:0),ntf
      common/z21g/punct1(0:127)
      character*(srec) auxlin
      common/z22g/auxlin
      common/z23g/tak(1:maxtak),ks
      common/z24g/iogp(1:4)
      common/z25g/ex(1:maxli),ey(1:maxli),ovm(1:maxn,1:maxdeg)
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z32g/acf0(0:127),acf1(0:127),sucpal(0:11,0:11)
      common/z33g/namep(0:0),namel(0:0)
      common/z34g/link(0:0),antiq(0:0),blok(0:0)
      common/z35g/stibs(0:mimal),stcbs(0:mcmal),imal(0:0),cmal(0:0)
      common/z36g/udkp(0:2),udkl(0:0),udkt(0:0),udki(0:0)
      common/z37g/drecp(0:0),drecl(0:0),drecii(0:0),irecc(0:0),
     :frecc(0:0),ndrec,ncom
      common/z38g/gmkp(0:0),gmkl(0:0),gmkd(0:0),gmko(0:0),
     :gmkvp(0:0),gmkvl(0:0)
      common/z39g/pmkp(0:0),pmkl(0:0),pmkd(0:0),pmkvpp(0:0),pmkvlp(0:0)
      common/z40g/vmkp(0:0),vmkl(0:0),vmkvpp(0:0),vmkvlp(0:0),vmks(0:0)
      common/z41g/nudk,ngmk,npmk,nvmk,tgmkd,tpmkd
      common/z43g/vmkr(0:0),vmkmao(0:0),vmkmio(0:0)
      common/z44g/sdiap,sdial,sgrap,sgral,stotp,stotl
      common/z45g/qdatp,qdatl,qvp,qvl,dprefp,dprefl
      integer lista(maxn+maxdeg),xlist(1:maxdeg)
      integer slist(2*maxli)
      i=srec
      if(i.lt.81)then
      auxlin(1:srec)='parameter "srec" is not properly set'
      call messag(1,0,0,0)
      endif
      i=maxdeg
      if(i.lt.3)then
      auxlin(1:srec)='parameter "maxdeg" is not properly set'
      call messag(1,0,0,0)
      endif
      i=maxleg
      if(i.lt.1)then
      auxlin(1:srec)='parameter "maxleg" is not properly set'
      call messag(1,0,0,0)
      endif
      i=maxrho
      if(i.lt.0)then
      auxlin(1:srec)='parameter "maxrho" is not properly set'
      call messag(1,0,0,0)
      endif
      i=maxleg
      j=maxrho
      if(i.lt.4.and.j.lt.1)then
      auxlin(1:srec)='parameters "maxdeg"/"maxrho" are too small'
      call messag(1,0,0,0)
      endif
      i=maxi
      j=max(maxleg,maxrho)
      if(i.lt.j)then
      auxlin(1:srec)='parameter "maxi" is not properly set'
      call messag(1,0,0,0)
      endif
      i=maxn
      j=2*max(maxleg,maxrho)-2
      if(i.lt.j)then
      auxlin(1:srec)='parameter "maxn" is not properly set'
      call messag(1,0,0,0)
      endif
      i=maxli
      j=2*max(maxleg,maxrho)+maxrho-3
      if(i.lt.j)then
      auxlin(1:srec)='parameter "maxli" is not properly set'
      call messag(1,0,0,0)
      endif
      i=sibuff
      if(i.lt.1)then
      auxlin(1:srec)='parameter "sibuff" is not properly set'
      call messag(1,0,0,0)
      endif
      i=scbuff
      if(i.lt.1)then
      auxlin(1:srec)='parameter "scbuff" is not properly set'
      call messag(1,0,0,0)
      endif
      i=mimal
      if(i.lt.1)then
      auxlin(1:srec)='parameter "mimal" is not properly set'
      call messag(1,0,0,0)
      endif
      i=mcmal
      if(i.lt.1)then
      auxlin(1:srec)='parameter "mcmal" is not properly set'
      call messag(1,0,0,0)
      endif
      call init0
      call inputs
      if(mflag(16).ne.0)then
      call iki
      call prepos(1)
      endif
      if(npart.eq.1)then
      mflag(12)=1
      endif
      do 06 i=1,npart
      if(stib(tpc(0)+i).eq.1)then
      mflag(10)=1
      elseif(stib(tpc(0)+i).eq.5)then
      mflag(11)=1
      endif
   06 continue
      if(nloop.eq.0.or.mflag(1).eq.1)then
      mflag(10)=0
      elseif(mflag(2).eq.1.or.mflag(3).eq.1)then
      mflag(10)=0
      endif
      if(nleg.eq.1)then
      if(stib(tpc(0)+leg(1)).eq.1)then
      mflag(20)=1
      endif
      endif
      if(mflag(20).eq.1)then
      goto 94
      endif
      hh=0
      gg=nleg-2+2*nloop
      rho(1)=nleg
      rho(2)=0
      m=gg
      do 44 i=nrho,3,-1
      rho(i)=m/(i-2)
      m=m-rho(i)*(i-2)
   44 continue
      goto 50
  182 continue
      do 20 i=4,nrho
      if(rho(i).gt.0)then
      goto 30
      endif
   20 continue
      goto 94
   30 continue
      rho(i)=rho(i)-1
      m=i-2+rho(3)
      do 91 j=i-1,3,-1
      rho(j)=m/(j-2)
      m=m-rho(j)*(j-2)
   91 continue
   50 continue
      do 97 i=3,nrho
      if(nv(i).eq.0.and.rho(i).gt.0)then
      goto 182
      endif
   97 continue
      nili=-nleg
      ii=0
      do 52 i=3,nrho
      nili=nili+i*rho(i)
      ii=ii+rho(i)
   52 continue
      nili=nili/2
      cntr10=0
      if(zpro(nili).eq.0)then
      goto 41
      endif
      do 81 i=nloop,nili
      if(zcho(i).ne.0.and.zbri(nili-i).ne.0)then
      do 79 j=0,nili-i
      if(rbri(j).ne.0.and.sbri(nili-i-j).ne.0)then
      goto 89
      endif
   79 continue
      endif
   81 continue
      goto 41
   89 continue
      if(mflag(21).ne.0)then
      do 71 i=1,ntf
      if(abs(stib(tftype(0)+i)).eq.6)then
      i1=0
      i2=0
      jj=stib(stib(tfo(0)+i)+1)
      if(stib(vmks(0)+jj).ne.1)then
      auxlin(1:srec)='  internal inconsistency'
      call messag(1,0,0,0)
      endif
      do 64 j=nrho,3,-1
      if(rho(j).gt.0)then
      i1=i1+rho(j)*stib(stib(vmkmio(0)+jj)+j)
      i2=i2+rho(j)*stib(stib(vmkmao(0)+jj)+j)
      endif
   64 continue
      if(stib(tftype(0)+i).gt.0)then
      if(i1.gt.stib(tfb(0)+i).or.i2.lt.stib(tfa(0)+i))then
      goto 41
      endif
      else
      if(i1.ge.stib(tfa(0)+i).and.i2.le.stib(tfb(0)+i))then
      goto 41
      endif
      endif
      endif
   71 continue
      endif
      cntr10=1
   40 continue
      call gen10(cntr10)
   41 continue
      if(cntr10.eq.0)then
      i0=stcbs(cmal(0))
      jj=2
      call vaocb(jj)
      ii=i0+jj
      stcb(i0+1:ii)='  '
      do 321 i=1,nrho
      if(nv(i).gt.0)then
      if(rho(i).gt.0)then
      jj=jj+wztos(i)+1
      call vaocb(jj)
      call karat(i,ii,stcb,scbuff,0)
      ii=i0+jj
      stcb(ii:ii)='^'
      jj=jj+wztos(gg/(i-2))+2
      call vaocb(jj)
      call karat(rho(i),ii,stcb,scbuff,0)
      i1=wztos(rho(i))
      do 307 j=ii+i1+1,i0+jj
      stcb(j:j)=' '
  307 continue
      ii=i0+jj
      else
      k=wztos(i)+1
      kk=k+wztos(gg/(i-2))
      jj=jj+kk+2
      call vaocb(jj)
      do 309 i1=ii+1,i0+jj
      stcb(i1:i1)=' '
  309 continue
      if(mod(kk,2).eq.1)then
      k=ii+(kk+1)/2
      stcb(k:k)='-'
      else
      if(kk.lt.2*k)then
      k=ii+1+kk/2
      stcb(k:k)='-'
      elseif(kk.gt.2*k)then
      k=ii-1+kk/2
      stcb(k:k)='-'
      else
      k=ii+kk/2
      stcb(k:k)='-'
      endif
      endif
      ii=i0+jj
      endif
      endif
  321 continue
      jj=jj+5
      call vaocb(jj)
      stcb(ii+1:i0+jj)='---  '
      ii=i0+jj
      jj=jj+wztos(ndiag-hh)
      call vaocb(jj)
      call karat(ndiag-hh,ii,stcb,scbuff,0)
      ii=i0+jj
      if(mflag(8).eq.0)then
      kk=sdial
      if(hh.eq.ndiag-1)then
      kk=kk-1
      endif
      jj=jj+sdial
      call vaocb(jj)
      stcb(ii+1:ii+sdial)=stcb(sdiap:sdiap-1+kk)
      else
      kk=sgral
      if(hh.eq.ndiag-1)then
      kk=kk-1
      endif
      jj=jj+sdial
      call vaocb(jj)
      stcb(ii+1:ii+jj)=stcb(sgrap:sgrap-1+kk)
      endif
      ii=i0+jj
      print *,stcb(i0+1:ii)
      hh=ndiag
      goto 182
      endif
      cntr10=0
      do 16 i=1,rho(1)
      pmap(i,1)=leg(p1(i))
      pmap(vmap(i,1),lmap(i,1))=stib(link(0)+leg(p1(i)))
   16 continue
      vind=rho(1)
   57 continue
      vind=vind+1
      vv=vlis(vind)
      vfo(vv)=stib(dpntro(degree(vv))+pmap(vv,1))
  100 continue
      do 49 j=1,rdeg(vv)
      aux=stib(vfo(vv)+j)-pmap(vv,j)
      if(aux.gt.0)then
      goto 104
      elseif(aux.lt.0)then
      vfo(vv)=vfo(vv)+degree(vv)+1
      goto 100
      endif
   49 continue
      goto 163
   58 continue
      vfo(vv)=vfo(vv)+degree(vv)+1
      do 56 j=1,rdeg(vv)
      if(stib(vfo(vv)+j).ne.pmap(vv,j))then
      goto 104
      endif
   56 continue
  163 continue
      do 66 j=rdeg(vv)+1,rdeg(vv)+g(vv,vv),2
      if(stib(vfo(vv)+j).ne.stib(link(0)+stib(vfo(vv)+j+1)))then
      goto 58
      endif
   66 continue
      do 76 j=rdeg(vv)+1,degree(vv)
      pmap(vv,j)=stib(vfo(vv)+j)
   76 continue
      goto 152
  104 continue
      vind=vind-1
      if(vind.eq.rho(1))then
      goto 40
      else
      vv=vlis(vind)
      goto 58
      endif
  152 continue
      do 86 i=rdeg(vv)+1,rdeg(vv)+g(vv,vv),2
      if(pmap(vv,i).gt.pmap(vv,i+1))then
      goto 58
      endif
   86 continue
      do 96 i=rdeg(vv)+1,rdeg(vv)+g(vv,vv)-2,2
      if(pmap(vv,i).gt.pmap(vv,i+2))then
      goto 58
      endif
   96 continue
      do 106 i=rdeg(vv)+g(vv,vv)+1,degree(vv)-1
      if(vmap(vv,i).eq.vmap(vv,i+1))then
      if(pmap(vv,i).gt.pmap(vv,i+1))then
      goto 58
      endif
      endif
  106 continue
      if(mflag(11).eq.1)then
      do 133 i=rdeg(vv)+1,degree(vv)
      if(stib(tpc(0)+pmap(vv,i)).eq.5)then
      goto 58
      endif
  133 continue
      endif
      if(mflag(10).eq.1)then
      do 136 i=1,ntadp
      if(xtail(i).eq.vv)then
      do 116 j=1,degree(vv)
      if(xhead(i).eq.vmap(vv,j))then
      if(stib(tpc(0)+pmap(vv,j)).eq.1)then
      goto 58
      endif
      endif
  116 continue
      elseif(xhead(i).eq.vv)then
      do 124 j=1,degree(vv)
      if(xtail(i).eq.vmap(vv,j))then
      if(stib(tpc(0)+pmap(vv,j)).eq.1)then
      goto 58
      endif
      endif
  124 continue
      endif
  136 continue
      endif
      do 146 i=rdeg(vv)+g(vv,vv)+1,degree(vv)
      pmap(vmap(vv,i),lmap(vv,i))=stib(link(0)+pmap(vv,i))
  146 continue
      if(vind.lt.n)then
      goto 57
      endif
      if(mflag(9).eq.0.or.nloop.eq.0)then
      goto 333
      endif
      do 14 i=1,n
      lista(i)=0
   14 continue
      do 74 i=1,n
      if(lista(i).eq.0)then
      ii=i
      kk=1
  820 continue
      do 24 j=1,degree(ii)
      if(j.ne.lista(ii))then
      jj=pmap(ii,j)
      if(stib(antiq(0)+jj).eq.1)then
      k=ii
      ii=vmap(k,j)
      lista(ii)=lmap(k,j)
      goto 830
      endif
      endif
   24 continue
      goto 74
  830 continue
      if(ii.gt.rho(1))then
      if(ii.ne.i)then
      kk=1-kk
      goto 820
      elseif(kk.eq.1)then
      goto 58
      endif
      endif
      endif
   74 continue
  333 continue
      dsym=1
      if(mflag(12).eq.0)then
      jk=psym(0)-rho(1)
      do 206 ii=1,nsym-1
      do 196 i=rho(1)+1,n
      ij=stib(jk+i)
      j=rdeg(i)+1
  400 continue
      if(j.le.degree(i))then
      kk=vmap(i,j)
      jj=1
  410 continue
      if(jj.le.degree(i).and.vmap(ij,jj).ne.stib(jk+kk))then
      jj=jj+1
      goto 410
      endif
      do 156 k=1,g(i,kk)
      lista(k)=pmap(ij,jj+k-1)
  156 continue
      do 176 k=1,g(i,kk)-1
      do 166 ik=k+1,g(i,kk)
      if(lista(k).gt.lista(ik))then
      aux=lista(k)
      lista(k)=lista(ik)
      lista(ik)=aux
      endif
  166 continue
  176 continue
      do 186 k=1,g(i,kk)
      aux=lista(k)-pmap(i,j)
      if(aux.lt.0)then
      goto 58
      elseif(aux.gt.0)then
      goto 339
      endif
      j=j+1
  186 continue
      goto 400
      endif
  196 continue
      dsym=dsym+1
  339 continue
      jk=jk+(n-rho(1))
  206 continue
      else
      dsym=nsym
      endif
      do 314 i=1,ntf
      if(stib(tfnarg(0)+i).eq.0)then
      auxlin(1:srec)='main_1'
      call messag(1,0,0,0)
      endif
      if(abs(stib(tftype(0)+i)).eq.1)then
      ii=0
      do 319 j=nleg+1,n
      do 315 jj=rdeg(j)+1,rdeg(j)+g(j,j),2
      do 311 ij=1,stib(tfnarg(0)+i)
      if(stib(stib(tfo(0)+i)+ij).eq.pmap(j,jj))then
      ii=ii+1
      elseif(stib(link(0)+stib(stib(tfo(0)+i)+ij)).eq.
     :pmap(j,jj))then
      ii=ii+1
      endif
  311 continue
  315 continue
      do 317 jj=rdeg(j)+g(j,j)+1,degree(j)
      do 308 ij=1,stib(tfnarg(0)+i)
      if(stib(stib(tfo(0)+i)+ij).eq.pmap(j,jj))then
      ii=ii+1
      elseif(stib(link(0)+stib(stib(tfo(0)+i)+ij)).eq.
     :pmap(j,jj))then
      ii=ii+1
      endif
  308 continue
  317 continue
  319 continue
      elseif(abs(stib(tftype(0)+i)).eq.3)then
      ii=0
      do 519 j=nleg+1,n
      do 515 jj=rdeg(j)+1,rdeg(j)+g(j,j),2
      if(flow(amap(j,jj),0).eq.1)then
      do 511 ij=1,stib(tfnarg(0)+i)
      if(stib(stib(tfo(0)+i)+ij).eq.pmap(j,jj))then
      ii=ii+1
      else
      k=stib(link(0)+stib(stib(tfo(0)+i)+ij))
      if(k.eq.pmap(j,jj))then
      ii=ii+1
      endif
      endif
  511 continue
      endif
  515 continue
      do 517 jj=rdeg(j)+g(j,j)+1,degree(j)
      if(flow(amap(j,jj),0).eq.1)then
      do 508 ij=1,stib(tfnarg(0)+i)
      if(stib(stib(tfo(0)+i)+ij).eq.pmap(j,jj))then
      ii=ii+1
      else
      k=stib(link(0)+stib(stib(tfo(0)+i)+ij))
      if(k.eq.pmap(j,jj))then
      ii=ii+1
      endif
      endif
  508 continue
      endif
  517 continue
  519 continue
      elseif(abs(stib(tftype(0)+i)).lt.6)then
      if(abs(stib(tftype(0)+i)).eq.2)then
      i1=2
      i2=3
      elseif(abs(stib(tftype(0)+i)).eq.4)then
      i1=2
      i2=2
      elseif(abs(stib(tftype(0)+i)).eq.5)then
      i1=3
      i2=3
      else
      auxlin(1:srec)='main_4'
      call messag(1,0,0,0)
      endif
      ii=0
      do 419 j=nleg+1,n
      do 415 jj=rdeg(j)+1,rdeg(j)+g(j,j),2
      ij=flow(amap(j,jj),0)
      if(ij.eq.i1.or.ij.eq.i2)then
      do 411 ij=1,stib(tfnarg(0)+i)
      if(stib(stib(tfo(0)+i)+ij).eq.pmap(j,jj))then
      ii=ii+1
      else
      k=stib(link(0)+stib(stib(tfo(0)+i)+ij))
      if(k.eq.pmap(j,jj))then
      ii=ii+1
      endif
      endif
  411 continue
      endif
  415 continue
      do 417 jj=rdeg(j)+g(j,j)+1,degree(j)
      ij=flow(amap(j,jj),0)
      if(ij.eq.i1.or.ij.eq.i2)then
      do 408 ij=1,stib(tfnarg(0)+i)
      if(stib(stib(tfo(0)+i)+ij).eq.pmap(j,jj))then
      ii=ii+1
      else
      k=stib(link(0)+stib(stib(tfo(0)+i)+ij))
      if(k.eq.pmap(j,jj))then
      ii=ii+1
      endif
      endif
  408 continue
      endif
  417 continue
  419 continue
      elseif(abs(stib(tftype(0)+i)).eq.6)then
      jj=stib(stib(tfo(0)+i)+1)
      ii=0
      i1=stib(vmkvpp(0)+jj)
      i2=stib(vmkvlp(0)+jj)
      do 191 k=rho(1)+1,n
      kk=stib(vfo(k))
      call stoz(stcb,stib(i1+kk),stib(i1+kk)-1+stib(i2+kk),ij)
      ii=ii+ij
  191 continue
      else
      jj=stib(stib(tfo(0)+i)+1)
      ii=0
      i1=stib(pmkvpp(0)+jj)
      i2=stib(pmkvlp(0)+jj)
      do 819 j=nleg+1,n
      do 815 k=rdeg(j)+1,rdeg(j)+g(j,j),2
      kk=pmap(j,k)
      call stoz(stcb,stib(i1+kk),stib(i1+kk)-1+stib(i2+kk),
     :ij)
      ii=ii+ij
  815 continue
      do 817 k=rdeg(j)+g(j,j)+1,degree(j)
      kk=pmap(j,k)
      call stoz(stcb,stib(i1+kk),stib(i1+kk)-1+stib(i2+kk),
     :ij)
      ii=ii+ij
  817 continue
  819 continue
      endif
      if(stib(tftype(0)+i).gt.0)then
      if(ii.lt.stib(tfa(0)+i).or.ii.gt.stib(tfb(0)+i))then
      goto 58
      endif
      elseif(stib(tftype(0)+i).lt.0)then
      if(ii.ge.stib(tfa(0)+i).and.ii.le.stib(tfb(0)+i))then
      goto 58
      endif
      endif
  314 continue
      i=ndiag+1
      if(ndiag.lt.i)then
      ndiag=i
      else
      auxlin(1:srec)='too many diagrams, arithmetic overflow'
      call messag(1,0,0,0)
      endif
      if(mflag(16).eq.0)then
      goto 58
      endif
      if(nloop.eq.0)then
      goto 231
      endif
      do 236 i=rho(1)+1,n
      j=rdeg(i)+1
  420 continue
      if(j.le.degree(i))then
      ii=vmap(i,j)
      aux=g(i,ii)
      k=j+aux
      if(i.ne.ii)then
  430 continue
      if(j.lt.k)then
      kk=1
      aux=pmap(i,j)
  440 continue
      if(j+kk.lt.k.and.aux.eq.pmap(i,j+kk))then
      kk=kk+1
      dsym=dsym*kk
      goto 440
      endif
      j=j+kk
      goto 430
      endif
      else
  450 continue
      if(j.lt.k)then
      kk=1
      aux=pmap(i,j)
  460 continue
      if(j+kk+kk.lt.k.and.aux.eq.pmap(i,j+kk+kk))then
      kk=kk+1
      dsym=dsym*kk
      goto 460
      endif
      if(aux.eq.stib(link(0)+aux))then
      do 226 ii=1,kk
      dsym=dsym+dsym
  226 continue
      endif
      j=j+kk+kk
      goto 450
      endif
      endif
      goto 420
      endif
  236 continue
  231 continue
      do 300 i=nleg+1,n
      do 250 j=1,degree(i)
      if(vmap(i,j).gt.nleg)then
      k=pmap(i,j)
      if(k.le.stib(link(0)+k))then
      if(k.eq.stib(link(0)+k))then
      if(i.gt.vmap(i,j))then
      goto 250
      endif
      if(i.eq.vmap(i,j).and.mod(j-rdeg(i),2).eq.0)then
      goto 250
      endif
      endif
      k=amap(i,j)-nleg
      ex(k)=i
      ey(k)=j
      endif
      endif
  250 continue
  300 continue
      do 130 i=nleg+1,n
      do 114 j=1,degree(i)
      xlist(j)=0
  114 continue
      do 125 j=1,degree(i)
      k=stib(stib(vparto(0)+stib(vfo(i)))+j)
      do 142 ik=1,degree(i)
      if(xlist(ik).eq.0.and.pmap(i,ik).eq.k)then
      xlist(ik)=1
      goto 148
      endif
  142 continue
      auxlin(1:srec)='main_2'
      call messag(1,0,0,0)
  148 continue
      ovm(i,j)=ik
  125 continue
  130 continue
      dis=1
      if(mflag(13).eq.0)then
      goto 495
      endif
      nf=0
      nf1=0
      do 256 i=rho(1)+1,n
      do 246 j=1,degree(i)
      k=ovm(i,j)
      ii=pmap(i,k)
      if(stib(antiq(0)+ii).eq.1)then
      nf=nf+1
      ij=vmap(i,k)
      if(ij.le.nleg)then
      ij=p1(ij)
      if(ij.le.incom)then
      jj=1-2*ij
      else
      jj=2*(incom-ij)
      endif
      elseif(ii.lt.stib(link(0)+ii))then
      jj=2*(amap(i,k)-nleg)-1
      elseif(ii.gt.stib(link(0)+ii))then
      jj=2*(amap(i,k)-nleg)
      elseif(i.lt.ij)then
      jj=2*(amap(i,k)-nleg)-1
      elseif(i.gt.ij)then
      jj=2*(amap(i,k)-nleg)
      elseif(mod(k-rdeg(i),2).eq.1)then
      jj=2*(amap(i,k)-nleg)-1
      else
      jj=2*(amap(i,k)-nleg)
      endif
      if(ij.lt.0)then
      nf1=nf1+1
      endif
      slist(nf)=jj
      endif
  246 continue
  256 continue
  266 continue
      ii=0
      do 276 i=1,nf
      if(slist(i).gt.ii)then
      ii=slist(i)
      endif
  276 continue
      if(ii.gt.0)then
      i=nf
      k1=0
      k2=0
  293 continue
      if(slist(i).gt.ii-2)then
      k2=k1
      k1=slist(i)
      if(i.ne.nf)then
      slist(i)=slist(nf)
      dis=-dis
      endif
      nf=nf-1
      endif
      if(i.gt.1.and.k2.eq.0)then
      i=i-1
      goto 293
      endif
      if(k1.gt.k2)then
      dis=-dis
      endif
      goto 266
      endif
      do 296 i=1,incom
      k=0
      do 291 j=nf,1,-1
      if(slist(j).eq.1-2*i)then
      if(j.ne.nf)then
      slist(j)=slist(nf)
      dis=-dis
      endif
      nf=nf-1
      goto 292
      endif
  291 continue
  292 continue
  296 continue
      do 306 i=rho(1),incom+1,-1
      k=0
      do 301 j=nf,1,-1
      if(slist(j).eq.2*(incom-i))then
      if(j.ne.nf)then
      slist(j)=slist(nf)
      dis=-dis
      endif
      nf=nf-1
      goto 302
      endif
  301 continue
  302 continue
  306 continue
      if(nf.ne.0)then
      auxlin(1:srec)='main_3'
      call messag(1,0,0,0)
      endif
  495 continue
      call compac
      goto 58
   94 continue
      if(mflag(16).ne.0)then
      call prepos(3)
      endif
      i0=stcbs(cmal(0))
      jj=stotl
      call vaocb(jj)
      stcb(i0+1:i0+stotl)=stcb(stotp:stotp-1+stotl)
      ii=i0+jj
      jj=jj+wztos(ndiag)
      call vaocb(jj)
      call karat(ndiag,ii,stcb,scbuff,0)
      if(mflag(8).eq.0)then
      ii=i0+jj
      kk=sdial
      if(ndiag.eq.1)then
      kk=kk-1
      endif
      jj=jj+kk
      call vaocb(jj)
      stcb(ii+1:ii+kk)=stcb(sdiap:sdiap+kk)
      else
      ii=i0+jj
      kk=sgral
      if(ndiag.eq.1)then
      kk=kk-1
      endif
      jj=jj+kk
      call vaocb(jj)
      stcb(ii+1:ii+kk)=stcb(sgrap:sgrap+kk)
      endif
      ii=i0+jj
      print *,' '
      print *,stcb(i0+1:ii)
      print *,' '
      stop
      end
      subroutine prepos(what)
      implicit integer(a-z)
      save
      parameter ( maxleg=10, maxrho=6, maxi=10 )
      parameter ( srec=81 )
      parameter ( lfeed=10, dquote=34, rquote=39 )
      parameter ( eoia=-63 )
      parameter ( sibuff=1048576, scbuff=524288, mimal=2, mcmal=2 )
      common/z2g/dsym,dis,ndiag
      common/z6in/dunit,munit,ounit,sunit,funit
      common/z9in/lofile,ofilea,ofileb
      common/z24g/iogp(1:4)
      character*(srec) auxlin
      common/z22g/auxlin
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z35g/stibs(0:mimal),stcbs(0:mcmal),imal(0:0),cmal(0:0)
      common/z36g/udkp(0:2),udkl(0:0),udkt(0:0),udki(0:0)
      common/z37g/drecp(0:0),drecl(0:0),drecii(0:0),irecc(0:0),
     :frecc(0:0),ndrec,ncom
      common/z38g/gmkp(0:0),gmkl(0:0),gmkd(0:0),gmko(0:0),
     :gmkvp(0:0),gmkvl(0:0)
      common/z39g/pmkp(0:0),pmkl(0:0),pmkd(0:0),pmkvpp(0:0),pmkvlp(0:0)
      common/z40g/vmkp(0:0),vmkl(0:0),vmkvpp(0:0),vmkvlp(0:0),vmks(0:0)
      common/z41g/nudk,ngmk,npmk,nvmk,tgmkd,tpmkd
      common/z45g/qdatp,qdatl,qvp,qvl,dprefp,dprefl
      logical lunit
      integer prea
      data prea/0/
      if(what.ne.1.and.what.ne.3)then
      auxlin(1:srec)='internal error in subroutine prepos'
      call messag(1,0,0,0)
      endif
      if(prea.eq.0)then
      ios=1
      lunit=.true.
      inquire(file=stcb(ofilea:ofileb),exist=lunit,iostat=ios)
      if(ios.ne.0)then
      auxlin(1:srec)='system/filesystem error'
      call messag(1,0,0,0)
      elseif(lunit)then
      auxlin(1:srec)='output file already exists'
      call messag(1,0,0,0)
      endif
      ios=1
      ounit=0
      lunit=.true.
   10 continue
      ounit=ounit+1
      inquire(unit=ounit,opened=lunit,iostat=ios)
      if(ios.ne.0)then
      ounit=0
      auxlin(1:srec)='system/filesystem error'
      call messag(1,0,0,0)
      elseif(lunit)then
      if(ounit.ge.99)then
      ounit=0
      auxlin(1:srec)='no logical unit number available'
      call messag(1,0,0,0)
      endif
      goto 10
      endif
      ios=1
      open(unit=ounit,file=stcb(ofilea:ofileb),status='new',
     :access='sequential',iostat=ios)
      if(ios.ne.0)then
      ounit=0
      auxlin(1:srec)='output file could not be opened'
      call messag(1,0,0,0)
      endif
      prea=1
      endif
      ki=iogp(what)
      klin=0
      lup=0
   20 continue
      if(ki.lt.iogp(what).or.ki.ge.iogp(what+1))then
      goto 80
      endif
      if(stib(ki).gt.0)then
      klin=klin+1
      call vaocb(klin)
      i=stcbs(1)+klin
      stcb(i:i)=char(stib(ki))
      ki=ki+1
      goto 20
      endif
      xstep=1
      if(stib(ki).eq.-1)then
      if(stib(ki+1).eq.11)then
      if(stib(ki+4).eq.0)then
      stib(ki+4)=ncom
      stib(ki+5)=0
      endif
      stib(ki+5)=stib(ki+5)+1
      if(stib(ki+4).lt.stib(ki+5))then
      stib(ki+4)=0
      stib(ki+5)=0
      xstep=0
      icom=0
      else
      icom=stib(ki+5)
      lup=11
      endif
      elseif(stib(ki+1).eq.12)then
      if(stib(ki+4).eq.0)then
      stib(ki+4)=stib(frecc(0)+icom)
      stib(ki+5)=stib(irecc(0)+icom)-1
      endif
      stib(ki+5)=stib(ki+5)+1
      if(stib(ki+4).lt.stib(ki+5))then
      stib(ki+4)=0
      stib(ki+5)=0
      xstep=0
      iline=0
      else
      iline=stib(ki+5)
      lup=12
      endif
      elseif(stib(ki+1).eq.19)then
      xstep=0
      lup=0
      else
      goto 80
      endif
      elseif(stib(ki).eq.-2)then
      if(stib(ki+1).eq.21)then
      if(klin.lt.1)then
      goto 80
      endif
      klin=klin-1
      elseif(stib(ki+1).eq.22)then
      klin=klin+1
      call vaocb(klin)
      i=stcbs(1)+klin
      stcb(i:i)=char(lfeed)
      else
      goto 80
      endif
      elseif(stib(ki).eq.-3)then
      if(stib(ki+1).eq.71)then
      jj=ndiag
      i=stcbs(1)+klin
      klin=klin+wztos(jj)
      call vaocb(klin)
      call karat(jj,i,stcb,scbuff,1)
      elseif(stib(ki+1).eq.81)then
      klin=klin+qvl
      call vaocb(klin)
      i=stcbs(1)+klin
      stcb(i-qvl+1:i)=stcb(qvp:qvp-1+qvl)
      else
      goto 80
      endif
      elseif(stib(ki).eq.-4)then
      if(stib(ki+1).eq.82)then
      if(lup.eq.11)then
      do 70 i1=stib(irecc(0)+icom),stib(frecc(0)+icom)
      j=stib(drecp(0)+i1)
      jj=stib(drecl(0)+i1)
      klin=klin+jj
      call vaocb(klin)
      i=stcbs(1)+klin
      stcb(i-jj+1:i)=stcb(j:j-1+jj)
   70 continue
      elseif(lup.eq.12)then
      jj=stib(drecl(0)+iline)
      klin=klin+jj
      call vaocb(klin)
      j=stib(drecp(0)+iline)
      i=stcbs(1)+klin
      stcb(i-jj+1:i)=stcb(j:j-1+jj)
      else
      goto 80
      endif
      else
      goto 80
      endif
      elseif(stib(ki).eq.-6)then
      ii=stib(ki+1)
      if(ii.lt.1.or.ii.gt.nudk)then
      goto 80
      endif
      if(stib(udkt(0)+ii).eq.1)then
      ij=stib(udki(0)+ii)
      if(ij.lt.1.or.ij.gt.ngmk)then
      goto 80
      elseif(stib(gmkd(0)+ij).ne.1)then
      goto 80
      endif
      ij=stib(gmko(0)+ij)+1
      j=stib(gmkvp(0)+ij)
      jj=stib(gmkvl(0)+ij)
      if(jj.gt.0)then
      klin=klin+jj
      call vaocb(klin)
      i=stcbs(1)+klin
      stcb(i-jj+1:i)=stcb(j:j-1+jj)
      endif
      if(jj.lt.0)then
      goto 80
      endif
      else
      goto 80
      endif
      elseif(stib(ki).eq.eoia)then
      goto 30
      else
      goto 80
      endif
      if(xstep.eq.1)then
      ki=ki+stib(ki+2)
      else
      ki=stib(ki+3)
      endif
      goto 20
   30 continue
      if(klin.gt.0)then
      i=stcbs(1)+klin
      if(stcb(i:i).ne.char(lfeed))then
      goto 80
      endif
      endif
      j=stcbs(1)+1
      do 29 i=stcbs(1)+1,stcbs(1)+klin
      if(stcb(i:i).eq.char(lfeed))then
      ios=1
      if(j.lt.i)then
      write(unit=ounit,fmt=404,iostat=ios)stcb(j:i-1)
      else
      write(unit=ounit,fmt=404,iostat=ios)
      endif
      if(ios.ne.0)then
      auxlin(1:srec)='run-time error while writing to file'
      call messag(1,0,0,0)
      endif
      j=i+1
      endif
   29 continue
  404 format(a)
      if(what.eq.3)then
      ios=1
      close(unit=ounit,status='keep',iostat=ios)
      if(ios.ne.0)then
      auxlin(1:srec)='system/filesystem error'
      call messag(1,0,0,0)
      endif
      ounit=0
      endif
      return
   80 continue
      if(what.eq.1)then
      auxlin(1:srec)='run-time error while processing <prologue>'
      else
      auxlin(1:srec)='run-time error while processing <epilogue>'
      endif
      call messag(1,0,0,0)
      end
      subroutine umpi(xx,situ)
      implicit integer(a-z)
      save
      parameter ( maxleg=10, maxrho=6, maxi=10 )
      parameter ( maxdeg=6 )
      parameter ( maxn=maxi+maxi-2, maxli=maxn+maxrho-1 )
      parameter ( srec=81 )
      parameter ( lfeed=10, dquote=34, rquote=39 )
      common/z1g/nrho,rho(1:maxdeg),g(1:maxn,1:maxn)
      common/z4g/n,nli
      common/z8g/degree(1:maxn)
      common/z17g/xtail(1:maxn),xhead(1:maxn),ntadp
      character*(srec) auxlin
      common/z22g/auxlin
      integer stack(1:maxn),lista(1:maxn)
      ntadp=0
      situ=-1
      if(xx.eq.5)then
      do 55 ii=rho(1)+1,n
      if(g(ii,ii).ne.0)then
      return
      endif
   55 continue
      endif
      do 66 ii=rho(1)+1,n
      ij=ii+1
      if(xx.eq.5)then
      if(degree(ii).lt.3)then
      goto 66
      endif
      ij=n
      endif
      do 56 jj=ij,n
      if(g(ii,jj).eq.1.or.xx.eq.5)then
      if(xx.ne.5)then
      g(ii,jj)=0
      k=1
      stack(1)=1
      lista(1)=1
      do 06 i=2,n
      lista(i)=0
   06 continue
      else
      k=rho(1)
      do 02 i=1,rho(1)
      lista(i)=1
      stack(i)=i
   02 continue
      do 04 i=rho(1)+1,n
      lista(i)=0
   04 continue
      endif
      do 46 j=1,n
      if(k.eq.n)then
      goto 20
      endif
      if(j.gt.k)then
      if(xx.ne.5)then
      g(ii,jj)=1
      endif
      if(xx.eq.1.or.xx.eq.5)then
      return
      endif
      aux=0
      do 16 kk=1,rho(1)
      aux=aux+lista(kk)
   16 continue
      if(xx.eq.2)then
      if(aux.eq.0.or.aux.eq.rho(1))then
      return
      endif
      elseif(xx.eq.3)then
      if(aux.eq.0.or.aux.eq.rho(1))then
      ntadp=ntadp+1
      xtail(ntadp)=ii
      xhead(ntadp)=jj
      endif
      elseif(xx.eq.4)then
      if(aux.eq.1.or.aux.eq.rho(1)-1)then
      return
      endif
      endif
      goto 20
      endif
      aux=stack(j)
      if(xx.ne.5.or.aux.ne.ii)then
      do 26 i=1,aux-1
      if(lista(i).eq.0)then
      if(g(i,aux).gt.0)then
      k=k+1
      stack(k)=i
      lista(i)=1
      endif
      endif
   26 continue
      do 36 i=aux+1,n
      if(lista(i).eq.0)then
      if(g(aux,i).gt.0)then
      k=k+1
      stack(k)=i
      lista(i)=1
      endif
      endif
   36 continue
      endif
   46 continue
      auxlin(1:srec)='umpi_1'
      call messag(1,0,0,0)
   20 continue
      if(xx.ne.5)then
      g(ii,jj)=1
      endif
      endif
   56 continue
   66 continue
      situ=1
      end
      subroutine karat(m,ic,c,lc,so)
      implicit integer(a-z)
      save
      parameter ( srec=81 )
      parameter ( lfeed=10, dquote=34, rquote=39 )
      character*(*) c
      character*(srec) auxlin
      character*1 c1
      common/z22g/auxlin
      integer m,ic,xm,jc,ii,jj,lc
      jc=ic
      if(m.lt.0)then
      xm=-m
      if(jc.ge.lc)then
      goto 30
      endif
      jc=jc+1
      c(jc:jc)='-'
      else
      xm=m
      endif
      ii=jc+1
   10 continue
      if(jc.ge.lc)then
      goto 30
      endif
      jc=jc+1
      if(xm.gt.9)then
      ym=xm/10
      c(jc:jc)=char((xm-10*ym)+48)
      xm=ym
      goto 10
      else
      c(jc:jc)=char(xm+48)
      endif
      jj=jc
   20 continue
      if(ii.lt.jj)then
      c1(1:1)=c(ii:ii)
      c(ii:ii)=c(jj:jj)
      c(jj:jj)=c1(1:1)
      ii=ii+1
      jj=jj-1
      goto 20
      endif
      return
   30 continue
      if(so.eq.0)then
      c(lc-2:lc)='...'
      endif
      auxlin(1:srec)='string too short in subroutine karat'
      call messag(so,0,0,0)
      end
      subroutine stoz(s,ia,ib,n)
      implicit integer(a-z)
      save
      parameter ( srec=81 )
      parameter ( zwidth=4 )
      character*(*) s
      character*(srec) auxlin
      common/z22g/auxlin
      if(ia.lt.1.or.ia.gt.ib)then
      goto 30
      endif
      if(s(ia:ia).eq.'+'.or.s(ia:ia).eq.'-')then
      k=1
      else
      k=0
      endif
      if(ia-1+k.lt.ib-zwidth)then
      auxlin(1:srec)='integer too large (in absolute value)'
      call messag(1,0,0,0)
      endif
      n=0
      i=ia+k
   10 continue
      j=ichar(s(i:i))-48
      if(j.lt.0.or.j.gt.9)then
      goto 30
      endif
      n=10*n+j
      if(i.lt.ib)then
      i=i+1
      goto 10
      endif
      if(s(ia:ia).eq.'-')then
      n=-n
      endif
      return
   30 continue
      auxlin(1:srec)='wrong input for subroutine stoz'
      call messag(1,0,0,0)
      end
      integer function wztos(n)
      implicit integer(a-z)
      save
      m=n
      if(m.lt.0)then
      wztos=2
   10 continue
      if(m.lt.-9)then
      m=m/10
      wztos=wztos+1
      goto 10
      endif
      else
      wztos=1
   20 continue
      if(m.gt.9)then
      m=m/10
      wztos=wztos+1
      goto 20
      endif
      endif
      end
      subroutine spack(s,ind,m,uc,nos)
      implicit integer(a-z)
      save
      parameter ( srec=81 )
      parameter ( sibuff=1048576, scbuff=524288, mimal=2, mcmal=2 )
      parameter ( lfeed=10, dquote=34, rquote=39 )
      parameter ( esse=115 )
      character*(srec) auxlin
      common/z22g/auxlin
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z35g/stibs(0:mimal),stcbs(0:mcmal),imal(0:0),cmal(0:0)
      character*(*) s
      if(cmal(0).ne.1.or.imal(0).ne.1)then
      goto 08
      elseif(m.ne.0.and.m.ne.1)then
      goto 08
      elseif(nos.ne.0.and.nos.ne.1)then
      goto 08
      elseif(uc.ne.0.and.uc.ne.1)then
      goto 08
      endif
      inos=nos
      iuc=uc
      i1=0
      sl=1
   03 continue
      if(s(sl:sl).ne.';')then
      if(sl.lt.srec)then
      if(i1.eq.0.and.s(sl:sl).eq.',')then
      i1=sl
      endif
      sl=sl+1
      goto 03
      else
      goto 08
      endif
      endif
      if(i1.eq.0)then
      i1=sl
      endif
      if(i1.lt.3.or.sl.lt.3)then
      goto 08
      elseif(mod(i1,2).eq.0)then
      goto 08
      endif
      sp=stcbs(cmal(0))+1
      i1=(i1-1)/2
      call aocb(i1+1)
      ij=0
      j=1
      k=sp
      do 04 i=1,i1
      ii=ichar(s(j:j))
      jj=ichar(s(j+1:j+1))
      if(ii.lt.48.or.ii.gt.57)then
      goto 08
      elseif(jj.lt.48.or.jj.gt.57)then
      goto 08
      endif
      kk=10*ii+jj-496
      if(iuc.eq.1.and.(kk.gt.64.and.kk.lt.91))then
      goto 08
      endif
      if(kk.gt.96.and.kk.lt.123)then
      ij=1
      endif
      stcb(k:k)=char(kk)
      j=j+2
      k=k+1
   04 continue
      stcb(k:k)=char(lfeed)
      if(uic.eq.1.and.ij.eq.0)then
      goto 08
      endif
      if(inos.eq.1)then
      if(stcb(k-1:k-1).ne.char(esse))then
      inos=0
      endif
      endif
      j0=stibs(imal(0))
      if(m.eq.0)then
      call aoib(2)
      elseif(m.eq.1)then
      call vaoib(2)
      endif
      stib(j0+1)=sp
      stib(j0+2)=i1
      if(m.eq.1)then
      if(nos.ne.0.or.uc.ne.0)then
      goto 08
      endif
      return
      endif
      j2=i1+i1+1
      j1=j2+1
   11 continue
      if(j2.lt.sl)then
      j2=j2+1
      if(s(j2:j2).eq.','.or.s(j2:j2).eq.';')then
      j2=j2-1
      if(j2.lt.j1)then
      goto 08
      endif
      call stoz(s,j1,j2,i3)
      call aoib(1)
      stib(stibs(imal(0)))=i3
      j2=j2+1
      j1=j2+1
      endif
      goto 11
      endif
      ind=ind+1
      if(iuc.eq.1)then
      ii=stcbs(cmal(0))
      call aocb(i1+1)
      i3=ii+1
      i4=ii-i1
      sp=i3
      do 07 i=1,i1
      j=ichar(stcb(i4:i4))
      if(j.gt.96.and.j.lt.123)then
      j=j-32
      endif
      stcb(i3:i3)=char(j)
      i3=i3+1
      i4=i4+1
   07 continue
      stcb(i3:i3)=stcb(i4:i4)
      j1=stibs(imal(0))
      jj=j1-j0
      if(jj.lt.2)then
      goto 08
      endif
      call aoib(jj)
      ii=j1+1
      stib(ii)=stib(ii-jj)+i1+1
      do 17 i=2,jj
      ii=ii+1
      stib(ii)=stib(ii-jj)
   17 continue
      ind=ind+1
      endif
      if(inos.eq.1)then
      j1=stibs(imal(0))
      jj=j1-j0
      if(jj.lt.4)then
      goto 08
      endif
      call aoib(jj)
      ii=j1
      do 19 i=1,jj
      ii=ii+1
      stib(ii)=stib(ii-jj)
   19 continue
      ii=j1+2
      stib(ii)=stib(ii)-1
      ind=ind+1
      if(uc.eq.1)then
      ii=ii+jj/2
      stib(ii)=stib(ii)-1
      ind=ind+1
      endif
      endif
      return
   08 continue
      auxlin(1:srec)='invalid input for subroutine spack'
      call messag(1,0,0,0)
      end
      integer function stds(s,ia,ib,iw,ilf)
      implicit integer(a-z)
      save
      parameter ( sibuff=1048576, scbuff=524288, mimal=2, mcmal=2 )
      parameter ( srec=81 )
      parameter ( lfeed=10, dquote=34, rquote=39 )
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z32g/acf0(0:127),acf1(0:127),sucpal(0:11,0:11)
      common/z35g/stibs(0:mimal),stcbs(0:mcmal),imal(0:0),cmal(0:0)
      character*(*) s
      sc=0
      sl=0
      sf=0
      itl=1
      stds=-1
      if(ia.lt.1.or.ia.ge.ib)then
      goto 30
      elseif(iw.ne.0.and.iw.ne.1.and.iw.ne.2)then
      goto 30
      elseif(ilf.ne.0.and.ilf.ne.1)then
      goto 30
      endif
      ii=0
      ll=0
      mm=0
      do 10 i=ia,ib
      if(s(i:i).eq.char(rquote))then
      ii=ii+1
      mm=1-mm
      ll=0
      if(ii.gt.1)then
      kk=mm
      else
      kk=0
      endif
      elseif(s(i:i).eq.char(lfeed))then
      sf=sf+1
      if(mm.ne.0.or.ll.eq.1)then
      goto 30
      endif
      ll=1
      if(ii.gt.0)then
      ii=0
      sc=sc+1
      endif
      kk=0
      elseif(s(i:i).eq.' ')then
      if(ii.gt.0)then
      ii=0
      if(mm.eq.0)then
      sc=sc+1
      endif
      endif
      kk=mm
      else
      if(mm.ne.1.or.acf0(ichar(s(i:i))).lt.0)then
      goto 30
      endif
      ii=0
      kk=1
      endif
      if(kk.ne.0)then
      if(iw.eq.2)then
      if(s(i:i).ne.' ')then
      itl=0
      endif
      endif
      if(iw.eq.1.or.(iw.eq.2.and.itl.eq.0))then
      sl=sl+1
      call vaocb(sl)
      j=stcbs(1)+sl
      stcb(j:j)=s(i:i)
      endif
      endif
   10 continue
      if(mm.ne.0)then
      goto 30
      endif
      if(s(ib:ib).eq.char(rquote))then
      sc=sc+1
      endif
      if(iw.eq.2.and.sl.gt.1)then
      do 15 i=stcbs(1)+sl,stcbs(1)+1,-1
      if(stcb(i:i).eq.' ')then
      sl=sl-1
      else
      goto 20
      endif
   15 continue
      endif
   20 continue
      stds=sl
   30 continue
      return
      end
      subroutine mstr0(s,ia,ib,pp,pl,ind)
      implicit integer(a-z)
      save
      parameter ( eoia=-63 )
      parameter ( sibuff=1048576, scbuff=524288, mimal=2, mcmal=2 )
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      character*(*) s
      ind=0
      ls=ib-ia+1
      i=0
   10 continue
      i=i+1
      ii=stib(pl+i)
      if(ii.ne.eoia)then
      if(ii.eq.ls)then
      jj=stib(pp+i)
      if(stcb(jj:jj-1+ls).eq.s(ia:ib))then
      ind=i
      return
      endif
      endif
      goto 10
      endif
      end
      subroutine mstr1(s,ia,ib,pp,kp,kl,ind)
      implicit integer(a-z)
      save
      parameter ( eoia=-63 )
      parameter ( sibuff=1048576, scbuff=524288, mimal=2, mcmal=2 )
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      character*(*) s
      ind=0
      lsm1=ib-ia
      ls=lsm1+1
      p0=pp
   10 continue
      p0=stib(p0)
      if(p0.eq.eoia)then
      ind=0
      return
      endif
      ind=ind+1
      if(stib(p0+kl).eq.ls)then
      jj=stib(p0+kp)
      if(stcb(jj:jj+lsm1).eq.s(ia:ib))then
      return
      endif
      endif
      goto 10
      end
      integer function stdw(s,ia,ib)
      implicit integer(a-z)
      save
      parameter ( srec=81 )
      character*(srec) auxlin
      common/z22g/auxlin
      common/z32g/acf0(0:127),acf1(0:127),sucpal(0:11,0:11)
      character*(*) s
      if(ia.gt.ib.or.ia.lt.1)then
      auxlin(1:srec)='input for function stdw is invalid'
      call messag(1,0,0,0)
      endif
      if(acf1(ichar(s(ia:ia))).lt.2)then
      stdw=0
      else
      stdw=1
      endif
      if(stdw.eq.1)then
      do 20 i=ia+1,ib
      if(acf1(ichar(s(i:i))).lt.0)then
      stdw=0
      goto 30
      endif
   20 continue
      endif
   30 continue
      end
      integer function stdun(s,ia,ib)
      implicit integer(a-z)
      save
      parameter ( srec=81 )
      character*(srec) auxlin
      common/z22g/auxlin
      common/z32g/acf0(0:127),acf1(0:127),sucpal(0:11,0:11)
      character*(*) s
      if(ia.gt.ib.or.ia.lt.1)then
      auxlin(1:srec)='input for function stdun is invalid'
      call messag(1,0,0,0)
      endif
      do 20 i=ia,ib
      if(acf1(ichar(s(i:i))).ne.1)then
      stdun=0
      goto 30
      endif
   20 continue
      stdun=1
   30 continue
      end
      integer function stdz(s,ia,ib)
      implicit integer(a-z)
      save
      parameter ( srec=81 )
      character*(srec) auxlin
      common/z22g/auxlin
      common/z32g/acf0(0:127),acf1(0:127),sucpal(0:11,0:11)
      character*(*) s
      if(ia.gt.ib.or.ia.lt.1)then
      auxlin(1:srec)='input for function stdz is invalid'
      call messag(1,0,0,0)
      endif
      stdz=0
      if(s(ia:ia).eq.'+'.or.s(ia:ia).eq.'-')then
      if(ia.lt.ib)then
      stdz=1
      endif
      elseif(acf1(ichar(s(ia:ia))).eq.1)then
      stdz=1
      endif
      if(stdz.eq.1)then
      do 20 i=ia+1,ib
      if(acf1(ichar(s(i:i))).ne.1)then
      stdz=0
      goto 30
      endif
   20 continue
      endif
   30 continue
      end
      integer function stdq(s,ia,ib)
      implicit integer(a-z)
      save
      parameter ( srec=81 )
      character*(srec) auxlin
      common/z22g/auxlin
      common/z32g/acf0(0:127),acf1(0:127),sucpal(0:11,0:11)
      character*(*) s
      xsla=0
      xnum=0
      xden=0
      if(ia.gt.ib.or.ia.lt.1)then
      auxlin(1:srec)='input for function stdq is invalid'
      call messag(1,0,0,0)
      endif
      if(s(ia:ia).eq.'+')then
      stdq=1
      elseif(s(ia:ia).eq.'-')then
      stdq=1
      elseif(acf1(ichar(s(ia:ia))).eq.1)then
      stdq=1
      xnum=1
      else
      stdq=0
      endif
      if(stdq.eq.1)then
      do 10 i=ia+1,ib
      if(acf1(ichar(s(i:i))).eq.1)then
      if(xsla.eq.0)then
      xnum=1
      elseif(s(i:i).ne.'0')then
      xden=1
      endif
      elseif(s(i:i).eq.'/')then
      if(xsla.ne.0)then
      stdq=0
      goto 20
      endif
      xsla=1
      else
      stdq=0
      goto 20
      endif
   10 continue
      if(xnum.ne.1)then
      stdq=0
      elseif(xsla.ne.0)then
      if(xden.ne.1)then
      stdq=0
      endif
      endif
      endif
   20 continue
      end
      subroutine compac
      implicit integer(a-z)
      save
      parameter ( maxleg=10, maxrho=6, maxi=10 )
      parameter ( maxdeg=6 )
      parameter ( maxn=maxi+maxi-2, maxli=maxn+maxrho-1 )
      parameter ( srec=81 )
      parameter ( lfeed=10, dquote=34, rquote=39 )
      parameter ( eoia=-63 )
      parameter ( sibuff=1048576, scbuff=524288, mimal=2, mcmal=2 )
      common/z1in/leg(1:maxleg),nleg,nloop,incom
      common/z3in/momep(0:maxleg),momel(0:maxleg)
      common/z6in/dunit,munit,ounit,sunit,funit
      common/z2g/dsym,dis,ndiag
      common/z3g/vparto(0:0),vdeg(0:0),nvert
      common/z4g/n,nli
      common/z6g/p1(1:maxleg),invp1(1:maxleg)
      common/z7g/lmap(1:maxn,1:maxdeg),vmap(1:maxn,1:maxdeg),
     :pmap(1:maxn,1:maxdeg),vlis(1:maxn),invlis(1:maxn)
      common/z8g/degree(1:maxn)
      common/z9g/tpc(0:0)
      common/z11g/npart,nblok,nprop
      common/z16g/rdeg(1:maxn),amap(1:maxn,1:maxdeg)
      common/z18g/eg(1:maxn,1:maxn),flow(1:maxli,0:maxleg+maxrho)
      common/z19g/vfo(1:maxn)
      character*(srec) auxlin
      common/z22g/auxlin
      common/z24g/iogp(1:4)
      common/z25g/ex(1:maxli),ey(1:maxli),ovm(1:maxn,1:maxdeg)
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z33g/namep(0:0),namel(0:0)
      common/z34g/link(0:0),antiq(0:0),blok(0:0)
      common/z35g/stibs(0:mimal),stcbs(0:mcmal),imal(0:0),cmal(0:0)
      common/z36g/udkp(0:2),udkl(0:0),udkt(0:0),udki(0:0)
      common/z38g/gmkp(0:0),gmkl(0:0),gmkd(0:0),gmko(0:0),
     :gmkvp(0:0),gmkvl(0:0)
      common/z39g/pmkp(0:0),pmkl(0:0),pmkd(0:0),pmkvpp(0:0),pmkvlp(0:0)
      common/z40g/vmkp(0:0),vmkl(0:0),vmkvpp(0:0),vmkvlp(0:0),vmks(0:0)
      common/z41g/nudk,ngmk,npmk,nvmk,tgmkd,tpmkd
      ig=iogp(2)
      klin=0
      lupt=0
   10 continue
      if(ig.lt.iogp(2).or.ig.ge.iogp(3))then
      goto 80
      endif
      k=stib(ig)
      if(k.gt.0)then
      klin=klin+1
      call vaocb(klin)
      ip=stcbs(1)+klin
      stcb(ip:ip)=char(k)
      ig=ig+1
      goto 10
      endif
      xstep=1
      ii=stib(ig+1)
      if(k.gt.-3)then
      if(k.eq.-1)then
      lupt=lupty(ii)
      if(lupt.gt.0)then
      if(stib(ig+4).eq.0)then
      if(lupt.eq.3)then
      if(ii.eq.13)then
      stib(ig+4)=incom
      stib(ig+5)=0
      elseif(ii.eq.14)then
      stib(ig+4)=nleg
      stib(ig+5)=incom
      endif
      else
      if(lupt.eq.4)then
      stib(ig+4)=nli-nleg
      elseif(lupt.eq.5)then
      stib(ig+4)=n-nleg
      elseif(lupt.eq.6)then
      stib(ig+4)=degree(nleg+lupi)
      endif
      stib(ig+5)=0
      endif
      endif
      stib(ig+5)=stib(ig+5)+1
      if(stib(ig+4).lt.stib(ig+5))then
      stib(ig+4)=0
      stib(ig+5)=0
      xstep=0
      lupt=0
      else
      if(lupt.lt.6)then
      lupi=stib(ig+5)
      else
      lupj=stib(ig+5)
      endif
      endif
      elseif(lupt.lt.0)then
      xstep=0
      lupt=0
      else
      goto 80
      endif
      elseif(k.eq.-2)then
      if(ii.eq.21)then
      if(klin.gt.0)then
      klin=klin-1
      else
      goto 80
      endif
      elseif(ii.eq.22)then
      klin=klin+1
      call vaocb(klin)
      ip=stcbs(1)+klin
      stcb(ip:ip)=char(lfeed)
      else
      goto 80
      endif
      else
      goto 80
      endif
      elseif(k.gt.-5)then
      if(k.eq.-3)then
      if(ii.gt.70.and.ii.lt.80)then
      if(ii.lt.74)then
      if(ii.eq.71)then
      jj=ndiag
      elseif(ii.eq.72)then
      jj=dsym
      if(dsym.gt.1)then
      klin=klin+2
      call vaocb(klin)
      ip=stcbs(1)+klin
      stcb(ip-1:ip)='1/'
      endif
      elseif(ii.eq.73)then
      jj=dsym
      endif
      elseif(ii.lt.77)then
      if(ii.eq.74)then
      jj=nli-nleg
      elseif(ii.eq.75)then
      jj=nleg
      elseif(ii.eq.76)then
      jj=nloop
      endif
      else
      if(ii.eq.77)then
      jj=n-nleg
      elseif(ii.eq.78)then
      jj=incom
      elseif(ii.eq.79)then
      jj=nleg-incom
      endif
      endif
      ip=stcbs(1)+klin
      klin=klin+wztos(jj)
      call vaocb(klin)
      call karat(jj,ip,stcb,scbuff,1)
      else
      if(ii.eq.61)then
      klin=klin+1
      call vaocb(klin)
      ip=stcbs(1)+klin
      stcb(ip:ip)=char(44-dis)
      elseif(ii.eq.62)then
      if(dis.lt.0)then
      klin=klin+1
      call vaocb(klin)
      ip=stcbs(1)+klin
      stcb(ip:ip)='-'
      endif
      else
      goto 80
      endif
      endif
      elseif(k.eq.-4)then
      if(lupt.eq.4)then
      if(ii.lt.40)then
      if(ii.eq.31.or.ii.eq.33)then
      i=lupi
      j=pmap(ex(i),ey(i))
      if(ii.eq.33)then
      j=stib(link(0)+j)
      endif
      il=stib(namel(0)+j)
      klin=klin+il
      call vaocb(klin)
      ia=stib(namep(0)+j)
      ip=stcbs(1)+klin
      stcb(ip-il+1:ip)=stcb(ia:ia-1+il)
      elseif(ii.eq.32.or.ii.eq.34)then
      i=lupi
      if(ii.eq.32)then
      ik=1
      else
      ik=2
      endif
      j=ey(i)
      i=ex(i)
      k0=klin
      other=0
      if(vmap(i,j).eq.i)then
      if(mod(j-rdeg(i)+ik,2).ne.0)then
      klin=klin+1
      call vaocb(klin)
      ip=stcbs(1)+klin
      stcb(ip:ip)='-'
      endif
      klin=klin+momel(0)
      call vaocb(klin)
      ip=stcbs(1)+klin
      stcb(ip-momel(0)+1:ip)=
     :stcb(momep(0):momep(0)-1+momel(0))
      ij=eg(i,i)+(j-1-rdeg(i))/2
      do 372 k=nleg+1,nleg+nloop
      if(flow(ij,k).ne.0)then
      ip=stcbs(1)+klin
      klin=klin+wztos(k-nleg)
      call vaocb(klin)
      call karat(k-nleg,ip,stcb,scbuff,1)
      goto 70
      endif
  372 continue
      else
      if(vmap(i,j).lt.i)then
      ij=-1
      else
      ij=1
      endif
      if(ik.eq.2)then
      ij=-ij
      endif
      do 371 k=nleg+1,nleg+nloop
      jj=flow(amap(i,j),k)
      if(jj.ne.0)then
      if(other.eq.1.or.ij.eq.jj)then
      klin=klin+1
      call vaocb(klin)
      ip=stcbs(1)+klin
      stcb(ip:ip)=char(44+ij*jj)
      endif
      klin=klin+momel(0)
      call vaocb(klin)
      ip=stcbs(1)+klin
      stcb(ip-momel(0)+1:ip)=
     :stcb(momep(0):momep(0)-1+momel(0))
      klin=klin+wztos(k-nleg)
      call vaocb(klin)
      call karat(k-nleg,ip,stcb,scbuff,1)
      other=1
      endif
  371 continue
      do 363 k=1,nleg
      jj=flow(amap(i,j),invp1(k))
      if(jj.ne.0)then
      if(k.gt.incom)then
      jj=-jj
      endif
      if(other.eq.1.or.ij.eq.jj)then
      klin=klin+1
      call vaocb(klin)
      ip=stcbs(1)+klin
      stcb(ip:ip)=char(44+ij*jj)
      endif
      klin=klin+momel(k)
      call vaocb(klin)
      ip=stcbs(1)+klin
      stcb(ip-momel(k)+1:ip)=
     :stcb(momep(k):momep(k)-1+momel(k))
      other=1
      endif
  363 continue
      if(k0.eq.klin)then
      klin=klin+1
      call vaocb(klin)
      ip=stcbs(1)+klin
      stcb(ip:ip)='0'
      endif
      endif
      elseif(ii.eq.35)then
      i=lupi
      j=pmap(ex(i),ey(i))
      klin=klin+1
      call vaocb(klin)
      ia=stib(antiq(0)+j)
      ip=stcbs(1)+klin
      stcb(ip:ip)=char(43+ia+ia)
      else
      goto 80
      endif
      elseif(ii.lt.50)then
      if(ii.lt.44)then
      if(ii.eq.40)then
      jj=3
      elseif(ii.eq.41)then
      jj=lupi
      elseif(ii.eq.42)then
      jj=2*lupi-1
      elseif(ii.eq.43)then
      i=lupi
      j=ey(i)
      i=ex(i)
      jj=0
      do 740 k=1,degree(i)
      if(ovm(i,k).eq.j)then
      jj=k
      endif
  740 continue
      else
      goto 80
      endif
      elseif(ii.lt.47)then
      if(ii.eq.44)then
      i=lupi
      jj=ex(i)-nleg
      elseif(ii.eq.45)then
      jj=2*lupi
      elseif(ii.eq.46)then
      i=lupi
      j=lmap(ex(i),ey(i))
      i=vmap(ex(i),ey(i))
      jj=0
      do 750 k=1,degree(i)
      if(ovm(i,k).eq.j)then
      jj=k
      endif
  750 continue
      else
      goto 80
      endif
      elseif(ii.lt.50)then
      if(ii.eq.47)then
      i=lupi
      jj=vmap(ex(i),ey(i))-nleg
      elseif(ii.eq.48)then
      i=lupi
      jj=degree(ex(i))
      elseif(ii.eq.49)then
      i=lupi
      jj=degree(vmap(ex(i),ey(i)))
      else
      goto 80
      endif
      else
      goto 80
      endif
      ip=stcbs(1)+klin
      klin=klin+wztos(jj)
      call vaocb(klin)
      call karat(jj,ip,stcb,scbuff,1)
      else
      goto 80
      endif
      elseif(lupt.eq.3)then
      if(ii.lt.40)then
      if(ii.eq.31.or.ii.eq.33)then
      i=lupi
      j=stib(link(0)+leg(i))
      if(ii.eq.33)then
      j=stib(link(0)+j)
      endif
      if(lupi.gt.incom)then
      j=stib(link(0)+j)
      endif
      il=stib(namel(0)+j)
      klin=klin+il
      call vaocb(klin)
      ia=stib(namep(0)+j)
      ip=stcbs(1)+klin
      stcb(ip-il+1:ip)=stcb(ia:ia-1+il)
      elseif(ii.eq.32.or.ii.eq.34)then
      i=lupi
      if(ii.eq.34)then
      klin=klin+1
      call vaocb(klin)
      ip=stcbs(1)+klin
      stcb(ip:ip)='-'
      endif
      klin=klin+momel(i)
      call vaocb(klin)
      ip=stcbs(1)+klin
      stcb(ip-momel(i)+1:ip)=
     :stcb(momep(i):momep(i)-1+momel(i))
      elseif(ii.eq.35)then
      i=lupi
      j=leg(i)
      klin=klin+1
      call vaocb(klin)
      ia=stib(antiq(0)+j)
      ip=stcbs(1)+klin
      stcb(ip:ip)=char(43+ia+ia)
      else
      goto 80
      endif
      elseif(ii.lt.50)then
      if(ii.eq.40)then
      i=lupi
      if(i.le.incom)then
      jj=1
      else
      jj=2
      endif
      elseif(ii.eq.42)then
      i=lupi
      if(i.le.incom)then
      jj=1-2*i
      else
      jj=2*(incom-i)
      endif
      elseif(ii.eq.43)then
      i=lupi
      j=lmap(invp1(i),1)
      i=vmap(invp1(i),1)
      jj=0
      do 744 k=1,degree(i)
      if(ovm(i,k).eq.j)then
      jj=k
      endif
  744 continue
      elseif(ii.eq.44)then
      i=lupi
      jj=vmap(invp1(i),1)-nleg
      elseif(ii.eq.48)then
      i=lupi
      j=lmap(invp1(i),1)
      jj=degree(vmap(invp1(i),1))
      else
      goto 80
      endif
      ip=stcbs(1)+klin
      klin=klin+wztos(jj)
      call vaocb(klin)
      call karat(jj,ip,stcb,scbuff,1)
      elseif(ii.lt.60)then
      if(ii.eq.51)then
      i=lupi
      elseif(ii.eq.52)then
      i=lupi
      elseif(ii.eq.53)then
      i=lupi-incom
      else
      goto 80
      endif
      ip=stcbs(1)+klin
      klin=klin+wztos(i)
      call vaocb(klin)
      call karat(i,ip,stcb,scbuff,1)
      else
      goto 80
      endif
      elseif(lupt.eq.5.or.lupt.eq.6)then
      if(ii.lt.40)then
      if(ii.eq.31.or.ii.eq.33)then
      i=nleg+lupi
      j=lupj
      j=pmap(i,ovm(i,j))
      if(ii.eq.33)then
      j=stib(link(0)+j)
      endif
      il=stib(namel(0)+j)
      klin=klin+il
      call vaocb(klin)
      ia=stib(namep(0)+j)
      ip=stcbs(1)+klin
      stcb(ip-il+1:ip)=stcb(ia:ia-1+il)
      elseif(ii.eq.32.or.ii.eq.34)then
      i=nleg+lupi
      j=lupj
      j=ovm(i,j)
      ij=vmap(i,j)
      if(ii.eq.32)then
      kk=0
      else
      kk=1
      endif
      k0=klin
      other=0
      if(ij.le.nleg)then
      ij=p1(ij)
      if(ij.gt.incom)then
      kk=1-kk
      endif
      if(kk.eq.1)then
      klin=klin+1
      call vaocb(klin)
      ip=stcbs(1)+klin
      stcb(ip:ip)='-'
      endif
      klin=klin+momel(ij)
      call vaocb(klin)
      ip=stcbs(1)+klin
      stcb(ip-momel(ij)+1:ip)=
     :stcb(momep(ij):momep(ij)-1+momel(ij))
      elseif(vmap(i,j).eq.i)then
      if(mod(j-rdeg(i),2).eq.kk)then
      klin=klin+1
      call vaocb(klin)
      ip=stcbs(1)+klin
      stcb(ip:ip)='-'
      endif
      klin=klin+momel(0)
      call vaocb(klin)
      ip=stcbs(1)+klin
      stcb(ip-momel(0)+1:ip)=
     :stcb(momep(0):momep(0)-1+momel(0))
      ij=eg(i,i)+(j-1-rdeg(i))/2
      do 172 k=nleg+1,nleg+nloop
      if(flow(ij,k).ne.0)then
      ip=stcbs(1)+klin
      klin=klin+wztos(k-nleg)
      call vaocb(klin)
      call karat(k-nleg,ip,stcb,scbuff,1)
      goto 70
      endif
  172 continue
      else
      ij=1
      if(vmap(i,j).lt.i)then
      ij=-ij
      endif
      if(kk.eq.1)then
      ij=-ij
      endif
      do 171 k=nleg+1,nleg+nloop
      jj=flow(amap(i,j),k)
      if(jj.ne.0)then
      if(other.eq.1.or.ij.eq.jj)then
      klin=klin+1
      call vaocb(klin)
      ip=stcbs(1)+klin
      stcb(ip:ip)=char(44+ij*jj)
      endif
      klin=klin+momel(0)
      call vaocb(klin)
      ip=stcbs(1)+klin
      stcb(ip-momel(0)+1:ip)=
     :stcb(momep(0):momep(0)-1+momel(0))
      klin=klin+wztos(k-nleg)
      call karat(k-nleg,ip,stcb,scbuff,1)
      other=1
      endif
  171 continue
      do 163 k=1,nleg
      jj=flow(amap(i,j),invp1(k))
      if(jj.ne.0)then
      if(k.gt.incom)then
      jj=-jj
      endif
      if(other.eq.1.or.ij.eq.jj)then
      klin=klin+1
      call vaocb(klin)
      ip=stcbs(1)+klin
      stcb(ip:ip)=char(44+ij*jj)
      endif
      klin=klin+momel(k)
      call vaocb(klin)
      ip=stcbs(1)+klin
      stcb(ip-momel(k)+1:ip)=
     :stcb(momep(k):momep(k)-1+momel(k))
      other=1
      endif
  163 continue
      if(k0.eq.klin)then
      klin=klin+1
      call vaocb(klin)
      ip=stcbs(1)+klin
      stcb(ip:ip)='0'
      endif
      endif
      elseif(ii.eq.35)then
      i=nleg+lupi
      j=lupj
      j=pmap(i,ovm(i,j))
      klin=klin+1
      call vaocb(klin)
      ia=stib(antiq(0)+j)
      ip=stcbs(1)+klin
      stcb(ip:ip)=char(43+ia+ia)
      endif
      elseif(ii.lt.50)then
      if(ii.eq.41)then
      i=nleg+lupi
      j=lupj
      k=stib(stib(vparto(0)+stib(vfo(i)))+j)
      j=ovm(i,j)
      ij=vmap(i,j)
      if(ij.le.nleg)then
      ij=p1(ij)
      if(ij.le.incom)then
      jj=1-2*ij
      else
      jj=2*(incom-ij)
      endif
      else
      jj=amap(i,j)-nleg
      endif
      elseif(ii.lt.46)then
      if(ii.eq.40)then
      i=nleg+lupi
      j=lupj
      k=stib(stib(vparto(0)+stib(vfo(i)))+j)
      j=ovm(i,j)
      ij=vmap(i,j)
      if(ij.le.nleg)then
      if(p1(ij).le.incom)then
      jj=1
      else
      jj=2
      endif
      else
      jj=3
      endif
      elseif(ii.eq.42.or.ii.eq.45)then
      i=nleg+lupi
      j=lupj
      k=stib(stib(vparto(0)+stib(vfo(i)))+j)
      j=ovm(i,j)
      ij=vmap(i,j)
      if(ij.le.nleg)then
      ij=p1(ij)
      if(ij.le.incom)then
      jj=1-2*ij
      else
      jj=2*(incom-ij)
      endif
      elseif(k.lt.stib(link(0)+k))then
      jj=2*(amap(i,j)-nleg)-1
      elseif(k.gt.stib(link(0)+k))then
      jj=2*(amap(i,j)-nleg)
      elseif(i.lt.ij)then
      jj=2*(amap(i,j)-nleg)-1
      elseif(i.gt.ij)then
      jj=2*(amap(i,j)-nleg)
      elseif(mod(j-rdeg(i),2).eq.1)then
      jj=2*(amap(i,j)-nleg)-1
      else
      jj=2*(amap(i,j)-nleg)
      endif
      if(ii.eq.45)then
      if(jj.gt.0)then
      jj=jj+mod(jj,2)-mod(jj+1,2)
      else
      jj=0
      endif
      endif
      elseif(ii.eq.43)then
      jj=lupj
      elseif(ii.eq.44)then
      jj=lupi
      else
      goto 80
      endif
      elseif(ii.lt.50)then
      if(ii.eq.46)then
      i=nleg+lupi
      j=lupj
      ij=vmap(i,ovm(i,j))
      j=lmap(i,ovm(i,j))
      jj=0
      if(ij.gt.nleg)then
      do 754 k=1,degree(ij)
      if(ovm(ij,k).eq.j)then
      jj=k
      endif
  754 continue
      endif
      elseif(ii.eq.47)then
      i=nleg+lupi
      j=lupj
      jj=vmap(i,ovm(i,j))-nleg
      if(jj.lt.0)then
      jj=0
      endif
      elseif(ii.eq.48)then
      jj=degree(nleg+lupi)
      elseif(ii.eq.49)then
      i=nleg+lupi
      j=lupj
      jj=vmap(i,ovm(i,j))
      if(jj.le.nleg)then
      jj=0
      else
      jj=degree(jj)
      endif
      else
      goto 80
      endif
      else
      goto 80
      endif
      ip=stcbs(1)+klin
      klin=klin+wztos(jj)
      call vaocb(klin)
      call karat(jj,ip,stcb,scbuff,1)
      elseif(ii.lt.60)then
      else
      goto 80
      endif
      endif
      else
      goto 80
      endif
      elseif(k.eq.-6)then
      if(ii.lt.1.or.ii.gt.nudk)then
      goto 80
      endif
      if(stib(udkt(0)+ii).eq.1)then
      ij=stib(udki(0)+ii)
      if(ij.lt.1.or.ij.gt.ngmk)then
      goto 80
      elseif(stib(gmkd(0)+ij).ne.1)then
      goto 80
      endif
      ij=stib(gmko(0)+ij)+1
      i=stib(gmkvp(0)+ij)
      j=stib(gmkvl(0)+ij)
      if(j.gt.0)then
      klin=klin+j
      call vaocb(klin)
      ip=stcbs(1)+klin
      stcb(ip-j+1:ip)=stcb(i:i-1+j)
      endif
      if(j.lt.0.or.stib(ig+3).ne.0)then
      goto 80
      endif
      elseif(stib(udkt(0)+ii).eq.2)then
      if(lupt.eq.3)then
      i=lupi
      j=leg(i)
      if(stib(ig+3).eq.0)then
      j=stib(link(0)+j)
      endif
      if(lupi.gt.incom)then
      j=stib(link(0)+j)
      endif
      elseif(lupt.eq.4)then
      i=lupi
      j=pmap(ex(i),ey(i))
      if(stib(ig+3).ne.0)then
      j=stib(link(0)+j)
      endif
      elseif(lupt.eq.5)then
      goto 80
      elseif(lupt.eq.6)then
      i=nleg+lupi
      j=lupj
      j=pmap(i,ovm(i,j))
      if(stib(ig+3).ne.0)then
      j=stib(link(0)+j)
      endif
      else
      goto 80
      endif
      ij=stib(udki(0)+ii)
      if(ij.lt.1.or.ij.gt.npmk)then
      goto 80
      elseif(j.lt.1.or.j.gt.npart)then
      goto 80
      endif
      il=stib(stib(pmkvlp(0)+ij)+j)
      if(il.gt.0)then
      klin=klin+il
      call vaocb(klin)
      ia=stib(stib(pmkvpp(0)+ij)+j)
      ip=stcbs(1)+klin
      stcb(ip-il+1:ip)=stcb(ia:ia-1+il)
      endif
      if(il.lt.0)then
      goto 80
      endif
      elseif(stib(udkt(0)+ii).eq.3)then
      if(lupt.eq.3)then
      i=lupi
      jj=vmap(invp1(i),1)
      elseif(lupt.eq.4)then
      i=lupi
      if(stib(ig+3).eq.0)then
      jj=ex(i)
      else
      jj=vmap(ex(i),ey(i))
      endif
      elseif(lupt.eq.5)then
      jj=nleg+lupi
      elseif(lupt.eq.6)then
      jj=nleg+lupi
      else
      goto 80
      endif
      j=stib(vfo(jj))
      ij=stib(udki(0)+ii)
      if(ij.lt.1.or.ij.gt.nvmk)then
      goto 80
      elseif(j.lt.1.or.j.gt.nvert)then
      goto 80
      endif
      il=stib(stib(vmkvlp(0)+ij)+j)
      if(il.gt.0)then
      klin=klin+il
      call vaocb(klin)
      ia=stib(stib(vmkvpp(0)+ij)+j)
      ip=stcbs(1)+klin
      stcb(ip-il+1:ip)=stcb(ia:ia-1+il)
      endif
      if(il.lt.0)then
      goto 80
      endif
      else
      goto 80
      endif
      elseif(k.eq.eoia)then
      goto 20
      else
      goto 80
      endif
   70 continue
      if(xstep.eq.1)then
      ig=ig+stib(ig+2)
      else
      ig=stib(ig+3)
      endif
      goto 10
   20 continue
      j=stcbs(1)+1
      do 30 i=stcbs(1)+1,stcbs(1)+klin
      if(stcb(i:i).eq.char(lfeed))then
      ios=1
      if(i.gt.j)then
      write(unit=ounit,fmt=404,iostat=ios)stcb(j:i-1)
      else
      write(unit=ounit,fmt=404,iostat=ios)
      endif
      if(ios.ne.0)then
      auxlin(1:srec)='run-time error while writing to file'
      call messag(1,0,0,0)
      endif
      j=i+1
      endif
   30 continue
  404 format(a)
      return
   80 continue
      auxlin(1:srec)='run-time error while processing <diagram>'
      call messag(1,0,0,0)
      end
      subroutine style
      implicit integer(a-z)
      save
      parameter ( maxtak=2 )
      parameter ( maxleg=10, maxrho=6, maxi=10 )
      parameter ( srec=81 )
      parameter ( lfeed=10, dquote=34, rquote=39 )
      parameter ( eoia=-63 )
      parameter ( sibuff=1048576, scbuff=524288, mimal=2, mcmal=2 )
      common/z6in/dunit,munit,ounit,sunit,funit
      common/z7in/lsfile,sfilea,sfileb
      character*(srec) auxlin
      common/z22g/auxlin
      common/z23g/tak(1:maxtak),ks
      common/z24g/iogp(1:4)
      common/z26g/kes(0:0),kloo(15:32,1:4),kle(0:0),
     :pstke(0:0),wstke(0:0)
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z35g/stibs(0:mimal),stcbs(0:mcmal),imal(0:0),cmal(0:0)
      common/z36g/udkp(0:2),udkl(0:0),udkt(0:0),udki(0:0)
      integer lg(srec)
      character*(srec) line
      logical lunit
      ks=0
      level=0
      bevel=0
      do 6 i=1,srec
      lg(i)=0
    6 continue
      call aoib(2)
      udkp(1)=stibs(1)-1
      stib(udkp(1))=eoia
      stib(stibs(1))=eoia
      ios=1
      lunit=.false.
      inquire(file=stcb(sfilea:sfileb),exist=lunit,iostat=ios)
      if(ios.ne.0)then
      auxlin(1:srec)='system/filesystem error'
      call messag(1,0,0,0)
      elseif(.not.lunit)then
      auxlin(1:srec)='style file does not exist'
      call messag(1,0,0,0)
      endif
      ios=1
      sunit=0
      lunit=.true.
   11 continue
      sunit=sunit+1
      inquire(unit=sunit,opened=lunit,iostat=ios)
      if(ios.ne.0)then
      sunit=0
      auxlin(1:srec)='system/filesystem error'
      call messag(1,0,0,0)
      elseif(lunit)then
      if(sunit.ge.99)then
      sunit=0
      auxlin(1:srec)='no logical unit number available'
      call messag(1,0,0,0)
      endif
      goto 11
      endif
      ios=1
      open(unit=sunit,file=stcb(sfilea:sfileb),status='old',
     :access='sequential',iostat=ios)
      if(ios.ne.0)then
      sunit=0
      auxlin(1:srec)='style file could not be opened'
      call messag(1,0,0,0)
      endif
      nline=0
   70 continue
      nline=nline+1
      ios=1
      read(unit=sunit,fmt=80,iostat=ios)line(1:srec)
   80 format(a)
      if(ios.ne.0)then
      if(ios.gt.0)then
      auxlin(1:srec)='run-time error while reading style file'
      call messag(1,0,0,0)
      else
      goto 444
      endif
      endif
      if(nline.eq.1)then
      i=srec
      call vaocb(i)
      stcb(stcbs(1)+1:stcbs(1)+srec)=line(1:srec)
      call bomf(sunit)
      line(1:srec)=stcb(stcbs(1)+1:stcbs(1)+srec)
      endif
      if(line(srec:srec).ne.' ')then
      goto 300
      endif
      do 75 i=1,srec-1
      j=ichar(line(i:i))
      if(j.lt.32.or.j.gt.126)then
      goto 300
      endif
   75 continue
      if(level.eq.0.or.level.eq.4)then
      if(line(1:1).eq.'#'.or.
     :line(1:1).eq.'%'.or.
     :line(1:1).eq.'*')then
      goto 70
      endif
      endif
      iw=srec-1
   16 continue
      if(line(iw:iw).eq.' ')then
      iw=iw-1
      if(iw.gt.0)then
      goto 16
      endif
      endif
      if(level.eq.0.or.level.eq.4)then
      if(iw.eq.0)then
      goto 70
      elseif(level.eq.4)then
      goto 300
      endif
      endif
      if(level.eq.0)then
      if(line(1:1).ne.'<')then
      goto 300
      endif
      endif
      j1=0
      j2=0
      k1=0
      k2=0
      nlg=0
      do 27 i=1,iw
      if(j2.eq.0.and.line(i:i).eq.'<')then
      j1=j1+1
      if(j1.eq.2)then
      call aoib(1)
      stib(stibs(1))=ichar(line(i:i))
      j1=0
      elseif(line(i+1:i+1).ne.'<')then
      if(j2.ne.0)then
      goto 300
      endif
      nlg=nlg+1
      if(mod(nlg,2).ne.1)then
      goto 300
      endif
      lg(nlg)=i
      endif
      elseif(j1.eq.0.and.line(i:i).eq.'[')then
      j2=j2+1
      if(j2.eq.2)then
      call aoib(1)
      stib(stibs(1))=ichar(line(i:i))
      j2=0
      elseif(line(i+1:i+1).ne.'[')then
      if(j1.ne.0)then
      goto 300
      endif
      nlg=nlg+1
      if(mod(nlg,2).ne.1)then
      goto 300
      endif
      lg(nlg)=i
      endif
      elseif(j2.eq.0.and.line(i:i).eq.'>')then
      if(j1.eq.0)then
      k1=k1+1
      if(k1.eq.2)then
      call aoib(1)
      stib(stibs(1))=ichar(line(i:i))
      k1=0
      elseif(line(i+1:i+1).ne.'>')then
      goto 300
      endif
      else
      j1=0
      if(j2.ne.0)then
      goto 300
      endif
      nlg=nlg+1
      if(mod(nlg,2).ne.0)then
      goto 300
      endif
      lg(nlg)=i
      ii=lg(nlg)-lg(nlg-1)-1
      if(ii.lt.1.or.ii.gt.srec-2)then
      goto 300
      endif
      ii=styki(line,lg(nlg-1)+1,lg(nlg)-1,level)
      if(ii.eq.0)then
      goto 300
      endif
      if(ii.gt.0.and.ii.lt.5)then
      if(nlg.ne.2.or.lg(1).ne.1.or.lg(2).ne.iw)then
      goto 300
      endif
      if(ks.ne.0)then
      goto 300
      endif
      level=ii
      else
      call ktoc(ii,nline)
      endif
      endif
      elseif(j1.eq.0.and.line(i:i).eq.']')then
      if(j2.eq.0)then
      k2=k2+1
      if(k2.eq.2)then
      call aoib(1)
      stib(stibs(1))=ichar(line(i:i))
      k2=0
      elseif(line(i+1:i+1).ne.']')then
      goto 300
      endif
      else
      j2=0
      if(j1.ne.0)then
      goto 300
      endif
      nlg=nlg+1
      if(mod(nlg,2).ne.0)then
      goto 300
      endif
      lg(nlg)=i
      lw=lg(nlg)-lg(nlg-1)-1
      if(lw.lt.1)then
      goto 300
      elseif(lw.gt.srec)then
      goto 300
      endif
      ii=udstyk(line,lg(nlg-1)+1,lg(nlg)-1)
      if(ii.eq.0)then
      goto 300
      endif
      endif
      elseif(j1.eq.0.and.j2.eq.0)then
      call aoib(1)
      stib(stibs(1))=ichar(line(i:i))
      endif
   27 continue
      if(mod(nlg,2).ne.0)then
      goto 300
      endif
      if(level.gt.1.or.bevel.eq.level)then
      call aoib(1)
      endif
      if(bevel.eq.level)then
      stib(stibs(1))=lfeed
      else
      bevel=level
      if(level.gt.1)then
      stib(stibs(1))=eoia
      endif
      iogp(level)=stibs(1)+1
      endif
      goto 70
  444 continue
      ios=1
      close(unit=sunit,status='keep',iostat=ios)
      if(ios.ne.0)then
      auxlin(1:srec)='system/filesystem error'
      call messag(1,0,0,0)
      endif
      sunit=0
      if(level.lt.4)then
      auxlin(1:srec)='style file is incomplete'
      call messag(1,0,0,0)
      endif
      if(ks.eq.0)then
      call rsfki
      return
      endif
  300 continue
      auxlin(1:srec)='wrong syntax or line too long'
      call messag(1,nline,nline,4)
      end
      subroutine ktoc(ii,nlin)
      implicit integer(a-z)
      save
      parameter ( maxtak=2 )
      parameter ( srec=81 )
      parameter ( lfeed=10, dquote=34, rquote=39 )
      parameter ( sibuff=1048576, scbuff=524288, mimal=2, mcmal=2 )
      character*(srec) auxlin
      common/z22g/auxlin
      common/z23g/tak(1:maxtak),ks
      common/z30g/stib(1:sibuff)
      common/z35g/stibs(0:mimal),stcbs(0:mcmal),imal(0:0),cmal(0:0)
      integer astak(1:maxtak)
      if(ii.lt.10.or.ii.gt.82)then
      auxlin(1:srec)='fatal error in subroutine ktoc'
      call messag(1,0,0,0)
      endif
      if(ii.le.30)then
      if(ii.lt.18)then
      if(ks.ge.maxtak)then
      auxlin(1:srec)='too many nested loops'
      call messag(1,nlin,nlin,4)
      endif
      call aoib(6)
      stib(stibs(1)-5)=-1
      stib(stibs(1)-4)=ii
      stib(stibs(1)-3)=6
      stib(stibs(1)-2)=0
      stib(stibs(1)-1)=0
      stib(stibs(1))=0
      ks=ks+1
      tak(ks)=ii
      astak(ks)=stibs(1)-5
      if(ii.eq.11.or.ii.eq.12)then
      if(ii.eq.11)then
      if(ks.ne.1)then
      goto 111
      endif
      else
      if(ks.ne.2)then
      goto 111
      elseif(tak(ks-1).ne.11)then
      goto 111
      endif
      endif
      elseif(ii.lt.17)then
      if(ks.ne.1)then
      goto 111
      endif
      elseif(ii.eq.17)then
      if(ks.ne.2)then
      goto 111
      elseif(tak(ks-1).ne.16)then
      goto 111
      endif
      endif
      elseif(ii.eq.19)then
      if(ks.lt.1)then
      goto 111
      endif
      call aoib(4)
      stib(stibs(1)-3)=-1
      stib(stibs(1)-2)=ii
      stib(stibs(1)-1)=4
      stib(stibs(1))=astak(ks)
      stib(astak(ks)+3)=stibs(1)+1
      ks=ks-1
      elseif(ii.eq.21.or.ii.eq.22)then
      call aoib(3)
      stib(stibs(1)-2)=-2
      stib(stibs(1)-1)=ii
      stib(stibs(1))=3
      else
      goto 111
      endif
      endif
      if(ii.gt.30)then
      if(ii.lt.60.or.ii.eq.82)then
      call aoib(3)
      stib(stibs(1)-2)=-4
      stib(stibs(1)-1)=ii
      stib(stibs(1))=3
      elseif(ii.lt.90)then
      call aoib(3)
      stib(stibs(1)-2)=-3
      stib(stibs(1)-1)=ii
      stib(stibs(1))=3
      endif
      endif
      if(ii.eq.82)then
      if(ks.ne.1.and.ks.ne.2)then
      goto 111
      endif
      endif
      if(ii.lt.31.or.ii.gt.60)then
      return
      endif
      j=0
      if(ks.gt.0)then
      j=tak(ks)
      endif
      if(styk(auxlin,0,0,ii,j).eq.0)then
      goto 111
      endif
      return
  111 continue
      auxlin(1:srec)='wrong syntax'
      call messag(1,nlin,nlin,4)
      end
      integer function styk(s,s1,s2,k,pl)
      implicit integer(a-z)
      save
      parameter ( srec=81 )
      parameter ( sibuff=1048576, scbuff=524288, mimal=2, mcmal=2 )
      character*(srec) s
      common/z26g/kes(0:0),kloo(15:32,1:4),kle(0:0),
     :pstke(0:0),wstke(0:0)
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      if(k.lt.stib(kes(0)+15).or.k.gt.stib(kes(0)+32))then
      goto 200
      elseif(lupty(pl).eq.0)then
      goto 200
      endif
      do 110 i=15,32
      if(stib(kes(0)+i).eq.k)then
      styk=kloo(i,lupty(pl)-2)
      goto 300
      endif
  110 continue
  200 continue
      styk=0
  300 continue
      if(k.eq.52)then
      if(pl.ne.13)then
      styk=0
      endif
      elseif(k.eq.53.or.k.eq.51)then
      if(pl.ne.14)then
      styk=0
      endif
      endif
      end
      integer function styki(s,s1,s2,level)
      implicit integer(a-z)
      save
      parameter ( srec=81 )
      parameter ( sibuff=1048576, scbuff=524288, mimal=2, mcmal=2 )
      character*(srec) s
      common/z4in/mflag(1:21)
      character*(srec) auxlin
      common/z22g/auxlin
      common/z26g/kes(0:0),kloo(15:32,1:4),kle(0:0),
     :pstke(0:0),wstke(0:0)
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      if(s1.gt.s2)then
      goto 100
      endif
      styki=0
      call mstr0(s,s1,s2,pstke(0),wstke(0),jk)
      if(jk.ne.0)then
      styki=stib(kes(0)+jk)
      if(styki.le.4)then
      if(styki.ne.level+1)then
      styki=0
      endif
      elseif(level.le.3)then
      k=1
      do 10 i=2,level
      k=k*2
   10 continue
      if(mod(stib(kle(0)+jk),2*k).lt.k)then
      styki=0
      endif
      endif
      endif
      if(styki.eq.32.or.styki.eq.34)then
      mflag(17)=1
      endif
      return
  100 continue
      auxlin(1:srec)='wrong input for function styki'
      call messag(1,0,0,0)
      end
      integer function udstyk(s,s1,s2)
      implicit integer(a-z)
      save
      parameter ( srec=81 )
      parameter ( lfeed=10, dquote=34, rquote=39 )
      parameter ( eoia=-63 )
      parameter ( sibuff=1048576, scbuff=524288, mimal=2, mcmal=2 )
      parameter ( maxtak=2 )
      common/z23g/tak(1:maxtak),ks
      common/z24g/iogp(1:4)
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z35g/stibs(0:mimal),stcbs(0:mcmal),imal(0:0),cmal(0:0)
      common/z36g/udkp(0:2),udkl(0:0),udkt(0:0),udki(0:0)
      common/z38g/gmkp(0:0),gmkl(0:0),gmkd(0:0),gmko(0:0),
     :gmkvp(0:0),gmkvl(0:0)
      common/z39g/pmkp(0:0),pmkl(0:0),pmkd(0:0),pmkvpp(0:0),pmkvlp(0:0)
      common/z40g/vmkp(0:0),vmkl(0:0),vmkvpp(0:0),vmkvlp(0:0),vmks(0:0)
      common/z41g/nudk,ngmk,npmk,nvmk,tgmkd,tpmkd
      common/z45g/qdatp,qdatl,qvp,qvl,dprefp,dprefl
      character*(srec) s,sloc
      integer ludkpa
      data ini/0/
      data kstep1/4/
      data kstep2/8/
      if(ini.eq.0)then
      nudk=0
      ludkpa=udkp(1)
      ini=1
      endif
      kl=0
      if(ks.gt.0)then
      kl=lupty(tak(ks))
      endif
      udstyk=0
      if(kl.lt.3.or.kl.gt.6)then
      if(kl.ne.0)then
      goto 70
      endif
      kl=2
      endif
      kd=0
      i2=s1-1+dprefl
      if(s2.ge.i2)then
      if(s(s1:i2).eq.stcb(dprefp:dprefp-1+dprefl))then
      if(s2.gt.i2)then
      kd=1
      else
      goto 70
      endif
      endif
      endif
      i1=s1+kd*dprefl
      if(stdw(s,i1,s2).eq.0)then
      goto 70
      endif
      lsloc=s2-i1+1
      sloc(1:lsloc)=s(i1:s2)
      if(nudk.gt.0)then
      call mstr1(sloc,1,lsloc,udkp(1),1,2,udstyk)
      endif
      j=stibs(1)
      call aoib(kstep1)
      stib(j+1)=-6
      if(udstyk.eq.0)then
      nudk=nudk+1
      stib(j+2)=nudk
      else
      stib(j+2)=udstyk
      endif
      stib(j+3)=kstep1
      stib(j+4)=kd
      if(udstyk.ne.0)then
      ii=udkp(1)
      do 111 i=1,udstyk
      ii=stib(ii)
  111 continue
      ii=ii+1+kl
      if(stib(ii).eq.0)then
      stib(ii)=1+kd
      elseif(stib(ii).eq.2-kd)then
      stib(ii)=3
      endif
      return
      endif
      udstyk=stib(j+2)
      j=stibs(1)
      call aoib(kstep2)
      stib(j-1)=stib(j-1)+kstep2
      stib(ludkpa)=j+1
      ludkpa=j+1
      stib(ludkpa)=eoia
      call aocb(lsloc+1)
      stib(j+2)=stcbs(1)-lsloc
      stcb(stib(j+2):stcbs(1))=sloc(1:lsloc)//char(lfeed)
      stib(j+3)=lsloc
      stib(j+4)=0
      stib(j+5)=0
      stib(j+6)=0
      stib(j+7)=0
      stib(j+8)=0
      stib(j+2+kl)=kd+1
   70 continue
      end
      subroutine stpa(s,ia,ib,ds)
      implicit integer(a-z)
      save
      parameter ( srec=81 )
      parameter ( lfeed=10, dquote=34, rquote=39 )
      parameter ( sibuff=1048576, scbuff=524288, mimal=2, mcmal=2 )
      common/z21g/punct1(0:127)
      character*(srec) auxlin
      common/z22g/auxlin
      common/z30g/stib(1:sibuff)
      common/z35g/stibs(0:mimal),stcbs(0:mcmal),imal(0:0),cmal(0:0)
      character*(*) s
      if(ia.lt.1.or.ib.lt.ia.or.ds.lt.0)then
      auxlin(1:srec)='input for subroutine stpa is invalid'
      call messag(1,0,0,0)
      endif
      i1=0
      i2=0
      it=0
      plic=0
      sts=ds+1
      call aoib(sts)
      call aoib(-sts)
      st0=stibs(1)+sts
      s0=st0
      stib(s0)=0
      do 100 i=ia,ib+1
      if(i.le.ib)then
      ic=ichar(s(i:i))
      else
      ic=lfeed
      endif
      if(char(ic).eq.char(rquote))then
      plic=1-plic
      if(plic.eq.1)then
      if(it.eq.0)then
      it=2
      i1=i
      endif
      endif
      i2=i
      elseif(plic.eq.1)then
      if(ic.eq.lfeed)then
      stib(s0)=-1
      endif
      i2=i
      elseif(punct1(ic).ne.0)then
      if(it.eq.2)then
      if(stds(s,i1,i2,0,1).ge.0)then
      it=2
      elseif(stdw(s,i1,i2).eq.1)then
      it=3
      elseif(stdz(s,i1,i2).eq.1)then
      it=4
      elseif(stdq(s,i1,i2).eq.1)then
      it=5
      else
      stib(s0)=-1
      goto 90
      endif
      sts=sts+3
      call aoib(sts)
      call aoib(-sts)
      st0=st0+3
      stib(st0-2)=it
      stib(st0-1)=i1
      stib(st0)=i2
      stib(s0)=stib(s0)+1
      i1=0
      it=0
      endif
      sts=sts+3
      call aoib(sts)
      call aoib(-sts)
      st0=st0+3
      stib(st0-2)=1
      stib(st0-1)=punct1(ic)
      stib(st0)=i
      if(stib(s0).ge.0)then
      stib(s0)=stib(s0)+1
      endif
      elseif(ic.eq.lfeed.or.ic.eq.32)then
      if(it.eq.2)then
      if(stds(s,i1,i2,0,1).ge.0)then
      it=2
      elseif(stdw(s,i1,i2).eq.1)then
      it=3
      elseif(stdz(s,i1,i2).eq.1)then
      it=4
      elseif(stdq(s,i1,i2).eq.1)then
      it=5
      else
      stib(s0)=-1
      goto 90
      endif
      sts=sts+3
      call aoib(sts)
      call aoib(-sts)
      st0=st0+3
      stib(st0-2)=it
      stib(st0-1)=i1
      stib(st0)=i2
      stib(s0)=stib(s0)+1
      i1=0
      it=0
      endif
      else
      if(it.eq.0)then
      it=2
      i1=i
      elseif(it.ne.2)then
      stib(s0)=-1
      endif
      i2=i
      endif
   90 continue
      if(stib(s0).eq.-1)then
      return
      endif
  100 continue
      if(plic.ne.0)then
      stib(s0)=-1
      endif
      end
      subroutine hsort(x,ia,ib)
      implicit integer(a-z)
      save
      parameter ( srec=81 )
      character*(srec) auxlin
      common/z22g/auxlin
      integer x(1)
      if(ia.lt.1.or.ib.lt.ia)then
      auxlin(1:srec)='wrong input for subroutine hsort'
      call messag(1,0,0,0)
      elseif(ia.eq.ib)then
      return
      endif
      n=ib-(ia-1)
      i0=ia-1
      hl=1+n/2
      hr=n
   20 continue
      if(hl.gt.1)then
      hl=hl-1
      xtmp=x(i0+hl)
      else
      xtmp=x(i0+hr)
      x(i0+hr)=x(i0+1)
      hr=hr-1
      if(hr.eq.1)then
      x(i0+1)=xtmp
      goto 500
      endif
      endif
      hj=hl
   40 continue
      hi=hj
      hj=hi+hi
      if(hj.le.hr)then
      if(hj.lt.hr)then
      if(x(i0+hj).lt.x(i0+hj+1))then
      hj=hj+1
      endif
      endif
      if(xtmp.lt.x(i0+hj))then
      x(i0+hi)=x(i0+hj)
      goto 40
      endif
      endif
      x(i0+hi)=xtmp
      goto 20
  500 continue
      end
      subroutine xipht(ia,ib,ixipht)
      implicit integer(a-z)
      save
      parameter ( sibuff=1048576, scbuff=524288, mimal=2, mcmal=2 )
      parameter ( srec=81 )
      character*(srec) auxlin
      common/z22g/auxlin
      common/z30g/stib(1:sibuff)
      jj=0
      if(ib.gt.sibuff)then
      jj=1
      elseif(ia.lt.1)then
      jj=1
      elseif(ib.lt.ia)then
      jj=1
      endif
      if(ixipht.lt.0)then
      if(ia+ixipht.lt.1)then
      jj=1
      endif
      i1=ia
      i2=ib
      i3=1
      elseif(ixipht.gt.0)then
      if(ib+ixipht.gt.sibuff)then
      jj=1
      endif
      i1=ib
      i2=ia
      i3=-1
      endif
      if(jj.ne.0)then
      auxlin(1:srec)='inconsistency found in subroutine xipht'
      call messag(1,0,0,0)
      endif
      if(ixipht.ne.0)then
      do 10 i=i1,i2,i3
      stib(i+ixipht)=stib(i)
   10 continue
      endif
      end
      subroutine cxipht(ia,ib,ixipht)
      implicit integer(a-z)
      save
      parameter ( sibuff=1048576, scbuff=524288, mimal=2, mcmal=2 )
      parameter ( srec=81 )
      character*(srec) auxlin
      common/z22g/auxlin
      character*(scbuff) stcb
      common/z31g/stcb
      jj=0
      if(ib.gt.scbuff)then
      jj=1
      elseif(ia.lt.1)then
      jj=1
      elseif(ib.lt.ia)then
      jj=1
      endif
      if(ixipht.lt.0)then
      if(ia+ixipht.lt.1)then
      jj=1
      endif
      i1=ia
      i2=ib
      i3=1
      elseif(ixipht.gt.0)then
      if(ib+ixipht.gt.scbuff)then
      jj=1
      endif
      i1=ib
      i2=ia
      i3=-1
      endif
      if(jj.ne.0)then
      auxlin(1:srec)='inconsistency found in subroutine cxipht'
      call messag(1,0,0,0)
      endif
      if(ixipht.ne.0)then
      do 10 i=i1,i2,i3
      stcb(i+ixipht:i+ixipht)=stcb(i:i)
   10 continue
      endif
      end
      subroutine trm(ia,ib)
      implicit integer(a-z)
      save
      parameter ( srec=81 )
      parameter ( sibuff=1048576, scbuff=524288, mimal=2, mcmal=2 )
      character*(srec) auxlin
      common/z22g/auxlin
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z35g/stibs(0:mimal),stcbs(0:mcmal),imal(0:0),cmal(0:0)
      if(ia.lt.1.or.ib.lt.1)then
      auxlin(1:srec)='input for subroutine trm is invalid'
      call messag(1,0,0,0)
      endif
      ab=ia*ib
      st0=stibs(1)-ab+1
      do 20 i1=1,ab-2
      ix=0
      k=i1
   10 continue
      q=k/ia
      kk=q+ib*(k-q*ia)
      if(ix.eq.1)then
      ytmp=stib(st0+kk)
      stib(st0+kk)=xtmp
      xtmp=ytmp
      endif
      if(kk.gt.i1)then
      k=kk
      goto 10
      elseif(kk.eq.i1)then
      if(ix.eq.0.and.kk.ne.k)then
      ix=1
      k=i1
      xtmp=stib(st0+i1)
      goto 10
      endif
      endif
   20 continue
      end
      subroutine bomf(junit)
      implicit integer(a-z)
      save
      parameter ( srec=81 )
      parameter ( sibuff=1048576, scbuff=524288, mimal=2, mcmal=2 )
      character*(scbuff) stcb
      common/z31g/stcb
      character*(srec) auxlin
      common/z22g/auxlin
      common/z35g/stibs(0:mimal),stcbs(0:mcmal),imal(0:0),cmal(0:0)
      common/z56g/prom,lrom
      ii=stcbs(1)
      if(stcb(ii+1:ii+lrom).eq.stcb(prom:prom-1+lrom))then
      rewind(unit=junit,iostat=ios)
      if(ios.ne.0)then
      goto 11
      endif
      i=srec+lrom
      call vaocb(i)
      read(unit=junit,fmt=80,iostat=ios)
     :stcb(stcbs(1)+1:stcbs(1)+srec+lrom)
   80 format(a)
      if(ios.ne.0)then
      goto 11
      endif
      do 7 i=ii+1,ii+srec
      stcb(i:i)=stcb(i+lrom:i+lrom)
    7 continue
      endif
      return
   11 continue
      auxlin(1:srec)='run-time error while reading file'
      call messag(1,0,0,0)
      end
      subroutine inputs
      implicit integer(a-z)
      save
      parameter ( maxleg=10, maxrho=6, maxi=10 )
      parameter ( maxdeg=6 )
      parameter ( maxn=maxi+maxi-2, maxli=maxn+maxrho-1 )
      parameter ( naux=5 )
      parameter ( srec=81 )
      parameter ( lfeed=10, dquote=34, rquote=39 )
      parameter ( eoia=-63 )
      parameter ( sibuff=1048576, scbuff=524288, mimal=2, mcmal=2 )
      parameter ( nmarks=6 )
      common/z1in/leg(1:maxleg),nleg,nloop,incom
      common/z2in/nv(1:maxdeg)
      common/z3in/momep(0:maxleg),momel(0:maxleg)
      common/z4in/mflag(1:21)
      common/z6in/dunit,munit,ounit,sunit,funit
      common/z7in/lsfile,sfilea,sfileb
      common/z8in/lmfile,mfilea,mfileb
      common/z9in/lofile,ofilea,ofileb
      common/z10in/lffile,ffilea,ffileb
      common/z1g/nrho,rho(1:maxdeg),g(1:maxn,1:maxn)
      common/z3g/vparto(0:0),vdeg(0:0),nvert
      common/z9g/tpc(0:0)
      common/z11g/npart,nblok,nprop
      common/z14g/zcho(0:maxli),zbri(0:maxli),zpro(0:maxli),
     :rbri(0:maxli),sbri(0:maxli)
      common/z15g/dpntro(1:maxdeg)
      common/z20g/tftype(0:0),tfnarg(0:0),tfa(0:0),tfb(0:0),tfo(0:0),ntf
      common/z21g/punct1(0:127)
      character*(srec) auxlin
      common/z22g/auxlin
      common/z27g/pkey(0:0),wkey(0:0),prevl(0:0),nextl(0:0)
      common/z28g/popt3(0:0),wopt3(0:0),fopt3(0:0),vopt3(0:0)
      common/z29g/popt5(0:0),wopt5(0:0),copt5(0:0)
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z32g/acf0(0:127),acf1(0:127),sucpal(0:11,0:11)
      common/z33g/namep(0:0),namel(0:0)
      common/z34g/link(0:0),antiq(0:0),blok(0:0)
      common/z35g/stibs(0:mimal),stcbs(0:mcmal),imal(0:0),cmal(0:0)
      common/z37g/drecp(0:0),drecl(0:0),drecii(0:0),irecc(0:0),
     :frecc(0:0),ndrec,ncom
      common/z39g/pmkp(0:0),pmkl(0:0),pmkd(0:0),pmkvpp(0:0),pmkvlp(0:0)
      common/z40g/vmkp(0:0),vmkl(0:0),vmkvpp(0:0),vmkvlp(0:0),vmks(0:0)
      common/z41g/nudk,ngmk,npmk,nvmk,tgmkd,tpmkd
      common/z42g/pmkr(0:0),pmkvma(0:0),pmkvmi(0:0)
      common/z43g/vmkr(0:0),vmkmao(0:0),vmkmio(0:0)
      common/z45g/qdatp,qdatl,qvp,qvl,dprefp,dprefl
      common/z46g/popt1(0:0),wopt1(0:0),copt1(0:0)
      common/z47g/popt7(0:0),wopt7(0:0),copt7(0:0)
      integer dreci(0:0),comii(0:0),comfi(0:0)
      integer dunit
      logical lunit
      character*(srec) xtlin
      datls=4
      incom=0
      outgo=0
      nleg=0
      ios=1
      lunit=.false.
      inquire(file=stcb(qdatp:qdatp-1+qdatl),exist=lunit,iostat=ios)
      if(ios.ne.0)then
      auxlin(1:srec)='system/filesystem error'
      call messag(1,0,0,0)
      elseif(.not.lunit)then
      auxlin(1:srec)='file '//stcb(qdatp:qdatp-1+qdatl)//
     :' does not exist'
      call messag(1,0,0,0)
      endif
      dunit=0
      ios=1
      lunit=.true.
   11 continue
      dunit=dunit+1
      inquire(unit=dunit,opened=lunit,iostat=ios)
      if(ios.ne.0)then
      dunit=0
      auxlin(1:srec)='system/filesystem error'
      call messag(1,0,0,0)
      elseif(lunit)then
      if(dunit.ge.99)then
      dunit=0
      auxlin(1:srec)='no logical unit number available'
      call messag(1,0,0,0)
      endif
      goto 11
      endif
      ios=1
      open(unit=dunit,file=stcb(qdatp:qdatp-1+qdatl),status='old',
     :access='sequential',iostat=ios)
      if(ios.ne.0)then
      dunit=0
      auxlin(1:srec)='file '//stcb(qdatp:qdatp-1+qdatl)//
     :' could not be opened'
      call messag(1,0,0,0)
      endif
      ndrec=0
      newc=1
      nline=0
      dbl=0
      level=0
      phase=1
   70 continue
      i=srec
      call vaocb(i)
      ios=1
      read(unit=dunit,fmt=80,iostat=ios)stcb(stcbs(1)+1:stcbs(1)+srec)
   80 format(a)
      if(ios.ne.0)then
      if(ios.gt.0)then
      auxlin(1:srec)='run-time error while reading '//
     :stcb(qdatp:qdatp-1+qdatl)
      call messag(1,0,0,0)
      else
      ios=1
      close(unit=dunit,status='keep',iostat=ios)
      if(ios.ne.0)then
      auxlin(1:srec)='system/filesystem error'
      call messag(1,0,0,0)
      endif
      dunit=0
      if(newc.eq.0)then
      goto 521
      endif
      goto 404
      endif
      endif
      nline=nline+1
      if(nline.eq.1)then
      call bomf(dunit)
      endif
      if(newc.eq.1)then
      dbl=nline
      nplic=0
      endif
      i=stcbs(1)+srec
      if(stcb(i:i).ne.' ')then
      goto 521
      endif
      do 75 i=stcbs(1)+1,stcbs(1)-1+srec
      j=ichar(stcb(i:i))
      if(j.lt.32.or.j.gt.126)then
      goto 521
      endif
   75 continue
      ii=stcbs(1)+1
      if(stcb(ii:ii).eq.'#'.or.
     :stcb(ii:ii).eq.'%'.or.
     :stcb(ii:ii).eq.'*')then
      if(dbl.eq.nline)then
      goto 70
      else
      goto 521
      endif
      endif
      top=stcbs(1)-1+srec
  465 continue
      if(top.le.stcbs(1))then
      goto 477
      elseif(stcb(top:top).ne.' ')then
      goto 477
      endif
      top=top-1
      goto 465
  477 continue
      if(top.eq.stcbs(1))then
      if(dbl.eq.nline)then
      goto 70
      else
      goto 521
      endif
      endif
      print *,stcb(stcbs(1)+1:top)
      do 342 i=top,stcbs(1)+1,-1
      if(stcb(i:i).ne.' ')then
      if(stcb(i:i).eq.char(rquote))then
      nplic=nplic+1
      endif
      bot=i
      endif
  342 continue
      if(mod(nplic,2).ne.0)then
      goto 521
      endif
      call aoib(datls)
      ndrec=ndrec+1
      stib(stibs(1)-3)=stcbs(1)+1
      stib(stibs(1)-2)=top-stcbs(1)+1
      stib(stibs(1)-1)=nline
      stib(stibs(1))=nline-dbl
      if(stcb(top:top).eq.';')then
      newc=1
      else
      newc=0
      endif
      i=top-stcbs(1)+1
      call aocb(i)
      stcb(stcbs(1):stcbs(1))=char(lfeed)
      if(newc.eq.0)then
      goto 70
      endif
      jj=stibs(1)-datls+1
      bot=stib(jj-datls*stib(stibs(1)))
      call stpa(stcb,bot,top,0)
      if(stib(stibs(1)+1).lt.3)then
      goto 521
      elseif(stib(stibs(1)+2).ne.3)then
      goto 521
      elseif(stib(stibs(1)+5).ne.1.or.stib(stibs(1)+6).ne.1)then
      goto 521
      endif
      j=stib(stibs(1)+3)
      jj=stib(stibs(1)+4)
      if(level.ge.0)then
      call mstr0(stcb,j,jj,pkey(0),wkey(0),ij)
      if(ij.eq.0)then
      goto 521
      elseif(level.ne.stib(prevl(0)+ij))then
      goto 521
      endif
      level=stib(nextl(0)+ij)
      else
      call mstr0(stcb,j,jj,popt1(0),wopt1(0),ij)
      if(ij.eq.0)then
      goto 521
      endif
      if(stib(stibs(1)+8).ne.3)then
      goto 521
      endif
      i1=stib(stibs(1)+9)
      i2=stib(stibs(1)+10)
      call mstr0(stcb,i1,i2,popt5(0),wopt5(0),ij)
      if(ij.eq.0)then
      goto 521
      elseif(stib(copt5(0)+ij).eq.6)then
      mflag(21)=1
      endif
      endif
      goto 70
  404 continue
      if(level.ne.-1)then
      auxlin(1:srec)='missing data in input file '//
     :stcb(qdatp:qdatp-1+qdatl)
      call messag(1,0,0,0)
      endif
      i=stibs(1)
      call aoib(datls)
      do 424 j=i+1,stibs(1)
      stib(j)=eoia
  424 continue
      call trm(datls,ndrec+1)
      i=ndrec+1
      drecii(0)=stibs(1)-i
      dreci(0)=drecii(0)-i
      drecl(0)=dreci(0)-i
      drecp(0)=drecl(0)-i
      ncom=0
      i=0
  556 continue
      i=i+1
      ii=stib(drecp(0)+i)
      jj=ii-1
      j=i
  558 continue
      jj=jj+stib(drecl(0)+j)
      if(j.lt.ndrec)then
      if(stib(drecii(0)+j+1).ne.0)then
      j=j+1
      goto 558
      endif
      endif
      i=j
      ncom=ncom+1
      call aoib(2)
      stib(stibs(1)-1)=ii
      stib(stibs(1))=jj
      if(i.lt.ndrec)then
      goto 556
      endif
      call aoib(2)
      i=stibs(1)
      stib(stibs(1)-1)=eoia
      stib(stibs(1))=eoia
      call trm(2,ncom+1)
      i=ncom+1
      comfi(0)=stibs(1)-i
      comii(0)=comfi(0)-i
      i=ncom+1
      irecc(0)=stibs(1)
      call aoib(2*i)
      frecc(0)=irecc(0)+i
      stib(frecc(0))=eoia
      stib(stibs(1))=eoia
      i=0
      i1=0
  137 continue
      if(i.lt.ndrec)then
      i=i+1
      else
      goto 143
      endif
      if(stib(drecii(0)+i).eq.0)then
      if(i1.gt.0)then
      stib(frecc(0)+i1)=i-1
      endif
      i1=i1+1
      stib(irecc(0)+i1)=i
      endif
      if(i1.le.ncom)then
      goto 137
      endif
  143 continue
      stib(frecc(0)+ncom)=ndrec
      if(i.ne.ndrec.or.i1.ne.ncom)then
      auxlin(1:srec)='internal inconsistency in routine inputs'
      call messag(1,0,0,0)
      endif
      phase=2
      icom=1
      call stpa(stcb,stib(comii(0)+icom),stib(comfi(0)+icom),0)
      j=stibs(1)+8
      lofile=0
      ii=stib(stibs(1)+1)-3
      do 257 i=1,ii
      if(stib(j).ne.2)then
      goto 521
      endif
      if(i.eq.1)then
      ofilea=stib(j+1)
      endif
      if(i.eq.ii)then
      ofileb=stib(j+2)
      lofile=ofileb-ofilea+1
      endif
      j=j+3
  257 continue
      if(lofile.gt.0)then
      k=stds(stcb,ofilea,ofileb,2,1)
      if(k.lt.0)then
      goto 521
      elseif(k.eq.0)then
      lofile=0
      elseif(k.lt.lofile-2)then
      lofile=k
      ofilea=stcbs(1)+1
      call aocb(k)
      ofileb=stcbs(1)
      else
      lofile=lofile-2
      ofilea=ofilea+1
      ofileb=ofileb-1
      endif
      endif
      if(lofile.gt.0)then
      mflag(16)=1
      endif
      icom=2
      call stpa(stcb,stib(comii(0)+icom),stib(comfi(0)+icom),0)
      j=stibs(1)+8
      lsfile=0
      ii=stib(stibs(1)+1)-3
      do 259 i=1,ii
      if(stib(j).ne.2)then
      goto 521
      endif
      if(i.eq.1)then
      sfilea=stib(j+1)
      endif
      if(i.eq.ii)then
      sfileb=stib(j+2)
      lsfile=sfileb-sfilea+1
      endif
      j=j+3
  259 continue
      if(lsfile.gt.0)then
      k=stds(stcb,sfilea,sfileb,2,1)
      if(k.lt.0)then
      goto 521
      elseif(k.eq.0)then
      lsfile=0
      elseif(k.lt.lsfile-2)then
      lsfile=k
      sfilea=stcbs(1)+1
      call aocb(k)
      sfileb=stcbs(1)
      else
      lsfile=lsfile-2
      sfilea=sfilea+1
      sfileb=sfileb-1
      endif
      endif
      if(mflag(16).ne.0)then
      if(lsfile.eq.0)then
      auxlin(1:srec)='name of style file is missing in file '//
     :stcb(qdatp:qdatp-1+qdatl)
      call messag(1,0,0,0)
      endif
      call style
      endif
      lffile=0
      icom=3
      call stpa(stcb,stib(comii(0)+icom),stib(comfi(0)+icom),0)
      if(stib(stibs(1)+1).lt.4)then
      goto 521
      endif
      j=stibs(1)+8
      lmfile=0
      ii=stib(stibs(1)+1)-3
      do 261 i=1,ii
      if(stib(j).ne.2)then
      goto 521
      endif
      if(i.eq.1)then
      mfilea=stib(j+1)
      endif
      if(i.eq.ii)then
      mfileb=stib(j+2)
      lmfile=mfileb-mfilea+1
      endif
      j=j+3
  261 continue
      if(lmfile.gt.0)then
      k=stds(stcb,mfilea,mfileb,2,1)
      if(k.lt.0)then
      goto 521
      elseif(k.eq.0)then
      lmfile=0
      elseif(k.lt.lmfile-2)then
      lmfile=k
      mfilea=stcbs(1)+1
      call aocb(k)
      mfileb=stcbs(1)
      else
      lmfile=lmfile-2
      mfilea=mfilea+1
      mfileb=mfileb-1
      endif
      endif
      if(lmfile.eq.0)then
      goto 521
      endif
      call model
      icom=4
  263 continue
      call stpa(stcb,stib(comii(0)+icom),stib(comfi(0)+icom),0)
      if(stib(stibs(1)+1).eq.3)then
      nphi=0
      else
      if(stib(stibs(1)+11).eq.1.and.stib(stibs(1)+12).eq.5)then
      if(mod(stib(stibs(1)+1),5).ne.2)then
      goto 521
      endif
      nphi=(stib(stibs(1)+1)-2)/5
      if(icom.eq.4)then
      mflag(15)=1
      elseif(incom.eq.0)then
      mflag(15)=1
      elseif(mflag(15).ne.1)then
      goto 521
      endif
      j=stibs(1)+8
      do 267 i=1,nphi
      if(stib(j).ne.3)then
      goto 521
      elseif(stib(j+3).ne.1)then
      goto 521
      elseif(stib(j+4).ne.5)then
      goto 521
      elseif(stib(j+6).ne.3)then
      goto 521
      elseif(stib(j+9).ne.1)then
      goto 521
      elseif(stib(j+10).ne.6)then
      goto 521
      elseif(stib(j+12).ne.1)then
      goto 521
      elseif(stib(j+13).ne.2)then
      if(i.ne.nphi)then
      goto 521
      endif
      endif
      nleg=nleg+1
      if(nleg.gt.maxleg)then
      auxlin(1:srec)='parameters maxleg/maxrho are '//
     :'too small'
      call messag(1,0,0,0)
      endif
      momep(nleg)=stib(j+7)
      momel(nleg)=stib(j+8)-momep(nleg)+1
      j=j+15
  267 continue
      else
      if(mod(stib(stibs(1)+1),2).ne.0)then
      goto 521
      endif
      nphi=(stib(stibs(1)+1)-2)/2
      if(icom.eq.4)then
      mflag(15)=0
      elseif(mflag(15).ne.0)then
      goto 521
      endif
      j=stibs(1)+8
      do 265 i=1,nphi
      if(stib(j).ne.3)then
      goto 521
      elseif(stib(j+3).ne.1)then
      goto 521
      elseif(stib(j+4).ne.2)then
      if(i.ne.nphi)then
      goto 521
      endif
      endif
      nleg=nleg+1
      if(nleg.gt.maxleg)then
      auxlin(1:srec)='parameters maxleg/maxrho are '//
     :'too small'
      call messag(1,0,0,0)
      endif
      call aocb(1)
      k=stcbs(1)
      if(icom.eq.4)then
      stcb(k:k)=char(112)
      else
      stcb(k:k)=char(113)
      endif
      momep(nleg)=k
      call karat(i,k,stcb,scbuff,1)
      k=k+wztos(i)
      call aocb(k-stcbs(1))
      momel(nleg)=stcbs(1)-momep(nleg)+1
      j=j+6
  265 continue
      endif
      endif
      if(icom.eq.4)then
      ij=0
      else
      ij=incom
      endif
      ii=6+9*mflag(15)
      jj=stibs(1)+8-ii
      do 273 i=1,nphi
      jj=jj+ii
      do 271 i1=1,npart
      ik=stib(namel(0)+i1)
      j=stib(jj+1)
      k=stib(jj+2)
      if(ik-1.eq.k-j)then
      jk=stib(namep(0)+i1)
      if(stcb(j:k).eq.stcb(jk:jk-1+ik))then
      ij=ij+1
      if(icom.eq.4)then
      leg(ij)=stib(link(0)+i1)
      else
      leg(ij)=i1
      endif
      goto 273
      endif
      endif
  271 continue
      auxlin(1:srec)='unknown external particle(s)'
      call messag(1,0,0,0)
  273 continue
      if(icom.eq.4)then
      icom=5
      incom=nphi
      goto 263
      else
      outgo=nphi
      endif
      if(nleg.ne.incom+outgo)then
      auxlin(1:srec)='internal error in subroutine inputs'
      call messag(1,0,0,0)
      endif
      if(nleg.eq.0)then
      auxlin(1:srec)='no external particles found'
      call messag(1,0,0,0)
      endif
      icom=6
      call stpa(stcb,stib(comii(0)+icom),stib(comfi(0)+icom),0)
      if(stib(stibs(1)+1).ne.4)then
      goto 521
      elseif(stib(stibs(1)+8).ne.4)then
      goto 521
      endif
      call stoz(stcb,stib(stibs(1)+9),stib(stibs(1)+10),nloop)
      if(nloop.lt.0)then
      auxlin(1:srec)='check number of loops in file '//
     :stcb(qdatp:qdatp-1+qdatl)
      call messag(1,0,0,0)
      elseif(nloop.gt.maxrho.or.nloop+nleg.gt.max(maxleg,maxrho))then
      auxlin(1:srec)='parameters maxleg/maxrho are too small'
      call messag(1,0,0,0)
      endif
      if(nleg+nloop.lt.2)then
      auxlin(1:srec)='no diagrams exist; check nloop and/or nleg'
      call messag(1,0,0,0)
      endif
      if(nleg.eq.2.and.nloop.eq.0)then
      auxlin(1:srec)='case legs=2, loops=0 not accepted'
      call messag(1,0,0,0)
      endif
      icom=7
      call stpa(stcb,stib(comii(0)+icom),stib(comfi(0)+icom),0)
      if(stib(stibs(1)+1).eq.3)then
      call aocb(1)
      stcb(stcbs(1):stcbs(1))='k'
      momep(0)=stcbs(1)
      momel(0)=1
      elseif(stib(stibs(1)+1).eq.4)then
      if(stib(stibs(1)+8).ne.3)then
      goto 521
      endif
      momep(0)=stib(stibs(1)+9)
      momel(0)=stib(stibs(1)+10)-stib(stibs(1)+9)+1
      else
      goto 521
      endif
      if(mflag(16).ne.0)then
      call vaocb(momel(0))
      k=stcbs(1)+momel(0)
      stcb(stcbs(1)+1:k)=stcb(momep(0):momep(0)-1+momel(0))
      do 289 i=1,nloop
      ik=k
      call karat(i,ik,stcb,scbuff,1)
      ik=ik+wztos(i)
      do 287 j=1,nleg
      if(ik-momel(j).eq.k-momel(0))then
      if(stcb(stcbs(1)+1:ik).eq.
     :stcb(momep(j):momep(j)-1+momel(j)))then
      auxlin(1:srec)='conflict between names of '//
     :'external and internal momenta'
      call messag(1,0,0,0)
      endif
      endif
  287 continue
  289 continue
      endif
      icom=8
  293 continue
      call stpa(stcb,stib(comii(0)+icom),stib(comfi(0)+icom),0)
      if(stib(stibs(1)+1).eq.3)then
      ii=0
      else
      if(mod(stib(stibs(1)+1),2).ne.0)then
      goto 521
      endif
      ii=-1+stib(stibs(1)+1)/2
      endif
      do 299 i=1,ii
      ik=stibs(1)+6*i
      if(stib(ik+5).ne.1.or.stib(ik+6).ne.2)then
      if(i.ne.ii)then
      goto 521
      endif
      endif
      i1=stib(ik+3)
      i2=stib(ik+4)
      call mstr0(stcb,i1,i2,popt3(0),wopt3(0),ij)
      if(ij.ne.0)then
      if(mflag(stib(fopt3(0)+ij)).ne.0)then
      auxlin(1:srec)='check options in file '//
     :stcb(qdatp:qdatp-1+qdatl)
      call messag(0,0,0,0)
      if(mflag(stib(fopt3(0)+ij)).ne.stib(vopt3(0)+ij))then
      mflag(20)=1
      endif
      endif
      mflag(stib(fopt3(0)+ij))=stib(vopt3(0)+ij)
      else
      auxlin(1:srec)='unknown option in file '//
     :stcb(qdatp:qdatp-1+qdatl)
      goto 522
      endif
  299 continue
      ipass=1
      nphia=0
      ntfphi=0
      ntf=0
      do 55 i=0,maxli
      zpro(i)=1
      zbri(i)=1
      rbri(i)=1
      sbri(i)=1
      if(i.lt.nloop)then
      zcho(i)=0
      else
      zcho(i)=1
      endif
   55 continue
  544 continue
      do 567 jcom=icom+1,ncom
      call stpa(stcb,stib(comii(0)+jcom),stib(comfi(0)+jcom),0)
      ii=stib(stibs(1)+1)
      ik=stibs(1)+3*ii
      if(ii.lt.9)then
      goto 521
      elseif(mod(ii,2).ne.1)then
      goto 521
      elseif(stib(stibs(1)+8).ne.3)then
      goto 521
      elseif(stib(stibs(1)+11).ne.1.or.stib(stibs(1)+12).ne.5)then
      goto 521
      elseif(stib(ik-13).ne.4)then
      goto 521
      elseif(stib(ik-10).ne.1)then
      goto 521
      elseif(stib(ik-9).ne.2)then
      goto 521
      elseif(stib(ik-7).ne.4)then
      goto 521
      elseif(stib(ik-4).ne.1)then
      goto 521
      elseif(stib(ik-3).ne.6)then
      goto 521
      endif
      i1=stib(stibs(1)+3)
      i2=stib(stibs(1)+4)
      call mstr0(stcb,i1,i2,popt1(0),wopt1(0),ij)
      if(ij.le.0)then
      goto 521
      endif
      tfv=stib(copt1(0)+ij)
      if(ipass.eq.1)then
      if(stib(stibs(1)+1).gt.9)then
      ntfphi=ntfphi+1
      nphia=nphia+(stib(stibs(1)+1)-9)/2
      endif
      goto 567
      endif
      ntf=ntf+1
      stib(tftype(0)+ntf)=0
      i1=stib(stibs(1)+9)
      i2=stib(stibs(1)+10)
      call mstr0(stcb,i1,i2,popt5(0),wopt5(0),ij)
      if(ij.le.0)then
      goto 521
      endif
      stib(tftype(0)+ntf)=stib(copt5(0)+ij)*tfv
      call stoz(stcb,stib(ik-12),stib(ik-11),stib(tfa(0)+ntf))
      call stoz(stcb,stib(ik-6),stib(ik-5),stib(tfb(0)+ntf))
      if(stib(tfa(0)+ntf).gt.stib(tfb(0)+ntf))then
      goto 521
      elseif(abs(stib(tftype(0)+ntf)).lt.6)then
      if(stib(tfa(0)+ntf).lt.0)then
      goto 521
      endif
      endif
      stib(tfnarg(0)+ntf)=(stib(stibs(1)+1)-9)/2
      stib(stib(tfo(0)+ntf))=eoia
      if(ntf.lt.ntfphi)then
      stib(tfo(0)+ntf+1)=stib(tfo(0)+ntf)+stib(tfnarg(0)+ntf)+1
      endif
      if(abs(stib(tftype(0)+ntf)).gt.5)then
      goto 347
      endif
      if(stib(tfnarg(0)+ntf).gt.0)then
      do 307 i=1,stib(tfnarg(0)+ntf)
      j=stibs(1)+6*i
      if(stib(j+8).ne.3)then
      goto 521
      elseif(stib(j+11).ne.1)then
      goto 521
      elseif(stib(j+12).ne.2)then
      goto 521
      endif
      i1=stib(j+9)
      i2=stib(j+10)
      call mstr0(stcb,i1,i2,namep(0),namel(0),jk)
      if(jk.eq.0)then
      auxlin(1:srec)='unknown field in file '//
     :stcb(qdatp:qdatp-1+qdatl)
      goto 522
      endif
      stib(stib(tfo(0)+ntf)+i)=jk
  307 continue
      i4=stib(tfo(0)+ntf)
      call hsort(stib,i4+1,i4+stib(tfnarg(0)+ntf))
      tfcl=stib(tfnarg(0)+ntf)
      tfch=1
      tfcn=1
      i3=1
      do 308 i1=2,stib(tfnarg(0)+ntf)
      if(stib(i4+i1).eq.stib(i4+i1-1).or.
     :stib(i4+i1).eq.stib(link(0)+stib(i4+i1-1)))then
      i3=i3+1
      if(i3.gt.tfch)then
      tfch=i3
      endif
      else
      if(i3.lt.tfcl)then
      tfcl=i3
      endif
      tfcn=tfcn+1
      i3=1
      endif
  308 continue
      if(i3.lt.tfcl)then
      tfcl=i3
      endif
      if(tfcn.ne.nprop)then
      tfcl=0
      endif
      i1=stib(tftype(0)+ntf)
      if(i1.gt.0)then
      do 311 i2=0,maxli
      if(i2*tfch.lt.stib(tfa(0)+ntf).or.
     :(i2*tfcl.gt.stib(tfb(0)+ntf)))then
      if(i1.eq.1)then
      zpro(i2)=0
      elseif(i1.eq.2)then
      zbri(i2)=0
      elseif(i1.eq.3)then
      zcho(i2)=0
      elseif(i1.eq.4)then
      rbri(i2)=0
      elseif(i1.eq.5)then
      sbri(i2)=0
      endif
      endif
  311 continue
      elseif(i1.lt.0)then
      do 315 i2=0,maxli
      if(i2*tfch.le.stib(tfb(0)+ntf).and.
     :i2*tfcl.ge.stib(tfa(0)+ntf))then
      if(i1.eq.-1)then
      zpro(i2)=0
      elseif(i1.eq.-2)then
      zbri(i2)=0
      elseif(i1.eq.-3)then
      zcho(i2)=0
      elseif(i1.eq.-4)then
      rbri(i2)=0
      elseif(i1.eq.-5)then
      sbri(i2)=0
      endif
      endif
  315 continue
      endif
      else
      if(stib(tftype(0)+ntf).gt.0)then
      do 302 i=0,maxli
      if(i.lt.stib(tfa(0)+ntf).or.i.gt.stib(tfb(0)+ntf))then
      if(stib(tftype(0)+ntf).eq.1)then
      zpro(i)=0
      elseif(stib(tftype(0)+ntf).eq.2)then
      zbri(i)=0
      elseif(stib(tftype(0)+ntf).eq.3)then
      zcho(i)=0
      elseif(stib(tftype(0)+ntf).eq.4)then
      rbri(i)=0
      elseif(stib(tftype(0)+ntf).eq.5)then
      sbri(i)=0
      endif
      endif
  302 continue
      else
      do 304 i=0,maxli
      if(i.ge.stib(tfa(0)+ntf).and.i.le.stib(tfb(0)+ntf))then
      if(stib(tftype(0)+ntf).eq.-1)then
      zpro(i)=0
      elseif(stib(tftype(0)+ntf).eq.-2)then
      zbri(i)=0
      elseif(stib(tftype(0)+ntf).eq.-3)then
      zcho(i)=0
      elseif(stib(tftype(0)+ntf).eq.-4)then
      rbri(i)=0
      elseif(stib(tftype(0)+ntf).eq.-5)then
      sbri(i)=0
      endif
      endif
  304 continue
      endif
      endif
      goto 741
  347 continue
      if(stib(tfnarg(0)+ntf).ne.1)then
      goto 521
      endif
      i1=stib(stibs(1)+15)
      i2=stib(stibs(1)+16)
      if(abs(stib(tftype(0)+ntf)).eq.6)then
      call mstr0(stcb,i1,i2,vmkp(0),vmkl(0),kid)
      if(kid.eq.0)then
      auxlin(1:srec)='wrong argument in vsum'
      goto 522
      endif
      if(stib(vmkr(0)+kid).ne.4)then
      auxlin(1:srec)='invalid keyword in vsum'
      goto 522
      endif
      stib(vmks(0)+kid)=1
      stib(stib(tfo(0)+ntf)+1)=kid
      elseif(abs(stib(tftype(0)+ntf)).eq.7)then
      call mstr0(stcb,i1,i2,pmkp(0),pmkl(0),kid)
      if(kid.eq.0)then
      auxlin(1:srec)='wrong argument in psum'
      goto 522
      endif
      if(stib(pmkr(0)+kid).ne.4.or.stib(pmkd(0)+kid).ne.1)then
      auxlin(1:srec)='invalid keyword in psum'
      goto 522
      endif
      stib(stib(tfo(0)+ntf)+1)=kid
      jj=stib(pmkvmi(0)+kid)
      kk=stib(pmkvma(0)+kid)
      do 373 i=0,maxli
      if(stib(tftype(0)+ntf).gt.0)then
      if(i*jj.gt.stib(tfb(0)+ntf).or.
     :i*kk.lt.stib(tfa(0)+ntf))then
      zpro(i)=0
      endif
      else
      if(i*jj.ge.stib(tfa(0)+ntf).and.
     :i*kk.le.stib(tfb(0)+ntf))then
      zpro(i)=0
      endif
      endif
  373 continue
      else
      auxlin(1:srec)='in subroutine inputs'
      goto 522
      endif
  741 continue
      i1=stib(tftype(0)+ntf)
      if(abs(i1).gt.1.and.abs(i1).lt.6)then
      mflag(17)=1
      endif
      if(stib(tfnarg(0)+ntf).eq.0)then
      ntf=ntf-1
      endif
  567 continue
      if(ipass.eq.1)then
      i=ntfphi+1
      j=5*i+nphia+ntfphi
      i1=stibs(1)
      call aoib(j)
      tftype(0)=i1
      tfnarg(0)=tftype(0)+i
      tfa(0)=tfnarg(0)+i
      tfb(0)=tfa(0)+i
      tfo(0)=tfb(0)+i
      stib(tfo(0)+1)=tfo(0)+i
      ipass=2
      goto 544
      else
      endif
      stib(tfnarg(0))=eoia
      stib(tfa(0))=eoia
      stib(tfb(0))=eoia
      stib(tfo(0))=eoia
      stib(stibs(1))=eoia
      print *,' '
      print *,' -------------------------------------------------------'
      print *,' '
      jj=0
      nrb=0
      nrf=0
      do 549 i=1,npart
      if(i.le.stib(link(0)+i))then
      if(stib(antiq(0)+i).eq.0)then
      jj=jj+1
      endif
      endif
      if(i.eq.stib(link(0)+i))then
      if(stib(antiq(0)+i).eq.0)then
      nrb=nrb+1
      else
      nrf=nrf+1
      endif
      endif
  549 continue
      ii=2
      if(ii.lt.srec)then
      xtlin(1:ii)='  '
      endif
      i1=ii+wztos(nprop)
      if(i1.lt.srec)then
      call karat(nprop,ii,xtlin,srec,0)
      endif
      ii=i1+6
      if(ii.lt.srec)then
      xtlin(ii-5:ii)=char(80)//'  ---'
      endif
      if(jj.gt.0)then
      ii=ii+2
      if(ii.lt.srec)then
      xtlin(ii-1:ii)='  '
      endif
      i1=ii+wztos(jj)
      if(i1.lt.srec)then
      call karat(jj,ii,xtlin,srec,0)
      endif
      ii=i1+1
      if(ii.lt.srec)then
      xtlin(ii:ii)='+'
      endif
      endif
      if(jj.lt.nprop)then
      ii=ii+2
      if(ii.lt.srec)then
      xtlin(ii-1:ii)='  '
      endif
      i1=ii+wztos(nprop-jj)
      if(i1.lt.srec)then
      call karat(nprop-jj,ii,xtlin,srec,0)
      endif
      ii=i1+1
      if(ii.lt.srec)then
      xtlin(ii:ii)='-'
      endif
      endif
      ii=ii+5
      if(ii.lt.srec)then
      xtlin(ii-4:ii)='  ---'
      endif
      if(nrb.gt.0)then
      ii=ii+2
      if(ii.lt.srec)then
      xtlin(ii-1:ii)='  '
      endif
      i1=ii+wztos(nrb)
      if(i1.lt.srec)then
      call karat(nrb,ii,xtlin,srec,0)
      endif
      ii=i1+2
      if(ii.lt.srec)then
      xtlin(ii-1:ii)=char(78)//'+'
      endif
      endif
      if(nrb.lt.jj)then
      ii=ii+2
      if(ii.lt.srec)then
      xtlin(ii-1:ii)='  '
      endif
      i1=ii+wztos(jj-nrb)
      if(i1.lt.srec)then
      call karat(jj-nrb,ii,xtlin,srec,0)
      endif
      ii=i1+2
      if(ii.lt.srec)then
      xtlin(ii-1:ii)=char(67)//'+'
      endif
      endif
      if(nrf.gt.0)then
      ii=ii+2
      if(ii.lt.srec)then
      xtlin(ii-1:ii)='  '
      endif
      i1=ii+wztos(nrf)
      if(i1.lt.srec)then
      call karat(nrf,ii,xtlin,srec,0)
      endif
      ii=i1+2
      if(ii.lt.srec)then
      xtlin(ii-1:ii)=char(78)//'-'
      endif
      endif
      if(nrf.lt.nprop-jj)then
      ii=ii+2
      if(ii.lt.srec)then
      xtlin(ii-1:ii)='  '
      endif
      i1=ii+wztos(nprop-jj-nrf)
      if(i1.lt.srec)then
      call karat(nprop-jj-nrf,ii,xtlin,srec,0)
      endif
      ii=i1+2
      if(ii.lt.srec)then
      xtlin(ii-1:ii)=char(67)//'-'
      endif
      endif
      if(ii.ge.srec)then
      xtlin(srec-3:srec-1)='...'
      endif
      print *,xtlin(1:min(ii,srec-1))
      print *,' '
      ii=2
      if(ii.lt.srec)then
      xtlin(1:ii)='  '
      endif
      i1=ii+wztos(nvert)
      if(i1.lt.srec)then
      call karat(nvert,ii,xtlin,srec,0)
      endif
      ii=i1+6
      if(ii.lt.srec)then
      xtlin(ii-5:ii)=char(86)//'  ---'
      endif
      do 321 i=1,nrho
      if(nv(i).gt.0)then
      ii=ii+2
      if(ii.lt.srec)then
      xtlin(ii-1:ii)='  '
      endif
      i1=ii+wztos(i)
      if(i1.lt.srec)then
      call karat(i,ii,xtlin,srec,0)
      endif
      ii=i1+1
      if(ii.lt.srec)then
      xtlin(ii:ii)='^'
      endif
      i1=ii+wztos(nv(i))
      if(i1.lt.srec)then
      call karat(nv(i),ii,xtlin,srec,0)
      endif
      ii=i1
      endif
  321 continue
      if(ii.ge.srec)then
      xtlin(srec-3:srec-1)='...'
      endif
      print *,xtlin(1:min(ii,srec-1))
      print *,' '
      print *,' -------------------------------------------------------'
      print *,' '
      call vsig
      if(mflag(1).eq.1)then
      do 305 i=1,maxli
      zbri(i)=0
  305 continue
      elseif(mflag(1).eq.-1)then
      zbri(0)=0
      endif
      if(mflag(2).eq.1)then
      do 317 i=1,maxli
      sbri(i)=0
  317 continue
      elseif(mflag(2).eq.-1)then
      sbri(0)=0
      zcho(0)=0
      endif
      if(mflag(3).eq.-1)then
      zcho(0)=0
      endif
      if(mflag(4).eq.-1)then
      zbri(0)=0
      zcho(0)=0
      endif
      if(mflag(5).eq.-1)then
      zcho(0)=0
      endif
      if(mflag(6).eq.-1)then
      zcho(0)=0
      endif
      if(mflag(5).ne.0)then
      mflag(17)=1
      endif
      if(nloop.eq.0)then
      if(mflag(2).ne.-1)then
      mflag(2)=0
      else
      mflag(20)=1
      endif
      if(mflag(3).ne.-1)then
      mflag(3)=0
      else
      mflag(20)=1
      endif
      if(mflag(4).ne.-1)then
      mflag(4)=0
      else
      mflag(20)=1
      endif
      if(mflag(5).ne.-1)then
      mflag(5)=0
      else
      mflag(20)=1
      endif
      if(mflag(6).ne.-1)then
      mflag(6)=0
      else
      mflag(20)=1
      endif
      mflag(9)=0
      endif
      if(mflag(1).eq.1.and.nleg.ne.1)then
      if(mflag(2).eq.-1)then
      mflag(20)=1
      endif
      if(mflag(4).eq.-1)then
      mflag(20)=1
      endif
      endif
      ii=0
      do 124 i=1,nleg
      if(stib(antiq(0)+leg(i)).eq.1)then
      ii=1-ii
      endif
  124 continue
      if(ii.ne.0)then
      auxlin(1:srec)='odd number of external anticommuting fields'
      call messag(0,0,0,0)
      mflag(20)=1
      endif
      j=stib(blok(0)+stib(link(0)+leg(1)))
      do 146 i=2,nleg
      if(stib(blok(0)+stib(link(0)+leg(i))).ne.j)then
      auxlin(1:srec)='external particles cannot be connected'
      call messag(0,0,0,0)
      mflag(20)=1
      endif
  146 continue
      if(nleg.ne.2.or.nloop.ne.0)then
      do 156 j=1,nleg
      jj=stib(link(0)+leg(j))
      do 151 k=1,nrho
      if(nv(k).gt.0)then
      if(stib(stib(dpntro(k)+jj)+1).eq.jj)then
      goto 400
      endif
      endif
  151 continue
      if(j.gt.incom)then
      ii=stib(namep(0)+leg(j))
      jj=stib(namel(0)+leg(j))
      auxlin(1:srec)='no vertex contains outgoing particle '//
     :stcb(ii:ii-1+jj)
      else
      ii=stib(namep(0)+jj)
      jj=stib(namel(0)+jj)
      auxlin(1:srec)='no vertex contains incoming particle '//
     :stcb(ii:ii-1+jj)
      endif
      call messag(0,0,0,0)
      mflag(20)=1
  400 continue
  156 continue
      endif
      if(mflag(8).eq.1)then
      if(npart.gt.1)then
      auxlin(1:srec)='option "topol" does not apply here'
      call messag(1,0,0,0)
      endif
      endif
      do 102 i=1,npart
      if(stib(antiq(0)+i).eq.1)then
      mflag(13)=1
      endif
  102 continue
      mflag(14)=0
      do 168 i=1,nvert
      k=stib(vparto(0)+i)
      ii=0
      do 268 j=1,stib(vdeg(0)+i)
      ii=ii+stib(antiq(0)+stib(k+j))
  268 continue
      if(ii.ne.0.and.ii.ne.2)then
      mflag(14)=1
      endif
  168 continue
      if(mflag(9).ne.0.and.mflag(14).ne.0)then
      auxlin(1:srec)='option "floop" does not apply here'
      call messag(1,0,0,0)
      endif
      return
  521 continue
      auxlin(1:srec)='wrong syntax or line too long in file '//
     :stcb(qdatp:qdatp-1+qdatl)
  522 continue
      if(phase.eq.1)then
      call messag(1,dbl,nline,0)
      else
      i=1
      j=0
  523 continue
      if(stib(drecii(0)+i).eq.0)then
      j=j+1
      nline=stib(dreci(0)+i)
      endif
      if(i.lt.ndrec.and.j.lt.jcom)then
      i=i+1
      goto 523
      endif
  524 continue
      if(i.lt.ndrec)then
      if(stib(drecii(0)+i+1).ne.0)then
      i=i+1
      goto 524
      endif
      endif
      call messag(1,nline,nline+stib(drecii(0)+i),0)
      endif
      end
      subroutine model
      implicit integer(a-z)
      save
      parameter ( maxleg=10, maxrho=6, maxi=10 )
      parameter ( maxdeg=6 )
      parameter ( maxn=maxi+maxi-2, maxli=maxn+maxrho-1 )
      parameter ( srec=81 )
      parameter ( lfeed=10, dquote=34, rquote=39 )
      parameter ( eoia=-63 )
      parameter ( sibuff=1048576, scbuff=524288, mimal=2, mcmal=2 )
      common/z1in/leg(1:maxleg),nleg,nloop,incom
      common/z4in/mflag(1:21)
      common/z6in/dunit,munit,ounit,sunit,funit
      common/z8in/lmfile,mfilea,mfileb
      common/z10in/lffile,ffilea,ffileb
      common/z3g/vparto(0:0),vdeg(0:0),nvert
      common/z9g/tpc(0:0)
      common/z11g/npart,nblok,nprop
      common/z15g/dpntro(1:maxdeg)
      character*(srec) auxlin
      common/z22g/auxlin
      common/z27g/pkey(0:0),wkey(0:0),prevl(0:0),nextl(0:0)
      common/z28g/popt3(0:0),wopt3(0:0),fopt3(0:0),vopt3(0:0)
      common/z29g/popt5(0:0),wopt5(0:0),copt5(0:0)
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z32g/acf0(0:127),acf1(0:127),sucpal(0:11,0:11)
      common/z33g/namep(0:0),namel(0:0)
      common/z34g/link(0:0),antiq(0:0),blok(0:0)
      common/z35g/stibs(0:mimal),stcbs(0:mcmal),imal(0:0),cmal(0:0)
      common/z36g/udkp(0:2),udkl(0:0),udkt(0:0),udki(0:0)
      common/z38g/gmkp(0:0),gmkl(0:0),gmkd(0:0),gmko(0:0),
     :gmkvp(0:0),gmkvl(0:0)
      common/z39g/pmkp(0:0),pmkl(0:0),pmkd(0:0),pmkvpp(0:0),pmkvlp(0:0)
      common/z40g/vmkp(0:0),vmkl(0:0),vmkvpp(0:0),vmkvlp(0:0),vmks(0:0)
      common/z41g/nudk,ngmk,npmk,nvmk,tgmkd,tpmkd
      common/z46g/popt1(0:0),wopt1(0:0),copt1(0:0)
      common/z47g/popt7(0:0),wopt7(0:0),copt7(0:0)
      integer apmkd(0:0),avmkd(0:0),llp(0:0)
      integer munit
      logical lunit
      llp(0)=-1
      if(lffile.gt.0)then
      xfile=0
      else
      xfile=1
      endif
   05 continue
      xfile=xfile+1
      if(xfile.eq.1)then
      ios=1
      lunit=.false.
      inquire(file=stcb(ffilea:ffileb),exist=lunit,iostat=ios)
      if(ios.ne.0)then
      auxlin(1:srec)='system/filesystem error'
      call messag(1,0,0,0)
      elseif(.not.lunit)then
      auxlin(1:srec)='library file does not exist'
      call messag(1,0,0,0)
      endif
      funit=0
      ios=1
      lunit=.true.
   11 continue
      funit=funit+1
      inquire(unit=funit,opened=lunit,iostat=ios)
      if(ios.ne.0)then
      funit=0
      auxlin(1:srec)='system/filesystem error'
      call messag(1,0,0,0)
      elseif(lunit)then
      if(funit.ge.99)then
      funit=0
      auxlin(1:srec)='no logical unit number available'
      call messag(1,0,0,0)
      endif
      goto 11
      endif
      ios=1
      open(unit=funit,file=stcb(ffilea:ffileb),status='old',
     :access='sequential',iostat=ios)
      if(ios.ne.0)then
      funit=0
      auxlin(1:srec)='library file could not be opened'
      call messag(1,0,0,0)
      endif
      elseif(xfile.eq.2)then
      ios=1
      lunit=.false.
      inquire(file=stcb(mfilea:mfileb),exist=lunit,iostat=ios)
      if(ios.ne.0)then
      auxlin(1:srec)='system/filesystem error'
      call messag(1,0,0,0)
      elseif(.not.lunit)then
      auxlin(1:srec)='model file does not exist'
      call messag(1,0,0,0)
      endif
      munit=0
      ios=1
      lunit=.true.
   12 continue
      munit=munit+1
      inquire(unit=munit,opened=lunit,iostat=ios)
      if(ios.ne.0)then
      munit=0
      auxlin(1:srec)='system/filesystem error'
      call messag(1,0,0,0)
      elseif(lunit)then
      if(munit.ge.99)then
      munit=0
      auxlin(1:srec)='no logical unit number available'
      call messag(1,0,0,0)
      endif
      goto 12
      endif
      ios=1
      open(unit=munit,file=stcb(mfilea:mfileb),status='old',
     :access='sequential',iostat=ios)
      if(ios.ne.0)then
      munit=0
      auxlin(1:srec)='model file could not be opened'
      call messag(1,0,0,0)
      endif
      endif
      nline=0
      dbl=0
      newc=1
      level=1
      if(xfile.eq.1)then
      lastp=0
      lastkp=0
      ngmk=0
      else
      npmk=0
      nvmk=0
      npart=0
      nprop=0
      nvert=0
      endif
   65 continue
      dip=0
      nest=0
      top=0
   70 continue
      nline=nline+1
      if(newc.eq.1)then
      dbl=nline
      endif
      i=srec
      call vaocb(i)
      ios=1
   80 format(a)
      if(xfile.eq.1)then
      read(unit=funit,fmt=80,iostat=ios)stcb(stcbs(1)+1:stcbs(1)+srec)
      if(ios.ne.0)then
      if(ios.gt.0)then
      auxlin(1:srec)='run-time error while reading library file'
      call messag(1,0,0,0)
      else
      ios=1
      close(unit=funit,status='keep',iostat=ios)
      if(ios.ne.0)then
      auxlin(1:srec)='system/filesystem error'
      call messag(1,0,0,0)
      endif
      funit=0
      goto 05
      endif
      endif
      elseif(xfile.eq.2)then
      read(unit=munit,fmt=80,iostat=ios)stcb(stcbs(1)+1:stcbs(1)+srec)
      if(ios.ne.0)then
      if(ios.gt.0)then
      auxlin(1:srec)='run-time error while reading model file'
      call messag(1,0,0,0)
      else
      ios=1
      close(unit=munit,status='keep',iostat=ios)
      if(ios.ne.0)then
      auxlin(1:srec)='system/filesystem error'
      call messag(1,0,0,0)
      endif
      munit=0
      goto 777
      endif
      endif
      endif
      if(nline.eq.1)then
      call bomf(munit)
      endif
      i=stcbs(1)+srec
      if(stcb(i:i).ne.' ')then
      goto 521
      endif
      do 75 i=stcbs(1)+1,stcbs(1)-1+srec
      j=ichar(stcb(i:i))
      if(j.lt.32.or.j.gt.126)then
      goto 521
      endif
   75 continue
      ii=stcbs(1)+1
      if(stcb(ii:ii).eq.'#'.or.
     :stcb(ii:ii).eq.'%'.or.
     :stcb(ii:ii).eq.'*')then
      if(top.eq.0)then
      goto 70
      else
      goto 521
      endif
      endif
      ibot=0
      itop=0
      do 16 i=stcbs(1)+1,stcbs(1)-1+srec
      if(stcb(i:i).ne.' ')then
      itop=i
      if(ibot.eq.0)then
      ibot=i
      endif
      endif
   16 continue
      if(itop.eq.0)then
      if(top.eq.0)then
      goto 70
      else
      goto 521
      endif
      endif
      ii=0
      quote=0
      do 26 i=ibot,itop
      if(stcb(i:i).eq.char(rquote))then
      if(nest.ne.1)then
      goto 521
      endif
      if(quote.eq.0)then
      quote=rquote
      endif
      if(quote.eq.rquote)then
      ii=1-ii
      endif
      elseif(stcb(i:i).eq.char(dquote))then
      if(nest.ne.1)then
      goto 521
      endif
      if(quote.eq.0)then
      quote=dquote
      endif
      if(quote.eq.dquote)then
      ii=1-ii
      endif
      elseif(ii.eq.0)then
      if(stcb(i:i).eq.'[')then
      if(nest.ne.0)then
      goto 521
      endif
      nest=1
      elseif(stcb(i:i).eq.']')then
      if(nest.ne.1)then
      goto 521
      endif
      nest=2
      elseif(stcb(i:i).eq.',')then
      if(nest.ne.1)then
      goto 521
      endif
      quote=0
      elseif(stcb(i:i).ne.' ')then
      if(nest.ne.1)then
      goto 521
      endif
      endif
      else
      if(nest.ne.1)then
      goto 521
      endif
      endif
   26 continue
      if(dbl.eq.nline)then
      jbot=ibot
      endif
      i=itop-stcbs(1)+1
      call aocb(i)
      stcb(stcbs(1):stcbs(1))=char(lfeed)
      if(stcb(itop:itop).ne.']')then
      newc=0
      goto 70
      endif
      if(ii.ne.0)then
      goto 521
      elseif(nest.ne.2)then
      goto 521
      endif
      newc=1
      ibot=jbot
      top=stcbs(1)
   90 continue
      pal=0
      mpal=11
      vel=0
      vin=0
      idin=0
      kin=0
      kid=0
      knf=0
      id1=0
      id2=0
      is1=0
      is2=0
      plm=0
      ix=ibot
  280 continue
      ic=ichar(stcb(ix:ix))
      if(sucpal(pal,acf0(ic)).lt.0)then
      goto 521
      endif
      if(pal.lt.3)then
      if(pal.eq.1)then
      if(acf0(ic).eq.0)then
      if(id1.eq.0)then
      id1=ix
      endif
      id2=ix
      elseif(acf0(ic).gt.2)then
      if(stcb(ix:ix).eq.',')then
      if(level.eq.1)then
      if(ngmk.eq.0)then
      call aoib(1)
      llp(0)=stibs(1)
      stib(llp(0))=eoia
      endif
      call rgmki(llp(0))
      if(xfile.eq.1)then
      goto 521
      endif
      level=2
      fpass=1
      lastkp=0
      off1=0
      endif
      elseif(stcb(ix:ix).eq.';')then
      if(level.lt.2)then
      goto 521
      endif
      elseif(stcb(ix:ix).eq.'=')then
      if(level.ne.1.or.kin.ne.0)then
      goto 521
      endif
      kin=1
      goto 283
      elseif(stcb(ix:ix).eq.']')then
      if(level.lt.2)then
      goto 521
      endif
      endif
      if(level.eq.2)then
      idin=idin+1
      if(idin.gt.4)then
      goto 521
      endif
      goto 281
      elseif(level.eq.3)then
      idin=idin+1
      goto 282
      endif
      endif
      elseif(pal.eq.2)then
      if(stcb(ix:ix).eq.'=')then
      kin=kin+1
      goto 283
      elseif(acf0(ic).eq.0)then
      if(id1.eq.0)then
      id1=ix
      endif
      id2=ix
      endif
      endif
      elseif(pal.lt.6)then
      if(pal.eq.3)then
      if(vel.ne.0.or.vin.ne.0)then
      goto 521
      endif
      if(stcb(ix:ix).eq.char(rquote))then
      plm=1-plm
      if(is1.eq.0)then
      is1=ix
      endif
      is2=ix
      else
      if(acf0(ic).eq.0)then
      if(is1.eq.0)then
      is1=ix
      endif
      is2=ix
      elseif(stcb(ix:ix).eq.'(')then
      if(is1.ne.0)then
      goto 521
      endif
      vel=1
      elseif(stcb(ix:ix).eq.','.or.stcb(ix:ix).eq.']')then
      goto 285
      endif
      endif
      elseif(pal.eq.4)then
      elseif(pal.eq.5)then
      if(vel.ne.0.or.vin.ne.0.or.plm.ne.1)then
      goto 521
      endif
      if(stcb(ix:ix).eq.char(rquote))then
      plm=1-plm
      endif
      is2=ix
      endif
      elseif(pal.lt.9)then
      if(pal.eq.6)then
      if(vel.ne.1.or.plm.ne.0)then
      goto 521
      endif
      if(stcb(ix:ix).eq.char(rquote))then
      plm=1-plm
      if(is1.eq.0)then
      is1=ix
      endif
      is2=ix
      else
      if(acf0(ic).eq.0)then
      if(is1.eq.0)then
      is1=ix
      endif
      is2=ix
      elseif(stcb(ix:ix).eq.','.or.stcb(ix:ix).eq.')')then
      vin=vin+1
      goto 285
      endif
      endif
      elseif(pal.eq.7)then
      elseif(pal.eq.8)then
      if(vel.ne.1.or.plm.ne.1)then
      goto 521
      endif
      if(stcb(ix:ix).eq.char(rquote))then
      plm=1-plm
      endif
      is2=ix
      endif
      elseif(pal.lt.12)then
      if(pal.eq.9)then
      if(level.gt.2.or.level.lt.1)then
      goto 521
      endif
      if(stcb(ix:ix).eq.','.or.stcb(ix:ix).eq.']')then
      vin=0
      vel=0
      endif
      elseif(pal.eq.10)then
      elseif(pal.eq.11)then
      endif
      else
      goto 521
      endif
  287 continue
      pal=sucpal(pal,acf0(ic))
      if(ix.lt.top)then
      ix=ix+1
      goto 280
      endif
      if(vel.ne.0)then
      goto 521
      elseif(pal.ne.mpal)then
      goto 521
      elseif(plm.ne.0)then
      goto 521
      elseif(id1.ne.0)then
      goto 521
      elseif(is1.ne.0)then
      goto 521
      endif
      if(level.eq.1)then
      goto 65
      elseif(level.eq.2)then
      if(fpass.eq.1)then
      else
      do 303 i=1,npmk
      if(stib(apmkd(0)+i).eq.3)then
      if(stib(lastp+1).eq.0)then
      auxlin(1:srec)='wrong dimension for M-function'
      call messag(1,dbl,nline,4-xfile)
      endif
      endif
      if(stib(apmkd(0)+i).ne.stib(pmkd(0)+i))then
      ii=1
      if(stib(apmkd(0)+i).eq.2)then
      if(stib(pmkd(0)+i).eq.3)then
      ii=stib(lastp+1)
      endif
      endif
      if(ii.ne.0)then
      if(stib(apmkd(0)+i).eq.0)then
      auxlin(1:srec)='missing M-function'
      else
      auxlin(1:srec)='wrong dimension for M-function'
      endif
      call messag(1,dbl,nline,2)
      endif
      endif
  303 continue
      nprop=nprop+1
      endif
      elseif(level.eq.3)then
      if(fpass.eq.0)then
      do 313 i=1,nvmk
      if(stib(avmkd(0)+i).eq.0)then
      auxlin(1:srec)='missing M-function'
      call messag(1,dbl,nline,2)
      endif
  313 continue
      j=stib(lastp+1)
      if(j.lt.3)then
      auxlin(1:srec)='vertex of wrong degree'
      call messag(1,dbl,nline,2)
      endif
      ii=0
      do 331 i=1,j
      if(stib(antiq(0)+stib(lastp+1+i)).ne.0)then
      ii=1-ii
      endif
  331 continue
      if(ii.ne.0)then
      auxlin(1:srec)='odd vertex'
      call messag(1,dbl,nline,2)
      endif
      endif
      endif
      if(level.eq.2)then
      if(fpass.eq.1)then
      i=stibs(1)
      call aoib(4)
      do 424 j=i+1,stibs(1)
      stib(j)=eoia
  424 continue
      call trm(4,npmk+1)
      i=npmk+1
      pmkd(0)=stibs(1)-i
      pmkl(0)=pmkd(0)-i
      pmkp(0)=pmkl(0)-i
      i=pmkp(0)-llp(0)+1
      call xipht(pmkp(0)+1,stibs(1),-i)
      call aoib(-i)
      pmkd(0)=pmkd(0)-i
      pmkl(0)=pmkl(0)-i
      pmkp(0)=pmkp(0)-i
      psize=6+2*npmk
      apmkd(0)=stibs(1)
      call aoib(npmk+1)
      stib(stibs(1))=eoia
      endif
      do 703 i=1,npmk
      stib(apmkd(0)+i)=0
  703 continue
      if(fpass.eq.1)then
      fpass=0
      goto 90
      endif
      elseif(level.eq.3)then
      if(fpass.eq.1)then
      i=stibs(1)
      call aoib(3)
      do 427 j=i+1,stibs(1)
      stib(j)=eoia
  427 continue
      call trm(3,nvmk+1)
      i=nvmk+1
      vmkl(0)=stibs(1)-i
      vmkp(0)=vmkl(0)-i
      avmkd(0)=vmkp(0)-i
      call xipht(vmkp(0),stibs(1),-1)
      call aoib(-1)
      vmkl(0)=vmkl(0)-1
      vmkp(0)=vmkp(0)-1
      avmkd(0)=avmkd(0)-1
      else
      nvert=nvert+1
      endif
      do 705 i=1,nvmk
      stib(avmkd(0)+i)=0
  705 continue
      if(fpass.eq.1)then
      fpass=0
      goto 90
      endif
      endif
      goto 65
  281 continue
      if(idin.eq.3)then
      if(stdw(stcb,id1,id2).eq.0)then
      if(id1.ne.id2)then
      goto 521
      elseif(stcb(id1:id1).ne.'+'.and.stcb(id1:id1).ne.'-')then
      goto 521
      endif
      else
      if(knf.ne.2)then
      auxlin(1:srec)='unknown field in vertex'
      call messag(1,dbl,nline,2)
      endif
      knf=0
      call rpi(npart)
      level=3
      fpass=1
      lastp=0
      lastkp=0
      goto 90
      endif
      if(knf.ne.0)then
      auxlin(1:srec)='field appeared in earlier propagator'
      call messag(1,dbl,nline,2)
      endif
      if(fpass.ne.1)then
      if(stcb(id1:id1).eq.'+')then
      stib(off1+3)=0
      elseif(stcb(id1:id1).eq.'-')then
      stib(off1+3)=1
      else
      goto 521
      endif
      stib(off1+4)=0
      if(off2.ne.0)then
      stib(off2+3)=stib(off1+3)
      stib(off2+4)=0
      endif
      stib(lastp)=off1+1
      stib(off1+1)=eoia
      lastp=off1+1
      if(off2.ne.0)then
      stib(lastp)=off2+1
      stib(off2+1)=eoia
      lastp=off2+1
      endif
      endif
      goto 284
      elseif(idin.eq.4)then
      if(fpass.ne.1)then
      call mstr0(stcb,id1,id2,popt7(0),wopt7(0),i)
      if(i.eq.0)then
      goto 521
      endif
      stib(off1+4)=stib(copt7(0)+i)
      if(off2.ne.0)then
      stib(off2+4)=stib(off1+4)
      endif
      endif
      goto 284
      endif
      if(id1.le.0)then
      goto 521
      elseif(id1.gt.id2)then
      goto 521
      elseif(stdw(stcb,id1,id2).eq.0)then
      auxlin(1:srec)='unacceptable field name'
      call messag(1,dbl,nline,2)
      endif
      if(fpass.eq.1)then
      goto 284
      endif
      kk=id2-id1+1
      if(kk.lt.1)then
      goto 521
      endif
      newid=1
      if(npart.gt.0)then
      call mstr1(stcb,id1,id2,llp(0),4,5,ij)
      if(ij.ne.0)then
      newid=0
      endif
      endif
      if(newid.eq.0)then
      if(idin.eq.1)then
      knf=1
      elseif(idin.eq.2)then
      if(knf.eq.1)then
      knf=2
      elseif(ij.ne.npart)then
      knf=1
      endif
      endif
      endif
      if(knf.ne.0)then
      goto 284
      endif
      if(idin.eq.1)then
      if(nprop.eq.0)then
      call aoib(1)
      lastp=stibs(1)
      stib(lastp)=eoia
      llp(0)=stibs(1)
      endif
      off1=stibs(1)
      off2=0
      call aoib(psize)
      do 278 i1=stibs(1)-psize+1,stibs(1)
      stib(i1)=0
  278 continue
      stib(off1+5)=id1
      stib(off1+6)=id2-id1+1
      npart=npart+1
      elseif(idin.eq.2)then
      newid=1
      i1=stib(off1+5)
      i2=stib(off1+6)
      if(id2-id1+1.eq.i2)then
      if(stcb(id1:id2).eq.stcb(i1:i1-1+i2))then
      newid=0
      endif
      endif
      if(newid.eq.0)then
      stib(off1+2)=0
      goto 284
      endif
      off2=stibs(1)
      call aoib(psize)
      do 279 i1=stibs(1)-psize+1,stibs(1)
      stib(i1)=0
  279 continue
      stib(off2+5)=id1
      stib(off2+6)=id2-id1+1
      stib(off1+2)=1
      stib(off2+2)=-1
      npart=npart+1
      endif
  284 continue
      id1=0
      goto 287
  282 continue
      if(id1.le.0)then
      goto 521
      elseif(id1.gt.id2)then
      goto 521
      elseif(stdw(stcb,id1,id2).eq.0)then
      auxlin(1:srec)='unacceptable field name'
      call messag(1,dbl,nline,2)
      endif
      if(idin.gt.maxdeg)then
      auxlin(1:srec)='vertex degree too high'
      call messag(1,dbl,nline,2)
      endif
      if(fpass.eq.1)then
      goto 286
      endif
      if(idin.eq.1)then
      if(nvert.eq.0)then
      call aoib(1)
      lastp=stibs(1)
      stib(lastp)=eoia
      llp(0)=stibs(1)
      endif
      i=stibs(1)
      ii=2+2*nvmk
      call aoib(ii)
      stib(lastp)=i+1
      lastp=i+1
      stib(lastp)=eoia
      do 378 i1=lastp+1,stibs(1)
      stib(i1)=0
  378 continue
      endif
      call mstr0(stcb,id1,id2,namep(0),namel(0),ij)
      if(ij.eq.0)then
      auxlin(1:srec)='unknown field in vertex'
      call messag(1,dbl,nline,2)
      endif
      call aoib(1)
      stib(stibs(1))=0
      stib(lastp+1+idin)=ij
      stib(lastp+1)=stib(lastp+1)+1
  286 continue
      id1=0
      goto 287
  283 continue
      if(id1.le.0)then
      goto 521
      elseif(id1.gt.id2)then
      goto 521
      elseif(stdw(stcb,id1,id2).eq.0)then
      goto 521
      endif
      if(level.eq.1)then
      if(kin.gt.1)then
      goto 521
      endif
      if(ngmk.eq.0)then
      call aoib(1)
      llp(0)=stibs(1)
      lastkp=stibs(1)
      stib(lastkp)=eoia
      else
      call mstr1(stcb,id1,id2,llp(0),1,2,ij)
      if(ij.ne.0)then
      auxlin(1:srec)='multiply defined M-function'
      call messag(1,dbl,nline,4-xfile)
      endif
      endif
      call aoib(4)
      ii=stibs(1)-3
      call mstr0(stcb,id1,id2,udkp(0),udkl(0),ij)
      if(ij.eq.0)then
      stib(stibs(1)-2)=id1
      else
      stib(stibs(1)-2)=stib(udkp(0)+ij)
      endif
      stib(stibs(1)-1)=id2-id1+1
      stib(lastkp)=ii
      lastkp=ii
      stib(lastkp)=eoia
      ngmk=ngmk+1
      stib(stibs(1))=0
      elseif(level.eq.2)then
      if(fpass.eq.1)then
      call mstr0(stcb,id1,id2,gmkp(0),gmkl(0),ij)
      if(ij.ne.0)then
      auxlin(1:srec)='multiply defined M-function'
      call messag(1,dbl,nline,2)
      endif
      if(npmk.eq.0)then
      call aoib(1)
      stib(stibs(1))=eoia
      lastkp=stibs(1)
      llp(0)=stibs(1)
      else
      call mstr1(stcb,id1,id2,llp(0),1,2,ij)
      if(ij.ne.0)then
      auxlin(1:srec)='multiply defined M-function'
      call messag(1,dbl,nline,2)
      endif
      endif
      call mstr0(stcb,id1,id2,udkp(0),udkl(0),ij)
      call aoib(4)
      ii=stibs(1)-3
      if(ij.eq.0)then
      stib(stibs(1)-2)=id1
      stib(stibs(1)-1)=id2-id1+1
      else
      stib(stibs(1)-2)=stib(udkp(0)+ij)
      stib(stibs(1)-1)=stib(udkl(0)+ij)
      endif
      stib(stibs(1))=0
      stib(lastkp)=ii
      lastkp=ii
      stib(lastkp)=eoia
      npmk=npmk+1
      else
      call mstr0(stcb,id1,id2,pmkp(0),pmkl(0),kid)
      if(kid.eq.0)then
      auxlin(1:srec)='unexpected M-function'
      call messag(1,dbl,nline,2)
      endif
      if(stib(apmkd(0)+kid).ne.0)then
      auxlin(1:srec)='multiply defined M-function'
      call messag(1,dbl,nline,2)
      endif
      endif
      elseif(level.eq.3)then
      if(fpass.eq.1)then
      call mstr0(stcb,id1,id2,gmkp(0),gmkl(0),ij)
      if(ij.ne.0)then
      auxlin(1:srec)='multiply defined M-function'
      call messag(1,dbl,nline,2)
      endif
      call mstr0(stcb,id1,id2,pmkp(0),pmkl(0),ij)
      if(ij.ne.0)then
      auxlin(1:srec)='multiply defined M-function'
      call messag(1,dbl,nline,2)
      endif
      if(nvmk.eq.0)then
      call aoib(1)
      lastkp=stibs(1)
      stib(lastkp)=eoia
      llp(0)=stibs(1)
      else
      call mstr1(stcb,id1,id2,llp(0),1,2,ij)
      if(ij.ne.0)then
      auxlin(1:srec)='multiply defined M-function'
      call messag(1,dbl,nline,2)
      endif
      endif
      call mstr0(stcb,id1,id2,udkp(0),udkl(0),ij)
      call aoib(3)
      ii=stibs(1)-2
      if(ij.eq.0)then
      stib(stibs(1)-1)=id1
      stib(stibs(1))=id2-id1+1
      else
      stib(stibs(1)-1)=stib(udkp(0)+ij)
      stib(stibs(1))=stib(udkl(0)+ij)
      endif
      stib(lastkp)=ii
      lastkp=ii
      stib(lastkp)=eoia
      nvmk=nvmk+1
      else
      call mstr0(stcb,id1,id2,vmkp(0),vmkl(0),kid)
      if(kid.eq.0)then
      auxlin(1:srec)='unexpected M-function'
      call messag(1,dbl,nline,2)
      endif
      if(stib(avmkd(0)+kid).ne.0)then
      auxlin(1:srec)='multiply defined M-function'
      call messag(1,dbl,nline,2)
      endif
      endif
      endif
      id1=0
      goto 287
  285 continue
      if(is1.eq.0)then
      goto 521
      elseif(is1.gt.is2)then
      goto 521
      endif
      if(level.eq.1)then
      call aoib(2)
      i=lastkp+3
      if(vel.eq.0)then
      stib(i)=1
      else
      stib(i)=vin+1
      endif
      if(vel.ne.0)then
      auxlin(1:srec)='wrong dimension for M-keyword'
      call messag(1,dbl,nline,2)
      endif
      elseif(level.eq.2)then
      if(fpass.eq.1)then
      i=lastkp+3
      if(vel.eq.0)then
      stib(i)=1
      else
      stib(i)=3
      endif
      goto 184
      else
      i=apmkd(0)+kid
      if(vel.eq.0)then
      stib(i)=1
      else
      stib(i)=vin+1
      endif
      if(stib(i).gt.stib(pmkd(0)+kid))then
      auxlin(1:srec)='wrong dimension for M-function'
      call messag(1,dbl,nline,2)
      endif
      endif
      if(vin.gt.2)then
      auxlin(1:srec)='wrong dimension for M-function'
      call messag(1,dbl,nline,2)
      endif
      elseif(level.eq.3)then
      if(vel.ne.0)then
      auxlin(1:srec)='wrong dimension for M-function'
      call messag(1,dbl,nline,2)
      endif
      if(fpass.eq.1)then
      goto 184
      else
      stib(avmkd(0)+kid)=1
      endif
      endif
      if(stcb(is1:is1).eq.char(rquote))then
      kk=stds(stcb,is1,is2,1,1)
      if(kk.lt.0)then
      goto 521
      elseif(kk.gt.0)then
      stcb(is1+1:is1+kk)=stcb(stcbs(1)+1:stcbs(1)+kk)
      endif
      stcb(is1:is1)=char(172)
      stcb(is1+kk+1:is1+kk+1)=char(172)
      is1=is1+1
      i=-(is2-is1-kk)
      if(i.lt.0)then
      call cxipht(is2+1,stcbs(1),i)
      call aocb(i)
      ix=ix+i
      top=top+i
      if(level.eq.2)then
      if(nprop.eq.0)then
      do 217 j=1,npmk
      if(stib(pmkp(0)+j).gt.is2)then
      stib(pmkp(0)+j)=stib(pmkp(0)+j)+i
      endif
 217  continue
      endif
      elseif(level.eq.3)then
      if(nvert.eq.0)then
      do 218 j=1,nvmk
      if(stib(vmkp(0)+j).gt.is2)then
      stib(vmkp(0)+j)=stib(vmkp(0)+j)+i
      endif
 218  continue
      endif
      endif
      endif
      else
      if(stdw(stcb,is1,is2).eq.0)then
      if(stdq(stcb,is1,is2).eq.0)then
      goto 521
      endif
      endif
      kk=is2-is1+1
      endif
      if(level.eq.1)then
      stib(stibs(1)-1)=is1
      stib(stibs(1))=kk
      elseif(level.eq.2)then
      j=5+2*kid
      if(vel.ne.0)then
      if(vin.eq.1)then
      i=off1+j
      elseif(vin.eq.2)then
      i=off2+j
      endif
      stib(i)=is1
      stib(i+1)=kk
      else
      stib(off1+j)=is1
      stib(off1+j+1)=kk
      if(off2.ne.0)then
      stib(off2+j)=is1
      stib(off2+j+1)=kk
      endif
      endif
      elseif(level.eq.3)then
      j=lastp+stib(lastp+1)+2*kid
      stib(j)=is1
      stib(j+1)=kk
      endif
  184 continue
      is1=0
      goto 287
  777 continue
      call rvi(llp(0))
      return
  521 continue
      auxlin(1:srec)='wrong syntax or line too long'
      call messag(1,dbl,nline,4-xfile)
      end
      subroutine iki
      implicit integer(a-z)
      save
      parameter ( sibuff=1048576, scbuff=524288, mimal=2, mcmal=2 )
      parameter ( srec=81 )
      parameter ( eoia=-63 )
      character*(srec) auxlin
      common/z22g/auxlin
      common/z24g/iogp(1:4)
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z35g/stibs(0:mimal),stcbs(0:mcmal),imal(0:0),cmal(0:0)
      common/z36g/udkp(0:2),udkl(0:0),udkt(0:0),udki(0:0)
      common/z38g/gmkp(0:0),gmkl(0:0),gmkd(0:0),gmko(0:0),
     :gmkvp(0:0),gmkvl(0:0)
      common/z39g/pmkp(0:0),pmkl(0:0),pmkd(0:0),pmkvpp(0:0),pmkvlp(0:0)
      common/z40g/vmkp(0:0),vmkl(0:0),vmkvpp(0:0),vmkvlp(0:0),vmks(0:0)
      common/z41g/nudk,ngmk,npmk,nvmk,tgmkd,tpmkd
      do 100 i=1,nudk
      i1=stib(udkp(0)+i)
      i2=i1-1+stib(udkl(0)+i)
      call mstr0(stcb,i1,i2,gmkp(0),gmkl(0),jj)
      if(jj.ne.0)then
      stib(udkt(0)+i)=1
      stib(udki(0)+i)=jj
      goto 100
      endif
      call mstr0(stcb,i1,i2,pmkp(0),pmkl(0),jj)
      if(jj.ne.0)then
      stib(udkt(0)+i)=2
      stib(udki(0)+i)=jj
      goto 100
      endif
      call mstr0(stcb,i1,i2,vmkp(0),vmkl(0),jj)
      if(jj.eq.0)then
      auxlin(1:srec)='function used in style file was not found'//
     :' in model file'
      call messag(1,0,0,0)
      endif
      stib(udkt(0)+i)=3
      stib(udki(0)+i)=jj
  100 continue
      j=udkp(1)
      do 500 i1=1,nudk
      j=stib(j)
      if(stib(udkt(0)+i1).eq.1)then
      do 200 i2=3,7
      if(stib(j+i2).gt.1)then
      auxlin(1:srec)='global M-functions cannot be '//
     :'prefixed with "dual-"'
      call messag(1,0,0,0)
      endif
  200 continue
      elseif(stib(udkt(0)+i1).eq.2)then
      if(stib(pmkd(0)+stib(udki(0)+i1)).eq.3)then
      if(stib(j+3).ne.0)then
      auxlin(1:srec)='field M-functions cannot be '//
     :'used outside a loop (in style file)'
      call messag(1,0,0,0)
      elseif(stib(j+6).ne.0)then
      auxlin(1:srec)='field M-functions cannot be '//
     :'used in the plain vertex loop (in style file)'
      call messag(1,0,0,0)
      endif
      else
      if(stib(j+3).ne.0)then
      auxlin(1:srec)='propagator M-functions cannot be '//
     :'used outside a loop (in style file)'
      call messag(1,0,0,0)
      elseif(stib(j+6).ne.0)then
      auxlin(1:srec)='propagator M-functions cannot be '//
     :'used in the plain vertex loop (in style file)'
      call messag(1,0,0,0)
      endif
      do 300 i2=4,7
      if(stib(j+i2).gt.1)then
      auxlin(1:srec)='propagator M-functions cannot '//
     :'be prefixed with "dual-"'
      call messag(1,0,0,0)
      endif
  300 continue
      endif
      elseif(stib(udkt(0)+i1).eq.3)then
      if(stib(j+3).ne.0)then
      auxlin(1:srec)='vertex M-functions cannot be '//
     :'used outside a loop (in style file)'
      call messag(1,0,0,0)
      else
      do 400 i2=4,7
      if(i2.ne.5.and.stib(j+i2).gt.1)then
      auxlin(1:srec)='vertex M-functions cannot be pre'//
     :'fixed with "dual-" except in the propagator_loop'
      call messag(1,0,0,0)
      endif
  400 continue
      endif
      endif
  500 continue
      end
      subroutine rsfki
      implicit integer(a-z)
      save
      parameter ( sibuff=1048576, scbuff=524288, mimal=2, mcmal=2 )
      parameter ( srec=81 )
      parameter ( lfeed=10, dquote=34, rquote=39 )
      parameter ( eoia=-63 )
      common/z24g/iogp(1:4)
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z35g/stibs(0:mimal),stcbs(0:mcmal),imal(0:0),cmal(0:0)
      common/z36g/udkp(0:2),udkl(0:0),udkt(0:0),udki(0:0)
      common/z38g/gmkp(0:0),gmkl(0:0),gmkd(0:0),gmko(0:0),
     :gmkvp(0:0),gmkvl(0:0)
      common/z39g/pmkp(0:0),pmkl(0:0),pmkd(0:0),pmkvpp(0:0),pmkvlp(0:0)
      common/z40g/vmkp(0:0),vmkl(0:0),vmkvpp(0:0),vmkvlp(0:0),vmks(0:0)
      common/z41g/nudk,ngmk,npmk,nvmk,tgmkd,tpmkd
      kstep=8
      udkp(2)=stibs(1)
      i=2*(nudk+1)
      call aoib(i)
      udkp(0)=udkp(2)+nudk+1
      j=udkp(1)
      do 10 i=1,nudk
      j=stib(j)
      stib(udkp(2)+i)=j+1
      stib(udkp(0)+i)=stib(j+1)
   10 continue
      stib(udkp(0))=eoia
      stib(stibs(1))=eoia
      udkl(0)=stibs(1)
      i=nudk+1
      call aoib(i)
      j=udkp(1)
      do 20 i=1,nudk
      j=stib(j)
      stib(udkl(0)+i)=stib(j+2)
   20 continue
      stib(stibs(1))=eoia
      udkt(0)=stibs(1)
      i=2*(nudk+1)
      call aoib(i)
      i=nudk+1
      udki(0)=udkt(0)+i
      do 30 i=1,nudk
      stib(udkt(0)+i)=0
      stib(udki(0)+i)=0
   30 continue
      stib(udki(0))=eoia
      stib(stibs(1))=eoia
      end
      subroutine rgmki(llpp)
      implicit integer(a-z)
      save
      parameter ( sibuff=1048576, scbuff=524288, mimal=2, mcmal=2 )
      parameter ( srec=81 )
      parameter ( eoia=-63 )
      common/z30g/stib(1:sibuff)
      common/z35g/stibs(0:mimal),stcbs(0:mcmal),imal(0:0),cmal(0:0)
      common/z38g/gmkp(0:0),gmkl(0:0),gmkd(0:0),gmko(0:0),
     :gmkvp(0:0),gmkvl(0:0)
      common/z39g/pmkp(0:0),pmkl(0:0),pmkd(0:0),pmkvpp(0:0),pmkvlp(0:0)
      common/z40g/vmkp(0:0),vmkl(0:0),vmkvpp(0:0),vmkvlp(0:0),vmks(0:0)
      common/z41g/nudk,ngmk,npmk,nvmk,tgmkd,tpmkd
      ii=ngmk+1
      gmkp(0)=stibs(1)
      call aoib(4*ii)
      gmkl(0)=gmkp(0)+ii
      gmkd(0)=gmkl(0)+ii
      gmko(0)=gmkd(0)+ii
      stib(gmkl(0))=eoia
      stib(gmkd(0))=eoia
      stib(gmko(0))=eoia
      stib(stibs(1))=eoia
      i=0
      tgmkd=0
      jj=llpp
   10 continue
      jj=stib(jj)
      if(jj.ne.eoia)then
      i=i+1
      stib(gmkp(0)+i)=stib(jj+1)
      stib(gmkl(0)+i)=stib(jj+2)
      stib(gmkd(0)+i)=stib(jj+3)
      stib(gmko(0)+i)=tgmkd
      tgmkd=tgmkd+max(1,stib(jj+3)-1)
      goto 10
      endif
      ii=tgmkd+1
      gmkvp(0)=stibs(1)
      jj=ii+ii
      call aoib(jj)
      gmkvl(0)=gmkvp(0)+ii
      stib(gmkvl(0))=eoia
      stib(stibs(1))=eoia
      i=0
      j=0
      jj=llpp
   20 continue
      jj=stib(jj)
      if(jj.ne.eoia)then
      i=i+1
      do 30 k=1,max(1,stib(gmkd(0)+i)-1)
      j=j+1
      stib(gmkvp(0)+j)=stib(jj+2+2*k)
      stib(gmkvl(0)+j)=stib(jj+3+2*k)
   30 continue
      goto 20
      endif
      ii=gmkp(0)-llpp+1
      call xipht(gmkp(0)+1,stibs(1),-ii)
      call aoib(-ii)
      gmkp(0)=gmkp(0)-ii
      gmkl(0)=gmkl(0)-ii
      gmkd(0)=gmkd(0)-ii
      gmko(0)=gmko(0)-ii
      gmkvp(0)=gmkvp(0)-ii
      gmkvl(0)=gmkvl(0)-ii
      end
      subroutine rpi(nphi)
      implicit integer(a-z)
      save
      parameter ( sibuff=1048576, scbuff=524288, mimal=2, mcmal=2 )
      parameter ( srec=81 )
      parameter ( lfeed=10, dquote=34, rquote=39 )
      parameter ( eoia=-63 )
      common/z9g/tpc(0:0)
      common/z11g/npart,nblok,nprop
      character*(srec) auxlin
      common/z22g/auxlin
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z33g/namep(0:0),namel(0:0)
      common/z34g/link(0:0),antiq(0:0),blok(0:0)
      common/z35g/stibs(0:mimal),stcbs(0:mcmal),imal(0:0),cmal(0:0)
      common/z38g/gmkp(0:0),gmkl(0:0),gmkd(0:0),gmko(0:0),
     :gmkvp(0:0),gmkvl(0:0)
      common/z39g/pmkp(0:0),pmkl(0:0),pmkd(0:0),pmkvpp(0:0),pmkvlp(0:0)
      common/z40g/vmkp(0:0),vmkl(0:0),vmkvpp(0:0),vmkvlp(0:0),vmks(0:0)
      common/z41g/nudk,ngmk,npmk,nvmk,tgmkd,tpmkd
      common/z42g/pmkr(0:0),pmkvma(0:0),pmkvmi(0:0)
      common/z43g/vmkr(0:0),vmkmao(0:0),vmkmio(0:0)
      if(nphi.lt.1)then
      auxlin(1:srec)='no propagators were found'
      call messag(1,0,0,0)
      endif
      psize=6+2*npmk
      i=stibs(1)
      call aoib(psize)
      do 424 j=i+1,stibs(1)
      stib(j)=eoia
  424 continue
      call trm(psize,nphi+1)
      ii=nphi+1
      blok(0)=stibs(1)-psize*ii
      link(0)=blok(0)+ii
      antiq(0)=link(0)+ii
      tpc(0)=antiq(0)+ii
      namep(0)=tpc(0)+ii
      namel(0)=namep(0)+ii
      do 11 i=1,nphi
      stib(link(0)+i)=stib(link(0)+i)+i
   11 continue
      stib(blok(0))=eoia
      pmkvpp(0)=stibs(1)
      call aoib(npmk+npmk+2)
      pmkvlp(0)=pmkvpp(0)+npmk+1
      stib(pmkvlp(0))=eoia
      stib(stibs(1))=eoia
      jj=namel(0)
      do 22 i=1,npmk
      jj=jj+ii
      stib(pmkvpp(0)+i)=jj
      jj=jj+ii
      stib(pmkvlp(0)+i)=jj
   22 continue
      pmkr(0)=stibs(1)
      call aoib(npmk+1)
      stib(stibs(1))=eoia
      do 81 i=1,npmk
      ii=stib(pmkvpp(0)+i)
      jj=stib(pmkvlp(0)+i)
      kk=4
      do 80 j=1,nphi
      i1=stib(ii+j)
      i2=stib(ii+j)-1+stib(jj+j)
      if(i1.gt.i2)then
      kk=1
      endif
      if(kk.eq.4)then
      if(stdz(stcb,i1,i2).eq.1)then
      goto 80
      else
      kk=3
      endif
      endif
      if(kk.eq.3)then
      if(stdq(stcb,i1,i2).eq.1)then
      goto 80
      elseif(j.eq.1)then
      kk=2
      else
      kk=1
      endif
      endif
      if(kk.eq.2)then
      if(stdw(stcb,i1,i2).eq.1)then
      goto 80
      else
      kk=1
      endif
      endif
   80 continue
      stib(pmkr(0)+i)=kk
   81 continue
      pmkvmi(0)=stibs(1)
      call aoib(npmk+1)
      stib(stibs(1))=eoia
      pmkvma(0)=stibs(1)
      call aoib(npmk+1)
      stib(stibs(1))=eoia
      do 910 i=1,npmk
      if(stib(pmkr(0)+i).eq.4)then
      ii=stib(pmkvpp(0)+i)
      jj=stib(pmkvlp(0)+i)
      do 901 j=1,npart
      i1=stib(ii+j)
      i2=stib(ii+j)-1+stib(jj+j)
      call stoz(stcb,i1,i2,ij)
      if(j.gt.1)then
      if(ij.lt.k)then
      k=ij
      elseif(ij.gt.kk)then
      kk=ij
      endif
      else
      k=ij
      kk=ij
      endif
  901 continue
      else
      k=0
      kk=0
      endif
      stib(pmkvmi(0)+i)=k
      stib(pmkvma(0)+i)=kk
  910 continue
      end
      subroutine rvi(llp)
      implicit integer(a-z)
      save
      parameter ( maxleg=10, maxrho=6, maxi=10 )
      parameter ( maxdeg=6 )
      parameter ( maxn=maxi+maxi-2, maxli=maxn+maxrho-1 )
      parameter ( sibuff=1048576, scbuff=524288, mimal=2, mcmal=2 )
      parameter ( srec=81 )
      parameter ( lfeed=10, dquote=34, rquote=39 )
      parameter ( eoia=-63 )
      common/z2in/nv(1:maxdeg)
      common/z4in/mflag(1:21)
      common/z1g/nrho,rho(1:maxdeg),g(1:maxn,1:maxn)
      common/z3g/vparto(0:0),vdeg(0:0),nvert
      common/z9g/tpc(0:0)
      common/z11g/npart,nblok,nprop
      common/z15g/dpntro(1:maxdeg)
      character*(srec) auxlin
      common/z22g/auxlin
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z33g/namep(0:0),namel(0:0)
      common/z34g/link(0:0),antiq(0:0),blok(0:0)
      common/z35g/stibs(0:mimal),stcbs(0:mcmal),imal(0:0),cmal(0:0)
      common/z38g/gmkp(0:0),gmkl(0:0),gmkd(0:0),gmko(0:0),
     :gmkvp(0:0),gmkvl(0:0)
      common/z39g/pmkp(0:0),pmkl(0:0),pmkd(0:0),pmkvpp(0:0),pmkvlp(0:0)
      common/z40g/vmkp(0:0),vmkl(0:0),vmkvpp(0:0),vmkvlp(0:0),vmks(0:0)
      common/z41g/nudk,ngmk,npmk,nvmk,tgmkd,tpmkd
      common/z42g/pmkr(0:0),pmkvma(0:0),pmkvmi(0:0)
      common/z43g/vmkr(0:0),vmkmao(0:0),vmkmio(0:0)
      integer xtmp(1:maxdeg+1),nrot(1:maxdeg),rotvpo(1:maxdeg)
      if(nvert.eq.0)then
      auxlin(1:srec)='no vertices were found'
      call messag(1,0,0,0)
      endif
      ii=nvert+1
      jj=nvmk+1
      i=2*jj*(ii+1)
      j=stibs(1)
      if(stib(stibs(1)).ne.eoia)then
      i=i+1
      j=j+1
      endif
      call aoib(i)
      do 05 i=1,maxdeg
      nv(i)=0
      nrot(i)=0
   05 continue
      vdeg(0)=j
      vparto(0)=vdeg(0)+ii
      vmkvpp(0)=vparto(0)+ii
      vmkvlp(0)=vmkvpp(0)+jj
      j=vmkvlp(0)+jj
      do 07 i=1,nvmk
      stib(vmkvpp(0)+i)=j
      j=j+ii
      stib(vmkvlp(0)+i)=j
      j=j+ii
   07 continue
      i=0
      jj=llp
   10 continue
      jj=stib(jj)
      if(jj.ne.eoia)then
      i=i+1
      k=stib(jj+1)
      stib(vdeg(0)+i)=k
      nv(k)=nv(k)+1
      goto 10
      endif
      if(i.ne.nvert.or.j.ne.stibs(1))then
      auxlin(1:srec)='internal run-time error'
      call messag(1,0,0,0)
      endif
      nrho=0
      do 17 i=maxdeg,3,-1
      if(nrho.eq.0.and.nv(i).gt.0)then
      nrho=i
      endif
   17 continue
      if(nrho.lt.3)then
      auxlin(1:srec)='internal run-time error'
      call messag(1,0,0,0)
      endif
      i=0
      jj=llp
   20 continue
      jj=stib(jj)
      if(jj.ne.eoia)then
      i=i+1
      do 30 i1=1,nvmk
      i2=jj+stib(vdeg(0)+i)+i1+i1
      stib(stib(vmkvpp(0)+i1)+i)=stib(i2)
      stib(stib(vmkvlp(0)+i1)+i)=stib(i2+1)
   30 continue
      goto 20
      endif
      if(i.ne.nvert)then
      auxlin(1:srec)='internal run-time error'
      call messag(1,0,0,0)
      endif
      stib(vdeg(0))=eoia
      stib(vparto(0))=eoia
      stib(vmkvpp(0))=eoia
      stib(vmkvlp(0))=eoia
      j=vmkvlp(0)+nvmk+1
      stib(j)=eoia
      do 40 i1=1,nvmk+nvmk
      j=j+ii
      stib(j)=eoia
   40 continue
      if(j.ne.stibs(1))then
      auxlin(1:srec)='internal run-time error'
      call messag(1,0,0,0)
      endif
      kk=llp
      if(llp.gt.1)then
      if(stib(llp-1).eq.eoia)then
      kk=llp-1
      endif
      endif
      ij=stib(llp)
      if(stib(llp).ne.eoia)then
      stib(llp)=eoia
      endif
      i=0
   50 continue
      jj=ij
      if(jj.ne.eoia)then
      i=i+1
      ij=stib(ij)
      k=jj-kk+1
      call xipht(jj+2,jj+1+stib(vdeg(0)+i),-k)
      stib(vparto(0)+i)=kk
      kk=kk+stib(vdeg(0)+i)
      goto 50
      endif
      if(i.ne.nvert)then
      auxlin(1:srec)='internal run-time error'
      call messag(1,0,0,0)
      endif
      k=vdeg(0)-kk-1
      call xipht(vdeg(0),stibs(1),-k)
      call aoib(-k)
      vdeg(0)=vdeg(0)-k
      vparto(0)=vparto(0)-k
      vmkvpp(0)=vmkvpp(0)-k
      vmkvlp(0)=vmkvlp(0)-k
      do 60 i=1,nvmk
      stib(vmkvpp(0)+i)=stib(vmkvpp(0)+i)-k
      stib(vmkvlp(0)+i)=stib(vmkvlp(0)+i)-k
   60 continue
      vmkr(0)=stibs(1)
      call aoib(nvmk+1)
      stib(stibs(1))=eoia
      ij=1
      do 91 i=1,nvmk
      ii=stib(vmkvpp(0)+i)
      jj=stib(vmkvlp(0)+i)
      kk=4
      do 90 j=1,nvert
      i1=stib(ii+j)
      i2=stib(ii+j)-1+stib(jj+j)
      if(i1.gt.i2)then
      kk=1
      endif
      if(kk.eq.4)then
      if(stdz(stcb,i1,i2).eq.1)then
      goto 90
      else
      kk=3
      endif
      endif
      if(kk.eq.3)then
      if(stdq(stcb,i1,i2).eq.1)then
      goto 90
      elseif(j.eq.1)then
      kk=2
      else
      kk=1
      endif
      endif
      if(kk.eq.2)then
      if(stdw(stcb,i1,i2).eq.1)then
      goto 90
      else
      kk=1
      endif
      endif
   90 continue
      stib(vmkr(0)+i)=kk
      if(kk.ne.4)then
      ij=0
      endif
   91 continue
      vmks(0)=stibs(1)
      ij=stibs(1)
      call aoib(nvmk+1)
      do 981 i=ij+1,ij+nvmk
      stib(i)=0
  981 continue
      stib(stibs(1))=eoia
      j=0
      do 388 i=1,nrho
      if(nv(i).gt.0)then
      dpntro(i)=stibs(1)
      call aoib(npart+1)
      stib(stibs(1))=eoia
      else
      dpntro(i)=-1
      endif
  388 continue
      nvrot=0
      do 403 ideg=1,nrho
      if(nv(ideg).gt.0)then
      rotvpo(ideg)=stibs(1)
      do 401 i1=1,nvert
      if(ideg.eq.stib(vdeg(0)+i1))then
      do 166 i2=1,ideg
      xtmp(i2)=stib(stib(vparto(0)+i1)+i2)
  166 continue
      do 186 i2=1,ideg-1
      do 176 i3=i2+1,ideg
      if(xtmp(i2).gt.xtmp(i3))then
      aux=xtmp(i2)
      xtmp(i2)=xtmp(i3)
      xtmp(i3)=aux
      endif
  176 continue
  186 continue
      goto 830
  805 continue
      do 196 i=ideg-1,1,-1
      if(xtmp(i).lt.xtmp(i+1))then
      goto 810
      endif
  196 continue
      goto 401
  810 continue
      do 206 ii=i+2,ideg
      if(xtmp(ii).le.xtmp(i))then
      goto 820
      endif
  206 continue
      ii=ideg+1
  820 continue
      aux=xtmp(i)
      xtmp(i)=xtmp(ii-1)
      xtmp(ii-1)=aux
      aux=ideg-i
      aux=aux/2
      i=i+1
      do 216 ii=0,aux-1
      aux=xtmp(i+ii)
      xtmp(i+ii)=xtmp(ideg-ii)
      xtmp(ideg-ii)=aux
  216 continue
  830 continue
      j=stibs(1)+1
      call aoib(ideg+1)
      stib(j)=i1
      do 256 i2=1,ideg
      stib(j+i2)=xtmp(i2)
  256 continue
      nvrot=nvrot+1
      nrot(ideg)=nrot(ideg)+1
      goto 805
      endif
  401 continue
      j=stibs(1)+1
      call aoib(ideg+1)
      stib(j)=0
      do 257 i2=1,ideg
      stib(j+i2)=npart+1
  257 continue
      endif
  403 continue
      ij=0
      do 790 ideg=1,nrho
      if(nv(ideg).gt.0.and.nrot(ideg).gt.1)then
      h1=1+nrot(ideg)/2
      h2=nrot(ideg)
  720 continue
      if(h1.gt.1)then
      h1=h1-1
      ii=rotvpo(ideg)+(ideg+1)*(h1-1)
      do 705 i=1,ideg+1
      xtmp(i)=stib(ii+i)
  705 continue
      else
      ii=rotvpo(ideg)+(ideg+1)*(h2-1)
      do 715 i=1,ideg+1
      xtmp(i)=stib(ii+i)
  715 continue
      ii=rotvpo(ideg)+(ideg+1)*(h2-1)
      do 725 i=1,ideg+1
      stib(ii+i)=stib(rotvpo(ideg)+i)
  725 continue
      h2=h2-1
      if(h2.eq.1)then
      do 735 i=1,ideg+1
      stib(rotvpo(ideg)+i)=xtmp(i)
  735 continue
      goto 790
      endif
      endif
      hj=h1
  740 continue
      hi=hj
      hj=hi+hi
      if(hj.le.h2)then
      if(hj.lt.h2)then
      kk=0
      jj=rotvpo(ideg)+(ideg+1)*(hj-1)
      jk=jj+ideg+1
      i=1
  745 continue
      if(i.le.ideg)then
      i=i+1
      if(stib(jj+i).lt.stib(jk+i))then
      kk=-1
      elseif(stib(jj+i).gt.stib(jk+i))then
      kk=1
      else
      goto 745
      endif
      endif
      if(kk.eq.-1)then
      hj=hj+1
      elseif(kk.eq.0)then
      ij=1
      endif
      endif
      kk=0
      jj=rotvpo(ideg)+(ideg+1)*(hj-1)
      i=1
  755 continue
      if(i.le.ideg)then
      i=i+1
      if(xtmp(i).lt.stib(jj+i))then
      kk=-1
      elseif(xtmp(i).gt.stib(jj+i))then
      kk=1
      else
      goto 755
      endif
      endif
      if(kk.eq.-1)then
      ii=rotvpo(ideg)+(ideg+1)*(hi-1)
      do 765 i=1,ideg+1
      stib(ii+i)=stib(jj+i)
  765 continue
      goto 740
      elseif(kk.eq.0)then
      ij=1
      endif
      endif
      ii=rotvpo(ideg)+(ideg+1)*(hi-1)
      do 775 i=1,ideg+1
      stib(ii+i)=xtmp(i)
  775 continue
      goto 720
      endif
  790 continue
      if(ij.ne.0)then
      auxlin(1:srec)='repeated vertices found'
      call messag(0,0,0,0)
      endif
      do 870 ideg=1,nrho
      if(nrot(ideg).gt.0)then
      ii=rotvpo(ideg)
      do 850 j=1,nrot(ideg)
      k=0
      kk=0
      i1=ii+2
      i2=i1+ideg+1
  840 continue
      if(k.lt.ideg)then
      if(stib(i1).lt.stib(i2))then
      kk=-1
      elseif(stib(i1).gt.stib(i2))then
      kk=1
      else
      k=k+1
      i1=i1+1
      i2=i2+1
      goto 840
      endif
      endif
      ii=ii+ideg+1
  850 continue
      endif
  870 continue
      do 970 ideg=1,nrho
      if(nrot(ideg).gt.0)then
      ii=rotvpo(ideg)+2
      jj=0
      do 950 j=1,nrot(ideg)+1
      if(stib(ii).ne.jj)then
      if(stib(ii).le.npart)then
      kk=stib(ii)
      else
      kk=npart
      endif
      do 930 k=jj+1,kk
      stib(dpntro(ideg)+k)=ii-1
  930 continue
      jj=stib(ii)
      endif
      ii=ii+ideg+1
  950 continue
      endif
  970 continue
      call aoib(1)
      stib(stibs(1))=eoia
      bloka=stibs(1)
      call aoib(npart)
      do 304 i1=1,npart
      stib(blok(0)+i1)=0
  304 continue
      nblok=0
      k=0
      jj=0
  310 continue
      if(k.lt.npart)then
      nblok=nblok+1
      i=nblok
  320 continue
      if(stib(blok(0)+i).ne.0)then
      i=i+1
      goto 320
      endif
      k=k+1
      jj=k
      stib(bloka+k)=i
      stib(blok(0)+i)=nblok
      if(stib(blok(0)+stib(link(0)+i)).eq.0)then
      k=k+1
      stib(bloka+k)=stib(link(0)+i)
      stib(blok(0)+stib(bloka+k))=nblok
      endif
  330 continue
      if(jj.le.k)then
      do 376 i=1,nrho
      if(dpntro(i).gt.0)then
      j=stib(dpntro(i)+stib(bloka+jj))
  366 continue
      if(stib(j+1).eq.stib(bloka+jj))then
      aux=stib(j+2)
      if(stib(blok(0)+aux).eq.0)then
      k=k+1
      stib(bloka+k)=aux
      stib(blok(0)+aux)=nblok
      endif
      if(stib(blok(0)+stib(link(0)+aux)).eq.0)then
      k=k+1
      stib(bloka+k)=stib(link(0)+aux)
      stib(blok(0)+stib(bloka+k))=nblok
      endif
      j=j+i+1
      goto 366
      endif
      endif
  376 continue
      jj=jj+1
      goto 330
      endif
      goto 310
      endif
      if(nblok.gt.1)then
      auxlin(1:srec)='model splits into disjoint components'
      call messag(0,0,0,0)
      endif
      i=-npart
      call aoib(i)
      j=0
      do 386 i=1,npart
      do 381 k=1,nrho
      if(nv(k).gt.0)then
      if(stib(stib(dpntro(k)+i)+1).eq.i)then
      goto 340
      endif
      endif
  381 continue
      j=j+1
  340 continue
  386 continue
      if(j.ne.0)then
      auxlin(1:srec)='model contains at least one non-interacting'//
     :' field'
      call messag(0,0,0,0)
      endif
      end
      subroutine vsig
      implicit integer(a-z)
      save
      parameter ( maxleg=10, maxrho=6, maxi=10 )
      parameter ( maxn=maxi+maxi-2, maxli=maxn+maxrho-1 )
      parameter ( maxdeg=6 )
      parameter ( sibuff=1048576, scbuff=524288, mimal=2, mcmal=2 )
      parameter ( srec=81 )
      parameter ( eoia=-63 )
      common/z2in/nv(1:maxdeg)
      common/z4in/mflag(1:21)
      common/z1g/nrho,rho(1:maxdeg),g(1:maxn,1:maxn)
      common/z3g/vparto(0:0),vdeg(0:0),nvert
      common/z9g/tpc(0:0)
      common/z11g/npart,nblok,nprop
      common/z15g/dpntro(1:maxdeg)
      character*(srec) auxlin
      common/z22g/auxlin
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z35g/stibs(0:mimal),stcbs(0:mcmal),imal(0:0),cmal(0:0)
      common/z40g/vmkp(0:0),vmkl(0:0),vmkvpp(0:0),vmkvlp(0:0),vmks(0:0)
      common/z41g/nudk,ngmk,npmk,nvmk,tgmkd,tpmkd
      common/z43g/vmkr(0:0),vmkmao(0:0),vmkmio(0:0)
      if(mflag(21).eq.0)then
      vmkmio(0)=0
      vmkmao(0)=0
      goto 981
      endif
      vmkmio(0)=stibs(1)
      call aoib(nvmk+1)
      stib(stibs(1))=eoia
      vmkmao(0)=stibs(1)
      call aoib(nvmk+1)
      stib(stibs(1))=eoia
      k0=0
      do 910 i=1,nvmk
      if(stib(vmks(0)+i).ne.1)then
      stib(vmkmio(0)+i)=0
      stib(vmkmao(0)+i)=0
      else
      if(k0.eq.0)then
      k0=1
      k1=0
      k2=nrho
      else
      k1=-nrho
      k2=0
      endif
      stib(vmkmio(0)+i)=stibs(1)+k1
      stib(vmkmao(0)+i)=stibs(1)+k2
      call aoib(2*nrho+k2)
      k2=stibs(1)-nrho
      do 900 j=k2+1,stibs(1)
      stib(j)=0
  900 continue
      ii=stib(vmkvpp(0)+i)
      jj=stib(vmkvlp(0)+i)
      do 901 j=1,nvert
      i1=stib(ii+j)
      i2=stib(ii+j)-1+stib(jj+j)
      call stoz(stcb,i1,i2,ij)
      i3=stib(vdeg(0)+j)
      k1=k2+i3
      if(stib(k1).ne.0)then
      i1=stib(vmkmao(0)+i)+i3
      if(ij.gt.stib(i1))then
      stib(i1)=ij
      else
      i1=stib(vmkmio(0)+i)+i3
      if(ij.lt.stib(i1))then
      stib(i1)=ij
      endif
      endif
      else
      stib(k1)=1
      i1=k1-nrho
      stib(i1)=ij
      stib(i1-nrho)=ij
      endif
  901 continue
      endif
  910 continue
      if(k0.ne.0)then
      k0=-nrho
      endif
      call aoib(k0+1)
      stib(stibs(1))=eoia
  981 continue
      end
      subroutine gen10(cntr)
      implicit integer(a-z)
      save
      parameter ( maxleg=10, maxrho=6, maxi=10 )
      parameter ( maxdeg=6 )
      parameter ( maxn=maxi+maxi-2, maxli=maxn+maxrho-1 )
      parameter ( ipar1= maxleg*maxleg-maxleg )
      parameter ( ipar2= ipar1/2 )
      parameter ( srec=81 )
      parameter ( lfeed=10, dquote=34, rquote=39 )
      common/z1in/leg(1:maxleg),nleg,nloop,incom
      common/z4in/mflag(1:21)
      common/z1g/nrho,rho(1:maxdeg),g(1:maxn,1:maxn)
      common/z4g/n,nli
      common/z6g/p1(1:maxleg),invp1(1:maxleg)
      common/z7g/lmap(1:maxn,1:maxdeg),vmap(1:maxn,1:maxdeg),
     :pmap(1:maxn,1:maxdeg),vlis(1:maxn),invlis(1:maxn)
      common/z8g/degree(1:maxn)
      common/z10g/p1l(1:ipar2),p1r(1:ipar2),ns1
      common/z13g/xn(1:maxn)
      common/z16g/rdeg(1:maxn),amap(1:maxn,1:maxdeg)
      common/z18g/eg(1:maxn,1:maxn),flow(1:maxli,0:maxleg+maxrho)
      character*(srec) auxlin
      common/z22g/auxlin
      integer vaux(1:maxn)
      if(cntr.eq.0)then
      if(mflag(8).eq.0)then
      goto 29
      else
      goto 27
      endif
      endif
      cntr22=1
   27 continue
      call gen22(cntr22)
      if(cntr22.eq.0)then
      cntr=0
      return
      endif
      do 06 i=1,n
      vaux(i)=xn(i)
   06 continue
      do 16 i=1,rho(1)
      vlis(i)=i
      invlis(i)=i
   16 continue
      do 26 i=rho(1)+1,n
      invlis(i)=0
   26 continue
      aux=rho(1)
      jj=rho(1)+1
   10 continue
      if(aux.lt.n)then
   20 continue
      if(invlis(jj).ne.0)then
      jj=jj+1
      goto 20
      endif
      ii=jj
      do 36 i=ii+1,n
      if(invlis(i).eq.0.and.vaux(i).gt.vaux(ii))then
      ii=i
      endif
   36 continue
      if(vaux(ii).eq.0)then
      auxlin(1:srec)='gen10_1'
      call messag(1,0,0,0)
      endif
      aux=aux+1
      vlis(aux)=ii
      invlis(ii)=aux
      rdeg(ii)=vaux(ii)
      do 46 i=rho(1)+1,n
      vaux(i)=vaux(i)+g(i,ii)
   46 continue
      goto 10
      endif
      do 56 i=1,n
      vaux(i)=0
   56 continue
      do 76 i=1,n
      ii=vlis(i)
      kk=0
      j=1
   30 continue
      if(kk.lt.degree(ii))then
      jj=vlis(j)
      aux=1
      do 66 k=1,g(ii,jj)
      kk=kk+1
      vmap(ii,kk)=jj
      vaux(jj)=vaux(jj)+1
      if(ii.ne.jj)then
      lmap(ii,kk)=vaux(jj)
      else
      lmap(ii,kk)=vaux(jj)+aux
      aux=-aux
      endif
   66 continue
      j=j+1
      goto 30
      endif
      if(kk.ne.degree(ii))then
      auxlin(1:srec)='gen10_2'
      call messag(1,0,0,0)
      endif
   76 continue
      do 96 i=rho(1)+1,n
      do 86 j=1,degree(i)-1
      if(invlis(vmap(i,j)).gt.invlis(vmap(i,j+1)))then
      auxlin(1:srec)='gen10_3'
      call messag(1,0,0,0)
      endif
   86 continue
   96 continue
      do 116 i=1,n
      do 106 j=1,degree(i)
      if(vmap(vmap(i,j),lmap(i,j)).ne.i)then
      auxlin(1:srec)='gen10_4'
      call messag(1,0,0,0)
      endif
      if(lmap(vmap(i,j),lmap(i,j)).ne.j)then
      auxlin(1:srec)='gen10_5'
      call messag(1,0,0,0)
      endif
  106 continue
  116 continue
      if(mflag(16)+mflag(17).ne.0)then
      do 122 i=rho(1)+1,n
      do 118 j=1,degree(i)
      if(vmap(i,j).gt.rho(1))then
      ii=min(i,vmap(i,j))
      jj=max(i,vmap(i,j))
      amap(i,j)=eg(ii,jj)
      if(g(i,vmap(i,j)).gt.1)then
      if(i.eq.vmap(i,j))then
      amap(i,j)=amap(i,j)+(j-1-rdeg(i))/2
      else
      ij=0
      do 117 k=j-1,1,-1
      if(vmap(i,k).eq.vmap(i,j))then
      ij=ij+1
      else
      goto 123
      endif
  117 continue
  123 continue
      amap(i,j)=eg(ii,jj)+ij
      endif
      endif
      endif
  118 continue
  122 continue
      endif
      do 84 i=1,rho(1)
      p1(i)=i
   84 continue
      goto 121
   29 continue
      do 17 i=rho(1)-1,1,-1
      if(p1(i).lt.p1(i+1))then
      goto 38
      endif
   17 continue
      goto 27
   38 continue
      do 75 ii=i+2,rho(1)
      if(p1(ii).lt.p1(i))then
      goto 61
      endif
   75 continue
      ii=rho(1)+1
   61 continue
      aux=p1(i)
      p1(i)=p1(ii-1)
      p1(ii-1)=aux
      aux=rho(1)-i
      aux=aux/2
      i=i+1
      do 33 ii=0,aux-1
      aux=p1(i+ii)
      p1(i+ii)=p1(rho(1)-ii)
      p1(rho(1)-ii)=aux
   33 continue
  121 continue
      do 124 i=1,ns1
      aux=p1(p1r(i))-p1(p1l(i))
      if(aux.lt.0)then
      goto 29
      endif
  124 continue
      do 136 i=1,rho(1)
      invp1(p1(i))=i
  136 continue
      cntr=1
      if(mflag(16).eq.0)then
      return
      endif
      do 114 i=1,rho(1)
      amap(i,1)=i
      amap(vmap(i,1),lmap(i,1))=amap(i,1)
  114 continue
      do 222 i=rho(1)+1,n
      do 225 j=1,degree(i)
      if(amap(i,j).gt.nli)then
      auxlin(1:srec)='gen10_6'
      call messag(1,0,0,0)
      endif
  225 continue
  222 continue
      do 228 ii=1,nli
      s=0
      do 231 i=1,n
      do 234 j=1,degree(i)
      if(amap(i,j).eq.ii)then
      s=s+1
      endif
  234 continue
  231 continue
      if(s.ne.2)then
      auxlin(1:srec)='gen10_7'
      call messag(1,0,0,0)
      endif
  228 continue
      end
      subroutine gen22(cntr22)
      implicit integer(a-z)
      save
      parameter ( maxleg=10, maxrho=6, maxi=10 )
      parameter ( maxdeg=6 )
      parameter ( maxn=maxi+maxi-2, maxli=maxn+maxrho-1 )
      parameter ( ipar1= maxleg*maxleg-maxleg )
      parameter ( ipar2= ipar1/2 )
      parameter ( sibuff=1048576, scbuff=524288, mimal=2, mcmal=2 )
      parameter ( srec=81 )
      parameter ( lfeed=10, dquote=34, rquote=39 )
      parameter ( eoia=-63 )
      common/z4in/mflag(1:21)
      common/z1g/nrho,rho(1:maxdeg),g(1:maxn,1:maxn)
      common/z4g/n,nli
      common/z5g/psym(0:0),psyms,nsym
      common/z7g/lmap(1:maxn,1:maxdeg),vmap(1:maxn,1:maxdeg),
     :pmap(1:maxn,1:maxdeg),vlis(1:maxn),invlis(1:maxn)
      common/z8g/degree(1:maxn)
      common/z10g/p1l(1:ipar2),p1r(1:ipar2),ns1
      common/z13g/xn(1:maxn)
      common/z14g/zcho(0:maxli),zbri(0:maxli),zpro(0:maxli),
     :rbri(0:maxli),sbri(0:maxli)
      common/z17g/xtail(1:maxn),xhead(1:maxn),ntadp
      common/z18g/eg(1:maxn,1:maxn),flow(1:maxli,0:maxleg+maxrho)
      character*(srec) auxlin
      common/z22g/auxlin
      common/z30g/stib(1:sibuff)
      common/z35g/stibs(0:mimal),stcbs(0:mcmal),imal(0:0),cmal(0:0)
      integer xc(1:maxdeg)
      integer xl(1:maxn),xt(1:maxn)
      integer bound(1:maxdeg),degr(1:maxdeg),xs(1:maxn)
      integer xg(1:maxn,1:maxn),ds(1:maxn,1:maxn)
      integer p1sym(1:maxleg,1:maxleg)
      integer lista(1:maxn),stack(1:maxn),listb(1:maxn)
      integer uset(0:maxn),xset(1:maxn),xp(1:maxn),a1(1:maxn)
      integer str(1:maxn),dta(1:maxn),lps(1:maxn)
      integer head(1:maxli),tail(1:maxli),intree(1:maxli)
      integer emul(1:maxdeg,1:maxli),nemul(1:maxdeg)
      if(cntr22.ne.1.and.cntr22.ne.-1)then
      goto 937
      endif
      if(nrho.lt.1.or.nrho.gt.nrho)then
      goto 937
      endif
      do 04 i=1,nrho
      if(rho(i).lt.0.or.rho(i).gt.maxn)then
      goto 937
      endif
   04 continue
      if(cntr22.eq.-1)then
      goto 05
      endif
      n=0
      nli=0
      do 16 i=1,nrho
      nli=nli+i*rho(i)
      do 06 j=1,rho(i)
      n=n+1
      degree(n)=i
   06 continue
   16 continue
      if(n.lt.1.or.n.gt.maxn)then
      goto 937
      endif
      if(mod(nli,2).ne.0)then
      auxlin(1:srec)='gen22_1'
      call messag(1,0,0,0)
      endif
      nli=nli/2
      ni=n-rho(1)
      loop=nli-n+1
      if(loop.lt.0)then
      goto 394
      endif
      do 26 i=rho(1)+1,n
      if(mflag(6).eq.1)then
      str(i)=1
      elseif(degree(i).ne.degree(n))then
      str(i)=min(degree(i),loop+1)
      elseif(n.gt.2)then
      str(i)=min(degree(i)-1,loop+1)
      else
      str(i)=min(degree(i),loop+1)
      endif
   26 continue
      j=1
      xl(1)=1
      do 46 i=2,n
      if(degree(i).ne.degree(i-1))then
      xt(j)=i-1
      j=j+1
      xl(j)=i
      endif
   46 continue
      xt(j)=n
      cntr22=-1
      do 821 i=1,n
      xn(i)=0
  821 continue
      if(rho(1).eq.0)then
      xc(1)=0
      goto 59
      endif
      ubound=0
      n1=1
      do 811 i=2,nrho
      if(rho(i).gt.0)then
      n1=n1+1
      degr(n1)=i
      bound(n1)=i*rho(i)
      ubound=ubound+bound(n1)
      endif
  811 continue
      xc(1)=max(0,rho(1)-ubound)
  810 continue
      aux=rho(1)-xc(1)
      do 831 i=n1,2,-1
      xc(i)=min(aux,bound(i))
      aux=aux-xc(i)
  831 continue
  813 continue
      do 861 i=n1,2,-1
      aux=xc(i)
      do 851 j=xl(i),xt(i)
      xn(j)=min(aux,degr(i))
      aux=aux-xn(j)
      xs(j)=aux
  851 continue
  861 continue
      goto 212
  814 continue
      do 921 i=n1,2,-1
      do 871 j=xt(i)-1,xl(i),-1
      if(xn(j).gt.1)then
      goto 815
      endif
  871 continue
      goto 820
  815 continue
      xn(j)=xn(j)-1
      xs(j)=xs(j)+1
      do 881 jj=j+1,xt(i)
      xn(jj)=min(xn(j),xs(jj-1))
      xs(jj)=xs(jj-1)-xn(jj)
  881 continue
      if(xs(xt(i)).gt.0)then
      do 891 j=j-1,xl(i),-1
      if(xn(j).gt.1)then
      goto 815
      endif
  891 continue
      goto 820
      endif
      do 911 ii=i+1,n1
      aux=xc(ii)
      do 901 jj=xl(ii),xt(ii)
      xn(jj)=min(aux,degr(ii))
      aux=aux-xn(jj)
      xs(jj)=aux
  901 continue
  911 continue
      goto 212
  820 continue
  921 continue
      do 931 i=n1,3,-1
      if(xc(i).gt.0)then
      goto 825
      endif
  931 continue
      goto 850
  825 continue
      do 941 j=i-1,2,-1
      if(xc(j).lt.bound(j))then
      goto 830
      endif
  941 continue
      goto 850
  830 continue
      xc(j)=xc(j)+1
      aux=-1
      do 951 jj=j+1,n1
      aux=aux+xc(jj)
  951 continue
      do 961 i=n1,j+1,-1
      xc(i)=min(aux,bound(i))
      aux=aux-xc(i)
  961 continue
      if(aux.ne.0)then
      auxlin(1:srec)='gen22_2'
      call messag(1,0,0,0)
      endif
      goto 813
  850 continue
      xc(1)=xc(1)+2
      if(xc(1).le.rho(1))then
      goto 810
      endif
      goto 394
  212 continue
      if(rho(1).lt.n)then
      if(xc(1).gt.0)then
      goto 394
      endif
      endif
      if(ni.gt.1)then
      j=0
      if(mflag(1).gt.0)then
      j=1
      endif
      do 56 i=rho(1)+1,n
      if(xn(i)+j.ge.degree(i))then
      goto 814
      endif
   56 continue
      endif
      if(mflag(6).eq.1.or.mflag(3).eq.1)then
      do 267 i=rho(1)+1,n
      dta(i)=0
  267 continue
      else
      do 36 i=rho(1)+1,n
      ii=0
      if(ni.gt.1)then
      if(mflag(1).eq.1)then
      ii=2
      elseif(mflag(2).eq.1.and.xn(i).eq.0)then
      ii=2
      else
      ii=1
      endif
      endif
      j=max(0,min(degree(i)-xn(i)-ii,loop+loop))
      dta(i)=j-mod(j,2)
   36 continue
      endif
   59 continue
      do 76 i=1,n
      do 66 j=1,n
      xg(i,j)=0
   66 continue
   76 continue
      do 86 i=2,xc(1),2
      xg(i-1,i)=1
   86 continue
      limin=rho(1)+1
      limax=max(limin,n-1)
      line=limin
      aux=xc(1)
      do 106 i=limin,n
      if(xn(i).gt.0)then
      a1(i)=aux+1
      do 96 j=1,xn(i)
      aux=aux+1
      xg(aux,i)=1
      vmap(aux,1)=i
   96 continue
      else
      a1(i)=0
      endif
  106 continue
      dsum=-1
  704 continue
      dsum=dsum+1
      if(dsum.gt.loop)then
      goto 814
      endif
      aux=2*dsum
      do 701 i=n,rho(1)+1,-1
      xg(i,i)=min(aux,dta(i))
      aux=aux-xg(i,i)
  701 continue
      if(aux.eq.0)then
      goto 308
      else
      goto 814
      endif
  180 continue
      do 803 i=n,rho(1)+1,-1
      if(xg(i,i).gt.0)then
      goto 807
      endif
  803 continue
      goto 704
  807 continue
      do 819 j=i-1,rho(1)+1,-1
      if(xg(j,j).lt.dta(j))then
      xg(j,j)=xg(j,j)+2
      goto 833
      endif
  819 continue
      goto 704
  833 continue
      aux=xg(i,i)-2
      do 844 k=j+1,i-1
      aux=aux+xg(k,k)
  844 continue
      do 855 k=n,j+1,-1
      xg(k,k)=min(aux,dta(k))
      aux=aux-xg(k,k)
  855 continue
      if(aux.ne.0)then
      auxlin(1:srec)='gen22_3'
      call messag(1,0,0,0)
      endif
  308 continue
      do 116 i=limin,n
      ds(limin,i)=degree(i)-xn(i)-xg(i,i)
  116 continue
      uset(0)=0
      xset(1)=1
      jj=1
      do 136 i=2,n
      if(degree(i-1).ne.degree(i))then
      uset(jj)=i-1
      jj=jj+1
      elseif(xn(i-1).ne.xn(i))then
      uset(jj)=i-1
      jj=jj+1
      elseif(xg(i-1,i-1).gt.xg(i,i))then
      uset(jj)=i-1
      jj=jj+1
      elseif(xg(i-1,i-1).lt.xg(i,i))then
      goto 180
      endif
      xset(i)=jj
  136 continue
      uset(jj)=n
      lps(line)=dsum
   10 continue
      aux=ds(line,line)
      stren=min(str(line),loop-lps(line)+1)
      do 146 i=n,line+1,-1
      xg(line,i)=min(aux,stren,ds(line,i))
      aux=aux-xg(line,i)
  146 continue
      if(aux.gt.0)then
      goto 15
      endif
      goto 28
   05 continue
      line=limax
      goto 17
   15 continue
      if(line.lt.limin.or.line.gt.n)then
      auxlin(1:srec)='gen22_4'
      call messag(1,0,0,0)
      endif
      if(line.eq.limin)then
      goto 180
      endif
      line=line-1
   17 continue
      if(line.lt.limin.or.line.gt.n)then
      auxlin(1:srec)='gen22_5'
      call messag(1,0,0,0)
      endif
      do 166 col=n,line+1,-1
      if(xg(line,col).gt.0)then
      goto 20
      endif
  166 continue
      goto 15
   20 continue
      if(line.lt.limin.or.line.gt.n)then
      auxlin(1:srec)='gen22_6'
      call messag(1,0,0,0)
      endif
      if(line.gt.limax)then
      line=limax
      endif
      aux=xg(line,col)-1
   22 continue
      stren=min(str(line),loop-lps(line)+1)
      do 176 i=col-1,line+1,-1
      if(min(ds(line,i),stren).gt.xg(line,i))then
      goto 25
      endif
      aux=aux+xg(line,i)
  176 continue
      goto 15
   25 continue
      xg(line,i)=xg(line,i)+1
      do 196 j=n,i+1,-1
      xg(line,j)=min(aux,stren,ds(line,j))
      aux=aux-xg(line,j)
  196 continue
      goto 28
   38 continue
      aux=-1
      do 206 i=col,n
      aux=aux+xg(line,i)
  206 continue
      goto 22
   28 continue
      if(line.eq.n)then
      goto 200
      endif
      msum=0
      do 257 i=line+1,n
      ii=xg(line,i)-1
      if(ii.gt.0)then
      msum=msum+ii
      endif
  257 continue
      if(lps(line)+msum.gt.loop)then
      goto 17
      endif
      if(line.gt.limin)then
      if(xset(line).eq.xset(line-1))then
      do 211 i=limin,line-2
      aux=xg(i,line-1)-xg(i,line)
      if(aux.gt.0)then
      goto 35
      elseif(aux.lt.0)then
      auxlin(1:srec)='gen22_7'
      call messag(1,0,0,0)
      endif
  211 continue
      do 216 col=line+1,n
      aux=xg(line-1,col)-xg(line,col)
      if(aux.lt.0)then
      goto 38
      elseif(aux.gt.0)then
      goto 35
      endif
  216 continue
      endif
      endif
   35 continue
      do 226 col=line+2,n
      if(xset(col).eq.xset(col-1))then
      do 219 i=limin,line
      aux=xg(i,col-1)-xg(i,col)
      if(aux.lt.0)then
      goto 38
      elseif(aux.gt.0)then
      goto 224
      endif
  219 continue
      endif
  224 continue
  226 continue
      do 236 i=line+1,n
      ds(line+1,i)=ds(line,i)-xg(line,i)
      if(ds(line+1,i).lt.0)then
      auxlin(1:srec)='gen22_8'
      call messag(1,0,0,0)
      endif
  236 continue
      line=line+1
      lps(line)=lps(line-1)+msum
      goto 10
  200 continue
      k=1
      stack(1)=limin
      lista(limin)=1
      do 276 i=limin+1,n
      lista(i)=0
  276 continue
      do 316 j=1,ni
      if(k.lt.ni)then
      if(j.gt.k)then
      do 289 i=n,limin,-1
      if(lista(i).gt.0)then
      line=max(limin,i-1)
      goto 17
      endif
  289 continue
      auxlin(1:srec)='gen22_9'
      call messag(1,0,0,0)
      endif
      aux=stack(j)
      do 296 i=limin+1,aux-1
      if(lista(i).eq.0)then
      if(xg(i,aux).gt.0)then
      k=k+1
      stack(k)=i
      lista(i)=1
      endif
      endif
  296 continue
      do 306 i=aux+1,n
      if(lista(i).eq.0)then
      if(xg(aux,i).gt.0)then
      k=k+1
      stack(k)=i
      lista(i)=1
      endif
      endif
  306 continue
      endif
  316 continue
      if(k.ne.ni)then
      auxlin(1:srec)='gen22_10'
      call messag(1,0,0,0)
      endif
      do 346 i=1,n
      do 326 j=1,i-1
      g(i,j)=xg(j,i)
  326 continue
      do 336 j=i,n
      g(i,j)=xg(i,j)
  336 continue
  346 continue
      if(mflag(6).eq.-1)then
      do 337 i=rho(1)+1,n
      do 327 j=i,n
      if(g(i,j).gt.1)then
      goto 347
      endif
  327 continue
  337 continue
      goto 05
      endif
  347 continue
      if(mflag(1).ne.0)then
      call umpi(1,aux)
      if(aux.ne.mflag(1))then
      goto 05
      endif
      endif
      if(mflag(2).ne.0)then
      call umpi(2,aux)
      if(aux.ne.mflag(2))then
      goto 05
      endif
      endif
      if(mflag(3).ne.0)then
      call umpi(5,aux)
      if(aux.ne.mflag(3))then
      goto 05
      endif
      endif
      if(mflag(4).ne.0)then
      call umpi(4,aux)
      if(aux.ne.mflag(4))then
      goto 05
      endif
      endif
      ntadp=0
      if(mflag(10).eq.1)then
      call umpi(3,aux)
      if(aux.ne.1)then
      auxlin(1:srec)='gen22_11'
      call messag(1,0,0,0)
      endif
      endif
      nsym=0
      do 396 i=1,rho(1)-1
      do 386 j=i+1,rho(1)
      p1sym(i,j)=0
  386 continue
  396 continue
      do 87 i=1,n
      xp(i)=i
   87 continue
      goto 93
   77 continue
      do 98 i=xset(n),1,-1
      do 18 j=uset(i)-1,uset(i-1)+1,-1
      if(xp(j).lt.xp(j+1))then
      goto 102
      endif
   18 continue
   98 continue
      goto 603
  102 continue
      do 47 ii=j+2,uset(i)
      if(xp(ii).lt.xp(j))then
      goto 202
      endif
   47 continue
      ii=uset(i)+1
  202 continue
      aux=xp(j)
      xp(j)=xp(ii-1)
      xp(ii-1)=aux
      aux=uset(i)-j
      aux=aux/2
      j=j+1
      do 43 ii=0,aux-1
      aux=xp(j+ii)
      xp(j+ii)=xp(uset(i)-ii)
      xp(uset(i)-ii)=aux
   43 continue
      do 57 ii=uset(i)+1,n
      xp(ii)=ii
   57 continue
   93 continue
      if(rho(1).gt.0)then
      if(xp(rho(1)).ne.rho(1))then
      goto 603
      endif
      endif
      do 416 i=1+rho(1),n-1
      do 406 j=i,n
      aux=g(xp(i),xp(j))-g(i,j)
      if(aux.gt.0)then
      goto 05
      endif
      if(aux.lt.0)then
      goto 77
      endif
  406 continue
  416 continue
      if(rho(1).eq.0)then
      goto 114
      endif
      i=1
  110 continue
      j=vmap(i,1)
      k=xn(j)
      if(xp(j).eq.j)then
      do 436 j=i,i+k-2
      do 426 jj=j+1,i+k-1
      p1sym(j,jj)=1
  426 continue
  436 continue
      i=i+k
      if(i.le.rho(1))then
      goto 110
      endif
      goto 114
      else
      ii=a1(xp(j))
      do 446 j=ii,ii+k-1
      p1sym(i,j)=1
  446 continue
      goto 77
      endif
  114 continue
      if(mflag(12).ne.1)then
      if(psym(0).lt.0)then
      psyms=0
      psym(0)=stibs(1)
      call aoib(1)
      stib(stibs(1))=eoia
      endif
      jj=n-rho(1)
      ii=nsym*jj-psyms
      if(ii.gt.0)then
      if(psym(0).eq.stibs(1)-1-psyms)then
      call aoib(ii)
      psyms=psyms+ii
      stib(stibs(1))=eoia
      else
      auxlin(1:srec)='internal inconsistency in routine gen22'
      call messag(1,0,0,0)
      endif
      endif
      ii=nsym*jj
      if(ii.le.psyms)then
      if(nsym.gt.0)then
      ii=psym(0)+(ii-n)
      do 456 j=rho(1)+1,n
      stib(ii+j)=xp(j)
  456 continue
      endif
      else
      auxlin(1:srec)='internal inconsistency in routine gen22'
      call messag(1,0,0,0)
      endif
      endif
      nsym=nsym+1
      goto 77
  603 continue
      ns1=0
      do 476 j=2,rho(1)
      do 466 i=1,j-1
      if(p1sym(i,j).eq.1)then
      ns1=ns1+1
      p1l(ns1)=i
      p1r(ns1)=j
      endif
  466 continue
  476 continue
      if(mflag(16).eq.0.and.mflag(17).eq.0)then
      return
      endif
      do 486 i=1,rho(1)
      tail(i)=i
      head(i)=vmap(i,1)
  486 continue
      do 493 i=1,degree(n)
      nemul(i)=0
  493 continue
      aux=rho(1)+1
      do 516 i=rho(1)+1,n
      do 506 j=i,n
      if(g(i,j).eq.0)then
      eg(i,j)=0
      else
      eg(i,j)=aux
      if(i.eq.j)then
      jj=g(i,j)/2
      else
      jj=g(i,j)
      nemul(jj)=nemul(jj)+1
      emul(jj,nemul(jj))=aux
      endif
      aux=aux+jj
      do 496 ii=1,jj
      tail(aux-ii)=i
      head(aux-ii)=j
  496 continue
      endif
  506 continue
  516 continue
      if(mflag(17).eq.0)then
      return
      endif
      do 526 i=rho(1)+1,nli
      intree(i)=0
  526 continue
      do 536 i=1,n
      listb(i)=0
  536 continue
      aux=0
      ntree=0
      do 553 i=rho(1)+1,n
      lista(i)=0
  553 continue
      ii=degree(n)
  571 continue
      if(ii.lt.1)then
      goto 582
      endif
      jj=nemul(ii)
  584 continue
      if(jj.lt.1)then
      ii=ii-1
      goto 571
      endif
      i=tail(emul(ii,jj))
      j=head(emul(ii,jj))
      if(lista(i).ne.0.and.lista(i).eq.lista(j))then
      jj=jj-1
      goto 584
      endif
      if(lista(i).eq.0.and.lista(j).eq.0)then
      ntree=ntree+1
      lista(i)=ntree
      lista(j)=ntree
      elseif(lista(i).eq.0)then
      lista(i)=lista(j)
      elseif(lista(j).eq.0)then
      lista(j)=lista(i)
      elseif(lista(i).ne.lista(j))then
      ij=lista(j)
      do 579 k=rho(1)+1,n
      if(lista(k).eq.ij)then
      lista(k)=lista(i)
      endif
  579 continue
      else
      auxlin(1:srec)='gen22_12'
      call messag(1,0,0,0)
      endif
      intree(eg(i,j))=1
      listb(i)=listb(i)+1
      listb(j)=listb(j)+1
      aux=aux+1
      if(aux+rho(1).lt.n-1)then
      jj=jj-1
      goto 584
      endif
  582 continue
      do 586 i=1,rho(1)
      intree(i)=1
  586 continue
      do 616 i=1,nli
      do 606 j=1,rho(1)+loop
      flow(i,j)=0
  606 continue
  616 continue
      do 596 i=1,rho(1)
      flow(i,i)=1
  596 continue
      aux=0
      do 626 i=rho(1)+1,nli
      if(intree(i).eq.0)then
      aux=aux+1
      flow(i,rho(1)+aux)=1
      endif
  626 continue
      i=rho(1)+1
  732 continue
      if(listb(i).eq.1)then
      goto 741
      endif
      i=i+1
      if(i.le.n)then
      goto 732
      else
      goto 837
      endif
  741 continue
      do 636 j=rho(1)+1,i-1
      if(xg(j,i).gt.0)then
      if(intree(eg(j,i)).ne.0)then
      goto 750
      endif
      endif
  636 continue
      do 646 j=i+1,n
      if(xg(i,j).gt.0)then
      if(intree(eg(i,j)).ne.0)then
      goto 750
      endif
      endif
  646 continue
      auxlin(1:srec)='gen22_13'
      call messag(1,0,0,0)
  750 continue
      j=eg(min(i,j),max(i,j))
      if(rho(1).gt.1.and.xn(i).gt.0)then
      do 656 k=1,rho(1)
      if(vmap(k,1).eq.i)then
      if(tail(j).eq.i)then
      flow(j,k)=flow(j,k)+1
      else
      flow(j,k)=flow(j,k)-1
      endif
      endif
  656 continue
      endif
      do 696 ii=rho(1)+1,n
      if(ii.ne.i)then
      aux=eg(min(i,ii),max(i,ii))
      do 686 k=0,g(i,ii)-1
      if(intree(aux+k).eq.0)then
      if(head(j).eq.head(aux+k).or.
     :tail(j).eq.tail(aux+k))then
      do 666 ij=1,rho(1)+loop
      flow(j,ij)=flow(j,ij)-flow(aux+k,ij)
  666 continue
      else
      do 676 ij=1,rho(1)+loop
      flow(j,ij)=flow(j,ij)+flow(aux+k,ij)
  676 continue
      endif
      endif
  686 continue
      endif
  696 continue
      intree(j)=0
      listb(tail(j))=listb(tail(j))-1
      listb(head(j))=listb(head(j))-1
      i=rho(1)+1
      goto 732
  837 continue
      do 736 i=rho(1)+1,nli
      ii=0
      jj=0
      do 706 j=1,rho(1)
      if(flow(i,j).eq.1)then
      ii=ii+1
      elseif(flow(i,j).eq.-1)then
      jj=jj+1
      endif
  706 continue
      if(ii.gt.0.and.jj.gt.0)then
      auxlin(1:srec)='gen22_14'
      call messag(1,0,0,0)
      endif
      if(2*(ii+jj).gt.rho(1))then
      if(ii.gt.jj)then
      do 716 j=1,rho(1)
      flow(i,j)=flow(i,j)-1
  716 continue
      else
      do 726 j=1,rho(1)
      flow(i,j)=flow(i,j)+1
  726 continue
      endif
      endif
  736 continue
      nbri=0
      nrbri=0
      nsbri=0
      do 60 i=rho(1)+1,nli
      do 70 j=rho(1)+1,rho(1)+loop
      if(flow(i,j).ne.0)then
      flow(i,0)=1
      goto 60
      endif
   70 continue
      nbri=nbri+1
      do 84 j=1,rho(1)
      if(flow(i,j).ne.0)then
      flow(i,0)=2
      nrbri=nrbri+1
      goto 60
      endif
   84 continue
      nsbri=nsbri+1
      flow(i,0)=3
   60 continue
      if(zbri(nbri).eq.0)then
      goto 05
      elseif(rbri(nrbri).eq.0)then
      goto 05
      elseif(sbri(nsbri).eq.0)then
      goto 05
      endif
      if(zcho(nli-rho(1)-nbri).eq.0)then
      goto 05
      endif
      if(mflag(5).ne.0)then
      do 1020 i=rho(1)+1,nli
      do 1015 j=1,i-1
      do 1003 k=rho(1)+1,rho(1)+loop
      if(flow(i,k)+flow(j,k).ne.0)then
      goto 1005
      endif
 1003 continue
      do 1004 k=2,rho(1)
      if(flow(i,k-1)+flow(j,k-1).ne.flow(i,k)+flow(j,k))then
      goto 1005
      endif
 1004 continue
      if(mflag(5).eq.1)then
      goto 05
      else
      goto 743
      endif
 1005 continue
      do 1007 k=rho(1)+1,rho(1)+loop
      if(flow(i,k)-flow(j,k).ne.0)then
      goto 1010
      endif
 1007 continue
      do 1008 k=2,rho(1)
      if(flow(i,k-1)-flow(j,k-1).ne.flow(i,k)-flow(j,k))then
      goto 1010
      endif
 1008 continue
      if(mflag(5).eq.1)then
      goto 05
      else
      goto 743
      endif
 1010 continue
 1015 continue
 1020 continue
      if(mflag(5).eq.-1)then
      goto 05
      endif
      endif
  743 continue
      do 756 i=1,nli
      s=0
      do 746 j=1,rho(1)+loop
      s=s+abs(flow(i,j))
      if(abs(flow(i,j)).gt.1)then
      auxlin(1:srec)='gen22_15'
      call messag(1,0,0,0)
      endif
  746 continue
  756 continue
      do 766 i=1,rho(1)
      eg(i,vmap(i,1))=i
  766 continue
      do 806 i=rho(1)+1,n
      aux=0
      do 796 ii=1,rho(1)
      s=0
      do 786 j=1,n
      if(j.ne.i)then
      k=eg(min(i,j),max(i,j))
      do 776 jj=0,g(i,j)-1
      if(head(k).eq.i)then
      s=s+flow(eg(min(i,j),max(i,j))+jj,ii)
      else
      s=s-flow(eg(min(i,j),max(i,j))+jj,ii)
      endif
  776 continue
      endif
  786 continue
      if(s.ne.aux)then
      if(ii.eq.1)then
      aux=s
      else
      auxlin(1:srec)='gen22_16'
      call messag(1,0,0,0)
      endif
      endif
  796 continue
  806 continue
      do 846 i=rho(1)+1,n
      do 836 ii=rho(1)+1,rho(1)+loop
      s=0
      do 826 j=1,n
      if(j.ne.i)then
      do 816 jj=0,g(i,j)-1
      k=eg(min(i,j),max(i,j))
      if(head(k).eq.i)then
      s=s+flow(k+jj,ii)
      else
      s=s-flow(k+jj,ii)
      endif
  816 continue
      endif
  826 continue
      if(s.ne.0)then
      auxlin(1:srec)='gen22_17'
      call messag(1,0,0,0)
      endif
  836 continue
  846 continue
      return
  394 continue
      cntr22=0
      return
  937 continue
      auxlin(1:srec)='gen22_18'
      call messag(1,0,0,0)
      end
      subroutine init0
      implicit integer(a-z)
      save
      parameter ( srec=81 )
      parameter ( lfeed=10, dquote=34, rquote=39 )
      parameter ( sibuff=1048576, scbuff=524288, mimal=2, mcmal=2 )
      parameter ( eoia=-63 )
      common/z4in/mflag(1:21)
      common/z6in/dunit,munit,ounit,sunit,funit
      common/z2g/dsym,dis,ndiag
      common/z5g/psym(0:0),psyms,nsym
      common/z21g/punct1(0:127)
      character*(srec) auxlin
      common/z22g/auxlin
      common/z26g/kes(0:0),kloo(15:32,1:4),kle(0:0),
     :pstke(0:0),wstke(0:0)
      common/z27g/pkey(0:0),wkey(0:0),prevl(0:0),nextl(0:0)
      common/z28g/popt3(0:0),wopt3(0:0),fopt3(0:0),vopt3(0:0)
      common/z29g/popt5(0:0),wopt5(0:0),copt5(0:0)
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z32g/acf0(0:127),acf1(0:127),sucpal(0:11,0:11)
      common/z35g/stibs(0:mimal),stcbs(0:mcmal),imal(0:0),cmal(0:0)
      common/z44g/sdiap,sdial,sgrap,sgral,stotp,stotl
      common/z45g/qdatp,qdatl,qvp,qvl,dprefp,dprefl
      common/z46g/popt1(0:0),wopt1(0:0),copt1(0:0)
      common/z47g/popt7(0:0),wopt7(0:0),copt7(0:0)
      common/z56g/prom,lrom
      imal(0)=1
      do 01 i=0,mimal
      stibs(i)=0
   01 continue
      cmal(0)=1
      do 03 i=0,mcmal
      stcbs(i)=0
   03 continue
      i=0
      call spack('8171826570131914171420;',i,1,uc,nos)
      qvp=stib(stibs(1)+1)
      qvl=stib(stibs(1)+2)
      print *,' '
      print *,' -------------------------------------------------------'
      print *,'                       '//stcb(qvp:qvp-1+qvl)
      print *,' -------------------------------------------------------'
      print *,' '
      dunit=0
      munit=0
      ounit=0
      sunit=0
      funit=0
      nos=0
      uc=0
      call spack('817182657014686584;',i,1,uc,nos)
      qdatp=stib(stibs(1)+1)
      qdatl=stib(stibs(1)+2)
      call spack('6885657613;',i,1,uc,nos)
      dprefp=stib(stibs(1)+1)
      dprefl=stib(stibs(1)+2)
      call spack('006873657182657783;',i,1,uc,nos)
      sdiap=stib(stibs(1)+1)
      sdial=stib(stibs(1)+2)
      call spack('00718265807283;',i,1,uc,nos)
      sgrap=stib(stibs(1)+1)
      sgral=stib(stibs(1)+2)
      call spack('0000847984657600290000;',i,1,uc,nos)
      stotp=stib(stibs(1)+1)
      stotl=stib(stibs(1)+2)
      ndiag=0
      psym(0)=-1
      nos=0
      uc=1
      nkey=0
      call spack('798584808584,0,1;',nkey,0,uc,nos)
      call spack('8384897669,1,2;',nkey,0,uc,nos)
      call spack('7779686976,2,4;',nkey,0,uc,nos)
      call spack('7378,4,5;',nkey,0,uc,nos)
      call spack('798584,5,6;',nkey,0,uc,nos)
      call spack('7679798083,6,7;',nkey,0,uc,nos)
      call spack('76797980637779776978848577,7,8;',nkey,0,uc,nos)
      call spack('79808473797883,8,-1;',nkey,0,uc,nos)
      kk=4
      call aoib(kk)
      j=stibs(1)
      do 05 i=1,kk
      stib(j)=eoia
      j=j-1
   05 continue
      i=nkey+1
      call trm(kk,i)
      nextl(0)=stibs(1)-i
      prevl(0)=nextl(0)-i
      wkey(0)=prevl(0)-i
      pkey(0)=wkey(0)-i
      nos=0
      uc=1
      nopt1=0
      call spack('84828569,1;',nopt1,0,uc,nos)
      call spack('7065768369,-1;',nopt1,0,uc,nos)
      kk=3
      call aoib(kk)
      j=stibs(1)
      do 08 i=1,kk
      stib(j)=eoia
      j=j-1
   08 continue
      i=nopt1+1
      call trm(kk,i)
      copt1(0)=stibs(1)-i
      wopt1(0)=copt1(0)-i
      popt1(0)=wopt1(0)-i
      nos=1
      uc=1
      nopt3=0
      call spack('7978698073,1,1;',nopt3,0,uc,nos)
      call spack('7978698082,1,-1;',nopt3,0,uc,nos)
      call spack('78798465688079766983,2,1;',nopt3,0,uc,nos)
      call spack('8465688079766983,2,-1;',nopt3,0,uc,nos)
      call spack('7879837865737683,3,1;',nopt3,0,uc,nos)
      call spack('837865737683,3,-1;',nopt3,0,uc,nos)
      call spack('79788372697676,4,1;',nopt3,0,uc,nos)
      call spack('7970708372697676,4,-1;',nopt3,0,uc,nos)
      call spack('7879837371776583,5,1;',nopt3,0,uc,nos)
      call spack('837371776583,5,-1;',nopt3,0,uc,nos)
      call spack('837377807669,6,1;',nopt3,0,uc,nos)
      call spack('787984837377807669,6,-1;',nopt3,0,uc,nos)
      call spack('8479807976,8,1;',nopt3,0,uc,nos)
      call spack('707679798083,9,1;',nopt3,0,uc,nos)
      kk=4
      call aoib(kk)
      j=stibs(1)
      do 11 i=1,kk
      stib(j)=eoia
      j=j-1
   11 continue
      i=nopt3+1
      call trm(kk,i)
      vopt3(0)=stibs(1)-i
      fopt3(0)=vopt3(0)-i
      wopt3(0)=fopt3(0)-i
      popt3(0)=wopt3(0)-i
      nos=0
      uc=1
      nopt5=0
      call spack('7380827980,1;',nopt5,0,uc,nos)
      call spack('668273687169,2;',nopt5,0,uc,nos)
      call spack('6772798268,3;',nopt5,0,uc,nos)
      call spack('82668273687169,4;',nopt5,0,uc,nos)
      call spack('83668273687169,5;',nopt5,0,uc,nos)
      call spack('86838577,6;',nopt5,0,uc,nos)
      call spack('80838577,7;',nopt5,0,uc,nos)
      kk=3
      call aoib(kk)
      j=stibs(1)
      do 25 i=1,kk
      stib(j)=eoia
      j=j-1
   25 continue
      i=nopt5+1
      call trm(kk,i)
      copt5(0)=stibs(1)-i
      wopt5(0)=copt5(0)-i
      popt5(0)=wopt5(0)-i
      nopt7=0
      nos=1
      uc=1
      call spack('6988846982786576,5;',nopt7,0,uc,nos)
      call spack('78798465688079766983,1;',nopt7,0,uc,nos)
      kk=3
      call aoib(kk)
      j=stibs(1)
      do 15 i=1,kk
      stib(j)=eoia
      j=j-1
   15 continue
      i=nopt7+1
      call trm(kk,i)
      copt7(0)=stibs(1)-i
      wopt7(0)=copt7(0)-i
      popt7(0)=wopt7(0)-i
      nos=0
      uc=0
      nstke=0
      call spack('8082797679718569,1,0;',nstke,0,uc,nos)
      call spack('68736571826577,2,0;',nstke,0,uc,nos)
      call spack('6980737679718569,3,0;',nstke,0,uc,nos)
      call spack('69887384,4,0;',nstke,0,uc,nos)
      call spack('677977776578686376797980,11,5;',nstke,0,uc,nos)
      call spack('6779777765786863767378696376797980,12,5;',nstke,0,
     :uc,nos)
      call spack('73786376797980,13,2;',nstke,0,uc,nos)
      call spack('7985846376797980,14,2;',nstke,0,uc,nos)
      call spack('808279806571658479826376797980,15,2;',nstke,0,uc,nos)
      call spack('8669828469886376797980,16,2;',nstke,0,uc,nos)
      call spack('8265896376797980,17,2;',nstke,0,uc,nos)
      call spack('697868,19,7;',nstke,0,uc,nos)
      call spack('66656775,21,7;',nstke,0,uc,nos)
      call spack('78698776737869,22,7;',nstke,0,uc,nos)
      call spack('7073697668,31,2;',nstke,0,uc,nos)
      kloo(15,1)=1
      kloo(15,2)=1
      kloo(15,3)=0
      kloo(15,4)=1
      call spack('7779776978848577,32,2;',nstke,0,uc,nos)
      kloo(16,1)=1
      kloo(16,2)=1
      kloo(16,3)=0
      kloo(16,4)=1
      call spack('68856576137073697668,33,2;',nstke,0,uc,nos)
      kloo(17,1)=1
      kloo(17,2)=1
      kloo(17,3)=0
      kloo(17,4)=1
      call spack('68856576137779776978848577,34,2;',nstke,0,uc,nos)
      kloo(18,1)=1
      kloo(18,2)=1
      kloo(18,3)=0
      kloo(18,4)=1
      call spack('70736976686383737178,35,2;',nstke,0,uc,nos)
      kloo(19,1)=1
      kloo(19,2)=1
      kloo(19,3)=0
      kloo(19,4)=1
      call spack('70736976686384898069,40,2;',nstke,0,uc,nos)
      kloo(20,1)=1
      kloo(20,2)=1
      kloo(20,3)=0
      kloo(20,4)=1
      call spack('80827980657165847982637378686988,41,2;',nstke,0,
     :uc,nos)
      kloo(21,1)=0
      kloo(21,2)=1
      kloo(21,3)=0
      kloo(21,4)=1
      call spack('7073697668637378686988,42,2;',nstke,0,uc,nos)
      kloo(22,1)=1
      kloo(22,2)=1
      kloo(22,3)=0
      kloo(22,4)=1
      call spack('826589637378686988,43,2;',nstke,0,uc,nos)
      kloo(23,1)=1
      kloo(23,2)=1
      kloo(23,3)=0
      kloo(23,4)=1
      call spack('866982846988637378686988,44,2;',nstke,0,uc,nos)
      kloo(24,1)=1
      kloo(24,2)=1
      kloo(24,3)=1
      kloo(24,4)=1
      call spack('68856576137073697668637378686988,45,2;',nstke,0,
     :uc,nos)
      kloo(25,1)=0
      kloo(25,2)=1
      kloo(25,3)=0
      kloo(25,4)=1
      call spack('6885657613826589637378686988,46,2;',nstke,0,uc,nos)
      kloo(26,1)=0
      kloo(26,2)=1
      kloo(26,3)=0
      kloo(26,4)=1
      call spack('6885657613866982846988637378686988,47,2;',nstke,0,
     :uc,nos)
      kloo(27,1)=0
      kloo(27,2)=1
      kloo(27,3)=0
      kloo(27,4)=1
      call spack('86698284698863686971826969,48,2;',nstke,0,uc,nos)
      kloo(28,1)=1
      kloo(28,2)=1
      kloo(28,3)=1
      kloo(28,4)=1
      call spack('688565761386698284698863686971826969,49,2;',nstke,0,
     :uc,nos)
      kloo(29,1)=0
      kloo(29,2)=1
      kloo(29,3)=0
      kloo(29,4)=1
      call spack('766971637378686988,51,2;',nstke,0,uc,nos)
      kloo(30,1)=1
      kloo(30,2)=0
      kloo(30,3)=0
      kloo(30,4)=0
      call spack('7378637378686988,52,2;',nstke,0,uc,nos)
      kloo(31,1)=1
      kloo(31,2)=0
      kloo(31,3)=0
      kloo(31,4)=0
      call spack('798584637378686988,53,2;',nstke,0,uc,nos)
      kloo(32,1)=1
      kloo(32,2)=0
      kloo(32,3)=0
      kloo(32,4)=0
      call spack('83737178,61,2;',nstke,0,uc,nos)
      call spack('7773788583,62,2;',nstke,0,uc,nos)
      call spack('68736571826577637378686988,71,6;',nstke,0,uc,nos)
      call spack('838977776984828963706567847982,72,2;',nstke,0,uc,nos)
      call spack('838977776984828963788577666982,73,2;',nstke,0,uc,nos)
      call spack('8082798065716584798283,74,2;',nstke,0,uc,nos)
      call spack('76697183,75,2;',nstke,0,uc,nos)
      call spack('7679798083,76,2;',nstke,0,uc,nos)
      call spack('8669828473676983,77,2;',nstke,0,uc,nos)
      call spack('76697183637378,78,2;',nstke,0,uc,nos)
      call spack('7669718363798584,79,2;',nstke,0,uc,nos)
      call spack('80827971826577,81,5;',nstke,0,uc,nos)
      call spack('677977776578686368658465,82,5;',nstke,0,uc,nos)
      kk=4
      call aoib(kk)
      j=stibs(1)
      do 45 i=1,kk
      stib(j)=eoia
      j=j-1
   45 continue
      i=nstke+1
      call trm(kk,i)
      kle(0)=stibs(1)-i
      kes(0)=kle(0)-i
      wstke(0)=kes(0)-i
      pstke(0)=wstke(0)-i
      do 102 i=0,31
      acf0(i)=-1
  102 continue
      do 103 i=32,126
      acf0(i)=0
  103 continue
      acf0(127)=-1
      acf0(lfeed)=1
      acf0(ichar(' '))=2
      acf0(dquote)=3
      acf0(rquote)=4
      acf0(ichar('('))=5
      acf0(ichar(')'))=6
      acf0(ichar(','))=7
      acf0(ichar(';'))=8
      acf0(ichar('='))=9
      acf0(ichar('['))=10
      acf0(ichar(']'))=11
      do 104 i=0,47
      acf1(i)=-1
  104 continue
      do 105 i=48,57
      acf1(i)=1
  105 continue
      do 106 i=58,64
      acf1(i)=-1
  106 continue
      do 107 i=65,90
      acf1(i)=2
  107 continue
      do 108 i=91,96
      acf1(i)=-1
  108 continue
      acf1(95)=0
      do 109 i=97,122
      acf1(i)=3
  109 continue
      do 110 i=123,127
      acf1(i)=-1
  110 continue
      do 02 i=0,127
      punct1(i)=0
   02 continue
      punct1(ichar('='))=1
      punct1(ichar(','))=2
      punct1(ichar(';'))=3
      punct1(rquote)=4
      punct1(ichar('['))=5
      punct1(ichar(']'))=6
      punct1(ichar('('))=-7
      punct1(ichar(')'))=7
      prom=stcbs(cmal(0))+1
      lrom=3
      call aocb(lrom)
      stcb(prom:stcbs(cmal(0)))=char(239)//char(187)//char(191)
      do 220 i=0,11
      do 210 j=0,11
      sucpal(i,j)=-1
  210 continue
  220 continue
      sucpal(0,2)=0
      sucpal(0,10)=1
      sucpal(1,0)=1
      sucpal(1,1)=1
      sucpal(1,2)=1
      sucpal(1,7)=1
      sucpal(1,8)=2
      sucpal(1,9)=3
      sucpal(1,11)=10
      sucpal(2,0)=2
      sucpal(2,1)=2
      sucpal(2,2)=2
      sucpal(2,9)=3
      sucpal(3,0)=3
      sucpal(3,1)=3
      sucpal(3,2)=3
      sucpal(3,3)=-1
      sucpal(3,4)=5
      sucpal(3,5)=6
      sucpal(3,7)=2
      sucpal(3,11)=10
      sucpal(4,0)=4
      sucpal(4,2)=4
      sucpal(4,3)=3
      sucpal(4,4)=4
      sucpal(4,5)=4
      sucpal(4,6)=4
      sucpal(4,7)=4
      sucpal(4,8)=4
      sucpal(4,9)=4
      sucpal(4,10)=4
      sucpal(4,11)=4
      sucpal(5,0)=5
      sucpal(5,2)=5
      sucpal(5,3)=5
      sucpal(5,4)=3
      sucpal(5,5)=5
      sucpal(5,6)=5
      sucpal(5,7)=5
      sucpal(5,8)=5
      sucpal(5,9)=5
      sucpal(5,10)=5
      sucpal(5,11)=5
      sucpal(6,0)=6
      sucpal(6,1)=6
      sucpal(6,2)=6
      sucpal(6,3)=-1
      sucpal(6,4)=8
      sucpal(6,6)=9
      sucpal(6,7)=6
      sucpal(7,0)=7
      sucpal(7,2)=7
      sucpal(7,3)=6
      sucpal(7,4)=7
      sucpal(7,5)=7
      sucpal(7,6)=7
      sucpal(7,7)=7
      sucpal(7,8)=7
      sucpal(7,9)=7
      sucpal(7,10)=7
      sucpal(7,11)=7
      sucpal(8,0)=8
      sucpal(8,2)=8
      sucpal(8,3)=8
      sucpal(8,4)=6
      sucpal(8,5)=8
      sucpal(8,6)=8
      sucpal(8,7)=8
      sucpal(8,8)=8
      sucpal(8,9)=8
      sucpal(8,10)=8
      sucpal(8,11)=8
      sucpal(9,1)=9
      sucpal(9,2)=9
      sucpal(9,7)=2
      sucpal(9,11)=10
      sucpal(10,1)=11
      do 101 i=1,21
      mflag(i)=0
  101 continue
      end
      integer function lupty(lcode)
      implicit integer(a-z)
      save
      parameter ( sibuff=1048576, scbuff=524288, mimal=2, mcmal=2 )
      parameter ( srec=81 )
      common/z35g/stibs(0:mimal),stcbs(0:mcmal),imal(0:0),cmal(0:0)
      character*(srec) auxlin
      common/z22g/auxlin
      integer loopt(11:19)
      data loopt(11)/1/
      data loopt(12)/2/
      data loopt(13)/3/
      data loopt(14)/3/
      data loopt(15)/4/
      data loopt(16)/5/
      data loopt(17)/6/
      data loopt(18)/0/
      data loopt(19)/-1/
      if(lcode.lt.11.or.lcode.gt.19.or.loopt(lcode).eq.0)then
      lupty=0
      auxlin(1:srec)='wrong input for function lupty'
      call messag(1,0,0,0)
      else
      lupty=loopt(lcode)
      endif
      end
      subroutine aocb(delta)
      implicit integer(a-z)
      save
      parameter ( sibuff=1048576, scbuff=524288, mimal=2, mcmal=2 )
      parameter ( srec=81 )
      common/z35g/stibs(0:mimal),stcbs(0:mcmal),imal(0:0),cmal(0:0)
      character*(srec) auxlin
      common/z22g/auxlin
      integer delta
      if(delta.gt.0)then
      if(delta.gt.scbuff.or.stcbs(cmal(0)).gt.scbuff-delta)then
      auxlin(1:srec)='internal buffer size (scbuff) is too small'
      call messag(1,0,0,0)
      endif
      elseif(delta.lt.0)then
      if(delta.lt.stcbs(cmal(0)-1)-stcbs(cmal(0)))then
      auxlin(1:srec)='internal error in subroutine aocb'
      call messag(1,0,0,0)
      endif
      endif
      stcbs(cmal(0))=stcbs(cmal(0))+delta
      end
      subroutine vaocb(delta)
      implicit integer(a-z)
      save
      parameter ( sibuff=1048576, scbuff=524288, mimal=2, mcmal=2 )
      parameter ( srec=81 )
      common/z35g/stibs(0:mimal),stcbs(0:mcmal),imal(0:0),cmal(0:0)
      character*(srec) auxlin
      common/z22g/auxlin
      integer delta
      if(delta.gt.0)then
      if(delta.gt.scbuff.or.stcbs(cmal(0)).gt.scbuff-delta)then
      auxlin(1:srec)='internal buffer size (scbuff) is too small'
      call messag(1,0,0,0)
      endif
      elseif(delta.lt.0)then
      auxlin(1:srec)='internal error in subroutine vaocb'
      call messag(1,0,0,0)
      endif
      end
      subroutine aoib(delta)
      implicit integer(a-z)
      save
      parameter ( sibuff=1048576, scbuff=524288, mimal=2, mcmal=2 )
      parameter ( srec=81 )
      common/z35g/stibs(0:mimal),stcbs(0:mcmal),imal(0:0),cmal(0:0)
      character*(srec) auxlin
      common/z22g/auxlin
      integer delta
      if(delta.gt.0)then
      if(delta.gt.sibuff.or.stibs(imal(0)).gt.sibuff-delta)then
      auxlin(1:srec)='internal buffer size (sibuff) is too small'
      call messag(1,0,0,0)
      endif
      elseif(delta.lt.0)then
      if(delta.lt.stibs(imal(0)-1)-stibs(imal(0)))then
      auxlin(1:srec)='internal error in subroutine aoib'
      call messag(1,0,0,0)
      endif
      endif
      stibs(imal(0))=stibs(imal(0))+delta
      end
      subroutine vaoib(delta)
      implicit integer(a-z)
      save
      parameter ( sibuff=1048576, scbuff=524288, mimal=2, mcmal=2 )
      parameter ( srec=81 )
      common/z35g/stibs(0:mimal),stcbs(0:mcmal),imal(0:0),cmal(0:0)
      character*(srec) auxlin
      common/z22g/auxlin
      integer delta
      if(delta.gt.0)then
      if(delta.gt.sibuff.or.stibs(imal(0)).gt.sibuff-delta)then
      auxlin(1:srec)='internal buffer size (sibuff) is too small'
      call messag(1,0,0,0)
      endif
      elseif(delta.lt.0)then
      if(delta.lt.stibs(imal(0)-1)-stibs(imal(0)))then
      auxlin(1:srec)='internal error in subroutine aoib'
      call messag(1,0,0,0)
      endif
      endif
      end
      subroutine eput(istop,nl1,nl2,nf1)
      implicit integer(a-z)
      save
      parameter ( srec=81 )
      parameter ( sibuff=1048576, scbuff=524288, mimal=2, mcmal=2 )
      character*(srec) auxlin
      common/z22g/auxlin
      character*(scbuff) stcb
      common/z31g/stcb
      common/z35g/stibs(0:mimal),stcbs(0:mcmal),imal(0:0),cmal(0:0)
      common/z45g/qdatp,qdatl,qvp,qvl,dprefp,dprefl
      character*(srec) lin1
      character*(2*srec) lin2
      slin1=56
      slin2=2*srec
      lin1(1:1)=' '
      do 10 i=2,slin1
      lin1(i:i)='='
   10 continue
      do 20 j2=srec,1,-1
      if(auxlin(j2:j2).ne.' ')then
      goto 30
      endif
   20 continue
      j1=0
      j2=0
      goto 50
   30 continue
      do 40 j1=1,j2
      if(auxlin(j1:j1).ne.' ')then
      goto 50
      endif
   40 continue
   50 continue
      if(istop.eq.0)then
      i2=11
      lin2(1:i2)='  warning: '
      else
      i2=9
      lin2(1:i2)='  error: '
      endif
      if(j1.gt.0)then
      if(j2-j1+1.le.slin2-i2)then
      lin2(i2+1:i2+1+(j2-j1))=auxlin(j1:j2)
      i2=i2+(j2-j1)+1
      endif
      endif
      if(nf1.gt.0)then
      if(nf1.eq.1)then
      i=7+qdatl
      if(i2+i.le.slin2)then
      lin2(i2+1:i2+i)=', file '//stcb(qdatp:qdatp-1+qdatl)
      endif
      elseif(nf1.eq.2)then
      i=12
      if(i2+i.le.slin2)then
      lin2(i2+1:i2+i)=', model file'
      endif
      elseif(nf1.eq.3)then
      i=11
      if(i2+i.le.slin2)then
      lin2(i2+1:i2+i)=', library file'
      endif
      elseif(nf1.eq.4)then
      i=12
      if(i2+i.le.slin2)then
      lin2(i2+1:i2+i)=', style file'
      endif
      elseif(nf1.eq.5)then
      i=13
      if(i2+i.le.slin2)then
      lin2(i2+1:i2+i)=', output file'
      endif
      endif
      i2=i2+i
      endif
      if(nl1.gt.0)then
      if(nl1.eq.nl2)then
      i=7
      if(i2+i.le.slin2)then
      lin2(i2+1:i2+i)=', line '
      endif
      i2=i2+i
      i3=i2+wztos(nl1)
      if(i3.le.slin2)then
      call karat(nl1,i2,lin2,srec,0)
      endif
      i2=i3
      else
      i=8
      if(i2+i.le.slin2)then
      lin2(i2+1:i2+i)=', lines '
      endif
      i2=i2+i
      i3=i2+wztos(nl1)
      if(i3.le.slin2)then
      call karat(nl1,i2,lin2,srec,0)
      endif
      i2=i3
      i=1
      if(i2+i.le.slin2)then
      lin2(i2+1:i2+1)='-'
      endif
      i2=i2+1
      i3=i2+wztos(nl2)
      if(i3.le.slin2)then
      call karat(nl2,i2,lin2,srec,0)
      endif
      i2=i3
      endif
      endif
      print *,' '
      print *,lin1(1:slin1)
      if(i2.ge.slin1)then
      j=slin1
   70 continue
      if(j.gt.0.and.lin2(j:j).ne.' ')then
      j=j-1
      goto 70
      endif
      if(j.gt.1)then
      print *,lin2(1:j-1)
      endif
      if(istop.eq.0)then
      print *,'          '//lin2(j:i2)
      else
      print *,'        '//lin2(j:i2)
      endif
      else
      print *,lin2(1:i2)
      endif
      print *,lin1(1:slin1)
      print *,' '
      end
      subroutine messag(istop,nl1,nl2,nf1)
      implicit integer(a-z)
      save
      parameter ( srec=81 )
      parameter ( lfeed=10, dquote=34, rquote=39 )
      parameter ( nfiles=5 )
      common/z6in/dunit,munit,ounit,sunit,funit
      character*(srec) auxlin
      common/z22g/auxlin
      integer xunit(nfiles)
      logical lunit
      call eput(istop,nl1,nl2,nf1)
      if(istop.eq.0)then
      return
      endif
      xunit(1)=dunit
      xunit(2)=munit
      xunit(3)=funit
      xunit(4)=sunit
      xunit(5)=ounit
      ioi=5
      do 10 i=1,nfiles
      if(xunit(i).ne.0)then
      ios=1
      lunit=.false.
      inquire(unit=xunit(i),opened=lunit,iostat=ios)
      if(ios.ne.0)then
      auxlin(1:srec)='system/filesystem error'
      call eput(istop,0,0,0)
      elseif(lunit)then
      ios=1
      if(i.ne.ioi)then
      close(unit=xunit(i),status='keep',iostat=ios)
      else
      close(unit=xunit(ioi),status='delete',iostat=ios)
      endif
      if(ios.ne.0)then
      auxlin(1:srec)='system/filesystem error'
      call eput(istop,0,0,0)
      endif
      endif
      endif
   10 continue
      stop
      end
