*
*  -------------------------------------------------------
*
*       qgraf-3.4.2
*       a module for Feynman diagram generation
*
*       copyright 1990-2019 by P. Nogueira
*
*       not to be redistributed without explicit permission
*
*       reference:
*        [1] Automatic Feynman graph generation
*            P. Nogueira
*            J. Comput. Phys. 105 (1993) 279-289
*            https://doi.org/10.1006/jcph.1993.1074
*
*       documentation:
*         [2] files 'qgraf-3.0.pdf' and 'qgraf-3.4.pdf'
*
*  -------------------------------------------------------
*
      program qgraf
      implicit integer(a-z)
      save
      parameter ( maxleg=10, maxrho=7, maxi=10, maxdeg=6 )
      parameter ( maxn=maxi+maxi-2, maxli=maxn+maxrho )
      parameter ( srec=81, ssrec=62 )
      parameter ( sibuff=1048574, scbuff=131064, swbuff=65528 )
      parameter ( ipar1= maxleg*maxleg-maxleg )
      parameter ( ipar2= ipar1/2 )
      parameter ( maxtak=2, maxpot=6 )
      parameter ( sxbuff=2040 )
      parameter ( nfiles=5 )
      common/z1in/leg(1:maxleg),nleg,nloop,incom
      common/z1g/g(1:maxn,1:maxn),rho(1:maxdeg),nrho,rhop1,rhop2
      common/z2g/dis,dsym
      common/z3g/nivd(0:0),dpntro(0:0),vparto(0:0),vval(0:0),nvert
      common/z4g/n,nli
      common/z5g/psym(0:0),psyms,nsym
      common/z6g/p1(1:maxleg),invp1(1:maxleg)
      common/z7g/lmap(1:maxn,1:maxdeg),vmap(1:maxn,1:maxdeg),
     :pmap(1:maxn,1:maxdeg),vlis(1:maxn),invlis(1:maxn)
      common/z8g/vdeg(1:maxn),xn(1:maxn)
      common/z9g/tpc(0:0)
      common/z11g/nphi,nblok,nprop,npprop
      common/z12g/cflag(1:7),dflag(1:20),jflag(1:10),mflag(1:6)
      common/z14g/zcho(0:maxli),zbri(0:maxli),zpro(0:maxli),
     :rbri(0:maxli),sbri(0:maxli)
      common/z16g/rdeg(1:maxn),amap(1:maxn,1:maxdeg)
      common/z17g/xtail(1:maxn),xhead(1:maxn),ntadp
      common/z18g/eg(1:maxn,1:maxn),flow(1:maxli,0:maxleg+maxrho),
     :net(-3:3)
      common/z19g/vfo(1:maxn)
      common/z20g/tftyp(0:0),tfnarg(0:0),tfa(0:0),tfb(0:0),tfc(0:0),
     :tfo(0:0),tf2(0:0),ntf
      character*(srec) mlin
      common/z22g/mlin
      common/z25g/ex(1:maxli),ey(1:maxli),ovm(1:maxn,1:maxdeg)
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z34g/link(0:0),antiq(0:0),blok(0:0)
      common/z35g/stibs(1:1),stcbs(1:1),stwbs(1:1)
      common/z39g/pmkp(0:0),pmkl(0:0),pmkd(0:0),pmkvpp(0:0),pmkvlp(0:0)
      common/z40g/vmkp(0:0),vmkl(0:0),vmkvpp(0:0),vmkvlp(0:0),vmks(0:0)
      common/z43g/vmkr(0:0),vmkmao(0:0),vmkmio(0:0)
      common/z2in/momep(0:maxleg),momel(0:maxleg),kpqs(1:4)
      common/z3in/lmfile,mfilea,mfileb
      common/z4in/llfile,lfilea,lfileb
      common/z5in/lsfile,sfilea,sfileb
      common/z6in/lofile,ofilea,ofileb
      common/z7in/aunit(1:nfiles)
      common/z10g/p1l(1:ipar2),p1r(1:ipar2),ns1
      common/z13g/kla(0:0),bbc(0:0)
      common/z15g/pten(1:10),iref,wiref,wsint
      common/z21g/punct1(0:0)
      common/z23g/tak(1:maxtak),pot(0:maxpot),ks
      common/z24g/iogp(1:4)
      common/z26g/kc(0:0),kse(0:0),pske(0:0),wske(0:0),klo(0:0),nske
      character*(sxbuff) stxb
      common/z27g/stxb
      common/z28g/wera(0:0),werb(0:0),nwer
      common/z29g/pkey(0:0),wkey(0:0),ikey(0:0),pokey(0:0),wokey(0:0),
     :cokey(0:0),popt1(0:0),wopt1(0:0),fopt1(0:0),vopt1(0:0),
     :popt5(0:0),wopt5(0:0),copt5(0:0),popt0(0:0),wopt0(0:0),
     :copt0(0:0),vopt0(0:0)
      common/z32g/aaf(0:4)
      common/z33g/namep(0:0),namel(0:0)
      common/z36g/udkp(0:2),udkl(0:0),udkt(0:0),udki(0:0)
      common/z37g/drecp(0:0),drecl(0:0),drecii(0:0),irecc(0:0),
     :frecc(0:0),ndrec,ncom
      common/z38g/gmkp(0:0),gmkl(0:0),gmkd(0:0),gmko(0:0),
     :gmkvp(0:0),gmkvl(0:0)
      common/z41g/nudk,ngmk,npmk0,npmk,nvmk0,nvmk
      common/z42g/pmkr(0:0),pmkvma(0:0),pmkvmi(0:0)
      common/z44g/xtstrp(0:0),xtstrl(0:0)
      common/z45g/qdatp,qdatl,qvp,qvl,dprefp,dprefl
      common/z46g/popt3(0:0),wopt3(0:0),copt3(0:0),popt9(0:0),
     :wopt9(0:0),copt9(0:0)
      common/z47g/ndiagp,ndiagl,hhp,hhl,doffp,doffl,noffp,noffl
      common/z48g/acomma,ascol,albra,arbra,alpar,arpar
      common/z49g/popt7(0:0),wopt7(0:0),copt7(0:0)
      common/z50g/tfta(0:0),tftb(0:0),tftic(0:0),ntft
      common/z51g/aeq,agt,alt,aminus,anine,aplus,aslash,atimes,azero
      common/z52g/acf(0:127),abo(0:6),acr,alf,aspace,dquote,squote
      common/z53g/pmkvp(0:0),pmkvl(0:0)
      character*(swbuff) stwb
      common/z54g/stwb
      common/z55g/aopp(0:0),aopl(0:0),aopna(0:0),aopnb(0:0)
      common/z56g/acz(0:127)
      integer xli(2*maxn+2*maxrho+maxdeg)
      call init0
      jflag(6)=1
      i=maxdeg
      if(i.lt.3)then
      mlin(1:srec)="parameter 'maxdeg' is properly set"
      call mput(1,0,0,0)
      endif
      i=maxleg
      if(i.lt.3)then
      mlin(1:srec)="parameter 'maxleg' is properly set"
      call mput(1,0,0,0)
      endif
      i=maxrho
      if(i.lt.3)then
      mlin(1:srec)="parameter 'maxrho' is properly set"
      call mput(1,0,0,0)
      endif
      i=maxi
      j=max(maxleg,maxrho)
      if(i.lt.j)then
      mlin(1:srec)="parameter 'maxi' is properly set"
      call mput(1,0,0,0)
      endif
      i=maxn
      j=2*max(maxleg,maxrho)-2
      if(i.lt.j)then
      mlin(1:srec)="parameter 'maxn' is properly set"
      call mput(1,0,0,0)
      endif
      i=maxli
      j=2*max(maxleg,maxrho)+maxrho-3
      if(i.lt.j)then
      mlin(1:srec)="parameter 'maxli' is properly set"
      call mput(1,0,0,0)
      endif
      call inputs
      if(cflag(2).ne.0)then
      call iki
      call prepos(1)
      endif
      call rflag
      if(jflag(1).ne.0)then
      goto 94
      endif
      rho(1)=nleg
      rho(2)=0
      rhop1=rho(1)+1
      jj=rhop1-3+2*nloop
      do i1=nrho,3,-1
      rho(i1)=jj/(i1-2)
      jj=jj-rho(i1)*(i1-2)
      enddo
      goto 50
   18 continue
      do i1=4,nrho
      if(rho(i1).gt.0)then
      goto 30
      endif
      enddo
      goto 94
   30 continue
      rho(i1)=rho(i1)-1
      jj=i1-2+rho(3)
      do i2=i1-1,3,-1
      rho(i2)=jj/(i2-2)
      jj=jj-rho(i2)*(i2-2)
      enddo
   50 continue
      do i1=3,nrho
      if(stib(nivd(0)+i1).eq.0)then
      if(rho(i1).gt.0)then
      goto 18
      endif
      endif
      enddo
      jflag(2)=0
      nili=-rho(1)
      do i1=3,nrho
      nili=nili+i1*rho(i1)
      enddo
      nili=nili/2
      call sdiag(1,-1)
      cntr10=0
      if(npprop.eq.0)then
      if(nili.gt.0)then
      goto 41
      endif
      endif
      if(nloop.eq.0)then
      if((dflag(1).gt.0).or.(dflag(2).gt.0))then
      if(nili.ne.0)then
      goto 41
      endif
      elseif((dflag(1).lt.0).or.(dflag(2).lt.0))then
      if(nili.eq.0)then
      goto 41
      endif
      endif
      endif
      if(zpro(nili).eq.0)then
      goto 41
      endif
      do i1=nloop,nili
      if(zcho(i1).ne.0)then
      if(zbri(nili-i1).ne.0)then
      do i2=0,nili-i1
      if(rbri(i2).ne.0)then
      if(sbri(nili-i1-i2).ne.0)then
      goto 89
      endif
      endif
      enddo
      endif
      endif
      enddo
      goto 41
   89 continue
      if(dflag(15).ne.0)then
      do i1=1,ntf
      if(abs(stib(tftyp(0)+i1)).eq.6)then
      j1=0
      j2=0
      jj=stib(stib(tfo(0)+i1)+1)
      if(stib(vmks(0)+jj).ne.1)then
      mlin(1:srec)='main_1'
      call mput(1,0,0,0)
      endif
      do i2=nrho,3,-1
      if(rho(i2).gt.0)then
      j1=j1+rho(i2)*stib(stib(vmkmio(0)+jj)+i2)
      j2=j2+rho(i2)*stib(stib(vmkmao(0)+jj)+i2)
      endif
      enddo
      if(stib(tftyp(0)+i1).gt.0)then
      if(j1.gt.stib(tfb(0)+i1))then
      goto 41
      elseif(j2.lt.stib(tfa(0)+i1))then
      goto 41
      endif
      else
      if(j1.ge.stib(tfa(0)+i1))then
      if(j2.le.stib(tfb(0)+i1))then
      goto 41
      endif
      endif
      endif
      endif
      enddo
      endif
      cntr10=1
   40 continue
      call gen10(cntr10)
      goto 42
   41 continue
      jflag(2)=1
   42 continue
      if(cntr10.eq.0)then
      if(cflag(1).ge.0)then
      call spp(3)
      endif
      goto 18
      endif
      cntr10=0
      do i1=1,rho(1)
      pmap(i1,1)=leg(p1(i1))
      pmap(vmap(i1,1),lmap(i1,1))=stib(link(0)+leg(p1(i1)))
      enddo
      vind=rho(1)
   57 continue
      vind=vind+1
      vv=vlis(vind)
      vfo(vv)=stib(stib(dpntro(0)+vdeg(vv))+pmap(vv,1))
  100 continue
      do i1=1,rdeg(vv)
      ii=stib(vfo(vv)+i1)-pmap(vv,i1)
      if(ii.gt.0)then
      goto 104
      elseif(ii.lt.0)then
      vfo(vv)=vfo(vv)+vdeg(vv)+1
      goto 100
      endif
      enddo
      goto 163
   58 continue
      vfo(vv)=vfo(vv)+vdeg(vv)+1
      do i1=1,rdeg(vv)
      if(stib(vfo(vv)+i1).ne.pmap(vv,i1))then
      goto 104
      endif
      enddo
      goto 163
  104 continue
      if(vind.ne.rhop1)then
      vind=vind-1
      vv=vlis(vind)
      goto 58
      endif
      goto 40
  163 continue
      sdeg=rdeg(vv)+g(vv,vv)
      do i1=rdeg(vv)+1,sdeg,2
      j1=vfo(vv)+i1
      j2=stib(j1)
      j3=stib(j1+1)
      if(j2.gt.j3)then
      goto 58
      elseif(j2.ne.stib(link(0)+j3))then
      goto 58
      elseif(i1+1.ne.sdeg)then
      if(j2.gt.stib(j1+2))then
      goto 58
      endif
      endif
      enddo
      do i1=rdeg(vv)+1,vdeg(vv)
      pmap(vv,i1)=stib(vfo(vv)+i1)
      enddo
      do i1=sdeg+1,vdeg(vv)-1
      if(vmap(vv,i1).eq.vmap(vv,i1+1))then
      if(pmap(vv,i1).gt.pmap(vv,i1+1))then
      goto 58
      endif
      endif
      enddo
      if(mflag(5).ne.0)then
      do i1=rdeg(vv)+1,vdeg(vv)
      if(stib(tpc(0)+pmap(vv,i1)).eq.5)then
      goto 58
      endif
      enddo
      endif
      if(mflag(6).ne.0)then
      do i1=1,ntadp
      if(xtail(i1).eq.vv)then
      jj=xhead(i1)
      elseif(xhead(i1).eq.vv)then
      jj=xtail(i1)
      else
      jj=0
      endif
      if(jj.ne.0)then
      do i2=1,vdeg(vv)
      if(jj.eq.vmap(vv,i2))then
      if(stib(tpc(0)+pmap(vv,i2)).eq.1)then
      goto 58
      endif
      endif
      enddo
      endif
      enddo
      endif
      do i1=sdeg+1,vdeg(vv)
      pmap(vmap(vv,i1),lmap(vv,i1))=stib(link(0)+pmap(vv,i1))
      enddo
      if(vind.lt.n)then
      goto 57
      endif
      if(dflag(14).eq.0)then
      goto 333
      endif
      if(nloop.eq.0)then
      if(dflag(14).gt.0)then
      goto 333
      elseif(dflag(14).lt.0)then
      goto 58
      endif
      endif
      do i1=1,n
      xli(i1)=0
      enddo
      do i1=1,n
      if(xli(i1).eq.0)then
      ii=i1
      kk=1
  820 continue
      do i2=1,vdeg(ii)
      if(i2.ne.xli(ii))then
      jj=pmap(ii,i2)
      if(stib(antiq(0)+jj).ne.0)then
      k=ii
      ii=vmap(k,i2)
      xli(ii)=lmap(k,i2)
      goto 830
      endif
      endif
      enddo
      goto 74
  830 continue
      if(ii.gt.rho(1))then
      if(ii.ne.i1)then
      kk=-kk
      goto 820
      endif
      if(kk.gt.0)then
      if(dflag(14).gt.0)then
      goto 58
      else
      goto 333
      endif
      endif
      endif
      endif
   74 continue
      enddo
      if(dflag(14).lt.0)then
      goto 58
      endif
  333 continue
      if(jflag(3).ne.0)then
      dsym=nsym
      else
      dsym=1
      jk=psym(0)-rho(1)
      do 206 i1=2,nsym
      if(mflag(1).eq.0)then
      do i2=rhop1,n
      j1=stib(jk+i2)
      j2=rdeg(i2)+1
  400 continue
      if(j2.le.vdeg(i2))then
      j3=vmap(i2,j2)
      j4=1
  410 continue
      if(j4.le.vdeg(i2))then
      if(vmap(j1,j4).ne.stib(jk+j3))then
      j4=j4+1
      goto 410
      endif
      endif
      do i3=1,g(i2,j3)
      xli(i3)=pmap(j1,j4+i3-1)
      enddo
      do i3=1,g(i2,j3)-1
      do i4=i3+1,g(i2,j3)
      if(xli(i3).gt.xli(i4))then
      ii=xli(i3)
      xli(i3)=xli(i4)
      xli(i4)=ii
      endif
      enddo
      enddo
      do i3=1,g(i2,j3)
      ii=xli(i3)-pmap(i2,j2)
      if(ii.lt.0)then
      goto 58
      elseif(ii.gt.0)then
      goto 339
      endif
      j2=j2+1
      enddo
      goto 400
      endif
      enddo
      endif
      if(mflag(4).ne.0)then
      do i2=rhop1,n
      j1=stib(jk+i2)
      if(i2.ne.j1)then
      ii=stib(vfo(j1))-stib(vfo(i2))
      if(ii.lt.0)then
      goto 58
      elseif(ii.gt.0)then
      goto 339
      endif
      endif
      enddo
      endif
      dsym=dsym+1
  339 continue
      jk=jk+(n-rho(1))
  206 continue
      endif
      do i=1,ntf
      if(stib(tfnarg(0)+i).eq.0)then
      mlin(1:srec)='main_2'
      call mput(1,0,0,0)
      endif
      if(abs(stib(tftyp(0)+i)).gt.10)then
      goto 314
      endif
      if(abs(stib(tftyp(0)+i)).eq.1)then
      ii=0
      do 319 j=rhop1,n
      do jj=rdeg(j)+1,rdeg(j)+g(j,j),2
      do ij=1,stib(tfnarg(0)+i)
      if(stib(stib(tfo(0)+i)+ij).eq.pmap(j,jj))then
      ii=ii+1
      elseif(stib(link(0)+stib(stib(tfo(0)+i)+ij)).eq.
     :pmap(j,jj))then
      ii=ii+1
      endif
      enddo
      enddo
      do jj=rdeg(j)+g(j,j)+1,vdeg(j)
      do ij=1,stib(tfnarg(0)+i)
      if(stib(stib(tfo(0)+i)+ij).eq.pmap(j,jj))then
      ii=ii+1
      elseif(stib(link(0)+stib(stib(tfo(0)+i)+ij)).eq.
     :pmap(j,jj))then
      ii=ii+1
      endif
      enddo
      enddo
  319 continue
      elseif(abs(stib(tftyp(0)+i)).eq.3)then
      ii=0
      do 519 j=rhop1,n
      do 515 jj=rdeg(j)+1,rdeg(j)+g(j,j),2
      if(flow(amap(j,jj),0).lt.0)then
      do ij=1,stib(tfnarg(0)+i)
      if(stib(stib(tfo(0)+i)+ij).eq.pmap(j,jj))then
      ii=ii+1
      else
      k=stib(link(0)+stib(stib(tfo(0)+i)+ij))
      if(k.eq.pmap(j,jj))then
      ii=ii+1
      endif
      endif
      enddo
      endif
  515 continue
      do 517 jj=rdeg(j)+g(j,j)+1,vdeg(j)
      if(flow(amap(j,jj),0).lt.0)then
      do ij=1,stib(tfnarg(0)+i)
      if(stib(stib(tfo(0)+i)+ij).eq.pmap(j,jj))then
      ii=ii+1
      else
      k=stib(link(0)+stib(stib(tfo(0)+i)+ij))
      if(k.eq.pmap(j,jj))then
      ii=ii+1
      endif
      endif
      enddo
      endif
  517 continue
  519 continue
      elseif(abs(stib(tftyp(0)+i)).lt.6)then
      if(abs(stib(tftyp(0)+i)).eq.2)then
      i1=1
      i2=2
      elseif(abs(stib(tftyp(0)+i)).eq.4)then
      i1=1
      i2=1
      elseif(abs(stib(tftyp(0)+i)).eq.5)then
      i1=2
      i2=2
      else
      mlin(1:srec)='main_3'
      call mput(1,0,0,0)
      endif
      ii=0
      do 419 j=rhop1,n
      do 415 jj=rdeg(j)+1,rdeg(j)+g(j,j),2
      ij=flow(amap(j,jj),0)
      if((ij.eq.i1).or.(ij.eq.i2))then
      do ij=1,stib(tfnarg(0)+i)
      if(stib(stib(tfo(0)+i)+ij).eq.pmap(j,jj))then
      ii=ii+1
      else
      k=stib(link(0)+stib(stib(tfo(0)+i)+ij))
      if(k.eq.pmap(j,jj))then
      ii=ii+1
      endif
      endif
      enddo
      endif
  415 continue
      do 417 jj=rdeg(j)+g(j,j)+1,vdeg(j)
      ij=flow(amap(j,jj),0)
      if((ij.eq.i1).or.(ij.eq.i2))then
      do ij=1,stib(tfnarg(0)+i)
      if(stib(stib(tfo(0)+i)+ij).eq.pmap(j,jj))then
      ii=ii+1
      else
      k=stib(link(0)+stib(stib(tfo(0)+i)+ij))
      if(k.eq.pmap(j,jj))then
      ii=ii+1
      endif
      endif
      enddo
      endif
  417 continue
  419 continue
      elseif(abs(stib(tftyp(0)+i)).eq.6)then
      jj=stib(stib(tfo(0)+i)+1)
      ii=0
      i1=stib(vmkvpp(0)+jj)
      i2=stib(vmkvlp(0)+jj)
      do k=rhop1,n
      kk=stib(vfo(k))
      ij=stoz(stcb,stib(i1+kk),stib(i1+kk)-1+stib(i2+kk))
      ii=ii+ij
      enddo
      elseif(abs(stib(tftyp(0)+i)).eq.7)then
      jj=stib(stib(tfo(0)+i)+1)
      ii=0
      i1=stib(pmkvpp(0)+jj)
      i2=stib(pmkvlp(0)+jj)
      do j=rhop1,n
      do k=rdeg(j)+1,rdeg(j)+g(j,j),2
      kk=pmap(j,k)
      ij=stoz(stcb,stib(i1+kk),stib(i1+kk)-1+stib(i2+kk))
      ii=ii+ij
      enddo
      do k=rdeg(j)+g(j,j)+1,vdeg(j)
      kk=pmap(j,k)
      ij=stoz(stcb,stib(i1+kk),stib(i1+kk)-1+stib(i2+kk))
      ii=ii+ij
      enddo
      enddo
      endif
      if(stib(tftyp(0)+i).gt.0)then
      if((ii.lt.stib(tfa(0)+i)).or.(ii.gt.stib(tfb(0)+i)))then
      goto 58
      endif
      elseif(stib(tftyp(0)+i).lt.0)then
      if((ii.ge.stib(tfa(0)+i)).and.(ii.le.stib(tfb(0)+i)))then
      goto 58
      endif
      endif
  314 continue
      enddo
      call sdiag(1,0)
      call sdiag(2,0)
      call sdiag(3,0)
      if(cflag(2).eq.0)then
      goto 58
      endif
      if(nloop.eq.0)then
      goto 231
      endif
      do 236 i=rhop1,n
      j=rdeg(i)+1
  420 continue
      if(j.le.vdeg(i))then
      ii=vmap(i,j)
      aux=g(i,ii)
      k=j+aux
      if(i.ne.ii)then
  430 continue
      if(j.lt.k)then
      kk=1
      aux=pmap(i,j)
  440 continue
      if(j+kk.lt.k)then
      if(aux.eq.pmap(i,j+kk))then
      kk=kk+1
      dsym=dsym*kk
      goto 440
      endif
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
      if(j+kk+kk.lt.k)then
      if(aux.eq.pmap(i,j+kk+kk))then
      kk=kk+1
      dsym=dsym*kk
      goto 460
      endif
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
      do i=rhop1,n
      do j=1,vdeg(i)
      if(vmap(i,j).gt.nleg)then
      k=pmap(i,j)
      if(k.le.stib(link(0)+k))then
      if(k.eq.stib(link(0)+k))then
      if(i.gt.vmap(i,j))then
      goto 250
      endif
      if((i.eq.vmap(i,j)).and.(mod(j-rdeg(i),2).eq.0))then
      goto 250
      endif
      endif
      k=amap(i,j)-nleg
      ex(k)=i
      ey(k)=j
      endif
      endif
  250 continue
      enddo
      enddo
      do i=rhop1,n
      do j=1,vdeg(i)
      xli(j)=0
      enddo
      do j=1,vdeg(i)
      k=stib(stib(vparto(0)+stib(vfo(i)))+j)
      do i3=1,vdeg(i)
      if((xli(i3).eq.0).and.(pmap(i,i3).eq.k))then
      xli(i3)=1
      goto 14
      endif
      enddo
      mlin(1:srec)='main_4'
      call mput(1,0,0,0)
   14 continue
      ovm(i,j)=i3
      enddo
      enddo
      dis=1
      if(mflag(2).gt.0)then
      goto 495
      endif
      nf=0
      nf1=0
      do i=rhop1,n
      do j=1,vdeg(i)
      k=ovm(i,j)
      ii=pmap(i,k)
      if(stib(antiq(0)+ii).ne.0)then
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
      elseif(mod(k-rdeg(i),2).ne.0)then
      jj=2*(amap(i,k)-nleg)-1
      else
      jj=2*(amap(i,k)-nleg)
      endif
      if(ij.lt.0)then
      nf1=nf1+1
      endif
      xli(nf)=jj
      endif
      enddo
      enddo
  266 continue
      ii=0
      do i1=1,nf
      if(xli(i1).gt.ii)then
      ii=xli(i1)
      endif
      enddo
      if(ii.gt.0)then
      j1=0
      j2=0
      j3=nf
  293 continue
      if(xli(j3).gt.ii-2)then
      j2=j1
      j1=xli(j3)
      if(j3.ne.nf)then
      xli(j3)=xli(nf)
      dis=-dis
      endif
      nf=nf-1
      endif
      if((j3.gt.1).and.(j2.eq.0))then
      j3=j3-1
      goto 293
      endif
      if(j1.gt.j2)then
      dis=-dis
      endif
      goto 266
      endif
      do i1=1,incom
      do i2=nf,1,-1
      if(xli(i2).eq.1-2*i1)then
      if(i2.ne.nf)then
      xli(i2)=xli(nf)
      dis=-dis
      endif
      nf=nf-1
      goto 292
      endif
      enddo
  292 continue
      enddo
      do i1=rho(1),incom+1,-1
      do i2=nf,1,-1
      if(xli(i2).eq.2*(incom-i1))then
      if(i2.ne.nf)then
      xli(i2)=xli(nf)
      dis=-dis
      endif
      nf=nf-1
      goto 302
      endif
      enddo
  302 continue
      enddo
      if(nf.ne.0)then
      mlin(1:srec)='main_5'
      call mput(1,0,0,0)
      endif
  495 continue
      call compac
      goto 58
   94 continue
      if(cflag(2).ne.0)then
      call prepos(3)
      endif
      if(cflag(1).ge.0)then
      call spp(4)
      endif
      stop
      end
      subroutine ctxb(ii,jj,kk)
      implicit integer(a-z)
      save
      parameter ( srec=81, ssrec=62 )
      parameter ( sibuff=1048574, scbuff=131064, swbuff=65528 )
      parameter ( sxbuff=2040 )
      character*(srec) mlin
      common/z22g/mlin
      character*(sxbuff) stxb
      common/z27g/stxb
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z35g/stibs(1:1),stcbs(1:1),stwbs(1:1)
      j1=0
      if(ii.lt.0)then
      j1=1
      elseif(kk.lt.0)then
      j1=1
      elseif(jj.le.0)then
      j1=1
      elseif(kk+jj.gt.scbuff)then
      j1=1
      endif
      if(j1.ne.0)then
      mlin(1:srec)='ctxb_1'
      call mput(1,0,0,0)
      endif
      if(ii+jj.gt.sxbuff)then
      call uput(6)
      endif
      stxb(ii+1:ii+jj)=stcb(kk+1:kk+jj)
      ii=ii+jj
      return
      end
      subroutine spp(what)
      implicit integer(a-z)
      save
      parameter ( maxleg=10, maxrho=7, maxi=10, maxdeg=6 )
      parameter ( maxn=maxi+maxi-2, maxli=maxn+maxrho )
      parameter ( srec=81, ssrec=62 )
      parameter ( sibuff=1048574, scbuff=131064, swbuff=65528 )
      parameter ( sxbuff=2040 )
      common/z1g/g(1:maxn,1:maxn),rho(1:maxdeg),nrho,rhop1,rhop2
      common/z3g/nivd(0:0),dpntro(0:0),vparto(0:0),vval(0:0),nvert
      common/z9g/tpc(0:0)
      common/z11g/nphi,nblok,nprop,npprop
      common/z12g/cflag(1:7),dflag(1:20),jflag(1:10),mflag(1:6)
      character*(srec) mlin
      common/z22g/mlin
      character*(sxbuff) stxb
      common/z27g/stxb
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z34g/link(0:0),antiq(0:0),blok(0:0)
      common/z35g/stibs(1:1),stcbs(1:1),stwbs(1:1)
      common/z44g/xtstrp(0:0),xtstrl(0:0)
      common/z45g/qdatp,qdatl,qvp,qvl,dprefp,dprefl
      common/z47g/ndiagp,ndiagl,hhp,hhl,doffp,doffl,noffp,noffl
      common/z51g/aeq,agt,alt,aminus,anine,aplus,aslash,atimes,azero
      integer nps(0:4)
      if(what.ne.3)then
      if(jflag(7).eq.0)then
      write(unit=*,fmt='(a)')
      jflag(7)=10
      endif
      endif
      if(what.eq.0)then
      goto 01
      elseif(what.eq.1)then
      goto 11
      elseif(what.eq.2)then
      goto 21
      elseif(what.eq.3)then
      goto 31
      elseif(what.eq.4)then
      goto 41
      endif
      mlin(1:srec)='spp_1'
      call mput(1,0,0,0)
   01 continue
      call hrul(1)
      kk=1+(ssrec-qvl)/2
      write(unit=*,fmt='(a,/,a,/,a,/)')stxb(1:ssrec),stcb(1:kk)//
     :stcb(qvp:qvp-1+qvl),stxb(1:ssrec)
      jflag(6)=2
      goto 90
   11 continue
      call hrul(1)
      write(unit=*,fmt='(a,/)')stxb(1:ssrec)
      goto 90
   21 continue
      if(cflag(1).le.0)then
      goto 90
      endif
      do i1=1,2
      do i2=0,4
      nps(i2)=0
      enddo
      do i2=1,nphi
      jj=0
      if((i1.eq.1).and.(stib(tpc(0)+i2).eq.5))then
      jj=1
      elseif((i1.eq.2).and.(stib(tpc(0)+i2).ne.5))then
      jj=1
      endif
      if(jj.eq.1)then
      if(i2.lt.stib(link(0)+i2))then
      if(stib(antiq(0)+i2).eq.0)then
      nps(2)=nps(2)+1
      else
      nps(4)=nps(4)+1
      endif
      elseif(i2.eq.stib(link(0)+i2))then
      if(stib(antiq(0)+i2).eq.0)then
      nps(1)=nps(1)+1
      else
      nps(3)=nps(3)+1
      endif
      endif
      endif
      enddo
      nps(0)=nps(1)+nps(2)+nps(3)+nps(4)
      if((i1.eq.2).or.(nps(0).gt.0))then
      ii=0
      call ctxb(ii,2,0)
      call ctxb(ii,stib(xtstrl(0)+i1),stib(xtstrp(0)+i1)-1)
      j1=stcbs(1)
      jj=2
      call vaocb(jj)
      call dkar(nps(0),jj)
      call ctxb(ii,jj,j1)
      call ctxb(ii,stib(xtstrl(0)+14),stib(xtstrp(0)+14)-1)
      if(nps(0).gt.0)then
      do i2=1,4
      if(nps(i2).gt.0)then
      call ctxb(ii,2,0)
      call dkar(nps(i2),jj)
      call ctxb(ii,jj,j1)
      j2=stib(xtstrp(0)+2+i2)
      stcb(j1+1:j1+2)=stcb(j2:j2+1)
      call ctxb(ii,2,j1)
      endif
      enddo
      endif
      if(ii.ge.srec)then
      j2=stib(xtstrp(0)+18)
      j3=stib(xtstrl(0)+18)
      ii=srec-1
      stxb(srec-j3:ii)=stcb(j2+1:j2+j3)
      endif
      write(unit=*,fmt='(a)')stxb(1:ii)
      endif
      enddo
      j1=stcbs(1)
      ii=0
      call ctxb(ii,2,0)
      call ctxb(ii,stib(xtstrl(0)+7),stib(xtstrp(0)+7)-1)
      call dkar(nvert,jj)
      call ctxb(ii,jj,j1)
      call ctxb(ii,stib(xtstrl(0)+14),stib(xtstrp(0)+14)-1)
      do i1=1,nrho
      if(stib(nivd(0)+i1).gt.0)then
      call ctxb(ii,2,0)
      call dkar(i1,jj)
      call ctxb(ii,jj,j1)
      stcb(j1+1:j1+1)='^'
      call ctxb(ii,1,j1)
      call dkar(stib(nivd(0)+i1),jj)
      call ctxb(ii,jj,j1)
      endif
      enddo
      if(ii.ge.srec)then
      j2=stib(xtstrp(0)+18)
      j3=stib(xtstrl(0)+18)
      ii=srec-1
      stxb(srec-3:ii)=stcb(j2+1:j2+j3)
      endif
      write(unit=*,fmt='(/,a,/)')stxb(1:ii)
      call hrul(1)
      write(unit=*,fmt='(a,/)')stxb(1:ssrec)
      goto 90
   31 continue
      gg=0
      do i1=3,nrho
      gg=gg+(i1-2)*rho(i1)
      enddo
      call vaocb(1)
      j1=stcbs(1)
      j2=j1+1
      ii=0
      do i1=1,nrho
      if(stib(nivd(0)+i1).gt.0)then
      kk=wztos(gg/(i1-2))+2
      call dkar(i1,jj)
      if(rho(i1).gt.0)then
      call ctxb(ii,jj,j1)
      stcb(j2:j2)='^'
      call ctxb(ii,1,j1)
      call dkar(rho(i1),jj)
      call ctxb(ii,jj,j1)
      kk=kk-jj
      else
      call ctxb(ii,jj,0)
      stcb(j2:j2)=char(aminus)
      call ctxb(ii,1,j1)
      endif
      call ctxb(ii,kk,0)
      endif
      enddo
      ii=ii-2
      j2=ii
      jj=5
      call ctxb(ii,jj,0)
      call ctxb(ii,stib(xtstrl(0)+15),stib(xtstrp(0)+15)-1)
      call ctxb(ii,jj,0)
      call ctxb(ii,hhl,hhp)
      if(jflag(2).ne.0)then
      call ctxb(ii,stib(xtstrl(0)+16),stib(xtstrp(0)+16)-1)
      endif
      i2=stib(xtstrl(0)+8)
      j1=stib(xtstrl(0)+15)
      if(j2.lt.i2)then
      j4=(i2-j2)/2
      else
      j4=0
      endif
      if(jflag(6).lt.3)then
      if(j2.lt.i2)then
      j3=0
      else
      j3=(j2-i2)/2
      endif
      jj=j2-i2+2*jj+j1+j4-j3-3
      i1=stib(xtstrp(0)+8)
      i3=stib(xtstrp(0)+9+dflag(13))
      i4=i3-1+stib(xtstrl(0)+9+dflag(13))
      write(unit=*,fmt='(a,/)')stcb(1:j3+4)//stcb(i1:i1-1+i2)//
     :stcb(1:jj)//stcb(i3:i4)
      jflag(6)=3
      endif
      write(unit=*,fmt='(a)')stcb(1:j4+4)//stxb(1:ii)
      jflag(7)=0
      goto 90
   41 continue
      ii=0
      call ctxb(ii,stib(xtstrl(0)+11),stib(xtstrp(0)+11)-1)
      call ctxb(ii,ndiagl,ndiagp)
      jj=0
      if(ndiagl.eq.1)then
      j1=ndiagp+1
      if(ichar(stcb(j1:j1)).eq.azero+1)then
      jj=1
      endif
      endif
      j1=12+dflag(13)
      call ctxb(ii,stib(xtstrl(0)+j1)-jj,stib(xtstrp(0)+j1)-1)
      if(jflag(1).ne.0)then
      call ctxb(ii,stib(xtstrl(0)+16),stib(xtstrp(0)+16)-1)
      endif
      if(jflag(8).gt.0)then
      call ctxb(ii,stib(xtstrl(0)+17),stib(xtstrp(0)+17)-1)
      endif
      write(unit=*,fmt='(/,a,/)')stxb(1:ii)
   90 return
      end
      subroutine rflag
      implicit integer(a-z)
      save
      parameter ( maxleg=10, maxrho=7, maxi=10, maxdeg=6 )
      parameter ( maxn=maxi+maxi-2, maxli=maxn+maxrho )
      parameter ( srec=81, ssrec=62 )
      parameter ( sibuff=1048574, scbuff=131064, swbuff=65528 )
      common/z1in/leg(1:maxleg),nleg,nloop,incom
      common/z1g/g(1:maxn,1:maxn),rho(1:maxdeg),nrho,rhop1,rhop2
      common/z9g/tpc(0:0)
      common/z12g/cflag(1:7),dflag(1:20),jflag(1:10),mflag(1:6)
      common/z14g/zcho(0:maxli),zbri(0:maxli),zpro(0:maxli),
     :rbri(0:maxli),sbri(0:maxli)
      common/z20g/tftyp(0:0),tfnarg(0:0),tfa(0:0),tfb(0:0),tfc(0:0),
     :tfo(0:0),tf2(0:0),ntf
      character*(srec) mlin
      common/z22g/mlin
      common/z30g/stib(1:sibuff)
      common/z50g/tfta(0:0),tftb(0:0),tftic(0:0),ntft
      if(dflag(13).ne.0)then
      if(mflag(1).eq.0)then
      mlin(1:srec)="option 'topol' does not apply here"
      call mput(1,0,0,0)
      endif
      ii=0
      if(dflag(17).ne.0)then
      ii=1
      elseif(dflag(18).ne.0)then
      ii=1
      endif
      if(ii.ne.0)then
      mlin(1:srec)="*link* statements incompatible with option "//
     :"'topol'"
      call mput(1,0,0,0)
      endif
      endif
      if(dflag(14).ne.0)then
      if((mflag(2).gt.0).or.(mflag(3).eq.0))then
      mlin(1:srec)="option 'floop' does not apply here"
      call mput(1,0,0,0)
      endif
      endif
      if(jflag(4).eq.0)then
      if(dflag(2).ne.0)then
      jflag(4)=1
      elseif(dflag(5).ne.0)then
      jflag(4)=1
      elseif(dflag(7).ne.0)then
      jflag(4)=1
      elseif(dflag(9).ne.0)then
      jflag(4)=1
      elseif(dflag(18).ne.0)then
      jflag(4)=1
      endif
      endif
      if(dflag(20).ne.0)then
      goto 90
      endif
      if(dflag(18).ne.0)then
      if((dflag(1).gt.0).or.(dflag(2).gt.0))then
      ii=stib(tftic(0)+12)
      do i1=stib(tfta(0)+ii),stib(tftb(0)+ii)
      styp=stib(tftyp(0)+stib(tf2(0)+i1))
      if(styp.gt.0)then
      jflag(1)=1
      endif
      enddo
      endif
      endif
      if(dflag(1).gt.0)then
      if(nleg.ne.1)then
      if(dflag(2).lt.0)then
      jflag(1)=1
      elseif(dflag(3).lt.0)then
      jflag(1)=1
      elseif(dflag(4).lt.0)then
      jflag(1)=1
      endif
      endif
      endif
      if(dflag(1).lt.0)then
      if(dflag(2).gt.0)then
      if(dflag(3).gt.0)then
      jflag(1)=1
      endif
      endif
      endif
      if(dflag(3).lt.0)then
      if(dflag(8).gt.0)then
      jflag(1)=1
      endif
      endif
      if(dflag(4).lt.0)then
      if(dflag(5).gt.0)then
      jflag(1)=1
      endif
      endif
      if(dflag(1).gt.0)then
      if(dflag(9).gt.0)then
      if(dflag(10).lt.0)then
      jflag(1)=1
      endif
      endif
      endif
      if(dflag(4).lt.0)then
      if(dflag(10).gt.0)then
      if(dflag(12).gt.0)then
      jflag(1)=1
      endif
      endif
      endif
      if(nloop.gt.0)then
      if(dflag(1).lt.0)then
      if(dflag(10).gt.0)then
      if(dflag(8).gt.0)then
      jflag(1)=1
      elseif(dflag(12).gt.0)then
      jflag(1)=1
      endif
      endif
      endif
      endif
      if(nloop.gt.1)then
      if(dflag(1).lt.0)then
      if(dflag(10).gt.0)then
      if(dflag(9).gt.0)then
      jflag(1)=1
      endif
      endif
      endif
      endif
      if(nloop.eq.0)then
      if(dflag(1).gt.0)then
      if(nrho.lt.nleg)then
      jflag(1)=1
      endif
      endif
      if(dflag(3).lt.0)then
      jflag(1)=1
      else
      dflag(3)=0
      endif
      if(dflag(4).lt.0)then
      jflag(1)=1
      else
      dflag(4)=0
      endif
      if(dflag(5).lt.0)then
      jflag(1)=1
      else
      dflag(5)=0
      endif
      if(dflag(7).lt.0)then
      jflag(1)=1
      else
      dflag(7)=0
      endif
      if(dflag(8).lt.0)then
      jflag(1)=1
      else
      dflag(8)=0
      endif
      if(dflag(9).lt.0)then
      jflag(1)=1
      else
      dflag(9)=0
      endif
      if(dflag(11).lt.0)then
      jflag(1)=1
      else
      dflag(11)=0
      endif
      if(dflag(12).lt.0)then
      jflag(1)=1
      else
      dflag(12)=0
      endif
      if(dflag(14).lt.0)then
      jflag(1)=1
      else
      dflag(14)=0
      endif
      elseif(nloop.eq.1)then
      if(dflag(9).lt.0)then
      jflag(1)=1
      else
      dflag(9)=0
      endif
      endif
      if(mflag(6).ne.0)then
      if(nloop.eq.0)then
      mflag(6)=0
      elseif(dflag(1).gt.0)then
      mflag(6)=0
      elseif(dflag(3).gt.0)then
      mflag(6)=0
      elseif(dflag(8).gt.0)then
      mflag(6)=0
      endif
      endif
      if(nleg.eq.1)then
      if(stib(tpc(0)+leg(1)).eq.1)then
      jflag(1)=1
      endif
      endif
      if(dflag(1).gt.0)then
      do i1=1,maxli
      rbri(i1)=0
      sbri(i1)=0
      zbri(i1)=0
      enddo
      elseif(dflag(1).lt.0)then
      zbri(0)=0
      endif
      if(dflag(2).gt.0)then
      do i1=1,maxli
      rbri(i1)=0
      enddo
      elseif(dflag(2).lt.0)then
      rbri(0)=0
      endif
      if(dflag(3).gt.0)then
      do i1=1,maxli
      sbri(i1)=0
      enddo
      elseif(dflag(3).lt.0)then
      sbri(0)=0
      endif
      if(dflag(4).lt.0)then
      zbri(0)=0
      endif
      if(dflag(8).gt.0)then
      do i1=1,maxli
      sbri(i1)=0
      enddo
      endif
      if(nloop.eq.0)then
      do i1=1,maxli
      zcho(i1)=0
      enddo
      endif
      if(zcho(0).ne.0)then
      if(dflag(3).lt.0)then
      zcho(0)=0
      elseif(dflag(4).lt.0)then
      zcho(0)=0
      elseif(dflag(5).lt.0)then
      zcho(0)=0
      elseif(dflag(7).lt.0)then
      zcho(0)=0
      elseif(dflag(8).lt.0)then
      zcho(0)=0
      elseif(dflag(12).lt.0)then
      zcho(0)=0
      endif
      endif
   90 return
      end
      subroutine prepos(sec)
      implicit integer(a-z)
      save
      parameter ( srec=81, ssrec=62 )
      parameter ( eoa=-2047, nap=-8191 )
      parameter ( nfiles=5 )
      parameter ( sibuff=1048574, scbuff=131064, swbuff=65528 )
      common/z6in/lofile,ofilea,ofileb
      common/z7in/aunit(1:nfiles)
      common/z12g/cflag(1:7),dflag(1:20),jflag(1:10),mflag(1:6)
      common/z24g/iogp(1:4)
      character*(srec) mlin
      common/z22g/mlin
      common/z26g/kc(0:0),kse(0:0),pske(0:0),wske(0:0),klo(0:0),nske
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z35g/stibs(1:1),stcbs(1:1),stwbs(1:1)
      common/z36g/udkp(0:2),udkl(0:0),udkt(0:0),udki(0:0)
      common/z37g/drecp(0:0),drecl(0:0),drecii(0:0),irecc(0:0),
     :frecc(0:0),ndrec,ncom
      common/z38g/gmkp(0:0),gmkl(0:0),gmkd(0:0),gmko(0:0),
     :gmkvp(0:0),gmkvl(0:0)
      common/z41g/nudk,ngmk,npmk0,npmk,nvmk0,nvmk
      common/z45g/qdatp,qdatl,qvp,qvl,dprefp,dprefl
      common/z47g/ndiagp,ndiagl,hhp,hhl,doffp,doffl,noffp,noffl
      common/z52g/acf(0:127),abo(0:6),acr,alf,aspace,dquote,squote
      character*(swbuff) stwb
      common/z54g/stwb
      if(abs(sec-2).ne.1)then
      mlin(1:srec)='prepos_1'
      call mput(1,0,0,0)
      endif
      if(sec.eq.1)then
      call qopen(ofilea,lofile,nfiles,1)
      elseif(stwbs(1).gt.0)then
      call qout
      endif
      stwbs(1)=0
      ig=iogp(sec)
      if((ig.lt.iogp(sec)).or.(ig.ge.iogp(sec+1)))then
      goto 80
      endif
      kma=-1
      lupt=0
   20 continue
      igk=stib(ig)
      if(igk.gt.0)then
      stwbs(1)=stwbs(1)+1
      if(cflag(3).gt.0)then
      if(igk.eq.alf)then
      stwb(stwbs(1):stwbs(1))=char(acr)
      stwbs(1)=stwbs(1)+1
      endif
      endif
      if(stwbs(1).ge.swbuff)then
      goto 85
      endif
      stwb(stwbs(1):stwbs(1))=char(igk)
      ig=ig+1
      goto 20
      elseif(igk.eq.eoa)then
      goto 30
      endif
      igc=stib(ig+1)
      if(igk.eq.-1)then
      if(igc.eq.1)then
      if(stib(ig+4).eq.0)then
      stib(ig+4)=ncom
      stib(ig+5)=1
      else
      stib(ig+5)=stib(ig+5)+1
      endif
      if(stib(ig+4).lt.stib(ig+5))then
      stib(ig+4)=0
      ig=ig+1
      lupt=0
      lupi=0
      else
      lupt=1
      lupi=stib(ig+5)
      endif
      elseif(igc.eq.2)then
      if(stib(ig+4).eq.0)then
      stib(ig+4)=stib(frecc(0)+lupi)
      stib(ig+5)=stib(irecc(0)+lupi)
      else
      stib(ig+5)=stib(ig+5)+1
      endif
      if(stib(ig+4).lt.stib(ig+5))then
      stib(ig+4)=0
      ig=ig+1
      lupt=1
      lupj=0
      else
      lupt=2
      lupj=stib(ig+5)
      endif
      elseif(igc.ne.-1)then
      goto 80
      endif
      elseif(igk.eq.-2)then
      if(igc.eq.22)then
      if(jflag(10).eq.0)then
      if(kma.ge.0)then
      stwbs(1)=kma
      kma=-1
      endif
      elseif(stwbs(1).gt.0)then
      stwbs(1)=stwbs(1)-1
      if(ichar(stwb(kli:kli)).eq.acr)then
      stwbs(1)=stwbs(1)-1
      endif
      else
      goto 80
      endif
      else
      goto 80
      endif
      elseif(igk.eq.-3)then
      ip=stwbs(1)+1
      if(igc.eq.71)then
      stwbs(1)=stwbs(1)+ndiagl
      if(stwbs(1).ge.swbuff)then
      goto 85
      endif
      stwb(ip:stwbs(1))=stcb(ndiagp+1:ndiagp+ndiagl)
      elseif(igc.eq.81)then
      stwbs(1)=stwbs(1)+qvl
      if(stwbs(1).ge.swbuff)then
      goto 85
      endif
      stwb(ip:stwbs(1))=stcb(qvp:qvp-1+qvl)
      else
      goto 80
      endif
      elseif(igk.eq.-4)then
      if(igc.eq.82)then
      if(lupt.eq.1)then
      k1=stib(irecc(0)+lupi)
      k2=stib(frecc(0)+lupi)
      elseif(lupt.eq.2)then
      k1=lupj
      k2=lupj
      else
      goto 80
      endif
      do i1=k1,k2
      j2=stib(drecl(0)+i1)
      ip=stwbs(1)+1
      stwbs(1)=stwbs(1)+j2
      if(stwbs(1).ge.swbuff)then
      goto 85
      endif
      j1=stib(drecp(0)+i1)
      stwb(ip:stwbs(1))=stcb(j1:j1+j2-1)
      if(cflag(3).gt.0)then
      if(ichar(stwb(stwbs(1):stwbs(1))).eq.alf)then
      stwbs(1)=stwbs(1)+1
      stwb(stwbs(1)-1:stwbs(1))=char(acr)//char(alf)
      else
      goto 80
      endif
      endif
      enddo
      else
      goto 80
      endif
      elseif(igk.eq.-6)then
      if((igc.le.0).or.(igc.gt.nudk))then
      goto 80
      endif
      if(stib(udkt(0)+igc).eq.1)then
      ij=stib(udki(0)+igc)
      if((ij.le.0).or.(ij.gt.ngmk))then
      goto 80
      endif
      j1=stib(gmkvp(0)+ij)
      j2=stib(gmkvl(0)+ij)
      if(j2.gt.0)then
      ip=stwbs(1)+1
      stwbs(1)=stwbs(1)+j2
      if(stwbs(1).ge.swbuff)then
      goto 85
      endif
      stwb(ip:stwbs(1))=stcb(j1:j1+j2-1)
      elseif(j2.lt.0)then
      goto 80
      endif
      else
      goto 80
      endif
      else
      goto 80
      endif
      ig=ig+stib(ig+2)
      goto 20
   30 continue
      if(lupt.ne.0)then
      goto 80
      endif
      if(stwbs(1).gt.0)then
      call qout
      stwbs(1)=0
      endif
      if(sec.ne.1)then
      endfile(unit=aunit(nfiles))
      call qclose(nfiles,0)
      endif
      goto 90
   80 continue
      j1=stib(pske(0)+sec)
      j2=j1-1+stib(wske(0)+sec)
      mlin(1:srec)='run-time error while processing <'//stcb(j1:j2)//'>'
      call mput(1,0,0,0)
   85 continue
      call uput(21)
   90 return
      end
      subroutine umpi(xx,situ)
      implicit integer(a-z)
      save
      parameter ( maxleg=10, maxrho=7, maxi=10, maxdeg=6 )
      parameter ( maxn=maxi+maxi-2, maxli=maxn+maxrho )
      parameter ( srec=81, ssrec=62 )
      common/z1g/g(1:maxn,1:maxn),rho(1:maxdeg),nrho,rhop1,rhop2
      common/z4g/n,nli
      common/z17g/xtail(1:maxn),xhead(1:maxn),ntadp
      character*(srec) mlin
      common/z22g/mlin
      integer aa(1:maxn),bb(1:maxn)
      ntadp=0
      situ=-1
      do i1=rhop1,n
      do i2=i1+1,n
      if(g(i1,i2).eq.1)then
      g(i1,i2)=0
      kk=1
      bb(1)=1
      aa(1)=1
      do i4=2,n
      aa(i4)=0
      enddo
      do i3=1,n
      if(kk.eq.n)then
      goto 20
      endif
      if(i3.gt.kk)then
      g(i1,i2)=1
      if(xx.eq.1)then
      goto 90
      endif
      ii=0
      do i4=1,rho(1)
      ii=ii+aa(i4)
      enddo
      if(xx.eq.2)then
      if((ii.eq.0).or.(ii.eq.rho(1)))then
      goto 90
      endif
      elseif(xx.eq.3)then
      if((ii.eq.0).or.(ii.eq.rho(1)))then
      ntadp=ntadp+1
      xtail(ntadp)=i1
      xhead(ntadp)=i2
      endif
      elseif(xx.eq.4)then
      if((ii.eq.1).or.(ii.eq.rho(1)-1))then
      goto 90
      endif
      endif
      goto 20
      endif
      jj=bb(i3)
      do i4=1,jj-1
      if(aa(i4).eq.0)then
      if(g(i4,jj).gt.0)then
      kk=kk+1
      bb(kk)=i4
      aa(i4)=1
      endif
      endif
      enddo
      do i4=jj+1,n
      if(aa(i4).eq.0)then
      if(g(jj,i4).gt.0)then
      kk=kk+1
      bb(kk)=i4
      aa(i4)=1
      endif
      endif
      enddo
      enddo
      mlin(1:srec)='umpi_1'
      call mput(1,0,0,0)
   20 continue
      g(i1,i2)=1
      endif
      enddo
      enddo
      situ=1
   90 return
      end
      subroutine umvi(xx,situ)
      implicit integer(a-z)
      save
      parameter ( maxleg=10, maxrho=7, maxi=10, maxdeg=6 )
      parameter ( maxn=maxi+maxi-2, maxli=maxn+maxrho )
      parameter ( srec=81, ssrec=62 )
      common/z1in/leg(1:maxleg),nleg,nloop,incom
      common/z1g/g(1:maxn,1:maxn),rho(1:maxdeg),nrho,rhop1,rhop2
      common/z4g/n,nli
      character*(srec) mlin
      common/z22g/mlin
      integer aa(1:maxn),bb(1:maxn),cc(1:maxn)
      situ=1
      if(xx.eq.1)then
      if(nloop.eq.0)then
      goto 90
      endif
      ii=0
      if(rho(1).ne.1)then
      ii=1
      elseif(nloop.gt.1)then
      ii=1
      endif
      if(ii.ne.0)then
      do i1=rhop1,n
      if(g(i1,i1).ne.0)then
      goto 80
      endif
      enddo
      endif
      elseif(xx.eq.2)then
      if(n-1-rhop1.le.0)then
      goto 90
      endif
      elseif(xx.eq.3)then
      if(nloop.eq.0)then
      goto 90
      elseif(rho(1).eq.0)then
      goto 90
      endif
      goto 30
      else
      mlin(1:srec)='umvi_0'
      call mput(1,0,0,0)
      endif
      do i1=rhop1,n
      if(xx.eq.1)then
      kk=rho(1)
      do i2=1,rho(1)
      aa(i2)=1
      bb(i2)=i2
      enddo
      do i2=rhop1,n
      aa(i2)=0
      enddo
      k1=n
      elseif(xx.eq.2)then
      kk=1
      do i2=rhop1,n
      aa(i2)=0
      enddo
      if(i1.ne.rhop1)then
      bb(1)=rhop1
      elseif(i1.ne.n)then
      bb(1)=n
      else
      goto 90
      endif
      aa(bb(1))=1
      k1=n-rho(1)
      endif
      do i2=1,k1
      if(kk.eq.k1)then
      goto 20
      endif
      if(i2.gt.kk)then
      goto 80
      endif
      j1=bb(i2)
      if(j1.ne.i1)then
      do i3=rhop1,j1-1
      if(aa(i3).eq.0)then
      if(g(i3,j1).gt.0)then
      kk=kk+1
      bb(kk)=i3
      aa(i3)=1
      endif
      endif
      enddo
      do i3=j1+1,n
      if(aa(i3).eq.0)then
      if(g(j1,i3).gt.0)then
      kk=kk+1
      bb(kk)=i3
      aa(i3)=1
      endif
      endif
      enddo
      endif
      enddo
      mlin(1:srec)='umvi_1'
      call mput(1,0,0,0)
   20 continue
      enddo
      goto 90
   30 continue
      do i1=rhop1,n
      do i2=1,n
      aa(i2)=0
      enddo
      jj=0
      kk=1
      aa(i1)=-1
      bb(1)=i1
      do i2=1,n
      if(kk.eq.n)then
      goto 45
      elseif(aa(i2).ne.0)then
      goto 40
      endif
      jj=jj+1
      aa(i2)=jj
      kk=kk+1
      bb(kk)=i2
      ii=kk
   35 continue
      j1=bb(ii)
      do i3=i2+1,j1-1
      if(aa(i3).eq.0)then
      if(g(i3,j1).gt.0)then
      kk=kk+1
      bb(kk)=i3
      aa(i3)=jj
      endif
      endif
      enddo
      do i3=j1+1,n
      if(aa(i3).eq.0)then
      if(g(j1,i3).gt.0)then
      kk=kk+1
      bb(kk)=i3
      aa(i3)=jj
      endif
      endif
      enddo
      if(ii.lt.kk)then
      if(kk.lt.n)then
      ii=ii+1
      goto 35
      endif
      endif
   40 continue
      enddo
   45 continue
      if(kk.ne.n)then
      mlin(1:srec)='umvi_2'
      call mput(1,0,0,0)
      endif
      do i2=1,jj
      cc(i2)=0
      enddo
      do i2=1,rho(1)
      cc(aa(i2))=cc(aa(i2))+1
      enddo
      do i2=1,jj
      if(cc(i2).eq.1)then
      goto 50
      endif
      enddo
      goto 70
   50 continue
      j0=0
      do i2=1,jj
      bb(i2)=cc(i2)
      if(cc(i2).eq.0)then
      j0=j0+1
      endif
      enddo
      do i2=rhop1,n
      if(i2.ne.i1)then
      bb(aa(i2))=bb(aa(i2))+1
      endif
      enddo
      do i2=1,jj
      if(cc(i2).eq.1)then
      j1=bb(i2)
      j2=g(i1,i1)
      j3=j0
      if(j1.eq.1)then
      if(j2.gt.1)then
      j2=j2-2
      elseif(j0.gt.0)then
      i4=0
      do i3=1,jj
      if(cc(i3).eq.0)then
      if(i4.eq.0)then
      i4=i3
      elseif(bb(i3).lt.bb(i4))then
      i4=i3
      endif
      endif
      enddo
      j3=j3-1
      j1=j1+bb(i4)
      else
      goto 55
      endif
      endif
      if(0.lt.j2+j3)then
      goto 80
      elseif(j1+2.lt.n)then
      goto 80
      endif
      endif
   55 continue
      enddo
   70 continue
      enddo
      goto 90
   80 continue
      situ=-1
   90 return
      end
      subroutine dkar(im,iw)
      implicit integer(a-z)
      save
      parameter ( sibuff=1048574, scbuff=131064, swbuff=65528 )
      common/z15g/pten(1:10),iref,wiref,wsint
      character*(scbuff) stcb
      common/z31g/stcb
      common/z35g/stibs(1:1),stcbs(1:1),stwbs(1:1)
      common/z51g/aeq,agt,alt,aminus,anine,aplus,aslash,atimes,azero
      jc=stcbs(1)
      if(im.lt.0)then
      if(im+iref.lt.0)then
      goto 80
      endif
      xm=-im
      if(jc.ge.scbuff)then
      call vaocb(scbuff+1)
      endif
      jc=jc+1
      stcb(jc:jc)=char(aminus)
      else
      if(im-iref.gt.0)then
      goto 80
      endif
      xm=im
      endif
      ii=jc+1
   10 continue
      if(jc.ge.scbuff)then
      call vaocb(scbuff+1)
      endif
      jc=jc+1
      if(xm.lt.10)then
      stcb(jc:jc)=char(xm+azero)
      else
      ym=xm/10
      stcb(jc:jc)=char((xm-10*ym)+azero)
      xm=ym
      goto 10
      endif
      jj=jc
   20 continue
      if(ii.lt.jj)then
      kk=ichar(stcb(ii:ii))
      stcb(ii:ii)=stcb(jj:jj)
      stcb(jj:jj)=char(kk)
      ii=ii+1
      jj=jj-1
      goto 20
      endif
      goto 90
   80 continue
      call wput(-16,0,0)
   90 continue
      iw=jc-stcbs(1)
      return
      end
      subroutine sdiag(ii,jj)
      implicit integer(a-z)
      save
      parameter ( srec=81, ssrec=62 )
      parameter ( sibuff=1048574, scbuff=131064, swbuff=65528 )
      common/z15g/pten(1:10),iref,wiref,wsint
      character*(srec) mlin
      common/z22g/mlin
      character*(scbuff) stcb
      common/z31g/stcb
      common/z47g/ndiagp,ndiagl,hhp,hhl,doffp,doffl,noffp,noffl
      common/z51g/aeq,agt,alt,aminus,anine,aplus,aslash,atimes,azero
      if(ii.eq.1)then
      xp=hhp
      xl=hhl
      elseif(ii.eq.2)then
      xp=noffp
      xl=noffl
      elseif(ii.eq.3)then
      xp=ndiagp
      xl=ndiagl
      else
      mlin(1:srec)='sdiag_1'
      call mput(1,0,0,0)
      endif
      if(jj.eq.0)then
      do i1=xp+xl,xp+1,-1
      j1=ichar(stcb(i1:i1))
      if(j1.lt.anine)then
      stcb(i1:i1)=char(j1+1)
      goto 90
      else
      stcb(i1:i1)=char(azero)
      endif
      enddo
      if(xl.ge.wsint)then
      mlin(1:srec)='sdiag_2'
      call mput(1,0,0,0)
      endif
      j1=xp+1
      stcb(j1:j1)=char(azero+1)
      j1=j1+xl
      stcb(j1:j1)=char(azero)
      xl=xl+1
      elseif(jj.eq.-1)then
      j1=xp+1
      stcb(j1:j1)=char(azero)
      xl=1
      else
      mlin(1:srec)='sdiag_3'
      call mput(1,0,0,0)
      endif
      if(ii.eq.1)then
      hhl=xl
      elseif(ii.eq.2)then
      noffl=xl
      elseif(ii.eq.3)then
      ndiagl=xl
      endif
   90 return
      end
      integer function stoz(s,j1,j2)
      implicit integer(a-z)
      save
      parameter ( srec=81, ssrec=62 )
      common/z15g/pten(1:10),iref,wiref,wsint
      character*(srec) mlin
      common/z22g/mlin
      common/z51g/aeq,agt,alt,aminus,anine,aplus,aslash,atimes,azero
      common/z52g/acf(0:127),abo(0:6),acr,alf,aspace,dquote,squote
      character*(*) s
      stoz=0
      if((j1.le.0).or.(j1.gt.j2))then
      goto 80
      endif
      is=ichar(s(j1:j1))
      if((is.eq.aplus).or.(is.eq.aminus))then
      i1=j1+1
      else
      i1=j1
      endif
   10 continue
      jj=ichar(s(i1:i1))
      if(acf(jj).ne.1)then
      goto 80
      endif
      stoz=10*stoz+(jj-azero)
      if(stoz-iref.gt.0)then
      call wput(-16,0,0)
      endif
      if(i1.lt.j2)then
      i1=i1+1
      goto 10
      endif
      if(is.eq.aminus)then
      stoz=-stoz
      endif
      goto 90
   80 continue
      mlin(1:srec)='stoz_1'
      call mput(1,0,0,0)
   90 return
      end
      integer function sigz(s,ia,ib)
      implicit integer(a-z)
      save
      parameter ( srec=81, ssrec=62 )
      character*(srec) mlin
      common/z22g/mlin
      common/z51g/aeq,agt,alt,aminus,anine,aplus,aslash,atimes,azero
      common/z52g/acf(0:127),abo(0:6),acr,alf,aspace,dquote,squote
      character*(*) s
      if(ia.le.0)then
      goto 30
      endif
      jj=ichar(s(ia:ia))
      if(jj.eq.aplus)then
      ii=1
      elseif(jj.eq.aminus)then
      ii=-1
      else
      ii=0
      endif
      ja=ia+ii*ii
      if(ja.gt.ib)then
      goto 30
      endif
      sigz=0
      do i1=ja,ib
      jj=ichar(s(i1:i1))
      if(acf(jj).ne.1)then
      goto 30
      elseif(sigz.eq.0)then
      if(jj.ne.azero)then
      if(ii.lt.0)then
      sigz=-1
      else
      sigz=1
      endif
      endif
      endif
      enddo
      goto 90
   30 continue
      mlin(1:srec)='sigz_1'
      call mput(1,0,0,0)
   90 return
      end
      integer function wztos(nn)
      implicit integer(a-z)
      save
      common/z15g/pten(1:10),iref,wiref,wsint
      if(nn.ge.0)then
      if(nn-iref.gt.0)then
      goto 80
      endif
      mm=nn
      wztos=0
      else
      if(nn+iref.lt.0)then
      goto 80
      endif
      mm=-nn
      wztos=1
      endif
      do i1=1,9
      if(mm.lt.pten(i1))then
      wztos=wztos+i1
      goto 90
      endif
      enddo
   80 continue
      call wput(-16,0,0)
   90 return
      end
      subroutine spak(s,ind,m,uc,nos)
      implicit integer(a-z)
      save
      parameter ( maxleg=10, maxrho=7, maxi=10, maxdeg=6 )
      parameter ( srec=81, ssrec=62 )
      parameter ( sibuff=1048574, scbuff=131064, swbuff=65528 )
      common/z2in/momep(0:maxleg),momel(0:maxleg),kpqs(1:4)
      character*(srec) mlin
      common/z22g/mlin
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z35g/stibs(1:1),stcbs(1:1),stwbs(1:1)
      common/z48g/acomma,ascol,albra,arbra,alpar,arpar
      common/z51g/aeq,agt,alt,aminus,anine,aplus,aslash,atimes,azero
      common/z52g/acf(0:127),abo(0:6),acr,alf,aspace,dquote,squote
      character*(*) s
      if((m.ne.0).and.(m.ne.1))then
      goto 08
      elseif((nos.ne.0).and.(nos.ne.1))then
      goto 08
      elseif((uc.ne.0).and.(uc.ne.1))then
      goto 08
      endif
      inos=nos
      iuc=uc
      i1=0
      sl=1
      ii=2*srec+2
   03 continue
      if(ichar(s(sl:sl)).ne.ascol)then
      if(sl.lt.ii)then
      if((i1.eq.0).and.(ichar(s(sl:sl)).eq.acomma))then
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
      if((i1.lt.3).or.(sl.lt.3))then
      goto 08
      elseif(mod(i1,2).eq.0)then
      goto 08
      endif
      sp=stcbs(1)+1
      i1=(i1-1)/2
      call aocb(i1+1)
      ij=0
      j2=1
      k=sp
      do i2=1,i1
      ii=ichar(s(j2:j2))
      jj=ichar(s(j2+1:j2+1))
      if((ii.lt.azero).or.(ii.gt.anine))then
      goto 08
      elseif((jj.lt.azero).or.(jj.gt.anine))then
      goto 08
      endif
      kk=10*ii+jj-496
      if((iuc.eq.1).and.(kk.ge.abo(3)).and.(kk.le.abo(4)))then
      goto 08
      endif
      if((kk.ge.abo(5)).and.(kk.le.abo(6)))then
      ij=1
      endif
      stcb(k:k)=char(kk)
      j2=j2+2
      k=k+1
      enddo
      stcb(k:k)=char(alf)
      if((uic.eq.1).and.(ij.eq.0))then
      goto 08
      endif
      if(inos.eq.1)then
      if(ichar(stcb(k-1:k-1)).ne.kpqs(4))then
      inos=0
      endif
      endif
      j0=stibs(1)
      if(m.eq.0)then
      call aoib(2)
      else
      call vaoib(2)
      endif
      stib(j0+1)=sp
      stib(j0+2)=i1
      if(m.ne.0)then
      if((nos.ne.0).or.(uc.ne.0))then
      goto 08
      endif
      goto 90
      endif
      j2=i1+i1+1
      j1=j2+1
   11 continue
      if(j2.lt.sl)then
      j2=j2+1
      j3=ichar(s(j2:j2))
      if((j3.eq.acomma).or.(j3.eq.ascol))then
      j2=j2-1
      if(j2.lt.j1)then
      goto 08
      endif
      call aoib(1)
      stib(stibs(1))=stoz(s,j1,j2)
      j2=j2+1
      j1=j2+1
      endif
      goto 11
      endif
      ind=ind+1
      if(iuc.eq.1)then
      ii=stcbs(1)
      call aocb(i1+1)
      i3=ii+1
      do i2=ii-i1,ii-1
      j2=ichar(stcb(i2:i2))
      if((j2.ge.abo(5)).and.(j2.le.abo(6)))then
      j2=j2+abo(0)
      endif
      stcb(i3:i3)=char(j2)
      i3=i3+1
      enddo
      stcb(i3:i3)=char(alf)
      j1=stibs(1)
      jj=j1-j0
      if(jj.lt.2)then
      goto 08
      endif
      call aoib(jj)
      ii=j1+1
      stib(ii)=stib(ii-jj)+i1+1
      do i2=2,jj
      ii=ii+1
      stib(ii)=stib(ii-jj)
      enddo
      ind=ind+1
      endif
      if(inos.eq.1)then
      j1=stibs(1)
      jj=j1-j0
      if(jj.lt.4)then
      goto 08
      endif
      call aoib(jj)
      ii=j1
      do i2=1,jj
      ii=ii+1
      stib(ii)=stib(ii-jj)
      enddo
      ii=j1+2
      stib(ii)=stib(ii)-1
      ind=ind+1
      if(uc.eq.1)then
      ii=ii+jj/2
      stib(ii)=stib(ii)-1
      ind=ind+1
      endif
      endif
      goto 90
   08 continue
      mlin(1:srec)='spak_1'
      call mput(1,0,0,0)
   90 return
      end
      integer function stds(s,ia,ib,iw)
      implicit integer(a-z)
      save
      parameter ( sibuff=1048574, scbuff=131064, swbuff=65528 )
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z35g/stibs(1:1),stcbs(1:1),stwbs(1:1)
      common/z52g/acf(0:127),abo(0:6),acr,alf,aspace,dquote,squote
      character*(*) s
      stds=-1
      if((ia.le.0).or.(ia.ge.ib))then
      goto 90
      elseif((iw.lt.0).or.(iw.gt.2))then
      goto 90
      endif
      sc=0
      sl=0
      sf=0
      itl=1
      ii=0
      ll=0
      mm=0
      do i1=ia,ib
      jj=ichar(s(i1:i1))
      if(jj.eq.squote)then
      ii=ii+1
      mm=1-mm
      ll=0
      if(ii.gt.1)then
      kk=mm
      else
      kk=0
      endif
      elseif(jj.eq.alf)then
      if((mm.ne.0).or.(ll.eq.1))then
      goto 90
      endif
      sf=sf+1
      ll=1
      if(ii.gt.0)then
      ii=0
      sc=sc+1
      endif
      kk=0
      elseif(jj.eq.aspace)then
      if(ii.gt.0)then
      ii=0
      if(mm.eq.0)then
      sc=sc+1
      endif
      endif
      kk=mm
      elseif((jj.lt.abo(1)).or.(jj.gt.abo(2)))then
      goto 90
      else
      if(mm.ne.1)then
      goto 90
      endif
      ii=0
      kk=1
      endif
      if(kk.ne.0)then
      if(iw.eq.2)then
      if(jj.ne.aspace)then
      itl=0
      endif
      endif
      if((iw.eq.1).or.((iw.eq.2).and.(itl.eq.0)))then
      sl=sl+1
      call vaocb(sl)
      j1=stcbs(1)+sl
      stcb(j1:j1)=char(jj)
      endif
      endif
      enddo
      if(mm.ne.0)then
      goto 90
      endif
      if(ichar(s(ib:ib)).eq.squote)then
      sc=sc+1
      endif
      if((iw.eq.2).and.(sl.gt.1))then
      do i1=stcbs(1)+sl,stcbs(1)+1,-1
      if(ichar(stcb(i1:i1)).eq.aspace)then
      sl=sl-1
      else
      goto 20
      endif
      enddo
      endif
   20 continue
      stds=sl
   90 return
      end
      subroutine mstr0(s,ia,ib,pp,pl,ind)
      implicit integer(a-z)
      save
      parameter ( srec=81, ssrec=62 )
      parameter ( eoa=-2047, nap=-8191 )
      parameter ( sibuff=1048574, scbuff=131064, swbuff=65528 )
      character*(srec) mlin
      common/z22g/mlin
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      character*(*) s
      ind=0
      if(pp.eq.nap)then
      goto 90
      endif
      ls=ib-ia+1
   10 continue
      ind=ind+1
      ii=stib(pl+ind)
      if(ii.ne.eoa)then
      if(ii.eq.ls)then
      jj=stib(pp+ind)
      if(stcb(jj:jj-1+ls).eq.s(ia:ib))then
      goto 90
      endif
      endif
      goto 10
      endif
      ind=0
   90 continue
      if(ind.lt.0)then
      mlin(1:srec)='mstr0_1'
      call mput(1,0,0,0)
      endif
      return
      end
      subroutine mstr1(s,ia,ib,pp,kp,kl,ind)
      implicit integer(a-z)
      save
      parameter ( srec=81, ssrec=62 )
      parameter ( eoa=-2047, nap=-8191 )
      parameter ( sibuff=1048574, scbuff=131064, swbuff=65528 )
      character*(srec) mlin
      common/z22g/mlin
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      character*(*) s
      ind=0
      if(pp.eq.nap)then
      goto 90
      endif
      ls=ib-ia+1
      p1=pp
   10 continue
      p1=stib(p1)
      if(p1.eq.eoa)then
      ind=0
      goto 90
      endif
      ind=ind+1
      if(stib(p1+kl).eq.ls)then
      jj=stib(p1+kp)
      if(stcb(jj:jj-1+ls).eq.s(ia:ib))then
      goto 90
      endif
      endif
      goto 10
   90 continue
      if(ind.lt.0)then
      mlin(1:srec)='mstr1_1'
      call mput(1,0,0,0)
      endif
      return
      end
      integer function stdw(s,ia,ib)
      implicit integer(a-z)
      save
      parameter ( srec=81, ssrec=62 )
      character*(srec) mlin
      common/z22g/mlin
      common/z52g/acf(0:127),abo(0:6),acr,alf,aspace,dquote,squote
      character*(*) s
      if((ia.le.0).or.(ia.gt.ib))then
      mlin(1:srec)='stdw_1'
      call mput(1,0,0,0)
      endif
      if(acf(ichar(s(ia:ia))).lt.2)then
      stdw=0
      goto 90
      else
      do i1=ia+1,ib
      if(acf(ichar(s(i1:i1))).lt.0)then
      stdw=0
      goto 90
      endif
      enddo
      endif
      stdw=1
   90 return
      end
      integer function stdz(s,ia,ib)
      implicit integer(a-z)
      save
      parameter ( srec=81, ssrec=62 )
      character*(srec) mlin
      common/z22g/mlin
      common/z51g/aeq,agt,alt,aminus,anine,aplus,aslash,atimes,azero
      common/z52g/acf(0:127),abo(0:6),acr,alf,aspace,dquote,squote
      character*(*) s
      if((ia.le.0).or.(ia.gt.ib))then
      mlin(1:srec)='stdz_1'
      call mput(1,0,0,0)
      endif
      stdz=0
      jj=ichar(s(ia:ia))
      if((jj.eq.aplus).or.(jj.eq.aminus))then
      if(ia.lt.ib)then
      stdz=1
      endif
      elseif(acf(jj).eq.1)then
      stdz=1
      endif
      if(stdz.eq.1)then
      do i1=ia+1,ib
      if(acf(ichar(s(i1:i1))).ne.1)then
      stdz=0
      goto 90
      endif
      enddo
      endif
   90 return
      end
      integer function stdq(s,ia,ib)
      implicit integer(a-z)
      save
      parameter ( srec=81, ssrec=62 )
      character*(srec) mlin
      common/z22g/mlin
      common/z51g/aeq,agt,alt,aminus,anine,aplus,aslash,atimes,azero
      common/z52g/acf(0:127),abo(0:6),acr,alf,aspace,dquote,squote
      character*(*) s
      if((ia.le.0).or.(ia.gt.ib))then
      mlin(1:srec)='stdq_1'
      call mput(1,0,0,0)
      endif
      xsla=0
      xnum=0
      xden=0
      jj=ichar(s(ia:ia))
      if((jj.eq.aplus).or.(jj.eq.aminus))then
      stdq=1
      elseif(acf(jj).eq.1)then
      stdq=1
      xnum=1
      else
      stdq=0
      endif
      if(stdq.eq.1)then
      do i1=ia+1,ib
      jj=ichar(s(i1:i1))
      if(acf(jj).eq.1)then
      if(xsla.eq.0)then
      xnum=1
      elseif(jj.ne.azero)then
      xden=1
      endif
      elseif(jj.eq.aslash)then
      if(xsla.ne.0)then
      stdq=0
      goto 90
      endif
      xsla=1
      else
      stdq=0
      goto 90
      endif
      enddo
      if(xnum.ne.1)then
      stdq=0
      elseif(xsla.ne.0)then
      if(xden.ne.1)then
      stdq=0
      endif
      endif
      endif
   90 return
      end
      subroutine compac
      implicit integer(a-z)
      save
      parameter ( maxleg=10, maxrho=7, maxi=10, maxdeg=6 )
      parameter ( maxn=maxi+maxi-2, maxli=maxn+maxrho )
      parameter ( srec=81, ssrec=62 )
      parameter ( eoa=-2047, nap=-8191 )
      parameter ( nfiles=5 )
      parameter ( sibuff=1048574, scbuff=131064, swbuff=65528 )
      parameter ( svbuff=32768 )
      common/z1in/leg(1:maxleg),nleg,nloop,incom
      common/z2in/momep(0:maxleg),momel(0:maxleg),kpqs(1:4)
      common/z7in/aunit(1:nfiles)
      common/z2g/dis,dsym
      common/z3g/nivd(0:0),dpntro(0:0),vparto(0:0),vval(0:0),nvert
      common/z4g/n,nli
      common/z6g/p1(1:maxleg),invp1(1:maxleg)
      common/z7g/lmap(1:maxn,1:maxdeg),vmap(1:maxn,1:maxdeg),
     :pmap(1:maxn,1:maxdeg),vlis(1:maxn),invlis(1:maxn)
      common/z8g/vdeg(1:maxn),xn(1:maxn)
      common/z11g/nphi,nblok,nprop,npprop
      common/z12g/cflag(1:7),dflag(1:20),jflag(1:10),mflag(1:6)
      common/z16g/rdeg(1:maxn),amap(1:maxn,1:maxdeg)
      common/z18g/eg(1:maxn,1:maxn),flow(1:maxli,0:maxleg+maxrho),
     :net(-3:3)
      common/z19g/vfo(1:maxn)
      character*(srec) mlin
      common/z22g/mlin
      common/z24g/iogp(1:4)
      common/z25g/ex(1:maxli),ey(1:maxli),ovm(1:maxn,1:maxdeg)
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z33g/namep(0:0),namel(0:0)
      common/z34g/link(0:0),antiq(0:0),blok(0:0)
      common/z35g/stibs(1:1),stcbs(1:1),stwbs(1:1)
      common/z36g/udkp(0:2),udkl(0:0),udkt(0:0),udki(0:0)
      common/z38g/gmkp(0:0),gmkl(0:0),gmkd(0:0),gmko(0:0),
     :gmkvp(0:0),gmkvl(0:0)
      common/z39g/pmkp(0:0),pmkl(0:0),pmkd(0:0),pmkvpp(0:0),pmkvlp(0:0)
      common/z40g/vmkp(0:0),vmkl(0:0),vmkvpp(0:0),vmkvlp(0:0),vmks(0:0)
      common/z41g/nudk,ngmk,npmk0,npmk,nvmk0,nvmk
      common/z47g/ndiagp,ndiagl,hhp,hhl,doffp,doffl,noffp,noffl
      common/z51g/aeq,agt,alt,aminus,anine,aplus,aslash,atimes,azero
      common/z52g/acf(0:127),abo(0:6),acr,alf,aspace,dquote,squote
      character*(swbuff) stwb
      common/z54g/stwb
   01 continue
      ig=iogp(2)
      if((ig.lt.iogp(2)).or.(ig.ge.iogp(3)))then
      goto 80
      endif
      kli=stwbs(1)
      kma=-1
      lupt=0
   10 continue
      igk=stib(ig)
      if(igk.gt.0)then
      kli=kli+1
      if(cflag(3).gt.0)then
      if(igk.eq.alf)then
      stwb(kli:kli)=char(acr)
      kli=kli+1
      endif
      endif
      if(kli.ge.swbuff)then
      goto 75
      endif
      stwb(kli:kli)=char(igk)
      ig=ig+1
      goto 10
      elseif(igk.eq.eoa)then
      goto 20
      endif
      igc=stib(ig+1)
      if(igk.gt.-3)then
      if(igk.eq.-1)then
      lupt=igc
      if(lupt.gt.0)then
      if(stib(ig+4).eq.0)then
      if(lupt.eq.1)then
      stib(ig+4)=incom
      stib(ig+5)=0
      elseif(lupt.eq.2)then
      stib(ig+4)=nleg
      stib(ig+5)=incom
      elseif(lupt.eq.3)then
      stib(ig+4)=nli-nleg
      stib(ig+5)=0
      elseif(lupt.eq.4)then
      stib(ig+4)=n-nleg
      stib(ig+5)=0
      elseif(lupt.eq.5)then
      stib(ig+4)=vdeg(nleg+lupi)
      stib(ig+5)=0
      endif
      endif
      stib(ig+5)=stib(ig+5)+1
      if(stib(ig+4).lt.stib(ig+5))then
      stib(ig+4)=0
      if(lupt.eq.5)then
      lupt=4
      else
      lupt=0
      endif
      ig=ig+1
      else
      if(lupt.eq.5)then
      lupj=stib(ig+5)
      else
      lupi=stib(ig+5)
      endif
      endif
      elseif(lupt.ne.-1)then
      goto 80
      endif
      elseif(igk.eq.-2)then
      if(igc.eq.22)then
      if(jflag(10).eq.0)then
      if(kma.ge.0)then
      kli=kma
      kma=-1
      endif
      elseif(kli.gt.0)then
      kli=kli-1
      if(ichar(stwb(kli:kli)).eq.acr)then
      kli=kli-1
      endif
      else
      goto 80
      endif
      else
      goto 80
      endif
      else
      goto 80
      endif
      elseif(igk.gt.-5)then
      if(igk.eq.-3)then
      if((igc.gt.70).and.(igc.lt.80))then
      if(igc.lt.74)then
      if(igc.eq.72)then
      jj=dsym
      if(dsym.gt.1)then
      kli=kli+2
      if(kli.ge.swbuff)then
      goto 75
      endif
      stwb(kli-1:kli)=char(azero+1)//char(aslash)
      endif
      elseif(igc.eq.73)then
      jj=dsym
      endif
      elseif(igc.lt.77)then
      if(igc.eq.74)then
      jj=nli-nleg
      elseif(igc.eq.75)then
      jj=nleg
      elseif(igc.eq.76)then
      jj=nloop
      endif
      else
      if(igc.eq.77)then
      jj=n-nleg
      elseif(igc.eq.78)then
      jj=incom
      elseif(igc.eq.79)then
      jj=nleg-incom
      endif
      endif
      if(igc.ne.71)then
      call dkar(jj,jk)
      ip=kli+1
      kli=kli+jk
      if(kli.ge.swbuff)then
      goto 75
      endif
      stwb(ip:kli)=stcb(stcbs(1)+1:stcbs(1)+jk)
      else
      ip=kli+1
      kli=kli+noffl
      if(kli.ge.swbuff)then
      goto 75
      endif
      stwb(ip:kli)=stcb(noffp+1:noffp+noffl)
      endif
      else
      if(igc.eq.61)then
      kli=kli+1
      if(kli.ge.swbuff)then
      goto 75
      endif
      if(dis.gt.0)then
      stwb(kli:kli)=char(aplus)
      elseif(dis.lt.0)then
      stwb(kli:kli)=char(aminus)
      else
      goto 80
      endif
      elseif(igc.eq.62)then
      if(dis.lt.0)then
      kli=kli+1
      if(kli.ge.swbuff)then
      goto 75
      endif
      stwb(kli:kli)=char(aminus)
      endif
      else
      goto 80
      endif
      endif
      elseif(igk.eq.-4)then
      if(lupt.eq.3)then
      if(igc.lt.40)then
      if((igc.eq.31).or.(igc.eq.33))then
      i=lupi
      j=pmap(ex(i),ey(i))
      if(igc.eq.33)then
      j=stib(link(0)+j)
      endif
      il=stib(namel(0)+j)
      ip=kli+1
      kli=kli+il
      if(kli.ge.swbuff)then
      goto 75
      endif
      ia=stib(namep(0)+j)
      stwb(ip:kli)=stcb(ia:ia-1+il)
      elseif((igc.eq.32).or.(igc.eq.34))then
      i=lupi
      if(igc.eq.32)then
      ik=1
      else
      ik=0
      endif
      j=ey(i)
      i=ex(i)
      k0=kli
      other=0
      if(vmap(i,j).eq.i)then
      if(mod(j-rdeg(i)+ik,2).ne.0)then
      kli=kli+1
      if(kli.ge.swbuff)then
      goto 75
      endif
      stwb(kli:kli)=char(aminus)
      endif
      ip=kli+1
      kli=kli+momel(0)
      if(kli.ge.swbuff)then
      goto 75
      endif
      stwb(ip:kli)=stcb(momep(0):momep(0)-1+momel(0))
      do k=nleg+1,nleg+nloop
      if(flow(amap(i,j),k).ne.0)then
      call dkar(k-nleg,jk)
      ip=kli+1
      kli=kli+jk
      if(kli.ge.swbuff)then
      goto 75
      endif
      stwb(ip:kli)=stcb(stcbs(1)+1:stcbs(1)+jk)
      goto 70
      endif
      enddo
      else
      if(vmap(i,j).lt.i)then
      ij=-1
      else
      ij=1
      endif
      if(ik.eq.0)then
      ij=-ij
      endif
      do k=nleg+1,nleg+nloop
      jj=flow(amap(i,j),k)
      if(jj.ne.0)then
      if((other.eq.1).or.(ij.eq.jj))then
      kli=kli+1
      if(kli.ge.swbuff)then
      goto 75
      endif
      if(ij.eq.jj)then
      stwb(kli:kli)=char(aminus)
      else
      stwb(kli:kli)=char(aplus)
      endif
      endif
      ip=kli+1
      kli=kli+momel(0)
      if(kli.ge.swbuff)then
      goto 75
      endif
      stwb(ip:kli)=stcb(momep(0):momep(0)-1+momel(0))
      call dkar(k-nleg,jk)
      ip=kli+1
      kli=kli+jk
      if(kli.ge.swbuff)then
      goto 75
      endif
      stwb(ip:kli)=stcb(stcbs(1)+1:stcbs(1)+jk)
      other=1
      endif
      enddo
      do k=1,nleg
      jj=flow(amap(i,j),invp1(k))
      if(jj.ne.0)then
      if(k.gt.incom)then
      jj=-jj
      endif
      if((other.eq.1).or.(ij.eq.jj))then
      kli=kli+1
      if(kli.ge.swbuff)then
      goto 75
      endif
      if(ij.eq.jj)then
      stwb(kli:kli)=char(aminus)
      else
      stwb(kli:kli)=char(aplus)
      endif
      endif
      ip=kli+1
      kli=kli+momel(k)
      if(kli.ge.swbuff)then
      goto 75
      endif
      stwb(ip:kli)=stcb(momep(k):momep(k)-1+momel(k))
      other=1
      endif
      enddo
      if(k0.eq.kli)then
      kli=kli+1
      if(kli.ge.swbuff)then
      goto 75
      endif
      stwb(kli:kli)=char(azero)
      endif
      endif
      elseif(igc.eq.35)then
      i=lupi
      j=pmap(ex(i),ey(i))
      kli=kli+1
      if(kli.ge.swbuff)then
      goto 75
      endif
      if(stib(antiq(0)+j).eq.0)then
      stwb(kli:kli)=char(aplus)
      else
      stwb(kli:kli)=char(aminus)
      endif
      else
      goto 80
      endif
      elseif(igc.lt.50)then
      if(igc.lt.44)then
      if(igc.eq.40)then
      jj=3
      elseif(igc.eq.41)then
      jj=lupi
      elseif(igc.eq.42)then
      jj=2*lupi-1
      elseif(igc.eq.43)then
      i=lupi
      j=ey(i)
      i=ex(i)
      jj=0
      do k=1,vdeg(i)
      if(ovm(i,k).eq.j)then
      jj=k
      endif
      enddo
      else
      goto 80
      endif
      elseif(igc.lt.47)then
      if(igc.eq.44)then
      i=lupi
      jj=ex(i)-nleg
      elseif(igc.eq.45)then
      jj=2*lupi
      elseif(igc.eq.46)then
      i=lupi
      j=lmap(ex(i),ey(i))
      i=vmap(ex(i),ey(i))
      jj=0
      do k=1,vdeg(i)
      if(ovm(i,k).eq.j)then
      jj=k
      endif
      enddo
      else
      goto 80
      endif
      else
      if(igc.eq.47)then
      i=lupi
      jj=vmap(ex(i),ey(i))-nleg
      elseif(igc.eq.48)then
      i=lupi
      jj=vdeg(ex(i))
      elseif(igc.eq.49)then
      i=lupi
      jj=vdeg(vmap(ex(i),ey(i)))
      else
      goto 80
      endif
      endif
      call dkar(jj,jk)
      ip=kli+1
      kli=kli+jk
      if(kli.ge.swbuff)then
      goto 75
      endif
      stwb(ip:kli)=stcb(stcbs(1)+1:stcbs(1)+jk)
      else
      goto 80
      endif
      elseif((lupt.eq.1).or.(lupt.eq.2))then
      if(igc.lt.40)then
      if((igc.eq.31).or.(igc.eq.33))then
      i=lupi
      j=stib(link(0)+leg(i))
      if(igc.eq.33)then
      j=stib(link(0)+j)
      endif
      if(lupi.gt.incom)then
      j=stib(link(0)+j)
      endif
      il=stib(namel(0)+j)
      ip=kli+1
      kli=kli+il
      if(kli.ge.swbuff)then
      goto 75
      endif
      ia=stib(namep(0)+j)
      stwb(ip:kli)=stcb(ia:ia-1+il)
      elseif((igc.eq.32).or.(igc.eq.34))then
      i=lupi
      if(igc.eq.34)then
      kli=kli+1
      if(kli.ge.swbuff)then
      goto 75
      endif
      stwb(kli:kli)=char(aminus)
      endif
      ip=kli+1
      kli=kli+momel(i)
      if(kli.ge.swbuff)then
      goto 75
      endif
      stwb(ip:kli)=stcb(momep(i):momep(i)-1+momel(i))
      elseif(igc.eq.35)then
      i=lupi
      j=leg(i)
      kli=kli+1
      if(kli.ge.swbuff)then
      goto 75
      endif
      if(stib(antiq(0)+j).eq.0)then
      stwb(kli:kli)=char(aplus)
      else
      stwb(kli:kli)=char(aminus)
      endif
      else
      goto 80
      endif
      elseif(igc.lt.50)then
      if(igc.eq.40)then
      i=lupi
      if(i.le.incom)then
      jj=1
      else
      jj=2
      endif
      elseif(igc.eq.42)then
      i=lupi
      if(i.le.incom)then
      jj=1-2*i
      else
      jj=2*(incom-i)
      endif
      elseif(igc.eq.43)then
      i=lupi
      j=lmap(invp1(i),1)
      i=vmap(invp1(i),1)
      jj=0
      do k=1,vdeg(i)
      if(ovm(i,k).eq.j)then
      jj=k
      endif
      enddo
      elseif(igc.eq.44)then
      i=lupi
      jj=vmap(invp1(i),1)-nleg
      elseif(igc.eq.48)then
      i=lupi
      j=lmap(invp1(i),1)
      jj=vdeg(vmap(invp1(i),1))
      else
      goto 80
      endif
      call dkar(jj,jk)
      ip=kli+1
      kli=kli+jk
      if(kli.ge.swbuff)then
      goto 75
      endif
      stwb(ip:kli)=stcb(stcbs(1)+1:stcbs(1)+jk)
      elseif(igc.lt.60)then
      if(igc.eq.51)then
      i=lupi
      elseif(igc.eq.52)then
      i=lupi
      elseif(igc.eq.53)then
      i=lupi-incom
      else
      goto 80
      endif
      call dkar(i,jk)
      ip=kli+1
      kli=kli+jk
      if(kli.ge.swbuff)then
      goto 75
      endif
      stwb(ip:kli)=stcb(stcbs(1)+1:stcbs(1)+jk)
      else
      goto 80
      endif
      elseif((lupt.eq.4).or.(lupt.eq.5))then
      if(igc.lt.40)then
      if((igc.eq.31).or.(igc.eq.33))then
      i=nleg+lupi
      j=lupj
      j=pmap(i,ovm(i,j))
      if(igc.eq.33)then
      j=stib(link(0)+j)
      endif
      il=stib(namel(0)+j)
      ip=kli+1
      kli=kli+il
      if(kli.ge.swbuff)then
      goto 75
      endif
      ia=stib(namep(0)+j)
      stwb(ip:kli)=stcb(ia:ia-1+il)
      elseif((igc.eq.32).or.(igc.eq.34))then
      i=nleg+lupi
      j=lupj
      j=ovm(i,j)
      ij=vmap(i,j)
      if(igc.eq.32)then
      kk=0
      else
      kk=1
      endif
      k0=kli
      other=0
      if(ij.le.nleg)then
      ij=p1(ij)
      if(ij.gt.incom)then
      kk=1-kk
      endif
      if(kk.eq.1)then
      kli=kli+1
      if(kli.ge.swbuff)then
      goto 75
      endif
      stwb(kli:kli)=char(aminus)
      endif
      ip=kli+1
      kli=kli+momel(ij)
      if(kli.ge.swbuff)then
      goto 75
      endif
      stwb(ip:kli)=stcb(momep(ij):momep(ij)-1+momel(ij))
      elseif(vmap(i,j).eq.i)then
      if(mod(j-rdeg(i)+kk,2).eq.0)then
      kli=kli+1
      if(kli.ge.swbuff)then
      goto 75
      endif
      stwb(kli:kli)=char(aminus)
      endif
      ip=kli+1
      kli=kli+momel(0)
      if(kli.ge.swbuff)then
      goto 75
      endif
      stwb(ip:kli)=stcb(momep(0):momep(0)-1+momel(0))
      do k=nleg+1,nleg+nloop
      if(flow(amap(i,j),k).ne.0)then
      call dkar(k-nleg,jk)
      ip=kli+1
      kli=kli+jk
      if(kli.ge.swbuff)then
      goto 75
      endif
      stwb(ip:kli)=stcb(stcbs(1)+1:stcbs(1)+jk)
      goto 70
      endif
      enddo
      else
      ij=1
      if(vmap(i,j).lt.i)then
      ij=-ij
      endif
      if(kk.eq.1)then
      ij=-ij
      endif
      do k=nleg+1,nleg+nloop
      jj=flow(amap(i,j),k)
      if(jj.ne.0)then
      if((other.eq.1).or.(ij.eq.jj))then
      kli=kli+1
      if(kli.ge.swbuff)then
      goto 75
      endif
      if(ij.eq.jj)then
      stwb(kli:kli)=char(aminus)
      else
      stwb(kli:kli)=char(aplus)
      endif
      endif
      ip=kli+1
      kli=kli+momel(0)
      if(kli.ge.swbuff)then
      goto 75
      endif
      stwb(ip:kli)=stcb(momep(0):momep(0)-1+momel(0))
      call dkar(k-nleg,jk)
      ip=kli+1
      kli=kli+jk
      if(kli.ge.swbuff)then
      goto 75
      endif
      stwb(ip:kli)=stcb(stcbs(1)+1:stcbs(1)+jk)
      other=1
      endif
      enddo
      do k=1,nleg
      jj=flow(amap(i,j),invp1(k))
      if(jj.ne.0)then
      if(k.gt.incom)then
      jj=-jj
      endif
      if((other.eq.1).or.(ij.eq.jj))then
      kli=kli+1
      if(kli.ge.swbuff)then
      goto 75
      endif
      if(ij.eq.jj)then
      stwb(kli:kli)=char(aminus)
      else
      stwb(kli:kli)=char(aplus)
      endif
      endif
      ip=kli+1
      kli=kli+momel(k)
      if(kli.ge.swbuff)then
      goto 75
      endif
      stwb(ip:kli)=stcb(momep(k):momep(k)-1+momel(k))
      other=1
      endif
      enddo
      if(k0.eq.kli)then
      kli=kli+1
      if(kli.ge.swbuff)then
      goto 75
      endif
      stwb(kli:kli)=char(azero)
      endif
      endif
      elseif(igc.eq.35)then
      i=nleg+lupi
      j=lupj
      j=pmap(i,ovm(i,j))
      kli=kli+1
      if(kli.ge.swbuff)then
      goto 75
      endif
      if(stib(antiq(0)+j).eq.0)then
      stwb(kli:kli)=char(aplus)
      else
      stwb(kli:kli)=char(aminus)
      endif
      endif
      elseif(igc.lt.50)then
      if(igc.eq.41)then
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
      elseif(igc.lt.46)then
      if(igc.eq.40)then
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
      elseif((igc.eq.42).or.(igc.eq.45))then
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
      elseif(mod(j-rdeg(i),2).ne.0)then
      jj=2*(amap(i,j)-nleg)-1
      else
      jj=2*(amap(i,j)-nleg)
      endif
      if(igc.eq.45)then
      if(jj.gt.0)then
      jj=jj-1+2*mod(jj,2)
      else
      jj=0
      endif
      endif
      elseif(igc.eq.43)then
      jj=lupj
      elseif(igc.eq.44)then
      jj=lupi
      else
      goto 80
      endif
      else
      if(igc.eq.46)then
      i=nleg+lupi
      j=lupj
      ij=vmap(i,ovm(i,j))
      j=lmap(i,ovm(i,j))
      jj=0
      if(ij.gt.nleg)then
      do k=1,vdeg(ij)
      if(ovm(ij,k).eq.j)then
      jj=k
      endif
      enddo
      endif
      elseif(igc.eq.47)then
      i=nleg+lupi
      j=lupj
      jj=vmap(i,ovm(i,j))-nleg
      if(jj.lt.0)then
      jj=0
      endif
      elseif(igc.eq.48)then
      jj=vdeg(nleg+lupi)
      elseif(igc.eq.49)then
      i=nleg+lupi
      j=lupj
      jj=vmap(i,ovm(i,j))
      if(jj.le.nleg)then
      jj=0
      else
      jj=vdeg(jj)
      endif
      else
      goto 80
      endif
      endif
      call dkar(jj,jk)
      ip=kli+1
      kli=kli+jk
      if(kli.ge.swbuff)then
      goto 75
      endif
      stwb(ip:kli)=stcb(stcbs(1)+1:stcbs(1)+jk)
      else
      goto 80
      endif
      else
      goto 80
      endif
      else
      goto 80
      endif
      elseif(igk.eq.-6)then
      if((igc.le.0).or.(igc.gt.nudk))then
      goto 80
      endif
      igut=stib(udkt(0)+igc)
      if(igut.eq.1)then
      ij=stib(udki(0)+igc)
      if((ij.le.0).or.(ij.gt.ngmk))then
      goto 80
      endif
      j1=stib(gmkvp(0)+ij)
      j2=stib(gmkvl(0)+ij)
      if(j2.gt.0)then
      ip=kli+1
      kli=kli+j2
      if(kli.ge.swbuff)then
      goto 75
      endif
      stwb(ip:kli)=stcb(j1:j1-1+j2)
      elseif((j2.lt.0).or.(stib(ig+3).ne.0))then
      goto 80
      endif
      elseif(igut.eq.2)then
      if(lupt.eq.1)then
      i=lupi
      j=stib(link(0)+leg(i))
      elseif(lupt.eq.2)then
      i=lupi
      j=leg(i)
      elseif(lupt.eq.3)then
      i=lupi
      j=pmap(ex(i),ey(i))
      elseif(lupt.eq.4)then
      goto 80
      elseif(lupt.eq.5)then
      i=nleg+lupi
      j=lupj
      j=pmap(i,ovm(i,j))
      else
      goto 80
      endif
      if(stib(ig+3).ne.0)then
      j=stib(link(0)+j)
      endif
      ij=stib(udki(0)+igc)
      if((ij.le.0).or.(ij.gt.npmk))then
      goto 80
      elseif((j.le.0).or.(j.gt.nphi))then
      goto 80
      endif
      il=stib(stib(pmkvlp(0)+ij)+j)
      if(il.gt.0)then
      ip=kli+1
      kli=kli+il
      if(kli.ge.swbuff)then
      goto 75
      endif
      ia=stib(stib(pmkvpp(0)+ij)+j)
      stwb(ip:kli)=stcb(ia:ia-1+il)
      elseif(il.lt.0)then
      goto 80
      endif
      elseif(igut.eq.3)then
      if((lupt.eq.1).or.(lupt.eq.2))then
      i=lupi
      jj=vmap(invp1(i),1)
      elseif(lupt.eq.3)then
      i=lupi
      if(stib(ig+3).eq.0)then
      jj=ex(i)
      else
      jj=vmap(ex(i),ey(i))
      endif
      elseif(lupt.eq.4)then
      jj=nleg+lupi
      elseif(lupt.eq.5)then
      jj=nleg+lupi
      else
      goto 80
      endif
      j=stib(vfo(jj))
      ij=stib(udki(0)+igc)
      if((ij.le.0).or.(ij.gt.nvmk))then
      goto 80
      elseif((j.le.0).or.(j.gt.nvert))then
      goto 80
      endif
      il=stib(stib(vmkvlp(0)+ij)+j)
      if(il.gt.0)then
      ip=kli+1
      kli=kli+il
      if(kli.ge.swbuff)then
      goto 75
      endif
      ia=stib(stib(vmkvpp(0)+ij)+j)
      stwb(ip:kli)=stcb(ia:ia-1+il)
      elseif(il.lt.0)then
      goto 80
      endif
      else
      goto 80
      endif
      else
      goto 80
      endif
   70 continue
      ig=ig+stib(ig+2)
      goto 10
   20 continue
      if(lupt.ne.0)then
      goto 80
      endif
      if(kli.gt.0)then
      if(kli.lt.svbuff)then
      stwbs(1)=kli
      elseif(kli.lt.swbuff)then
      stwbs(1)=kli
      call qout
      stwbs(1)=0
      else
      goto 80
      endif
      endif
      goto 90
   75 continue
      if(stwbs(1).gt.0)then
      call qout
      stwbs(1)=0
      goto 01
      endif
   80 continue
      mlin(1:srec)='run-time error while processing <diagram>'
      call mput(1,0,0,0)
   90 return
      end
      subroutine style
      implicit integer(a-z)
      save
      parameter ( maxtak=2, maxpot=6 )
      parameter ( maxleg=10, maxrho=7, maxi=10, maxdeg=6 )
      parameter ( srec=81, ssrec=62 )
      parameter ( eoa=-2047, nap=-8191 )
      parameter ( nfiles=5 )
      parameter ( sibuff=1048574, scbuff=131064, swbuff=65528 )
      parameter ( sxbuff=2040 )
      common/z5in/lsfile,sfilea,sfileb
      common/z7in/aunit(1:nfiles)
      common/z12g/cflag(1:7),dflag(1:20),jflag(1:10),mflag(1:6)
      character*(srec) mlin
      common/z22g/mlin
      common/z23g/tak(1:maxtak),pot(0:maxpot),ks
      common/z24g/iogp(1:4)
      common/z26g/kc(0:0),kse(0:0),pske(0:0),wske(0:0),klo(0:0),nske
      character*(sxbuff) stxb
      common/z27g/stxb
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z35g/stibs(1:1),stcbs(1:1),stwbs(1:1)
      common/z36g/udkp(0:2),udkl(0:0),udkt(0:0),udki(0:0)
      common/z48g/acomma,ascol,albra,arbra,alpar,arpar
      common/z51g/aeq,agt,alt,aminus,anine,aplus,aslash,atimes,azero
      common/z52g/acf(0:127),abo(0:6),acr,alf,aspace,dquote,squote
      character*(srec) styb
      level=0
      bevel=0
      ks=0
      nudk=0
      call aoib(2)
      udkp(1)=stibs(1)-1
      stib(udkp(1))=eoa
      stib(stibs(1))=eoa
      call qopen(sfilea,lsfile,4,0)
      nlin=0
      clin=0
      jflag(9)=-1
      jflag(10)=1
      icc=0
   20 continue
      slin=0
      call qrlin(4,nlin,slin,qc)
      if(slin.eq.-1)then
      goto 70
      endif
      if(slin.gt.0)then
      styb(1:slin+1)=stcb(stcbs(1)+1:stcbs(1)+slin)//char(alf)
      endif
      if(abs(level-2).eq.2)then
      if(slin.eq.0)then
      goto 20
      endif
      if(qc.ne.0)then
      if(jflag(9).lt.0)then
      jflag(9)=1
      goto 20
      elseif(jflag(9).gt.0)then
      goto 20
      elseif(clin.eq.0)then
      goto 80
      endif
      endif
      if(slin.gt.2)then
      j1=ichar(styb(1:1))
      j2=ichar(styb(slin:slin))
      if((j1.eq.alt).and.(j2.eq.agt))then
      ii=styki(styb,2,slin-1,level)
      if(ii.eq.24)then
      if(jflag(9).lt.0)then
      jflag(9)=0
      endif
      if(abs(level-2).eq.2)then
      if(icc.ne.0)then
      goto 80
      else
      icc=1
      endif
      endif
      if(clin.eq.0)then
      clin=nlin
      goto 20
      else
      goto 80
      endif
      elseif(ii.eq.25)then
      if(clin.ne.0)then
      clin=0
      goto 20
      else
      goto 80
      endif
      elseif((level.eq.0).and.(ii.eq.1))then
      level=1
      bevel=1
      icc=0
      iogp(1)=stibs(1)+1
      goto 20
      endif
      endif
      endif
      if(clin.ne.0)then
      goto 20
      else
      goto 80
      endif
      elseif(clin.gt.0)then
      if(nlin-1.gt.clin)then
      goto 80
      endif
      endif
      j1=0
      j2=0
      k1=0
      k2=0
      kp1=0
      kp2=0
      lgm=0
      do i1=1,slin
      jj=ichar(styb(i1:i1))
      if((j2.eq.0).and.(jj.eq.alt))then
      j1=j1+1
      if(j1.eq.2)then
      call aoib(1)
      stib(stibs(1))=jj
      j1=0
      elseif(ichar(styb(i1+1:i1+1)).ne.alt)then
      lgm=1-lgm
      if(lgm.eq.0)then
      goto 80
      endif
      kp1=i1+1
      endif
      elseif((j1.eq.0).and.(jj.eq.albra))then
      j2=j2+1
      if(j2.eq.2)then
      call aoib(1)
      stib(stibs(1))=jj
      j2=0
      elseif(ichar(styb(i1+1:i1+1)).ne.albra)then
      if(j1.ne.0)then
      goto 80
      endif
      lgm=1-lgm
      if(lgm.eq.0)then
      goto 80
      endif
      kp1=i1+1
      endif
      elseif((j2.eq.0).and.(jj.eq.agt))then
      if(j1.eq.0)then
      k1=k1+1
      if(k1.eq.2)then
      call aoib(1)
      stib(stibs(1))=jj
      k1=0
      elseif(ichar(styb(i1+1:i1+1)).ne.agt)then
      goto 80
      endif
      else
      j1=0
      lgm=1-lgm
      if(lgm.ne.0)then
      goto 80
      endif
      kp2=i1-1
      if(kp2-kp1.lt.0)then
      goto 80
      endif
      ii=styki(styb,kp1,kp2,level)
      if(ii.le.0)then
      goto 80
      elseif(ii.lt.5)then
      if((kp1.ne.2).or.(kp2+1.ne.slin))then
      goto 80
      elseif((ks.ne.0).or.(clin.ne.0))then
      goto 80
      endif
      level=ii
      icc=0
      elseif(ii.eq.24)then
      if(clin.ne.0)then
      goto 80
      endif
      if(abs(level-2).eq.2)then
      if((kp1.ne.2).or.(kp2+1.ne.slin))then
      goto 80
      endif
      endif
      clin=nlin
      elseif(ii.eq.25)then
      if(clin.eq.0)then
      goto 80
      endif
      if(abs(level-2).eq.2)then
      if((kp1.ne.2).or.(kp2+1.ne.slin))then
      goto 80
      endif
      endif
      clin=0
      elseif(clin.eq.0)then
      call ktoc(level,ii,nlin)
      endif
      endif
      elseif((j1.eq.0).and.(jj.eq.arbra))then
      if(j2.eq.0)then
      k2=k2+1
      if(k2.eq.2)then
      call aoib(1)
      stib(stibs(1))=jj
      k2=0
      elseif(ichar(styb(i1+1:i1+1)).ne.arbra)then
      goto 80
      endif
      else
      j2=0
      if(j1.ne.0)then
      goto 80
      endif
      lgm=1-lgm
      if(lgm.ne.0)then
      goto 80
      endif
      kp2=i1-1
      if(kp2-kp1.lt.0)then
      goto 80
      endif
      if(clin.eq.0)then
      ii=udstyk(styb,kp1,kp2)
      endif
      if(ii.eq.0)then
      goto 80
      endif
      endif
      elseif((j1.eq.0).and.(j2.eq.0))then
      if(clin.eq.0)then
      call aoib(1)
      stib(stibs(1))=jj
      elseif(jj.ne.aspace)then
      if(abs(level-2).lt.2)then
      goto 80
      endif
      endif
      endif
      enddo
      if(lgm.ne.0)then
      goto 80
      endif
      if(clin.ne.0)then
      goto 20
      endif
      if(bevel.eq.level)then
      call aoib(1)
      stib(stibs(1))=alf
      else
      bevel=level
      if(level.gt.1)then
      call aoib(3)
      do i1=stibs(1)-2,stibs(1)
      stib(i1)=eoa
      enddo
      endif
      iogp(level)=stibs(1)+1
      endif
      goto 20
   70 continue
      if(level.lt.4)then
      mlin(1:srec)='is incomplete'
      call mput(1,0,0,-4)
      endif
      if(ks.eq.0)then
      call rsfki
      goto 90
      endif
   80 continue
      mlin(1:srec)='wrong syntax,'
      if(clin.eq.0)then
      ii=nlin
      else
      ii=clin
      endif
      call mput(1,ii,nlin,4)
   90 return
      end
      subroutine ktoc(isec,kco,nlin)
      implicit integer(a-z)
      save
      parameter ( maxtak=2, maxpot=6 )
      parameter ( srec=81, ssrec=62 )
      parameter ( sibuff=1048574, scbuff=131064, swbuff=65528 )
      common/z12g/cflag(1:7),dflag(1:20),jflag(1:10),mflag(1:6)
      common/z13g/kla(0:0),bbc(0:0)
      character*(srec) mlin
      common/z22g/mlin
      common/z23g/tak(1:maxtak),pot(0:maxpot),ks
      common/z26g/kc(0:0),kse(0:0),pske(0:0),wske(0:0),klo(0:0),nske
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z35g/stibs(1:1),stcbs(1:1),stwbs(1:1)
      common/z52g/acf(0:127),abo(0:6),acr,alf,aspace,dquote,squote
      integer astak(1:maxtak)
      ka=-1
      do i1=1,nske
      if(stib(kc(0)+i1).eq.kco)then
      ik=i1
      ka=stib(kla(0)+i1)
      j1=stib(pske(0)+i1)-1
      j2=stib(wske(0)+i1)
      goto 10
      endif
      enddo
   10 continue
      ii=0
      if(max(3,5,2).gt.maxpot)then
      ii=1
      elseif(abs(isec-2).gt.1)then
      ii=1
      elseif(abs(ka-3).gt.2)then
      ii=1
      endif
      if(ii.ne.0)then
      mlin(1:srec)='ktoc_1'
      call mput(1,0,0,0)
      endif
      if(mod(stib(kse(0)+ik),pot(isec)).lt.pot(isec-1))then
      goto 80
      endif
      kk=stib(klo(0)+i1)
      if(ks.eq.0)then
      if(kk.ne.0)then
      goto 80
      endif
      elseif(kk.ne.0)then
      if(mod(kk,pot(tak(ks))).lt.pot(tak(ks)-1))then
      goto 80
      endif
      elseif(ks.ne.0)then
      if(ka.eq.1)then
      goto 80
      endif
      endif
      if(ka.eq.1)then
      lupt=lupty(kco)
      if(lupt.gt.0)then
      if(ks.ge.maxtak)then
      mlin(1:srec)='too many nested loops,'
      call mput(1,nlin,nlin,4)
      endif
      if((lupt.eq.2).and.(abs(isec-2).eq.1))then
      if(tak(1).ne.1)then
      goto 80
      endif
      elseif(lupt.eq.5)then
      if(tak(1).ne.4)then
      goto 80
      endif
      else
      if(ks.ne.0)then
      goto 80
      endif
      endif
      ii=stibs(1)+1
      call aoib(6)
      stib(ii)=-ka
      stib(ii+1)=lupt
      stib(ii+2)=6
      stib(ii+3)=0
      stib(ii+4)=0
      stib(stibs(1))=0
      ks=ks+1
      astak(ks)=ii
      tak(ks)=lupt
      elseif(lupt.eq.-1)then
      if(ks.le.0)then
      goto 80
      endif
      call aoib(3)
      ii=stibs(1)-astak(ks)
      stib(stibs(1)-2)=-ka
      stib(stibs(1)-1)=-1
      stib(stibs(1))=2-ii
      stib(astak(ks)+3)=ii
      ks=ks-1
      endif
      elseif(ka.eq.5)then
      call aoib(1)
      stib(stibs(1))=alf
      else
      call aoib(3)
      stib(stibs(1)-2)=-ka
      stib(stibs(1)-1)=kco
      stib(stibs(1))=3
      endif
      if(stib(bbc(0)+ik).ne.0)then
      jflag(10)=0
      endif
      goto 90
   80 continue
      mlin(1:srec)='wrong syntax,'
      call mput(1,nlin,nlin,4)
   90 return
      end
      integer function styki(s,s1,s2,level)
      implicit integer(a-z)
      save
      parameter ( srec=81, ssrec=62 )
      parameter ( sibuff=1048574, scbuff=131064, swbuff=65528 )
      parameter ( maxtak=2, maxpot=6 )
      common/z12g/cflag(1:7),dflag(1:20),jflag(1:10),mflag(1:6)
      character*(srec) mlin
      common/z22g/mlin
      common/z23g/tak(1:maxtak),pot(0:maxpot),ks
      common/z26g/kc(0:0),kse(0:0),pske(0:0),wske(0:0),klo(0:0),nske
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      character*(*) s
      if(s1.gt.s2)then
      mlin(1:srec)='styki_1'
      call mput(1,0,0,0)
      endif
      styki=0
      call mstr0(s,s1,s2,pske(0),wske(0),ij)
      if(ij.ne.0)then
      styki=stib(kc(0)+ij)
      if(styki.le.4)then
      if(styki.ne.level+1)then
      styki=0
      endif
      elseif(abs(level-2)-1.le.0)then
      if(mod(stib(kse(0)+ij),pot(level)).lt.pot(level-1))then
      styki=0
      endif
      endif
      endif
      if(abs(styki-33).eq.1)then
      jflag(4)=1
      endif
      return
      end
      integer function udstyk(s,s1,s2)
      implicit integer(a-z)
      save
      parameter ( srec=81, ssrec=62 )
      parameter ( eoa=-2047, nap=-8191 )
      parameter ( sibuff=1048574, scbuff=131064, swbuff=65528 )
      parameter ( maxtak=2, maxpot=6 )
      common/z23g/tak(1:maxtak),pot(0:maxpot),ks
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z35g/stibs(1:1),stcbs(1:1),stwbs(1:1)
      common/z36g/udkp(0:2),udkl(0:0),udkt(0:0),udki(0:0)
      common/z41g/nudk,ngmk,npmk0,npmk,nvmk0,nvmk
      common/z45g/qdatp,qdatl,qvp,qvl,dprefp,dprefl
      common/z52g/acf(0:127),abo(0:6),acr,alf,aspace,dquote,squote
      character*(*) s
      udstyk=0
      jj=ks
      if(jj.eq.0)then
      kt=4
      else
      kt=max(tak(jj)+3,5)
      if(kt.gt.8)then
      goto 90
      endif
      endif
      kd=0
      j1=s1-1+dprefl
      if(s2.ge.j1)then
      if(s(s1:j1).eq.stcb(dprefp:dprefp-1+dprefl))then
      if(s2.gt.j1)then
      kd=1
      else
      goto 90
      endif
      endif
      endif
      j1=s1+kd*dprefl
      jj=s2-j1+1
      if(stdw(s,j1,s2).eq.0)then
      goto 90
      endif
      if(nudk.gt.0)then
      call mstr1(s,j1,s2,udkp(1),1,2,udstyk)
      endif
      j2=stibs(1)
      kstep=4
      call aoib(kstep)
      stib(j2+1)=-6
      if(udstyk.eq.0)then
      nudk=nudk+1
      stib(j2+2)=nudk
      else
      stib(j2+2)=udstyk
      endif
      stib(j2+3)=kstep
      stib(j2+4)=kd
      ii=udkp(1)
      if(udstyk.ne.0)then
      do i1=1,udstyk
      ii=stib(ii)
      enddo
      j2=ii-1+kt
      if(stib(j2).eq.0)then
      stib(j2)=1+kd
      elseif(stib(j2).eq.2-kd)then
      stib(j2)=3
      endif
      goto 90
      endif
      udstyk=stib(j2+2)
      do i1=2,nudk
      ii=stib(ii)
      enddo
      j3=stcbs(1)+1
      call aocb(jj+1)
      stcb(j3:stcbs(1))=s(j1:s2)//char(alf)
      j2=stibs(1)
      kstep=8
      call aoib(kstep)
      stib(j2-1)=stib(j2-1)+kstep
      stib(ii)=j2+1
      stib(j2+1)=eoa
      stib(j2+2)=j3
      stib(j2+3)=jj
      stib(j2+4)=0
      stib(j2+5)=0
      stib(j2+6)=0
      stib(j2+7)=0
      stib(j2+8)=0
      stib(j2+kt)=kd+1
   90 return
      end
      integer function lupty(ll)
      implicit integer(a-z)
      save
      parameter ( sibuff=1048574, scbuff=131064, swbuff=65528 )
      parameter ( srec=81, ssrec=62 )
      parameter ( lup1=11, lup2=19 )
      common/z13g/kla(0:0),bbc(0:0)
      character*(srec) mlin
      common/z22g/mlin
      common/z35g/stibs(1:1),stcbs(1:1),stwbs(1:1)
      if((ll.lt.lup1).or.(ll.gt.lup2))then
      mlin(1:srec)='lupty_1'
      call mput(1,0,0,0)
      endif
      if(ll.ne.lup2)then
      if(ll.lt.13)then
      lupty=ll-10
      else
      lupty=ll-12
      endif
      else
      lupty=-1
      endif
      return
      end
      subroutine stpa(s,ia,ib,iu)
      implicit integer(a-z)
      save
      parameter ( srec=81, ssrec=62 )
      parameter ( sibuff=1048574, scbuff=131064, swbuff=65528 )
      common/z21g/punct1(0:0)
      character*(srec) mlin
      common/z22g/mlin
      common/z30g/stib(1:sibuff)
      common/z35g/stibs(1:1),stcbs(1:1),stwbs(1:1)
      common/z52g/acf(0:127),abo(0:6),acr,alf,aspace,dquote,squote
      character*(*) s
      if((ia.le.0).or.(ib.lt.ia))then
      mlin(1:srec)='stpa_1'
      call mput(1,0,0,0)
      endif
      j1=0
      j2=0
      it=0
      bplic=0
      iu=0
      sts=1
      call vaoib(sts)
      st0=stibs(1)+sts
      s0=st0
      stib(s0)=0
      do i1=ia,ib+1
      if(i1.le.ib)then
      ic=ichar(s(i1:i1))
      else
      ic=alf
      endif
      if(ic.eq.squote)then
      bplic=1-bplic
      if(bplic.ne.0)then
      if(it.eq.0)then
      it=2
      j1=i1
      endif
      endif
      j2=i1
      elseif(bplic.ne.0)then
      if(ic.eq.alf)then
      goto 80
      endif
      j2=i1
      elseif(stib(punct1(0)+ic).ne.0)then
      if(it.eq.2)then
      if(stds(s,j1,j2,0).ge.0)then
      it=2
      elseif(stdw(s,j1,j2).ne.0)then
      it=3
      elseif(stdz(s,j1,j2).ne.0)then
      it=4
      elseif(stdq(s,j1,j2).ne.0)then
      it=5
      else
      it=0
      iu=1
      endif
      sts=sts+3
      call vaoib(sts)
      st0=st0+3
      stib(st0-2)=it
      stib(st0-1)=j1
      stib(st0)=j2
      stib(s0)=stib(s0)+1
      j1=0
      it=0
      endif
      sts=sts+3
      call vaoib(sts)
      st0=st0+3
      stib(st0-2)=1
      stib(st0-1)=ic
      stib(st0)=i1
      if(stib(s0).ge.0)then
      stib(s0)=stib(s0)+1
      endif
      elseif((ic.eq.aspace).or.(ic.eq.alf))then
      if(it.eq.2)then
      if(stds(s,j1,j2,0).ge.0)then
      it=2
      elseif(stdw(s,j1,j2).ne.0)then
      it=3
      elseif(stdz(s,j1,j2).ne.0)then
      it=4
      elseif(stdq(s,j1,j2).ne.0)then
      it=5
      else
      it=0
      iu=1
      endif
      sts=sts+3
      call vaoib(sts)
      st0=st0+3
      stib(st0-2)=it
      stib(st0-1)=j1
      stib(st0)=j2
      stib(s0)=stib(s0)+1
      j1=0
      it=0
      endif
      else
      if(it.eq.0)then
      it=2
      j1=i1
      elseif(it.ne.2)then
      goto 80
      endif
      j2=i1
      endif
      enddo
      if(bplic.eq.0)then
      goto 90
      endif
   80 continue
      iu=-1
   90 return
      end
      subroutine hsort(x,ia,ib)
      implicit integer(a-z)
      save
      parameter ( srec=81, ssrec=62 )
      character*(srec) mlin
      common/z22g/mlin
      integer x(*)
      if((ia.le.0).or.(ia.gt.ib))then
      mlin(1:srec)='hsort_1'
      call mput(1,0,0,0)
      elseif(ia.eq.ib)then
      goto 90
      endif
      i0=ia-1
      hr=ib-i0
      hl=1+hr/2
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
      goto 90
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
   90 return
      end
      subroutine xipht(ia,ib,ixipht)
      implicit integer(a-z)
      save
      parameter ( sibuff=1048574, scbuff=131064, swbuff=65528 )
      parameter ( srec=81, ssrec=62 )
      character*(srec) mlin
      common/z22g/mlin
      common/z30g/stib(1:sibuff)
      common/z35g/stibs(1:1),stcbs(1:1),stwbs(1:1)
      ii=0
      if(ia.le.0)then
      ii=1
      elseif(ib.lt.ia)then
      ii=1
      elseif(sibuff.lt.ib)then
      ii=1
      elseif(stibs(1).gt.sibuff)then
      ii=1
      elseif(ixipht.lt.0)then
      if(ia+ixipht.le.0)then
      ii=1
      endif
      endif
      if(ii.ne.0)then
      mlin(1:srec)='xipht_1'
      call mput(1,0,0,0)
      endif
      if(ixipht.lt.0)then
      j1=ia
      j2=ib
      j3=1
      elseif(ixipht.gt.0)then
      jj=ib-(stibs(1)-ixipht)
      if(jj.gt.0)then
      call vaoib(jj)
      endif
      j1=ib
      j2=ia
      j3=-1
      endif
      if(ixipht.ne.0)then
      do i1=j1,j2,j3
      stib(i1+ixipht)=stib(i1)
      enddo
      endif
      return
      end
      subroutine cxipht(ia,ib,ixipht)
      implicit integer(a-z)
      save
      parameter ( sibuff=1048574, scbuff=131064, swbuff=65528 )
      parameter ( srec=81, ssrec=62 )
      character*(srec) mlin
      common/z22g/mlin
      character*(scbuff) stcb
      common/z31g/stcb
      common/z35g/stibs(1:1),stcbs(1:1),stwbs(1:1)
      ii=0
      if(ia.le.0)then
      ii=1
      elseif(ib.lt.ia)then
      ii=1
      elseif(scbuff.lt.ib)then
      ii=1
      elseif(stcbs(1).gt.scbuff)then
      ii=1
      elseif(ixipht.lt.0)then
      if(ia+ixipht.le.0)then
      ii=1
      endif
      endif
      if(ii.ne.0)then
      mlin(1:srec)='cxipht_1'
      call mput(1,0,0,0)
      endif
      if(ixipht.lt.0)then
      j1=ia
      j2=ib
      j3=1
      elseif(ixipht.gt.0)then
      jj=ib-(stcbs(1)-ixipht)
      if(jj.gt.0)then
      call vaocb(jj)
      endif
      j1=ib
      j2=ia
      j3=-1
      endif
      if(ixipht.ne.0)then
      do i1=j1,j2,j3
      stcb(i1+ixipht:i1+ixipht)=stcb(i1:i1)
      enddo
      endif
      return
      end
      subroutine trm(ia,ib)
      implicit integer(a-z)
      save
      parameter ( srec=81, ssrec=62 )
      parameter ( sibuff=1048574, scbuff=131064, swbuff=65528 )
      character*(srec) mlin
      common/z22g/mlin
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z35g/stibs(1:1),stcbs(1:1),stwbs(1:1)
      if((ia.le.0).or.(ib.le.0))then
      mlin(1:srec)='trm_1'
      call mput(1,0,0,0)
      endif
      ab=ia*ib
      st0=stibs(1)-ab+1
      do i1=1,ab-2
      ix=0
      jj=i1
   10 continue
      ii=jj/ia
      kk=ii+ib*(jj-ii*ia)
      if(ix.eq.1)then
      ytmp=stib(st0+kk)
      stib(st0+kk)=xtmp
      xtmp=ytmp
      endif
      if(kk.gt.i1)then
      jj=kk
      goto 10
      elseif(kk.eq.i1)then
      if((ix.eq.0).and.(kk.ne.jj))then
      ix=1
      jj=i1
      xtmp=stib(st0+i1)
      goto 10
      endif
      endif
      enddo
      return
      end
      subroutine qrlin(junit,nlin,slin,qc)
      implicit integer(a-z)
      save
      parameter ( srec=81, ssrec=62 )
      parameter ( srecx=87 )
      parameter ( sibuff=1048574, scbuff=131064, swbuff=65528 )
      parameter ( sxbuff=2040 )
      parameter ( nfiles=5 )
      parameter ( lun1=7, lun2=63 )
      common/z7in/aunit(1:nfiles)
      character*(scbuff) stcb
      common/z31g/stcb
      character*(srec) mlin
      common/z22g/mlin
      character*(sxbuff) stxb
      common/z27g/stxb
      common/z35g/stibs(1:1),stcbs(1:1),stwbs(1:1)
      common/z52g/acf(0:127),abo(0:6),acr,alf,aspace,dquote,squote
      ii=0
      if(srecx.le.srec)then
      ii=1
      elseif((junit.le.0).or.(junit.ge.nfiles))then
      ii=1
      endif
      iunit=aunit(junit)
      if((iunit.lt.lun1).or.(iunit.gt.lun2))then
      ii=1
      endif
      if(ii.ne.0)then
      mlin(1:srec)='qrlin_1'
      call mput(1,0,0,0)
      endif
      do i1=1,srecx
      stxb(i1:i1)=char(aspace)
      enddo
      ios=1
      read(unit=iunit,fmt='(a)',iostat=ios)stxb(1:srecx)
      if(ios.ne.0)then
      if(ios.gt.0)then
      call uput(1)
      endif
      call qclose(junit,0)
      slin=-1
      goto 90
      endif
      nlin=nlin+1
      slin=0
      do i1=srecx,1,-1
      if(ichar(stxb(i1:i1)).ne.aspace)then
      slin=i1
      goto 10
      endif
      enddo
   10 continue
      j1=0
      if(nlin.eq.1)then
      if(slin.gt.2)then
      if(ichar(stxb(1:1)).eq.239)then
      if(ichar(stxb(2:2)).eq.187)then
      if(ichar(stxb(3:3)).eq.191)then
      j1=3
      slin=slin-j1
      endif
      endif
      endif
      endif
      endif
      if(slin.ge.srec)then
      mlin(1:srec)='line too long,'
      call mput(1,nlin,nlin,junit)
      endif
      jj=stcbs(1)
      ii=slin+1
      call vaocb(ii)
      i1=jj+ii
      stcb(i1:i1)=char(alf)
      qc=0
      if(slin.gt.0)then
      do i1=jj+1,jj+slin
      j1=j1+1
      j2=ichar(stxb(j1:j1))
      if((j2.lt.abo(1)).or.(j2.gt.abo(2)))then
      mlin(1:srec)='wrong syntax,'
      call mput(1,nlin,nlin,junit)
      endif
      stcb(i1:i1)=char(j2)
      enddo
      i1=jj+1
      if(acf(ichar(stcb(i1:i1))).eq.-1)then
      qc=1
      endif
      endif
   90 return
      end
      subroutine qout
      implicit integer(a-z)
      save
      parameter ( sibuff=1048574, scbuff=131064, swbuff=65528 )
      parameter ( lun1=7, lun2=63 )
      parameter ( nfiles=5 )
      parameter ( srec=81, ssrec=62 )
      common/z7in/aunit(1:nfiles)
      common/z12g/cflag(1:7),dflag(1:20),jflag(1:10),mflag(1:6)
      character*(srec) mlin
      common/z22g/mlin
      common/z35g/stibs(1:1),stcbs(1:1),stwbs(1:1)
      common/z52g/acf(0:127),abo(0:6),acr,alf,aspace,dquote,squote
      character*(swbuff) stwb
      common/z54g/stwb
      iunit=aunit(nfiles)
      ii=0
      if((iunit.lt.lun1).or.(iunit.gt.lun2))then
      ii=1
      elseif(stwbs(1).le.0)then
      ii=1
      elseif(ichar(stwb(stwbs(1):stwbs(1))).ne.alf)then
      ii=1
      endif
      if(ii.ne.0)then
      mlin(1:srec)='qout_1'
      call mput(1,0,0,0)
      endif
      j1=1
      if(cflag(3).eq.0)then
      do i1=1,stwbs(1)
      if(ichar(stwb(i1:i1)).eq.alf)then
      ios=1
      if(j1.lt.i1)then
      write(unit=iunit,fmt='(a)',iostat=ios)stwb(j1:i1-1)
      else
      write(unit=iunit,fmt='(a)',iostat=ios)
      endif
      if(ios.ne.0)then
      call uput(1)
      endif
      j1=i1+1
      endif
      enddo
      else
      write(unit=iunit,fmt='(a)',iostat=ios)stwb(1:stwbs(1))
      backspace(unit=iunit)
      if(ios.ne.0)then
      call uput(1)
      endif
      endif
      return
      end
      subroutine qopen(j1,j2,junit,jstat)
      implicit integer(a-z)
      save
      parameter ( srec=81, ssrec=62 )
      parameter ( sibuff=1048574, scbuff=131064, swbuff=65528 )
      parameter ( lun1=7, lun2=63 )
      parameter ( nfiles=5 )
      common/z7in/aunit(1:nfiles)
      character*(srec) mlin
      common/z22g/mlin
      character*(scbuff) stcb
      common/z31g/stcb
      common/z35g/stibs(1:1),stcbs(1:1),stwbs(1:1)
      character*(12) stact
      logical lun,loun
      ii=0
      if(j1.le.0)then
      ii=1
      elseif(j2.le.0)then
      ii=1
      elseif(junit.le.0)then
      ii=1
      elseif(junit.gt.nfiles)then
      ii=1
      elseif(aunit(junit).ne.0)then
      ii=1
      elseif(jstat.eq.0)then
      if(junit.eq.nfiles)then
      ii=1
      endif
      elseif(jstat.eq.1)then
      if(junit.ne.nfiles)then
      ii=1
      endif
      else
      ii=1
      endif
      if(ii.ne.0)then
      mlin(1:srec)='qopen_1'
      call mput(1,0,0,0)
      endif
      iunit=lun1-1+junit
   11 continue
      if(iunit.gt.lun2)then
      mlin(1:srec)='no logical unit number available'
      call mput(1,0,0,0)
      endif
      do i1=1,nfiles
      if(aunit(i1).eq.iunit)then
      iunit=iunit+1
      goto 11
      endif
      enddo
      ios=1
      loun=.true.
      inquire(unit=iunit,opened=loun,iostat=ios)
      if(ios.ne.0)then
      call uput(1)
      elseif(loun)then
      iunit=iunit+1
      goto 11
      endif
      if(jstat.eq.0)then
      jj=7
      stact(1:jj)='oldread'
      else
      jj=12
      stact(1:jj)='newreadwrite'
      endif
      lun=.true.
      inquire(file=stcb(j1:j1-1+j2),exist=lun,opened=loun,iostat=ios)
      if(ios.ne.0)then
      call uput(1)
      elseif(loun)then
      mlin(1:srec)='already opened'
      call mput(1,0,0,-junit)
      elseif(jstat.eq.0)then
      if(.not.lun)then
      mlin(1:srec)='could not be found'
      call mput(1,0,0,-junit)
      endif
      else
      if(lun)then
      mlin(1:srec)='already exists'
      call mput(1,0,0,-junit)
      endif
      endif
      open(unit=iunit,file=stcb(j1:j1-1+j2),access='sequential',
     :status=stact(1:3),action=stact(4:jj),iostat=ios)
      if(ios.ne.0)then
      mlin(1:srec)='could not be opened'
      call mput(1,0,0,-junit)
      endif
      aunit(junit)=iunit
      return
      end
      subroutine qclose(junit,istop)
      implicit integer(a-z)
      save
      parameter ( srec=81, ssrec=62 )
      parameter ( nfiles=5 )
      parameter ( lun1=7, lun2=63 )
      common/z7in/aunit(1:nfiles)
      common/z12g/cflag(1:7),dflag(1:20),jflag(1:10),mflag(1:6)
      character*(srec) mlin
      common/z22g/mlin
      logical lun
      character*(6) fstat
      ii=0
      if((istop.lt.0).or.(istop.gt.1))then
      ii=1
      elseif((junit.lt.0).or.(junit.gt.nfiles))then
      ii=1
      endif
      if(ii.ne.0)then
      mlin(1:srec)='qclose_1'
      call mput(-1,0,0,0)
      endif
      if(junit.eq.0)then
      j1=1
      j2=nfiles
      else
      j1=junit
      j2=junit
      endif
      do i1=j1,j2
      iunit=aunit(i1)
      if(iunit.gt.0)then
      if((iunit.lt.lun1).or.(iunit.gt.lun2))then
      mlin(1:srec)='qclose_1'
      call mput(-1,0,0,0)
      endif
      ios=1
      lun=.false.
      inquire(unit=iunit,opened=lun,iostat=ios)
      if((ios.eq.0).and.(lun))then
      jj=0
      if(istop.ne.0)then
      if(cflag(2).ne.0)then
      if(i1.eq.nfiles)then
      jj=6
      fstat(1:jj)='delete'
      endif
      endif
      endif
      if(jj.eq.0)then
      jj=4
      fstat(1:jj)='keep'
      endif
      close(unit=iunit,status=fstat(1:jj),iostat=ios)
      endif
      if(ios.ne.0)then
      call uput(1)
      endif
      aunit(i1)=0
      endif
      enddo
      return
      end
      subroutine inputs
      implicit integer(a-z)
      save
      parameter ( maxleg=10, maxrho=7, maxi=10, maxdeg=6 )
      parameter ( maxn=maxi+maxi-2, maxli=maxn+maxrho )
      parameter ( srec=81, ssrec=62 )
      parameter ( eoa=-2047, nap=-8191 )
      parameter ( nfiles=5 )
      parameter ( sibuff=1048574, scbuff=131064, swbuff=65528 )
      parameter ( sxbuff=2040 )
      common/z1in/leg(1:maxleg),nleg,nloop,incom
      common/z2in/momep(0:maxleg),momel(0:maxleg),kpqs(1:4)
      common/z3in/lmfile,mfilea,mfileb
      common/z4in/llfile,lfilea,lfileb
      common/z5in/lsfile,sfilea,sfileb
      common/z6in/lofile,ofilea,ofileb
      common/z7in/aunit(1:nfiles)
      common/z1g/g(1:maxn,1:maxn),rho(1:maxdeg),nrho,rhop1,rhop2
      common/z3g/nivd(0:0),dpntro(0:0),vparto(0:0),vval(0:0),nvert
      common/z11g/nphi,nblok,nprop,npprop
      common/z12g/cflag(1:7),dflag(1:20),jflag(1:10),mflag(1:6)
      common/z14g/zcho(0:maxli),zbri(0:maxli),zpro(0:maxli),
     :rbri(0:maxli),sbri(0:maxli)
      common/z15g/pten(1:10),iref,wiref,wsint
      common/z20g/tftyp(0:0),tfnarg(0:0),tfa(0:0),tfb(0:0),tfc(0:0),
     :tfo(0:0),tf2(0:0),ntf
      character*(srec) mlin
      common/z22g/mlin
      character*(sxbuff) stxb
      common/z27g/stxb
      common/z29g/pkey(0:0),wkey(0:0),ikey(0:0),pokey(0:0),wokey(0:0),
     :cokey(0:0),popt1(0:0),wopt1(0:0),fopt1(0:0),vopt1(0:0),
     :popt5(0:0),wopt5(0:0),copt5(0:0),popt0(0:0),wopt0(0:0),
     :copt0(0:0),vopt0(0:0)
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z33g/namep(0:0),namel(0:0)
      common/z34g/link(0:0),antiq(0:0),blok(0:0)
      common/z35g/stibs(1:1),stcbs(1:1),stwbs(1:1)
      common/z37g/drecp(0:0),drecl(0:0),drecii(0:0),irecc(0:0),
     :frecc(0:0),ndrec,ncom
      common/z39g/pmkp(0:0),pmkl(0:0),pmkd(0:0),pmkvpp(0:0),pmkvlp(0:0)
      common/z40g/vmkp(0:0),vmkl(0:0),vmkvpp(0:0),vmkvlp(0:0),vmks(0:0)
      common/z41g/nudk,ngmk,npmk0,npmk,nvmk0,nvmk
      common/z42g/pmkr(0:0),pmkvma(0:0),pmkvmi(0:0)
      common/z43g/vmkr(0:0),vmkmao(0:0),vmkmio(0:0)
      common/z44g/xtstrp(0:0),xtstrl(0:0)
      common/z45g/qdatp,qdatl,qvp,qvl,dprefp,dprefl
      common/z46g/popt3(0:0),wopt3(0:0),copt3(0:0),popt9(0:0),
     :wopt9(0:0),copt9(0:0)
      common/z47g/ndiagp,ndiagl,hhp,hhl,doffp,doffl,noffp,noffl
      common/z48g/acomma,ascol,albra,arbra,alpar,arpar
      common/z49g/popt7(0:0),wopt7(0:0),copt7(0:0)
      common/z50g/tfta(0:0),tftb(0:0),tftic(0:0),ntft
      common/z51g/aeq,agt,alt,aminus,anine,aplus,aslash,atimes,azero
      common/z52g/acf(0:127),abo(0:6),acr,alf,aspace,dquote,squote
      integer dreci(0:0),comii(0:0),comfi(0:0),ist1(0:9)
      integer xli(1:maxleg)
      datls=4
      incom=0
      outgo=0
      nleg=0
      level=0
      phase=1
      ncom=0
      call qopen(qdatp,qdatl,1,0)
      ndrec=0
      newc=1
      dbl=0
      nlin=0
      tsl=0
      do i1=0,9
      ist1(i1)=0
      enddo
   70 continue
      slin=0
      call qrlin(1,nlin,slin,qc)
      if(slin.eq.-1)then
      if(newc.eq.0)then
      goto 521
      endif
      if(level.eq.0)then
      mlin(1:srec)='is incomplete'
      call mput(1,0,0,-1)
      endif
      goto 404
      endif
      if(newc.ne.0)then
      dbl=nlin
      endif
      if((slin.eq.0).or.(qc.ne.0))then
      if(newc.eq.0)then
      goto 521
      else
      goto 70
      endif
      endif
      top=stcbs(1)+slin
      jflag(7)=0
      kk=0
      do i1=top,stcbs(1)+1,-1
      j1=ichar(stcb(i1:i1))
      if(j1.ne.aspace)then
      if(j1.eq.squote)then
      kk=1-kk
      endif
      bot=i1
      endif
      enddo
      if(kk.ne.0)then
      goto 521
      endif
      ii=slin+1
      tsl=tsl+ii
      if(tsl.gt.sxbuff)then
      mlin(1:srec)='statement too long,'
      goto 522
      endif
      if(ichar(stcb(top:top)).eq.ascol)then
      newc=1
      tsl=0
      else
      newc=0
      endif
      call aoib(datls)
      ndrec=ndrec+1
      stib(stibs(1)-3)=stcbs(1)+1
      stib(stibs(1)-2)=ii
      stib(stibs(1)-1)=nlin
      stib(stibs(1))=nlin-dbl
      call aocb(ii)
      if(newc.eq.0)then
      goto 70
      endif
      jj=stibs(1)-datls+1
      bot=stib(jj-datls*stib(stibs(1)))
      call stpa(stcb,bot,top,ii)
      jj=stib(stibs(1)+1)-3
      if(ii.ne.0)then
      goto 521
      elseif(jj.lt.0)then
      goto 521
      elseif(stib(stibs(1)+2).ne.3)then
      goto 521
      elseif(stib(stibs(1)+5).ne.1)then
      goto 521
      elseif(stib(stibs(1)+6).ne.aeq)then
      goto 521
      endif
      id1=stib(stibs(1)+3)
      id2=stib(stibs(1)+4)
      ncom=ncom+1
   16 continue
      if(level.eq.0)then
      level=level+1
      call mstr0(stcb,id1,id2,pokey(0),wokey(0),ij)
      if(ij.le.0)then
      goto 16
      endif
      ij=stib(cokey(0)+ij)
      if(ij.ne.1)then
      goto 521
      endif
      ist1(0)=1
      if(jj.gt.0)then
      if(mod(jj,2).eq.0)then
      goto 521
      endif
      j2=stibs(1)+5+3*jj
      do i1=stibs(1)+8,j2,6
      if(stib(i1).ne.3)then
      goto 521
      elseif(stib(i1+3).ne.1)then
      goto 521
      elseif(stib(i1+4).ne.acomma)then
      if(i1.ne.j2)then
      goto 521
      endif
      endif
      call mstr0(stcb,stib(i1+1),stib(i1+2),popt0(0),wopt0(0),ij)
      if(ij.le.0)then
      goto 521
      endif
      ik=stib(vopt0(0)+ij)
      ij=stib(copt0(0)+ij)
      if(cflag(ij).eq.0)then
      cflag(ij)=ik
      elseif(cflag(ij).eq.ik)then
      call wput(11,0,0)
      else
      call wput(-11,0,0)
      endif
      enddo
      endif
      elseif(level.eq.1)then
      call mstr0(stcb,id1,id2,pkey(0),wkey(0),ij)
      if(ij.le.0)then
      goto 521
      endif
      ij=stib(ikey(0)+ij)
      if((ij+ist1(0).ne.ncom).or.(ist1(ij).ne.0))then
      goto 521
      endif
      ist1(ij)=ncom
      if(ncom-ist1(0).eq.8)then
      level=level+1
      endif
      elseif(level.eq.2)then
      level=level+1
      call mstr0(stcb,id1,id2,pokey(0),wokey(0),ij)
      if(ij.le.0)then
      goto 16
      endif
      ij=stib(cokey(0)+ij)
      if(ij.eq.2)then
      ist1(9)=ncom
      if(cflag(6).eq.0)then
      cflag(6)=1
      else
      goto 521
      endif
      else
      goto 521
      endif
      else
      call mstr0(stcb,id1,id2,popt3(0),wopt3(0),ij)
      if(ij.eq.0)then
      goto 521
      endif
      if(stib(stibs(1)+8).ne.3)then
      goto 521
      endif
      endif
      if(cflag(1).ge.0)then
      if(ncom.eq.1)then
      call spp(0)
      endif
      j1=bot
      do i1=bot,top+1
      if(ichar(stcb(i1:i1)).eq.alf)then
      if(j1.lt.i1)then
      write(unit=*,fmt='(1x,a)')stcb(j1:i1-1)
      else
      write(unit=*,fmt='(1x,a)')
      endif
      j1=i1+1
      endif
      enddo
      endif
      goto 70
  404 continue
      if(cflag(1).ge.0)then
      call spp(1)
      endif
      ii=stibs(1)
      call aoib(datls)
      do i1=ii+1,stibs(1)
      stib(i1)=eoa
      enddo
      call trm(datls,ndrec+1)
      ii=ndrec+1
      drecii(0)=stibs(1)-ii
      dreci(0)=drecii(0)-ii
      drecl(0)=dreci(0)-ii
      drecp(0)=drecl(0)-ii
      ncom=0
      j1=0
   56 continue
      j1=j1+1
      ii=stib(drecp(0)+j1)
      jj=ii-1
      j2=j1
   58 continue
      jj=jj+stib(drecl(0)+j2)
      if(j2.lt.ndrec)then
      if(stib(drecii(0)+j2+1).ne.0)then
      j2=j2+1
      goto 58
      endif
      endif
      j1=j2
      ncom=ncom+1
      call aoib(2)
      stib(stibs(1)-1)=ii
      stib(stibs(1))=jj
      if(j1.lt.ndrec)then
      goto 56
      endif
      call aoib(2)
      stib(stibs(1)-1)=eoa
      stib(stibs(1))=eoa
      call trm(2,ncom+1)
      ii=ncom+1
      comfi(0)=stibs(1)-ii
      comii(0)=comfi(0)-ii
      irecc(0)=stibs(1)
      call aoib(2*ii)
      frecc(0)=irecc(0)+ii
      stib(frecc(0))=eoa
      stib(stibs(1))=eoa
      j1=0
      j2=0
  137 continue
      if(j1.lt.ndrec)then
      j1=j1+1
      else
      goto 143
      endif
      if(stib(drecii(0)+j1).eq.0)then
      if(j2.gt.0)then
      stib(frecc(0)+j2)=j1-1
      endif
      j2=j2+1
      stib(irecc(0)+j2)=j1
      endif
      if(j2.le.ncom)then
      goto 137
      endif
  143 continue
      stib(frecc(0)+ncom)=ndrec
      if((j1.ne.ndrec).or.(j2.ne.ncom))then
      mlin(1:srec)='inputs_1'
      call mput(1,0,0,0)
      endif
      phase=2
      icom=ist1(1)
      call stpa(stcb,stib(comii(0)+icom),stib(comfi(0)+icom),ii)
      jj=stib(stibs(1)+1)-3
      lofile=0
      j1=stibs(1)+5
      do i1=1,jj
      j1=j1+3
      if(stib(j1).ne.2)then
      goto 521
      endif
      if(i1.eq.1)then
      ofilea=stib(j1+1)
      endif
      if(i1.eq.jj)then
      ofileb=stib(j1+2)
      lofile=ofileb-ofilea+1
      endif
      enddo
      if(lofile.gt.0)then
      ii=stds(stcb,ofilea,ofileb,2)
      if(ii.lt.0)then
      goto 521
      elseif(ii.eq.0)then
      lofile=0
      elseif(ii.lt.lofile-2)then
      lofile=ii
      ofilea=stcbs(1)+1
      call aocb(ii)
      ofileb=stcbs(1)
      else
      lofile=lofile-2
      ofilea=ofilea+1
      ofileb=ofileb-1
      endif
      endif
      if(cflag(2).eq.0)then
      if(lofile.gt.0)then
      cflag(2)=1
      endif
      else
      cflag(2)=0
      endif
      if(cflag(1).eq.0)then
      if(cflag(2).eq.0)then
      cflag(1)=1
      endif
      elseif(cflag(1).gt.0)then
      cflag(1)=cflag(1)-1
      endif
      icom=ist1(2)
      call stpa(stcb,stib(comii(0)+icom),stib(comfi(0)+icom),ii)
      jj=stib(stibs(1)+1)-3
      lsfile=0
      j1=stibs(1)+5
      do i1=1,jj
      j1=j1+3
      if(stib(j1).ne.2)then
      goto 521
      endif
      if(i1.eq.1)then
      sfilea=stib(j1+1)
      endif
      if(i1.eq.jj)then
      sfileb=stib(j1+2)
      lsfile=sfileb-sfilea+1
      endif
      enddo
      if(lsfile.gt.0)then
      ii=stds(stcb,sfilea,sfileb,2)
      if(ii.lt.0)then
      goto 521
      elseif(ii.eq.0)then
      lsfile=0
      elseif(ii.lt.lsfile-2)then
      lsfile=ii
      sfilea=stcbs(1)+1
      call aocb(ii)
      sfileb=stcbs(1)
      else
      lsfile=lsfile-2
      sfilea=sfilea+1
      sfileb=sfileb-1
      endif
      endif
      if(cflag(2).ne.0)then
      if(lsfile.le.0)then
      mlin(1:srec)='name of style-file is missing,'
      call mput(1,0,0,1)
      endif
      call style
      endif
      llfile=0
      icom=ist1(3)
      call stpa(stcb,stib(comii(0)+icom),stib(comfi(0)+icom),ii)
      jj=stib(stibs(1)+1)-3
      lmfile=0
      j1=stibs(1)+5
      do i1=1,jj
      j1=j1+3
      if(stib(j1).ne.2)then
      goto 521
      endif
      if(i1.eq.1)then
      mfilea=stib(j1+1)
      endif
      if(i1.eq.jj)then
      mfileb=stib(j1+2)
      lmfile=mfileb-mfilea+1
      endif
      enddo
      if(lmfile.gt.0)then
      ii=stds(stcb,mfilea,mfileb,2)
      if(ii.lt.0)then
      goto 521
      elseif(ii.eq.0)then
      lmfile=0
      elseif(ii.lt.lmfile-2)then
      lmfile=ii
      mfilea=stcbs(1)+1
      call aocb(ii)
      mfileb=stcbs(1)
      else
      lmfile=lmfile-2
      mfilea=mfilea+1
      mfileb=mfileb-1
      endif
      endif
      if(lmfile.le.0)then
      goto 521
      endif
      call model
      icom=ist1(4)
   26 continue
      call stpa(stcb,stib(comii(0)+icom),stib(comfi(0)+icom),ii)
      if(stib(stibs(1)+1).eq.3)then
      nphia=0
      else
      if((stib(stibs(1)+11).eq.1).and.
     :(stib(stibs(1)+12).eq.albra))then
      j=stib(stibs(1)+1)-2
      nphia=j/5
      if((nphia.lt.0).or.(j.ne.5*nphia))then
      goto 521
      endif
      if(icom.eq.ist1(4))then
      cflag(7)=1
      elseif(incom.eq.0)then
      cflag(7)=1
      elseif(cflag(7).eq.0)then
      goto 521
      endif
      j=stibs(1)+8
      do i1=1,nphia
      if(stib(j).ne.3)then
      goto 521
      elseif(stib(j+3).ne.1)then
      goto 521
      elseif(stib(j+4).ne.albra)then
      goto 521
      elseif(stib(j+6).ne.3)then
      goto 521
      elseif(stib(j+9).ne.1)then
      goto 521
      elseif(stib(j+10).ne.arbra)then
      goto 521
      elseif(stib(j+12).ne.1)then
      goto 521
      elseif(stib(j+13).ne.acomma)then
      if(i1.ne.nphia)then
      goto 521
      endif
      endif
      nleg=nleg+1
      if(nleg.gt.maxleg)then
      mlin(1:srec)='too many legs'
      call mput(1,0,0,0)
      endif
      momep(nleg)=stib(j+7)
      momel(nleg)=stib(j+8)-momep(nleg)+1
      j=j+15
      enddo
      else
      jj=stib(stibs(1)+1)-2
      nphia=jj/2
      if(jj.ne.2*nphia)then
      goto 521
      endif
      if(icom.eq.ist1(4))then
      cflag(7)=0
      elseif(cflag(7).ne.0)then
      goto 521
      endif
      j1=stibs(1)+2
      do i1=1,nphia
      j1=j1+6
      if(stib(j1).ne.3)then
      goto 521
      elseif(stib(j1+3).ne.1)then
      goto 521
      elseif(stib(j1+4).ne.acomma)then
      if(i1.ne.nphia)then
      goto 521
      endif
      endif
      nleg=nleg+1
      if(nleg.gt.maxleg)then
      mlin(1:srec)='too many legs'
      call mput(1,0,0,0)
      endif
      call aocb(1)
      j2=stcbs(1)
      if(icom.eq.ist1(4))then
      stcb(j2:j2)=char(kpqs(2))
      else
      stcb(j2:j2)=char(kpqs(3))
      endif
      momep(nleg)=j2
      call dkar(i1,ii)
      momel(nleg)=ii+1
      j2=j2+ii
      call aocb(ii)
      enddo
      endif
      endif
      if(icom.eq.ist1(4))then
      ij=0
      else
      ij=incom
      endif
      ii=6+9*cflag(7)
      jj=stibs(1)+8-ii
      do i1=1,nphia
      jj=jj+ii
      call mstr0(stcb,stib(jj+1),stib(jj+2),namep(0),namel(0),jk)
      if(jk.eq.0)then
      mlin(1:srec)='unknown external particle(s)'
      call mput(1,0,0,0)
      endif
      ij=ij+1
      if(icom.eq.ist1(4))then
      leg(ij)=stib(link(0)+jk)
      else
      leg(ij)=jk
      endif
      enddo
      if(icom.eq.ist1(4))then
      icom=ist1(5)
      incom=nphia
      goto 26
      else
      outgo=nphia
      endif
      if(nleg.ne.incom+outgo)then
      mlin(1:srec)='inputs_2'
      call mput(1,0,0,0)
      endif
      if(nleg.eq.0)then
      mlin(1:srec)='no external particles listed'
      call mput(1,0,0,0)
      endif
      icom=ist1(6)
      call stpa(stcb,stib(comii(0)+icom),stib(comfi(0)+icom),ii)
      if(stib(stibs(1)+1).ne.4)then
      goto 521
      elseif(stib(stibs(1)+8).ne.4)then
      goto 521
      endif
      nloop=stoz(stcb,stib(stibs(1)+9),stib(stibs(1)+10))
      if(nloop.lt.0)then
      mlin(1:srec)='check loops in'
      call mput(1,0,0,1)
      elseif((nloop.gt.maxrho).or.
     :(nloop+nleg.gt.max(maxleg,maxrho)))then
      mlin(1:srec)='too many legs and/or loops'
      call mput(1,0,0,0)
      endif
      if(nleg+nloop.lt.2)then
      mlin(1:srec)='check legs, loops in'
      call mput(1,0,0,1)
      endif
      if((nleg.eq.2).and.(nloop.eq.0))then
      mlin(1:srec)='case legs=2, loops=0 not accepted'
      call mput(1,0,0,0)
      endif
      icom=ist1(7)
      call stpa(stcb,stib(comii(0)+icom),stib(comfi(0)+icom),ii)
      if(stib(stibs(1)+1).eq.3)then
      call aocb(1)
      stcb(stcbs(1):stcbs(1))=char(kpqs(1))
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
      j1=stcbs(1)+1
      ii=momel(0)
      call aocb(ii)
      k=stcbs(1)
      j2=momep(0)-1+momel(0)
      stcb(j1:k)=stcb(momep(0):j2)
      do i1=1,nloop
      call dkar(i1,jk)
      ik=k+jk
      do i2=1,nleg
      if(momel(i2).eq.momel(0)+jk)then
      j2=momep(i2)-1+momel(i2)
      if(stcb(j1:ik).eq.stcb(momep(i2):j2))then
      mlin(1:srec)='conflict between names of external and'
     ://' internal momenta'
      if(cflag(2).gt.0)then
      call mput(1,0,0,0)
      elseif(cflag(1).gt.0)then
      call mput(0,0,0,0)
      else
      jflag(8)=min(jflag(8)+1,iref)
      endif
      goto 55
      endif
      endif
      enddo
      enddo
   55 continue
      call aocb(-ii)
      icom=ist1(8)
      call stpa(stcb,stib(comii(0)+icom),stib(comfi(0)+icom),ii)
      jj=stib(stibs(1)+1)-2
      if(jj.eq.1)then
      ii=0
      else
      ii=jj/2
      if(jj.ne.2*ii)then
      goto 521
      endif
      endif
      do i1=1,ii
      jj=stibs(1)+6*i1
      if((stib(jj+5).ne.1).or.(stib(jj+6).ne.acomma))then
      if(i1.ne.ii)then
      goto 521
      endif
      endif
      j1=stib(jj+3)
      j2=stib(jj+4)
      call mstr0(stcb,j1,j2,popt1(0),wopt1(0),ij)
      if(ij.ne.0)then
      j1=stib(fopt1(0)+ij)
      j2=stib(vopt1(0)+ij)
      if(dflag(j1).ne.0)then
      call wput(12,0,0)
      if(dflag(j1)*j2.lt.0)then
      jflag(1)=1
      elseif(dflag(j1).lt.j2)then
      dflag(j1)=j2
      endif
      endif
      dflag(j1)=j2
      else
      mlin(1:srec)='unknown option,'
      goto 522
      endif
      enddo
      if(cflag(6).ne.0)then
      icom=ist1(9)
      call stpa(stcb,stib(comii(0)+icom),stib(comfi(0)+icom),ii)
      jj=stibs(1)+1
      if(stib(jj).ne.4)then
      goto 521
      endif
      jj=stibs(1)+8
      if(stib(jj).ne.4)then
      goto 521
      endif
      jj=jj+1
      j1=stib(jj)
      j2=stib(jj+1)
      if(j2.lt.j1)then
      goto 521
      elseif(stdz(stcb,j1,j2).eq.0)then
      goto 521
      endif
      kk=sigz(stcb,j1,j2)
      if(kk.lt.0)then
      goto 521
      endif
      if(kk.gt.0)then
      if(ichar(stcb(j1:j1)).eq.aplus)then
      j1=j1+1
      endif
      j3=j1-1
      ii=0
      do i1=j1,j2-1
      if(ii.eq.0)then
      if(ichar(stcb(i1:i1)).eq.azero)then
      j3=i1
      else
      ii=1
      endif
      endif
      enddo
      noffl=j2-j3
      if(noffl.ge.wsint)then
      ii=stoz(stcb,j1,j2)
      endif
      stcb(noffp+1:noffp+noffl)=stcb(j3+1:j2)
      doffl=noffl
      stcb(doffp+1:doffp+doffl)=stcb(j3+1:j2)
      else
      cflag(6)=0
      endif
      endif
      ipass=1
      nphib=0
      ntf1=0
      ntf=0
      do i1=0,maxli
      zpro(i1)=1
      zbri(i1)=1
      rbri(i1)=1
      sbri(i1)=1
      if(i1.lt.nloop)then
      zcho(i1)=0
      else
      zcho(i1)=1
      endif
      enddo
      jcom=1
      do i1=0,9
      if(ist1(i1).gt.0)then
      jcom=jcom+1
      endif
      enddo
   54 continue
      do 71 icom=jcom,ncom
      call stpa(stcb,stib(comii(0)+icom),stib(comfi(0)+icom),ii)
      jj=stib(stibs(1)+1)
      ik=stibs(1)+3*jj
      if(ii.ne.0)then
      goto 521
      elseif(jj.lt.7)then
      goto 521
      elseif(mod(jj,2).eq.0)then
      goto 521
      elseif(stib(stibs(1)+8).ne.3)then
      goto 521
      elseif(stib(stibs(1)+11).ne.1)then
      goto 521
      elseif(stib(stibs(1)+12).ne.albra)then
      goto 521
      elseif(stib(ik-7).ne.4)then
      goto 521
      elseif(stib(ik-4).ne.1)then
      goto 521
      elseif(stib(ik-3).ne.arbra)then
      goto 521
      endif
      j1=stib(stibs(1)+3)
      j2=stib(stibs(1)+4)
      call mstr0(stcb,j1,j2,popt3(0),wopt3(0),ij)
      if(ij.eq.0)then
      goto 521
      endif
      tfv=stib(copt3(0)+ij)
      j1=stib(stibs(1)+9)
      j2=stib(stibs(1)+10)
      call mstr0(stcb,j1,j2,popt5(0),wopt5(0),ij)
      if(ij.eq.0)then
      goto 521
      endif
      atyp=stib(copt5(0)+ij)
      if(atyp.lt.0)then
      goto 521
      endif
      styp=atyp*tfv
      if(atyp.lt.12)then
      if(jj.lt.9)then
      goto 521
      elseif(stib(ik-13).ne.4)then
      goto 521
      elseif(stib(ik-10).ne.1)then
      goto 521
      elseif(stib(ik-9).ne.acomma)then
      goto 521
      endif
      endif
      if(atyp.eq.6)then
      dflag(15)=1
      elseif(atyp.eq.7)then
      dflag(16)=1
      elseif(atyp.eq.11)then
      dflag(17)=1
      elseif(atyp.eq.12)then
      dflag(18)=1
      endif
      j1=stib(stibs(1)+1)
      if(atyp.lt.11)then
      j3=max(0,(j1-9)/2)
      elseif(atyp.eq.11)then
      j3=max(0,(j1-11)/2)
      elseif(atyp.eq.12)then
      j3=max(0,(j1-5)/2)
      endif
      if(atyp.ge.11)then
      if(j3.le.0)then
      goto 521
      endif
      endif
      if(ipass.eq.1)then
      if(j3.ne.0)then
      ntf1=ntf1+1
      nphib=nphib+j3
      endif
      goto 71
      endif
      if(atyp.lt.12)then
      j1=stoz(stcb,stib(ik-12),stib(ik-11))
      j2=stoz(stcb,stib(ik-6),stib(ik-5))
      if(j1.gt.j2)then
      goto 521
      endif
      if(atyp.lt.6)then
      if(j1.lt.0)then
      goto 521
      endif
      if(atyp.gt.1)then
      jflag(4)=1
      endif
      elseif(atyp.eq.11)then
      if(j1.le.0)then
      goto 521
      endif
      endif
      else
      j1=0
      j2=0
      endif
      if(j3.eq.0)then
      goto 73
      endif
      ntf=ntf+1
      stib(tftyp(0)+ntf)=styp
      stib(tfa(0)+ntf)=j1
      stib(tfb(0)+ntf)=j2
      stib(tfc(0)+ntf)=0
      stib(tfnarg(0)+ntf)=j3
      stib(stib(tfo(0)+ntf))=eoa
      if(ntf.lt.ntf1)then
      stib(tfo(0)+ntf+1)=stib(tfo(0)+ntf)+stib(tfnarg(0)+ntf)+1
      endif
      if(atyp.gt.5)then
      goto 37
      endif
      do i1=1,stib(tfnarg(0)+ntf)
      j=stibs(1)+6*i1
      if(stib(j+8).ne.3)then
      goto 521
      elseif(stib(j+11).ne.1)then
      goto 521
      elseif(stib(j+12).ne.acomma)then
      goto 521
      endif
      j1=stib(j+9)
      j2=stib(j+10)
      call mstr0(stcb,j1,j2,namep(0),namel(0),jk)
      if(jk.eq.0)then
      mlin(1:srec)='unknown field,'
      goto 522
      endif
      stib(stib(tfo(0)+ntf)+i1)=jk
      enddo
      i4=stib(tfo(0)+ntf)
      call hsort(stib,i4+1,i4+stib(tfnarg(0)+ntf))
      tfcl=stib(tfnarg(0)+ntf)
      tfch=1
      tfcn=1
      i3=1
      do i1=2,stib(tfnarg(0)+ntf)
      if((stib(i4+i1).eq.stib(i4+i1-1)).or.
     :(stib(i4+i1).eq.stib(link(0)+stib(i4+i1-1))))then
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
      enddo
      if(i3.lt.tfcl)then
      tfcl=i3
      endif
      if(tfcn.ne.nprop)then
      tfcl=0
      endif
      if(styp.gt.0)then
      do i2=0,maxli
      if((i2*tfch.lt.stib(tfa(0)+ntf)).or.
     :(i2*tfcl.gt.stib(tfb(0)+ntf)))then
      if(styp.eq.1)then
      zpro(i2)=0
      elseif(styp.eq.2)then
      zbri(i2)=0
      elseif(styp.eq.3)then
      zcho(i2)=0
      elseif(styp.eq.4)then
      rbri(i2)=0
      elseif(styp.eq.5)then
      sbri(i2)=0
      endif
      endif
      enddo
      elseif(styp.lt.0)then
      do i2=0,maxli
      if(i2*tfch.le.stib(tfb(0)+ntf))then
      if(i2*tfcl.ge.stib(tfa(0)+ntf))then
      if(styp.eq.-1)then
      zpro(i2)=0
      elseif(styp.eq.-2)then
      zbri(i2)=0
      elseif(styp.eq.-3)then
      zcho(i2)=0
      elseif(styp.eq.-4)then
      rbri(i2)=0
      elseif(styp.eq.-5)then
      sbri(i2)=0
      endif
      endif
      endif
      enddo
      endif
      goto 71
   73 continue
      if(styp.gt.0)then
      do i2=0,maxli
      if((i2.lt.j1).or.(i2.gt.j2))then
      if(styp.eq.1)then
      zpro(i2)=0
      elseif(styp.eq.2)then
      zbri(i2)=0
      elseif(styp.eq.3)then
      zcho(i2)=0
      elseif(styp.eq.4)then
      rbri(i2)=0
      elseif(styp.eq.5)then
      sbri(i2)=0
      endif
      endif
      enddo
      else
      do i2=0,maxli
      if(i2.ge.j1)then
      if(i2.le.j2)then
      if(styp.eq.-1)then
      zpro(i2)=0
      elseif(styp.eq.-2)then
      zbri(i2)=0
      elseif(styp.eq.-3)then
      zcho(i2)=0
      elseif(styp.eq.-4)then
      rbri(i2)=0
      elseif(styp.eq.-5)then
      sbri(i2)=0
      endif
      endif
      endif
      enddo
      endif
      goto 71
   37 continue
      if(atyp.eq.6)then
      if(stib(tfnarg(0)+ntf).ne.1)then
      goto 521
      endif
      j1=stib(stibs(1)+15)
      j2=stib(stibs(1)+16)
      call mstr0(stcb,j1,j2,vmkp(0),vmkl(0),kid)
      if(kid.eq.0)then
      mlin(1:srec)='wrong argument in vsum,'
      goto 522
      endif
      if(stib(vmkr(0)+kid).ne.4)then
      mlin(1:srec)='invalid keyword in vsum,'
      goto 522
      endif
      stib(vmks(0)+kid)=1
      stib(stib(tfo(0)+ntf)+1)=kid
      elseif(atyp.eq.7)then
      if(stib(tfnarg(0)+ntf).ne.1)then
      goto 521
      endif
      j1=stib(stibs(1)+15)
      j2=stib(stibs(1)+16)
      call mstr0(stcb,j1,j2,pmkp(0),pmkl(0),kid)
      if(kid.eq.0)then
      mlin(1:srec)='wrong argument in psum,'
      goto 522
      endif
      if((stib(pmkr(0)+kid).ne.4).or.(stib(pmkd(0)+kid).ne.1))then
      mlin(1:srec)='invalid keyword in psum,'
      goto 522
      endif
      stib(stib(tfo(0)+ntf)+1)=kid
      j1=stib(pmkvmi(0)+kid)
      j2=stib(pmkvma(0)+kid)
      do i2=0,maxli
      if(styp.gt.0)then
      if((i2*j1.gt.stib(tfb(0)+ntf)).or.
     :(i2*j2.lt.stib(tfa(0)+ntf)))then
      zpro(i2)=0
      endif
      else
      if((i2*j1.ge.stib(tfa(0)+ntf)).and.
     :(i2*j2.le.stib(tfb(0)+ntf)))then
      zpro(i2)=0
      endif
      endif
      enddo
      elseif((atyp.eq.11).or.(atyp.eq.12))then
      if(nleg.lt.2)then
      call wput(-15,0,0)
      endif
      j3=stib(tfnarg(0)+ntf)
      if(j3.le.0)then
      goto 521
      elseif(j3.gt.nleg)then
      goto 521
      endif
      if(atyp.eq.11)then
      if(j3.lt.stib(tfb(0)+ntf))then
      goto 521
      endif
      else
      if(j3.eq.nleg)then
      goto 521
      endif
      endif
      do i1=1,nleg
      xli(i1)=0
      enddo
      do i1=1,j3
      j1=stibs(1)+9+6*i1
      if(stib(j1-1).ne.4)then
      goto 521
      endif
      j2=stoz(stcb,stib(j1),stib(j1+1))
      if(j2.ge.0)then
      goto 521
      elseif(j2+2*nleg.lt.0)then
      goto 521
      endif
      j1=(-j2)/2
      if(j2+2*j1.eq.0)then
      j1=j1+incom
      if(j1.gt.nleg)then
      goto 521
      endif
      else
      j1=j1+1
      if(j1.gt.incom)then
      goto 521
      endif
      endif
      if(xli(j1).ne.0)then
      goto 521
      endif
      xli(j1)=1
      stib(stib(tfo(0)+ntf)+i1)=j1
      enddo
      if(atyp.eq.11)then
      j1=stibs(1)+15+6*j3
      call mstr0(stcb,stib(j1),stib(j1+1),popt7(0),wopt7(0),ij)
      if(ij.eq.0)then
      mlin(1:srec)='invalid elink statement,'
      goto 522
      endif
      j1=stib(copt7(0)+ij)
      stib(tfc(0)+ntf)=j1
      j2=0
      if(j1.ne.0)then
      if(j3.eq.nleg)then
      if(stib(tfa(0)+ntf).eq.1)then
      if(stib(tfb(0)+ntf).eq.nleg)then
      j2=styp
      endif
      endif
      endif
      else
      if(stib(tfa(0)+ntf).eq.1)then
      if(stib(tfb(0)+ntf).eq.j3)then
      j2=styp
      endif
      endif
      endif
      if(j2.ne.0)then
      if(j2.gt.0)then
      call wput(14,0,0)
      else
      jflag(1)=1
      call wput(13,0,0)
      endif
      endif
      endif
      else
      mlin(1:srec)='inputs_3'
      call mput(1,0,0,0)
      endif
   71 continue
      if(ipass.eq.1)then
      ii=ntf1+1
      jj=7*ii+nphib+ntf1
      tftyp(0)=stibs(1)
      call aoib(jj)
      tfnarg(0)=tftyp(0)+ii
      tfa(0)=tfnarg(0)+ii
      tfb(0)=tfa(0)+ii
      tfc(0)=tfb(0)+ii
      tf2(0)=tfc(0)+ii
      tfo(0)=tf2(0)+ii
      stib(tfo(0)+1)=tfo(0)+ii
      stib(tfnarg(0))=eoa
      stib(tfa(0))=eoa
      stib(tfb(0))=eoa
      stib(tfc(0))=eoa
      stib(tf2(0))=eoa
      stib(tfo(0))=eoa
      stib(stibs(1))=eoa
      ipass=2
      goto 54
      endif
      if(cflag(1).ge.0)then
      call spp(2)
      endif
      call vsig
      do i1=1,ntft
      stib(tfta(0)+i1)=0
      stib(tftb(0)+i1)=0
      enddo
      j1=tftic(0)
      j3=1
      jj=0
   29 continue
      j1=j1+1
      j2=stib(j1)
      if(j2.ge.0)then
      if(j2.gt.0)then
      j3=j2
      jj=jj+1
      endif
      goto 29
      endif
      j1=j1-tftic(0)-1
      ii=0
      do i1=1,j1
      j2=stib(tftic(0)+i1)
      if(j2.gt.0)then
      do i2=1,ntf
      if(abs(stib(tftyp(0)+i2)).eq.i1)then
      ii=ii+1
      j3=stib(tfta(0)+j2)
      if(j3.eq.0)then
      stib(tfta(0)+j2)=ii
      endif
      stib(tftb(0)+j2)=ii
      stib(tf2(0)+ii)=i2
      endif
      enddo
      endif
      enddo
      do i1=1,ntft
      if(stib(tfta(0)+i1).eq.0)then
      stib(tfta(0)+i1)=1
      endif
      enddo
      ii=1
      do i1=1,nleg
      if(stib(antiq(0)+leg(i1)).ne.0)then
      ii=-ii
      endif
      enddo
      if(ii.lt.0)then
      call wput(6,0,0)
      jflag(1)=1
      endif
      if(nleg.gt.1)then
      ii=0
      jj=stib(blok(0)+leg(1))
      do i1=2,nleg
      if(ii.eq.0)then
      if(stib(blok(0)+leg(i1)).ne.jj)then
      ii=1
      call wput(5,0,0)
      jflag(1)=1
      endif
      endif
      enddo
      endif
      if((nleg.ne.2).or.(nloop.ne.0))then
      do i1=1,nleg
      jj=stib(link(0)+leg(i1))
      do i2=1,nrho
      if(stib(nivd(0)+i2).gt.0)then
      if(stib(stib(stib(dpntro(0)+i2)+jj)+1).eq.jj)then
      goto 40
      endif
      endif
      enddo
      jflag(1)=1
      kk=0
      if(i1.gt.incom)then
      ii=stib(namep(0)+leg(i1))
      do i2=incom+1,i1-1
      if(leg(i2).eq.leg(i1))then
      kk=1
      endif
      enddo
      if(kk.eq.0)then
      j1=ii-1+stib(namel(0)+leg(i1))
      call wput(4,ii,j1)
      endif
      else
      ii=stib(namep(0)+jj)
      do i2=1,i1-1
      if(leg(i2).eq.leg(i1))then
      kk=1
      endif
      enddo
      if(kk.eq.0)then
      j1=ii-1+stib(namel(0)+jj)
      call wput(3,ii,j1)
      endif
      endif
   40 continue
      enddo
      endif
      if(nphi.eq.1)then
      mflag(1)=1
      endif
      ii=0
      jj=0
      do i1=1,nphi
      if(stib(antiq(0)+i1).eq.0)then
      ii=1
      else
      jj=1
      endif
      enddo
      if(ii.eq.0)then
      mflag(2)=-1
      elseif(jj.eq.0)then
      mflag(2)=1
      endif
      mflag(3)=1
      j1=0
   02 continue
      if(j1.lt.nvert)then
      j1=j1+1
      j2=stib(vparto(0)+j1)
      ii=-1
      do i2=1,stib(vval(0)+j1)
      ii=ii+stib(antiq(0)+stib(j2+i2))
      enddo
      if(abs(ii).eq.1)then
      goto 02
      endif
      mflag(3)=0
      endif
      if(mflag(1).eq.0)then
      jflag(3)=0
      elseif(mflag(4).ne.0)then
      jflag(3)=0
      else
      jflag(3)=1
      endif
      goto 90
  521 continue
      mlin(1:srec)='wrong syntax,'
  522 continue
      if(phase.eq.1)then
      call mput(1,dbl,nlin,1)
      else
      j1=1
      j2=0
  523 continue
      if(stib(drecii(0)+j1).eq.0)then
      j2=j2+1
      nlin=stib(dreci(0)+j1)
      endif
      if((j1.lt.ndrec).and.(j2.lt.icom))then
      j1=j1+1
      goto 523
      endif
  524 continue
      if(j1.lt.ndrec)then
      if(stib(drecii(0)+j1+1).ne.0)then
      j1=j1+1
      goto 524
      endif
      endif
      call mput(1,nlin,nlin+stib(drecii(0)+j1),1)
      endif
   90 return
      end
      subroutine iki
      implicit integer(a-z)
      save
      parameter ( sibuff=1048574, scbuff=131064, swbuff=65528 )
      parameter ( srec=81, ssrec=62 )
      character*(srec) mlin
      common/z22g/mlin
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z35g/stibs(1:1),stcbs(1:1),stwbs(1:1)
      common/z36g/udkp(0:2),udkl(0:0),udkt(0:0),udki(0:0)
      common/z38g/gmkp(0:0),gmkl(0:0),gmkd(0:0),gmko(0:0),
     :gmkvp(0:0),gmkvl(0:0)
      common/z39g/pmkp(0:0),pmkl(0:0),pmkd(0:0),pmkvpp(0:0),pmkvlp(0:0)
      common/z40g/vmkp(0:0),vmkl(0:0),vmkvpp(0:0),vmkvlp(0:0),vmks(0:0)
      common/z41g/nudk,ngmk,npmk0,npmk,nvmk0,nvmk
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
      mlin(1:srec)='function used in style-file was not found in'
      call mput(1,0,0,2)
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
      mlin(1:srec)="global M-functions cannot have 'dual-' as prefix"
      call mput(1,0,0,4)
      endif
  200 continue
      elseif(stib(udkt(0)+i1).eq.2)then
      if(stib(pmkd(0)+stib(udki(0)+i1)).eq.3)then
      if(stib(j+3).ne.0)then
      mlin(1:srec)='field M-functions cannot be used outside a loop,'
      call mput(1,0,0,4)
      elseif(stib(j+6).ne.0)then
      mlin(1:srec)='field M-functions cannot be used in the '//
     :'(strict) vertex loop,'
      call mput(1,0,0,4)
      endif
      else
      if(stib(j+3).ne.0)then
      mlin(1:srec)='propagator M-functions cannot be used '//
     :'outside a loop,'
      call mput(1,0,0,4)
      elseif(stib(j+6).ne.0)then
      mlin(1:srec)='propagator M-functions cannot be used in the '//
     :'(strict) vertex loop,'
      call mput(1,0,0,4)
      endif
      do 300 i2=4,7
      if(stib(j+i2).gt.1)then
      mlin(1:srec)='propagator M-functions cannot be prefixed with '//
     :'"dual-",'
      call mput(1,0,0,4)
      endif
  300 continue
      endif
      elseif(stib(udkt(0)+i1).eq.3)then
      if(stib(j+3).ne.0)then
      mlin(1:srec)='vertex M-functions cannot be used outside a loop,'
      call mput(1,0,0,4)
      else
      do 400 i2=4,7
      if((i2.ne.5).and.(stib(j+i2).gt.1))then
      mlin(1:srec)='vertex M-functions cannot be prefixed with '//
     :'"dual-" except in the propagator_loop,'
      call mput(1,0,0,4)
      endif
  400 continue
      endif
      endif
  500 continue
      return
      end
      subroutine rsfki
      implicit integer(a-z)
      save
      parameter ( sibuff=1048574, scbuff=131064, swbuff=65528 )
      parameter ( eoa=-2047, nap=-8191 )
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z35g/stibs(1:1),stcbs(1:1),stwbs(1:1)
      common/z36g/udkp(0:2),udkl(0:0),udkt(0:0),udki(0:0)
      common/z41g/nudk,ngmk,npmk0,npmk,nvmk0,nvmk
      ii=nudk+1
      udkp(2)=stibs(1)
      call aoib(ii)
      stib(stibs(1))=eoa
      udkp(0)=stibs(1)
      call aoib(ii)
      stib(stibs(1))=eoa
      udkl(0)=stibs(1)
      call aoib(ii)
      stib(stibs(1))=eoa
      udkt(0)=stibs(1)
      call aoib(ii)
      stib(stibs(1))=eoa
      udki(0)=stibs(1)
      call aoib(ii)
      stib(stibs(1))=eoa
      j1=udkp(1)
      do i1=1,nudk
      stib(udkt(0)+i1)=0
      stib(udki(0)+i1)=0
      j1=stib(j1)
      stib(udkp(2)+i1)=j1+1
      stib(udkp(0)+i1)=stib(j1+1)
      stib(udkl(0)+i1)=stib(j1+2)
      enddo
      return
      end
      subroutine rvmki(llp,llq)
      implicit integer(a-z)
      save
      parameter ( srec=81, ssrec=62 )
      parameter ( sibuff=1048574, scbuff=131064, swbuff=65528 )
      parameter ( eoa=-2047, nap=-8191 )
      character*(srec) mlin
      common/z22g/mlin
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z32g/aaf(0:4)
      common/z35g/stibs(1:1),stcbs(1:1),stwbs(1:1)
      common/z40g/vmkp(0:0),vmkl(0:0),vmkvpp(0:0),vmkvlp(0:0),vmks(0:0)
      common/z41g/nudk,ngmk,npmk0,npmk,nvmk0,nvmk
      common/z43g/vmkr(0:0),vmkmao(0:0),vmkmio(0:0)
      common/z57g/vmkvp(0:0),vmkvl(0:0)
      ii=0
      if(llp.gt.llq)then
      ii=1
      elseif(llq.lt.stibs(1))then
      ii=1
      elseif(aaf(3).ne.5)then
      ii=1
      elseif(aaf(4).ne.5)then
      ii=1
      endif
      if(ii.ne.0)then
      mlin(1:srec)='rvmki_1'
      call mput(1,0,0,0)
      endif
      if(nvmk.eq.0)then
      vmkp(0)=nap
      vmkl(0)=nap
      vmkr(0)=nap
      vmkvp(0)=nap
      vmkvl(0)=nap
      goto 90
      endif
      jj=aaf(3)
      j1=stibs(1)+1
      call xipht(j1,llq,jj)
      llq=llq+jj
      call aoib(jj)
      do i1=stibs(1)-jj+1,stibs(1)
      stib(i1)=eoa
      enddo
      ii=nvmk+1
      call trm(jj,ii)
      vmkvl(0)=stibs(1)-ii
      vmkvp(0)=vmkvl(0)-ii
      vmkr(0)=vmkvp(0)-ii
      vmkl(0)=vmkr(0)-ii
      vmkp(0)=vmkl(0)-ii
   90 continue
      do i1=1,nvmk
      j1=stib(vmkp(0)+i1)
      j2=j1-1+stib(vmkl(0)+i1)
      j3=stib(vmkvp(0)+i1)
      j4=j3-1+stib(vmkvl(0)+i1)
      if(stcb(j3:j3).eq."'")then
      j3=j3+1
      j4=j4+1
      endif
      enddo
      return
      end
      subroutine rpmki(llp,llq)
      implicit integer(a-z)
      save
      parameter ( srec=81, ssrec=62 )
      parameter ( sibuff=1048574, scbuff=131064, swbuff=65528 )
      parameter ( eoa=-2047, nap=-8191 )
      character*(srec) mlin
      common/z22g/mlin
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z32g/aaf(0:4)
      common/z35g/stibs(1:1),stcbs(1:1),stwbs(1:1)
      common/z39g/pmkp(0:0),pmkl(0:0),pmkd(0:0),pmkvpp(0:0),pmkvlp(0:0)
      common/z41g/nudk,ngmk,npmk0,npmk,nvmk0,nvmk
      common/z42g/pmkr(0:0),pmkvma(0:0),pmkvmi(0:0)
      common/z53g/pmkvp(0:0),pmkvl(0:0)
      ii=0
      if(llp.gt.llq)then
      ii=1
      elseif(llq.lt.stibs(1))then
      ii=1
      elseif(aaf(1).ne.6)then
      ii=1
      elseif(aaf(2).ne.6)then
      ii=1
      endif
      if(ii.ne.0)then
      mlin(1:srec)='rpmki_1'
      call mput(1,0,0,0)
      endif
      if(npmk.eq.0)then
      pmkp(0)=nap
      pmkl(0)=nap
      pmkr(0)=nap
      pmkvp(0)=nap
      pmkvl(0)=nap
      pmkd(0)=nap
      goto 90
      endif
      jj=aaf(1)
      j1=stibs(1)+1
      call xipht(j1,llq,jj)
      llq=llq+jj
      call aoib(jj)
      do i1=stibs(1)-jj+1,stibs(1)
      stib(i1)=eoa
      enddo
      ii=npmk+1
      call trm(jj,ii)
      pmkd(0)=stibs(1)-ii
      pmkvl(0)=pmkd(0)-ii
      pmkvp(0)=pmkvl(0)-ii
      pmkr(0)=pmkvp(0)-ii
      pmkl(0)=pmkr(0)-ii
      pmkp(0)=pmkl(0)-ii
   90 continue
      do i1=1,npmk
      j1=stib(pmkp(0)+i1)
      j2=j1-1+stib(pmkl(0)+i1)
      j3=stib(pmkvp(0)+i1)
      j4=j3-1+stib(pmkvl(0)+i1)
      if(stcb(j3:j3).eq."'")then
      j3=j3+1
      j4=j4+1
      endif
      enddo
      return
      end
      subroutine vsig
      implicit integer(a-z)
      save
      parameter ( maxleg=10, maxrho=7, maxi=10, maxdeg=6 )
      parameter ( maxn=maxi+maxi-2, maxli=maxn+maxrho )
      parameter ( sibuff=1048574, scbuff=131064, swbuff=65528 )
      parameter ( srec=81, ssrec=62 )
      parameter ( eoa=-2047, nap=-8191 )
      common/z1g/g(1:maxn,1:maxn),rho(1:maxdeg),nrho,rhop1,rhop2
      common/z3g/nivd(0:0),dpntro(0:0),vparto(0:0),vval(0:0),nvert
      common/z11g/nphi,nblok,nprop,npprop
      common/z12g/cflag(1:7),dflag(1:20),jflag(1:10),mflag(1:6)
      character*(srec) mlin
      common/z22g/mlin
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z35g/stibs(1:1),stcbs(1:1),stwbs(1:1)
      common/z40g/vmkp(0:0),vmkl(0:0),vmkvpp(0:0),vmkvlp(0:0),vmks(0:0)
      common/z41g/nudk,ngmk,npmk0,npmk,nvmk0,nvmk
      common/z43g/vmkr(0:0),vmkmao(0:0),vmkmio(0:0)
      if(dflag(15).eq.0)then
      vmkmio(0)=nap
      vmkmao(0)=nap
      goto 91
      endif
      ii=nvmk+1
      vmkmio(0)=stibs(1)
      call aoib(ii)
      stib(stibs(1))=eoa
      vmkmao(0)=stibs(1)
      call aoib(ii)
      stib(stibs(1))=eoa
      do i1=1,nvmk
      if(stib(vmks(0)+i1).ne.1)then
      stib(vmkmio(0)+i1)=nap
      stib(vmkmao(0)+i1)=nap
      else
      ii=nrho
      p1=stibs(1)
      stib(vmkmio(0)+i1)=p1
      call aoib(ii)
      p2=stibs(1)
      stib(vmkmao(0)+i1)=p2
      call aoib(ii)
      call vaoib(ii)
      do i2=stibs(1)+1,stibs(1)+ii
      stib(i2)=0
      enddo
      ii=stib(vmkvpp(0)+i1)
      jj=stib(vmkvlp(0)+i1)
      do i2=1,nvert
      j1=stib(ii+i2)
      j2=j1-1+stib(jj+i2)
      ij=stoz(stcb,j1,j2)
      j1=stib(vval(0)+i2)
      j2=stibs(1)+j1
      if(stib(j2).ne.0)then
      if(ij.gt.stib(p2+j1))then
      stib(p2+j1)=ij
      elseif(ij.lt.stib(p1+j1))then
      stib(p1+j1)=ij
      endif
      else
      stib(j2)=1
      stib(p1+j1)=ij
      stib(p2+j1)=ij
      endif
      enddo
      endif
      enddo
      if(stib(stibs(1)).ne.eoa)then
      call aoib(1)
      stib(stibs(1))=eoa
      endif
   91 continue
      return
      end
      subroutine gen10(cntr)
      implicit integer(a-z)
      save
      parameter ( maxleg=10, maxrho=7, maxi=10, maxdeg=6 )
      parameter ( maxn=maxi+maxi-2, maxli=maxn+maxrho )
      parameter ( ipar1= maxleg*maxleg-maxleg )
      parameter ( ipar2= ipar1/2 )
      parameter ( sibuff=1048574, scbuff=131064, swbuff=65528 )
      parameter ( srec=81, ssrec=62 )
      common/z1in/leg(1:maxleg),nleg,nloop,incom
      common/z1g/g(1:maxn,1:maxn),rho(1:maxdeg),nrho,rhop1,rhop2
      common/z4g/n,nli
      common/z6g/p1(1:maxleg),invp1(1:maxleg)
      common/z7g/lmap(1:maxn,1:maxdeg),vmap(1:maxn,1:maxdeg),
     :pmap(1:maxn,1:maxdeg),vlis(1:maxn),invlis(1:maxn)
      common/z8g/vdeg(1:maxn),xn(1:maxn)
      common/z10g/p1l(1:ipar2),p1r(1:ipar2),ns1
      common/z12g/cflag(1:7),dflag(1:20),jflag(1:10),mflag(1:6)
      common/z16g/rdeg(1:maxn),amap(1:maxn,1:maxdeg)
      common/z18g/eg(1:maxn,1:maxn),flow(1:maxli,0:maxleg+maxrho),
     :net(-3:3)
      common/z20g/tftyp(0:0),tfnarg(0:0),tfa(0:0),tfb(0:0),tfc(0:0),
     :tfo(0:0),tf2(0:0),ntf
      character*(srec) mlin
      common/z22g/mlin
      common/z30g/stib(1:sibuff)
      common/z50g/tfta(0:0),tftb(0:0),tftic(0:0),ntft
      integer vaux(1:maxli),xli(1:maxleg),xlj(1:maxn)
      if(cntr.eq.0)then
      if(dflag(13).eq.0)then
      goto 29
      else
      goto 27
      endif
      endif
      cntr21=1
   27 continue
      call gen21(cntr21)
      if(cntr21.eq.0)then
      cntr=0
      goto 90
      endif
      do i1=1,n
      vaux(i1)=xn(i1)
      enddo
      do i1=1,rho(1)
      vlis(i1)=i1
      invlis(i1)=i1
      enddo
      do i1=rhop1,n
      invlis(i1)=0
      enddo
      aux=rho(1)
      jj=rhop1
   10 continue
      if(aux.lt.n)then
   20 continue
      if(invlis(jj).ne.0)then
      jj=jj+1
      goto 20
      endif
      ii=jj
      do i1=ii+1,n
      if(invlis(i1).eq.0)then
      if(vaux(i1).gt.vaux(ii))then
      ii=i1
      endif
      endif
      enddo
      if(vaux(ii).eq.0)then
      mlin(1:srec)='gen10_1'
      call mput(1,0,0,0)
      endif
      aux=aux+1
      vlis(aux)=ii
      invlis(ii)=aux
      rdeg(ii)=vaux(ii)
      do i1=rhop1,n
      vaux(i1)=vaux(i1)+g(i1,ii)
      enddo
      goto 10
      endif
      do i1=1,n
      vaux(i1)=0
      enddo
      do i1=1,n
      ii=vlis(i1)
      kk=0
      j1=1
   30 continue
      if(kk.lt.vdeg(ii))then
      jj=vlis(j1)
      aux=1
      do i3=1,g(ii,jj)
      kk=kk+1
      vmap(ii,kk)=jj
      vaux(jj)=vaux(jj)+1
      if(ii.ne.jj)then
      lmap(ii,kk)=vaux(jj)
      else
      lmap(ii,kk)=vaux(jj)+aux
      aux=-aux
      endif
      enddo
      j1=j1+1
      goto 30
      endif
      if(kk.ne.vdeg(ii))then
      mlin(1:srec)='gen10_2'
      call mput(1,0,0,0)
      endif
      enddo
      do i1=rhop1,n
      do i2=1,vdeg(i1)-1
      if(invlis(vmap(i1,i2)).gt.invlis(vmap(i1,i2+1)))then
      mlin(1:srec)='gen10_3'
      call mput(1,0,0,0)
      endif
      enddo
      enddo
      do i1=1,n
      do i2=1,vdeg(i1)
      j1=vmap(i1,i2)
      j2=lmap(i1,i2)
      if(vmap(j1,j2).ne.i1)then
      mlin(1:srec)='gen10_4'
      call mput(1,0,0,0)
      endif
      if(lmap(j1,j2).ne.i2)then
      mlin(1:srec)='gen10_5'
      call mput(1,0,0,0)
      endif
      enddo
      enddo
      if(cflag(2)+jflag(4).ne.0)then
      do i1=rhop1,n
      do i2=1,vdeg(i1)
      j3=vmap(i1,i2)
      if(j3.gt.rho(1))then
      j1=min(i1,j3)
      j2=max(i1,j3)
      ii=eg(j1,j2)
      jj=g(i1,j3)-1
      if(jj.eq.0)then
      amap(i1,i2)=ii
      elseif(i1.eq.j3)then
      amap(i1,i2)=ii+(i2-1-rdeg(i1))/2
      else
      kk=0
      do i3=i2-1,1,-1
      if(vmap(i1,i3).eq.j3)then
      kk=kk+1
      else
      goto 37
      endif
      enddo
   37 continue
      amap(i1,i2)=ii+kk
      endif
      endif
      enddo
      enddo
      endif
      if(dflag(9).ne.0)then
      call ccyc(1,j1)
      if(j1.ne.dflag(9))then
      goto 27
      endif
      endif
      do i1=1,rho(1)
      p1(i1)=i1
      enddo
      goto 121
   29 continue
      do i1=rho(1)-1,1,-1
      if(p1(i1).lt.p1(i1+1))then
      j1=i1
      goto 38
      endif
      enddo
      goto 27
   38 continue
      do i1=j1+2,rho(1)
      if(p1(i1).lt.p1(j1))then
      j2=i1-1
      goto 61
      endif
      enddo
      j2=rho(1)
   61 continue
      ii=p1(j1)
      p1(j1)=p1(j2)
      p1(j2)=ii
      j1=j1+1
      j2=rho(1)
   33 continue
      if(j1.lt.j2)then
      ii=p1(j1)
      p1(j1)=p1(j2)
      p1(j2)=ii
      j1=j1+1
      j2=j2-1
      goto 33
      endif
  121 continue
      do i1=1,ns1
      if(p1(p1r(i1)).lt.p1(p1l(i1)))then
      do i2=p1r(i1)+1,rho(1)-1
      do i3=i2+1,rho(1)
      if(p1(i2).lt.p1(i3))then
      ii=p1(i2)
      p1(i2)=p1(i3)
      p1(i3)=ii
      endif
      enddo
      enddo
      goto 29
      endif
      enddo
      do i1=1,rho(1)
      invp1(p1(i1))=i1
      enddo
      if(dflag(17).ne.0)then
      ii=stib(tftic(0)+11)
      do i1=stib(tfta(0)+ii),stib(tftb(0)+ii)
      i3=stib(tf2(0)+i1)
      styp=stib(tftyp(0)+i3)
      j1=stib(tfnarg(0)+i3)
      j2=stib(tfo(0)+i3)
      do i2=rhop1,n
      xlj(i2)=0
      enddo
      ii=0
      do i2=1,j1
      j3=vmap(invp1(stib(j2+i2)),1)
      if(xlj(j3).eq.0)then
      xlj(j3)=1
      ii=ii+1
      endif
      enddo
      if((ii.lt.stib(tfa(0)+i3)).or.(ii.gt.stib(tfb(0)+i3)))then
      if(styp.gt.0)then
      goto 29
      else
      goto 55
      endif
      endif
      if(stib(tfc(0)+i3).eq.0)then
      if(styp.gt.0)then
      goto 55
      else
      goto 29
      endif
      endif
      do i2=1,rho(1)
      xli(i2)=0
      enddo
      do i2=1,j1
      xli(invp1(stib(j2+i2)))=1
      enddo
      do i2=1,rho(1)
      if(xli(i2).eq.0)then
      if(xlj(vmap(i2,1)).ne.0)then
      if(styp.gt.0)then
      goto 29
      else
      goto 55
      endif
      endif
      endif
      enddo
      if(styp.lt.0)then
      goto 29
      endif
   55 continue
      enddo
      endif
      if(dflag(18).ne.0)then
      ii=stib(tftic(0)+12)
      do i2=stib(tfta(0)+ii),stib(tftb(0)+ii)
      j3=stib(tf2(0)+i2)
      styp=stib(tftyp(0)+j3)
      j1=stib(tfnarg(0)+j3)
      j2=stib(tfo(0)+j3)
      do i3=1,rho(1)
      xli(i3)=0
      enddo
      do i3=1,j1
      xli(invp1(stib(j2+i3)))=1
      enddo
      do i3=rhop1,nli
      if(flow(i3,0).eq.1)then
      kk=0
      do i4=1,rho(1)
      if(flow(i3,i4).eq.0)then
      if(xli(i4).ne.0)then
      kk=kk+1
      endif
      else
      if(xli(i4).eq.0)then
      kk=kk+1
      endif
      endif
      enddo
      if(kk.gt.0)then
      kk=kk-rho(1)
      endif
      if(kk.eq.0)then
      if(styp.gt.0)then
      goto 57
      else
      goto 29
      endif
      endif
      endif
      enddo
      if(styp.gt.0)then
      goto 29
      endif
   57 continue
      enddo
      endif
      cntr=1
      if(cflag(2).eq.0)then
      goto 90
      endif
      do i1=1,rho(1)
      amap(i1,1)=i1
      amap(vmap(i1,1),lmap(i1,1))=i1
      enddo
      do i1=1,nli
      vaux(i1)=-2
      enddo
      do i1=1,n
      do i2=1,vdeg(i1)
      ii=amap(i1,i2)
      if(ii.gt.0)then
      if(ii.le.nli)then
      vaux(ii)=vaux(ii)+1
      endif
      endif
      enddo
      enddo
      do i1=1,nli
      if(vaux(i1).ne.0)then
      mlin(1:srec)='gen10_6'
      call mput(1,0,0,0)
      endif
      enddo
   90 return
      end
      subroutine gen21(cntr21)
      implicit integer(a-z)
      save
      parameter ( maxleg=10, maxrho=7, maxi=10, maxdeg=6 )
      parameter ( maxn=maxi+maxi-2, maxli=maxn+maxrho )
      parameter ( ipar1= maxleg*maxleg-maxleg )
      parameter ( ipar2= ipar1/2 )
      parameter ( sibuff=1048574, scbuff=131064, swbuff=65528 )
      parameter ( srec=81, ssrec=62 )
      parameter ( eoa=-2047, nap=-8191 )
      common/z1g/g(1:maxn,1:maxn),rho(1:maxdeg),nrho,rhop1,rhop2
      common/z4g/n,nli
      common/z5g/psym(0:0),psyms,nsym
      common/z7g/lmap(1:maxn,1:maxdeg),vmap(1:maxn,1:maxdeg),
     :pmap(1:maxn,1:maxdeg),vlis(1:maxn),invlis(1:maxn)
      common/z8g/vdeg(1:maxn),xn(1:maxn)
      common/z10g/p1l(1:ipar2),p1r(1:ipar2),ns1
      common/z12g/cflag(1:7),dflag(1:20),jflag(1:10),mflag(1:6)
      common/z14g/zcho(0:maxli),zbri(0:maxli),zpro(0:maxli),
     :rbri(0:maxli),sbri(0:maxli)
      common/z17g/xtail(1:maxn),xhead(1:maxn),ntadp
      common/z18g/eg(1:maxn,1:maxn),flow(1:maxli,0:maxleg+maxrho),
     :net(-3:3)
      character*(srec) mlin
      common/z22g/mlin
      common/z30g/stib(1:sibuff)
      common/z35g/stibs(1:1),stcbs(1:1),stwbs(1:1)
      integer xc(1:maxdeg)
      integer xl(1:maxn),xt(1:maxn)
      integer bound(1:maxdeg),degr(1:maxdeg),xs(1:maxn)
      integer xg(1:maxn,1:maxn),ds(1:maxn,1:maxn)
      integer p1s(1:maxleg,1:maxleg)
      integer aa(1:maxn),bb(1:maxn)
      integer uset(0:maxn),xset(1:maxn),xp(1:maxn),a1(1:maxn)
      integer str(1:maxn),dta(1:maxn),lps(1:maxn)
      integer head(1:maxli),tail(1:maxli),intree(1:maxli)
      integer emul(1:maxdeg,1:maxli),nemul(1:maxdeg)
      if(abs(cntr21).ne.1)then
      goto 97
      elseif(nrho.le.0)then
      goto 97
      endif
      do i1=1,nrho
      if((rho(i1).lt.0).or.(rho(i1).gt.maxn))then
      goto 97
      endif
      enddo
      if(cntr21.eq.-1)then
      goto 05
      endif
      n=0
      jj=0
      do i1=1,nrho
      jj=jj+i1*rho(i1)
      do i2=1,rho(i1)
      n=n+1
      vdeg(n)=i1
      enddo
      enddo
      if((n.le.0).or.(n.gt.maxn))then
      goto 97
      endif
      nli=jj/2
      if(jj.ne.nli+nli)then
      goto 97
      endif
      loop=nli-n+1
      rhop2=rho(1)+loop
      if(loop.lt.0)then
      goto 97
      endif
      j1=0
      if(dflag(12).gt.0)then
      j1=1
      endif
      do i1=rhop1,n
      if(j1.ne.0)then
      str(i1)=1
      elseif(vdeg(i1).ne.vdeg(n))then
      str(i1)=min(vdeg(i1),loop+1)
      elseif(n.gt.2)then
      str(i1)=min(vdeg(i1)-1,loop+1)
      else
      str(i1)=min(vdeg(i1),loop+1)
      endif
      enddo
      j1=1
      xl(1)=1
      do i1=2,n
      if(vdeg(i1).ne.vdeg(i1-1))then
      xt(j1)=i1-1
      j1=j1+1
      xl(j1)=i1
      endif
      enddo
      xt(j1)=n
      cntr21=-1
      do i1=1,n
      xn(i1)=0
      enddo
      if(rho(1).eq.0)then
      xc(1)=0
      goto 59
      endif
      ubo=0
      n1=1
      do i1=2,nrho
      if(rho(i1).gt.0)then
      n1=n1+1
      degr(n1)=i1
      bound(n1)=i1*rho(i1)
      ubo=ubo+bound(n1)
      endif
      enddo
      xc(1)=max(0,rho(1)-ubo)
   61 continue
      ii=rho(1)-xc(1)
      do i1=n1,2,-1
      xc(i1)=min(ii,bound(i1))
      ii=ii-xc(i1)
      enddo
   43 continue
      do i1=n1,2,-1
      ii=xc(i1)
      do i2=xl(i1),xt(i1)
      xn(i2)=min(ii,degr(i1))
      ii=ii-xn(i2)
      xs(i2)=ii
      enddo
      enddo
      goto 12
   76 continue
      do i1=n1,2,-1
      do i2=xt(i1)-1,xl(i1),-1
      if(xn(i2).gt.1)then
      j1=i2
      goto 89
      endif
      enddo
      goto 82
   89 continue
      xn(j1)=xn(j1)-1
      xs(j1)=xs(j1)+1
      do i2=j1+1,xt(i1)
      xn(i2)=min(xn(j1),xs(i2-1))
      xs(i2)=xs(i2-1)-xn(i2)
      enddo
      if(xs(xt(i1)).gt.0)then
      do i2=j1-1,xl(i1),-1
      if(xn(i2).gt.1)then
      j1=i2
      goto 89
      endif
      enddo
      goto 82
      endif
      do i2=i1+1,n1
      ii=xc(i2)
      do i3=xl(i2),xt(i2)
      xn(i3)=min(ii,degr(i2))
      ii=ii-xn(i3)
      xs(i3)=ii
      enddo
      enddo
      goto 12
   82 continue
      enddo
      do i1=n1,3,-1
      if(xc(i1).gt.0)then
      j1=i1
      goto 45
      endif
      enddo
      goto 85
   45 continue
      do i1=j1-1,2,-1
      if(xc(i1).lt.bound(i1))then
      j1=i1
      goto 55
      endif
      enddo
      goto 85
   55 continue
      xc(j1)=xc(j1)+1
      ii=-1
      do i1=j1+1,n1
      ii=ii+xc(i1)
      enddo
      do i1=n1,j1+1,-1
      xc(i1)=min(ii,bound(i1))
      ii=ii-xc(i1)
      enddo
      if(ii.ne.0)then
      mlin(1:srec)='gen21_2'
      call mput(1,0,0,0)
      endif
      goto 43
   85 continue
      xc(1)=xc(1)+2
      if(xc(1).le.rho(1))then
      goto 61
      endif
      goto 80
   12 continue
      if(rho(1).lt.n)then
      if(xc(1).gt.0)then
      goto 80
      endif
      endif
      if(n.gt.rhop1)then
      j1=0
      if(dflag(1).gt.0)then
      j1=1
      endif
      do i1=rhop1,n
      if(xn(i1)+j1.ge.vdeg(i1))then
      goto 76
      endif
      enddo
      endif
      ii=1
      if(loop.eq.0)then
      ii=0
      elseif((dflag(8).gt.0).and.(rhop1.ne.2))then
      ii=0
      elseif(dflag(11).gt.0)then
      ii=0
      elseif(dflag(12).gt.0)then
      ii=0
      else
      if(dflag(9).gt.0)then
      if(loop.ne.1)then
      ii=0
      endif
      endif
      endif
      if(ii.eq.0)then
      do i1=rhop1,n
      dta(i1)=0
      enddo
      else
      do i1=rhop1,n
      ii=0
      if(n.gt.rhop1)then
      if(dflag(1).gt.0)then
      ii=2
      elseif((dflag(3).gt.0).and.(xn(i1).eq.0))then
      ii=2
      else
      ii=1
      endif
      endif
      j=max(0,min(vdeg(i1)-xn(i1)-ii,loop+loop))/2
      dta(i1)=j+j
      enddo
      endif
   59 continue
      do i1=1,n
      do i2=i1+1,n
      xg(i1,i2)=0
      enddo
      enddo
      do i1=2,xc(1),2
      xg(i1-1,i1)=1
      enddo
      limin=rhop1
      limax=max(limin,n-1)
      lin=limin
      j1=xc(1)
      do i1=limin,n
      if(xn(i1).gt.0)then
      a1(i1)=j1+1
      do i2=1,xn(i1)
      j1=j1+1
      xg(j1,i1)=1
      vmap(j1,1)=i1
      enddo
      else
      a1(i1)=0
      endif
      enddo
      dsum=-1
   70 continue
      dsum=dsum+1
      if(dsum.gt.loop)then
      goto 76
      endif
      ii=2*dsum
      do i1=n,rhop1,-1
      jj=min(ii,dta(i1))
      xg(i1,i1)=jj
      ii=ii-jj
      enddo
      if(ii.eq.0)then
      goto 19
      else
      goto 76
      endif
   32 continue
      do i1=n,rhop1,-1
      if(xg(i1,i1).gt.0)then
      j1=i1
      goto 33
      endif
      enddo
      goto 70
   33 continue
      do i1=j1-1,rhop1,-1
      if(xg(i1,i1).lt.dta(i1))then
      xg(i1,i1)=xg(i1,i1)+2
      j2=i1
      goto 64
      endif
      enddo
      goto 70
   64 continue
      ii=-2
      do i1=j2+1,j1
      ii=ii+xg(i1,i1)
      enddo
      do i1=n,j2+1,-1
      jj=min(ii,dta(i1))
      xg(i1,i1)=jj
      ii=ii-jj
      enddo
      if(ii.ne.0)then
      mlin(1:srec)='gen21_3'
      call mput(1,0,0,0)
      endif
   19 continue
      do i1=limin,n
      ds(limin,i1)=vdeg(i1)-xn(i1)-xg(i1,i1)
      enddo
      uset(0)=0
      xset(1)=1
      jj=1
      do i1=2,n
      if(vdeg(i1-1).ne.vdeg(i1))then
      uset(jj)=i1-1
      jj=jj+1
      elseif(xn(i1-1).ne.xn(i1))then
      uset(jj)=i1-1
      jj=jj+1
      else
      ii=xg(i1-1,i1-1)-xg(i1,i1)
      if(ii.gt.0)then
      uset(jj)=i1-1
      jj=jj+1
      elseif(ii.lt.0)then
      goto 32
      endif
      endif
      xset(i1)=jj
      enddo
      uset(jj)=n
      lps(lin)=loop-dsum+1
   10 continue
      ii=ds(lin,lin)
      bond=min(str(lin),lps(lin))
      do i1=n,lin+1,-1
      jj=min(ii,bond,ds(lin,i1))
      xg(lin,i1)=jj
      ii=ii-jj
      enddo
      if(ii.gt.0)then
      goto 15
      endif
      goto 28
   05 continue
      lin=limax
      goto 17
   15 continue
      if(lin.eq.limin)then
      goto 32
      endif
      lin=lin-1
   17 continue
      if((lin.lt.limin).or.(lin.gt.n))then
      mlin(1:srec)='gen21_5'
      call mput(1,0,0,0)
      endif
      do col=n,lin+1,-1
      aux=xg(lin,col)-1
      if(aux.ge.0)then
      goto 23
      endif
      enddo
      goto 15
   23 continue
      bond=min(str(lin),lps(lin))
      do i1=col-1,lin+1,-1
      if(min(ds(lin,i1),bond).gt.xg(lin,i1))then
      xg(lin,i1)=xg(lin,i1)+1
      do i2=n,i1+1,-1
      xg(lin,i2)=min(aux,bond,ds(lin,i2))
      aux=aux-xg(lin,i2)
      enddo
      goto 28
      endif
      aux=aux+xg(lin,i1)
      enddo
      goto 15
   38 continue
      aux=-1
      do i1=col,n
      aux=aux+xg(lin,i1)
      enddo
      goto 23
   28 continue
      if(lin.eq.n)then
      goto 200
      endif
      msum=0
      do i1=lin+1,n
      ii=xg(lin,i1)-1
      if(ii.gt.0)then
      msum=msum+ii
      endif
      enddo
      if(msum.ge.lps(lin))then
      goto 17
      endif
      if(lin.gt.limin)then
      if(xset(lin).eq.xset(lin-1))then
      do i1=limin,lin-2
      ii=xg(i1,lin-1)-xg(i1,lin)
      if(ii.gt.0)then
      goto 35
      elseif(ii.lt.0)then
      mlin(1:srec)='gen21_7'
      call mput(1,0,0,0)
      endif
      enddo
      do col=lin+1,n
      ii=xg(lin-1,col)-xg(lin,col)
      if(ii.lt.0)then
      goto 38
      elseif(ii.gt.0)then
      goto 35
      endif
      enddo
      endif
      endif
   35 continue
      do col=lin+2,n
      if(xset(col).eq.xset(col-1))then
      do i1=limin,lin
      ii=xg(i1,col-1)-xg(i1,col)
      if(ii.lt.0)then
      goto 38
      elseif(ii.gt.0)then
      goto 24
      endif
      enddo
      endif
   24 continue
      enddo
      do i1=lin+1,n
      ds(lin+1,i1)=ds(lin,i1)-xg(lin,i1)
      if(ds(lin+1,i1).lt.0)then
      mlin(1:srec)='gen21_8'
      call mput(1,0,0,0)
      endif
      enddo
      lin=lin+1
      lps(lin)=lps(lin-1)-msum
      goto 10
  200 continue
      j1=rhop1
      j2=n-1
      lin=n
      do i1=j1,j2
      aa(i1)=0
      enddo
      aa(n)=1
      i1=n-rho(1)
      bb(i1)=lin
      kk=i1-1
   21 continue
      if(kk.gt.0)then
      ii=bb(i1)
      i1=i1-1
      do i2=j2,ii+1,-1
      if(aa(i2).eq.0)then
      if(xg(ii,i2).ne.0)then
      bb(kk)=i2
      kk=kk-1
      aa(i2)=-aa(ii)
      endif
      endif
      enddo
      if(aa(j2).ne.0)then
      j2=j2-1
      endif
      do i2=j1,ii-1
      if(aa(i2).eq.0)then
      if(xg(i2,ii).ne.0)then
      bb(kk)=i2
      kk=kk-1
      aa(i2)=-aa(ii)
      endif
      endif
      enddo
      if(aa(j1).ne.0)then
      j1=j1+1
      endif
      if(i1.eq.kk)then
      do i2=j2,j1,-1
      if(aa(i2).eq.0)then
      aa(i2)=1
      bb(kk)=i2
      kk=kk-1
      lin=i2
      j2=lin-1
      goto 21
      endif
      enddo
      mlin(1:srec)='gen21_9'
      call mput(1,0,0,0)
      endif
      goto 21
      endif
      if(lin.ne.n)then
      goto 17
      endif
      if(dflag(11).ne.0)then
      do i1=rhop1,n
      do i2=i1+1,n
      if(xg(i1,i2).ne.0)then
      if(aa(i1)+aa(i2).ne.0)then
      if(dflag(11).gt.0)then
      goto 05
      else
      goto 62
      endif
      endif
      endif
      enddo
      if(xg(i1,i1).ne.0)then
      if(dflag(11).gt.0)then
      goto 05
      else
      goto 62
      endif
      endif
      enddo
      if(dflag(11).lt.0)then
      goto 05
      endif
      endif
   62 continue
      if(dflag(12).lt.0)then
      do i1=rhop1,n
      do i2=i1,n
      if(xg(i1,i2)-1.gt.0)then
      goto 34
      endif
      enddo
      enddo
      goto 05
      endif
   34 continue
      do i1=1,n
      do i2=1,i1-1
      j1=xg(i2,i1)
      g(i1,i2)=j1
      g(i2,i1)=j1
      enddo
      g(i1,i1)=xg(i1,i1)
      enddo
      if(dflag(1).ne.0)then
      call umpi(1,ii)
      if(ii.ne.dflag(1))then
      goto 05
      endif
      endif
      if(dflag(3).ne.0)then
      call umpi(2,ii)
      if(ii.ne.dflag(3))then
      goto 05
      endif
      endif
      if(dflag(4).ne.0)then
      call umpi(4,ii)
      if(ii.ne.dflag(4))then
      goto 05
      endif
      endif
      if(dflag(5).ne.0)then
      call umvi(3,ii)
      if(ii.ne.dflag(5))then
      goto 05
      endif
      endif
      if(dflag(8).ne.0)then
      if(rho(1).ne.1)then
      call umvi(1,ii)
      else
      call umvi(3,ii)
      endif
      if(ii.ne.dflag(8))then
      goto 05
      endif
      endif
      if(dflag(10).ne.0)then
      call umvi(2,ii)
      if(ii.ne.dflag(10))then
      goto 05
      endif
      endif
      ntadp=0
      if(mflag(6).ne.0)then
      call umpi(3,ii)
      if(ii.ne.1)then
      mlin(1:srec)='gen21_11'
      call mput(1,0,0,0)
      endif
      endif
      nsym=0
      do i1=1,rho(1)-1
      do i2=i1+1,rho(1)
      p1s(i1,i2)=0
      enddo
      enddo
      do i1=1,n
      xp(i1)=i1
      enddo
      goto 93
   77 continue
      do i1=xset(n),1,-1
      do i2=uset(i1)-1,uset(i1-1)+1,-1
      if(xp(i2).lt.xp(i2+1))then
      goto 102
      endif
      enddo
      enddo
      goto 63
  102 continue
      j1=uset(i1)
      do i1=i2+2,j1
      if(xp(i1).lt.xp(i2))then
      goto 202
      endif
      enddo
      i1=j1+1
  202 continue
      i1=i1-1
      ii=xp(i1)
      xp(i1)=xp(i2)
      xp(i2)=ii
      i1=i2+1
      i2=j1
  204 continue
      if(i1.lt.i2)then
      ii=xp(i1)
      xp(i1)=xp(i2)
      xp(i2)=ii
      i1=i1+1
      i2=i2-1
      goto 204
      endif
      do i1=j1+1,n
      xp(i1)=i1
      enddo
   93 continue
      if(rho(1).gt.0)then
      if(xp(rho(1)).ne.rho(1))then
      goto 63
      endif
      endif
      do i1=rhop1,n-1
      do i2=i1,n
      ii=g(xp(i1),xp(i2))-g(i1,i2)
      if(ii.gt.0)then
      goto 05
      endif
      if(ii.lt.0)then
      j1=xset(i2)
      j3=uset(j1)
      do i3=i2+1,j3-1
      do i4=i3+1,j3
      if(xp(i3).lt.xp(i4))then
      ii=xp(i3)
      xp(i3)=xp(i4)
      xp(i4)=ii
      endif
      enddo
      enddo
      do i3=j1+1,xset(n)
      j2=j3+1
      j3=uset(i3)
      j4=j2+j3
      do i4=j2,j3
      xp(i4)=j4-i4
      enddo
      enddo
      goto 77
      endif
      enddo
      enddo
      if(rho(1).eq.0)then
      goto 114
      endif
      i1=1
  110 continue
      j1=vmap(i1,1)
      if(xp(j1).eq.j1)then
      i1=i1+xn(j1)
      if(i1.le.rho(1))then
      goto 110
      endif
      goto 114
      else
      p1s(i1,a1(xp(j1)))=1
      goto 77
      endif
  114 continue
      if(jflag(3).eq.0)then
      if(psym(0).lt.0)then
      psyms=0
      psym(0)=stibs(1)
      call aoib(1)
      stib(stibs(1))=eoa
      endif
      jj=n-rho(1)
      ii=nsym*jj-psyms
      if(ii.gt.0)then
      if(psym(0).eq.stibs(1)-1-psyms)then
      call aoib(ii)
      psyms=psyms+ii
      stib(stibs(1))=eoa
      else
      mlin(1:srec)='gen21_13'
      call mput(1,0,0,0)
      endif
      endif
      ii=nsym*jj
      if(ii.le.psyms)then
      if(nsym.gt.0)then
      ii=psym(0)+(ii-n)
      do i1=rhop1,n
      stib(ii+i1)=xp(i1)
      enddo
      endif
      else
      mlin(1:srec)='gen21_14'
      call mput(1,0,0,0)
      endif
      endif
      nsym=nsym+1
      goto 77
   63 continue
      ns1=0
      do i1=2,rho(1)
      j1=i1-1
      if(vmap(j1,1).eq.vmap(i1,1))then
      p1s(j1,i1)=1
      endif
      do i2=1,j1
      if(p1s(i2,i1).eq.1)then
      ns1=ns1+1
      p1l(ns1)=i2
      p1r(ns1)=i1
      endif
      enddo
      enddo
      if(cflag(2)+jflag(4).eq.0)then
      goto 90
      endif
      do i1=1,rho(1)
      tail(i1)=i1
      head(i1)=vmap(i1,1)
      enddo
      do i1=1,vdeg(n)
      nemul(i1)=0
      enddo
      ii=rhop1
      do i1=rhop1,n
      do i2=i1,n
      jj=g(i1,i2)
      if(jj.eq.0)then
      eg(i1,i2)=0
      else
      eg(i1,i2)=ii
      if(i1.eq.i2)then
      jj=jj/2
      else
      nemul(jj)=nemul(jj)+1
      emul(jj,nemul(jj))=ii
      endif
      ii=ii+jj
      do i3=ii-jj,ii-1
      tail(i3)=i1
      head(i3)=i2
      enddo
      endif
      enddo
      enddo
      if(jflag(4).eq.0)then
      goto 90
      endif
      do i1=rhop1,nli
      intree(i1)=0
      enddo
      do i1=1,n
      bb(i1)=0
      enddo
      do i1=rhop1,n
      aa(i1)=0
      enddo
      kk=0
      ntree=0
      ii=vdeg(n)
   51 continue
      if(ii.le.0)then
      goto 52
      endif
      jj=nemul(ii)
   58 continue
      if(jj.le.0)then
      ii=ii-1
      goto 51
      endif
      j1=tail(emul(ii,jj))
      j2=head(emul(ii,jj))
      if(aa(j1).eq.aa(j2))then
      if(aa(j1).eq.0)then
      ntree=ntree+1
      aa(j1)=ntree
      aa(j2)=ntree
      else
      jj=jj-1
      goto 58
      endif
      else
      if(aa(j1).eq.0)then
      aa(j1)=aa(j2)
      elseif(aa(j2).eq.0)then
      aa(j2)=aa(j1)
      else
      j3=aa(j2)
      do i1=rhop1,n
      if(aa(i1).eq.j3)then
      aa(i1)=aa(j1)
      endif
      enddo
      endif
      endif
      intree(eg(j1,j2))=1
      bb(j1)=bb(j1)+1
      bb(j2)=bb(j2)+1
      kk=kk+1
      if(kk.lt.n-rhop1)then
      jj=jj-1
      goto 58
      endif
   52 continue
      do i1=1,rho(1)
      intree(i1)=1
      enddo
      do i1=1,nli
      do i2=1,rhop2
      flow(i1,i2)=0
      enddo
      enddo
      do i1=1,rho(1)
      flow(i1,i1)=1
      enddo
      ii=rhop1
      do i1=rhop1,nli
      if(intree(i1).eq.0)then
      flow(i1,ii)=1
      ii=ii+1
      endif
      enddo
      ii=rhop1
   73 continue
      if(bb(ii).eq.1)then
      goto 71
      elseif(ii.lt.n)then
      ii=ii+1
      goto 73
      else
      goto 83
      endif
   71 continue
      do i1=rhop1,ii-1
      if(xg(i1,ii).ne.0)then
      j1=eg(i1,ii)
      if(intree(j1).ne.0)then
      goto 75
      endif
      endif
      enddo
      do i1=ii+1,n
      if(xg(ii,i1).ne.0)then
      j1=eg(ii,i1)
      if(intree(j1).ne.0)then
      goto 75
      endif
      endif
      enddo
      mlin(1:srec)='gen21_16'
      call mput(1,0,0,0)
   75 continue
      if(rho(1).gt.1)then
      if(xn(ii).gt.0)then
      do i1=1,rho(1)
      if(vmap(i1,1).eq.ii)then
      if(tail(j1).eq.ii)then
      flow(j1,i1)=flow(j1,i1)+1
      else
      flow(j1,i1)=flow(j1,i1)-1
      endif
      endif
      enddo
      endif
      endif
      do i1=rhop1,n
      if(i1.ne.ii)then
      if(i1.lt.ii)then
      jj=eg(i1,ii)
      else
      jj=eg(ii,i1)
      endif
      do i2=jj,jj+g(ii,i1)-1
      if(intree(i2).eq.0)then
      if((head(j1).eq.head(i2)).or.(tail(j1).eq.tail(i2)))then
      do i3=1,rhop2
      flow(j1,i3)=flow(j1,i3)-flow(i2,i3)
      enddo
      else
      do i3=1,rhop2
      flow(j1,i3)=flow(j1,i3)+flow(i2,i3)
      enddo
      endif
      endif
      enddo
      endif
      enddo
      intree(j1)=0
      bb(tail(j1))=bb(tail(j1))-1
      bb(head(j1))=bb(head(j1))-1
      ii=rhop1
      goto 73
   83 continue
      do i1=rhop1,nli
      ii=0
      jj=0
      do i2=1,rho(1)
      kk=flow(i1,i2)
      if(kk.gt.0)then
      ii=ii+1
      elseif(kk.lt.0)then
      jj=jj+1
      endif
      enddo
      if((ii.gt.0).and.(jj.gt.0))then
      mlin(1:srec)='gen21_17'
      call mput(1,0,0,0)
      endif
      if(2*(ii+jj).gt.rho(1))then
      if(jj.eq.0)then
      do i2=1,rho(1)
      flow(i1,i2)=flow(i1,i2)-1
      enddo
      else
      do i2=1,rho(1)
      flow(i1,i2)=flow(i1,i2)+1
      enddo
      endif
      endif
      enddo
      do i1=1,nli
      do i2=1,rhop2
      if(abs(flow(i1,i2)).gt.1)then
      mlin(1:srec)='gen21_18'
      call mput(1,0,0,0)
      endif
      enddo
      enddo
      do i1=-2,2
      net(i1)=0
      enddo
      do i1=rhop1,nli
      do i2=rhop1,rhop2
      if(flow(i1,i2).ne.0)then
      if(head(i1).ne.tail(i1))then
      kk=-1
      else
      kk=-2
      endif
      flow(i1,0)=kk
      net(kk)=net(kk)+1
      goto 60
      endif
      enddo
      do i2=1,rho(1)
      if(flow(i1,i2).ne.0)then
      if(dflag(2).gt.0)then
      goto 05
      endif
      flow(i1,0)=1
      net(1)=net(1)+1
      goto 60
      endif
      enddo
      net(2)=net(2)+1
      flow(i1,0)=2
   60 continue
      enddo
      net(3)=net(1)+net(2)
      net(-3)=net(-1)+net(-2)
      net(0)=net(-3)+net(3)
      if(dflag(2).lt.0)then
      if(net(1).eq.0)then
      goto 05
      endif
      endif
      if(zbri(net(3)).eq.0)then
      goto 05
      elseif(zcho(nli-rho(1)-net(3)).eq.0)then
      goto 05
      elseif(rbri(net(1)).eq.0)then
      goto 05
      elseif(sbri(net(2)).eq.0)then
      goto 05
      endif
      if(dflag(7).ne.0)then
      call gsig(ii,head,tail)
      if(dflag(7).ne.ii)then
      goto 05
      endif
      endif
      do i1=1,rho(1)
      eg(i1,vmap(i1,1))=i1
      enddo
      do i1=rhop1,n
      ii=0
      do i2=1,rho(1)
      jj=0
      do i3=1,n
      if(i3.ne.i1)then
      j2=g(i1,i3)
      if(j2.gt.0)then
      if(i1.lt.i3)then
      j1=eg(i1,i3)
      else
      j1=eg(i3,i1)
      endif
      j3=head(j1)-i1
      do i4=j1,j1+j2-1
      if(j3.eq.0)then
      jj=jj+flow(i4,i2)
      else
      jj=jj-flow(i4,i2)
      endif
      enddo
      endif
      endif
      enddo
      if(jj.ne.ii)then
      if(i2.eq.1)then
      ii=jj
      else
      mlin(1:srec)='gen21_19'
      call mput(1,0,0,0)
      endif
      endif
      enddo
      enddo
      do i1=rhop1,n
      do i2=rhop1,rhop2
      ii=0
      do i3=1,n
      if(i3.ne.i1)then
      j2=g(i1,i3)
      if(j2.gt.0)then
      if(i1.lt.i3)then
      j1=eg(i1,i3)
      else
      j1=eg(i3,i1)
      endif
      j3=head(j1)-i1
      do i4=j1,j1+j2-1
      if(j3.eq.0)then
      ii=ii+flow(i4,i2)
      else
      ii=ii-flow(i4,i2)
      endif
      enddo
      endif
      endif
      enddo
      if(ii.ne.0)then
      mlin(1:srec)='gen21_20'
      call mput(1,0,0,0)
      endif
      enddo
      enddo
      goto 90
   97 continue
      mlin(1:srec)='gen21_1'
      call mput(1,0,0,0)
   80 continue
      cntr21=0
   90 return
      end
      subroutine gsig(situ,head,tail)
      implicit integer(a-z)
      save
      parameter ( maxleg=10, maxrho=7, maxi=10, maxdeg=6 )
      parameter ( maxn=maxi+maxi-2, maxli=maxn+maxrho )
      parameter ( srec=81, ssrec=62 )
      common/z1g/g(1:maxn,1:maxn),rho(1:maxdeg),nrho,rhop1,rhop2
      common/z4g/n,nli
      common/z8g/vdeg(1:maxn),xn(1:maxn)
      common/z18g/eg(1:maxn,1:maxn),flow(1:maxli,0:maxleg+maxrho),
     :net(-3:3)
      character*(srec) mlin
      common/z22g/mlin
      integer head(1:maxli),tail(1:maxli)
      integer aa(1:maxli),bb(1:maxli)
      situ=1
      do i1=rhop1,nli
      f1=flow(i1,0)
      if(abs(f1).eq.1)then
      do i2=1,i1-1
      ii=1
      if(i2.lt.rhop1)then
      if(f1.eq.1)then
      ii=0
      endif
      elseif(f1.eq.flow(i2,0))then
      ii=0
      endif
      if(ii.eq.0)then
      do i3=-1,1,2
      if(f1.lt.0)then
      do i4=rhop1,rhop2
      ii=flow(i1,i4)+i3*flow(i2,i4)
      if(ii.ne.0)then
      goto 10
      endif
      enddo
      endif
      ii=0
      jj=0
      do i4=1,rho(1)
      if(flow(i1,i4)+i3*flow(i2,i4).ne.0)then
      if(jj.ne.0)then
      goto 10
      endif
      ii=ii+1
      else
      if(ii.ne.0)then
      goto 10
      endif
      jj=jj+1
      endif
      enddo
      goto 80
  10  continue
      enddo
      endif
      enddo
      endif
      enddo
      if(rho(1)-1.ne.0)then
      if(net(2)-1.le.0)then
      goto 90
      endif
      else
      if(net(2).ne.0)then
      goto 80
      else
      goto 90
      endif
      endif
      rr=0
      do i3=rhop1,n
      if(xn(i3).eq.0)then
      aa(i3)=0
      else
      aa(i3)=1
      rr=rr+1
      bb(rr)=i3
      endif
      enddo
      if(rr.le.0)then
      mlin(1:srec)='gsig_1'
      call mput(1,0,0,0)
      endif
      ss=1
   30 continue
      tt=bb(ss)
      do i3=rhop1,n
      if(aa(i3).eq.0)then
      if(i3.ne.tt)then
      if(i3.lt.tt)then
      j3=i3
      j4=tt
      else
      j3=tt
      j4=i3
      endif
      if(g(j3,j4).ne.0)then
      if(flow(eg(j3,j4),0).ne.2)then
      aa(i3)=1
      rr=rr+1
      bb(rr)=i3
      endif
      endif
      endif
      endif
      enddo
      if(ss.lt.rr)then
      ss=ss+1
      goto 30
      endif
      do i1=rhop1,nli
      if(flow(i1,0).eq.2)then
      if(aa(head(i1)).eq.0)then
      if(aa(tail(i1)).eq.0)then
      goto 80
      endif
      endif
      endif
      enddo
      goto 90
   80 continue
      situ=-1
   90 return
      end
      subroutine init0
      implicit integer(a-z)
      save
      parameter ( sibuff=1048574, scbuff=131064, swbuff=65528 )
      parameter ( sxbuff=2040 )
      parameter ( srec=81, ssrec=62 )
      common/z12g/cflag(1:7),dflag(1:20),jflag(1:10),mflag(1:6)
      common/z15g/pten(1:10),iref,wiref,wsint
      character*(sxbuff) stxb
      common/z27g/stxb
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z35g/stibs(1:1),stcbs(1:1),stwbs(1:1)
      common/z48g/acomma,ascol,albra,arbra,alpar,arpar
      common/z51g/aeq,agt,alt,aminus,anine,aplus,aslash,atimes,azero
      common/z52g/acf(0:127),abo(0:6),acr,alf,aspace,dquote,squote
      common/z56g/acz(0:127)
      integer lsize(1:3),usize(1:3)
      wiref=4
      iref=1
      do i1=1,wiref
      iref=10*iref
      enddo
      iref=iref-1
      do i1=1,7
      cflag(i1)=0
      enddo
      do i1=1,20
      dflag(i1)=0
      enddo
      do i1=1,10
      jflag(i1)=0
      enddo
      do i1=1,6
      mflag(i1)=0
      enddo
      abo(1)=aspace
      abo(2)=126
      aeq=ichar('=')
      aminus=ichar('-')
      anine=ichar('9')
      aplus=ichar('+')
      atimes=ichar('*')
      azero=ichar('0')
      ii=max(aplus,aminus,azero,anine,aeq)
      jj=min(aplus,aminus,azero,anine,aeq)
      aspace=ichar(' ')
      dquote=ichar('"')
      squote=ichar("'")
      ii=max(ii,aspace,dquote,squote)
      jj=min(jj,aspace,dquote,squote)
      acomma=ichar(',')
      ascol=ichar(';')
      albra=ichar('[')
      arbra=ichar(']')
      alpar=ichar('(')
      arpar=ichar(')')
      ii=max(ii,acomma,ascol,albra,arbra,alpar,arpar)
      jj=min(jj,acomma,ascol,albra,arbra,alpar,arpar)
      alt=ichar('<')
      agt=ichar('>')
      aslash=ichar('/')
      ausco=ichar('_')
      ii=max(ii,alt,agt,aslash,ausco)
      jj=min(jj,alt,agt,aslash,ausco)
      abo(3)=65
      abo(4)=90
      abo(5)=97
      abo(6)=122
      ii=max(ii,abo(3),abo(4),abo(5),abo(6))
      jj=min(jj,abo(3),abo(4),abo(5),abo(6))
      if((jj.lt.abo(1)).or.(ii.gt.abo(2)))then
      call uput(2)
      endif
      ii=0
      if(anine-azero.ne.9)then
      ii=1
      elseif(abo(4)-abo(3).ne.25)then
      ii=1
      elseif(abo(6)-abo(5).ne.25)then
      ii=1
      elseif((abo(4).ge.abo(5)).and.(abo(6).ge.abo(3)))then
      ii=1
      endif
      if(ii.ne.0)then
      call uput(2)
      endif
      abo(0)=abo(3)-abo(5)
      alf=10
      acr=13
      do i1=0,127
      acz(i1)=-1
      acf(i1)=-2
      enddo
      acf(ichar('%'))=-1
      acf(ichar('*'))=-1
      acf(ichar('#'))=-1
      acf(ausco)=0
      do i1=azero,anine
      acz(i1)=1
      acf(i1)=1
      enddo
      acz(aminus)=2
      acz(aplus)=2
      acz(atimes)=3
      acz(alpar)=6
      acz(arpar)=6
      acz(acomma)=7
      acz(aspace)=0
      do i1=abo(3),abo(4)
      acf(i1)=2
      enddo
      do i1=abo(5),abo(6)
      acz(i1)=5
      acf(i1)=3
      enddo
      lsize(1)=2040
      lsize(2)=8184
      lsize(3)=8190
      usize(1)=32760
      usize(2)=134217720
      usize(3)=134217726
      if((sxbuff.lt.lsize(1)).or.(sxbuff.gt.usize(1)))then
      call uput(7)
      endif
      if((scbuff.lt.lsize(2)).or.(scbuff.gt.usize(2)))then
      call uput(8)
      endif
      if((sibuff.lt.lsize(3)).or.(sibuff.gt.usize(3)))then
      call uput(9)
      endif
      if((srec.lt.81).or.(srec.gt.lsize(1)))then
      call uput(10)
      endif
      if((ssrec.ge.min(srec,128)).or.(ssrec.le.48))then
      call uput(11)
      endif
      stcbs(1)=0
      stibs(1)=0
      call init1
      return
      end
      subroutine init1
      implicit integer(a-z)
      save
      parameter ( maxleg=10, maxrho=7, maxi=10, maxdeg=6 )
      parameter ( srec=81, ssrec=62 )
      parameter ( sibuff=1048574, scbuff=131064, swbuff=65528 )
      parameter ( sxbuff=2040 )
      parameter ( nfiles=5 )
      parameter ( eoa=-2047, nap=-8191 )
      parameter ( maxtak=2, maxpot=6 )
      common/z2in/momep(0:maxleg),momel(0:maxleg),kpqs(1:4)
      common/z7in/aunit(1:nfiles)
      common/z5g/psym(0:0),psyms,nsym
      common/z13g/kla(0:0),bbc(0:0)
      common/z15g/pten(1:10),iref,wiref,wsint
      common/z21g/punct1(0:0)
      character*(srec) mlin
      common/z22g/mlin
      common/z23g/tak(1:maxtak),pot(0:maxpot),ks
      common/z26g/kc(0:0),kse(0:0),pske(0:0),wske(0:0),klo(0:0),nske
      character*(sxbuff) stxb
      common/z27g/stxb
      common/z28g/wera(0:0),werb(0:0),nwer
      common/z29g/pkey(0:0),wkey(0:0),ikey(0:0),pokey(0:0),wokey(0:0),
     :cokey(0:0),popt1(0:0),wopt1(0:0),fopt1(0:0),vopt1(0:0),
     :popt5(0:0),wopt5(0:0),copt5(0:0),popt0(0:0),wopt0(0:0),
     :copt0(0:0),vopt0(0:0)
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z32g/aaf(0:4)
      common/z35g/stibs(1:1),stcbs(1:1),stwbs(1:1)
      common/z44g/xtstrp(0:0),xtstrl(0:0)
      common/z45g/qdatp,qdatl,qvp,qvl,dprefp,dprefl
      common/z46g/popt3(0:0),wopt3(0:0),copt3(0:0),popt9(0:0),
     :wopt9(0:0),copt9(0:0)
      common/z47g/ndiagp,ndiagl,hhp,hhl,doffp,doffl,noffp,noffl
      common/z48g/acomma,ascol,albra,arbra,alpar,arpar
      common/z49g/popt7(0:0),wopt7(0:0),copt7(0:0)
      common/z50g/tfta(0:0),tftb(0:0),tftic(0:0),ntft
      common/z51g/aeq,agt,alt,aminus,anine,aplus,aslash,atimes,azero
      common/z52g/acf(0:127),abo(0:6),acr,alf,aspace,dquote,squote
      common/z55g/aopp(0:0),aopl(0:0),aopna(0:0),aopnb(0:0)
      integer weri(0:0),aopc(0:0)
      kpqs(1)=107
      kpqs(2)=112
      kpqs(3)=113
      kpqs(4)=115
      j1=2*srec
      do i1=1,j1
      stxb(i1:i1)=char(azero)
      enddo
      j1=j1+1
      stxb(j1:j1)=char(ascol)
      ii=0
      call spak(stxb(1:j1),ii,0,0,0)
      call spak('8171826570131914201418;',ii,1,0,0)
      qvp=stib(stibs(1)+1)
      qvl=stib(stibs(1)+2)
      wsint=15
      jj=wsint+1
      ndiagp=stcbs(1)
      call aocb(jj)
      ii=stcbs(1)
      stcb(ii:ii)=char(alf)
      call sdiag(3,-1)
      hhp=stcbs(1)
      call aocb(jj)
      ii=stcbs(1)
      stcb(ii:ii)=char(alf)
      noffp=stcbs(1)
      noffl=1
      call aocb(jj)
      ii=noffp+1
      stcb(ii:ii)=char(azero)
      ii=stcbs(1)
      stcb(ii:ii)=char(alf)
      doffp=stcbs(1)
      doffl=1
      call aocb(jj)
      ii=doffp+1
      stcb(ii:ii)=char(azero)
      ii=stcbs(1)
      stcb(ii:ii)=char(alf)
      psym(0)=nap
      do i1=1,nfiles
      aunit(i1)=0
      enddo
      nxts=0
      call spak('78797813808279806571658479828326000008;',nxts,0,0,0)
      call spak('808279806571658479828326000008;',nxts,0,0,0)
      call spak('4611;',nxts,0,0,0)
      call spak('3511;',nxts,0,0,0)
      call spak('4613;',nxts,0,0,0)
      call spak('3513;',nxts,0,0,0)
      call spak('866982847367698326000008;',nxts,0,0,0)
      call spak('861368697182696983;',nxts,0,0,0)
      call spak('036873657182657783;',nxts,0,0,0)
      call spak('03718265807283;',nxts,0,0,0)
      call spak('00000000000000847984657600290000;',nxts,0,0,0)
      call spak('00677978786967846968006873657182657783;',nxts,0,0,0)
      call spak('0067797878696784696800718265807283;',nxts,0,0,0)
      call spak('0900;',nxts,0,0,0)
      call spak('14141414;',nxts,0,0,0)
      call spak('0000001010;',nxts,0,0,0)
      call spak('00000008870109;',nxts,0,0,0)
      call spak('141414;',nxts,0,0,0)
      call spak('000000876582787378712600;',nxts,0,0,0)
      call spak('00000069828279822600;',nxts,0,0,0)
      call spak('7073766900;',nxts,0,0,0)
      call spak('77796869761370737669;',nxts,0,0,0)
      call spak('767366826582891370737669;',nxts,0,0,0)
      call spak('83848976691370737669;',nxts,0,0,0)
      call spak('7985848085841370737669;',nxts,0,0,0)
      call spak('12007673786900;',nxts,0,0,0)
      call spak('1200767378698300;',nxts,0,0,0)
      kk=2
      call aoib(kk)
      jj=stibs(1)
      do i1=1,kk
      stib(jj)=eoa
      jj=jj-1
      enddo
      ii=nxts+1
      call trm(kk,ii)
      xtstrl(0)=stibs(1)-ii
      xtstrp(0)=xtstrl(0)-ii
      call spak('817182657014686584;',ii,1,0,0)
      qdatp=stib(stibs(1)+1)
      qdatl=stib(stibs(1)+2)
      uc=1
      nos=0
      nkey=0
      call spak('798584808584,1;',nkey,0,uc,0)
      call spak('8384897669,2;',nkey,0,uc,0)
      call spak('7779686976,3;',nkey,0,uc,0)
      call spak('7378,4;',nkey,0,uc,0)
      call spak('798584,5;',nkey,0,uc,0)
      call spak('7679798083,6;',nkey,0,uc,0)
      call spak('76797980637779776978848577,7;',nkey,0,uc,0)
      call spak('79808473797883,8;',nkey,0,uc,0)
      kk=3
      call aoib(kk)
      jj=stibs(1)
      do i1=1,kk
      stib(jj)=eoa
      jj=jj-1
      enddo
      ii=nkey+1
      call trm(kk,ii)
      ikey(0)=stibs(1)-ii
      wkey(0)=ikey(0)-ii
      pkey(0)=wkey(0)-ii
      uc=1
      nos=0
      nokey=0
      call spak('677978707371,1;',nokey,0,uc,0)
      call spak('737868698863797070836984,2;',nokey,0,uc,0)
      kk=3
      call aoib(kk)
      jj=stibs(1)
      do i1=1,kk
      stib(jj)=eoa
      jj=jj-1
      enddo
      ii=nokey+1
      call trm(kk,ii)
      cokey(0)=stibs(1)-ii
      wokey(0)=cokey(0)-ii
      pokey(0)=wokey(0)-ii
      uc=1
      nos=0
      nopt0=0
      call spak('73787079,1,1;',nopt0,0,uc,0)
      call spak('86698266798369,1,2;',nopt0,0,uc,0)
      call spak('787973787079,1,-1;',nopt0,0,uc,0)
      call spak('787976738384,2,-1;',nopt0,0,uc,0)
      call spak('7670,3,-1;',nopt0,0,uc,0)
      call spak('7667797780,4,1;',nopt0,0,uc,0)
      call spak('8367797780,5,1;',nopt0,0,uc,0)
      kk=4
      call aoib(kk)
      jj=stibs(1)
      do i1=1,kk
      stib(jj)=eoa
      jj=jj-1
      enddo
      ii=nopt0+1
      call trm(kk,ii)
      vopt0(0)=stibs(1)-ii
      copt0(0)=vopt0(0)-ii
      wopt0(0)=copt0(0)-ii
      popt0(0)=wopt0(0)-ii
      uc=1
      nos=1
      nopt1=0
      call spak('787966827368716983,1,1;',nopt1,0,uc,nos)
      call spak('7978698073,1,1;',nopt1,0,uc,nos)
      call spak('66827368716983,1,-1;',nopt1,0,uc,nos)
      call spak('7978698082,1,-1;',nopt1,0,uc,nos)
      call spak('78798266827368716983,2,1;',nopt1,0,uc,nos)
      call spak('8266827368716983,2,-1;',nopt1,0,uc,nos)
      call spak('78798366827368716983,3,1;',nopt1,0,uc,nos)
      call spak('78798465688079766983,3,1;',nopt1,0,uc,nos)
      call spak('8366827368716983,3,-1;',nopt1,0,uc,nos)
      call spak('8465688079766983,3,-1;',nopt1,0,uc,nos)
      call spak('79788372697676,4,1;',nopt1,0,uc,nos)
      call spak('7970708372697676,4,-1;',nopt1,0,uc,nos)
      call spak('7978837269767688,5,1;',nopt1,0,uc,nos)
      call spak('797070837269767688,5,-1;',nopt1,0,uc,nos)
      call spak('7879837371776583,7,1;',nopt1,0,uc,nos)
      call spak('837371776583,7,-1;',nopt1,0,uc,nos)
      call spak('7879837865737683,8,1;',nopt1,0,uc,nos)
      call spak('837865737683,8,-1;',nopt1,0,uc,nos)
      call spak('6789677673,9,1;',nopt1,0,uc,nos)
      call spak('6789677682,9,-1;',nopt1,0,uc,nos)
      call spak('7978698673,10,1;',nopt1,0,uc,nos)
      call spak('7978698682,10,-1;',nopt1,0,uc,nos)
      call spak('667380658284,11,1;',nopt1,0,uc,nos)
      call spak('787978667380658284,11,-1;',nopt1,0,uc,nos)
      call spak('837377807669,12,1;',nopt1,0,uc,nos)
      call spak('787984837377807669,12,-1;',nopt1,0,uc,nos)
      call spak('8479807976,13,1;',nopt1,0,uc,nos)
      call spak('707679798083,14,1;',nopt1,0,uc,nos)
      call spak('787984707679798083,14,-1;',nopt1,0,uc,nos)
      call spak('7879798084,20,1;',nopt1,0,uc,nos)
      kk=4
      call aoib(kk)
      jj=stibs(1)
      do i1=1,kk
      stib(jj)=eoa
      jj=jj-1
      enddo
      ii=nopt1+1
      call trm(kk,ii)
      vopt1(0)=stibs(1)-ii
      fopt1(0)=vopt1(0)-ii
      wopt1(0)=fopt1(0)-ii
      popt1(0)=wopt1(0)-ii
      uc=1
      nos=0
      nopt3=0
      call spak('84828569,1;',nopt3,0,uc,0)
      call spak('7065768369,-1;',nopt3,0,uc,0)
      kk=3
      call aoib(kk)
      jj=stibs(1)
      do i1=1,kk
      stib(jj)=eoa
      jj=jj-1
      enddo
      ii=nopt3+1
      call trm(kk,ii)
      copt3(0)=stibs(1)-ii
      wopt3(0)=copt3(0)-ii
      popt3(0)=wopt3(0)-ii
      uc=1
      nos=0
      nopt5=0
      call spak('7380827980,1;',nopt5,0,uc,nos)
      call spak('668273687169,2;',nopt5,0,uc,nos)
      call spak('6772798268,3;',nopt5,0,uc,nos)
      call spak('82668273687169,4;',nopt5,0,uc,nos)
      call spak('83668273687169,5;',nopt5,0,uc,nos)
      call spak('86838577,6;',nopt5,0,uc,nos)
      call spak('80838577,7;',nopt5,0,uc,nos)
      call spak('6976737875,11;',nopt5,0,uc,nos)
      call spak('8076737875,12;',nopt5,0,uc,nos)
      kk=3
      call aoib(kk)
      jj=stibs(1)
      do i1=1,kk
      stib(jj)=eoa
      jj=jj-1
      enddo
      ii=nopt5+1
      call trm(kk,ii)
      copt5(0)=stibs(1)-ii
      wopt5(0)=copt5(0)-ii
      popt5(0)=wopt5(0)-ii
      kk=0
      do i1=1,nopt5
      j1=stib(copt5(0)+i1)
      if(j1.ne.-1)then
      if(j1.gt.kk)then
      ntft=ntft+1
      kk=j1
      elseif(j1.lt.kk)then
      mlin(1:srec)='init1_1'
      call mput(1,0,0,0)
      endif
      endif
      enddo
      ii=ntft+1
      jj=2*ii
      call aoib(jj)
      tftb(0)=stibs(1)-ii
      tfta(0)=tftb(0)-ii
      stib(tftb(0))=eoa
      stib(stibs(1))=eoa
      jj=kk+1
      tftic(0)=stibs(1)
      call aoib(jj)
      stib(stibs(1))=eoa
      do i1=tftic(0)+1,tftic(0)+kk
      stib(tftic(0)+i1)=0
      enddo
      ii=0
      do i1=1,nopt5
      j1=tftic(0)+stib(copt5(0)+i1)
      if(stib(j1).eq.0)then
      ii=ii+1
      stib(j1)=ii
      endif
      enddo
      uc=1
      nos=0
      nopt7=0
      call spak('69886776,1;',nopt7,0,uc,0)
      call spak('73786776,0;',nopt7,0,uc,0)
      kk=3
      call aoib(kk)
      jj=stibs(1)
      do i1=1,kk
      stib(jj)=eoa
      jj=jj-1
      enddo
      ii=nopt7+1
      call trm(kk,ii)
      copt7(0)=stibs(1)-ii
      wopt7(0)=copt7(0)-ii
      popt7(0)=wopt7(0)-ii
      pot(0)=1
      do i1=1,maxpot
      pot(i1)=2*pot(i1-1)
      enddo
      call spak('6885657613;',ii,1,0,0)
      dprefp=stib(stibs(1)+1)
      dprefl=stib(stibs(1)+2)
      uc=0
      nos=0
      nske=0
      ii=stibs(1)
      call spak('8082797679718569,1,0,0,0,0;',nske,0,uc,0)
      call spak('68736571826577,2,0,0,0,0;',nske,0,uc,0)
      call spak('6980737679718569,3,0,0,0,0;',nske,0,uc,0)
      call spak('69887384,4,0,0,0,0;',nske,0,uc,0)
      call spak('677977776578686376797980,11,1,5,0,0;',nske,0,uc,0)
      call spak('6779777765786863767378696376797980,12,1,5,1,0;',nske,
     :0,uc,0)
      call spak('73786376797980,13,1,2,0,0;',nske,0,uc,0)
      call spak('736376797980,13,1,2,0,0;',nske,0,uc,0)
      call spak('7985846376797980,14,1,2,0,0;',nske,0,uc,0)
      call spak('796376797980,14,1,2,0,0;',nske,0,uc,0)
      call spak('808279806571658479826376797980,15,1,2,0,0;',nske,0,uc
     :,0)
      call spak('806376797980,15,1,2,0,0;',nske,0,uc,0)
      call spak('8669828469886376797980,16,1,2,0,0;',nske,0,uc,0)
      call spak('866376797980,16,1,2,0,0;',nske,0,uc,0)
      call spak('8265896376797980,17,1,2,8,0;',nske,0,uc,0)
      call spak('826376797980,17,1,2,8,0;',nske,0,uc,0)
      call spak('697868,19,1,7,31,0;',nske,0,uc,0)
      call spak('66656775,22,2,7,0,0;',nske,0,uc,0)
      call spak('78698776737869,29,5,7,0,0;',nske,0,uc,0)
      call spak('7073697668,31,4,2,23,0;',nske,0,uc,0)
      call spak('7779776978848577,32,4,2,23,0;',nske,0,uc,0)
      call spak('68856576137073697668,33,4,2,23,0;',nske,0,uc,0)
      call spak('68856576137779776978848577,34,4,2,23,0;',nske,0,uc,0)
      call spak('70736976686383737178,35,4,2,23,0;',nske,0,uc,0)
      call spak('70736976686384898069,40,4,2,23,0;',nske,0,uc,0)
      call spak('80827980657165847982637378686988,41,4,2,20,0;',nske,0
     :,uc,0)
      call spak('80637378686988,41,4,2,20,0;',nske,0,uc,0)
      call spak('7073697668637378686988,42,4,2,23,0;',nske,0,uc,0)
      call spak('70637378686988,42,4,2,23,0;',nske,0,uc,0)
      call spak('826589637378686988,43,4,2,23,0;',nske,0,uc,0)
      call spak('82637378686988,43,4,2,23,0;',nske,0,uc,0)
      call spak('866982846988637378686988,44,4,2,31,0;',nske,0,uc,0)
      call spak('86637378686988,44,4,2,31,0;',nske,0,uc,0)
      call spak('68856576137073697668637378686988,45,4,2,20,0;',nske,0
     :,uc,0)
      call spak('688565761370637378686988,45,4,2,20,0;',nske,0,uc,0)
      call spak('6885657613826589637378686988,46,4,2,20,0;',nske,0,uc,0)
      call spak('688565761382637378686988,46,4,2,20,0;',nske,0,uc,0)
      call spak('6885657613866982846988637378686988,47,4,2,20,0;',nske
     :,0,uc,0)
      call spak('688565761386637378686988,47,4,2,20,0;',nske,0,uc,0)
      call spak('86698284698863686971826969,48,4,2,31,0;',nske,0,uc,0)
      call spak('8663686971826969,48,4,2,31,0;',nske,0,uc,0)
      call spak('688565761386698284698863686971826969,49,4,2,20,0;',
     :nske,0,uc,0)
      call spak('68856576138663686971826969,49,4,2,20,0;',nske,0,uc,0)
      call spak('766971637378686988,51,4,2,2,0;',nske,0,uc,0)
      call spak('76637378686988,51,4,2,2,0;',nske,0,uc,0)
      call spak('7378637378686988,52,4,2,1,0;',nske,0,uc,0)
      call spak('73637378686988,52,4,2,1,0;',nske,0,uc,0)
      call spak('798584637378686988,53,4,2,2,0;',nske,0,uc,0)
      call spak('79637378686988,53,4,2,2,0;',nske,0,uc,0)
      call spak('83737178,61,3,2,0,0;',nske,0,uc,0)
      call spak('7773788583,62,3,2,0,0;',nske,0,uc,0)
      call spak('68736571826577637378686988,71,3,6,0,0;',nske,0,uc,0)
      call spak('838977776984828963706567847982,72,3,2,0,0;',nske,0,uc
     :,0)
      call spak('838977776984828963788577666982,73,3,2,0,0;',nske,0,uc
     :,0)
      call spak('8082798065716584798283,74,3,2,0,0;',nske,0,uc,0)
      call spak('76697183,75,3,2,0,0;',nske,0,uc,0)
      call spak('7679798083,76,3,2,0,0;',nske,0,uc,0)
      call spak('8669828473676983,77,3,2,0,0;',nske,0,uc,0)
      call spak('76697183637378,78,3,2,0,0;',nske,0,uc,0)
      call spak('7669718363798584,79,3,2,0,0;',nske,0,uc,0)
      call spak('80827971826577,81,3,5,0,0;',nske,0,uc,0)
      call spak('677977776578686368658465,82,4,5,3,0;',nske,0,uc,0)
      kk=7
      jj=stibs(1)
      if(jj-ii.ne.kk*nske)then
      mlin(1:srec)='init1_2'
      call mput(1,0,0,0)
      endif
      call aoib(kk)
      do i1=jj+1,stibs(1)
      stib(i1)=eoa
      enddo
      ii=nske+1
      call trm(kk,ii)
      bbc(0)=stibs(1)-ii
      klo(0)=bbc(0)-ii
      kse(0)=klo(0)-ii
      kla(0)=kse(0)-ii
      kc(0)=kla(0)-ii
      wske(0)=kc(0)-ii
      pske(0)=wske(0)-ii
      uc=1
      nos=1
      nopt9=0
      call spak('6988846982786576,5;',nopt9,0,uc,nos)
      call spak('78798465688079766983,1;',nopt9,0,uc,nos)
      kk=3
      call aoib(kk)
      jj=stibs(1)
      do i1=1,kk
      stib(jj)=eoa
      jj=jj-1
      enddo
      ii=nopt9+1
      call trm(kk,ii)
      copt9(0)=stibs(1)-ii
      wopt9(0)=copt9(0)-ii
      popt9(0)=wopt9(0)-ii
      aaf(0)=5
      aaf(1)=6
      aaf(2)=6
      aaf(3)=5
      aaf(4)=5
      punct1(0)=stibs(1)+1
      ii=128
      call aoib(ii)
      do i1=0,127
      stib(punct1(0)+i1)=0
      enddo
      stib(punct1(0)+squote)=1
      stib(punct1(0)+acomma)=1
      stib(punct1(0)+ascol)=1
      stib(punct1(0)+aeq)=1
      stib(punct1(0)+albra)=1
      stib(punct1(0)+arbra)=1
      stib(punct1(0)+alpar)=1
      stib(punct1(0)+arpar)=1
      uc=0
      nos=0
      nwer=0
      call spak('7779686976006765780066690083807673840073788479006873'
     ://'8374797378840083856613777968697683,1;',nwer,0,0,0)
      call spak('6584007669658384007978690070736976680073830078798400'
     ://'806582840079700065788900866982846988,2;',nwer,0,0,0)
      call spak('7879008669828469880077658900767378751385800087738472'
     ://'00737867797773787100707369766800,3;',nwer,0,0,0)
      call spak('7879008669828469880077658900767378751385800087738472'
     ://'00798584717973787100707369766800,4;',nwer,0,0,0)
      call spak('6988846982786576007073697668830067657878798400666900'
     ://'677978786967846968,5;',nwer,0,0,0)
      call spak('7968680078857766698200797000698884698278657600657884'
     ://'7367797777858473787100707369766883,6;',nwer,0,0,0)
      call spak('7779686976137073766900726583006885807673676584690086'
     ://'69828473676983,7;',nwer,0,0,0)
      call spak('6977808489138384827378710067797883846578840068698469'
     ://'678469682600,8;',nwer,0,0,0)
      call spak('6977808489138384827378710070857867847379780068698469'
     ://'678469682600,9;',nwer,0,0,0)
      call spak('788576760070857867847379780068698469678469682600,10;'
     :,nwer,0,0,0)
      call spak('6772696775006779787073710079808473797883,11;',nwer,0,
     :0,0)
      call spak('67726967750079808473797883,12;',nwer,0,0,0)
      call spak('6779788482656873678479828900697673787500838465846977'
     ://'697884,13;',nwer,0,0,0)
      call spak('8482738673657600697673787500838465846977697884,14;',
     :nwer,0,0,0)
      call spak('1076737875100083846584697769788483008269818573826900'
     ://'6873657182657783008773847200658400766965838400180076697183,'
     ://'15;',nwer,0,0,0)
      call spak('7665827169007885776669820070798578680079820071697869'
     ://'8265846968,16;',nwer,0,0,0)
      call spak('7376761370798277696800658273847277698473676576006988'
     ://'8082698383737978,17;',nwer,0,0,0)
      kk=3
      call aoib(kk)
      jj=stibs(1)
      do i1=1,kk
      stib(jj)=eoa
      jj=jj-1
      enddo
      ii=nwer+1
      call trm(kk,ii)
      weri(0)=stibs(1)-ii
      werb(0)=weri(0)-ii
      wera(0)=werb(0)-ii
      do i1=1,nwer
      stib(werb(0)+i1)=stib(werb(0)+i1)-1+stib(wera(0)+i1)
      if(stib(weri(0)+i1).ne.i1)then
      mlin(1:srec)='init1_3'
      call mput(1,0,0,0)
      endif
      enddo
      call aoib(-ii)
      weri(0)=nap
      uc=0
      nos=0
      naop=0
      call spak('656683,1,1,1;',naop,0,0,0)
      call spak('8472698465,1,0,2;',naop,0,0,0)
      call spak('6869768465,2,2,3;',naop,0,0,0)
      call spak('777968,2,2,4;',naop,0,0,0)
      call spak('777378,1,0,5;',naop,0,0,0)
      call spak('776588,1,0,6;',naop,0,0,0)
      call spak('8079876982,2,2,7;',naop,0,0,0)
      kk=5
      call aoib(kk)
      jj=stibs(1)
      do i1=1,kk
      stib(jj)=eoa
      jj=jj-1
      enddo
      ii=naop+1
      call trm(kk,ii)
      aopc(0)=stibs(1)-ii
      aopnb(0)=aopc(0)-ii
      aopna(0)=aopnb(0)-ii
      aopl(0)=aopna(0)-ii
      aopp(0)=aopl(0)-ii
      do i1=1,naop
      if(stib(aopc(0)+i1).ne.i1)then
      mlin(1:srec)='init1_4'
      call mput(1,0,0,0)
      endif
      enddo
      call aoib(-ii)
      aopc(0)=nap
      wiref=5
      iref=46340
      pten(1)=10
      do i1=2,9
      pten(i1)=10*pten(i1-1)
      enddo
      return
      end
      subroutine aocb(delta)
      implicit integer(a-z)
      save
      parameter ( sibuff=1048574, scbuff=131064, swbuff=65528 )
      parameter ( srec=81, ssrec=62 )
      common/z35g/stibs(1:1),stcbs(1:1),stwbs(1:1)
      character*(srec) mlin
      common/z22g/mlin
      integer delta
      if(delta.gt.0)then
      if((delta.gt.scbuff).or.(stcbs(1).gt.scbuff-delta))then
      call uput(5)
      endif
      elseif(delta.lt.0)then
      if(delta.lt.-stcbs(1))then
      mlin(1:srec)='aocb_1'
      call mput(1,0,0,0)
      endif
      endif
      stcbs(1)=stcbs(1)+delta
      return
      end
      subroutine vaocb(delta)
      implicit integer(a-z)
      save
      parameter ( sibuff=1048574, scbuff=131064, swbuff=65528 )
      parameter ( srec=81, ssrec=62 )
      common/z35g/stibs(1:1),stcbs(1:1),stwbs(1:1)
      character*(srec) mlin
      common/z22g/mlin
      integer delta
      if((delta.lt.0).or.(stcbs(1).lt.0))then
      mlin(1:srec)='vaocb_1'
      call mput(1,0,0,0)
      elseif(stcbs(1).gt.scbuff-delta)then
      call uput(5)
      endif
      return
      end
      subroutine aoib(delta)
      implicit integer(a-z)
      save
      parameter ( sibuff=1048574, scbuff=131064, swbuff=65528 )
      parameter ( srec=81, ssrec=62 )
      common/z35g/stibs(1:1),stcbs(1:1),stwbs(1:1)
      character*(srec) mlin
      common/z22g/mlin
      integer delta
      if(delta.gt.0)then
      if((delta.gt.sibuff).or.(stibs(1).gt.sibuff-delta))then
      call uput(4)
      endif
      elseif(delta.lt.0)then
      if(delta.lt.-stibs(1))then
      mlin(1:srec)='aoib_1'
      call mput(1,0,0,0)
      endif
      endif
      stibs(1)=stibs(1)+delta
      return
      end
      subroutine vaoib(delta)
      implicit integer(a-z)
      save
      parameter ( sibuff=1048574, scbuff=131064, swbuff=65528 )
      parameter ( srec=81, ssrec=62 )
      common/z35g/stibs(1:1),stcbs(1:1),stwbs(1:1)
      character*(srec) mlin
      common/z22g/mlin
      integer delta
      if((delta.lt.0).or.(stibs(1).lt.0))then
      mlin(1:srec)='vaoib_1'
      call mput(1,0,0,0)
      elseif(stibs(1).gt.sibuff-delta)then
      call uput(4)
      endif
      return
      end
      subroutine mput(istop,nl1,nl2,nf1)
      implicit integer(a-z)
      save
      parameter ( srec=81, ssrec=62 )
      parameter ( sibuff=1048574, scbuff=131064, swbuff=65528 )
      parameter ( sxbuff=2040 )
      common/z12g/cflag(1:7),dflag(1:20),jflag(1:10),mflag(1:6)
      common/z15g/pten(1:10),iref,wiref,wsint
      character*(srec) mlin
      common/z22g/mlin
      character*(sxbuff) stxb
      common/z27g/stxb
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z35g/stibs(1:1),stcbs(1:1),stwbs(1:1)
      common/z44g/xtstrp(0:0),xtstrl(0:0)
      common/z45g/qdatp,qdatl,qvp,qvl,dprefp,dprefl
      common/z51g/aeq,agt,alt,aminus,anine,aplus,aslash,atimes,azero
      common/z52g/acf(0:127),abo(0:6),acr,alf,aspace,dquote,squote
      if(istop.ne.0)then
      if(jflag(6).eq.1)then
      call spp(0)
      jflag(6)=2
      endif
      endif
      if(jflag(7).eq.0)then
      write(unit=*,fmt='(a)')
      endif
      call hrul(2)
      write(unit=*,fmt='(a)')stxb(1:ssrec-1)
      slin1=ssrec-1
      slin2=2*srec
      slin3=(scbuff-wiref)-stcbs(1)
      stab=0
      do j2=srec-1,1,-1
      if(ichar(mlin(j2:j2)).ne.aspace)then
      goto 30
      endif
      enddo
      j1=0
      j2=0
      goto 50
   30 continue
      do j1=1,j2
      if(ichar(mlin(j1:j1)).ne.aspace)then
      goto 50
      endif
      enddo
   50 continue
      if(jflag(6).eq.0)then
      write(unit=*,fmt='(a)')'   error: '//mlin(j1:j2)
      goto 70
      endif
      if(istop.eq.0)then
      j3=19
      else
      j3=20
      endif
      i1=stib(xtstrp(0)+j3)
      i2=stib(xtstrl(0)+j3)
      stab=i2-1
      if(i2.ge.slin2)then
      write(unit=*,fmt='(/,a,/)')' error: mput_1'
      stop
      endif
      stxb(1:i2)=stcb(i1:i1-1+i2)
      if(nf1.ge.0)then
      if(j1.gt.0)then
      j3=min(slin2-i2,j2-j1+1)
      if(j3.gt.0)then
      stxb(i2+1:i2+j3)=mlin(j1:j1-1+j3)
      i2=i2+j3
      endif
      endif
      endif
      if(nf1.ne.0)then
      j4=abs(nf1)
      if(nf1.gt.0)then
      if(i2.lt.slin2)then
      i2=i2+1
      stxb(i2:i2)=stcb(1:1)
      endif
      endif
      j3=20+j4
      i1=stib(xtstrl(0)+j3)
      if(i2+i1.le.slin2)then
      j3=stib(xtstrp(0)+j3)
      stxb(i2+1:i2+i1)=stcb(j3:j3-1+i1)
      endif
      i2=i2+i1
      if(j4.eq.1)then
      if(i2+qdatl.le.slin2)then
      stxb(i2+1:i2+qdatl)=stcb(qdatp:qdatp-1+qdatl)
      endif
      i2=i2+qdatl
      endif
      endif
      if(nf1.lt.0)then
      if(j1.gt.0)then
      if(j2-j1+2.le.slin2-i2)then
      stxb(i2+1:i2+2+(j2-j1))=stcb(1:1)//mlin(j1:j2)
      i2=i2+(j2-j1)+2
      endif
      endif
      endif
      if(nl1.gt.0)then
      if(nl1.eq.nl2)then
      i1=stib(xtstrl(0)+26)
      if(i2+i1.le.slin2)then
      j3=stib(xtstrp(0)+26)
      stxb(i2+1:i2+i1)=stcb(j3:j3-1+i1)
      endif
      i2=i2+i1
      if(i2+wztos(nl1).le.slin2)then
      if(slin3.ge.0)then
      call dkar(nl1,jj)
      j1=stcbs(1)
      call ctxb(i2,jj,j1)
      else
      i1=stib(xtstrl(0)+18)
      j3=stib(xtstrp(0)+18)
      stxb(i2+1:i2+i1)=stcb(j3:j3-1+i1)
      i2=i2+i1
      endif
      endif
      else
      i1=stib(xtstrl(0)+27)
      if(i2+i1.le.slin2)then
      j3=stib(xtstrp(0)+27)
      stxb(i2+1:i2+i1)=stcb(j3:j3-1+i1)
      endif
      i2=i2+i1
      if(i2+wztos(nl1).le.slin2)then
      if(slin3.ge.0)then
      call dkar(nl1,jj)
      j1=stcbs(1)
      call ctxb(i2,jj,j1)
      else
      i1=stib(xtstrl(0)+18)
      j3=stib(xtstrp(0)+18)
      stxb(i2+1:i2+i1)=stcb(j3:j3-1+i1)
      i2=i2+i1
      endif
      endif
      if(i2+1.le.slin2)then
      stxb(i2+1:i2+1)=char(aminus)
      endif
      i2=i2+1
      if(i2+wztos(nl2).le.slin2)then
      if(slin3.ge.0)then
      call dkar(nl2,jj)
      j1=stcbs(1)
      call ctxb(i2,jj,j1)
      else
      i1=stib(xtstrl(0)+18)
      j3=stib(xtstrp(0)+18)
      stxb(i2+1:i2+i1)=stcb(j3:j3-1+i1)
      i2=i2+i1
      endif
      endif
      endif
      endif
      if(i2.ge.slin1)then
      j1=slin1
   60 continue
      if(j1.gt.0)then
      if(ichar(stxb(j1:j1)).ne.aspace)then
      j1=j1-1
      goto 60
      endif
      endif
      if(j1.gt.1)then
      write(unit=*,fmt='(a)')stxb(1:j1-1)
      endif
      write(unit=*,fmt='(a)')stcb(1:stab)//stxb(j1:i2)
      else
      write(unit=*,fmt='(a)')stxb(1:i2)
      endif
   70 continue
      call hrul(2)
      write(unit=*,fmt='(a,/)')stxb(1:ssrec-1)
      jflag(7)=10
      if(istop.ne.0)then
      if(istop.gt.0)then
      call qclose(0,istop)
      endif
      stop
      endif
      return
      end
      subroutine hrul(ii)
      implicit integer(a-z)
      save
      parameter ( srec=81, ssrec=62 )
      parameter ( sxbuff=2040 )
      character*(sxbuff) stxb
      common/z27g/stxb
      common/z51g/aeq,agt,alt,aminus,anine,aplus,aslash,atimes,azero
      common/z52g/acf(0:127),abo(0:6),acr,alf,aspace,dquote,squote
      stxb(1:1)=char(aspace)
      j1=ssrec
      if(ii.eq.1)then
      do i1=2,j1
      stxb(i1:i1)=char(aminus)
      enddo
      elseif(ii.eq.2)then
      do i1=3,j1-1
      stxb(i1:i1)=char(aeq)
      enddo
      stxb(2:2)=char(aspace)
      stxb(j1:j1)=char(aspace)
      endif
      return
      end
      subroutine uput(ind)
      implicit integer(a-z)
      save
      parameter ( srec=81, ssrec=62 )
      parameter ( sxbuff=2040 )
      common/z12g/cflag(1:7),dflag(1:20),jflag(1:10),mflag(1:6)
      character*(srec) mlin
      common/z22g/mlin
      character*(sxbuff) stxb
      common/z27g/stxb
      common/z52g/acf(0:127),abo(0:6),acr,alf,aspace,dquote,squote
      if(ind.eq.1)then
      mlin(1:srec)='system/filesystem I/O error'
      elseif(ind.eq.2)then
      mlin(1:srec)='invalid character set'
      elseif(ind.eq.3)then
      mlin(1:srec)='largest available integer is too small'
      elseif(ind.eq.4)then
      mlin(1:srec)="parameter 'sibuff' is too small"
      elseif(ind.eq.5)then
      mlin(1:srec)="parameter 'scbuff' is too small"
      elseif(ind.eq.6)then
      mlin(1:srec)="parameter 'sxbuff' is too small"
      elseif(ind.eq.7)then
      mlin(1:srec)="parameter 'sxbuff' not properly set"
      elseif(ind.eq.8)then
      mlin(1:srec)="parameter 'scbuff' not properly set"
      elseif(ind.eq.9)then
      mlin(1:srec)="parameter 'sibuff' not properly set"
      elseif(ind.eq.10)then
      mlin(1:srec)="parameter 'srec' not properly set"
      elseif(ind.eq.11)then
      mlin(1:srec)="parameter 'ssrec' not properly set"
      else
      mlin(1:srec)='uput_1'
      endif
      ii=sxbuff-max(80,ssrec)
      do i1=srec,1,-1
      if(ichar(mlin(i1:i1)).ne.aspace)then
      if(ii.gt.0)then
      call hrul(2)
      write(unit=*,fmt='(/,a,/,a,/,a,/)')stxb(1:ssrec-1),
     :'   error: '//mlin(1:i1),stxb(1:ssrec-1)
      else
      write(unit=*,fmt='(//,a,/)')'   error: '//mlin(1:i1)
      endif
      goto 90
      endif
      enddo
   90 stop
      end
      subroutine wput(sind,j1,j2)
      implicit integer(a-z)
      save
      parameter ( sibuff=1048574, scbuff=131064, swbuff=65528 )
      parameter ( srec=81, ssrec=62 )
      common/z12g/cflag(1:7),dflag(1:20),jflag(1:10),mflag(1:6)
      common/z15g/pten(1:10),iref,wiref,wsint
      character*(srec) mlin
      common/z22g/mlin
      common/z28g/wera(0:0),werb(0:0),nwer
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z35g/stibs(1:1),stcbs(1:1),stwbs(1:1)
      if(sind.gt.0)then
      if(cflag(1).le.0)then
      jflag(8)=min(jflag(8)+1,iref)
      goto 90
      endif
      endif
      aind=abs(sind)
      ii=0
      if(aind.gt.nwer)then
      ii=1
      elseif(j1.lt.0)then
      ii=1
      elseif(j1.eq.0)then
      if(j2.ne.0)then
      ii=1
      endif
      else
      jj=j2-j1
      if(jj.lt.0)then
      ii=1
      elseif(jj.ge.srec-1)then
      ii=1
      endif
      endif
      if(ii.ne.0)then
      mlin(1:srec)='wput_1'
      call mput(1,0,0,0)
      endif
      i1=stib(wera(0)+aind)
      i2=stib(werb(0)+aind)
      ll=i2-i1
      ii=0
      if(i1.le.0)then
      ii=1
      elseif(i1.ge.stcbs(1))then
      ii=1
      elseif(ll.lt.0)then
      ii=1
      elseif(ll.ge.srec)then
      ii=1
      endif
      if(ii.ne.0)then
      mlin(1:srec)='wput_2'
      call mput(1,0,0,0)
      endif
      ll=ll+1
      if(j1.ne.0)then
      jj=(ll-j1)+(j2+1)
      else
      jj=0
      endif
      kk=srec-6
      if(jj.gt.0)then
      if(jj.lt.srec)then
      mlin(1:srec)=stcb(i1:i2)//stcb(j1:j2)
      elseif(ll.le.kk)then
      mlin(1:srec)=stcb(i1:i2)//stcb(j1:j1+(kk-ll))//' ...'
      else
      mlin(1:srec)=stcb(i1:i1+kk)//' ...'
      endif
      else
      if(ll.lt.srec)then
      mlin(1:srec)=stcb(i1:i2)
      else
      mlin(1:srec)=stcb(i1:i1+kk)//' ...'
      endif
      endif
      if(sind.gt.0)then
      call mput(0,0,0,0)
      else
      call mput(1,0,0,0)
      endif
   90 return
      end
      subroutine ccyc(xx,situ)
      implicit integer(a-z)
      save
      parameter ( maxleg=10, maxrho=7, maxi=10, maxdeg=6 )
      parameter ( maxn=maxi+maxi-2, maxli=maxn+maxrho )
      parameter ( srec=81, ssrec=62 )
      common/z1in/leg(1:maxleg),nleg,nloop,incom
      common/z1g/g(1:maxn,1:maxn),rho(1:maxdeg),nrho,rhop1,rhop2
      common/z4g/n,nli
      common/z18g/eg(1:maxn,1:maxn),flow(1:maxli,0:maxleg+maxrho),
     :net(-3:3)
      character*(srec) mlin
      common/z22g/mlin
      integer bb(1:maxli),cc(1:maxli)
      if(xx.eq.1)then
      if(nloop.le.1)then
      situ=1
      goto 90
      endif
      else
      mlin(1:srec)='ccyc_2'
      call mput(1,0,0,0)
      endif
      j1=0
      do i1=rhop1,nli
      if(flow(i1,0).lt.0)then
      j1=j1+1
      bb(j1)=j1
      cc(j1)=i1
      endif
      enddo
      jj=1-net(-3)
      if(jj.ge.0)then
      goto 20
      endif
      k1=1
      k2=net(-3)
      do i1=1,net(-3)
      j1=cc(i1)
      do i2=i1+1,net(-3)
      if(bb(i1).ne.bb(i2))then
      kk=0
      j2=cc(i2)
      do i3=rhop1,rhop2
      if(flow(j1,i3).ne.0)then
      if(flow(j2,i3).ne.0)then
      kk=1
      goto 15
      endif
      endif
      enddo
   15 continue
      if(kk.ne.0)then
      if(bb(k1).eq.1)then
      k1=k1+1
      endif
      if(bb(k2).eq.1)then
      k2=k2-1
      endif
      j2=min(bb(i1),bb(i2))
      j3=max(bb(i1),bb(i2))
      do i3=k1,k2
      if(bb(i3).eq.j3)then
      bb(i3)=j2
      endif
      enddo
      jj=jj+1
      if(jj.eq.0)then
      goto 20
      endif
      endif
      endif
      enddo
      enddo
   20 continue
      if(xx.eq.1)then
      if(jj.lt.0)then
      situ=-1
      else
      situ=1
      endif
      endif
   90 return
      end
      subroutine model
      implicit integer(a-z)
      save
      parameter ( maxleg=10, maxrho=7, maxi=10, maxdeg=6 )
      parameter ( maxn=maxi+maxi-2, maxli=maxn+maxrho )
      parameter ( srec=81, ssrec=62 )
      parameter ( eoa=-2047, nap=-8191 )
      parameter ( nfiles=5 )
      parameter ( sibuff=1048574, scbuff=131064, swbuff=65528 )
      common/z1in/leg(1:maxleg),nleg,nloop,incom
      common/z3in/lmfile,mfilea,mfileb
      common/z4in/llfile,lfilea,lfileb
      common/z7in/aunit(1:nfiles)
      common/z1g/g(1:maxn,1:maxn),rho(1:maxdeg),nrho,rhop1,rhop2
      common/z3g/nivd(0:0),dpntro(0:0),vparto(0:0),vval(0:0),nvert
      common/z11g/nphi,nblok,nprop,npprop
      character*(srec) mlin
      common/z22g/mlin
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z33g/namep(0:0),namel(0:0)
      common/z34g/link(0:0),antiq(0:0),blok(0:0)
      common/z35g/stibs(1:1),stcbs(1:1),stwbs(1:1)
      common/z36g/udkp(0:2),udkl(0:0),udkt(0:0),udki(0:0)
      common/z38g/gmkp(0:0),gmkl(0:0),gmkd(0:0),gmko(0:0),
     :gmkvp(0:0),gmkvl(0:0)
      common/z39g/pmkp(0:0),pmkl(0:0),pmkd(0:0),pmkvpp(0:0),pmkvlp(0:0)
      common/z40g/vmkp(0:0),vmkl(0:0),vmkvpp(0:0),vmkvlp(0:0),vmks(0:0)
      common/z41g/nudk,ngmk,npmk0,npmk,nvmk0,nvmk
      common/z46g/popt3(0:0),wopt3(0:0),copt3(0:0),popt9(0:0),
     :wopt9(0:0),copt9(0:0)
      common/z51g/aeq,agt,alt,aminus,anine,aplus,aslash,atimes,azero
      common/z52g/acf(0:127),abo(0:6),acr,alf,aspace,dquote,squote
      integer apmkd(0:0),avmkd(0:0),llp(0:0)
      integer sucpal(0:11,0:11),acf0(0:127)
      do i1=0,31
      acf0(i1)=-1
      enddo
      do i1=32,126
      acf0(i1)=0
      enddo
      acf0(127)=-1
      acf0(alf)=1
      acf0(aspace)=2
      acf0(dquote)=3
      acf0(squote)=4
      acf0(ichar('('))=5
      acf0(ichar(')'))=6
      acf0(ichar(','))=7
      acf0(ichar(';'))=8
      acf0(aeq)=9
      acf0(ichar('['))=10
      acf0(ichar(']'))=11
      do i1=0,11
      do i2=0,11
      sucpal(i1,i2)=-1
      enddo
      enddo
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
      if(llfile.gt.0)then
      xfile=0
      else
      xfile=1
      endif
   05 continue
      xfile=xfile+1
      if(xfile.eq.1)then
      call qopen(lfilea,llfile,3,0)
      elseif(xfile.eq.2)then
      call qopen(mfilea,lmfile,2,0)
      endif
      nlin=0
      dbl=0
      newc=1
      level=1
      if(xfile.eq.1)then
      lastp=0
      lastkp=0
      else
      llp(0)=nap
      ngmk=0
      npmk=0
      nvmk=0
      nrho=0
      nphi=0
      nprop=0
      nvert=0
      endif
   65 continue
      dip=0
      nest=0
   70 continue
      if(xfile.eq.1)then
      goto 05
      elseif(xfile.eq.2)then
      call qrlin(2,nlin,slin,qc)
      if(slin.eq.-1)then
      goto 777
      endif
      endif
      if(newc.ne.0)then
      dbl=nlin
      endif
      if((slin.eq.0).or.(qc.ne.0))then
      if(newc.ne.0)then
      goto 70
      else
      goto 521
      endif
      endif
      ibot=stcbs(1)+1
      itop=stcbs(1)+slin
      do i1=ibot,itop
      if(ichar(stcb(i1:i1)).ne.aspace)then
      ibot=i1
      goto 16
      endif
      enddo
   16 continue
      ii=0
      quote=0
      do i1=ibot,itop
      jj=ichar(stcb(i1:i1))
      if(jj.eq.squote)then
      if(nest.ne.1)then
      goto 521
      endif
      if(quote.eq.0)then
      quote=squote
      endif
      if(quote.eq.squote)then
      ii=1-ii
      endif
      elseif(jj.eq.dquote)then
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
      if(stcb(i1:i1).eq.'[')then
      if(nest.ne.0)then
      goto 521
      endif
      nest=1
      elseif(stcb(i1:i1).eq.']')then
      if(nest.ne.1)then
      goto 521
      endif
      nest=2
      elseif(stcb(i1:i1).eq.',')then
      if(nest.ne.1)then
      goto 521
      endif
      quote=0
      elseif(jj.ne.aspace)then
      if(nest.ne.1)then
      goto 521
      endif
      endif
      else
      if(nest.ne.1)then
      goto 521
      endif
      endif
      enddo
      if(ii.ne.0)then
      goto 521
      endif
      if(dbl.eq.nlin)then
      jbot=ibot
      endif
      jj=slin+1
      call aocb(jj)
      if(stcb(itop:itop).ne.']')then
      newc=0
      goto 70
      endif
      if(nest.ne.2)then
      goto 521
      endif
      newc=1
      ibot=jbot
      top=stcbs(1)
   50 continue
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
      stib(stibs(1))=eoa
      llp(0)=stibs(1)
      endif
      call rgmki(llp(0))
      if(xfile.eq.1)then
      goto 521
      endif
      level=2
      fpass=1
      call aoib(1)
      stib(stibs(1))=eoa
      llp(0)=stibs(1)
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
      if(ichar(stcb(ix:ix)).eq.squote)then
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
      if(ichar(stcb(ix:ix)).eq.squote)then
      plm=1-plm
      endif
      is2=ix
      endif
      elseif(pal.lt.9)then
      if(pal.eq.6)then
      if(vel.ne.1.or.plm.ne.0)then
      goto 521
      endif
      if(ichar(stcb(ix:ix)).eq.squote)then
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
      if(ichar(stcb(ix:ix)).eq.squote)then
      plm=1-plm
      endif
      is2=ix
      endif
      elseif(pal.lt.12)then
      if(pal.eq.9)then
      if(level.gt.2.or.level.le.0)then
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
      if(fpass.eq.0)then
      do i1=1,npmk
      if(stib(apmkd(0)+i1).eq.3)then
      if(stib(lastp+1).eq.0)then
      mlin(1:srec)='wrong dimension for M-function,'
      call mput(1,dbl,nlin,4-xfile)
      endif
      endif
      if(stib(apmkd(0)+i1).ne.stib(pmkd(0)+i1))then
      ii=1
      if(stib(apmkd(0)+i1).eq.2)then
      if(stib(pmkd(0)+i1).eq.3)then
      ii=stib(lastp+1)
      endif
      endif
      if(ii.ne.0)then
      if(stib(apmkd(0)+i1).eq.0)then
      mlin(1:srec)='missing M-function,'
      else
      mlin(1:srec)='wrong dimension for M-function,'
      endif
      call mput(1,dbl,nlin,2)
      endif
      endif
      enddo
      nprop=nprop+1
      endif
      elseif(level.eq.3)then
      if(fpass.eq.0)then
      do i1=1,nvmk
      if(stib(avmkd(0)+i1).eq.0)then
      mlin(1:srec)='missing M-function,'
      call mput(1,dbl,nlin,2)
      endif
      enddo
      jj=stib(lastp+1)
      if(jj.lt.3)then
      mlin(1:srec)='vertex degree is too small,'
      call mput(1,dbl,nlin,2)
      endif
      ii=0
      do i1=1,jj
      if(stib(antiq(0)+stib(lastp+1+i1)).ne.0)then
      ii=1-ii
      endif
      enddo
      if(ii.ne.0)then
      mlin(1:srec)='odd vertex,'
      call mput(1,dbl,nlin,2)
      endif
      endif
      endif
      if(level.eq.2)then
      if(fpass.eq.1)then
      ii=stibs(1)+1
      call aoib(4)
      do i1=ii,stibs(1)
      stib(i1)=eoa
      enddo
      call trm(4,npmk+1)
      ii=npmk+1
      pmkd(0)=stibs(1)-ii
      pmkl(0)=pmkd(0)-ii
      pmkp(0)=pmkl(0)-ii
      jj=pmkp(0)-llp(0)+1
      call xipht(pmkp(0)+1,stibs(1),-jj)
      call aoib(-jj)
      pmkd(0)=pmkd(0)-jj
      pmkl(0)=pmkl(0)-jj
      pmkp(0)=pmkp(0)-jj
      psize=6+2*npmk
      ii=npmk+1
      apmkd(0)=stibs(1)
      call aoib(ii)
      stib(stibs(1))=eoa
      endif
      do i1=1,npmk
      stib(apmkd(0)+i1)=0
      enddo
      if(fpass.eq.1)then
      fpass=0
      if(idin.lt.3)then
      goto 521
      endif
      goto 50
      endif
      elseif(level.eq.3)then
      if(fpass.eq.1)then
      ii=stibs(1)+1
      call aoib(3)
      do i1=ii,stibs(1)
      stib(i1)=eoa
      enddo
      call trm(3,nvmk+1)
      ii=nvmk+1
      vmkl(0)=stibs(1)-ii
      vmkp(0)=vmkl(0)-ii
      avmkd(0)=vmkp(0)-ii
      call xipht(vmkp(0),stibs(1),-1)
      call aoib(-1)
      vmkl(0)=vmkl(0)-1
      vmkp(0)=vmkp(0)-1
      avmkd(0)=avmkd(0)-1
      else
      nvert=nvert+1
      endif
      do i1=1,nvmk
      stib(avmkd(0)+i1)=0
      enddo
      if(fpass.eq.1)then
      fpass=0
      goto 50
      endif
      endif
      goto 65
  281 continue
      if(idin.eq.3)then
      if(stdw(stcb,id1,id2).eq.0)then
      if(id1.ne.id2)then
      goto 521
      endif
      ii=ichar(stcb(id1:id1))
      if((ii.ne.aplus).and.(ii.ne.aminus))then
      goto 521
      endif
      else
      if(knf.ne.2)then
      mlin(1:srec)='unknown field in vertex,'
      call mput(1,dbl,nlin,2)
      endif
      knf=0
      call rpi
      level=3
      fpass=1
      lastp=0
      lastkp=0
      goto 50
      endif
      if(knf.ne.0)then
      mlin(1:srec)='field appeared in previous propagator,'
      call mput(1,dbl,nlin,2)
      endif
      if(fpass.eq.0)then
      ii=ichar(stcb(id1:id1))
      if(ii.eq.aplus)then
      stib(off1+3)=0
      elseif(ii.eq.aminus)then
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
      stib(off1+1)=eoa
      lastp=off1+1
      if(off2.ne.0)then
      stib(lastp)=off2+1
      stib(off2+1)=eoa
      lastp=off2+1
      endif
      endif
      goto 284
      elseif(idin.eq.4)then
      if(fpass.eq.0)then
      call mstr0(stcb,id1,id2,popt9(0),wopt9(0),ij)
      if(ij.eq.0)then
      goto 521
      endif
      stib(off1+4)=stib(copt9(0)+ij)
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
      mlin(1:srec)='unacceptable field name,'
      call mput(1,dbl,nlin,2)
      endif
      if(fpass.eq.1)then
      goto 284
      endif
      kk=id2-id1+1
      if(kk.le.0)then
      goto 521
      endif
      newid=1
      if(nphi.gt.0)then
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
      elseif(ij.ne.nphi)then
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
      stib(stibs(1))=eoa
      llp(0)=stibs(1)
      lastp=stibs(1)
      endif
      off1=stibs(1)
      off2=0
      call aoib(psize)
      do i1=stibs(1)-psize+1,stibs(1)
      stib(i1)=0
      enddo
      stib(off1+5)=id1
      stib(off1+6)=id2-id1+1
      nphi=nphi+1
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
      do i1=stibs(1)-psize+1,stibs(1)
      stib(i1)=0
      enddo
      stib(off2+5)=id1
      stib(off2+6)=id2-id1+1
      stib(off1+2)=1
      stib(off2+2)=-1
      nphi=nphi+1
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
      mlin(1:srec)='unacceptable field name,'
      call mput(1,dbl,nlin,2)
      endif
      if(idin.gt.nrho)then
      if(idin.gt.maxdeg)then
      mlin(1:srec)='vertex degree is too large,'
      call mput(1,dbl,nlin,2)
      endif
      nrho=idin
      endif
      if(fpass.eq.1)then
      goto 286
      endif
      if(idin.eq.1)then
      if(nvert.eq.0)then
      call aoib(1)
      stib(stibs(1))=eoa
      llp(0)=stibs(1)
      lastp=stibs(1)
      endif
      ii=stibs(1)+1
      jj=2+2*nvmk
      call aoib(jj)
      stib(lastp)=ii
      lastp=ii
      stib(lastp)=eoa
      do i1=lastp+1,stibs(1)
      stib(i1)=0
      enddo
      endif
      call mstr0(stcb,id1,id2,namep(0),namel(0),ij)
      if(ij.eq.0)then
      mlin(1:srec)='unknown field in vertex,'
      call mput(1,dbl,nlin,2)
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
      stib(stibs(1))=eoa
      llp(0)=stibs(1)
      lastkp=stibs(1)
      else
      call mstr1(stcb,id1,id2,llp(0),1,2,ij)
      if(ij.ne.0)then
      mlin(1:srec)='multiply defined M-function,'
      call mput(1,dbl,nlin,4-xfile)
      endif
      endif
      ii=stibs(1)+1
      call aoib(4)
      call mstr0(stcb,id1,id2,udkp(0),udkl(0),ij)
      if(ij.eq.0)then
      stib(stibs(1)-2)=id1
      else
      stib(stibs(1)-2)=stib(udkp(0)+ij)
      endif
      stib(stibs(1)-1)=id2-id1+1
      stib(lastkp)=ii
      lastkp=ii
      stib(lastkp)=eoa
      ngmk=ngmk+1
      stib(stibs(1))=0
      elseif(level.eq.2)then
      if(fpass.eq.1)then
      call mstr0(stcb,id1,id2,gmkp(0),gmkl(0),ij)
      if(ij.ne.0)then
      mlin(1:srec)='multiply defined M-function,'
      call mput(1,dbl,nlin,2)
      endif
      if(npmk.eq.0)then
      call aoib(1)
      stib(stibs(1))=eoa
      llp(0)=stibs(1)
      lastkp=stibs(1)
      else
      call mstr1(stcb,id1,id2,llp(0),1,2,ij)
      if(ij.ne.0)then
      mlin(1:srec)='multiply defined M-function,'
      call mput(1,dbl,nlin,2)
      endif
      endif
      call mstr0(stcb,id1,id2,udkp(0),udkl(0),ij)
      ii=stibs(1)+1
      call aoib(4)
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
      stib(lastkp)=eoa
      npmk=npmk+1
      else
      call mstr0(stcb,id1,id2,pmkp(0),pmkl(0),kid)
      if(kid.eq.0)then
      mlin(1:srec)='unexpected M-function,'
      call mput(1,dbl,nlin,2)
      endif
      if(stib(apmkd(0)+kid).ne.0)then
      mlin(1:srec)='multiply defined M-function,'
      call mput(1,dbl,nlin,2)
      endif
      endif
      elseif(level.eq.3)then
      if(fpass.eq.1)then
      call mstr0(stcb,id1,id2,gmkp(0),gmkl(0),ij)
      if(ij.ne.0)then
      mlin(1:srec)='multiply defined M-function,'
      call mput(1,dbl,nlin,2)
      endif
      call mstr0(stcb,id1,id2,pmkp(0),pmkl(0),ij)
      if(ij.ne.0)then
      mlin(1:srec)='multiply defined M-function,'
      call mput(1,dbl,nlin,2)
      endif
      if(nvmk.eq.0)then
      call aoib(1)
      stib(stibs(1))=eoa
      llp(0)=stibs(1)
      lastkp=stibs(1)
      else
      call mstr1(stcb,id1,id2,llp(0),1,2,ij)
      if(ij.ne.0)then
      mlin(1:srec)='multiply defined M-function,'
      call mput(1,dbl,nlin,2)
      endif
      endif
      call mstr0(stcb,id1,id2,udkp(0),udkl(0),ij)
      ii=stibs(1)+1
      call aoib(3)
      if(ij.eq.0)then
      stib(stibs(1)-1)=id1
      stib(stibs(1))=id2-id1+1
      else
      stib(stibs(1)-1)=stib(udkp(0)+ij)
      stib(stibs(1))=stib(udkl(0)+ij)
      endif
      stib(lastkp)=ii
      lastkp=ii
      stib(lastkp)=eoa
      nvmk=nvmk+1
      else
      call mstr0(stcb,id1,id2,vmkp(0),vmkl(0),kid)
      if(kid.eq.0)then
      mlin(1:srec)='unexpected M-function,'
      call mput(1,dbl,nlin,2)
      endif
      if(stib(avmkd(0)+kid).ne.0)then
      mlin(1:srec)='multiply defined M-function,'
      call mput(1,dbl,nlin,2)
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
      j1=lastkp+3
      if(vel.eq.0)then
      stib(j1)=1
      else
      stib(j1)=vin+1
      endif
      if(vel.ne.0)then
      mlin(1:srec)='wrong dimension for M-keyword,'
      call mput(1,dbl,nlin,2)
      endif
      elseif(level.eq.2)then
      if(fpass.eq.1)then
      j1=lastkp+3
      if(vel.eq.0)then
      stib(j1)=1
      else
      stib(j1)=3
      endif
      goto 184
      else
      j1=apmkd(0)+kid
      if(vel.eq.0)then
      stib(j1)=1
      else
      stib(j1)=vin+1
      endif
      if(stib(j1).gt.stib(pmkd(0)+kid))then
      mlin(1:srec)='wrong dimension for M-function,'
      call mput(1,dbl,nlin,2)
      endif
      endif
      if(vin.gt.2)then
      mlin(1:srec)='wrong dimension for M-function,'
      call mput(1,dbl,nlin,2)
      endif
      elseif(level.eq.3)then
      if(vel.ne.0)then
      mlin(1:srec)='wrong dimension for M-function,'
      call mput(1,dbl,nlin,2)
      endif
      if(fpass.eq.1)then
      goto 184
      else
      stib(avmkd(0)+kid)=1
      endif
      endif
      if(ichar(stcb(is1:is1)).eq.squote)then
      kk=stds(stcb,is1,is2,1)
      if(kk.lt.0)then
      goto 521
      elseif(kk.gt.0)then
      stcb(is1+1:is1+kk)=stcb(stcbs(1)+1:stcbs(1)+kk)
      endif
      stcb(is1:is1)=char(172)
      stcb(is1+kk+1:is1+kk+1)=char(172)
      is1=is1+1
      ii=-(is2-is1-kk)
      if(ii.lt.0)then
      call cxipht(is2+1,stcbs(1),ii)
      call aocb(ii)
      ix=ix+ii
      top=top+ii
      if(level.eq.2)then
      if(nprop.eq.0)then
      do i1=1,npmk
      if(stib(pmkp(0)+i1).gt.is2)then
      stib(pmkp(0)+i1)=stib(pmkp(0)+i1)+ii
      endif
      enddo
      endif
      elseif(level.eq.3)then
      if(nvert.eq.0)then
      do i1=1,nvmk
      if(stib(vmkp(0)+i1).gt.is2)then
      stib(vmkp(0)+i1)=stib(vmkp(0)+i1)+ii
      endif
      enddo
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
      jj=5+2*kid
      if(vel.ne.0)then
      if(vin.eq.1)then
      j1=off1+jj
      elseif(vin.eq.2)then
      j1=off2+jj
      endif
      stib(j1)=is1
      stib(j1+1)=kk
      else
      stib(off1+jj)=is1
      stib(off1+jj+1)=kk
      if(off2.ne.0)then
      stib(off2+jj)=is1
      stib(off2+jj+1)=kk
      endif
      endif
      elseif(level.eq.3)then
      jj=lastp+stib(lastp+1)+2*kid
      stib(jj)=is1
      stib(jj+1)=kk
      endif
  184 continue
      is1=0
      goto 287
  521 continue
      mlin(1:srec)='wrong syntax,'
      call mput(1,dbl,nlin,4-xfile)
  777 continue
      call rvi(llp(0))
      return
      end
      subroutine rgmki(llpp)
      implicit integer(a-z)
      save
      parameter ( sibuff=1048574, scbuff=131064, swbuff=65528 )
      parameter ( eoa=-2047, nap=-8191 )
      common/z30g/stib(1:sibuff)
      common/z35g/stibs(1:1),stcbs(1:1),stwbs(1:1)
      common/z38g/gmkp(0:0),gmkl(0:0),gmkd(0:0),gmko(0:0),
     :gmkvp(0:0),gmkvl(0:0)
      common/z41g/nudk,ngmk,npmk0,npmk,nvmk0,nvmk
      if(ngmk.eq.0)then
      gmkp(0)=nap
      gmkl(0)=nap
      gmkd(0)=nap
      gmko(0)=nap
      gmkvp(0)=nap
      gmkvl(0)=nap
      call aoib(-1)
      goto 90
      endif
      ii=ngmk+1
      gmkp(0)=stibs(1)
      call aoib(ii)
      stib(stibs(1))=eoa
      gmkl(0)=stibs(1)
      call aoib(ii)
      stib(stibs(1))=eoa
      gmkd(0)=stibs(1)
      call aoib(ii)
      stib(stibs(1))=eoa
      gmko(0)=stibs(1)
      call aoib(ii)
      stib(stibs(1))=eoa
      j1=0
      tgmkd=0
      jj=llpp
   10 continue
      jj=stib(jj)
      if(jj.ne.eoa)then
      j1=j1+1
      stib(gmkp(0)+j1)=stib(jj+1)
      stib(gmkl(0)+j1)=stib(jj+2)
      stib(gmkd(0)+j1)=stib(jj+3)
      stib(gmko(0)+j1)=tgmkd
      tgmkd=tgmkd+max(1,stib(jj+3)-1)
      goto 10
      endif
      ii=tgmkd+1
      gmkvp(0)=stibs(1)
      call aoib(ii)
      stib(stibs(1))=eoa
      gmkvl(0)=stibs(1)
      call aoib(ii)
      stib(stibs(1))=eoa
      j1=gmkd(0)
      j2=0
      jj=llpp
   20 continue
      jj=stib(jj)
      if(jj.ne.eoa)then
      j1=j1+1
      j3=jj+3
      do i1=1,max(1,stib(j1)-1)
      j2=j2+1
      j3=j3+1
      stib(gmkvp(0)+j2)=stib(j3)
      j3=j3+1
      stib(gmkvl(0)+j2)=stib(j3)
      if(stib(j3).eq.0)then
      i2=stib(gmkp(0)+j2)
      i3=i2-1+stib(gmkl(0)+j2)
      call wput(8,i2,i3)
      endif
      enddo
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
   90 continue
      return
      end
      subroutine rpi
      implicit integer(a-z)
      save
      parameter ( sibuff=1048574, scbuff=131064, swbuff=65528 )
      parameter ( srec=81, ssrec=62 )
      parameter ( eoa=-2047, nap=-8191 )
      common/z9g/tpc(0:0)
      common/z11g/nphi,nblok,nprop,npprop
      common/z12g/cflag(1:7),dflag(1:20),jflag(1:10),mflag(1:6)
      character*(srec) mlin
      common/z22g/mlin
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z33g/namep(0:0),namel(0:0)
      common/z34g/link(0:0),antiq(0:0),blok(0:0)
      common/z35g/stibs(1:1),stcbs(1:1),stwbs(1:1)
      common/z39g/pmkp(0:0),pmkl(0:0),pmkd(0:0),pmkvpp(0:0),pmkvlp(0:0)
      common/z41g/nudk,ngmk,npmk0,npmk,nvmk0,nvmk
      common/z42g/pmkr(0:0),pmkvma(0:0),pmkvmi(0:0)
      if(nphi.eq.0)then
      mlin(1:srec)='no propagators listed in'
      call mput(1,0,0,2)
      endif
      ii=stibs(1)
      psize=6+2*npmk
      call aoib(psize)
      do i1=ii+1,stibs(1)
      stib(i1)=eoa
      enddo
      jj=nphi+1
      call trm(psize,jj)
      blok(0)=stibs(1)-psize*jj
      link(0)=blok(0)+jj
      antiq(0)=link(0)+jj
      tpc(0)=antiq(0)+jj
      namep(0)=tpc(0)+jj
      namel(0)=namep(0)+jj
      if(stib(blok(0)).ne.eoa)then
      stib(blok(0))=eoa
      endif
      npprop=0
      do i1=1,nphi
      j1=stib(tpc(0)+i1)
      j2=stib(link(0)+i1)
      stib(link(0)+i1)=j2+i1
      if(j1.eq.5)then
      mflag(5)=1
      else
      if(j1.eq.1)then
      mflag(6)=1
      endif
      if(j2.ge.0)then
      npprop=npprop+1
      endif
      endif
      enddo
      if(npmk.eq.0)then
      pmkvpp(0)=nap
      pmkvlp(0)=nap
      pmkr(0)=nap
      pmkvmi(0)=nap
      pmkvma(0)=nap
      goto 90
      endif
      ii=npmk+1
      pmkvpp(0)=stibs(1)
      call aoib(ii)
      stib(stibs(1))=eoa
      pmkvlp(0)=stibs(1)
      call aoib(ii)
      stib(stibs(1))=eoa
      pmkr(0)=stibs(1)
      call aoib(ii)
      stib(stibs(1))=eoa
      ii=nphi+1
      jj=namel(0)
      do i1=1,npmk
      jj=jj+ii
      stib(pmkvpp(0)+i1)=jj
      jj=jj+ii
      stib(pmkvlp(0)+i1)=jj
      enddo
      do i1=1,npmk
      ii=stib(pmkvpp(0)+i1)
      jj=stib(pmkvlp(0)+i1)
      ij=0
      do i2=4,2,-2
      if(ij.eq.0)then
      do i3=1,nphi
      j1=stib(ii+i3)
      j2=stib(ii+i3)-1+stib(jj+i3)
      if(j1.gt.j2)then
      ij=1
      goto 40
      elseif(i2.eq.4)then
      if(stdz(stcb,j1,j2).eq.0)then
      goto 40
      endif
      elseif(i2.eq.2)then
      if(stdw(stcb,j1,j2).eq.0)then
      goto 40
      endif
      endif
      enddo
      ij=i2
   40 continue
      endif
      enddo
      if(ij.eq.0)then
      ij=1
      endif
      stib(pmkr(0)+i1)=ij
      enddo
      ii=npmk+1
      pmkvmi(0)=stibs(1)
      call aoib(ii)
      stib(stibs(1))=eoa
      pmkvma(0)=stibs(1)
      call aoib(ii)
      stib(stibs(1))=eoa
      do i1=1,npmk
      if(stib(pmkr(0)+i1).eq.4)then
      ii=stib(pmkvpp(0)+i1)
      jj=stib(pmkvlp(0)+i1)
      do i2=1,nphi
      j1=stib(ii+i2)
      j2=stib(ii+i2)-1+stib(jj+i2)
      ij=stoz(stcb,j1,j2)
      if(i2.gt.1)then
      if(ij.lt.j3)then
      j3=ij
      elseif(ij.gt.j4)then
      j4=ij
      endif
      else
      j3=ij
      j4=ij
      endif
      enddo
      if(j3.eq.0)then
      if(j4.eq.0)then
      j1=stib(pmkp(0)+i1)
      j2=j1-1+stib(pmkl(0)+i1)
      call wput(10,j1,j2)
      endif
      endif
      else
      j3=0
      j4=0
      endif
      stib(pmkvmi(0)+i1)=j3
      stib(pmkvma(0)+i1)=j4
      enddo
      do i1=1,npmk
      if(stib(pmkr(0)+i1).eq.1)then
      jj=stib(pmkvlp(0)+i1)
      do i2=1,nphi
      if(stib(jj+i2).gt.0)then
      goto 85
      endif
      enddo
      j1=stib(pmkp(0)+i1)
      j2=j1-1+stib(pmkl(0)+i1)
      call wput(9,j1,j2)
      endif
   85 continue
      enddo
   90 continue
      return
      end
      subroutine rvi(llp)
      implicit integer(a-z)
      save
      parameter ( maxleg=10, maxrho=7, maxi=10, maxdeg=6 )
      parameter ( maxn=maxi+maxi-2, maxli=maxn+maxrho )
      parameter ( sibuff=1048574, scbuff=131064, swbuff=65528 )
      parameter ( srec=81, ssrec=62 )
      parameter ( eoa=-2047, nap=-8191 )
      common/z1g/g(1:maxn,1:maxn),rho(1:maxdeg),nrho,rhop1,rhop2
      common/z3g/nivd(0:0),dpntro(0:0),vparto(0:0),vval(0:0),nvert
      common/z11g/nphi,nblok,nprop,npprop
      common/z12g/cflag(1:7),dflag(1:20),jflag(1:10),mflag(1:6)
      character*(srec) mlin
      common/z22g/mlin
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z34g/link(0:0),antiq(0:0),blok(0:0)
      common/z35g/stibs(1:1),stcbs(1:1),stwbs(1:1)
      common/z40g/vmkp(0:0),vmkl(0:0),vmkvpp(0:0),vmkvlp(0:0),vmks(0:0)
      common/z41g/nudk,ngmk,npmk0,npmk,nvmk0,nvmk
      common/z43g/vmkr(0:0),vmkmao(0:0),vmkmio(0:0)
      integer xtmp(1:maxdeg+1),nrot(1:maxdeg),rotvpo(1:maxdeg)
      if(nvert.eq.0)then
      mlin(1:srec)='no vertices listed in'
      call mput(1,0,0,2)
      endif
      ii=nrho+1
      nivd(0)=stibs(1)
      call aoib(ii)
      stib(stibs(1))=eoa
      ii=nvert+1
      vval(0)=stibs(1)
      call aoib(ii)
      stib(stibs(1))=eoa
      vparto(0)=stibs(1)
      call aoib(ii)
      stib(stibs(1))=eoa
      jj=nvmk+1
      vmkvpp(0)=stibs(1)
      call aoib(jj)
      stib(stibs(1))=eoa
      vmkvlp(0)=stibs(1)
      call aoib(jj)
      stib(stibs(1))=eoa
      j1=stibs(1)
      call aoib(2*nvmk*ii)
      do i1=1,nvmk
      stib(vmkvpp(0)+i1)=j1
      j1=j1+ii
      stib(j1)=eoa
      stib(vmkvlp(0)+i1)=j1
      j1=j1+ii
      stib(j1)=eoa
      enddo
      if(j1.ne.stibs(1))then
      mlin(1:srec)='rvi_4'
      call mput(1,0,0,0)
      endif
      do i1=1,nrho
      nrot(i1)=0
      stib(nivd(0)+i1)=0
      enddo
      j1=0
      jj=llp
   10 continue
      jj=stib(jj)
      if(jj.ne.eoa)then
      j1=j1+1
      j2=stib(jj+1)
      stib(vval(0)+j1)=j2
      j2=j2+nivd(0)
      stib(j2)=stib(j2)+1
      goto 10
      endif
      if((j1.ne.nvert).or.(nrho.lt.3))then
      mlin(1:srec)='rvi_1'
      call mput(1,0,0,0)
      endif
      j1=0
      jj=llp
   30 continue
      jj=stib(jj)
      if(jj.ne.eoa)then
      j1=j1+1
      j2=jj+stib(vval(0)+j1)
      do i1=1,nvmk
      j2=j2+2
      stib(stib(vmkvpp(0)+i1)+j1)=stib(j2)
      stib(stib(vmkvlp(0)+i1)+j1)=stib(j2+1)
      enddo
      goto 30
      endif
      if(j1.ne.nvert)then
      mlin(1:srec)='rvi_3'
      call mput(1,0,0,0)
      endif
      kk=llp
      if(llp.gt.1)then
      if(stib(llp-1).eq.eoa)then
      kk=llp-1
      endif
      endif
      ij=stib(llp)
      if(stib(llp).ne.eoa)then
      stib(llp)=eoa
      endif
      j1=0
   50 continue
      jj=ij
      if(jj.ne.eoa)then
      j1=j1+1
      ij=stib(ij)
      j2=jj-kk+1
      call xipht(jj+2,jj+1+stib(vval(0)+j1),-j2)
      stib(vparto(0)+j1)=kk
      kk=kk+stib(vval(0)+j1)
      goto 50
      endif
      if(j1.ne.nvert)then
      mlin(1:srec)='rvi_5'
      call mput(1,0,0,0)
      endif
      kk=nivd(0)-kk-1
      call xipht(nivd(0),stibs(1),-kk)
      call aoib(-kk)
      nivd(0)=nivd(0)-kk
      vval(0)=vval(0)-kk
      vparto(0)=vparto(0)-kk
      vmkvpp(0)=vmkvpp(0)-kk
      vmkvlp(0)=vmkvlp(0)-kk
      do i1=1,nvmk
      stib(vmkvpp(0)+i1)=stib(vmkvpp(0)+i1)-kk
      stib(vmkvlp(0)+i1)=stib(vmkvlp(0)+i1)-kk
      enddo
      stib(nivd(0))=eoa
      jj=nvmk+1
      vmkr(0)=stibs(1)
      call aoib(jj)
      stib(stibs(1))=eoa
      vmks(0)=stibs(1)
      call aoib(jj)
      stib(stibs(1))=eoa
      do i1=vmks(0)+1,vmks(0)+nvmk
      stib(i1)=0
      enddo
      do i1=1,nvmk
      ii=stib(vmkvpp(0)+i1)
      jj=stib(vmkvlp(0)+i1)
      ij=0
      do i2=4,2,-2
      if(ij.eq.0)then
      do i3=1,nvert
      j1=stib(ii+i3)
      j2=stib(ii+i3)-1+stib(jj+i3)
      if(j1.gt.j2)then
      ij=1
      goto 90
      elseif(i2.eq.4)then
      if(stdz(stcb,j1,j2).eq.0)then
      goto 90
      endif
      elseif(i2.eq.2)then
      if(stdw(stcb,j1,j2).eq.0)then
      goto 90
      endif
      endif
      enddo
      ij=i2
   90 continue
      endif
      enddo
      if(ij.eq.0)then
      ij=1
      endif
      stib(vmkr(0)+i1)=ij
      enddo
      do i1=1,nvmk
      if(stib(vmkr(0)+i1).eq.1)then
      jj=stib(vmkvlp(0)+i1)
      do i2=1,nvert
      if(stib(jj+i2).gt.0)then
      goto 23
      endif
      enddo
      j1=stib(vmkp(0)+i1)
      j2=j1-1+stib(vmkl(0)+i1)
      call wput(9,j1,j2)
      elseif(stib(vmkr(0)+i1).eq.4)then
      ii=stib(vmkvpp(0)+i1)
      jj=stib(vmkvlp(0)+i1)
      do i2=1,nvert
      j1=stib(ii+i2)
      j2=j1-1+stib(jj+i2)
      if(sigz(stcb,j1,j2).ne.0)then
      goto 23
      endif
      enddo
      j1=stib(vmkp(0)+i1)
      j2=j1-1+stib(vmkl(0)+i1)
      call wput(10,j1,j2)
      endif
   23 continue
      enddo
      nvrot=0
      do 403 i1=1,nrho
      if(stib(nivd(0)+i1).gt.0)then
      rotvpo(i1)=stibs(1)
      do i2=1,nvert
      if(i1.eq.stib(vval(0)+i2))then
      j1=stib(vparto(0)+i2)
      do i3=1,i1
      xtmp(i3)=stib(j1+i3)
      enddo
      do i3=1,i1-1
      do i4=i3+1,i1
      if(xtmp(i3).gt.xtmp(i4))then
      ii=xtmp(i3)
      xtmp(i3)=xtmp(i4)
      xtmp(i4)=ii
      endif
      enddo
      enddo
      goto 830
  805 continue
      do i3=i1-1,1,-1
      if(xtmp(i3).lt.xtmp(i3+1))then
      j1=i3
      goto 810
      endif
      enddo
      goto 401
  810 continue
      do i3=j1+2,i1
      if(xtmp(i3).le.xtmp(j1))then
      j2=i3-1
      goto 820
      endif
      enddo
      j2=i1
  820 continue
      ii=xtmp(j1)
      xtmp(j1)=xtmp(j2)
      xtmp(j2)=ii
      jj=(i1-j1)/2
      j1=j1+1
      do i3=0,jj-1
      ii=xtmp(j1+i3)
      xtmp(j1+i3)=xtmp(i1-i3)
      xtmp(i1-i3)=ii
      enddo
  830 continue
      j1=stibs(1)+1
      call aoib(i1+1)
      stib(j1)=i2
      do i3=1,i1
      stib(j1+i3)=xtmp(i3)
      enddo
      nvrot=nvrot+1
      nrot(i1)=nrot(i1)+1
      goto 805
      endif
  401 continue
      enddo
      j1=stibs(1)+1
      call aoib(i1+1)
      stib(j1)=0
      do i2=1,i1
      stib(j1+i2)=nphi+1
      enddo
      endif
  403 continue
      call aoib(1)
      stib(stibs(1))=eoa
      do 790 i1=1,nrho
      if((stib(nivd(0)+i1).gt.0).and.(nrot(i1).gt.1))then
      h2=nrot(i1)
      h1=(h2/2)+1
  720 continue
      if(h1.gt.1)then
      h1=h1-1
      ii=rotvpo(i1)+(i1+1)*(h1-1)
      do i2=1,i1+1
      xtmp(i2)=stib(ii+i2)
      enddo
      else
      ii=rotvpo(i1)+(i1+1)*(h2-1)
      do i2=1,i1+1
      xtmp(i2)=stib(ii+i2)
      enddo
      do i2=1,i1+1
      stib(ii+i2)=stib(rotvpo(i1)+i2)
      enddo
      h2=h2-1
      if(h2.eq.1)then
      do i2=1,i1+1
      stib(rotvpo(i1)+i2)=xtmp(i2)
      enddo
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
      jj=rotvpo(i1)+(i1+1)*(hj-1)
      jk=jj+i1+1
      j1=1
  745 continue
      if(j1.le.i1)then
      j1=j1+1
      if(stib(jj+j1).lt.stib(jk+j1))then
      kk=-1
      elseif(stib(jj+j1).gt.stib(jk+j1))then
      kk=1
      else
      goto 745
      endif
      endif
      if(kk.eq.-1)then
      hj=hj+1
      elseif(kk.eq.0)then
      mflag(4)=1
      endif
      endif
      kk=0
      jj=rotvpo(i1)+(i1+1)*(hj-1)
      j1=1
  755 continue
      if(j1.le.i1)then
      j1=j1+1
      if(xtmp(j1).lt.stib(jj+j1))then
      kk=-1
      elseif(xtmp(j1).gt.stib(jj+j1))then
      kk=1
      else
      goto 755
      endif
      endif
      if(kk.eq.-1)then
      ii=rotvpo(i1)+(i1+1)*(hi-1)
      do i2=1,i1+1
      stib(ii+i2)=stib(jj+i2)
      enddo
      goto 740
      elseif(kk.eq.0)then
      mflag(4)=1
      endif
      endif
      ii=rotvpo(i1)+(i1+1)*(hi-1)
      do i2=1,i1+1
      stib(ii+i2)=xtmp(i2)
      enddo
      goto 720
      endif
  790 continue
      if(mflag(4).ne.0)then
      call wput(7,0,0)
      endif
      do i1=1,nrho
      if(nrot(i1).gt.0)then
      ii=rotvpo(i1)
      do i2=1,nrot(i1)
      j1=ii+2
      j2=j1+i1+1
      jj=i1
   84 continue
      if(jj.gt.0)then
      kk=stib(j2)-stib(j1)
      if(kk.eq.0)then
      jj=jj-1
      j1=j1+1
      j2=j2+1
      goto 84
      elseif(kk.lt.0)then
      mlin(1:srec)='rvi_2'
      call mput(1,0,0,0)
      endif
      endif
      ii=ii+i1+1
      enddo
      endif
      enddo
      ii=nrho+1
      dpntro(0)=stibs(1)
      call aoib(ii)
      stib(stibs(1))=eoa
      do i1=1,nrho
      if(nrot(i1).gt.0)then
      stib(dpntro(0)+i1)=stibs(1)
      call aoib(nphi+1)
      stib(stibs(1))=nap
      ii=rotvpo(i1)+2
      jj=0
      do i2=1,nrot(i1)+1
      if(stib(ii).ne.jj)then
      if(stib(ii).le.nphi)then
      kk=stib(ii)
      else
      kk=nphi
      endif
      do i3=jj+1,kk
      stib(stib(dpntro(0)+i1)+i3)=ii-1
      enddo
      jj=stib(ii)
      endif
      ii=ii+i1+1
      enddo
      else
      stib(dpntro(0)+i1)=nap
      endif
      enddo
      stib(stibs(1))=eoa
      bloka=stibs(1)
      ii=nphi
      call vaoib(ii)
      do i1=blok(0)+1,blok(0)+nphi
      stib(i1)=0
      enddo
      nblok=0
      j1=0
      j2=0
   31 continue
      if(j1.lt.nphi)then
      nblok=nblok+1
      j3=nblok
   32 continue
      if(stib(blok(0)+j3).gt.0)then
      j3=j3+1
      goto 32
      endif
      j1=j1+1
      j2=j1
      stib(bloka+j1)=j3
      stib(blok(0)+j3)=nblok
      ii=stib(link(0)+j3)
      jj=blok(0)+ii
      if(stib(jj).eq.0)then
      j1=j1+1
      stib(bloka+j1)=ii
      stib(jj)=nblok
      endif
   33 continue
      if(j2.le.j1)then
      jj=stib(bloka+j2)
      do i1=2,nrho
      if(stib(dpntro(0)+i1).gt.0)then
      j3=stib(stib(dpntro(0)+i1)+jj)
   36 continue
      if(stib(j3+1).eq.jj)then
      ii=stib(j3+2)
      if(stib(blok(0)+ii).eq.0)then
      j1=j1+1
      stib(bloka+j1)=ii
      stib(blok(0)+ii)=nblok
      endif
      ii=stib(link(0)+ii)
      if(stib(blok(0)+ii).eq.0)then
      j1=j1+1
      stib(bloka+j1)=ii
      stib(blok(0)+ii)=nblok
      endif
      j3=j3+i1+1
      goto 36
      endif
      endif
      enddo
      j2=j2+1
      goto 33
      endif
      goto 31
      endif
      if(nblok.gt.1)then
      call wput(1,0,0)
      endif
      j1=0
   34 continue
      if(j1.lt.nphi)then
      j1=j1+1
      do i1=1,nrho
      if(stib(nivd(0)+i1).gt.0)then
      if(stib(stib(stib(dpntro(0)+i1)+j1)+1).eq.j1)then
      goto 34
      endif
      endif
      enddo
      call wput(2,0,0)
      endif
      return
      end
