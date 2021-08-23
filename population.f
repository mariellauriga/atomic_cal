      program population

c***********************************************************************
c     This program reads oscillator strength Aij and collision strength
c     gamma(i,j) from data files for a specific temperature(tem). Then  
c     for a given electrons density(de), it calculates the population
c     of excited energy levels relative to the ground level population.
c     Also, the emissivity ratios are found out.
cc   
cc    ---------         AKP   12.01.92     ___________________
cc
cc    The aij(i,j) are read in file aij.in such that i,j are in observed
cc    order. gamma(i,j) however are in stgfj order that is reset to
cc    the observed order using the correspondence defined in file
cc    ei.dat.
cc
c***********************************************************************

c
c     alpha and U are the power law order and ionization parameter
c     tt--radiation temperature, w--dilution coefficient
c***********************************************************************
       implicit double precision (a-h,o-z)
      character*4 dete
      integer nt,ne
      dimension tem(290),de(290)


      COMMON /BLK1/gamma(290,290),ei(290),gi(290),indo(290),inda(290)
      COMMON /QIJ/qij(290,290)
      COMMON /NQIJ/rqij(290,290)
      COMMON /AIJ/aij(290,290),wave(290,290)
      COMMON /C/cij(290,290)
      COMMON /CT/ctail(290,290) 
      COMMON /P/para(290)
      COMMON /IND/index(290) 

c**********************************************************************
c     set initial temperature 'tem' and number of energy level 'ni'
c**********************************************************************

       open(unit=11,file='temp.in',status='unknown')
       read(11,*) nt
       read(11,*) (tem(i),i=1,nt) 
       read(11,*) ni

       open(unit=12,file='den.in', status='unknown')
       read(12,*) ne  
       do i=1,ne
       read(12,*) d
       de(i)=d
c       write(6,*)i, de(i)
       enddo
       close(unit=12)

c********************************************************************
c     readin the energy level in the array ei()and gi() and ind()
c********************************************************************
       call getei(ni)
       call reada(ni)
       close(unit=11) 
c**********************************************************************
c     Do loop for each temperature and electron density 'de'
c**********************************************************************

       do 5 j=1,nt

       do 100 i=1,ne

c**********************************************************************
c     setup, read data , prepare for finding CIJ
c**********************************************************************
          call readg(ni,tem(j))   
          call gqij(tem(j),ni)

c*********************************************************************
c     Get coefficient matrix and constants for linear equation
c*********************************************************************
          call nuqij(ni,de(i))
          call getcij(ni)
          call gctail(ni)
          call gpara(ni)

c*********************************************************************
c     Solve equation for population and flux ratios.
c*********************************************************************
          call cmpctail(ni-1)
          call slu(ni-1,tem(j),de(i))
          call getrr(tem(j),de(i))

 100   continue  
 5     continue
       stop
       end

c-------------bellow are subrutine-------------------

c************************************************************************
c     The program GETEI reads the energy values from a data file 'ei.dat' 
c     and store them in 1-D COMMON block EI.
c
c     called by: main
c     call subtoutines: none
c************************************************************************

      subroutine getei(n)

      implicit double precision (a-h,o-z)
       COMMON /BLK1/gamma(290,290),ei(290),gi(290),indo(290),inda(290)

      open(unit=15,file='d415t006',status='unknown')
      open(unit=16,file='ei.out',status='unknown')

      do 100 i=1,n
         read(15,*) i1,i2,i3,i4,i5,i6,en,gw
         indo(i)=i    
         inda(i)=i
         gi(i)=gw 
         ei(i)=en
c         ei(i)=en/13.6 
         write(16,1000)i,ei(i),gi(i)
 100     continue

c     do i=1,n
c        write(16,1000) i,ei(i),gi(i)
c     enddo

      if(ei(1).ne.0.) then
         print *,' e1 .ne. 0', ei(1)
      endif

 1000    format(i5,2(3x,e14.6))

         close(unit=15)
         close(unit=16)

         return
         end



c**********************************************************************
c     This subroutine reads oscilator strength Aij from a data file
c     'aij.in' and store them in a 2-D COMMON block AIJ.
c
c     called by: main
c     call subroutine: none
c**********************************************************************

      subroutine reada(n)

       implicit double precision (a-h,o-z)

      COMMON /BLK1/gamma(290,290),ei(290),gi(290),indo(290),inda(290)
      COMMON /AIJ/aij(290,290),wave(290,290)

      do 400 i=1,n      
         do 450 j=1,n
            aij(i,j)=0.
 450        continue
 400        continue

      open(unit=15,file='d415t050',status='old')
      open(unit=16,file='aij.out',status='unknown')
c      open(unit=99,file='paver.out',status='unknown') 
c      read(15,*)
   50 read(15,*,end=500) i1,i2,i3,i4,i5,i6,ww,r,aa,il,iu

      if (iu.gt.n) goto 50
      aij(iu,il)=aa
      wave(iu,il)=ww
      write(99,*)iu,il,aa,ww
      
      goto 50
 500  continue
      close(15)

         do 300 j=1,n-1
           do i=j+1,n 
               write(16,1201) i,j,aij(i,j)
           enddo
  300    continue 

c1100 format(i2,i2,e9.3)      
 1200 format(9(e8.2,1x))
 1201 format(2i4,3x,e9.3)

  301    continue
      close(unit=15)
      close(unit=16)
      close(unit=17)
      close(unit=18)

      return
      end


c**********************************************************************
c     The subroutine READG reads gamma(i,j) from a data file 'gamma.in' 
c     and store them in a 2-D COMMON block gamma.
c
c     Called by: main
c     Call subroutines:none
c**********************************************************************

      subroutine readg(n,ti)
      implicit double precision(a-h,o-z)

       COMMON /BLK1/gamma(290,290),ei(290),gi(290),indo(290),inda(290)
      COMMON /AIJ/aij(290,290),wave(290,290)
      real temp(30)
      real gam(30)       

      do 500 i=1,n
         do 550 j=1,n
            gamma(i,j)=0
 550        continue
 500        continue
      do i=2,n
       gamma(1,i)=.001
      enddo 

      open(unit=15,file='upsilon.dat',status='unknown')
      open(unit=16,file='upsilon.out',status='unknown')
      open(unit=24,status='unknown')

      read(15,*)ntemp
      read(15,*)(temp(k),k=1,ntemp)
      it=0
      do i=1,ntemp-1
       if (temp(i).le.ti .and. temp(i+1).gt.ti) it=i
      enddo
 50   read(15,*,end=200)i,j,(gam(k),k=1,ntemp)
      rm=(gam(it+1)-gam(it))/(temp(it+1)-temp(it))
      i1=i
      i2=j
      if (i2.lt.i1) then
       ii=i2
       i2=i1
       i1=ii
      endif
      gamma(i1,i2)=gam(it)+rm*(ti-temp(it))
      write(16,*)i,i1,j,i2,gamma(i1,i2)
      goto 50
 200  continue
 

      write(96,*)ti
      do 410 i=1,n-1
         write(96,2100) (gamma(i,j),j=i+1,n)
 410     continue
            
 1105 format(2(f9.6,2x))
 1100 format(1x,f6.0,1x,9(f7.4,1x))      
 1200 format(8x,9(f7.4,1x))
 1300 format(8x,1(f7.4,1x))
 1400 format(8x,2(f7.4,1x))
 1500 format(8x,3(f7.4,1x))
 1600 format(8x,4(f7.4,1x))
 1700 format(8x,5(f7.4,1x))
 1800 format(8x,6(f7.4,1x))
 1900 format(8x,7(f7.4,1x))
 2000 format(8x,8(f7.4,1x))
 2100 format(8(1pe10.4,1x))

      close(unit=14)
      close(unit=15)
      close(unit=16)

      return
      end



c*********************************************************************
c     This subroutine get matrix QIJ from data of gamma, ei and gi and
c     store in 2-D COMMON QIJ.
c
c     The equation of calculation is :
c     Qij(i,j)=(8.63E-06)*gamma(i,j)*exp(-(ei(j)-ei(i))/kT)/sqrt(T)
c     for j>i
c     Qij(i,j)=Qij(j,i)*gi(j)/gi(i) for i>j
c
c     called by:main
c     call subroutine:none
c*********************************************************************

         subroutine gqij(ti,n)

         implicit double precision (a-h,o-z)    
         COMMON /BLK1/gamma(290,290),ei(290),gi(290),indo(290),inda(290)
         COMMON /QIJ/qij(290,290)

         open(unit=16, file='qij.out', status='unknown')

         do 10 i=1,n
            do 20 j=1,n
               qij(i,j)=0
 20            continue
 10            continue

         do 100 i=1,n-1
            do 150 j=i+1,n
               deltae=ei(j)-ei(i)
               para=157885.*deltae/ti
               if(deltae.gt.0.) then
                  pexp=exp(-para)
                  psq=sqrt(ti)*gi(i)
                  pij=gamma(i,j)*(8.63e-06)/psq
                  qij(i,j)=pij*pexp
c     qij(j,i)=qij(i,j)*gi(i)/gi(j)/pexp
                  qij(j,i)=pij*gi(i)/gi(j)
               else
                  pexp=exp(para)
                  psq=sqrt(ti)*gi(j)
                  pij=gamma(i,j)*(8.63e-06)/psq
                  qij(j,i)=pij*pexp
                  qij(i,j)=pij*gi(j)/gi(i)
c     qij(i,j)=qij(j,i)*gi(j)/gi(i)/pexp
               endif
 150           continue           
 100           continue


         do 200 i=1,n-1
            do 250 j=i+1,n
              qlow=gi(i)*qij(i,j)
              qhigh=gi(j)*qij(j,i) 
         if(abs(qlow-qhigh).gt.1.e-04) then
             write(6,6000) i,j,qlow,qhigh
 6000 format(' i,j,qlow,qhigh =', 2i5,5x,2e12.4)
         stop
         endif
 250           continue
 200           continue


 1100       format(9(D8.3,1x))

            close(unit=16)

            return
            end
            


c*********************************************************************
c     This subroutine multiply the qij elements by electron density
c     ne, and store the new elements in 2-D COMMON block NQIJ.
c
c     called by: main
c     call subroutines: none
c*********************************************************************

      subroutine nuqij(n,x)

      implicit double precision(a-h,o-z)
      COMMON /QIJ/qij(290,290)
      COMMON /NQIJ/rqij(290,290)

      open(unit=16,file='nqij.out',status='unknown')

      do 100 i=1,n
         do 150 j=1,n
            rqij(i,j)=x*qij(i,j)
 150        continue
 100        continue

      do 200 i=1,n
         write(16,1100) (rqij(i,j),j=1,n)
 200     continue

 1100    format(9(e8.3,1x))

         close(unit=16)

         return
         end

c**********************************************************************
c     This subroutine get matrix Cij by Cij=ne*Qij+Aij+u*Bij. The result is
c     stored in 2-D COMMON block CIJ.
c
c     called by: main
c     call subroutine: none
c**********************************************************************

         subroutine getcij(n)
         implicit double precision(a-h,o-z)

         COMMON /NQIJ/rqij(290,290)
         COMMON /AIJ/aij(290,290),wave(290,290)
         COMMON /C/cij(290,290)
c this is the added block         
         COMMON/UBIJ/ubij(290,290) 

         open(unit=16,file='cij.out',status='unknown')

         do 100 i=1,n
            sumcu=0.
            sumcd=0.
            sumc=0. 
            do 150 j=1,n
c               cij(i,j)=aij(i,j)+rqij(i,j)+ubij(i,j)
             cij(i,j)=aij(i,j)+rqij(i,j)
             if (j.ne.i) sumc=sumc+cij(i,j)
             if (j.lt.i) sumcd=sumcd+cij(i,j)
             if (j.gt.i) sumcu=sumcu+cij(i,j)
 150           continue
            if (sumc.eq.0.) print*, i
c           if (i.gt.1)write(16,1101) i,sumcu,sumcd,
c    #          cij(1,i)/(sumcd+sumcu)
 100           continue

               do 200 i=1,n
                write(16,*) i 
                  write(16,1100) (cij(i,j),j=1,n)
 200           continue
               
 1100          format(9(D8.3,1x))
 1101          format(i3,3(2x,d9.4))
               
               close(unit=16)
               
               return
               end

c*********************************************************************
c     This subroutine get Ctail matrix by rearrange:
c     Ctail(i,i)=Sum(Cij(i,j)) i.ne.j
c     Ctail(i,j)=-CIJ(j+1,i+1) for i.ne.j
c     Ctail will be used as coefficient matrix for solving population.
c     It is stored in a 2-D COMMON block CT/ctail.
c
c     called by:main
c     call subroutines: none
c*********************************************************************

            subroutine gctail(n)

            implicit double precision(a-h,o-z)
            COMMON /C/cij(290,290) 
            COMMON /CT/ctail(290,290)

            open(unit=26,file='ctail.out',status='unknown')

            do 100 i=1,n-1
               do 150 j=1,n-1
                  if(j.eq.i) then
                     ctail(i,j)=cij(i+1,1)
                     do 10 k=2,n
                        ctail(i,j)=ctail(i,j)+cij(i+1,k)
 10                     continue
                     ctail(i,j)=ctail(i,j)-cij(i+1,i+1)
                  else
                     ctail(i,j)=-cij(j+1,i+1)
                  end if
 150           continue
 100        continue
   
            do 200 i=1,n-1
               write(26,1100) (ctail(i,j), j=1,n-1)
 200        continue

 1100       format(9(D8.2,1x))

            return
            end

c*********************************************************************
c     This subroutine get constant term of linear equation by:
c     para(i)=CIJ(1,i+1). It is stored in 1-D COMMOn P/para.
c     
c     called by: main
c     call subroutines: none
c********************************************************************

            subroutine gpara(n)

            implicit double precision(a-h,o-z)
            COMMON /C/cij(290,290)
            COMMON /P/para(290)

            open(unit=36,file='para.out',status='unknown')

            do 100 i=1,n-1
               para(i)=cij(1,i+1)
 100        continue

            do 200 i=1,n-1
               write(36,1100) para(i)
 200        continue
 1100       format(D10.3)
            close(unit=36)

            return
            end

c********************************************************************
c     This subroutine composite the coefficient matrix Ctail in a 
c     upper teiangle and a lower triangle matrix, which is easily to
c     find the solutions of linear equation. The new matrix will
c     replace the old one to store in 2-D COMMON block Ctail.
c
c     For upper triangle matrix(j>i), new
c     Ctail(i,j)=Ctail(i,j)-Sum(k=1,i-1)Ctail(i,k)*Ctail(k,j)
c
c     For lower triangle matrix with 1 in diagonal(i>j), new
c     Ctail(i,j)=1/Ctail(j,j)*(Ctail(i,j)-Sum(k=1,i-1)Ctail(i,k)*Ctail(k,j)
c
c     called by: main
c     call subroutines: none
c********************************************************************

            subroutine cmpctail(n)

            implicit double precision(a-h,o-z)
            parameter (tiny=1.0e-30)
            dimension vv(290)
            COMMON /CT/ctail(290,290) 
            COMMON /IND/index(290) 

            d=1
            do 40 i=1,n
               aamax=0.
               do 30 j=1,n
                  if(abs(ctail(i,j)).gt.aamax) aamax=abs(ctail(i,j))
 30            continue
               if(aamax.eq.0.) pause 'Singular matrix.'
               vv(i)=1/aamax
 40         continue

            do 100 j=1,n
               if(j.gt.1) then
                  do 120 i=1,j-1
                     sum=ctail(i,j)
                     if(i.gt.1) then
                        do 140 k=1,i-1
                           sum=sum-ctail(i,k)*ctail(k,j)
 140                    continue
                        ctail(i,j)=sum
                     end if
 120              continue
               end if
               aamax=0.
               do 200 i=j,n
                  sum=ctail(i,j)
                  if (j.gt.1) then
                     do 220 k=1,j-1
                        sum=sum-ctail(i,k)*ctail(k,j)
 220                 continue
                     ctail(i,j)=sum
                  end if
                  dum=vv(i)*abs(sum)
                  if(dum.ge.aamax) then
                     imax=i
                     aamax=dum
                  end if
 200           continue
               if(j.ne.imax) then
                  do 240 k=1,n
                     dum=ctail(imax,k)
                     ctail(imax,k)=ctail(j,k)
                     ctail(j,k)=dum
 240              continue
                  d=-d
                  vv(imax)=vv(j)
               end if
               index(j)=imax
               if(j.ne.n) then
                  if(ctail(j,j).eq.0.) ctail(j,j)=tiny
                  dum=1./ctail(j,j)
                  do 260 i=j+1,n
                     ctail(i,j)=ctail(i,j)*dum
 260              continue
               end if
 100        continue
            if(ctail(n,n).eq.0) ctail(n,n)=tiny

            open(unit=60,file='lud.out',status='unknown')

            do 700 i=1,n
               write(60,1100) (ctail(i,j), j=1,n)
 700        continue

 1100       format(8(D10.2,1x))

            close(unit=60)

            return
            end

c********************************************************************
c     This subroutine slove the linear equation with coefficient
c     matrix are upper and lower triangle matrixes. The solution will
c     replace the old constant terms and stored in 1-D COMMON block
c     P/para. The solution are electron populations of each exicting
c     states relatively to the ground states.
c
c     called by: main
c     call subroutines: none
c********************************************************************

            subroutine slu(n,ti,x1)

            implicit double precision(a-h,o-z)
            character*9 id(290)
            common /IDS/id

            COMMON /C/cij(290,290) 
            COMMON /CT/ctail(290,290) 
            COMMON /P/para(290) 
            COMMON /IND/index(290)
        COMMON /BLK1/gamma(290,290),ei(290),gi(290),indo(290),inda(290)
            COMMON /AIJ/aij(290,290),wave(290,290) 
         COMMON /NQIJ/rqij(290,290)
         COMMON/UBIJ/ubij(290,290) 

            ii=0
            do 100 i=1,n
               ll=index(i)
               sum=para(ll)
               para(ll)=para(i)
               if(ii.ne.0) then
                  do 120 j=ii,i-1
                     sum=sum-ctail(i,j)*para(j)
 120              continue
               else if(sum.ne.0.) then
                  ii=i
               end if
               para(i)=sum
 100        continue
            do 200 i=n,1,-1
               sum=para(i)
               if(i.lt.n) then
                  do 220 j=i+1,n
                     sum=sum-ctail(i,j)*para(j)
 220              continue
               end if
               para(i)=sum/ctail(i,i)
 200        continue
            
           open(unit=41,file='pop2.out',status='unknown')
           open(unit=42,file='pop3.out',status='unknown')
           open(unit=43,file='pop4.out',status='unknown')
           open(unit=44,file='pop5.out',status='unknown')

           open(unit=46,file='pop.out',status='unknown')
            open(unit=57,file='line.out',status='unknown')

c          write(46,*) '----------'
           total=1.
           do 699 i=1,n
              total=total+para(i)
 699       continue

c           xnormal=para(3)*aij(4,1)*(ei(4)-ei(1))
           xnormal=1.
             
            do 700 i=1,n
c           do 700 i=1,52
               do 777 j=1,i
                 enrx=ei(i+1)-ei(j)
                 strength=para(i)*enrx*aij(i+1,j)
c                 if (strength/xnormal.gt.0.0001) then
                if(wave(i+1,j).gt.2.e3.and.wave(i+1,j).lt.11.e3) then
                 write(57,1102)i+1,j,wave(i+1,j),strength/xnormal,
     *                 id(i+1),id(j)
c                 write(57,1102)i+1,j,wave(i+1,j),strength/1.0,
c     #                 id(i+1),id(j)

                  endif
c                 endif
777            continue
700         continue

               write(46,1103)1,1.,1.
           do 399 i=1,n
            boltz=gi(i+1)/gi(1)*exp(-ei(i+1)*157889.32/ti)
              write(46,1103)i+1,para(i)/total,boltz
 399       continue

c          write(46,1100)23,ti,x1,para(23)/total
           write(46,*)ti,x1,xnormal/total*2.1799e-18
           
 1100       format(1x,i3,2x,f10.0, 1x, d10.3,1x,D10.3)
 1101       format(f10.0, 1x, d10.3,1x,D10.3)
 1102       format(i3,2x,i3,2x,d12.5,2x,d10.3,2(2x,a9))
 1103       format(i3,2(2x,d10.3))
 1104       format(2i5,3(2x,d10.3))

c            close(unit=46)

            return
            end

c*******************************************************************
c     This subroutine get the flux ratio of two transitions by:
c     Bigg(ij,nm)=(Ni*Aij*(ei(j)-ei(i)))/(Nn*Anm*(ei(m)-ei(n)))
c
c     called by: main
c     call subroutines: none
c*******************************************************************
            subroutine getrr(ti,x2)
            implicit double precision(a-h,o-z)

        COMMON /UBIJ/ubij(290,290) 
        COMMON /BLK1/gamma(290,290),ei(290),gi(290),indo(290),inda(290)
            COMMON /P/para(290) 
            COMMON /AIJ/aij(290,290),wave(290,290)

            open(unit=80, file='nratios.in',status='old')
            open(unit=81, file='nratios.out',status='unknown')

            read(80,*) nr
            do 100 i=1,nr
               read(80,*) nin, nout
c               write(6,*) nin, nout
               sg1=0.
               do 150 j=1, nin
                  read(80, *) nup, ndown
c                  write(6,*)nup, ndown
             dele=(ei(nup)-ei(ndown))
            sg1=(aij(nup,ndown))*dele*para(nup-1)+sg1
 150           continue

               sg2=0.
               do 160 j=1, nout
                  read(80, *) nup, ndown
             dele=(ei(nup)-ei(ndown))
             sg2=(aij(nup,ndown))*dele*para(nup-1)+sg2
 160           continue


            bigg1=sg2/sg1

c            x1=dlog10(x2)
c            bigg2=dlog10(bigg1)
            
c            write(81,1100)i, ti, x2, w, bigg1
             write(81,1100)i, ti, x2, bigg1
 100        continue
            close(unit=80)
c 1100       format(i3,1x,f10.0,2(1x,d10.3),1x,f15.6)
1100         format(i3,1x,f10.0,1x,d10.3,1x,f15.6)
c2200       format(i2)
c2201       format(i2,1x,i2)


            return
            end