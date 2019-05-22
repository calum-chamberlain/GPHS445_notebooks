c Sample main driver program for using the focal mechanism subroutines.  
c Uses P polarities and S/P amplitude ratios.
c
c   program modified for seisan by j.havskov, oct 2010
c   program modified for VUW GPHS445 class by Calum Chamberlain, May 2019
c
c   - all include files included so no include files needed
c   - input to use focmec format, assume free surface correction included so
c     only z used for P
c   - only use sh, sv commented out in calculations
c   - variables declare originally still there
c   - only input of one set of observations so error calculation will give
c     too small value
c   - maximum solutions out fixed to 1000
c   - How it works:
c      First search within critria given to find minimum number of pol misfits
c         and minimum average amp error
c      If no solutions within critria (min number pol, min av amp. error) + extra
c         pol and amp error, then increase limits to min error found+ extra
c      This does not ensure a solution since both minima were found independently
c         so if no solutions, search for minimum amplitude error within given
c         polarity error and increase amp error to this + extra, in this way
c         solutions will always be found. if too fe solutions is found, the errro 
c         limits of both pol and amplitude can be increased.
c      Extra for pol is half the given max polarity error, at least 1
c      Extra for amp is half the given max av amp error
c
c   changes
c
c nov 6  2010 jh: run focmec first to prepare
c nov 12 2010 jh: rename program to hash_foc, hash is used by linux
c jan 11 2011 jh: general output unit 13 has been diabled, problem with format
c                 in gfortran. most of it undefined or written elsewhere
c jan 18 2011 jh: hash.out to hash_foc.out,change output 
c                 format a bit in hash_foc.out
c jan 31 2011 jh: format change
c feb 23 2011 jh: change hash_foc to hash_seisan
c may 20 2019 cc: Change to inline input
c
c      include 'param.inc'
c npick0 = maximum number of picks per event
c nmc0 = maximum number of trials of location/take-off angles
c nmax0 = maximum number of acceptable mechanisms output
      parameter (npick0=1000,nmc0=500,nmax0=1000) 

c      include 'rot.inc'
c dang0 : minimum grid spacing (degrees)
c ncoor : number of test mechanisms

c choose one, based on your computer's memory and speed
c this is grid resolution
c
c      parameter (dang0=10.0,ncoor=4032) 
c      parameter (dang0=7.0,ncoor=10426) 
c      parameter (dang0=5.0,ncoor=31032) 
c      parameter (dang0=3.0,ncoor=141180) 
       parameter (dang0=2.0,ncoor=472410) 


       character*80 text

c variables for storing earthquake input information  
      integer icusp,icusp2          ! event ID
      real qlat,qlon,qdep           ! location
      real qmag                     ! magnitude
      integer iyr,imon,idy,ihr,imn  ! origin time, year, month, day, hour, minute
      real qsec                     ! origin time, seconds
      real seh, sez                 ! location error, horizontal & vertical 
      real rms                      ! residual travel time error 
      real terr                     ! origin time error 
      character*1 evtype            ! event type
      character*1 magtype           ! magnitude type
      character*1 locqual           ! location quality
      character*1 cns,cew           ! north/south and east/west codes
c
c variables for polarity information, input to HASH subroutines
      character*4 sname(npick0)                        ! station name
      character*3 scomp(npick0)                        ! station component
      character*2 snet(npick0)                         ! network code
      character*1 pickpol,pickonset                    ! polarity pick : U, u, +, D, d, or - ; onset, I or E
      integer p_pol(npick0),spol                       ! polarity pick (+/-1), and reversal (+/-1)
      real sp_ratio(npick0),spin                       ! S/P ratio (log10)
      real p_azi_mc(npick0,nmc0),p_the_mc(npick0,nmc0) ! azimuth and takeoff angle for each trial
      integer index(nmc0)                              ! index into velocity models, for each trial
      real qdep2(nmc0)                                 ! new depth, for each trail
      integer nmc                                      ! number of trails with different azimuth and take-off angles
      integer npol,nppl,nspr                           ! number of observations, P polarities, and S/P ratios
c
c variables for set of acceptable mechanisms, output of HASH subroutines
      integer nout2                                    ! number of acceptable mechanisms returned
      integer nmult                                    ! number of solutions (1 if no mulitples)
      real str_avg(5),dip_avg(5),rak_avg(5)            ! solution(s)
      real f1norm(3,nmax0),f2norm(3,nmax0)             ! normal vectors to the fault planes
      real strike2(nmax0),dip2(nmax0),rake2(nmax0)     ! strike, dip and rake
      real var_est(2,5),var_avg(5)                     ! variance of each plane, average
      real mfrac(5),stdr(5),mavg(5)                    ! fraction misfit polarities, station distribution       
      real prob(5)                                     ! probability true mechanism is "close" to preferred solution(s)
      character*1 qual(5),mflag                        ! solution quality rating, multiple flag
c
c control parameters
      integer npolmin                                  ! minimum number of observations
      real delmax                                      ! maximum station distance
      real dang,dang2                                  ! grid angle to sample focal sphere
      integer maxout                                   ! max number of acceptable mechanisms output
      real badfrac                                     ! assumed rate of polarity error (fraction)
      real cangle                                      ! definition of "close" == 45 degrees
      real ratmin                                      ! minimum allowed signal to noise ratio
      real qbadfac                                     ! assumed noise in amplitude ratios, log10 (0.3 for factor of 2)
c
c file names
c
      character*100 outfile1,corfile,fpfile
      character*100 stfile,plfile,ampfile
      
      degrad=180./3.1415927
      rad=1./degrad
      
c
      open (13,file='hash_seisan.out',status='unknown')

c   get arguments from command line
      if (command_argument_count() .ne. 5)then
        print *,'Grid angle for focal mech. search, enter for def 2' 
 
        read(5,'(a)') text
        if(text.eq.' ')then
           dang2=2.0
        else
           read(text,*) dang2
           if(dang2.lt.2) then
             write(6,*)' Angle set to 2'
             dang2=2
           endif
        endif
c
c   maxim number of output solutions now hardwired to 1000
c
   
        maxout=1000


         write(6,*) ' Max number of polarity errors'
         read(5,*) nmismax


         write(6,*) ' Max average error in amp rat, log10, def 0.2'
         read(5,'(a)') text
         if(text.eq.' ') then
            qmismax=0.2
         else             
            read(text,*) qmismax
         endif


        print *,' Enter angle for computing mechanisms probability'
        print *,'    default is 60'
        read(5,'(a)') text
        if(text.eq.' ') then
           cangle=60
        else
           read(text,*) cangle
        endif
       

        print *,' Enter probability threshold for multiples, def is 0.1'
        read(5,'(a)') text
        if(text.eq.' ') then
           prob_max=0.1
        else
           read(text,*) prob_max
        endif       
      else
        call get_command_argument(1, text)
        read(text,*) dang2
        ! print *,dang2
        call get_command_argument(2, text)
        read(text,*) nmismax
        ! print *,nmismax
        call get_command_argument(3, text)
        read(text,*) qmismax
        ! print *,qmismax
        call get_command_argument(4, text)
        read(text,*) cangle
        ! print *,cangle
        call get_command_argument(5, text)
        read(text,*) prob_max
        ! print *,prob_max
      endif


      terr=-9                  ! set parameters not given in input file
      rms=-9
      nppick=-9
      nspick=-9
      evtype='L'
      magtype='X'
      locqual='X'


      k=1


      icusp=1

c
c   start focmec to prepare input for hash
c
      call systemc('focmec p',8)

c
c   read focmec input file
c
       k=1
       npol=0
       nspr=0
       nppl=0
       open(12,file='focmec.dat',status='old')
       read(12,'(a)') pickpol             ! read one line
 221   continue
       read(12,'(6x,f6.1,f8.1,a1,f8.3)',end=222)
     * p_azi_mc(k,1),p_the_mc(k,1),pickpol,sp_ratio(k)
       p_pol(k)=0

       p_the_mc(k,1)=180.0-p_the_mc(k,1)
       if(pickpol.eq.'C') p_pol(k)=1
       if(pickpol.eq.'D') p_pol(k)=-1
       if(p_pol(k).ne.0) nppl=nppl+1
c
c   only use ratio sh to p identified as H in input file, discard sv
c
       if(sp_ratio(k).ne.0.0) then
          if(pickpol.eq.'H') then
             nspr=nspr+1
          else
             k=k-1
          endif
       endif
       k=k+1
       goto 221
 222   continue
       npol=k-1
       close(12)

c
cc print data
c      do k=1,npol
c        print *,k,p_azi_mc(k,1),p_the_mc(k,1),p_pol(k),sp_ratio(k)
c      end do


      if (nppl.lt.1) then
        print *,'*** warning - no p-wave polarity data for event'
           stop
      else
        write(6,'(a,i5)')
     * ' Number of polarities is                      : ',nppl
      endif   
      if (nspr.lt.1) then
        print *,'*** warning - no s/p amplitude ratios for event'
      else
        write(6,'(a,i5)')
     *' Number of amplitude ratios is                : ',nspr
      endif
      
       nmc=1    !  no trial set


       nextra=nmismax/2
       if(nextra.eq.0) nextra=1

c
c  assume total amplitude error to be maximum average error times 
c  number of amplitude observertions
c
       qmismax=qmismax*nspr
       qextra=qmismax/2.0

  

c
c jh: added nspr to call in order to calculate average amplitude error
c

      call FOCALAMP_MC(p_azi_mc,p_the_mc,sp_ratio,p_pol,npol,nmc,
     &    dang2,nmax0,nextra,nmismax,qextra,qmismax,nf2,strike2,dip2,
     &    rake2,f1norm,f2norm,nspr)

      nout2=min(nmax0,nf2)  ! number mechs returned from sub
      nout1=min(maxout,nf2)  ! number mechs to return

      write(6,'(a,i5)') 
     *' Number of solutions found                      ',nf2

c      write(6,*) 'Finished focalamp'
c       
c find the probable mechanism from the set of acceptable solutions          
c
      call MECH_PROB(nout2,f1norm,f2norm,cangle,prob_max,nmult,
     &        str_avg,dip_avg,rak_avg,prob,var_est)           

c      write(6,*)'Finished mech_prob'

      do 390 imult=1,nmult
      
      var_avg(imult)=(var_est(1,imult)+var_est(2,imult))/2.

c      print *,imult,'  mech = ',
c     &          str_avg(imult),dip_avg(imult),rak_avg(imult)
c
c find misfit for prefered solutions
c
      call GET_MISF_AMP(npol,p_azi_mc,p_the_mc,sp_ratio,
     &      p_pol,str_avg(imult),dip_avg(imult),rak_avg(imult),
     &      mfrac(imult),mavg(imult),stdr(imult))
c      
c solution quality rating, completely ad-hoc - make up your own!
c
      if ((prob(imult).gt.0.8).and.(var_avg(imult).le.25)) then
        qual(imult)='A'
      else if ((prob(imult).gt.0.6).and.(var_avg(imult).le.35)) then
        qual(imult)='B'
      else if ((prob(imult).gt.0.5).and.(var_avg(imult).le.45)) then
        qual(imult)='C'
      else
        qual(imult)='D'
      end if

390   continue

400   continue

c       
c output prefered mechanisms  ** YOU MAY WISH TO CHANGE THE OUTPUT FORMAT **
c
      if (nmult.gt.1) then
        mflag='*'
      else
        mflag=' '
      end if

      do i=1,nmult
cx      write (13,*) icusp,iyr,imon,idy,ihr,imn,qsec,evtype,
cx     &   qmag,magtype,qlat,qlon,qdep,locqual,rms,seh,sez,terr,
cx     &   nppick+nspick,nppick,nspick,
cx     &   nint(str_avg(i)),nint(dip_avg(i)),nint(rak_avg(i)),
cx     &   nint(var_est(1,i)),nint(var_est(2,i)),nppl,nint(mfrac(i)*100.),
cx     &   qual(i),nint(100*prob(i)),nint(100*stdr(i)),nspr,
cx     &   nint(mavg(i)*100.),mflag
c
         write(13,'(a,3f8.1)')'Strike,dip,rake            ',
     *   str_avg(i),dip_avg(i),rak_avg(i)

         write(6,*) 
         write(6,'(a,3f8.1)')' Strike,dip,rake            ',
     *   str_avg(i),dip_avg(i),rak_avg(i)

c
c   write  solution in fps.out
c
         call add_fps(str_avg(i),dip_avg(i),rak_avg(i),'HASH   ','H')

         write(13,'(a,2f8.1)')'Fault+aux plane uncertainty',
     *   var_est(1,i),var_est(2,i)


         write(6,'(a,2f8.1)')' Fault+aux plane uncertainty',
     *   var_est(1,i),var_est(2,i)

         write(13,'(a,f5.2)')
     *  'Weighted fraction of pol misfits',mfrac(i)
         write(13,'(a,f5.2)')
     *  'Average amplitrude error        ',mavg(i)
         write(13,'(a,f5.2)')
     *  'Station dist ratio              ',stdr(i)

cx      write (13,*) icusp,iyr,imon,idy,ihr,imn,qsec,evtype,
cx     &   qmag,magtype,qlat,qlon,qdep,locqual,rms,seh,sez,terr,
cx     &   nppick+nspick,nppick,nspick,
cx     &   nint(str_avg(i)),nint(dip_avg(i)),nint(rak_avg(i)),
cx     &   nint(var_est(1,i)),nint(var_est(2,i)),nppl,nint(mfrac(i)*100.),
cx     &   qual(i),nint(100*prob(i)),nint(100*stdr(i)),nspr,
cx     &   nint(mavg(i)*100.),mflag
c
      end do
411   format(i16,1x,i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,f6.3,1x,a1,1x,
     &  f5.3,1x,a1,1x,f9.5,1x,f10.5,1x,f7.3,1x,a1,1x,f7.3,1x,f7.3,
     &  1x,f7.3,1x,f7.3,3x,i4,1x,i4,1x,i4,1x,i4,1x,i3,1x,i4,3x,i2,
     &  1x,i2,1x,i3,1x,i2,1x,a1,1x,i3,1x,i2,1x,i3,1x,i3,1x,a1)

      
505   continue
      close(11)
      close(12)
      close(13)
      stop
      end




c subroutine FOCALAMP_MC performs grid search to find focal mechanisms, using both
c            P-polarity and S/P amplitude ratio information
c  Inputs:  
c           p_azi_mc(npsta,nmc)  =  azimuth to station from event (deg. E of N)
c           p_the_mc(npsta,nmc)  =  takeoff angle (from vert, <90 upgoing, >90 downgoing)
c           sp_amp(npsta)  =  amplitude ratios
c           p_pol(nspta)   =  P polarities
c           npsta  =  number of observations
c           nmc    =  number of trials
c           dang   =  desired angle spacing for grid search
c           maxout =  maximum number of fault planes to return:
c                     if more are found, a random selection will be returned
c           nextra =  number of polarity additional misfits allowed above minimum
c           ntotal =  total number of allowed polarity misfits
c           qextra =  additional amplitude misfit allowed above minimum
c           qtotal =  total allowed amplitude misfit
c           nspr   =  number of ampltude ratios , jh added
c  Outputs: 
c           nf     =  number of fault planes found
c           strike(min(maxout,nf)) = strike
c           dip(min(maxout,nf))    = dip
c           rake(min(maxout,nf))   = rake
c           faults(3,min(maxout,nf)) = fault normal vector
c           slips(3,min(maxout,nf))  = slip vector
c

      subroutine FOCALAMP_MC(p_azi_mc,p_the_mc,sp_amp,p_pol,npsta,nmc,
     &    dang,maxout,nextra,ntotal,qextra,qtotal,nf,strike,dip,rake,
     &    faults,slips,nspr)
     
c      include 'param.inc'
c npick0 = maximum number of picks per event
c nmc0 = maximum number of trials of location/take-off angles
c nmax0 = maximum number of acceptable mechanisms output
      parameter (npick0=1000,nmc0=500,nmax0=1000) 

c      include 'rot.inc'
c dang0 : minimum grid spacing (degrees)
c ncoor : number of test mechanisms

c choose one, based on your computer's memory and speed
c      parameter (dang0=10.0,ncoor=4032) 
c      parameter (dang0=7.0,ncoor=10426) 
c      parameter (dang0=5.0,ncoor=31032) 
c      parameter (dang0=3.0,ncoor=141180) 
       parameter (dang0=2.0,ncoor=472410) 

      parameter (ntab=180) 

c input and output arrays
      dimension p_azi_mc(npick0,nmc0),p_the_mc(npick0,nmc0)
      real sp_amp(npsta)
      real p_a1(npick0),p_a2(npick0),p_a3(npick0)
      real faultnorm(3),slip(3),faults(3,nmax0),slips(3,nmax0)
      real strike(nmax0),dip(nmax0),rake(nmax0)
      integer p_pol(npsta)
      save dangold,nrot,b1,b2,b3    
      save amptable,phitable,thetable

c coordinate transformation arrays
      real b1(3,ncoor),bb1(3)
      real b2(3,ncoor),bb2(3)
      real b3(3,ncoor),bb3(3)
c P and S amplitude arrays
      real amptable(2,ntab,2*ntab)
      real phitable(2*ntab+1,2*ntab+1)
      real thetable(2*ntab+1)
c misfit arrays
      real qmis(ncoor)
      integer nmis(ncoor)
      integer irotgood(ncoor),irotgood2(ncoor)
      
      pi=3.1415927
      degrad=180./pi

      if (maxout.gt.nmax0) then
        maxout=nmax0
      end if

c Set up array with direction cosines for all coordinate transformations
      if (dang.eq.dangold) go to 8
      irot=0
       do 5 ithe=0,int(90.1/dang)
          the=real(ithe)*dang
         rthe=the/degrad
         costhe=cos(rthe)
         sinthe=sin(rthe)
         fnumang=360./dang
         numphi=nint(fnumang*sin(rthe))
         if (numphi.ne.0) then
            dphi=360./float(numphi)
         else
            dphi=10000.
         end if
         do 4 iphi=0,int(359.9/dphi)
            phi=real(iphi)*dphi
            rphi=phi/degrad
            cosphi=cos(rphi)
            sinphi=sin(rphi)
            bb3(3)=costhe
            bb3(1)=sinthe*cosphi
            bb3(2)=sinthe*sinphi
            bb1(3)=-sinthe
            bb1(1)=costhe*cosphi
            bb1(2)=costhe*sinphi
            call CROSS(bb3,bb1,bb2)
            do 3 izeta=0,int(179.9/dang)
               zeta=real(izeta)*dang
               rzeta=zeta/degrad
               coszeta=cos(rzeta)
               sinzeta=sin(rzeta)
               irot=irot+1
               if (irot.gt.ncoor) then
                  print *,'***FOCAL error: # of rotations too big'
                  return
               end if
               b3(3,irot)=bb3(3)
               b3(1,irot)=bb3(1)
               b3(2,irot)=bb3(2)
               b1(1,irot)=bb1(1)*coszeta+bb2(1)*sinzeta
               b1(2,irot)=bb1(2)*coszeta+bb2(2)*sinzeta                
               b1(3,irot)=bb1(3)*coszeta+bb2(3)*sinzeta
               b2(1,irot)=bb2(1)*coszeta-bb1(1)*sinzeta
               b2(2,irot)=bb2(2)*coszeta-bb1(2)*sinzeta                
               b2(3,irot)=bb2(3)*coszeta-bb1(3)*sinzeta
3           continue
4        continue
5     continue
      nrot=irot
      dangold=dang
      
      astep=1./real(ntab)
      do 150 i=1,2*ntab+1
        bbb3=-1.+real(i-1)*astep
        thetable(i)=aacos(bbb3)
        do 140 j=1,2*ntab+1
          bbb1=-1.+real(j-1)*astep
          phitable(i,j)=atan2(bbb3,bbb1)
          if (phitable(i,j).lt.0.) then
            phitable(i,j)=phitable(i,j)+2.*pi
          end if
140     continue
150   continue

      do 250 i=1,2*ntab
        phi=real(i-1)*pi*astep
        do 240 j=1,ntab
          theta=real(j-1)*pi*astep
          amptable(1,j,i)=abs(sin(2*theta)*cos(phi))                
          s1=cos(2*theta)*cos(phi)  
          s2=-cos(theta)*sin(phi)         
          s2=0.0                        ! jens, null out sv
          amptable(2,j,i)=sqrt(s1*s1+s2*s2)
240     continue
250   continue

8     continue

      do irot=1,nrot
        irotgood(irot)=0
      end do

c loop over multiple trials
      do 430 im=1,nmc 

c  Convert data to Cartesian coordinates
      do 40 i=1,npsta
        call TO_CAR(p_the_mc(i,im),p_azi_mc(i,im),1.,
     &               p_a1(i),p_a2(i),p_a3(i))
40    continue

c  find misfit for each solution and minimum misfit
         nmis0min=1e5
         qmis0min=1.0e5
         do 420 irot=1,nrot  
           qmis(irot)=0.
           nmis(irot)=0
           do 400 ista=1,npsta
             p_b1= b1(1,irot)*p_a1(ista)
     &              +b1(2,irot)*p_a2(ista)
     &              +b1(3,irot)*p_a3(ista) 
             p_b3= b3(1,irot)*p_a1(ista)
     &              +b3(2,irot)*p_a2(ista)
     &              +b3(3,irot)*p_a3(ista) 
             if (sp_amp(ista).ne.0.) then          ! jh: only test for non zero amprat
               p_proj1=p_a1(ista)-p_b3*b3(1,irot)
               p_proj2=p_a2(ista)-p_b3*b3(2,irot)
               p_proj3=p_a3(ista)-p_b3*b3(3,irot)
               plen=sqrt(p_proj1*p_proj1+p_proj2*p_proj2+
     &                    p_proj3*p_proj3)
               p_proj1=p_proj1/plen
               p_proj2=p_proj2/plen
               p_proj3=p_proj3/plen
               pp_b1=b1(1,irot)*p_proj1+b1(2,irot)*p_proj2
     &                +b1(3,irot)*p_proj3
               pp_b2=b2(1,irot)*p_proj1+b2(2,irot)*p_proj2
     &              +b2(3,irot)*p_proj3
               i=nint((p_b3+1.)/astep)+1
               theta=thetable(i)
               i=nint((pp_b2+1.)/astep)+1
               j=nint((pp_b1+1.)/astep)+1
               phi=phitable(i,j)
               i=nint(phi/(pi*astep))+1
               if (i.gt.2*ntab) i=1
               j=nint(theta/(pi*astep))+1
               if (j.gt.ntab) j=1
               p_amp=amptable(1,j,i)
               s_amp=amptable(2,j,i)
               if (p_amp.eq.0.0) then
                 sp_ratio=4.0
               else if (s_amp.eq.0.0) then 
                 sp_ratio=-2.0
               else
                 sp_ratio=log10(4.9*s_amp/p_amp) ! jh: theoretical ratio
               end if
               qmis(irot)=qmis(irot)+abs(sp_amp(ista)-sp_ratio) ! jh: sum amp misfit
                                                                ! from all stations
             end if
             if (p_pol(ista).ne.0) then   ! jh: only count if a polarity
               prod=p_b1*p_b3
               ipol=-1
               if (prod.gt.0.) ipol=1 
               if (ipol.ne.p_pol(ista)) then
                  nmis(irot)=nmis(irot)+1   ! jens, number of misfits this angle                    
               end if
             end if
400         continue                        ! jens, end, station loop

            if (nmis(irot).lt.nmis0min) nmis0min=nmis(irot)  ! jh: minimum misifts
            if (qmis(irot).lt.qmis0min) qmis0min=qmis(irot)
420      continue

        write(6,'(a,i6)')  
     *' Minimum number of polarity misfits overall   :', 
     *  nmis0min  ! jens

        if(nspr.gt.0)
     *   write(6,'(a,f5.2)')
     * ' Minimum average amplitude error overall      : ',qmis0min/nspr
 
         nmis0max=ntotal
c
c jh: if number of polarity misfits + nextra is larger than the given minimum
c     misfits, set new limit to minimum number of misfits found + nextra
c 
         if (nmis0max.lt.nmis0min+nextra) then   ! jh: min pol err nmis0max is >
            nmis0max=nmis0min+nextra             !     allowed, inc. to min+extra
            write(6,'(a,i6)')
     *      ' New number of pol. misfits inc. extra is     :',
     *      nmis0max  ! jens
         end if
c
         qmis0max=qtotal
c
c jh : if the min errro + extra error larger than errr limit, increase to min
c      error+extra 
c
         if (qmis0max.lt.qmis0min+qextra) then
            qmis0max=qmis0min+qextra             ! jh: max amp error
            write(6,'(a,f6.2)')
     *     ' New average amp limit inc. extra             :'
     *     ,qmis0max/nspr  
         end if

c
c loop over rotations - find those meeting fit criteria. 
c
c jh: this means sum of wrong polarities under max level and sum of amplitude
c     errors under max level
c
c
c jh: enter here from below to try with new higher errror limits for amplitude
c
425      nadd=0
         do irot=1,nrot        
            if ((nmis(irot).le.nmis0max).and. 
     &            (qmis(irot).le.qmis0max)) then
              irotgood(irot)=1
              nadd=nadd+1
            end if
         end do
       
c         print *,im,nmis0max,qmis0max,nadd

         if (nadd.eq.0) then  ! if there are none that meet criteria
           qmis0min=1.0e5     ! loosen the amplitude criteria
c
c   jh:loook for smallest amplitude error that fit the criterias for minimum
c   number of polarities
c
           do irot=1,nrot        
             if ((nmis(irot).le.nmis0max).and.
     &           (qmis(irot).lt.qmis0min)) then
               qmis0min=qmis(irot)
             end if
           end do
c
c    jh: add qextra to the amplitude error known to give a solution, then at least
c    one solution will be found
c    
           qmis0max=qtotal
           if (qmis0max.lt.qmis0min+qextra) then
              qmis0max=qmis0min+qextra
           end if
           if(nspr.gt.0)
     *     write(6,'(a,f5.2)')
     *     ' Minimum average amplitude error for pol ok   : ',
     *     qmis0min/nspr
           write(6,'(a,f6.2)')
     *     ' New average amp limit is                     :'
     *     ,qmis0max/nspr  ! jens add
           goto 425   ! jh: try again
         end if

430     continue

        nfault=0
        do irot=1,nrot
          if (irotgood(irot).gt.0) then
            nfault=nfault+1
            irotgood2(nfault)=irot
          end if
        end do

c  Select output solutions  
        nf=0      
        if (nfault.le.maxout) then
          do i=1,nfault
            irot=irotgood2(i)
            nf=nf+1
            faultnorm(1)=b3(1,irot)
            faultnorm(2)=b3(2,irot)
            faultnorm(3)=b3(3,irot)
            slip(1)=b1(1,irot)
            slip(2)=b1(2,irot)
            slip(3)=b1(3,irot)
            do m=1,3
              faults(m,nf)=faultnorm(m)
              slips(m,nf)=slip(m)
            end do
            call FPCOOR(s1,d1,r1,faultnorm,slip,2)
            strike(nf)=s1
            dip(nf)=d1
            rake(nf)=r1
          end do
        else
          do 441 i=1,99999
            fran=rand(0)
            iscr=nint(fran*float(nfault)+0.5)
            if (iscr.lt.1) iscr=1
            if (iscr.gt.nfault) iscr=nfault
            if (irotgood2(iscr).le.0) goto 441
            irot=irotgood2(iscr)
            irotgood2(iscr)=-1
            nf=nf+1
            faultnorm(1)=b3(1,irot)
            faultnorm(2)=b3(2,irot)
            faultnorm(3)=b3(3,irot)
            slip(1)=b1(1,irot)
            slip(2)=b1(2,irot)
            slip(3)=b1(3,irot)
            do m=1,3
              faults(m,nf)=faultnorm(m)
              slips(m,nf)=slip(m)
            end do
            call FPCOOR(s1,d1,r1,faultnorm,slip,2)
            strike(nf)=s1
            dip(nf)=d1
            rake(nf)=r1
            if (nf.eq.maxout) go to 445
441       continue
445       continue
        end if

c
c   jh added
c

c        if(nspr.gt.0)
c     *   write(6,'(a,f6.2)')
c     * ' Minimum average amplitude error      : ',qmis0min/nspr

         
      return 
      end

c ------------------------------------------------------------------- c
      

c subroutine GET_MISF_AMP finds the percent of misfit polarities and the
c                         average S/P ratio misfit for a given mechanism  
c    Inputs:    npol   = number of polarity observations
c               p_azi_mc(npol) = azimuths
c               p_the_mc(npol) = takeoff angles
c               sp_ratio(npol) = S/P ratio
c               p_pol(npol)  = polarity observations
c               str_avg,dip_avg,rak_avg = mechanism
c    Outputs:   mfrac = weighted fraction misfit polarities
c               mavg = average S/P misfit (log10)
c               stdr = station distribution ratio

      subroutine GET_MISF_AMP(npol,p_azi_mc,p_the_mc,sp_ratio,p_pol,
     &          str_avg,dip_avg,rak_avg,mfrac,mavg,stdr)

      dimension p_azi_mc(npol),p_the_mc(npol)
      real str_avg,dip_avg,rak_avg,M(3,3),a(3),b(3),sp_ratio(npol)
      real strike,dip,rake,mfrac,mavg,qcount,azi,toff,pol,wt,wo
      integer k,npol,p_pol(npol)
      real bb1(3),bb2(3),bb3(3)
      
      rad=3.14159265/180.

      strike=str_avg*rad
      dip=dip_avg*rad
      rake=rak_avg*rad
      
      M(1,1)=-sin(dip)*cos(rake)*sin(2*strike)-sin(2*dip)*sin(rake)*
     & sin(strike)*sin(strike)
      M(2,2)=sin(dip)*cos(rake)*sin(2*strike)-sin(2*dip)*sin(rake)*
     & cos(strike)*cos(strike)
      M(3,3)=sin(2*dip)*sin(rake)
      M(1,2)=sin(dip)*cos(rake)*cos(2*strike)+0.5*sin(2*dip)*sin(rake)*
     & sin(2*strike)
      M(2,1)=M(1,2)
      M(1,3)=-cos(dip)*cos(rake)*cos(strike)-cos(2*dip)*sin(rake)*
     & sin(strike)
      M(3,1)=M(1,3)
      M(2,3)=-cos(dip)*cos(rake)*sin(strike)+cos(2*dip)*sin(rake)*
     & cos(strike)
      M(3,2)=M(2,3)
      call FPCOOR(strike,dip,rake,bb3,bb1,1)
      call CROSS(bb3,bb1,bb2)
      
      mfrac=0.
      qcount=0.
      stdr=0.
      scount=0.
      mavg=0.
      acount=0.
      
      do 600 k=1,npol
          call TO_CAR(p_the_mc(k),p_azi_mc(k),1.,p_a1,
     &                p_a2,p_a3)
          p_b1= bb1(1)*p_a1
     &              +bb1(2)*p_a2
     &              +bb1(3)*p_a3 
          p_b3= bb3(1)*p_a1
     &              +bb3(2)*p_a2
     &              +bb3(3)*p_a3
          p_proj1=p_a1-p_b3*bb3(1)
          p_proj2=p_a2-p_b3*bb3(2)
          p_proj3=p_a3-p_b3*bb3(3)
          plen=sqrt(p_proj1*p_proj1+p_proj2*p_proj2+
     &                    p_proj3*p_proj3)
          p_proj1=p_proj1/plen
          p_proj2=p_proj2/plen
          p_proj3=p_proj3/plen
          pp_b1=bb1(1)*p_proj1+bb1(2)*p_proj2
     &              +bb1(3)*p_proj3
          pp_b2=bb2(1)*p_proj1+bb2(2)*p_proj2
     &              +bb2(3)*p_proj3
          phi=atan2(pp_b2,pp_b1)
          theta=aacos(p_b3)
          p_amp=abs(sin(2*theta)*cos(phi))     
          wt=sqrt(p_amp)
          if (p_pol(k).ne.0) then
            azi=rad*p_azi_mc(k)
            toff=rad*p_the_mc(k)        
            a(1)=sin(toff)*cos(azi)
            a(2)=sin(toff)*sin(azi)
            a(3)=-cos(toff)
            do 615 in=1,3
              b(in)=0
              do 610 jn=1,3
                 b(in)=b(in)+M(in,jn)*a(jn)
610           continue
615         continue
            if ((a(1)*b(1)+a(2)*b(2)+a(3)*b(3)).lt.0) then
              pol=-1
            else
             pol=1
            end if
            if ((pol*p_pol(k)).lt.0) then
              mfrac=mfrac+wt
            end if
            qcount=qcount+wt
            stdr=stdr+wt
            scount=scount+1.0
          end if
          if (sp_ratio(k).ne.0.) then
            s1=cos(2*theta)*cos(phi)  
            s2=-cos(theta)*sin(phi)
c            s2=0.0                     ! jens   null out sv
            s_amp=sqrt(s1*s1+s2*s2)
            sp_rat=log10(4.9*s_amp/p_amp)

            mavg=mavg+abs(sp_ratio(k)-sp_rat)

            write(17,*) k,sp_ratio(k),sp_rat

            acount=acount+1.0
            stdr=stdr+wt
            scount=scount+1.0
          end if
600    continue
       mfrac=mfrac/qcount
       if (qcount.eq.0.0) mfrac=0.0

       write(17,*)mavg,acount

       mavg=mavg/acount

       if (acount.eq.0.0) mavg=0.0
       stdr=stdr/scount
       if (scount.eq.0.0) stdr=0.0
       
       return 
       end




c subroutine GET_MISF finds the percent of misfit polarities for a given mechanism  
c    Inputs:    npol   = number of polarity observations
c               p_azi_mc(npol) = azimuths
c               p_the_mc(npol) = takeoff angles
c               p_pol(npol)  = polarity observations
c               p_qual(npol) = quality of polarity observations
c               str_avg,dip_avg,rak_avg = mechanism
c    Outputs:   mfrac = weighted fraction misfit polarities (like FPFIT)
c               stdr = station distribution ratio (like FPFIT)

      subroutine GET_MISF(npol,p_azi_mc,p_the_mc,p_pol,p_qual,str_avg,
     &                  dip_avg,rak_avg,mfrac,stdr)

      dimension p_azi_mc(npol),p_the_mc(npol)
      real str_avg,dip_avg,rak_avg,M(3,3),a(3),b(3)
      real strike,dip,rake,mfrac,qcount,azi,toff,pol,wt,wo
      integer k,npol,p_pol(npol),p_qual(npol)
      real bb1(3),bb2(3),bb3(3)
      
      rad=3.14159265/180.

      strike=str_avg*rad
      dip=dip_avg*rad
      rake=rak_avg*rad
      
      M(1,1)=-sin(dip)*cos(rake)*sin(2*strike)-sin(2*dip)*sin(rake)*
     & sin(strike)*sin(strike)
      M(2,2)=sin(dip)*cos(rake)*sin(2*strike)-sin(2*dip)*sin(rake)*
     & cos(strike)*cos(strike)
      M(3,3)=sin(2*dip)*sin(rake)
      M(1,2)=sin(dip)*cos(rake)*cos(2*strike)+0.5*sin(2*dip)*sin(rake)*
     & sin(2*strike)
      M(2,1)=M(1,2)
      M(1,3)=-cos(dip)*cos(rake)*cos(strike)-cos(2*dip)*sin(rake)*
     & sin(strike)
      M(3,1)=M(1,3)
      M(2,3)=-cos(dip)*cos(rake)*sin(strike)+cos(2*dip)*sin(rake)*
     & cos(strike)
      M(3,2)=M(2,3)
      call FPCOOR(strike,dip,rake,bb3,bb1,1)
      call CROSS(bb3,bb1,bb2)
      mfrac=0.
      qcount=0.
      scount=0.
      
      do 600 k=1,npol
          call TO_CAR(p_the_mc(k),p_azi_mc(k),1.,p_a1,
     &                p_a2,p_a3)
          p_b1= bb1(1)*p_a1
     &              +bb1(2)*p_a2
     &              +bb1(3)*p_a3 
          p_b3= bb3(1)*p_a1
     &              +bb3(2)*p_a2
     &              +bb3(3)*p_a3
          p_proj1=p_a1-p_b3*bb3(1)
          p_proj2=p_a2-p_b3*bb3(2)
          p_proj3=p_a3-p_b3*bb3(3)
          plen=sqrt(p_proj1*p_proj1+p_proj2*p_proj2+
     &                    p_proj3*p_proj3)
          pp_b1=bb1(1)*p_proj1+bb1(2)*p_proj2
     &              +bb1(3)*p_proj3
          pp_b2=bb2(1)*p_proj1+bb2(2)*p_proj2
     &              +bb2(3)*p_proj3
          phi=atan2(pp_b2,pp_b1)
          theta=aacos(p_b3)
          p_amp=abs(sin(2*theta)*cos(phi))     
          wt=sqrt(p_amp)
         azi=rad*p_azi_mc(k)
         toff=rad*p_the_mc(k)        
         a(1)=sin(toff)*cos(azi)
         a(2)=sin(toff)*sin(azi)
         a(3)=-cos(toff)
         do in=1,3
           b(in)=0
           do jn=1,3
              b(in)=b(in)+M(in,jn)*a(jn)
           end do
         end do
         if ((a(1)*b(1)+a(2)*b(2)+a(3)*b(3)).lt.0) then
           pol=-1
         else
           pol=1
         end if
         if (p_qual(k).eq.0) then
           wo=1
         else
           wo=0.5
         end if
         if ((pol*p_pol(k)).lt.0) then
           mfrac=mfrac+wt*wo
         end if
         qcount=qcount+wt*wo
         scount=scount+wo
600   continue
      mfrac=mfrac/qcount
      stdr=qcount/scount
      
      return 
      end

c --------------------------------------------------------------- c


c subroutine GET_GAP finds the maximum azimuthal and takeoff angle gaps  
c    Inputs:    npol   = number of polarity observations
c               p_azi_mc(npol) = azimuths
c               p_the_mc(npol) = takeoff angles
c    Outputs:   magap  = maximum azimuthal gap
c               mpgap  = maximum takeoff angle gap

      subroutine GET_GAP(npol,p_azi_mc,p_the_mc,magap,mpgap)

c      include 'param.inc'
c npick0 = maximum number of picks per event
c nmc0 = maximum number of trials of location/take-off angles
c nmax0 = maximum number of acceptable mechanisms output
      parameter (npick0=1000,nmc0=500,nmax0=1000) 
      dimension p_azi_mc(npol),p_the_mc(npol)
      real p2_azi(npick0),p2_the(npick0)

      do 403 k=1,npol
        if (p_the_mc(k).gt.90) then
          p2_the(k)=180.-p_the_mc(k)
          p2_azi(k)=p_azi_mc(k)-180.
          if (p2_azi(k).lt.0) p2_azi(k)=p2_azi(k)+360.
        else
          p2_the(k)=p_the_mc(k)
          p2_azi(k)=p_azi_mc(k)
        end if
403   continue
      call sort(npol,p2_azi)
      call sort(npol,p2_the)
      magap=0
      mpgap=0
      do 405 k=2,npol
        if (p2_azi(k)-p2_azi(k-1).gt.magap) then
          magap=p2_azi(k)-p2_azi(k-1)
        end if
        if (p2_the(k)-p2_the(k-1).gt.mpgap) then
          mpgap=p2_the(k)-p2_the(k-1)
        end if
405   continue
      if (p2_azi(1)-p2_azi(npol)+360.gt.magap) then
        magap=p2_azi(1)-p2_azi(npol)+360
      end if
      if (90.-p2_the(npol).gt.mpgap) then
        mpgap=90.-p2_the(npol)
      end if
      if (p2_the(1).gt.mpgap) then
        mpgap=p2_the(1)
      end if
      
      return 
      end

c --------------------------------------------------------------- c

      SUBROUTINE SORT(N,RA)   ! modified from numerical recipies
      DIMENSION RA(N)
      if (n.eq.0) then
         print *,'***n=0 in SORT'
         return
      end if
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(L)
        ELSE
          RRA=RA(IR)
          RA(IR)=RA(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(1)=RRA
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(RA(J).LT.RA(J+1))J=J+1
          ENDIF
          IF(RRA.LT.RA(J))THEN
            RA(I)=RA(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)=RRA
      GO TO 10
      END

c GETSTAT_TRI finds station locations for TriNet stations
c
c   inputs:
c     stlfile - file with the stations and locations, in alphabetical order!
c     snam - name of the station of interest, 4 characters
c     scom - station component, 3 characters
c     snet - network, 2 characters
c   outputs:
c     flat,flon,felev - station lat, lon, and elevation
c
c   input file "stlfile" format:
c     columns  format   value
c     -------------------------
c     1-4        a4     station name
c     6-8        a3     station component
c     42-50      f9.5   station latitude (degrees, signed)
c     52-61      f10.5  station longitude (degrees, signed)
c     63-67      i5     station elevation (meters)
c     91-92      a2     network code
c
c
      subroutine GETSTAT_TRI(stlfile,snam,scom,snet,flat,
     &                       flon,felev)
      parameter(nsta0=20000)
      character stlfile*100
      character*4 snam,stname(nsta0)
      character*3 scom,scompt(nsta0),scom2
      character*2 snet,snetwk(nsta0)
      real slat(nsta0),slon(nsta0),selev(nsta0)
      logical firstcall
      save firstcall,stname,slat,slon,selev,nsta,scompt,snetwk
      data firstcall/.true./
      
c read in station list - in alphabetical order!
      if (firstcall) then
         firstcall=.false.
         open (19,file=stlfile)
         do i=1,nsta0
           read (19,11,end=12) stname(i),scompt(i),slat(i),
     &                         slon(i),ntemp,snetwk(i)
           selev(i)=real(ntemp)/1000.0
         end do
11       format (a4,1x,a3,33x,f9.5,1x,f10.5,1x,i5,23x,a2)
12       nsta=i-1
         close (19)
      end if  
       
      scom2=scom                             ! short-period stations are
      if (scom(1:1).eq."V") scom2(1:1)="E"   ! interchangably called V and E   
      if (scom(1:1).eq."E") scom2(1:1)="V"           

c binary search for station name
      i1=1
      i2=nsta
      do 30 it=1,30
         i=(i1+i2)/2
         if (snam.eq.stname(i)) then
           goto 40
         else if (i1.eq.i2) then
           goto 999
         else if (snam.lt.stname(i)) then
            i2=i-1
         else 
            i1=i+1
         end if
30    continue
      goto 999
      
c search for proper component/network
40    i1=i
45    continue
        if (i1.gt.nsta) goto 50
        if (snam.ne.stname(i1)) goto 50
        if ((scom(1:2).eq.scompt(i1)(1:2)).and.
     &             snet.eq.snetwk(i1)) goto 900
        if ((scom2(1:2).eq.scompt(i1)(1:2)).and.
     &             snet.eq.snetwk(i1)) goto 900
        if ((scom(1:2).eq.'XX').and.
     &                 (snet.eq.'XX')) goto 900
        i1=i1+1
      goto 45
50    i1=i-1
55    continue
        if (i1.lt.1) goto 999
        if (snam.ne.stname(i1)) goto 999
        if ((scom(1:2).eq.scompt(i1)(1:2)).and.
     &             snet.eq.snetwk(i1)) goto 900
        if ((scom2(1:2).eq.scompt(i1)(1:2)).and.
     &             snet.eq.snetwk(i1)) goto 900
        if ((scom(1:2).eq.'XX').and.
     &                 (snet.eq.'XX')) goto 900
        i1=i1-1
      goto 55
      
900   flat=slat(i1)
      flon=slon(i1)
      felev=selev(i1)
      return
999   print *,'***station not found ',snam,' ',scom,' ',snet
      flat=999.
      flon=999.
      felev=999.
      return
      end

c --------------------------------------------------------------- c


c CHECK_POL determines whether the polarity of a given station
c     was reversed at a given time.
c
c   inputs:
c     polfile - file with the stations and times of reversal;
c               must be in alphabetical order; if a station is
c               reversed during multiple time periods, it is 
c               listed on multiple lines (same as FPFIT)
c     snam - name of the station of interest, 4 characters
c     evyr  - year of interest, 4 digits
c     evmon - month of interest
c     evdy  - day of interest
c     evhr  - hour of interest (not implemented)
c   output:
c     stpol - station polarity: 1=correct, -1=reversed 
c
c   input file "polfile" format:
c     columns  format   value
c     -------------------------
c     1-4        a4     station name
c     6-9        i4     begining of reversal: year
c     10-11      i2                           month
c     12-13      i2                           day
c     15-18      i4     end of reversal: year
c     19-20      i2                      month
c     21-22      i2                      day

      subroutine CHECK_POL(polfile,snam,evyr,evmon,evdy,evhr,stpol)
      parameter(nsta0=300)
      character polfile*100,polfileold*100
      character*4 snam,statname(nsta0)
      integer begtime(nsta0),endtime(nsta0)
      integer evyr,evmon,evdy,evhr,stpol,nstat(nsta0),nth(nsta0)
      integer i,i1,i2,count,itemp,nrev,evtime
      save polfileold,statname,begtime,endtime,nrev,nstat,nth

c read in polarity reversal file - in alphabetical order
      if (polfile.ne.polfileold) then
         print *,'Reading polarity reversal file ',polfile
         open (19,file=polfile)
         do i=1,300
           read (19,11,end=12) statname(i),begtime(i),endtime(i)  
           if (endtime(i).eq.0) then
             endtime(i)=99999999
           end if
         end do
11       format (a4,1x,i8,1x,i8)
12       nrev=i-1
         close (19)
         polfileold=polfile
         nth(1)=1
         nstat(1)=1
         count=1
         do 17 i=2,nrev
           if (statname(i).eq.statname(i-1)) then
             count=count+1
           else
             count=1
           end if
           nth(i)=count
           do itemp=i-count+1,i
             nstat(itemp)=count
           end do
17       continue
      end if
      
      evtime=evyr*10000+evmon*100+evdy
      stpol=1

c binary search for a station
      i1=1
      i2=nrev
20    continue
         i=(i1+i2)/2 
         if (snam.eq.statname(i)) then
            do 25 itemp=i-nth(i)+1,i+nstat(i)-nth(i)
              if ((evtime.ge.begtime(itemp)).and.
     &                (evtime.le.endtime(itemp))) then
                stpol=-1
                goto 30
              end if
25          continue
            goto 30
         else if (i1.ge.i2) then
            goto 30
         else if (snam.lt.statname(i)) then
            i2=i-1
         else
            i1=i+1
         end if
      goto 20
30    continue
      return
      end

c ------------------------------------------------------------------- c
      
c GET_COR reads a file of station amplitude corrections
c
c   inputs:
c     stlfile - file with the stations and locations
c     snam - name of the station of interest, 4 characters
c     scom - station component, 3 characters
c     snet - network, 2 characters
c   outputs:
c     qcor - corrections to be subtracted from log(S/P)
c
c   input file format:
c     columns  format   variable
c     -------------------------
c     1-4        a4     station name
c     7-9        a3     station component (vertical - Z optional)
c     11-12      a2     network code
c     14-20      f7.4   correction to be subtracted from log(S/P)
c
c
      subroutine GET_COR(stlfile,snam,scom,snet,qcor)
      parameter(nsta0=10000)
      character stlfile*100
      character*4 snam,stname(nsta0)
      character*3 scom,scom2,scompt(nsta0)
      character*2 snet,snetwk(nsta0)
      real corr_val(nsta0)
      logical firstcall
      save firstcall,stname,corr_val,nsta,scompt,snetwk
      data firstcall/.true./
      
c read in station list - in alphabetical order!
      if (firstcall) then
         firstcall=.false.
         open (19,file=stlfile)
         do 10 i=1,nsta0
           read (19,11,end=12) stname(i),scompt(i),snetwk(i),
     &                              corr_val(i)
10       continue
11       format (a4,2x,a3,a2,1x,f7.4)
12       nsta=i-1
         close (19)
      end if  
      
      scom2=scom                             ! short-period stations are
      if (scom(1:1).eq."V") scom2(1:1)="E"   ! called both V and E     
      if (scom(1:1).eq."E") scom2(1:1)="V"           

c binary search for station name
      i1=1
      i2=nsta
      do 30 it=1,30
         i=(i1+i2)/2
         if (snam.eq.stname(i)) then
           goto 40
         else if (i1.eq.i2) then
           goto 999
         else if (snam.lt.stname(i)) then
            i2=i-1
         else 
            i1=i+1
         end if
30    continue
      print *,'station not found'
      goto 999
      
c search for proper component/network
40    i1=i
45    continue
        if (i1.gt.nsta) goto 50
        if (snam.ne.stname(i1)) goto 50
        if (scom(1:2).eq.scompt(i1)(1:2)) goto 900
        if (scom2(1:2).eq.scompt(i1)(1:2)) goto 900
        i1=i1+1
      goto 45
50    i1=i-1
55    continue
        if (i1.lt.1) goto 999
        if (snam.ne.stname(i1)) goto 999
        if (scom(1:2).eq.scompt(i1)(1:2)) goto 900
        if (scom2(1:2).eq.scompt(i1)(1:2)) goto 900
        i1=i1-1
      goto 55

900   qcor=corr_val(i1)
      return
999   print *,'GET_COR ***station not found ',snam,' ',scom,' ',snet,
     &  ' in file ',stlfile
      qcor=-999.
      return
      end

c subroutine MECH_PROB determines the average focal mechanism of a set
c   of mechanisms, after removing outliers, and returns the probability
c   of the mechanism being within a given cutoff angle of the result,
c   also checks outliers for multiple solutions and returns any
c
c  Inputs:  nf     =  number of fault planes
c           norm1(3,nf) = normal to fault plane
c           norm2(3,nf) = slip vector
c           cangle  =  cutoff angle 
c           prob_max = cutoff percent for multiples
c  Output:  nsltn = number of solutions, up to 5
c           str_avg(5)   = average strike
c           dip_avg(5)   = average dip
c           rak_avg(5)   = average rake
c           prob(5)      = percent of mechs within cutoff angle
c                          of average mechanism
c           rms_diff(2,5)  = RMS angular difference of all planes to average 
c                          plane (1=fault plane, 2=auxiliary plane)
c

      subroutine MECH_PROB(nf,norm1in,norm2in,cangle,prob_max,nsltn,
     &             str_avg,dip_avg,rak_avg,prob,rms_diff)

c      include 'param.inc'
c npick0 = maximum number of picks per event
c nmc0 = maximum number of trials of location/take-off angles
c nmax0 = maximum number of acceptable mechanisms output
      parameter (npick0=1000,nmc0=500,nmax0=1000) 
      
      integer nf
      real str_avg(5),dip_avg(5),rak_avg(5),rota(nmax0)
      real norm1(3,nmax0),norm2(3,nmax0),temp1(3),temp2(3)
      real norm1in(3,nf),norm2in(3,nf),ln_norm1,ln_norm2
      real norm1_avg(3),norm2_avg(3),slipol,rms_diff(2,5)
      real stv(2),udv(3),dd_rad,di_rad,a,b,prob(5)

      pi=3.1415927
      degrad=180./3.1415927

c if there is only one mechanism, return that mechanism
   
      if (nf.le.1) then
        do i=1,3
          norm1_avg(i)=norm1in(i,1)
          norm2_avg(i)=norm2in(i,1)
        end do
        call fpcoor(str_avg(1),dip_avg(1),rak_avg(1),
     &              norm1_avg,norm2_avg,2)
        prob(1)=1.
        rms_diff(1,1)=0.
        rms_diff(2,1)=0.
        nsltn=1
        return
      end if
      
c otherwise find the prefered mechanism and any multiples

      do j=1,nf
        do i=1,3
          norm1(i,j)=norm1in(i,j)
          norm2(i,j)=norm2in(i,j)
        end do
      end do

      nfault=nf
      nc=nf
      do 380 imult=1,5
        if (nc.lt.1) goto 385

c find the average mechanism by repeatedly throwing out the mechanism
c with the largest angular difference from the average, and finding the average of 
c the remaining mechanisms - stop when all mechanisms are within cangle of average
   
      do 77 icount=1,nf 
        call MECH_AVG(nc,norm1,norm2,norm1_avg,norm2_avg) 
        do i=1,nc      
          do j=1,3
            temp1(j)=norm1(j,i)
            temp2(j)=norm2(j,i)
          end do
          call MECH_ROT(norm1_avg,temp1,norm2_avg,temp2,rota(i))
        end do
        maxrot=0.
        do i=1,nc
          if (abs(rota(i)).gt.maxrot) then
            maxrot=abs(rota(i))
            imax=i
          end if
        end do
        if (maxrot.le.cangle) goto 78
        nc=nc-1
        do i=1,3
          temp1(i)=norm1(i,imax)
          temp2(i)=norm2(i,imax)
        end do
        do j=imax,nc            
          do i=1,3
            norm1(i,j)=norm1(i,j+1)
            norm2(i,j)=norm2(i,j+1)
          end do
        end do
        do i=1,3
          norm1(i,nc+1)=temp1(i)
          norm2(i,nc+1)=temp2(i)
        end do
77    continue

78    continue
      a=nc
      b=nfault
      prob(imult)=a/b 

      if ((imult.gt.1).and.(prob(imult).lt.prob_max)) goto 385

      do j=1,nf-nc   ! set up for next round
        do i=1,3
          norm1(i,j)=norm1(i,j+nc)
          norm2(i,j)=norm2(i,j+nc)
        end do
      end do
      nc=nf-nc 
      nf=nc

c determine the RMS observed angular difference between the average 
c normal vectors and the normal vectors of each mechanism

      rms_diff(1,imult)=0.
      rms_diff(2,imult)=0.
      do 80 i=1,nfault
        do j=1,3
          temp1(j)=norm1in(j,i)
          temp2(j)=norm2in(j,i)
        end do
        call MECH_ROT(norm1_avg,temp1,norm2_avg,temp2,rota)
        d11=temp1(1)*norm1_avg(1)+temp1(2)*norm1_avg(2)+
     &                            temp1(3)*norm1_avg(3)
        d22=temp2(1)*norm2_avg(1)+temp2(2)*norm2_avg(2)+
     &                            temp2(3)*norm2_avg(3)
        if (d11.ge.1.) d11=1.
        if (d11.le.-1.) d11=-1.
        if (d22.ge.1.) d22=1.
        if (d22.le.-1.) d22=-1.
        a11=aacos(d11)
        a22=aacos(d22)
        rms_diff(1,imult)=rms_diff(1,imult)+a11*a11
        rms_diff(2,imult)=rms_diff(2,imult)+a22*a22
80    continue
      rms_diff(1,imult)=degrad*sqrt(rms_diff(1,imult)/nfault)
      rms_diff(2,imult)=degrad*sqrt(rms_diff(2,imult)/nfault)

      call fpcoor(str_avg(imult),dip_avg(imult),rak_avg(imult),
     &            norm1_avg,norm2_avg,2)

380   continue

385   nsltn=imult-1  ! only use ones with certain probability
          
      return 
      end

c ------------------------------------------------------------ c

c subroutine MECH_AVG determines the average focal mechanism of a set
c   of mechanisms
c
c  Inputs:  nf     =  number of fault planes
c           norm1(3,nf) = normal to fault plane
c           norm2(3,nf) = slip vector
c  Output:  norm1_avg(3)    = normal to avg plane 1
c           norm2_avg(3)    = normal to avg plane 2
c
c    Written  10/4/2000 by Jeanne Hardebeck                              
c    Modified 5/14/2001 by Jeanne Hardebeck                              
c
      subroutine MECH_AVG(nf,norm1,norm2,norm1_avg,norm2_avg)
           
      real dot1,fract1
      real misf,maxmisf,avang1,avang2
      real norm1(3,nf),norm2(3,nf),temp1(3),temp2(3)
      real norm1_avg(3),norm2_avg(3),ln_norm1,ln_norm2
      real theta1,theta2,ref1(3),ref2(3)
      integer nf

      pi=3.1415927
      degrad=180./3.1415927
      
c if there is only one mechanism, return that mechanism
   
      if (nf.le.1) then
        do 5 i=1,3
          norm1_avg(i)=norm1(i,1)
          norm2_avg(i)=norm2(i,1)
5       continue
        goto 120
      end if
            
c find the average normal vector for each plane - determine which
c nodal plane of each event corresponds to which by finding the
c minimum focal mechanism rotation

      do j=1,3
        norm1_avg(j)=norm1(j,1)
        norm2_avg(j)=norm2(j,1)
        ref1(j)=norm1(j,1)
        ref2(j)=norm2(j,1)
      end do
      do 50 i=2,nf
        do j=1,3
          temp1(j)=norm1(j,i)
          temp2(j)=norm2(j,i)
        end do
        call MECH_ROT(ref1,temp1,ref2,temp2,rota)
        do j=1,3
          norm1_avg(j)=norm1_avg(j)+temp1(j)
          norm2_avg(j)=norm2_avg(j)+temp2(j)
        end do
50    continue
      ln_norm1=0
      ln_norm2=0
      do 60 j=1,3
        ln_norm1=ln_norm1+norm1_avg(j)*norm1_avg(j)
        ln_norm2=ln_norm2+norm2_avg(j)*norm2_avg(j)
60    continue
      ln_norm1=sqrt(ln_norm1)
      ln_norm2=sqrt(ln_norm2)
      do 70 i=1,3
        norm1_avg(i)=norm1_avg(i)/ln_norm1
        norm2_avg(i)=norm2_avg(i)/ln_norm2
70    continue

c determine the RMS observed angular difference between the average 
c normal vectors and the normal vectors of each mechanism

      avang1=0.
      avang2=0.
      do 80 i=1,nf
        do j=1,3
          temp1(j)=norm1(j,i)
          temp2(j)=norm2(j,i)
        end do
        call MECH_ROT(norm1_avg,temp1,norm2_avg,temp2,rota)
        d11=temp1(1)*norm1_avg(1)+temp1(2)*norm1_avg(2)+
     &                            temp1(3)*norm1_avg(3)
        d22=temp2(1)*norm2_avg(1)+temp2(2)*norm2_avg(2)+
     &                            temp2(3)*norm2_avg(3)
        if (d11.ge.1.) d11=1.
        if (d11.le.-1.) d11=-1.
        if (d22.ge.1.) d22=1.
        if (d22.le.-1.) d22=-1.
        a11=aacos(d11)
        a22=aacos(d22)
        avang1=avang1+a11*a11
        avang2=avang2+a22*a22
80    continue
      avang1=sqrt(avang1/nf)
      avang2=sqrt(avang2/nf)

c the average normal vectors may not be exactly orthogonal (although
c usually they are very close) - find the misfit from orthogonal and 
c adjust the vectors to make them orthogonal - adjust the more poorly 
c constrained plane more
 
      if ((avang1+avang2).lt.0.0001) goto 120

      maxmisf=0.01
      fract1=avang1/(avang1+avang2)
90    do 115 icount=1,100  
        dot1=norm1_avg(1)*norm2_avg(1)+norm1_avg(2)
     &     *norm2_avg(2)+norm1_avg(3)*norm2_avg(3)
        misf=90.-aacos(dot1)*degrad
        if (abs(misf).le.maxmisf) goto 120
        theta1=misf*fract1/degrad
        theta2=misf*(1.-fract1)/degrad
        do j=1,3
          temp=norm1_avg(j)
          norm1_avg(j)=norm1_avg(j)-norm2_avg(j)*sin(theta1)
          norm2_avg(j)=norm2_avg(j)-temp*sin(theta2)
        end do
        ln_norm1=0
        ln_norm2=0
        do j=1,3
          ln_norm1=ln_norm1+norm1_avg(j)*norm1_avg(j)
          ln_norm2=ln_norm2+norm2_avg(j)*norm2_avg(j)
        end do
        ln_norm1=sqrt(ln_norm1)
        ln_norm2=sqrt(ln_norm2)
        do i=1,3
          norm1_avg(i)=norm1_avg(i)/ln_norm1
          norm2_avg(i)=norm2_avg(i)/ln_norm2
        end do
115   continue

120   continue      
      return 
      end
      

c ------------------------------------------------------------ c


c subroutine MECH_ROT finds the minimum rotation angle between two
c mechanisms, given NORMAL and SLIP (doesn't work for P & T axes!)
c
c  Inputs:  norm1(3) = normal to fault plane 1
c           slip1(3) = slip vector 1
c           norm2(3) = normal to fault plane 2
c           slip2(3) = slip vector 2
c  Output:  rota  = rotation angle
c
c Edited 10/23/07 JLH - Does NOT assume that the normal and slip vectors
c                       have been matched!  
c Tries 4 different combinations, then CHANGES norm2 & slip2 to the best combo
c     (1) norm1 & slip1 <=> norm2 & slip2
c     (2) norm1 & slip1 <=> -norm2 & -slip2
c     (3) norm1 & slip1 <=> slip2 & norm2
c     (4) norm1 & slip1 <=> -slip2 & -norm2
c
c
      subroutine MECH_ROT(norm1,norm2,slip1,slip2,rota)
      
      real norm1(3),norm2(3),slip1(3),slip2(3),B1(3),B2(3)
      real norm2_temp(3),slip2_temp(3),rotemp(4)
      real rota,phi(3),n(3,3),scale(3),R(3),qdot(3)
      real theta(3),n1(3),n2(3),phi1
      
      pi=3.1415927
      degrad=180./3.1415927

      do 200 iter=1,4    ! iteration over the 4 possibilities

      if (iter.lt.3) then
        do i=1,3
          norm2_temp(i)=norm2(i)
          slip2_temp(i)=slip2(i)
        end do
      else
        do i=1,3
          norm2_temp(i)=slip2(i)
          slip2_temp(i)=norm2(i)
        end do
      end if
      if ((iter.eq.2).or.(iter.eq.4)) then
        do i=1,3
          norm2_temp(i)=-norm2_temp(i)
          slip2_temp(i)=-slip2_temp(i)
        end do
      end if

      call cross(norm1,slip1,B1)
      call cross(norm2_temp,slip2_temp,B2)

      phi1=norm1(1)*norm2_temp(1)+norm1(2)*norm2_temp(2)+
     &     norm1(3)*norm2_temp(3)
      phi(1)=aacos(phi1)
      phi1=slip1(1)*slip2_temp(1)+slip1(2)*slip2_temp(2)+
     &     slip1(3)*slip2_temp(3)
      phi(2)=aacos(phi1)
      phi1=B1(1)*B2(1)+B1(2)*B2(2)+B1(3)*B2(3)
      phi(3)=aacos(phi1)

c if the mechanisms are very close, rotation = 0
      if ((phi(1).lt.1e-4).and.(phi(2).lt.1e-4).and.
     &    (phi(3).lt.1e-4)) then
        rotemp(iter)=0.0
c if one vector is the same, it is the rotation axis
      else if (phi(1).lt.1e-4) then
        rotemp(iter)=degrad*phi(2)
      else if (phi(2).lt.1e-4) then
        rotemp(iter)=degrad*phi(3)
      else if (phi(3).lt.1e-4) then
        rotemp(iter)=degrad*phi(1)
      else
c find difference vectors - the rotation axis must be orthogonal
c to all three of these vectors
        do i=1,3
          n(i,1)=norm1(i)-norm2_temp(i)
          n(i,2)=slip1(i)-slip2_temp(i)
          n(i,3)=B1(i)-B2(i)
        end do
        do j=1,3
          scale(j)=sqrt(n(1,j)*n(1,j)+n(2,j)*n(2,j)+n(3,j)*n(3,j))
          do i=1,3
            n(i,j)=n(i,j)/scale(j)
          end do
        end do
        qdot(3)=n(1,1)*n(1,2)+n(2,1)*n(2,2)+n(3,1)*n(3,2)
        qdot(2)=n(1,1)*n(1,3)+n(2,1)*n(2,3)+n(3,1)*n(3,3)
        qdot(1)=n(1,2)*n(1,3)+n(2,2)*n(2,3)+n(3,2)*n(3,3)
c use the two largest difference vectors, as long as they aren't orthogonal 
        iout=0
        do i=1,3
          if (qdot(i).gt.0.9999) iout=i
        end do
        if (iout.eq.0) then
          qmins=10000.
          do i=1,3
            if (scale(i).lt.qmins) then
              qmins=scale(i)
              iout=i
            end if
          end do
        end if
        k=1
        do j=1,3
          if (j.ne.iout) then
            if (k.eq.1) then
              do i=1,3
                n1(i)=n(i,j)
              end do
              k=2
            else
              do i=1,3
                n2(i)=n(i,j)
              end do
            end if
          end if
        end do
c  find rotation axis by taking cross product
         call CROSS(n1,n2,R)
         scaleR=sqrt(R(1)*R(1)+R(2)*R(2)+R(3)*R(3))
         do i=1,3
           R(i)=R(i)/scaleR
         end do
c find rotation using axis furthest from rotation axis
         theta(1)=aacos(norm1(1)*R(1)+norm1(2)*R(2)+norm1(3)*R(3))
         theta(2)=aacos(slip1(1)*R(1)+slip1(2)*R(2)+slip1(3)*R(3))
         theta(3)=aacos(B1(1)*R(1)+B1(2)*R(2)+B1(3)*R(3))
         qmindif=1000.
         do i=1,3
           if (abs(theta(i)-pi/2.0).lt.qmindif) then
             qmindif=abs(theta(i)-pi/2.0)
             iuse=i
           end if
         end do
         rotemp(iter)=(cos(phi(iuse))-cos(theta(iuse))*
     &     cos(theta(iuse)))/(sin(theta(iuse))*sin(theta(iuse)))
         if (rotemp(iter).gt.1.0) then
           rotemp(iter)=1.0
         end if
         if (rotemp(iter).lt.-1.0) then
           rotemp(iter)=-1.0
         end if
         rotemp(iter)=degrad*aacos(rotemp(iter))
       end if

200    continue

c find the minimum rotation for the 4 combos, and change norm2 and slip2
       rota=180.0
       do iter=1,4
         if (abs(rotemp(iter)).lt.rota) then
           rota=abs(rotemp(iter))
           irot=iter
         end if
       end do
       if (irot.ge.3) then
         do i=1,3
           qtemp=slip2(i)
           slip2(i)=norm2(i)
           norm2(i)=qtemp
         end do
       end if
       if ((irot.eq.2).or.(irot.eq.4)) then
         do i=1,3
           norm2(i)=-norm2(i)
           slip2(i)=-slip2(i)
         end do
       end if

       return
       end      

c
c  cross product of two vectors, sets v3 = v1 x v2
c
      subroutine CROSS(v1,v2,v3)
      real v1(3),v2(3),v3(3)
      v3(1)=v1(2)*v2(3)-v1(3)*v2(2)   
      v3(2)=v1(3)*v2(1)-v1(1)*v2(3)
      v3(3)=v1(1)*v2(2)-v1(2)*v2(1)
      return
      end   

c ------------------------------------------------------------ c

c
c  transforms spherical co-ordinates to cartesian
c
      subroutine TO_CAR(the,phi,r,x,y,z)
      degrad=3.1415927/180.
      z=-r*cos(the*degrad)
      x=r*sin(the*degrad)*cos(phi*degrad)
      y=r*sin(the*degrad)*sin(phi*degrad)
      return
      end


c ------------------------------------------------------------ c

c subroutine FPCOOR gets fault normal vector,fnorm, and slip 
c vector, slip, from (strike,dip,rake) or vice versa.
c   idir = 1 compute fnorm,slip
c   idir = c compute strike,dip,rake
c Reference:  Aki and Richards, p. 115
c   uses (x,y,z) coordinate system with x=north, y=east, z=down
      subroutine FPCOOR(strike,dip,rake,fnorm,slip,idir)
      real fnorm(3),slip(3),phi,del,lam,a,clam,slam
      degrad=180./3.1415927
      pi=3.1415927
      phi=strike/degrad
      del=dip/degrad
      lam=rake/degrad
      if (idir.eq.1) then
         fnorm(1)=-sin(del)*sin(phi)
         fnorm(2)= sin(del)*cos(phi)
         fnorm(3)=-cos(del)
         slip(1)= cos(lam)*cos(phi)+cos(del)*sin(lam)*sin(phi)
         slip(2)= cos(lam)*sin(phi)-cos(del)*sin(lam)*cos(phi)
         slip(3)=-sin(lam)*sin(del)
      else
         if ((1.-abs(fnorm(3))).le.1e-7) then
           print *,'***FPCOOR warning, horz fault, strike undefined'
           del=0.
           phi=atan2(-slip(1),slip(2))
           clam=cos(phi)*slip(1)+sin(phi)*slip(2)
           slam=sin(phi)*slip(1)-cos(phi)*slip(2)
           lam=atan2(slam,clam)
         else
           phi=atan2(-fnorm(1),fnorm(2))
           a=sqrt(fnorm(1)*fnorm(1)+fnorm(2)*fnorm(2))
           del=atan2(a,-fnorm(3))
           clam=cos(phi)*slip(1)+sin(phi)*slip(2)
           slam=-slip(3)/sin(del)
           lam=atan2(slam,clam)
           if (del.gt.(0.5*pi)) then
             del=pi-del
             phi=phi+pi
             lam=-lam
           end if
         end if
         strike=phi*degrad
         if (strike.lt.0.) strike=strike+360.
         dip=del*degrad
         rake=lam*degrad
         if (rake.le.-180.) rake=rake+360.
         if (rake.gt.180.) rake=rake-360.
      end if
      return
      end


c ------------------------------------------------------------ c

c normally-distributed random numbers, from numerical recipes      
      subroutine RAN_NORM(fran)
      save jran,ifirst
      im=120050                          !overflow at 2**28
      ia=2311
      ic=25367
      if (ifirst.ne.12345) then
         jran=314159
         ifirst=12345
      end if
      fran=0
      do 10 i=1,12
        jran=mod(jran*ia+ic,im)
        fran=fran+(float(jran)/float(im))
10    continue
      fran=fran-6.
      return
      end

c MK_TABLES creates tables of takeoff angles given 1D velocity models.

c   output:
c     ntab - number of tables (max nindex)
c
c   you are prompted for the names of the 1D velocity model files,
c   velocity file format (free format):
c     depth(km) P_velocity(km/s)
      
      subroutine MK_TABLE(ntab)

c      include 'vel.inc'
c nx0 = maximum source-station distance bins for look-up tables
c nd0 = maximum source depth bins for look-up tables
c nindex = maximum number of velocity model look-up tables
c dep1,dep2,dep3 = minimum source depth, maximum, and interval
c del1,del2,del3 = minimum source-station distance, maximum, and interval
c pmin = minimum ray parameter for ray tracing
c nump = number of rays traced
      parameter(nx0=101,nd0=14,nindex=10)
      parameter(dep1=0,dep2=39,dep3=3)
      parameter(del1=0,del2=200,del3=2)
      parameter(pmin=0,nump=9000)

      real table(nx0,nd0,nindex),delttab(nx0),deptab(nd0)
      integer ndel,ndep
      
      common /angtable/ table,delttab,deptab,ndel,ndep
c  common block:
c    table(nx0,nd0,nindex)  =  takeoff angle table
c        delttab(nx0)  =  list of ranges for tables
c         deptab(nd0)  =  list of source depths for tables
c           ndel       =  number of distance points in table
c           ndep       =  number of source depths in table

      parameter (nray0=10001)
      real z(1000),alpha(1000),slow(1000)
      character vmodel*100
      real deltab(nray0),tttab(nray0),ptab(nray0),tt(nray0,nd0)
      real depxcor(nray0,nd0),depucor(nray0,nd0),deptcor(nray0,nd0)
      real xsave(20000),tsave(20000),psave(20000),usave(20000)
c
      degrad=180./3.14159265

      if (ntab.ne.1) then
        print *,'Enter number of velocity models (max ',nindex,')'
        read *,ntab
      end if
      
      do 300 itab=1,ntab
      
      print *,'Enter file name for velocity model ',itab
      read (*,'(a)') vmodel

c set up table
      qtempdep2=dep2+dep3/20.
      ndep=int((qtempdep2-dep1)/dep3)+1
      do idep=1,ndep
         dep=dep1+dep3*real(idep-1)
         deptab(idep)=dep
      end do

c read velocity model
      open (7,file=vmodel,status='old')
      do i=1,1000
         read (7,*,end=30) z(i),alpha(i)
      end do
      print *,'***1000 point maximum exceeded in model'
30    close (7)
38    z(i)=z(i-1)           
      alpha(i)=alpha(i-1)
      npts=i
      npts_old=npts
      do i=npts_old,2,-1
        do idep=ndep,1,-1
          if ((z(i-1).le.(deptab(idep)-0.1)).and.
     &        (z(i).ge.(deptab(idep)+0.1))) then
            npts=npts+1
            do j=npts,i+1,-1
              z(j)=z(j-1)
              alpha(j)=alpha(j-1)
            end do
            z(i)=deptab(idep)
            frac=(z(i)-z(i-1))/(z(i+1)-z(i-1))
            alpha(i)=alpha(i-1)+frac*(alpha(i+1)-alpha(i-1))
          end if
        end do
      end do
      do i=1,npts
         slow(i)=1./alpha(i)
      end do
      pmax=slow(1)
      plongcut=slow(npts)
      pstep=(pmax-pmin)/float(nump)


c do P-wave ray tracing
       npmax=int((pmax+pstep/2.-pmin)/pstep)+1
       do 200 np=1,npmax
         p=pmin+pstep*real(np-1)
         ptab(np)=p
         x=0.
         t=0.
         imth=3
         do 70 idep=1,ndep
            if (deptab(idep).eq.0.) then
               depxcor(np,idep)=0.
               deptcor(np,idep)=0.
               depucor(np,idep)=slow(1)
            else
               depxcor(np,idep)=-999.
               deptcor(np,idep)=-999.
               depucor(np,idep)=-999.
            end if
70       continue
         do 100 i=1,npts-1
           if (z(i).ge.9999) then 
              deltab(np)=-999.
              tttab(np)=-999.
              go to 200
           end if
           h=z(i+1)-z(i)
           if (h.eq.0.) go to 100          !skip if interface
           call LAYERTRACE(p,h,slow(i),slow(i+1),imth,dx,dt,irtr)
           x=x+dx
           t=t+dt
           if (irtr.eq.0.or.irtr.eq.2) go to 105   !ray has turned
           xdeg=x                       ! actually in km 
           tmin=t                       ! actually in s
           do 80 idep=1,ndep
            if (abs(z(i+1)-deptab(idep)).lt.0.1) then
               depxcor(np,idep)=xdeg
               deptcor(np,idep)=tmin
               depucor(np,idep)=slow(i+1)            
            end if
80         continue
100      continue
105      xdeg=2.*x                     ! actually in km 
         tmin=2.*t                     ! actually in s
110      deltab(np)=xdeg
         tttab(np)=tmin
200   continue    !  end loop on ray parameter p

c create table
      do 250 idep=1,ndep
         icount=0
         xold=-999.
         if (deptab(idep).eq.0.) then
            i2=np
            go to 223
         end if
         do 220 i=1,np                    !upgoing rays from source
            x2=depxcor(i,idep)
            if (x2.eq.-999.) go to 221
            if (x2.le.xold) go to 221     !stop when heads inward
            t2=deptcor(i,idep)
            icount=icount+1
            xsave(icount)=x2
            tsave(icount)=t2
            psave(icount)=-ptab(i)
            usave(icount)=depucor(i,idep)
            xold=x2
220      continue
221      continue
         i2=i-1
223      do 225 i=i2,1,-1                 !downgoing rays from source
            if (depxcor(i,idep).eq.-999.) go to 225
            if (deltab(i).eq.-999.) go to 225
            x2=deltab(i)-depxcor(i,idep)
            t2=tttab(i)-deptcor(i,idep)
            icount=icount+1
            xsave(icount)=x2
            tsave(icount)=t2
            psave(icount)=ptab(i)
            usave(icount)=depucor(i,idep)
            xold=x2
225      continue
226      ncount=icount

         ndel=int((del2-del1)/del3)+1
         do 240 idel=1,ndel
            del=del1+del3*real(idel-1)
            delttab(idel)=del
            tt(idel,idep)=999.
            do 230 i=2,ncount
               x1=xsave(i-1)
               x2=xsave(i)
               if (x1.gt.del.or.x2.lt.del) go to 230
               if (psave(i).gt.0..and.psave(i).lt.plongcut) go to 230
               frac=(del-x1)/(x2-x1)
               t1=tsave(i-1)+frac*(tsave(i)-tsave(i-1))
               if (t1.lt.tt(idel,idep)) then
                  tt(idel,idep)=t1
                  scr1=psave(i)/usave(i)
                  angle=asin(scr1)*degrad
                  if (angle.lt.0.) then
                     angle=-angle
                  else
                     angle=180.-angle
                  end if
                  table(idel,idep,itab)=angle
               end if
230         continue
240      continue


250   continue
      if (delttab(1).eq.0.) then
         do idep=1,ndep
            table(1,idep,itab)=0.        !straight up at zero range
         end do
      end if
300   continue

c      do idel=1,ndel
c        do idep=1,ndep
c          print *,idel,delttab(idel),idep,deptab(idep),
c     &     table(idel,idep,1)
c        end do
c      end do
      
999   return
      end


c ------------------------------------------------------------ c


c subroutine GET_TTS obtains the takeoff angle for a velocity model
c at a specified range and earthquake depth by interpolating
c from a table of takeoff angles.  
c    Inputs:    ip     =  index number for model (up to nindex)
c               del    =  range
c               qdep   =  earthquake depth
c    Returns:   tt     =  takeoff angle (degrees)
c               iflag  = -1 if outside depth range
c                      =  0 for interpolation
c                      =  1 for extrapolation in range
c
      subroutine GET_TTS(ip,del,qdep,tt,iflag)

c      include 'vel.inc'

c nx0 = maximum source-station distance bins for look-up tables
c nd0 = maximum source depth bins for look-up tables
c nindex = maximum number of velocity model look-up tables
c dep1,dep2,dep3 = minimum source depth, maximum, and interval
c del1,del2,del3 = minimum source-station distance, maximum, and interval
c pmin = minimum ray parameter for ray tracing
c nump = number of rays traced
      parameter(nx0=101,nd0=14,nindex=10)
      parameter(dep1=0,dep2=39,dep3=3)
      parameter(del1=0,del2=200,del3=2)
      parameter(pmin=0,nump=9000)

      real t(nx0,nd0,nindex),x(nx0),d(nd0)
      integer nx,nd
      
      common /angtable/ t,x,d,nx,nd
c  common block:
c    t(nx0,nd0,nindex)  =  takeoff angle tables
c          x(nx0)  =  list of ranges for tables
c          d(nd0)  =  list of source depths for tables
c           nx     =  number of distance points in table
c           nd     =  number of source depths in table

c
c check if outside depth range
      if (qdep.lt.d(1).or.qdep.gt.d(nd0)) then
         iflag=-1
         tt=999
         print *,'*** event outside of velocity table depth range, 
     &           event depth=',qdep,' table range=',d(1),d(nd0)
         return
      end if
c first check to see if interpolation alone will work
      do 30 id=2,nd
         if (d(id).lt.qdep) go to 30
         id1=id-1
         id2=id
         go to 32
30    continue
      id1=nd-1
      id2=nd
32    do 35 ix=2,nx
         if (x(ix).lt.del) go to 35
         ix1=ix-1
         ix2=ix
         go to 37
35    continue
      ix1=nx-1
      ix2=nx
37    if (t(ix1,id1,ip).eq.0.) go to 50
      if (t(ix1,id2,ip).eq.0.) go to 50
      if (t(ix2,id1,ip).eq.0.) go to 50
      if (t(ix2,id2,ip).eq.0.) go to 50
      if (x(ix2).lt.del) go to 50
      iflag=0
      xfrac=(del-x(ix1))/(x(ix2)-x(ix1))
      t1=t(ix1,id1,ip)+xfrac*(t(ix2,id1,ip)-t(ix1,id1,ip))
      t2=t(ix1,id2,ip)+xfrac*(t(ix2,id2,ip)-t(ix1,id2,ip))
      dfrac=(qdep-d(id1))/(d(id2)-d(id1))
      tt=t1+dfrac*(t2-t1)
      return
c extrapolate to get tt
50    iflag=1
      xoffmin1=999.
      xoffmin2=999.
      ixbest1=999
      ixbest2=999
      do 60 ix=2,nx
         if (t(ix-1,id1,ip).eq.0) go to 55
         if (t(ix,id1,ip).eq.0) go to 55
         xoff=abs((x(ix-1)+x(ix))/2.-del)
         if (xoff.lt.xoffmin1) then
            xoffmin1=xoff
            ixbest1=ix
         end if
55       if (t(ix-1,id2,ip).eq.0) go to 60
         if (t(ix,id2,ip).eq.0) go to 60
         xoff=abs((x(ix-1)+x(ix))/2.-del)
         if (xoff.lt.xoffmin2) then
            xoffmin2=xoff
            ixbest2=ix
         end if
60    continue
      if (ixbest1.eq.999.or.ixbest2.eq.999) then
         iflag=-1
         tt=999
         return
      end if
      xfrac1=(del-x(ixbest1-1))/(x(ixbest1)-x(ixbest1-1))
      t1=t(ixbest1-1,id1,ip)
      t2=t(ixbest1,id1,ip)
      tt1=t1+xfrac1*(t2-t1)
      xfrac2=(del-x(ixbest2-1))/(x(ixbest2)-x(ixbest2-1))
      t1=t(ixbest2-1,id2,ip)
      t2=t(ixbest2,id2,ip)
      tt2=t1+xfrac2*(t2-t1)
      dfrac=(qdep-d(id1))/(d(id2)-d(id1))
      tt=tt1+dfrac*(tt2-tt1)
999   return
      end


c ------------------------------------------------------------ c

c LAYERTRACE calculates the travel time and range offset
c for ray tracing through a single layer.
c
c Input:    p     =  horizontal slowness
c           h     =  layer thickness
c           utop  =  slowness at top of layer
c           ubot  =  slowness at bottom of layer
c           imth  =  interpolation method
c                    imth = 1,  v(z) = 1/sqrt(a - 2*b*z)
c                         = 2,  v(z) = a - b*z
c                         = 3,  v(z) = a*exp(-b*z)
c
c Returns:  dx    =  range offset
c           dt    =  travel time
c           irtr  =  return code
c                 = -1, zero thickness layer
c                 =  0,  ray turned above layer
c                 =  1,  ray passed through layer
c                 =  2,  ray turned within layer, 1 segment counted
c
c Note:  This version does calculation in double precision,
c        but all i/o is still single precision
c
      subroutine LAYERTRACE(p1,h1,utop1,ubot1,imth,dx1,dt1,irtr)
      implicit real*8 (a-h,o-z)
      real*4 p1,h1,utop1,ubot1,dx1,dt1
      p=dble(p1)
      h=dble(h1)
      utop=dble(utop1)
      ubot=dble(ubot1)
c
      if (h.eq.0.) then      !check for zero thickness layer
         dx1=0.
         dt1=0.
         irtr=-1
         return         
      end if
c
      u=utop
      y=u-p
      if (y.le.0.) then   !complex vertical slowness
         dx1=0.
         dt1=0.
         irtr=0
         return
      end if
c
      q=y*(u+p)
      qs=dsqrt(q)
c
c special function needed for integral at top of layer
      if (imth.eq.2) then
         y=u+qs
         if (p.ne.0.) y=y/p
         qr=dlog(y)
      else if (imth.eq.3) then
         qr=atan2(qs,p)
      end if      
c
      if (imth.eq.1) then
          b=-(utop**2-ubot**2)/(2.*h)
      else if (imth.eq.2) then
          vtop=1./utop
          vbot=1./ubot
          b=-(vtop-vbot)/h
      else
          b=-dlog(ubot/utop)/h
      end if  
c
      if (b.eq.0.) then     !constant velocity layer
         b=1./h
         etau=qs
         ex=p/qs
         irtr=1
         go to 160
      end if
c
c integral at upper limit, 1/b factor omitted until end
      if (imth.eq.1) then
         etau=-q*qs/3.
         ex=-qs*p
      else if (imth.eq.2) then
         ex=qs/u                 !*** - in some versions (wrongly)
         etau=qr-ex
         if (p.ne.0.) ex=ex/p
      else
         etau=qs-p*qr
         ex=qr
      end if
c
c check lower limit to see if we have turning point
      u=ubot
      if (u.le.p) then   !if turning point,
         irtr=2          !then no contribution
         go to 160       !from bottom point
      end if 
      irtr=1
      q=(u-p)*(u+p)
      qs=dsqrt(q)
c
      if (imth.eq.1) then
         etau=etau+q*qs/3.
         ex=ex+qs*p
      else if (imth.eq.2) then
         y=u+qs
         z=qs/u
         etau=etau+z
         if (p.ne.0.) then
            y=y/p
            z=z/p
         end if
         qr=dlog(y)
         etau=etau-qr
         ex=ex-z
      else
         qr=atan2(qs,p)
         etau=etau-qs+p*qr
         ex=ex-qr
      end if      
c
160   dx=ex/b
      dtau=etau/b
      dt=dtau+p*dx     !convert tau to t
c
      dx1=sngl(dx)
      dt1=sngl(dt)
      return
      end


      function aacos(x)
      real x
      if(x.lt.-1.0) then
c         write(6,*) x
         x=-1.0
      endif
      if(x.gt.1.0) then
c         write(6,*) x
         x=1.0
      endif
      aacos=acos(x)
      return
      end  





