C ---------------------------------------------------------------------------------
C Add shear fracture
C 2024-12-11, lizt, Tsinghua
C ---------------------------------------------------------------------------------

      module vars_module
        parameter (NodeNum = 70023, NumEle= 65850, NPT=1)
        real*8,save :: allD(NodeNum), allH(NumEle), allDP(NumEle), allH_glo(NumEle),allD_glo(NodeNum)
        real*8,save :: allG(NumEle), allG_glo(NumEle)
        integer, save :: NUMPROCESSES = 1
      end module

      subroutine vexternaldb(lOp, i_Array, niArray, r_Array, nrArray)
C
      use vars_module
      include 'vaba_param.inc'
      include 'mpif.h'

C     Contents of i_Array
      parameter( i_int_nTotalNodes     = 1,
     *           i_int_nTotalElements  = 2,
     *           i_int_kStep           = 3,
     *           i_int_kInc            = 4,
     *           i_int_iStatus         = 5,
     *           i_int_lWriteRestart   = 6  )

C     Possible values for the lOp argument
      parameter( j_int_StartAnalysis    = 0,      
     *           j_int_StartStep        = 1,      
     *           j_int_SetupIncrement   = 2,      
     *           j_int_StartIncrement   = 3,      
     *           j_int_EndIncrement     = 4,      
     *           j_int_EndStep          = 5,      
     *           j_int_EndAnalysis      = 6 )     


C     Possible values for i_Array(i_int_iStatus)
      parameter( j_int_Continue          = 0,      
     *           j_int_TerminateStep     = 1,      
     *           j_int_TerminateAnalysis = 2)

C     Contents of r_Array
      parameter( i_flt_TotalTime   = 1,
     *           i_flt_StepTime    = 2,
     *           i_flt_dTime       = 3 )
C
      dimension i_Array(niArray),      
     *   r_Array(nrArray)

      
      if (lOp == 0)then
          call VGETNUMCPUS( NUMPROCESSES )
      endif
  
      if(lOp == j_int_EndIncrement)then
        if(NUMPROCESSES > 1)then
          call MPI_Allreduce(allH, allH_glo, NumEle, MPI_doUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD,ierr)
          call MPI_Allreduce(allD, allD_glo, NodeNum, MPI_doUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD,ierr)
          call MPI_Allreduce(allG, allG_glo, NumEle, MPI_doUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD,ierr)
        else
          allH_glo = allH
          allD_glo = allD
          allG_glo = allG
        endif
      endif
      
      if(lOp == j_int_StartIncrement)then
C initialize allG to zero
        allG(:) = 0.0d0
      endif


      return
      end subroutine

      
C ---------------------------------------------------------------------------------
C Vumat interface
C ---------------------------------------------------------------------------------

      subroutine vumat(
     1  jblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2  stepTime, totalTime, dt, cmname, coordMp, charLength,
     3  props, density, strainInc, relSpinInc,
     4  tempOld, stretchOld, defgradOld, fieldOld,
     5  stressOld, stateOld, enerInternOld, enerInelasOld,
     6  tempNew, stretchNew, defgradNew, fieldNew,
     7  stressNew, stateNew, enerInternNew, enerInelasNew )
C
      include 'vaba_param.inc'
C
      dimension jblock(*), props(nprops),density(*), coordMp(*),
     1     charLength(*), strainInc(*),
     2     relSpinInc(*), tempOld(*),
     3     stretchOld(*),
     4     defgradOld(*),
     5     fieldOld(*), stressOld(*),
     6     stateOld(*), enerInternOld(*),
     7     enerInelasOld(*), tempNew(*),
     8     stretchNew(*),
     9     defgradNew(*),
     1     fieldNew(*),
     2     stressNew(*), stateNew(*),
     3     enerInternNew(*), enerInelasNew(*)
C
      character*80 cmname
      parameter (
     *     i_umt_nblock = 1,
     *     i_umt_npt    = 2,
     *     i_umt_layer  = 3,
     *     i_umt_kspt   = 4,
     *     i_umt_noel   = 5 )
!      Local variables
!
      call vumat_pfm_elas_plas(
     1  jblock(i_umt_nblock), ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2  stepTime, totalTime, dt, cmname, coordMp, charLength,
     3  props, density, strainInc, relSpinInc,
     4  tempOld, stretchOld, defgradOld, fieldOld,
     5  stressOld, stateOld, enerInternOld, enerInelasOld,
     6  tempNew, stretchNew, defgradNew, fieldNew,
     7  stressNew, stateNew, enerInternNew, enerInelasNew,
     8  jblock(i_umt_noel) )

      return
      end subroutine

      subroutine vumat_pfm_elas_plas(
     1  nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2  stepTime, totalTime, dt, cmname, coordMp, charLength,
     3  props, density, strainInc, relSpinInc,
     4  tempOld, stretchOld, defgradOld, fieldOld,
     5  stressOld, stateOld, enerInternOld, enerInelasOld,
     6  tempNew, stretchNew, defgradNew, fieldNew,
     7  stressNew, stateNew, enerInternNew, enerInelasNew,
c Read only extra arguments -
     8  nElem )
C
      include 'vaba_param.inc'
C
      dimension props(nprops), density(nblock), coordMp(nblock,*),
     1  charLength(nblock), strainInc(nblock,ndir+nshr),
     2  relSpinInc(nblock,nshr), tempOld(nblock),
     3  stretchOld(nblock,ndir+nshr),
     4  defgradOld(nblock,ndir+nshr+nshr),
     5  fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     6  stateOld(nblock,nstatev), enerInternOld(nblock),
     7  enerInelasOld(nblock), tempNew(nblock),
     8  stretchNew(nblock,ndir+nshr),
     8  defgradNew(nblock,ndir+nshr+nshr),
     9  fieldNew(nblock,nfieldv),
     1  stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     2  enerInternNew(nblock), enerInelasNew(nblock)
C
      character*80 cmname
c
      parameter ( zero = 0.d0, one = 1.d0, two = 2.d0,d_thresh = 0.95d0,
     *     third = 1.d0 / 3.d0, half = 0.5d0, op5 = 1.5d0)

      parameter (nblkLocal = 144)
      dimension dYield(nblock), dHard(nblock), gdk(nblock)
      dimension e_tensor(3,3), e_prin(3), e_voigt(6), e_dir(3,3), alphai(3)
      dimension s_tensor(3,3), s_prin(3,3)
      real*8 phi, phiev, phiplus, phie
C --- Parameters read in
      real*8 E, nu, w0, w0s
      real*8 twomu, alamda, thremu, alphaT, module_K
      integer nvalue
      real*8 triaxiality, trix_cr, tri_delta
      real*8 phase, gd, trace, factor
      real*8 de_th_expan, de_v_3, de_v
      real*8 s11, s22, s33, s12, s13, s23, smean, vmises
      real*8 e_tr, e1plus, e2plus, e3plus
      real*8 dep11, dep22, dep33, dep12, dep13, dep23
      real*8 yieldOld, yieldNew, hard, sigdif, facyld, deqps
      real*8 plasticWorkInc
*
* --- exit if nblkLocal is smaller than nblock
*
      if (nblkLocal .lt. nblock) then
        call xplb_abqerr (-2,'Change nblkLocal to be greater '//
     *       'than or equal to %i',nblock,zero,' ')
        call xplb_exit
      end if

      E      = props(1)
      nu    = props(2)
C w0T default value is 0.0
      w0s    = props(3)
      twomu  = E / ( one + nu )
      alamda = nu * twomu / ( one - two * nu )
      thremu = op5 * twomu
      module_K = E/3.0d0/(one - two * nu)
      alphaT = zero
      
      nvalue = nprops/2-1

      if ( stepTime .eq. zero ) then
        do k = 1, nblock
          stateNew(k,16) = one
          trace = strainInc(k,1) + strainInc(k,2) + strainInc(k,3)
          stressNew(k,1) = stressOld(k,1) 
     *         + twomu * strainInc(k,1) + alamda * trace
          stressNew(k,2) = stressOld(k,2) 
     *         + twomu * strainInc(k,2) + alamda * trace
          stressNew(k,3) = stressOld(k,3) 
     *         + twomu * strainInc(k,3) + alamda * trace
          stressNew(k,4) = stressOld(k,4) + twomu * strainInc(k,4)
        end do
        if ( nshr .gt. 1 ) then
          do k = 1, nblock
            stressNew(k,5) = stressOld(k,5) + twomu * strainInc(k,5)
            stressNew(k,6) = stressOld(k,6) + twomu * strainInc(k,6)
          end do
        end if
        return
      end if


      do k=1 , nblock
C          gdk(k) = (one-fieldOLd(k,1))**2
            gdk(k) = one
      enddo
C      print*, 'nprops = ', props(4:nprops)
      call xhard (nblock, nvalue, stateOld(:,17), dYield, dHard, props(4:nprops), gdk)

      if ( nshr .gt. 1 ) then         ! 3D case
      
        do k = 1, nblock
          twomu  = e / ( one + xnu )
          alamda = xnu * twomu / ( one - two * xnu )
          thremu = op5 * twomu

C          phase = fieldOLd(k,1)
          phase = zero
          if (phase >= d_thresh)stateNew(k,16) = zero
          if (phase <= zero)phase = zero   
          if (phase >= d_thresh)phase = d_thresh
           
          de_th_expan = alphaT*(tempNew(k)-tempOld(k))
          de_v_3 = de_th_expan
          de_v = 3.*de_v_3

          gd = (one-phase)**2
                 
          trace = strainInc(k,1) + strainInc(k,2) + strainInc(k,3) - de_v
          
           s11 = stressOld(k,1) + twomu * strainInc(k,1) + alamda * trace
           s22 = stressOld(k,2) + twomu * strainInc(k,2) + alamda * trace
           s33 = stressOld(k,3) + twomu * strainInc(k,3) + alamda * trace
           s12 = stressOld(k,4) + twomu * strainInc(k,4)
           s23 = stressOld(k,5) + twomu * strainInc(k,5)
           s13 = stressOld(k,6) + twomu * strainInc(k,6)

c     
            smean = third * ( s11 + s22 + s33 ) 
            s11 = s11 - smean
            s22 = s22 - smean
            s33 = s33 - smean
c
          vmises = sqrt( op5 * ( s11 * s11 + s22 * s22 + s33 * s33 +
     *         two * s12 * s12 + two * s13 * s13 + two * s23 * s23 ) )
c

          yieldOld = dYield(k)
          hard = dHard(k)
          sigdif = vmises - yieldOld!*gd
          facyld = half + sign(half,sigdif) ! first value, second sign
          deqps = facyld * sigdif / ( thremu + hard )

            trix_cr = 0.7d0
            tri_delta = 0.2d0

          if (vmises > 0)then
            dep11 = sqrt(op5)*deqps*s11/vmises
            dep22 = sqrt(op5)*deqps*s22/vmises
            dep33 = sqrt(op5)*deqps*s33/vmises
            dep12 = sqrt(op5)*deqps*s12/vmises
            dep23 = sqrt(op5)*deqps*s23/vmises
            dep13 = sqrt(op5)*deqps*s13/vmises
c Stress triaxiality
            triaxiality = smean / vmises
          else
            dep11 = 0.d0
            dep22 = 0.d0
            dep33 = 0.d0
            dep12 = 0.d0
            dep23 = 0.d0
            dep13 = 0.d0
            triaxiality = 1.0d5
          endif
            stateNew(k,20) = triaxiality
            stateNew(k,21) = tanh((triaxiality-trix_cr)/tri_delta)

c
c Update the stress
c
          yieldNew = yieldOld + hard * deqps 
          factor = yieldNew / ( yieldNew + thremu * deqps )

          stressNew(k,1) = s11 * factor + smean
          stressNew(k,2) = s22 * factor + smean
          stressNew(k,3) = s33 * factor + smean
          stressNew(k,4) = s12 * factor
          stressNew(k,5) = s23 * factor
          stressNew(k,6) = s13 * factor

c
c Update the state variables
c
          stateNew(k,17) = stateOld(k,17) + deqps

!c Update the dissipated inelastic specific energy
          plasticWorkInc = half * ( yieldOld + yieldNew ) * deqps
          
          stateNew(k,18) = stateOld(k,18) + plasticWorkInc
          stateNew(k,19) = stateOld(k,19) + plasticWorkInc/gd
          

          !stateNew(k,1:6) = stateOld(k,1:6) + strainInc(k,1:6)
          stateNew(k,1) = stateOld(k,1) + strainInc(k,1) - de_v_3 - dep11
          stateNew(k,2) = stateOld(k,2) + strainInc(k,2) - de_v_3 - dep22
          stateNew(k,3) = stateOld(k,3) + strainInc(k,3) - de_v_3 - dep33
          stateNew(k,4) = stateOld(k,4) + strainInc(k,4) - dep12
          stateNew(k,5) = stateOld(k,5) + strainInc(k,5) - dep23
          stateNew(k,6) = stateOld(k,6) + strainInc(k,6) - dep13

          e_voigt = stateNew(k,1:6)

          call voigt_convection(e_voigt, e_tensor,.false.,.true.)

          call eig3(e_tensor, e_prin, e_dir)
           
          stateNew(k,7:9) = e_prin
          stateNew(k,10) = de_th_expan

          e_tr = e_prin(1) + e_prin(2) + e_prin(3)
          e1plus = max(e_prin(1),0.0)
          e2plus = max(e_prin(2),0.0)
          e3plus = max(e_prin(3),0.0)           
           
          phie = half*(alamda*e_tr**2 + twomu*(e_prin(1)**2 + e_prin(2)**2 + e_prin(3)**2) )
          if(e_tr < 0.0) then
              phiev = phie - half*module_K*e_tr**2
          else
              phiev = phie
          endif
          phieplus = half*(alamda*max(e_tr,0.0)**2 + twomu*(e1plus**2 + e2plus**2 + e3plus**2) )
          
          
          alpha=zero
          if (e_tr .gt. zero) alpha=one

          do K1=1,3
              alphai(K1)=zero
              if (e_prin(K1) .gt. zero) alphai(K1)=one
          enddo
          
          alpha = zero
          alphai = zero

          sigma1 = twomu*e_prin(1)*(1-alphai(1)*phase)**2 + alamda*e_tr*(1-alpha*phase)**2
          sigma2 = twomu*e_prin(2)*(1-alphai(2)*phase)**2 + alamda*e_tr*(1-alpha*phase)**2
          sigma3 = twomu*e_prin(3)*(1-alphai(3)*phase)**2 + alamda*e_tr*(1-alpha*phase)**2

          s_prin = 0.d0
          s_prin(1,1) = sigma1
          s_prin(2,2) = sigma2
          s_prin(3,3) = sigma3

          s_tensor = matmul(matmul(e_dir,s_prin),transpose(e_dir))
          
          call voigt_convection(stressNew(k,:), s_tensor,.true.,.true.)

c
c Update the state variables
c
c Update the specific internal energy
          stressPower = half * (
     *      ( stressOld(k,1) + stressNew(k,1) ) * strainInc(k,1) +
     *      ( stressOld(k,2) + stressNew(k,2) ) * strainInc(k,2) +
     *      ( stressOld(k,3) + stressNew(k,3) ) * strainInc(k,3) ) +
     *      ( stressOld(k,4) + stressNew(k,4) ) * strainInc(k,4) +    
     *      ( stressOld(k,5) + stressNew(k,5) ) * strainInc(k,5) +    
     *      ( stressOld(k,6) + stressNew(k,6) ) * strainInc(k,6)     

          enerInternNew(k) = enerInternOld(k) + stressPower / density(k)

c     
c Update the dissipated inelastic specific energy
          stateNew(k,11) = phie
          stateNew(k,12) = phiev
          stateNew(k,13) = phieplus
          ! Select model A, B or C by modifying them here
          !stateNew(k,14) = max( stateOld(k,14), stateNew(k,13) ) ! H
C w0 = 0.5*w0s - 0.5*w0s*tanh((triaxiality-trix_cr)/tri_delta)
          w0 = half*w0s - half*w0s*stateNew(k,21)
          if (phase < d_thresh)then
              stateNew(k,14) = max( stateOld(k,14),  max(stateNew(k,13) + 0.1*stateNew(k,19) - 0.1*w0, 0.0) ) ! H
          else
              stateNew(k,14) = stateOld(k,14)
          endif
          
C          stateNew(k,15) = fieldOLd(k,1)  ! d
          stateNew(k,15) = zero  ! d
          
          enerInelasNew(k) = stateNew(k,18)/density(k)

        end do

      else if ( nshr .eq. 1 ) then     ! plane strain/axisymmetric case

        write(*,*)"We have not implement 2D cases!"
        call xplb_exit

      end if

      return
      end subroutine


      subroutine xhard(nblock, nvalue, peeqOld, syield, hard, table_one, gdk)
      include 'vaba_param.inc'

c
      dimension peeqOld(nblock), syield(nblock), hard(nblock), gdk(nblock),
     *     table_one(2*nvalue), table(2, nvalue)
c
      parameter(zero=0.d0)
c
c    set yield stress to last value of table, hardening to zero
c
C Reshape table_one to table, table(1, k) is the yield stress, table(2, k) is the strain
      do k = 1, nvalue
C      print*, 'sigma_y = ', table_one(2*k-1)
C      print*, 'strain = ', table_one(2*k)
      table(1, k) = table_one(2*k-1)
      table(2, k) = table_one(2*k)
      end do

      do k = 1,nblock
        syield(k) = table(1, nvalue)*gdk(k)
        hard(k) = zero
c C
c         print*, 'yield stress: ', syield(k)
c         print*, 'hardening: ', hard(k)
c         print*, 'table 1: ', table(1, 1)
c         print*, 'table 2: ', table(2, 1)
      end do
c    if more than one entry, search table
c
      if(nvalue.gt.1) then
        do k = 1, nblock
c equivalent plastic strain
          eqplas = peeqOld(k)
          interval_search: do k1 = 1, nvalue-1
            eqpl1 = table(2,k1+1)
            if(eqplas.lt.eqpl1) then
              eqpl0 = table(2, k1)
c
c          yield stress and hardening
c
              deqpl = eqpl1-eqpl0
              syiel0 = table(1, k1)*gdk(k)
              syiel1 = table(1, k1+1)*gdk(k)
              dsyiel = syiel1-syiel0
              hard(k) = dsyiel/deqpl
              syield(k) = syiel0+(eqplas-eqpl0)*hard(k)
              exit interval_search
            endif
          end do interval_search
        end do
      endif

      return
      end


C***********************************************************************
C
C     VUEL for crack propagation by PFM
C     doF 11 represent the phase d
C
C***********************************************************************
      subroutine vuel(
     *     nblock,
c          to be defined
     *     rhs,amass,dtimeStable,
     *     svars,nsvars,
     *     energy,
c          
     *     nnode,ndofel,
     *     props,nprops,
     *     jprops,njprops,
     *     coords,ncrd,
     *     u,du,v,a,
     *     jtype,jelem,
     *     time,period,dtimeCur,dtimePrev,kstep,kinc,lflags,
     *     dMassScaleFactor,
     *     predef,npredef,
     *     jdltyp,adlmag)

      use vars_module
C     
      include 'vaba_param.inc'

      parameter ( zero = 0.d0, half = 0.5d0, one = 1.d0, two=2.d0 )
      parameter (scaleTemp = 0.9d0)

c     operation code
      parameter ( jMassCalc            = 1,
     *            jIntForceAndDtStable = 2,
     *            jExternForce         = 3)

c     flags
      parameter (iProcedure = 1,
     *           iNlgeom    = 2,
     *           iOpCode    = 3,
     *           nFlags     = 3)

c     time
      parameter (iStepTime  = 1,
     *           iTotalTime = 2,
     *           nTime      = 2)

c     procedure flags
      parameter ( jDynExplicit = 17 )

c     energies 
      parameter ( iElPd = 1,
     *            iElCd = 2,
     *            iElIe = 3,
     *            iElTs = 4,
     *            iElDd = 5,
     *            iElBv = 6,
     *            iElDe = 7,
     *            iElHe = 8,
     *            iElKe = 9,
     *            iElTh = 10,
     *            iElDmd = 11,
     *            iElDc = 12,
     *            nElEnergy = 12)


c     predefined variables
      parameter ( iPredValueNew = 1,
     *            iPredValueOld = 2,
     *            nPred         = 2)    

      parameter (factorStable = 0.99d0)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      dimension rhs(nblock,ndofel), amass(nblock,ndofel,ndofel),
     *     dtimeStable(nblock),
     *     svars(nblock,nsvars), energy(nblock,nElEnergy),
     *     props(nprops), jprops(njprops),
     *     jelem(nblock), time(nTime), lflags(nFlags),
     *     coords(nblock,nnode,ncrd), u(nblock,ndofel),
     *     du(nblock,ndofel), v(nblock,ndofel), a(nblock, ndofel),
     *     dMassScaleFactor(nblock),
     *     predef(nblock, nnode, npredef, nPred), adlmag(nblock)
      
!     Declaration of variables for user element
      parameter(mone=-1.d0,three=3.d0,d_thresh=0.95d0)
      
      integer i,j,L,k,K1,K2,K3,K4,IX,IY,IZ

      dimension Awt(NPT),XII(NPT,3),XI(3),dNdxi(nnode,3),
     1 VJ(3,3),dNdx(nnode,3),VJ_inv(3,3),AN(nnode),BP(3,8),dp(3)
   
      dimension UD(8),DUD(8)
      dimension sg(4,NPT),ss(3),e_coord(3,8),shp(4,8),bn(8),cvn(8),vn(8)
      dimension cmass(8,8)
      dimension ID(8),IZ(16),IU(24)
      real*8 det, hist, eta, Gc, alc, WF, WO
      real*8 GcT, GcS, Gf, w0s
      data ID/1,2,3,4,5,6,7,8/
      logical,save :: firstcall = .true.
      
C      
!     begin the VUEL      
      if (jtype .eq. 2) then 
         ! determine the gauss locations and weights
         call make_quadrature(sg)
        
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                     -----------amass------------
c
         if ( lflags(iOpCode).eq.jMassCalc ) then
             call random_seed()
             call random_number(svars(1:nblock,1))
             do kblock=1,nblock
             !     the coordinate of element nodes
                do i=1,nnode
                   do j=1,ncrd
                      e_coord(j,i)=coords(kblock,i,j)
                   enddo
                enddo
                           
               cmass=0.0
               do i=1,NPT
                  do j=1,3
                     ss(j)=sg(j,i)
                  enddo 
                  call shapef3D(ss,nnode,e_coord,shp,dj)
                  we=sg(4,i)*dj
                  !calculate the mass matrix of phase field              
                  call cal_cmass(shp,we,cmass,props,nprops)          
               enddo
               
               !assemble mass matrix 
               do i=1,8
                   do j=1,8
                     amass(kblock,ID(i),ID(i))=amass(kblock,ID(i),ID(i))
     1                 +cmass(i,j)
                   enddo
               enddo
           enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                         --------rhs----------
c   
      else if ( lflags(iOpCode) == jIntForceAndDtStable) then
C
              !Material parameters  
              eta=props(1)
              alc=props(2)
              Gc = props(3)
              w0s = props(4)
              GcT = Gc
              GcS = Gc*10.0d0
              
              do kblock = 1, nblock
C calculate the local Gf according to stress trixiality
                 Gf = GcS - half*(GcS-GcT) - half*(GcS-GcT)*allG_glo(jelem(kblock) - NumEle)
                 WF = Gf/2.0/alc
                 WO = half*w0s - half*w0s*allG_glo(jelem(kblock) - NumEle)
C w0 = w0/(1-chi), chi - work-heat convertion factor
                 WO = WO*10.0d0

                 
              !  the coordinate of element nodes
                 do i = 1, nnode
                    do j = 1, ncrd
                       e_coord(j,i) = coords(kblock,i,j)
                    end do
                 end do
c                  
                !  seperate displacements 
                 do i = 1, nnode
                    UD(i) = U(kblock,ID(i))
                    DUD(i) = DU(kblock,ID(i))  
                 end do
C
                 !Local coordinates and weights
                 if(NPT==1)then
                   XII(1,1) = zero
                   XII(1,2) = zero
                   XII(1,3) = zero
                   Awt(1) = 8.d0
                 else                  
                   XII(1,1) = mone/three**half
                   XII(1,2) = mone/three**half
                   XII(1,3) = mone/three**half
                   XII(2,1) = one/three**half
                   XII(2,2) = mone/three**half
                   XII(2,3) = mone/three**half
                   XII(3,1) = one/three**half
                   XII(3,2) = one/three**half
                   XII(3,3) = mone/three**half
                   XII(4,1) = mone/three**half
                   XII(4,2) = one/three**half
                   XII(4,3) = mone/three**half
                   XII(5,1) = mone/three**half
                   XII(5,2) = mone/three**half
                   XII(5,3) = one/three**half
                   XII(6,1) = one/three**half
                   XII(6,2) = mone/three**half
                   XII(6,3) = one/three**half
                   XII(7,1) = one/three**half
                   XII(7,2) = one/three**half
                   XII(7,3) = one/three**half
                   XII(8,1) = mone/three**half
                   XII(8,2) = one/three**half
                   XII(8,3) = one/three**half
                   do i = 1, NPT
                     Awt(i) = one
                   end do 

                  end if  

                  hist = allH_glo(jelem(kblock) - NumEle)

                  !hist = hist - WO
                  if(hist < zero) hist = zero 
                  
                  !if(hist > zero)write(*,*)"Hello",hist
C
                 !Calculating properties at each integration point
                 do INPT = 1, NPT
c                    
                    !Local coordinates of the integration point
                    XI(1) = XII(INPT,1)
                    XI(2) = XII(INPT,2) 
                    XI(3) = XII(INPT,3) 
                    !Shape functions and local derivatives
                    call shapefun(AN,dNdxi,XI)
                     !  Jacobian
                    do i = 1,3
                       do j = 1,3
                          VJ(i,j) = zero
                          do k = 1, nnode
                             VJ(i,j) = VJ(i,j) + e_coord(i,k)
     1                        *dNdxi(k,j)
                          end do
                       end do
                     end do
C        
                     det = VJ(1,1)*VJ(2,2)*VJ(3,3)+
     1                  VJ(1,2)*VJ(2,3)*VJ(3,1)+VJ(1,3)*
     2                  VJ(2,1)*VJ(3,2)-VJ(3,1)*VJ(2,2)*
     3                  VJ(1,3)-VJ(3,2)*VJ(2,3)*VJ(1,1)-
     4                  VJ(3,3)*VJ(2,1)*VJ(1,2)
C     
C
                    !Inverse of Jacobian
                      VJ_inv(1,1)=(VJ(2,2)*VJ(3,3)-
     1                          VJ(2,3)*VJ(3,2))/det
                      VJ_inv(1,2)=-(VJ(1,2)*VJ(3,3)-
     1                          VJ(3,2)*VJ(1,3))/det
                      VJ_inv(1,3)=(VJ(1,2)*VJ(2,3)-
     1                          VJ(1,3)*VJ(2,2))/det
                      VJ_inv(2,1)=-(VJ(2,1)*VJ(3,3)-
     1                          VJ(2,3)*VJ(3,1))/det
                      VJ_inv(2,2)=(VJ(1,1)*VJ(3,3)-
     1                          VJ(1,3)*VJ(3,1))/det
                      VJ_inv(2,3)=-(VJ(1,1)*VJ(2,3)-
     1                          VJ(1,3)*VJ(2,1))/det
                      VJ_inv(3,1)=(VJ(2,1)*VJ(3,2)-
     1                          VJ(2,2)*VJ(3,1))/det
                      VJ_inv(3,2)=-(VJ(1,1)*VJ(3,2)-
     1                          VJ(1,2)*VJ(3,1))/det
                      VJ_inv(3,3)=(VJ(1,1)*VJ(2,2)-
     1                          VJ(1,2)*VJ(2,1))/det
C        
                      !Derivatives of shape functions respect to global ccordinates
                      do k = 1, nnode
                         do i = 1,3
                            dNdx(k,i) = zero
                            do j = 1,3
                               dNdx(k,i) = dNdx(k,i) + dNdxi(k,j)*VJ_inv(j,i)
                            end do
                         end do
                      end do
C
                      !Calculating B matrix (B=LN)
                      do i=1,nnode
                          BP(1,i) = dNdx(i,1)
                          BP(2,i) = dNdx(i,2)
                          BP(3,i) = dNdx(i,3)
                      end do    
C
                      ! Nodal phase-field value
                      phase = zero
                      do i=1, nnode
                         phase = phase + AN(i)*UD(i)
                      end do
C       
                      do i=1,3
                         dp(i) = zero
                      end do
                      do i=1,3
                         do j=1, nnode
                            dp(i) = dp(i) + BP(i,j)*UD(j)
                         end do
                      end do     

                      !if(phase >= one) then
                      !    phase = one
                      !    dp = zero
                      !end if
C                  
C        
                  !Internal forces (residual vector)
                  if(phase < d_thresh)then
                    do i = 1, nnode
                       do j = 1, 3
                          rhs(kblock,ID(i)) = rhs(kblock,ID(i)) + BP(j,i)
     1                    *dp(j)*(WF-WO)*two*alc*alc*Awt(INPT)*det
                       end do
                       rhs(kblock,ID(i)) = rhs(kblock,ID(i)) + AN(i)*
     1                     Awt(INPT)*det*((two*(WF-WO)+two*hist)*phase-two*hist)
                    end do
                  !else
                  !  do i = 1, nnode
                  !     do j = 1, 3
                  !        rhs(kblock,ID(i)) = zero
                  !     end do
                  !     rhs(kblock,ID(i)) = zero
                  !  end do
                  end if
C
                  end do       

                  !dtimeStable(kblock) = factorStable*eta/(500*(WF-WO))

             end do
          end if
      end if
c
      return
      end subroutine
      

C
C User subroutine vufield
C
      subroutine VUField( 
C Write only - 
     *     rUserField, 
C Read only - 
     *     nBlock, nField, kField, nComp,
     *     kStep, kInc, jNodeUid, time, 
     *     coords, U, V, A )

      use vars_module
*
      include 'vaba_param.inc'

      dimension rUserField(nBlock,nComp,nField)
      dimension jNodeUid(nBlock), time(4), coords(3,nBlock)
      dimension U(8,nBlock), V(8,nBlock), A(8,nBlock)
*
      parameter ( i_ufld_Current   = 1, 
     *            i_ufld_Increment = 2, 
     *            i_ufld_Period    = 3, 
     *            i_ufld_Total     = 4 )
*
      parameter ( i_ufld_CoordX = 1,
     *            i_ufld_CoordY = 2,
     *            i_ufld_CoordZ = 3 )
*
      parameter ( i_ufld_SpaDisplX = 1,
     *            i_ufld_SpaDisplY = 2,
     *            i_ufld_SpaDisplZ = 3,
     *            i_ufld_RotDisplX = 4,
     *            i_ufld_RotDisplY = 5,
     *            i_ufld_RotDisplZ = 6, 
     *            i_ufld_AcoPress  = 7,
     *            i_ufld_Temp      = 8 )
*
      parameter ( i_ufld_SpaVelX   = 1,
     *            i_ufld_SpaVelY   = 2,
     *            i_ufld_SpaVelZ   = 3,
     *            i_ufld_RotVelX   = 4,
     *            i_ufld_RotVelY   = 5,
     *            i_ufld_RotVelZ   = 6,
     *            i_ufld_DAcoPress = 7,
     *            i_ufld_DTemp     = 8 )
*
      parameter ( i_ufld_SpaAccelX  = 1,
     *            i_ufld_SpaAccelY  = 2,
     *            i_ufld_SpaAccelZ  = 3,
     *            i_ufld_RotAccelX  = 4,
     *            i_ufld_RotAccelY  = 5,
     *            i_ufld_RotAccelZ  = 6, 
     *            i_ufld_DDAcoPress = 7,
     *            i_ufld_DDTemp     = 8 )

      parameter (oneHundred = 100.d0, twoHundred = 200.d0)
*
      if (kField .eq. 1) then

         do kComp = 1, nComp
            do kNod = 1, nBlock
                if(JNODEUID(kNod) > NodeNum)then
                    allD(JNODEUID(kNod)-NodeNum) = U(i_ufld_Temp,kNod)
                    rUserField(kNod,kComp,1) = U(i_ufld_Temp,kNod) 
                else
                    rUserField(kNod,kComp,1) = allD_glo(JNODEUID(kNod))
                    !j = JNODEUID(kNod)
                    !if( j==4765 .or. j==4766 .or. j==496112 .or. j==560551 .or. j==745 .or. j==746 .or. j==110576 .or. j==46137)then
                    !    rUserField(kNod,kComp,1) = 0.99
                    !endif
                
                end if
            end do
         end do
         
      end if
*
      return
      end subroutine


c
c User subroutine VUSDFLD for user-defined fields
c
      subroutine vusdfld(
c Read only -
     *   nblock, nstatev, nfieldv, nprops, ndir, nshr, 
     *   jElemUid, kIntPt, kLayer, kSecPt, 
     *   stepTime, totalTime, dt, cmname, 
     *   coordMp, direct, T, charLength, props, 
     *   stateOld, 
c Write only -
     *   stateNew, field )
c
      use vars_module
      include 'vaba_param.inc'
c
      dimension props(nprops),
     *          jElemUid(nblock), coordMp(nblock, *), 
     *          direct(nblock, 3, 3), T(nblock,3,3), 
     *          charLength(nblock),
     *          stateOld(nblock, nstatev), 
     *          stateNew(nblock, nstatev),
     *          field(nblock, nfieldv)
      character*80 cmname
c
      character*3 cData(maxblk)
      dimension jData(maxblk)
      dimension eqps(maxblk)
c
      parameter ( zero = 0.d0 )
c
      do k = 1, nblock
        if(jElemUid(k) <= NumEle)then
            allH(jElemUid(k)) = stateOld(k,14)  
            allG(jElemUid(k)) = stateOld(k,21)  
        end if
      end do

      do k = 1, nblock
        if(jElemUid(k) <= NumEle)then
            stateNew(k,15) = allDP(jElemUid(k))
        end if
      end do
c
      return
      end subroutine



      subroutine eig3(u, d, v)
      !=======================================================================
      !  eig3= compute eigenvalues/vectors for 3x3 symmetric matrix
      !
      !  arguments description
      !  ---------------------
      !  inoutput:
      !  --------
      !  u(3,3) : matrix with initial values (only upper half used)
      !
      !  output:
      !  ------
      !  d(3) : eigenvalues associated with columns of v
      !  v(3,3) : matrix of eigenvectors (by column)
      !
      ! ======================================================================
      include 'vaba_param.inc'
      
      dimension :: u(3,3)
      dimension :: d(3)
      dimension v(3,3)
      ! ====================================
      ! local variable
      ! ==============
      integer :: rot, its
      !real :: g, h, aij, sm, thresh, t, c, s, tau

      !real :: a(3), b(3), z(3)
      dimension a(3), b(3), z(3)

      ! loop index
      integer :: i, j, k
      ! ====================================

      ! initialize: do not initialize v
      d(:)= 0.0d0
      v(:,:)= 0.0d0

      ! copy u to v
      v(1:3,1:3)= u(1:3,1:3)

      ! move array into 1d arrays
      a(1) = v(1,2)
      a(2) = v(2,3)
      a(3) = v(1,3)

      do i = 1, 3
         d(i) = v(i,i)
         b(i) = d(i)
         z(i) = 0.0d0

         do j = 1,3
            v(i,j) = 0.0d0
         end do ! j

         v(i,i) = 1.0d0

      end do ! i

      ! check for diagonal case
      sm = abs(a(1)) + abs(a(2)) + abs(a(3))
      g = abs(d(1)) + abs(d(2)) + abs(d(3))

      if (sm < 1.0e-13*g) return

      rot = 0
      do its = 1, 50

         ! set convergence test and threshold
         sm = abs(a(1)) + abs(a(2)) + abs(a(3))
         if ( sm == 0.0d0 ) return

         if( its < 4 ) then
            thresh = 0.011d0*sm
         else
            thresh = 0.0d0
         end if

         ! perform sweeps for rotations
         do i = 1, 3
            j = mod(i,3) + 1
            k = mod(j,3) + 1

            aij = a(i)
            g = 100.0d0 * abs(aij)

          if((abs(d(i))+g/=abs(d(i))).or.(abs(d(j))+g/=abs(d(j)))) then

               if( abs(aij) > thresh ) then

                  a(i) = 0.0d0
                  h = d(j) - d(i)

                  if( abs(h)+g == abs(h) ) then
                     t=aij / h
                  else
                     t=sign(2.0,h/aij)/(abs(h/aij)+sqrt(4.0+(h/aij)**2))
                  end if

                  ! set rotation parameters
                  c = 1.0d0/sqrt(1.0d0+t*t)
                  s = t * c
                  tau = s / (1.0d0+c)

                  ! rotate diagonal terms
                  h = t * aij
                  z(i) = z(i) - h
                  z(j) = z(j) + h
                  d(i) = d(i) - h
                  d(j) = d(j) + h

                  ! rotate off-diagonal terms
                  h = a(j)
                  g = a(k)
                  a(j) = h + s*(g - h*tau)
                  a(k) = g - s*(h + g*tau)

                  ! rotate eigenvectors
                  do k = 1, 3
                     g = v(k,i)
                     h = v(k,j)
                     v(k,i) = g - s*(h + g*tau)
                     v(k,j) = h + s*(g - h*tau)
                  end do ! k

                  rot = rot + 1

               end if

            else

               a(i) = 0.0d0

            end if

         end do ! i

         ! update diagonal terms
         do i = 1, 3
            b(i) = b(i) + z(i)
            d(i) = b(i)
            z(i) = 0.0d0
         end do ! i

      end do ! its

      return
      end subroutine


c ----------------------------------------------------------------------
c Convert between voigt notation and tensor form. it only works for 3D. 
c ----------------------------------------------------------------------
      subroutine voigt_convection(A_voigt, A, convert2voigt, kinetic)
      ! Checked.
        include 'vaba_param.inc'
        integer i,j
        dimension A(3,3), A_voigt(6)
        
        logical convert2voigt, kinetic

        if (convert2voigt) then 
        ! Convect to voigt notation.
            A_voigt(:) = 0
            do i = 1, 3
                A_voigt(i) = A(i,i)        
            end do

            if (kinetic) then
                A_voigt(4) = A(1,2)
                A_voigt(5) = A(2,3)
                A_voigt(6) = A(1,3)
            else  
                A_voigt(4) = 2.0*A(1,2)
                A_voigt(5) = 2.0*A(1,3)
                A_voigt(6) = 2.0*A(2,3)
            end if
        else 
        ! Convect to tensor form.
            A(:,:) = 0
            do i = 1, 3
                A(i,i) = A_voigt(i)        
            end do
            if (kinetic) then
                A(1,2) =  A_voigt(4) 
                A(2,3) =  A_voigt(5) 
                A(1,3) =  A_voigt(6) 
            else 
                A(1,2) = 0.5 * A_voigt(4) 
                A(1,3) = 0.5 * A_voigt(5) 
                A(2,3) = 0.5 * A_voigt(6) 
            end if 
            do i = 1,2
                do j = i+1,3
                    A(j,i) = A(i,j)
                end do
            end do
        end if
        return
      end subroutine



C      
C************************************************************************ 
C      
      subroutine shapefun(AN, dNdxi, XI)
      include 'vaba_param.inc'
      
      dimension AN(8), dNdxi(8,3)
      dimension XI(3)
      parameter(zero=0.d0,one=1.d0,mone=-1.d0,eight=8.d0)

C     Values of shape functions as a function of local coord.
      AN(1) = one/eight*(one-XI(1))*(one-XI(2))*(one-XI(3))
      AN(2) = one/eight*(one+XI(1))*(one-XI(2))*(one-XI(3))
      AN(3) = one/eight*(one+XI(1))*(one+XI(2))*(one-XI(3))
      AN(4) = one/eight*(one-XI(1))*(one+XI(2))*(one-XI(3))
      AN(5) = one/eight*(one-XI(1))*(one-XI(2))*(one+XI(3))
      AN(6) = one/eight*(one+XI(1))*(one-XI(2))*(one+XI(3))
      AN(7) = one/eight*(one+XI(1))*(one+XI(2))*(one+XI(3))
      AN(8) = one/eight*(one-XI(1))*(one+XI(2))*(one+XI(3))
      
C     Derivatives of shape functions respect to local coordinates
      do i=1,8
        do j=1,3
            dNdxi(i,j) =  zero
        end do
      end do
      dNdxi(1,1) =  mone/eight*(one-XI(2))*(one-XI(3))
      dNdxi(1,2) =  mone/eight*(one-XI(1))*(one-XI(3))
      dNdxi(1,3) =  mone/eight*(one-XI(1))*(one-XI(2))
      dNdxi(2,1) =  one/eight*(one-XI(2))*(one-XI(3))
      dNdxi(2,2) =  mone/eight*(one+XI(1))*(one-XI(3))
      dNdxi(2,3) =  mone/eight*(one+XI(1))*(one-XI(2))
      dNdxi(3,1) =  one/eight*(one+XI(2))*(one-XI(3))
      dNdxi(3,2) =  one/eight*(one+XI(1))*(one-XI(3))
      dNdxi(3,3) =  mone/eight*(one+XI(1))*(one+XI(2))
      dNdxi(4,1) =  mone/eight*(one+XI(2))*(one-XI(3))
      dNdxi(4,2) =  one/eight*(one-XI(1))*(one-XI(3))
      dNdxi(4,3) =  mone/eight*(one-XI(1))*(one+XI(2))
      dNdxi(5,1) =  mone/eight*(one-XI(2))*(one+XI(3))
      dNdxi(5,2) =  mone/eight*(one-XI(1))*(one+XI(3))
      dNdxi(5,3) =  one/eight*(one-XI(1))*(one-XI(2))
      dNdxi(6,1) =  one/eight*(one-XI(2))*(one+XI(3))
      dNdxi(6,2) =  mone/eight*(one+XI(1))*(one+XI(3))
      dNdxi(6,3) =  one/eight*(one+XI(1))*(one-XI(2))
      dNdxi(7,1) =  one/eight*(one+XI(2))*(one+XI(3))
      dNdxi(7,2) =  one/eight*(one+XI(1))*(one+XI(3))
      dNdxi(7,3) =  one/eight*(one+XI(1))*(one+XI(2))
      dNdxi(8,1) =  mone/eight*(one+XI(2))*(one+XI(3))
      dNdxi(8,2) =  one/eight*(one-XI(1))*(one+XI(3))
      dNdxi(8,3) =  one/eight*(one-XI(1))*(one+XI(2))
      
      return
      end subroutine          
  
c----------------------------------------------------------------------------------
c     
c     determine the shape function, J matrix(inverse and determinant) and dN/dX matrix
c      
c-------------------------------------------------------------------------------------         
      subroutine shapef3D(ss,nnode,e_coord,shp,dj)

      include 'vaba_param.inc'
        
  	dimension ss(3)
  	dimension e_coord(3,nnode) !the coordinate of element node
  	dimension shp(4,nnode)

c     loacal variables
      dimension bjacob(3,3),binvjacob(3,3)      
      dimension bn(nnode),dndxsi(nnode),dndeta(nnode),dndzta(nnode)

      xsi      = ss(1)
      eta      = ss(2)
      zta      = ss(3)

      call make_N(xsi,eta,zta,nnode,bn,dndxsi,dndeta,dndzta)
        !
        !---- Jacobian:
        !--------------
      bjacob = 0.0
      do i = 1, nnode
       bjacob(1,1) = bjacob(1,1) + dndxsi(i)*e_coord(1,i)   !dX/dxsi
       bjacob(1,2) = bjacob(1,2) + dndeta(i)*e_coord(1,i)   !dX/deta
       bjacob(1,3) = bjacob(1,3) + dndzta(i)*e_coord(1,i)   !dX/zeta
       
       bjacob(2,1) = bjacob(2,1) + dndxsi(i)*e_coord(2,i)   !dY/dxsi
       bjacob(2,2) = bjacob(2,2) + dndeta(i)*e_coord(2,i)   !dY/deta
       bjacob(2,3) = bjacob(2,3) + dndzta(i)*e_coord(2,i)   !dY/zeta

       bjacob(3,1) = bjacob(3,1) + dndxsi(i)*e_coord(3,i)   !dZ/dxsi
       bjacob(3,2) = bjacob(3,2) + dndeta(i)*e_coord(3,i)   !dZ/deta
       bjacob(3,3) = bjacob(3,3) + dndzta(i)*e_coord(3,i)   !dZ/zeta
      end do
  
      binvjacob = bjacob

      call inverseJ(binvjacob,dj) !dj is the determinant of the Jacobi matrix
      shp = 0.0
      do i = 1, nnode
        shp(1,i)=dndxsi(i)*binvjacob(1,1)+dndeta(i)*binvjacob(2,1)
     *           +dndzta(i)*binvjacob(3,1) !dNi/dX
        shp(2,i)=dndxsi(i)*binvjacob(1,2)+dndeta(i)*binvjacob(2,2)
     *           +dndzta(i)*binvjacob(3,2) !dNi/dY
        shp(3,i)=dndxsi(i)*binvjacob(1,3)+dndeta(i)*binvjacob(2,3)
     *           +dndzta(i)*binvjacob(3,3) !dNi/dZ
      end do

      do i = 1, nnode
          shp(4,i) = bn(i)
      end do

      return
      end subroutine
c-----------------------------------------------------------------------------
c      
!     determine the shape function  
c      
c------------------------------------------------------------------------------
      subroutine make_N(xsi,eta,zta,nnode,bn,dndxsi,dndeta,dndzta)

      include 'vaba_param.inc'
      
      dimension bn(nnode),dndxsi(nnode),dndeta(nnode),dndzta(nnode)

!c      8/------/|7      
!c     5|------|6|
!c      | 4    | /3
!c	  |      |/
!C     1-------/2
!c            
!c      z/|\ /y
!c        | /
!c        |/
!c         ------>x
      ! loacal variables
      dimension bNodeNat(3,8)
      bNodeNat(1:3,1) = (/-1,-1,-1/)
      bNodeNat(1:3,2) = (/ 1,-1,-1/)
      bNodeNat(1:3,3) = (/ 1, 1,-1/)
      bNodeNat(1:3,4) = (/-1, 1,-1/)
      bNodeNat(1:3,5) = (/-1,-1, 1/)
      bNodeNat(1:3,6) = (/ 1,-1, 1/)
      bNodeNat(1:3,7) = (/ 1, 1, 1/)
      bNodeNat(1:3,8) = (/-1, 1, 1/)
      do i = 1, nnode
        bn(i) = 1/8.0d0*(1+bNodeNat(1,i)*xsi)*(1+bNodeNat(2,i)*eta)
     *                                     *(1+bNodeNat(3,i)*zta)

        dndxsi(i)=1/8.0d0*bNodeNat(1,i)*(1+bNodeNat(2,i)*eta)
     *                                     *(1+bNodeNat(3,i)*zta)

        dndeta(i)=1/8.0d0*bNodeNat(2,i)*(1+bNodeNat(1,i)*xsi)
     *                                     *(1+bNodeNat(3,i)*zta)
    
        dndzta(i)=1/8.0d0*bNodeNat(3,i)*(1+bNodeNat(2,i)*eta)
     *                                     *(1+bNodeNat(1,i)*xsi)
      end do
      
      return
      end subroutine
c----------------------------------------------------------------------
c
!     return the inverse of a 3 by 3 Matrix bM and its determinant
c      
c----------------------------------------------------------------------     
      subroutine inverseJ(bM,DJ)
  
      include 'vaba_param.inc'
  
      dimension bM(3,3)
          
      !local variables
      dimension AdjM(3,3)
    
      AdjM(1:3,1:3) = 0.0
      
      AdjM(1,1) =  bM(2,2)*bM(3,3)-bM(2,3)*bM(3,2)
      AdjM(2,1) = -bM(2,1)*bM(3,3)+bM(2,3)*bM(3,1)
      AdjM(3,1) =  bM(2,1)*bM(3,2)-bM(2,2)*bM(3,1)
    
      AdjM(1,2) = -bM(1,2)*bM(3,3)+bM(1,3)*bM(3,2)
      AdjM(2,2) =  bM(1,1)*bM(3,3)-bM(1,3)*bM(3,1)
      AdjM(3,2) = -bM(1,1)*bM(3,2)+bM(1,2)*bM(3,1)
    
      AdjM(1,3) =  bM(1,2)*bM(2,3)-bM(2,2)*bM(1,3)
      AdjM(2,3) = -bM(1,1)*bM(2,3)+bM(2,1)*bM(1,3)
      AdjM(3,3) =  bM(1,1)*bM(2,2)-bM(2,1)*bM(1,2)
    
      
      DJ = bM(1,1)*AdjM(1,1)+bM(1,2)*AdjM(2,1)+bM(1,3)*AdjM(3,1)

      if(DJ==0) then
        write(*,*) "warning: zero Jacob value!"
      else
        bM(:,:) = AdjM(:,:)/DJ
      end if

      return
      end subroutine
c--------------------------------------------------------------------  
c      
!     determine the gauss point locations and weights     
c      
c---------------------------------------------------------------------- 
      subroutine make_quadrature(xg)

      use vars_module

      include 'vaba_param.inc'
      
      dimension xg(4,NPT)
      
      !local variables
      dimension r1pt(2), r1wt(2)
 
      xg(:,:) = 0.d0
      if(NPT==1)then
        xg(4,:) = 8.d0
        return
      else
        r1pt(1) = -0.577350269189625764509148780502
        r1pt(2) =  0.577350269189625764509148780502
        r1wt(1) = 1.000000000000000 
        r1wt(2) = 1.000000000000000

        kk = 1
        do i = 1, 2
      	  do j = 1, 2
      	    do k = 1, 2
             xg(1,kk)=r1pt(k)
             xg(2,kk)=r1pt(j)
             xg(3,kk)=r1pt(i)
             xg(4,kk)=r1wt(i)*r1wt(j)*r1wt(k)
             kk=kk+1
            end do
          enddo
        end do
      end if
  
      return 
      end subroutine
c----------------------------------------------------------------------
c      
!     determine the mass matrix of phase field   
c     
c-----------------------------------------------------------------------       
      subroutine cal_cmass(shp,we,cmass,props,nprops)      
      
      include 'vaba_param.inc'
       
      dimension shp(4,8),cmass(8,8)
      dimension props(nprops)
      
      eta=props(1)
      
      do i = 1, 8
          do j = 1, 8
              cmass(i,j) = cmass(i,j)+we*eta*shp(4,i)*shp(4,j)
          end do
      end do
      
      return 
      end subroutine

      
      
      SUBROUTINE VUAMP(
     *     ampName, time, ampValueOld, dt, nprops, props, nSvars, 
     *     svars, lFlagsInfo, nSensor, sensorValues, sensorNames,	
     *     jSensorLookUpTable,
     *     AmpValueNew,
     *     lFlagsDefine,
     *     AmpDerivative, AmpSecDerivative, AmpIncIntegral)

      INCLUDE 'VABA_PARAM.INC'

C     time indices
      parameter (iStepTime        = 1,
     *           iTotalTime       = 2,
     *           nTime            = 2)
C     flags passed in for information
      parameter (iInitialization   = 1,
     *           iRegularInc       = 2,
     *           ikStep            = 3,
     *           nFlagsInfo        = 3)
C     optional flags to be defined
      parameter (iComputeDeriv     = 1,
     *           iComputeSecDeriv  = 2,
     *           iComputeInteg     = 3,
     *           iStopAnalysis     = 4,
     *           iConcludeStep     = 5,
     *           nFlagsDefine      = 5)
      dimension time(nTime), lFlagsInfo(nFlagsInfo),
     *          lFlagsDefine(nFlagsDefine),
     *          sensorValues(nSensor),
     *          props(nprops),
     *          sVars(nSvars)

      character*80 sensorNames(nSensor)
      character*80 ampName
      dimension jSensorLookUpTable(*)

      erate = props(1)
      AmpValueNew = erate*exp(erate*time(iStepTime))

      RETURN
      END
      

C
C User subroutine VDLOAD
      subroutine vdload (
C Read only (unmodifiable) variables -
     *     nblock, ndim, stepTime, totalTime, 
     *     amplitude, curCoords, velocity, dircos, 
     *     jltyp, sname,
C Write only (modifiable) variable -
     *     value )
C
      include 'vaba_param.inc'
      parameter ( 
     *     d0 = 120.00d0,
     *     pa = 0.0d0,
     *     ps = 364.0d0,
     *     theta = 0.124651928275d-3,
c     *     theta = 0.1569262d-3,
     *     timea = 0.01875d-3,
     *     timed = 0.195d-3,
     *     zero = 0.0d0,
     *     one = 1.0d0
     *     )
C
      dimension curCoords(nblock,ndim), 
     *     velocity(nblock,ndim),
     *     dircos(nblock,ndim,ndim), 
     *     value(nblock)
      character*80 sname
c    
c     local variables
c 
      real xx,yy,zz,dd,dd0,xc,yc,zc,tt,Tint,press
c
      xc = zero
      yc = -5.0
      zc = zero
      time = totalTime
      dd0 = d0 * d0
      Tint = timed-timea
      tt = time - timea
      do k = 1, nblock
         press = zero
         if ((time .gt. timea) .and. (time .le. timed)) then
                xx = curCoords(k,1) - xc
                yy = curCoords(k,2) - yc
                zz = curCoords(k,3) - zc
                xx = xx * xx
                yy = yy * yy
                zz = zz * zz
                dd = xx + yy + zz
                press = (ps-pa) * (one - (tt/Tint))
                press = press * (exp(-tt/theta))
                press = press * (exp(-dd/dd0))
         end if
         value(k) = amplitude * press
c         write(*,*) time, value(k)
      end do
*     
      return
      end
c
c===================================================


      















