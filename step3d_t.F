#include "cppdefs.h"
      MODULE step3d_t_mod
#if !defined TS_FIXED && (defined NONLINEAR && defined SOLVE3D)
!
!svn $Id: step3d_t.F 719 2014-03-13 22:25:13Z arango $
!=======================================================================
!  Copyright (c) 2002-2014 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This routine time-steps tracer equations.  Notice that advective    !
!  and diffusive terms are time-stepped differently. It applies the    !
!  corrector time-step for horizontal/vertical advection,  vertical    !
!  diffusion, nudging if necessary, and lateral boundary conditions.   !
!                                                                      !
!  Notice that at input the tracer arrays have:                        !
!                                                                      !
!    t(:,:,:,nnew,:)   m Tunits  n+1     horizontal/vertical diffusion !
!                                        terms plus source/sink terms  !
!                                        (biology, sediment), if any   !
!                                                                      !
!    t(:,:,:,3   ,:)   Tunits    n+1/2   advective terms and vertical  !
!                                        diffusion predictor step      !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: step3d_t
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE step3d_t (ng, tile)
!***********************************************************************
!
      USE mod_param
# ifdef DIAGNOSTICS_TS
      USE mod_diags
# endif
      USE mod_grid
      USE mod_mixing
      USE mod_ocean
      USE mod_stepping
      USE netcdf  
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 35)
# endif
      CALL step3d_t_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    nrhs(ng), nstp(ng), nnew(ng),                 &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
     &                    GRID(ng) % umask,                             &
     &                    GRID(ng) % vmask,                             &
# endif
# ifdef TS_MPDATA
#  ifdef WET_DRY
     &                    GRID(ng) % rmask_wet,                         &
     &                    GRID(ng) % umask_wet,                         &
     &                    GRID(ng) % vmask_wet,                         &
#  endif
     &                    GRID(ng) % omn,                               &
     &                    GRID(ng) % om_u,                              &
     &                    GRID(ng) % om_v,                              &
     &                    GRID(ng) % on_u,                              &
     &                    GRID(ng) % on_v,                              &
# endif
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    GRID(ng) % Hz,                                &
     &                    GRID(ng) % Huon,                              &
     &                    GRID(ng) % Hvom,                              &
     &                    GRID(ng) % z_r,                               &
     &                    MIXING(ng) % Akt,                             &
     &                    OCEAN(ng) % W,                                &
# if defined FLOATS && defined FLOAT_VWALK
     &                    MIXING(ng) % dAktdz,                          &
# endif
# ifdef DIAGNOSTICS_TS
     &                    DIAGS(ng) % DiaTwrk,                          &
# endif
     &                    OCEAN(ng) % t,                                &
     &                    OCEAN(ng) % u,                                &
     &                    OCEAN(ng) % v,                                &
     &                    OCEAN(ng) % pden,                             &
     &                    GRID(ng)  % lonr,                             &
     &                    GRID(ng)  % latr)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 35)
# endif
      RETURN
      END SUBROUTINE step3d_t
!
!***********************************************************************
      SUBROUTINE step3d_t_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          nrhs, nstp, nnew,                       &
# ifdef MASKING
     &                          rmask, umask, vmask,                    &
# endif
# ifdef TS_MPDATA
#  ifdef WET_DRY
     &                          rmask_wet, umask_wet, vmask_wet,        &
#  endif
     &                          omn, om_u, om_v, on_u, on_v,            &
# endif
     &                          pm, pn,                                 &
     &                          Hz, Huon, Hvom,                         &
     &                          z_r,                                    &
     &                          Akt,                                    &
     &                          W,                                      &
# if defined FLOATS && defined FLOAT_VWALK
     &                          dAktdz,                                 &
# endif
# ifdef DIAGNOSTICS_TS
     &                          DiaTwrk,                                &
# endif
     &                          t, u, v, pden, lonr, latr)
!***********************************************************************
!
      USE mod_param
      USE mod_clima
      USE mod_ncparam
      USE netcdf
# if defined NESTING && !defined ONE_WAY
      USE mod_nesting
# endif
      USE mod_scalars
      USE mod_sources
!
      USE exchange_3d_mod, ONLY : exchange_r3d_tile
# ifdef DISTRIBUTE
#  if defined FLOATS && defined FLOAT_VWALK
      USE mp_exchange_mod, ONLY : mp_exchange3d
#  endif
      USE mp_exchange_mod, ONLY : mp_exchange4d
# endif
# ifdef TS_MPDATA
      USE mpdata_adiff_mod
# endif
# ifdef NESTING
      USE nesting_mod, ONLY : bry_fluxes
# endif
      USE t3dbc_mod, ONLY : t3dbc_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs, nstp, nnew
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef TS_MPDATA
#   ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:,LBj:)
      real(r8), intent(in) :: umask_wet(LBi:,LBj:)
      real(r8), intent(in) :: vmask_wet(LBi:,LBj:)
#   endif
      real(r8), intent(in) :: omn(LBi:,LBj:)
      real(r8), intent(in) :: om_u(LBi:,LBj:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: on_v(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: Huon(LBi:,LBj:,:)
      real(r8), intent(in) :: Hvom(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      !!! zhousy
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
      real(r8), intent(in) :: pden(LBi:,LBj:,:)
      real(r8), intent(in) :: lonr(LBi:,LBj:)
      real(r8), intent(in) :: latr(LBi:,LBj:)
#  ifdef SUN
      real(r8), intent(in) :: Akt(LBi:UBi,LBj:UBj,0:N(ng),NAT)
#  else
      real(r8), intent(in) :: Akt(LBi:,LBj:,0:,:)
#  endif
      real(r8), intent(in) :: W(LBi:,LBj:,0:)
#  ifdef DIAGNOSTICS_TS
      real(r8), intent(inout) :: DiaTwrk(LBi:,LBj:,:,:,:)
#  endif
#  ifdef SUN
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
#  else
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
#  endif
#  if defined FLOATS && defined FLOAT_VWALK
      real(r8), intent(out) :: dAktdz(LBi:,LBj:,:)
#  endif

# else

#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef TS_MPDATA
#   ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask_wet(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(in) :: omn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_v(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: Huon(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: Hvom(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: Akt(LBi:UBi,LBj:UBj,0:N(ng),NAT)
      real(r8), intent(in) :: W(LBi:UBi,LBj:UBj,0:N(ng))

#  ifdef DIAGNOSTICS_TS
      real(r8), intent(inout) :: DiaTwrk(LBi:UBi,LBj:UBj,N(ng),NT(ng),  &
     &                                   NDT)
#  endif
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))

#  if defined FLOATS && defined FLOAT_VWALK
      real(r8), intent(out) :: dAktdz(LBi:UBi,LBj:UBj,N(ng))
#  endif
# endif
!
!  Local variable declarations.
!
# ifdef NESTING
      integer :: ILB, IUB, JLB, JUB
      integer :: dg, cr, rg
# endif
      integer :: i, ic, is, itrc, j, k, ltrc
# ifdef DIAGNOSTICS_TS
      integer :: idiag
# endif
      real(r8), parameter :: eps = 1.0E-16_r8

      real(r8) :: cff, cff1, cff2, cff3

      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: CF
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: BC
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: DC
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FC

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FX
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: curv
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: grad

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: oHz
# ifdef TS_MPDATA
#  ifdef DIAGNOSTICS_TS
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,3) :: Dhadv
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng),NT(ng)) :: Dvadv
#  endif
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng),NT(ng)) :: Ta
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: Ua
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: Va
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: Wa
      !!! zhousy
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: dudx
      real(r8), dimension(IminS:ImaxS+1,JminS:JmaxS+1,0:N(ng)) :: dude_p
      real(r8), dimension(IminS:ImaxS+1,JminS:JmaxS+1,0:N(ng)) :: dvdx_p
      real(r8), dimension(IminS:ImaxS+1,JminS:JmaxS+1,0:N(ng)) :: dude_r
      real(r8), dimension(IminS:ImaxS+1,JminS:JmaxS+1,0:N(ng)) :: dvdx_r
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: dvde
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: dtdx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: dtde
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: dbdx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: dbde
      !real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: wtemp
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: k_rho
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: t_mld
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: k_mld
  
      INTEGER :: fidA, status, w1_hs_id, w2_hs_id, w3_hs_id, w4_hs_id
      INTEGER :: b1_hs_id, b2_hs_id, b3_hs_id, b4_hs_id, input_id
      INTEGER :: fidB, w1_hl_id, w2_hl_id, w3_hl_id, w4_hl_id
      INTEGER :: b1_hl_id, b2_hl_id, b3_hl_id, b4_hl_id

      CHARACTER(LEN=100) :: file_in_hs
      CHARACTER(LEN=100) :: file_in_hl

      REAL*8 :: b1_hs(1,64), b2_hs(1,32), b3_hs(1,16), b4_hs(1,1)
      REAL*8 :: input(1,10), w1_hs(10,64), w2_hs(64,32), w3_hs(32,16) 
      REAL*8 :: w4_hs(16,1),wtemp_01(1,1),wtemp_02(1,1),tao(1,1)
      REAL*8 :: wb1_hs(1,64),wb2_hs(1,32),wb3_hs(1,16),wb4_hs(1,1) 

      REAL*8 :: b1_hl(1,32), b2_hl(1,16), b3_hl(1,8), b4_hl(1,1)
      REAL*8 :: w1_hl(5,32), w2_hl(32,16), w3_hl(16,8), w4_hl(8,1)
      REAL*8 :: wb1_hl(1,32),wb2_hl(1,16) ,wb3_hl(1,8) ,wb4_hl(1,1)
# endif

# include "set_bounds.h"

# ifdef NESTING
!
!  Notice that the trace flux boundary arrays are dimensioned with the
!  global dimensions of grid to facilitate processing.
!
      ILB=BOUNDS(ng)%LBi(-1)
      IUB=BOUNDS(ng)%UBi(-1)
      JLB=BOUNDS(ng)%LBj(-1)
      JUB=BOUNDS(ng)%UBj(-1)
# endif

!
!-----------------------------------------------------------------------
!  Compute wtemp. zhousy
!-----------------------------------------------------------------------
!
!
      file_in_hs='/home/zsy/project/kpp/nn/20240102_test01.nc'
      file_in_hl='/home/zsy/project/kpp/nn/hl_5_0624.nc' 
      ! Read hs dimension IDs
      !print *, 'ttt',TRIM(file_in_hs)
      status = NF90_OPEN(TRIM(file_in_hs),0,fidA)
      if (status /= nf90_noerr) print *, '1t',nf90_strerror(status)
      ! Read hs dimension IDs
      status = NF90_INQ_VARID(fidA,"w1",w1_hs_id)
      if (status /= nf90_noerr) print *, '2t',nf90_strerror(status)
      status = NF90_INQ_VARID(fidA,"b1",b1_hs_id)
      if (status /= nf90_noerr) print *, '3t',nf90_strerror(status)
      status = NF90_INQ_VARID(fidA,"w2",w2_hs_id)
      !if (status /= nf90_noerr) print *, '4',nf90_strerror(status)
      status = NF90_INQ_VARID(fidA,"b2",b2_hs_id)
      !if (status /= nf90_noerr) print *, '5',nf90_strerror(status)
      status = NF90_INQ_VARID(fidA,"w3",w3_hs_id)
      !if (status /= nf90_noerr) print *, '6',nf90_strerror(status)
      status = NF90_INQ_VARID(fidA,"b3",b3_hs_id)
      !if (status /= nf90_noerr) print *, '7',nf90_strerror(status)
      status = NF90_INQ_VARID(fidA,"w4",w4_hs_id)
      !if (status /= nf90_noerr) print *, '8',nf90_strerror(status)
      status = NF90_INQ_VARID(fidA,"b4",b4_hs_id)
      !if (status /= nf90_noerr) print *, '9',nf90_strerror(status)

      ! Read hs variable values
      status = NF90_GET_VAR(fidA,w1_hs_id,w1_hs)
      !if (status /= nf90_noerr) print *, nf90_strerror(status)
      status = NF90_GET_VAR(fidA,b1_hs_id,b1_hs)
      status = NF90_GET_VAR(fidA,w2_hs_id,w2_hs)
      status = NF90_GET_VAR(fidA,b2_hs_id,b2_hs)
      status = NF90_GET_VAR(fidA,w3_hs_id,w3_hs)
      status = NF90_GET_VAR(fidA,b3_hs_id,b3_hs)
      status = NF90_GET_VAR(fidA,w4_hs_id,w4_hs)
      status = NF90_GET_VAR(fidA,b4_hs_id,b4_hs)

      status = NF90_CLOSE(fidA)
  
      ! Read hl dimension IDs
      status = NF90_OPEN(TRIM(file_in_hl),0,fidB)

      ! Read hs dimension IDs
      status = NF90_INQ_VARID(fidB,"w1",w1_hl_id)
      status = NF90_INQ_VARID(fidB,"b1",b1_hl_id)
      status = NF90_INQ_VARID(fidB,"w2",w2_hl_id)
      status = NF90_INQ_VARID(fidB,"b2",b2_hl_id)
      status = NF90_INQ_VARID(fidB,"w3",w3_hl_id)
      status = NF90_INQ_VARID(fidB,"b3",b3_hl_id)
      status = NF90_INQ_VARID(fidB,"w4",w4_hl_id)
      status = NF90_INQ_VARID(fidB,"b4",b4_hl_id)

 
      ! Read hs variable values
      status = NF90_GET_VAR(fidB,w1_hl_id,w1_hl)
      status = NF90_GET_VAR(fidB,b1_hl_id,b1_hl)
      status = NF90_GET_VAR(fidB,w2_hl_id,w2_hl)
      status = NF90_GET_VAR(fidB,b2_hl_id,b2_hl)
      status = NF90_GET_VAR(fidB,w3_hl_id,w3_hl)
      status = NF90_GET_VAR(fidB,b3_hl_id,b3_hl)
      status = NF90_GET_VAR(fidB,w4_hl_id,w4_hl)
      status = NF90_GET_VAR(fidB,b4_hl_id,b4_hl)

      status = NF90_CLOSE(fidB)  

      DO j=Jstr,Jend        !!! RHO-points
        DO i=Istr,Iend
          DO k=1,N(ng)
            dudx(i,j,k)= pn(i, j)*                                      &
     &        (u(i+1,j,k,nstp)-u(i,j,k,nstp))
            dvde(i,j,k)= pm(i, j)*                                      &
     &        (v(i,j+1,k,nstp)-v(i,j,k,nstp))                        
          END DO
        END DO
      END DO
  
      DO j=Jstr,Jend+1      !!! PSI-points
        DO i=Istr,Iend+1
          DO k=1,N(ng)
            dude_p(i,j,k)= pn(i,j)*                                     &
     &        (u(i,j+1,k,nstp)-u(i,j,k,nstp))                        
            dvdx_p(i,j,k)= pm(i,j)*                                     &
     &        (v(i+1,j,k,nstp)-v(i,j,k,nstp))
          END DO
        END DO
      END DO
  
      DO j=Jstr,Jend        !!! PSI ==> RHO 
        DO i=Istr,Iend
          DO k=1,N(ng)
            dude_r(i,j,k)=0.25*                                         &
     &          (dude_p(i,j,k)+dude_p(i+1,j,k)+                         &
     &           dude_p(i,j+1,k)+dude_p(i+1,j+1,k))             
            dvdx_r(i,j,k)=0.25*                                         &
     &          (dvdx_p(i,j,k)+dvdx_p(i+1,j,k)+                         &
     &           dvdx_p(i,j+1,k)+dvdx_p(i+1,j+1,k))            
          END DO
        END DO
      END DO
  
      DO j=Jstr,Jend        !!! RHO-points
        DO i=Istr,Iend
          DO k=1,N(ng)
            dtdx(i,j,k)= 0.5*pm(i,j)*                                   &
     &        (t(i+1,j,k,nrhs,1)-t(i-1,j,k,nrhs,1))                
            dtde(i,j,k)= 0.5*pn(i,j)*                                   &
     &        (t(i,j+1,k,nrhs,1)-t(i,j-1,k,nrhs,1))
            !print *, "Hz: ", Hz(i,j,k)
            !print *, "dt: ", dt(ng)                
          END DO
        END DO
      END DO
  
      DO j=Jstr,Jend        !!! RHO-points
        DO i=Istr,Iend
          DO k=1,N(ng)
            
            dbdx(i,j,k)= -9.8/1.025*0.5*pm(i,j)*                        &
     &        (pden(i+1,j,k)-pden(i-1,j,k))                
            dbde(i,j,k)= -9.8/1.025*0.5*pn(i,j)*                        &
     &        (pden(i,j+1,k)-pden(i,j-1,k))
            !print *, "dbdx: ", dbdx(i,j,k)
            !print *, "dbde: ", dbde(i,j,k)                
          END DO
        END DO
      END DO
 
  
      DO j=Jstr,Jend        !!! RHO-points
        DO i=Istr,Iend
          !print *, z_r(i,j,N(ng))
          !k_rho(i,j)=50
          IF (z_r(i,j,1) < -500) THEN
            DO k=1,N(ng)
              IF (z_r(i,j,k) > -10) THEN
                k_rho(i,j)=k
                EXIT
              END IF
            END DO
          END IF

        END DO
      END DO
      

 
      DO j=Jstr,Jend        !!! RHO-points
        DO i=Istr,Iend
          
          IF (z_r(i,j,1) < -500) THEN
          k = k_rho(i,j)
          t_mld(i,j)=t(i,j,k,nstp,1)-(t(i,j,k,nstp,1)-t(i,j,k-1,nstp,1))&
     &       /(z_r(i,j,k)-z_r(i,j,k-1))*(z_r(i,j,k)+10)-0.2
          END IF

        END DO
      END DO

      DO j=Jstr,Jend        !!! RHO-points
        DO i=Istr,Iend

          IF (z_r(i,j,1) < -500) THEN
            !print *, lonr(i,j),"--",latr(i,j)
            DO k=N(ng)-1,1,-1
              k_mld(i,j)=999;
              IF (t(i,j,k,nstp,1) < t_mld(i,j)) THEN
                k_mld(i,j)=k
                EXIT
              END IF
            END DO
          END IF

        END DO
      END DO

      DO j=Jstr,Jend
        DO i=Istr,Iend
          IF (z_r(i,j,1) < -500 .and. k_mld(i,j)<100) THEN
            cff = z_r(i,j,k_mld(i,j))*1.2
            DO k=N(ng)-1,1,-1
               IF (z_r(i,j,k)<cff) THEN
                 k_mld(i,j)=k
                 EXIT
               END IF
            END DO
          END IF
        END DO
      END DO
      !print *, "5,5: ", k_rho(5,5)
      !print *, "2,2: ", k_rho(2,2)
      !DO j=Jstr,Jend        !!! RHO-points
      !  DO i=Istr,Iend
      !    print *, k_mld(i,j), k_rho(i,j),z_r(i,j,k_mld(i,j)) 
      !  END DO
      !END DO

!
!-----------------------------------------------------------------------
!  Time-step horizontal advection term.
!-----------------------------------------------------------------------
!
!  Compute inverse thickness.
!
# ifdef TS_MPDATA
#  define I_RANGE Istrm2,Iendp2
#  define J_RANGE Jstrm2,Jendp2
# else
#  define I_RANGE Istr,Iend
#  define J_RANGE Jstr,Jend
# endif

      DO k=1,N(ng)
        DO j=J_RANGE
          DO i=I_RANGE
            oHz(i,j,k)=1.0_r8/Hz(i,j,k)
          END DO
        END DO
      END DO
# undef I_RANGE
# undef J_RANGE
# ifdef TS_MPDATA
!
!  The MPDATA algorithm requires a three-point footprint, so exchange
!  boundary data on t(:,:,:,nnew,:) so other processes computed earlier
!  (horizontal diffusion, biology, or sediment) are accounted.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        DO itrc=1,NT(ng)
          CALL exchange_r3d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            t(:,:,:,nnew,itrc))
        END DO
      END IF

#  ifdef DISTRIBUTE
      CALL mp_exchange4d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),      &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    t(:,:,:,nnew,:))
#  endif
# endif
!
!  Compute horizontal tracer advection fluxes.
!
      T_LOOP : DO itrc=1,NT(ng)
        K_LOOP : DO k=1,N(ng)

# if defined TS_C2HADVECTION
!
!  Second-order, centered differences horizontal advective fluxes.
!
          DO j=Jstr,Jend
            DO i=Istr,Iend+1
              FX(i,j)=Huon(i,j,k)*                                      &
     &                0.5_r8*(t(i-1,j,k,3,itrc)+                        &
     &                        t(i  ,j,k,3,itrc))
            END DO
          END DO
          DO j=Jstr,Jend+1
            DO i=Istr,Iend
              FE(i,j)=Hvom(i,j,k)*                                      &
     &                0.5_r8*(t(i,j-1,k,3,itrc)+                        &
     &                        t(i,j  ,k,3,itrc))
            END DO
          END DO

# elif defined TS_MPDATA
!
!  First-order, upstream differences horizontal advective fluxes.
!
          DO j=JstrVm2,Jendp2i
            DO i=IstrUm2,Iendp3
              cff1=MAX(Huon(i,j,k),0.0_r8)
              cff2=MIN(Huon(i,j,k),0.0_r8)
              FX(i,j)=cff1*t(i-1,j,k,3,itrc)+                           &
     &                cff2*t(i  ,j,k,3,itrc)
            END DO
          END DO
          DO j=JstrVm2,Jendp3
            DO i=IstrUm2,Iendp2i
              cff1=MAX(Hvom(i,j,k),0.0_r8)
              cff2=MIN(Hvom(i,j,k),0.0_r8)
              FE(i,j)=cff1*t(i,j-1,k,3,itrc)+                           &
     &                cff2*t(i,j  ,k,3,itrc)
            END DO
          END DO

# else
!
#  if defined TS_U3HADVECTION
!  Third-order, uptream-biased horizontal advective fluxes.
#  elif defined TS_A4HADVECTION
!  Fourth-order, Akima horizontal advective fluxes.
#  else
!  Fourth-order, centered differences horizontal advective fluxes.
#  endif
!
          DO j=Jstr,Jend
            DO i=Istrm1,Iendp2
              FX(i,j)=t(i  ,j,k,3,itrc)-                                &
     &                t(i-1,j,k,3,itrc)
#  ifdef MASKING
              FX(i,j)=FX(i,j)*umask(i,j)
#  endif
            END DO
          END DO
          IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
            IF (DOMAIN(ng)%Western_Edge(tile)) THEN
              DO j=Jstr,Jend
                FX(Istr-1,j)=FX(Istr,j)
              END DO
            END IF
          END IF
          IF (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) THEN
            IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
              DO j=Jstr,Jend
                FX(Iend+2,j)=FX(Iend+1,j)
              END DO
            END IF
          END IF
!
          DO j=Jstr,Jend
            DO i=Istr-1,Iend+1
#  if defined TS_U3HADVECTION
              curv(i,j)=FX(i+1,j)-FX(i,j)
#  elif defined TS_A4HADVECTION
              cff=2.0_r8*FX(i+1,j)*FX(i,j)
              IF (cff.gt.eps) THEN
                grad(i,j)=cff/(FX(i+1,j)+FX(i,j))
              ELSE
                grad(i,j)=0.0_r8
              END IF
#  else
              grad(i,j)=0.5_r8*(FX(i+1,j)+FX(i,j))
#  endif
            END DO
          END DO
!
          cff1=1.0_r8/6.0_r8
          cff2=1.0_r8/3.0_r8
          DO j=Jstr,Jend
            DO i=Istr,Iend+1
#  ifdef TS_U3HADVECTION
              FX(i,j)=Huon(i,j,k)*0.5_r8*                               &
     &                (t(i-1,j,k,3,itrc)+                               &
     &                 t(i  ,j,k,3,itrc))-                              &
     &                cff1*(curv(i-1,j)*MAX(Huon(i,j,k),0.0_r8)+        &
     &                      curv(i  ,j)*MIN(Huon(i,j,k),0.0_r8))
#  else
              FX(i,j)=Huon(i,j,k)*0.5_r8*                               &
     &                (t(i-1,j,k,3,itrc)+                               &
     &                 t(i  ,j,k,3,itrc)-                               &
     &                 cff2*(grad(i  ,j)-                               &
     &                       grad(i-1,j)))
#  endif
            END DO
          END DO
!
          DO j=Jstrm1,Jendp2
            DO i=Istr,Iend
              FE(i,j)=t(i,j  ,k,3,itrc)-                                &
     &                t(i,j-1,k,3,itrc)
#  ifdef MASKING
              FE(i,j)=FE(i,j)*vmask(i,j)
#  endif
            END DO
          END DO
          IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) THEN
            IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
              DO i=Istr,Iend
                FE(i,Jstr-1)=FE(i,Jstr)
              END DO
            END IF
          END IF
          IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
            IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
              DO i=Istr,Iend
                FE(i,Jend+2)=FE(i,Jend+1)
              END DO
            END IF
          END IF
!
          DO j=Jstr-1,Jend+1
            DO i=Istr,Iend
#  if defined TS_U3HADVECTION
              curv(i,j)=FE(i,j+1)-FE(i,j)
#  elif defined TS_A4HADVECTION
              cff=2.0_r8*FE(i,j+1)*FE(i,j)
              IF (cff.gt.eps) THEN
                grad(i,j)=cff/(FE(i,j+1)+FE(i,j))
              ELSE
                grad(i,j)=0.0_r8
              END IF
#  else
              grad(i,j)=0.5_r8*(FE(i,j+1)+FE(i,j))
#  endif
            END DO
          END DO
!
          cff1=1.0_r8/6.0_r8
          cff2=1.0_r8/3.0_r8
          DO j=Jstr,Jend+1
            DO i=Istr,Iend
#  ifdef TS_U3HADVECTION
              FE(i,j)=Hvom(i,j,k)*0.5_r8*                               &
     &                (t(i,j-1,k,3,itrc)+                               &
     &                 t(i,j  ,k,3,itrc))-                              &
     &                cff1*(curv(i,j-1)*MAX(Hvom(i,j,k),0.0_r8)+        &
     &                      curv(i,j  )*MIN(Hvom(i,j,k),0.0_r8))
#  else
              FE(i,j)=Hvom(i,j,k)*0.5_r8*                               &
     &                (t(i,j-1,k,3,itrc)+                               &
     &                 t(i,j  ,k,3,itrc)-                               &
     &                 cff2*(grad(i,j  )-                               &
     &                       grad(i,j-1)))
#  endif
            END DO
          END DO
# endif
!
!  Apply tracers point sources to the horizontal advection terms,
!  if any.
!
          IF (.not.LwSrc(ng).and.ANY(LtracerSrc(:,ng))) THEN
            DO is=1,Nsrc(ng)
              i=SOURCES(ng)%Isrc(is)
              j=SOURCES(ng)%Jsrc(is)
              IF (INT(SOURCES(ng)%Dsrc(is)).eq.0) THEN
# ifdef TS_MPDATA
                IF (((IstrUm2.le.i).and.(i.le.Iendp3)).and.             &
     &              ((JstrVm2.le.j).and.(j.le.Jendp2i))) THEN
# else
                IF (((Istr.le.i).and.(i.le.Iend+1)).and.                &
     &              ((Jstr.le.j).and.(j.le.Jend))) THEN
# endif
                  IF (LtracerSrc(itrc,ng)) THEN
                    FX(i,j)=Huon(i,j,k)*                                &
     &                      SOURCES(ng)%Tsrc(is,k,itrc)
# ifdef MASKING
                  ELSE
                    IF ((rmask(i  ,j).eq.0.0_r8).and.                   &
     &                  (rmask(i-1,j).eq.1.0_r8)) THEN
                      FX(i,j)=Huon(i,j,k)*t(i-1,j,k,3,itrc)
                    ELSE IF ((rmask(i  ,j).eq.1.0_r8).and.              &
     &                       (rmask(i-1,j).eq.0.0_r8)) THEN
                      FX(i,j)=Huon(i,j,k)*t(i  ,j,k,3,itrc)
                    END IF
# endif
                  END IF
                END IF
              ELSE IF (INT(SOURCES(ng)%Dsrc(is)).eq.1) THEN
# ifdef TS_MPDATA
                IF (((IstrUm2.le.i).and.(i.le.Iendp2i)).and.            &
     &              ((JstrVm2.le.j).and.(j.le.Jendp3))) THEN
# else
                IF (((Istr.le.i).and.(i.le.Iend)).and.                  &
     &              ((Jstr.le.j).and.(j.le.Jend+1))) THEN
# endif
                  IF (LtracerSrc(itrc,ng)) THEN
                    FE(i,j)=Hvom(i,j,k)*                                &
     &                      SOURCES(ng)%Tsrc(is,k,itrc)
# ifdef MASKING
                  ELSE
                    IF ((rmask(i,j  ).eq.0.0_r8).and.                   &
     &                  (rmask(i,j-1).eq.1.0_r8)) THEN
                      FE(i,j)=Hvom(i,j,k)*t(i,j-1,k,3,itrc)
                    ELSE IF ((rmask(i,j  ).eq.1.0_r8).and.              &
     &                       (rmask(i,j-1).eq.0.0_r8)) THEN
                      FE(i,j)=Hvom(i,j,k)*t(i,j  ,k,3,itrc)
                    END IF
# endif
                  END IF
                END IF
              END IF
            END DO
          END IF

# if defined NESTING && !defined ONE_WAY
!
!  If refinement grids, extract tracer horizontal advection fluxes
!  (Hz*u*T/n, Hz*v*T/m) at the grid contact boundary (physical
!  domain period) to be used in two-way nesting.
!
         IF (RefinedGrid(ng)) THEN
           DO cr=1,Ncontact
             dg=Rcontact(cr)%donor_grid
             rg=Rcontact(cr)%receiver_grid
             IF (ng.eq.rg) THEN
               CALL bry_fluxes (dg, rg, cr, iNLM, tile,                 &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          ILB, IUB, JLB, JUB,                     &
     &                          dt(ng), FX, FE,                         &
     &                          BRY_CONTACT(iwest, cr)%Tflux(:,k,itrc), &
     &                          BRY_CONTACT(ieast, cr)%Tflux(:,k,itrc), &
     &                          BRY_CONTACT(isouth,cr)%Tflux(:,k,itrc), &
     &                          BRY_CONTACT(inorth,cr)%Tflux(:,k,itrc))
             END IF
           END DO
         END IF
# endif
!
# ifdef TS_MPDATA
!  Time-step horizontal advection for intermediate diffusive tracer, Ta.
!  Advective fluxes have units of Tunits m3/s.  The new tracer has
!  units of m Tunits.
# else
!  Time-step horizontal advection term.  Advective fluxes have units
!  of Tunits m3/s.  The new tracer has units of m Tunits.
# endif
!
# ifdef TS_MPDATA
#  define I_RANGE IstrUm2,Iendp2i
#  define J_RANGE JstrVm2,Jendp2i
# else
#  define I_RANGE Istr,Iend
#  define J_RANGE Jstr,Jend
# endif
          DO j=J_RANGE
            DO i=I_RANGE
              cff=dt(ng)*pm(i,j)*pn(i,j)
              cff1=cff*(FX(i+1,j)-FX(i,j))
              cff2=cff*(FE(i,j+1)-FE(i,j))
              cff3=cff1+cff2
# ifdef TS_MPDATA
              Ta(i,j,k,itrc)=t(i,j,k,nnew,itrc)-cff3
# else
              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)-cff3
# endif
# ifdef DIAGNOSTICS_TS
#  ifdef TS_MPDATA
              Dhadv(i,j,iTxadv)=-cff1
              Dhadv(i,j,iTyadv)=-cff2
              Dhadv(i,j,iThadv)=-cff3
#  else
              DiaTwrk(i,j,k,itrc,iTxadv)=-cff1
              DiaTwrk(i,j,k,itrc,iTyadv)=-cff2
              DiaTwrk(i,j,k,itrc,iThadv)=-cff3
#  endif
# endif
            END DO
          END DO
# if defined DIAGNOSTICS_TS && defined TS_MPDATA
          DO j=Jstr,Jend
            DO i=Istr,Iend
              DiaTwrk(i,j,k,itrc,iTxadv)=Dhadv(i,j,iTxadv)
              DiaTwrk(i,j,k,itrc,iTyadv)=Dhadv(i,j,iTyadv)
              DiaTwrk(i,j,k,itrc,iThadv)=Dhadv(i,j,iThadv)
            END DO
          END DO
# endif
        END DO K_LOOP
      END DO T_LOOP
!
!-----------------------------------------------------------------------
!  Time-step vertical advection term.
!-----------------------------------------------------------------------
!
      DO j=J_RANGE
        DO itrc=1,NT(ng)

# if defined TS_SVADVECTION
!
!  Build conservative parabolic splines for the vertical derivatives
!  "FC" of the tracer.  Then, the interfacial "FC" values are
!  converted to vertical advective flux.
!
          DO i=Istr,Iend
#  ifdef NEUMANN
            FC(i,0)=1.5_r8*t(i,j,1,3,itrc)
            CF(i,1)=0.5_r8
#  else
            FC(i,0)=2.0_r8*t(i,j,1,3,itrc)
            CF(i,1)=1.0_r8
#  endif
          END DO
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              cff=1.0_r8/(2.0_r8*Hz(i,j,k)+                             &
     &                    Hz(i,j,k+1)*(2.0_r8-CF(i,k)))
              CF(i,k+1)=cff*Hz(i,j,k)
              FC(i,k)=cff*(3.0_r8*(Hz(i,j,k  )*t(i,j,k+1,3,itrc)+       &
     &                             Hz(i,j,k+1)*t(i,j,k  ,3,itrc))-      &
     &                     Hz(i,j,k+1)*FC(i,k-1))
            END DO
          END DO
          DO i=Istr,Iend
#  ifdef NEUMANN
            FC(i,N(ng))=(3.0_r8*t(i,j,N(ng),3,itrc)-FC(i,N(ng)-1))/     &
     &                  (2.0_r8-CF(i,N(ng)))
#  else
            FC(i,N(ng))=(2.0_r8*t(i,j,N(ng),3,itrc)-FC(i,N(ng)-1))/     &
     &                  (1.0_r8-CF(i,N(ng)))
#  endif
          END DO
          DO k=N(ng)-1,0,-1
            DO i=Istr,Iend
              FC(i,k)=FC(i,k)-CF(i,k+1)*FC(i,k+1)
              FC(i,k+1)=W(i,j,k+1)*FC(i,k+1)
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,N(ng))=0.0_r8
            FC(i,0)=0.0_r8
          END DO

# elif defined TS_A4VADVECTION
!
!  Fourth-order, Akima vertical advective flux.
!
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              FC(i,k)=t(i,j,k+1,3,itrc)-                                &
     &                t(i,j,k  ,3,itrc)
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,0)=FC(i,1)
            FC(i,N(ng))=FC(i,N(ng)-1)
          END DO
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff=2.0_r8*FC(i,k)*FC(i,k-1)
              IF (cff.gt.eps) THEN
                CF(i,k)=cff/(FC(i,k)+FC(i,k-1))
              ELSE
                CF(i,k)=0.0_r8
              END IF
            END DO
          END DO
          cff1=1.0_r8/3.0_r8
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              FC(i,k)=W(i,j,k)*                                         &
     &                0.5_r8*(t(i,j,k  ,3,itrc)+                        &
     &                        t(i,j,k+1,3,itrc)-                        &
     &                        cff1*(CF(i,k+1)-CF(i,k)))
            END DO
          END DO
          DO i=Istr,Iend
#  ifdef SED_MORPH
            FC(i,0)=W(i,j,0)*t(i,j,1,3,itrc)
#  else
            FC(i,0)=0.0_r8
#  endif
            FC(i,N(ng))=0.0_r8
          END DO
# elif defined TS_C2VADVECTION
!
!  Second-order, central differences vertical advective flux.
!
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              FC(i,k)=W(i,j,k)*                                         &
     &                0.5_r8*(t(i,j,k  ,3,itrc)+                        &
     &                        t(i,j,k+1,3,itrc))
            END DO
          END DO
          DO i=Istr,Iend
#  ifdef SED_MORPH
            FC(i,0)=W(i,j,0)*t(i,j,1,3,itrc)
#  else
            FC(i,0)=0.0_r8
#  endif
            FC(i,N(ng))=0.0_r8
          END DO
# elif defined TS_MPDATA
!
!  First_order, upstream differences vertical advective flux.
!
          DO i=I_RANGE
            DO k=1,N(ng)-1
              cff1=MAX(W(i,j,k),0.0_r8)
              cff2=MIN(W(i,j,k),0.0_r8)
              FC(i,k)=cff1*t(i,j,k  ,3,itrc)+                           &
     &                cff2*t(i,j,k+1,3,itrc)
            END DO
#  ifdef SED_MORPH
            FC(i,0)=W(i,j,0)*t(i,j,1,3,itrc)
#  else
            FC(i,0)=0.0_r8
#  endif
            FC(i,N(ng))=0.0_r8
          END DO
# else
!
!  Fourth-order, central differences vertical advective flux.
!
          cff1=0.5_r8
          cff2=7.0_r8/12.0_r8
          cff3=1.0_r8/12.0_r8
          DO k=2,N(ng)-2
            DO i=Istr,Iend
              FC(i,k)=W(i,j,k)*                                         &
     &                (cff2*(t(i,j,k  ,3,itrc)+                         &
     &                       t(i,j,k+1,3,itrc))-                        &
     &                 cff3*(t(i,j,k-1,3,itrc)+                         &
     &                       t(i,j,k+2,3,itrc)))
            END DO
          END DO
          DO i=Istr,Iend
#  ifdef SED_MORPH
            FC(i,0)=W(i,j,0)*2.0_r8*                                    &
     &              (cff2*t(i,j,1,3,itrc)-                              &
     &               cff3*t(i,j,2,3,itrc))
#  else
            FC(i,0)=0.0_r8
#  endif
            FC(i,1)=W(i,j,1)*                                           &
     &              (cff1*t(i,j,1,3,itrc)+                              &
     &               cff2*t(i,j,2,3,itrc)-                              &
     &               cff3*t(i,j,3,3,itrc))
            FC(i,N(ng)-1)=W(i,j,N(ng)-1)*                               &
     &                    (cff1*t(i,j,N(ng)  ,3,itrc)+                  &
     &                     cff2*t(i,j,N(ng)-1,3,itrc)-                  &
     &                     cff3*t(i,j,N(ng)-2,3,itrc))
            FC(i,N(ng))=0.0_r8
          END DO
# endif
!
# ifdef TS_MPDATA
!  Time-step vertical advection for intermediate diffusive tracer, Ta
!  (Tunits).
# else
#  ifdef SPLINES
!  Time-step vertical advection term (Tunits).
#  else
!  Time-step vertical advection term (m Tunits).
#  endif
# endif
# ifdef DIAGNOSTICS_TS
!  Convert units of tracer diagnostic terms to Tunits.
# endif
!
          DO i=I_RANGE
            CF(i,0)=dt(ng)*pm(i,j)*pn(i,j)
          END DO
!
!  Apply mass point sources (volume vertical influx), if any.
!
          IF (LwSrc(ng)) THEN
            DO is=1,Nsrc(ng)
              i=SOURCES(ng)%Isrc(is)
#  ifdef TS_MPDATA
              IF (((IstrUm2.le.i).and.(i.le.Iendp2i)).and.              &
     &            ((JstrVm2.le.j).and.(j.le.Jendp2i))) THEN
#  else
              IF (((IstrR.le.i).and.(i.le.IendR)).and.                  &
     &            ((JstrR.le.j).and.(j.le.JendR))) THEN
#  endif
                IF (j.eq.SOURCES(ng)%Jsrc(is)) THEN
                  DO k=1,N(ng)-1
                    FC(i,k)=FC(i,k)+0.5_r8*                             &
     &                      (SOURCES(ng)%Qsrc(is,k  )*                  &
                             SOURCES(ng)%Tsrc(is,k  ,itrc)+             &
     &                       SOURCES(ng)%Qsrc(is,k+1)*                  &
     &                       SOURCES(ng)%Tsrc(is,k+1,itrc))
                  END DO
                END IF
              END IF
            END DO
          END IF
!
          DO k=1,N(ng)
            DO i=I_RANGE
              cff1=CF(i,0)*(FC(i,k)-FC(i,k-1))
# ifdef TS_MPDATA
              Ta(i,j,k,itrc)=(Ta(i,j,k,itrc)-cff1)*oHz(i,j,k)
#  ifdef DIAGNOSTICS_TS
              Dvadv(i,j,k,itrc)=-cff1
#  endif
# else
              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)-cff1
#  ifdef SPLINES
              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)*oHz(i,j,k)
#  endif
#  ifdef DIAGNOSTICS_TS
              DiaTwrk(i,j,k,itrc,iTvadv)=-cff1
              DO idiag=1,NDT
                DiaTwrk(i,j,k,itrc,idiag)=DiaTwrk(i,j,k,itrc,idiag)*    &
     &                                    oHz(i,j,k)
              END DO
#  endif
# endif
            END DO
          END DO
        END DO
# undef I_RANGE
# undef J_RANGE
# ifdef TS_MPDATA
      END DO
!
!-----------------------------------------------------------------------
!  Compute anti-diffusive velocities to corrected advected tracers
!  using MPDATA recursive method.  Notice that pipelined J-loop ended.
!-----------------------------------------------------------------------
!
      DO itrc=1,NT(ng)
        CALL mpdata_adiff_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
#  ifdef MASKING
     &                          rmask, umask, vmask,                    &
#  endif
#  ifdef WET_DRY
     &                          rmask_wet, umask_wet, vmask_wet,        &
#  endif
     &                          pm, pn, omn, om_u, on_v,                &
     &                          z_r, oHz,                               &
     &                          Huon, Hvom, W,                          &
     &                          t(:,:,:,3,itrc),                        &
     &                          Ta(:,:,:,itrc),  Ua, Va, Wa)
!
!  Compute anti-diffusive corrected advection fluxes.
!
        DO k=1,N(ng)
          DO j=Jstr,Jend
            DO i=Istr,Iend+1
              cff1=MAX(Ua(i,j,k),0.0_r8)
              cff2=MIN(Ua(i,j,k),0.0_r8)
              FX(i,j)=(cff1*Ta(i-1,j,k,itrc)+                           &
     &                 cff2*Ta(i  ,j,k,itrc))*                          &
     &                0.5_r8*(Hz(i,j,k)+Hz(i-1,j,k))*on_u(i,j)
            END DO
          END DO
          DO j=Jstr,Jend+1
            DO i=Istr,Iend
              cff1=MAX(Va(i,j,k),0.0_r8)
              cff2=MIN(Va(i,j,k),0.0_r8)
              FE(i,j)=(cff1*Ta(i,j-1,k,itrc)+                           &
     &                 cff2*Ta(i,j  ,k,itrc))*                          &
     &                0.5_r8*(Hz(i,j,k)+Hz(i,j-1,k))*om_v(i,j)
            END DO
          END DO
!
!  Time-step corrected horizontal advection (Tunits m).
!
          DO j=Jstr,Jend
            DO i=Istr,Iend
              cff=dt(ng)*pm(i,j)*pn(i,j)
              cff1=cff*(FX(i+1,j)-FX(i,j))
              cff2=cff*(FE(i,j+1)-FE(i,j))
              cff3=cff1+cff2
              t(i,j,k,nnew,itrc)=Ta(i,j,k,itrc)*Hz(i,j,k)-cff3
#  ifdef DIAGNOSTICS_TS
              DiaTwrk(i,j,k,itrc,iTxadv)=DiaTwrk(i,j,k,itrc,iTxadv)-    &
     &                                   cff1
              DiaTwrk(i,j,k,itrc,iTyadv)=DiaTwrk(i,j,k,itrc,iTyadv)-    &
     &                                   cff2
              DiaTwrk(i,j,k,itrc,iThadv)=DiaTwrk(i,j,k,itrc,iThadv)-    &
     &                                   cff3
#  endif
            END DO
          END DO
        END DO
!
!  Compute anti-diffusive corrected vertical advection flux.
!
        DO j=Jstr,Jend
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              cff1=MAX(Wa(i,j,k),0.0_r8)
              cff2=MIN(Wa(i,j,k),0.0_r8)
              FC(i,k)=cff1*Ta(i,j,k  ,itrc)+                            &
     &                cff2*Ta(i,j,k+1,itrc)
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,0)=0.0_r8
            FC(i,N(ng))=0.0_r8
          END DO
!
!  Time-step corrected vertical advection (Tunits).
#  ifdef DIAGNOSTICS_TS
!  Convert units of tracer diagnostic terms to Tunits.
#  endif
!
          DO i=Istr,Iend
            CF(i,0)=dt(ng)*pm(i,j)*pn(i,j)
          END DO
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff1=CF(i,0)*(FC(i,k)-FC(i,k-1))
              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)-cff1
#  ifdef DIAGNOSTICS_TS
              DiaTwrk(i,j,k,itrc,iTvadv)=Dvadv(i,j,k,itrc)-             &
     &                                   cff1
              DO idiag=1,NDT
                DiaTwrk(i,j,k,itrc,idiag)=DiaTwrk(i,j,k,itrc,idiag)*    &
     &                                    oHz(i,j,k)
              END DO
#  endif
            END DO
          END DO
        END DO
      END DO

!    zhousy !!!!
      DO j=Jstr,Jend            
        DO i=Istr,Iend
          IF (z_r(i,j,1)< -500 .and. latr(i,j)>27 .and. latr(i,j)<47    &
     &    .and. k_mld(i,j)<100 .and. z_r(i,j,k_mld(i,j))/z_r(i,j,1)<0.7 &
     &    .and. lonr(i,j)> -72 .and. lonr(i,j)< -35) THEN
            !print *, z_r(i,j,k_mld(i,j))
            DO k=k_mld(i,j),N(ng)-1
              input(1,1) = (t(i,j,k,nstp,1)-16.0148)/2.0323
              input(1,2) = (dtdx(i,j,k)*1000+0.0074944)/0.025787
              input(1,3) = (dtde(i,j,k)*1000+0.0006884)/0.027229175
              input(1,4) = (u(i,j,k,nstp)-0.0867963)/0.46863
              input(1,5) = (dudx(i,j,k)*1000-0.000032812)/0.009736
              input(1,6) = (dude_r(i,j,k)*1000-0.0002731)/0.005895
              input(1,7) = (v(i,j,k,nstp)+0.00485236)/0.548829
              input(1,8) = (dvdx_r(i,j,k)*1000+0.00024663)/0.0057913
              input(1,9) = (dvde(i,j,k)*1000+0.000814153)/0.0105653
              input(1,10) = z_r(i,j,k)/z_r(i,j,k_mld(i,j))
              wb1_hs = matmul(input, w1_hs)+b1_hs
              where (wb1_hs<0) wb1_hs = 0.3*wb1_hs
              wb2_hs = matmul(wb1_hs, w2_hs)+b2_hs
              where (wb2_hs<0) wb2_hs = 0.3*wb2_hs
              wb3_hs = matmul(wb2_hs, w3_hs)+b3_hs
              where (wb3_hs<0) wb3_hs = 0.3*wb3_hs
              wb4_hs = matmul(wb3_hs, w4_hs)+b4_hs
              wtemp_01(1,1)=wb4_hs(1,1)*0.000031706+0.00003252756
              !print *, "wtemp_01:",wb4_hs(1,1),wtemp_01(1,1)
              input(1,1) = (t(i,j,k-1,nstp,1)-16.0148)/2.0323
              input(1,2) = (dtdx(i,j,k-1)*1000+0.0074944)/0.025787
              input(1,3) = (dtde(i,j,k-1)*1000+0.0006884)/0.027229175
              input(1,4) = (u(i,j,k-1,nstp)-0.0867963)/0.46863
              input(1,5) = (dudx(i,j,k-1)*1000-0.000032812)/0.009736
              input(1,6) = (dude_r(i,j,k-1)*1000-0.0002731)/0.005895
              input(1,7) = (v(i,j,k-1,nstp)+0.00485236)/0.548829
              input(1,8) = (dvdx_r(i,j,k-1)*1000+0.00024663)/0.0057913
              input(1,9) = (dvde(i,j,k-1)*1000+0.000814153)/0.0105653
              input(1,10) = z_r(i,j,k-1)/z_r(i,j,k_mld(i,j))
              wb1_hs = matmul(input, w1_hs)+b1_hs
              where (wb1_hs<0) wb1_hs = 0.3*wb1_hs
              wb2_hs = matmul(wb1_hs, w2_hs)+b2_hs
              where (wb2_hs<0) wb2_hs = 0.3*wb2_hs
              wb3_hs = matmul(wb2_hs, w3_hs)+b3_hs
              where (wb3_hs<0) wb3_hs = 0.3*wb3_hs
              wb4_hs = matmul(wb3_hs, w4_hs)+b4_hs
              !wtemp_02(1,1)=wb4_hs(1,1)
              wtemp_02(1,1)=wb4_hs(1,1)*0.000031706+0.00003252756
              !print *, "wtemp_02:",wb4_hs(1,1), wtemp_02(1,1) 
              !cff=abs(dvdx_r(i,j,k)*1000-dude_r(i,j,k)*1000)*10000
              tao(1,1)=sqrt((dudx(i,j,k)*1000-dvde(i,j,k)*1000)**2+     &
     &           (dvdx_r(i,j,k)*1000+dude_r(i,j,k)*1000)**2)/1e-04

              !print *, tao(1,1)
              
              !cff = sqrt(dbdx(i,j,k)**2+dbde(i,j,k)**2)
              
              !print *, wtemp_01(1,1)
              !print *, z_r(i,j,k),"m:  ",wtemp_01(1,1)
              !print *,"nnew:", nnew,"nstp:",nstp,"nhrs:",nrhs
              !print *, dt(ng), oHz(i,j,k)
              !print *, -dt(ng)*oHz(i,j,k)*                              &
              !                      (wtemp_01(1,1)-wtemp_02(1,1))
              !print *, "OHz: ", oHz(i,j,k), 1/Hz(i,j,k)
              !print *, latr(i,j), "--", lonr(i,j)

             
              !IF (z_r(i,j,k_mld(i,j))/z_r(i,j,1)>0.5) THEN              
              !  print *,                                                &
     !&          z_r(i,j,k_mld(i,j))/z_r(i,j,1),lonr(i,j),latr(i,j) 
              !END IF
              IF (tao(1,1)>200                                          &
     &        .and. abs(wtemp_01(1,1)) < 1E-03                          &
     &        .and. abs(wtemp_02(1,1)) < 1E-03 ) THEN
                !print *,lonr(i,j),latr(i,j),z_r(i,j,k) 
                !print *, -dt(ng)*oHz(i,j,k)*                            &
                !                    (wtemp_01(1,1)-wtemp_02(1,1))
                t(i,j,k,nnew,1)=t(i,j,k,nnew,1)-dt(ng)*                 &
     &              oHz(i,j,k)*(wtemp_01(1,1)-wtemp_02(1,1))
              END IF
              
              !print *, "t:",t(i,j,k,nnew,1),"H",t(i,j,k,nrhs,1)
              !print *, "U:",u(i,j,1,nstp), "H:",u(i,j,1,nrhs)
              !print *, "V:",v(i,j,1,nstp), "H:",v(i,j,1,nrhs)
              !print *, "dtdx, dtde:", pm(i,j), pn(i,j)
              !print *, "dudx, dude:", dudx(i,j,k), dude_r(i,j,k)
              !print *, "dvdx, dvde:", dvdx_r(i,j,k), dvde(i,j,k)
            END DO
          END IF
        END DO
      END DO


!
!  Start pipelined J-loop.
!
      DO j=Jstr,Jend
# endif /* TS_MPDATA */
!
!-----------------------------------------------------------------------
!  Time-step vertical diffusion term.
!-----------------------------------------------------------------------
!
        DO itrc=1,NT(ng)
          ltrc=MIN(NAT,itrc)

# if defined SPLINES && !defined TS_MPDATA
!
!  Use conservative, parabolic spline reconstruction of vertical
!  diffusion derivatives.  Then, time step vertical diffusion term
!  implicitly.
!
          cff1=1.0_r8/6.0_r8
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              FC(i,k)=cff1*Hz(i,j,k  )-                                 &
     &                dt(ng)*Akt(i,j,k-1,ltrc)*oHz(i,j,k  )
              CF(i,k)=cff1*Hz(i,j,k+1)-                                 &
     &                dt(ng)*Akt(i,j,k+1,ltrc)*oHz(i,j,k+1)
            END DO
          END DO
          DO i=Istr,Iend
            CF(i,0)=0.0_r8
            DC(i,0)=0.0_r8
          END DO
!
!  LU decomposition and forward substitution.
!
          cff1=1.0_r8/3.0_r8
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              BC(i,k)=cff1*(Hz(i,j,k)+Hz(i,j,k+1))+                     &
     &                dt(ng)*Akt(i,j,k,ltrc)*(oHz(i,j,k)+oHz(i,j,k+1))
              cff=1.0_r8/(BC(i,k)-FC(i,k)*CF(i,k-1))
              CF(i,k)=cff*CF(i,k)
              DC(i,k)=cff*(t(i,j,k+1,nnew,itrc)-t(i,j,k,nnew,itrc)-     &
     &                     FC(i,k)*DC(i,k-1))
            END DO
          END DO
!
!  Backward substitution.
!
          DO i=Istr,Iend
            DC(i,N(ng))=0.0_r8
          END DO
          DO k=N(ng)-1,1,-1
            DO i=Istr,Iend
              DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
            END DO
          END DO
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              DC(i,k)=DC(i,k)*Akt(i,j,k,ltrc)
              cff1=dt(ng)*oHz(i,j,k)*(DC(i,k)-DC(i,k-1))
              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)+cff1
#  ifdef DIAGNOSTICS_TS
              DiaTwrk(i,j,k,itrc,iTvdif)=DiaTwrk(i,j,k,itrc,iTvdif)+    &
     &                                   cff1
#  endif
            END DO
          END DO
# else
!
!  Compute off-diagonal coefficients FC [lambda*dt*Akt/Hz] for the
!  implicit vertical diffusion terms at future time step, located
!  at horizontal RHO-points and vertical W-points.
!  Also set FC at the top and bottom levels.
!
          cff=-dt(ng)*lambda
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              cff1=1.0_r8/(z_r(i,j,k+1)-z_r(i,j,k))
              FC(i,k)=cff*cff1*Akt(i,j,k,ltrc)
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,0)=0.0_r8
            FC(i,N(ng))=0.0_r8
          END DO
!
!  Compute diagonal matrix coefficients BC and load right-hand-side
!  terms for the tracer equation into DC.
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              BC(i,k)=Hz(i,j,k)-FC(i,k)-FC(i,k-1)
              DC(i,k)=t(i,j,k,nnew,itrc)
            END DO
          END DO
!
!  Solve the tridiagonal system.
!
          DO i=Istr,Iend
            cff=1.0_r8/BC(i,1)
            CF(i,1)=cff*FC(i,1)
            DC(i,1)=cff*DC(i,1)
          END DO
          DO k=2,N(ng)-1
            DO i=Istr,Iend
              cff=1.0_r8/(BC(i,k)-FC(i,k-1)*CF(i,k-1))
              CF(i,k)=cff*FC(i,k)
              DC(i,k)=cff*(DC(i,k)-FC(i,k-1)*DC(i,k-1))
            END DO
          END DO
!
!  Compute new solution by back substitution.
!
          DO i=Istr,Iend
#  ifdef DIAGNOSTICS_TS
             cff1=t(i,j,N(ng),nnew,itrc)*oHz(i,j,N(ng))
#  endif
             DC(i,N(ng))=(DC(i,N(ng))-FC(i,N(ng)-1)*DC(i,N(ng)-1))/     &
     &                   (BC(i,N(ng))-FC(i,N(ng)-1)*CF(i,N(ng)-1))
             t(i,j,N(ng),nnew,itrc)=DC(i,N(ng))
#  ifdef DIAGNOSTICS_TS
             DiaTwrk(i,j,N(ng),itrc,iTvdif)=                            &
     &                              DiaTwrk(i,j,N(ng),itrc,iTvdif)+     &
     &                              t(i,j,N(ng),nnew,itrc)-cff1
#  endif
          END DO
          DO k=N(ng)-1,1,-1
            DO i=Istr,Iend
#  ifdef DIAGNOSTICS_TS
              cff1=t(i,j,k,nnew,itrc)*oHz(i,j,k)
#  endif
              DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
              t(i,j,k,nnew,itrc)=DC(i,k)
#  ifdef DIAGNOSTICS_TS
              DiaTwrk(i,j,k,itrc,iTvdif)=DiaTwrk(i,j,k,itrc,iTvdif)+    &
     &                                   t(i,j,k,nnew,itrc)-cff1
#  endif
            END DO
          END DO
# endif
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Apply lateral boundary conditions and, if appropriate, nudge
!  to tracer data and apply Land/Sea mask.
!-----------------------------------------------------------------------
!
!  Initialize tracer counter index. The "tclm" array is only allocated
!  to the NTCLM fields that need to be processed. This is done to
!  reduce memory.
!
      ic=0
!
      DO itrc=1,NT(ng)
!
!  Set compact reduced memory tracer index for nudging coefficients and
!  climatology arrays.
!
        IF (LtracerCLM(itrc,ng).and.LnudgeTCLM(itrc,ng)) THEN
          ic=ic+1
        END IF
!
!  Set lateral boundary conditions.
!
        CALL t3dbc_tile (ng, tile, itrc, ic,                            &
     &                   LBi, UBi, LBj, UBj, N(ng), NT(ng),             &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nstp, nnew,                                    &
     &                   t)
!
!  Nudge towards tracer climatology.
!
        IF (LtracerCLM(itrc,ng).and.LnudgeTCLM(itrc,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)+                  &
     &                             dt(ng)*                              &
     &                             CLIMA(ng)%Tnudgcof(i,j,k,ic)*        &
     &                             (CLIMA(ng)%tclm(i,j,k,ic)-           &
     &                              t(i,j,k,nnew,itrc))
              END DO
            END DO
          END DO
        END IF

# ifdef MASKING
!
!  Apply Land/Sea mask.
!
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)*rmask(i,j)
            END DO
          END DO
        END DO
# endif
# ifdef DIAGNOSTICS_TS
!
!  Compute time-rate-of-change diagnostic term.
!
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              DiaTwrk(i,j,k,itrc,iTrate)=t(i,j,k,nnew,itrc)-            &
     &                                   DiaTwrk(i,j,k,itrc,iTrate)
!!            DiaTwrk(i,j,k,itrc,iTrate)=t(i,j,k,nnew,itrc)-            &
!!   &                                   t(i,j,k,nstp,itrc)
            END DO
          END DO
        END DO
# endif
!
!  Apply periodic boundary conditions.
!
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          CALL exchange_r3d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            t(:,:,:,nnew,itrc))
        END IF
      END DO
# ifdef DISTRIBUTE
!
!  Exchange boundary data.
!
      CALL mp_exchange4d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),      &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    t(:,:,:,nnew,:))
# endif
# if defined FLOATS && defined FLOAT_VWALK
!
!-----------------------------------------------------------------------
!  Compute vertical gradient in vertical T-diffusion coefficient for
!  floats random walk.
!-----------------------------------------------------------------------
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          DO k=1,N(ng)
            dAktdz(i,j,k)=(Akt(i,j,k,1)-Akt(i,j,k-1,1))/Hz(i,j,k)
          END DO
        END DO
      END DO
!
!  Apply periodic boundary conditions.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          dAktdz)
      END IF

#  ifdef DISTRIBUTE
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    dAktdz)
#  endif
# endif

      RETURN
      END SUBROUTINE step3d_t_tile
#endif
      END MODULE step3d_t_mod
