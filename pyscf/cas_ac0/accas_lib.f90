! Copyright 2023 Quantinuum
!
! You may not use this file except in compliance with the Licence.
! You may obtain a copy of the Licence in the LICENCE file accompanying
! these documents
!
! This code is derived from GAMMCOR
! https://github.com/pernalk/GAMMCOR
! Caution: GAMMCOR is released under GPL v 3.0
! hence this code can be distributed only under GPL
! Isolate in a separate executable before distribution!

module accas_lib
implicit none
contains

function get_two_el_size(NBas) result(two_el_size)

    integer(8), intent(in) :: nbas ! Cono
    integer(8)    ::  NInte1 ! Cono
    integer(8) :: two_el_size
    NInte1 = NBas*(NBas+1)/2
    two_el_size = NInte1*(NInte1+1)/2
end function get_two_el_size

function get_two_el(eri_matrix, two_el_size, nbas) result(twoel)
    implicit none
    integer(4), intent(in) :: nbas
    integer(8), intent(in) :: two_el_size
    real(8), dimension(nbas,nbas,nbas,nbas), intent(in) :: eri_matrix
    integer(4)    ::  I,J,K,L
    real(8), dimension(two_el_size) :: twoel

    ! print *, 'two_el_size=', two_el_size



    Do I=1,nbas
      Do J=I,nbas
        Do K=1,nbas
          Do L=K,nbas
            ! print *, 'inner loop'
            twoel(NAddr3(I,J,K,L)) = eri_matrix(I,J,K,L)
          end do
        end do
       end do
    end do
end function get_two_el


subroutine get_rdm2_act(NAct, NRDM2Act, rdm2, RDM2Act)
      integer, intent(in) :: NAct
      integer, intent(in) :: NRDM2Act
      real(8), dimension(NAct,NAct,NAct,NAct), intent(in) :: rdm2
      real(8), dimension(NRDM2Act), intent(out) :: RDM2Act
      integer :: i,j,k,l

      RDM2Act(1:NRDM2Act)=0.0D0

      do I = 1, NAct
          do J = 1, NAct
              do K = 1, NAct
                  do L = 1, NAct
                     RDM2Act(NAddrRDM(J,L,I,K,NAct))=rdm2(I,J,K,L) * 0.5
                     !print *, 'RDM2Act(', NAddrRDM(J,L,I,K,NAct), ')=', rdm2(I,J,K,L)
                     !print *, 'I=', I, 'J=', J, 'K=', K, 'L=', L
                  end do
              end do
          end do
      end do
end subroutine get_rdm2_act


Integer Function GetNRDM2Act(NAct)
    integer, intent(in) :: NAct
    GetNRDM2Act=NAct**2*(NAct**2+1)/2
    Return
End



SUBROUTINE diag8(a, nmax, n, eigenvalues)
implicit none
integer, intent(in) :: nmax,n
real(8), dimension(nmax,n), intent(inout) :: a
real(8), dimension(n), intent(out)  :: eigenvalues
real(8), allocatable  :: work(:,:)
integer , allocatable  :: iwork(:)
integer                        :: lwork, liwork, info
integer                        :: i,j



     lwork = 1
      liwork = 1
      allocate (work(lwork,1),iwork(liwork))

      ! Determine optimal size for work arrays
      lwork = -1
      liwork = -1
      CALL dsyevd( 'V', 'L', n, a, nmax, eigenvalues, work, lwork,  &
    iwork, liwork, info )
      if (info < 0) then
         deallocate (work,iwork)
         PRINT *, 'Diag8: ',  &
      ' DSYEVD: the ',-info,'-th argument had an illegal value'
         stop 2
      endif

      lwork  = max(int(work(1,1)), 2*n*n + 6*n+ 1)
      liwork = max(iwork(1), 5*n + 3)

!     /!\ liwork becomes negative when > 2147483648 (integer*4 overflow)
      if ((liwork < 0) .or. (lwork < 0)) then
         deallocate (work,iwork)
         PRINT *, 'Diag8: ', ' Required work space too large'
         stop 3
      endif

      ! Allocate temporary arrays
      deallocate (work,iwork)
      allocate (work(lwork,1),iwork(liwork))

      ! Diagonalize
      CALL dsyevd( 'V', 'L', n, a, nmax, eigenvalues, work, lwork,  &
    iwork, liwork, info)

      deallocate(work,iwork)

      if (info < 0) then
         deallocate (work,iwork)
         PRINT *, 'Diag8:',  &
      ': DSYEVD: the ',-info,'-th argument had an illegal value'
         stop 2
      else if( info > 0 ) then
         deallocate (work,iwork)
         write(*,*)'DSYEVD Failed'
         stop 1
      end if

      ! Transpose result
      allocate(work(size(A,1),size(A,2)))
      work = A
      do j=1,n
        do i=1,n
          A(i,j) = work(j,i)
        enddo
      enddo
      deallocate(work)

! Suggestion :
!   Remove the transposition, and use the eigenvectors as column vectors
!   outside of this routine

END SUBROUTINE diag8

function Naddr3(IAddr1,IAddr2,IAddr3,IAddr4) result(NAddr3_res)
! POINTER FOR TWO-ELECTRON INTEGRALS
! parameter(Zero=0.0D0,One=1.0D0,Two=2.0D0)
 integer(4)    :: IAddr1, IAddr2, IAddr3, IAddr4
 integer(8) :: IAddr12, IAddr34
 integer(8) :: NAddr3_res

! CHANGE THE ORDER IF NECESSARY

 IAddr12 = Max(IAddr1,IAddr2)*(Max(IAddr1,IAddr2)-1)/2 + &
          Min(IAddr2,IAddr1)
 IAddr34 = Max(IAddr3,IAddr4)*(Max(IAddr3,IAddr4)-1)/2 + &
          Min(IAddr3,IAddr4)

! GET THE POSITION OF THE ELEMEMT (12|34)

 NAddr3_res = Max(IAddr12,IAddr34)*(Max(IAddr12,IAddr34)-1)/2 + &
          Min(IAddr12,IAddr34)

end function NAddr3

SUBROUTINE accas(etot,ecorr,enuc,twono,ure,occ,xone,  &
    nbasis,ninte1,ninte2,rdm2act, n_electron, naccas, n_act_ele)
!     A ROUTINE FOR COMPUTING ELECTRONIC ENERGY USING ERPA TRANSITION
!     DENSITY MATRIX ELEMENTS
REAL(8), INTENT(OUT)                     :: etot, ecorr
REAL(8), INTENT(IN )                     :: enuc
REAL(8), INTENT(IN )                     :: twono(ninte2)
! twono=TwoEl - 2-electron integrals in NO representation
REAL(8), INTENT(IN )                     :: ure(nbasis,nbasis)

REAL(8), INTENT(IN)                         :: occ(nbasis)
REAL(8), INTENT(IN )                     :: xone(ninte1)

!XOne(I) 1-electron integrals (in NO representation)

INTEGER, INTENT(IN)                      :: nbasis
!NBasis - number of basis functions
INTEGER, INTENT(IN )                  :: ninte1
 !NInte1 - dimension of 1-electron integrals
INTEGER, INTENT(IN )                  :: ninte2
!NInte2 - dimension of 2-electron integrals

REAL(8), dimension(:), INTENT(IN )                     :: rdm2act
integer, intent(in) :: naccas , n_act_ele, n_electron
!naccas number of active CAS orbitals

real(8) :: ThrAct, ThrSelAct, ThrVirt, ThrQVirt


!


REAL(8), PARAMETER :: zero=0.d0
REAL(8), PARAMETER :: half=0.5D0
REAL(8), PARAMETER :: one=1.d0
REAL(8), PARAMETER :: two=2.d0
integer :: i,j, icount, ij, ind1,ind2,ind, ndimx
integer ::  iflcore
integer :: nele ! number of electron pairs
real(8) :: xele !half of number of electrons


integer, dimension(nbasis*(nbasis-1)/2) :: indx
integer, dimension(2,nbasis*(nbasis-1)/2) :: indn
integer, dimension(nbasis) :: indaux
integer, dimension(nbasis,nbasis) :: ipair

!System%NELE = (System%ZNucl - System%Charge)/2
nele= n_electron /2
XELE = n_electron /2.0d0




iflcore = 1

! transform 1-electron integrals to MOs: (mattr is in matvec)
!this is already done!!!
!!Call MatTr(XOne,UNOAO,NBasis)

!     CONSTRUCT LOOK-UP TABLES

DO i=1,nele
  indaux(i)=0
END DO
DO i=1+nele,nbasis
  indaux(i)=2
END DO

icount=0

DO i=1,nbasis
  IF(occ(i) < one.AND.occ(i) /= zero) THEN
    indaux(i)=1
    !WRITE(6,'(X," Active Orbital: ",I4,E14.4)') i, occ(i)
    icount=icount+1
  END IF
END DO

!WRITE(6,'(/,X," In ACCAS: Active Orbitals ",I4,/)')icount

ipair(1:nbasis,1:nbasis)=0

! ThrAct for active geminals selection
 ThrAct    = 0.992d0
 ThrSelAct = 1.d-8
! ThrVirt for reduction of virtual orbs
 ThrVirt   = 1.d-6
! ThrQVirt for quasi-virtual orbs
ThrQVirt  = 1.d-7
!WRITE(*,'(2x,a,4x,2e15.5)') 'Threshold for quasi-degeneracy ', thrselact

!WRITE(*,'(2x,a,2e15.5)') 'Threshold for quasi-virtual orbital', thrqvirt

ij=0
ind=0
DO i=1,nbasis
  DO j=1,i-1
    ij=ij+1
    IF(indaux(i)+indaux(j) /= 0.AND.indaux(i)+indaux(j) /= 4) THEN
      IF((indaux(i) == 1).AND.(indaux(j) == 1)  &
            .AND.(ABS(occ(i)-occ(j))/occ(i) < thrselact) ) THEN
        !WRITE(6,'(2X,"Discarding nearly degenerate pair ",2I4)')i,j
      ELSE
!     If IFlCore=0 do not include core (inactive) orbitals
        IF((iflcore == 1).OR.  &
              (iflcore == 0.AND.occ(i) /= one.AND.occ(j) /= one)) THEN

          IF(ABS(occ(i)+occ(j)-two) > 1.d-10.AND.  &
                ABS(occ(i)+occ(j)) > thrqvirt) THEN
            ind=ind+1
            indx(ind)=ind
            indn(1,ind)=i
            indn(2,ind)=j
            ipair(i,j)=1
            ipair(j,i)=1
          END IF
        END IF
      END IF
    END IF
  END DO
END DO

ndimx=ind
!WRITE(6,'(/2X,"Number of pairs reduced to:",I6)')ind
!WRITE(6,'(2X,"Accepted pairs read:")')
DO i=1,ind
  ind1=indn(1,i)
  ind2=indn(2,i)
  !WRITE(6,'(2X,3I5,2E14.4)')i,ind1,ind2,occ(ind1),occ(ind2)
END DO
!WRITE(6,'(/," *** ADIABATIC CONNECTION CALCULATIONS ***",/)')

!IF(iflcore == 0) THEN
  !WRITE(6,'(/," *** IFlCore=0: Inactive orbitals (n_p=1)  &
  !  excluded from ERPA correlation ***",/)')
!ELSE
  !WRITE(6,'(/," *** IFlCore=1: Inactive orbitals (n_p=1)  &
  !included in ERPA correlation ***",/)')
!END IF





CALL ac0cas(ecorr,etot,twono,occ,ure,xone,  &
          indn,ndimx,nbasis,ninte1,ninte2, rdm2act, naccas, n_electron, n_act_ele)
    !WRITE (6,'(/,X,''ECASSCF+ENUC, AC0-CORR, AC0-CASSCF '',4X,3F15.8)')  &
     !  etot+enuc,ecorr,etot+enuc+ecorr
etot=etot+enuc+ecorr
END SUBROUTINE accas



SUBROUTINE ab0element(abpl,abmin,ir,is,ipp,iqq,occ,hno,igfact,  &
    twono,auxi,auxio,wmat,rdm2act,c,ind1,ind2,nact,nrdm2act,ninte1, ninte2,nbasis, igem) !added igem

!     FOR A GIVEN SET OF INDICES IR,IS,IPP,IQQ RETURNS
!     VALUES OF ABPLUS AND ABMIN MATRICES FOR ALPHA=0

real(8), intent(out) :: abpl, abmin
INTEGER, INTENT(IN ) :: nbasis,ninte1, ninte2, nrdm2act, nact
real(8), dimension(nbasis), intent(in) :: occ, c
integer, dimension(nbasis), intent(in) :: ind1, ind2
real(8), dimension(ninte1), intent(in) :: hno, auxi, auxio
real(8), dimension(ninte2), intent(in) :: twono
integer, dimension(ninte2), intent(in) :: igfact
real(8), dimension(nbasis, nbasis), intent(in) :: wmat
real(8), dimension(nrdm2act), intent(in) :: rdm2act
integer, dimension(nbasis), intent(in) :: igem
integer, intent(in) :: IR,IS,IPP,IQQ
real(8) :: arspq, auxpqrs, auxtwopqrs
integer :: ip, iq, it, itt, ipq, iu, iuu, ipr, iqp, iqs

REAL(8), PARAMETER :: zero=0.d0
REAL(8), PARAMETER :: half=0.5D0
REAL(8), PARAMETER :: one=1.d0
REAL(8), PARAMETER :: two=2.d0
real(8), PARAMETER :: three=3.d0
real(8), PARAMETER :: four=4.d0


abpl=zero
abmin=zero

DO ip=iqq,ipp,ipp-iqq
  DO iq=iqq,ipp,ipp-iqq
    IF(ip /= iq) THEN

      IF(ip > iq) ipq=(ip**2-3*ip+2)/2+iq
      IF(iq > ip) iqp=(iq**2-3*iq+2)/2+ip

      iqs=(MAX(is,iq)*(MAX(is,iq)-1))/2+MIN(is,iq)
      ipr=(MAX(ir,ip)*(MAX(ir,ip)-1))/2+MIN(ir,ip)

      arspq=zero

      IF(ip == ir) arspq=arspq+(occ(ip)-occ(is))*hno(iqs)
      IF(is == iq) arspq=arspq+(occ(iq)-occ(ir))*hno(ipr)

      auxtwopqrs=one
      IF(igfact(naddr3(ip,iq,ir,is)) == 0) auxtwopqrs=zero

      auxpqrs=zero
      IF(auxtwopqrs == one)  &
          auxpqrs=two*twono(naddr3(ip,iq,ir,is))-twono(naddr3(ip,ir,iq,is))

!     T1+T2

      IF(occ(ip)*occ(ir) /= zero) THEN

        IF(occ(ip) == one.OR.occ(ir) == one) THEN

          IF(auxtwopqrs == one) arspq=arspq+occ(ip)*occ(ir)*  &
              (two*twono(naddr3(ip,iq,ir,is))-twono(naddr3(ip,ir,iq,is)))
          IF(ip == ir) arspq=arspq+occ(ip)*auxi(iqs)

        ELSE

          IF(igem(is) == igem(iq)) THEN
            DO itt=1,nact
              DO iuu=1,nact
                it=ind1(itt)
                iu=ind1(iuu)

                IF(igfact(naddr3(is,iq,it,iu)) == 1) arspq=arspq  &
                    +twono(naddr3(is,iq,it,iu))*  &
                    frdm2(ip,iu,ir,it,rdm2act,occ,ind2,nact)  &
                    +twono(naddr3(is,iu,it,iq))*  &
                    frdm2(ip,iu,it,ir,rdm2act,occ,ind2,nact)

              END DO
            END DO
          END IF

          IF(ip == ir) arspq=arspq+occ(ip)*auxio(iqs)

        END IF
      END IF

!     T3+T4

      IF(occ(iq)*occ(is) /= zero) THEN

        IF(occ(iq) == one.OR.occ(is) == one) THEN

          IF(auxtwopqrs == one) arspq=arspq+occ(iq)*occ(is)*  &
              (two*twono(naddr3(ip,iq,ir,is))-twono(naddr3(ip,ir,iq,is)))
          IF(iq == is) arspq=arspq+occ(iq)*auxi(ipr)

        ELSE

          IF(igem(ip) == igem(ir)) THEN
            DO itt=1,nact
              DO iuu=1,nact
                it=ind1(itt)
                iu=ind1(iuu)

                IF(igfact(naddr3(iu,it,ip,ir)) == 1) arspq=arspq  &
                    +twono(naddr3(iu,it,ip,ir))*  &
                    frdm2(is,it,iq,iu,rdm2act,occ,ind2,nact)  &
                    +twono(naddr3(iu,ir,ip,it))*  &
                    frdm2(is,it,iu,iq,rdm2act,occ,ind2,nact)

              END DO
            END DO
          END IF

          IF(iq == is) arspq=arspq+occ(iq)*auxio(ipr)

        END IF

      END IF

!     T5

      IF(occ(ir)*occ(iq) /= zero) THEN

        IF(occ(iq) == one.OR.occ(ir) == one) THEN

          arspq=arspq-occ(iq)*occ(ir)*auxpqrs

        ELSE

          IF(igem(ip) == igem(is)) THEN
            DO itt=1,nact
              DO iuu=1,nact
                it=ind1(itt)
                iu=ind1(iuu)

                IF(igfact(naddr3(ip,it,is,iu)) == 1) arspq=arspq  &
                    -twono(naddr3(ip,it,is,iu))*  &
                    frdm2(it,iu,iq,ir,rdm2act,occ,ind2,nact)

              END DO
            END DO
          END IF

        END IF
      END IF

!     T6

      IF(occ(ip)*occ(is) /= zero) THEN

        IF(occ(ip) == one.OR.occ(is) == one) THEN

          arspq=arspq-occ(ip)*occ(is)*auxpqrs

        ELSE

          IF(igem(ir) == igem(iq)) THEN
            DO itt=1,nact
              DO iuu=1,nact
                it=ind1(itt)
                iu=ind1(iuu)

                IF(igfact(naddr3(it,iq,iu,ir)) == 1) arspq=arspq  &
                    -twono(naddr3(it,iq,iu,ir))*  &
                    frdm2(is,ip,iu,it,rdm2act,occ,ind2,nact)

              END DO
            END DO
          END IF

        END IF
      END IF

      IF(is == iq) arspq=arspq-half*wmat(ip,ir)
      IF(ip == ir) arspq=arspq-half*wmat(iq,is)

      abpl=abpl+arspq
      IF(ip > iq) THEN
        abmin=abmin+arspq
      ELSE
        abmin=abmin-arspq
      END IF

!     If(IP.Ne.IQ)
    END IF
!     end of IP,IQ LOOPS
  END DO
END DO

ip=ipp
iq=iqq

IF((c(ip)+c(iq))*(c(ir)+c(is)) /= zero) abpl=abpl/(c(ip)+c(iq))/(c(ir)+c(is))
IF((c(ip)-c(iq))*(c(ir)-c(is)) /= zero)  &
    abmin=abmin/(c(ip)-c(iq))/(c(ir)-c(is))

RETURN
END SUBROUTINE ab0element


REAL*8 FUNCTION frdm2(ip,iq,ir,is,rdm2act,occ,ind2,nact)

!     FOR A GIVEN SET OF INDICES AND THE KNOWN PART OF ACTIVE RDM2
!     RETURNS THE ELEMENT OF RDM2_PQRS FOR CAS

INTEGER, INTENT(IN )                  :: ip,iq,ir, is
REAL(8), INTENT(IN), dimension(:)           :: rdm2act,occ
integer, intent(in), dimension(:) :: ind2
INTEGER, INTENT(IN )                  ::  nact
REAL(8) :: rdm2
REAL, PARAMETER :: zero=0.d0
REAL, PARAMETER :: half=0.5D0
REAL, PARAMETER :: one=1.d0
REAL, PARAMETER :: two=2.d0
REAL, PARAMETER :: four=4.d0

rdm2=zero
IF(ip == ir.AND.iq == is.AND. (occ(ip) == one.OR.occ(iq) == one) )  &
    rdm2=rdm2+two*occ(ip)*occ(iq)
IF(ip == is.AND.iq == ir.AND. (occ(ip) == one.OR.occ(iq) == one) )  &
    rdm2=rdm2-occ(ip)*occ(iq)

!     ACTIVE PART

IF(ind2(ip)*ind2(iq)*ind2(ir)*ind2(is) /= zero) THEN
  rdm2=rdm2act(naddrrdm(ind2(ip),ind2(iq),ind2(ir),ind2(is),nact))
END IF

frdm2=rdm2

RETURN
END FUNCTION frdm2

INTEGER FUNCTION naddrrdm(ind1,ind2,ind3,ind4,nbasis)

!     A POINTER FOR 2-RDM ELEMENT

INTEGER, INTENT(IN)                  :: ind1,ind2,ind3,ind4,nbasis
integer :: ind12, ind34
ind12=(ind1-1)*nbasis+ind2
ind34=(ind3-1)*nbasis+ind4
naddrrdm=MAX(ind12,ind34)*(MAX(ind12,ind34)-1)/2+ MIN(ind12,ind34)
RETURN
END FUNCTION naddrrdm

 !NDim=NBasis*(NBasis-1)/2



SUBROUTINE erpasymm0(eigy,eigx,eig,aplsqrt,abmin,ndimx)

!     ALMOST THE SAME AS ERPASYMM1 BUT BOTH Y AND X ARE RETURNED

!     A SYMMETRIZED PROBLEM A+^(1/2) A- A+^(1/2) [A+^(-1/2)] Y = om^2 [A+^(-1/2)] Y IS SOLVED

!     ABPLUS IS CHANGED AND TURNS INTO ABPLUS^(1/2)
!     THIS ALLOWS ONE TO GET RID OF ONE LOCAL BIG ARRAY (COMPARING WITH ERPASYMM)

INTEGER, INTENT(IN)  :: ndimx
real(8), dimension(ndimx,ndimx), INTENT(INOUT) :: aplsqrt, abmin
real(8), INTENT(INOUT), dimension(ndimx*ndimx)  :: eigx, eigy
real(8), dimension(ndimx),intent(out) :: eig

real(8), dimension(ndimx,ndimx) :: hlpab

REAL(8), PARAMETER :: zero=0.d0
REAL(8), PARAMETER :: half=0.5D0
REAL(8), PARAMETER :: one=1.d0
REAL(8), PARAMETER :: two=2.d0
real(8), PARAMETER :: three=3.d0
real(8), PARAMETER :: four=4.d0
real(8), PARAMETER :: small=1.d-6

integer :: i,j,k, nu, noneg
real(8) :: sqrteig,sumnu



!     SYMMETRIZE A+,A-

DO i=1,ndimx
  DO j=i+1,ndimx
    aplsqrt(i,j)=half*(aplsqrt(i,j)+aplsqrt(j,i))
    aplsqrt(j,i)=aplsqrt(i,j)
    abmin(i,j)=half*(abmin(i,j)+abmin(j,i))
    abmin(j,i)=abmin(i,j)
  END DO
END DO

!     FIND A+^(1/2)

noneg=0
hlpab=aplsqrt
! CALL cpym(hlpab,aplsqrt,ndimx)
CALL diag8(hlpab,ndimx,ndimx,eig)

DO i=1,ndimx

  IF(eig(i) < zero) noneg=noneg+1

  DO j=i,ndimx

    aplsqrt(i,j)=zero

    DO k=1,ndimx

      sqrteig=SQRT(ABS(eig(k)))

      aplsqrt(i,j)=aplsqrt(i,j)+hlpab(k,i)*sqrteig*hlpab(k,j)

    END DO

    aplsqrt(j,i)=aplsqrt(i,j)

  END DO
END DO

IF(noneg /= 0) THEN
  WRITE(6,*)"The ERPA A+ matrix is not nonnegative definite"
  WRITE(6,*)"The number of negative eigenvalues of A+ is",noneg
END IF



CALL dgemm('N','N',ndimx,ndimx,ndimx,1D0,abmin,ndimx,  &
    aplsqrt,ndimx,0D0,hlpab,ndimx)
CALL dgemm('N','N',ndimx,ndimx,ndimx,1D0,aplsqrt,ndimx,  &
    hlpab,ndimx,0D0,eigy,ndimx)

CALL diag8(eigy,ndimx,ndimx,eig)

!     COMPUTE Y's

DO i=1,ndimx
  DO j=i,ndimx
    hlpab(i,j)=eigy(ndimx*(i-1)+j)
    hlpab(j,i)=eigy(ndimx*(j-1)+i)
  END DO
END DO



CALL dgemm('N','N',ndimx,ndimx,ndimx,1D0,aplsqrt,ndimx,  &
    hlpab,ndimx,0D0,eigy,ndimx)

!     IMPOSE THE NORMALIZATION 2 Y*X = 1 ON THE EIGENVECTORS CORRESPONDING TO POSITIVE
!     OMEGA'S

!     SINCE X = Om^-1 ABMIN.Y THEN THE NORMALIZATION READS 2 Om^-1 Y^T AMIN Y = 1

DO nu=1,ndimx
  sumnu=zero

  IF(eig(nu) > small) THEN

    eig(nu)=SQRT(eig(nu))

    DO i=1,ndimx
      DO j=1,ndimx
        sumnu=sumnu+two/eig(nu)*abmin(i,j)*  &
            eigy((nu-1)*ndimx+i)*eigy((nu-1)*ndimx+j)
      END DO
    END DO

    IF(sumnu > zero) THEN
      sumnu=one/SQRT(sumnu)
    ELSE
      sumnu=zero
    END IF

    DO i=1,ndimx
      eigy((nu-1)*ndimx+i)=eigy((nu-1)*ndimx+i)*sumnu
    END DO

  END IF
!     enddo NU
END DO

!     COMPUTE EigX

DO nu=1,ndimx

  IF(eig(nu) > small) THEN

    DO i=1,ndimx
      eigx((nu-1)*ndimx+i)=zero
      DO j=1,ndimx
        eigx((nu-1)*ndimx+i)=eigx((nu-1)*ndimx+i)  &
            +one/eig(nu)*abmin(i,j)*eigy((nu-1)*ndimx+j)
      END DO
    END DO

  ELSE

    DO i=1,ndimx
      eigy((nu-1)*ndimx+i)=zero
      eigx((nu-1)*ndimx+i)=zero
    END DO

  END IF
!     enddo NU
END DO

RETURN
END SUBROUTINE erpasymm0



SUBROUTINE ab1_cas(abplus,abmin,ure,occ,xone,twono,  &
    rdm2act,igfact,c,ind1,ind2,  &
    indblock,noeig,ndimx,nbasis,ninte1,ninte2, naccas, ninaccas, igem)
!added naccas, ninaccas, igem

!     COMPUTE THE ALPHA-DEPENDENT PARTS OF A+B AND A-B MATRICES

!     RDM2 IS IN NO REPRESENTATION. IT IS PARTIALLY SPIN-SUMMED
!     THE FOLLOWING SYMMETRY IS ASSUMED
!     RDM2(ij,kl) = RDM2(kl,ij)
!     ONLY ELEMENTS ij >= kl ARE STORED BUT ONLY FOR THE ACTIVE ELEMENTS
!     SIZE: NRDM2 = NBasis**2*(NBasis**2+1)/2
implicit none
INTEGER, INTENT(IN )   :: ndimx, noeig
integer, dimension(nbasis), intent(in) :: igem
real(8), dimension((NBasis*(NBasis-1)/2)*(NBasis*(NBasis-1)/2)), intent(inout) :: abplus, abmin
INTEGER, INTENT(IN ) :: nbasis,ninte1, ninte2, naccas, ninaccas
Integer, dimension(nbasis), INTENT(IN ) :: ind1,ind2
real(8), INTENT(IN ), dimension(ninte2) ::  twono
integer, dimension(2,ndimx), intent(in) :: indblock
real(8), INTENT(IN ) , dimension(nbasis) :: occ
real(8), INTENT(IN ), dimension(nbasis, nbasis) :: ure
real(8), INTENT(IN ), dimension(ninte1)   :: xone
real(8), dimension(:), intent(in) :: RDM2Act
integer, dimension(ninte2), intent(in) :: igfact
Real(8), dimension(nbasis), intent(inout) :: c

REAL(8), PARAMETER :: zero=0.d0
REAL(8), PARAMETER :: half=0.5D0
REAL(8), PARAMETER :: one=1.d0
REAL(8), PARAMETER :: two=2.d0
real(8), PARAMETER :: three=3.d0
real(8), PARAMETER :: four=4.d0

real(8), dimension(ninte1) :: hno,auxi,auxio
real(8), dimension(nbasis,nbasis) :: wmat
real(8) :: arspq, auxpqrs, auxtwopqrs, aux
integer :: i, ia, iab, icol, ib, icolend, icount, ifunsr, iqs, iw, j, nact, noccup
integer :: ij, inactive, ip, ipp,ipq, ipr, iq, iqq, ir, is, it, itt, irow, iu, iuu
real :: start_time, end_time
ifunsr=0 ! no short-range functional


DO i=1,noeig*noeig
  abplus(i)=zero
  abmin(i)=zero
END DO

!     ONE-ELECTRON MATRIX IN A NO REPRESENTATION

ij=0
DO i=1,nbasis
  DO j=1,i
    ij=ij+1
    hno(ij)=zero

    DO ia=1,nbasis
      DO ib=1,nbasis
        iab=(MAX(ia,ib)*(MAX(ia,ib)-1))/2+MIN(ia,ib)
        hno(ij)=hno(ij)+ure(i,ia)*ure(j,ib)*xone(iab)
      END DO
    END DO

  END DO
END DO

!     CONSTRUCT ONE-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN

ij=0
DO i=1,nbasis
  DO j=1,i
    ij=ij+1

    IF(igem(i) == igem(j)) THEN

      aux=zero

      DO it=1,nbasis
        IF(igem(it) /= igem(i)) aux=aux+occ(it)*  &
            (two*twono(naddr3(it,it,i,j))-twono(naddr3(it,i,it,j)))
      END DO

      hno(ij)=-aux

    END IF

  END DO
END DO

nact=naccas
inactive=ninaccas
noccup=inactive+nact

!     AUXILIARY MATRIX AuxI

ipq=0
DO ip=1,nbasis
  DO iq=1,ip
    ipq=ipq+1
    auxi(ipq)=zero
    auxio(ipq)=zero
    DO it=1,noccup
      IF(igfact(naddr3(it,it,ip,iq)) == 0) THEN
        auxi(ipq)=auxi(ipq)+occ(it)*  &
            (two*twono(naddr3(ip,iq,it,it))-twono(naddr3(ip,it,iq,it)))
        IF(it <= inactive) auxio(ipq)=auxio(ipq)+occ(it)*  &
            (two*twono(naddr3(ip,iq,it,it))-twono(naddr3(ip,it,iq,it)))
      END IF
    END DO
  END DO
END DO

!     AUXILIARY MATRIX WMAT

DO i=1,nbasis
  DO j=1,nbasis
    wmat(i,j)=zero
  END DO
END DO

DO ip=1,nbasis
  DO it=1,noccup
    DO iw=1,noccup
      DO iu=1,noccup

        IF(igfact(naddr3(it,iw,ip,iu)) == 0) THEN

          DO ir=1,noccup
            wmat(ip,ir)=wmat(ip,ir) +twono(naddr3(it,iw,ip,iu))  &
                *frdm2(iw,iu,it,ir,rdm2act,occ,ind2,nact)  &
                +twono(naddr3(it,iu,ip,iw))  &
                *frdm2(iw,iu,ir,it,rdm2act,occ,ind2,nact)
          END DO

        END IF

      END DO
    END DO
  END DO
END DO

icount=0
CALL cpu_time(start_time)

!     HAP-TEST!
!      IFunSR=4

DO irow=1,noeig

  ir=indblock(1,irow)
  is=indblock(2,irow)

  IF(ifunsr == 4) THEN
    icolend=noeig
  ELSE
    icolend=irow
  END IF

  DO icol=1,icolend

    ipp=indblock(1,icol)
    iqq=indblock(2,icol)

    IF( .NOT.(igem(ir) == igem(is).AND.igem(ir) == igem(ipp)  &
          .AND.igem(ir) == igem(iqq)) ) THEN
!C
      IF( (occ(ir)*occ(is) == zero.AND.occ(ipp)*occ(iqq) == zero  &
            .AND.ABS(twono(naddr3(ir,is,ipp,iqq))) < 1.d-25) .OR.  &
            ((occ(ir) == one.OR.occ(is) == one) .AND.  &
            (occ(ipp) == one.OR.occ(iqq) == one)  &
            .AND.ABS(twono(naddr3(ir,is,ipp,iqq))) < 1.d-25)) THEN

        icount=icount+1

      ELSE

        DO ip=iqq,ipp,ipp-iqq
          DO iq=iqq,ipp,ipp-iqq
            IF(ip /= iq) THEN

              iqs=(MAX(is,iq)*(MAX(is,iq)-1))/2+MIN(is,iq)
              ipr=(MAX(ir,ip)*(MAX(ir,ip)-1))/2+MIN(ir,ip)

              arspq=zero

              IF(ip == ir) arspq=arspq+(occ(ip)-occ(is))*hno(iqs)
              IF(is == iq) arspq=arspq+(occ(iq)-occ(ir))*hno(ipr)

              auxtwopqrs=zero
              IF(igfact(naddr3(ip,iq,ir,is)) == 0) auxtwopqrs=one
              auxpqrs=auxtwopqrs*  &
                  (two*twono(naddr3(ip,iq,ir,is))-twono(naddr3(ip,ir,iq,is)))

!     T1+T2

              IF(occ(ip)*occ(ir) /= zero) THEN

                IF(occ(ip) == one.OR.occ(ir) == one) THEN

                  IF(auxtwopqrs == one) arspq=arspq+occ(ip)*occ(ir)*  &
                      (two*twono(naddr3(ip,iq,ir,is))-twono(naddr3(ip,ir,iq,is)))
                  IF(ip == ir) arspq=arspq+occ(ip)*auxi(iqs)

                ELSE

                  DO itt=1,nact
                    DO iuu=1,nact
                      it=ind1(itt)
                      iu=ind1(iuu)

                      IF(igfact(naddr3(is,iq,it,iu)) == 0) arspq=arspq  &
                          +twono(naddr3(is,iq,it,iu))*  &
                          frdm2(ip,iu,ir,it,rdm2act,occ,ind2,nact)  &
                          +twono(naddr3(is,iu,it,iq))*  &
                          frdm2(ip,iu,it,ir,rdm2act,occ,ind2,nact)

                    END DO
                  END DO

                  IF(ip == ir) arspq=arspq+occ(ip)*auxio(iqs)

                END IF
              END IF

!     T3+T4

              IF(occ(iq)*occ(is) /= zero) THEN

                IF(occ(iq) == one.OR.occ(is) == one) THEN

                  IF(auxtwopqrs == one) arspq=arspq+occ(iq)*occ(is)*  &
                      (two*twono(naddr3(ip,iq,ir,is))-twono(naddr3(ip,ir,iq,is)))
                  IF(iq == is) arspq=arspq+occ(iq)*auxi(ipr)

                ELSE

                  DO itt=1,nact
                    DO iuu=1,nact
                      it=ind1(itt)
                      iu=ind1(iuu)

                      IF(igfact(naddr3(iu,it,ip,ir)) == 0) arspq=arspq  &
                          +twono(naddr3(iu,it,ip,ir))*  &
                          frdm2(is,it,iq,iu,rdm2act,occ,ind2,nact)  &
                          +twono(naddr3(iu,ir,ip,it))*  &
                          frdm2(is,it,iu,iq,rdm2act,occ,ind2,nact)

                    END DO
                  END DO

                  IF(iq == is) arspq=arspq+occ(iq)*auxio(ipr)

                END IF

              END IF

!     T5

              IF(occ(ir)*occ(iq) /= zero) THEN
                IF(occ(iq) == one.OR.occ(ir) == one) THEN

                  arspq=arspq-occ(iq)*occ(ir)*auxpqrs

                ELSE

                  DO itt=1,nact
                    DO iuu=1,nact
                      it=ind1(itt)
                      iu=ind1(iuu)

                      IF(igfact(naddr3(ip,it,is,iu)) == 0) arspq=arspq  &
                          -twono(naddr3(ip,it,is,iu))*  &
                          frdm2(it,iu,iq,ir,rdm2act,occ,ind2,nact)

                    END DO
                  END DO

                END IF
              END IF

!     T6

              IF(occ(ip)*occ(is) /= zero) THEN
                IF(occ(ip) == one.OR.occ(is) == one) THEN

                  arspq=arspq-occ(ip)*occ(is)*auxpqrs

                ELSE

                  DO itt=1,nact
                    DO iuu=1,nact
                      it=ind1(itt)
                      iu=ind1(iuu)

                      IF(igfact(naddr3(it,iq,iu,ir)) == 0) arspq=arspq  &
                          -twono(naddr3(it,iq,iu,ir))*  &
                          frdm2(is,ip,iu,it,rdm2act,occ,ind2,nact)

                    END DO
                  END DO

                END IF
              END IF

              IF(is == iq) arspq=arspq-half*wmat(ip,ir)
              IF(ip == ir) arspq=arspq-half*wmat(iq,is)

              IF(ir > is.AND.ip > iq) THEN
                abplus(irow+(icol-1)*noeig)=abplus(irow+(icol-1)*noeig)+arspq
                abmin(irow+(icol-1)*noeig)=abmin(irow+(icol-1)*noeig)+arspq
              END IF

              IF(ir > is.AND.iq > ip) THEN
                abplus(irow+(icol-1)*noeig)=abplus(irow+(icol-1)*noeig)+arspq
                abmin(irow+(icol-1)*noeig)=abmin(irow+(icol-1)*noeig)-arspq
              END IF

!     If(IP.Ne.IQ)
            END IF

!     end of IP,IQ LOOPS
          END DO
        END DO
        
!     end icount
      END IF
!     If IGem ....
    END IF

  END DO
END DO
CALL cpu_time(end_time)
!WRITE(6,'(X,"TIME SPENT ON CONSTRUCTING AB(1)" ,F10.2)')end_time-start_time
!WRITE(*,*)'icount',icount

IF(ifunsr == 4) THEN

!     POSTCAS: DIVIDE BY C'c AND SYMMETRIZE

  DO i=1,noeig
    ip=indblock(1,i)
    iq=indblock(2,i)
    DO j=1,noeig
      ir=indblock(1,j)
      is=indblock(2,j)

      IF((c(ip)+c(iq))*(c(ir)+c(is)) /= zero)  &
          abplus(i+(j-1)*noeig)=abplus(i+(j-1)*noeig) /(c(ip)+c(iq))/(c(ir)+c(is))
      IF((c(ip)-c(iq))*(c(ir)-c(is)) /= zero)  &
          abmin(i+(j-1)*noeig)=abmin(i+(j-1)*noeig) /(c(ip)-c(iq))/(c(ir)-c(is))

    END DO
  END DO

  DO i=1,noeig
    DO j=i+1,noeig
      abplus((j-1)*noeig+i)=  &
          half*(abplus((j-1)*noeig+i)+abplus((i-1)*noeig+j))
      abplus((i-1)*noeig+j)=abplus((j-1)*noeig+i)
      abmin((j-1)*noeig+i)= half*(abmin((j-1)*noeig+i)+abmin((i-1)*noeig+j))
      abmin((i-1)*noeig+j)=abmin((j-1)*noeig+i)
    END DO
  END DO

ELSE
!     DIVIDE BY C'c AND COPY TRIANGLE

  DO i=1,noeig
    ip=indblock(1,i)
    iq=indblock(2,i)
    DO j=1,i
      ir=indblock(1,j)
      is=indblock(2,j)

      IF((c(ip)+c(iq))*(c(ir)+c(is)) /= zero)  &
          abplus(i+(j-1)*noeig)=abplus(i+(j-1)*noeig) /(c(ip)+c(iq))/(c(ir)+c(is))
      IF((c(ip)-c(iq))*(c(ir)-c(is)) /= zero)  &
          abmin(i+(j-1)*noeig)=abmin(i+(j-1)*noeig) /(c(ip)-c(iq))/(c(ir)-c(is))

      abplus(j+(i-1)*noeig)=abplus(i+(j-1)*noeig)
      abmin(j+(i-1)*noeig)=abmin(i+(j-1)*noeig)

    END DO
  END DO

!     end IFunSR
END IF

RETURN
END SUBROUTINE ab1_cas


SUBROUTINE ac0cas(ecorr,etot,twono,occ,ure,xone,  &
    indn,ndimx,nbasis,ninte1,ninte2,rdm2act, naccas, n_electron, n_act_ele)
real(8), INTENT(OUT)  :: ecorr, etot
real(8), INTENT(IN ), dimension(ninte2) ::  twono
real(8), INTENT(IN ) , dimension(nbasis) :: occ
real(8), INTENT(IN ), dimension(nbasis, nbasis) :: ure
real(8), INTENT(IN ), dimension(ninte1)   :: xone


INTEGER, INTENT(IN )                  :: nbasis,ninte1, ninte2
INTEGER, INTENT(IN ), dimension(2,NBasis*(NBasis-1)/2) :: indn
INTEGER, INTENT(IN )                  :: ndimx
real(8), dimension(:), intent(in) :: RDM2Act
integer, intent(in) :: naccas, n_act_ele, n_electron

Real(8), dimension(nbasis) :: c
Integer, dimension(nbasis) :: ind1,ind2
real(8), dimension(nbasis,nbasis) :: wmat
real(8), dimension(ninte1) :: hno, auxi, auxio
integer,dimension(ninte2) :: igfact
integer, dimension(nbasis,nbasis) :: ipair
real(8), dimension(ndimx*ndimx) :: eigx, xmaux
integer, dimension(2,ndimx) :: ieigaddy, ieigaddind, indblock
integer :: nrdm2act, nact,inactive, nele, noccup, ngem
integer :: i, ia, iab, ii, ij, ib, ip, iq, ir, is, j, it, k, l, kl, nadd
integer :: ipp, ipq, iqq, iu, iw, ndimb, nfree1, nfree2, noeig, icol, irow, irr, iss
real(8) :: aux, abp, abm
real(8), dimension(ndimx*ndimx)  :: abplus, abmin, eigy !  Cono
real(8), dimension(nbasis) :: cicoef !needed in several places downstream, the dimension should probably be nbasis(!) rather than 1000
integer, dimension(nbasis) :: igem
integer :: nost, istart, mu, nu
real(8) :: eall, eintra, sumy

real(8), dimension(ndimx) :: eig ! Cono
REAL(8), PARAMETER :: zero=0.d0
REAL(8), PARAMETER :: half=0.5D0
REAL(8), PARAMETER :: one=1.d0
REAL(8), PARAMETER :: two=2.d0
real(8), PARAMETER :: three=3.d0
real(8), PARAMETER :: four=4.d0

nost=1


nact=naccas !number of active orbitals (argument)
!n_act_ele !number of active electrons is an argument now!!

NELE=Int(n_electron)/2 !def copied from the original code
!print *, 'NELE=',NELE ! NELE  - half of the number of the electrons
!print *, 'number of active orbitals nact=',nact
!print *, 'number of active electrons =',n_act_ele
INActive=NELE-(n_act_ele/2)
!print *, 'number of inactive electron pairs =',INActive

noccup=inactive+nact


If(INActive == 0) Then
  NGem=2
  IGem(1:noccup)=1
  IGem(noccup+1:NBasis)=2
Else
  NGem=3
  IGem(1:INActive)=1
  IGem(INActive+1:noccup)=2
  IGem(noccup+1:NBasis)=3
EndIf
!print *, 'igem=',igem


nrdm2act = nact**2*(nact**2+1)/2



!     A ROUTINE FOR COMPUTING AC INTEGRAND

ipair(1:nbasis,1:nbasis)=0
DO ii=1,ndimx
  i=indn(1,ii)
  j=indn(2,ii)
  ipair(i,j)=1
  ipair(j,i)=1
END DO



!     AUXILIARY STUFF LATER NEEDED TO GET A+ AND A- MATRICES FOR ALPHA=0

!     ONE-ELECTRON MATRIX IN A NO REPRESENTATION

! sets hno=xone, because ure=1
ij=0
DO i=1,nbasis
  DO j=1,i
    ij=ij+1
    hno(ij)=zero

    DO ia=1,nbasis
      DO ib=1,nbasis
        iab=(MAX(ia,ib)*(MAX(ia,ib)-1))/2+MIN(ia,ib)
        hno(ij)=hno(ij)+ure(i,ia)*ure(j,ib)*xone(iab)
      END DO
    END DO

  END DO
END DO


!     COMPUTE THE ENERGY

!print *, 'inactive=',inactive


ind2(1:nbasis)=0
DO i=1,nact
  ind1(i)=inactive+i
  ind2(inactive+i)=i
END DO



!     COMPUTE THE ENERGY FOR CHECKING

etot=zero
DO i=1,nbasis
  ii=(i*(i+1))/2
  etot=etot+two*occ(i)*hno(ii)
END DO

!print *, '1-electron etot=',etot

DO ip=1,noccup
  DO iq=1,noccup
    DO ir=1,noccup
      DO is=1,noccup
        !print *, 'p=',ip,' q=',iq,' r=',ir,' s=',is
        !print *, 'frdm2=',frdm2(ip,iq,ir,is,rdm2act,occ,ind2,nact,nbasis) ! seems correct
        etot=etot+frdm2(ip,iq,ir,is,rdm2act,occ,ind2,nact)  &
            *twono(naddr3(ip,ir,iq,is))
      END DO
    END DO
  END DO
END DO

!WRITE(6,'(/,X,''CASSCF ENERGY (W/O ENUC)'',5X,F15.8)')etot



DO i=1,nbasis
  c(i)=SQRT(occ(i))
  IF(occ(i) < half) c(i)=-c(i)
  cicoef(i)=c(i)
  ! this is needed in many places and it's a global variable
  ! defined in commons.inc as CICoef(1000)
  ! so it has to be passed as argument downstream!!
END DO

!     CONSTRUCT ONE-ELECTRON PART OF THE AC ALPHA=0 HAMILTONIAN

ij=0
DO i=1,nbasis
  DO j=1,i
    ij=ij+1
    !igem is defined in commons.inc so we need to pass it as variable downstream
    IF(igem(i) /= igem(j)) THEN

      hno(ij)=zero

    ELSE

      aux=zero

      DO it=1,nbasis
        IF(igem(it) /= igem(i)) aux=aux+occ(it)*  &
            (two*twono(naddr3(it,it,i,j))-twono(naddr3(it,i,it,j)))
      END DO

      hno(ij)=hno(ij)+aux

    END IF

  END DO
END DO

!     CONSTRUCT TWO-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN

nadd=0
ij=0
DO i=1,nbasis
  DO j=1,i
    ij=ij+1
    kl=0
    DO k=1,nbasis
      DO l=1,k
        kl=kl+1

        IF(ij >= kl) THEN
          nadd=nadd+1

          igfact(nadd)=1
          IF(.NOT.(  &
              igem(i) == igem(j).AND.igem(j) == igem(k).AND.igem(k) == igem(l)))  &
              igfact(nadd)=0

        END IF

      END DO
    END DO
  END DO
END DO


!     AUXILIARY MATRIX AuxI AND AuxIO

ipq=0
DO ip=1,nbasis
  DO iq=1,ip
    ipq=ipq+1
    auxi(ipq)=zero
    auxio(ipq)=zero
    DO it=1,noccup
      IF(igfact(naddr3(it,it,ip,iq)) == 1) THEN
        auxi(ipq)=auxi(ipq)+occ(it)*  &
            (two*twono(naddr3(ip,iq,it,it))-twono(naddr3(ip,it,iq,it)))
        IF(it <= inactive) auxio(ipq)=auxio(ipq)+occ(it)*  &
            (two*twono(naddr3(ip,iq,it,it))-twono(naddr3(ip,it,iq,it)))
      END IF
    END DO
  END DO
END DO

!     AUXILIARY MATRIX WMAT

DO i=1,nbasis
  DO j=1,nbasis
    wmat(i,j)=zero
  END DO
END DO

DO ip=1,nbasis
  DO ir=1,noccup
    DO it=1,noccup
      DO iw=1,noccup
        DO iu=1,noccup
          IF(igfact(naddr3(it,iw,ip,iu)) == 1) wmat(ip,ir)=wmat(ip,ir)  &
              +twono(naddr3(it,iw,ip,iu))  &
              *frdm2(iw,iu,it,ir,rdm2act,occ,ind2,nact)  &
              +twono(naddr3(it,iu,ip,iw))  &
              *frdm2(iw,iu,ir,it,rdm2act,occ,ind2,nact)

        END DO
      END DO
    END DO
  END DO
END DO
!     FIND THE 0TH-ORDER SOLUTION FOR THE ACTIVE-ACTIVE BLOCK

!WRITE(6,'(" *** ACTIVE-ACTIVE BLOCK ***")')

nfree1=1
nfree2=1
noeig=0

ndimb=0
DO iqq=1,nact
  DO ipp=iqq+1,nact
    ip=ind1(ipp)
    iq=ind1(iqq)
    IF(ipair(ip,iq) == 1) THEN
      ndimb=ndimb+1
      indblock(1,nfree1-1+ndimb)=ip
      indblock(2,nfree1-1+ndimb)=iq
    END IF
  END DO
END DO

DO i=1,ndimb
  ieigaddy(1,nfree1-1+i)=nfree2+(i-1)*ndimb
  ieigaddy(2,nfree1-1+i)=ieigaddy(1,nfree1-1+i)+ndimb-1
  ieigaddind(1,nfree1-1+i)=nfree1
  ieigaddind(2,nfree1-1+i)=nfree1+ndimb-1
END DO


irow=0
DO iqq=1,nact
  DO ipp=iqq+1,nact
    ip=ind1(ipp)
    iq=ind1(iqq)
    IF(ipair(ip,iq) == 1) THEN

      irow=irow+1

      icol=0
      DO iss=1,nact
        DO irr=iss+1,nact
          ir=ind1(irr)
          is=ind1(iss)
          IF(ipair(ir,is) == 1) THEN

            icol=icol+1

            IF(irow >= icol) THEN

              CALL ab0element(abp,abm,ip,iq,ir,is,occ,hno,igfact,  &
                  twono,auxi,auxio,wmat,rdm2act,c,ind1,ind2,nact,nrdm2act,  &
                  ninte1,ninte2,nbasis,igem)

              abplus((icol-1)*ndimb+irow)=abp
              abplus((irow-1)*ndimb+icol)=abp
              abmin((icol-1)*ndimb+irow)=abm
              abmin((irow-1)*ndimb+icol)=abm
              !print *, 'abm=',abm

            END IF

          END IF
        END DO
      END DO

    END IF
  END DO
END DO
!print *, abmin
IF(ndimb /= 0) THEN
      !Print *, 'ACT-KA',norm2(ABPLUS(1:NDimB**2)),norm2(ABMIN)
  !IF(nost == 1) THEN
    CALL erpasymm0(eigy(nfree2),eigx(nfree2),eig(nfree1),abplus,abmin, ndimb)
  !ELSE
  !  CALL erpavecyx(eigy(nfree2),eigx(nfree2),eig(nfree1),abplus,abmin, ndimb)
  !END IF
ELSE
  !print *, 'ndimb == 0, whatever this means'
END IF

noeig=noeig+ndimb
nfree1=noeig+1
nfree2=nfree2+ndimb*ndimb

!     FIND THE 0TH-ORDER SOLUTION FOR THE ACTIVE-INACTIVE BLOCKS

!WRITE(6,'(" *** ACTIVE-INACTIVE BLOCKS ***")')

DO iq=1,inactive

  ndimb=0
  DO ipp=1,nact
    ip=ind1(ipp)
    IF(ipair(ip,iq) == 1) THEN
      ndimb=ndimb+1
      indblock(1,nfree1-1+ndimb)=ip
      indblock(2,nfree1-1+ndimb)=iq
    END IF
  END DO

  DO i=1,ndimb
    ieigaddy(1,nfree1-1+i)=nfree2+(i-1)*ndimb
    ieigaddy(2,nfree1-1+i)=ieigaddy(1,nfree1-1+i)+ndimb-1
    ieigaddind(1,nfree1-1+i)=nfree1
    ieigaddind(2,nfree1-1+i)=nfree1+ndimb-1
  END DO

  irow=0
  DO ipp=1,nact
    ip=ind1(ipp)

    IF(ipair(ip,iq) == 1) THEN

      irow=irow+1

      icol=0
      is=iq
      DO irr=1,nact
        ir=ind1(irr)

        IF(ipair(ir,is) == 1) THEN

          icol=icol+1

          IF(irow >= icol) THEN

            CALL ab0element(abp,abm,ip,iq,ir,is,occ,hno,igfact,  &
                twono,auxi,auxio,wmat,rdm2act,c,ind1,ind2,nact,nrdm2act,  &
                ninte1,ninte2,nbasis,igem)

            abplus((icol-1)*ndimb+irow)=abp
            abplus((irow-1)*ndimb+icol)=abp
            abmin((icol-1)*ndimb+irow)=abm
            abmin((irow-1)*ndimb+icol)=abm

          END IF

        END IF

      END DO

    END IF

  END DO

  IF(ndimb /= 0) THEN
     !Print*, 'AI-KA',norm2(ABPLUS(1:NDimB**2)),norm2(ABMIN)
    !IF(nost == 1) THEN
      CALL erpasymm0(eigy(nfree2),eigx(nfree2),eig(nfree1),abplus,abmin,  &
          ndimb)
    !ELSE
     ! CALL erpavecyx(eigy(nfree2),eigx(nfree2),eig(nfree1),abplus,abmin,  &
    !      ndimb)
    !END IF
  END IF

  noeig=noeig+ndimb
  nfree1=noeig+1
  nfree2=nfree2+ndimb*ndimb

!     Do IP
END DO

!     FIND THE 0TH-ORDER SOLUTION FOR THE VIRTUAL-ACTIVE BLOCKS

!WRITE(6,'(" *** VIRTUAL-ACTIVE BLOCKS ***")')

DO ip=noccup+1,nbasis

  ndimb=0
  DO iqq=1,nact
    iq=ind1(iqq)
    IF(ipair(ip,iq) == 1) THEN
      ndimb=ndimb+1
      indblock(1,nfree1-1+ndimb)=ip
      indblock(2,nfree1-1+ndimb)=iq
    END IF
  END DO

  DO i=1,ndimb
    ieigaddy(1,nfree1-1+i)=nfree2+(i-1)*ndimb
    ieigaddy(2,nfree1-1+i)=ieigaddy(1,nfree1-1+i)+ndimb-1
    ieigaddind(1,nfree1-1+i)=nfree1
    ieigaddind(2,nfree1-1+i)=nfree1+ndimb-1
  END DO

  irow=0
  DO iqq=1,nact
    iq=ind1(iqq)

    IF(ipair(ip,iq) == 1) THEN

      irow=irow+1

      icol=0
      ir=ip
      DO iss=1,nact
        is=ind1(iss)

        IF(ipair(ir,is) == 1) THEN

          icol=icol+1

          IF(irow >= icol) THEN

            CALL ab0element(abp,abm,ip,iq,ir,is,occ,hno,igfact,  &
                twono,auxi,auxio,wmat,rdm2act,c,ind1,ind2,nact,nrdm2act,  &
                ninte1,ninte2,nbasis,igem)

            abplus((icol-1)*ndimb+irow)=abp
            abplus((irow-1)*ndimb+icol)=abp
            abmin((icol-1)*ndimb+irow)=abm
            abmin((irow-1)*ndimb+icol)=abm

          END IF

        END IF

      END DO

    END IF

  END DO

  IF(ndimb /= 0) THEN
     !Print*, 'AV-KA',norm2(ABPLUS(1:NDimB**2)),norm2(ABMIN(1:NDimB**2))

    CALL erpasymm0(eigy(nfree2),eigx(nfree2),eig(nfree1),abplus,abmin,  &
          ndimb)
  END IF

  noeig=noeig+ndimb
  nfree1=noeig+1
  nfree2=nfree2+ndimb*ndimb

!     Do IP
END DO

DO ip=noccup+1,nbasis
  DO iq=1,inactive

    ndimb=0

    IF(ipair(ip,iq) == 1) THEN

      ndimb=1
      indblock(1,nfree1)=ip
      indblock(2,nfree1)=iq

      ieigaddy(1,nfree1)=nfree2
      ieigaddy(2,nfree1)=ieigaddy(1,nfree1)
      ieigaddind(1,nfree1)=nfree1
      ieigaddind(2,nfree1)=nfree1

      CALL ab0element(abp,abm,ip,iq,ip,iq,occ,hno,igfact,  &
          twono,auxi,auxio,wmat,rdm2act,c,ind1,ind2,nact,nrdm2act,  &
          ninte1,ninte2,nbasis,igem)

      eig(nfree1)=abp



      eigy(nfree2)=one/SQRT(two)
      eigx(nfree2)=one/SQRT(two)

      noeig=noeig+ndimb
      nfree1=noeig+1
      nfree2=nfree2+ndimb*ndimb

    END IF

  END DO
END DO

!WRITE(6,'(/," *** DONE WITH 0TH-ORDER IN AC0-CASSCF ***")')
!PRINT*, 'NoEig,NDimX',noeig,ndimx
!      Print*, 'Eig,Y,X',norm2(Eig(1:NoEig)),
!     $ norm2(EigY(1:NoEig**2)),norm2(EigX(1:NoEig**2))

!     DONE 0TH-ORDER CALCULATIONS

!WRITE(6,'(/, " *** COMPUTING ABPLUS(1) AND ABMIN(1) MATRICES ***"  &
!    )')
CALL ab1_cas(abplus,abmin,ure,occ,xone,twono,  &
    rdm2act,igfact,c,ind1,ind2,  &
    indblock,noeig,ndimx,nbasis,ninte1,ninte2,naccas, INActive, igem)
        ! added naccas, INActive, igem to arguments!!
!WRITE(6,'(/," *** DONE WITH COMPUTING AB(1) MATRICES ***")')

!PRINT*, 'AB1-Ka',norm2(abplus),norm2(abmin)
! ----------------------------------------------------------------
!     1ST-ORDER PART

DO nu=1,noeig
  DO mu=1,noeig

    xmaux(nu+(mu-1)*noeig)=zero

    istart=ieigaddy(1,mu)
    ii=0
    DO i=ieigaddind(1,mu),ieigaddind(2,mu)
      xmaux(nu+(mu-1)*noeig)=xmaux(nu+(mu-1)*noeig)  &
          +abplus(nu+(i-1)*noeig)*eigx(istart+ii)
      ii=ii+1
    END DO

  END DO
END DO

DO nu=1,noeig
  DO mu=1,noeig

    abplus(nu+(mu-1)*noeig)=zero

    istart=ieigaddy(1,nu)
    ii=0
    DO i=ieigaddind(1,nu),ieigaddind(2,nu)
      abplus(nu+(mu-1)*noeig)=abplus(nu+(mu-1)*noeig)  &
          +eigx(istart+ii)*xmaux(i+(mu-1)*noeig)
      ii=ii+1
    END DO

  END DO
END DO

DO nu=1,noeig
  DO mu=1,noeig

    xmaux(nu+(mu-1)*noeig)=zero

    istart=ieigaddy(1,mu)
    ii=0
    DO i=ieigaddind(1,mu),ieigaddind(2,mu)
      xmaux(nu+(mu-1)*noeig)=xmaux(nu+(mu-1)*noeig)  &
          +abmin(nu+(i-1)*noeig)*eigy(istart+ii)
      ii=ii+1
    END DO

  END DO
END DO

DO nu=1,noeig
  DO mu=1,noeig

    abmin(nu+(mu-1)*noeig)=zero

    istart=ieigaddy(1,nu)
    ii=0
    DO i=ieigaddind(1,nu),ieigaddind(2,nu)
      abmin(nu+(mu-1)*noeig)=abmin(nu+(mu-1)*noeig)  &
          +eigy(istart+ii)*xmaux(i+(mu-1)*noeig)
      ii=ii+1
    END DO

  END DO
END DO

!PRINT*, 'AB-KA',norm2(abplus),norm2(abmin)

xmaux(1:noeig*noeig)=zero

DO mu=1,noeig
! KP 03/04/2018

  IF(eig(mu) /= zero) THEN

    DO nu=1,noeig

! KP 03/04/2018

      IF(eig(nu) /= zero) THEN

        istart=ieigaddy(1,nu)
        ii=0
        DO i=ieigaddind(1,nu),ieigaddind(2,nu)

          xmaux(mu+(i-1)*noeig)=xmaux(mu+(i-1)*noeig)+two*  &
              (abplus(mu+(nu-1)*noeig)-abmin(mu+(nu-1)*noeig))/  &
              (eig(mu)+eig(nu))*eigy(istart+ii)

          ii=ii+1
        END DO

      END IF

    END DO

  END IF

END DO

!      Print*, 'FIRST',norm2(XMAux)

abplus(1:noeig*noeig)=zero

DO mu=1,noeig

  istart=ieigaddy(1,mu)
  ii=0
  DO i=ieigaddind(1,mu),ieigaddind(2,mu)

    DO j=1,noeig
      abplus(i+(j-1)*noeig)=abplus(i+(j-1)*noeig)+xmaux(mu+(j-1)*noeig)  &
          *eigy(istart+ii)
    END DO

    ii=ii+1
  END DO
END DO
!PRINT*, 'ABPLUS-Ka',norm2(abplus(1:noeig*noeig))

!     FINALLY THE ENERGY CORRECTION
!     TESTY Z sapt.f90 -- remove later!
!      Call check_loop(ABPLUS,Occ,IndN,IndBlock,
!     $ NAct,INActive,NDimX,NDimX,NBasis)

eall=zero
eintra=zero

DO i=1,noeig

  ip=indblock(1,i)
  ir=indblock(2,i)

  DO j=1,noeig

    iq=indblock(1,j)
    is=indblock(2,j)

    IF(ip > ir.AND.iq > is) THEN

      sumy=abplus(i+(j-1)*noeig)
      aux=(c(is)+c(iq))*(c(ip)+c(ir))*sumy

      eall=eall+aux*twono(naddr3(ip,ir,iq,is))

      IF(igem(ip) == igem(ir).AND.igem(iq) == igem(is) .AND. &
         igem(ip) == igem(iq)) eintra=eintra+aux*twono(naddr3(ip,ir,iq,is))

!     endinf of If(IP.Gt.IR.And.IQ.Gt.IS)
    END IF

  END DO
END DO

ecorr=eall-eintra
!PRINT*, 'EAll,EIntra',eall,eintra
!print *, ecorr

end subroutine ac0cas

function trrdm2(rdm2,ure,nbasis) result(out_rdm)
!     TRANSFORM RDM2 WITH URe
implicit none
INTEGER, INTENT(IN)                      :: nbasis
REAL(8), INTENT(IN), dimension(nbasis,nbasis,nbasis,nbasis) :: rdm2
REAL(8), INTENT(IN), dimension(nbasis,nbasis)    :: ure
real(8), dimension(nbasis,nbasis,nbasis,nbasis) :: out_rdm


integer :: i,ia,ib,ic,id,j,k,l

REAL, PARAMETER :: zero=0.d0
REAL, PARAMETER :: half=0.5D0
REAL, PARAMETER :: one=1.d0



!     LOCAL ARRAYS

real(8), DIMENSION(nbasis,nbasis,nbasis,nbasis) :: aux1, aux2

!WRITE(6,'(/,X,"FCI RDM2 TRANSFORMATION TO NO IN PROCESS...")')

!     COPY RDM2 TO Aux1
aux1=rdm2

!     FIRST INDEX

DO ia=1,nbasis
  DO j=1,nbasis
    DO k=1,nbasis
      DO l=1,nbasis

        aux2(ia,j,k,l)=zero
        DO i=1,nbasis
          aux2(ia,j,k,l)=aux2(ia,j,k,l)+ure(ia,i)*aux1(i,j,k,l)
        END DO

      END DO
    END DO
  END DO
END DO

!     SECOND INDEX

DO ia=1,nbasis
  DO ib=1,nbasis
    DO k=1,nbasis
      DO l=1,nbasis

        aux1(ia,ib,k,l)=zero
        DO j=1,nbasis
          aux1(ia,ib,k,l)=aux1(ia,ib,k,l)+ure(ib,j)*aux2(ia,j,k,l)
        END DO

      END DO
    END DO
  END DO
END DO

!     THIRD INDEX

DO ia=1,nbasis
  DO ib=1,nbasis
    DO ic=1,nbasis
      DO l=1,nbasis

        aux2(ia,ib,ic,l)=zero
        DO k=1,nbasis
          aux2(ia,ib,ic,l)=aux2(ia,ib,ic,l)+ure(ic,k)*aux1(ia,ib,k,l)
        END DO

      END DO
    END DO
  END DO
END DO

!     FOURTH INDEX

DO ia=1,nbasis
  DO ib=1,nbasis
    DO ic=1,nbasis
      DO id=1,nbasis

        out_rdm(ia,ib,ic,id)=zero
        DO l=1,nbasis
          out_rdm(ia,ib,ic,id)=  &
              out_rdm(ia,ib,ic,id)+ure(id,l)*aux2(ia,ib,ic,l)
        END DO

      END DO
    END DO
  END DO
END DO

!WRITE(6,'(/,X,"DONE WITH FCI RDM2 TRANSFORMATION")')

RETURN
END function trrdm2



end module accas_lib