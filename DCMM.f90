!----------------------------------------------------------------------
! Simple Fortran 90 DCM code for many players.
! This is a stochastic Montecarlo iterative procedure
! Prizes are integers (but can be changed to real numbers)
! Code inspired by  http://arxiv.org/abs/2102.07738  article
! By E. Besalu 2021
! No exhaustive checks made here (in order to maintain a clear code)
! http://iqcc.udg.edu/~emili/Poker
!----------------------------------------------------------------------
program Test_DCM_Montecarlo_Code
implicit none

! Tournament Variables
integer N ! Total number of players in contest
double precision, allocatable :: DCM(:) ! Cummulated DCM values
integer, allocatable :: S(:) ! Players' stacks
integer, allocatable :: P(:) ! Prizes to be shared (some can be zero)

! Montecarlo items
integer Iterations ! Number of iterations to do

! Auxiliar
integer i

! Enters number of players and prizes
write(*,*) "Enter the number of players (e.g. 3):"
read(*,*) N

allocate(S(N),P(N))
write(*,*) "Enter the stacks (e.g. 1000,800,500):"
read(*,*) S(1:N)
write(*,*) "Enter the prizes (complete with zeroes, e.g. 500,200,0):"
read(*,*) P(1:N)
write(*,*) "For the proposed example the exact results are 310.681696 257.845768 131.472535"
write(*,*)
write(*,*) "Enter the number of iterations to do (e.g. 1000000):"
read(*,*) Iterations
write(*,*)

allocate(DCM(N))     ! DCM values
DCM(:)=0.0d0
call DCM_Iterative_Montecarlo(Iterations,N,S,P,DCM)

write(*,'("   Prizes:",99(i8,1x))') (P(i),i=1,N)
write(*,'("   Player:",99(1x,"--",i2," ---"))') (i,i=1,N)
write(*,'("   Stacks:",99(i8,1x))') (S(i),i=1,N)
write(*,'("      DCM:",99f9.3)')  DCM(:)
END
!----------------------------------------------------------------------

!----------------------------------------------------------------------
! DCM (iterative,Montecarlo) simple calling code. E.Besalu. Feb.2021
!----------------------------------------------------------------------
subroutine DCM_Iterative_Montecarlo(Iterations,n,So,Po,DCM)
implicit none

! Parameter
integer MaxDeep ! Maximal tree deep search to reach (prune tree there)
parameter (MaxDeep=50)

! Input variables
integer Iterations ! Number of iterations to do
integer n ! Number of players
integer So(n) ! Stacks
integer Po(n) ! Prizes (as many as players, insert zeroes if necessary)

! Output result
double precision DCM(n) ! Expected prizes distribution

! Working info
integer wo(n) ! Identities
integer w(n),ww(n) ! Dynamic identities
integer S(n),SS(n) ! Dynamic stacks (pre-showdown & post-showdown)
integer P(n),Ps(n) ! Dynamic prizes
integer Iter ! Iteration number
integer winner ! Who wins the hand inside a node (randomly chosen)
integer nr ! Number of players remaining in table
integer e ! Efective stack
integer L ! Number of players in a list (of losers or survival ones)
integer Deep ! Deep in path

! Auxiliar
integer i,j,k
double precision x

CALL RANDOM_SEED() ! PLEASE, USE YOUR PREFERRED RANDOM SEED GENERATOR!

! The prizes are to be sorted from lower to higher
Ps(:)=Po(:) ! In order to avoid outside collateral effects
do i=1,n-1
   do j=i+1,n
      if (Ps(i)>Ps(j)) then
         k=Ps(i); Ps(i)=Ps(j); Ps(j)=k
      end if
   end do
end do

DCM(:)=0.0d0
do i=1,n
   wo(i)=i
end do

Iter=0
DO WHILE (Iter<Iterations)

   Iter=Iter+1 ! Path number

   ! Initiation of a new iteration or path
   S(:)=So(:) ! The pre-showdown stacks
   w(:)=wo(:) !
   P(:)=Ps(:)
   nr=n
   Deep=0
   do while (Deep<MaxDeep)
      
      Deep=Deep+1

      ! Stacks. Entering stacks in node are S(:) & w(:)
      SS(1:nr)=S(1:nr) ! Will be the post-showdown stacks
      ww(1:nr)=w(1:nr) !

      ! The winner is chosen randomly
      call random_number(x) ! Please, use your own or effective
      winner=1+INT(nr*x)    ! random number generator

      ! New stacks: post-showdown stacks
      do i=1,nr
         if (i/=winner) then
            e=MIN(S(i),S(winner)) ! Stacks just before the hand showdown
            SS(i)=SS(i)-e
            SS(winner)=SS(winner)+e
         end if
      end do

      ! Losers. Share prizes, if any
      L=0
      do i=1,nr
         if (SS(i)==0) then
            L=L+1
            S(L)=S(i) ! The pre-showdown stacks, S(:), are the ones used to resolve ties
            w(L)=w(i)
         end if
      end do
      if (L>0) then ! Prizes are to be shared among losers
         call Solve_Bankruptcy_and_Ties(n,DCM,L,w,S,P)
         P(1:nr-L)=P(L+1:nr) ! Redefines initial part of prizes vector
      end if

      ! Winner/s
      L=0
      do i=1,nr
         if (SS(i)>0) then
            L=L+1
            S(L)=SS(i) ! The post-showdown become the new pre-showdown
            w(L)=ww(i)
         end if
      end do
      if (L==1) then ! Winner of the path! Path ends
         DCM(w(1))=DCM(w(1))+P(1)
         EXIT
      else if (Deep==MaxDeep) then ! Must finish!
         call Solve_Bankruptcy_and_Ties(n,DCM,L,w,S,P)
         EXIT
      end if
      nr=L

   end do ! Go to next node

   if (MOD(Iter,100000)==0) write(*,'("+",i13,99f9.2)') Iter,DCM(:)/Iter

END DO

DCM(:)=DCM(:)/Iter
write(*,'("+",i13,99f9.2)') Iter,DCM(:)
write(*,*)

END
!----------------------------------------------------------------------

!----------------------------------------------------------------------
! Solving ties when bankruptcy is met (http://arxiv.org/abs/2102.07738)
! This is almost the same as article Algorithm 2.
! The probability term is not needed.
! Player stacks and prizes are considered to be integer values.
! No optimized code (e.g. in sorting code, null prizes, ...)
!----------------------------------------------------------------------
subroutine Solve_Bankruptcy_and_Ties(N,DCM,L,w,S,P)
implicit none

! Input variables
integer N ! Total number of players in contest
double precision DCM(N) ! Cummulated DCM values
integer L ! Number of players that get into bankruptcy in this node.
integer w(L) ! Players' identities
integer S(L) ! Players' stacks
integer P(L) ! Prizes to be shared in this node (some can be zero)

! Working variables
integer PL ! Prize level
double precision SharedMoney

! Auxiliar variables
integer i,j,k,m

! Sort stacks and prizes from lower to higher. Stacks drag players' identities.
! A little bit redundant, and a shell sort code could be used here.
do i=1,L-1
   do j=i+1,L
      if (S(i)>S(j)) then
         k=S(i); S(i)=S(j); S(j)=k
         k=w(i); w(i)=w(j); w(j)=k
      end if
      if (P(i)>P(j)) then
         k=P(i); P(i)=P(j); P(j)=k
      end if
   end do
end do

! Let's go
i=0
PL=1
do while (i<L)
   i=i+1
   if (i<L) then
      k=i+1
      do
         if (S(i)==S(k)) then
            if (k==L) Exit
            k=k+1
         else
            k=k-1
            Exit
         end if
      end do
   else
      k=i
   end if
   ! Update the DCM values of involved players.
   SharedMoney=SUM(P(PL:PL+k-i))/(k-i+1) ! Shared among tied players
   do m=i,k
      DCM(w(m))=DCM(w(m))+SharedMoney
   end do
   PL=PL+k-i+1
   i=k
end do
End subroutine Solve_Bankruptcy_and_Ties
!----------------------------------------------------------------------

