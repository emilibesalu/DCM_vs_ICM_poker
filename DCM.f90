!----------------------------------------------------------------------
! Simple Fortran 90 DCM code.
! Code and procedure from  http://arxiv.org/abs/2102.07738  article
! By E. Besalu 2021
! No exhaustive checks made here (in order to maintain a clear code)
! http://iqcc.udg.edu/~emili/Poker/DCM/index.html
!----------------------------------------------------------------------
program Test_DCM_Code
implicit none

! Contest Variables
integer N ! Total number of players in contest
double precision, allocatable :: DCM(:) ! Cummulated DCM values
integer, allocatable :: S(:) ! Players' stacks
integer, allocatable :: P(:) ! Prizes to be shared (some can be zero)
double precision, allocatable :: ProbWin(:) ! Probability a player has to win

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

write(*,*) "For the proposed example the result must be 310.681696 257.845768 131.472535"
write(*,*)

allocate(ProbWin(N)) ! Probability a player has to win
allocate(DCM(N))     ! DCM values

call DCM_Calculation(N,S,P,ProbWin,DCM)

write(*,'("   Prizes:",99(i8,1x))') (P(i),i=1,N)
write(*,'("   Player:",99(1x,"--",i2," ---"))') (i,i=1,N)
write(*,'("   Stacks:",99(i8,1x))') (S(i),i=1,N)
write(*,'(" Prob.Win:",99f9.3)')  ProbWin(:)*100
write(*,'("      DCM:",99f9.3)')  DCM(:)
END
!----------------------------------------------------------------------

!----------------------------------------------------------------------
! DCM calling code
!----------------------------------------------------------------------
subroutine DCM_Calculation(N,S,P,ProbWin,DCM)
implicit none

! Input variables
integer N ! Number of players
integer S(N) ! Stacks
integer P(N) ! Prizes (as many as players, insert zeroes if necessary)

! Output results
double precision ProbWin(N) ! Probability each player has to win
double precision DCM(N) ! Expected prizes distribution

! Needed for the recursive routine New_Hand()
integer w(N) ! Player identities
double precision Prob ! Node probability
integer Deep ! Actual tree deep

! Auxiliar
integer i,j,k,np

! Identities initialization
do i=1,N
   w(i)=i
end do

! Initial Winning probabilities and DCM values to be updated
ProbWin(:)=0.0d0
DCM(:)=0.0d0

Prob=1.0d0/N ! Node probability

! The prizes for the N podium positions. Sorted from lower to higher.
do i=1,N-1
   do j=i+1,N
      if (P(i)>P(j)) then
         k=P(i); P(i)=P(j); P(j)=k
      end if
   end do
end do

Deep=1 ! Actual tree deep

! Here, at the first call, n=N (i.e., the set of involved players in the hand are all the players)
np=N
call New_Hand(Deep,Prob,np,w,S,P,N,ProbWin,DCM)

END
!----------------------------------------------------------------------

!----------------------------------------------------------------------
! Solving ties when bankruptcy is met (http://arxiv.org/abs/2102.07738)
! This corresponds to the article Algorithm 2
! Player stacks and prizes are considered to be integer values.
! No optimized code (e.g. in sorting code, null prizes, ...)
!----------------------------------------------------------------------
subroutine Solve_Bankruptcy(Prob,N,DCM,L,w,S,P)
implicit none

! Input variables
integer N ! Total number of players in contest
double precision DCM(N) ! Cummulated DCM values
double precision prob ! Probability attached to the actual node.
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
   SharedMoney=Prob*SUM(P(PL:PL+k-i))/(k-i+1) ! Shared among tied players
   do m=i,k
      DCM(w(m))=DCM(w(m))+SharedMoney
   end do
   PL=PL+k-i+1
   i=k
end do
End subroutine Solve_Bankruptcy
!----------------------------------------------------------------------

!----------------------------------------------------------------------
! New_Hand() routine. No optimized code. For basic concepts exhibition.
! Specific MD & ProbMin values are set as fixed parameters.
! Algorithm taken from http://arxiv.org/abs/2102.07738 article.
!----------------------------------------------------------------------
recursive subroutine New_Hand(Deep,Prob,n,w,S,P,Ntp,ProbWin,DCM)
implicit none

integer MD ! Maximal acceptable tree deep search.
parameter (MD=90)
double precision ProbMin ! Minimal acceptable probability for a node.
parameter (ProbMin=1.0d-15)

! Variables
integer Deep
double precision Prob
integer n ! Number of involved players in the hand.
integer w(n) ! Player identities
integer S(n) ! Stacks
integer P(n) ! Prizes still to be distributed (sorted from lower to higher)
integer Ntp ! Total number of players
double precision ProbWin(Ntp)
double precision DCM(Ntp)

! Auxiliar variables for new calls (avoiding collateral effects)
integer L,np ! List of players
integer wp(n) ! Player identities
integer Sp(n),Spp(n) ! Stacks
integer Ppp(n) ! Prizes still to be distributed (sorted from lower to higher)
integer Deepp
double precision Probp

! Auxiliar
integer i,j,k,e

! Go ahead!
if (Deep==MD .or. Prob<ProbMin) then ! Tree search termination forced
   call Solve_Bankruptcy(Prob,Ntp,DCM,n,w,S,P)
   return
end if
if (n==1) then ! Assign first prize to winner
   DCM(w(1))=DCM(w(1))+Prob*P(1) ! DCM value update
   ProbWin(w(1))=ProbWin(w(1))+Prob ! Probability to win value update
   return
end if

! One player, one at a time, wins the all-in situation
do j=1,n ! Winnig player is the j-th one

   Sp(:)=S(:)

   ! Player j takes all
   do k=1,n
      if (k/=j) then
         e=MIN(S(j),S(k)) ! Comparison among original stacks
         Sp(j)=Sp(j)+e
         Sp(k)=Sp(k)-e
      end if
   end do

   ! List of bankrupted players
   L=0
   do k=1,n
      if (Sp(k)==0) then ! Bankrupcy
         L=L+1
         wp(L)=w(k)
         Spp(L)=S(k) ! The evaluating stack is the original one when entering the node
      end if
   end do
   if (L>0) then ! There are some player/s eliminated
      Ppp(1:L)=P(1:L) ! Sorted list of prizes
      Call Solve_Bankruptcy(Prob,Ntp,DCM,L,wp,Spp,Ppp)
   end if

   ! List of still remaining players in play
   np=0
   do k=1,n
      if (Sp(k)>0) then ! Alive!
         np=np+1
         wp(np)=w(k)
         Spp(np)=Sp(k) ! The continuation stack is the actual one
      end if
   end do
   Ppp(1:np)=P(n-np+1:n) ! List of remaining prizes (still sorted)
   Probp=Prob/np ! Next node probability term
   Deepp=Deep+1 ! Tree search deep will be incremented
   call New_Hand(Deepp,Probp,np,wp,Spp,Ppp,Ntp,ProbWin,DCM)

end do

END
!----------------------------------------------------------------------
