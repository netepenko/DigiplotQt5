
! f90 version

subroutine find_peaks(N, delta, val, results, min_pos, max_pos)
      
implicit none

! values can be returned via array pos      

integer, intent(in):: N
real*8, intent(in):: delta
real*8, dimension(:), intent(in) ::val
integer, dimension(2), intent(inout):: results
integer, dimension(:), intent(inout):: min_pos
integer, dimension(:), intent(inout):: max_pos

real*8:: inf = 1.0D30

integer :: i, max_i, min_i
integer :: n_max, n_min

logical :: look_for_max

real*8 :: max_val, min_val, this
real*8 ::frac, rn

max_val = -inf
min_val = inf

n_max = 0
n_min = 0

look_for_max = .True.
rn = real(N)


print *, 'Searching : ', N, ' points'
do i = 1, N
   this = val(i)
   frac = real(i)/rn
   if (this .gt. max_val) then
      max_val = this
      max_i = i-1 ! make sure the indices start at 0
   endif
   if (this .lt. min_val) then
      min_val = this
      min_i = i-1
   endif
   if (look_for_max) then
      if (this .lt. (max_val-delta)) then
         n_max = n_max + 1
         max_pos(n_max) = max_i
         min_val = this
         min_i = i-1
         look_for_max = .False.
      endif
   else
      if (this .gt. min_val+delta) then
         n_min = n_min + 1
         min_pos(n_min) = min_i
         max_val = this
         max_i = i-1
         look_for_max = .True.
      endif
   endif
enddo
print *, 'PD-found : ', n_max, ' peaks'
!do i=1, n_max
!   print *, i, max_pos(i)
!enddo
results(1) = n_min
results(2) = n_max

return

end subroutine find_peaks
