      module utils


      CONTAINS
      ! gets the first non blank characters of a string
      character(len=256) function string_trim(strin,delimiter)
      implicit none

      character(*), intent(in) :: strin
      character(len=1), intent(in) :: delimiter

      integer :: i,j

      do i=1,len_trim(strin)
        if (strin(i:i) /= delimiter)then
          string_trim=trim(strin(1:i))
          cycle
        else
          exit
        end if
      end do
      end function string_trim
      

      end module utils


