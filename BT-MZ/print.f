c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine print_zone(u, nx, nxmax, ny, nz, zone, step)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      integer nx, nxmax, ny, nz, zone, step
      double precision u      (5,0:nxmax-1,0:ny-1,0:nz-1)

      integer i, j, k, m

      write(*, 80) zone, step
  80  format('------ Printing U of zone',i6,' at step ',i6,' -------')
      do k = 0, nz-1
         do j = 0, ny-1
            do i = 0, nx-1
                write(*, 100) u(1,i,j,k),u(2,i,j,k),u(3,i,j,k),
     $                        u(4,i,j,k),u(5,i,j,k)
            enddo
         enddo
      enddo
      write(*, 110)
  100 format(F20.13, F20.13, F20.13, F20.13, F20.13)
  110 format('------ End of U ------- ')
  120 format(' ')
      write(*, 120)
      write(*, 120)

      return
      end




