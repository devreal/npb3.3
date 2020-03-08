c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine  adi(rho_i, us, vs, ws, qs, square, rhs,
     $                forcing, u, nx, nxmax, ny, nz)

      integer nx, nxmax, ny, nz
      double precision rho_i  (  0:nxmax-1,0:ny-1,0:nz-1),
     $                 us     (  0:nxmax-1,0:ny-1,0:nz-1),
     $                 vs     (  0:nxmax-1,0:ny-1,0:nz-1),
     $                 ws     (  0:nxmax-1,0:ny-1,0:nz-1),
     $                 qs     (  0:nxmax-1,0:ny-1,0:nz-1),
     $                 square (  0:nxmax-1,0:ny-1,0:nz-1),
     $                 rhs    (5,0:nxmax-1,0:ny-1,0:nz-1),
     $                 forcing(5,0:nxmax-1,0:ny-1,0:nz-1),
     $                 u      (5,0:nxmax-1,0:ny-1,0:nz-1)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      call compute_rhs(rho_i, us, vs, ws, qs, square, rhs,
     $                 forcing, u, nx, nxmax, ny, nz)
c            call print_zone(forcing,
c     $              nx, nxmax, ny, nz, 4, 10000)
c            call print_zone(rhs,
c     $              nx, nxmax, ny, nz, 4, 10001)

      call x_solve(rho_i, qs, square, u, rhs, nx, nxmax, ny, nz)
c            call print_zone(rhs,
c     $              nx, nxmax, ny, nz, 4, 10002)

      call y_solve(rho_i, qs, square, u, rhs, nx, nxmax, ny, nz)
c            call print_zone(rhs,
c     $              nx, nxmax, ny, nz, 4, 10003)

      call z_solve(rho_i, qs, square, u, rhs, nx, nxmax, ny, nz)
c            call print_zone(rhs,
c     $              nx, nxmax, ny, nz, 4, 10004)

      call add(u, rhs, nx, nxmax, ny, nz)

c            call print_zone(u,
c     $              nx, nxmax, ny, nz, 4, 10005)
      return
      end

