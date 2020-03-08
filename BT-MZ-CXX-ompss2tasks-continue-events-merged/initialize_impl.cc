#include "NArray.h"

#include "exact_solution.h"

void initialize_1(size_type k, double* u_ptr, size_type nx, size_type nxmax, size_type ny, size_type nz)
{
  ArrayViewT4<double, size_type> u{u_ptr, 5, nxmax, ny, nz};

  for (size_type j = 0; j <= ny-1; ++j) {
    // do i = 0, nx-1
    for (size_type i = 0; i <= nx-1; ++i) {
      // do m = 1, 5
      for (size_type m = 1; m <= 5; ++m) {
        u(m,i,j,k) = 1.0;
      } // end do
    } // end do
  } //end do
}

void initialize_2(size_type k, double* u_ptr, size_type nx, size_type nxmax, size_type ny, size_type nz)
{
  ArrayViewT4<double, size_type> u{u_ptr, 5, nxmax, ny, nz};
  StaticArrayT3<double, size_type, 5, 3, 2, 1, 1, 1> Pface{};
  // do j = 0, ny-1
  double zeta = dble(k) * dnzm1;
  for (size_type j = 0; j <= ny-1; ++j) {
    double eta = dble(j) * dnym1;
    // do i = 0, nx-1
    for (size_type i = 0; i <= nx-1; ++i) {
      double xi = dble(i) * dnxm1;

      // do ix = 1, 2
      for  (size_type ix = 1; ix <= 2; ++ix) {
        exact_solution(dble(ix-1), eta, zeta, &Pface(1,1,ix));
      } // enddo

      // do iy = 1, 2
      for  (size_type iy = 1; iy <= 2; ++iy) {
        exact_solution(xi, dble(iy-1) , zeta, &Pface(1,2,iy));
      } // enddo

      // do iz = 1, 2
      for  (size_type iz = 1; iz <= 2; ++iz) {
        exact_solution(xi, eta, dble(iz-1), &Pface(1,3,iz));
      } // enddo

      // do m = 1, 5
      for (int m = 1; m <= 5; ++m) {
        double Pxi   = xi   * Pface(m,1,2) +
                    (1.0-xi)   * Pface(m,1,1);
        double Peta  = eta  * Pface(m,2,2) +
                    (1.0-eta)  * Pface(m,2,1);
        double Pzeta = zeta * Pface(m,3,2) +
                    (1.0-zeta) * Pface(m,3,1);

        u(m,i,j,k) = Pxi + Peta + Pzeta -
                    Pxi*Peta - Pxi*Pzeta - Peta*Pzeta +
                    Pxi*Peta*Pzeta;

      } // enddo
    } // enddo
  } // enddo
}


//---------------------------------------------------------------------
//     now store the exact values on the boundaries
//---------------------------------------------------------------------

//---------------------------------------------------------------------
//     west face
//---------------------------------------------------------------------
void initialize_3(size_type k, double* u_ptr, size_type nx, size_type nxmax, size_type ny, size_type nz)
{
  ArrayViewT4<double, size_type> u{u_ptr, 5, nxmax, ny, nz};
  StaticArrayT3<double, size_type, 5, 3, 2, 1, 1, 1> Pface{};
  StaticArrayT1<double, size_type, 5, 1> temp{};
  size_type i = 0;
  double xi = 0.0;
  double zeta = dble(k) * dnzm1;
  //do j = 0, ny-1
  for (size_type j = 0; j <= ny-1; ++j) {
    double eta = dble(j) * dnym1;
    exact_solution(xi, eta, zeta, temp.begin());
    //do m = 1, 5
    for (int m = 1; m <= 5; ++m) {
      u(m,i,j,k) = temp(m);
    } // enddo
  } // enddo
}

//---------------------------------------------------------------------
//     east face
//---------------------------------------------------------------------

void initialize_4(size_type k, double* u_ptr, size_type nx, size_type nxmax, size_type ny, size_type nz)
{
  ArrayViewT4<double, size_type> u{u_ptr, 5, nxmax, ny, nz};
  StaticArrayT1<double, size_type, 5, 1> temp{};
  size_type i = nx-1;
  double xi = 1.0;
  double zeta = dble(k) * dnzm1;
  //do j = 0, ny-1
  for (size_type j = 0; j <= ny-1; ++j) {
    double eta = dble(j) * dnym1;
    exact_solution(xi, eta, zeta, temp.begin());
    // do m = 1, 5
    for (int m = 1; m <= 5; ++m) {
      u(m,i,j,k) = temp(m);
    } // enddo
  } // enddo
}

//---------------------------------------------------------------------
//     south face
//---------------------------------------------------------------------
void initialize_5(size_type k, double* u_ptr, size_type nx, size_type nxmax, size_type ny, size_type nz)
{
  ArrayViewT4<double, size_type> u{u_ptr, 5, nxmax, ny, nz};
  StaticArrayT1<double, size_type, 5, 1> temp{};
  size_type j = 0;
  double eta = 0.0;
  double zeta = dble(k) * dnzm1;
  // do i = 0, nx-1
  for (size_type i = 0; i <= nx-1; ++i) {
    double xi = dble(i) * dnxm1;
    exact_solution(xi, eta, zeta, temp.begin());
    // do m = 1, 5
    for (int m = 1; m <= 5; ++m) {
        u(m,i,j,k) = temp(m);
    } // enddo
  } // enddo
}


//---------------------------------------------------------------------
//     north face
//---------------------------------------------------------------------
void initialize_6(size_type k, double* u_ptr, size_type nx, size_type nxmax, size_type ny, size_type nz)
{
  ArrayViewT4<double, size_type> u{u_ptr, 5, nxmax, ny, nz};
  StaticArrayT1<double, size_type, 5, 1> temp{};
  size_type j = ny-1;
  double eta = 1.0;
  double zeta = dble(k) * dnzm1;
  // do i = 0, nx-1
  for (size_type i = 0; i <= nx-1; ++i) {
    double  xi = dble(i) * dnxm1;
    exact_solution(xi, eta, zeta, temp.begin());
    // do m = 1, 5
    for (int m = 1; m <= 5; ++m) {
      u(m,i,j,k) = temp(m);
    }// enddo
  }// enddo
}

//---------------------------------------------------------------------
//     bottom face
//---------------------------------------------------------------------
void initialize_7(size_type j, double* u_ptr, size_type nx, size_type nxmax, size_type ny, size_type nz)
{
  ArrayViewT4<double, size_type> u{u_ptr, 5, nxmax, ny, nz};
  StaticArrayT1<double, size_type, 5, 1> temp{};
  size_type k = 0;
  double zeta = 0.0;
  double  eta = dble(j) * dnym1;
  // do i =0, nx-1
  for (size_type i = 0; i <= nx-1; ++i) {
    double xi = dble(i) *dnxm1;
    exact_solution(xi, eta, zeta, temp.begin());
    //do m = 1, 5
    for (int m = 1; m <= 5; ++m) {
      u(m,i,j,k) = temp(m);
    } // enddo
  } // enddo
}

//---------------------------------------------------------------------
//     top face
//---------------------------------------------------------------------
void initialize_8(size_type j, double* u_ptr, size_type nx, size_type nxmax, size_type ny, size_type nz)
{
  ArrayViewT4<double, size_type> u{u_ptr, 5, nxmax, ny, nz};
  StaticArrayT1<double, size_type, 5, 1> temp{};
  size_type k = nz-1;
  double zeta = 1.0;
  double eta = dble(j) * dnym1;
  //do i =0, nx-1
  for (size_type i = 0; i <= nx-1; ++i) {
    double xi = dble(i) * dnxm1;
    exact_solution(xi, eta, zeta, temp.begin());
    // do m = 1, 5
    for (int m = 1; m <= 5; ++m) {
      u(m,i,j,k) = temp(m);
    } // enddo
  } // enddo
}
