
/*
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine  set_constants

c---------------------------------------------------------------------
c---------------------------------------------------------------------
*/


#include "header.h"

double  tx1, tx2, tx3, ty1, ty2, ty3, tz1, tz2, tz3,
        dx1, dx2, dx3, dx4, dx5, dy1, dy2, dy3, dy4,
        dy5, dz1, dz2, dz3, dz4, dz5, dssp, dt,
        dxmax, dymax, dzmax, xxcon1, xxcon2,
        xxcon3, xxcon4, xxcon5, dx1tx1, dx2tx1, dx3tx1,
        dx4tx1, dx5tx1, yycon1, yycon2, yycon3, yycon4,
        yycon5, dy1ty1, dy2ty1, dy3ty1, dy4ty1, dy5ty1,
        zzcon1, zzcon2, zzcon3, zzcon4, zzcon5, dz1tz1,
        dz2tz1, dz3tz1, dz4tz1, dz5tz1, dnxm1, dnym1,
        dnzm1, c1c2, c1c5, c3c4, c1345, conz1, c1, c2,
        c3, c4, c5, c4dssp, c5dssp, dtdssp, dttx1,
        dttx2, dtty1, dtty2, dttz1, dttz2, c2dttx1,
        c2dtty1, c2dttz1, comz1, comz4, comz5, comz6,
        c3c4tx3, c3c4ty3, c3c4tz3, c2iv, con43, con16;


StaticArrayT2<double, size_type, 5, 13, 1, 1> ce/*(5,13)*/;

void set_constants()
{

  ce(1,1)  = 2.0;
  ce(1,2)  = 0.0;
  ce(1,3)  = 0.0;
  ce(1,4)  = 4.0;
  ce(1,5)  = 5.0;
  ce(1,6)  = 3.0;
  ce(1,7)  = 0.5;
  ce(1,8)  = 0.02;
  ce(1,9)  = 0.01;
  ce(1,10) = 0.03;
  ce(1,11) = 0.5;
  ce(1,12) = 0.4;
  ce(1,13) = 0.3;

  ce(2,1)  = 1.0;
  ce(2,2)  = 0.0;
  ce(2,3)  = 0.0;
  ce(2,4)  = 0.0;
  ce(2,5)  = 1.0;
  ce(2,6)  = 2.0;
  ce(2,7)  = 3.0;
  ce(2,8)  = 0.01;
  ce(2,9)  = 0.03;
  ce(2,10) = 0.02;
  ce(2,11) = 0.4;
  ce(2,12) = 0.3;
  ce(2,13) = 0.5;

  ce(3,1)  = 2.0;
  ce(3,2)  = 2.0;
  ce(3,3)  = 0.0;
  ce(3,4)  = 0.0;
  ce(3,5)  = 0.0;
  ce(3,6)  = 2.0;
  ce(3,7)  = 3.0;
  ce(3,8)  = 0.04;
  ce(3,9)  = 0.03;
  ce(3,10) = 0.05;
  ce(3,11) = 0.3;
  ce(3,12) = 0.5;
  ce(3,13) = 0.4;

  ce(4,1)  = 2.0;
  ce(4,2)  = 2.0;
  ce(4,3)  = 0.0;
  ce(4,4)  = 0.0;
  ce(4,5)  = 0.0;
  ce(4,6)  = 2.0;
  ce(4,7)  = 3.0;
  ce(4,8)  = 0.03;
  ce(4,9)  = 0.05;
  ce(4,10) = 0.04;
  ce(4,11) = 0.2;
  ce(4,12) = 0.1;
  ce(4,13) = 0.3;

  ce(5,1)  = 5.0;
  ce(5,2)  = 4.0;
  ce(5,3)  = 3.0;
  ce(5,4)  = 2.0;
  ce(5,5)  = 0.1;
  ce(5,6)  = 0.4;
  ce(5,7)  = 0.3;
  ce(5,8)  = 0.05;
  ce(5,9)  = 0.04;
  ce(5,10) = 0.03;
  ce(5,11) = 0.1;
  ce(5,12) = 0.3;
  ce(5,13) = 0.2;

  c1 = 1.4;
  c2 = 0.4;
  c3 = 0.1;
  c4 = 1.0;
  c5 = 1.4;

  /**
   * The following three settings are based on average cell size,
   * not actual cell size.
   */

  dnxm1 = 1.0 / (dble(gx_size-1)/dble(x_zones));
  dnym1 = 1.0 / (dble(gy_size-1)/dble(y_zones));
  dnzm1 = 1.0 / (dble(gz_size-1));

  c1c2 = c1 * c2;
  c1c5 = c1 * c5;
  c3c4 = c3 * c4;
  c1345 = c1c5 * c3c4;

  conz1 = (1.0-c1c5);

  tx1 = 1.0 / (dnxm1 * dnxm1);
  tx2 = 1.0 / (2.0 * dnxm1);
  tx3 = 1.0 / dnxm1;

  ty1 = 1.0 / (dnym1 * dnym1);
  ty2 = 1.0 / (2.0 * dnym1);
  ty3 = 1.0 / dnym1;

  tz1 = 1.0 / (dnzm1 * dnzm1);
  tz2 = 1.0 / (2.0 * dnzm1);
  tz3 = 1.0 / dnzm1;

  dx1 = 0.75;
  dx2 = 0.75;
  dx3 = 0.75;
  dx4 = 0.75;
  dx5 = 0.75;

  dy1 = 0.75;
  dy2 = 0.75;
  dy3 = 0.75;
  dy4 = 0.75;
  dy5 = 0.75;

  dz1 = 1.0;
  dz2 = 1.0;
  dz3 = 1.0;
  dz4 = 1.0;
  dz5 = 1.0;

  dxmax = dmax1(dx3, dx4);
  dymax = dmax1(dy2, dy4);
  dzmax = dmax1(dz2, dz3);

  dssp = 0.25 * dmax1(dx1, dmax1(dy1, dz1) );

  c4dssp = 4.0 * dssp;
  c5dssp = 5.0 * dssp;

  dttx1 = dt*tx1;
  dttx2 = dt*tx2;
  dtty1 = dt*ty1;
  dtty2 = dt*ty2;
  dttz1 = dt*tz1;
  dttz2 = dt*tz2;

  c2dttx1 = 2.0*dttx1;
  c2dtty1 = 2.0*dtty1;
  c2dttz1 = 2.0*dttz1;

  dtdssp = dt*dssp;

  comz1  = dtdssp;
  comz4  = 4.0*dtdssp;
  comz5  = 5.0*dtdssp;
  comz6  = 6.0*dtdssp;

  c3c4tx3 = c3c4*tx3;
  c3c4ty3 = c3c4*ty3;
  c3c4tz3 = c3c4*tz3;

  dx1tx1 = dx1*tx1;
  dx2tx1 = dx2*tx1;
  dx3tx1 = dx3*tx1;
  dx4tx1 = dx4*tx1;
  dx5tx1 = dx5*tx1;

  dy1ty1 = dy1*ty1;
  dy2ty1 = dy2*ty1;
  dy3ty1 = dy3*ty1;
  dy4ty1 = dy4*ty1;
  dy5ty1 = dy5*ty1;

  dz1tz1 = dz1*tz1;
  dz2tz1 = dz2*tz1;
  dz3tz1 = dz3*tz1;
  dz4tz1 = dz4*tz1;
  dz5tz1 = dz5*tz1;

  c2iv  = 2.5;
  con43 = 4.0/3.0;
  con16 = 1.0/6.0;

  xxcon1 = c3c4tx3*con43*tx3;
  xxcon2 = c3c4tx3*tx3;
  xxcon3 = c3c4tx3*conz1*tx3;
  xxcon4 = c3c4tx3*con16*tx3;
  xxcon5 = c3c4tx3*c1c5*tx3;

  yycon1 = c3c4ty3*con43*ty3;
  yycon2 = c3c4ty3*ty3;
  yycon3 = c3c4ty3*conz1*ty3;
  yycon4 = c3c4ty3*con16*ty3;
  yycon5 = c3c4ty3*c1c5*ty3;

  zzcon1 = c3c4tz3*con43*tz3;
  zzcon2 = c3c4tz3*tz3;
  zzcon3 = c3c4tz3*conz1*tz3;
  zzcon4 = c3c4tz3*con16*tz3;
  zzcon5 = c3c4tz3*c1c5*tz3;

  return;

}
