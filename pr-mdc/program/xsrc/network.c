#define XR_INCLUDED
#include  <stdlib.h>
#include  <math.h>
#include  "../md.h"
#include  "../prototypes.h"

/**** Network Resolver
*  
*   The number [i] atom and [i+1] atom will be connected. 
*           net.pos[0]   <->   net.pos[1]
*           net.pos[2]   <->   net.pos[3]
*           .............................
*   Total number of connection is (net.counter/2).
*   The "net.length" is defined in network_resolver structure of 
*   "program/md.h".  see also init_mem in "init.c", and "xsrc/xnet.c".
***/
void   network(void)
{
  int i, j;
  static int cnt=0;
  double dr2, dx, dy, dz;
  double rdx, rdy, rdz;
  double rxi, ryi, rzi;

  if(cnt==0){
    net.length=10.0; /* initialize the network-length [A] (see xnet.c) */
    cnt=1;
  }

  net.counter = 0;
  rdx = cell.hm[0][0]*0.5;
  rdy = cell.hm[1][1]*0.5;
  rdz = cell.hm[2][2]*0.5;

  for(i=0; i<sys.N; i++) {
    rxi = cell.rx[i];
    ryi = cell.ry[i];
    rzi = cell.rz[i];

    for(j=i+1; j<sys.N; j++) {
      dx = rxi - cell.rx[j];
      dy = ryi - cell.ry[j];
      dz = rzi - cell.rz[j];

      /* neglect the PBC effects */
      if(dx > rdx) continue; 
      if(dy > rdy) continue; 
      if(dz > rdz) continue; 
      if(dx < -rdx) continue; 
      if(dy < -rdy) continue; 
      if(dz < -rdz) continue; 

      dr2 = dx*dx + dy*dy + dz*dz;
      /* check the network distances */
      if(dr2 > net.length) continue; 
      if(net.counter>10*sys.N) continue;
      net.pos[net.counter]=(int)i;
      net.pos[net.counter+1]=(int)j;
      net.counter +=2;
    }
  }
}
