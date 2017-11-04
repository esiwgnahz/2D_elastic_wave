#include <stdio.h>
#include <math.h>
#include "petscvec.h"

void assemble_body_force( double dim_x, double dim_y, int nNode, double *node2xy, int *node2DOF, Vec f, int *node_load ) {

	int i0;
	double x, y, x0, y0, r, dist_min;

	printf("----- force vector assembly -----\n");

	x0 = .3*dim_x;
	y0 = .3*dim_y;
	dist_min = sqrt( dim_x*dim_x + dim_y*dim_y );

	for( i0=0; i0<nNode; i0++ ) {
		x = node2xy[2*i0  ] - x0;
		y = node2xy[2*i0+1] - y0;
		r = sqrt( x*x + y*y );
		if( r<dist_min ) {
			dist_min = r;
			*node_load = i0;
		}
	}

	VecSet( f, 0. );
	VecSetValue( f, node2DOF[(*node_load)*2], 1., INSERT_VALUES );
	VecAssemblyBegin( f );
	VecAssemblyEnd  ( f );

	printf(" loading node=%i, coordinate=%f %f\n",(*node_load),node2xy[2*(*node_load)],node2xy[2*(*node_load)+1]);
}

double loading_fw_time_signal( double t, double offset ) {

	/* derivative of Gaussian */
	if( t<2*offset ) {
		t -= offset;
		return -2*(512*t)*exp(-(512*t)*(512*t));
	} else {
		return 0.;
	}
}