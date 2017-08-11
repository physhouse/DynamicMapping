#ifndef GEOM_H
#define GEOM_H

void wrap_coord_into_box(double &x, const double L);

void wrap_coord_relative(double &x, const double x_ref, const double L);

double calc_sq_distance(const double *R, const double *r, const double L);

#endif
