#include "geom.h"

void wrap_coord_into_box(double &x, const double L) {
    if (x > L) x -= L;
    else if (x < 0) x += L;
}

void wrap_coord_relative(double &x, const double x_ref, const double L) {
    if (x - x_ref > 0.5 * L) x -= L;
    else if (x - x_ref < -0.5 * L) x += L;
}

double calc_sq_distance(const double *R, const double *r, const double L)
{
    double sq_distance = 0.0;
    for (int dim = 0; dim < 3; dim++)
    {
        double r_dim = R[dim] - r[dim];
        if (r_dim > 0.5 * L) sq_distance += (r_dim - L) * (r_dim - L);
        else if (r_dim < -0.5 * L) sq_distance += (r_dim + L) * (r_dim + L);
        else sq_distance += r_dim * r_dim;
    }
    return sq_distance;
}