#pragma once
#ifndef SLOPE_INTERCEPT_H
#define SLOPE_INTERCEPT_H

#include "init_KARL.h"
#include "ip_bound.h"

//Used in polynomial kernel function
void poly_rotate_up(double x_min, double x_max, double deg, double& m, double& c);
void poly_rotate_down(double x_min, double x_max, double deg, double& m, double& c);
void poly_tangent(double*q, double*a_G, double b_G, model& our_model, double deg, double& m, double& c);
void poly_chord(double x_min, double x_max, double deg, double& m, double& c);

//Used in sigmoid kernel function
void sigmoid_chord(double x_min, double x_max, double& m, double& c);
void sigmoid_tangent(double*q, double*a_G, double b_G, model& our_model, double& m, double& c);

void sigmoid_rotate_up(double x_min, double x_max, double& m, double& c);
void sigmoid_rotate_down(double x_min, double x_max, double& m, double& c);

#endif