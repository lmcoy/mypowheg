#ifndef MATH_DILOG_H_H
#define MATH_DILOG_H_H

#include <cmath>
#include <cassert>

namespace Math {

namespace {
// order of Li2 expansion (max=22)
const int n_dilog = 18;

double dilog_expansion1(double x) {
    // expansion of Li2(x) at x = 0
    static const double c[22] = {
        0.0,         1.0,         0.25,        1. / 9.0,    1.0 / 16.0,
        1.0 / 25.0,  1.0 / 36.0,  1.0 / 49.0,  1.0 / 64.0,  1.0 / 81.0,
        0.01,        1.0 / 121.0, 1.0 / 144.0, 1.0 / 169.0, 1.0 / 196.0,
        1.0 / 225.0, 1.0 / 256.0, 1.0 / 289.0, 1.0 / 324.0, 1.0 / 361.0,
        1.0 / 400.0, 1.0 / 441.0
    };
    double res = 0.0;
    double xx = x;
    for (int i = 1; i < n_dilog; i++) {
        res += c[i] * xx;
        xx *= x;
    }
    return res;
}

double dilog_expansion2(double x) {
    // expansion of Li2(x) at x = 1/2
    static const double c[22] = {
        0.5822405264650125059026563, 1.386294361119890618834464,
        0.6137056388801093811655358, 0.5150591481598541584459523,
        0.5607446110935520956644048, 0.7028086222503166469369522,
        0.9619856295828055884384129, 1.398500825477095181724625,
        2.124052126843654860553334,  3.335018441166835803460740,
        5.374744583677473331548447,  8.845918938768230306275552,
        14.81278497589460807485846,  25.16614055732277483615874,
        43.28485984409286871087003,  75.23969019578854983494738,
        131.9922475495631357261403,  233.4263575537635092213829,
        415.7632854049826525164728,  745.2438335016702958167416,
        1343.447242662615911632402,  2434.271918737874454033521
    };
    double res = c[0];
    double xx = x - 0.5;
    for (int i = 1; i < n_dilog; i++) {
        res += c[i] * xx;
        xx *= (x - 0.5);
    }
    return res;
}

} // end namespace

/**
 * Dilog returns Li2(x) for x in [0,1].
 */
double Dilog(double x) {
    assert(x >= 0.0 && x <= 1.0 && "implementation only valid for x in [0,1]");
    if (x == 1.0) {
        return M_PI * M_PI / 6.0;
    }
    if (x == 0.0) {
        return 0.0;
    }
    if (x < 1.0 / 3.0) {
        return dilog_expansion1(x);
    }
    if (x < 2.0 / 3.0) {
        return dilog_expansion2(x);
    }
    return -dilog_expansion1(1.0 - x) + M_PI * M_PI / 6.0 -
           log(x) * log(1.0 - x);
}

} // end namespace Math
#endif
