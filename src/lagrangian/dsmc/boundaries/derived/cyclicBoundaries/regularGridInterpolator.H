#ifndef REGULARGRIDINTERPOLATOR_H
#define REGULARGRIDINTERPOLATOR_H

#include <limits>      // for std::numeric_limits
#include <cstdint>     // for SIZE_MAX
#include <vector>
#include <iostream>

class RegularGridInterpolator {
public:
    RegularGridInterpolator(const std::vector<double>& x_vals,
                            const std::vector<double>& y_vals,
                            const std::vector<std::vector<double>>& z_vals);

    double interpolate(double x, double y);

private:
    size_t findIndex(double value, const std::vector<double>& values);

    std::vector<double> x_vals_;
    std::vector<double> y_vals_;
    std::vector<std::vector<double>> z_vals_;
};

#endif // REGULARGRIDINTERPOLATOR_H
