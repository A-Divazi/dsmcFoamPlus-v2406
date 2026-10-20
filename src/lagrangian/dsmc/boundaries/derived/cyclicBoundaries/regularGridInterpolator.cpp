#include <iostream>
#include <vector>

class RegularGridInterpolator {
public:
    RegularGridInterpolator(const std::vector<double>& x_vals,
                            const std::vector<double>& y_vals,
                            const std::vector<std::vector<double>>& z_vals)
            : x_vals_(x_vals), y_vals_(y_vals), z_vals_(z_vals) {}

    double interpolate(double x, double y) {
        size_t x_index = findIndex(x, x_vals_);
        size_t y_index = findIndex(y, y_vals_);

        if (x_index == SIZE_MAX || y_index == SIZE_MAX) {
            return 0.0; // Return default value or handle out-of-bounds case
        }

        double x1 = x_vals_[x_index];
        double x2 = x_vals_[x_index + 1];
        double y1 = y_vals_[y_index];
        double y2 = y_vals_[y_index + 1];

        double z11 = z_vals_[x_index][y_index];
        double z12 = z_vals_[x_index][y_index + 1];
        double z21 = z_vals_[x_index + 1][y_index];
        double z22 = z_vals_[x_index + 1][y_index + 1];

        double a = (x2 - x) / (x2 - x1);
        double b = (x - x1) / (x2 - x1);
        double c = (y2 - y) / (y2 - y1);
        double d = (y - y1) / (y2 - y1);

        return a * c * z11 + a * d * z12 + b * c * z21 + b * d * z22;
    }

private:
    size_t findIndex(double value, const std::vector<double>& values) {
        for (size_t i = 0; i < values.size() - 1; ++i) {
            if (value >= values[i] && value <= values[i + 1]) {
                return i;
            }
        }
        return SIZE_MAX; // Value is out of bounds
    }

    std::vector<double> x_vals_;
    std::vector<double> y_vals_;
    std::vector<std::vector<double>> z_vals_;
};

int main() {

    #include "grid.h"

    RegularGridInterpolator interp(x_vals, y_vals, z_vals);

    double x = 76;
    double y = 0.0;

    double interpolatedValue = interp.interpolate(x, y);

    std::cout << "Interpolated value at (" << x << ", " << y << "): " << interpolatedValue << std::endl;

    return 0;
}
