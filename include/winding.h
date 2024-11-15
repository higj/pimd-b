#pragma once

#include "common.h"

class WindingProbability {
public:
    WindingProbability(double diff_, int max_wind_, double beta_half_k_, double size_)
        : diff(diff_), max_wind(max_wind_), beta_half_k(beta_half_k_), size(size_),
        shift(getWindingShift(diff_)), denominator(computeDenominator()) {
    }

    [[nodiscard]] double getWindingShift(const double diff_) const {
        // Start from the zero winding
        double diff_shift = std::min(std::numeric_limits<double>::max(), diff_ * diff_);

        for (int wind_idx = 1; wind_idx <= max_wind; ++wind_idx) {
            const double diff_plus_wind = diff_ + wind_idx * size;
            const double diff_minus_wind = diff_ - wind_idx * size;

            diff_shift = std::min(diff_shift, diff_plus_wind * diff_plus_wind);
            diff_shift = std::min(diff_shift, diff_minus_wind * diff_minus_wind);
        }

        return diff_shift;
    }

    [[nodiscard]] double getLogWindingWeight() const {
        double result = 0.0;
        double weight = exp(-beta_half_k * (diff * diff - shift)); // Zero winding contribution

        for (int wind_num = 1; wind_num <= max_wind; ++wind_num) {
            const double diff_plus = diff + wind_num * size;
            const double diff_minus = diff - wind_num * size;
            weight += exp(-beta_half_k * (diff_plus * diff_plus - shift)) + exp(-beta_half_k * (diff_minus * diff_minus - shift));
        }

        result += log(weight) - beta_half_k * shift;

        return result;
    }

    // Get the probability for a specific winding number
    [[nodiscard]] double getProbability(const int winding) const {
        // Important note: The winding probability will be wrong if winding_number is greater than max_wind
        const double diff_val = diff + winding * size;
        return exp(-beta_half_k * (diff_val * diff_val - shift)) / denominator;
    }

    [[nodiscard]] double getExpectation() const {
        double wind_mean = 0.0;

        for (int wind_idx = 1; wind_idx <= max_wind; ++wind_idx) {
            wind_mean += wind_idx * getProbability(wind_idx);
            wind_mean -= wind_idx * getProbability(-wind_idx);
        }

        return wind_mean;
    }

    [[nodiscard]] double getSquaredExpectation() const {
        double wind_squared_mean = 0.0;

        for (int wind_idx = 1; wind_idx <= max_wind; ++wind_idx) {
            const double wind_squared = wind_idx * wind_idx;
            wind_squared_mean += wind_squared * getProbability(wind_idx);
            wind_squared_mean += wind_squared * getProbability(-wind_idx);
        }

        return wind_squared_mean;
    }

    [[nodiscard]] double getDiffSquaredExpectation() const {
        double diff_squared_mean = diff * diff * getProbability(0);

        for (int wind_idx = 1; wind_idx <= max_wind; ++wind_idx) {
            const double diff_plus = diff + wind_idx * size;
            const double diff_minus = diff - wind_idx * size;

            diff_squared_mean += diff_plus * diff_plus * getProbability(wind_idx);
            diff_squared_mean += diff_minus * diff_minus * getProbability(-wind_idx);
        }

        return diff_squared_mean;
    }

    [[nodiscard]] int getMinimumImageWindingNumber() const {
        return -static_cast<int>(std::floor(diff / size + 0.5));
    }

private:
    double diff;
    int max_wind;
    double beta_half_k;
    double size;
    double shift;
    double denominator;

    // Precompute the denominator when the object is created
    [[nodiscard]] double computeDenominator() const {
        // Start with the zero winding contribution
        double denom = exp(-beta_half_k * (diff * diff - shift));

        // Add the contribution of the nonzero winding numbers
        for (int wind_idx = 1; wind_idx <= max_wind; ++wind_idx) {
            const double diff_plus = diff + wind_idx * size;
            const double diff_minus = diff - wind_idx * size;

            denom += exp(-beta_half_k * (diff_plus * diff_plus - shift));
            denom += exp(-beta_half_k * (diff_minus * diff_minus - shift));
        }

        return denom;
    }
};
