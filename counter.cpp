#include "counter.h"

void Counter::increment (int value) {
    assert(0 <= value);
    if ( value > limit ) {
        ++counts[limit + 1];
    } else {
        ++counts[value];
    }
}

void Counter::nullify () {
    for ( auto & c : counts ) {
        c = 0;
    }
}

Counter & Counter::operator+=(Counter const & rhs) {
    assert(rhs.limit == limit);
    for ( int i = 0; i <= limit + 1; ++i ) {
        counts[i] += rhs.counts[i];
    }
    return *this;

}

std::vector<double> Counter::getDistribution() const {
    double sum = 0.0;
    for ( auto const & c : counts ) {
        sum += c;
    }

    std::vector<double> distr(limit + 2);
    for ( int i = 0; i <= limit + 1; ++i ) {
        distr[i] = counts[i] / sum;
    }
    return distr;
}
