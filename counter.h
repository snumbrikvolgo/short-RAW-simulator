#ifndef IITP_COUNTER_H
#define IITP_COUNTER_H
#pragma once
#include <limits>
#include <cassert>
#include <cstdint>
#include <vector>

class Counter {
 public:
    explicit Counter(int limit)
    : limit(limit) {
        counts = std::vector<int64_t>(limit + 2, 0);
    }

    void increment(int value);

    void nullify();

    Counter & operator+=(Counter const & rhs);

    std::vector<int64_t> getCounts() const {
        return counts;
    }

    std::vector<double> getDistribution() const;

 private:
    int limit;
    std::vector<int64_t> counts;
};

#endif //IITP_COUNTER_H
