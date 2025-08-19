#ifndef TRIPLET_HPP
#define TRIPLET_HPP

#include <string>

template <typename T1, typename T2, typename T3>
struct Triplet
{
    T1 first;
    T2 second;
    T3 third;

    Triplet(T1 f, T2 s, T3 t) : first(f), second(s), third(t) {}
    Triplet(const Triplet &other) : first(other.first), second(other.second), third(other.third) {}
    Triplet& operator=(const Triplet &other) {
        if (this != &other) {
            first = other.first;
            second = other.second;
            third = other.third;
        }
        return *this;
    }

    std::string join(const std::string &sep) const
    {
        return first + sep + second + sep + third;
    }

    // Comparison operators for std::set compatibility
    bool operator<(const Triplet& other) const {
        if (first != other.first) return first < other.first;
        if (second != other.second) return second < other.second;
        return third < other.third;
    }

    bool operator==(const Triplet& other) const {
        return first == other.first && second == other.second && third == other.third;
    }

    bool operator!=(const Triplet& other) const {
        return !(*this == other);
    }
};

#endif // TRIPLET_HPP