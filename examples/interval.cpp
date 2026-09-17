/// @file interval.cpp
/// @brief Example to demonstrating the interval class.
///
/// This program creates a real interval with open/closed boundaries, inspects its
/// properties, bisects it and maps it linearly onto another interval.
///
/// @author Roland Abel
/// @date September 17, 2026
///
/// Copyright (c) 2026 Roland Abel

#include <iostream>
#include <xpolynomial.h>

using namespace std;
using namespace xmath;

auto main() -> int {
    auto I = interval(1., 5.); // half-open: (1,5]

    cout << "I = [" << I.lower() << ", " << I.upper() << "]" << endl
         << "length = " << I.length() << endl
         << "is_empty = " << boolalpha << I.is_empty() << endl
         << "is_closed = " << I.is_closed() << endl
         << "is_half_open = " << I.is_half_open() << endl << endl;

    auto [left, right] = I.bisect();
    cout << "bisect:" << endl
         << "left  = [" << left.lower() << ", " << left.upper() << "]" << endl
         << "right = [" << right.lower() << ", " << right.upper() << "]" << endl << endl;

    auto map = I.linear_transform(interval(0., 1.));
    cout << "linear transform onto [0, 1]:" << endl
         << "map(I.lower) = " << map(I.lower()) << endl
         << "map(I.upper) = " << map(I.upper()) << endl;

    return 0;
}