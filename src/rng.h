//
// Created by Jerry Vinokurov on 8/31/20.
//

#include <random>

#ifndef SLEIPNIR_RNG_H
#define SLEIPNIR_RNG_H

namespace Sleipnir {

    static std::random_device rng;
    static std::mt19937 g(rng());

}

#endif //SLEIPNIR_RNG_H
