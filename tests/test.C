#include <iostream>
#include <vector>
#include <unistd.h>
#include "TBenchmark.h"
#include "TStopwatch.h"
#include "../include/my_func.h"


void test() {
    std::vector<int> teste = {1, 2, 3};
    TBenchmark benchmark;
    TStopwatch stopwatch, stopwatch1, stopwatch2;

    for (int i = 0; i < 2; i++) {
        std::cout << "teste" << std::endl;
        if (i > 1)
        std::cout << "teste3" << std::endl;

    }
    
}