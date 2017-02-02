#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>

#include "phasespace.h"

int main() {

    double S = 14000.0 * 14000.0;
    srand48(time(NULL));
    std::array<double, 7> v = {drand48(), drand48(), drand48(), drand48(),
                               drand48(), drand48(), drand48()};
    Phasespace::Phasespace ps = GenPhasespace(S, v);

    printf("ps.S = %.16g;\n", ps.S);
    printf("ps.X1 = %.16g;\n", ps.X1);
    printf("ps.X2 = %.16g; \n", ps.X2);
    printf("ps.Jacobian = %.16g;\n", ps.Jacobian);
    printf("ps.N = %d;\n", ps.N);
    for (int i = 0; i < ps.N + 2; i++) {
        std::cout << "ps.Momenta[" << i << "].Set" << ps.Momenta[i].ToString(16)
                  << ";\n";
    }
    return 0;
}
