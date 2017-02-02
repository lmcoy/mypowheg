#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cstring>

#include "powheg/color.h"

#define QUARK 2
#define ANTIQUARK 4
#define GLUON 8
#define PHOTON 16

namespace {

static constexpr int newcol = Powheg::NewColorID;

int get_type(int pdg) {
    if (pdg == 0 || pdg == 21) {
        return GLUON;
    }
    if (pdg == 22) {
        return PHOTON;
    }
    if (pdg > 0 && pdg <= 6) {
        return QUARK;
    }
    if (pdg < 0 && pdg >= -6) {
        return ANTIQUARK;
    }
    return 0;
}

bool IsValidISRSplitting(int born, int real, int rad) {
    auto isquark = [](int pdg)->bool {
        return pdg != 0 && pdg >= -6 && pdg <= 6;
    };
    auto isgluon = [](int pdg)->bool {
        return pdg == 0 || pdg == 21;
    };
    auto isphoton = [](int pdg)->bool {
        return pdg == 22;
    };

    if (isquark(born) && born == real && (isgluon(rad) || isphoton(rad))) {
        return true;
    }
    if (isquark(born) && rad == -born && isgluon(real)) {
        return true;
    }
    if (isgluon(born) && isquark(rad) && rad == real) {
        return true;
    }
    if (isgluon(born) && isgluon(real) && real == rad) {
        return true;
    }
    return false;
}

int getmother(int i, int j) {
    if (i < j) {
        return i;
    }
    return j;
}

int select_radiating_particle(int *pdgs_b, int *pdgs, int i, double y) {
    bool rad1 = false;
    bool rad2 = false;
    if (IsValidISRSplitting(pdgs_b[0], pdgs[0], pdgs[i])) {
        rad1 = true;
    }
    if (IsValidISRSplitting(pdgs_b[1], pdgs[1], pdgs[i])) {
        rad2 = true;
    }
    if (rad1 && !rad2) {
        return -1;
    }
    if (rad2 && !rad1) {
        return -2;
    }
    if (rad1 && rad2) {
        if (y > 0.0) {
            return -1;
        }
        return -2;
    }
    assert(0);
    return 0;
}

void color_fsr_splitting(int i, int j, int *pdgs_b, int *pdgs, int *color1,
                         int *color2) {
    int mother = getmother(i, j);

    int type_m = get_type(pdgs_b[mother]);
    int type_i = get_type(pdgs[i]);
    int type_j = get_type(pdgs[j]);

    int col1 = color1[mother];
    int col2 = color2[mother];

    if (type_m == GLUON) {
        if (type_i == GLUON && type_j == GLUON) {
            double rnd = drand48();
            if (rnd > 0.5) {
                color1[i] = col1;
                color2[i] = newcol;
                color1[j] = newcol;
                color2[j] = col2;
                return;
            } else {
                color1[i] = newcol;
                color2[i] = col2;
                color1[j] = col1;
                color2[j] = newcol;
                return;
            }
        }
        if (type_i == QUARK && type_j == ANTIQUARK) {
            color1[i] = col1;
            color2[i] = 0;
            color1[j] = 0;
            color2[j] = col2;
            return;
        }
        if (type_i == ANTIQUARK && type_j == QUARK) {
            color1[i] = 0;
            color2[i] = col2;
            color1[j] = col1;
            color2[j] = 0;
            return;
        }
    }
    if (type_m == QUARK) {
        if (type_i == GLUON) {
            color1[i] = col1;
            color2[i] = newcol;
            color1[j] = newcol;
            color2[j] = col2;
            return;
        }
        if (type_i == PHOTON) {
            color1[i] = 0;
            color2[i] = 0;
            color1[j] = col1;
            color2[j] = col2;
            return;
        }
        if (type_j == GLUON) {
            assert(0 && "not implemented");
            return;
        }
        if (type_j == PHOTON) {
            assert(0 && "not implemented");
            return;
        }
    }
    if (type_m == ANTIQUARK) {
        if (type_i == GLUON) {
            color1[i] = newcol;
            color2[i] = col2;
            color1[j] = col1;
            color2[j] = newcol;
            return;
        }
        if (type_i == PHOTON) {
            color1[i] = 0;
            color2[i] = 0;
            color1[j] = col1;
            color2[j] = col2;
            return;
        }
        if (type_j == GLUON) {
            assert(0 && "not implemented");
            return;
        }
        if (type_j == PHOTON) {
            assert(0 && "not implemented");
            return;
        }
    }
}

static void array_insert(size_t n, int * array, size_t i, int value) {
    size_t l = n - i - 1;
    assert(l < n);
    if (l == 0) {
        array[n-1] = value;
        return;
    }
    memmove(&array[i+1], &array[i], sizeof(int)*l);
    array[i] = value;
}

void color_isr_splitting(int i, int j, int *pdgs_b, int *pdgs, int *color1,
                         int *color2, double y) {
    int jp = select_radiating_particle(pdgs_b, pdgs, i, y);

    int mother = 0;
    if (jp == -1) {
        mother = 0;
    }
    if (jp == -2) {
        mother = 1;
    }

    int type_m = get_type(pdgs_b[mother]);
    int type_i = get_type(pdgs[i]);

    int col1 = color1[mother];
    int col2 = color2[mother];
    size_t ncol = 10;

    if (type_m == GLUON) {
        if (type_i == QUARK) {
            color1[mother] = col1;
            color2[mother] = 0;
            array_insert(ncol, color1, i, col2);
            array_insert(ncol, color2, i, 0);
            return;
        }
        if (type_i == ANTIQUARK) {
            color1[mother] = 0;
            color2[mother] = col2;
            array_insert(ncol, color1, i, 0);
            array_insert(ncol, color2, i, col1);
            return;
        }
        double rnd = drand48();
        if (type_i == GLUON && rnd > 0.5) {
            color1[mother] = col1;
            color2[mother] = newcol;
            array_insert(ncol, color1, i, col2);
            array_insert(ncol, color2, i, newcol);
            return;
        }
        if (type_i == GLUON && rnd <= 0.5) {
            color1[mother] = newcol;
            color2[mother] = col2;
            array_insert(ncol, color1, i, newcol);
            array_insert(ncol, color2, i, col1);
            return;
        }
    }
    if (type_m == QUARK) {
        if (type_i == GLUON) {
            color1[mother] = newcol;
            color2[mother] = col2;
            array_insert(ncol, color1, i, newcol);
            array_insert(ncol, color2, i, col1);
            return;
        }
        if (type_i == PHOTON) {
            color1[mother] = col1;
            color2[mother] = col2;
            array_insert(ncol, color1, i, 0);
            array_insert(ncol, color2, i, 0);
            return;
        }
        if (type_i == ANTIQUARK) {
            color1[mother] = col1;
            color2[mother] = newcol;
            array_insert(ncol, color1, i, col2);
            array_insert(ncol, color2, i, newcol);
            return;
        }
    }
    if (type_m == ANTIQUARK) {
        if (type_i == GLUON) {
            color1[mother] = col1;
            color2[mother] = newcol;
            array_insert(ncol, color1, i, col2);
            array_insert(ncol, color2, i, newcol);
            return;
        }
        if (type_i == PHOTON) {
            color1[mother] = col1;
            color2[mother] = col2;
            array_insert(ncol, color1, i, 0);
            array_insert(ncol, color2, i, 0);
            return;
        }
        if (type_i == QUARK) {
            color1[mother] = newcol;
            color2[mother] = col2;
            array_insert(ncol, color1, i, newcol);
            array_insert(ncol, color2, i, col1);
            return;
        }
    }
}

} // namespace

namespace Powheg {
void Color(int i, int j, int *pdgs_b, int *pdgs, int *color1, int *color2,
           double y) {
    if (i == 0 && j == 0) {
        // born process, nothing to do
        return;
    }
    if (j >= 2) {
        color_fsr_splitting(i, j, pdgs_b, pdgs, color1, color2);
    } else {
        color_isr_splitting(i, j, pdgs_b, pdgs, color1, color2, y);
    }
}
}

