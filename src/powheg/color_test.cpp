#include "gtest/gtest.h"

#include "powheg/color.h"

void color_reset(int *color, int *color_b) {
    for (int i = 0; i < 5; i++) {
        color[i] = color_b[i];
    }
    color[5] = 0;
}

int newcol = Powheg::NewColorID;;

TEST(WJ, Born_3_Real_32_54) {
    int pdgs_b[] = {2, -1, -13, 14, 0};
    int pdgs_r[] = {2, -1, -13, 14, 21, 21};

    int color1_b[] = {501, 0, 0, 0, 501};
    int color2_b[] = {0, 502, 0, 0, 502};

    int color1[10] = {0};
    int color2[10] = {0};

    color_reset(color1, color1_b);
    color_reset(color2, color2_b);

    Powheg::Color(5, 4, pdgs_b, pdgs_r, color1, color2, 0.5);

    EXPECT_EQ(color1[0], 501);
    EXPECT_EQ(color1[1], 0);
    EXPECT_EQ(color1[2], 0);
    EXPECT_EQ(color1[3], 0);

    EXPECT_EQ(color2[0], 0);
    EXPECT_EQ(color2[1], 502);
    EXPECT_EQ(color2[2], 0);
    EXPECT_EQ(color2[3], 0);

    bool col1 = color1[4] == 501 && color1[5] == newcol &&
                color2[4] == newcol && color2[5] == 502;
    bool col2 = color1[5] == 501 && color1[4] == newcol &&
                color2[5] == newcol && color2[4] == 502;
    EXPECT_TRUE(col1 ^ col2);
}

TEST(WJ, Born_3_Real_32_45) {
    int pdgs_b[] = {2, -1, -13, 14, 0};
    int pdgs_r[] = {2, -1, -13, 14, 21, 21};

    int color1_b[] = {501, 0, 0, 0, 501};
    int color2_b[] = {0, 502, 0, 0, 502};

    int color1[10] = {0};
    int color2[10] = {0};

    color_reset(color1, color1_b);
    color_reset(color2, color2_b);

    Powheg::Color(4, 5, pdgs_b, pdgs_r, color1, color2, 0.5);

    EXPECT_EQ(color1[0], 501);
    EXPECT_EQ(color1[1], 0);
    EXPECT_EQ(color1[2], 0);
    EXPECT_EQ(color1[3], 0);

    EXPECT_EQ(color2[0], 0);
    EXPECT_EQ(color2[1], 502);
    EXPECT_EQ(color2[2], 0);
    EXPECT_EQ(color2[3], 0);

    bool col1 = color1[4] == 501 && color1[5] == newcol &&
                color2[4] == newcol && color2[5] == 502;
    bool col2 = color1[5] == 501 && color1[4] == newcol &&
                color2[5] == newcol && color2[4] == 502;
    EXPECT_TRUE(col1 ^ col2);
}

TEST(WJ, Born_3_Real_32_50_0) {
    int pdgs_b[] = {2, -1, -13, 14, 0};
    int pdgs_r[] = {2, -1, -13, 14, 21, 21};

    int color1_b[] = {501, 0, 0, 0, 501};
    int color2_b[] = {0, 502, 0, 0, 502};

    int color1[10] = {0};
    int color2[10] = {0};

    color_reset(color1, color1_b);
    color_reset(color2, color2_b);

    Powheg::Color(5, 0, pdgs_b, pdgs_r, color1, color2, 0.5);

    EXPECT_EQ(color1[0], newcol);
    EXPECT_EQ(color1[1], 0);
    EXPECT_EQ(color1[2], 0);
    EXPECT_EQ(color1[3], 0);
    EXPECT_EQ(color1[4], 501);
    EXPECT_EQ(color1[5], newcol);

    EXPECT_EQ(color2[0], 0);
    EXPECT_EQ(color2[1], 502);
    EXPECT_EQ(color2[2], 0);
    EXPECT_EQ(color2[3], 0);
    EXPECT_EQ(color2[4], 502);
    EXPECT_EQ(color2[5], 501);
}


TEST(WJ, Born_3_Real_32_50_1) {
    int pdgs_b[] = {2, -1, -13, 14, 0};
    int pdgs_r[] = {2, -1, -13, 14, 21, 21};

    int color1_b[] = {501, 0, 0, 0, 501};
    int color2_b[] = {0, 502, 0, 0, 502};

    int color1[10] = {0};
    int color2[10] = {0};

    color_reset(color1, color1_b);
    color_reset(color2, color2_b);

    Powheg::Color(5, 0, pdgs_b, pdgs_r, color1, color2, -0.5);

    EXPECT_EQ(color1[0], 501);
    EXPECT_EQ(color1[1], 0);
    EXPECT_EQ(color1[2], 0);
    EXPECT_EQ(color1[3], 0);
    EXPECT_EQ(color1[4], 501);
    EXPECT_EQ(color1[5], 502);

    EXPECT_EQ(color2[0], 0);
    EXPECT_EQ(color2[1], newcol);
    EXPECT_EQ(color2[2], 0);
    EXPECT_EQ(color2[3], 0);
    EXPECT_EQ(color2[4], 502);
    EXPECT_EQ(color2[5], newcol);
}

TEST(WJ, Born_3_Real_32_40_0) {
    int pdgs_b[] = {2, -1, -13, 14, 0};
    int pdgs_r[] = {2, -1, -13, 14, 21, 21};

    int color1_b[] = {501, 0, 0, 0, 501};
    int color2_b[] = {0, 502, 0, 0, 502};

    int color1[10] = {0};
    int color2[10] = {0};

    color_reset(color1, color1_b);
    color_reset(color2, color2_b);

    Powheg::Color(4, 0, pdgs_b, pdgs_r, color1, color2, 0.5);

    EXPECT_EQ(color1[0], newcol);
    EXPECT_EQ(color1[1], 0);
    EXPECT_EQ(color1[2], 0);
    EXPECT_EQ(color1[3], 0);
    EXPECT_EQ(color1[4], newcol);
    EXPECT_EQ(color1[5], 501);

    EXPECT_EQ(color2[0], 0);
    EXPECT_EQ(color2[1], 502);
    EXPECT_EQ(color2[2], 0);
    EXPECT_EQ(color2[3], 0);
    EXPECT_EQ(color2[4], 501);
    EXPECT_EQ(color2[5], 502);
}

TEST(WJ, Born_3_Real_32_40_1) {
    int pdgs_b[] = {2, -1, -13, 14, 0};
    int pdgs_r[] = {2, -1, -13, 14, 21, 21};

    int color1_b[] = {501, 0, 0, 0, 501};
    int color2_b[] = {0, 502, 0, 0, 502};

    int color1[10] = {0};
    int color2[10] = {0};

    color_reset(color1, color1_b);
    color_reset(color2, color2_b);

    Powheg::Color(4, 0, pdgs_b, pdgs_r, color1, color2, -0.5);

    EXPECT_EQ(color1[0], 501);
    EXPECT_EQ(color1[1], 0);
    EXPECT_EQ(color1[2], 0);
    EXPECT_EQ(color1[3], 0);
    EXPECT_EQ(color1[4], 502);
    EXPECT_EQ(color1[5], 501);

    EXPECT_EQ(color2[0], 0);
    EXPECT_EQ(color2[1], newcol);
    EXPECT_EQ(color2[2], 0);
    EXPECT_EQ(color2[3], 0);
    EXPECT_EQ(color2[4], newcol);
    EXPECT_EQ(color2[5], 502);
}

TEST(WJ, Born_3_Real_25_54) {
    int pdgs_b[] = {2, -1, -13, 14, 0};
    int pdgs_r[] = {2, -1, -13, 14, -2, 2};

    int color1_b[] = {501, 0, 0, 0, 501};
    int color2_b[] = {0, 502, 0, 0, 502};

    int color1[10] = {0};
    int color2[10] = {0};

    color_reset(color1, color1_b);
    color_reset(color2, color2_b);

    Powheg::Color(5, 4, pdgs_b, pdgs_r, color1, color2, -0.5);

    EXPECT_EQ(color1[0], 501);
    EXPECT_EQ(color1[1], 0);
    EXPECT_EQ(color1[2], 0);
    EXPECT_EQ(color1[3], 0);
    EXPECT_EQ(color1[4], 0);
    EXPECT_EQ(color1[5], 501);

    EXPECT_EQ(color2[0], 0);
    EXPECT_EQ(color2[1], 502);
    EXPECT_EQ(color2[2], 0);
    EXPECT_EQ(color2[3], 0);
    EXPECT_EQ(color2[4], 502);
    EXPECT_EQ(color2[5], 0);
}

TEST(WJ, Born_3_Real_40_52) {
    int pdgs_b[] = {2, -1, -13, 14, 0};
    int pdgs_r[] = {2, -1, -13, 14, 21, 22};

    int color1_b[] = {501, 0, 0, 0, 501};
    int color2_b[] = {0, 502, 0, 0, 502};

    int color1[10] = {0};
    int color2[10] = {0};

    color_reset(color1, color1_b);
    color_reset(color2, color2_b);

    Powheg::Color(5, 2, pdgs_b, pdgs_r, color1, color2, -0.5);

    EXPECT_EQ(color1[0], 501);
    EXPECT_EQ(color1[1], 0);
    EXPECT_EQ(color1[2], 0);
    EXPECT_EQ(color1[3], 0);
    EXPECT_EQ(color1[4], 501);
    EXPECT_EQ(color1[5], 0);

    EXPECT_EQ(color2[0], 0);
    EXPECT_EQ(color2[1], 502);
    EXPECT_EQ(color2[2], 0);
    EXPECT_EQ(color2[3], 0);
    EXPECT_EQ(color2[4], 502);
    EXPECT_EQ(color2[5], 0);
}

TEST(WJ, Born_3_Real_40_50) {
    int pdgs_b[] = {2, -1, -13, 14, 0};
    int pdgs_r[] = {2, -1, -13, 14, 21, 22};

    int color1_b[] = {501, 0, 0, 0, 501};
    int color2_b[] = {0, 502, 0, 0, 502};

    int color1[10] = {0};
    int color2[10] = {0};

    color_reset(color1, color1_b);
    color_reset(color2, color2_b);

    Powheg::Color(5, 0, pdgs_b, pdgs_r, color1, color2, -0.5);

    EXPECT_EQ(color1[0], 501);
    EXPECT_EQ(color1[1], 0);
    EXPECT_EQ(color1[2], 0);
    EXPECT_EQ(color1[3], 0);
    EXPECT_EQ(color1[4], 501);
    EXPECT_EQ(color1[5], 0);

    EXPECT_EQ(color2[0], 0);
    EXPECT_EQ(color2[1], 502);
    EXPECT_EQ(color2[2], 0);
    EXPECT_EQ(color2[3], 0);
    EXPECT_EQ(color2[4], 502);
    EXPECT_EQ(color2[5], 0);
}

TEST(WJ, Born_1_Real_34_50) {
    int pdgs_b[] = {21, -1, -13, 14, -2};
    int pdgs_r[] = {-2, -1, -13, 14, -2, -2};

    int color1_b[] = {501, 0, 0, 0, 0};
    int color2_b[] = {502, 501, 0, 0, 502};

    int color1[10] = {0};
    int color2[10] = {0};

    color_reset(color1, color1_b);
    color_reset(color2, color2_b);

    Powheg::Color(5, 0, pdgs_b, pdgs_r, color1, color2, -0.5);

    EXPECT_EQ(color1[0], 0);
    EXPECT_EQ(color1[1], 0);
    EXPECT_EQ(color1[2], 0);
    EXPECT_EQ(color1[3], 0);
    EXPECT_EQ(color1[4], 0);
    EXPECT_EQ(color1[5], 0);

    EXPECT_EQ(color2[0], 502);
    EXPECT_EQ(color2[1], 501);
    EXPECT_EQ(color2[2], 0);
    EXPECT_EQ(color2[3], 0);
    EXPECT_EQ(color2[4], 502);
    EXPECT_EQ(color2[5], 501);
}

TEST(WJ, Born_1_Real_34_40) {
    int pdgs_b[] = {21, -1, -13, 14, -2};
    int pdgs_r[] = {-2, -1, -13, 14, -2, -2};

    int color1_b[] = {501, 0, 0, 0, 0};
    int color2_b[] = {502, 501, 0, 0, 502};

    int color1[10] = {0};
    int color2[10] = {0};

    color_reset(color1, color1_b);
    color_reset(color2, color2_b);

    Powheg::Color(4, 0, pdgs_b, pdgs_r, color1, color2, -0.5);

    EXPECT_EQ(color1[0], 0);
    EXPECT_EQ(color1[1], 0);
    EXPECT_EQ(color1[2], 0);
    EXPECT_EQ(color1[3], 0);
    EXPECT_EQ(color1[4], 0);
    EXPECT_EQ(color1[5], 0);

    EXPECT_EQ(color2[0], 502);
    EXPECT_EQ(color2[1], 501);
    EXPECT_EQ(color2[2], 0);
    EXPECT_EQ(color2[3], 0);
    EXPECT_EQ(color2[4], 501);
    EXPECT_EQ(color2[5], 502);
}

TEST(WJ, Born_1_Real_30_54) {
    int pdgs_b[] = {21, -1, -13, 14, -2};
    int pdgs_r[] = {21, -1, -13, 14, -2, 21};

    int color1_b[] = {501, 0, 0, 0, 0};
    int color2_b[] = {502, 501, 0, 0, 502};

    int color1[10] = {0};
    int color2[10] = {0};

    color_reset(color1, color1_b);
    color_reset(color2, color2_b);

    Powheg::Color(5, 4, pdgs_b, pdgs_r, color1, color2, -0.5);

    EXPECT_EQ(color1[0], 501);
    EXPECT_EQ(color1[1], 0);
    EXPECT_EQ(color1[2], 0);
    EXPECT_EQ(color1[3], 0);
    EXPECT_EQ(color1[4], 0);
    EXPECT_EQ(color1[5], newcol);

    EXPECT_EQ(color2[0], 502);
    EXPECT_EQ(color2[1], 501);
    EXPECT_EQ(color2[2], 0);
    EXPECT_EQ(color2[3], 0);
    EXPECT_EQ(color2[4], newcol);
    EXPECT_EQ(color2[5], 502);
}

TEST(WJ, Born_1_Real_30_50_0) {
    int pdgs_b[] = {21, -1, -13, 14, -2};
    int pdgs_r[] = {21, -1, -13, 14, -2, 21};

    int color1_b[] = {501, 0, 0, 0, 0};
    int color2_b[] = {502, 501, 0, 0, 502};

    int color1[10] = {0};
    int color2[10] = {0};

    color_reset(color1, color1_b);
    color_reset(color2, color2_b);

    Powheg::Color(5, 0, pdgs_b, pdgs_r, color1, color2, 0.5);

    EXPECT_EQ(color1[1], 0);
    EXPECT_EQ(color1[2], 0);
    EXPECT_EQ(color1[3], 0);
    EXPECT_EQ(color1[4], 0);

    EXPECT_EQ(color2[1], 501);
    EXPECT_EQ(color2[2], 0);
    EXPECT_EQ(color2[3], 0);
    EXPECT_EQ(color2[4], 502);

    bool c1 = color1[0] == 501 && color2[0] == newcol && color1[5] == 502 && color2[5] == newcol;
    bool c2 = color1[0] == newcol && color2[0] == 502 && color1[5] == newcol && color2[5] == 501;

    EXPECT_TRUE(c1 ^ c2);
}

TEST(WJ, Born_1_Real_30_50_1) {
    int pdgs_b[] = {21, -1, -13, 14, -2};
    int pdgs_r[] = {21, -1, -13, 14, -2, 21};

    int color1_b[] = {501, 0, 0, 0, 0};
    int color2_b[] = {502, 501, 0, 0, 502};

    int color1[10] = {0};
    int color2[10] = {0};

    color_reset(color1, color1_b);
    color_reset(color2, color2_b);

    Powheg::Color(5, 0, pdgs_b, pdgs_r, color1, color2, -0.5);

    EXPECT_EQ(color1[0], 501);
    EXPECT_EQ(color1[1], 0);
    EXPECT_EQ(color1[2], 0);
    EXPECT_EQ(color1[3], 0);
    EXPECT_EQ(color1[4], 0);
    EXPECT_EQ(color1[5], 501);

    EXPECT_EQ(color2[0], 502);
    EXPECT_EQ(color2[1], newcol);
    EXPECT_EQ(color2[2], 0);
    EXPECT_EQ(color2[3], 0);
    EXPECT_EQ(color2[4], 502);
    EXPECT_EQ(color2[5], newcol);
}
