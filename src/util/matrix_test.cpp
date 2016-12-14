#include "gtest/gtest.h"

#include "util/matrix.h"


TEST(Matrix,Index) {
    Util::Matrix<4,4> m;
    for (int i = 0; i < 4; i++) {
        for(int j = 0; j < 4; j++) {
            m[i][j] = (double)(i*(j+1));
        }
    }

    EXPECT_DOUBLE_EQ( m[0][0], 0.0 );
    EXPECT_DOUBLE_EQ( m[1][0], 1.0 );
    EXPECT_DOUBLE_EQ( m[2][0], 2.0 );
    EXPECT_DOUBLE_EQ( m[3][0], 3.0 );

    EXPECT_DOUBLE_EQ( m[0][1], 0.0 );
    EXPECT_DOUBLE_EQ( m[1][1], 2.0 );
    EXPECT_DOUBLE_EQ( m[2][1], 4.0 );
    EXPECT_DOUBLE_EQ( m[3][1], 6.0 );
    
    EXPECT_DOUBLE_EQ( m[0][2], 0.0 );
    EXPECT_DOUBLE_EQ( m[1][2], 3.0 );
    EXPECT_DOUBLE_EQ( m[2][2], 6.0 );
    EXPECT_DOUBLE_EQ( m[3][2], 9.0 );

    EXPECT_DOUBLE_EQ( m[0][3], 0.0 );
    EXPECT_DOUBLE_EQ( m[1][3], 4.0 );
    EXPECT_DOUBLE_EQ( m[2][3], 8.0 );
    EXPECT_DOUBLE_EQ( m[3][3], 12.0 );
}

TEST(Matrix, Copy) {
    Util::Matrix<4,4> m;
    for (int i = 0; i < 4; i++) {
        for(int j = 0; j < 4; j++) {
            m[i][j] = (double)(i*(j+1));
        }
    }
    Util::Matrix<4,4> m2 = m;
    for (int i = 0; i < 4; i++) {
        for(int j = 0; j < 4; j++) {
            m2[i][j] = 0.0;
        }
    }
    EXPECT_DOUBLE_EQ( m[0][0], 0.0 );
    EXPECT_DOUBLE_EQ( m[1][0], 1.0 );
    EXPECT_DOUBLE_EQ( m[2][0], 2.0 );
    EXPECT_DOUBLE_EQ( m[3][0], 3.0 );

    EXPECT_DOUBLE_EQ( m[0][1], 0.0 );
    EXPECT_DOUBLE_EQ( m[1][1], 2.0 );
    EXPECT_DOUBLE_EQ( m[2][1], 4.0 );
    EXPECT_DOUBLE_EQ( m[3][1], 6.0 );
    
    EXPECT_DOUBLE_EQ( m[0][2], 0.0 );
    EXPECT_DOUBLE_EQ( m[1][2], 3.0 );
    EXPECT_DOUBLE_EQ( m[2][2], 6.0 );
    EXPECT_DOUBLE_EQ( m[3][2], 9.0 );

    EXPECT_DOUBLE_EQ( m[0][3], 0.0 );
    EXPECT_DOUBLE_EQ( m[1][3], 4.0 );
    EXPECT_DOUBLE_EQ( m[2][3], 8.0 );
    EXPECT_DOUBLE_EQ( m[3][3], 12.0 );

    EXPECT_DOUBLE_EQ( m2[0][0], 0.0 );
    EXPECT_DOUBLE_EQ( m2[1][0], 0.0 );
    EXPECT_DOUBLE_EQ( m2[2][0], 0.0 );
    EXPECT_DOUBLE_EQ( m2[3][0], 0.0 );

    EXPECT_DOUBLE_EQ( m2[0][1], 0.0 );
    EXPECT_DOUBLE_EQ( m2[1][1], 0.0 );
    EXPECT_DOUBLE_EQ( m2[2][1], 0.0 );
    EXPECT_DOUBLE_EQ( m2[3][1], 0.0 );
    
    EXPECT_DOUBLE_EQ( m2[0][2], 0.0 );
    EXPECT_DOUBLE_EQ( m2[1][2], 0.0 );
    EXPECT_DOUBLE_EQ( m2[2][2], 0.0 );
    EXPECT_DOUBLE_EQ( m2[3][2], 0.0 );

    EXPECT_DOUBLE_EQ( m2[0][3], 0.0 );
    EXPECT_DOUBLE_EQ( m2[1][3], 0.0 );
    EXPECT_DOUBLE_EQ( m2[2][3], 0.0 );
    EXPECT_DOUBLE_EQ( m2[3][3], 0.0 );
}
