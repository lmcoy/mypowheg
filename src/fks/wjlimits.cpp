
// bornpdgs = 2, -1, -13, 14, 0
// realpdgs = 2, -1, -13, 14, 21, 21
TEST(QCDCollinearFSR, Wj_born_3_real_32) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 3.907985046680551e-14;
    ps.X2 = 0.0009853946746503084;
    ps.Jacobian = 1.571531091663275e-12;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(4.033619322210627e-05, 0, 0, 4.033619322210627e-05);
    ps.Momenta[1].Set(4.033619322210627e-05, 0, 0, -4.033619322210627e-05);
    ps.Momenta[2].Set(1.197745634288066e-05, 1.068563839979435e-05,
                      -5.06149328043127e-06, 1.912558026989396e-06);
    ps.Momenta[3].Set(3.003797300357577e-05, 8.75773984736797e-06,
                      -1.710411577839625e-05, 2.30874649619352e-05);
    ps.Momenta[4].Set(3.865695709775614e-05, -1.944337824716231e-05,
                      2.216560905882752e-05, -2.500002298892459e-05);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.305891697171142e-18;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(3.166715203862643e-19, -5.349805401075655e-36);
    spin[1][2] =
        std::complex<double>(7.845810300233236e-20, 3.930067047385589e-19);
    spin[1][3] =
        std::complex<double>(-1.76723537985724e-19, 3.484489985706078e-19);
    spin[2][1] =
        std::complex<double>(7.845810300233236e-20, -3.930067047385589e-19);
    spin[2][2] =
        std::complex<double>(5.071815226714396e-19, 2.139922160430262e-35);
    spin[2][3] =
        std::complex<double>(3.886595138121028e-19, 3.056548394890632e-19);
    spin[3][1] =
        std::complex<double>(-1.76723537985724e-19, -3.484489985706078e-19);
    spin[3][2] =
        std::complex<double>(3.886595138121028e-19, -3.056548394890632e-19);
    spin[3][3] =
        std::complex<double>(4.820386541134384e-19, -1.069961080215131e-35);

    // radiation variables
    double phi = 2.855289411123452;
    double xi = 0.4718482340746312;
    int bornpdgs[] = {2, -1, -13, 14, 0};
    int realpdgs[] = {2, -1, -13, 14, 21, 21};

    // FKS g = real matrix element * (1.0 - y) * xi * xi
    // real = 3.53656943764855
    double fks_g = 7.873844971986796e-09;
    double limit = FKS::QCD::CollinearLimitFSR(
        realpdgs[4], bornpdgs[4], ps, 4, xi, phi, alphas, born,
        spin); // limit = 7.873862492601418e-09
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, -1, -13, 14, 0
// realpdgs = 2, -1, -13, 14, 21, 21
TEST(QCDSoftCollinearFSR, Wj_born_3_real_32) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.2331784335639391;
    ps.X2 = 0.8312917879809874;
    ps.Jacobian = 563.4942164290239;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2861.767136122338, 0, 0, 2861.767136122338);
    ps.Momenta[1].Set(2861.767136122338, 0, 0, -2861.767136122338);
    ps.Momenta[2].Set(2765.702666419686, 220.616275966976, -1173.152791443982,
                      -2494.825089208985);
    ps.Momenta[3].Set(2762.463004523391, -38.96371823281797, 1239.966177166542,
                      2468.231666498278);
    ps.Momenta[4].Set(195.3686013015997, -181.652557734158, -66.8133857225595,
                      26.59342271070611);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 2.888278107118312e-06;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(4.916177823892196e-08, 1.176434648964306e-23);
    spin[1][2] =
        std::complex<double>(1.491703712983412e-08, 1.117961972459648e-09);
    spin[1][3] =
        std::complex<double>(3.732885617411153e-07, 2.808770623535729e-09);
    spin[2][1] =
        std::complex<double>(1.491703712983412e-08, -1.117961972459648e-09);
    spin[2][2] = std::complex<double>(4.551662769748373e-09, 0);
    spin[2][3] =
        std::complex<double>(1.133299004625928e-07, -7.636499218472382e-09);
    spin[3][1] =
        std::complex<double>(3.732885617411153e-07, -2.808770623535729e-09);
    spin[3][2] =
        std::complex<double>(1.133299004625928e-07, 7.636499218472382e-09);
    spin[3][3] = std::complex<double>(2.834564666109642e-06, 0);

    // radiation variables
    double phi = 1.873628098368453;
    int bornpdgs[] = {2, -1, -13, 14, 0};
    int realpdgs[] = {2, -1, -13, 14, 21, 21};

    // FKS g = real matrix element * (1.0 - y) * xi * xi
    // real = 336542913.139859
    double fks_g = 1.568488721588071e-12;
    double limit = FKS::QCD::SoftCollinearLimitFSR(realpdgs[4], bornpdgs[4], ps,
                                                   4, phi, alphas, born, spin);
    // limit = 1.568857688225015e-12
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, -1, -13, 14, 0
// realpdgs = 2, -1, -13, 14, 21, 21
TEST(QCDSoftLimit, Wj_born_3_real_32_4) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.8759808671584857;
    ps.X2 = 0.5315568646241715;
    ps.Jacobian = 1581.038953306218;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(4435.427986783298, 0, 0, 4435.427986783298);
    ps.Momenta[1].Set(4435.427986783298, 0, 0, -4435.427986783298);
    ps.Momenta[2].Set(4139.169990888839, 2449.678499761794, -2018.125682349901,
                      -2656.872633665238);
    ps.Momenta[3].Set(4378.009132554647, -2580.700036498459, 2346.457112223858,
                      2645.957351667411);
    ps.Momenta[4].Set(353.6768501231129, 131.0215367366658, -328.3314298739572,
                      10.91528199782754);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    // double born = 6.332522298326254e-07;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, -1.055420383054376e-07);
    ColorCorr.Set(0, 4, 9.498783447489381e-07);
    ColorCorr.Set(1, 0, -1.055420383054376e-07);
    ColorCorr.Set(1, 4, 9.498783447489381e-07);
    ColorCorr.Set(4, 0, 9.498783447489381e-07);
    ColorCorr.Set(4, 1, 9.498783447489381e-07);

    // radiation variables
    double y = -0.8464508702323386;
    double phi = 5.122516878283343;
    int bornpdgs[] = {2, -1, -13, 14, 0};
    int realpdgs[] = {2, -1, -13, 14, 21, 21};

    // FKS g = real matrix element * (1.0 - y) * xi * xi
    // real = 103732.6230425458
    double fks_g = 1.217854259935875e-13;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 4, alphas, ColorCorr, y, phi);
    // limit = 1.217854199521461e-13
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, -1, -13, 14, 0
// realpdgs = 2, -1, -13, 14, 21, 21
TEST(QCDCollinearISR1, Wj_born_3_real_32) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.9848910176245091;
    ps.X2 = 0.1183517093587625;
    ps.Jacobian = 526.6002472768062;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2219.191152992676, 0, 0, 2219.191152992676);
    ps.Momenta[1].Set(2219.191152992676, 0, 0, -2219.191152992676);
    ps.Momenta[2].Set(2061.000998049534, 1733.614369304639, 170.1335497582736,
                      -1101.52662598246);
    ps.Momenta[3].Set(2141.938221278446, -1889.597652688503, -284.8398601213539,
                      967.5672115226005);
    ps.Momenta[4].Set(235.4430866573712, 155.9832833838641, 114.7063103630803,
                      133.959414459859);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.788691288438134e-06;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(3.631546644747231e-07, 0);
    spin[1][2] =
        std::complex<double>(2.800591976659514e-07, 8.60282298169047e-09);
    spin[1][3] =
        std::complex<double>(-6.626679770613571e-07, -7.366395911144522e-09);
    spin[2][1] =
        std::complex<double>(2.800591976659514e-07, -8.60282298169047e-09);
    spin[2][2] =
        std::complex<double>(2.161810667477059e-07, -4.705738595857224e-23);
    spin[2][3] =
        std::complex<double>(-5.112136836884315e-07, 1.001718752254131e-08);
    spin[3][1] =
        std::complex<double>(-6.626679770613571e-07, 7.366395911144522e-09);
    spin[3][2] =
        std::complex<double>(-5.112136836884315e-07, -1.001718752254131e-08);
    spin[3][3] =
        std::complex<double>(1.209355557215705e-06, 2.941086622410765e-24);

    // radiation variables
    double phi = 4.278365629852333;
    double xi = 0.00581081945588911;
    int bornpdgs[] = {2, -1, -13, 14, 0};
    int realpdgs[] = {2, -1, -13, 14, 21, 21};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 2.113807580765473
    double fks_g = 1.427480585968142e-12;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, 1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 1.427851289566259e-12
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, -1, -13, 14, 0
// realpdgs = 2, -1, -13, 14, 21, 21
TEST(QCDSoftCollinearISR1, Wj_born_3_real_32) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.7527066611667799;
    ps.X2 = 0.006229453387625483;
    ps.Jacobian = 74.99298422253216;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(445.0934534425409, 0, 0, 445.0934534425409);
    ps.Momenta[1].Set(445.0934534425409, 0, 0, -445.0934534425409);
    ps.Momenta[2].Set(437.7239657504573, 160.2498704310209, 313.815996844995,
                      259.6955320049435);
    ps.Momenta[3].Set(285.2888288289449, -78.08010427724273, -239.0002489327759,
                      -134.8039101091547);
    ps.Momenta[4].Set(167.1741123056796, -82.16976615377818, -74.81574791221908,
                      -124.8916218957888);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 5.476594358144457e-07;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(2.166181761697318e-07, -1.176434648964306e-23);
    spin[1][2] =
        std::complex<double>(1.386938280722629e-07, -3.406425856854789e-08);
    spin[1][3] =
        std::complex<double>(-2.256031824302052e-07, 2.040603639536097e-08);
    spin[2][1] =
        std::complex<double>(1.386938280722629e-07, 3.406425856854789e-08);
    spin[2][2] =
        std::complex<double>(9.415808044280495e-08, 1.176434648964306e-23);
    spin[2][3] =
        std::complex<double>(-1.476555941067437e-07, -2.241184891581423e-08);
    spin[3][1] =
        std::complex<double>(-2.256031824302052e-07, -2.040603639536097e-08);
    spin[3][2] =
        std::complex<double>(-1.476555941067437e-07, 2.241184891581423e-08);
    spin[3][3] =
        std::complex<double>(2.368831792019089e-07, 2.941086622410765e-24);

    // radiation variables
    double phi = 1.985103300474277;
    int bornpdgs[] = {2, -1, -13, 14, 0};
    int realpdgs[] = {2, -1, -13, 14, 21, 21};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 89390556.24951422
    double fks_g = 1.093318345214027e-11;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, 1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 1.093122327552658e-11
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, -1, -13, 14, 0
// realpdgs = 2, -1, -13, 14, 21, 21
TEST(QCDSoftLimit, Wj_born_3_real_32_0) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.7017025611364041;
    ps.X2 = 0.3672466319664665;
    ps.Jacobian = 1357.26386145169;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3299.657007747561, 0, 0, 3299.657007747561);
    ps.Momenta[1].Set(3299.657007747561, 0, 0, -3299.657007747561);
    ps.Momenta[2].Set(2894.081526465613, -2758.617212940952, -780.9010226269794,
                      -394.882953741446);
    ps.Momenta[3].Set(3297.105717379912, 3147.972753194641, 843.2019966692561,
                      500.1840156746654);
    ps.Momenta[4].Set(408.1267716495974, -389.3555402536896, -62.30097404227676,
                      -105.3010619332193);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    // double born = 2.552036435384857e-07;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, -4.253394058974761e-08);
    ColorCorr.Set(0, 4, 3.828054653077285e-07);
    ColorCorr.Set(1, 0, -4.253394058974761e-08);
    ColorCorr.Set(1, 4, 3.828054653077285e-07);
    ColorCorr.Set(4, 0, 3.828054653077285e-07);
    ColorCorr.Set(4, 1, 3.828054653077285e-07);

    // radiation variables
    double y = -0.6577163407772204;
    double phi = 5.219939256185132;
    int bornpdgs[] = {2, -1, -13, 14, 0};
    int realpdgs[] = {2, -1, -13, 14, 21, 21};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 2390.730716928965
    double fks_g = 6.86453115373543e-14;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 0, alphas, ColorCorr, y, phi);
    // limit = 6.864531008826912e-14
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, -1, -13, 14, 0
// realpdgs = 2, -1, -13, 14, 21, 21
TEST(QCDCollinearISR2, Wj_born_3_real_32) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.8619053990967949;
    ps.X2 = 0.7692211391682449;
    ps.Jacobian = 26747.3477965829;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(5292.596223695325, 0, 0, 5292.596223695325);
    ps.Momenta[1].Set(5292.596223695325, 0, 0, -5292.596223695325);
    ps.Momenta[2].Set(4287.125027438549, 540.0532645737458, -406.6525429092141,
                      -4233.487590823752);
    ps.Momenta[3].Set(1283.753046029903, 368.5092458245906, -1145.539689313042,
                      -447.1707046994748);
    ps.Momenta[4].Set(5014.314373922201, -908.5625103983358, 1552.192232222255,
                      4680.658295523222);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 2.312271205551725e-07;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(1.028340238169638e-07, 0);
    spin[1][2] =
        std::complex<double>(-1.315140684877203e-08, 1.058597153250035e-07);
    spin[1][3] =
        std::complex<double>(2.432235877928692e-08, -3.510502528883323e-08);
    spin[2][1] =
        std::complex<double>(-1.315140684877203e-08, -1.058597153250035e-07);
    spin[2][2] = std::complex<double>(1.1065636069096e-07, 0);
    spin[2][3] =
        std::complex<double>(-3.9248500346199e-08, -2.054842772815292e-08);
    spin[3][1] =
        std::complex<double>(2.432235877928692e-08, 3.510502528883323e-08);
    spin[3][2] =
        std::complex<double>(-3.9248500346199e-08, 2.054842772815292e-08);
    spin[3][3] = std::complex<double>(1.773673604724873e-08, 0);

    // radiation variables
    double phi = 5.090376302500245;
    double xi = 0.01927459206107722;
    int bornpdgs[] = {2, -1, -13, 14, 0};
    int realpdgs[] = {2, -1, -13, 14, 21, 21};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 0.004309024985371799
    double fks_g = 3.201690873564955e-14;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, -1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 3.201780286365613e-14
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, -1, -13, 14, 0
// realpdgs = 2, -1, -13, 14, 21, 21
TEST(QCDSoftCollinearISR2, Wj_born_3_real_32) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.4642112321783465;
    ps.X2 = 0.569185048047558;
    ps.Jacobian = 5488.268509468664;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3341.164977634618, 0, 0, 3341.164977634618);
    ps.Momenta[1].Set(3341.164977634618, 0, 0, -3341.164977634618);
    ps.Momenta[2].Set(2136.04721269626, -1359.158486561564, -15.47942307520213,
                      1647.770096444404);
    ps.Momenta[3].Set(2916.472560417077, 1472.445119804218, 1616.706619591913,
                      -1929.760936229711);
    ps.Momenta[4].Set(1629.810182155899, -113.2866332426541, -1601.227196516711,
                      281.9908397853069);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.234289136552575e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(5.331465857591936e-10, 0);
    spin[1][2] =
        std::complex<double>(-1.203000665920395e-10, -6.477000382566862e-11);
    spin[1][3] =
        std::complex<double>(-4.689136595541281e-10, -3.677831936778958e-10);
    spin[2][1] =
        std::complex<double>(-1.203000665920395e-10, 6.477000382566862e-11);
    spin[2][2] =
        std::complex<double>(3.501337139962051e-11, -7.180387261745032e-28);
    spin[2][3] =
        std::complex<double>(1.5048706204798e-10, 2.602061710270543e-11);
    spin[3][1] =
        std::complex<double>(-4.689136595541281e-10, 3.677831936778958e-10);
    spin[3][2] =
        std::complex<double>(1.5048706204798e-10, -2.602061710270543e-11);
    spin[3][3] =
        std::complex<double>(6.661291793937612e-10, 4.59544784751682e-26);

    // radiation variables
    double phi = 2.468210925794011;
    int bornpdgs[] = {2, -1, -13, 14, 0};
    int realpdgs[] = {2, -1, -13, 14, 21, 21};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 1176.451382885885
    double fks_g = 4.367024606878663e-16;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, -1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 4.372015261917824e-16
    EXPECT_NEAR(limit / fks_g, 1.0, 5e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, -1, -13, 14, 0
// realpdgs = 2, -1, -13, 14, -2, 2
TEST(QCDCollinearFSR, Wj_born_3_real_25) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.1769277311465025;
    ps.X2 = 0.1118690230642088;
    ps.Jacobian = 124.1629753575585;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(914.4632007001608, 0, 0, 914.4632007001608);
    ps.Momenta[1].Set(914.4632007001608, 0, 0, -914.4632007001608);
    ps.Momenta[2].Set(800.7141710118716, 270.0810555264717, 395.0125079815961,
                      642.0004093791725);
    ps.Momenta[3].Set(893.4942777556722, -371.3228613798716, -365.5212738417503,
                      -725.8412742245367);
    ps.Momenta[4].Set(134.7179526327778, 101.2418058533999, -29.4912341398458,
                      83.84086484536422);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.012617844128507e-07;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(2.948411994749914e-08, 1.470543311205382e-24);
    spin[1][2] =
        std::complex<double>(-1.846066674870543e-08, -1.724824527733284e-09);
    spin[1][3] =
        std::complex<double>(-4.209705373932526e-08, -6.067113464460378e-10);
    spin[2][1] =
        std::complex<double>(-1.846066674870543e-08, 1.724824527733284e-09);
    spin[2][2] =
        std::complex<double>(1.16595386625195e-08, 2.941086622410765e-24);
    spin[2][3] =
        std::complex<double>(2.63933993004135e-08, -2.082807116673139e-09);
    spin[3][1] =
        std::complex<double>(-4.209705373932526e-08, 6.067113464460378e-10);
    spin[3][2] =
        std::complex<double>(2.63933993004135e-08, 2.082807116673139e-09);
    spin[3][3] =
        std::complex<double>(6.011812580283209e-08, 1.470543311205382e-24);

    // radiation variables
    double phi = 0.02569947906114568;
    double xi = 0.1015508507984796;
    int bornpdgs[] = {2, -1, -13, 14, 0};
    int realpdgs[] = {2, -1, -13, 14, -2, 2};

    // FKS g = real matrix element * (1.0 - y) * xi * xi
    // real = 0.0003439061808661828
    double fks_g = 3.546558403414632e-14;
    double limit = FKS::QCD::CollinearLimitFSR(
        realpdgs[4], bornpdgs[4], ps, 4, xi, phi, alphas, born,
        spin); // limit = 3.545838867546777e-14
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, -1, -13, 14, 0
// realpdgs = 2, -1, -13, 14, -2, 2
TEST(QCDSoftCollinearFSR, Wj_born_3_real_25) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.0005027845529248509;
    ps.X2 = 0.1929519171122678;
    ps.Jacobian = 1.707022785428134;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(64.02194571285122, 0, 0, 64.02194571285122);
    ps.Momenta[1].Set(64.02194571285122, 0, 0, -64.02194571285122);
    ps.Momenta[2].Set(49.87471505361805, -1.735571411463013, -49.84418872930428,
                      -0.1784473950588898);
    ps.Momenta[3].Set(51.71403749357381, -23.51703535535106, 45.55356123924987,
                      6.794393306531273);
    ps.Momenta[4].Set(26.4551388785106, 25.25260676681408, 4.290627490054411,
                      -6.615945911472384);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 0.0005230268022084378;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(4.052728902155148e-05, 1.505836350674312e-21);
    spin[1][2] =
        std::complex<double>(-2.813225842771274e-05, 1.012582832660559e-05);
    spin[1][3] =
        std::complex<double>(0.0001364452889576454, 6.566885213249359e-06);
    spin[2][1] =
        std::complex<double>(-2.813225842771274e-05, -1.012582832660559e-05);
    spin[2][2] = std::complex<double>(2.205813379395091e-05, 0);
    spin[2][3] =
        std::complex<double>(-9.307355781642977e-05, -3.864958455549571e-05);
    spin[3][1] =
        std::complex<double>(0.0001364452889576454, -6.566885213249359e-06);
    spin[3][2] =
        std::complex<double>(-9.307355781642977e-05, 3.864958455549571e-05);
    spin[3][3] = std::complex<double>(0.0004604413793929354, 0);

    // radiation variables
    double phi = 5.361774725498351;
    int bornpdgs[] = {2, -1, -13, 14, 0};
    int realpdgs[] = {2, -1, -13, 14, -2, 2};

    // FKS g = real matrix element * (1.0 - y) * xi * xi
    // real = 553971.0223624608
    double fks_g = 9.459090960047918e-14;
    double limit = FKS::QCD::SoftCollinearLimitFSR(realpdgs[4], bornpdgs[4], ps,
                                                   4, phi, alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, -1, -13, 14, 0
// realpdgs = 2, -1, -13, 14, -2, 2
TEST(QCDSoftLimit, Wj_born_3_real_25_4) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.5259719694901115;
    ps.X2 = 0.7284167167451798;
    ps.Jacobian = 11969.60113036257;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(4023.320301523392, 0, 0, 4023.320301523392);
    ps.Momenta[1].Set(4023.320301523392, 0, 0, -4023.320301523392);
    ps.Momenta[2].Set(2095.697212156931, -486.4816727167104, 319.3622302012809,
                      -2013.278458900791);
    ps.Momenta[3].Set(2999.090798808033, 2869.191188176078, -804.3072995882264,
                      339.6723612789986);
    ps.Momenta[4].Set(2951.852592081821, -2382.709515459367, 484.9450693869456,
                      1673.606097621793);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    // double born = 1.650469392179905e-08;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, -2.750782320299841e-09);
    ColorCorr.Set(0, 4, 2.475704088269857e-08);
    ColorCorr.Set(1, 0, -2.750782320299841e-09);
    ColorCorr.Set(1, 4, 2.475704088269857e-08);
    ColorCorr.Set(4, 0, 2.475704088269857e-08);
    ColorCorr.Set(4, 1, 2.475704088269857e-08);

    // radiation variables
    double y = 0.9321860003070057;
    double phi = 4.75537087283174;
    int bornpdgs[] = {2, -1, -13, 14, 0};
    int realpdgs[] = {2, -1, -13, 14, -2, 2};

    // FKS g = real matrix element * (1.0 - y) * xi * xi
    // real = 2.187371122383598e-06
    double fks_g = 7.984761511680043e-24;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 4, alphas, ColorCorr, y, phi);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, -1, -13, 14, 0
// realpdgs = 2, -1, -13, 14, -1, 1
TEST(QCDCollinearFSR, Wj_born_3_real_4) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.4660568837601637;
    ps.X2 = 0.8689664427537878;
    ps.Jacobian = 14933.69419952127;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(4136.512326705176, 0, 0, 4136.512326705176);
    ps.Momenta[1].Set(4136.512326705176, 0, 0, -4136.512326705176);
    ps.Momenta[2].Set(1604.82086137663, -1381.846595934496, -744.7698302174558,
                      333.5684073945481);
    ps.Momenta[3].Set(3086.146502118359, -1186.991288417874, 451.5627662647746,
                      -2812.728743032424);
    ps.Momenta[4].Set(3582.057289915365, 2568.837884352369, 293.2070639526812,
                      2479.160335637875);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.656133860092891e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(3.834606151250774e-10, 0);
    spin[1][2] =
        std::complex<double>(-2.80451855631159e-10, 5.224315164503088e-10);
    spin[1][3] =
        std::complex<double>(-3.641626873190428e-10, -6.178729501790344e-11);
    spin[2][1] =
        std::complex<double>(-2.80451855631159e-10, -5.224315164503088e-10);
    spin[2][2] = std::complex<double>(9.168814705855736e-10, 0);
    spin[2][3] =
        std::complex<double>(1.821581367721567e-10, 5.413291960771512e-10);
    spin[3][1] =
        std::complex<double>(-3.641626873190428e-10, 6.178729501790344e-11);
    spin[3][2] =
        std::complex<double>(1.821581367721567e-10, -5.413291960771512e-10);
    spin[3][3] = std::complex<double>(3.557917743822397e-10, 0);

    // radiation variables
    double phi = 3.835635259657122;
    double xi = 0.6944761176878219;
    int bornpdgs[] = {2, -1, -13, 14, 0};
    int realpdgs[] = {2, -1, -13, 14, -1, 1};

    // FKS g = real matrix element * (1.0 - y) * xi * xi
    // real = 4.872696786723906e-08
    double fks_g = 2.350087434214367e-16;
    double limit = FKS::QCD::CollinearLimitFSR(
        realpdgs[4], bornpdgs[4], ps, 4, xi, phi, alphas, born,
        spin); // limit = 2.350066517650932e-16
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, -1, -13, 14, 0
// realpdgs = 2, -1, -13, 14, -1, 1
TEST(QCDSoftCollinearFSR, Wj_born_3_real_4) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.9082825574177029;
    ps.X2 = 0.7744147172895239;
    ps.Jacobian = 18987.44939303779;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(5451.432545825834, 0, 0, 5451.432545825834);
    ps.Momenta[1].Set(5451.432545825834, 0, 0, -5451.432545825834);
    ps.Momenta[2].Set(4030.776147759332, 1083.689411855779, 3348.794952333068,
                      -1964.267288127);
    ps.Momenta[3].Set(3416.233379265683, -2271.667285411504, -2299.083587860817,
                      -1106.522888163272);
    ps.Momenta[4].Set(3455.855564626655, 1187.977873555725, -1049.711364472251,
                      3070.790176290271);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 2.264600628170428e-08;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(7.288122975144015e-09, 0);
    spin[1][2] =
        std::complex<double>(-5.316559478258201e-09, 7.460524468366369e-09);
    spin[1][3] =
        std::complex<double>(-4.63691132281991e-09, 2.550287342923688e-09);
    spin[2][1] =
        std::complex<double>(-5.316559478258201e-09, -7.460524468366369e-09);
    spin[2][2] = std::complex<double>(1.151534219649909e-08, 0);
    spin[2][3] =
        std::complex<double>(5.993161217969855e-09, 2.88620761586758e-09);
    spin[3][1] =
        std::complex<double>(-4.63691132281991e-09, -2.550287342923688e-09);
    spin[3][2] =
        std::complex<double>(5.993161217969855e-09, -2.88620761586758e-09);
    spin[3][3] =
        std::complex<double>(3.842541110061177e-09, -4.59544784751682e-26);

    // radiation variables
    double phi = 0.1955849808608072;
    int bornpdgs[] = {2, -1, -13, 14, 0};
    int realpdgs[] = {2, -1, -13, 14, -1, 1};

    // FKS g = real matrix element * (1.0 - y) * xi * xi
    // real = 0.001405774562732521
    double fks_g = 5.64944209364513e-22;
    double limit = FKS::QCD::SoftCollinearLimitFSR(realpdgs[4], bornpdgs[4], ps,
                                                   4, phi, alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, -1, -13, 14, 0
// realpdgs = 2, -1, -13, 14, -1, 1
TEST(QCDSoftLimit, Wj_born_3_real_4_4) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.6084302189408284;
    ps.X2 = 0.792261168114706;
    ps.Jacobian = 16780.14589096991;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(4512.871106060548, 0, 0, 4512.871106060548);
    ps.Momenta[1].Set(4512.871106060548, 0, 0, -4512.871106060548);
    ps.Momenta[2].Set(1065.645638572622, 701.6073323527532, 84.51020137348621,
                      -797.6251024519711);
    ps.Momenta[3].Set(4270.809933599997, -3665.553439739725, 2104.439615779153,
                      612.2657698156643);
    ps.Momenta[4].Set(3689.286639948478, 2963.946107386972, -2188.949817152638,
                      185.3593326363067);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 8.125048563578259e-09;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, -1.354174760596376e-09);
    ColorCorr.Set(0, 4, 1.218757284536739e-08);
    ColorCorr.Set(1, 0, -1.354174760596376e-09);
    ColorCorr.Set(1, 4, 1.218757284536739e-08);
    ColorCorr.Set(4, 0, 1.218757284536739e-08);
    ColorCorr.Set(4, 1, 1.218757284536739e-08);

    // radiation variables
    double y = 0.861024528044517;
    double phi = 5.726395709196769;
    int bornpdgs[] = {2, -1, -13, 14, 0};
    int realpdgs[] = {2, -1, -13, 14, -1, 1};

    // FKS g = real matrix element * (1.0 - y) * xi * xi
    // real = 3.184672012894753e-07
    double fks_g = 2.957888561183157e-24;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 4, alphas, ColorCorr, y, phi);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, -1, -13, 14, 0
// realpdgs = 2, -1, -13, 14, -5, 5
TEST(QCDCollinearFSR, Wj_born_3_real_13) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.3377183506057939;
    ps.X2 = 0.532812323335655;
    ps.Jacobian = 2282.756959431699;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2757.260612196078, 0, 0, 2757.260612196078);
    ps.Momenta[1].Set(2757.260612196078, 0, 0, -2757.260612196078);
    ps.Momenta[2].Set(2204.137025271802, 1663.94578545695, 1211.162299898735,
                      789.0439357463487);
    ps.Momenta[3].Set(2488.933591920535, -1695.200144047135, -937.2997228617485,
                      -1562.867917047477);
    ps.Momenta[4].Set(821.4506071998196, 31.25435859018495, -273.8625770369866,
                      773.8239813011284);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.058259782700706e-07;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(4.435996566399695e-09, 2.941086622410765e-24);
    spin[1][2] =
        std::complex<double>(-1.190773097810147e-09, 1.995299036928533e-08);
    spin[1][3] =
        std::complex<double>(-6.005919018459998e-10, 7.061524964551637e-09);
    spin[2][1] =
        std::complex<double>(-1.190773097810147e-09, -1.995299036928533e-08);
    spin[2][2] =
        std::complex<double>(9.00676452893501e-08, -1.470543311205382e-24);
    spin[2][3] =
        std::complex<double>(3.192376418013851e-08, 8.058911729507052e-10);
    spin[3][1] =
        std::complex<double>(-6.005919018459998e-10, -7.061524964551637e-09);
    spin[3][2] =
        std::complex<double>(3.192376418013851e-08, -8.058911729507052e-10);
    spin[3][3] = std::complex<double>(1.132233641432076e-08, 0);

    // radiation variables
    double phi = 0.1242209071389916;
    double xi = 0.2548511635502116;
    int bornpdgs[] = {2, -1, -13, 14, 0};
    int realpdgs[] = {2, -1, -13, 14, -5, 5};

    // FKS g = real matrix element * (1.0 - y) * xi * xi
    // real = 5.145127951656996e-05
    double fks_g = 3.341715115971917e-14;
    double limit = FKS::QCD::CollinearLimitFSR(
        realpdgs[4], bornpdgs[4], ps, 4, xi, phi, alphas, born,
        spin); // limit = 3.341100019871197e-14
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, -1, -13, 14, 0
// realpdgs = 2, -1, -13, 14, -5, 5
TEST(QCDSoftCollinearFSR, Wj_born_3_real_13) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.4190535227529644;
    ps.X2 = 0.296369894962794;
    ps.Jacobian = 3288.203976573942;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2290.683817129305, 0, 0, 2290.683817129305);
    ps.Momenta[1].Set(2290.683817129305, 0, 0, -2290.683817129305);
    ps.Momenta[2].Set(1190.13451245862, -815.230630355271, 834.6207349601417,
                      -235.0051187745864);
    ps.Momenta[3].Set(1966.960261212058, -91.07349019087042, -1908.657641322614,
                      466.5450640534651);
    ps.Momenta[4].Set(1424.272860587933, 906.3041205461415, 1074.036906362473,
                      -231.5399452788787);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 3.047681725987262e-08;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(1.666077977546348e-09, 0);
    spin[1][2] =
        std::complex<double>(-2.104346804247505e-09, 1.211981018937017e-09);
    spin[1][3] =
        std::complex<double>(-3.239928187953755e-09, 5.621977419841323e-09);
    spin[2][1] =
        std::complex<double>(-2.104346804247505e-09, -1.211981018937017e-09);
    spin[2][2] = std::complex<double>(3.539554295949057e-09, 0);
    spin[2][3] =
        std::complex<double>(8.181887422896964e-09, -4.743990891780385e-09);
    spin[3][1] =
        std::complex<double>(-3.239928187953755e-09, -5.621977419841323e-09);
    spin[3][2] =
        std::complex<double>(8.181887422896964e-09, 4.743990891780385e-09);
    spin[3][3] = std::complex<double>(2.527118498637722e-08, 0);

    // radiation variables
    double phi = 0.2331650596184928;
    int bornpdgs[] = {2, -1, -13, 14, 0};
    int realpdgs[] = {2, -1, -13, 14, -5, 5};

    // FKS g = real matrix element * (1.0 - y) * xi * xi
    // real = 0.01114006599257405
    double fks_g = 4.306692519671516e-21;
    double limit = FKS::QCD::SoftCollinearLimitFSR(realpdgs[4], bornpdgs[4], ps,
                                                   4, phi, alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, -1, -13, 14, 0
// realpdgs = 2, -1, -13, 14, -5, 5
TEST(QCDSoftLimit, Wj_born_3_real_13_4) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.2210590367059631;
    ps.X2 = 0.5254345453825877;
    ps.Jacobian = 2500.0776786712;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2215.27070596228, 0, 0, 2215.27070596228);
    ps.Momenta[1].Set(2215.27070596228, 0, 0, -2215.27070596228);
    ps.Momenta[2].Set(2129.538759670396, 474.0736521226813, 911.0130286695054,
                      -1865.541412806482);
    ps.Momenta[3].Set(1181.239184705242, -581.6406645587474, -197.9780458398671,
                      1008.873055532851);
    ps.Momenta[4].Set(1119.763467548921, 107.5670124360661, -713.0349828296382,
                      856.6683572736307);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 2.469949118028816e-07;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, -4.116581863381359e-08);
    ColorCorr.Set(0, 4, 3.704923677043223e-07);
    ColorCorr.Set(1, 0, -4.116581863381359e-08);
    ColorCorr.Set(1, 4, 3.704923677043223e-07);
    ColorCorr.Set(4, 0, 3.704923677043223e-07);
    ColorCorr.Set(4, 1, 3.704923677043223e-07);

    // radiation variables
    double y = -0.3161686048191328;
    double phi = 2.870233676175648;
    int bornpdgs[] = {2, -1, -13, 14, 0};
    int realpdgs[] = {2, -1, -13, 14, -5, 5};

    // FKS g = real matrix element * (1.0 - y) * xi * xi
    // real = 1.109649363968894e-05
    double fks_g = 3.731610234532009e-22;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 4, alphas, ColorCorr, y, phi);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, -1, -13, 14, 0
// realpdgs = 2, -1, -13, 14, -4, 4
TEST(QCDCollinearFSR, Wj_born_3_real_13_n1) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.5881217415065016;
    ps.X2 = 0.5851499257200174;
    ps.Jacobian = 13629.45607979168;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3813.120686436137, 0, 0, 3813.120686436137);
    ps.Momenta[1].Set(3813.120686436137, 0, 0, -3813.120686436137);
    ps.Momenta[2].Set(1528.27950050781, 1291.472343481101, 276.7614125218596,
                      -768.8566434874195);
    ps.Momenta[3].Set(2551.481424764558, 949.6583133911068, -1837.995263354962,
                      -1493.311742607838);
    ps.Momenta[4].Set(3546.480447599906, -2241.130656872208, 1561.233850833102,
                      2262.168386095258);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.268681374882911e-08;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(6.549454913482073e-09, 0);
    spin[1][2] =
        std::complex<double>(4.031077912588189e-09, 2.629633827797739e-09);
    spin[1][3] =
        std::complex<double>(3.706500785470525e-09, -1.814839855639632e-09);
    spin[2][1] =
        std::complex<double>(4.031077912588189e-09, -2.629633827797739e-09);
    spin[2][2] = std::complex<double>(3.536868871021651e-09, 0);
    spin[2][3] =
        std::complex<double>(1.552622211650454e-09, -2.605178740915161e-09);
    spin[3][1] =
        std::complex<double>(3.706500785470525e-09, 1.814839855639632e-09);
    spin[3][2] =
        std::complex<double>(1.552622211650454e-09, 2.605178740915161e-09);
    spin[3][3] =
        std::complex<double>(2.600489964325383e-09, 1.838179139006728e-25);

    // radiation variables
    double phi = 1.305040092370744;
    double xi = 0.8308848363637686;
    int bornpdgs[] = {2, -1, -13, 14, 0};
    int realpdgs[] = {2, -1, -13, 14, -4, 4};

    // FKS g = real matrix element * (1.0 - y) * xi * xi
    // real = 5.429409566509146e-07
    double fks_g = 3.748299390849632e-15;
    double limit = FKS::QCD::CollinearLimitFSR(
        realpdgs[4], bornpdgs[4], ps, 4, xi, phi, alphas, born,
        spin); // limit = 3.748243604955768e-15
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, -1, -13, 14, 0
// realpdgs = 2, -1, -13, 14, -4, 4
TEST(QCDSoftCollinearFSR, Wj_born_3_real_13_n1) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.1157478839406245;
    ps.X2 = 0.1589145687876723;
    ps.Jacobian = 5.76156974313868;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(881.5597308041824, 0, 0, 881.5597308041824);
    ps.Momenta[1].Set(881.5597308041824, 0, 0, -881.5597308041824);
    ps.Momenta[2].Set(880.8611323549894, -403.7466105790157, 742.3500173546864,
                      248.6392178904355);
    ps.Momenta[3].Set(875.7736473668966, 401.4591491473317, -736.757766892409,
                      -250.9940754518491);
    ps.Momenta[4].Set(6.484681886479051, 2.287461431683957, -5.592250462277333,
                      2.354857561413608);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 0.000422470843301392;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(8.052376671603524e-06, 3.011672701348623e-21);
    spin[1][2] =
        std::complex<double>(-1.957172522298483e-05, 9.083042554943614e-07);
    spin[1][3] =
        std::complex<double>(-5.430030783633779e-05, 2.157015768557819e-06);
    spin[2][1] =
        std::complex<double>(-1.957172522298483e-05, -9.083042554943614e-07);
    spin[2][2] = std::complex<double>(4.767256432232007e-05, 0);
    spin[2][3] =
        std::complex<double>(0.000132223065874788, 8.823085466920937e-07);
    spin[3][1] =
        std::complex<double>(-5.430030783633779e-05, -2.157015768557819e-06);
    spin[3][2] =
        std::complex<double>(0.000132223065874788, -8.823085466920937e-07);
    spin[3][3] = std::complex<double>(0.0003667459023074685, 0);

    // radiation variables
    double phi = 4.531816487336044;
    int bornpdgs[] = {2, -1, -13, 14, 0};
    int realpdgs[] = {2, -1, -13, 14, -4, 4};

    // FKS g = real matrix element * (1.0 - y) * xi * xi
    // real = 7450854.121422853
    double fks_g = 4.031623106099267e-16;
    double limit = FKS::QCD::SoftCollinearLimitFSR(realpdgs[4], bornpdgs[4], ps,
                                                   4, phi, alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, -1, -13, 14, 0
// realpdgs = 2, -1, -13, 14, -4, 4
TEST(QCDSoftLimit, Wj_born_3_real_13_4_n1) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.65798240964833;
    ps.X2 = 0.35877448545466;
    ps.Jacobian = 3916.707797483481;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3158.139237656633, 0, 0, 3158.139237656633);
    ps.Momenta[1].Set(3158.139237656633, 0, 0, -3158.139237656633);
    ps.Momenta[2].Set(3090.927385314439, -581.0004639643798, 2657.326620752739,
                      1467.952925947816);
    ps.Momenta[3].Set(1994.828770307596, 99.73436962828421, -1544.307335967821,
                      -1258.773105224979);
    ps.Momenta[4].Set(1230.522319691231, 481.2660943360955, -1113.019284784918,
                      -209.1798207228369);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 5.323128382333341e-09;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, -8.871880637222234e-10);
    ColorCorr.Set(0, 4, 7.984692573500011e-09);
    ColorCorr.Set(1, 0, -8.871880637222234e-10);
    ColorCorr.Set(1, 4, 7.984692573500011e-09);
    ColorCorr.Set(4, 0, 7.984692573500011e-09);
    ColorCorr.Set(4, 1, 7.984692573500011e-09);

    // radiation variables
    double y = 0.2020824686384302;
    double phi = 3.230129730183163;
    int bornpdgs[] = {2, -1, -13, 14, 0};
    int realpdgs[] = {2, -1, -13, 14, -4, 4};

    // FKS g = real matrix element * (1.0 - y) * xi * xi
    // real = 3.266568834310361e-07
    double fks_g = 3.957001961125758e-24;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 4, alphas, ColorCorr, y, phi);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, -1, -13, 14, 0
// realpdgs = 2, -1, -13, 14, -3, 3
TEST(QCDCollinearFSR, Wj_born_3_real_13_n2) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.6302145965532411;
    ps.X2 = 0.1808932032536248;
    ps.Jacobian = 3502.463651054462;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2194.667387737969, 0, 0, 2194.667387737969);
    ps.Momenta[1].Set(2194.667387737969, 0, 0, -2194.667387737969);
    ps.Momenta[2].Set(1350.25569845314, -386.2549167041066, 1126.621547062445,
                      636.1772396287882);
    ps.Momenta[3].Set(1455.628462004993, 1368.056780775653, 106.5622013234637,
                      -485.7153088197632);
    ps.Momenta[4].Set(1583.450615017805, -981.8018640715461, -1233.183748385909,
                      -150.4619308090251);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 6.542579119910052e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(2.945755578158335e-09, 0);
    spin[1][2] =
        std::complex<double>(-2.051258779869677e-09, 9.231684190900513e-11);
    spin[1][3] =
        std::complex<double>(-2.409708055162517e-09, -7.566274640526394e-10);
    spin[2][1] =
        std::complex<double>(-2.051258779869677e-09, -9.231684190900513e-11);
    spin[2][2] =
        std::complex<double>(1.431274547200694e-09, 9.19089569503364e-26);
    spin[2][3] =
        std::complex<double>(1.654273485289793e-09, 6.02390564736943e-10);
    spin[3][1] =
        std::complex<double>(-2.409708055162517e-09, 7.566274640526394e-10);
    spin[3][2] =
        std::complex<double>(1.654273485289793e-09, -6.02390564736943e-10);
    spin[3][3] = std::complex<double>(2.165548994551022e-09, 0);

    // radiation variables
    double phi = 5.368582249099881;
    double xi = 0.4486798985911715;
    int bornpdgs[] = {2, -1, -13, 14, 0};
    int realpdgs[] = {2, -1, -13, 14, -3, 3};

    // FKS g = real matrix element * (1.0 - y) * xi * xi
    // real = 7.379073184363611e-07
    double fks_g = 1.485508174154791e-15;
    double limit = FKS::QCD::CollinearLimitFSR(
        realpdgs[4], bornpdgs[4], ps, 4, xi, phi, alphas, born,
        spin); // limit = 1.485502556348801e-15
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, -1, -13, 14, 0
// realpdgs = 2, -1, -13, 14, -3, 3
TEST(QCDSoftCollinearFSR, Wj_born_3_real_13_n2) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.1814570983640813;
    ps.X2 = 0.3079279642937109;
    ps.Jacobian = 2323.800608574036;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1536.472894252962, 0, 0, 1536.472894252962);
    ps.Momenta[1].Set(1536.472894252962, 0, 0, -1536.472894252962);
    ps.Momenta[2].Set(1509.675469400844, -220.7397326071563, 116.5831445142879,
                      1488.892999437615);
    ps.Momenta[3].Set(62.6405608830928, -52.31114395332521, -30.61188909726854,
                      -15.82075636432899);
    ps.Momenta[4].Set(1500.629758221989, 273.0508765604813, -85.97125541701929,
                      -1473.072243073285);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 9.254516066150693e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(1.389066456321072e-09, -2.757268708510092e-25);
    spin[1][2] =
        std::complex<double>(-1.225470406331042e-09, 3.04696424310509e-09);
    spin[1][3] =
        std::complex<double>(3.290001865738801e-10, -1.778265407024445e-10);
    spin[2][1] =
        std::complex<double>(-1.225470406331042e-09, -3.04696424310509e-09);
    spin[2][2] =
        std::complex<double>(7.76476083377619e-09, 1.838179139006728e-25);
    spin[2][3] =
        std::complex<double>(-6.803210163334829e-10, -5.647898542250231e-10);
    spin[3][1] =
        std::complex<double>(3.290001865738801e-10, 1.778265407024445e-10);
    spin[3][2] =
        std::complex<double>(-6.803210163334829e-10, 5.647898542250231e-10);
    spin[3][3] =
        std::complex<double>(1.006887760534314e-10, 1.148861961879205e-26);

    // radiation variables
    double phi = 3.515865250032069;
    int bornpdgs[] = {2, -1, -13, 14, 0};
    int realpdgs[] = {2, -1, -13, 14, -3, 3};

    // FKS g = real matrix element * (1.0 - y) * xi * xi
    // real = 0.00304746157151132
    double fks_g = 2.906936471483654e-21;
    double limit = FKS::QCD::SoftCollinearLimitFSR(realpdgs[4], bornpdgs[4], ps,
                                                   4, phi, alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, -1, -13, 14, 0
// realpdgs = 2, -1, -13, 14, -3, 3
TEST(QCDSoftLimit, Wj_born_3_real_13_4_n2) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.546866177274719;
    ps.X2 = 0.4391279598723585;
    ps.Jacobian = 5626.752210117113;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3185.293340444623, 0, 0, 3185.293340444623);
    ps.Momenta[1].Set(3185.293340444623, 0, 0, -3185.293340444623);
    ps.Momenta[2].Set(1728.695510037313, 1458.123953085152, -127.4712573039479,
                      -919.7900751931625);
    ps.Momenta[3].Set(2889.189692485466, -1949.805754201536, -1412.521441425507,
                      1597.015271548788);
    ps.Momenta[4].Set(1752.701478366466, 491.6818011163838, 1539.992698729455,
                      -677.2251963556257);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 3.415354236127075e-08;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, -5.692257060211792e-09);
    ColorCorr.Set(0, 4, 5.123031354190612e-08);
    ColorCorr.Set(1, 0, -5.692257060211792e-09);
    ColorCorr.Set(1, 4, 5.123031354190612e-08);
    ColorCorr.Set(4, 0, 5.123031354190612e-08);
    ColorCorr.Set(4, 1, 5.123031354190612e-08);

    // radiation variables
    double y = 0.149126683298384;
    double phi = 5.74798468761835;
    int bornpdgs[] = {2, -1, -13, 14, 0};
    int realpdgs[] = {2, -1, -13, 14, -3, 3};

    // FKS g = real matrix element * (1.0 - y) * xi * xi
    // real = 9.687614198352607e-07
    double fks_g = 2.495736773625638e-23;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 4, alphas, ColorCorr, y, phi);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, -1, -13, 14, 0
// realpdgs = 2, 21, -13, 14, 21, 1
TEST(QCDCollinearISR1, Wj_born_3_real_36) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.4156099401429323;
    ps.X2 = 0.6541588924411883;
    ps.Jacobian = 7459.059113145995;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3389.205826746692, 0, 0, 3389.205826746692);
    ps.Momenta[1].Set(3389.205826746692, 0, 0, -3389.205826746692);
    ps.Momenta[2].Set(2168.732560059794, 650.1187987636935, 57.84135416594256,
                      -2068.187816013287);
    ps.Momenta[3].Set(2426.015663500963, 1023.567086748272, -1408.315193737033,
                      1689.529737402422);
    ps.Momenta[4].Set(2183.663429932627, -1673.685885511966, 1350.47383957109,
                      378.6580786108654);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 2.770883025361118e-08;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(8.121794527129341e-10, 0);
    spin[1][2] =
        std::complex<double>(-2.992373768055832e-10, 6.98075615086533e-11);
    spin[1][3] =
        std::complex<double>(4.657092071495251e-09, -2.489667881048109e-10);
    spin[2][1] =
        std::complex<double>(-2.992373768055832e-10, -6.98075615086533e-11);
    spin[2][2] = std::complex<double>(1.162502978940086e-10, 0);
    spin[2][3] =
        std::complex<double>(-1.737246337219591e-09, -3.085525887303465e-10);
    spin[3][1] =
        std::complex<double>(4.657092071495251e-09, 2.489667881048109e-10);
    spin[3][2] =
        std::complex<double>(-1.737246337219591e-09, 3.085525887303465e-10);
    spin[3][3] = std::complex<double>(2.678040050300424e-08, 0);

    // radiation variables
    double phi = 5.817831403638403;
    double xi = 0.05884105830074772;
    int bornpdgs[] = {2, -1, -13, 14, 0};
    int realpdgs[] = {2, 21, -13, 14, 21, 1};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 2.905298401524096e-14
    double fks_g = 2.01178558075263e-24;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, 1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, -1, -13, 14, 0
// realpdgs = 2, 21, -13, 14, 21, 1
TEST(QCDSoftCollinearISR1, Wj_born_3_real_36) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.9212841669384311;
    ps.X2 = 0.6059274985997618;
    ps.Jacobian = 3300.653715176364;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(4856.467554214795, 0, 0, 4856.467554214795);
    ps.Momenta[1].Set(4856.467554214795, 0, 0, -4856.467554214795);
    ps.Momenta[2].Set(4799.916273292397, 2273.223304421934, 2452.821895765281,
                      3443.155062802808);
    ps.Momenta[3].Set(4238.678607876011, -2210.32058729132, -2303.092381277186,
                      -2788.663609251723);
    ps.Momenta[4].Set(674.3402272611841, -62.90271713061387, -149.7295144880947,
                      -654.4914535510841);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.104048920051832e-07;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(1.458115233470621e-08, 0);
    spin[1][2] =
        std::complex<double>(3.593498025040806e-08, 3.556254108365529e-09);
    spin[1][3] =
        std::complex<double>(-9.622312426864621e-09, -8.135724281085699e-10);
    spin[2][1] =
        std::complex<double>(3.593498025040806e-08, -3.556254108365529e-09);
    spin[2][2] =
        std::complex<double>(8.942844289314239e-08, -2.941086622410765e-24);
    spin[2][3] =
        std::complex<double>(-2.391243636366382e-08, 3.417891020721178e-10);
    spin[3][1] =
        std::complex<double>(-9.622312426864621e-09, 8.135724281085699e-10);
    spin[3][2] =
        std::complex<double>(-2.391243636366382e-08, -3.417891020721178e-10);
    spin[3][3] =
        std::complex<double>(6.39529677733462e-09, -9.19089569503364e-26);

    // radiation variables
    double phi = 2.457722797226403;
    int bornpdgs[] = {2, -1, -13, 14, 0};
    int realpdgs[] = {2, 21, -13, 14, 21, 1};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 3.498034982927538e-08
    double fks_g = 4.334894538381191e-28;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, 1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                   << ", fks_g = " << fks_g;
}

// bornpdgs = 2, -1, -13, 14, 0
// realpdgs = 2, 21, -13, 14, 21, 1
TEST(QCDSoftLimit, Wj_born_3_real_36_0) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.865675992138069;
    ps.X2 = 0.01505738249782951;
    ps.Jacobian = 42.02365381820108;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(742.1057296714955, 0, 0, 742.1057296714955);
    ps.Momenta[1].Set(742.1057296714955, 0, 0, -742.1057296714955);
    ps.Momenta[2].Set(721.8146266363908, 611.2001120384724, 328.2844793098095,
                      199.1985916483303);
    ps.Momenta[3].Set(706.2108964850942, -592.5092444285054, -297.9684886397867,
                      -242.6549100240159);
    ps.Momenta[4].Set(56.18593622150621, -18.69086760996692, -30.31599067002272,
                      43.45631837568558);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    // double born = 1.140699613452073e-05;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, -1.901166022420121e-06);
    ColorCorr.Set(0, 4, 1.711049420178109e-05);
    ColorCorr.Set(1, 0, -1.901166022420121e-06);
    ColorCorr.Set(1, 4, 1.711049420178109e-05);
    ColorCorr.Set(4, 0, 1.711049420178109e-05);
    ColorCorr.Set(4, 1, 1.711049420178109e-05);

    // radiation variables
    double y = -0.9934957403925395;
    double phi = 3.309165884967592;
    int bornpdgs[] = {2, -1, -13, 14, 0};
    int realpdgs[] = {2, 21, -13, 14, 21, 1};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 0.2389466587674635
    double fks_g = 3.015057436857883e-19;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 0, alphas, ColorCorr, y, phi);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, -1, -13, 14, 0
// realpdgs = 2, 21, -13, 14, 21, 1
TEST(QCDCollinearISR2, Wj_born_3_real_36) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.7897562615032676;
    ps.X2 = 0.654018747873458;
    ps.Jacobian = 9987.781389479822;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(4671.485385164893, 0, 0, 4671.485385164893);
    ps.Momenta[1].Set(4671.485385164893, 0, 0, -4671.485385164893);
    ps.Momenta[2].Set(4651.572316316498, 43.18952091604409, 3966.676834024006,
                      -2429.142765179181);
    ps.Momenta[3].Set(2570.04209627395, 204.5600614760369, -2282.612111095846,
                      1163.165297001112);
    ps.Momenta[4].Set(2121.356357739337, -247.749582392081, -1684.064722928161,
                      1265.977468178069);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 2.784438233586024e-08;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(2.488437993457381e-09, -9.19089569503364e-26);
    spin[1][2] =
        std::complex<double>(6.129869492654234e-10, 4.694045027932738e-09);
    spin[1][3] =
        std::complex<double>(1.302407990668116e-09, 6.244246709030648e-09);
    spin[2][1] =
        std::complex<double>(6.129869492654234e-10, -4.694045027932738e-09);
    spin[2][2] = std::complex<double>(9.005573690463588e-09, 0);
    spin[2][3] =
        std::complex<double>(1.209961204489418e-08, -9.186164245668733e-10);
    spin[3][1] =
        std::complex<double>(1.302407990668116e-09, -6.244246709030648e-09);
    spin[3][2] =
        std::complex<double>(1.209961204489418e-08, 9.186164245668733e-10);
    spin[3][3] = std::complex<double>(1.635037065193927e-08, 0);

    // radiation variables
    double phi = 2.937933679372612;
    double xi = 0.05084917383206267;
    int bornpdgs[] = {2, -1, -13, 14, 0};
    int realpdgs[] = {2, 21, -13, 14, 21, 1};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 1.68081725626762e-06
    double fks_g = 8.691971544726295e-17;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, -1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 8.691975957307947e-17
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, -1, -13, 14, 0
// realpdgs = 2, 21, -13, 14, 21, 1
TEST(QCDSoftCollinearISR2, Wj_born_3_real_36) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.1008759035480331;
    ps.X2 = 0.5042142095686017;
    ps.Jacobian = 306.6960045936707;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1465.934668672714, 0, 0, 1465.934668672714);
    ps.Momenta[1].Set(1465.934668672714, 0, 0, -1465.934668672714);
    ps.Momenta[2].Set(1454.120541700396, 490.8126442387927, 80.72945317890498,
                      1366.401205151625);
    ps.Momenta[3].Set(1270.16516718696, -508.8113254911511, -43.43989707145674,
                      -1162.989063718759);
    ps.Momenta[4].Set(207.583628458072, 17.99868125235842, -37.28955610744823,
                      -203.4121414328662);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 9.11114772328275e-08;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(2.253832264816311e-08, 0);
    spin[1][2] =
        std::complex<double>(-3.724507767129041e-08, 8.825085807710947e-09);
    spin[1][3] =
        std::complex<double>(8.82205204771028e-09, -1.617816567199814e-09);
    spin[2][1] =
        std::complex<double>(-3.724507767129041e-08, -8.825085807710947e-09);
    spin[2][2] = std::complex<double>(6.500385912140277e-08, 0);
    spin[2][3] =
        std::complex<double>(-1.521210735670047e-08, -7.808772148939036e-10);
    spin[3][1] =
        std::complex<double>(8.82205204771028e-09, 1.617816567199814e-09);
    spin[3][2] =
        std::complex<double>(-1.521210735670047e-08, 7.808772148939036e-10);
    spin[3][3] =
        std::complex<double>(3.569295463261627e-09, 4.59544784751682e-26);

    // radiation variables
    double phi = 0.238523706891768;
    int bornpdgs[] = {2, -1, -13, 14, 0};
    int realpdgs[] = {2, 21, -13, 14, 21, 1};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 0.06342246060528269
    double fks_g = 3.117893979389817e-20;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, -1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, -1, -13, 14, 0
// realpdgs = 21, -1, -13, 14, 21, -2
TEST(QCDCollinearISR1, Wj_born_3_real_11) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.1902569095273172;
    ps.X2 = 0.2766655764547252;
    ps.Jacobian = 1238.799680780106;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1491.286679830456, 0, 0, 1491.286679830456);
    ps.Momenta[1].Set(1491.286679830456, 0, 0, -1491.286679830456);
    ps.Momenta[2].Set(1475.481056059834, -860.3964115168108, -1186.210374956355,
                      -172.2420047094807);
    ps.Momenta[3].Set(682.8791475771393, 402.5662586556916, 546.6149498069358,
                      -74.00293396346001);
    ps.Momenta[4].Set(824.2131560239394, 457.8301528611191, 639.5954251494195,
                      246.2449386729407);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 7.542875504463212e-08;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(1.065416131877432e-08, 7.352716556026912e-25);
    spin[1][2] =
        std::complex<double>(-3.596830471925205e-09, 8.560200837260397e-09);
    spin[1][3] =
        std::complex<double>(-1.046632675656004e-08, -2.223422468448739e-08);
    spin[2][1] =
        std::complex<double>(-3.596830471925205e-09, -8.560200837260397e-09);
    spin[2][2] = std::complex<double>(8.092070810499188e-09, 0);
    spin[2][3] =
        std::complex<double>(-1.4330909884367e-08, 1.591552735648358e-08);
    spin[3][1] =
        std::complex<double>(-1.046632675656004e-08, 2.223422468448739e-08);
    spin[3][2] =
        std::complex<double>(-1.4330909884367e-08, -1.591552735648358e-08);
    spin[3][3] = std::complex<double>(5.668252291535861e-08, 0);

    // radiation variables
    double phi = 2.47630618431996;
    double xi = 0.436400513455364;
    int bornpdgs[] = {2, -1, -13, 14, 0};
    int realpdgs[] = {21, -1, -13, 14, 21, -2};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 2.927743944618705e-06
    double fks_g = 1.115150780361718e-14;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, 1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 1.115147342610844e-14
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, -1, -13, 14, 0
// realpdgs = 21, -1, -13, 14, 21, -2
TEST(QCDSoftCollinearISR1, Wj_born_3_real_11) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.7219923085429976;
    ps.X2 = 0.7247018112167467;
    ps.Jacobian = 15459.31475740804;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(4701.747642974852, 0, 0, 4701.747642974852);
    ps.Momenta[1].Set(4701.747642974852, 0, 0, -4701.747642974852);
    ps.Momenta[2].Set(4105.566676325526, -851.5859431995074, 2907.111920659326,
                      2771.133233148501);
    ps.Momenta[3].Set(2035.57884827286, -853.6377408907234, -1836.624650675997,
                      -204.1904684164158);
    ps.Momenta[4].Set(3262.349761351321, 1705.223684090231, -1070.487269983329,
                      -2566.942764732085);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 2.808038371907531e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(1.946743779127105e-09, 2.29772392375841e-26);
    spin[1][2] =
        std::complex<double>(1.429663573902393e-10, 3.384912507419753e-10);
    spin[1][3] =
        std::complex<double>(1.233603638112759e-09, -1.411603639545265e-10);
    spin[2][1] =
        std::complex<double>(1.429663573902393e-10, -3.384912507419753e-10);
    spin[2][2] =
        std::complex<double>(6.9354635993669e-11, -3.446585885637615e-26);
    spin[2][3] =
        std::complex<double>(6.604991978544941e-11, -2.248602133062303e-10);
    spin[3][1] =
        std::complex<double>(1.233603638112759e-09, 1.411603639545265e-10);
    spin[3][2] =
        std::complex<double>(6.604991978544941e-11, 2.248602133062303e-10);
    spin[3][3] = std::complex<double>(7.91939956786757e-10, 0);

    // radiation variables
    double phi = 0.1645029168349104;
    int bornpdgs[] = {2, -1, -13, 14, 0};
    int realpdgs[] = {21, -1, -13, 14, 21, -2};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 0.0003388463867299689
    double fks_g = 5.237772537373555e-23;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, 1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, -1, -13, 14, 0
// realpdgs = 21, -1, -13, 14, 21, -2
TEST(QCDSoftLimit, Wj_born_3_real_11_0) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.8492578810742337;
    ps.X2 = 0.04181605797614196;
    ps.Jacobian = 1427.335558487138;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1224.911449635908, 0, 0, 1224.911449635908);
    ps.Momenta[1].Set(1224.911449635908, 0, 0, -1224.911449635908);
    ps.Momenta[2].Set(1218.831306693684, 1201.292503767818, -152.3460958011003,
                      138.6965812918319);
    ps.Momenta[3].Set(74.82347922790504, -67.10171684036197, 7.132816652295347,
                      32.32700987486079);
    ps.Momenta[4].Set(1156.168113350227, -1134.190786927456, 145.2132791488049,
                      -171.0235911666927);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    // double born = 1.668228470885324e-07;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, -2.780380784808873e-08);
    ColorCorr.Set(0, 4, 2.502342706327986e-07);
    ColorCorr.Set(1, 0, -2.780380784808873e-08);
    ColorCorr.Set(1, 4, 2.502342706327986e-07);
    ColorCorr.Set(4, 0, 2.502342706327986e-07);
    ColorCorr.Set(4, 1, 2.502342706327986e-07);

    // radiation variables
    double y = 0.8253119627582848;
    double phi = 4.715005871866179;
    int bornpdgs[] = {2, -1, -13, 14, 0};
    int realpdgs[] = {21, -1, -13, 14, 21, -2};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 0.0003685439192660372
    double fks_g = 3.159750810315479e-22;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 0, alphas, ColorCorr, y, phi);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, -1, -13, 14, 0
// realpdgs = 21, -1, -13, 14, 21, -2
TEST(QCDCollinearISR2, Wj_born_3_real_11) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.7138128889650979;
    ps.X2 = 0.9711323391765916;
    ps.Jacobian = 28655.92544834277;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(5411.837625071725, 0, 0, 5411.837625071725);
    ps.Momenta[1].Set(5411.837625071725, 0, 0, -5411.837625071725);
    ps.Momenta[2].Set(919.8598489074177, -860.2499984403164, 283.792905342672,
                      -159.9176934944287);
    ps.Momenta[3].Set(4650.066921480655, -3742.440136800921, -815.6733265212483,
                      2636.653413156882);
    ps.Momenta[4].Set(5253.74847975538, 4602.690135241235, 531.8804211785761,
                      -2476.735719662452);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 4.447182581424717e-08;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(5.748424060368606e-09, -9.19089569503364e-26);
    spin[1][2] =
        std::complex<double>(6.112325679602035e-10, 1.003204109240439e-08);
    spin[1][3] =
        std::complex<double>(1.081395852569079e-08, 2.154386598113115e-09);
    spin[2][1] =
        std::complex<double>(6.112325679602035e-10, -1.003204109240439e-08);
    spin[2][2] = std::complex<double>(1.75727212660348e-08, 0);
    spin[2][3] =
        std::complex<double>(4.909647970344672e-09, -1.86432392466312e-08);
    spin[3][1] =
        std::complex<double>(1.081395852569079e-08, -2.154386598113115e-09);
    spin[3][2] =
        std::complex<double>(4.909647970344672e-09, 1.86432392466312e-08);
    spin[3][3] = std::complex<double>(2.115068048784376e-08, 0);

    // radiation variables
    double phi = 5.895190952510198;
    double xi = 0.00355673523031769;
    int bornpdgs[] = {2, -1, -13, 14, 0};
    int realpdgs[] = {21, -1, -13, 14, 21, -2};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 3.591538337733647e-13
    double fks_g = 9.086854530117874e-26;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, -1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, -1, -13, 14, 0
// realpdgs = 21, -1, -13, 14, 21, -2
TEST(QCDSoftCollinearISR2, Wj_born_3_real_11) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.1288037846300476;
    ps.X2 = 0.5926116736081681;
    ps.Jacobian = 2852.076091270048;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1795.819858564513, 0, 0, 1795.819858564513);
    ps.Momenta[1].Set(1795.819858564513, 0, 0, -1795.819858564513);
    ps.Momenta[2].Set(1468.573486688077, -1373.791742082022, -111.1627887365797,
                      -506.9981948594357);
    ps.Momenta[3].Set(547.2777017206961, -40.27449933023058, -509.9876786180147,
                      194.4310036033216);
    ps.Momenta[4].Set(1575.788528720254, 1414.066241412253, 621.1504673545945,
                      312.5671912561142);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 7.143002851882041e-08;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(7.987958794234536e-09, 0);
    spin[1][2] =
        std::complex<double>(-9.968593421912583e-09, 5.33388300716509e-09);
    spin[1][3] =
        std::complex<double>(-1.632770984478817e-08, -1.059978147226774e-08);
    spin[2][1] =
        std::complex<double>(-9.968593421912583e-09, -5.33388300716509e-09);
    spin[2][2] = std::complex<double>(1.600198073602761e-08, 0);
    spin[2][3] =
        std::complex<double>(1.329830428716424e-08, 2.41306960777416e-08);
    spin[3][1] =
        std::complex<double>(-1.632770984478817e-08, 1.059978147226774e-08);
    spin[3][2] =
        std::complex<double>(1.329830428716424e-08, -2.41306960777416e-08);
    spin[3][3] = std::complex<double>(4.744008898855826e-08, 0);

    // radiation variables
    double phi = 0.1067914552613785;
    int bornpdgs[] = {2, -1, -13, 14, 0};
    int realpdgs[] = {21, -1, -13, 14, 21, -2};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 5.152371603003015e-08
    double fks_g = 1.710229773529958e-26;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, -1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, -1, -13, 14, -2
// realpdgs = 2, -1, -13, 14, -2, 2
TEST(QCDCollinearISR1, Wj_born_1_real_25) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.9734254526432728;
    ps.X2 = 0.7623585453721233;
    ps.Jacobian = 4361.616641055772;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(5599.436731623105, 0, 0, 5599.436731623105);
    ps.Momenta[1].Set(5599.436731623105, 0, 0, -5599.436731623105);
    ps.Momenta[2].Set(5222.049507857027, 248.4348906194253, 2170.597140209088,
                      -4743.056949115081);
    ps.Momenta[3].Set(5203.960547426938, -760.4259737213775, -1610.526671654505,
                      4889.699536531404);
    ps.Momenta[4].Set(772.8634079622444, 511.9910831019521, -560.0704685545834,
                      -146.6425874163231);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 2.071585714247723e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(1.106204613089639e-09, 0);
    spin[1][2] =
        std::complex<double>(1.451662687696058e-10, -1.023149931291665e-09);
    spin[2][1] =
        std::complex<double>(1.451662687696058e-10, 1.023149931291665e-09);
    spin[2][2] = std::complex<double>(9.653811011580838e-10, 0);

    // radiation variables
    double phi = 3.54593237037814;
    double xi = 0.01140313936010562;
    int bornpdgs[] = {0, -1, -13, 14, -2};
    int realpdgs[] = {2, -1, -13, 14, -2, 2};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 1.163075626731978e-06
    double fks_g = 3.024731395494817e-18;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, 1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 3.024146735537628e-18
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, -1, -13, 14, -2
// realpdgs = 2, -1, -13, 14, -2, 2
TEST(QCDSoftCollinearISR1, Wj_born_1_real_25) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.576261904360134;
    ps.X2 = 0.4682228102207588;
    ps.Jacobian = 6156.513825134336;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3376.366598866706, 0, 0, 3376.366598866706);
    ps.Momenta[1].Set(3376.366598866706, 0, 0, -3376.366598866706);
    ps.Momenta[2].Set(2331.092300103302, 2088.270030667285, 607.3969060872988,
                      -839.1594539140563);
    ps.Momenta[3].Set(2612.448058323814, -1592.319803195826, -2068.389527741685,
                      -105.6752729771761);
    ps.Momenta[4].Set(1809.192839306296, -495.9502274714592, 1460.992621654386,
                      944.8347268912322);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.570177951405508e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(5.599218718100605e-10, 0);
    spin[1][2] =
        std::complex<double>(6.604620968360021e-10, -3.597975732791537e-10);
    spin[2][1] =
        std::complex<double>(6.604620968360021e-10, 3.597975732791537e-10);
    spin[2][2] = std::complex<double>(1.010256079595447e-09, 0);

    // radiation variables
    double phi = 1.133787925250638;
    int bornpdgs[] = {0, -1, -13, 14, -2};
    int realpdgs[] = {2, -1, -13, 14, -2, 2};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 0.0006419192650374367
    double fks_g = 2.305183760882653e-22;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, 1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, -1, -13, 14, -2
// realpdgs = 2, -1, -13, 14, -2, 2
TEST(QCDSoftLimit, Wj_born_1_real_25_0) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.8930769661806472;
    ps.X2 = 0.7308787987043388;
    ps.Jacobian = 19392.6936406977;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(5251.465091109193, 0, 0, 5251.465091109193);
    ps.Momenta[1].Set(5251.465091109193, 0, 0, -5251.465091109193);
    ps.Momenta[2].Set(3140.196906464638, -2967.868082764888, -448.7451924519471,
                      922.6176927237054);
    ps.Momenta[3].Set(3698.718214553306, 1822.183861443121, 2729.077972897383,
                      1706.54499606982);
    ps.Momenta[4].Set(3664.015061200442, 1145.684221321767, -2280.332780445436,
                      -2629.162688793525);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    // double born = 3.967602122273385e-10;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, 5.951403183410078e-10);
    ColorCorr.Set(0, 4, 5.951403183410078e-10);
    ColorCorr.Set(1, 0, 5.951403183410078e-10);
    ColorCorr.Set(1, 4, -6.612670203788975e-11);
    ColorCorr.Set(4, 0, 5.951403183410078e-10);
    ColorCorr.Set(4, 1, -6.612670203788975e-11);

    // radiation variables
    double y = 0.00537805788361112;
    double phi = 5.640160255134591;
    int bornpdgs[] = {0, -1, -13, 14, -2};
    int realpdgs[] = {2, -1, -13, 14, -2, 2};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 4.558915255867973e-07
    double fks_g = 1.850087672999746e-24;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 0, alphas, ColorCorr, y, phi);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, -1, -13, 14, -2
// realpdgs = 2, -1, -13, 14, -2, 2
TEST(QCDCollinearISR2, Wj_born_1_real_25) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.7248083484442525;
    ps.X2 = 0.2612821929102509;
    ps.Jacobian = 4394.69471637505;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2828.654184761633, 0, 0, 2828.654184761633);
    ps.Momenta[1].Set(2828.654184761633, 0, 0, -2828.654184761633);
    ps.Momenta[2].Set(1529.871026874466, -111.2910198955404, -35.60443383912468,
                      -1525.402239427863);
    ps.Momenta[3].Set(2585.919870087374, -1202.758272169533, 368.6980431342324,
                      2259.295435800605);
    ps.Momenta[4].Set(1541.517472561425, 1314.049292065074, -333.0936092951077,
                      -733.8931963727423);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 2.847627012466848e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(1.165956693082959e-09, 0);
    spin[1][2] =
        std::complex<double>(1.153549231359039e-10, -1.395509944841997e-09);
    spin[2][1] =
        std::complex<double>(1.153549231359039e-10, 1.395509944841997e-09);
    spin[2][2] = std::complex<double>(1.681670319383889e-09, 0);

    // radiation variables
    double phi = 3.274309217448327;
    double xi = 0.6670187310700969;
    int bornpdgs[] = {0, -1, -13, 14, -2};
    int realpdgs[] = {2, -1, -13, 14, -2, 2};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 3.56627925693797e-15
    double fks_g = 3.173375048516869e-23;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, -1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, -1, -13, 14, -2
// realpdgs = 2, -1, -13, 14, -2, 2
TEST(QCDSoftCollinearISR2, Wj_born_1_real_25) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.0847015779963165;
    ps.X2 = 0.7954019687676492;
    ps.Jacobian = 376.0741751215844;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1687.145112343886, 0, 0, 1687.145112343886);
    ps.Momenta[1].Set(1687.145112343886, 0, 0, -1687.145112343886);
    ps.Momenta[2].Set(1492.617411530664, 25.916868466718, 275.3550763397022,
                      -1466.770137092835);
    ps.Momenta[3].Set(1660.50563855768, -120.9625378131744, -186.4060200209988,
                      1645.569760243416);
    ps.Momenta[4].Set(221.1671745994281, 95.04566934645636, -88.94905631870346,
                      -178.7996231505819);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.844281179204461e-08;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(9.315213228064565e-09, 0);
    spin[1][2] =
        std::complex<double>(2.819411072717112e-10, -9.216617389025134e-09);
    spin[2][1] =
        std::complex<double>(2.819411072717112e-10, 9.216617389025134e-09);
    spin[2][2] = std::complex<double>(9.127598563980048e-09, 0);

    // radiation variables
    double phi = 5.868292873842881;
    int bornpdgs[] = {0, -1, -13, 14, -2};
    int realpdgs[] = {2, -1, -13, 14, -2, 2};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 0.000360755657877787
    double fks_g = 3.020273138256386e-23;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, -1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, -1, -13, 14, -2
// realpdgs = 1, -1, -13, 14, -2, 1
TEST(QCDCollinearISR1, Wj_born_1_real_7) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.09842615279202604;
    ps.X2 = 0.5708980042871517;
    ps.Jacobian = 1178.227740234687;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1540.805691802868, 0, 0, 1540.805691802868);
    ps.Momenta[1].Set(1540.805691802868, 0, 0, -1540.805691802868);
    ps.Momenta[2].Set(1013.395669379556, -623.3687420369502, -688.9960267315038,
                      404.5573745672036);
    ps.Momenta[3].Set(1309.496703766136, 60.98111385338137, 1137.150261336611,
                      -646.4922304789552);
    ps.Momenta[4].Set(758.7190104600444, 562.3876281835687, -448.1542346051068,
                      241.9348559117516);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.920388400474834e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(6.726267942530981e-10, 0);
    spin[1][2] =
        std::complex<double>(-8.540907829368556e-10, -3.313711268160394e-10);
    spin[2][1] =
        std::complex<double>(-8.540907829368556e-10, 3.313711268160394e-10);
    spin[2][2] = std::complex<double>(1.247761606221736e-09, 0);

    // radiation variables
    double phi = 4.797306396660884;
    double xi = 0.73139755122204;
    int bornpdgs[] = {0, -1, -13, 14, -2};
    int realpdgs[] = {1, -1, -13, 14, -2, 1};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 1.779653274655055e-06
    double fks_g = 1.90402390828038e-14;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, 1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 1.904481874489658e-14
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, -1, -13, 14, -2
// realpdgs = 1, -1, -13, 14, -2, 1
TEST(QCDSoftCollinearISR1, Wj_born_1_real_7) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.7687275807334721;
    ps.X2 = 0.5423603406341329;
    ps.Jacobian = 9213.240089012847;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(4197.044274829208, 0, 0, 4197.044274829208);
    ps.Momenta[1].Set(4197.044274829208, 0, 0, -4197.044274829208);
    ps.Momenta[2].Set(4098.722965675915, -4021.904250156021, -696.4031822349274,
                      372.6107348281516);
    ps.Momenta[3].Set(2117.312629638866, 2101.21901661352, -260.5595673587057,
                      0.357362288779342);
    ps.Momenta[4].Set(2178.052954343636, 1920.685233542501, 956.962749593633,
                      -372.9680971169309);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 6.754263692291483e-10;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(7.578878429365016e-11, -4.308232357047019e-27);
    spin[1][2] =
        std::complex<double>(-7.728352675048865e-11, 1.986782828394898e-10);
    spin[2][1] =
        std::complex<double>(-7.728352675048865e-11, -1.986782828394898e-10);
    spin[2][2] =
        std::complex<double>(5.996375849354981e-10, 4.308232357047019e-27);

    // radiation variables
    double phi = 0.05891201425936118;
    int bornpdgs[] = {0, -1, -13, 14, -2};
    int realpdgs[] = {1, -1, -13, 14, -2, 1};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 0.0003278401049788915
    double fks_g = 3.507033622813431e-23;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, 1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, -1, -13, 14, -2
// realpdgs = 1, -1, -13, 14, -2, 1
TEST(QCDSoftLimit, Wj_born_1_real_7_0) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.763973485661932;
    ps.X2 = 0.1011929509771008;
    ps.Jacobian = 114.2715984066875;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1807.289103915334, 0, 0, 1807.289103915334);
    ps.Momenta[1].Set(1807.289103915334, 0, 0, -1807.289103915334);
    ps.Momenta[2].Set(1803.844398550767, -45.79477403909833, -1094.321817771175,
                      1433.25406401087);
    ps.Momenta[3].Set(1747.998756401476, 19.09976228192524, 1070.511269733383,
                      -1381.716495101546);
    ps.Momenta[4].Set(62.73505287842593, 26.69501175717309, 23.81054803779275,
                      -51.53756890932333);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    // double born = 6.407070699665886e-10;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, 9.610606049498829e-10);
    ColorCorr.Set(0, 4, 9.610606049498829e-10);
    ColorCorr.Set(1, 0, 9.610606049498829e-10);
    ColorCorr.Set(1, 4, -1.067845116610981e-10);
    ColorCorr.Set(4, 0, 9.610606049498829e-10);
    ColorCorr.Set(4, 1, -1.067845116610981e-10);

    // radiation variables
    double y = 0.1709205110464822;
    double phi = 1.512946721959911;
    int bornpdgs[] = {0, -1, -13, 14, -2};
    int realpdgs[] = {1, -1, -13, 14, -2, 1};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 1.266537072184461e-07
    double fks_g = 1.677192837724106e-24;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 0, alphas, ColorCorr, y, phi);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, -1, -13, 14, -2
// realpdgs = 1, -1, -13, 14, -2, 1
TEST(QCDCollinearISR2, Wj_born_1_real_7) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.8844603043821344;
    ps.X2 = 0.3960932036041953;
    ps.Jacobian = 14074.83960064421;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3847.257234269791, 0, 0, 3847.257234269791);
    ps.Momenta[1].Set(3847.257234269791, 0, 0, -3847.257234269791);
    ps.Momenta[2].Set(772.3631480670093, 591.8626106748341, 449.819885744174,
                      -209.5365194095027);
    ps.Momenta[3].Set(3292.275016898646, 1210.436533831412, -1005.2676355022,
                      -2891.946604879314);
    ps.Momenta[4].Set(3629.876303573927, -1802.299144506246, 555.4477497580261,
                      3101.483124288816);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 2.311557134869073e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(1.796503234581331e-09, 0);
    spin[1][2] =
        std::complex<double>(-3.656906582305252e-10, 8.89700140684256e-10);
    spin[2][1] =
        std::complex<double>(-3.656906582305252e-10, -8.89700140684256e-10);
    spin[2][2] = std::complex<double>(5.150539002877427e-10, 0);

    // radiation variables
    double phi = 5.358860804963292;
    double xi = 0.2983735361576609;
    int bornpdgs[] = {0, -1, -13, 14, -2};
    int realpdgs[] = {1, -1, -13, 14, -2, 1};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 3.100696248966356e-15
    double fks_g = 5.520899251897253e-24;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, -1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, -1, -13, 14, -2
// realpdgs = 1, -1, -13, 14, -2, 1
TEST(QCDSoftCollinearISR2, Wj_born_1_real_7) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.2183642600838276;
    ps.X2 = 0.9205110028751662;
    ps.Jacobian = 8228.67659337329;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2914.195128293315, 0, 0, 2914.195128293315);
    ps.Momenta[1].Set(2914.195128293315, 0, 0, -2914.195128293315);
    ps.Momenta[2].Set(2145.972630374488, -2038.314300683268, -660.7554472156289,
                      -117.7946557418987);
    ps.Momenta[3].Set(880.7864443906442, -559.9153821202408, -252.9448925996016,
                      631.1088707942819);
    ps.Momenta[4].Set(2801.631181821499, 2598.229682803508, 913.7003398152301,
                      -513.3142150523829);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.686890618748467e-08;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(3.632568072464284e-09, 2.757268708510092e-25);
    spin[1][2] =
        std::complex<double>(-2.701720232642749e-09, 6.386126135419964e-09);
    spin[2][1] =
        std::complex<double>(-2.701720232642749e-09, -6.386126135419964e-09);
    spin[2][2] =
        std::complex<double>(1.323633811502039e-08, -2.757268708510092e-25);

    // radiation variables
    double phi = 0.2674598209866651;
    int bornpdgs[] = {0, -1, -13, 14, -2};
    int realpdgs[] = {1, -1, -13, 14, -2, 1};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 2.470263536052514e-08
    double fks_g = 3.121673795408823e-28;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, -1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, -1, -13, 14, -2
// realpdgs = -2, -1, -13, 14, -2, -2
TEST(QCDCollinearISR1, Wj_born_1_real_34) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.5009024151073369;
    ps.X2 = 0.7247493488632131;
    ps.Jacobian = 7818.446446013481;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3916.371093360089, 0, 0, 3916.371093360089);
    ps.Momenta[1].Set(3916.371093360089, 0, 0, -3916.371093360089);
    ps.Momenta[2].Set(2180.56189431423, 1490.088492297333, -1520.973754659539,
                      -470.239617317793);
    ps.Momenta[3].Set(3671.400367837865, -3329.105201977276, 1108.999930857231,
                      1079.980725979175);
    ps.Momenta[4].Set(1980.779924568085, 1839.016709679943, 411.9738238023085,
                      -609.7411086613818);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 8.345014605241818e-10;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(5.023011217290512e-10, 0);
    spin[1][2] =
        std::complex<double>(-4.008020596233463e-12, -4.08470976430964e-10);
    spin[2][1] =
        std::complex<double>(-4.008020596233463e-12, 4.08470976430964e-10);
    spin[2][2] = std::complex<double>(3.322003387951306e-10, 0);

    // radiation variables
    double phi = 5.802762078289265;
    double xi = 0.2544797115980877;
    int bornpdgs[] = {0, -1, -13, 14, -2};
    int realpdgs[] = {-2, -1, -13, 14, -2, -2};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 3.198271898522339e-08
    double fks_g = 4.142396874788714e-17;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, 1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 4.142549549202445e-17
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, -1, -13, 14, -2
// realpdgs = -2, -1, -13, 14, -2, -2
TEST(QCDSoftCollinearISR1, Wj_born_1_real_34) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.01534709501823883;
    ps.X2 = 0.2899615153127755;
    ps.Jacobian = 144.1365325569509;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(433.6073427325314, 0, 0, 433.6073427325314);
    ps.Momenta[1].Set(433.6073427325314, 0, 0, -433.6073427325314);
    ps.Momenta[2].Set(152.1171288323945, -84.52853790464192, 109.3134431164243,
                      63.60124462577021);
    ps.Momenta[3].Set(385.2775176980138, -154.2606337869226, -286.4370692140778,
                      -206.3885362302129);
    ps.Momenta[4].Set(329.8200389346545, 238.7891716915645, 177.1236260976535,
                      142.7872916044427);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 3.427027966538917e-08;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(2.325916596240543e-08, 0);
    spin[1][2] =
        std::complex<double>(1.295857605704046e-08, -9.390667048858312e-09);
    spin[2][1] =
        std::complex<double>(1.295857605704046e-08, 9.390667048858312e-09);
    spin[2][2] = std::complex<double>(1.101111370298374e-08, 0);

    // radiation variables
    double phi = 3.232799374501395;
    int bornpdgs[] = {0, -1, -13, 14, -2};
    int realpdgs[] = {-2, -1, -13, 14, -2, -2};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 0.1829976029786621
    double fks_g = 3.548474863568298e-19;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, 1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, -1, -13, 14, -2
// realpdgs = -2, -1, -13, 14, -2, -2
TEST(QCDSoftLimit, Wj_born_1_real_34_0) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.768635332322436;
    ps.X2 = 0.7215874788060006;
    ps.Jacobian = 4638.997288751528;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(4840.809842774948, 0, 0, 4840.809842774948);
    ps.Momenta[1].Set(4840.809842774948, 0, 0, -4840.809842774948);
    ps.Momenta[2].Set(4523.508520346027, 3541.153268913659, 1233.980941042538,
                      2529.75372256846);
    ps.Momenta[3].Set(4207.274999349526, -2708.183079603762, -1201.62252624749,
                      -2987.140879154185);
    ps.Momenta[4].Set(950.8361658543434, -832.970189309897, -32.35841479504828,
                      457.3871565857249);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    // double born = 1.934523869725366e-10;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, 2.901785804588049e-10);
    ColorCorr.Set(0, 4, 2.901785804588049e-10);
    ColorCorr.Set(1, 0, 2.901785804588049e-10);
    ColorCorr.Set(1, 4, -3.224206449542277e-11);
    ColorCorr.Set(4, 0, 2.901785804588049e-10);
    ColorCorr.Set(4, 1, -3.224206449542277e-11);

    // radiation variables
    double y = 0.7234951429580079;
    double phi = 3.463127559995238;
    int bornpdgs[] = {0, -1, -13, 14, -2};
    int realpdgs[] = {-2, -1, -13, 14, -2, -2};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 1.120144256254337e-08
    double fks_g = 3.706042849071611e-26;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 0, alphas, ColorCorr, y, phi);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, -1, -13, 14, -2
// realpdgs = -2, -1, -13, 14, -2, -2
TEST(QCDCollinearISR2, Wj_born_1_real_34) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.8286544853750932;
    ps.X2 = 0.1240575109469759;
    ps.Jacobian = 1268.707763360116;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2084.066780271517, 0, 0, 2084.066780271517);
    ps.Momenta[1].Set(2084.066780271517, 0, 0, -2084.066780271517);
    ps.Momenta[2].Set(1670.322046703015, -778.9642671110435, -1410.102078726807,
                      -441.3644048138033);
    ps.Momenta[3].Set(1893.793956482275, 732.057087955726, 1400.603925726202,
                      1043.482924068679);
    ps.Momenta[4].Set(604.017557357744, 46.90717915531733, 9.49815300060439,
                      -602.1185192548752);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 3.184606241447028e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(1.786198877607122e-09, -3.446585885637615e-26);
    spin[1][2] =
        std::complex<double>(-3.564550185059997e-10, -1.539731627105926e-09);
    spin[2][1] =
        std::complex<double>(-3.564550185059997e-10, 1.539731627105926e-09);
    spin[2][2] =
        std::complex<double>(1.398407363839906e-09, 3.446585885637615e-26);

    // radiation variables
    double phi = 1.717167876535634;
    double xi = 0.2691194450242542;
    int bornpdgs[] = {0, -1, -13, 14, -2};
    int realpdgs[] = {-2, -1, -13, 14, -2, -2};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 7.381042710979451e-16
    double fks_g = 1.069148105884351e-24;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, -1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, -1, -13, 14, -2
// realpdgs = -2, -1, -13, 14, -2, -2
TEST(QCDSoftCollinearISR2, Wj_born_1_real_34) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.1242907805022959;
    ps.X2 = 0.1343285368944578;
    ps.Jacobian = 636.3341151646703;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(839.8794525621029, 0, 0, 839.8794525621029);
    ps.Momenta[1].Set(839.8794525621029, 0, 0, -839.8794525621029);
    ps.Momenta[2].Set(586.2965134490188, -43.07179754266373, 581.9349144008806,
                      -56.92255563622311);
    ps.Momenta[3].Set(341.722067009464, -274.9160374383094, 87.43781698855034,
                      183.1659673648885);
    ps.Momenta[4].Set(751.7403246657232, 317.9878349809731, -669.3727313894308,
                      -126.2434117286654);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 6.827154353419213e-08;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(5.843177199634649e-08, 0);
    spin[1][2] =
        std::complex<double>(-1.422402952577316e-08, 1.930368542650521e-08);
    spin[2][1] =
        std::complex<double>(-1.422402952577316e-08, -1.930368542650521e-08);
    spin[2][2] = std::complex<double>(9.83977153784564e-09, 0);

    // radiation variables
    double phi = 2.217433148733355;
    int bornpdgs[] = {0, -1, -13, 14, -2};
    int realpdgs[] = {-2, -1, -13, 14, -2, -2};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 5.526175234883479e-08
    double fks_g = 8.282489224783473e-26;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, -1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, -1, -13, 14, -2
// realpdgs = -1, -1, -13, 14, -2, -1
TEST(QCDCollinearISR1, Wj_born_1_real_2) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.7619300578602797;
    ps.X2 = 0.4212026093633554;
    ps.Jacobian = 4281.598822246106;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3682.276840502614, 0, 0, 3682.276840502614);
    ps.Momenta[1].Set(3682.276840502614, 0, 0, -3682.276840502614);
    ps.Momenta[2].Set(2739.708053658417, -1256.249867067813, 1035.108137266909,
                      -2203.721315170483);
    ps.Momenta[3].Set(3471.155616806589, 1810.484387387334, -2045.074350450566,
                      2142.134099500738);
    ps.Momenta[4].Set(1153.690010540222, -554.2345203195214, 1009.966213183657,
                      61.58721566974481);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 2.380340427671212e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(1.333244684650514e-09, 0);
    spin[1][2] =
        std::complex<double>(1.059726114710597e-10, -1.176777225867551e-09);
    spin[2][1] =
        std::complex<double>(1.059726114710597e-10, 1.176777225867551e-09);
    spin[2][2] = std::complex<double>(1.047095743020698e-09, 0);

    // radiation variables
    double phi = 4.599000879144474;
    double xi = 0.1778902608803425;
    int bornpdgs[] = {0, -1, -13, 14, -2};
    int realpdgs[] = {-1, -1, -13, 14, -2, -1};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 2.366540849564519e-07
    double fks_g = 1.497781095733908e-16;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, 1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 1.497699990701111e-16
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, -1, -13, 14, -2
// realpdgs = -1, -1, -13, 14, -2, -1
TEST(QCDSoftCollinearISR1, Wj_born_1_real_2) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.1737601750924291;
    ps.X2 = 0.5566972737985196;
    ps.Jacobian = 1295.381203661664;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2021.613023362271, 0, 0, 2021.613023362271);
    ps.Momenta[1].Set(2021.613023362271, 0, 0, -2021.613023362271);
    ps.Momenta[2].Set(1575.993997500351, -1359.570293543777, -279.7406484009716,
                      -746.3718019208843);
    ps.Momenta[3].Set(1831.463323272595, 1758.137051156029, 432.9074795337366,
                      275.3236786199623);
    ps.Momenta[4].Set(635.7687259515969, -398.5667576122516, -153.1668311327651,
                      471.0481233009219);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.382887772091806e-08;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(5.332202009979181e-09, -2.757268708510092e-25);
    spin[1][2] =
        std::complex<double>(-1.827793745595405e-09, -6.478052278854494e-09);
    spin[2][1] =
        std::complex<double>(-1.827793745595405e-09, 6.478052278854494e-09);
    spin[2][2] =
        std::complex<double>(8.49667571093888e-09, 2.757268708510092e-25);

    // radiation variables
    double phi = 4.729447831185459;
    int bornpdgs[] = {0, -1, -13, 14, -2};
    int realpdgs[] = {-1, -1, -13, 14, -2, -1};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 0.008097002776105236
    double fks_g = 1.105519913962897e-20;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, 1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, -1, -13, 14, -2
// realpdgs = -1, -1, -13, 14, -2, -1
TEST(QCDSoftLimit, Wj_born_1_real_2_0) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.6869607166203728;
    ps.X2 = 0.3361660105968483;
    ps.Jacobian = 6286.304843806281;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3123.605711304279, 0, 0, 3123.605711304279);
    ps.Momenta[1].Set(3123.605711304279, 0, 0, -3123.605711304279);
    ps.Momenta[2].Set(2738.055415693766, 2376.13982808871, 1051.305021896009,
                      -863.5188056532104);
    ps.Momenta[3].Set(1512.336430580968, -1427.807557310962, 497.4174925126259,
                      33.20988844509515);
    ps.Momenta[4].Set(1996.819576333825, -948.3322707777482, -1548.722514408635,
                      830.3089172081151);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    // double born = 2.80312473380606e-09;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, 4.20468710070909e-09);
    ColorCorr.Set(0, 4, 4.20468710070909e-09);
    ColorCorr.Set(1, 0, 4.20468710070909e-09);
    ColorCorr.Set(1, 4, -4.671874556343433e-10);
    ColorCorr.Set(4, 0, 4.20468710070909e-09);
    ColorCorr.Set(4, 1, -4.671874556343433e-10);

    // radiation variables
    double y = 0.3350240336923704;
    double phi = 1.378817748427015;
    int bornpdgs[] = {0, -1, -13, 14, -2};
    int realpdgs[] = {-1, -1, -13, 14, -2, -1};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 2.831872924497448e-07
    double fks_g = 4.693555925074885e-24;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 0, alphas, ColorCorr, y, phi);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, -1, -13, 14, -2
// realpdgs = -1, -1, -13, 14, -2, -1
TEST(QCDCollinearISR2, Wj_born_1_real_2) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.08881691613774834;
    ps.X2 = 0.617571035736546;
    ps.Jacobian = 2290.858745323132;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1522.315471282929, 0, 0, 1522.315471282929);
    ps.Momenta[1].Set(1522.315471282929, 0, 0, -1522.315471282929);
    ps.Momenta[2].Set(463.5533415929864, -279.7138708682052, -284.7889975329686,
                      235.6630599608599);
    ps.Momenta[3].Set(1087.962641653964, -55.99714785677401, -899.994352180111,
                      608.7176645296909);
    ps.Momenta[4].Set(1493.114959318906, 335.7110187249794, 1184.78334971308,
                      -844.3807244905514);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 6.955317749823015e-08;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(9.015241605931842e-09, 0);
    spin[1][2] =
        std::complex<double>(1.150815003618604e-08, -2.033043534059498e-08);
    spin[2][1] =
        std::complex<double>(1.150815003618604e-08, 2.033043534059498e-08);
    spin[2][2] = std::complex<double>(6.053793589229831e-08, 0);

    // radiation variables
    double phi = 4.977235951205702;
    double xi = 0.1151145236251334;
    int bornpdgs[] = {0, -1, -13, 14, -2};
    int realpdgs[] = {-1, -1, -13, 14, -2, -1};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 6.544621399787415e-05
    double fks_g = 1.734501839403498e-14;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, -1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, -1, -13, 14, -2
// realpdgs = -1, -1, -13, 14, -2, -1
TEST(QCDSoftCollinearISR2, Wj_born_1_real_2) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.3971001323140619;
    ps.X2 = 0.7270586265187902;
    ps.Jacobian = 3721.642040426628;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3492.593877679843, 0, 0, 3492.593877679843);
    ps.Momenta[1].Set(3492.593877679843, 0, 0, -3492.593877679843);
    ps.Momenta[2].Set(2466.319798360556, 926.8574637141095, 1924.141298654915,
                      1233.429711233387);
    ps.Momenta[3].Set(3461.597475728877, -1656.757868688546, -2636.30423192476,
                      -1512.517915677933);
    ps.Momenta[4].Set(1057.270481270255, 729.9004049744369, 712.1629332698446,
                      279.0882044445462);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 4.627072997507311e-10;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(2.449136072801405e-10, 0);
    spin[1][2] =
        std::complex<double>(5.681449453562776e-12, -2.308860324877634e-10);
    spin[2][1] =
        std::complex<double>(5.681449453562776e-12, 2.308860324877634e-10);
    spin[2][2] = std::complex<double>(2.177936924705906e-10, 0);

    // radiation variables
    double phi = 1.073161930260761;
    int bornpdgs[] = {0, -1, -13, 14, -2};
    int realpdgs[] = {-1, -1, -13, 14, -2, -1};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 0.00104055967647148
    double fks_g = 1.550371909797342e-22;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, -1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-6) << "limit = " << limit
                                         << ", fks_g = " << fks_g;
}

// bornpdgs = 0, -1, -13, 14, -2
// realpdgs = 3, -1, -13, 14, -2, 3
TEST(QCDCollinearISR1, Wj_born_1_real_17) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.1700830414284162;
    ps.X2 = 0.5286458568785548;
    ps.Jacobian = 2630.736180521104;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1949.064806825164, 0, 0, 1949.064806825164);
    ps.Momenta[1].Set(1949.064806825164, 0, 0, -1949.064806825164);
    ps.Momenta[2].Set(1806.948450359055, 842.1714214482886, 277.1900261133375,
                      -1574.476322011241);
    ps.Momenta[3].Set(751.9652094578083, -92.25090945065659, -554.0170887115429,
                      500.0065113136415);
    ps.Momenta[4].Set(1339.215953833465, -749.920511997632, 276.8270625982055,
                      1074.4698106976);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 3.355485864693342e-08;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(4.136448637313496e-09, 0);
    spin[1][2] =
        std::complex<double>(5.48265692597866e-09, 9.572262795611091e-09);
    spin[2][1] =
        std::complex<double>(5.48265692597866e-09, -9.572262795611091e-09);
    spin[2][2] = std::complex<double>(2.941841000961993e-08, 0);

    // radiation variables
    double phi = 5.290991500304154;
    double xi = 0.4481032372571858;
    int bornpdgs[] = {0, -1, -13, 14, -2};
    int realpdgs[] = {3, -1, -13, 14, -2, 3};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 8.502867646558198e-06
    double fks_g = 3.414692316137656e-14;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, 1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 3.414376885801069e-14
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, -1, -13, 14, -2
// realpdgs = 3, -1, -13, 14, -2, 3
TEST(QCDSoftCollinearISR1, Wj_born_1_real_17) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.9697532130532771;
    ps.X2 = 0.3404124738139771;
    ps.Jacobian = 8888.108387908564;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3734.62244582112, 0, 0, 3734.62244582112);
    ps.Momenta[1].Set(3734.62244582112, 0, 0, -3734.62244582112);
    ps.Momenta[2].Set(2766.939118300098, -1147.895320073962, 1179.143709685032,
                      -2224.389473640315);
    ps.Momenta[3].Set(2340.945593346594, -901.4952342635661, -1783.590966341345,
                      1219.071728158135);
    ps.Momenta[4].Set(2361.360179995547, 2049.390554337529, 604.4472566563131,
                      1005.31774548218);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 2.214505710398975e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(3.647583123156395e-11, 0);
    spin[1][2] =
        std::complex<double>(9.425252755866224e-11, -2.65634921158002e-10);
    spin[2][1] =
        std::complex<double>(9.425252755866224e-11, 2.65634921158002e-10);
    spin[2][2] = std::complex<double>(2.178029879167411e-09, 0);

    // radiation variables
    double phi = 4.321895517127354;
    int bornpdgs[] = {0, -1, -13, 14, -2};
    int realpdgs[] = {3, -1, -13, 14, -2, 3};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 0.01042325131065879
    double fks_g = 1.9071809922491e-23;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, 1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, -1, -13, 14, -2
// realpdgs = 3, -1, -13, 14, -2, 3
TEST(QCDSoftLimit, Wj_born_1_real_17_0) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.6887023855403562;
    ps.X2 = 0.808772521904217;
    ps.Jacobian = 11275.49292210567;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(4851.123646072686, 0, 0, 4851.123646072686);
    ps.Momenta[1].Set(4851.123646072686, 0, 0, -4851.123646072686);
    ps.Momenta[2].Set(2706.110189398457, 610.092400241746, -2598.65713082466,
                      -444.7479474349083);
    ps.Momenta[3].Set(4689.9592931907, 307.1931577179431, 4679.488943313354,
                      61.10454193560722);
    ps.Momenta[4].Set(2306.177809556214, -917.2855579596891, -2080.831812488693,
                      383.643405499301);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    // double born = 5.538977686269973e-10;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, 8.30846652940496e-10);
    ColorCorr.Set(0, 4, 8.30846652940496e-10);
    ColorCorr.Set(1, 0, 8.30846652940496e-10);
    ColorCorr.Set(1, 4, -9.231629477116622e-11);
    ColorCorr.Set(4, 0, 8.30846652940496e-10);
    ColorCorr.Set(4, 1, -9.231629477116622e-11);

    // radiation variables
    double y = -0.366143369701426;
    double phi = 0.4144822879749155;
    int bornpdgs[] = {0, -1, -13, 14, -2};
    int realpdgs[] = {3, -1, -13, 14, -2, 3};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 1.273381095530638e-08
    double fks_g = 7.890188714287815e-26;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 0, alphas, ColorCorr, y, phi);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, -1, -13, 14, -2
// realpdgs = 3, -1, -13, 14, -2, 3
TEST(QCDCollinearISR2, Wj_born_1_real_17) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.1484057281277096;
    ps.X2 = 0.4726213565781237;
    ps.Jacobian = 2271.1325540867;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1721.453753171597, 0, 0, 1721.453753171597);
    ps.Momenta[1].Set(1721.453753171597, 0, 0, -1721.453753171597);
    ps.Momenta[2].Set(1203.519468900981, 506.3921286968242, 371.31822856572,
                      -1026.717535231776);
    ps.Momenta[3].Set(930.3667204088309, -216.3879189966118, -894.3108547192063,
                      -137.7192727544049);
    ps.Momenta[4].Set(1309.021317033382, -290.0042097002125, 522.9926261534863,
                      1164.436807986181);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 5.005569104956846e-08;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(2.679184372266386e-08, 0);
    spin[1][2] =
        std::complex<double>(1.996817364554767e-08, 1.498510604672502e-08);
    spin[2][1] =
        std::complex<double>(1.996817364554767e-08, -1.498510604672502e-08);
    spin[2][2] = std::complex<double>(2.32638473269046e-08, 0);

    // radiation variables
    double phi = 3.006345778771419;
    double xi = 0.46152893644024;
    int bornpdgs[] = {0, -1, -13, 14, -2};
    int realpdgs[] = {3, -1, -13, 14, -2, 3};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 7.219205376942111e-15
    double fks_g = 3.075510845158954e-23;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, -1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, -1, -13, 14, -2
// realpdgs = 3, -1, -13, 14, -2, 3
TEST(QCDSoftCollinearISR2, Wj_born_1_real_17) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.2167577146343156;
    ps.X2 = 0.7727723998172991;
    ps.Jacobian = 6400.13205384654;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2660.274426847332, 0, 0, 2660.274426847332);
    ps.Momenta[1].Set(2660.274426847332, 0, 0, -2660.274426847332);
    ps.Momenta[2].Set(296.4318570945322, -116.9671868384885, 155.5448421556499,
                      -223.5985804569813);
    ps.Momenta[3].Set(2637.063832468389, -113.7878988643132, -2078.559644057908,
                      1618.872378132005);
    ps.Momenta[4].Set(2387.053164131743, 230.7550857028018, 1923.014801902258,
                      -1395.273797675023);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 6.995629836866172e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(3.460370575704956e-09, -4.308232357047019e-27);
    spin[1][2] =
        std::complex<double>(-6.647446342196777e-11, -3.496982738093438e-09);
    spin[2][1] =
        std::complex<double>(-6.647446342196777e-11, 3.496982738093438e-09);
    spin[2][2] =
        std::complex<double>(3.535259261161216e-09, 4.308232357047019e-27);

    // radiation variables
    double phi = 1.622661156071877;
    int bornpdgs[] = {0, -1, -13, 14, -2};
    int realpdgs[] = {3, -1, -13, 14, -2, 3};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 4.300443749755643e-09
    double fks_g = 4.44084482969026e-28;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, -1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, -1, -13, 14, -2
// realpdgs = 4, -1, -13, 14, -2, 4
TEST(QCDCollinearISR1, Wj_born_1_real_17_n1) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.8795619480533752;
    ps.X2 = 0.215200237904007;
    ps.Jacobian = 2699.479973370911;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2827.925385323811, 0, 0, 2827.925385323811);
    ps.Momenta[1].Set(2827.925385323811, 0, 0, -2827.925385323811);
    ps.Momenta[2].Set(2162.321069475402, -1613.677269528281, -1179.555376306846,
                      -824.8194902708178);
    ps.Momenta[3].Set(2546.394975407223, 1371.65697662884, 2063.676580304305,
                      586.4498965368716);
    ps.Momenta[4].Set(947.1347257649966, 242.0202928994406, -884.1212039974596,
                      238.3695937339463);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 2.673912299555396e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(1.20969331274332e-09, 0);
    spin[1][2] =
        std::complex<double>(4.478072811591321e-10, -1.253285504455881e-09);
    spin[2][1] =
        std::complex<double>(4.478072811591321e-10, 1.253285504455881e-09);
    spin[2][2] = std::complex<double>(1.464218986812076e-09, 0);

    // radiation variables
    double phi = 5.070508956808587;
    double xi = 0.03107407911653776;
    int bornpdgs[] = {0, -1, -13, 14, -2};
    int realpdgs[] = {4, -1, -13, 14, -2, 4};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 2.17748220247513e-06
    double fks_g = 4.205146628521165e-17;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, 1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 4.205198461422234e-17
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, -1, -13, 14, -2
// realpdgs = 4, -1, -13, 14, -2, 4
TEST(QCDSoftCollinearISR1, Wj_born_1_real_17_n1) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.4050845717699509;
    ps.X2 = 0.438403598409252;
    ps.Jacobian = 1842.57098262014;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2739.19697325507, 0, 0, 2739.19697325507);
    ps.Momenta[1].Set(2739.19697325507, 0, 0, -2739.19697325507);
    ps.Momenta[2].Set(2116.720355482532, -1379.663718192918, 1361.615734839009,
                      850.3150467062264);
    ps.Momenta[3].Set(2694.25152538686, 1924.471239919917, -1429.276418117284,
                      -1229.866110348742);
    ps.Momenta[4].Set(667.4220656407479, -544.8075217269989, 67.66068327827448,
                      379.5510636425162);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.353261653275337e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(6.967278309308959e-10, 0);
    spin[1][2] =
        std::complex<double>(-8.028048700883726e-11, -6.71550764559491e-10);
    spin[2][1] =
        std::complex<double>(-8.028048700883726e-11, 6.71550764559491e-10);
    spin[2][2] = std::complex<double>(6.565338223444416e-10, 0);

    // radiation variables
    double phi = 5.659920318754485;
    int bornpdgs[] = {0, -1, -13, 14, -2};
    int realpdgs[] = {4, -1, -13, 14, -2, 4};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 0.0005992518862185885
    double fks_g = 4.241797746338795e-22;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, 1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, -1, -13, 14, -2
// realpdgs = 4, -1, -13, 14, -2, 4
TEST(QCDSoftLimit, Wj_born_1_real_17_0_n1) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.6123538297690736;
    ps.X2 = 0.8341650906937232;
    ps.Jacobian = 5267.632044308142;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(4645.586824150115, 0, 0, 4645.586824150115);
    ps.Momenta[1].Set(4645.586824150115, 0, 0, -4645.586824150115);
    ps.Momenta[2].Set(4490.167951848332, -63.52233074346157, 4400.909293945158,
                      -888.577816384077);
    ps.Momenta[3].Set(3675.948903941561, -644.766261841811, -3560.338209348711,
                      648.7438994174148);
    ps.Momenta[4].Set(1125.056792510337, 708.2885925852726, -840.5710845964481,
                      239.8339169666622);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    // double born = 9.392983780866885e-10;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, 1.408947567130033e-09);
    ColorCorr.Set(0, 4, 1.408947567130033e-09);
    ColorCorr.Set(1, 0, 1.408947567130033e-09);
    ColorCorr.Set(1, 4, -1.565497296811147e-10);
    ColorCorr.Set(4, 0, 1.408947567130033e-09);
    ColorCorr.Set(4, 1, -1.565497296811147e-10);

    // radiation variables
    double y = -0.7895144773721441;
    double phi = 3.081639218169141;
    int bornpdgs[] = {0, -1, -13, 14, -2};
    int realpdgs[] = {4, -1, -13, 14, -2, 4};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 2.62003514149644e-08
    double fks_g = 3.32419042788146e-26;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 0, alphas, ColorCorr, y, phi);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, -1, -13, 14, -2
// realpdgs = 4, -1, -13, 14, -2, 4
TEST(QCDCollinearISR2, Wj_born_1_real_17_n1) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.8450100862902339;
    ps.X2 = 0.7041264128369136;
    ps.Jacobian = 22684.6419894367;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(5013.830188267336, 0, 0, 5013.830188267336);
    ps.Momenta[1].Set(5013.830188267336, 0, 0, -5013.830188267336);
    ps.Momenta[2].Set(3538.73779568208, -2559.168235413514, -2097.011006932431,
                      -1255.335798200454);
    ps.Momenta[3].Set(1999.795445664387, -1492.854996854404, 1242.03518386619,
                      -477.4037964923842);
    ps.Momenta[4].Set(4489.127135188204, 4052.023232267919, 854.9758230662419,
                      1732.739594692839);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 2.993139878832116e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(1.981318896702269e-09, 0);
    spin[1][2] =
        std::complex<double>(-5.195813934313038e-10, 1.317108654429732e-09);
    spin[2][1] =
        std::complex<double>(-5.195813934313038e-10, -1.317108654429732e-09);
    spin[2][2] = std::complex<double>(1.011820982129847e-09, 0);

    // radiation variables
    double phi = 1.973038413273342;
    double xi = 0.0835490036432065;
    int bornpdgs[] = {0, -1, -13, 14, -2};
    int realpdgs[] = {4, -1, -13, 14, -2, 4};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 1.433973822401102e-15
    double fks_g = 2.001952500338224e-25;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, -1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, -1, -13, 14, -2
// realpdgs = 4, -1, -13, 14, -2, 4
TEST(QCDSoftCollinearISR2, Wj_born_1_real_17_n1) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.2183907262440492;
    ps.X2 = 0.6075509908348735;
    ps.Jacobian = 2375.895060475542;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2367.673534192564, 0, 0, 2367.673534192564);
    ps.Momenta[1].Set(2367.673534192564, 0, 0, -2367.673534192564);
    ps.Momenta[2].Set(1506.988676586175, -124.9477473908146, 1246.836826607111,
                      -837.1384948723942);
    ps.Momenta[3].Set(2232.712212785211, 630.6933103330391, -1312.860618155376,
                      1692.225390046512);
    ps.Momenta[4].Set(995.6461790137442, -505.7455629422245, 66.02379154826453,
                      -855.0868951741178);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 3.060370402923478e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(1.687585905198639e-09, -8.616464714094038e-27);
    spin[1][2] =
        std::complex<double>(-1.382195160184583e-10, -1.515779381911251e-09);
    spin[2][1] =
        std::complex<double>(-1.382195160184583e-10, 1.515779381911251e-09);
    spin[2][2] =
        std::complex<double>(1.37278449772484e-09, 8.616464714094038e-27);

    // radiation variables
    double phi = 0.5841226502461656;
    int bornpdgs[] = {0, -1, -13, 14, -2};
    int realpdgs[] = {4, -1, -13, 14, -2, 4};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 1.375140052197762e-09
    double fks_g = 4.235878874994145e-28;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, -1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, -1, -13, 14, -2
// realpdgs = 5, -1, -13, 14, -2, 5
TEST(QCDCollinearISR1, Wj_born_1_real_17_n2) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.05162608603655272;
    ps.X2 = 0.02245220413280435;
    ps.Jacobian = 3.280705598045086;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(221.2979791840054, 0, 0, 221.2979791840054);
    ps.Momenta[1].Set(221.2979791840054, 0, 0, -221.2979791840054);
    ps.Momenta[2].Set(215.7672936460447, -180.0676376620092, 114.9100357441471,
                      30.44428615423016);
    ps.Momenta[3].Set(212.1194527022739, 181.7484056353065, -107.9468429988108,
                      -17.5971119882317);
    ps.Momenta[4].Set(14.70921201969244, -1.680767973297372, -6.963192745336302,
                      -12.84717416599846);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 4.137674826121312e-07;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(2.038622499547445e-07, 0);
    spin[1][2] =
        std::complex<double>(3.85843983152527e-08, -2.03231388356865e-07);
    spin[2][1] =
        std::complex<double>(3.85843983152527e-08, 2.03231388356865e-07);
    spin[2][2] = std::complex<double>(2.099052326573867e-07, 0);

    // radiation variables
    double phi = 4.39674838042675;
    double xi = 0.9373343384288418;
    int bornpdgs[] = {0, -1, -13, 14, -2};
    int realpdgs[] = {5, -1, -13, 14, -2, 5};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 0.05989627014105456
    double fks_g = 1.052492061761634e-09;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, 1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 1.052932138536187e-09
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, -1, -13, 14, -2
// realpdgs = 5, -1, -13, 14, -2, 5
TEST(QCDSoftCollinearISR1, Wj_born_1_real_17_n2) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.687136647855862;
    ps.X2 = 0.5350253333101378;
    ps.Jacobian = 6637.212978211919;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3941.142025930719, 0, 0, 3941.142025930719);
    ps.Momenta[1].Set(3941.142025930719, 0, 0, -3941.142025930719);
    ps.Momenta[2].Set(3357.179157155109, 3013.345596186435, -1217.932437164968,
                      -840.8571755547482);
    ps.Momenta[3].Set(2854.155626625255, -1672.28328532339, 2204.652763331288,
                      699.413431212417);
    ps.Momenta[4].Set(1670.949268081074, -1341.062310863045, -986.7203261663201,
                      141.4437443423312);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 9.050571489727599e-10;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(4.990177868295062e-11, 0);
    spin[1][2] =
        std::complex<double>(-1.477777522122539e-11, -2.060470610983074e-10);
    spin[2][1] =
        std::complex<double>(-1.477777522122539e-11, 2.060470610983074e-10);
    spin[2][2] = std::complex<double>(8.551553702898093e-10, 0);

    // radiation variables
    double phi = 3.582651406787277;
    int bornpdgs[] = {0, -1, -13, 14, -2};
    int realpdgs[] = {5, -1, -13, 14, -2, 5};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 0.0003683527980111953
    double fks_g = 7.211133012669063e-23;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, 1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, -1, -13, 14, -2
// realpdgs = 5, -1, -13, 14, -2, 5
TEST(QCDSoftLimit, Wj_born_1_real_17_0_n2) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.2900628152662144;
    ps.X2 = 0.8321158066012657;
    ps.Jacobian = 7524.740201549803;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3193.384929814182, 0, 0, 3193.384929814182);
    ps.Momenta[1].Set(3193.384929814182, 0, 0, -3193.384929814182);
    ps.Momenta[2].Set(1921.404078099025, 1631.969628506218, 479.194305989647,
                      -893.7793799791648);
    ps.Momenta[3].Set(2127.390931498796, 422.1877916366664, -565.4659038343769,
                      2006.937456827452);
    ps.Momenta[4].Set(2337.974850030542, -2054.157420142885, 86.2715978447299,
                      -1113.158076848288);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    // double born = 2.023348906395561e-09;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, 3.035023359593342e-09);
    ColorCorr.Set(0, 4, 3.035023359593342e-09);
    ColorCorr.Set(1, 0, 3.035023359593342e-09);
    ColorCorr.Set(1, 4, -3.372248177325935e-10);
    ColorCorr.Set(4, 0, 3.035023359593342e-09);
    ColorCorr.Set(4, 1, -3.372248177325935e-10);

    // radiation variables
    double y = 0.2350634844343418;
    double phi = 4.420701294662219;
    int bornpdgs[] = {0, -1, -13, 14, -2};
    int realpdgs[] = {5, -1, -13, 14, -2, 5};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 1.349778875960534e-07
    double fks_g = 1.840767377857709e-24;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 0, alphas, ColorCorr, y, phi);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, -1, -13, 14, -2
// realpdgs = 5, -1, -13, 14, -2, 5
TEST(QCDCollinearISR2, Wj_born_1_real_17_n2) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.258865561149431;
    ps.X2 = 0.8031918762725851;
    ps.Jacobian = 7718.701348494738;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2963.876809339941, 0, 0, 2963.876809339941);
    ps.Momenta[1].Set(2963.876809339941, 0, 0, -2963.876809339941);
    ps.Momenta[2].Set(2124.246200590715, -1261.800899168695, 1653.235884704746,
                      -432.5407739227344);
    ps.Momenta[3].Set(1219.559919382379, 351.0670715002116, 714.4380636187371,
                      923.9353665300021);
    ps.Momenta[4].Set(2583.947498706788, 910.7338276684837, -2367.673948323484,
                      -491.3945926072677);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 4.605913041493761e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(2.196064415235913e-09, 0);
    spin[1][2] =
        std::complex<double>(2.051624463010718e-09, 1.04068240929017e-09);
    spin[2][1] =
        std::complex<double>(2.051624463010718e-09, -1.04068240929017e-09);
    spin[2][2] = std::complex<double>(2.409848626257848e-09, 0);

    // radiation variables
    double phi = 3.411255429663737;
    double xi = 0.13659361407889;
    int bornpdgs[] = {0, -1, -13, 14, -2};
    int realpdgs[] = {5, -1, -13, 14, -2, 5};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 3.296554548098746e-15
    double fks_g = 1.230130124111919e-24;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, -1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, -1, -13, 14, -2
// realpdgs = 5, -1, -13, 14, -2, 5
TEST(QCDSoftCollinearISR2, Wj_born_1_real_17_n2) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.5080577166616287;
    ps.X2 = 0.9133773894697015;
    ps.Jacobian = 11089.19352350855;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(4427.87152110336, 0, 0, 4427.87152110336);
    ps.Momenta[1].Set(4427.87152110336, 0, 0, -4427.87152110336);
    ps.Momenta[2].Set(2105.400095896536, 1570.483179149309, -1383.417296444333,
                      229.0168808382319);
    ps.Momenta[3].Set(4265.467938261232, -2422.574081777605, 2858.749404452063,
                      -2037.867363981007);
    ps.Momenta[4].Set(2484.875008048953, 852.0909026282957, -1475.332108007729,
                      1808.850483142775);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 5.535120956696211e-10;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(3.501452082096654e-10, -6.462348535570529e-27);
    spin[1][2] =
        std::complex<double>(2.931540257219637e-11, -2.652330076077408e-10);
    spin[2][1] =
        std::complex<double>(2.931540257219637e-11, 2.652330076077408e-10);
    spin[2][2] =
        std::complex<double>(2.033668874599558e-10, 6.462348535570529e-27);

    // radiation variables
    double phi = 0.5797228001696154;
    int bornpdgs[] = {0, -1, -13, 14, -2};
    int realpdgs[] = {5, -1, -13, 14, -2, 5};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 3.221854301541872e-10
    double fks_g = 4.835023915833429e-30;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, -1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, -1, -13, 14, -2
// realpdgs = -5, -1, -13, 14, -2, -5
TEST(QCDCollinearISR1, Wj_born_1_real_12) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.6748118914964145;
    ps.X2 = 0.1311214121132096;
    ps.Jacobian = 392.0732448499118;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1933.488213884526, 0, 0, 1933.488213884526);
    ps.Momenta[1].Set(1933.488213884526, 0, 0, -1933.488213884526);
    ps.Momenta[2].Set(1781.585205729765, 323.7884919443782, 1719.660459995432,
                      334.6262991575139);
    ps.Momenta[3].Set(1884.192473845867, -492.701599454986, -1799.895928603237,
                      -260.386364061941);
    ps.Momenta[4].Set(201.1987481934198, 168.9131075106078, 80.23546860780529,
                      -74.23993509557286);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 4.175749747471624e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(2.276978726967185e-09, 0);
    spin[1][2] =
        std::complex<double>(1.674781208771169e-10, -2.072537647450065e-09);
    spin[2][1] =
        std::complex<double>(1.674781208771169e-10, 2.072537647450065e-09);
    spin[2][2] = std::complex<double>(1.89877102050444e-09, 0);

    // radiation variables
    double phi = 5.057954048050618;
    double xi = 0.2264450094201777;
    int bornpdgs[] = {0, -1, -13, 14, -2};
    int realpdgs[] = {-5, -1, -13, 14, -2, -5};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 1.256318043123879e-06
    double fks_g = 1.288413005801933e-15;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, 1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 1.288449317461593e-15
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, -1, -13, 14, -2
// realpdgs = -5, -1, -13, 14, -2, -5
TEST(QCDSoftCollinearISR1, Wj_born_1_real_12) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.8787333952442857;
    ps.X2 = 0.7438494503614521;
    ps.Jacobian = 11765.92934171706;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(5255.141878871463, 0, 0, 5255.141878871463);
    ps.Momenta[1].Set(5255.141878871463, 0, 0, -5255.141878871463);
    ps.Momenta[2].Set(3915.477271856347, 2706.581142585379, 2804.098112834322,
                      377.3785349272539);
    ps.Momenta[3].Set(4373.331707977937, -3557.543077208958, -1559.869938081246,
                      -2009.159838857586);
    ps.Momenta[4].Set(2221.474777908642, 850.9619346235792, -1244.228174753076,
                      1631.781303930332);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 4.952273557739867e-10;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(4.847264703341293e-10, 0);
    spin[1][2] =
        std::complex<double>(-1.057189823756099e-11, -7.055701313335032e-11);
    spin[2][1] =
        std::complex<double>(-1.057189823756099e-11, 7.055701313335032e-11);
    spin[2][2] = std::complex<double>(1.050088543985737e-11, 0);

    // radiation variables
    double phi = 5.845006233936573;
    int bornpdgs[] = {0, -1, -13, 14, -2};
    int realpdgs[] = {-5, -1, -13, 14, -2, -5};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 0.0002920271268798801
    double fks_g = 8.588865834043109e-24;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, 1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, -1, -13, 14, -2
// realpdgs = -5, -1, -13, 14, -2, -5
TEST(QCDSoftLimit, Wj_born_1_real_12_0) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.4448861122385459;
    ps.X2 = 0.1026127905725716;
    ps.Jacobian = 111.5123933726894;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1388.796234457987, 0, 0, 1388.796234457987);
    ps.Momenta[1].Set(1388.796234457987, 0, 0, -1388.796234457987);
    ps.Momenta[2].Set(1349.665484196803, -419.823463734903, -920.1612002505623,
                      893.6713848407953);
    ps.Momenta[3].Set(1348.258931793139, 382.5138438570275, 877.3814585855381,
                      -949.5720523204175);
    ps.Momenta[4].Set(79.66805292603225, 37.3096198778756, 42.77974166502428,
                      55.9006674796222);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    // double born = 7.769106055759288e-09;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, 1.165365908363893e-08);
    ColorCorr.Set(0, 4, 1.165365908363893e-08);
    ColorCorr.Set(1, 0, 1.165365908363893e-08);
    ColorCorr.Set(1, 4, -1.294851009293215e-09);
    ColorCorr.Set(4, 0, 1.165365908363893e-08);
    ColorCorr.Set(4, 1, -1.294851009293215e-09);

    // radiation variables
    double y = 0.9633381779940393;
    double phi = 2.191402489098995;
    int bornpdgs[] = {0, -1, -13, 14, -2};
    int realpdgs[] = {-5, -1, -13, 14, -2, -5};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 3.86151460125633e-05
    double fks_g = 8.795751187025891e-23;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 0, alphas, ColorCorr, y, phi);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, -1, -13, 14, -2
// realpdgs = -5, -1, -13, 14, -2, -5
TEST(QCDCollinearISR2, Wj_born_1_real_12) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.753355357836714;
    ps.X2 = 0.2098973586353452;
    ps.Jacobian = 3889.750325294726;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2584.739525237856, 0, 0, 2584.739525237856);
    ps.Momenta[1].Set(2584.739525237856, 0, 0, -2584.739525237856);
    ps.Momenta[2].Set(2555.549098337114, -2062.371651507209, -380.7873443470451,
                      -1460.292903309534);
    ps.Momenta[3].Set(1120.776171641962, 1018.334848760255, 239.4942957875321,
                      402.2136807793527);
    ps.Momenta[4].Set(1493.153780496636, 1044.036802746954, 141.293048559513,
                      1058.079222530181);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 9.663570141658432e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(1.066154444449661e-09, 0);
    spin[1][2] =
        std::complex<double>(-1.213038420289421e-09, 2.773934164190304e-09);
    spin[2][1] =
        std::complex<double>(-1.213038420289421e-09, -2.773934164190304e-09);
    spin[2][2] = std::complex<double>(8.597415697208771e-09, 0);

    // radiation variables
    double phi = 3.944138175873924;
    double xi = 0.4987363700537236;
    int bornpdgs[] = {0, -1, -13, 14, -2};
    int realpdgs[] = {-5, -1, -13, 14, -2, -5};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 7.363292562070251e-16
    double fks_g = 3.663060839969304e-24;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, -1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, -1, -13, 14, -2
// realpdgs = -5, -1, -13, 14, -2, -5
TEST(QCDSoftCollinearISR2, Wj_born_1_real_12) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.7902718666192463;
    ps.X2 = 0.1125353395050652;
    ps.Jacobian = 427.1937966149412;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1938.411957319022, 0, 0, 1938.411957319022);
    ps.Momenta[1].Set(1938.411957319022, 0, 0, -1938.411957319022);
    ps.Momenta[2].Set(1752.256146971812, 218.7571493974348, -1733.48831383036,
                      -132.5336938351278);
    ps.Momenta[3].Set(1905.903180850201, -332.9883020920556, 1876.544176197006,
                      12.94913996624604);
    ps.Momenta[4].Set(218.6645868160317, 114.2311526946209, -143.0558623666453,
                      119.5845538688818);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.620454969484637e-08;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(7.756244414150807e-09, 0);
    spin[1][2] =
        std::complex<double>(4.254966699351759e-10, -8.083691806737177e-09);
    spin[2][1] =
        std::complex<double>(4.254966699351759e-10, 8.083691806737177e-09);
    spin[2][2] = std::complex<double>(8.448305280695559e-09, 0);

    // radiation variables
    double phi = 1.243689808415833;
    int bornpdgs[] = {0, -1, -13, 14, -2};
    int realpdgs[] = {-5, -1, -13, 14, -2, -5};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 4.803600652401992e-09
    double fks_g = 7.566569953898706e-27;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, -1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, -1, -13, 14, -2
// realpdgs = -4, -1, -13, 14, -2, -4
TEST(QCDCollinearISR1, Wj_born_1_real_12_n1) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.297003726451937;
    ps.X2 = 0.4765416578967319;
    ps.Jacobian = 4134.54558701851;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2445.370909833025, 0, 0, 2445.370909833025);
    ps.Momenta[1].Set(2445.370909833025, 0, 0, -2445.370909833025);
    ps.Momenta[2].Set(1715.825694052779, 1627.159633788413, -540.465613731689,
                      -65.62209166802791);
    ps.Momenta[3].Set(1497.338444087164, -239.2870758019091, 1427.072543857699,
                      385.0039818794809);
    ps.Momenta[4].Set(1677.577681526107, -1387.872557986504, -886.60693012601,
                      -319.3818902114531);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 2.192884043273137e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(2.794384476375967e-10, 0);
    spin[1][2] =
        std::complex<double>(7.311818109174863e-10, 7.964060950108475e-12);
    spin[2][1] =
        std::complex<double>(7.311818109174863e-10, -7.964060950108475e-12);
    spin[2][2] = std::complex<double>(1.91344559563554e-09, 0);

    // radiation variables
    double phi = 4.791201102945314;
    double xi = 0.03969538821707171;
    int bornpdgs[] = {0, -1, -13, 14, -2};
    int realpdgs[] = {-4, -1, -13, 14, -2, -4};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 2.00028855295633e-06
    double fks_g = 6.303804739046175e-17;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, 1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 6.303839899118012e-17
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, -1, -13, 14, -2
// realpdgs = -4, -1, -13, 14, -2, -4
TEST(QCDSoftCollinearISR1, Wj_born_1_real_12_n1) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.4349942926847881;
    ps.X2 = 0.4903070064922233;
    ps.Jacobian = 4238.501938364773;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3001.851373043983, 0, 0, 3001.851373043983);
    ps.Momenta[1].Set(3001.851373043983, 0, 0, -3001.851373043983);
    ps.Momenta[2].Set(1611.514538055187, -387.670913543706, 477.8924898954706,
                      1489.398918106536);
    ps.Momenta[3].Set(2991.237688059122, 926.684295308171, -628.9995958599108,
                      -2773.647171451796);
    ps.Momenta[4].Set(1400.950519973658, -539.013381764465, 151.1071059644401,
                      1284.24825334526);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 7.988053740876809e-11;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(4.691268539290636e-11, 0);
    spin[1][2] =
        std::complex<double>(-3.31278067256338e-12, -3.918719119706726e-11);
    spin[2][1] =
        std::complex<double>(-3.31278067256338e-12, 3.918719119706726e-11);
    spin[2][2] = std::complex<double>(3.296785201586174e-11, 0);

    // radiation variables
    double phi = 2.089634384026967;
    int bornpdgs[] = {0, -1, -13, 14, -2};
    int realpdgs[] = {-4, -1, -13, 14, -2, -4};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 3.101739110551611e-05
    double fks_g = 1.980345774014328e-23;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, 1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, -1, -13, 14, -2
// realpdgs = -4, -1, -13, 14, -2, -4
TEST(QCDSoftLimit, Wj_born_1_real_12_0_n1) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.8682591701325215;
    ps.X2 = 0.9456862415258342;
    ps.Jacobian = 10392.54701934665;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(5889.949638263703, 0, 0, 5889.949638263703);
    ps.Momenta[1].Set(5889.949638263703, 0, 0, -5889.949638263703);
    ps.Momenta[2].Set(5859.540647549627, -1399.56096260971, -3979.784333947149,
                      4066.541819223568);
    ps.Momenta[3].Set(4169.665557096763, 923.6197950470831, 2560.708427996474,
                      -3158.45051869173);
    ps.Momenta[4].Set(1750.693071881022, 475.9411675626272, 1419.075905950674,
                      -908.0913005318382);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    // double born = 9.472165281196034e-11;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, 1.420824792179405e-10);
    ColorCorr.Set(0, 4, 1.420824792179405e-10);
    ColorCorr.Set(1, 0, 1.420824792179405e-10);
    ColorCorr.Set(1, 4, -1.578694213532672e-11);
    ColorCorr.Set(4, 0, 1.420824792179405e-10);
    ColorCorr.Set(4, 1, -1.578694213532672e-11);

    // radiation variables
    double y = 0.7095042901714876;
    double phi = 1.172908209236108;
    int bornpdgs[] = {0, -1, -13, 14, -2};
    int realpdgs[] = {-4, -1, -13, 14, -2, -4};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 1.219351269759469e-08
    double fks_g = 1.406419811836095e-26;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 0, alphas, ColorCorr, y, phi);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, -1, -13, 14, -2
// realpdgs = -4, -1, -13, 14, -2, -4
TEST(QCDCollinearISR2, Wj_born_1_real_12_n1) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.06140691786948693;
    ps.X2 = 0.9182233416206103;
    ps.Jacobian = 1655.840658642875;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1543.462814573249, 0, 0, 1543.462814573249);
    ps.Momenta[1].Set(1543.462814573249, 0, 0, -1543.462814573249);
    ps.Momenta[2].Set(1308.57049589061, 119.7969607860818, -404.7940075065546,
                      1238.606976562218);
    ps.Momenta[3].Set(713.9132175197337, 85.81143033588802, -406.3256285740266,
                      -580.6961030840571);
    ps.Momenta[4].Set(1064.441915736155, -205.6083911219698, 811.1196360805811,
                      -657.9108734781609);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 3.980571530761365e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(2.282609115123752e-09, 8.616464714094038e-27);
    spin[1][2] =
        std::complex<double>(1.240926129335541e-10, 1.964786377825153e-09);
    spin[2][1] =
        std::complex<double>(1.240926129335541e-10, -1.964786377825153e-09);
    spin[2][2] =
        std::complex<double>(1.697962415637612e-09, -8.616464714094038e-27);

    // radiation variables
    double phi = 2.876126274670849;
    double xi = 0.0695508424722252;
    int bornpdgs[] = {0, -1, -13, 14, -2};
    int realpdgs[] = {-4, -1, -13, 14, -2, -4};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 2.36596673956425e-14
    double fks_g = 2.288987497166799e-24;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, -1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, -1, -13, 14, -2
// realpdgs = -4, -1, -13, 14, -2, -4
TEST(QCDSoftCollinearISR2, Wj_born_1_real_12_n1) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.7104920998101427;
    ps.X2 = 0.5003225953914949;
    ps.Jacobian = 10586.56080051896;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3875.413445156067, 0, 0, 3875.413445156067);
    ps.Momenta[1].Set(3875.413445156067, 0, 0, -3875.413445156067);
    ps.Momenta[2].Set(3410.00932662997, -2600.714206573691, 1105.91462174612,
                      1908.24581038117);
    ps.Momenta[3].Set(1630.398452686997, 35.62713708952379, -1251.160112887037,
                      -1044.76226652165);
    ps.Momenta[4].Set(2710.419110995168, 2565.087069484167, 145.2454911409166,
                      -863.4835438595203);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 9.225931043611038e-10;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(5.509201199822063e-10, 4.308232357047019e-27);
    spin[1][2] =
        std::complex<double>(7.724982065320791e-11, 4.458638697598339e-10);
    spin[2][1] =
        std::complex<double>(7.724982065320791e-11, -4.458638697598339e-10);
    spin[2][2] =
        std::complex<double>(3.716729843788975e-10, -4.308232357047019e-27);

    // radiation variables
    double phi = 3.82372195374913;
    int bornpdgs[] = {0, -1, -13, 14, -2};
    int realpdgs[] = {-4, -1, -13, 14, -2, -4};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 1.215304450625031e-10
    double fks_g = 6.068685268696779e-29;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, -1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, -1, -13, 14, -2
// realpdgs = -3, -1, -13, 14, -2, -3
TEST(QCDCollinearISR1, Wj_born_1_real_12_n2) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.7643694716688536;
    ps.X2 = 0.2181973208902477;
    ps.Jacobian = 5721.274064212278;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2654.542789264566, 0, 0, 2654.542789264566);
    ps.Momenta[1].Set(2654.542789264566, 0, 0, -2654.542789264566);
    ps.Momenta[2].Set(1299.025160539692, -326.0362453325363, 1171.988806820726,
                      -455.6412746136651);
    ps.Momenta[3].Set(1871.593095436332, 1490.159230280158, 538.4804108205467,
                      996.1551236922554);
    ps.Momenta[4].Set(2138.467322553109, -1164.122984947622, -1710.469217641272,
                      -540.5138490785903);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 3.088516358547716e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(2.753951497201484e-09, -3.446585885637615e-26);
    spin[1][2] =
        std::complex<double>(4.736291270817815e-10, -8.348957125265478e-10);
    spin[2][1] =
        std::complex<double>(4.736291270817815e-10, 8.348957125265478e-10);
    spin[2][2] =
        std::complex<double>(3.345648613462316e-10, 3.446585885637615e-26);

    // radiation variables
    double phi = 5.46673528427572;
    double xi = 0.1044125762631013;
    int bornpdgs[] = {0, -1, -13, 14, -2};
    int realpdgs[] = {-3, -1, -13, 14, -2, -3};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 8.681594797121329e-07
    double fks_g = 1.89293251194153e-16;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, 1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 1.892930111408727e-16
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, -1, -13, 14, -2
// realpdgs = -3, -1, -13, 14, -2, -3
TEST(QCDSoftCollinearISR1, Wj_born_1_real_12_n2) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.6661010372116642;
    ps.X2 = 0.3061556205219347;
    ps.Jacobian = 7727.739275577053;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2935.313756988135, 0, 0, 2935.313756988135);
    ps.Momenta[1].Set(2935.313756988135, 0, 0, -2935.313756988135);
    ps.Momenta[2].Set(1929.595594393482, -481.036193231171, 1212.073066381535,
                      -1422.259547500687);
    ps.Momenta[3].Set(1328.885436542485, -773.1669592584399, -711.3536162443191,
                      -813.7108756935096);
    ps.Momenta[4].Set(2612.146483040304, 1254.203152489611, -500.7194501372157,
                      2235.970423194196);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 2.85193169811916e-08;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(1.592117899881565e-08, 0);
    spin[1][2] =
        std::complex<double>(4.13720840704411e-09, 1.354476712482977e-08);
    spin[2][1] =
        std::complex<double>(4.13720840704411e-09, -1.354476712482977e-08);
    spin[2][2] = std::complex<double>(1.259813798237595e-08, 0);

    // radiation variables
    double phi = 1.840853682522748;
    int bornpdgs[] = {0, -1, -13, 14, -2};
    int realpdgs[] = {-3, -1, -13, 14, -2, -3};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 0.01958741761593887
    double fks_g = 4.367545752263893e-21;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, 1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, -1, -13, 14, -2
// realpdgs = -3, -1, -13, 14, -2, -3
TEST(QCDSoftLimit, Wj_born_1_real_12_0_n2) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.6331890441503383;
    ps.X2 = 0.547850533401153;
    ps.Jacobian = 6096.592108307173;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3828.345252628033, 0, 0, 3828.345252628033);
    ps.Momenta[1].Set(3828.345252628033, 0, 0, -3828.345252628033);
    ps.Momenta[2].Set(2833.494449946549, -771.9657287199534, -1459.295766321729,
                      2302.87115965163);
    ps.Momenta[3].Set(3243.128559883165, -474.43604222521, 882.8137989549928,
                      -3084.385367320721);
    ps.Momenta[4].Set(1580.067495426354, 1246.401770945163, 576.4819673667367,
                      781.514207669091);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    // double born = 1.491330305183105e-10;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, 2.236995457774657e-10);
    ColorCorr.Set(0, 4, 2.236995457774657e-10);
    ColorCorr.Set(1, 0, 2.236995457774657e-10);
    ColorCorr.Set(1, 4, -2.485550508638508e-11);
    ColorCorr.Set(4, 0, 2.236995457774657e-10);
    ColorCorr.Set(4, 1, -2.485550508638508e-11);

    // radiation variables
    double y = 0.9380640107562144;
    double phi = 0.671043924366573;
    int bornpdgs[] = {0, -1, -13, 14, -2};
    int realpdgs[] = {-3, -1, -13, 14, -2, -3};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 8.631611875051199e-08
    double fks_g = 1.467226216645736e-25;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 0, alphas, ColorCorr, y, phi);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, -1, -13, 14, -2
// realpdgs = -3, -1, -13, 14, -2, -3
TEST(QCDCollinearISR2, Wj_born_1_real_12_n2) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.09064131419906474;
    ps.X2 = 0.4400251678888907;
    ps.Jacobian = 1052.693381877836;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1298.121109063126, 0, 0, 1298.121109063126);
    ps.Momenta[1].Set(1298.121109063126, 0, 0, -1298.121109063126);
    ps.Momenta[2].Set(635.7542916589318, 393.8437642527149, -141.1432864178303,
                      -478.6952907865046);
    ps.Momenta[3].Set(1155.876349808953, -1067.907645507949, 419.8640562281604,
                      139.0595951416996);
    ps.Momenta[4].Set(804.6115766583675, 674.0638812552343, -278.7207698103301,
                      339.6356956448051);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.284645867893325e-08;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(4.605325322902048e-09, -1.378634354255046e-25);
    spin[1][2] =
        std::complex<double>(1.985932535067196e-09, -5.831738342904696e-09);
    spin[2][1] =
        std::complex<double>(1.985932535067196e-09, 5.831738342904696e-09);
    spin[2][2] =
        std::complex<double>(8.241133356031202e-09, 1.378634354255046e-25);

    // radiation variables
    double phi = 2.437796066458127;
    double xi = 0.1024430952320856;
    int bornpdgs[] = {0, -1, -13, 14, -2};
    int realpdgs[] = {-3, -1, -13, 14, -2, -3};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 3.816224348684821e-14
    double fks_g = 8.009940264165552e-24;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, -1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, -1, -13, 14, -2
// realpdgs = -3, -1, -13, 14, -2, -3
TEST(QCDSoftCollinearISR2, Wj_born_1_real_12_n2) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.9453050075532552;
    ps.X2 = 0.7142889861784667;
    ps.Jacobian = 15486.61200558397;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(5341.16891408647, 0, 0, 5341.16891408647);
    ps.Momenta[1].Set(5341.16891408647, 0, 0, -5341.16891408647);
    ps.Momenta[2].Set(2751.328545509963, 2225.133037981836, -1255.040376628315,
                      1021.501532868329);
    ps.Momenta[3].Set(5054.142832434458, -3473.986594542189, 3670.042421662479,
                      -81.02798768068766);
    ps.Momenta[4].Set(2876.86645022852, 1248.853556560353, -2415.002045034164,
                      -940.4735451876415);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 2.626405599851841e-10;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(1.05856936893669e-10, 0);
    spin[1][2] =
        std::complex<double>(1.570731850617501e-11, -1.278667831804056e-10);
    spin[2][1] =
        std::complex<double>(1.570731850617501e-11, 1.278667831804056e-10);
    spin[2][2] = std::complex<double>(1.567836230915151e-10, 0);

    // radiation variables
    double phi = 5.593169353643097;
    int bornpdgs[] = {0, -1, -13, 14, -2};
    int realpdgs[] = {-3, -1, -13, 14, -2, -3};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 3.185387495186363e-11
    double fks_g = 5.200515392040829e-30;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, -1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, -1, -13, 14, -2
// realpdgs = 0, -1, -13, 14, -2, 0
TEST(QCDCollinearFSR, Wj_born_1_real_30) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.6314141956981167;
    ps.X2 = 0.3578484001275228;
    ps.Jacobian = 8970.473785942677;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3089.726711113615, 0, 0, 3089.726711113615);
    ps.Momenta[1].Set(3089.726711113615, 0, 0, -3089.726711113615);
    ps.Momenta[2].Set(892.539049986984, 566.9598846544266, -226.1604480346336,
                      -651.1788515370265);
    ps.Momenta[3].Set(2406.235082024925, 1077.079508659552, 1887.944169467961,
                      -1032.246973822217);
    ps.Momenta[4].Set(2880.679290215321, -1644.039393313979, -1661.783721433328,
                      1683.425825359244);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 4.44288292783178e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(2.512202409422824e-12, 3.446585885637615e-26);
    spin[1][2] =
        std::complex<double>(-1.050056773460562e-10, -1.135419571046337e-11);
    spin[2][1] =
        std::complex<double>(-1.050056773460562e-10, 1.135419571046337e-11);
    spin[2][2] =
        std::complex<double>(4.440370725422357e-09, -3.446585885637615e-26);

    // radiation variables
    double phi = 3.315387883576041;
    double xi = 0.2796824711440804;
    int bornpdgs[] = {0, -1, -13, 14, -2};
    int realpdgs[] = {0, -1, -13, 14, -2, 0};

    // FKS g = real matrix element * (1.0 - y) * xi * xi
    // real = 2.504587852673515e-06
    double fks_g = 1.95914584965402e-15;
    double limit = FKS::QCD::CollinearLimitFSR(
        realpdgs[4], bornpdgs[4], ps, 4, xi, phi, alphas, born,
        spin); // limit = 1.958569765663863e-15
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, -1, -13, 14, -2
// realpdgs = 0, -1, -13, 14, -2, 0
TEST(QCDSoftCollinearFSR, Wj_born_1_real_30) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.7979718303876773;
    ps.X2 = 0.6979032807407819;
    ps.Jacobian = 18072.83983636562;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(4850.703808827706, 0, 0, 4850.703808827706);
    ps.Momenta[1].Set(4850.703808827706, 0, 0, -4850.703808827706);
    ps.Momenta[2].Set(3778.72391928797, 3410.251851634262, -571.2990309867716,
                      1523.992842443068);
    ps.Momenta[3].Set(2225.923841983469, -785.9753771757864, 1971.411770831518,
                      671.2043553276196);
    ps.Momenta[4].Set(3696.759856383973, -2624.276474458476, -1400.112739844746,
                      -2195.197197770688);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 6.322488018099813e-10;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(1.964618689823696e-10, 0);
    spin[1][2] =
        std::complex<double>(2.392022413997489e-10, 1.68516476971076e-10);
    spin[2][1] =
        std::complex<double>(2.392022413997489e-10, -1.68516476971076e-10);
    spin[2][2] = std::complex<double>(4.357869328276117e-10, 0);

    // radiation variables
    double phi = 2.524119393613116;
    int bornpdgs[] = {0, -1, -13, 14, -2};
    int realpdgs[] = {0, -1, -13, 14, -2, 0};

    // FKS g = real matrix element * (1.0 - y) * xi * xi
    // real = 183.1684735255642
    double fks_g = 1.063858068162393e-16;
    double limit = FKS::QCD::SoftCollinearLimitFSR(realpdgs[4], bornpdgs[4], ps,
                                                   4, phi, alphas, born, spin);
    // limit = 1.062524872142025e-16
    EXPECT_NEAR(limit / fks_g, 1.0, 5e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, -1, -13, 14, -2
// realpdgs = 0, -1, -13, 14, -2, 0
TEST(QCDSoftLimit, Wj_born_1_real_30_4) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.01056946004455384;
    ps.X2 = 0.4779955034617416;
    ps.Jacobian = 97.70088326052321;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(462.0103054663069, 0, 0, 462.0103054663069);
    ps.Momenta[1].Set(462.0103054663069, 0, 0, -462.0103054663069);
    ps.Momenta[2].Set(397.977038807958, 75.49535489712585, 294.4314228000875,
                      256.8974738591721);
    ps.Momenta[3].Set(316.2237939980747, 14.22293893412434, -308.552590306199,
                      -67.75319117169144);
    ps.Momenta[4].Set(209.8197781265813, -89.71829383125021, 14.12116750611148,
                      -189.1442826874807);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    // double born = 1.975758722141095e-08;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, 2.963638083211642e-08);
    ColorCorr.Set(0, 4, 2.963638083211642e-08);
    ColorCorr.Set(1, 0, 2.963638083211642e-08);
    ColorCorr.Set(1, 4, -3.292931203568491e-09);
    ColorCorr.Set(4, 0, 2.963638083211642e-08);
    ColorCorr.Set(4, 1, -3.292931203568491e-09);

    // radiation variables
    double y = 0.9238316669310862;
    double phi = 0.4689923217761751;
    int bornpdgs[] = {0, -1, -13, 14, -2};
    int realpdgs[] = {0, -1, -13, 14, -2, 0};

    // FKS g = real matrix element * (1.0 - y) * xi * xi
    // real = 1143112.692016057
    double fks_g = 1.795779539625652e-12;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 4, alphas, ColorCorr, y, phi);
    // limit = 1.795779119040393e-12
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, -1, -13, 14, -2
// realpdgs = 0, -1, -13, 14, -2, 0
TEST(QCDCollinearISR1, Wj_born_1_real_30) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.1512161607565439;
    ps.X2 = 0.1738797477916805;
    ps.Jacobian = 358.9352508245257;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1053.991142532661, 0, 0, 1053.991142532661);
    ps.Momenta[1].Set(1053.991142532661, 0, 0, -1053.991142532661);
    ps.Momenta[2].Set(785.8782867698321, 472.6624506164883, 40.42594028045433,
                      -626.5465926377251);
    ps.Momenta[3].Set(984.211346313301, -471.0532632609103, 226.2666164861278,
                      834.0169156839731);
    ps.Momenta[4].Set(337.8926519821887, -1.609187355578001, -266.6925567665821,
                      -207.470323046248);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 2.164488839805205e-08;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(1.134896917860204e-08, -1.378634354255046e-25);
    spin[1][2] =
        std::complex<double>(-1.565467261351439e-09, -1.069567118700605e-08);
    spin[2][1] =
        std::complex<double>(-1.565467261351439e-09, 1.069567118700605e-08);
    spin[2][2] =
        std::complex<double>(1.029591921945002e-08, 1.378634354255046e-25);

    // radiation variables
    double phi = 3.242184825341755;
    double xi = 0.2758395497171194;
    int bornpdgs[] = {0, -1, -13, 14, -2};
    int realpdgs[] = {0, -1, -13, 14, -2, 0};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 0.0002019254561228895
    double fks_g = 3.072798897971943e-13;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, 1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 3.072812117906379e-13
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, -1, -13, 14, -2
// realpdgs = 0, -1, -13, 14, -2, 0
TEST(QCDSoftCollinearISR1, Wj_born_1_real_30) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.8975567276565393;
    ps.X2 = 0.916574169833126;
    ps.Jacobian = 29905.82036661989;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(5895.601449758076, 0, 0, 5895.601449758076);
    ps.Momenta[1].Set(5895.601449758076, 0, 0, -5895.601449758076);
    ps.Momenta[2].Set(1851.176854732907, 1155.048052551133, 1392.876545678302,
                      390.65927903473);
    ps.Momenta[3].Set(4907.022923166541, -3355.120287211727, 2670.636619808476,
                      -2385.317981267938);
    ps.Momenta[4].Set(5033.003121616705, 2200.072234660594, -4063.513165486778,
                      1994.658702233208);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 3.171344987900771e-10;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(5.817767497175427e-11, -8.616464714094038e-27);
    spin[1][2] =
        std::complex<double>(-1.167657470412806e-10, -3.783207958481634e-11);
    spin[2][1] =
        std::complex<double>(-1.167657470412806e-10, 3.783207958481634e-11);
    spin[2][2] =
        std::complex<double>(2.589568238183228e-10, 8.616464714094038e-27);

    // radiation variables
    double phi = 4.19298333132849;
    int bornpdgs[] = {0, -1, -13, 14, -2};
    int realpdgs[] = {0, -1, -13, 14, -2, 0};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 7736.65481009928
    double fks_g = 1.623866401311156e-16;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, 1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 1.62353020098213e-16
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, -1, -13, 14, -2
// realpdgs = 0, -1, -13, 14, -2, 0
TEST(QCDSoftLimit, Wj_born_1_real_30_0) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.6651231505959281;
    ps.X2 = 0.2486478532211223;
    ps.Jacobian = 6063.997636074371;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2643.362629088449, 0, 0, 2643.362629088449);
    ps.Momenta[1].Set(2643.362629088449, 0, 0, -2643.362629088449);
    ps.Momenta[2].Set(2250.711241036382, 1810.172739930966, -1004.400418646511,
                      -883.2641400855484);
    ps.Momenta[3].Set(759.8588027329421, -178.1837316304195, -537.042832697539,
                      507.1695512559895);
    ps.Momenta[4].Set(2276.155214407575, -1631.989008300547, 1541.44325134405,
                      376.094588829559);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    // double born = 9.140538918809305e-09;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, 1.371080837821396e-08);
    ColorCorr.Set(0, 4, 1.371080837821396e-08);
    ColorCorr.Set(1, 0, 1.371080837821396e-08);
    ColorCorr.Set(1, 4, -1.523423153134884e-09);
    ColorCorr.Set(4, 0, 1.371080837821396e-08);
    ColorCorr.Set(4, 1, -1.523423153134884e-09);

    // radiation variables
    double y = 0.5130396182974195;
    double phi = 6.232253920401152;
    int bornpdgs[] = {0, -1, -13, 14, -2};
    int realpdgs[] = {0, -1, -13, 14, -2, 0};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 1248.429670299144
    double fks_g = 1.611650152696661e-14;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 0, alphas, ColorCorr, y, phi);
    // limit = 1.611650156643007e-14
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, -1, -13, 14, -2
// realpdgs = 0, -1, -13, 14, -2, 0
TEST(QCDCollinearISR2, Wj_born_1_real_30) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.6169127149285778;
    ps.X2 = 0.09970071674136705;
    ps.Jacobian = 2384.403471392025;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1612.034594374906, 0, 0, 1612.034594374906);
    ps.Momenta[1].Set(1612.034594374906, 0, 0, -1612.034594374906);
    ps.Momenta[2].Set(1052.197960163238, 597.4354649734764, 669.2471289964839,
                      -549.8178724755754);
    ps.Momenta[3].Set(704.2804161426682, 332.2063619790241, -321.5853905107028,
                      -531.2557521879929);
    ps.Momenta[4].Set(1467.590812443906, -929.6418269525007, -347.661738485781,
                      1081.073624663569);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 6.536363564260115e-08;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(3.84845412322164e-08, 0);
    spin[1][2] =
        std::complex<double>(-2.106738243879361e-09, 3.209347708013055e-08);
    spin[2][1] =
        std::complex<double>(-2.106738243879361e-09, -3.209347708013055e-08);
    spin[2][2] = std::complex<double>(2.687909441038475e-08, 0);

    // radiation variables
    double phi = 3.672161958025353;
    double xi = 0.504890492412355;
    int bornpdgs[] = {0, -1, -13, 14, -2};
    int realpdgs[] = {0, -1, -13, 14, -2, 0};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 2.428798103466256e-05
    double fks_g = 1.238271267194259e-13;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, -1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 1.238408660404547e-13
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, -1, -13, 14, -2
// realpdgs = 0, -1, -13, 14, -2, 0
TEST(QCDSoftCollinearISR2, Wj_born_1_real_30) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.7356747476553629;
    ps.X2 = 0.2931963724412228;
    ps.Jacobian = 749.6784274987947;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3018.808592609375, 0, 0, 3018.808592609375);
    ps.Momenta[1].Set(3018.808592609375, 0, 0, -3018.808592609375);
    ps.Momenta[2].Set(2972.374771550083, -279.5034443642897, -2245.666457881347,
                      -1927.140826998003);
    ps.Momenta[3].Set(2818.843361526707, 144.6979003041179, 2039.823021302996,
                      1940.119186089629);
    ps.Momenta[4].Set(246.3990521419589, 134.8055440601718, 205.8434365783511,
                      -12.97835909162577);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 9.789604391065043e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(5.436480550330135e-09, 0);
    spin[1][2] =
        std::complex<double>(-6.115123967837757e-10, -4.82615019263942e-09);
    spin[2][1] =
        std::complex<double>(-6.115123967837757e-10, 4.82615019263942e-09);
    spin[2][2] = std::complex<double>(4.353123840734908e-09, 0);

    // radiation variables
    double phi = 0.2893190495448469;
    int bornpdgs[] = {0, -1, -13, 14, -2};
    int realpdgs[] = {0, -1, -13, 14, -2, 0};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 8502.08835535764
    double fks_g = 8.494801065558901e-15;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, -1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 8.495414749215598e-15
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, 0, -13, 14, 1
// realpdgs = 2, 1, -13, 14, 1, 1
TEST(QCDCollinearISR1, Wj_born_2_real_33) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.07758841763459046;
    ps.X2 = 0.7282667730951857;
    ps.Jacobian = 1109.818667595705;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1545.101634627275, 0, 0, 1545.101634627275);
    ps.Momenta[1].Set(1545.101634627275, 0, 0, -1545.101634627275);
    ps.Momenta[2].Set(1247.210285014386, 690.9430239198134, -510.3136376720612,
                      -904.2738655673504);
    ps.Momenta[3].Set(1130.312986471911, -1033.384131256337, 344.5961501255657,
                      301.6258907542172);
    ps.Momenta[4].Set(712.6799977682537, 342.4411073365232, 165.7174875464955,
                      602.6479748131333);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 5.374552926304765e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(1.402507604845269e-09, -6.89317177127523e-26);
    spin[1][2] =
        std::complex<double>(5.100381348585423e-10, -2.304492323946217e-09);
    spin[2][1] =
        std::complex<double>(5.100381348585423e-10, 2.304492323946217e-09);
    spin[2][2] =
        std::complex<double>(3.972045321459496e-09, 6.89317177127523e-26);

    // radiation variables
    double phi = 6.281457887187575;
    double xi = 0.2917956953978194;
    int bornpdgs[] = {2, 0, -13, 14, 1};
    int realpdgs[] = {2, 1, -13, 14, 1, 1};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 2.38440386019391e-15
    double fks_g = 4.060388353205403e-24;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, 1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, 0, -13, 14, 1
// realpdgs = 2, 1, -13, 14, 1, 1
TEST(QCDSoftCollinearISR1, Wj_born_2_real_33) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.2551774098233288;
    ps.X2 = 0.4841680334500964;
    ps.Jacobian = 4592.702949175039;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2284.717589411409, 0, 0, 2284.717589411409);
    ps.Momenta[1].Set(2284.717589411409, 0, 0, -2284.717589411409);
    ps.Momenta[2].Set(1542.207267276503, -1333.735126174661, 546.321265754896,
                      548.7138990697146);
    ps.Momenta[3].Set(1032.721540868111, -660.1344550611727, -526.1500373539178,
                      -594.8969830202235);
    ps.Momenta[4].Set(1994.506370678205, 1993.869581235834, -20.17122840097813,
                      46.18308395050892);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 2.917744945729427e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(1.960297577981062e-09, 0);
    spin[1][2] =
        std::complex<double>(1.251768756262474e-09, 5.567375835046188e-10);
    spin[2][1] =
        std::complex<double>(1.251768756262474e-09, -5.567375835046188e-10);
    spin[2][2] = std::complex<double>(9.574473677483649e-10, 0);

    // radiation variables
    double phi = 2.05975036774104;
    int bornpdgs[] = {2, 0, -13, 14, 1};
    int realpdgs[] = {2, 1, -13, 14, 1, 1};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 3.709386194633922e-10
    double fks_g = 4.115643821108984e-28;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, 1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(limit, limit, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, 0, -13, 14, 1
// realpdgs = 2, 1, -13, 14, 1, 1
TEST(QCDSoftLimit, Wj_born_2_real_33_0) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.4576223094267;
    ps.X2 = 0.4698841959449425;
    ps.Jacobian = 1203.836890489981;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3014.132709587955, 0, 0, 3014.132709587955);
    ps.Momenta[1].Set(3014.132709587955, 0, 0, -3014.132709587955);
    ps.Momenta[2].Set(2789.016662632833, -1130.110450312729, 393.0222960629792,
                      -2519.324867764165);
    ps.Momenta[3].Set(2842.966274771802, 1390.036503943381, -96.74985701315913,
                      2478.082973260309);
    ps.Momenta[4].Set(396.2824817712752, -259.9260536306521, -296.27243904982,
                      41.24189450385583);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    // double born = 7.672096590598254e-09;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, 1.150814488589738e-08);
    ColorCorr.Set(0, 4, -1.278682765099709e-09);
    ColorCorr.Set(1, 0, 1.150814488589738e-08);
    ColorCorr.Set(1, 4, 1.150814488589738e-08);
    ColorCorr.Set(4, 0, -1.278682765099709e-09);
    ColorCorr.Set(4, 1, 1.150814488589738e-08);

    // radiation variables
    double y = -0.08250042069035857;
    double phi = 0.5671514978987908;
    int bornpdgs[] = {2, 0, -13, 14, 1};
    int realpdgs[] = {2, 1, -13, 14, 1, 1};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 1.203695056149864e-07
    double fks_g = 6.831019800522025e-24;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 0, alphas, ColorCorr, y, phi);
    // limit = 0
    EXPECT_NEAR(limit, limit, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, 0, -13, 14, 1
// realpdgs = 2, 1, -13, 14, 1, 1
TEST(QCDCollinearISR2, Wj_born_2_real_33) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.7256612093424586;
    ps.X2 = 0.8641161144178326;
    ps.Jacobian = 21974.44247894028;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(5147.14452481973, 0, 0, 5147.14452481973);
    ps.Momenta[1].Set(5147.14452481973, 0, 0, -5147.14452481973);
    ps.Momenta[2].Set(4152.581582854026, 2540.249940523745, -2604.575049452682,
                      -2001.812392232626);
    ps.Momenta[3].Set(1905.754797576284, -1671.222088971238, -884.6748102827214,
                      -237.2099447132537);
    ps.Momenta[4].Set(4235.952669209148, -869.0278515525072, 3489.249859735405,
                      2239.022336945881);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 7.530380207277413e-10;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(2.209850158040251e-10, 0);
    spin[1][2] =
        std::complex<double>(-1.561092902877277e-10, -3.052959730941573e-10);
    spin[2][1] =
        std::complex<double>(-1.561092902877277e-10, 3.052959730941573e-10);
    spin[2][2] = std::complex<double>(5.320530049237161e-10, 0);

    // radiation variables
    double phi = 1.666160651467829;
    double xi = 0.002086610016152125;
    int bornpdgs[] = {2, 0, -13, 14, 1};
    int realpdgs[] = {2, 1, -13, 14, 1, 1};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 1.352153093914062e-06
    double fks_g = 1.177439055375636e-19;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, -1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 1.177448120052351e-19
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, 0, -13, 14, 1
// realpdgs = 2, 1, -13, 14, 1, 1
TEST(QCDSoftCollinearISR2, Wj_born_2_real_33) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.8768982887213532;
    ps.X2 = 0.09064635782166874;
    ps.Jacobian = 340.5777028835869;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1832.580864034175, 0, 0, 1832.580864034175);
    ps.Momenta[1].Set(1832.580864034175, 0, 0, -1832.580864034175);
    ps.Momenta[2].Set(1790.417250474345, -1572.793088485339, 281.1418788698186,
                      808.006853655592);
    ps.Momenta[3].Set(1690.347982897397, 1513.941398462505, -119.2016292751949,
                      -742.3265567692102);
    ps.Momenta[4].Set(184.3964946966076, 58.85169002283416, -161.9402495946237,
                      -65.68029688638175);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 4.164608487211801e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(2.073430921031457e-09, -3.446585885637615e-26);
    spin[1][2] =
        std::complex<double>(-3.547521899490839e-10, -2.051843831974539e-09);
    spin[2][1] =
        std::complex<double>(-3.547521899490839e-10, 2.051843831974539e-09);
    spin[2][2] =
        std::complex<double>(2.091177566180344e-09, 3.446585885637615e-26);

    // radiation variables
    double phi = 3.165543581789232;
    int bornpdgs[] = {2, 0, -13, 14, 1};
    int realpdgs[] = {2, 1, -13, 14, 1, 1};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 0.001348193090835021
    double fks_g = 2.229706673533882e-21;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, -1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, 0, -13, 14, 1
// realpdgs = 2, -2, -13, 14, 1, -2
TEST(QCDCollinearISR1, Wj_born_2_real_10) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.9570144050547178;
    ps.X2 = 0.639628411520313;
    ps.Jacobian = 24292.70340703105;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(5085.532888167127, 0, 0, 5085.532888167127);
    ps.Momenta[1].Set(5085.532888167127, 0, 0, -5085.532888167127);
    ps.Momenta[2].Set(4632.596797865876, 4176.532107424096, -1285.997530729954,
                      -1537.446909070118);
    ps.Momenta[3].Set(798.8985950164358, 25.49877208262342, -636.0154930933584,
                      482.7764185280032);
    ps.Momenta[4].Set(4739.570383451943, -4202.030879506719, 1922.013023823313,
                      1054.670490542115);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 2.113078337061546e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(1.279399727613017e-09, 0);
    spin[1][2] =
        std::complex<double>(2.063045509691856e-10, -1.011951885266611e-09);
    spin[2][1] =
        std::complex<double>(2.063045509691856e-10, 1.011951885266611e-09);
    spin[2][2] = std::complex<double>(8.336786094485288e-10, 0);

    // radiation variables
    double phi = 3.611714834362836;
    double xi = 0.01897425061481637;
    int bornpdgs[] = {2, 0, -13, 14, 1};
    int realpdgs[] = {2, -2, -13, 14, 1, -2};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 4.55488462646533e-15
    double fks_g = 3.279719042257353e-26;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, 1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, 0, -13, 14, 1
// realpdgs = 2, -2, -13, 14, 1, -2
TEST(QCDSoftCollinearISR1, Wj_born_2_real_10) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.3505606338757339;
    ps.X2 = 0.2592452200420183;
    ps.Jacobian = 2705.010252432974;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1959.522741942161, 0, 0, 1959.522741942161);
    ps.Momenta[1].Set(1959.522741942161, 0, 0, -1959.522741942161);
    ps.Momenta[2].Set(598.6117302639652, 436.0201874817323, 196.9218805065208,
                      359.7835080931152);
    ps.Momenta[3].Set(1950.756625372606, -1254.296964707357, -403.1733660182417,
                      -1438.624959012788);
    ps.Momenta[4].Set(1369.677128247753, 818.2767772256249, 206.251485511721,
                      1078.841450919673);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 3.795779934211197e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(1.666795877911227e-09, 0);
    spin[1][2] =
        std::complex<double>(-2.363008843509199e-10, 1.868888370458355e-09);
    spin[2][1] =
        std::complex<double>(-2.363008843509199e-10, -1.868888370458355e-09);
    spin[2][2] = std::complex<double>(2.128984056299969e-09, 0);

    // radiation variables
    double phi = 0.99882179824752;
    int bornpdgs[] = {2, 0, -13, 14, 1};
    int realpdgs[] = {2, -2, -13, 14, 1, -2};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 1.504743489464664e-09
    double fks_g = 1.26931603059336e-27;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, 1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, 0, -13, 14, 1
// realpdgs = 2, -2, -13, 14, 1, -2
TEST(QCDSoftLimit, Wj_born_2_real_10_0) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.09717008156877327;
    ps.X2 = 0.3636706297094676;
    ps.Jacobian = 663.4322094062719;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1221.894625495904, 0, 0, 1221.894625495904);
    ps.Momenta[1].Set(1221.894625495904, 0, 0, -1221.894625495904);
    ps.Momenta[2].Set(1202.151457369618, 502.9119077552094, -1089.52294197829,
                      -72.02429032542943);
    ps.Momenta[3].Set(702.9186677204802, -298.0868994935503, 626.6027114397415,
                      -112.284886768526);
    ps.Momenta[4].Set(538.7191259017104, -204.8250082616591, 462.9202305385488,
                      184.3091770939554);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    // double born = 6.219873032043599e-09;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, 9.329809548065399e-09);
    ColorCorr.Set(0, 4, -1.0366455053406e-09);
    ColorCorr.Set(1, 0, 9.329809548065399e-09);
    ColorCorr.Set(1, 4, 9.329809548065399e-09);
    ColorCorr.Set(4, 0, -1.0366455053406e-09);
    ColorCorr.Set(4, 1, 9.329809548065399e-09);

    // radiation variables
    double y = 0.6868759154130473;
    double phi = 3.271674054060605;
    int bornpdgs[] = {2, 0, -13, 14, 1};
    int realpdgs[] = {2, -2, -13, 14, 1, -2};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 5.095895388204703e-07
    double fks_g = 2.471166612425786e-23;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 0, alphas, ColorCorr, y, phi);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, 0, -13, 14, 1
// realpdgs = 2, -2, -13, 14, 1, -2
TEST(QCDCollinearISR2, Wj_born_2_real_10) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.5615879445749989;
    ps.X2 = 0.4041197462546506;
    ps.Jacobian = 2189.500800100292;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3096.544179596142, 0, 0, 3096.544179596142);
    ps.Momenta[1].Set(3096.544179596142, 0, 0, -3096.544179596142);
    ps.Momenta[2].Set(2923.365731111519, 884.6463534495589, -1232.274454555497,
                      -2498.993336488823);
    ps.Momenta[3].Set(2568.158420631916, -690.7234617167079, 1567.098829904573,
                      1913.776379363174);
    ps.Momenta[4].Set(701.5642074488507, -193.9228917328509, -334.8243753490759,
                      585.2169571256493);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 2.810284378719605e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(1.583688723575221e-09, 0);
    spin[1][2] =
        std::complex<double>(9.121894122697477e-12, -1.393722532818543e-09);
    spin[2][1] =
        std::complex<double>(9.121894122697477e-12, 1.393722532818543e-09);
    spin[2][2] = std::complex<double>(1.226595655144384e-09, 0);

    // radiation variables
    double phi = 1.406800147439511;
    double xi = 0.2404867476163311;
    int bornpdgs[] = {2, 0, -13, 14, 1};
    int realpdgs[] = {2, -2, -13, 14, 1, -2};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 3.175732002285324e-07
    double fks_g = 3.673297800622966e-16;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, -1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 3.6731479455468e-16
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, 0, -13, 14, 1
// realpdgs = 2, -2, -13, 14, 1, -2
TEST(QCDSoftCollinearISR2, Wj_born_2_real_10) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.7211421976499643;
    ps.X2 = 0.2380741556027104;
    ps.Jacobian = 2608.10205090858;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2693.270272455713, 0, 0, 2693.270272455713);
    ps.Momenta[1].Set(2693.270272455713, 0, 0, -2693.270272455713);
    ps.Momenta[2].Set(1789.351161771322, 1108.061473422694, 490.2098725779096,
                      -1316.689649110926);
    ps.Momenta[3].Set(2636.364504888947, -1270.691748582451, -462.0971737778959,
                      2263.233634583336);
    ps.Momenta[4].Set(960.8248782511575, 162.6302751597562, -28.11269880001378,
                      -946.5439854724098);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 2.024993418878227e-07;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(1.72997758705216e-08, 4.411629933616147e-24);
    spin[1][2] =
        std::complex<double>(2.651606231074127e-08, -5.000809357444117e-08);
    spin[2][1] =
        std::complex<double>(2.651606231074127e-08, 5.000809357444117e-08);
    spin[2][2] =
        std::complex<double>(1.851995660173011e-07, -4.411629933616147e-24);

    // radiation variables
    double phi = 1.17250120764826;
    int bornpdgs[] = {2, 0, -13, 14, 1};
    int realpdgs[] = {2, -2, -13, 14, 1, -2};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 0.07242470429827345
    double fks_g = 8.408958092849543e-20;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, -1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, 0, -13, 14, 1
// realpdgs = 2, -1, -13, 14, 1, -1
TEST(QCDCollinearISR1, Wj_born_2_real_20) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.513874269274698;
    ps.X2 = 0.1800962752008246;
    ps.Jacobian = 6.374914179401963;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1977.398307576666, 0, 0, 1977.398307576666);
    ps.Momenta[1].Set(1977.398307576666, 0, 0, -1977.398307576666);
    ps.Momenta[2].Set(1976.682425953547, -1844.691670436846, -590.4707464634461,
                      394.6268511713582);
    ps.Momenta[3].Set(1974.915443015329, 1843.241309124744, 590.9579782755083,
                      -391.8177526721519);
    ps.Momenta[4].Set(3.19874618445585, 1.450361312102103, -0.4872318120622718,
                      -2.809098499206257);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 2.341689417108779e-06;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(1.169943923915003e-06, 1.378634354255046e-25);
    spin[1][2] =
        std::complex<double>(1.91826208227952e-09, -1.17084279064588e-06);
    spin[2][1] =
        std::complex<double>(1.91826208227952e-09, 1.17084279064588e-06);
    spin[2][2] =
        std::complex<double>(1.171745493193776e-06, -1.378634354255046e-25);

    // radiation variables
    double phi = 2.044818826838902;
    double xi = 0.4509060629924296;
    int bornpdgs[] = {2, 0, -13, 14, 1};
    int realpdgs[] = {2, -1, -13, 14, 1, -1};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 4.993200424563512e-13
    double fks_g = 2.030397846630014e-21;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, 1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, 0, -13, 14, 1
// realpdgs = 2, -1, -13, 14, 1, -1
TEST(QCDSoftCollinearISR1, Wj_born_2_real_20) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.09150016469774513;
    ps.X2 = 0.2433189101996831;
    ps.Jacobian = 119.300175200509;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(969.8670966157708, 0, 0, 969.8670966157708);
    ps.Momenta[1].Set(969.8670966157708, 0, 0, -969.8670966157708);
    ps.Momenta[2].Set(884.876166739425, 265.3910276653644, -679.8046798806457,
                      -500.4388375321736);
    ps.Momenta[3].Set(932.8106489680791, -215.421462027748, 783.0817637576703,
                      458.8161416038593);
    ps.Momenta[4].Set(122.0473775240376, -49.96956563761643, -103.2770838770247,
                      41.62269592831431);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 4.212854792926005e-08;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(2.605324750716318e-08, -1.378634354255046e-25);
    spin[1][2] =
        std::complex<double>(-1.02217178530161e-09, -2.043939689644066e-08);
    spin[2][1] =
        std::complex<double>(-1.02217178530161e-09, 2.043939689644066e-08);
    spin[2][2] =
        std::complex<double>(1.607530042209687e-08, 1.378634354255046e-25);

    // radiation variables
    double phi = 3.778923968502202;
    int bornpdgs[] = {2, 0, -13, 14, 1};
    int realpdgs[] = {2, -1, -13, 14, 1, -1};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 5.703214169686974e-05
    double fks_g = 9.414546440578125e-23;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, 1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, 0, -13, 14, 1
// realpdgs = 2, -1, -13, 14, 1, -1
TEST(QCDSoftLimit, Wj_born_2_real_20_0) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.5333736568728149;
    ps.X2 = 0.9069545246404971;
    ps.Jacobian = 2089.904625543698;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(4520.868696688702, 0, 0, 4520.868696688702);
    ps.Momenta[1].Set(4520.868696688702, 0, 0, -4520.868696688702);
    ps.Momenta[2].Set(4369.043256763806, -1312.188484430708, 534.8965486870966,
                      -4132.866564867944);
    ps.Momenta[3].Set(4214.020086805363, 1096.552052842391, -135.5339589292252,
                      4066.591869536848);
    ps.Momenta[4].Set(458.6740498082341, 215.636431588317, -399.3625897578713,
                      66.27469533109672);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    // double born = 4.670020636281969e-09;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, 7.005030954422954e-09);
    ColorCorr.Set(0, 4, -7.783367727136615e-10);
    ColorCorr.Set(1, 0, 7.005030954422954e-09);
    ColorCorr.Set(1, 4, 7.005030954422954e-09);
    ColorCorr.Set(4, 0, -7.783367727136615e-10);
    ColorCorr.Set(4, 1, 7.005030954422954e-09);

    // radiation variables
    double y = 0.141656823646322;
    double phi = 1.349074644577219;
    int bornpdgs[] = {2, 0, -13, 14, 1};
    int realpdgs[] = {2, -1, -13, 14, 1, -1};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 5.882496782223445e-05
    double fks_g = 2.38486889721045e-22;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 0, alphas, ColorCorr, y, phi);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                   << ", fks_g = " << fks_g;
}

// bornpdgs = 2, 0, -13, 14, 1
// realpdgs = 2, -1, -13, 14, 1, -1
TEST(QCDCollinearISR2, Wj_born_2_real_20) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.07798170239616908;
    ps.X2 = 0.2034921121910784;
    ps.Jacobian = 158.5058665936186;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(818.8106871023481, 0, 0, 818.8106871023481);
    ps.Momenta[1].Set(818.8106871023481, 0, 0, -818.8106871023481);
    ps.Momenta[2].Set(670.5278424807741, -236.7903336581669, -469.3957242396477,
                      416.1797442133775);
    ps.Momenta[3].Set(775.0226921037989, 124.0083647899126, 622.1270003326956,
                      -445.2416132786713);
    ps.Momenta[4].Set(192.070839620123, 112.7819688682542, -152.731276093048,
                      29.06186906529373);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 4.648264976924459e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(4.497960523844861e-09, 0);
    spin[1][2] =
        std::complex<double>(-7.420635408845817e-10, 3.541259631825389e-10);
    spin[2][1] =
        std::complex<double>(-7.420635408845817e-10, -3.541259631825389e-10);
    spin[2][2] = std::complex<double>(1.503044530795982e-10, 0);

    // radiation variables
    double phi = 2.193732841825843;
    double xi = 0.4733837199970981;
    int bornpdgs[] = {2, 0, -13, 14, 1};
    int realpdgs[] = {2, -1, -13, 14, 1, -1};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 6.752177739905208e-06
    double fks_g = 3.026220003063195e-14;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, -1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 3.025906200530137e-14
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, 0, -13, 14, 1
// realpdgs = 2, -1, -13, 14, 1, -1
TEST(QCDSoftCollinearISR2, Wj_born_2_real_20) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.7465045380221973;
    ps.X2 = 0.8341780021549852;
    ps.Jacobian = 26172.71069382949;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(5129.310022738457, 0, 0, 5129.310022738457);
    ps.Momenta[1].Set(5129.310022738457, 0, 0, -5129.310022738457);
    ps.Momenta[2].Set(744.6207341458374, 607.777158044571, -429.0977381726943,
                      -30.69356568123155);
    ps.Momenta[3].Set(4451.215860776849, 1383.125050443781, -4206.897642958047,
                      449.7776736792271);
    ps.Momenta[4].Set(5062.783450554231, -1990.90220848835, 4635.995381130738,
                      -419.0841079979952);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 2.706726382882004e-08;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(1.405952691181389e-08, -1.378634354255046e-25);
    spin[1][2] =
        std::complex<double>(-2.315529225586444e-09, 1.332369887275669e-08);
    spin[2][1] =
        std::complex<double>(-2.315529225586444e-09, -1.332369887275669e-08);
    spin[2][2] =
        std::complex<double>(1.300773691700615e-08, 1.378634354255046e-25);

    // radiation variables
    double phi = 5.306673450936237;
    int bornpdgs[] = {2, 0, -13, 14, 1};
    int realpdgs[] = {2, -1, -13, 14, 1, -1};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 0.0122688784036898
    double fks_g = 6.747133846677047e-22;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, -1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, 0, -13, 14, 1
// realpdgs = 2, 3, -13, 14, 1, 3
TEST(QCDCollinearISR1, Wj_born_2_real_27) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.06200600499994025;
    ps.X2 = 0.3505148727122034;
    ps.Jacobian = 76.17303288121829;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(958.2602144695518, 0, 0, 958.2602144695518);
    ps.Momenta[1].Set(958.2602144695518, 0, 0, -958.2602144695518);
    ps.Momenta[2].Set(891.516753943382, 559.8746063691214, -692.8734352535889,
                      35.62233042147643);
    ps.Momenta[3].Set(946.1326678123534, -638.5391222553529, 697.0464528582236,
                      -39.51021398135524);
    ps.Momenta[4].Set(78.87100718336835, 78.66451588623156, -4.173017604634571,
                      3.887883559878818);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 3.169277827470145e-08;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(1.317421675298259e-08, 2.757268708510092e-25);
    spin[1][2] =
        std::complex<double>(4.018065461065332e-09, -1.509379652013241e-08);
    spin[2][1] =
        std::complex<double>(4.018065461065332e-09, 1.509379652013241e-08);
    spin[2][2] =
        std::complex<double>(1.851856152171886e-08, -2.757268708510092e-25);

    // radiation variables
    double phi = 1.361240751104167;
    double xi = 0.7142612989418797;
    int bornpdgs[] = {2, 0, -13, 14, 1};
    int realpdgs[] = {2, 3, -13, 14, 1, 3};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 8.517724232540583e-16
    double fks_g = 8.690961164433785e-24;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, 1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, 0, -13, 14, 1
// realpdgs = 2, 3, -13, 14, 1, 3
TEST(QCDSoftCollinearISR1, Wj_born_2_real_27) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.4697779501367201;
    ps.X2 = 0.0226138155782678;
    ps.Jacobian = 451.1119849263805;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(669.9564828562958, 0, 0, 669.9564828562958);
    ps.Momenta[1].Set(669.9564828562958, 0, 0, -669.9564828562958);
    ps.Momenta[2].Set(177.5935268665352, 60.87525848398278, 38.75826732894026,
                      162.269714990387);
    ps.Momenta[3].Set(494.2257521018409, 274.3711053692979, 89.43220011210079,
                      401.2249645300864);
    ps.Momenta[4].Set(668.0936867442163, -335.2463638532803, -128.1904674410409,
                      -563.4946795204728);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 0.0003526664174311941;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(0.0001789955738969037, 0);
    spin[1][2] =
        std::complex<double>(-4.797455052561024e-06, 0.0001762478275968199);
    spin[2][1] =
        std::complex<double>(-4.797455052561024e-06, -0.0001762478275968199);
    spin[2][2] = std::complex<double>(0.0001736708435342904, 0);

    // radiation variables
    double phi = 5.583989594008592;
    int bornpdgs[] = {2, 0, -13, 14, 1};
    int realpdgs[] = {2, 3, -13, 14, 1, 3};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 0.001464919584276543
    double fks_g = 8.236817650221471e-22;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, 1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, 0, -13, 14, 1
// realpdgs = 2, 3, -13, 14, 1, 3
TEST(QCDSoftLimit, Wj_born_2_real_27_0) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.0671210107517588;
    ps.X2 = 0.9175221520130847;
    ps.Jacobian = 710.6763988401185;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1613.061329034932, 0, 0, 1613.061329034932);
    ps.Momenta[1].Set(1613.061329034932, 0, 0, -1613.061329034932);
    ps.Momenta[2].Set(1514.972025846343, -1304.854208442879, -36.01716277317736,
                      768.894334607591);
    ps.Momenta[3].Set(1274.010570038201, 1222.432285750048, 100.5185205249038,
                      -344.4680919279153);
    ps.Momenta[4].Set(437.1400621853207, 82.42192269283059, -64.50135775172649,
                      -424.4262426796757);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    // double born = 4.308884105271926e-08;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, 6.463326157907889e-08);
    ColorCorr.Set(0, 4, -7.181473508786543e-09);
    ColorCorr.Set(1, 0, 6.463326157907889e-08);
    ColorCorr.Set(1, 4, 6.463326157907889e-08);
    ColorCorr.Set(4, 0, -7.181473508786543e-09);
    ColorCorr.Set(4, 1, 6.463326157907889e-08);

    // radiation variables
    double y = 0.7795370323097686;
    double phi = 1.363887548596495;
    int bornpdgs[] = {2, 0, -13, 14, 1};
    int realpdgs[] = {2, 3, -13, 14, 1, 3};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 3.534210262871851e-06
    double fks_g = 3.757699473392314e-23;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 0, alphas, ColorCorr, y, phi);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, 0, -13, 14, 1
// realpdgs = 2, 3, -13, 14, 1, 3
TEST(QCDCollinearISR2, Wj_born_2_real_27) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.5458289826933544;
    ps.X2 = 0.7267882686312674;
    ps.Jacobian = 7767.237301207409;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(4093.978966719885, 0, 0, 4093.978966719885);
    ps.Momenta[1].Set(4093.978966719885, 0, 0, -4093.978966719885);
    ps.Momenta[2].Set(3889.03635908096, 3776.62224394592, 374.500881695833,
                      -849.3982095552196);
    ps.Momenta[3].Set(2416.484086443811, -2340.185084829986, -481.1202833199352,
                      -362.5636244214131);
    ps.Momenta[4].Set(1882.437487915, -1436.437159115933, 106.6194016241022,
                      1211.961833976633);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 5.547239651337023e-10;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(2.396563965613554e-10, 2.154116178523509e-27);
    spin[1][2] =
        std::complex<double>(-2.909915648838428e-11, -2.732420122328572e-10);
    spin[2][1] =
        std::complex<double>(-2.909915648838428e-11, 2.732420122328572e-10);
    spin[2][2] =
        std::complex<double>(3.150675685723469e-10, -2.154116178523509e-27);

    // radiation variables
    double phi = 2.47717784962372;
    double xi = 0.08808034849971352;
    int bornpdgs[] = {2, 0, -13, 14, 1};
    int realpdgs[] = {2, 3, -13, 14, 1, 3};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 8.309269005976028e-08
    double fks_g = 1.289290739130226e-17;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, -1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 1.289297952118176e-17
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, 0, -13, 14, 1
// realpdgs = 2, 3, -13, 14, 1, 3
TEST(QCDSoftCollinearISR2, Wj_born_2_real_27) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.3537826111161415;
    ps.X2 = 0.7988843001673089;
    ps.Jacobian = 5795.089518546054;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3455.600604601791, 0, 0, 3455.600604601791);
    ps.Momenta[1].Set(3455.600604601791, 0, 0, -3455.600604601791);
    ps.Momenta[2].Set(2512.656340042298, -1811.874998444238, -1021.928594468918,
                      -1409.330628693717);
    ps.Momenta[3].Set(2734.610454890277, 1474.764132740705, -358.6772900247103,
                      2274.756183505672);
    ps.Momenta[4].Set(1663.934414271008, 337.110865703533, 1380.605884493629,
                      -865.4255548119555);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 3.118652240883451e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(2.77552338920454e-09, 6.89317177127523e-26);
    spin[1][2] =
        std::complex<double>(6.598174638180063e-11, -9.736573126568304e-10);
    spin[2][1] =
        std::complex<double>(6.598174638180063e-11, 9.736573126568304e-10);
    spin[2][2] =
        std::complex<double>(3.431288516789109e-10, -6.89317177127523e-26);

    // radiation variables
    double phi = 4.943080590020089;
    int bornpdgs[] = {2, 0, -13, 14, 1};
    int realpdgs[] = {2, 3, -13, 14, 1, 3};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 0.002568767871483691
    double fks_g = 2.078006869720368e-22;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, -1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, 0, -13, 14, 1
// realpdgs = 2, 4, -13, 14, 1, 4
TEST(QCDCollinearISR1, Wj_born_2_real_27_n1) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.7276392483873941;
    ps.X2 = 0.3262402584369681;
    ps.Jacobian = 6616.932868440826;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3166.94259416044, 0, 0, 3166.94259416044);
    ps.Momenta[1].Set(3166.94259416044, 0, 0, -3166.94259416044);
    ps.Momenta[2].Set(2998.85688574049, -1867.05530714683, 2289.056926537426,
                      517.1706568232955);
    ps.Momenta[3].Set(1261.947958936347, 472.7977505583833, -1086.145474848767,
                      435.0436134413308);
    ps.Momenta[4].Set(2073.080343644042, 1394.257556588447, -1202.911451688659,
                      -952.2142702646266);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.435411833086016e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(5.596207946697954e-10, 8.616464714094038e-27);
    spin[1][2] =
        std::complex<double>(1.070822428542782e-10, -6.91841217439716e-10);
    spin[2][1] =
        std::complex<double>(1.070822428542782e-10, 6.91841217439716e-10);
    spin[2][2] =
        std::complex<double>(8.75791038416221e-10, -8.616464714094038e-27);

    // radiation variables
    double phi = 2.872984179582253;
    double xi = 0.118503254783488;
    int bornpdgs[] = {2, 0, -13, 14, 1};
    int realpdgs[] = {2, 4, -13, 14, 1, 4};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 6.457489489333408e-16
    double fks_g = 1.813653260086353e-25;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, 1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, 0, -13, 14, 1
// realpdgs = 2, 4, -13, 14, 1, 4
TEST(QCDSoftCollinearISR1, Wj_born_2_real_27_n1) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.4418277783606008;
    ps.X2 = 0.2241217401234294;
    ps.Jacobian = 324.1851755025829;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2045.4169854859, 0, 0, 2045.4169854859);
    ps.Momenta[1].Set(2045.4169854859, 0, 0, -2045.4169854859);
    ps.Momenta[2].Set(2021.098628596167, -908.1115930625563, -316.4164279037133,
                      1777.653972293637);
    ps.Momenta[3].Set(1912.478015289979, 807.1214614756066, 396.7804549738217,
                      -1687.806972359027);
    ps.Momenta[4].Set(157.2573270856539, 100.9901315869497, -80.36402707010842,
                      -89.84699993460991);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 3.256941637138867e-10;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(1.697083670660473e-10, 0);
    spin[1][2] =
        std::complex<double>(-6.219415430375467e-11, -1.503462071521662e-10);
    spin[2][1] =
        std::complex<double>(-6.219415430375467e-11, 1.503462071521662e-10);
    spin[2][2] = std::complex<double>(1.559857966478394e-10, 0);

    // radiation variables
    double phi = 1.464828631390083;
    int bornpdgs[] = {2, 0, -13, 14, 1};
    int realpdgs[] = {2, 4, -13, 14, 1, 4};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 1.378645067503217e-10
    double fks_g = 8.590511065591535e-29;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, 1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, 0, -13, 14, 1
// realpdgs = 2, 4, -13, 14, 1, 4
TEST(QCDSoftLimit, Wj_born_2_real_27_0_n1) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.6793099935594178;
    ps.X2 = 0.4586084234789922;
    ps.Jacobian = 5974.493043904044;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3628.009137211748, 0, 0, 3628.009137211748);
    ps.Momenta[1].Set(3628.009137211748, 0, 0, -3628.009137211748);
    ps.Momenta[2].Set(3044.605006820753, 2357.573623341302, 687.1908162019421,
                      -1799.787498626176);
    ps.Momenta[3].Set(2577.487642057629, -2475.731787994221, 405.9160400272714,
                      591.1233604953504);
    ps.Momenta[4].Set(1633.925625545114, 118.1581646529182, -1093.106856229213,
                      1208.664138130825);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 9.361519849356805e-10;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, 1.404227977403521e-09);
    ColorCorr.Set(0, 4, -1.560253308226134e-10);
    ColorCorr.Set(1, 0, 1.404227977403521e-09);
    ColorCorr.Set(1, 4, 1.404227977403521e-09);
    ColorCorr.Set(4, 0, -1.560253308226134e-10);
    ColorCorr.Set(4, 1, 1.404227977403521e-09);

    // radiation variables
    double y = -0.9399684909423272;
    double phi = 2.770055682743375;
    int bornpdgs[] = {2, 0, -13, 14, 1};
    int realpdgs[] = {2, 4, -13, 14, 1, 4};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 4.232384112354395e-07
    double fks_g = 1.509770925047145e-24;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 0, alphas, ColorCorr, y, phi);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, 0, -13, 14, 1
// realpdgs = 2, 4, -13, 14, 1, 4
TEST(QCDCollinearISR2, Wj_born_2_real_27_n1) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.5127006420170197;
    ps.X2 = 0.6433648946156652;
    ps.Jacobian = 9764.296904798532;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3733.137335874217, 0, 0, 3733.137335874217);
    ps.Momenta[1].Set(3733.137335874217, 0, 0, -3733.137335874217);
    ps.Momenta[2].Set(1842.197455409808, -854.9356566062997, 717.844936868426,
                      -1465.426604915406);
    ps.Momenta[3].Set(3028.902473555646, -1276.666271083322, 23.63902819428961,
                      2746.600557586732);
    ps.Momenta[4].Set(2595.17474278298, 2131.601927689621, -741.4839650627157,
                      -1281.173952671327);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 3.408607795513846e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(2.208999073689689e-10, 0);
    spin[1][2] =
        std::complex<double>(8.078572264164807e-10, 2.270045791106549e-10);
    spin[2][1] =
        std::complex<double>(8.078572264164807e-10, -2.270045791106549e-10);
    spin[2][2] = std::complex<double>(3.187707888144877e-09, 0);

    // radiation variables
    double phi = 2.064162574707082;
    double xi = 0.1619070668495389;
    int bornpdgs[] = {2, 0, -13, 14, 1};
    int realpdgs[] = {2, 4, -13, 14, 1, 4};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 3.754722230599892e-07
    double fks_g = 1.968518132603939e-16;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, -1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 1.968481133036204e-16
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, 0, -13, 14, 1
// realpdgs = 2, 4, -13, 14, 1, 4
TEST(QCDSoftCollinearISR2, Wj_born_2_real_27_n1) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.9784759605801447;
    ps.X2 = 0.1226292857322733;
    ps.Jacobian = 4345.000447520756;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2251.570428484983, 0, 0, 2251.570428484983);
    ps.Momenta[1].Set(2251.570428484983, 0, 0, -2251.570428484983);
    ps.Momenta[2].Set(470.5632531734356, -203.9441739808373, 161.9957398870325,
                      391.910614037045);
    ps.Momenta[3].Set(2117.863722712123, 894.8546120059972, 1309.664875903483,
                      -1403.338834414683);
    ps.Momenta[4].Set(1914.713881084407, -690.9104380251599, -1471.660615790515,
                      1011.428220377638);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 7.173353271610537e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(3.501947226768883e-09, 0);
    spin[1][2] =
        std::complex<double>(1.337358204987732e-10, 3.583180842127741e-09);
    spin[2][1] =
        std::complex<double>(1.337358204987732e-10, -3.583180842127741e-09);
    spin[2][2] = std::complex<double>(3.671406044841654e-09, 0);

    // radiation variables
    double phi = 2.661246249935623;
    int bornpdgs[] = {2, 0, -13, 14, 1};
    int realpdgs[] = {2, 4, -13, 14, 1, 4};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 0.003187113443949893
    double fks_g = 4.90674866067634e-21;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, -1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, 0, -13, 14, 1
// realpdgs = 2, 5, -13, 14, 1, 5
TEST(QCDCollinearISR1, Wj_born_2_real_27_n2) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.07267831941426195;
    ps.X2 = 0.8796507711947541;
    ps.Jacobian = 2203.254895349373;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1643.504655682493, 0, 0, 1643.504655682493);
    ps.Momenta[1].Set(1643.504655682493, 0, 0, -1643.504655682493);
    ps.Momenta[2].Set(796.200342761695, 115.6205158961081, 78.85685425498588,
                      -783.8038521561158);
    ps.Momenta[3].Set(1160.681073510786, 785.5689162188655, 788.0047863225122,
                      330.3187687820358);
    ps.Momenta[4].Set(1330.127895092505, -901.1894321149736, -866.8616405774983,
                      453.4850833740801);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.006051791202788e-08;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(5.076749100274659e-09, 0);
    spin[1][2] =
        std::complex<double>(-5.01773604090335e-09, 3.516658284005982e-10);
    spin[2][1] =
        std::complex<double>(-5.01773604090335e-09, -3.516658284005982e-10);
    spin[2][2] = std::complex<double>(4.983768811753226e-09, 0);

    // radiation variables
    double phi = 5.312132943495435;
    double xi = 0.01543284828919841;
    int bornpdgs[] = {2, 0, -13, 14, 1};
    int realpdgs[] = {2, 5, -13, 14, 1, 5};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 2.322370731072227e-13
    double fks_g = 1.106251108075512e-24;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, 1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, 0, -13, 14, 1
// realpdgs = 2, 5, -13, 14, 1, 5
TEST(QCDSoftCollinearISR1, Wj_born_2_real_27_n2) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.8881596180531695;
    ps.X2 = 0.1700863516417286;
    ps.Jacobian = 4893.538764534212;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2526.35048635475, 0, 0, 2526.35048635475);
    ps.Momenta[1].Set(2526.35048635475, 0, 0, -2526.35048635475);
    ps.Momenta[2].Set(2133.708442944184, -326.6925184200375, 755.1235759218166,
                      1968.698073089603);
    ps.Momenta[3].Set(997.1003115833382, 853.9830363575979, -495.6195558352606,
                      -138.8641812949128);
    ps.Momenta[4].Set(1921.892218181978, -527.2905179375604, -259.5040200865561,
                      -1829.83389179469);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 5.108562803981937e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(1.627821611270265e-09, 0);
    spin[1][2] =
        std::complex<double>(4.456238733743926e-10, 2.338256850778695e-09);
    spin[2][1] =
        std::complex<double>(4.456238733743926e-10, -2.338256850778695e-09);
    spin[2][2] = std::complex<double>(3.480741192711672e-09, 0);

    // radiation variables
    double phi = 5.176328074366501;
    int bornpdgs[] = {2, 0, -13, 14, 1};
    int realpdgs[] = {2, 5, -13, 14, 1, 5};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 7.074785794995595e-09
    double fks_g = 1.769867550660703e-28;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, 1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, 0, -13, 14, 1
// realpdgs = 2, 5, -13, 14, 1, 5
TEST(QCDSoftLimit, Wj_born_2_real_27_0_n2) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.6954046152875115;
    ps.X2 = 0.9129522644152175;
    ps.Jacobian = 26026.1278704867;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(5179.122413057738, 0, 0, 5179.122413057738);
    ps.Momenta[1].Set(5179.122413057738, 0, 0, -5179.122413057738);
    ps.Momenta[2].Set(4463.179833529996, 3161.770996388174, -2996.489532809864,
                      -971.7143986729329);
    ps.Momenta[3].Set(909.0568943201122, -218.8799505025352, -789.9876190708583,
                      -392.9320120500388);
    ps.Momenta[4].Set(4986.008098265362, -2942.891045885643, 3786.477151880727,
                      1364.646410722973);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 2.708721051088922e-09;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, 4.063081576633383e-09);
    ColorCorr.Set(0, 4, -4.514535085148203e-10);
    ColorCorr.Set(1, 0, 4.063081576633383e-09);
    ColorCorr.Set(1, 4, 4.063081576633383e-09);
    ColorCorr.Set(4, 0, -4.514535085148203e-10);
    ColorCorr.Set(4, 1, 4.063081576633383e-09);

    // radiation variables
    double y = -0.1057052689504943;
    double phi = 2.804323396422849;
    int bornpdgs[] = {2, 0, -13, 14, 1};
    int realpdgs[] = {2, 5, -13, 14, 1, 5};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 1.469678852563599e-07
    double fks_g = 3.353522527345009e-25;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 0, alphas, ColorCorr, y, phi);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, 0, -13, 14, 1
// realpdgs = 2, 5, -13, 14, 1, 5
TEST(QCDCollinearISR2, Wj_born_2_real_27_n2) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.3459991942756382;
    ps.X2 = 0.01963334585335019;
    ps.Jacobian = 213.6353923652897;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(535.7325806796719, 0, 0, 535.7325806796719);
    ps.Momenta[1].Set(535.7325806796719, 0, 0, -535.7325806796719);
    ps.Momenta[2].Set(190.9327334201298, 182.4398476031996, 48.68101620475002,
                      -28.30493524240101);
    ps.Momenta[3].Set(484.870065139994, -194.2492015706228, -443.9734683772604,
                      15.93069788500746);
    ps.Momenta[4].Set(395.6623627992199, 11.80935396742321, 395.2924521725104,
                      12.37423735739355);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 9.374871232011947e-08;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(5.535459059278309e-08, 0);
    spin[1][2] =
        std::complex<double>(-2.570602110039611e-08, 3.82686734624049e-08);
    spin[2][1] =
        std::complex<double>(-2.570602110039611e-08, -3.82686734624049e-08);
    spin[2][2] = std::complex<double>(3.839412172733638e-08, 0);

    // radiation variables
    double phi = 4.767779376789532;
    double xi = 0.4495253089061138;
    int bornpdgs[] = {2, 0, -13, 14, 1};
    int realpdgs[] = {2, 5, -13, 14, 1, 5};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 0.0002857629723107367
    double fks_g = 1.154899640596818e-12;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, -1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 1.154866383902436e-12
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, 0, -13, 14, 1
// realpdgs = 2, 5, -13, 14, 1, 5
TEST(QCDSoftCollinearISR2, Wj_born_2_real_27_n2) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.3149823629405368;
    ps.X2 = 0.09525777352192932;
    ps.Jacobian = 379.7261556037765;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1125.917808069734, 0, 0, 1125.917808069734);
    ps.Momenta[1].Set(1125.917808069734, 0, 0, -1125.917808069734);
    ps.Momenta[2].Set(811.8611391469356, -338.2231108896053, 422.1760138278459,
                      605.3850426509956);
    ps.Momenta[3].Set(1105.345596368789, 555.3802821809596, -416.9715864526596,
                      -859.9280933122399);
    ps.Momenta[4].Set(334.6288806237425, -217.1571712913542, -5.2044273751864,
                      254.5430506612443);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 2.087389411366352e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(1.052635747977662e-09, 0);
    spin[1][2] =
        std::complex<double>(6.191941314036664e-10, 8.401293496057851e-10);
    spin[2][1] =
        std::complex<double>(6.191941314036664e-10, -8.401293496057851e-10);
    spin[2][2] = std::complex<double>(1.034753663388691e-09, 0);

    // radiation variables
    double phi = 3.883770301804703;
    int bornpdgs[] = {2, 0, -13, 14, 1};
    int realpdgs[] = {2, 5, -13, 14, 1, 5};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 0.00359836315233666
    double fks_g = 5.890941743567575e-21;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, -1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, 0, -13, 14, 1
// realpdgs = 2, -5, -13, 14, 1, -5
TEST(QCDCollinearISR1, Wj_born_2_real_38) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.0934123978832595;
    ps.X2 = 0.282962511611899;
    ps.Jacobian = 667.4106567042554;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1056.769006903185, 0, 0, 1056.769006903185);
    ps.Momenta[1].Set(1056.769006903185, 0, 0, -1056.769006903185);
    ps.Momenta[2].Set(647.8290695935733, -128.0466212474549, 620.6412272269797,
                      134.5029117297242);
    ps.Momenta[3].Set(839.0768050498112, 266.1027773042136, -637.3246452237479,
                      476.5044525227513);
    ps.Momenta[4].Set(626.6321391629864, -138.0561560567588, 16.68341799676827,
                      -611.0073642524756);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 3.518101314060807e-07;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(3.621950598255754e-08, 0);
    spin[1][2] =
        std::complex<double>(7.470429040718595e-08, 7.648402147066537e-08);
    spin[2][1] =
        std::complex<double>(7.470429040718595e-08, -7.648402147066537e-08);
    spin[2][2] = std::complex<double>(3.155906254235232e-07, 0);

    // radiation variables
    double phi = 5.442676280466418;
    double xi = 0.2226345499626115;
    int bornpdgs[] = {2, 0, -13, 14, 1};
    int realpdgs[] = {2, -5, -13, 14, 1, -5};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 8.105236583946086e-14
    double fks_g = 8.034906280730616e-23;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, 1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, 0, -13, 14, 1
// realpdgs = 2, -5, -13, 14, 1, -5
TEST(QCDSoftCollinearISR1, Wj_born_2_real_38) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.5083319970640297;
    ps.X2 = 0.906233730316778;
    ps.Jacobian = 14083.49796105796;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(4411.712386580819, 0, 0, 4411.712386580819);
    ps.Momenta[1].Set(4411.712386580819, 0, 0, -4411.712386580819);
    ps.Momenta[2].Set(3276.071923222501, -1775.142352757268, 2323.106279663823,
                      1478.071069658249);
    ps.Momenta[3].Set(2379.952648709348, -433.8916637729579, -219.0244168687354,
                      -2329.794183833044);
    ps.Momenta[4].Set(3167.400201229791, 2209.034016530225, -2104.081862795088,
                      851.7231141747952);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 2.425485921670822e-10;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(1.815011271552214e-10, 0);
    spin[1][2] =
        std::complex<double>(1.015159789462784e-10, 2.783324861025818e-11);
    spin[2][1] =
        std::complex<double>(1.015159789462784e-10, -2.783324861025818e-11);
    spin[2][2] = std::complex<double>(6.104746501186084e-11, 0);

    // radiation variables
    double phi = 5.167794710696446;
    int bornpdgs[] = {2, 0, -13, 14, 1};
    int realpdgs[] = {2, -5, -13, 14, 1, -5};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 2.505605700508366e-11
    double fks_g = 1.211397648713153e-29;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, 1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, 0, -13, 14, 1
// realpdgs = 2, -5, -13, 14, 1, -5
TEST(QCDSoftLimit, Wj_born_2_real_38_0) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.6490760994986182;
    ps.X2 = 0.6787902850304022;
    ps.Jacobian = 3102.443915928419;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(4314.485109746005, 0, 0, 4314.485109746005);
    ps.Momenta[1].Set(4314.485109746005, 0, 0, -4314.485109746005);
    ps.Momenta[2].Set(3897.037426157165, 1778.828917147699, -65.97575226395895,
                      -3466.744234076405);
    ps.Momenta[3].Set(4018.4646913815, -2421.889338918865, -98.53378946468484,
                      3205.121151882512);
    ps.Momenta[4].Set(713.4681019533471, 643.0604217711664, 164.5095417286438,
                      261.6230821938935);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 2.363132770292018e-09;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, 3.544699155438027e-09);
    ColorCorr.Set(0, 4, -3.938554617153363e-10);
    ColorCorr.Set(1, 0, 3.544699155438027e-09);
    ColorCorr.Set(1, 4, 3.544699155438027e-09);
    ColorCorr.Set(4, 0, -3.938554617153363e-10);
    ColorCorr.Set(4, 1, 3.544699155438027e-09);

    // radiation variables
    double y = -0.4958220496580381;
    double phi = 2.961159446648732;
    int bornpdgs[] = {2, 0, -13, 14, 1};
    int realpdgs[] = {2, -5, -13, 14, 1, -5};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 1.226075905745939e-07
    double fks_g = 1.524404353700275e-24;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 0, alphas, ColorCorr, y, phi);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, 0, -13, 14, 1
// realpdgs = 2, -5, -13, 14, 1, -5
TEST(QCDCollinearISR2, Wj_born_2_real_38) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.226752738596371;
    ps.X2 = 0.6241126477553465;
    ps.Jacobian = 2799.001243182022;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2445.23790253758, 0, 0, 2445.23790253758);
    ps.Momenta[1].Set(2445.23790253758, 0, 0, -2445.23790253758);
    ps.Momenta[2].Set(1399.85732382194, -912.8763168606822, -1058.596574458777,
                      75.03765532678244);
    ps.Momenta[3].Set(2354.871558507868, 995.3404269541586, 2030.668244641167,
                      -656.5850834019532);
    ps.Momenta[4].Set(1135.746922745353, -82.46411009347636, -972.071670182391,
                      581.5474280751708);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.282430246989327e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(9.592649380387613e-10, 0);
    spin[1][2] =
        std::complex<double>(-4.847955798336384e-10, 2.738145281764792e-10);
    spin[2][1] =
        std::complex<double>(-4.847955798336384e-10, -2.738145281764792e-10);
    spin[2][2] = std::complex<double>(3.231653089505658e-10, 0);

    // radiation variables
    double phi = 0.4249526895543917;
    double xi = 0.1759873597875576;
    int bornpdgs[] = {2, 0, -13, 14, 1};
    int realpdgs[] = {2, -5, -13, 14, 1, -5};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 2.76715831277946e-07
    double fks_g = 1.714063684492073e-16;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, -1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 1.714064682298519e-16
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, 0, -13, 14, 1
// realpdgs = 2, -5, -13, 14, 1, -5
TEST(QCDSoftCollinearISR2, Wj_born_2_real_38) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.6984308920005908;
    ps.X2 = 0.4595016995089196;
    ps.Jacobian = 2519.795202805865;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3682.295504674447, 0, 0, 3682.295504674447);
    ps.Momenta[1].Set(3682.295504674447, 0, 0, -3682.295504674447);
    ps.Momenta[2].Set(3285.571035993529, -997.0637076715622, 121.07696918751,
                      3128.287289066685);
    ps.Momenta[3].Set(3400.056824315664, 517.063146164127, 306.6836825288706,
                      -3346.487297200588);
    ps.Momenta[4].Set(678.9631490397036, 480.000561507435, -427.7606517163807,
                      218.2000081339025);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 5.036960348690437e-11;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(1.562648667043509e-11, 5.385290446308774e-28);
    spin[1][2] =
        std::complex<double>(9.6750188917107e-13, 2.328039501203484e-11);
    spin[2][1] =
        std::complex<double>(9.6750188917107e-13, -2.328039501203484e-11);
    spin[2][2] =
        std::complex<double>(3.474311681646928e-11, -5.385290446308774e-28);

    // radiation variables
    double phi = 6.132872047855688;
    int bornpdgs[] = {2, 0, -13, 14, 1};
    int realpdgs[] = {2, -5, -13, 14, 1, -5};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 1.358914349624687e-05
    double fks_g = 7.939823450023886e-24;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, -1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, 0, -13, 14, 1
// realpdgs = 2, -4, -13, 14, 1, -4
TEST(QCDCollinearISR1, Wj_born_2_real_38_n1) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.6415063474526299;
    ps.X2 = 0.1636808977239284;
    ps.Jacobian = 3768.534862688064;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2106.264144705203, 0, 0, 2106.264144705203);
    ps.Momenta[1].Set(2106.264144705203, 0, 0, -2106.264144705203);
    ps.Momenta[2].Set(353.419251031058, 219.5138762125191, 61.85709747627357,
                      269.9861563882852);
    ps.Momenta[3].Set(2083.85994741459, -1193.775800417749, -1230.742292994344,
                      -1184.333157101804);
    ps.Momenta[4].Set(1775.249090964759, 974.2619242052297, 1168.885195518071,
                      914.3470007135187);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 8.794802780526728e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(4.244154038321255e-09, 0);
    spin[1][2] =
        std::complex<double>(-1.423957115454247e-10, 4.392422759428601e-09);
    spin[2][1] =
        std::complex<double>(-1.423957115454247e-10, -4.392422759428601e-09);
    spin[2][2] = std::complex<double>(4.550648742205473e-09, 0);

    // radiation variables
    double phi = 0.8957419042360475;
    double xi = 0.3082886521519465;
    int bornpdgs[] = {2, 0, -13, 14, 1};
    int realpdgs[] = {2, -4, -13, 14, 1, -4};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 6.339792297254159e-15
    double fks_g = 1.205091722260399e-23;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, 1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, 0, -13, 14, 1
// realpdgs = 2, -4, -13, 14, 1, -4
TEST(QCDSoftCollinearISR1, Wj_born_2_real_38_n1) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.7177036572377382;
    ps.X2 = 0.129070663933561;
    ps.Jacobian = 1624.615910118979;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1978.334425437225, 0, 0, 1978.334425437225);
    ps.Momenta[1].Set(1978.334425437225, 0, 0, -1978.334425437225);
    ps.Momenta[2].Set(1696.760667634687, -21.0174318518238, 1004.125313862124,
                      1367.584507389489);
    ps.Momenta[3].Set(1445.108990372528, -631.6543034229906, -954.3524040227217,
                      -882.3629207751056);
    ps.Momenta[4].Set(814.7991928672357, 652.6717352748145, -49.77290983940223,
                      -485.2215866143837);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 7.982560947903866e-10;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(5.35413410792959e-10, 0);
    spin[1][2] =
        std::complex<double>(-3.19961560515429e-10, 1.958420223893994e-10);
    spin[2][1] =
        std::complex<double>(-3.19961560515429e-10, -1.958420223893994e-10);
    spin[2][2] = std::complex<double>(2.628426839974276e-10, 0);

    // radiation variables
    double phi = 4.799568232670476;
    int bornpdgs[] = {2, 0, -13, 14, 1};
    int realpdgs[] = {2, -4, -13, 14, 1, -4};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 7.142297502087223e-10
    double fks_g = 1.138357285014694e-28;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, 1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, 0, -13, 14, 1
// realpdgs = 2, -4, -13, 14, 1, -4
TEST(QCDSoftLimit, Wj_born_2_real_38_0_n1) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.9272895720433034;
    ps.X2 = 0.3079451688788382;
    ps.Jacobian = 191.5994401112016;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3473.423531357638, 0, 0, 3473.423531357638);
    ps.Momenta[1].Set(3473.423531357638, 0, 0, -3473.423531357638);
    ps.Momenta[2].Set(3429.120386043398, 1368.00516381102, 3001.616846993668,
                      -936.8696801629314);
    ps.Momenta[3].Set(3462.995335577462, -1357.593762402024, -3030.63967256302,
                      982.0891230461116);
    ps.Momenta[4].Set(54.73134109441717, -10.41140140899546, 29.02282556935184,
                      -45.21944288318004);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.367089790591581e-07;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, 2.050634685887372e-07);
    ColorCorr.Set(0, 4, -2.278482984319302e-08);
    ColorCorr.Set(1, 0, 2.050634685887372e-07);
    ColorCorr.Set(1, 4, 2.050634685887372e-07);
    ColorCorr.Set(4, 0, -2.278482984319302e-08);
    ColorCorr.Set(4, 1, 2.050634685887372e-07);

    // radiation variables
    double y = 0.3818423360165966;
    double phi = 3.459869997331597;
    int bornpdgs[] = {2, 0, -13, 14, 1};
    int realpdgs[] = {2, -4, -13, 14, 1, -4};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 1.566061038546611e-05
    double fks_g = 1.433691043186376e-23;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 0, alphas, ColorCorr, y, phi);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, 0, -13, 14, 1
// realpdgs = 2, -4, -13, 14, 1, -4
TEST(QCDCollinearISR2, Wj_born_2_real_38_n1) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.3367513561198798;
    ps.X2 = 0.4822198464215006;
    ps.Jacobian = 4227.489546635713;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2619.332149705862, 0, 0, 2619.332149705862);
    ps.Momenta[1].Set(2619.332149705862, 0, 0, -2619.332149705862);
    ps.Momenta[2].Set(1932.413341540401, -759.2530397138045, 670.5933986354246,
                      -1645.61861862197);
    ps.Momenta[3].Set(1704.881411078636, 1519.701142015384, 359.6128069544015,
                      683.9646875921251);
    ps.Momenta[4].Set(1601.369546792686, -760.4481023015799, -1030.206205589826,
                      961.6539310298444);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 2.310216028863483e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(6.153441237511031e-10, 0);
    spin[1][2] =
        std::complex<double>(-5.534347191504407e-10, -8.582770409143036e-10);
    spin[2][1] =
        std::complex<double>(-5.534347191504407e-10, 8.582770409143036e-10);
    spin[2][2] = std::complex<double>(1.69487190511238e-09, 0);

    // radiation variables
    double phi = 1.673774329891087;
    double xi = 0.3512544806511668;
    int bornpdgs[] = {2, 0, -13, 14, 1};
    int realpdgs[] = {2, -4, -13, 14, 1, -4};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 4.422010257062609e-07
    double fks_g = 1.091172687262474e-15;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, -1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 1.091126855361873e-15
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, 0, -13, 14, 1
// realpdgs = 2, -4, -13, 14, 1, -4
TEST(QCDSoftCollinearISR2, Wj_born_2_real_38_n1) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.5818030077562604;
    ps.X2 = 0.1957110048928641;
    ps.Jacobian = 615.3642574790694;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2193.355162149189, 0, 0, 2193.355162149189);
    ps.Momenta[1].Set(2193.355162149189, 0, 0, -2193.355162149189);
    ps.Momenta[2].Set(2094.103099065282, -1704.183509828873, 1031.395485287236,
                      645.9486878005616);
    ps.Momenta[3].Set(2014.236921970512, 1715.879171378816, -987.4362665447089,
                      -371.3201672005355);
    ps.Momenta[4].Set(278.3703032625844, -11.69566154994361, -43.95921874252757,
                      -274.6285206000261);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.744069562787878e-07;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(1.044771923166144e-07, 0);
    spin[1][2] =
        std::complex<double>(-1.588267632017499e-08, -8.398693940750545e-08);
    spin[2][1] =
        std::complex<double>(-1.588267632017499e-08, 8.398693940750545e-08);
    spin[2][2] = std::complex<double>(6.992976396217334e-08, 0);

    // radiation variables
    double phi = 0.5375541506564148;
    int bornpdgs[] = {2, 0, -13, 14, 1};
    int realpdgs[] = {2, -4, -13, 14, 1, -4};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 0.08907970430792482
    double fks_g = 1.152479098549563e-19;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, -1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, 0, -13, 14, 1
// realpdgs = 2, -3, -13, 14, 1, -3
TEST(QCDCollinearISR1, Wj_born_2_real_38_n2) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.1450605395396316;
    ps.X2 = 0.639990346793752;
    ps.Jacobian = 3115.265054311012;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1980.499388160959, 0, 0, 1980.499388160959);
    ps.Momenta[1].Set(1980.499388160959, 0, 0, -1980.499388160959);
    ps.Momenta[2].Set(460.5449970500319, -137.0581499096358, -146.298536785089,
                      -414.6245241006656);
    ps.Momenta[3].Set(1939.752137512653, 691.3193173777195, -398.6769990342045,
                      1767.98546567736);
    ps.Momenta[4].Set(1560.701641759234, -554.2611674680838, 544.9755358192937,
                      -1353.360941576695);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 8.27363303021152e-08;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(3.831859087148119e-08, 0);
    spin[1][2] =
        std::complex<double>(2.178371798176433e-08, 3.503562209479774e-08);
    spin[2][1] =
        std::complex<double>(2.178371798176433e-08, -3.503562209479774e-08);
    spin[2][2] = std::complex<double>(4.441773943063401e-08, 0);

    // radiation variables
    double phi = 5.450306773385022;
    double xi = 0.7776783673974849;
    int bornpdgs[] = {2, 0, -13, 14, 1};
    int realpdgs[] = {2, -3, -13, 14, 1, -3};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 1.730675946925351e-15
    double fks_g = 2.093369007374619e-23;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, 1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, 0, -13, 14, 1
// realpdgs = 2, -3, -13, 14, 1, -3
TEST(QCDSoftCollinearISR1, Wj_born_2_real_38_n2) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.07101888248355337;
    ps.X2 = 0.1155125424478527;
    ps.Jacobian = 338.4998578316498;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(588.728208406521, 0, 0, 588.728208406521);
    ps.Momenta[1].Set(588.728208406521, 0, 0, -588.728208406521);
    ps.Momenta[2].Set(444.73066569509, 114.1760402737009, 86.82620128842009,
                      -420.9636654235599);
    ps.Momenta[3].Set(162.2420488837698, -78.94334734680311, 0.2040069996975831,
                      -141.7405683522879);
    ps.Momenta[4].Set(570.4837022341819, -35.23269292689784, -87.03020828811775,
                      562.7042337758483);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 4.526049051971873e-07;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(1.047078810798922e-07, 0);
    spin[1][2] =
        std::complex<double>(-5.994280677743967e-08, -1.812027045740892e-07);
    spin[2][1] =
        std::complex<double>(-5.994280677743967e-08, 1.812027045740892e-07);
    spin[2][2] = std::complex<double>(3.478970241172951e-07, 0);

    // radiation variables
    double phi = 1.119641633544825;
    int bornpdgs[] = {2, 0, -13, 14, 1};
    int realpdgs[] = {2, -3, -13, 14, 1, -3};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 1.389577083784841e-06
    double fks_g = 2.398426575300156e-24;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, 1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, 0, -13, 14, 1
// realpdgs = 2, -3, -13, 14, 1, -3
TEST(QCDSoftLimit, Wj_born_2_real_38_0_n2) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.7003098619995178;
    ps.X2 = 0.4453717181547319;
    ps.Jacobian = 5036.679924947356;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3630.110084247563, 0, 0, 3630.110084247563);
    ps.Momenta[1].Set(3630.110084247563, 0, 0, -3630.110084247563);
    ps.Momenta[2].Set(2930.137221442057, 2152.65834419302, 1987.771997654985,
                      22.99293358954048);
    ps.Momenta[3].Set(2953.430998291912, -2125.839724592252, -1589.822402191716,
                      1294.613786619956);
    ps.Momenta[4].Set(1376.651948761159, -26.81861960076849, -397.9495954632699,
                      -1317.606720209496);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.670249611419021e-08;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, 2.505374417128532e-08);
    ColorCorr.Set(0, 4, -2.783749352365036e-09);
    ColorCorr.Set(1, 0, 2.505374417128532e-08);
    ColorCorr.Set(1, 4, 2.505374417128532e-08);
    ColorCorr.Set(4, 0, -2.783749352365036e-09);
    ColorCorr.Set(4, 1, 2.505374417128532e-08);

    // radiation variables
    double y = 0.98401511229455;
    double phi = 1.203364360855731;
    int bornpdgs[] = {2, 0, -13, 14, 1};
    int realpdgs[] = {2, -3, -13, 14, 1, -3};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 4.185954576695086e-07
    double fks_g = 1.208688776938016e-25;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 0, alphas, ColorCorr, y, phi);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, 0, -13, 14, 1
// realpdgs = 2, -3, -13, 14, 1, -3
TEST(QCDCollinearISR2, Wj_born_2_real_38_n2) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.3905811283392175;
    ps.X2 = 0.8668161193378054;
    ps.Jacobian = 2745.369036061362;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3782.095352914682, 0, 0, 3782.095352914682);
    ps.Momenta[1].Set(3782.095352914682, 0, 0, -3782.095352914682);
    ps.Momenta[2].Set(3152.764943232693, 467.3136246686109, 2884.346690077547,
                      -1184.098363699075);
    ps.Momenta[3].Set(3691.201324206067, -1014.138151721009, -3352.666133554659,
                      1164.52600742112);
    ps.Momenta[4].Set(720.2244383906045, 546.8245270523981, 468.3194434771116,
                      19.57235627795556);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.770137293059626e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(1.075960708590328e-09, 0);
    spin[1][2] =
        std::complex<double>(-4.336715086632036e-10, -7.475531768952496e-10);
    spin[2][1] =
        std::complex<double>(-4.336715086632036e-10, 7.475531768952496e-10);
    spin[2][2] = std::complex<double>(6.941765844692979e-10, 0);

    // radiation variables
    double phi = 5.41936878608916;
    double xi = 0.09694249257990635;
    int bornpdgs[] = {2, 0, -13, 14, 1};
    int realpdgs[] = {2, -3, -13, 14, 1, -3};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 3.065011560263842e-07
    double fks_g = 5.760901855127569e-17;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, -1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 5.760891568094953e-17
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, 0, -13, 14, 1
// realpdgs = 2, -3, -13, 14, 1, -3
TEST(QCDSoftCollinearISR2, Wj_born_2_real_38_n2) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.2654624806503456;
    ps.X2 = 0.1669362634361562;
    ps.Jacobian = 366.6918555946463;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1368.328192337479, 0, 0, 1368.328192337479);
    ps.Momenta[1].Set(1368.328192337479, 0, 0, -1368.328192337479);
    ps.Momenta[2].Set(1276.436910059769, -858.1051294069015, -944.9541682951357,
                      -2.896906996125939);
    ps.Momenta[3].Set(1194.324220663109, 964.7401978313461, 703.4581740648601,
                      -28.86679913344696);
    ps.Momenta[4].Set(265.8952539520801, -106.6350684244446, 241.4959942302756,
                      31.7637061295729);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 7.661183046104932e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(2.762986116780574e-09, -3.446585885637615e-26);
    spin[1][2] =
        std::complex<double>(1.515878236701571e-10, -3.67568649978414e-09);
    spin[2][1] =
        std::complex<double>(1.515878236701571e-10, 3.67568649978414e-09);
    spin[2][2] =
        std::complex<double>(4.898196929324358e-09, 3.446585885637615e-26);

    // radiation variables
    double phi = 5.871895318482464;
    int bornpdgs[] = {2, 0, -13, 14, 1};
    int realpdgs[] = {2, -3, -13, 14, 1, -3};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 0.009708595318932128
    double fks_g = 1.347543801532584e-20;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, -1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, 0, -13, 14, 1
// realpdgs = 2, 0, -13, 14, 1, 0
TEST(QCDCollinearFSR, Wj_born_2_real_31) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.3770457874348772;
    ps.X2 = 0.4800333686035074;
    ps.Jacobian = 5462.163984588711;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2765.324598883525, 0, 0, 2765.324598883525);
    ps.Momenta[1].Set(2765.324598883525, 0, 0, -2765.324598883525);
    ps.Momenta[2].Set(1868.370860045532, 975.6181079518453, -1589.328174941782,
                      114.0829980319787);
    ps.Momenta[3].Set(1702.449183803549, -986.6739669806036, 293.3234589942759,
                      1356.01218826372);
    ps.Momenta[4].Set(1959.82915391797, 11.05585902875824, 1296.004715947506,
                      -1470.095186295699);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 5.87667211182889e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(5.290999813295833e-09, -1.378634354255046e-25);
    spin[1][2] =
        std::complex<double>(-1.759115445778105e-09, 6.56115128300073e-11);
    spin[2][1] =
        std::complex<double>(-1.759115445778105e-09, -6.56115128300073e-11);
    spin[2][2] =
        std::complex<double>(5.856722985330567e-10, 1.378634354255046e-25);

    // radiation variables
    double phi = 0.2313956609670529;
    double xi = 0.1636588267787018;
    int bornpdgs[] = {2, 0, -13, 14, 1};
    int realpdgs[] = {2, 0, -13, 14, 1, 0};

    // FKS g = real matrix element * (1.0 - y) * xi * xi
    // real = 1.17433065044995e-05
    double fks_g = 3.145352076760794e-15;
    double limit = FKS::QCD::CollinearLimitFSR(
        realpdgs[4], bornpdgs[4], ps, 4, xi, phi, alphas, born,
        spin); // limit = 3.1441340363137e-15
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, 0, -13, 14, 1
// realpdgs = 2, 0, -13, 14, 1, 0
TEST(QCDSoftCollinearFSR, Wj_born_2_real_31) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.9004719521828584;
    ps.X2 = 0.2441919490225786;
    ps.Jacobian = 849.4315419555834;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3047.99410171612, 0, 0, 3047.99410171612);
    ps.Momenta[1].Set(3047.99410171612, 0, 0, -3047.99410171612);
    ps.Momenta[2].Set(3039.525943171606, 734.5629498003329, 167.8508536884311,
                      -2944.649609530808);
    ps.Momenta[3].Set(2779.950340282629, -612.7685970573668, -222.2885368456735,
                      2702.448213617582);
    ps.Momenta[4].Set(276.5119199780028, -121.7943527429661, 54.43768315724241,
                      242.2013959132256);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 7.542021729236961e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(3.774633570642777e-09, -2.154116178523509e-27);
    spin[1][2] =
        std::complex<double>(2.103546729485068e-11, -3.770950453955523e-09);
    spin[2][1] =
        std::complex<double>(2.103546729485068e-11, 3.770950453955523e-09);
    spin[2][2] =
        std::complex<double>(3.767388158594184e-09, 2.154116178523509e-27);

    // radiation variables
    double phi = 3.604174131961644;
    int bornpdgs[] = {2, 0, -13, 14, 1};
    int realpdgs[] = {2, 0, -13, 14, 1, 0};

    // FKS g = real matrix element * (1.0 - y) * xi * xi
    // real = 389664.3173433932
    double fks_g = 3.206934642719195e-15;
    double limit = FKS::QCD::SoftCollinearLimitFSR(realpdgs[4], bornpdgs[4], ps,
                                                   4, phi, alphas, born, spin);
    // limit = 3.210111654961085e-15
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, 0, -13, 14, 1
// realpdgs = 2, 0, -13, 14, 1, 0
TEST(QCDSoftLimit, Wj_born_2_real_31_4) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.01202114845693458;
    ps.X2 = 0.3135983172075996;
    ps.Jacobian = 12.43908029074016;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(399.0921621826729, 0, 0, 399.0921621826729);
    ps.Momenta[1].Set(399.0921621826729, 0, 0, -399.0921621826729);
    ps.Momenta[2].Set(393.4892325588161, -186.4244541352915, -224.1388571527413,
                      264.2753710703055);
    ps.Momenta[3].Set(373.7697385821492, 194.6167506715562, 214.8524046114242,
                      -235.9376656460395);
    ps.Momenta[4].Set(30.92535322438049, -8.192296536264619, 9.28645254131713,
                      -28.33770542426595);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    // double born = 3.186866472774556e-07;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, 4.780299709161834e-07);
    ColorCorr.Set(0, 4, -5.311444121290927e-08);
    ColorCorr.Set(1, 0, 4.780299709161834e-07);
    ColorCorr.Set(1, 4, 4.780299709161834e-07);
    ColorCorr.Set(4, 0, -5.311444121290927e-08);
    ColorCorr.Set(4, 1, 4.780299709161834e-07);

    // radiation variables
    double y = -0.6656446808069134;
    double phi = 4.999983716412815;
    int bornpdgs[] = {2, 0, -13, 14, 1};
    int realpdgs[] = {2, 0, -13, 14, 1, 0};

    // FKS g = real matrix element * (1.0 - y) * xi * xi
    // real = 37435756.3182152
    double fks_g = 3.744138574125524e-11;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 4, alphas, ColorCorr, y, phi);
    // limit = 3.744139309923325e-11
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, 0, -13, 14, 1
// realpdgs = 2, 0, -13, 14, 1, 0
TEST(QCDCollinearISR1, Wj_born_2_real_31) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.6060544111572632;
    ps.X2 = 0.8044048071168959;
    ps.Jacobian = 3623.594779703089;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(4538.438905859372, 0, 0, 4538.438905859372);
    ps.Momenta[1].Set(4538.438905859372, 0, 0, -4538.438905859372);
    ps.Momenta[2].Set(4058.997080835908, 2220.105105539229, 470.8897287152479,
                      3365.241965443726);
    ps.Momenta[3].Set(4225.684579581853, -2050.797967765768, -1235.9516219834,
                      -3481.818698559741);
    ps.Momenta[4].Set(792.1961513009851, -169.3071377734603, 765.0618932681521,
                      116.5767331160149);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 5.61768439051389e-11;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(1.223652240905937e-11, -2.154116178523509e-27);
    spin[1][2] =
        std::complex<double>(-2.145879422057551e-11, 8.786175461885876e-12);
    spin[2][1] =
        std::complex<double>(-2.145879422057551e-11, -8.786175461885876e-12);
    spin[2][2] =
        std::complex<double>(4.394032149607953e-11, 2.154116178523509e-27);

    // radiation variables
    double phi = 0.5198639002936677;
    double xi = 0.2808949888078284;
    int bornpdgs[] = {2, 0, -13, 14, 1};
    int realpdgs[] = {2, 0, -13, 14, 1, 0};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 1.036788727620007e-08
    double fks_g = 1.636093973747185e-17;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, 1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 1.636148057834657e-17
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, 0, -13, 14, 1
// realpdgs = 2, 0, -13, 14, 1, 0
TEST(QCDSoftCollinearISR1, Wj_born_2_real_31) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.4019528323223192;
    ps.X2 = 0.173289429878146;
    ps.Jacobian = 2345.647528646842;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1715.48505811956, 0, 0, 1715.48505811956);
    ps.Momenta[1].Set(1715.48505811956, 0, 0, -1715.48505811956);
    ps.Momenta[2].Set(1147.741051823002, -819.8057540647067, -630.7959161387662,
                      497.3173632850788);
    ps.Momenta[3].Set(926.5553530653015, 663.6709937660273, -92.26234892838846,
                      639.9478832669072);
    ps.Momenta[4].Set(1356.673711350817, 156.1347602986792, 723.0582650671546,
                      -1137.265246551986);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 2.089148670256226e-08;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(1.850564889565343e-08, 0);
    spin[1][2] =
        std::complex<double>(4.028195510192175e-09, 5.28442217354917e-09);
    spin[2][1] =
        std::complex<double>(4.028195510192175e-09, -5.28442217354917e-09);
    spin[2][2] = std::complex<double>(2.385837806908828e-09, 0);

    // radiation variables
    double phi = 6.028970976122172;
    int bornpdgs[] = {2, 0, -13, 14, 1};
    int realpdgs[] = {2, 0, -13, 14, 1, 0};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 78486.74314012504
    double fks_g = 5.614321349615434e-14;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, 1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 5.61417656066481e-14
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, 0, -13, 14, 1
// realpdgs = 2, 0, -13, 14, 1, 0
TEST(QCDSoftLimit, Wj_born_2_real_31_0) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.1629217384040551;
    ps.X2 = 0.8488584038264122;
    ps.Jacobian = 2232.114141677898;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2417.244054235065, 0, 0, 2417.244054235065);
    ps.Momenta[1].Set(2417.244054235065, 0, 0, -2417.244054235065);
    ps.Momenta[2].Set(1592.884813106145, -1463.226748746994, 441.1818777377197,
                      -449.0078622154808);
    ps.Momenta[3].Set(2325.392309880836, 2176.642629133542, -818.1994888660812,
                      15.02851658788455);
    ps.Momenta[4].Set(916.2109854831507, -713.4158803865479, 377.0176111283616,
                      433.9793456275962);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    // double born = 1.636838418111114e-09;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, 2.455257627166672e-09);
    ColorCorr.Set(0, 4, -2.728064030185191e-10);
    ColorCorr.Set(1, 0, 2.455257627166672e-09);
    ColorCorr.Set(1, 4, 2.455257627166672e-09);
    ColorCorr.Set(4, 0, -2.728064030185191e-10);
    ColorCorr.Set(4, 1, 2.455257627166672e-09);

    // radiation variables
    double y = -0.8536953201371134;
    double phi = 4.241563089381033;
    int bornpdgs[] = {2, 0, -13, 14, 1};
    int realpdgs[] = {2, 0, -13, 14, 1, 0};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 6871.57671279481
    double fks_g = 4.896722652233856e-15;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 0, alphas, ColorCorr, y, phi);
    // limit = 4.896722482819271e-15
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, 0, -13, 14, 1
// realpdgs = 2, 0, -13, 14, 1, 0
TEST(QCDCollinearISR2, Wj_born_2_real_31) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.8412021518792265;
    ps.X2 = 0.408338042899139;
    ps.Jacobian = 9396.009290260135;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3809.548136733981, 0, 0, 3809.548136733981);
    ps.Momenta[1].Set(3809.548136733981, 0, 0, -3809.548136733981);
    ps.Momenta[2].Set(3453.322714473803, 2489.930318953436, 2247.996402634923,
                      819.876180162346);
    ps.Momenta[3].Set(1718.572971100881, -246.4219166942273, -1700.78534484746,
                      -9.925055332783245);
    ps.Momenta[4].Set(2447.20058789328, -2243.508402259208, -547.2110577874631,
                      -809.9511248295626);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 6.468808524832781e-10;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(3.557984158535375e-10, 0);
    spin[1][2] =
        std::complex<double>(1.598286601933609e-10, -2.793232342939709e-10);
    spin[2][1] =
        std::complex<double>(1.598286601933609e-10, 2.793232342939709e-10);
    spin[2][2] = std::complex<double>(2.910824366297405e-10, 0);

    // radiation variables
    double phi = 3.743058894295403;
    double xi = 0.5440971379468943;
    int bornpdgs[] = {2, 0, -13, 14, 1};
    int realpdgs[] = {2, 0, -13, 14, 1, 0};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 2.09338450817452e-07
    double fks_g = 1.239458197706557e-15;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, -1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 1.23967578310109e-15
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, 0, -13, 14, 1
// realpdgs = 2, 0, -13, 14, 1, 0
TEST(QCDSoftCollinearISR2, Wj_born_2_real_31) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.5493304195846385;
    ps.X2 = 0.5119538935535672;
    ps.Jacobian = 1350.307093689597;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3447.03431115021, 0, 0, 3447.03431115021);
    ps.Momenta[1].Set(3447.03431115021, 0, 0, -3447.03431115021);
    ps.Momenta[2].Set(3443.273652930334, 1047.592736319852, -3013.64616194475,
                      -1294.843511150967);
    ps.Momenta[3].Set(3062.120030198128, -868.1313862282127, 2708.617873831481,
                      1134.158802437602);
    ps.Momenta[4].Set(388.6749391719572, -179.4613500916395, 305.0282881132693,
                      160.6847087133652);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 3.162907884601914e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(1.580142269233051e-09, -2.692645223154387e-28);
    spin[1][2] =
        std::complex<double>(-4.846443563199304e-12, -1.581445972241161e-09);
    spin[2][1] =
        std::complex<double>(-4.846443563199304e-12, 1.581445972241161e-09);
    spin[2][2] =
        std::complex<double>(1.582765615368864e-09, 2.692645223154387e-28);

    // radiation variables
    double phi = 4.943958970847864;
    int bornpdgs[] = {2, 0, -13, 14, 1};
    int realpdgs[] = {2, 0, -13, 14, 1, 0};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 9940.630314217671
    double fks_g = 4.735498840304087e-15;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, -1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 4.736618954945576e-15
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, 0, -13, 14, 1
// realpdgs = 0, 0, -13, 14, 1, -2
TEST(QCDCollinearISR1, Wj_born_2_real_26) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.2716942749649824;
    ps.X2 = 0.04101078130802449;
    ps.Jacobian = 359.2010517480393;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(686.124017462599, 0, 0, 686.124017462599);
    ps.Momenta[1].Set(686.124017462599, 0, 0, -686.124017462599);
    ps.Momenta[2].Set(248.9582942573876, -11.77539822745925, -248.2521773345484,
                      14.57493481250191);
    ps.Momenta[3].Set(603.850716309844, -340.9757878675783, 352.6540302774654,
                      352.1453316530781);
    ps.Momenta[4].Set(519.4390243579666, 352.7511860950376, -104.401852942917,
                      -366.7202664655799);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.978067184767551e-07;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(8.633761917052647e-08, 0);
    spin[1][2] =
        std::complex<double>(4.784876912782129e-08, 8.564153161950816e-08);
    spin[2][1] =
        std::complex<double>(4.784876912782129e-08, -8.564153161950816e-08);
    spin[2][2] = std::complex<double>(1.114690993062286e-07, 0);

    // radiation variables
    double phi = 4.442151311420399;
    double xi = 0.7099248560683751;
    int bornpdgs[] = {2, 0, -13, 14, 1};
    int realpdgs[] = {0, 0, -13, 14, 1, -2};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 2.580338723158978e-05
    double fks_g = 2.600946861557899e-13;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, 1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 2.601467884992388e-13
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, 0, -13, 14, 1
// realpdgs = 0, 0, -13, 14, 1, -2
TEST(QCDSoftCollinearISR1, Wj_born_2_real_26) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.5373491643615758;
    ps.X2 = 0.2573734246602015;
    ps.Jacobian = 2377.893889899067;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2417.260727519692, 0, 0, 2417.260727519692);
    ps.Momenta[1].Set(2417.260727519692, 0, 0, -2417.260727519692);
    ps.Momenta[2].Set(1567.459984991357, -838.6172250836894, 284.3375934923987,
                      -1293.369277224796);
    ps.Momenta[3].Set(2291.01932413014, 690.4916223293042, 89.69906917690159,
                      2182.646315833045);
    ps.Momenta[4].Set(976.042145917888, 148.1256027543853, -374.0366626693003,
                      -889.2770386082493);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 4.265041352876175e-08;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(3.625860355682318e-08, 0);
    spin[1][2] =
        std::complex<double>(1.130227180313715e-08, -1.019886051826774e-08);
    spin[2][1] =
        std::complex<double>(1.130227180313715e-08, 1.019886051826774e-08);
    spin[2][2] = std::complex<double>(6.391809971938567e-09, 0);

    // radiation variables
    double phi = 1.01319551372975;
    int bornpdgs[] = {2, 0, -13, 14, 1};
    int realpdgs[] = {0, 0, -13, 14, 1, -2};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 0.01169872013829285
    double fks_g = 5.008125067496244e-21;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, 1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, 0, -13, 14, 1
// realpdgs = 0, 0, -13, 14, 1, -2
TEST(QCDSoftLimit, Wj_born_2_real_26_0) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.3248988346351176;
    ps.X2 = 0.06976190867322174;
    ps.Jacobian = 712.4571864598287;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(978.5806198577689, 0, 0, 978.5806198577689);
    ps.Momenta[1].Set(978.5806198577689, 0, 0, -978.5806198577689);
    ps.Momenta[2].Set(460.4310147229306, -197.648229271071, -413.7594154376089,
                      -41.65384641933773);
    ps.Momenta[3].Set(774.3568054347287, 659.2533567168592, -5.603014713787843,
                      406.1798616459816);
    ps.Momenta[4].Set(722.3734195578785, -461.6051274457882, 419.3624301513967,
                      -364.5260152266438);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 4.301890052484555e-08;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, 6.452835078726833e-08);
    ColorCorr.Set(0, 4, -7.169816754140925e-09);
    ColorCorr.Set(1, 0, 6.452835078726833e-08);
    ColorCorr.Set(1, 4, 6.452835078726833e-08);
    ColorCorr.Set(4, 0, -7.169816754140925e-09);
    ColorCorr.Set(4, 1, 6.452835078726833e-08);

    // radiation variables
    double y = -0.9924472175164922;
    double phi = 1.971817627765092;
    int bornpdgs[] = {2, 0, -13, 14, 1};
    int realpdgs[] = {0, 0, -13, 14, 1, -2};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 0.0001516186745517742
    double fks_g = 1.982304325817053e-22;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 0, alphas, ColorCorr, y, phi);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, 0, -13, 14, 1
// realpdgs = 0, 0, -13, 14, 1, -2
TEST(QCDCollinearISR2, Wj_born_2_real_26) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.46408041839587;
    ps.X2 = 0.6130212690402459;
    ps.Jacobian = 881.8299780992559;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3466.951370681388, 0, 0, 3466.951370681388);
    ps.Momenta[1].Set(3466.951370681388, 0, 0, -3466.951370681388);
    ps.Momenta[2].Set(3413.989122451107, -585.5016796032095, 1749.369510650855,
                      2872.667023273231);
    ps.Momenta[3].Set(3267.544213171275, 521.3698706870169, -1843.070114803864,
                      -2647.283739021143);
    ps.Momenta[4].Set(252.3694057403936, 64.1318089161926, 93.70060415300858,
                      -225.3832842520881);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 8.206778279567344e-10;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(5.53879510112962e-10, -4.308232357047019e-27);
    spin[1][2] =
        std::complex<double>(-5.75991855204986e-11, -3.800742761798622e-10);
    spin[2][1] =
        std::complex<double>(-5.75991855204986e-11, 3.800742761798622e-10);
    spin[2][2] =
        std::complex<double>(2.667983178437724e-10, 4.308232357047019e-27);

    // radiation variables
    double phi = 2.111650368988298;
    double xi = 0.09047085337046804;
    int bornpdgs[] = {2, 0, -13, 14, 1};
    int realpdgs[] = {0, 0, -13, 14, 1, -2};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 3.43092528014217e-07
    double fks_g = 5.616407738439771e-17;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, -1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-7) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 2, 0, -13, 14, 1
// realpdgs = 0, 0, -13, 14, 1, -2
TEST(QCDSoftCollinearISR2, Wj_born_2_real_26) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.898007467535944;
    ps.X2 = 0.6673236838687515;
    ps.Jacobian = 18141.83648105767;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(5031.779483513688, 0, 0, 5031.779483513688);
    ps.Momenta[1].Set(5031.779483513688, 0, 0, -5031.779483513688);
    ps.Momenta[2].Set(2886.135653303078, -2367.730364384875, 1073.74620148683,
                      -1253.276117074627);
    ps.Momenta[3].Set(3600.091337148059, 2788.049464943166, -419.8056617013511,
                      -2238.571201289366);
    ps.Momenta[4].Set(3577.331976576239, -420.3191005582914, -653.9405397854792,
                      3491.847318363993);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 3.841465337339535e-10;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(2.518961354981078e-11, 0);
    spin[1][2] =
        std::complex<double>(9.483737252792148e-11, 6.918010730627876e-12);
    spin[2][1] =
        std::complex<double>(9.483737252792148e-11, -6.918010730627876e-12);
    spin[2][2] = std::complex<double>(3.589569201841428e-10, 0);

    // radiation variables
    double phi = 4.362348938352135;
    int bornpdgs[] = {2, 0, -13, 14, 1};
    int realpdgs[] = {0, 0, -13, 14, 1, -2};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 0.001658989789085926
    double fks_g = 3.672126392856574e-22;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, -1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 2, -13, 14, 0
// realpdgs = -1, 2, -13, 14, 0, 0
TEST(QCDCollinearFSR, Wj_born_5_real_21) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.7140069660625166;
    ps.X2 = 0.3350350816829888;
    ps.Jacobian = 5492.404454075345;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3179.140512437651, 0, 0, 3179.140512437651);
    ps.Momenta[1].Set(3179.140512437651, 0, 0, -3179.140512437651);
    ps.Momenta[2].Set(1952.83162390325, -1455.285743039469, 92.11977892993116,
                      1298.925980859789);
    ps.Momenta[3].Set(2691.285347114948, 1563.65742687891, 1608.229897859524,
                      -1487.141172381694);
    ps.Momenta[4].Set(1714.164053857105, -108.3716838394409, -1700.349676789454,
                      188.2151915219045);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 3.111517755303003e-08;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(1.273496923343392e-09, -1.102907483404037e-24);
    spin[1][2] =
        std::complex<double>(4.526355007374315e-10, -4.195622926857606e-10);
    spin[1][3] =
        std::complex<double>(4.822403685856375e-09, -3.79035615028051e-09);
    spin[2][1] =
        std::complex<double>(4.526355007374315e-10, 4.195622926857606e-10);
    spin[2][2] =
        std::complex<double>(2.991066621280298e-10, 2.29772392375841e-26);
    spin[2][3] =
        std::complex<double>(2.962772468823241e-09, 2.415781200562685e-10);
    spin[3][1] =
        std::complex<double>(4.822403685856375e-09, 3.79035615028051e-09);
    spin[3][2] =
        std::complex<double>(2.962772468823241e-09, -2.415781200562685e-10);
    spin[3][3] = std::complex<double>(2.95425739675586e-08, 0);

    // radiation variables
    double phi = 5.610367429203538;
    double xi = 0.2436093611221547;
    int bornpdgs[] = {-1, 2, -13, 14, 0};
    int realpdgs[] = {-1, 2, -13, 14, 0, 0};

    // FKS g = real matrix element * (1.0 - y) * xi * xi
    // real = 4.580551744890834e-05
    double fks_g = 2.718352303384737e-14;
    double limit = FKS::QCD::CollinearLimitFSR(
        realpdgs[4], bornpdgs[4], ps, 4, xi, phi, alphas, born,
        spin); // limit = 2.718356152808246e-14
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 2, -13, 14, 0
// realpdgs = -1, 2, -13, 14, 0, 0
TEST(QCDSoftCollinearFSR, Wj_born_5_real_21) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.1932731337073754;
    ps.X2 = 0.1689652565384776;
    ps.Jacobian = 98.7745787899064;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1174.621124085093, 0, 0, 1174.621124085093);
    ps.Momenta[1].Set(1174.621124085093, 0, 0, -1174.621124085093);
    ps.Momenta[2].Set(1103.823207450129, 724.7366192517022, -650.7648760599072,
                      519.3145310000838);
    ps.Momenta[3].Set(1161.984289380631, -730.0197560726353, 677.9877209717906,
                      -598.0060992352306);
    ps.Momenta[4].Set(83.43475133942643, 5.283136820933232, -27.22284491188343,
                      78.69156823514678);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 8.743171681860516e-05;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(3.088257232799647e-06, 3.011672701348623e-21);
    spin[1][2] =
        std::complex<double>(-1.481910709639075e-05, -3.33023538333176e-06);
    spin[1][3] =
        std::complex<double>(-5.33391250333837e-06, -1.152073638812238e-06);
    spin[2][1] =
        std::complex<double>(-1.481910709639075e-05, 3.33023538333176e-06);
    spin[2][2] = std::complex<double>(7.470116167543346e-05, 0);
    spin[2][3] =
        std::complex<double>(2.683730362408645e-05, -2.235829018361916e-07);
    spin[3][1] =
        std::complex<double>(-5.33391250333837e-06, 1.152073638812238e-06);
    spin[3][2] =
        std::complex<double>(2.683730362408645e-05, 2.235829018361916e-07);
    spin[3][3] = std::complex<double>(9.642297910372053e-06, 0);

    // radiation variables
    double phi = 2.686186968259284;
    int bornpdgs[] = {-1, 2, -13, 14, 0};
    int realpdgs[] = {-1, 2, -13, 14, 0, 0};

    // FKS g = real matrix element * (1.0 - y) * xi * xi
    // real = 56067576106.0399
    double fks_g = 2.828851302447262e-10;
    double limit = FKS::QCD::SoftCollinearLimitFSR(realpdgs[4], bornpdgs[4], ps,
                                                   4, phi, alphas, born, spin);
    // limit = 2.818944456298923e-10
    EXPECT_NEAR(limit / fks_g, 1.0, 5e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 2, -13, 14, 0
// realpdgs = -1, 2, -13, 14, 0, 0
TEST(QCDSoftLimit, Wj_born_5_real_21_4) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.1141202038556379;
    ps.X2 = 0.2083632180375901;
    ps.Jacobian = 349.5735590901243;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1002.317133349126, 0, 0, 1002.317133349126);
    ps.Momenta[1].Set(1002.317133349126, 0, 0, -1002.317133349126);
    ps.Momenta[2].Set(969.9613350354595, -171.6861386815058, -898.9350507831126,
                      321.3481534441548);
    ps.Momenta[3].Set(688.6275811003942, 262.5473898790415, 566.4574809176009,
                      -290.5214894502097);
    ps.Momenta[4].Set(346.0453505623988, -90.86125119753572, 332.4775698655117,
                      -30.82666399394506);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    // double born = 5.064640633104903e-07;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, -8.441067721841505e-08);
    ColorCorr.Set(0, 4, 7.596960949657354e-07);
    ColorCorr.Set(1, 0, -8.441067721841505e-08);
    ColorCorr.Set(1, 4, 7.596960949657354e-07);
    ColorCorr.Set(4, 0, 7.596960949657354e-07);
    ColorCorr.Set(4, 1, 7.596960949657354e-07);

    // radiation variables
    double y = 0.01901059708843178;
    double phi = 6.250535619214387;
    int bornpdgs[] = {-1, 2, -13, 14, 0};
    int realpdgs[] = {-1, 2, -13, 14, 0, 0};

    // FKS g = real matrix element * (1.0 - y) * xi * xi
    // real = 25740747.56690777
    double fks_g = 3.009824688445307e-10;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 4, alphas, ColorCorr, y, phi);
    // limit = 3.009822609492805e-10
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 2, -13, 14, 0
// realpdgs = -1, 2, -13, 14, 0, 0
TEST(QCDCollinearISR1, Wj_born_5_real_21) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.2955827910734392;
    ps.X2 = 0.4354549847491072;
    ps.Jacobian = 2041.968105931397;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2331.978610678501, 0, 0, 2331.978610678501);
    ps.Momenta[1].Set(2331.978610678501, 0, 0, -2331.978610678501);
    ps.Momenta[2].Set(1846.849240724905, 300.7794429479623, -31.0933519597141,
                      1821.926740604567);
    ps.Momenta[3].Set(1948.299626731319, -1124.373147624833, 168.9834904897455,
                      -1582.119161240241);
    ps.Momenta[4].Set(868.8083539007786, 823.5937046768703, -137.8901385300314,
                      -239.807579364326);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.635437792018212e-07;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(1.209971093976014e-08, 7.352716556026912e-25);
    spin[1][2] =
        std::complex<double>(-2.069473579730214e-09, -8.672222075069684e-10);
    spin[1][3] =
        std::complex<double>(4.274512834065342e-08, 4.986555914806281e-10);
    spin[2][1] =
        std::complex<double>(-2.069473579730214e-09, 8.672222075069684e-10);
    spin[2][2] =
        std::complex<double>(4.161087218910514e-10, 2.205814966808074e-24);
    spin[2][3] =
        std::complex<double>(-7.346651453786801e-09, 2.978382720646271e-09);
    spin[3][1] =
        std::complex<double>(4.274512834065342e-08, -4.986555914806281e-10);
    spin[3][2] =
        std::complex<double>(-7.346651453786801e-09, -2.978382720646271e-09);
    spin[3][3] =
        std::complex<double>(1.5102795954017e-07, 7.352716556026912e-25);

    // radiation variables
    double phi = 2.068260883723013;
    double xi = 0.0329271967931884;
    int bornpdgs[] = {-1, 2, -13, 14, 0};
    int realpdgs[] = {-1, 2, -13, 14, 0, 0};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 0.005306147605041111
    double fks_g = 1.150585352403205e-13;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, 1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 1.150662756062159e-13
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 2, -13, 14, 0
// realpdgs = -1, 2, -13, 14, 0, 0
TEST(QCDSoftCollinearISR1, Wj_born_5_real_21) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.9815916754408427;
    ps.X2 = 0.9678971094992761;
    ps.Jacobian = 17546.95796835647;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(6335.682223864069, 0, 0, 6335.682223864069);
    ps.Momenta[1].Set(6335.682223864069, 0, 0, -6335.682223864069);
    ps.Momenta[2].Set(3940.551597966761, 2307.691973587672, 1139.090385215021,
                      -2984.120933472388);
    ps.Momenta[3].Set(5982.867788425233, -3387.972338755544, -3642.89013217317,
                      3323.50746810036);
    ps.Momenta[4].Set(2747.945061336145, 1080.280365167872, 2503.799746958149,
                      -339.3865346279721);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 7.087033554397144e-10;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(1.322161857044793e-10, 5.744309809396025e-27);
    spin[1][2] =
        std::complex<double>(-6.940139010801051e-11, 3.373958517906099e-11);
    spin[1][3] =
        std::complex<double>(-9.115456996829437e-11, 2.489113627516419e-10);
    spin[2][1] =
        std::complex<double>(-6.940139010801051e-11, -3.373958517906099e-11);
    spin[2][2] = std::complex<double>(4.503921002750293e-11, 0);
    spin[2][3] =
        std::complex<double>(1.113662439862222e-10, -1.073943945295363e-10);
    spin[3][1] =
        std::complex<double>(-9.115456996829437e-11, -2.489113627516419e-10);
    spin[3][2] =
        std::complex<double>(1.113662439862222e-10, 1.073943945295363e-10);
    spin[3][3] = std::complex<double>(5.314479597077322e-10, 0);

    // radiation variables
    double phi = 3.066114741923329;
    int bornpdgs[] = {-1, 2, -13, 14, 0};
    int realpdgs[] = {-1, 2, -13, 14, 0, 0};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 103331.5709707144
    double fks_g = 7.003123199838378e-17;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, 1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 6.981332453905548e-17
    EXPECT_NEAR(limit / fks_g, 1.0, 5e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 2, -13, 14, 0
// realpdgs = -1, 2, -13, 14, 0, 0
TEST(QCDSoftLimit, Wj_born_5_real_21_0) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.5270808780048775;
    ps.X2 = 0.380983867168986;
    ps.Jacobian = 3479.581881833191;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2912.763876244423, 0, 0, 2912.763876244423);
    ps.Momenta[1].Set(2912.763876244423, 0, 0, -2912.763876244423);
    ps.Momenta[2].Set(1734.943992153746, 1441.671664662246, -926.0693049230974,
                      -272.0461536203309);
    ps.Momenta[3].Set(2905.302611364137, -2271.287972700735, 1743.169873204466,
                      493.3487628611077);
    ps.Momenta[4].Set(1185.281148970962, 829.616308038489, -817.1005682813684,
                      -221.3026092407768);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    // double born = 1.750187876338254e-08;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, -2.91697979389709e-09);
    ColorCorr.Set(0, 4, 2.625281814507381e-08);
    ColorCorr.Set(1, 0, -2.91697979389709e-09);
    ColorCorr.Set(1, 4, 2.625281814507381e-08);
    ColorCorr.Set(4, 0, 2.625281814507381e-08);
    ColorCorr.Set(4, 1, 2.625281814507381e-08);

    // radiation variables
    double y = 0.769530977122308;
    double phi = 0.9240097150415104;
    int bornpdgs[] = {-1, 2, -13, 14, 0};
    int realpdgs[] = {-1, 2, -13, 14, 0, 0};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 689.7343626617319
    double fks_g = 7.543429493291522e-15;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 0, alphas, ColorCorr, y, phi);
    // limit = 7.54342898218582e-15
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 2, -13, 14, 0
// realpdgs = -1, 2, -13, 14, 0, 0
TEST(QCDCollinearISR2, Wj_born_5_real_21) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.7371694117050289;
    ps.X2 = 0.6937193831286663;
    ps.Jacobian = 1092.454975318548;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(4648.244074751172, 0, 0, 4648.244074751172);
    ps.Momenta[1].Set(4648.244074751172, 0, 0, -4648.244074751172);
    ps.Momenta[2].Set(4444.747974498943, -3569.730934443255, 1646.059789265835,
                      -2074.437943794284);
    ps.Momenta[3].Set(4618.547868019406, 3800.867137452328, -1616.053619773162,
                      2067.066547220868);
    ps.Momenta[4].Set(233.1923069839932, -231.1362030090729, -30.00616949267223,
                      7.371396573415715);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.631660372816551e-07;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(1.937318843299091e-10, 9.19089569503364e-26);
    spin[1][2] =
        std::complex<double>(-1.261630914620783e-10, -1.896022297759807e-10);
    spin[1][3] =
        std::complex<double>(5.561060327091111e-09, -7.717990188405259e-10);
    spin[2][1] =
        std::complex<double>(-1.261630914620783e-10, 1.896022297759807e-10);
    spin[2][2] = std::complex<double>(2.6772119293988e-10, 0);
    spin[2][3] =
        std::complex<double>(-2.866155716610238e-09, 5.945133874701711e-09);
    spin[3][1] =
        std::complex<double>(5.561060327091111e-09, 7.717990188405259e-10);
    spin[3][2] =
        std::complex<double>(-2.866155716610238e-09, -5.945133874701711e-09);
    spin[3][3] =
        std::complex<double>(1.627045842043853e-07, 2.941086622410765e-24);

    // radiation variables
    double phi = 5.24707316341762;
    double xi = 0.03351738318458638;
    int bornpdgs[] = {-1, 2, -13, 14, 0};
    int realpdgs[] = {-1, 2, -13, 14, 0, 0};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 0.001285196220628203
    double fks_g = 2.887617360008505e-14;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, -1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 2.887745404746219e-14
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 2, -13, 14, 0
// realpdgs = -1, 2, -13, 14, 0, 0
TEST(QCDSoftCollinearISR2, Wj_born_5_real_21) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.9883039101214415;
    ps.X2 = 0.5075506104788658;
    ps.Jacobian = 16248.4834393184;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(4603.607518664388, 0, 0, 4603.607518664388);
    ps.Momenta[1].Set(4603.607518664388, 0, 0, -4603.607518664388);
    ps.Momenta[2].Set(2372.989434298853, -1089.181754529702, 1698.924674916027,
                      1248.365695561626);
    ps.Momenta[3].Set(3332.241761968226, -1852.183268301521, -2749.689544574955,
                      335.3501293934068);
    ps.Momenta[4].Set(3501.983841061699, 2941.365022831222, 1050.764869658929,
                      -1583.715824955032);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 6.701415019497999e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(1.845252737595532e-09, 0);
    spin[1][2] =
        std::complex<double>(-8.866280611337678e-10, -2.83429070806174e-10);
    spin[1][3] =
        std::complex<double>(2.838845309689029e-09, -1.880497157068294e-10);
    spin[2][1] =
        std::complex<double>(-8.866280611337678e-10, 2.83429070806174e-10);
    spin[2][2] = std::complex<double>(4.69551590042292e-10, 0);
    spin[2][3] =
        std::complex<double>(-1.33515648366438e-09, 5.264002178841051e-10);
    spin[3][1] =
        std::complex<double>(2.838845309689029e-09, 1.880497157068294e-10);
    spin[3][2] =
        std::complex<double>(-1.33515648366438e-09, -5.264002178841051e-10);
    spin[3][3] = std::complex<double>(4.386610691860174e-09, 0);

    // radiation variables
    double phi = 3.859678493712229;
    int bornpdgs[] = {-1, 2, -13, 14, 0};
    int realpdgs[] = {-1, 2, -13, 14, 0, 0};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 2572.776043656853
    double fks_g = 1.247829635804024e-15;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, -1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 1.250346943812107e-15
    EXPECT_NEAR(limit / fks_g, 1.0, 5e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 2, -13, 14, 0
// realpdgs = -1, 2, -13, 14, -2, 2
TEST(QCDCollinearFSR, Wj_born_5_real_6) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.8345358254526793;
    ps.X2 = 0.0835927732879469;
    ps.Jacobian = 1155.123186746572;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1716.802021620234, 0, 0, 1716.802021620234);
    ps.Momenta[1].Set(1716.802021620234, 0, 0, -1716.802021620234);
    ps.Momenta[2].Set(1696.041634433575, -130.1707911254569, -1627.973769747801,
                      457.5086839427422);
    ps.Momenta[3].Set(1069.975699526718, -40.92717723178393, 1057.864223033983,
                      -155.2296665146015);
    ps.Momenta[4].Set(667.5867092801744, 171.0979683572409, 570.1095467138176,
                      -302.2790174281407);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.445079124117254e-07;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(1.037741249178095e-08, 0);
    spin[1][2] =
        std::complex<double>(2.503267217332804e-09, -1.671618063641133e-08);
    spin[1][3] =
        std::complex<double>(1.059514735748173e-08, -3.152734267331769e-08);
    spin[2][1] =
        std::complex<double>(2.503267217332804e-09, 1.671618063641133e-08);
    spin[2][2] =
        std::complex<double>(2.753066258634141e-08, 1.470543311205382e-24);
    spin[2][3] =
        std::complex<double>(5.334077647914849e-08, 9.461803104684725e-09);
    spin[3][1] =
        std::complex<double>(1.059514735748173e-08, 3.152734267331769e-08);
    spin[3][2] =
        std::complex<double>(5.334077647914849e-08, -9.461803104684725e-09);
    spin[3][3] = std::complex<double>(1.06599837333603e-07, 0);

    // radiation variables
    double phi = 3.144051437396133;
    double xi = 0.07185921169140896;
    int bornpdgs[] = {-1, 2, -13, 14, 0};
    int realpdgs[] = {-1, 2, -13, 14, -2, 2};

    // FKS g = real matrix element * (1.0 - y) * xi * xi
    // real = 7.03685994024635e-05
    double fks_g = 3.633655969720387e-15;
    double limit = FKS::QCD::CollinearLimitFSR(
        realpdgs[4], bornpdgs[4], ps, 4, xi, phi, alphas, born,
        spin); // limit = 3.63371468679435e-15
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 2, -13, 14, 0
// realpdgs = -1, 2, -13, 14, -2, 2
TEST(QCDSoftCollinearFSR, Wj_born_5_real_6) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.2527091354409343;
    ps.X2 = 0.09129826312121025;
    ps.Jacobian = 598.9988499160211;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(987.3135227430009, 0, 0, 987.3135227430009);
    ps.Momenta[1].Set(987.3135227430009, 0, 0, -987.3135227430009);
    ps.Momenta[2].Set(407.5502683094778, -6.365542758505995, -59.74517002554173,
                      403.0970301591615);
    ps.Momenta[3].Set(965.112804635671, -258.5451993832311, -73.66424192191896,
                      -926.9146050255574);
    ps.Momenta[4].Set(601.9639725408534, 264.910742141737, 133.4094119474607,
                      523.8175748663958);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.914484090125211e-06;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(1.059710267166041e-06, 0);
    spin[1][2] =
        std::complex<double>(3.887302047433551e-07, -5.745840355580673e-07);
    spin[1][3] =
        std::complex<double>(-6.349326889920235e-07, 1.463389584012184e-07);
    spin[2][1] =
        std::complex<double>(3.887302047433551e-07, 5.745840355580673e-07);
    spin[2][2] = std::complex<double>(4.541410996092567e-07, 0);
    spin[2][3] =
        std::complex<double>(-3.122566174173031e-07, -2.905849108276009e-07);
    spin[3][1] =
        std::complex<double>(-6.349326889920235e-07, -1.463389584012184e-07);
    spin[3][2] =
        std::complex<double>(-3.122566174173031e-07, 2.905849108276009e-07);
    spin[3][3] =
        std::complex<double>(4.006327233499134e-07, -2.941086622410765e-24);

    // radiation variables
    double phi = 5.176516702987613;
    int bornpdgs[] = {-1, 2, -13, 14, 0};
    int realpdgs[] = {-1, 2, -13, 14, -2, 2};

    // FKS g = real matrix element * (1.0 - y) * xi * xi
    // real = 3.915551918325301
    double fks_g = 1.455538895394903e-18;
    double limit = FKS::QCD::SoftCollinearLimitFSR(realpdgs[4], bornpdgs[4], ps,
                                                   4, phi, alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 2, -13, 14, 0
// realpdgs = -1, 2, -13, 14, -2, 2
TEST(QCDSoftLimit, Wj_born_5_real_6_4) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.8004996831027746;
    ps.X2 = 0.02541097273112314;
    ps.Jacobian = 127.522031654645;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(927.0530431888657, 0, 0, 927.0530431888657);
    ps.Momenta[1].Set(927.0530431888657, 0, 0, -927.0530431888657);
    ps.Momenta[2].Set(914.0910845053388, 292.977522119455, -344.4566207032034,
                      -794.403121065508);
    ps.Momenta[3].Set(803.5314752697032, -197.143866309025, 347.2995431460855,
                      697.2662009995572);
    ps.Momenta[4].Set(136.4835266026897, -95.83365581042999, -2.842922442882021,
                      97.13692006595078);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 5.72754300200579e-08;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, -9.545905003342982e-09);
    ColorCorr.Set(0, 4, 8.591314503008685e-08);
    ColorCorr.Set(1, 0, -9.545905003342982e-09);
    ColorCorr.Set(1, 4, 8.591314503008685e-08);
    ColorCorr.Set(4, 0, 8.591314503008685e-08);
    ColorCorr.Set(4, 1, 8.591314503008685e-08);

    // radiation variables
    double y = 0.7515627675309773;
    double phi = 5.931953683899809;
    int bornpdgs[] = {-1, 2, -13, 14, 0};
    int realpdgs[] = {-1, 2, -13, 14, -2, 2};

    // FKS g = real matrix element * (1.0 - y) * xi * xi
    // real = 0.0009128671287512891
    double fks_g = 4.915588904691237e-22;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 4, alphas, ColorCorr, y, phi);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 2, -13, 14, 0
// realpdgs = -1, 2, -13, 14, -1, 1
TEST(QCDCollinearFSR, Wj_born_5_real_5) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.1461904275777606;
    ps.X2 = 0.9583799014262482;
    ps.Jacobian = 3300.519292113687;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2432.99756060156, 0, 0, 2432.99756060156);
    ps.Momenta[1].Set(2432.99756060156, 0, 0, -2432.99756060156);
    ps.Momenta[2].Set(1593.59674941102, 629.7536295507271, -1321.966827058409,
                      -628.7803065921964);
    ps.Momenta[3].Set(1926.413438412829, 511.7505328943618, 1284.865649213374,
                      1341.007230126435);
    ps.Momenta[4].Set(1345.98493337927, -1141.504162445089, 37.1011778450353,
                      -712.2269235342383);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 5.542943378409787e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(8.754887354531755e-10, 0);
    spin[1][2] =
        std::complex<double>(-1.214853587050982e-09, -6.772804292394587e-10);
    spin[1][3] =
        std::complex<double>(-1.466451913246489e-09, -3.528075227974313e-11);
    spin[2][1] =
        std::complex<double>(-1.214853587050982e-09, 6.772804292394587e-10);
    spin[2][2] =
        std::complex<double>(2.209712060772588e-09, -1.838179139006728e-25);
    spin[2][3] =
        std::complex<double>(2.062184534171826e-09, -1.085494529304006e-09);
    spin[3][1] =
        std::complex<double>(-1.466451913246489e-09, 3.528075227974313e-11);
    spin[3][2] =
        std::complex<double>(2.062184534171826e-09, 1.085494529304006e-09);
    spin[3][3] = std::complex<double>(2.457742582184023e-09, 0);

    // radiation variables
    double phi = 5.714434372446001;
    double xi = 0.123092403674553;
    int bornpdgs[] = {-1, 2, -13, 14, 0};
    int realpdgs[] = {-1, 2, -13, 14, -1, 1};

    // FKS g = real matrix element * (1.0 - y) * xi * xi
    // real = 1.16598473872582e-06
    double fks_g = 1.76666975101289e-16;
    double limit = FKS::QCD::CollinearLimitFSR(
        realpdgs[4], bornpdgs[4], ps, 4, xi, phi, alphas, born,
        spin); // limit = 1.766676759755361e-16
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 2, -13, 14, 0
// realpdgs = -1, 2, -13, 14, -1, 1
TEST(QCDSoftCollinearFSR, Wj_born_5_real_5) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.7426813273517787;
    ps.X2 = 0.1251244713203477;
    ps.Jacobian = 3925.443452829919;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1981.46194936347, 0, 0, 1981.46194936347);
    ps.Momenta[1].Set(1981.46194936347, 0, 0, -1981.46194936347);
    ps.Momenta[2].Set(1668.596740400969, -464.9011664623634, 947.9545627326338,
                      1292.077449107852);
    ps.Momenta[3].Set(328.6934773091348, -65.86573670212169, 290.3987948561872,
                      139.1748781222083);
    ps.Momenta[4].Set(1965.633681016835, 530.7669031644854, -1238.353357588821,
                      -1431.252327230061);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 2.577228353120625e-06;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(1.133706967797891e-06, 1.176434648964306e-23);
    spin[1][2] =
        std::complex<double>(9.83908868323701e-08, -9.306308259680095e-07);
    spin[1][3] =
        std::complex<double>(3.352947919797093e-07, 8.052037967641305e-07);
    spin[2][1] =
        std::complex<double>(9.83908868323701e-08, 9.306308259680095e-07);
    spin[2][2] =
        std::complex<double>(7.724698936574614e-07, 5.88217324482153e-24);
    spin[2][3] =
        std::complex<double>(-6.318718530181322e-07, 3.451159743749711e-07);
    spin[3][1] =
        std::complex<double>(3.352947919797093e-07, -8.052037967641305e-07);
    spin[3][2] =
        std::complex<double>(-6.318718530181322e-07, -3.451159743749711e-07);
    spin[3][3] =
        std::complex<double>(6.710514916652724e-07, -2.941086622410765e-24);

    // radiation variables
    double phi = 5.507474678894449;
    int bornpdgs[] = {-1, 2, -13, 14, 0};
    int realpdgs[] = {-1, 2, -13, 14, -1, 1};

    // FKS g = real matrix element * (1.0 - y) * xi * xi
    // real = 0.4945199223154267
    double fks_g = 4.866508529970381e-19;
    double limit = FKS::QCD::SoftCollinearLimitFSR(realpdgs[4], bornpdgs[4], ps,
                                                   4, phi, alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 2, -13, 14, 0
// realpdgs = -1, 2, -13, 14, -1, 1
TEST(QCDSoftLimit, Wj_born_5_real_5_4) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.8263052424456951;
    ps.X2 = 0.381593116115873;
    ps.Jacobian = 10065.14337846077;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3649.924461663117, 0, 0, 3649.924461663117);
    ps.Momenta[1].Set(3649.924461663117, 0, 0, -3649.924461663117);
    ps.Momenta[2].Set(1436.637923700578, -1282.708522536517, 380.5802452098276,
                      -523.207460748527);
    ps.Momenta[3].Set(3127.087565204863, 1938.090172320792, 2226.841685299599,
                      1031.338563724097);
    ps.Momenta[4].Set(2736.123434420792, -655.3816497842745, -2607.421930509426,
                      -508.1311029755705);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 3.254208794448786e-09;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, -5.42368132408131e-10);
    ColorCorr.Set(0, 4, 4.881313191673179e-09);
    ColorCorr.Set(1, 0, -5.42368132408131e-10);
    ColorCorr.Set(1, 4, 4.881313191673179e-09);
    ColorCorr.Set(4, 0, 4.881313191673179e-09);
    ColorCorr.Set(4, 1, 4.881313191673179e-09);

    // radiation variables
    double y = -0.9447938609183169;
    double phi = 1.800695509456996;
    int bornpdgs[] = {-1, 2, -13, 14, 0};
    int realpdgs[] = {-1, 2, -13, 14, -1, 1};

    // FKS g = real matrix element * (1.0 - y) * xi * xi
    // real = 1.657152016334797e-08
    double fks_g = 1.8110879588259e-24;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 4, alphas, ColorCorr, y, phi);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 2, -13, 14, 0
// realpdgs = -1, 2, -13, 14, -5, 5
TEST(QCDCollinearFSR, Wj_born_5_real_37) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.0854831677857284;
    ps.X2 = 0.2379311944563369;
    ps.Jacobian = 802.357059386003;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(926.9991861784049, 0, 0, 926.9991861784049);
    ps.Momenta[1].Set(926.9991861784049, 0, 0, -926.9991861784049);
    ps.Momenta[2].Set(388.0671105427764, 60.64745471409614, 352.2157333361865,
                      -151.2020030031792);
    ps.Momenta[3].Set(607.1393974621262, 141.9716034008882, 32.43347307138774,
                      -589.4152878938432);
    ps.Momenta[4].Set(858.7918643519073, -202.6190581149844, -384.6492064075742,
                      740.6172908970224);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.96625513506708e-06;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(8.047329883118846e-07, 0);
    spin[1][2] =
        std::complex<double>(2.125742357931761e-07, -7.83901242702628e-07);
    spin[1][3] =
        std::complex<double>(3.305631049064897e-07, -4.071292887886442e-07);
    spin[2][1] =
        std::complex<double>(2.125742357931761e-07, 7.83901242702628e-07);
    spin[2][2] =
        std::complex<double>(8.19761304203061e-07, 2.352869297928612e-23);
    spin[2][3] =
        std::complex<double>(4.839100179624398e-07, 2.144607389049696e-07);
    spin[3][1] =
        std::complex<double>(3.305631049064897e-07, 4.071292887886442e-07);
    spin[3][2] =
        std::complex<double>(4.839100179624398e-07, -2.144607389049696e-07);
    spin[3][3] = std::complex<double>(3.417608425521339e-07, 0);

    // radiation variables
    double phi = 5.348267873542412;
    double xi = 0.7740139660756671;
    int bornpdgs[] = {-1, 2, -13, 14, 0};
    int realpdgs[] = {-1, 2, -13, 14, -5, 5};

    // FKS g = real matrix element * (1.0 - y) * xi * xi
    // real = 0.001033491689919175
    double fks_g = 6.191624145009706e-12;
    double limit = FKS::QCD::CollinearLimitFSR(
        realpdgs[4], bornpdgs[4], ps, 4, xi, phi, alphas, born,
        spin); // limit = 6.19151881319558e-12
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 2, -13, 14, 0
// realpdgs = -1, 2, -13, 14, -5, 5
TEST(QCDSoftCollinearFSR, Wj_born_5_real_37) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.02924533232778614;
    ps.X2 = 0.6658615786806266;
    ps.Jacobian = 452.8062446919553;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(907.054986319249, 0, 0, 907.054986319249);
    ps.Momenta[1].Set(907.054986319249, 0, 0, -907.054986319249);
    ps.Momenta[2].Set(431.1263275667946, 90.19199440070867, 411.2589355090781,
                      92.74374604895701);
    ps.Momenta[3].Set(887.672170737959, -412.5749652417981, -784.5409455377937,
                      -47.32108972787904);
    ps.Momenta[4].Set(495.3114743337446, 322.3829708410895, 373.2820100287157,
                      -45.42265632107797);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.912420617820632e-07;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(1.163680958777752e-08, 2.941086622410765e-24);
    spin[1][2] =
        std::complex<double>(-7.061111028844379e-09, 4.578735676789709e-09);
    spin[1][3] =
        std::complex<double>(2.456315016563854e-08, 3.762791072236642e-08);
    spin[2][1] =
        std::complex<double>(-7.061111028844379e-09, -4.578735676789709e-09);
    spin[2][2] =
        std::complex<double>(6.086213650343071e-09, 7.352716556026912e-25);
    spin[2][3] =
        std::complex<double>(-9.924311826899249e-11, -3.249713974773817e-08);
    spin[3][1] =
        std::complex<double>(2.456315016563854e-08, -3.762791072236642e-08);
    spin[3][2] =
        std::complex<double>(-9.924311826899249e-11, 3.249713974773817e-08);
    spin[3][3] =
        std::complex<double>(1.735190385439426e-07, 1.470543311205382e-24);

    // radiation variables
    double phi = 4.08065923819759;
    int bornpdgs[] = {-1, 2, -13, 14, 0};
    int realpdgs[] = {-1, 2, -13, 14, -5, 5};

    // FKS g = real matrix element * (1.0 - y) * xi * xi
    // real = 0.5779813101803727
    double fks_g = 1.723468470705832e-19;
    double limit = FKS::QCD::SoftCollinearLimitFSR(realpdgs[4], bornpdgs[4], ps,
                                                   4, phi, alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 2, -13, 14, 0
// realpdgs = -1, 2, -13, 14, -5, 5
TEST(QCDSoftLimit, Wj_born_5_real_37_4) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.1324368033482877;
    ps.X2 = 0.5016974430081547;
    ps.Jacobian = 409.6720488232385;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1675.477674157556, 0, 0, 1675.477674157556);
    ps.Momenta[1].Set(1675.477674157556, 0, 0, -1675.477674157556);
    ps.Momenta[2].Set(1483.675400482303, -707.2065451391578, -1017.22667425424,
                      816.3341777080742);
    ps.Momenta[3].Set(1624.676323279985, 949.2597340600943, 1033.350279450227,
                      -818.9421912924428);
    ps.Momenta[4].Set(242.6036245528225, -242.0531889209365, -16.12360519598683,
                      2.608013584368647);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.196903372850883e-06;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, -1.994838954751471e-07);
    ColorCorr.Set(0, 4, 1.795355059276324e-06);
    ColorCorr.Set(1, 0, -1.994838954751471e-07);
    ColorCorr.Set(1, 4, 1.795355059276324e-06);
    ColorCorr.Set(4, 0, 1.795355059276324e-06);
    ColorCorr.Set(4, 1, 1.795355059276324e-06);

    // radiation variables
    double y = -0.8460607800482975;
    double phi = 0.4401584486846244;
    int bornpdgs[] = {-1, 2, -13, 14, 0};
    int realpdgs[] = {-1, 2, -13, 14, -5, 5};

    // FKS g = real matrix element * (1.0 - y) * xi * xi
    // real = 0.000816733568986179
    double fks_g = 3.1611396452783e-21;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 4, alphas, ColorCorr, y, phi);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 2, -13, 14, 0
// realpdgs = -1, 2, -13, 14, -4, 4
TEST(QCDCollinearFSR, Wj_born_5_real_37_n1) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.8184538147407601;
    ps.X2 = 0.1903833306503984;
    ps.Jacobian = 1670.490825014996;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2565.812434031569, 0, 0, 2565.812434031569);
    ps.Momenta[1].Set(2565.812434031569, 0, 0, -2565.812434031569);
    ps.Momenta[2].Set(2240.427969939166, 1650.586172673372, -348.6982278098515,
                      1474.276880707104);
    ps.Momenta[3].Set(2245.217324243026, -2000.281040003195, 53.53359603758005,
                      -1018.337246780277);
    ps.Momenta[4].Set(645.9795738809484, 349.6948673298223, 295.1646317722714,
                      -455.939633926827);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 3.556752943789873e-07;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(9.405904873462333e-08, 1.176434648964306e-23);
    spin[1][2] =
        std::complex<double>(8.748379418103642e-08, -1.616236572258972e-08);
    spin[1][3] =
        std::complex<double>(1.287760135222926e-07, -1.046313672270615e-08);
    spin[2][1] =
        std::complex<double>(8.748379418103642e-08, 1.616236572258972e-08);
    spin[2][2] = std::complex<double>(8.414540032603259e-08, 0);
    spin[2][3] =
        std::complex<double>(1.215717515602221e-07, 1.239615053514771e-08);
    spin[3][1] =
        std::complex<double>(1.287760135222926e-07, 1.046313672270615e-08);
    spin[3][2] =
        std::complex<double>(1.215717515602221e-07, -1.239615053514771e-08);
    spin[3][3] = std::complex<double>(1.774708453183314e-07, 0);

    // radiation variables
    double phi = 0.2101347053417509;
    double xi = 0.0367859411573524;
    int bornpdgs[] = {-1, 2, -13, 14, 0};
    int realpdgs[] = {-1, 2, -13, 14, -4, 4};

    // FKS g = real matrix element * (1.0 - y) * xi * xi
    // real = 0.0002700987994585121
    double fks_g = 3.65499173848616e-15;
    double limit = FKS::QCD::CollinearLimitFSR(
        realpdgs[4], bornpdgs[4], ps, 4, xi, phi, alphas, born,
        spin); // limit = 3.654693365747288e-15
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 2, -13, 14, 0
// realpdgs = -1, 2, -13, 14, -4, 4
TEST(QCDSoftCollinearFSR, Wj_born_5_real_37_n1) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.4299293039017051;
    ps.X2 = 0.9841883145061558;
    ps.Jacobian = 8884.643427615691;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(4228.155806225874, 0, 0, 4228.155806225874);
    ps.Momenta[1].Set(4228.155806225874, 0, 0, -4228.155806225874);
    ps.Momenta[2].Set(3753.05227688335, 3115.413520065219, -2036.729941769335,
                      -480.9689556659688);
    ps.Momenta[3].Set(2618.343074597532, -1580.340954032455, 1264.093361682146,
                      1661.418339326907);
    ps.Momenta[4].Set(2084.916260970865, -1535.072566032765, 772.6365800871889,
                      -1180.449383660938);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 6.514397797224514e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(1.477206183441383e-09, 0);
    spin[1][2] =
        std::complex<double>(-1.132858668669715e-10, -1.553576584796747e-09);
    spin[1][3] =
        std::complex<double>(-1.995127892751404e-09, -1.016858593003147e-09);
    spin[2][1] =
        std::complex<double>(-1.132858668669715e-10, 1.553576584796747e-09);
    spin[2][2] =
        std::complex<double>(1.642583086680403e-09, -9.19089569503364e-26);
    spin[2][3] =
        std::complex<double>(1.222434290637076e-09, -2.020292295088669e-09);
    spin[3][1] =
        std::complex<double>(-1.995127892751404e-09, 1.016858593003147e-09);
    spin[3][2] =
        std::complex<double>(1.222434290637076e-09, 2.020292295088669e-09);
    spin[3][3] =
        std::complex<double>(3.394608527102728e-09, -4.59544784751682e-26);

    // radiation variables
    double phi = 3.567406027272053;
    int bornpdgs[] = {-1, 2, -13, 14, 0};
    int realpdgs[] = {-1, 2, -13, 14, -4, 4};

    // FKS g = real matrix element * (1.0 - y) * xi * xi
    // real = 0.001110181433684766
    double fks_g = 2.699412335195566e-22;
    double limit = FKS::QCD::SoftCollinearLimitFSR(realpdgs[4], bornpdgs[4], ps,
                                                   4, phi, alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 2, -13, 14, 0
// realpdgs = -1, 2, -13, 14, -4, 4
TEST(QCDSoftLimit, Wj_born_5_real_37_4_n1) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.5302299447420005;
    ps.X2 = 0.5857611968361098;
    ps.Jacobian = 6223.471231252225;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3622.47820794483, 0, 0, 3622.47820794483);
    ps.Momenta[1].Set(3622.47820794483, 0, 0, -3622.47820794483);
    ps.Momenta[2].Set(2275.796856623159, -2086.714295938104, 680.9105817257598,
                      601.0287509240142);
    ps.Momenta[3].Set(3264.543792216175, 1963.33665698685, -1925.842488164343,
                      -1758.887731894744);
    ps.Momenta[4].Set(1704.615767050326, 123.3776389512539, 1244.931906438584,
                      1157.85898097073);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 4.623876077451274e-08;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, -7.706460129085456e-09);
    ColorCorr.Set(0, 4, 6.935814116176912e-08);
    ColorCorr.Set(1, 0, -7.706460129085456e-09);
    ColorCorr.Set(1, 4, 6.935814116176912e-08);
    ColorCorr.Set(4, 0, 6.935814116176912e-08);
    ColorCorr.Set(4, 1, 6.935814116176912e-08);

    // radiation variables
    double y = 0.447843054683311;
    double phi = 2.310525984320454;
    int bornpdgs[] = {-1, 2, -13, 14, 0};
    int realpdgs[] = {-1, 2, -13, 14, -4, 4};

    // FKS g = real matrix element * (1.0 - y) * xi * xi
    // real = 2.13674422784234e-06
    double fks_g = 2.612501317096704e-23;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 4, alphas, ColorCorr, y, phi);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 2, -13, 14, 0
// realpdgs = -1, 2, -13, 14, -3, 3
TEST(QCDCollinearFSR, Wj_born_5_real_37_n2) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.6959044703539412;
    ps.X2 = 0.3854025141216511;
    ps.Jacobian = 597.6783835205159;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3366.242830895854, 0, 0, 3366.242830895854);
    ps.Momenta[1].Set(3366.242830895854, 0, 0, -3366.242830895854);
    ps.Momenta[2].Set(3364.959458640019, 642.6069515539424, -1766.260746940719,
                      -2791.116521738547);
    ps.Momenta[3].Set(3191.360361232654, -594.432435612748, 1655.62215086489,
                      2662.770423507071);
    ps.Momenta[4].Set(176.165841919036, -48.17451594119445, 110.6385960758293,
                      128.3460982314766);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 5.839207767127728e-08;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(5.134912076218002e-09, 2.941086622410765e-24);
    spin[1][2] =
        std::complex<double>(-1.14010116233487e-08, 1.743327387476161e-09);
    spin[1][3] =
        std::complex<double>(1.175543194778543e-08, -1.502806063516176e-09);
    spin[2][1] =
        std::complex<double>(-1.14010116233487e-08, -1.743327387476161e-09);
    spin[2][2] = std::complex<double>(2.590545941998501e-08, 0);
    spin[2][3] =
        std::complex<double>(-2.661071839489429e-08, -6.543553265423251e-10);
    spin[3][1] =
        std::complex<double>(1.175543194778543e-08, 1.502806063516176e-09);
    spin[3][2] =
        std::complex<double>(-2.661071839489429e-08, 6.543553265423251e-10);
    spin[3][3] =
        std::complex<double>(2.735170617507425e-08, 3.676358278013456e-25);

    // radiation variables
    double phi = 5.159677970581597;
    double xi = 0.02407774107011108;
    int bornpdgs[] = {-1, 2, -13, 14, 0};
    int realpdgs[] = {-1, 2, -13, 14, -3, 3};

    // FKS g = real matrix element * (1.0 - y) * xi * xi
    // real = 0.0004582738436271467
    double fks_g = 2.656785864742725e-15;
    double limit = FKS::QCD::CollinearLimitFSR(
        realpdgs[4], bornpdgs[4], ps, 4, xi, phi, alphas, born,
        spin); // limit = 2.656803443247798e-15
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 2, -13, 14, 0
// realpdgs = -1, 2, -13, 14, -3, 3
TEST(QCDSoftCollinearFSR, Wj_born_5_real_37_n2) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.8423033507041247;
    ps.X2 = 0.9716574387075454;
    ps.Jacobian = 26170.14176131137;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(5880.364007968215, 0, 0, 5880.364007968215);
    ps.Momenta[1].Set(5880.364007968215, 0, 0, -5880.364007968215);
    ps.Momenta[2].Set(4206.333498174005, -699.2515444748819, -243.973502546459,
                      4140.623830471136);
    ps.Momenta[3].Set(3138.675190921598, -2360.585707068023, 1629.504113188794,
                      -1274.218748407702);
    ps.Momenta[4].Set(4415.71932684083, 3059.837251542904, -1385.530610642335,
                      -2866.405082063433);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.376269285960142e-08;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(4.860192998806242e-09, 1.838179139006728e-25);
    spin[1][2] =
        std::complex<double>(-1.292740860872715e-09, -2.515346950166471e-09);
    spin[1][3] =
        std::complex<double>(5.813041473470945e-09, 1.21584008403051e-09);
    spin[2][1] =
        std::complex<double>(-1.292740860872715e-09, 2.515346950166471e-09);
    spin[2][2] =
        std::complex<double>(1.645644363309482e-09, -1.838179139006728e-25);
    spin[2][3] =
        std::complex<double>(-2.175431281958733e-09, 2.685088840665082e-09);
    spin[3][1] =
        std::complex<double>(5.813041473470945e-09, -1.21584008403051e-09);
    spin[3][2] =
        std::complex<double>(-2.175431281958733e-09, -2.685088840665082e-09);
    spin[3][3] = std::complex<double>(7.2568554974857e-09, 0);

    // radiation variables
    double phi = 2.942232973104738;
    int bornpdgs[] = {-1, 2, -13, 14, 0};
    int realpdgs[] = {-1, 2, -13, 14, -3, 3};

    // FKS g = real matrix element * (1.0 - y) * xi * xi
    // real = 0.0005232335399943431
    double fks_g = 2.950462338074091e-22;
    double limit = FKS::QCD::SoftCollinearLimitFSR(realpdgs[4], bornpdgs[4], ps,
                                                   4, phi, alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 2, -13, 14, 0
// realpdgs = -1, 2, -13, 14, -3, 3
TEST(QCDSoftLimit, Wj_born_5_real_37_4_n2) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.1651056264678097;
    ps.X2 = 0.5419133059978662;
    ps.Jacobian = 3505.084676480029;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1944.281754490929, 0, 0, 1944.281754490929);
    ps.Momenta[1].Set(1944.281754490929, 0, 0, -1944.281754490929);
    ps.Momenta[2].Set(1619.820232801839, 1139.973234065428, 1107.534892627497,
                      -312.4821176022931);
    ps.Momenta[3].Set(480.0374980667111, 2.306937220960776, 18.13324338269899,
                      -479.6893401726057);
    ps.Momenta[4].Set(1788.705778113309, -1142.280171286388, -1125.668136010196,
                      792.1714577748987);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 5.124598309411399e-08;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, -8.540997182352331e-09);
    ColorCorr.Set(0, 4, 7.686897464117098e-08);
    ColorCorr.Set(1, 0, -8.540997182352331e-09);
    ColorCorr.Set(1, 4, 7.686897464117098e-08);
    ColorCorr.Set(4, 0, 7.686897464117098e-08);
    ColorCorr.Set(4, 1, 7.686897464117098e-08);

    // radiation variables
    double y = 0.855779533771063;
    double phi = 3.683233146995973;
    int bornpdgs[] = {-1, 2, -13, 14, 0};
    int realpdgs[] = {-1, 2, -13, 14, -3, 3};

    // FKS g = real matrix element * (1.0 - y) * xi * xi
    // real = 8.234116876411465e-06
    double fks_g = 1.005086271019313e-22;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 4, alphas, ColorCorr, y, phi);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 2, -13, 14, 0
// realpdgs = 0, 2, -13, 14, 0, 1
TEST(QCDCollinearISR1, Wj_born_5_real_39) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.1054593776153858;
    ps.X2 = 0.7185739986916317;
    ps.Jacobian = 2375.059765133873;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1789.335209489305, 0, 0, 1789.335209489305);
    ps.Momenta[1].Set(1789.335209489305, 0, 0, -1789.335209489305);
    ps.Momenta[2].Set(1122.624196297594, 919.601141573731, -287.5704230288821,
                      -576.1267901496166);
    ps.Momenta[3].Set(1139.056362859153, -35.100855538161, 1134.78911936457,
                      92.03793936304808);
    ps.Momenta[4].Set(1316.989859821862, -884.5002860355698, -847.2186963356875,
                      484.0888507865684);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.656177330209444e-08;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(2.353425612451279e-09, 0);
    spin[1][2] =
        std::complex<double>(-3.872376682883033e-09, -1.740415482746988e-09);
    spin[1][3] =
        std::complex<double>(-2.477116123760881e-09, -3.04595434077752e-09);
    spin[2][1] =
        std::complex<double>(-3.872376682883033e-09, 1.740415482746988e-09);
    spin[2][2] =
        std::complex<double>(7.658770743107466e-09, -1.838179139006728e-25);
    spin[2][3] =
        std::complex<double>(6.328457009253844e-09, 3.179990594307563e-09);
    spin[3][1] =
        std::complex<double>(-2.477116123760881e-09, 3.04595434077752e-09);
    spin[3][2] =
        std::complex<double>(6.328457009253844e-09, -3.179990594307563e-09);
    spin[3][3] = std::complex<double>(6.549576946535691e-09, 0);

    // radiation variables
    double phi = 3.705894784624187;
    double xi = 0.359688065348204;
    int bornpdgs[] = {-1, 2, -13, 14, 0};
    int realpdgs[] = {0, 2, -13, 14, 0, 1};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 5.750363161251325e-07
    double fks_g = 1.487912267627172e-15;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, 1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 1.488099717312641e-15
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 2, -13, 14, 0
// realpdgs = 0, 2, -13, 14, 0, 1
TEST(QCDSoftCollinearISR1, Wj_born_5_real_39) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.7430949830168068;
    ps.X2 = 0.9735293733190389;
    ps.Jacobian = 11804.80183935059;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(5528.534842963723, 0, 0, 5528.534842963723);
    ps.Momenta[1].Set(5528.534842963723, 0, 0, -5528.534842963723);
    ps.Momenta[2].Set(5259.840890639833, -1441.351075644012, -5045.172247560678,
                      366.9744732678349);
    ps.Momenta[3].Set(3678.632314252216, 1136.078937694752, 3379.452668163121,
                      906.0684380509787);
    ps.Momenta[4].Set(2118.596481035398, 305.2721379492602, 1665.719579397558,
                      -1273.042911318814);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.091609770953933e-08;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(1.259632405252448e-09, 2.757268708510092e-25);
    spin[1][2] =
        std::complex<double>(5.468392678528722e-10, -1.99823791673971e-09);
    spin[1][3] =
        std::complex<double>(1.017571003388685e-09, -2.614604733755397e-09);
    spin[2][1] =
        std::complex<double>(5.468392678528722e-10, 1.99823791673971e-09);
    spin[2][2] =
        std::complex<double>(3.40733371010882e-09, -9.19089569503364e-26);
    spin[2][3] =
        std::complex<double>(4.589473940692994e-09, 4.791718767299592e-10);
    spin[3][1] =
        std::complex<double>(1.017571003388685e-09, 2.614604733755397e-09);
    spin[3][2] =
        std::complex<double>(4.589473940692994e-09, -4.791718767299592e-10);
    spin[3][3] =
        std::complex<double>(6.249131594178065e-09, -4.59544784751682e-26);

    // radiation variables
    double phi = 0.8781085119472924;
    int bornpdgs[] = {-1, 2, -13, 14, 0};
    int realpdgs[] = {0, 2, -13, 14, 0, 1};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 0.001031444829906781
    double fks_g = 1.36151155449609e-22;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, 1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 2, -13, 14, 0
// realpdgs = 0, 2, -13, 14, 0, 1
TEST(QCDSoftLimit, Wj_born_5_real_39_0) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.7807393511535814;
    ps.X2 = 0.5584780028258578;
    ps.Jacobian = 18342.99752120979;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(4292.096001710807, 0, 0, 4292.096001710807);
    ps.Momenta[1].Set(4292.096001710807, 0, 0, -4292.096001710807);
    ps.Momenta[2].Set(3896.917429531766, -2456.67948433659, -2892.567880741811,
                      -885.2922789277013);
    ps.Momenta[3].Set(446.9365679404309, -44.1362658293906, -314.9144848588965,
                      -314.0591552989092);
    ps.Momenta[4].Set(4240.338005949423, 2500.815750165978, 3207.482365600703,
                      1199.351434226609);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 3.995465300518171e-08;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, -6.659108834196952e-09);
    ColorCorr.Set(0, 4, 5.993197950777257e-08);
    ColorCorr.Set(1, 0, -6.659108834196952e-09);
    ColorCorr.Set(1, 4, 5.993197950777257e-08);
    ColorCorr.Set(4, 0, 5.993197950777257e-08);
    ColorCorr.Set(4, 1, 5.993197950777257e-08);

    // radiation variables
    double y = -0.1770079324563838;
    double phi = 4.9180273483863;
    int bornpdgs[] = {-1, 2, -13, 14, 0};
    int realpdgs[] = {0, 2, -13, 14, 0, 1};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 4.878052035082372e-07
    double fks_g = 9.570241081769632e-24;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 0, alphas, ColorCorr, y, phi);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 2, -13, 14, 0
// realpdgs = 0, 2, -13, 14, 0, 1
TEST(QCDCollinearISR2, Wj_born_5_real_39) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.02359616890737826;
    ps.X2 = 0.8088456116080067;
    ps.Jacobian = 501.4491805260784;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(897.980532428532, 0, 0, 897.980532428532);
    ps.Momenta[1].Set(897.980532428532, 0, 0, -897.980532428532);
    ps.Momenta[2].Set(351.1944855903818, -215.2283806589375, -215.6296592962831,
                      174.6944787329961);
    ps.Momenta[3].Set(890.7030031950753, 660.4458822506895, 528.4854956748545,
                      -279.0450812684898);
    ps.Momenta[4].Set(554.0635760716071, -445.2175015917518, -312.8558363785714,
                      104.3506025354937);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 2.725909811214387e-07;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(1.589232774064456e-08, 0);
    spin[1][2] =
        std::complex<double>(-1.186338972083903e-08, -1.703762572423354e-08);
    spin[1][3] =
        std::complex<double>(3.223758805441282e-08, -5.108088038156842e-08);
    spin[2][1] =
        std::complex<double>(-1.186338972083903e-08, 1.703762572423354e-08);
    spin[2][2] =
        std::complex<double>(2.712130740201388e-08, 2.941086622410765e-24);
    spin[2][3] =
        std::complex<double>(3.069719293596422e-08, 7.269195360341602e-08);
    spin[3][1] =
        std::complex<double>(3.223758805441282e-08, 5.108088038156842e-08);
    spin[3][2] =
        std::complex<double>(3.069719293596422e-08, -7.269195360341602e-08);
    spin[3][3] = std::complex<double>(2.295773459787803e-07, 0);

    // radiation variables
    double phi = 3.427906659716672;
    double xi = 0.02563239061323235;
    int bornpdgs[] = {-1, 2, -13, 14, 0};
    int realpdgs[] = {0, 2, -13, 14, 0, 1};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 9.639421621447566e-12
    double fks_g = 1.266657494944866e-22;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, -1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 2, -13, 14, 0
// realpdgs = 0, 2, -13, 14, 0, 1
TEST(QCDSoftCollinearISR2, Wj_born_5_real_39) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.1610318426268513;
    ps.X2 = 0.4135211503831684;
    ps.Jacobian = 1250.762807315179;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1677.328404421947, 0, 0, 1677.328404421947);
    ps.Momenta[1].Set(1677.328404421947, 0, 0, -1677.328404421947);
    ps.Momenta[2].Set(1205.149302821369, -439.1887156731451, -624.0667597593783,
                      -932.7587005649618);
    ps.Momenta[3].Set(1409.635742085479, 1074.702844179962, 715.971021464863,
                      565.2187350919467);
    ps.Momenta[4].Set(739.8717639370451, -635.5141285068171, -91.90426170548471,
                      367.5399654730152);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.82459390415921e-08;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(2.940056982150723e-09, 0);
    spin[1][2] =
        std::complex<double>(1.239869606935056e-09, -3.677675757626804e-09);
    spin[1][3] =
        std::complex<double>(5.393691129771528e-09, -9.196117621164275e-10);
    spin[2][1] =
        std::complex<double>(1.239869606935056e-09, 3.677675757626804e-09);
    spin[2][2] = std::complex<double>(5.123225744223001e-09, 0);
    spin[2][3] =
        std::complex<double>(3.424936198888527e-09, 6.359076899381284e-09);
    spin[3][1] =
        std::complex<double>(5.393691129771528e-09, 9.196117621164275e-10);
    spin[3][2] =
        std::complex<double>(3.424936198888527e-09, -6.359076899381284e-09);
    spin[3][3] = std::complex<double>(1.018265631521837e-08, 0);

    // radiation variables
    double phi = 0.1674793789834257;
    int bornpdgs[] = {-1, 2, -13, 14, 0};
    int realpdgs[] = {0, 2, -13, 14, 0, 1};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 9.221950259980867e-09
    double fks_g = 6.343918137711897e-27;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, -1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 2, -13, 14, 0
// realpdgs = -1, 0, -13, 14, 0, -2
TEST(QCDCollinearISR1, Wj_born_5_real_14) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.4422314567487042;
    ps.X2 = 0.9935184366274115;
    ps.Jacobian = 12568.14005629592;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(4308.500401405853, 0, 0, 4308.500401405853);
    ps.Momenta[1].Set(4308.500401405853, 0, 0, -4308.500401405853);
    ps.Momenta[2].Set(1812.968555103613, 160.0073402902268, -1672.491975252857,
                      681.1925025736321);
    ps.Momenta[3].Set(3909.726077959639, -826.3447941353115, 3547.582488959015,
                      1420.482583460456);
    ps.Momenta[4].Set(2894.306169748454, 666.3374538450847, -1875.090513706158,
                      -2101.675086034088);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 7.721505400494129e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(7.182921621092545e-10, 0);
    spin[1][2] =
        std::complex<double>(-1.555627922067752e-09, -5.532819036554256e-12);
    spin[1][3] =
        std::complex<double>(1.615648466552433e-09, 4.93631796771808e-12);
    spin[2][1] =
        std::complex<double>(-1.555627922067752e-09, 5.532819036554256e-12);
    spin[2][2] = std::complex<double>(3.369114925181699e-09, 0);
    spin[2][3] =
        std::complex<double>(-3.499098711015364e-09, 1.754183876423829e-12);
    spin[3][1] =
        std::complex<double>(1.615648466552433e-09, -4.93631796771808e-12);
    spin[3][2] =
        std::complex<double>(-3.499098711015364e-09, -1.754183876423829e-12);
    spin[3][3] =
        std::complex<double>(3.634098313203175e-09, 2.29772392375841e-26);

    // radiation variables
    double phi = 2.944719474708537;
    double xi = 0.3377038499667901;
    int bornpdgs[] = {-1, 2, -13, 14, 0};
    int realpdgs[] = {-1, 0, -13, 14, 0, -2};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 7.134551118178152e-15
    double fks_g = 1.627303929014716e-23;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, 1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 2, -13, 14, 0
// realpdgs = -1, 0, -13, 14, 0, -2
TEST(QCDSoftCollinearISR1, Wj_born_5_real_14) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.5621461918767032;
    ps.X2 = 0.4785177872975552;
    ps.Jacobian = 8209.429088638239;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3371.219544423249, 0, 0, 3371.219544423249);
    ps.Momenta[1].Set(3371.219544423249, 0, 0, -3371.219544423249);
    ps.Momenta[2].Set(2943.318153074536, -1184.094117077355, 2684.405733203433,
                      234.5394032256394);
    ps.Momenta[3].Set(1382.961903892202, 815.8285342183319, -432.6738974645971,
                      -1029.466234912592);
    ps.Momenta[4].Set(2416.159031879759, 368.2655828590227, -2251.731835738835,
                      794.9268316869521);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.425612877817179e-08;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(4.356071106537832e-10, 1.838179139006728e-25);
    spin[1][2] =
        std::complex<double>(5.822621539316016e-10, 6.303878869368309e-10);
    spin[1][3] =
        std::complex<double>(1.44752834648856e-09, 1.78565424803624e-09);
    spin[2][1] =
        std::complex<double>(5.822621539316016e-10, -6.303878869368309e-10);
    spin[2][2] = std::complex<double>(1.690555745962212e-09, 0);
    spin[2][3] =
        std::complex<double>(4.518970726093962e-09, 2.920396612068093e-10);
    spin[3][1] =
        std::complex<double>(1.44752834648856e-09, -1.78565424803624e-09);
    spin[3][2] =
        std::complex<double>(4.518970726093962e-09, -2.920396612068093e-10);
    spin[3][3] = std::complex<double>(1.212996592155579e-08, 0);

    // radiation variables
    double phi = 4.136703741819963;
    int bornpdgs[] = {-1, 2, -13, 14, 0};
    int realpdgs[] = {-1, 0, -13, 14, 0, -2};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 2.60450008185399e-09
    double fks_g = 9.986487336320121e-28;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, 1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 2, -13, 14, 0
// realpdgs = -1, 0, -13, 14, 0, -2
TEST(QCDSoftLimit, Wj_born_5_real_14_0) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.06465883461865829;
    ps.X2 = 0.2946621084256549;
    ps.Jacobian = 53.93179068752883;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(897.2003598369831, 0, 0, 897.2003598369831);
    ps.Momenta[1].Set(897.2003598369831, 0, 0, -897.2003598369831);
    ps.Momenta[2].Set(849.0326646026189, -66.7404825251069, 377.3574531603122,
                      -757.6302040566372);
    ps.Momenta[3].Set(885.7256706487515, 61.57671899811962, -350.9909813646522,
                      810.8782906973681);
    ps.Momenta[4].Set(59.64238442259636, 5.163763526987276, -26.36647179566004,
                      -53.2480866407309);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 3.874984666697015e-07;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, -6.458307777828358e-08);
    ColorCorr.Set(0, 4, 5.812477000045522e-07);
    ColorCorr.Set(1, 0, -6.458307777828358e-08);
    ColorCorr.Set(1, 4, 5.812477000045522e-07);
    ColorCorr.Set(4, 0, 5.812477000045522e-07);
    ColorCorr.Set(4, 1, 5.812477000045522e-07);

    // radiation variables
    double y = 0.9076100860985505;
    double phi = 0.235774136855068;
    int bornpdgs[] = {-1, 2, -13, 14, 0};
    int realpdgs[] = {-1, 0, -13, 14, 0, -2};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 3.008593779392085e-05
    double fks_g = 4.837048521678736e-22;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 0, alphas, ColorCorr, y, phi);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 2, -13, 14, 0
// realpdgs = -1, 0, -13, 14, 0, -2
TEST(QCDCollinearISR2, Wj_born_5_real_14) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.02147544489680442;
    ps.X2 = 0.6477513979637095;
    ps.Jacobian = 143.7247800086063;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(766.634961649253, 0, 0, 766.634961649253);
    ps.Momenta[1].Set(766.634961649253, 0, 0, -766.634961649253);
    ps.Momenta[2].Set(672.8011138185841, 505.2247041034775, 283.2065162351339,
                      -342.3498302629047);
    ps.Momenta[3].Set(674.4560960424808, -365.8220156837729, -373.5076538550368,
                      426.0954245725146);
    ps.Momenta[4].Set(186.0127134374413, -139.4026884197047, 90.30113761990283,
                      -83.74559430960983);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 2.102022389547997e-07;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(5.018066452946037e-08, 0);
    spin[1][2] =
        std::complex<double>(-3.608755062044814e-09, -1.315720691365558e-08);
    spin[1][3] =
        std::complex<double>(-8.742184338105787e-08, -1.418714335957866e-08);
    spin[2][1] =
        std::complex<double>(-3.608755062044814e-09, 1.315720691365558e-08);
    spin[2][2] =
        std::complex<double>(3.709301353657896e-09, 5.88217324482153e-24);
    spin[2][3] =
        std::complex<double>(1.000678658281915e-08, -2.190145082829085e-08);
    spin[3][1] =
        std::complex<double>(-8.742184338105787e-08, 1.418714335957866e-08);
    spin[3][2] =
        std::complex<double>(1.000678658281915e-08, 2.190145082829085e-08);
    spin[3][3] =
        std::complex<double>(1.563122730716814e-07, 5.88217324482153e-24);

    // radiation variables
    double phi = 0.1500719432461769;
    double xi = 0.1059822448291994;
    int bornpdgs[] = {-1, 2, -13, 14, 0};
    int realpdgs[] = {-1, 0, -13, 14, 0, -2};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 0.0002027935111123865
    double fks_g = 4.555649238606312e-14;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, -1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 4.555521399836432e-14
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 2, -13, 14, 0
// realpdgs = -1, 0, -13, 14, 0, -2
TEST(QCDSoftCollinearISR2, Wj_born_5_real_14) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.4542767039484268;
    ps.X2 = 0.7910359190160108;
    ps.Jacobian = 1767.286479684321;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3896.473184471567, 0, 0, 3896.473184471567);
    ps.Momenta[1].Set(3896.473184471567, 0, 0, -3896.473184471567);
    ps.Momenta[2].Set(3649.756706452038, -2946.154170925535, -2065.527210477879,
                      611.9615675883291);
    ps.Momenta[3].Set(3693.166507966026, 2862.343690452035, 2320.275704103435,
                      -250.5755579435073);
    ps.Momenta[4].Set(450.0231545250713, 83.81048047350068, -254.748493625556,
                      -361.3860096448219);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 5.23610905825663e-07;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(4.169851728490149e-08, 1.176434648964306e-23);
    spin[1][2] =
        std::complex<double>(-1.099792796961962e-07, -1.626990974907123e-08);
    spin[1][3] =
        std::complex<double>(8.719714587811618e-08, 1.146899683270301e-08);
    spin[2][1] =
        std::complex<double>(-1.099792796961962e-07, 1.626990974907123e-08);
    spin[2][2] = std::complex<double>(2.964170606184248e-07, 0);
    spin[2][3] =
        std::complex<double>(-2.344562702751008e-07, 3.773220094132365e-09);
    spin[3][1] =
        std::complex<double>(8.719714587811618e-08, -1.146899683270301e-08);
    spin[3][2] =
        std::complex<double>(-2.344562702751008e-07, -3.773220094132365e-09);
    spin[3][3] = std::complex<double>(1.854953279223367e-07, 0);

    // radiation variables
    double phi = 2.199626290277565;
    int bornpdgs[] = {-1, 2, -13, 14, 0};
    int realpdgs[] = {-1, 0, -13, 14, 0, -2};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 0.1223612779288412
    double fks_g = 1.0686056204265e-20;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, -1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, 2, -13, 14, 1
// realpdgs = 2, 2, -13, 14, 1, 2
TEST(QCDCollinearISR1, Wj_born_0_real_18) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.8094952977355341;
    ps.X2 = 0.6647995197594909;
    ps.Jacobian = 4507.359091088889;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(4768.325240474459, 0, 0, 4768.325240474459);
    ps.Momenta[1].Set(4768.325240474459, 0, 0, -4768.325240474459);
    ps.Momenta[2].Set(4730.656105732294, 3648.938102506375, 1917.042295570406,
                      2321.488046872857);
    ps.Momenta[3].Set(3868.095781290271, -2763.690080045053, -1680.628694067406,
                      -2121.242397128024);
    ps.Momenta[4].Set(937.8985939263536, -885.248022461322, -236.4136015030005,
                      -200.2456497448332);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.411138824743814e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(6.922195621900509e-10, 0);
    spin[1][2] =
        std::complex<double>(6.9994613910473e-12, 7.054083815175402e-10);
    spin[2][1] =
        std::complex<double>(6.9994613910473e-12, -7.054083815175402e-10);
    spin[2][2] = std::complex<double>(7.189192625537632e-10, 0);

    // radiation variables
    double phi = 4.334930338417124;
    double xi = 0.1484942332936416;
    int bornpdgs[] = {0, 2, -13, 14, 1};
    int realpdgs[] = {2, 2, -13, 14, 1, 2};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 9.977557759689719e-08
    double fks_g = 4.400210192826432e-17;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, 1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 4.400298615397618e-17
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, 2, -13, 14, 1
// realpdgs = 2, 2, -13, 14, 1, 2
TEST(QCDSoftCollinearISR1, Wj_born_0_real_18) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.9655009280283515;
    ps.X2 = 0.009215925405811731;
    ps.Jacobian = 146.3360177464005;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(613.1393369169376, 0, 0, 613.1393369169376);
    ps.Momenta[1].Set(613.1393369169376, 0, 0, -613.1393369169376);
    ps.Momenta[2].Set(496.9447807788736, -52.92870890786504, 456.0514571785618,
                      -190.1834254661001);
    ps.Momenta[3].Set(492.5284672682515, 113.2308521012882, -477.8452050646891,
                      -37.77598706472706);
    ps.Momenta[4].Set(236.8054257867501, -60.30214319342312, 21.7937478861274,
                      227.9594125308272);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 3.374594444506631e-07;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(8.495000534606217e-08, 8.823259867232295e-24);
    spin[1][2] =
        std::complex<double>(1.391542980959782e-07, 4.568106307084639e-08);
    spin[2][1] =
        std::complex<double>(1.391542980959782e-07, -4.568106307084639e-08);
    spin[2][2] =
        std::complex<double>(2.525094391046009e-07, -8.823259867232295e-24);

    // radiation variables
    double phi = 1.251920264736616;
    int bornpdgs[] = {0, 2, -13, 14, 1};
    int realpdgs[] = {2, 2, -13, 14, 1, 2};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 51.58906473937017
    double fks_g = 1.228012210843386e-19;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, 1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, 2, -13, 14, 1
// realpdgs = 2, 2, -13, 14, 1, 2
TEST(QCDSoftLimit, Wj_born_0_real_18_0) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.9486164207209917;
    ps.X2 = 0.5682542489466087;
    ps.Jacobian = 8820.420095942685;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(4772.325106184967, 0, 0, 4772.325106184967);
    ps.Momenta[1].Set(4772.325106184967, 0, 0, -4772.325106184967);
    ps.Momenta[2].Set(3212.576872870922, 1726.820438602698, -2681.572745191072,
                      384.5893253556411);
    ps.Momenta[3].Set(4498.244303321523, -2365.201074806806, 3691.429437916376,
                      1006.664985475134);
    ps.Momenta[4].Set(1833.82903617749, 638.3806362041084, -1009.856692725304,
                      -1391.254310830775);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 2.954059371683711e-10;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, 4.431089057525566e-10);
    ColorCorr.Set(0, 4, 4.431089057525566e-10);
    ColorCorr.Set(1, 0, 4.431089057525566e-10);
    ColorCorr.Set(1, 4, -4.923432286139518e-11);
    ColorCorr.Set(4, 0, 4.431089057525566e-10);
    ColorCorr.Set(4, 1, -4.923432286139518e-11);

    // radiation variables
    double y = 0.1211562391775018;
    double phi = 2.531432972592501;
    int bornpdgs[] = {0, 2, -13, 14, 1};
    int realpdgs[] = {2, 2, -13, 14, 1, 2};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 1.981037831110539e-07
    double fks_g = 1.574612746553541e-25;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 0, alphas, ColorCorr, y, phi);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, 2, -13, 14, 1
// realpdgs = 2, 2, -13, 14, 1, 2
TEST(QCDCollinearISR2, Wj_born_0_real_18) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.6602005291238875;
    ps.X2 = 0.8113832947023951;
    ps.Jacobian = 5023.078303806879;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(4757.341431985199, 0, 0, 4757.341431985199);
    ps.Momenta[1].Set(4757.341431985199, 0, 0, -4757.341431985199);
    ps.Momenta[2].Set(3901.64343316976, -2926.534905451256, 1471.692484871963,
                      -2119.513188622213);
    ps.Momenta[3].Set(4565.415951838151, 3832.847674225944, -1798.546520526437,
                      1708.078432912629);
    ps.Momenta[4].Set(1047.623478962487, -906.3127687746885, 326.8540356544738,
                      411.4347557095848);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 3.214232024547605e-10;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(5.159040243880552e-11, 0);
    spin[1][2] =
        std::complex<double>(1.178068940281126e-10, 6.505264288395642e-12);
    spin[2][1] =
        std::complex<double>(1.178068940281126e-10, -6.505264288395642e-12);
    spin[2][2] = std::complex<double>(2.69832800015955e-10, 0);

    // radiation variables
    double phi = 1.049563536032786;
    double xi = 0.1012471733008871;
    int bornpdgs[] = {0, 2, -13, 14, 1};
    int realpdgs[] = {2, 2, -13, 14, 1, 2};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 9.61713132967772e-08
    double fks_g = 1.971702360253809e-17;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, -1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, 2, -13, 14, 1
// realpdgs = 2, 2, -13, 14, 1, 2
TEST(QCDSoftCollinearISR2, Wj_born_0_real_18) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.5439424796741115;
    ps.X2 = 0.3574461279428007;
    ps.Jacobian = 3465.61504363866;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2866.125106653002, 0, 0, 2866.125106653002);
    ps.Momenta[1].Set(2866.125106653002, 0, 0, -2866.125106653002);
    ps.Momenta[2].Set(2669.402306851475, 653.5859658649679, 530.9103197472416,
                      2533.114346695429);
    ps.Momenta[3].Set(1863.114441521963, -958.9111846739779, 145.8752530553429,
                      -1590.724731894518);
    ps.Momenta[4].Set(1199.733464932567, 305.3252188090102, -676.7855728025845,
                      -942.3896148009112);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 2.611939787966089e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(1.262300165531638e-09, 0);
    spin[1][2] =
        std::complex<double>(1.914909305043981e-10, 1.291116393800971e-09);
    spin[2][1] =
        std::complex<double>(1.914909305043981e-10, -1.291116393800971e-09);
    spin[2][2] = std::complex<double>(1.349639622434451e-09, 0);

    // radiation variables
    double phi = 1.977350231206304;
    int bornpdgs[] = {0, 2, -13, 14, 1};
    int realpdgs[] = {2, 2, -13, 14, 1, 2};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 5.187874951581478e-05
    double fks_g = 4.283893471820721e-23;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, -1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, 2, -13, 14, 1
// realpdgs = 1, 2, -13, 14, 1, 1
TEST(QCDCollinearISR1, Wj_born_0_real_28) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.8598998470831454;
    ps.X2 = 0.4445033263135088;
    ps.Jacobian = 3547.84561638069;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(4018.600187033855, 0, 0, 4018.600187033855);
    ps.Momenta[1].Set(4018.600187033855, 0, 0, -4018.600187033855);
    ps.Momenta[2].Set(3639.968660267412, 120.2580281412154, -873.3225406150818,
                      3531.602694875251);
    ps.Momenta[3].Set(3521.261155451748, -825.6913954858255, 1258.541849584675,
                      -3183.329429568676);
    ps.Momenta[4].Set(875.9705583485501, 705.4333673446101, -385.2193089695929,
                      -348.2732653065755);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 2.354082778080956e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(1.103959728580227e-09, 0);
    spin[1][2] =
        std::complex<double>(2.347179630333071e-10, 1.151083394132765e-09);
    spin[2][1] =
        std::complex<double>(2.347179630333071e-10, -1.151083394132765e-09);
    spin[2][2] = std::complex<double>(1.250123049500729e-09, 0);

    // radiation variables
    double phi = 2.004608209161613;
    double xi = 0.09843460400744029;
    int bornpdgs[] = {0, 2, -13, 14, 1};
    int realpdgs[] = {1, 2, -13, 14, 1, 1};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 1.603757061510391e-07
    double fks_g = 3.107879516285408e-17;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, 1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 3.10781576863571e-17
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, 2, -13, 14, 1
// realpdgs = 1, 2, -13, 14, 1, 1
TEST(QCDSoftCollinearISR1, Wj_born_0_real_28) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.3005430860587346;
    ps.X2 = 0.8896548448279376;
    ps.Jacobian = 4212.783447890874;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3361.069566670561, 0, 0, 3361.069566670561);
    ps.Momenta[1].Set(3361.069566670561, 0, 0, -3361.069566670561);
    ps.Momenta[2].Set(2146.607527257493, -1320.265908489586, 14.86267215239633,
                      1692.513192839331);
    ps.Momenta[3].Set(3331.901490481967, 2178.935664837342, 405.5746369403184,
                      -2487.833620772605);
    ps.Momenta[4].Set(1243.630115601663, -858.6697563477566, -420.4373090927147,
                      795.3204279332739);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 5.413883363302129e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(9.748577067244423e-10, 1.378634354255046e-25);
    spin[1][2] =
        std::complex<double>(-1.727235109454526e-09, 1.159343455723927e-09);
    spin[2][1] =
        std::complex<double>(-1.727235109454526e-09, -1.159343455723927e-09);
    spin[2][2] =
        std::complex<double>(4.439025656577687e-09, -1.378634354255046e-25);

    // radiation variables
    double phi = 4.190127349394888;
    int bornpdgs[] = {0, 2, -13, 14, 1};
    int realpdgs[] = {1, 2, -13, 14, 1, 1};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 0.0006772085561210857
    double fks_g = 6.626350930033437e-22;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, 1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, 2, -13, 14, 1
// realpdgs = 1, 2, -13, 14, 1, 1
TEST(QCDSoftLimit, Wj_born_0_real_28_0) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.4031367411163167;
    ps.X2 = 0.6604288903321418;
    ps.Jacobian = 2994.313543916942;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3353.919067646776, 0, 0, 3353.919067646776);
    ps.Momenta[1].Set(3353.919067646776, 0, 0, -3353.919067646776);
    ps.Momenta[2].Set(2982.92315636685, 828.0501430945592, -1093.613755236533,
                      2648.805857677777);
    ps.Momenta[3].Set(2839.097404082006, -24.41215144385734, 733.0337036196614,
                      -2742.724868827479);
    ps.Momenta[4].Set(885.8175748446953, -803.637991650702, 360.5800516168717,
                      93.91901114970176);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 4.128256015968351e-09;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, 6.192384023952526e-09);
    ColorCorr.Set(0, 4, 6.192384023952526e-09);
    ColorCorr.Set(1, 0, 6.192384023952526e-09);
    ColorCorr.Set(1, 4, -6.880426693280584e-10);
    ColorCorr.Set(4, 0, 6.192384023952526e-09);
    ColorCorr.Set(4, 1, -6.880426693280584e-10);

    // radiation variables
    double y = 0.7835418735449267;
    double phi = 4.340078501999828;
    int bornpdgs[] = {0, 2, -13, 14, 1};
    int realpdgs[] = {1, 2, -13, 14, 1, 1};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 2.600579775216211e-07
    double fks_g = 4.170263031700361e-24;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 0, alphas, ColorCorr, y, phi);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, 2, -13, 14, 1
// realpdgs = 1, 2, -13, 14, 1, 1
TEST(QCDCollinearISR2, Wj_born_0_real_28) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.7618126084569461;
    ps.X2 = 0.001037398846868598;
    ps.Jacobian = 25.32714475658464;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(182.7301939614796, 0, 0, 182.7301939614796);
    ps.Momenta[1].Set(182.7301939614796, 0, 0, -182.7301939614796);
    ps.Momenta[2].Set(110.7084582469822, -80.34307211949408, -58.42905181217374,
                      48.86102121482828);
    ps.Momenta[3].Set(117.2288739278569, -53.46374486810424, 67.96908223564043,
                      -79.14821998703682);
    ps.Momenta[4].Set(137.5230557481202, 133.8068169875983, -9.540030423466694,
                      30.28719877220855);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.240867464456604e-06;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(3.391936952641257e-07, 0);
    spin[1][2] =
        std::complex<double>(5.49396500520889e-07, 6.328935858839013e-08);
    spin[2][1] =
        std::complex<double>(5.49396500520889e-07, -6.328935858839013e-08);
    spin[2][2] = std::complex<double>(9.016737691924783e-07, 0);

    // radiation variables
    double phi = 0.4616249524625016;
    double xi = 0.5150770595989516;
    int bornpdgs[] = {0, 2, -13, 14, 1};
    int realpdgs[] = {1, 2, -13, 14, 1, 1};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 1.753173613760565e-10
    double fks_g = 9.302492673934484e-19;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, -1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, 2, -13, 14, 1
// realpdgs = 1, 2, -13, 14, 1, 1
TEST(QCDSoftCollinearISR2, Wj_born_0_real_28) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.7698537516323576;
    ps.X2 = 0.5518784335492626;
    ps.Jacobian = 3453.243395105129;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(4236.811901202311, 0, 0, 4236.811901202311);
    ps.Momenta[1].Set(4236.811901202311, 0, 0, -4236.811901202311);
    ps.Momenta[2].Set(3766.914406237276, 1545.643667560278, 855.5970070482044,
                      -3326.948084713511);
    ps.Momenta[3].Set(3898.00909868064, -2014.66501066772, -1434.696023034948,
                      3012.847017303732);
    ps.Momenta[4].Set(808.700297486705, 469.021343107442, 579.0990159867433,
                      314.1010674097784);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 9.949839002862395e-11;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(3.659140507411931e-11, -2.154116178523509e-27);
    spin[1][2] =
        std::complex<double>(-3.586955685632674e-11, -3.186267188094168e-11);
    spin[2][1] =
        std::complex<double>(-3.586955685632674e-11, 3.186267188094168e-11);
    spin[2][2] =
        std::complex<double>(6.290698495450464e-11, 2.154116178523509e-27);

    // radiation variables
    double phi = 4.182971891760776;
    int bornpdgs[] = {0, 2, -13, 14, 1};
    int realpdgs[] = {1, 2, -13, 14, 1, 1};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 6.113813177930482e-12
    double fks_g = 2.455466254835358e-30;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, -1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, 2, -13, 14, 1
// realpdgs = -2, 2, -13, 14, 1, -2
TEST(QCDCollinearISR1, Wj_born_0_real_1) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.823860075853684;
    ps.X2 = 0.7029747531034687;
    ps.Jacobian = 10339.09498540151;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(4946.635948984502, 0, 0, 4946.635948984502);
    ps.Momenta[1].Set(4946.635948984502, 0, 0, -4946.635948984502);
    ps.Momenta[2].Set(3892.454018290578, -1908.17009560393, 1416.637266411637,
                      -3082.729963224661);
    ps.Momenta[3].Set(3926.992561317112, 2320.696562616803, -2759.013044424472,
                      1556.754656840341);
    ps.Momenta[4].Set(2073.825318361314, -412.5264670128726, 1342.375778012835,
                      1525.97530638432);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 3.725032212017281e-10;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(2.226334057502654e-10, 0);
    spin[1][2] =
        std::complex<double>(1.192457465656923e-10, -1.3837080385368e-10);
    spin[2][1] =
        std::complex<double>(1.192457465656923e-10, 1.3837080385368e-10);
    spin[2][2] = std::complex<double>(1.498698154514627e-10, 0);

    // radiation variables
    double phi = 0.2393375073893342;
    double xi = 0.1232513283188984;
    int bornpdgs[] = {0, 2, -13, 14, 1};
    int realpdgs[] = {-2, 2, -13, 14, 1, -2};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 3.149097945644448e-08
    double fks_g = 9.567520050673356e-18;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, 1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 9.567643960935375e-18
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, 2, -13, 14, 1
// realpdgs = -2, 2, -13, 14, 1, -2
TEST(QCDSoftCollinearISR1, Wj_born_0_real_1) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.9224623531296352;
    ps.X2 = 0.2134759176727101;
    ps.Jacobian = 2425.162748755032;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2884.444099503236, 0, 0, 2884.444099503236);
    ps.Momenta[1].Set(2884.444099503236, 0, 0, -2884.444099503236);
    ps.Momenta[2].Set(2099.3664242193, -313.7642493561609, 87.88178840243374,
                      2073.92578706001);
    ps.Momenta[3].Set(2835.306076081356, 834.6555228940779, -286.688420265979,
                      -2694.461069091863);
    ps.Momenta[4].Set(834.2156987058153, -520.8912735379171, 198.8066318635453,
                      620.5352820318528);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.534539726927635e-08;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(3.883376822721056e-09, -2.757268708510092e-25);
    spin[1][2] =
        std::complex<double>(3.349686304908998e-09, 5.769830690965965e-09);
    spin[2][1] =
        std::complex<double>(3.349686304908998e-09, -5.769830690965965e-09);
    spin[2][2] =
        std::complex<double>(1.146202044655529e-08, 2.757268708510092e-25);

    // radiation variables
    double phi = 1.426468042251543;
    int bornpdgs[] = {0, 2, -13, 14, 1};
    int realpdgs[] = {-2, 2, -13, 14, 1, -2};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 0.04698716419978819
    double fks_g = 5.649820688490024e-22;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, 1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, 2, -13, 14, 1
// realpdgs = -2, 2, -13, 14, 1, -2
TEST(QCDSoftLimit, Wj_born_0_real_1_0) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.1227951114665835;
    ps.X2 = 0.1072080701032334;
    ps.Jacobian = 123.3098343004067;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(745.7918525327643, 0, 0, 745.7918525327643);
    ps.Momenta[1].Set(745.7918525327643, 0, 0, -745.7918525327643);
    ps.Momenta[2].Set(589.7996773085964, 160.5584295888518, -222.7269931152095,
                      521.9936173739518);
    ps.Momenta[3].Set(737.7327278459082, -187.1973812494973, 206.5074562207826,
                      -683.0529911460142);
    ps.Momenta[4].Set(164.051299911024, 26.63895166064553, 16.21953689442694,
                      161.0593737720623);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 3.9018543280789e-06;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, 5.85278149211835e-06);
    ColorCorr.Set(0, 4, 5.85278149211835e-06);
    ColorCorr.Set(1, 0, 5.85278149211835e-06);
    ColorCorr.Set(1, 4, -6.503090546798167e-07);
    ColorCorr.Set(4, 0, 5.85278149211835e-06);
    ColorCorr.Set(4, 1, -6.503090546798167e-07);

    // radiation variables
    double y = 0.4318735874343744;
    double phi = 1.406961511628433;
    int bornpdgs[] = {0, 2, -13, 14, 1};
    int realpdgs[] = {-2, 2, -13, 14, 1, -2};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 0.002530841581497845
    double fks_g = 1.91569158994007e-19;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 0, alphas, ColorCorr, y, phi);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, 2, -13, 14, 1
// realpdgs = -2, 2, -13, 14, 1, -2
TEST(QCDCollinearISR2, Wj_born_0_real_1) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.5714548346548192;
    ps.X2 = 0.7263871307465024;
    ps.Jacobian = 12048.21973747107;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(4187.823628409097, 0, 0, 4187.823628409097);
    ps.Momenta[1].Set(4187.823628409097, 0, 0, -4187.823628409097);
    ps.Momenta[2].Set(2709.279609671949, -140.8034524959513, -2349.859619953406,
                      -1341.092896730813);
    ps.Momenta[3].Set(2811.841055140808, 1152.188840645152, -22.45220937811428,
                      2564.840520025015);
    ps.Momenta[4].Set(2854.526592005438, -1011.385388149201, 2372.31182933152,
                      -1223.747623294202);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 3.132001304048327e-10;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(2.551595645863285e-10, 0);
    spin[1][2] =
        std::complex<double>(-6.732602700663522e-11, -1.013746101846784e-10);
    spin[2][1] =
        std::complex<double>(-6.732602700663522e-11, 1.013746101846784e-10);
    spin[2][2] = std::complex<double>(5.804056581850416e-11, 0);

    // radiation variables
    double phi = 0.9085588144752423;
    double xi = 0.09543309053567034;
    int bornpdgs[] = {0, 2, -13, 14, 1};
    int realpdgs[] = {-2, 2, -13, 14, 1, -2};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 1.791522624642989e-16
    double fks_g = 3.263249418755988e-26;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, -1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, 2, -13, 14, 1
// realpdgs = -2, 2, -13, 14, 1, -2
TEST(QCDSoftCollinearISR2, Wj_born_0_real_1) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.2569977667652203;
    ps.X2 = 0.8686069697959127;
    ps.Jacobian = 9121.015946416446;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3071.069792938165, 0, 0, 3071.069792938165);
    ps.Momenta[1].Set(3071.069792938165, 0, 0, -3071.069792938165);
    ps.Momenta[2].Set(688.9926807730375, 439.5751814901412, 33.30540959375253,
                      -529.5047909778707);
    ps.Momenta[3].Set(2506.330186678197, -174.9390780699474, -1424.449732286217,
                      -2054.757962342068);
    ps.Momenta[4].Set(2946.816718425096, -264.6361034201936, 1391.144322692465,
                      2584.262753319938);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.516575543066736e-07;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(7.886626835087488e-08, 0);
    spin[1][2] =
        std::complex<double>(-4.061004288079662e-09, -7.565900697185688e-08);
    spin[2][1] =
        std::complex<double>(-4.061004288079662e-09, 7.565900697185688e-08);
    spin[2][2] = std::complex<double>(7.27912859557987e-08, 0);

    // radiation variables
    double phi = 5.9175192251276;
    int bornpdgs[] = {0, 2, -13, 14, 1};
    int realpdgs[] = {-2, 2, -13, 14, 1, -2};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 1.209798266470968e-07
    double fks_g = 4.177224333103751e-27;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, -1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, 2, -13, 14, 1
// realpdgs = -1, 2, -13, 14, 1, -1
TEST(QCDCollinearISR1, Wj_born_0_real_15) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.4429642537681282;
    ps.X2 = 0.4289578149620894;
    ps.Jacobian = 2303.152658964316;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2833.381078766779, 0, 0, 2833.381078766779);
    ps.Momenta[1].Set(2833.381078766779, 0, 0, -2833.381078766779);
    ps.Momenta[2].Set(2666.084394962163, -954.0080078166772, -1469.439664100395,
                      2009.632253834252);
    ps.Momenta[3].Set(2194.153694308497, 1036.513181408429, 693.6228381653793,
                      -1805.225253919744);
    ps.Momenta[4].Set(806.524068262899, -82.50517359175154, 775.816825935016,
                      -204.4069999145077);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 3.758215604428651e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(1.99450685914538e-09, 0);
    spin[1][2] =
        std::complex<double>(2.414024950105822e-10, 1.859960759103493e-09);
    spin[2][1] =
        std::complex<double>(2.414024950105822e-10, -1.859960759103493e-09);
    spin[2][2] = std::complex<double>(1.76370874528327e-09, 0);

    // radiation variables
    double phi = 0.7883091565661946;
    double xi = 0.4295448482681829;
    int bornpdgs[] = {0, 2, -13, 14, 1};
    int realpdgs[] = {-1, 2, -13, 14, 1, -1};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 4.889329223309314e-07
    double fks_g = 1.804248306546225e-15;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, 1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 1.804441735876655e-15
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, 2, -13, 14, 1
// realpdgs = -1, 2, -13, 14, 1, -1
TEST(QCDSoftCollinearISR1, Wj_born_0_real_15) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.1873690913883728;
    ps.X2 = 0.3734114257950623;
    ps.Jacobian = 201.1909543924212;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1719.317696538988, 0, 0, 1719.317696538988);
    ps.Momenta[1].Set(1719.317696538988, 0, 0, -1719.317696538988);
    ps.Momenta[2].Set(1646.735890691937, -1468.10192746194, -729.1715456245423,
                      -157.2408386107687);
    ps.Momenta[3].Set(1675.794237419175, 1457.37541542099, 786.0199121170634,
                      257.9068096217977);
    ps.Momenta[4].Set(116.1052649668652, 10.72651204095023, -56.84836649252112,
                      -100.665971011029);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 6.137131060724224e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(2.940666177371613e-09, 3.446585885637615e-26);
    spin[1][2] =
        std::complex<double>(-5.785172537630251e-10, 3.010822803940744e-09);
    spin[2][1] =
        std::complex<double>(-5.785172537630251e-10, -3.010822803940744e-09);
    spin[2][2] =
        std::complex<double>(3.196464883352611e-09, -3.446585885637615e-26);

    // radiation variables
    double phi = 0.5083474928913504;
    int bornpdgs[] = {0, 2, -13, 14, 1};
    int realpdgs[] = {-1, 2, -13, 14, 1, -1};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 0.005111664069897415
    double fks_g = 6.751169547982389e-21;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, 1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, 2, -13, 14, 1
// realpdgs = -1, 2, -13, 14, 1, -1
TEST(QCDSoftLimit, Wj_born_0_real_15_0) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.3629356357831206;
    ps.X2 = 0.326693044531595;
    ps.Jacobian = 3956.469909010913;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2238.195957802099, 0, 0, 2238.195957802099);
    ps.Momenta[1].Set(2238.195957802099, 0, 0, -2238.195957802099);
    ps.Momenta[2].Set(588.7843625696711, 466.1809125861379, 137.9632047733906,
                      -332.1272895676972);
    ps.Momenta[3].Set(2133.689265428791, -2082.686517231115, -336.8620581879544,
                      -318.7015941333757);
    ps.Momenta[4].Set(1753.918287605736, 1616.505604644977, 198.8988534145637,
                      650.8288837010729);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.153749130344337e-08;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, 1.730623695516506e-08);
    ColorCorr.Set(0, 4, 1.730623695516506e-08);
    ColorCorr.Set(1, 0, 1.730623695516506e-08);
    ColorCorr.Set(1, 4, -1.922915217240562e-09);
    ColorCorr.Set(4, 0, 1.730623695516506e-08);
    ColorCorr.Set(4, 1, -1.922915217240562e-09);

    // radiation variables
    double y = 0.1809662034160695;
    double phi = 1.132449041261327;
    int bornpdgs[] = {0, 2, -13, 14, 1};
    int realpdgs[] = {-1, 2, -13, 14, 1, -1};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 2.293364014761761e-06
    double fks_g = 1.526391814118044e-22;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 0, alphas, ColorCorr, y, phi);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, 2, -13, 14, 1
// realpdgs = -1, 2, -13, 14, 1, -1
TEST(QCDCollinearISR2, Wj_born_0_real_15) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.8881620883886221;
    ps.X2 = 0.4411674860080339;
    ps.Jacobian = 49.77445820502539;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(4068.75201485809, 0, 0, 4068.75201485809);
    ps.Momenta[1].Set(4068.75201485809, 0, 0, -4068.75201485809);
    ps.Momenta[2].Set(4064.981058186872, 792.0150649759867, 3979.720830099901,
                      242.0852220564084);
    ps.Momenta[3].Set(4060.385033471989, -794.6242255083807, -3975.156726249813,
                      -231.1448940153873);
    ps.Momenta[4].Set(12.13793805732038, 2.609160532394131, -4.564103850088213,
                      -10.94032804102112);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 3.366180034674636e-08;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(1.687295772160152e-08, 0);
    spin[1][2] =
        std::complex<double>(-4.129769925519987e-11, 1.683079695999944e-08);
    spin[2][1] =
        std::complex<double>(-4.129769925519987e-11, -1.683079695999944e-08);
    spin[2][2] = std::complex<double>(1.678884262514483e-08, 0);

    // radiation variables
    double phi = 2.252010730932418;
    double xi = 0.269960099641042;
    int bornpdgs[] = {0, 2, -13, 14, 1};
    int realpdgs[] = {-1, 2, -13, 14, 1, -1};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 2.764214124692448e-11
    double fks_g = 4.029033113828829e-20;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, -1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, 2, -13, 14, 1
// realpdgs = -1, 2, -13, 14, 1, -1
TEST(QCDSoftCollinearISR2, Wj_born_0_real_15) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.699641444785545;
    ps.X2 = 0.6185783285528572;
    ps.Jacobian = 7394.517905451413;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(4276.106084973673, 0, 0, 4276.106084973673);
    ps.Momenta[1].Set(4276.106084973673, 0, 0, -4276.106084973673);
    ps.Momenta[2].Set(3464.457933239186, 905.5592204834654, 107.5501091115976,
                      -3342.28428524735);
    ps.Momenta[3].Set(3371.976649995504, -1083.187414768729, 1514.320459470254,
                      2811.363565715279);
    ps.Momenta[4].Set(1715.777586712658, 177.6281942852637, -1621.870568581852,
                      530.920719532071);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.343679178236145e-10;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(5.721161847711828e-11, 0);
    spin[1][2] =
        std::complex<double>(1.293887568340965e-12, -6.64271225952063e-11);
    spin[2][1] =
        std::complex<double>(1.293887568340965e-12, 6.64271225952063e-11);
    spin[2][2] = std::complex<double>(7.715629934649626e-11, 0);

    // radiation variables
    double phi = 1.823982906608731;
    int bornpdgs[] = {0, 2, -13, 14, 1};
    int realpdgs[] = {-1, 2, -13, 14, 1, -1};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 5.356675585021842e-11
    double fks_g = 1.558605502113516e-29;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, -1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, 2, -13, 14, 1
// realpdgs = 3, 2, -13, 14, 1, 3
TEST(QCDCollinearISR1, Wj_born_0_real_3) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.6229221781105672;
    ps.X2 = 0.8263115314978151;
    ps.Jacobian = 16953.12728004708;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(4663.394542893247, 0, 0, 4663.394542893247);
    ps.Momenta[1].Set(4663.394542893247, 0, 0, -4663.394542893247);
    ps.Momenta[2].Set(1937.835422399733, 357.8802182324973, -296.5032954133751,
                      1881.27979564817);
    ps.Momenta[3].Set(3781.944390075022, -1319.468896794609, -3167.73169909275,
                      -1589.83681004514);
    ps.Momenta[4].Set(3607.00927331174, 961.5886785621112, 3464.234994506125,
                      -291.4429856030299);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.391715073827747e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(1.153979261250673e-09, 0);
    spin[1][2] =
        std::complex<double>(-4.667633345408758e-10, -2.376429820104153e-10);
    spin[2][1] =
        std::complex<double>(-4.667633345408758e-10, 2.376429820104153e-10);
    spin[2][2] = std::complex<double>(2.377358125770739e-10, 0);

    // radiation variables
    double phi = 6.069164112973635;
    double xi = 0.1069648499858544;
    int bornpdgs[] = {0, 2, -13, 14, 1};
    int realpdgs[] = {3, 2, -13, 14, 1, 3};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 1.588320804657634e-07
    double fks_g = 3.634547866526958e-17;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, 1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 3.634540945198449e-17
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, 2, -13, 14, 1
// realpdgs = 3, 2, -13, 14, 1, 3
TEST(QCDSoftCollinearISR1, Wj_born_0_real_3) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.06084664907834636;
    ps.X2 = 0.1819200683370923;
    ps.Jacobian = 144.3012529534145;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(683.8675471851764, 0, 0, 683.8675471851764);
    ps.Momenta[1].Set(683.8675471851764, 0, 0, -683.8675471851764);
    ps.Momenta[2].Set(516.1102809838302, -144.0219516674478, 455.7540366231662,
                      -194.7196900079629);
    ps.Momenta[3].Set(642.262888011182, 19.70730065700174, -546.6805737935958,
                      336.5317070559377);
    ps.Momenta[4].Set(209.3619253753404, 124.314651010446, 90.92653717042971,
                      -141.8120170479748);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 7.866544997168871e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(7.838229567341976e-09, 0);
    spin[1][2] =
        std::complex<double>(4.397295997938564e-10, -1.690583282370337e-10);
    spin[2][1] =
        std::complex<double>(4.397295997938564e-10, 1.690583282370337e-10);
    spin[2][2] = std::complex<double>(2.831542982689509e-11, 0);

    // radiation variables
    double phi = 2.709265898123502;
    int bornpdgs[] = {0, 2, -13, 14, 1};
    int realpdgs[] = {3, 2, -13, 14, 1, 3};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 0.03541282674101068
    double fks_g = 6.246886687693694e-20;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, 1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, 2, -13, 14, 1
// realpdgs = 3, 2, -13, 14, 1, 3
TEST(QCDSoftLimit, Wj_born_0_real_3_0) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.4585802405421298;
    ps.X2 = 0.5045498345172703;
    ps.Jacobian = 3042.245580895471;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3126.605298756674, 0, 0, 3126.605298756674);
    ps.Momenta[1].Set(3126.605298756674, 0, 0, -3126.605298756674);
    ps.Momenta[2].Set(2286.958799951692, 259.4685519572212, -1595.043700877153,
                      1618.237379222782);
    ps.Momenta[3].Set(3000.821762647558, -780.8678601530929, 1577.187901911888,
                      -2430.566797728855);
    ps.Momenta[4].Set(965.4300349140987, 521.3993081958718, 17.85579896526537,
                      812.3294185060738);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.571607435270763e-08;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, 2.357411152906145e-08);
    ColorCorr.Set(0, 4, 2.357411152906145e-08);
    ColorCorr.Set(1, 0, 2.357411152906145e-08);
    ColorCorr.Set(1, 4, -2.619345725451272e-09);
    ColorCorr.Set(4, 0, 2.357411152906145e-08);
    ColorCorr.Set(4, 1, -2.619345725451272e-09);

    // radiation variables
    double y = 0.7868489204702982;
    double phi = 5.823706889697633;
    int bornpdgs[] = {0, 2, -13, 14, 1};
    int realpdgs[] = {3, 2, -13, 14, 1, 3};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 2.54426863831299e-06
    double fks_g = 3.328852415954425e-23;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 0, alphas, ColorCorr, y, phi);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, 2, -13, 14, 1
// realpdgs = 3, 2, -13, 14, 1, 3
TEST(QCDCollinearISR2, Wj_born_0_real_3) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.5273237537487496;
    ps.X2 = 0.9264903289922728;
    ps.Jacobian = 8948.228667886624;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(4543.311031567291, 0, 0, 4543.311031567291);
    ps.Momenta[1].Set(4543.311031567291, 0, 0, -4543.311031567291);
    ps.Momenta[2].Set(3314.949853962658, 702.519773649945, 1647.029829982359,
                      -2789.74035370071);
    ps.Momenta[3].Set(3817.493844029289, -2577.133409680894, -1705.820852244841,
                      2240.941288377438);
    ps.Momenta[4].Set(1954.178365142634, 1874.613636030949, 58.7910222624819,
                      548.7990653232715);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 2.363331209594398e-10;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(1.264686076024208e-10, 0);
    spin[1][2] =
        std::complex<double>(-4.759755129622084e-11, -1.078373086634833e-10);
    spin[2][1] =
        std::complex<double>(-4.759755129622084e-11, 1.078373086634833e-10);
    spin[2][2] = std::complex<double>(1.09864513357019e-10, 0);

    // radiation variables
    double phi = 0.7215808522119772;
    double xi = 0.04311815940906569;
    int bornpdgs[] = {0, 2, -13, 14, 1};
    int realpdgs[] = {3, 2, -13, 14, 1, 3};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 2.596739531754066e-16
    double fks_g = 9.655589916734166e-27;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, -1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, 2, -13, 14, 1
// realpdgs = 3, 2, -13, 14, 1, 3
TEST(QCDSoftCollinearISR2, Wj_born_0_real_3) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.5607225376678819;
    ps.X2 = 0.3920563058067259;
    ps.Jacobian = 275.5625224184842;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3047.625400718098, 0, 0, 3047.625400718098);
    ps.Momenta[1].Set(3047.625400718098, 0, 0, -3047.625400718098);
    ps.Momenta[2].Set(3035.291625144264, -1854.328488534021, 77.00355833410813,
                      2401.776750300895);
    ps.Momenta[3].Set(2970.245601275101, 1792.403930491315, -21.34319027414367,
                      -2368.373186405025);
    ps.Momenta[4].Set(89.71357501683114, 61.92455804270531, -55.66036805996446,
                      -33.40356389587043);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 2.458354684719711e-08;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(1.237647798643299e-08, 0);
    spin[1][2] =
        std::complex<double>(7.25559960070423e-11, 1.229126741477713e-08);
    spin[2][1] =
        std::complex<double>(7.25559960070423e-11, -1.229126741477713e-08);
    spin[2][2] = std::complex<double>(1.220706886076412e-08, 0);

    // radiation variables
    double phi = 1.972445891022457;
    int bornpdgs[] = {0, 2, -13, 14, 1};
    int realpdgs[] = {3, 2, -13, 14, 1, 3};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 4.303608855635093e-09
    double fks_g = 3.181189861250256e-27;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, -1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, 2, -13, 14, 1
// realpdgs = 4, 2, -13, 14, 1, 4
TEST(QCDCollinearISR1, Wj_born_0_real_3_n1) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.03958128950842266;
    ps.X2 = 0.8658374536588269;
    ps.Jacobian = 368.829822494331;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1203.307185797275, 0, 0, 1203.307185797275);
    ps.Momenta[1].Set(1203.307185797275, 0, 0, -1203.307185797275);
    ps.Momenta[2].Set(965.2274891798608, 614.0791488481215, 689.4996308808392,
                      -281.3559379746787);
    ps.Momenta[3].Set(1137.263987534699, -627.8716540278587, -941.6156064694472,
                      111.8329694894093);
    ps.Momenta[4].Set(304.1228948799914, 13.79250517973721, 252.115975588608,
                      169.5229684852694);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.117491454918094e-08;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(9.700210468576565e-09, 0);
    spin[1][2] =
        std::complex<double>(-3.247465194822241e-09, 1.938790800769725e-09);
    spin[2][1] =
        std::complex<double>(-3.247465194822241e-09, -1.938790800769725e-09);
    spin[2][2] = std::complex<double>(1.474704080604379e-09, 0);

    // radiation variables
    double phi = 3.463534523321626;
    double xi = 0.396431027706174;
    int bornpdgs[] = {0, 2, -13, 14, 1};
    int realpdgs[] = {4, 2, -13, 14, 1, 4};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 8.592326284873486e-06
    double fks_g = 2.700698061216388e-14;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, 1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 2.700525592753531e-14
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, 2, -13, 14, 1
// realpdgs = 4, 2, -13, 14, 1, 4
TEST(QCDSoftCollinearISR1, Wj_born_0_real_3_n1) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.8300181466990288;
    ps.X2 = 0.3824230949274821;
    ps.Jacobian = 10296.33420927389;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3662.091626981021, 0, 0, 3662.091626981021);
    ps.Momenta[1].Set(3662.091626981021, 0, 0, -3662.091626981021);
    ps.Momenta[2].Set(3247.813792086601, -794.3547424396204, -424.3554864495637,
                      -3120.451472519928);
    ps.Momenta[3].Set(1286.698248334398, -580.0748999193379, -652.8388634391609,
                      944.9375170527644);
    ps.Momenta[4].Set(2789.671213541042, 1374.429642358959, 1077.194349888725,
                      2175.513955467163);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 2.752572600271915e-10;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(1.074308481980039e-10, 0);
    spin[1][2] =
        std::complex<double>(-2.82251300691993e-11, -1.312748102471069e-10);
    spin[2][1] =
        std::complex<double>(-2.82251300691993e-11, 1.312748102471069e-10);
    spin[2][2] = std::complex<double>(1.678264118291877e-10, 0);

    // radiation variables
    double phi = 5.030829519102077;
    int bornpdgs[] = {0, 2, -13, 14, 1};
    int realpdgs[] = {4, 2, -13, 14, 1, 4};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 0.0002386754446262016
    double fks_g = 1.379250138455744e-23;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, 1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, 2, -13, 14, 1
// realpdgs = 4, 2, -13, 14, 1, 4
TEST(QCDSoftLimit, Wj_born_0_real_3_0_n1) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.5241707423190469;
    ps.X2 = 0.6773339181790305;
    ps.Jacobian = 5905.282430867813;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3873.032637177594, 0, 0, 3873.032637177594);
    ps.Momenta[1].Set(3873.032637177594, 0, 0, -3873.032637177594);
    ps.Momenta[2].Set(3422.911029894348, 399.7069789190325, -2543.000569562907,
                      -2256.036868665612);
    ps.Momenta[3].Set(2810.327780409289, -1271.89879583118, 1312.886770142803,
                      2134.606290463745);
    ps.Momenta[4].Set(1512.82646405155, 872.1918169121478, 1230.113799420103,
                      121.4305782018662);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.240207185995922e-10;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, 1.860310778993883e-10);
    ColorCorr.Set(0, 4, 1.860310778993883e-10);
    ColorCorr.Set(1, 0, 1.860310778993883e-10);
    ColorCorr.Set(1, 4, -2.06701197665987e-11);
    ColorCorr.Set(4, 0, 1.860310778993883e-10);
    ColorCorr.Set(4, 1, -2.06701197665987e-11);

    // radiation variables
    double y = -0.1230508815466678;
    double phi = 5.964702871580374;
    int bornpdgs[] = {0, 2, -13, 14, 1};
    int realpdgs[] = {4, 2, -13, 14, 1, 4};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 2.90075202400488e-09
    double fks_g = 7.192989039978742e-26;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 0, alphas, ColorCorr, y, phi);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, 2, -13, 14, 1
// realpdgs = 4, 2, -13, 14, 1, 4
TEST(QCDCollinearISR2, Wj_born_0_real_3_n1) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.8298487745566341;
    ps.X2 = 0.08572268811672501;
    ps.Jacobian = 1907.827598496038;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1733.647213162661, 0, 0, 1733.647213162661);
    ps.Momenta[1].Set(1733.647213162661, 0, 0, -1733.647213162661);
    ps.Momenta[2].Set(946.3101717506145, -627.0306546226558, 444.6506122906768,
                      -551.9251147680156);
    ps.Momenta[3].Set(1429.096494653996, 161.3307287287803, 26.8324747316465,
                      1419.707436516396);
    ps.Momenta[4].Set(1091.887759920712, 465.6999258938756, -471.4830870223234,
                      -867.78232174838);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.965925028855484e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(7.983417251607995e-10, 0);
    spin[1][2] =
        std::complex<double>(4.249403768162993e-10, -8.669233790201887e-10);
    spin[2][1] =
        std::complex<double>(4.249403768162993e-10, 8.669233790201887e-10);
    spin[2][2] = std::complex<double>(1.167583303694685e-09, 0);

    // radiation variables
    double phi = 1.922852704216628;
    double xi = 0.8689583317451095;
    int bornpdgs[] = {0, 2, -13, 14, 1};
    int realpdgs[] = {4, 2, -13, 14, 1, 4};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 6.804581888943227e-16
    double fks_g = 1.02761241780497e-23;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, -1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, 2, -13, 14, 1
// realpdgs = 4, 2, -13, 14, 1, 4
TEST(QCDSoftCollinearISR2, Wj_born_0_real_3_n1) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.6499720950832035;
    ps.X2 = 0.25416505664748;
    ps.Jacobian = 4044.381462722153;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2641.913740447988, 0, 0, 2641.913740447988);
    ps.Momenta[1].Set(2641.913740447988, 0, 0, -2641.913740447988);
    ps.Momenta[2].Set(1611.07580467802, -347.4455536008892, 1347.372468379877,
                      812.0555813205265);
    ps.Momenta[3].Set(2153.838064904369, 365.0334627589693, -2054.969384014446,
                      531.8550664048246);
    ps.Momenta[4].Set(1518.913611313587, -17.58790915808008, 707.5969156345694,
                      -1343.910647725351);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.269720989269733e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(1.219544136159692e-09, -1.723292942818808e-26);
    spin[1][2] =
        std::complex<double>(2.25511420528384e-10, 1.016734291374463e-10);
    spin[2][1] =
        std::complex<double>(2.25511420528384e-10, -1.016734291374463e-10);
    spin[2][2] =
        std::complex<double>(5.017685311004127e-11, 1.723292942818808e-26);

    // radiation variables
    double phi = 1.438078335223723;
    int bornpdgs[] = {0, 2, -13, 14, 1};
    int realpdgs[] = {4, 2, -13, 14, 1, 4};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 2.411172649674404e-10
    double fks_g = 2.682525216458289e-28;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, -1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, 2, -13, 14, 1
// realpdgs = 5, 2, -13, 14, 1, 5
TEST(QCDCollinearISR1, Wj_born_0_real_3_n2) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.1451157861428776;
    ps.X2 = 0.4325778983934541;
    ps.Jacobian = 1045.874243315898;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1628.556571252944, 0, 0, 1628.556571252944);
    ps.Momenta[1].Set(1628.556571252944, 0, 0, -1628.556571252944);
    ps.Momenta[2].Set(1268.924491833581, 181.0432117727815, -682.0298101352927,
                      -1054.622235462975);
    ps.Momenta[3].Set(1350.987993678682, 48.40678677900793, 1156.554914080104,
                      696.5673497766353);
    ps.Momenta[4].Set(637.2006569936258, -229.4499985517893, -474.5251039448115,
                      358.0548856863398);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 2.162195102225687e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(1.6727643699298e-09, 0);
    spin[1][2] =
        std::complex<double>(5.53355118249776e-12, -9.048047691875414e-10);
    spin[2][1] =
        std::complex<double>(5.53355118249776e-12, 9.048047691875414e-10);
    spin[2][2] = std::complex<double>(4.894307322958867e-10, 0);

    // radiation variables
    double phi = 0.2040659019869994;
    double xi = 0.5093039956965585;
    int bornpdgs[] = {0, 2, -13, 14, 1};
    int realpdgs[] = {5, 2, -13, 14, 1, 5};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 1.143692103491527e-06
    double fks_g = 5.93325870146484e-15;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, 1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 5.932764688180686e-15
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, 2, -13, 14, 1
// realpdgs = 5, 2, -13, 14, 1, 5
TEST(QCDSoftCollinearISR1, Wj_born_0_real_3_n2) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.7092473433543098;
    ps.X2 = 0.06643094673729877;
    ps.Jacobian = 1813.024444200923;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1410.90390803201, 0, 0, 1410.90390803201);
    ps.Momenta[1].Set(1410.90390803201, 0, 0, -1410.90390803201);
    ps.Momenta[2].Set(227.54062467376, -36.36726111141945, 99.80444391708079,
                      201.2243304637413);
    ps.Momenta[3].Set(1319.280006791373, -1175.083252434572, -60.15156461405093,
                      -596.708367162481);
    ps.Momenta[4].Set(1274.987184598887, 1211.450513545991, -39.65287930302987,
                      395.4840366987397);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 6.952598014249023e-08;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(2.307127059939242e-08, 0);
    spin[1][2] =
        std::complex<double>(-3.468395679237131e-09, -3.255363890408008e-08);
    spin[2][1] =
        std::complex<double>(-3.468395679237131e-09, 3.255363890408008e-08);
    spin[2][2] = std::complex<double>(4.64547095430978e-08, 0);

    // radiation variables
    double phi = 5.069320585450784;
    int bornpdgs[] = {0, 2, -13, 14, 1};
    int realpdgs[] = {5, 2, -13, 14, 1, 5};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 0.2375686474018896
    double fks_g = 4.016674674084314e-20;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, 1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, 2, -13, 14, 1
// realpdgs = 5, 2, -13, 14, 1, 5
TEST(QCDSoftLimit, Wj_born_0_real_3_0_n2) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.1193777762495891;
    ps.X2 = 0.9585142249204281;
    ps.Jacobian = 2697.775817411306;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2198.742546207221, 0, 0, 2198.742546207221);
    ps.Momenta[1].Set(2198.742546207221, 0, 0, -2198.742546207221);
    ps.Momenta[2].Set(1434.253366441934, 636.2169085865851, -784.1339181579104,
                      -1018.550324123886);
    ps.Momenta[3].Set(1745.837959734726, -538.1956196925976, 1652.036455438887,
                      170.5028049485637);
    ps.Momenta[4].Set(1217.393766237783, -98.02128889398749, -867.9025372809768,
                      848.0475191753219);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 4.034568534644386e-09;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, 6.051852801966579e-09);
    ColorCorr.Set(0, 4, 6.051852801966579e-09);
    ColorCorr.Set(1, 0, 6.051852801966579e-09);
    ColorCorr.Set(1, 4, -6.724280891073977e-10);
    ColorCorr.Set(4, 0, 6.051852801966579e-09);
    ColorCorr.Set(4, 1, -6.724280891073977e-10);

    // radiation variables
    double y = -0.3508231545663492;
    double phi = 2.358287182525622;
    int bornpdgs[] = {0, 2, -13, 14, 1};
    int realpdgs[] = {5, 2, -13, 14, 1, 5};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 2.008635656299015e-06
    double fks_g = 6.513599740812322e-25;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 0, alphas, ColorCorr, y, phi);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, 2, -13, 14, 1
// realpdgs = 5, 2, -13, 14, 1, 5
TEST(QCDCollinearISR2, Wj_born_0_real_3_n2) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.5044398540994202;
    ps.X2 = 0.02237789042524696;
    ps.Jacobian = 29.91923987408232;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(690.6016693827169, 0, 0, 690.6016693827169);
    ps.Momenta[1].Set(690.6016693827169, 0, 0, -690.6016693827169);
    ps.Momenta[2].Set(658.0827219955127, -261.1163408913529, 331.102751900937,
                      -505.2346912001991);
    ps.Momenta[3].Set(680.135062776494, 265.677808085201, -311.4001659661417,
                      543.1656676791049);
    ps.Momenta[4].Set(42.98555399342727, -4.561467193848142, -19.70258593479522,
                      -37.93097647890584);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 3.126819944281042e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(1.659782032656309e-09, 0);
    spin[1][2] =
        std::complex<double>(1.122390326397582e-09, 1.08406785864608e-09);
    spin[2][1] =
        std::complex<double>(1.122390326397582e-09, -1.08406785864608e-09);
    spin[2][2] = std::complex<double>(1.467037911624733e-09, 0);

    // radiation variables
    double phi = 1.002851281689447;
    double xi = 0.8201955520812064;
    int bornpdgs[] = {0, 2, -13, 14, 1};
    int realpdgs[] = {5, 2, -13, 14, 1, 5};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 1.220042194595639e-15
    double fks_g = 1.641495384010762e-23;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, -1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, 2, -13, 14, 1
// realpdgs = 5, 2, -13, 14, 1, 5
TEST(QCDSoftCollinearISR2, Wj_born_0_real_3_n2) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.799906859924004;
    ps.X2 = 0.5329390354337775;
    ps.Jacobian = 16675.9313918488;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(4243.965385451575, 0, 0, 4243.965385451575);
    ps.Momenta[1].Set(4243.965385451575, 0, 0, -4243.965385451575);
    ps.Momenta[2].Set(1922.94979269235, 1494.126899000553, 67.07240474388823,
                      1208.644698586154);
    ps.Momenta[3].Set(2666.298528746143, -224.1823246339779, -1520.476649963524,
                      2178.770498843614);
    ps.Momenta[4].Set(3898.682449464655, -1269.944574366575, 1453.404245219636,
                      -3387.415197429769);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.949995969439254e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(6.989159130477895e-10, -6.89317177127523e-26);
    spin[1][2] =
        std::complex<double>(-8.60870330138215e-10, -3.651055115941704e-10);
    spin[2][1] =
        std::complex<double>(-8.60870330138215e-10, 3.651055115941704e-10);
    spin[2][2] =
        std::complex<double>(1.251080056391464e-09, 6.89317177127523e-26);

    // radiation variables
    double phi = 4.35698593352235;
    int bornpdgs[] = {0, 2, -13, 14, 1};
    int realpdgs[] = {5, 2, -13, 14, 1, 5};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 2.291484968875022e-10
    double fks_g = 9.997565726869643e-29;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, -1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, 2, -13, 14, 1
// realpdgs = -5, 2, -13, 14, 1, -5
TEST(QCDCollinearISR1, Wj_born_0_real_19) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.5367555725350321;
    ps.X2 = 0.3851596232486685;
    ps.Jacobian = 677.0147812268387;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2955.439096899247, 0, 0, 2955.439096899247);
    ps.Momenta[1].Set(2955.439096899247, 0, 0, -2955.439096899247);
    ps.Momenta[2].Set(2843.23038421345, -140.8731483462938, 405.4357262575247,
                      2810.646837591589);
    ps.Momenta[3].Set(2840.360214772935, -14.96078545916347, -568.537253522225,
                      -2782.838068576898);
    ps.Momenta[4].Set(227.2875948121098, 155.8339338054573, 163.1015272647003,
                      -27.80876901469024);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.481531786707404e-08;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(7.519886340096118e-09, -3.446585885637615e-26);
    spin[1][2] =
        std::complex<double>(-5.763087196125443e-10, 7.384354009967723e-09);
    spin[2][1] =
        std::complex<double>(-5.763087196125443e-10, -7.384354009967723e-09);
    spin[2][2] =
        std::complex<double>(7.295431526977921e-09, 3.446585885637615e-26);

    // radiation variables
    double phi = 3.224797783014643;
    double xi = 0.02209140955048874;
    int bornpdgs[] = {0, 2, -13, 14, 1};
    int realpdgs[] = {-5, 2, -13, 14, 1, -5};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 1.553205708273418e-05
    double fks_g = 1.51602313060463e-16;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, 1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 1.516016765408397e-16
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, 2, -13, 14, 1
// realpdgs = -5, 2, -13, 14, 1, -5
TEST(QCDSoftCollinearISR1, Wj_born_0_real_19) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.2265758612386399;
    ps.X2 = 0.04383045938981667;
    ps.Jacobian = 358.5183285247571;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(647.7511424767747, 0, 0, 647.7511424767747);
    ps.Momenta[1].Set(647.7511424767747, 0, 0, -647.7511424767747);
    ps.Momenta[2].Set(537.7231158226502, 167.1404386881056, -508.5727875498272,
                      50.63539087402395);
    ps.Momenta[3].Set(208.6142687122693, -161.9788153510034, -9.619741371924761,
                      131.1115443573819);
    ps.Momenta[4].Set(549.1649004186299, -5.161623337102142, 518.192528921752,
                      -181.7469352314058);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 3.862780700979991e-08;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(1.270864900921578e-08, 0);
    spin[1][2] =
        std::complex<double>(1.081721443039253e-08, 1.457344686787373e-08);
    spin[2][1] =
        std::complex<double>(1.081721443039253e-08, -1.457344686787373e-08);
    spin[2][2] = std::complex<double>(2.591915800058414e-08, 0);

    // radiation variables
    double phi = 0.412274377053483;
    int bornpdgs[] = {0, 2, -13, 14, 1};
    int realpdgs[] = {-5, 2, -13, 14, 1, -5};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 0.2354354172527713
    double fks_g = 2.816678542277666e-19;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, 1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, 2, -13, 14, 1
// realpdgs = -5, 2, -13, 14, 1, -5
TEST(QCDSoftLimit, Wj_born_0_real_19_0) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.4991851594486185;
    ps.X2 = 0.2892406960182257;
    ps.Jacobian = 87.52235791921585;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2469.868824471755, 0, 0, 2469.868824471755);
    ps.Momenta[1].Set(2469.868824471755, 0, 0, -2469.868824471755);
    ps.Momenta[2].Set(2461.785540628995, 2086.455504468357, 923.2367669287961,
                      -924.5135737859945);
    ps.Momenta[3].Set(2442.792444209704, -2085.836057393187, -897.5315075807456,
                      900.5332087526234);
    ps.Momenta[4].Set(35.15966410481195, -0.6194470751699876, -25.7052593480506,
                      23.98036503337099);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 4.004522054271675e-08;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, 6.006783081407513e-08);
    ColorCorr.Set(0, 4, 6.006783081407513e-08);
    ColorCorr.Set(1, 0, 6.006783081407513e-08);
    ColorCorr.Set(1, 4, -6.674203423786125e-09);
    ColorCorr.Set(4, 0, 6.006783081407513e-08);
    ColorCorr.Set(4, 1, -6.674203423786125e-09);

    // radiation variables
    double y = -0.7307661741954163;
    double phi = 2.729705340271881;
    int bornpdgs[] = {0, 2, -13, 14, 1};
    int realpdgs[] = {-5, 2, -13, 14, 1, -5};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 9.696164053278044e-07
    double fks_g = 2.7024224823942e-23;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 0, alphas, ColorCorr, y, phi);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, 2, -13, 14, 1
// realpdgs = -5, 2, -13, 14, 1, -5
TEST(QCDCollinearISR2, Wj_born_0_real_19) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.4199595581953979;
    ps.X2 = 0.9227273328477352;
    ps.Jacobian = 10802.78914878613;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(4046.259987733932, 0, 0, 4046.259987733932);
    ps.Momenta[1].Set(4046.259987733932, 0, 0, -4046.259987733932);
    ps.Momenta[2].Set(3884.766162969353, 2193.018847929061, 3018.732831200967,
                      -1081.354968273099);
    ps.Momenta[3].Set(1558.755317836453, -478.6509833616622, -999.3487780622738,
                      1096.318109307592);
    ps.Momenta[4].Set(2648.998494662059, -1714.367864567399, -2019.384053138693,
                      -14.96314103449324);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 3.676315573529796e-10;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(2.123041188490372e-10, 0);
    spin[1][2] =
        std::complex<double>(-3.058962217696e-11, 1.789998044125421e-10);
    spin[2][1] =
        std::complex<double>(-3.058962217696e-11, -1.789998044125421e-10);
    spin[2][2] = std::complex<double>(1.553274385039424e-10, 0);

    // radiation variables
    double phi = 0.1318222602254123;
    double xi = 0.06940492316162931;
    int bornpdgs[] = {0, 2, -13, 14, 1};
    int realpdgs[] = {-5, 2, -13, 14, 1, -5};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 2.391115326851988e-16
    double fks_g = 2.303621239984846e-26;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, -1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, 2, -13, 14, 1
// realpdgs = -5, 2, -13, 14, 1, -5
TEST(QCDSoftCollinearISR2, Wj_born_0_real_19) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.6031344727174321;
    ps.X2 = 0.6005094594826375;
    ps.Jacobian = 13814.01296135063;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3911.833476739665, 0, 0, 3911.833476739665);
    ps.Momenta[1].Set(3911.833476739665, 0, 0, -3911.833476739665);
    ps.Momenta[2].Set(2095.504957768069, -1346.933976266176, -1485.535827445184,
                      608.3528556578802);
    ps.Momenta[3].Set(2224.36370375025, 981.7531656985321, -1980.507024483182,
                      -248.0857395643927);
    ps.Momenta[4].Set(3503.798291961009, 365.180810567644, 3466.042851928367,
                      -360.2671160934876);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 2.394736063452918e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(1.027111171022898e-09, 0);
    spin[1][2] =
        std::complex<double>(1.134291877389902e-09, -3.436346048801768e-10);
    spin[2][1] =
        std::complex<double>(1.134291877389902e-09, 3.436346048801768e-10);
    spin[2][2] = std::complex<double>(1.36762489243002e-09, 0);

    // radiation variables
    double phi = 5.743426914415542;
    int bornpdgs[] = {0, 2, -13, 14, 1};
    int realpdgs[] = {-5, 2, -13, 14, 1, -5};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 3.872496106014531e-10
    double fks_g = 1.236044527526446e-28;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, -1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, 2, -13, 14, 1
// realpdgs = -4, 2, -13, 14, 1, -4
TEST(QCDCollinearISR1, Wj_born_0_real_19_n1) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.4365339044334426;
    ps.X2 = 0.8021165715072343;
    ps.Jacobian = 3155.415190101142;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3846.28172109989, 0, 0, 3846.28172109989);
    ps.Momenta[1].Set(3846.28172109989, 0, 0, -3846.28172109989);
    ps.Momenta[2].Set(3513.593247819825, 51.07699319285456, 1337.100274159706,
                      -3248.829252012713);
    ps.Momenta[3].Set(3364.987789388413, -758.0697402177634, -949.2074792121423,
                      3138.069191892846);
    ps.Momenta[4].Set(813.9824049915427, 706.9927470249088, -387.8927949475632,
                      110.7600601198674);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 4.572156082468887e-11;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(3.64860540129704e-11, 0);
    spin[1][2] =
        std::complex<double>(-7.813044156366166e-12, -1.661094643240381e-11);
    spin[2][1] =
        std::complex<double>(-7.813044156366166e-12, 1.661094643240381e-11);
    spin[2][2] = std::complex<double>(9.235506811718472e-12, 0);

    // radiation variables
    double phi = 0.0863628324046487;
    double xi = 0.1832355836469942;
    int bornpdgs[] = {0, 2, -13, 14, 1};
    int realpdgs[] = {-4, 2, -13, 14, 1, -4};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 5.010370831694503e-09
    double fks_g = 3.364491981050379e-18;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, 1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 3.364611101864092e-18
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, 2, -13, 14, 1
// realpdgs = -4, 2, -13, 14, 1, -4
TEST(QCDSoftCollinearISR1, Wj_born_0_real_19_n1) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.6581697044412103;
    ps.X2 = 0.5654048072951525;
    ps.Jacobian = 13355.60920842084;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3965.172165849114, 0, 0, 3965.172165849114);
    ps.Momenta[1].Set(3965.172165849114, 0, 0, -3965.172165849114);
    ps.Momenta[2].Set(2235.903525929745, 1369.889629958613, 1266.717940243869,
                      -1232.109020688207);
    ps.Momenta[3].Set(2352.480783618851, 491.0747221645143, -1742.034139894108,
                      -1502.773605700607);
    ps.Momenta[4].Set(3341.960022149632, -1860.964352123127, 475.3161996502388,
                      2734.882626388813);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 6.138338336589822e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(3.368590123153746e-09, 1.378634354255046e-25);
    spin[1][2] =
        std::complex<double>(1.615645041790443e-09, -2.592264950644794e-09);
    spin[2][1] =
        std::complex<double>(1.615645041790443e-09, 2.592264950644794e-09);
    spin[2][2] =
        std::complex<double>(2.769748213436075e-09, -1.378634354255046e-25);

    // radiation variables
    double phi = 0.3996473683041796;
    int bornpdgs[] = {0, 2, -13, 14, 1};
    int realpdgs[] = {-4, 2, -13, 14, 1, -4};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 0.002259516404627205
    double fks_g = 5.280398978815108e-22;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, 1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, 2, -13, 14, 1
// realpdgs = -4, 2, -13, 14, 1, -4
TEST(QCDSoftLimit, Wj_born_0_real_19_0_n1) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.1848955090061217;
    ps.X2 = 0.2741675282518159;
    ps.Jacobian = 2091.841518786698;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1463.472433328758, 0, 0, 1463.472433328758);
    ps.Momenta[1].Set(1463.472433328758, 0, 0, -1463.472433328758);
    ps.Momenta[2].Set(865.7673997485683, -303.9272215882445, 656.9445228156787,
                      474.9792925881375);
    ps.Momenta[3].Set(642.9567381659197, 160.3076961834113, 345.2348008490114,
                      518.1773267715113);
    ps.Momenta[4].Set(1418.220728743029, 143.6195254048331, -1002.17932366469,
                      -993.1566193596482);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 4.078708071004146e-08;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, 6.118062106506219e-08);
    ColorCorr.Set(0, 4, 6.118062106506219e-08);
    ColorCorr.Set(1, 0, 6.118062106506219e-08);
    ColorCorr.Set(1, 4, -6.79784678500691e-09);
    ColorCorr.Set(4, 0, 6.118062106506219e-08);
    ColorCorr.Set(4, 1, -6.79784678500691e-09);

    // radiation variables
    double y = 0.09731285498209985;
    double phi = 3.966582658529261;
    int bornpdgs[] = {0, 2, -13, 14, 1};
    int realpdgs[] = {-4, 2, -13, 14, 1, -4};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 4.452923858191712e-06
    double fks_g = 3.870046994787233e-22;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 0, alphas, ColorCorr, y, phi);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, 2, -13, 14, 1
// realpdgs = -4, 2, -13, 14, 1, -4
TEST(QCDCollinearISR2, Wj_born_0_real_19_n1) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.8881232617167321;
    ps.X2 = 0.6400433582662544;
    ps.Jacobian = 14580.97396082545;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(4900.661173561713, 0, 0, 4900.661173561713);
    ps.Momenta[1].Set(4900.661173561713, 0, 0, -4900.661173561713);
    ps.Momenta[2].Set(3523.722855054705, -2651.240527286832, -1939.12549730665,
                      -1275.671874504553);
    ps.Momenta[3].Set(3325.496847068133, 1416.975000063028, -109.2030859357069,
                      3006.523875687434);
    ps.Momenta[4].Set(2952.10264500059, 1234.265527223804, 2048.328583242357,
                      -1730.852001182881);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.621978160174791e-10;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(2.982222278410504e-11, 4.308232357047019e-27);
    spin[1][2] =
        std::complex<double>(-6.069690108535268e-11, -1.623639834503527e-11);
    spin[2][1] =
        std::complex<double>(-6.069690108535268e-11, 1.623639834503527e-11);
    spin[2][2] =
        std::complex<double>(1.323755932333741e-10, -4.308232357047019e-27);

    // radiation variables
    double phi = 0.4048459544894644;
    double xi = 0.1229472959326385;
    int bornpdgs[] = {0, 2, -13, 14, 1};
    int realpdgs[] = {-4, 2, -13, 14, 1, -4};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 4.577397726413229e-17
    double fks_g = 1.383842320031872e-26;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, -1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, 2, -13, 14, 1
// realpdgs = -4, 2, -13, 14, 1, -4
TEST(QCDSoftCollinearISR2, Wj_born_0_real_19_n1) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.7476324262441345;
    ps.X2 = 0.979440778511325;
    ps.Jacobian = 13883.77829364541;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(5562.198865254332, 0, 0, 5562.198865254332);
    ps.Momenta[1].Set(5562.198865254332, 0, 0, -5562.198865254332);
    ps.Momenta[2].Set(4735.547405585879, -4599.218015427701, -682.2910816370555,
                      898.3772910280901);
    ps.Momenta[3].Set(3912.222464212013, 2988.369378324803, -291.2684880165647,
                      -2508.046198949351);
    ps.Momenta[4].Set(2476.627860710775, 1610.848637102898, 973.5595696536202,
                      1609.668907921261);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.17146697539425e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(5.733950599617429e-10, 0);
    spin[1][2] =
        std::complex<double>(-3.162751639381401e-10, 4.928503854992501e-10);
    spin[2][1] =
        std::complex<double>(-3.162751639381401e-10, -4.928503854992501e-10);
    spin[2][2] = std::complex<double>(5.980719154325071e-10, 0);

    // radiation variables
    double phi = 1.710961452739064;
    int bornpdgs[] = {0, 2, -13, 14, 1};
    int realpdgs[] = {-4, 2, -13, 14, 1, -4};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 1.82066825684895e-09
    double fks_g = 1.53912665465774e-30;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, -1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, 2, -13, 14, 1
// realpdgs = -3, 2, -13, 14, 1, -3
TEST(QCDCollinearISR1, Wj_born_0_real_19_n2) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.02914244709976188;
    ps.X2 = 0.5267511478031111;
    ps.Jacobian = 637.3866418108926;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(805.3397032729235, 0, 0, 805.3397032729235);
    ps.Momenta[1].Set(805.3397032729235, 0, 0, -805.3397032729235);
    ps.Momenta[2].Set(575.394172899378, 72.19100906261019, 528.0703865169697,
                      -216.8146196659706);
    ps.Momenta[3].Set(250.0072104520985, -64.88967222821864, 237.6305539687647,
                      42.71598689779076);
    ps.Momenta[4].Set(785.2780231943702, -7.301336834391554, -765.7009404857349,
                      174.0986327681799);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.484716214069972e-07;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(8.489756090201895e-08, 0);
    spin[1][2] =
        std::complex<double>(-7.120145563281892e-08, -1.810070134767563e-08);
    spin[2][1] =
        std::complex<double>(-7.120145563281892e-08, 1.810070134767563e-08);
    spin[2][2] = std::complex<double>(6.357406050497826e-08, 0);

    // radiation variables
    double phi = 6.25725560532146;
    double xi = 0.506258231613151;
    int bornpdgs[] = {0, 2, -13, 14, 1};
    int realpdgs[] = {-3, 2, -13, 14, 1, -3};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 0.0002628751702915217
    double fks_g = 1.347484437323682e-12;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, 1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 1.347542810697811e-12
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, 2, -13, 14, 1
// realpdgs = -3, 2, -13, 14, 1, -3
TEST(QCDSoftCollinearISR1, Wj_born_0_real_19_n2) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.2806143011303668;
    ps.X2 = 0.7383029332514823;
    ps.Jacobian = 3745.93607674814;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2958.595237466174, 0, 0, 2958.595237466174);
    ps.Momenta[1].Set(2958.595237466174, 0, 0, -2958.595237466174);
    ps.Momenta[2].Set(2157.628063564697, -102.3554191299582, -648.6487523123743,
                      2055.270547879348);
    ps.Momenta[3].Set(2503.317242562832, 1334.372773032148, 862.4213898841371,
                      -1934.496282181279);
    ps.Momenta[4].Set(1256.245168804819, -1232.01735390219, -213.7726375717627,
                      -120.7742656980689);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 2.785521644414933e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(6.555937757655528e-10, 0);
    spin[1][2] =
        std::complex<double>(-5.593916615697355e-10, 1.040888285303651e-09);
    spin[2][1] =
        std::complex<double>(-5.593916615697355e-10, -1.040888285303651e-09);
    spin[2][2] = std::complex<double>(2.12992786864938e-09, 0);

    // radiation variables
    double phi = 1.461567424445045;
    int bornpdgs[] = {0, 2, -13, 14, 1};
    int realpdgs[] = {-3, 2, -13, 14, 1, -3};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 0.0008749389404291268
    double fks_g = 9.055895500549973e-22;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, 1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, 2, -13, 14, 1
// realpdgs = -3, 2, -13, 14, 1, -3
TEST(QCDSoftLimit, Wj_born_0_real_19_0_n2) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.4070154667237382;
    ps.X2 = 0.1719702131335481;
    ps.Jacobian = 749.8175208456482;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1719.671238844163, 0, 0, 1719.671238844163);
    ps.Momenta[1].Set(1719.671238844163, 0, 0, -1719.671238844163);
    ps.Momenta[2].Set(1429.989796462966, -128.3806097098807, -1026.85650071477,
                      986.8915664745198);
    ps.Momenta[3].Set(1576.729533527377, -277.3190125521074, 1174.481989680532,
                      -1014.772015326912);
    ps.Momenta[4].Set(432.6231476979833, 405.6996222619882, -147.6254889657622,
                      27.88044885239196);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.129075997889266e-08;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, 1.693613996833899e-08);
    ColorCorr.Set(0, 4, 1.693613996833899e-08);
    ColorCorr.Set(1, 0, 1.693613996833899e-08);
    ColorCorr.Set(1, 4, -1.881793329815443e-09);
    ColorCorr.Set(4, 0, 1.693613996833899e-08);
    ColorCorr.Set(4, 1, -1.881793329815443e-09);

    // radiation variables
    double y = 0.8484580724367916;
    double phi = 3.445467958691602;
    int bornpdgs[] = {0, 2, -13, 14, 1};
    int realpdgs[] = {-3, 2, -13, 14, 1, -3};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 7.960287281816456e-06
    double fks_g = 8.73120323516025e-23;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 0, alphas, ColorCorr, y, phi);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, 2, -13, 14, 1
// realpdgs = -3, 2, -13, 14, 1, -3
TEST(QCDCollinearISR2, Wj_born_0_real_19_n2) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.1074179441762588;
    ps.X2 = 0.3660360679235382;
    ps.Jacobian = 253.6527856807142;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1288.883652905649, 0, 0, 1288.883652905649);
    ps.Momenta[1].Set(1288.883652905649, 0, 0, -1288.883652905649);
    ps.Momenta[2].Set(1233.028792128426, -785.2875716406807, -948.1774956664428,
                      68.13859960462938);
    ps.Momenta[3].Set(1149.47301400331, 790.361282026086, 830.0611110939897,
                      87.26858338112046);
    ps.Momenta[4].Set(195.2654996795618, -5.073710385405374, 118.1163845724532,
                      -155.4071829857498);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 7.544572359371241e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(3.498104390294062e-09, 3.446585885637615e-26);
    spin[1][2] =
        std::complex<double>(-5.815997617654213e-10, 3.717083411080195e-09);
    spin[2][1] =
        std::complex<double>(-5.815997617654213e-10, -3.717083411080195e-09);
    spin[2][2] =
        std::complex<double>(4.046467969077179e-09, -3.446585885637615e-26);

    // radiation variables
    double phi = 4.974817671583712;
    double xi = 0.1068329360404132;
    int bornpdgs[] = {0, 2, -13, 14, 1};
    int realpdgs[] = {-3, 2, -13, 14, 1, -3};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 1.232127685988034e-14
    double fks_g = 2.81252272296076e-24;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, -1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, 2, -13, 14, 1
// realpdgs = -3, 2, -13, 14, 1, -3
TEST(QCDSoftCollinearISR2, Wj_born_0_real_19_n2) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.03872358237464724;
    ps.X2 = 0.8987647598908381;
    ps.Jacobian = 543.806636714748;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1212.618356630151, 0, 0, 1212.618356630151);
    ps.Momenta[1].Set(1212.618356630151, 0, 0, -1212.618356630151);
    ps.Momenta[2].Set(927.6919443852097, 32.42908957847723, 217.6423735325573,
                      -901.2172296785259);
    ps.Momenta[3].Set(1052.585803604775, -427.2357709647245, -21.54733808740215,
                      961.7391445604339);
    ps.Momenta[4].Set(444.9589652703172, 394.8066813862472, -196.0950354451552,
                      -60.52191488190807);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.609263886838498e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(8.700519809024351e-10, 4.308232357047019e-27);
    spin[1][2] =
        std::complex<double>(-3.585955486837056e-11, -8.011659474734172e-10);
    spin[2][1] =
        std::complex<double>(-3.585955486837056e-11, 8.011659474734172e-10);
    spin[2][2] =
        std::complex<double>(7.392119059360629e-10, -4.308232357047019e-27);

    // radiation variables
    double phi = 0.263880221102185;
    int bornpdgs[] = {0, 2, -13, 14, 1};
    int realpdgs[] = {-3, 2, -13, 14, 1, -3};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 1.068680529859501e-08
    double fks_g = 2.190491248773846e-28;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, -1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, 2, -13, 14, 1
// realpdgs = 0, 2, -13, 14, 1, 0
TEST(QCDCollinearFSR, Wj_born_0_real_23) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.9743712860102782;
    ps.X2 = 0.2328917387638185;
    ps.Jacobian = 70.07100557199294;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3096.3684731908, 0, 0, 3096.3684731908);
    ps.Momenta[1].Set(3096.3684731908, 0, 0, -3096.3684731908);
    ps.Momenta[2].Set(3094.235317601828, -2389.049171568186, 1888.842942053382,
                      -546.816785385171);
    ps.Momenta[3].Set(3076.048063529063, 2367.965319759478, -1883.942787942976,
                      552.78522530449);
    ps.Momenta[4].Set(22.45356525070805, 21.08385180870788, -4.900154110406025,
                      -5.968439919318976);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 2.179087539783625e-08;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(1.092419057150366e-08, 0);
    spin[1][2] =
        std::complex<double>(2.247465718998324e-11, 1.089537657968317e-08);
    spin[2][1] =
        std::complex<double>(2.247465718998324e-11, -1.089537657968317e-08);
    spin[2][2] = std::complex<double>(1.086668482633259e-08, 0);

    // radiation variables
    double phi = 2.078380818328196;
    double xi = 0.004721853436713573;
    int bornpdgs[] = {0, 2, -13, 14, 1};
    int realpdgs[] = {0, 2, -13, 14, 1, 0};

    // FKS g = real matrix element * (1.0 - y) * xi * xi
    // real = 0.06480824262611339
    double fks_g = 1.444958096108806e-14;
    double limit = FKS::QCD::CollinearLimitFSR(
        realpdgs[4], bornpdgs[4], ps, 4, xi, phi, alphas, born,
        spin); // limit = 1.444890064671212e-14
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, 2, -13, 14, 1
// realpdgs = 0, 2, -13, 14, 1, 0
TEST(QCDSoftCollinearFSR, Wj_born_0_real_23) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.7187146410578471;
    ps.X2 = 0.01475795139536729;
    ps.Jacobian = 250.2446287618149;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(669.4291822213391, 0, 0, 669.4291822213391);
    ps.Momenta[1].Set(669.4291822213391, 0, 0, -669.4291822213391);
    ps.Momenta[2].Set(592.7275202053801, -414.6324550774693, -215.1717012543475,
                      -364.838292102898);
    ps.Momenta[3].Set(375.2283785204203, 215.5494422182895, -89.59924077123694,
                      293.7801049418867);
    ps.Momenta[4].Set(370.902465716878, 199.0830128591799, 304.7709420255845,
                      71.05818716101122);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 5.376460375549586e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(2.170263331147216e-10, -3.446585885637615e-26);
    spin[1][2] =
        std::complex<double>(3.101905659698315e-10, 1.011689114283741e-09);
    spin[2][1] =
        std::complex<double>(3.101905659698315e-10, -1.011689114283741e-09);
    spin[2][2] =
        std::complex<double>(5.159434042434864e-09, 3.446585885637615e-26);

    // radiation variables
    double phi = 1.256733353057939;
    int bornpdgs[] = {0, 2, -13, 14, 1};
    int realpdgs[] = {0, 2, -13, 14, 1, 0};

    // FKS g = real matrix element * (1.0 - y) * xi * xi
    // real = 154472.0318525586
    double fks_g = 4.741983287785126e-14;
    double limit = FKS::QCD::SoftCollinearLimitFSR(realpdgs[4], bornpdgs[4], ps,
                                                   4, phi, alphas, born, spin);
    // limit = 4.744033771632552e-14
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, 2, -13, 14, 1
// realpdgs = 0, 2, -13, 14, 1, 0
TEST(QCDSoftLimit, Wj_born_0_real_23_4) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.05491719826602193;
    ps.X2 = 0.6685818979208236;
    ps.Jacobian = 439.2289118454129;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1245.503206041367, 0, 0, 1245.503206041367);
    ps.Momenta[1].Set(1245.503206041367, 0, 0, -1245.503206041367);
    ps.Momenta[2].Set(925.0089019889519, -211.6103148578804, -446.6251753602844,
                      781.9133559026538);
    ps.Momenta[3].Set(1216.096121466739, 449.7048675769544, 667.9976851984165,
                      -911.2817354109525);
    ps.Momenta[4].Set(349.901388627043, -238.094552719074, -221.3725098381321,
                      129.3683795082987);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    // double born = 3.059697767219673e-08;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, 4.589546650829509e-08);
    ColorCorr.Set(0, 4, 4.589546650829509e-08);
    ColorCorr.Set(1, 0, 4.589546650829509e-08);
    ColorCorr.Set(1, 4, -5.099496278699454e-09);
    ColorCorr.Set(4, 0, 4.589546650829509e-08);
    ColorCorr.Set(4, 1, -5.099496278699454e-09);

    // radiation variables
    double y = 0.708127652639682;
    double phi = 4.022286959599939;
    int bornpdgs[] = {0, 2, -13, 14, 1};
    int realpdgs[] = {0, 2, -13, 14, 1, 0};

    // FKS g = real matrix element * (1.0 - y) * xi * xi
    // real = 112748.5784737567
    double fks_g = 2.597201580188214e-13;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 4, alphas, ColorCorr, y, phi);
    // limit = 2.59720136674225e-13
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, 2, -13, 14, 1
// realpdgs = 0, 2, -13, 14, 1, 0
TEST(QCDCollinearISR1, Wj_born_0_real_23) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.9215670742351882;
    ps.X2 = 0.911077829010253;
    ps.Jacobian = 21362.79942626762;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(5955.998376606937, 0, 0, 5955.998376606937);
    ps.Momenta[1].Set(5955.998376606937, 0, 0, -5955.998376606937);
    ps.Momenta[2].Set(3487.013621582384, -432.8330176536454, 3281.608692023625,
                      1096.797140935789);
    ps.Momenta[3].Set(4866.18635410004, -2416.274797134338, -2789.506866650614,
                      -3171.756166301553);
    ps.Momenta[4].Set(3558.796777531452, 2849.107814787983, -492.1018253730119,
                      2074.959025365764);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.011770677908468e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(4.420204587491116e-11, 3.446585885637615e-26);
    spin[1][2] =
        std::complex<double>(-2.033623021072816e-10, -3.758040901718621e-11);
    spin[2][1] =
        std::complex<double>(-2.033623021072816e-10, 3.758040901718621e-11);
    spin[2][2] =
        std::complex<double>(9.675686320335573e-10, -3.446585885637615e-26);

    // radiation variables
    double phi = 3.547351312700729;
    double xi = 0.0383661602340967;
    int bornpdgs[] = {0, 2, -13, 14, 1};
    int realpdgs[] = {0, 2, -13, 14, 1, 0};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 1.660251570551556e-05
    double fks_g = 4.887655275817947e-16;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, 1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 4.888210899075372e-16
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, 2, -13, 14, 1
// realpdgs = 0, 2, -13, 14, 1, 0
TEST(QCDSoftCollinearISR1, Wj_born_0_real_23) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.8049266763906751;
    ps.X2 = 0.4210238785066807;
    ps.Jacobian = 9976.954820581077;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3783.945571558424, 0, 0, 3783.945571558424);
    ps.Momenta[1].Set(3783.945571558424, 0, 0, -3783.945571558424);
    ps.Momenta[2].Set(3144.184169022034, 2998.944019088991, 219.1253460014245,
                      -918.8106126068081);
    ps.Momenta[3].Set(1807.616730776586, -510.6404635195026, -174.3520638112886,
                      1725.203153325649);
    ps.Momenta[4].Set(2616.090243318229, -2488.303555569489, -44.77328219013582,
                      -806.392540718841);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 3.239851572030173e-10;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(2.719555366903029e-11, 0);
    spin[1][2] =
        std::complex<double>(-3.780957609865963e-11, 8.149719962402654e-11);
    spin[2][1] =
        std::complex<double>(-3.780957609865963e-11, -8.149719962402654e-11);
    spin[2][2] = std::complex<double>(2.96789603533987e-10, 0);

    // radiation variables
    double phi = 0.8219008549437927;
    int bornpdgs[] = {0, 2, -13, 14, 1};
    int realpdgs[] = {0, 2, -13, 14, 1, 0};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 5292.244560440449
    double fks_g = 4.027780940955125e-16;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, 1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 4.026321682920157e-16
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, 2, -13, 14, 1
// realpdgs = 0, 2, -13, 14, 1, 0
TEST(QCDSoftLimit, Wj_born_0_real_23_0) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.832835071459602;
    ps.X2 = 0.5616935147702549;
    ps.Jacobian = 5298.738467915391;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(4445.724684698254, 0, 0, 4445.724684698254);
    ps.Momenta[1].Set(4445.724684698254, 0, 0, -4445.724684698254);
    ps.Momenta[2].Set(3871.168472707985, 678.0190801422261, 3010.138698478335,
                      2337.798213485184);
    ps.Momenta[3].Set(3837.703655712955, -87.04268895328238, -2297.022642880415,
                      -3073.122174184918);
    ps.Momenta[4].Set(1182.577240975569, -590.9763911889437, -713.1160555979206,
                      735.3239606997344);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    // double born = 3.887439241223282e-09;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, 5.831158861834923e-09);
    ColorCorr.Set(0, 4, 5.831158861834923e-09);
    ColorCorr.Set(1, 0, 5.831158861834923e-09);
    ColorCorr.Set(1, 4, -6.479065402038803e-10);
    ColorCorr.Set(4, 0, 5.831158861834923e-09);
    ColorCorr.Set(4, 1, -6.479065402038803e-10);

    // radiation variables
    double y = -0.4334770480156394;
    double phi = 3.332784202005477;
    int bornpdgs[] = {0, 2, -13, 14, 1};
    int realpdgs[] = {0, 2, -13, 14, 1, 0};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 94.98347225982521
    double fks_g = 1.696705714148526e-15;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 0, alphas, ColorCorr, y, phi);
    // limit = 1.696705774830263e-15
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, 2, -13, 14, 1
// realpdgs = 0, 2, -13, 14, 1, 0
TEST(QCDCollinearISR2, Wj_born_0_real_23) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.1781111741840071;
    ps.X2 = 0.5937905623201445;
    ps.Jacobian = 1243.673185776798;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2113.856906956077, 0, 0, 2113.856906956077);
    ps.Momenta[1].Set(2113.856906956077, 0, 0, -2113.856906956077);
    ps.Momenta[2].Set(1696.429568267045, -1115.645204702526, -303.9986383426208,
                      1241.287188849061);
    ps.Momenta[3].Set(1947.529691267541, 839.5582722063161, 626.8729989001688,
                      -1641.658932052726);
    ps.Momenta[4].Set(583.7545543775685, 276.0869324962099, -322.874360557548,
                      400.3717432036648);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 2.0937781353547e-08;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(1.303187876796033e-08, -5.514537417020184e-25);
    spin[1][2] =
        std::complex<double>(4.872611687813634e-09, 8.904292188964356e-09);
    spin[2][1] =
        std::complex<double>(4.872611687813634e-09, -8.904292188964356e-09);
    spin[2][2] =
        std::complex<double>(7.905902585586679e-09, 5.514537417020184e-25);

    // radiation variables
    double phi = 2.849200847758474;
    double xi = 0.2060717881833477;
    int bornpdgs[] = {0, 2, -13, 14, 1};
    int realpdgs[] = {0, 2, -13, 14, 1, 0};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 3.556603727251765e-05
    double fks_g = 3.020664934658046e-14;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, -1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 3.020740955261091e-14
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = 0, 2, -13, 14, 1
// realpdgs = 0, 2, -13, 14, 1, 0
TEST(QCDSoftCollinearISR2, Wj_born_0_real_23) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.07788031208940893;
    ps.X2 = 0.3604197285868338;
    ps.Jacobian = 1065.896495298269;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1089.009017386143, 0, 0, 1089.009017386143);
    ps.Momenta[1].Set(1089.009017386143, 0, 0, -1089.009017386143);
    ps.Momenta[2].Set(183.8837871019377, 116.3725540775648, -43.25208899185179,
                      135.6463512751299);
    ps.Momenta[3].Set(1022.991313332313, -896.6018415712344, 327.5250841590013,
                      367.890858939278);
    ps.Momenta[4].Set(971.1429343380358, 780.2292874936696, -284.2729951671495,
                      -503.537210214408);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 4.727588137433381e-08;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(1.746959964128405e-08, 0);
    spin[1][2] =
        std::complex<double>(5.238492728828552e-09, -2.220950252936134e-08);
    spin[2][1] =
        std::complex<double>(5.238492728828552e-09, 2.220950252936134e-08);
    spin[2][2] = std::complex<double>(2.980628173304976e-08, 0);

    // radiation variables
    double phi = 4.129582704703264;
    int bornpdgs[] = {0, 2, -13, 14, 1};
    int realpdgs[] = {0, 2, -13, 14, 1, 0};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 385431.2523075554
    double fks_g = 3.153313266551151e-13;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, -1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 3.152588461720709e-13
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 0, -13, 14, -2
// realpdgs = -1, 2, -13, 14, -2, 2
TEST(QCDCollinearISR1, Wj_born_4_real_6) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.6583042717237895;
    ps.X2 = 0.666011429169977;
    ps.Jacobian = 15965.8548677505;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(4303.953140249933, 0, 0, 4303.953140249933);
    ps.Momenta[1].Set(4303.953140249933, 0, 0, -4303.953140249933);
    ps.Momenta[2].Set(1090.957581360972, 378.9436313548752, 965.6154376979928,
                      -337.9008657341301);
    ps.Momenta[3].Set(3836.301091060253, -986.7429944869664, 452.6732487563417,
                      3679.487906483912);
    ps.Momenta[4].Set(3680.647608078642, 607.7993631320912, -1418.288686454335,
                      -3341.587040749783);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 6.467443593539917e-10;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(2.340671077604805e-10, 8.616464714094038e-27);
    spin[1][2] =
        std::complex<double>(-5.610893235476996e-11, -3.05689316838265e-10);
    spin[2][1] =
        std::complex<double>(-5.610893235476996e-11, 3.05689316838265e-10);
    spin[2][2] =
        std::complex<double>(4.126772515935111e-10, -8.616464714094038e-27);

    // radiation variables
    double phi = 1.293662573573874;
    double xi = 0.2258622686082057;
    int bornpdgs[] = {-1, 0, -13, 14, -2};
    int realpdgs[] = {-1, 2, -13, 14, -2, 2};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 6.962806836550341e-16
    double fks_g = 7.103979744043131e-25;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, 1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 0, -13, 14, -2
// realpdgs = -1, 2, -13, 14, -2, 2
TEST(QCDSoftCollinearISR1, Wj_born_4_real_6) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.6889736966259257;
    ps.X2 = 0.7630228342169829;
    ps.Jacobian = 16604.59503280917;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(4712.848130281165, 0, 0, 4712.848130281165);
    ps.Momenta[1].Set(4712.848130281165, 0, 0, -4712.848130281165);
    ps.Momenta[2].Set(2780.649784494886, 705.7114505272868, 1431.465737886202,
                      -2277.035444138725);
    ps.Momenta[3].Set(3149.263689816548, -2735.664729710488, 1142.390254782761,
                      1062.518131794633);
    ps.Momenta[4].Set(3495.782786250896, 2029.953279183201, -2573.855992668963,
                      1214.517312344092);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 3.122829394636801e-10;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(1.195905433754963e-11, 0);
    spin[1][2] =
        std::complex<double>(3.339900925381544e-13, -5.99289421874398e-11);
    spin[2][1] =
        std::complex<double>(3.339900925381544e-13, 5.99289421874398e-11);
    spin[2][2] = std::complex<double>(3.003238851261305e-10, 0);

    // radiation variables
    double phi = 0.6481809581955211;
    int bornpdgs[] = {-1, 0, -13, 14, -2};
    int realpdgs[] = {-1, 2, -13, 14, -2, 2};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 3.31352427544632e-10
    double fks_g = 6.410834114906836e-29;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, 1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 0, -13, 14, -2
// realpdgs = -1, 2, -13, 14, -2, 2
TEST(QCDSoftLimit, Wj_born_4_real_6_0) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.5892361603354708;
    ps.X2 = 0.2809019799404311;
    ps.Jacobian = 4982.77733650028;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2644.450561616417, 0, 0, 2644.450561616417);
    ps.Momenta[1].Set(2644.450561616417, 0, 0, -2644.450561616417);
    ps.Momenta[2].Set(1972.76281513582, 997.5914781210447, 1497.519234727158,
                      -808.7277101627368);
    ps.Momenta[3].Set(1446.594592975565, 446.235772490678, -420.868906014979,
                      1310.106452058588);
    ps.Momenta[4].Set(1869.54371512145, -1443.827250611723, -1076.650328712179,
                      -501.3787418958511);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 2.187159079950194e-09;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, 3.280738619925292e-09);
    ColorCorr.Set(0, 4, -3.645265133250324e-10);
    ColorCorr.Set(1, 0, 3.280738619925292e-09);
    ColorCorr.Set(1, 4, 3.280738619925292e-09);
    ColorCorr.Set(4, 0, -3.645265133250324e-10);
    ColorCorr.Set(4, 1, 3.280738619925292e-09);

    // radiation variables
    double y = 0.05962930199182637;
    double phi = 4.277839347281578;
    int bornpdgs[] = {-1, 0, -13, 14, -2};
    int realpdgs[] = {-1, 2, -13, 14, -2, 2};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 5.253122299787856e-07
    double fks_g = 2.097907636548803e-23;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 0, alphas, ColorCorr, y, phi);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 0, -13, 14, -2
// realpdgs = -1, 2, -13, 14, -2, 2
TEST(QCDCollinearISR2, Wj_born_4_real_6) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.4840872684855277;
    ps.X2 = 0.3055641733494667;
    ps.Jacobian = 3248.032426999542;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2499.921683674266, 0, 0, 2499.921683674266);
    ps.Momenta[1].Set(2499.921683674266, 0, 0, -2499.921683674266);
    ps.Momenta[2].Set(1498.489323217201, 920.4135896108039, -683.0898813657195,
                      965.2446787377103);
    ps.Momenta[3].Set(2212.233441599615, -284.8049861301422, 676.9758516050063,
                      -2086.759836768659);
    ps.Momenta[4].Set(1289.120602531716, -635.6086034806616, 6.114029760713271,
                      1121.515158030949);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 2.902874296125871e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(1.215121572898997e-09, 3.446585885637615e-26);
    spin[1][2] =
        std::complex<double>(2.997280377770611e-10, 1.40035275808716e-09);
    spin[2][1] =
        std::complex<double>(2.997280377770611e-10, -1.40035275808716e-09);
    spin[2][2] =
        std::complex<double>(1.687752723226873e-09, -3.446585885637615e-26);

    // radiation variables
    double phi = 5.04250334475462;
    double xi = 0.2800551903277695;
    int bornpdgs[] = {-1, 0, -13, 14, -2};
    int realpdgs[] = {-1, 2, -13, 14, -2, 2};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 4.916953410142555e-07
    double fks_g = 7.712822567209646e-16;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, -1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 7.7128101042198e-16
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 0, -13, 14, -2
// realpdgs = -1, 2, -13, 14, -2, 2
TEST(QCDSoftCollinearISR2, Wj_born_4_real_6) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.05926845292215077;
    ps.X2 = 0.7412695882847018;
    ps.Jacobian = 1509.91234253284;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1362.427005989926, 0, 0, 1362.427005989926);
    ps.Momenta[1].Set(1362.427005989926, 0, 0, -1362.427005989926);
    ps.Momenta[2].Set(930.9595583864123, 679.7585004975987, 193.6056555555771,
                      605.9133027828116);
    ps.Momenta[3].Set(694.2859543551259, -164.4971580237395, -619.3144142930321,
                      267.2514315505148);
    ps.Momenta[4].Set(1099.608499238314, -515.261342473859, 425.708758737455,
                      -873.1647343333264);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 5.054661638085392e-08;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(1.136208975623299e-08, 0);
    spin[1][2] =
        std::complex<double>(4.61378979590624e-09, -2.058958601535339e-08);
    spin[2][1] =
        std::complex<double>(4.61378979590624e-09, 2.058958601535339e-08);
    spin[2][2] = std::complex<double>(3.918452662462093e-08, 0);

    // radiation variables
    double phi = 3.33109571451458;
    int bornpdgs[] = {-1, 0, -13, 14, -2};
    int realpdgs[] = {-1, 2, -13, 14, -2, 2};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 0.2081058712049388
    double fks_g = 2.78618178602076e-20;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, -1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 0, -13, 14, -2
// realpdgs = -1, 1, -13, 14, -2, 1
TEST(QCDCollinearISR1, Wj_born_4_real_24) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.5424122701911536;
    ps.X2 = 0.2717804181093726;
    ps.Jacobian = 5336.251510118984;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2495.670184291925, 0, 0, 2495.670184291925);
    ps.Momenta[1].Set(2495.670184291925, 0, 0, -2495.670184291925);
    ps.Momenta[2].Set(1129.213577583269, -701.001427526294, 856.7436383938066,
                      222.9588313484046);
    ps.Momenta[3].Set(1740.599139793119, 8.627453748842584, 551.1468094881309,
                      -1651.014272161634);
    ps.Momenta[4].Set(2121.527651207462, 692.3739737774515, -1407.890447881937,
                      1428.05544081323);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 4.701704712465863e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(2.644410871023415e-09, 0);
    spin[1][2] =
        std::complex<double>(1.722260688885461e-09, 1.572942566885297e-09);
    spin[2][1] =
        std::complex<double>(1.722260688885461e-09, -1.572942566885297e-09);
    spin[2][2] = std::complex<double>(2.057293841442448e-09, 0);

    // radiation variables
    double phi = 5.944974155543912;
    double xi = 0.07799509584280274;
    int bornpdgs[] = {-1, 0, -13, 14, -2};
    int realpdgs[] = {-1, 1, -13, 14, -2, 1};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 9.763976895560439e-15
    double fks_g = 1.187931314401134e-24;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, 1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 0, -13, 14, -2
// realpdgs = -1, 1, -13, 14, -2, 1
TEST(QCDSoftCollinearISR1, Wj_born_4_real_24) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.1957743877853453;
    ps.X2 = 0.4169797306076006;
    ps.Jacobian = 1365.451127178298;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1857.157626581788, 0, 0, 1857.157626581788);
    ps.Momenta[1].Set(1857.157626581788, 0, 0, -1857.157626581788);
    ps.Momenta[2].Set(1788.471354285457, -589.8231406677731, -1023.136870399386,
                      -1343.104386211962);
    ps.Momenta[3].Set(1196.341052821614, 593.9332832353317, 404.5352762709379,
                      956.4655665438362);
    ps.Momenta[4].Set(729.5028460565069, -4.110142567558642, 618.6015941284476,
                      386.6388196681257);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.170449595197742e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(5.650562111731183e-10, 0);
    spin[1][2] =
        std::complex<double>(-2.046695702354734e-10, -5.478974893772128e-10);
    spin[2][1] =
        std::complex<double>(-2.046695702354734e-10, 5.478974893772128e-10);
    spin[2][2] = std::complex<double>(6.053933840246239e-10, 0);

    // radiation variables
    double phi = 2.619407428176251;
    int bornpdgs[] = {-1, 0, -13, 14, -2};
    int realpdgs[] = {-1, 1, -13, 14, -2, 1};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 4.171362413556861e-10
    double fks_g = 5.395898375662348e-28;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, 1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 0, -13, 14, -2
// realpdgs = -1, 1, -13, 14, -2, 1
TEST(QCDSoftLimit, Wj_born_4_real_24_0) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.5831793475483131;
    ps.X2 = 0.8094651315103931;
    ps.Jacobian = 16154.8999277954;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(4465.946307516766, 0, 0, 4465.946307516766);
    ps.Momenta[1].Set(4465.946307516766, 0, 0, -4465.946307516766);
    ps.Momenta[2].Set(2348.998814581381, -1438.491337803678, -1449.792481314699,
                      -1160.44813028928);
    ps.Momenta[3].Set(2993.754030916136, 2503.580742106913, -774.1784324169301,
                      -1447.513184796075);
    ps.Momenta[4].Set(3589.139769536016, -1065.089404303235, 2223.970913731629,
                      2607.961315085354);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 7.373046184606927e-10;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, 1.105956927691039e-09);
    ColorCorr.Set(0, 4, -1.228841030767821e-10);
    ColorCorr.Set(1, 0, 1.105956927691039e-09);
    ColorCorr.Set(1, 4, 1.105956927691039e-09);
    ColorCorr.Set(4, 0, -1.228841030767821e-10);
    ColorCorr.Set(4, 1, 1.105956927691039e-09);

    // radiation variables
    double y = -0.5201622316242762;
    double phi = 4.214686406527258;
    int bornpdgs[] = {-1, 0, -13, 14, -2};
    int realpdgs[] = {-1, 1, -13, 14, -2, 1};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 6.265190738520101e-08
    double fks_g = 2.701235513322061e-25;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 0, alphas, ColorCorr, y, phi);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 0, -13, 14, -2
// realpdgs = -1, 1, -13, 14, -2, 1
TEST(QCDCollinearISR2, Wj_born_4_real_24) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.8494484015377388;
    ps.X2 = 0.2233193333053478;
    ps.Jacobian = 7247.009851093999;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2831.033573174756, 0, 0, 2831.033573174756);
    ps.Momenta[1].Set(2831.033573174756, 0, 0, -2831.033573174756);
    ps.Momenta[2].Set(1600.352880754225, 1175.974513309047, 1065.524036200256,
                      -207.0551020028527);
    ps.Momenta[3].Set(1521.832919904066, -422.4239331937287, 1324.896722177707,
                      618.2088096522418);
    ps.Momenta[4].Set(2539.881345691221, -753.5505801153188, -2390.420758377963,
                      -411.1537076493891);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 4.360062294785545e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(1.235241806853128e-09, 0);
    spin[1][2] =
        std::complex<double>(9.566516727231199e-10, -1.716020536788931e-09);
    spin[2][1] =
        std::complex<double>(9.566516727231199e-10, 1.716020536788931e-09);
    spin[2][2] = std::complex<double>(3.124820487932416e-09, 0);

    // radiation variables
    double phi = 0.9538729228871103;
    double xi = 0.6277148384956048;
    int bornpdgs[] = {-1, 0, -13, 14, -2};
    int realpdgs[] = {-1, 1, -13, 14, -2, 1};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 9.63819796765441e-07
    double fks_g = 7.595399609156579e-15;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, -1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 7.594396478017691e-15
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 0, -13, 14, -2
// realpdgs = -1, 1, -13, 14, -2, 1
TEST(QCDSoftCollinearISR2, Wj_born_4_real_24) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.282891613856993;
    ps.X2 = 0.3291447353940242;
    ps.Jacobian = 1004.274079722191;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1983.429872127894, 0, 0, 1983.429872127894);
    ps.Momenta[1].Set(1983.429872127894, 0, 0, -1983.429872127894);
    ps.Momenta[2].Set(1953.168826660653, -40.86646042857983, 854.1837726738326,
                      1756.009248367357);
    ps.Momenta[3].Set(1511.307842034205, -84.44641091844593, -812.9069866665512,
                      -1271.260173255572);
    ps.Momenta[4].Set(502.3830755609291, 125.3128713470258, -41.27678600728137,
                      -484.7490751117851);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 2.365528573337477e-07;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(6.380129927274095e-08, 0);
    spin[1][2] =
        std::complex<double>(4.705382669518386e-08, 9.384940729180009e-08);
    spin[2][1] =
        std::complex<double>(4.705382669518386e-08, -9.384940729180009e-08);
    spin[2][2] = std::complex<double>(1.727515580610067e-07, 0);

    // radiation variables
    double phi = 5.268591207941529;
    int bornpdgs[] = {-1, 0, -13, 14, -2};
    int realpdgs[] = {-1, 1, -13, 14, -2, 1};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 0.1772143707117628
    double fks_g = 1.595095422184376e-19;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, -1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 0, -13, 14, -2
// realpdgs = -1, -2, -13, 14, -2, -2
TEST(QCDCollinearISR1, Wj_born_4_real_9) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.9014421143360742;
    ps.X2 = 0.721440222568269;
    ps.Jacobian = 2270.857068694535;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(5241.824237139104, 0, 0, 5241.824237139104);
    ps.Momenta[1].Set(5241.824237139104, 0, 0, -5241.824237139104);
    ps.Momenta[2].Set(4823.018792039478, 4453.462494500362, -1117.842834313644,
                      -1476.01134013875);
    ps.Momenta[3].Set(5230.789569252193, -4882.526069769167, 1137.615218578989,
                      1492.625308066114);
    ps.Momenta[4].Set(429.8401129865361, 429.0635752688057, -19.77238426534533,
                      -16.61396792736481);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 7.192631517880075e-10;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(3.5944817762576e-10, 0);
    spin[1][2] =
        std::complex<double>(4.242267995754086e-12, 3.596065070425339e-10);
    spin[2][1] =
        std::complex<double>(4.242267995754086e-12, -3.596065070425339e-10);
    spin[2][2] = std::complex<double>(3.598149741622475e-10, 0);

    // radiation variables
    double phi = 6.26064770955903;
    double xi = 0.08482166753113742;
    int bornpdgs[] = {-1, 0, -13, 14, -2};
    int realpdgs[] = {-1, -2, -13, 14, -2, -2};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 1.025013534682163e-16
    double fks_g = 1.474936107826956e-26;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, 1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 0, -13, 14, -2
// realpdgs = -1, -2, -13, 14, -2, -2
TEST(QCDSoftCollinearISR1, Wj_born_4_real_9) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.4143059261034274;
    ps.X2 = 0.9416247086858789;
    ps.Jacobian = 10174.09338155796;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(4059.876777335781, 0, 0, 4059.876777335781);
    ps.Momenta[1].Set(4059.876777335781, 0, 0, -4059.876777335781);
    ps.Momenta[2].Set(3975.547479743667, 1017.971159028678, 3841.765623173489,
                      -97.72092756645125);
    ps.Momenta[3].Set(1657.740427530126, -25.43063071381073, -1612.983167656279,
                      -381.7615865259097);
    ps.Momenta[4].Set(2486.46564739777, -992.5405283148677, -2228.78245551721,
                      479.482514092361);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 9.96281942802874e-10;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(8.219465821598912e-10, 8.616464714094038e-27);
    spin[1][2] =
        std::complex<double>(-1.535268757400612e-10, -3.460113470041228e-10);
    spin[2][1] =
        std::complex<double>(-1.535268757400612e-10, 3.460113470041228e-10);
    spin[2][2] =
        std::complex<double>(1.743353606429827e-10, -8.616464714094038e-27);

    // radiation variables
    double phi = 5.034580117839345;
    int bornpdgs[] = {-1, 0, -13, 14, -2};
    int realpdgs[] = {-1, -2, -13, 14, -2, -2};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 5.101013740507819e-11
    double fks_g = 3.499679218839705e-29;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, 1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 0, -13, 14, -2
// realpdgs = -1, -2, -13, 14, -2, -2
TEST(QCDSoftLimit, Wj_born_4_real_9_0) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.545108810945532;
    ps.X2 = 0.3017069279213231;
    ps.Jacobian = 2595.327331036428;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2636.013310850046, 0, 0, 2636.013310850046);
    ps.Momenta[1].Set(2636.013310850046, 0, 0, -2636.013310850046);
    ps.Momenta[2].Set(1974.343798572266, 1186.265750305128, -921.3328569474838,
                      1281.387049771755);
    ps.Momenta[3].Set(2320.796256179683, -684.6570490502937, 1756.49129198643,
                      -1353.542880391759);
    ps.Momenta[4].Set(976.8865669481415, -501.6087012548343, -835.1584350389461,
                      72.15583062000381);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 3.233623971389079e-09;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, 4.850435957083619e-09);
    ColorCorr.Set(0, 4, -5.389373285648465e-10);
    ColorCorr.Set(1, 0, 4.850435957083619e-09);
    ColorCorr.Set(1, 4, 4.850435957083619e-09);
    ColorCorr.Set(4, 0, -5.389373285648465e-10);
    ColorCorr.Set(4, 1, 4.850435957083619e-09);

    // radiation variables
    double y = 0.4155377571031167;
    double phi = 0.7564062448911092;
    int bornpdgs[] = {-1, 0, -13, 14, -2};
    int realpdgs[] = {-1, -2, -13, 14, -2, -2};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 5.59357092960155e-08
    double fks_g = 1.562196666730791e-24;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 0, alphas, ColorCorr, y, phi);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 0, -13, 14, -2
// realpdgs = -1, -2, -13, 14, -2, -2
TEST(QCDCollinearISR2, Wj_born_4_real_9) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.9367326188745331;
    ps.X2 = 0.8411321378715826;
    ps.Jacobian = 2256.252699661956;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(5769.700790453285, 0, 0, 5769.700790453285);
    ps.Momenta[1].Set(5769.700790453285, 0, 0, -5769.700790453285);
    ps.Momenta[2].Set(5702.936569379263, 772.0294862537644, 181.6263389227119,
                      -5647.518734783233);
    ps.Momenta[3].Set(5448.462940669719, -763.8293436994968, 111.5529761504435,
                      5393.502487536914);
    ps.Momenta[4].Set(388.0020708575893, -8.200142554267606, -293.1793150731554,
                      254.0162472463188);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.016550444375706e-11;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(4.241263790021424e-12, -3.365806528942984e-29);
    spin[1][2] =
        std::complex<double>(1.927805104259745e-13, -5.008902379050893e-12);
    spin[2][1] =
        std::complex<double>(1.927805104259745e-13, 5.008902379050893e-12);
    spin[2][2] =
        std::complex<double>(5.924240653735636e-12, 3.365806528942984e-29);

    // radiation variables
    double phi = 6.079851959293325;
    double xi = 0.1535116869398976;
    int bornpdgs[] = {-1, 0, -13, 14, -2};
    int realpdgs[] = {-1, -2, -13, 14, -2, -2};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 2.258600328393997e-10
    double fks_g = 1.064516189578938e-19;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, -1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 1.064541520790057e-19
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 0, -13, 14, -2
// realpdgs = -1, -2, -13, 14, -2, -2
TEST(QCDSoftCollinearISR2, Wj_born_4_real_9) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.009178326586404495;
    ps.X2 = 0.9950359385777787;
    ps.Jacobian = 310.8789363705209;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(621.1757506538523, 0, 0, 621.1757506538523);
    ps.Momenta[1].Set(621.1757506538523, 0, 0, -621.1757506538523);
    ps.Momenta[2].Set(320.7205695798097, 56.31755226547173, 299.774678851548,
                      -99.12193994204637);
    ps.Momenta[3].Set(425.0655945248027, 349.0530092786346, -48.65153816745773,
                      237.6463426947569);
    ps.Momenta[4].Set(496.5653372030922, -405.3705615441064, -251.1231406840902,
                      -138.5244027527106);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 3.219378842375811e-08;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(1.945187195643222e-08, 0);
    spin[1][2] =
        std::complex<double>(1.20661763522366e-08, -1.0112443613623e-08);
    spin[2][1] =
        std::complex<double>(1.20661763522366e-08, 1.0112443613623e-08);
    spin[2][2] = std::complex<double>(1.27419164673259e-08, 0);

    // radiation variables
    double phi = 0.4238661813638395;
    int bornpdgs[] = {-1, 0, -13, 14, -2};
    int realpdgs[] = {-1, -2, -13, 14, -2, -2};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 15.56398947845902
    double fks_g = 7.670531069656586e-22;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, -1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 0, -13, 14, -2
// realpdgs = -1, 3, -13, 14, -2, 3
TEST(QCDCollinearISR1, Wj_born_4_real_16) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.105414401942884;
    ps.X2 = 0.4578377789482317;
    ps.Jacobian = 508.9548225601559;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1427.970199762758, 0, 0, 1427.970199762758);
    ps.Momenta[1].Set(1427.970199762758, 0, 0, -1427.970199762758);
    ps.Momenta[2].Set(1205.858493716329, -1072.769572750395, -326.0057133883234,
                      -443.8247688974996);
    ps.Momenta[3].Set(1296.443286769212, 1135.382212824755, 28.92926066795917,
                      625.1683969038468);
    ps.Momenta[4].Set(353.638619039975, -62.6126400743597, 297.0764527203642,
                      -181.3436280063471);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 3.880485788340332e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(2.307070775817247e-09, 0);
    spin[1][2] =
        std::complex<double>(-1.266154860874173e-09, 1.423668382000823e-09);
    spin[2][1] =
        std::complex<double>(-1.266154860874173e-09, -1.423668382000823e-09);
    spin[2][2] = std::complex<double>(1.573415012523085e-09, 0);

    // radiation variables
    double phi = 2.466165458456015;
    double xi = 0.7185739381721269;
    int bornpdgs[] = {-1, 0, -13, 14, -2};
    int realpdgs[] = {-1, 3, -13, 14, -2, 3};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 5.171618932271124e-16
    double fks_g = 5.340715401476299e-24;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, 1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 0, -13, 14, -2
// realpdgs = -1, 3, -13, 14, -2, 3
TEST(QCDSoftCollinearISR1, Wj_born_4_real_16) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.07088666960957468;
    ps.X2 = 0.2015495041129967;
    ps.Jacobian = 165.3528262348746;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(776.938262550079, 0, 0, 776.938262550079);
    ps.Momenta[1].Set(776.938262550079, 0, 0, -776.938262550079);
    ps.Momenta[2].Set(685.840406864241, 184.3115151187386, -29.47833590972903,
                      -659.9525413197041);
    ps.Momenta[3].Set(656.8697704427642, -169.043340861163, 228.847926403711,
                      592.0564760331129);
    ps.Momenta[4].Set(211.166347793153, -15.26817425757563, -199.369590493982,
                      67.89606528659121);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.838729511685992e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(9.452759604458022e-10, 0);
    spin[1][2] =
        std::complex<double>(4.618519736957511e-10, -7.945142655456995e-10);
    spin[2][1] =
        std::complex<double>(4.618519736957511e-10, 7.945142655456995e-10);
    spin[2][2] = std::complex<double>(8.934535512401895e-10, 0);

    // radiation variables
    double phi = 1.727539840592873;
    int bornpdgs[] = {-1, 0, -13, 14, -2};
    int realpdgs[] = {-1, 3, -13, 14, -2, 3};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 3.240981720732423e-09
    double fks_g = 5.595565385487124e-27;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, 1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 0, -13, 14, -2
// realpdgs = -1, 3, -13, 14, -2, 3
TEST(QCDSoftLimit, Wj_born_4_real_16_0) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.3574614665814302;
    ps.X2 = 0.8806686282421623;
    ps.Jacobian = 2113.777703789026;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3646.987174456651, 0, 0, 3646.987174456651);
    ps.Momenta[1].Set(3646.987174456651, 0, 0, -3646.987174456651);
    ps.Momenta[2].Set(3640.46887330842, -2508.782578795202, 2216.29487405096,
                      1430.75526248084);
    ps.Momenta[3].Set(3078.430230737571, 2114.973493688854, -1937.54671740733,
                      -1117.824907745673);
    ps.Momenta[4].Set(575.0752448673113, 393.8090851063478, -278.7481566436301,
                      -312.9303547351673);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 4.963432052467984e-09;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, 7.445148078701976e-09);
    ColorCorr.Set(0, 4, -8.272386754113307e-10);
    ColorCorr.Set(1, 0, 7.445148078701976e-09);
    ColorCorr.Set(1, 4, 7.445148078701976e-09);
    ColorCorr.Set(4, 0, -8.272386754113307e-10);
    ColorCorr.Set(4, 1, 7.445148078701976e-09);

    // radiation variables
    double y = -0.3412862061776636;
    double phi = 0.1029211653286697;
    int bornpdgs[] = {-1, 0, -13, 14, -2};
    int realpdgs[] = {-1, 3, -13, 14, -2, 3};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 6.485052386642057e-07
    double fks_g = 1.709226887126472e-24;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 0, alphas, ColorCorr, y, phi);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 0, -13, 14, -2
// realpdgs = -1, 3, -13, 14, -2, 3
TEST(QCDCollinearISR2, Wj_born_4_real_16) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.01708254629926031;
    ps.X2 = 0.5935527704281682;
    ps.Jacobian = 425.0986063430867;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(654.5145841079453, 0, 0, 654.5145841079453);
    ps.Momenta[1].Set(654.5145841079453, 0, 0, -654.5145841079453);
    ps.Momenta[2].Set(537.192724476805, 51.98442703908609, -27.77719014174031,
                      533.9495016235712);
    ps.Momenta[3].Set(127.4150505554119, -61.76693813638536, 15.98769945021262,
                      110.2897725429869);
    ps.Momenta[4].Set(644.4213931836741, 9.782511097299256, 11.78949069152768,
                      -644.2392741665577);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 0.006958169042061263;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(0.003508921325349249, -9.03501810404587e-21);
    spin[1][2] =
        std::complex<double>(-0.0001089943182086368, -0.003477248784338162);
    spin[2][1] =
        std::complex<double>(-0.0001089943182086368, 0.003477248784338162);
    spin[2][2] =
        std::complex<double>(0.003449247716712015, 9.03501810404587e-21);

    // radiation variables
    double phi = 3.48560265689846;
    double xi = 0.3993053900819165;
    int bornpdgs[] = {-1, 0, -13, 14, -2};
    int realpdgs[] = {-1, 3, -13, 14, -2, 3};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 15.41932911781844
    double fks_g = 4.917063523943525e-08;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, -1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 4.90490955235954e-08
    EXPECT_NEAR(limit / fks_g, 1.0, 5e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 0, -13, 14, -2
// realpdgs = -1, 3, -13, 14, -2, 3
TEST(QCDSoftCollinearISR2, Wj_born_4_real_16) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.2482573761847426;
    ps.X2 = 0.1848097008997911;
    ps.Jacobian = 1723.654316771208;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1392.280752324122, 0, 0, 1392.280752324122);
    ps.Momenta[1].Set(1392.280752324122, 0, 0, -1392.280752324122);
    ps.Momenta[2].Set(544.7781154554863, 22.36666786586761, -359.1428310905015,
                      -409.0224371890966);
    ps.Momenta[3].Set(1011.431074582938, -770.3441249322514, 332.2279977125811,
                      -564.9666409180712);
    ps.Momenta[4].Set(1228.352314609821, 747.9774570663836, 26.91483337792046,
                      973.9890781071676);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.249167030885816e-08;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(8.116767492937793e-09, 0);
    spin[1][2] =
        std::complex<double>(4.521589444878681e-09, 3.881404108436301e-09);
    spin[2][1] =
        std::complex<double>(4.521589444878681e-09, -3.881404108436301e-09);
    spin[2][2] = std::complex<double>(4.374902815920371e-09, 0);

    // radiation variables
    double phi = 2.009113931519106;
    int bornpdgs[] = {-1, 0, -13, 14, -2};
    int realpdgs[] = {-1, 3, -13, 14, -2, 3};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 0.01563350783078184
    double fks_g = 2.077803516937306e-20;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, -1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 0, -13, 14, -2
// realpdgs = -1, 4, -13, 14, -2, 4
TEST(QCDCollinearISR1, Wj_born_4_real_16_n1) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.2503840090096894;
    ps.X2 = 0.5325800109679975;
    ps.Jacobian = 1857.589115936028;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2373.608465328474, 0, 0, 2373.608465328474);
    ps.Momenta[1].Set(2373.608465328474, 0, 0, -2373.608465328474);
    ps.Momenta[2].Set(2013.777632964636, 1262.651255473962, 1323.469796576368,
                      -842.2825295759058);
    ps.Momenta[3].Set(1956.941600743892, -803.2017909582448, -1743.908886899531,
                      378.5090565843478);
    ps.Momenta[4].Set(776.4976969484196, -459.4494645157175, 420.439090323163,
                      463.773472991558);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 8.033580963499765e-10;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(3.561469371086886e-10, 0);
    spin[1][2] =
        std::complex<double>(-3.285425736034775e-10, 2.265671245671231e-10);
    spin[2][1] =
        std::complex<double>(-3.285425736034775e-10, -2.265671245671231e-10);
    spin[2][2] = std::complex<double>(4.472111592412879e-10, 0);

    // radiation variables
    double phi = 0.9393965695287076;
    double xi = 0.5487762255863173;
    int bornpdgs[] = {-1, 0, -13, 14, -2};
    int realpdgs[] = {-1, 4, -13, 14, -2, 4};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 7.552242067673308e-17
    double fks_g = 4.548796140044845e-25;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, 1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 0, -13, 14, -2
// realpdgs = -1, 4, -13, 14, -2, 4
TEST(QCDSoftCollinearISR1, Wj_born_4_real_16_n1) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.3962293511973947;
    ps.X2 = 0.5375553430555584;
    ps.Jacobian = 7111.812603158703;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2999.841229680325, 0, 0, 2999.841229680325);
    ps.Momenta[1].Set(2999.841229680325, 0, 0, -2999.841229680325);
    ps.Momenta[2].Set(915.1071676765896, -504.5388693013035, -737.967456655378,
                      -195.6161818838725);
    ps.Momenta[3].Set(2732.33529081729, -649.9177071217562, 1805.565858423489,
                      1945.043713213318);
    ps.Momenta[4].Set(2352.240000866772, 1154.456576423059, -1067.598401768111,
                      -1749.427531329445);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 7.263962257990994e-10;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(1.517238712240011e-10, -1.723292942818808e-26);
    spin[1][2] =
        std::complex<double>(-2.632686803139256e-10, 1.337202912325193e-10);
    spin[2][1] =
        std::complex<double>(-2.632686803139256e-10, -1.337202912325193e-10);
    spin[2][2] =
        std::complex<double>(5.746723545750983e-10, 1.723292942818808e-26);

    // radiation variables
    double phi = 1.312391711225357;
    int bornpdgs[] = {-1, 0, -13, 14, -2};
    int realpdgs[] = {-1, 4, -13, 14, -2, 4};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 1.321617560027602e-10
    double fks_g = 9.635624687307161e-29;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, 1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 0, -13, 14, -2
// realpdgs = -1, 4, -13, 14, -2, 4
TEST(QCDSoftLimit, Wj_born_4_real_16_0_n1) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.6984455399068246;
    ps.X2 = 0.007770514181707711;
    ps.Jacobian = 120.764910659211;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(478.8555326078015, 0, 0, 478.8555326078015);
    ps.Momenta[1].Set(478.8555326078015, 0, 0, -478.8555326078015);
    ps.Momenta[2].Set(450.6574061548927, 278.5428787814178, 42.32358117651219,
                      351.7309723052712);
    ps.Momenta[3].Set(256.8256898673023, -76.45102388666305, 53.26971431581921,
                      -239.3261654295401);
    ps.Momenta[4].Set(250.2279691934081, -202.0918548947547, -95.59329549233141,
                      -112.404806875731);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.777430243222923e-07;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, 2.666145364834385e-07);
    ColorCorr.Set(0, 4, -2.962383738704872e-08);
    ColorCorr.Set(1, 0, 2.666145364834385e-07);
    ColorCorr.Set(1, 4, 2.666145364834385e-07);
    ColorCorr.Set(4, 0, -2.962383738704872e-08);
    ColorCorr.Set(4, 1, 2.666145364834385e-07);

    // radiation variables
    double y = -0.7028647277729831;
    double phi = 4.887905443416265;
    int bornpdgs[] = {-1, 0, -13, 14, -2};
    int realpdgs[] = {-1, 4, -13, 14, -2, 4};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 0.0006123563932820063
    double fks_g = 2.19811322108183e-20;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 0, alphas, ColorCorr, y, phi);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 0, -13, 14, -2
// realpdgs = -1, 4, -13, 14, -2, 4
TEST(QCDCollinearISR2, Wj_born_4_real_16_n1) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.1018424091801684;
    ps.X2 = 0.1537450191595156;
    ps.Jacobian = 117.0071606939764;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(813.3513958402144, 0, 0, 813.3513958402144);
    ps.Momenta[1].Set(813.3513958402144, 0, 0, -813.3513958402144);
    ps.Momenta[2].Set(676.0027914429448, 7.183226247720075, 569.9708893131708,
                      -363.4024774748801);
    ps.Momenta[3].Set(807.963902000644, -11.12134511299776, -647.4404707533845,
                      483.221294491064);
    ps.Momenta[4].Set(142.7360982368397, 3.938118865277686, 77.46958144021372,
                      -119.818817016184);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 2.937101990726912e-08;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(1.563998977253301e-08, -1.723292942818808e-26);
    spin[1][2] =
        std::complex<double>(-1.912985006120946e-10, 1.465321042493006e-08);
    spin[2][1] =
        std::complex<double>(-1.912985006120946e-10, -1.465321042493006e-08);
    spin[2][2] =
        std::complex<double>(1.373103013473611e-08, 1.723292942818808e-26);

    // radiation variables
    double phi = 5.003269030179543;
    double xi = 0.2856744550309387;
    int bornpdgs[] = {-1, 0, -13, 14, -2};
    int realpdgs[] = {-1, 4, -13, 14, -2, 4};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 4.538301912134019e-05
    double fks_g = 7.40740677923335e-14;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, -1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 7.404854675969405e-14
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 0, -13, 14, -2
// realpdgs = -1, 4, -13, 14, -2, 4
TEST(QCDSoftCollinearISR2, Wj_born_4_real_16_n1) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.4003268213069582;
    ps.X2 = 0.6529788837733932;
    ps.Jacobian = 9171.866839839195;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3323.305522959973, 0, 0, 3323.305522959973);
    ps.Momenta[1].Set(3323.305522959973, 0, 0, -3323.305522959973);
    ps.Momenta[2].Set(840.0238843835167, 247.6279615855601, 83.81506288102787,
                      -798.3079319473163);
    ps.Momenta[3].Set(3068.249068250299, -219.1148660410865, 2452.065862872343,
                      1831.260228486037);
    ps.Momenta[4].Set(2738.338093286131, -28.51309554447363, -2535.880925753371,
                      -1032.952296538721);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 4.576616029569742e-10;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(3.472787140166504e-10, 0);
    spin[1][2] =
        std::complex<double>(1.094901917763898e-10, 1.623130482290558e-10);
    spin[2][1] =
        std::complex<double>(1.094901917763898e-10, -1.623130482290558e-10);
    spin[2][2] = std::complex<double>(1.103828889403238e-10, 0);

    // radiation variables
    double phi = 3.503285751322375;
    int bornpdgs[] = {-1, 0, -13, 14, -2};
    int realpdgs[] = {-1, 4, -13, 14, -2, 4};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 0.0002361198131164679
    double fks_g = 5.68688404467481e-23;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, -1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 0, -13, 14, -2
// realpdgs = -1, 5, -13, 14, -2, 5
TEST(QCDCollinearISR1, Wj_born_4_real_16_n2) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.6606998005811029;
    ps.X2 = 0.6335701822068671;
    ps.Jacobian = 17426.69555938027;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(4205.453249159284, 0, 0, 4205.453249159284);
    ps.Momenta[1].Set(4205.453249159284, 0, 0, -4205.453249159284);
    ps.Momenta[2].Set(1664.546982540554, -620.2176450506051, 1499.638657299769,
                      370.3115274804348);
    ps.Momenta[3].Set(2634.844968361158, -505.0410257881304, 2416.331507121734,
                      -921.2402603337491);
    ps.Momenta[4].Set(4111.514547416859, 1125.258670838735, -3915.9701644215,
                      550.9287328533139);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 8.239266035568307e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(5.941696527466915e-09, 0);
    spin[1][2] =
        std::complex<double>(3.659638503326439e-09, 5.084356329669805e-10);
    spin[2][1] =
        std::complex<double>(3.659638503326439e-09, -5.084356329669805e-10);
    spin[2][2] = std::complex<double>(2.297569508101393e-09, 0);

    // radiation variables
    double phi = 1.382572116274675;
    double xi = 0.1278783945333864;
    int bornpdgs[] = {-1, 0, -13, 14, -2};
    int realpdgs[] = {-1, 5, -13, 14, -2, 5};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 2.824293091469502e-15
    double fks_g = 9.237067336994935e-25;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, 1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 0, -13, 14, -2
// realpdgs = -1, 5, -13, 14, -2, 5
TEST(QCDSoftCollinearISR1, Wj_born_4_real_16_n2) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.8384620424809519;
    ps.X2 = 0.2087746608756262;
    ps.Jacobian = 5343.202501986206;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2719.530622614485, 0, 0, 2719.530622614485);
    ps.Momenta[1].Set(2719.530622614485, 0, 0, -2719.530622614485);
    ps.Momenta[2].Set(976.9129854401253, -265.0733376849438, -524.0289196208416,
                      780.6976355613714);
    ps.Momenta[3].Set(2512.719913496575, 2184.846438186455, 832.6884386730017,
                      -920.2376700214537);
    ps.Momenta[4].Set(1949.428346292271, -1919.773100501511, -308.6595190521601,
                      139.5400344600822);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 2.829930962744924e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(1.146880457146934e-09, -1.723292942818808e-26);
    spin[1][2] =
        std::complex<double>(2.026934291591542e-10, 1.374471937522663e-09);
    spin[2][1] =
        std::complex<double>(2.026934291591542e-10, -1.374471937522663e-09);
    spin[2][2] =
        std::complex<double>(1.683050505597991e-09, 1.723292942818808e-26);

    // radiation variables
    double phi = 6.271179221625208;
    int bornpdgs[] = {-1, 0, -13, 14, -2};
    int realpdgs[] = {-1, 5, -13, 14, -2, 5};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 2.341605491135885e-09
    double fks_g = 1.222061550989107e-28;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, 1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 0, -13, 14, -2
// realpdgs = -1, 5, -13, 14, -2, 5
TEST(QCDSoftLimit, Wj_born_4_real_16_0_n2) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.02840760192097846;
    ps.X2 = 0.9498566141527931;
    ps.Jacobian = 241.7260455524285;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1067.725633003337, 0, 0, 1067.725633003337);
    ps.Momenta[1].Set(1067.725633003337, 0, 0, -1067.725633003337);
    ps.Momenta[2].Set(899.8927290493507, 698.6921785020274, 499.5937370437001,
                      -268.4068952219545);
    ps.Momenta[3].Set(1010.930803612724, -897.952330175639, -414.9732393171696,
                      208.4704129506692);
    ps.Momenta[4].Set(224.6277333445996, 199.2601516736116, -84.62049772653057,
                      59.93648227128537);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 6.77678980955308e-09;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, 1.016518471432962e-08);
    ColorCorr.Set(0, 4, -1.129464968258847e-09);
    ColorCorr.Set(1, 0, 1.016518471432962e-08);
    ColorCorr.Set(1, 4, 1.016518471432962e-08);
    ColorCorr.Set(4, 0, -1.129464968258847e-09);
    ColorCorr.Set(4, 1, 1.016518471432962e-08);

    // radiation variables
    double y = -0.2970789400610485;
    double phi = 5.015425835321595;
    int bornpdgs[] = {-1, 0, -13, 14, -2};
    int realpdgs[] = {-1, 5, -13, 14, -2, 5};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 2.192269158689096e-05
    double fks_g = 1.162617691480065e-23;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 0, alphas, ColorCorr, y, phi);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 0, -13, 14, -2
// realpdgs = -1, 5, -13, 14, -2, 5
TEST(QCDCollinearISR2, Wj_born_4_real_16_n2) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.6220706643260989;
    ps.X2 = 0.5669438493993724;
    ps.Jacobian = 911.695417946175;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3860.141647605607, 0, 0, 3860.141647605607);
    ps.Momenta[1].Set(3860.141647605607, 0, 0, -3860.141647605607);
    ps.Momenta[2].Set(3842.585240275701, 781.6691305162797, -1920.276945458041,
                      -3235.272964051862);
    ps.Momenta[3].Set(3643.358215975167, -633.2948204364473, 1873.142333792818,
                      3060.120023410081);
    ps.Momenta[4].Set(234.3398389603459, -148.3743100798324, 47.13461166522322,
                      175.1529406417808);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 6.713909795486522e-11;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(5.547313215514883e-11, 0);
    spin[1][2] =
        std::complex<double>(2.352824330975261e-11, 9.673129239403131e-12);
    spin[2][1] =
        std::complex<double>(2.352824330975261e-11, -9.673129239403131e-12);
    spin[2][2] = std::complex<double>(1.16659657997164e-11, 0);

    // radiation variables
    double phi = 2.384907786246244;
    double xi = 0.280291594789909;
    int bornpdgs[] = {-1, 0, -13, 14, -2};
    int realpdgs[] = {-1, 5, -13, 14, -2, 5};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 3.123120435597224e-09
    double fks_g = 4.907257830706802e-18;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, -1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 4.907817717048684e-18
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 0, -13, 14, -2
// realpdgs = -1, 5, -13, 14, -2, 5
TEST(QCDSoftCollinearISR2, Wj_born_4_real_16_n2) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.3197095031516604;
    ps.X2 = 0.3209245791181417;
    ps.Jacobian = 4073.66176100953;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2082.057022387585, 0, 0, 2082.057022387585);
    ps.Momenta[1].Set(2082.057022387585, 0, 0, -2082.057022387585);
    ps.Momenta[2].Set(432.738191813068, -108.4172822265843, -392.6457045755033,
                      146.07322220945);
    ps.Momenta[3].Set(1790.079032663112, -633.7300490965694, 307.1720336903217,
                      1645.726134498299);
    ps.Momenta[4].Set(1941.296820298991, 742.1473313231538, 85.47367088518165,
                      -1791.799356707749);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.127723497232436e-08;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(8.03842955542751e-09, 0);
    spin[1][2] =
        std::complex<double>(-1.069086476072775e-09, -4.989184632195993e-09);
    spin[2][1] =
        std::complex<double>(-1.069086476072775e-09, 4.989184632195993e-09);
    spin[2][2] = std::complex<double>(3.238805416896852e-09, 0);

    // radiation variables
    double phi = 0.8545990265900053;
    int bornpdgs[] = {-1, 0, -13, 14, -2};
    int realpdgs[] = {-1, 5, -13, 14, -2, 5};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 0.007573351566471241
    double fks_g = 6.984803715275576e-21;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, -1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 0, -13, 14, -2
// realpdgs = -1, -5, -13, 14, -2, -5
TEST(QCDCollinearISR1, Wj_born_4_real_8) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.2629099199970391;
    ps.X2 = 0.5514966799636554;
    ps.Jacobian = 1061.526466456173;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2475.074605609459, 0, 0, 2475.074605609459);
    ps.Momenta[1].Set(2475.074605609459, 0, 0, -2475.074605609459);
    ps.Momenta[2].Set(2311.074025080603, -786.2833585683321, 2117.179458598994,
                      -490.2782572450422);
    ps.Momenta[3].Set(2213.533481686346, 1112.397218916618, -1844.512235998514,
                      509.9779535870215);
    ps.Momenta[4].Set(425.5417044519695, -326.113860348286, -272.66722260048,
                      -19.69969634197933);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.75350596162617e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(1.361590095807601e-09, 1.723292942818808e-26);
    spin[1][2] =
        std::complex<double>(-8.382561979150877e-11, 7.256734987272122e-10);
    spin[2][1] =
        std::complex<double>(-8.382561979150877e-11, -7.256734987272122e-10);
    spin[2][2] =
        std::complex<double>(3.919158658185687e-10, -1.723292942818808e-26);

    // radiation variables
    double phi = 5.220242399711965;
    double xi = 0.1609215230016346;
    int bornpdgs[] = {-1, 0, -13, 14, -2};
    int realpdgs[] = {-1, -5, -13, 14, -2, -5};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 3.233667419790289e-16
    double fks_g = 1.674763991963441e-25;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, 1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 0, -13, 14, -2
// realpdgs = -1, -5, -13, 14, -2, -5
TEST(QCDSoftCollinearISR1, Wj_born_4_real_8) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.05526406740080247;
    ps.X2 = 0.573846325595099;
    ps.Jacobian = 1345.118668797483;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1157.530870063621, 0, 0, 1157.530870063621);
    ps.Momenta[1].Set(1157.530870063621, 0, 0, -1157.530870063621);
    ps.Momenta[2].Set(225.6528517718633, 64.80652829977352, -20.88503363329276,
                      -215.1351639612673);
    ps.Momenta[3].Set(936.4133820420781, 108.1463370208141, 159.776751265181,
                      -916.3218766413949);
    ps.Momenta[4].Set(1152.995506313302, -172.9528653205874, -138.8917176318881,
                      1131.457040602661);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.534288605521144e-06;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(8.938490875617727e-07, 0);
    spin[1][2] =
        std::complex<double>(3.224850181680811e-07, 6.844411529296128e-07);
    spin[2][1] =
        std::complex<double>(3.224850181680811e-07, -6.844411529296128e-07);
    spin[2][2] = std::complex<double>(6.404395179593709e-07, 0);

    // radiation variables
    double phi = 3.040630653315087;
    int bornpdgs[] = {-1, 0, -13, 14, -2};
    int realpdgs[] = {-1, -5, -13, 14, -2, -5};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 1.19820574894323e-06
    double fks_g = 2.138859585435555e-24;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, 1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 0, -13, 14, -2
// realpdgs = -1, -5, -13, 14, -2, -5
TEST(QCDSoftLimit, Wj_born_4_real_8_0) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.9125537839328608;
    ps.X2 = 0.9618633464900554;
    ps.Jacobian = 36224.30986946823;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(6089.747411894938, 0, 0, 6089.747411894938);
    ps.Momenta[1].Set(6089.747411894938, 0, 0, -6089.747411894938);
    ps.Momenta[2].Set(6061.285743518443, -6035.390572285462, 536.1576039538439,
                      160.5631597539256);
    ps.Momenta[3].Set(216.1923382475237, 172.4847389446682, 94.15092112909588,
                      90.13182567031282);
    ps.Momenta[4].Set(5902.016742023913, 5862.905833340791, -630.3085250829397,
                      -250.6949854242384);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 9.890873599975163e-09;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, 1.483631039996274e-08);
    ColorCorr.Set(0, 4, -1.648478933329194e-09);
    ColorCorr.Set(1, 0, 1.483631039996274e-08);
    ColorCorr.Set(1, 4, 1.483631039996274e-08);
    ColorCorr.Set(4, 0, -1.648478933329194e-09);
    ColorCorr.Set(4, 1, 1.483631039996274e-08);

    // radiation variables
    double y = 0.2286408413598551;
    double phi = 4.36928089763517;
    int bornpdgs[] = {-1, 0, -13, 14, -2};
    int realpdgs[] = {-1, -5, -13, 14, -2, -5};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 4.476267175551617e-07
    double fks_g = 3.899859675388004e-25;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 0, alphas, ColorCorr, y, phi);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 0, -13, 14, -2
// realpdgs = -1, -5, -13, 14, -2, -5
TEST(QCDCollinearISR2, Wj_born_4_real_8) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.2908913616777404;
    ps.X2 = 0.9895209682467083;
    ps.Jacobian = 9835.500123126478;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3487.31573759358, 0, 0, 3487.31573759358);
    ps.Momenta[1].Set(3487.31573759358, 0, 0, -3487.31573759358);
    ps.Momenta[2].Set(2029.95152485244, -539.1921376260173, 482.4713615119043,
                      -1896.627643290554);
    ps.Momenta[3].Set(2146.31227911517, 1095.908831245091, 1805.052553179611,
                      383.9603017752481);
    ps.Momenta[4].Set(2798.367671219551, -556.7166936190737, -2287.523914691516,
                      1512.667341515305);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 8.0280161430429e-10;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(4.020745799126404e-10, 3.446585885637615e-26);
    spin[1][2] =
        std::complex<double>(-3.68436861221411e-10, -1.592998220525543e-10);
    spin[2][1] =
        std::complex<double>(-3.68436861221411e-10, 1.592998220525543e-10);
    spin[2][2] =
        std::complex<double>(4.007270343916497e-10, -3.446585885637615e-26);

    // radiation variables
    double phi = 5.284260356386956;
    double xi = 0.008074354031389298;
    int bornpdgs[] = {-1, 0, -13, 14, -2};
    int realpdgs[] = {-1, -5, -13, 14, -2, -5};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 1.651627497937747e-06
    double fks_g = 2.153563469509434e-18;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, -1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 2.153565735078309e-18
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 0, -13, 14, -2
// realpdgs = -1, -5, -13, 14, -2, -5
TEST(QCDSoftCollinearISR2, Wj_born_4_real_8) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.838025506219477;
    ps.X2 = 0.3193990552600141;
    ps.Jacobian = 4994.148048594314;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3362.860010094831, 0, 0, 3362.860010094831);
    ps.Momenta[1].Set(3362.860010094831, 0, 0, -3362.860010094831);
    ps.Momenta[2].Set(3174.566514126513, 1593.6690157239, -1216.478236505798,
                      -2461.355789202918);
    ps.Momenta[3].Set(2077.646448739257, -562.9791544016711, 279.7594780724168,
                      1980.253486828561);
    ps.Momenta[4].Set(1473.507057323893, -1030.689861322229, 936.7187584333815,
                      481.1023023743567);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 4.456984906849295e-10;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(2.134118475510764e-10, 0);
    spin[1][2] =
        std::complex<double>(1.397510547433242e-11, -2.222103015155113e-10);
    spin[2][1] =
        std::complex<double>(1.397510547433242e-11, 2.222103015155113e-10);
    spin[2][2] = std::complex<double>(2.322866431338531e-10, 0);

    // radiation variables
    double phi = 1.794990597140928;
    int bornpdgs[] = {-1, 0, -13, 14, -2};
    int realpdgs[] = {-1, -5, -13, 14, -2, -5};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 0.0001144713901400984
    double fks_g = 1.060503526861272e-22;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, -1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 0, -13, 14, -2
// realpdgs = -1, -4, -13, 14, -2, -4
TEST(QCDCollinearISR1, Wj_born_4_real_8_n1) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.9763541688829811;
    ps.X2 = 0.3467596375601261;
    ps.Jacobian = 9048.971303868997;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3782.085297714031, 0, 0, 3782.085297714031);
    ps.Momenta[1].Set(3782.085297714031, 0, 0, -3782.085297714031);
    ps.Momenta[2].Set(2484.257567821626, 1100.212608075974, 1190.323260168175,
                      1882.604158237277);
    ps.Momenta[3].Set(2705.985327228526, -1542.906995568247, -2212.724916768025,
                      213.6423106890941);
    ps.Momenta[4].Set(2373.927700377911, 442.6943874922721, 1022.40165659985,
                      -2096.246468926372);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 5.454556556199438e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(4.839489316212613e-09, 1.378634354255046e-25);
    spin[1][2] =
        std::complex<double>(-1.678909624549946e-09, 3.973333729529585e-10);
    spin[2][1] =
        std::complex<double>(-1.678909624549946e-09, -3.973333729529585e-10);
    spin[2][2] =
        std::complex<double>(6.150672399868249e-10, -1.378634354255046e-25);

    // radiation variables
    double phi = 2.230643492460946;
    double xi = 0.003200589121881023;
    int bornpdgs[] = {-1, 0, -13, 14, -2};
    int realpdgs[] = {-1, -4, -13, 14, -2, -4};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 1.076132176420803e-13
    double fks_g = 2.20473025630217e-26;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, 1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 0, -13, 14, -2
// realpdgs = -1, -4, -13, 14, -2, -4
TEST(QCDSoftCollinearISR1, Wj_born_4_real_8_n1) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.02636889638456452;
    ps.X2 = 0.7341231014770599;
    ps.Jacobian = 4.857341971806851;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(904.365067794178, 0, 0, 904.365067794178);
    ps.Momenta[1].Set(904.365067794178, 0, 0, -904.365067794178);
    ps.Momenta[2].Set(900.3948284586326, -87.9631599369585, -802.6725474686932,
                      -398.3592739589146);
    ps.Momenta[3].Set(903.0061994241013, 84.88804489004616, 806.7198424939893,
                      396.7585055919761);
    ps.Momenta[4].Set(5.32910770562231, 3.075115046912355, -4.047295025295973,
                      1.600768366938445);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.419889722934891e-07;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(7.134756770817116e-08, 0);
    spin[1][2] =
        std::complex<double>(-6.755212996754549e-10, 7.099039419262984e-08);
    spin[2][1] =
        std::complex<double>(-6.755212996754549e-10, -7.099039419262984e-08);
    spin[2][2] = std::complex<double>(7.064140458531798e-08, 0);

    // radiation variables
    double phi = 4.260821062617159;
    int bornpdgs[] = {-1, 0, -13, 14, -2};
    int realpdgs[] = {-1, -4, -13, 14, -2, -4};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 1.761785779506323e-07
    double fks_g = 3.340196221673558e-25;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, 1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 0, -13, 14, -2
// realpdgs = -1, -4, -13, 14, -2, -4
TEST(QCDSoftLimit, Wj_born_4_real_8_0_n1) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.559581307137524;
    ps.X2 = 0.2219330690607677;
    ps.Jacobian = 3320.055091984896;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2290.635385273553, 0, 0, 2290.635385273553);
    ps.Momenta[1].Set(2290.635385273553, 0, 0, -2290.635385273553);
    ps.Momenta[2].Set(1910.749804969382, -1358.788549468031, -1343.249609280253,
                      -18.41146907609334);
    ps.Momenta[3].Set(1232.421511597577, 648.0707504395348, 486.2532715581796,
                      928.6683156930476);
    ps.Momenta[4].Set(1438.099453980148, 710.7177990284963, 856.9963377220741,
                      -910.2568466169544);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 5.4529219276233e-09;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, 8.17938289143495e-09);
    ColorCorr.Set(0, 4, -9.0882032127055e-10);
    ColorCorr.Set(1, 0, 8.17938289143495e-09);
    ColorCorr.Set(1, 4, 8.17938289143495e-09);
    ColorCorr.Set(4, 0, -9.0882032127055e-10);
    ColorCorr.Set(4, 1, 8.17938289143495e-09);

    // radiation variables
    double y = -0.2651692436507673;
    double phi = 1.892305427654365;
    int bornpdgs[] = {-1, 0, -13, 14, -2};
    int realpdgs[] = {-1, -4, -13, 14, -2, -4};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 3.597751994775167e-07
    double fks_g = 2.020347220433514e-23;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 0, alphas, ColorCorr, y, phi);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 0, -13, 14, -2
// realpdgs = -1, -4, -13, 14, -2, -4
TEST(QCDCollinearISR2, Wj_born_4_real_8_n1) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.9460001630368104;
    ps.X2 = 0.4395103956514603;
    ps.Jacobian = 2226.771202977191;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(4191.249727238543, 0, 0, 4191.249727238543);
    ps.Momenta[1].Set(4191.249727238543, 0, 0, -4191.249727238543);
    ps.Momenta[2].Set(4001.034971167758, 3815.58676493136, -855.5793294452195,
                      847.0906036629646);
    ps.Momenta[3].Set(3854.317586195373, -3779.621071482395, 665.5824035113318,
                      -356.6912886478632);
    ps.Momenta[4].Set(527.146897113957, -35.96569344896517, 189.9969259338878,
                      -490.3993150151014);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 2.1624133269823e-08;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(1.344164304280689e-08, 0);
    spin[1][2] =
        std::complex<double>(2.114031683383543e-09, 1.027214597221325e-08);
    spin[2][1] =
        std::complex<double>(2.114031683383543e-09, -1.027214597221325e-08);
    spin[2][2] = std::complex<double>(8.182490227016108e-09, 0);

    // radiation variables
    double phi = 4.203379112539103;
    double xi = 0.544844702171661;
    int bornpdgs[] = {-1, 0, -13, 14, -2};
    int realpdgs[] = {-1, -4, -13, 14, -2, -4};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 1.315657071431018e-06
    double fks_g = 7.811207315973956e-15;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, -1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 7.802212466938883e-15
    EXPECT_NEAR(limit / fks_g, 1.0, 5e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 0, -13, 14, -2
// realpdgs = -1, -4, -13, 14, -2, -4
TEST(QCDSoftCollinearISR2, Wj_born_4_real_8_n1) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.9591322145283279;
    ps.X2 = 0.03121857552301321;
    ps.Jacobian = 971.0765521181272;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1124.758119487416, 0, 0, 1124.758119487416);
    ps.Momenta[1].Set(1124.758119487416, 0, 0, -1124.758119487416);
    ps.Momenta[2].Set(575.0603100282474, 549.8463839780482, -78.35769466987443,
                      149.0751014826806);
    ps.Momenta[3].Set(817.8246520580607, -5.253588668448572, 598.6384975610092,
                      -557.1727833965764);
    ps.Momenta[4].Set(856.6312768885235, -544.5927953095996, -520.2808028911346,
                      408.0976819138958);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.502201088791173e-08;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(1.953037853600106e-09, -1.378634354255046e-25);
    spin[1][2] =
        std::complex<double>(6.159114089984825e-10, 5.014464296408456e-09);
    spin[2][1] =
        std::complex<double>(6.159114089984825e-10, -5.014464296408456e-09);
    spin[2][2] =
        std::complex<double>(1.306897303431163e-08, 1.378634354255046e-25);

    // radiation variables
    double phi = 5.209321319619472;
    int bornpdgs[] = {-1, 0, -13, 14, -2};
    int realpdgs[] = {-1, -4, -13, 14, -2, -4};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 0.02423677777944543
    double fks_g = 4.549424786153564e-20;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, -1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 0, -13, 14, -2
// realpdgs = -1, -3, -13, 14, -2, -3
TEST(QCDCollinearISR1, Wj_born_4_real_8_n2) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.111445085731912;
    ps.X2 = 0.7519894639314444;
    ps.Jacobian = 1328.569305794359;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1881.697014457263, 0, 0, 1881.697014457263);
    ps.Momenta[1].Set(1881.697014457263, 0, 0, -1881.697014457263);
    ps.Momenta[2].Set(1256.03722948582, 1056.665183860886, -360.5957095590315,
                      -575.3772200202862);
    ps.Momenta[3].Set(1806.814900614327, -1384.473069010866, 977.6230721549452,
                      626.1529645828713);
    ps.Momenta[4].Set(700.5418988143797, 327.8078851499798, -617.0273625959138,
                      -50.77574456258507);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.467672763388727e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(5.623600476785931e-10, 0);
    spin[1][2] =
        std::complex<double>(-6.454804594441652e-11, 7.105949983890577e-10);
    spin[2][1] =
        std::complex<double>(-6.454804594441652e-11, -7.105949983890577e-10);
    spin[2][2] = std::complex<double>(9.053127157101336e-10, 0);

    // radiation variables
    double phi = 3.747100176432543;
    double xi = 0.8474051950249027;
    int bornpdgs[] = {-1, 0, -13, 14, -2};
    int realpdgs[] = {-1, -3, -13, 14, -2, -3};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 5.939150192764152e-17
    double fks_g = 8.529754816812433e-25;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, 1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 0, -13, 14, -2
// realpdgs = -1, -3, -13, 14, -2, -3
TEST(QCDSoftCollinearISR1, Wj_born_4_real_8_n2) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.3197132183525646;
    ps.X2 = 0.754757668122906;
    ps.Jacobian = 6801.758397386959;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3192.988981059238, 0, 0, 3192.988981059238);
    ps.Momenta[1].Set(3192.988981059238, 0, 0, -3192.988981059238);
    ps.Momenta[2].Set(2541.928657340648, 2031.082648621566, 1298.541700702987,
                      806.1600492541357);
    ps.Momenta[3].Set(1730.446492787214, -85.09280562949436, -1725.438679697936,
                      -100.3276703203643);
    ps.Momenta[4].Set(2113.602811990613, -1945.989842992071, 426.8969789949492,
                      -705.8323789337713);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 2.163289991408223e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(4.332703283068237e-10, 0);
    spin[1][2] =
        std::complex<double>(-7.253169207154113e-10, -4.727383546245399e-10);
    spin[2][1] =
        std::complex<double>(-7.253169207154113e-10, 4.727383546245399e-10);
    spin[2][2] = std::complex<double>(1.730019663101399e-09, 0);

    // radiation variables
    double phi = 1.360192760162521;
    int bornpdgs[] = {-1, 0, -13, 14, -2};
    int realpdgs[] = {-1, -3, -13, 14, -2, -3};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 3.08338657618497e-10
    double fks_g = 2.853922052780731e-28;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, 1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 0, -13, 14, -2
// realpdgs = -1, -3, -13, 14, -2, -3
TEST(QCDSoftLimit, Wj_born_4_real_8_0_n2) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.5469684255747289;
    ps.X2 = 0.4995550053430335;
    ps.Jacobian = 220.9239579448127;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3397.70870199745, 0, 0, 3397.70870199745);
    ps.Momenta[1].Set(3397.70870199745, 0, 0, -3397.70870199745);
    ps.Momenta[2].Set(3379.133112928631, -613.3338002792184, 1478.642895871856,
                      -2975.899432241198);
    ps.Momenta[3].Set(3351.769949515407, 611.8681218692482, -1414.149036288444,
                      2976.602375037549);
    ps.Momenta[4].Set(64.51434155086261, 1.465678409970122, -64.49385958341277,
                      -0.7029427963518189);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.708897514065822e-10;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, 2.563346271098734e-10);
    ColorCorr.Set(0, 4, -2.848162523443037e-11);
    ColorCorr.Set(1, 0, 2.563346271098734e-10);
    ColorCorr.Set(1, 4, 2.563346271098734e-10);
    ColorCorr.Set(4, 0, -2.848162523443037e-11);
    ColorCorr.Set(4, 1, 2.563346271098734e-10);

    // radiation variables
    double y = 0.723695446054812;
    double phi = 1.39028074096988;
    int bornpdgs[] = {-1, 0, -13, 14, -2};
    int realpdgs[] = {-1, -3, -13, 14, -2, -3};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 3.353387761392978e-09
    double fks_g = 4.094413392878634e-26;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 0, alphas, ColorCorr, y, phi);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 0, -13, 14, -2
// realpdgs = -1, -3, -13, 14, -2, -3
TEST(QCDCollinearISR2, Wj_born_4_real_8_n2) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.2053821103358793;
    ps.X2 = 0.3898885282692959;
    ps.Jacobian = 1190.590421685807;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1839.352179141951, 0, 0, 1839.352179141951);
    ps.Momenta[1].Set(1839.352179141951, 0, 0, -1839.352179141951);
    ps.Momenta[2].Set(1788.970128991037, -1710.733610543654, -29.44567026015767,
                      -522.434291257199);
    ps.Momenta[3].Set(1247.49461564931, 1174.620610040847, 299.7311589094484,
                      294.3984900098607);
    ps.Momenta[4].Set(642.2396136435553, 536.1130005028061, -270.2854886492907,
                      228.0358012473383);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.845123674008399e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(2.125291078453612e-11, 0);
    spin[1][2] =
        std::complex<double>(-3.810044893984298e-12, -1.968452340576895e-10);
    spin[2][1] =
        std::complex<double>(-3.810044893984298e-12, 1.968452340576895e-10);
    spin[2][2] = std::complex<double>(1.823870763223863e-09, 0);

    // radiation variables
    double phi = 1.403374297503039;
    double xi = 0.3721938036657526;
    int bornpdgs[] = {-1, 0, -13, 14, -2};
    int realpdgs[] = {-1, -3, -13, 14, -2, -3};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 8.418864286259971e-07
    double fks_g = 2.33250069283372e-15;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, -1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 2.332449090033888e-15
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 0, -13, 14, -2
// realpdgs = -1, -3, -13, 14, -2, -3
TEST(QCDSoftCollinearISR2, Wj_born_4_real_8_n2) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.5129015838501907;
    ps.X2 = 0.6596771672696562;
    ps.Jacobian = 819.5322390439873;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3780.907939995628, 0, 0, 3780.907939995628);
    ps.Momenta[1].Set(3780.907939995628, 0, 0, -3780.907939995628);
    ps.Momenta[2].Set(3663.136846441475, -2563.598577064426, 1197.602133735061,
                      -2326.431391789509);
    ps.Momenta[3].Set(3683.614133748998, 2607.636382387368, -1000.946113105552,
                      2401.53127447135);
    ps.Momenta[4].Set(215.064899800785, -44.03780532294231, -196.6560206295097,
                      -75.09988268184023);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 6.628729018072073e-10;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(3.546011938428524e-10, 0);
    spin[1][2] =
        std::complex<double>(9.557471522583516e-11, 3.165106467113483e-10);
    spin[2][1] =
        std::complex<double>(9.557471522583516e-11, -3.165106467113483e-10);
    spin[2][2] = std::complex<double>(3.08271707964355e-10, 0);

    // radiation variables
    double phi = 5.396494015841644;
    int bornpdgs[] = {-1, 0, -13, 14, -2};
    int realpdgs[] = {-1, -3, -13, 14, -2, -3};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 0.0002692600067741164
    double fks_g = 6.23712095488647e-23;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, -1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 0, -13, 14, -2
// realpdgs = -1, 0, -13, 14, -2, 0
TEST(QCDCollinearFSR, Wj_born_4_real_22) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.101039105598705;
    ps.X2 = 0.5309631434221416;
    ps.Jacobian = 325.3274194162649;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1505.533040887452, 0, 0, 1505.533040887452);
    ps.Momenta[1].Set(1505.533040887452, 0, 0, -1505.533040887452);
    ps.Momenta[2].Set(1488.695878832563, 1297.062374979603, 501.7098913178125,
                      -531.1608042948956);
    ps.Momenta[3].Set(1307.967639509255, -1145.454177020183, -504.6075532581789,
                      379.5856840625346);
    ps.Momenta[4].Set(214.4025634330849, -151.6081979594198, 2.897661940366407,
                      151.575120232361);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 2.677019096514131e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(5.581830837758498e-10, 0);
    spin[1][2] =
        std::complex<double>(-5.671045101694401e-10, 9.279498338547507e-10);
    spin[2][1] =
        std::complex<double>(-5.671045101694401e-10, -9.279498338547507e-10);
    spin[2][2] = std::complex<double>(2.118836012738281e-09, 0);

    // radiation variables
    double phi = 2.636642405893279;
    double xi = 0.1357439670858863;
    int bornpdgs[] = {-1, 0, -13, 14, -2};
    int realpdgs[] = {-1, 0, -13, 14, -2, 0};

    // FKS g = real matrix element * (1.0 - y) * xi * xi
    // real = 0.0002713212185383983
    double fks_g = 4.999480000957211e-14;
    double limit = FKS::QCD::CollinearLimitFSR(
        realpdgs[4], bornpdgs[4], ps, 4, xi, phi, alphas, born,
        spin); // limit = 4.999670142254868e-14
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 0, -13, 14, -2
// realpdgs = -1, 0, -13, 14, -2, 0
TEST(QCDSoftCollinearFSR, Wj_born_4_real_22) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.3260534059150153;
    ps.X2 = 0.9577103609031141;
    ps.Jacobian = 6815.773522351497;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(3632.242369868868, 0, 0, 3632.242369868868);
    ps.Momenta[1].Set(3632.242369868868, 0, 0, -3632.242369868868);
    ps.Momenta[2].Set(2862.198625327269, -2523.84492519382, -1322.836616082837,
                      269.2416229880135);
    ps.Momenta[3].Set(2540.456528346391, 1994.432828987972, 100.8467325009044,
                      -1570.346140064897);
    ps.Momenta[4].Set(1861.829586064075, 529.4120962058473, 1221.989883581933,
                      1301.104517076883);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 8.25037445370874e-10;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(2.647949713380562e-10, -1.723292942818808e-26);
    spin[1][2] =
        std::complex<double>(-3.082315290379702e-10, 2.309604173030451e-10);
    spin[2][1] =
        std::complex<double>(-3.082315290379702e-10, -2.309604173030451e-10);
    spin[2][2] =
        std::complex<double>(5.602424740328178e-10, 1.723292942818808e-26);

    // radiation variables
    double phi = 0.9434432790883504;
    int bornpdgs[] = {-1, 0, -13, 14, -2};
    int realpdgs[] = {-1, 0, -13, 14, -2, 0};

    // FKS g = real matrix element * (1.0 - y) * xi * xi
    // real = 942.0125820602134
    double fks_g = 2.475066761582029e-16;
    double limit = FKS::QCD::SoftCollinearLimitFSR(realpdgs[4], bornpdgs[4], ps,
                                                   4, phi, alphas, born, spin);
    // limit = 2.47277512539989e-16
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 0, -13, 14, -2
// realpdgs = -1, 0, -13, 14, -2, 0
TEST(QCDSoftLimit, Wj_born_4_real_22_4) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.8288442858497582;
    ps.X2 = 0.8561362700602047;
    ps.Jacobian = 17319.32196926665;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(5475.468421830074, 0, 0, 5475.468421830074);
    ps.Momenta[1].Set(5475.468421830074, 0, 0, -5475.468421830074);
    ps.Momenta[2].Set(2661.096752010761, 1158.562963850498, 1753.370029537264,
                      -1632.440296573404);
    ps.Momenta[3].Set(5151.433526084111, -547.0922025518221, -2043.581237699081,
                      4696.991933194548);
    ps.Momenta[4].Set(3138.406565565279, -611.4707612986757, 290.2112081618166,
                      -3064.551636621143);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    // double born = 3.40505202017623e-10;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, 5.107578030264345e-10);
    ColorCorr.Set(0, 4, -5.675086700293716e-11);
    ColorCorr.Set(1, 0, 5.107578030264345e-10);
    ColorCorr.Set(1, 4, 5.107578030264345e-10);
    ColorCorr.Set(4, 0, -5.675086700293716e-11);
    ColorCorr.Set(4, 1, 5.107578030264345e-10);

    // radiation variables
    double y = -0.3131695386675162;
    double phi = 5.552335715201948;
    int bornpdgs[] = {-1, 0, -13, 14, -2};
    int realpdgs[] = {-1, 0, -13, 14, -2, 0};

    // FKS g = real matrix element * (1.0 - y) * xi * xi
    // real = 2.869895600166518
    double fks_g = 1.238120026788435e-16;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 4, alphas, ColorCorr, y, phi);
    // limit = 1.238120035729216e-16
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 0, -13, 14, -2
// realpdgs = -1, 0, -13, 14, -2, 0
TEST(QCDCollinearISR1, Wj_born_4_real_22) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.7929555131779971;
    ps.X2 = 0.2386232433291084;
    ps.Jacobian = 754.3171851766655;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2827.444834411853, 0, 0, 2827.444834411853);
    ps.Momenta[1].Set(2827.444834411853, 0, 0, -2827.444834411853);
    ps.Momenta[2].Set(2686.588136466034, -982.0734593136114, 1385.210194281698,
                      2081.941462474364);
    ps.Momenta[3].Set(2703.59816978931, 740.9125848650511, -1401.899528360844,
                      -2189.787505137961);
    ps.Momenta[4].Set(264.7033625683626, 241.1608744485601, 16.68933407914686,
                      107.8460426635976);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 8.847574138582833e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(4.052025097310504e-09, 1.723292942818808e-26);
    spin[1][2] =
        std::complex<double>(2.47660219679255e-10, 4.401175920842885e-09);
    spin[2][1] =
        std::complex<double>(2.47660219679255e-10, -4.401175920842885e-09);
    spin[2][2] =
        std::complex<double>(4.795549041272328e-09, -1.723292942818808e-26);

    // radiation variables
    double phi = 1.752370732369184;
    double xi = 0.03963407065159093;
    int bornpdgs[] = {-1, 0, -13, 14, -2};
    int realpdgs[] = {-1, 0, -13, 14, -2, 0};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 0.0002677639029924763
    double fks_g = 8.412389713147909e-15;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, 1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 8.412368986925861e-15
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 0, -13, 14, -2
// realpdgs = -1, 0, -13, 14, -2, 0
TEST(QCDSoftCollinearISR1, Wj_born_4_real_22) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.05189452889716861;
    ps.X2 = 0.9600865375951315;
    ps.Jacobian = 1131.969274220911;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1450.872781997376, 0, 0, 1450.872781997376);
    ps.Momenta[1].Set(1450.872781997376, 0, 0, -1450.872781997376);
    ps.Momenta[2].Set(1402.481648383662, -20.75630256256945, 1122.795784845148,
                      -840.1509242337671);
    ps.Momenta[3].Set(725.149951754976, 263.6485578755991, -582.9383292328557,
                      341.3426354435782);
    ps.Momenta[4].Set(774.1139638561134, -242.8922553130296, -539.8574556122918,
                      498.808288790189);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 3.506456653787993e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] =
        std::complex<double>(2.513757498994862e-09, -3.446585885637615e-26);
    spin[1][2] =
        std::complex<double>(4.154798442352477e-10, -1.524067401279008e-09);
    spin[2][1] =
        std::complex<double>(4.154798442352477e-10, 1.524067401279008e-09);
    spin[2][2] =
        std::complex<double>(9.926991547931316e-10, 3.446585885637615e-26);

    // radiation variables
    double phi = 1.342181301574761;
    int bornpdgs[] = {-1, 0, -13, 14, -2};
    int realpdgs[] = {-1, 0, -13, 14, -2, 0};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 7330.042517953836
    double fks_g = 1.317800919191135e-14;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, 1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 1.317347667466452e-14
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 0, -13, 14, -2
// realpdgs = -1, 0, -13, 14, -2, 0
TEST(QCDSoftLimit, Wj_born_4_real_22_0) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.5929804417485762;
    ps.X2 = 0.3505650048734239;
    ps.Jacobian = 4018.57748449916;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2963.587958678317, 0, 0, 2963.587958678317);
    ps.Momenta[1].Set(2963.587958678317, 0, 0, -2963.587958678317);
    ps.Momenta[2].Set(2314.089907348445, -524.9490394432655, -1709.909581155203,
                      -1468.213141731529);
    ps.Momenta[3].Set(2267.677651405755, 338.2941759675261, 2229.290218949898,
                      241.213807411365);
    ps.Momenta[4].Set(1345.408358602434, 186.6548634757394, -519.3806377946951,
                      1226.999334320164);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    // double born = 4.651277001704906e-10;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, 6.976915502557359e-10);
    ColorCorr.Set(0, 4, -7.752128336174842e-11);
    ColorCorr.Set(1, 0, 6.976915502557359e-10);
    ColorCorr.Set(1, 4, 6.976915502557359e-10);
    ColorCorr.Set(4, 0, -7.752128336174842e-11);
    ColorCorr.Set(4, 1, 6.976915502557359e-10);

    // radiation variables
    double y = -0.9392535087036507;
    double phi = 4.185624224870872;
    int bornpdgs[] = {-1, 0, -13, 14, -2};
    int realpdgs[] = {-1, 0, -13, 14, -2, 0};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 186.5505041059561
    double fks_g = 9.657233958010885e-16;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 0, alphas, ColorCorr, y, phi);
    // limit = 9.65723661514197e-16
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 0, -13, 14, -2
// realpdgs = -1, 0, -13, 14, -2, 0
TEST(QCDCollinearISR2, Wj_born_4_real_22) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.7580320195852757;
    ps.X2 = 0.8521880209449364;
    ps.Jacobian = 18018.6129107097;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(5224.260744654992, 0, 0, 5224.260744654992);
    ps.Momenta[1].Set(5224.260744654992, 0, 0, -5224.260744654992);
    ps.Momenta[2].Set(2652.439292232177, 1170.710427445151, -485.9272841500037,
                      -2329.966945809945);
    ps.Momenta[3].Set(4373.955318735968, -2845.381915593232, 3163.826280577814,
                      1012.664875977638);
    ps.Momenta[4].Set(3422.12687834184, 1674.671488148082, -2677.89899642781,
                      1317.302069832307);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 1.712208592552505e-10;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(6.738585763929082e-11, 0);
    spin[1][2] =
        std::complex<double>(6.652096833828525e-11, 5.071460744406328e-11);
    spin[2][1] =
        std::complex<double>(6.652096833828525e-11, -5.071460744406328e-11);
    spin[2][2] = std::complex<double>(1.038350016159596e-10, 0);

    // radiation variables
    double phi = 2.106180044559471;
    double xi = 0.04405287944358757;
    int bornpdgs[] = {-1, 0, -13, 14, -2};
    int realpdgs[] = {-1, 0, -13, 14, -2, 0};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 2.757044648399734e-06
    double fks_g = 1.070095150536761e-16;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, -1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 1.070144858176234e-16
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 0, -13, 14, -2
// realpdgs = -1, 0, -13, 14, -2, 0
TEST(QCDSoftCollinearISR2, Wj_born_4_real_22) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.9493359847416052;
    ps.X2 = 0.2101528818290106;
    ps.Jacobian = 6803.455934810677;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2903.29391036912, 0, 0, 2903.29391036912);
    ps.Momenta[1].Set(2903.29391036912, 0, 0, -2903.29391036912);
    ps.Momenta[2].Set(1862.988162863636, 1352.858324927779, -1280.467046332922,
                      -30.05646185221997);
    ps.Momenta[3].Set(1618.518217447949, -402.6027211399805, -194.6605688073618,
                      1555.512626786797);
    ps.Momenta[4].Set(2325.081440426654, -950.2556037877983, 1475.127615140284,
                      -1525.456164934578);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 4.775663877345382e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(2.276951834406547e-09, 0);
    spin[1][2] =
        std::complex<double>(-3.7435981011628e-11, -2.384962372271103e-09);
    spin[2][1] =
        std::complex<double>(-3.7435981011628e-11, 2.384962372271103e-09);
    spin[2][2] = std::complex<double>(2.498712042938835e-09, 0);

    // radiation variables
    double phi = 3.987211750243199;
    int bornpdgs[] = {-1, 0, -13, 14, -2};
    int realpdgs[] = {-1, 0, -13, 14, -2, 0};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 8075.811915240508
    double fks_g = 1.007632839142356e-14;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, -1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 1.008149324687161e-14
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 0, -13, 14, -2
// realpdgs = 0, 0, -13, 14, -2, 1
TEST(QCDCollinearISR1, Wj_born_4_real_35) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.08286203418895255;
    ps.X2 = 0.4416075491771636;
    ps.Jacobian = 741.6453840626589;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1243.395801085176, 0, 0, 1243.395801085176);
    ps.Momenta[1].Set(1243.395801085176, 0, 0, -1243.395801085176);
    ps.Momenta[2].Set(1129.860337114778, -333.8733657364214, 517.5205766492321,
                      -947.251502918628);
    ps.Momenta[3].Set(765.1155434494499, 458.4076646512663, -488.6748732076824,
                      369.4063834212596);
    ps.Momenta[4].Set(591.8157216061244, -124.5342989148448, -28.84570344154966,
                      577.8451194973684);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 2.567413789251735e-09;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(1.445312990701573e-09, 0);
    spin[1][2] =
        std::complex<double>(8.923741132146776e-10, -9.085457077579946e-10);
    spin[2][1] =
        std::complex<double>(8.923741132146776e-10, 9.085457077579946e-10);
    spin[2][2] = std::complex<double>(1.122100798550161e-09, 0);

    // radiation variables
    double phi = 3.518227489020948;
    double xi = 0.8419599081621627;
    int bornpdgs[] = {-1, 0, -13, 14, -2};
    int realpdgs[] = {0, 0, -13, 14, -2, 1};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 1.073339873333072e-07
    double fks_g = 1.521773730222583e-15;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, 1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 1.521532954317785e-15
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-3) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 0, -13, 14, -2
// realpdgs = 0, 0, -13, 14, -2, 1
TEST(QCDSoftCollinearISR1, Wj_born_4_real_35) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.5519652887906794;
    ps.X2 = 0.3761384282095932;
    ps.Jacobian = 1800.895451607503;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(2961.71382773913, 0, 0, 2961.71382773913);
    ps.Momenta[1].Set(2961.71382773913, 0, 0, -2961.71382773913);
    ps.Momenta[2].Set(2754.720432983149, 1866.969042446565, -793.2239964887932,
                      -1863.788332894959);
    ps.Momenta[3].Set(2565.390997983722, -1839.262502097968, 1144.595488680132,
                      1374.13448694272);
    ps.Momenta[4].Set(603.3162245113881, -27.70654034859747, -351.3714921913386,
                      489.6538459522392);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 2.54721151147001e-10;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(4.332989175371216e-11, 0);
    spin[1][2] =
        std::complex<double>(8.209493351745589e-11, 4.919331537878479e-11);
    spin[2][1] =
        std::complex<double>(8.209493351745589e-11, -4.919331537878479e-11);
    spin[2][2] = std::complex<double>(2.113912593932889e-10, 0);

    // radiation variables
    double phi = 4.594552553629486;
    int bornpdgs[] = {-1, 0, -13, 14, -2};
    int realpdgs[] = {0, 0, -13, 14, -2, 1};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 4.806060171891003e-05
    double fks_g = 1.929490494538239e-23;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, 1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 0, -13, 14, -2
// realpdgs = 0, 0, -13, 14, -2, 1
TEST(QCDSoftLimit, Wj_born_4_real_35_0) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.1493751036273956;
    ps.X2 = 0.556697945610523;
    ps.Jacobian = 2264.00704642324;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1874.40000067951, 0, 0, 1874.40000067951);
    ps.Momenta[1].Set(1874.40000067951, 0, 0, -1874.40000067951);
    ps.Momenta[2].Set(679.5851732825192, 546.8859706034718, -400.5785267000682,
                      47.83917692992279);
    ps.Momenta[3].Set(1870.778142423937, -1384.768574397406, 1257.352481898571,
                      -35.93869496845063);
    ps.Momenta[4].Set(1198.436685652564, 837.8826037939341, -856.7739551985023,
                      -11.90048196147215);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 3.636940548773866e-09;

    // color correlated born
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, 5.455410823160799e-09);
    ColorCorr.Set(0, 4, -6.061567581289776e-10);
    ColorCorr.Set(1, 0, 5.455410823160799e-09);
    ColorCorr.Set(1, 4, 5.455410823160799e-09);
    ColorCorr.Set(4, 0, -6.061567581289776e-10);
    ColorCorr.Set(4, 1, 5.455410823160799e-09);

    // radiation variables
    double y = 0.09503683363795545;
    double phi = 5.51965364508433;
    int bornpdgs[] = {-1, 0, -13, 14, -2};
    int realpdgs[] = {0, 0, -13, 14, -2, 1};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 2.056914433497829e-07
    double fks_g = 1.061515106444762e-23;
    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 0, alphas, ColorCorr, y, phi);
    // limit = 0
    EXPECT_NEAR(fks_g / born, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 0, -13, 14, -2
// realpdgs = 0, 0, -13, 14, -2, 1
TEST(QCDCollinearISR2, Wj_born_4_real_35) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.0304772565431044;
    ps.X2 = 0.388137088232078;
    ps.Jacobian = 304.2465190966799;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(706.9584076199776, 0, 0, 706.9584076199776);
    ps.Momenta[1].Set(706.9584076199776, 0, 0, -706.9584076199776);
    ps.Momenta[2].Set(381.4769819557209, 233.9256353913204, -297.611014650765,
                      47.23524983901896);
    ps.Momenta[3].Set(605.4364158973041, -534.7374095268226, 169.7941665047141,
                      227.5502088948412);
    ps.Momenta[4].Set(427.0034173869303, 300.8117741355022, 127.8168481460509,
                      -274.7854587338602);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 2.469656075018611e-08;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(1.936460390994799e-08, 0);
    spin[1][2] =
        std::complex<double>(-7.536553769099545e-09, 6.815540298583856e-09);
    spin[2][1] =
        std::complex<double>(-7.536553769099545e-09, -6.815540298583856e-09);
    spin[2][2] = std::complex<double>(5.331956840238124e-09, 0);

    // radiation variables
    double phi = 0.9917645376242432;
    double xi = 0.07857370097766184;
    int bornpdgs[] = {-1, 0, -13, 14, -2};
    int realpdgs[] = {0, 0, -13, 14, -2, 1};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 3.454832350920523e-05
    double fks_g = 4.265907091850161e-15;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, -1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

// bornpdgs = -1, 0, -13, 14, -2
// realpdgs = 0, 0, -13, 14, -2, 1
TEST(QCDSoftCollinearISR2, Wj_born_4_real_35) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 0.6440682030956459;
    ps.X2 = 0.05155437352355818;
    ps.Jacobian = 144.7169539587829;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(1184.437211208317, 0, 0, 1184.437211208317);
    ps.Momenta[1].Set(1184.437211208317, 0, 0, -1184.437211208317);
    ps.Momenta[2].Set(1099.11604972025, 875.1741013426983, 212.2282868697479,
                      -630.1472346555043);
    ps.Momenta[3].Set(1148.529248928517, -847.6149541703132, -203.370432491073,
                      747.8695021857455);
    ps.Momenta[4].Set(121.2291237678663, -27.55914717238513, -8.857854378674887,
                      -117.7222675302411);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 8.929309411457764e-08;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[1][1] = std::complex<double>(3.176528413827715e-08, 0);
    spin[1][2] =
        std::complex<double>(-8.633749570936952e-09, 4.186699891437582e-08);
    spin[2][1] =
        std::complex<double>(-8.633749570936952e-09, -4.186699891437582e-08);
    spin[2][2] = std::complex<double>(5.75278099763005e-08, 0);

    // radiation variables
    double phi = 4.117200485935669;
    int bornpdgs[] = {-1, 0, -13, 14, -2};
    int realpdgs[] = {0, 0, -13, 14, -2, 1};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    // real = 0.01587120681942724
    double fks_g = 2.855386056038934e-20;
    double limit = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, -1,
                                                   phi, ps.X1 * ps.X2 * ps.S,
                                                   alphas, born, spin);
    // limit = 0
    EXPECT_NEAR(limit, 0.0, 1e-10) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

