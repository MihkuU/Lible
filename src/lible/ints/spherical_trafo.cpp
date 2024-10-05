#include <lible/ints/spherical_trafo.hpp>
#include <lible/ints/defs.hpp>
#include <lible/ints/utils.hpp>

#include <cassert>
#include <stdexcept>

namespace LI = lible::ints;

using std::tuple, std::vector;

arma::dmat LI::returnSphericalTrafo(const int l)
{
    int dim_cart = dimCartesians(l);
    int dim_sph = dimSphericals(l);

    arma::dmat trafo(dim_sph, dim_cart, arma::fill::zeros);

    switch (l)
    {
    case (0):
    {
        trafo(0, 0) = 1.00000000000000;

        return trafo;
    }
    case (1):
    {
        trafo(0, 2) = 1.00000000000000;
        trafo(1, 0) = 1.00000000000000;
        trafo(2, 1) = 1.00000000000000;

        return trafo;
    }
    case (2):
    {
        trafo(0, 0) = -0.50000000000000;
        trafo(0, 3) = -0.50000000000000;
        trafo(0, 5) = 1.00000000000000;
        trafo(1, 2) = 1.73205080756888;
        trafo(2, 4) = 1.73205080756888;
        trafo(3, 0) = 0.86602540378444;
        trafo(3, 3) = -0.86602540378444;
        trafo(4, 1) = 1.73205080756888;

        return trafo;
    }
    case (3):
    {
        trafo(0, 2) = -1.50000000000000;
        trafo(0, 7) = -1.50000000000000;
        trafo(0, 9) = 1.00000000000000;
        trafo(1, 0) = -0.61237243569579;
        trafo(1, 3) = -0.61237243569579;
        trafo(1, 5) = 2.44948974278318;
        trafo(2, 1) = -0.61237243569579;
        trafo(2, 6) = -0.61237243569579;
        trafo(2, 8) = 2.44948974278318;
        trafo(3, 2) = 1.93649167310371;
        trafo(3, 7) = -1.93649167310371;
        trafo(4, 4) = 3.87298334620742;
        trafo(5, 0) = 0.79056941504209;
        trafo(5, 3) = -2.37170824512628;
        trafo(6, 1) = 2.37170824512628;
        trafo(6, 6) = -0.79056941504209;

        return trafo;
    }
    case (4):
    {
        trafo(0, 0) = 0.37500000000000;
        trafo(0, 3) = 0.75000000000000;
        trafo(0, 5) = -3.00000000000000;
        trafo(0, 10) = 0.37500000000000;
        trafo(0, 12) = -3.00000000000000;
        trafo(0, 14) = 1.00000000000000;
        trafo(1, 2) = -2.37170824512628;
        trafo(1, 7) = -2.37170824512628;
        trafo(1, 9) = 3.16227766016838;
        trafo(2, 4) = -2.37170824512628;
        trafo(2, 11) = -2.37170824512628;
        trafo(2, 13) = 3.16227766016838;
        trafo(3, 0) = -0.55901699437495;
        trafo(3, 5) = 3.35410196624968;
        trafo(3, 10) = 0.55901699437495;
        trafo(3, 12) = -3.35410196624968;
        trafo(4, 1) = -1.11803398874989;
        trafo(4, 6) = -1.11803398874989;
        trafo(4, 8) = 6.70820393249937;
        trafo(5, 2) = 2.09165006633519;
        trafo(5, 7) = -6.27495019900557;
        trafo(6, 4) = 6.27495019900557;
        trafo(6, 11) = -2.09165006633519;
        trafo(7, 0) = 0.73950997288745;
        trafo(7, 3) = -4.43705983732471;
        trafo(7, 10) = 0.73950997288745;
        trafo(8, 1) = 2.95803989154981;
        trafo(8, 6) = -2.95803989154981;

        return trafo;
    }
    case (5):
    {
        trafo(0, 2) = 1.87500000000000;
        trafo(0, 7) = 3.75000000000000;
        trafo(0, 9) = -5.00000000000000;
        trafo(0, 16) = 1.87500000000000;
        trafo(0, 18) = -5.00000000000000;
        trafo(0, 20) = 1.00000000000000;
        trafo(1, 0) = 0.48412291827593;
        trafo(1, 3) = 0.96824583655185;
        trafo(1, 5) = -5.80947501931113;
        trafo(1, 10) = 0.48412291827593;
        trafo(1, 12) = -5.80947501931113;
        trafo(1, 14) = 3.87298334620742;
        trafo(2, 1) = 0.48412291827593;
        trafo(2, 6) = 0.96824583655185;
        trafo(2, 8) = -5.80947501931113;
        trafo(2, 15) = 0.48412291827593;
        trafo(2, 17) = -5.80947501931113;
        trafo(2, 19) = 3.87298334620742;
        trafo(3, 2) = -2.56173769148990;
        trafo(3, 9) = 5.12347538297980;
        trafo(3, 16) = 2.56173769148990;
        trafo(3, 18) = -5.12347538297980;
        trafo(4, 4) = -5.12347538297980;
        trafo(4, 11) = -5.12347538297980;
        trafo(4, 13) = 10.24695076595960;
        trafo(5, 0) = -0.52291251658380;
        trafo(5, 3) = 1.04582503316759;
        trafo(5, 5) = 4.18330013267038;
        trafo(5, 10) = 1.56873754975139;
        trafo(5, 12) = -12.54990039801113;
        trafo(6, 1) = -1.56873754975139;
        trafo(6, 6) = -1.04582503316759;
        trafo(6, 8) = 12.54990039801113;
        trafo(6, 15) = 0.52291251658380;
        trafo(6, 17) = -4.18330013267038;
        trafo(7, 2) = 2.21852991866236;
        trafo(7, 7) = -13.31117951197414;
        trafo(7, 16) = 2.21852991866236;
        trafo(8, 4) = 8.87411967464942;
        trafo(8, 11) = -8.87411967464942;
        trafo(9, 0) = 0.70156076002011;
        trafo(9, 3) = -7.01560760020114;
        trafo(9, 10) = 3.50780380010057;
        trafo(10, 1) = 3.50780380010057;
        trafo(10, 6) = -7.01560760020114;
        trafo(10, 15) = 0.70156076002011;

        return trafo;
    }
    case (6):
    {
        trafo(0, 0) = -0.31250000000000;
        trafo(0, 3) = -0.93750000000000;
        trafo(0, 5) = 5.62500000000000;
        trafo(0, 10) = -0.93750000000000;
        trafo(0, 12) = 11.25000000000000;
        trafo(0, 14) = -7.50000000000000;
        trafo(0, 21) = -0.31250000000000;
        trafo(0, 23) = 5.62500000000000;
        trafo(0, 25) = -7.50000000000000;
        trafo(0, 27) = 1.00000000000000;
        trafo(1, 2) = 2.86410980934740;
        trafo(1, 7) = 5.72821961869480;
        trafo(1, 9) = -11.45643923738960;
        trafo(1, 16) = 2.86410980934740;
        trafo(1, 18) = -11.45643923738960;
        trafo(1, 20) = 4.58257569495584;
        trafo(2, 4) = 2.86410980934740;
        trafo(2, 11) = 5.72821961869480;
        trafo(2, 13) = -11.45643923738960;
        trafo(2, 22) = 2.86410980934740;
        trafo(2, 24) = -11.45643923738960;
        trafo(2, 26) = 4.58257569495584;
        trafo(3, 0) = 0.45285552331842;
        trafo(3, 3) = 0.45285552331842;
        trafo(3, 5) = -7.24568837309472;
        trafo(3, 10) = -0.45285552331842;
        trafo(3, 14) = 7.24568837309472;
        trafo(3, 21) = -0.45285552331842;
        trafo(3, 23) = 7.24568837309472;
        trafo(3, 25) = -7.24568837309472;
        trafo(4, 1) = 0.90571104663684;
        trafo(4, 6) = 1.81142209327368;
        trafo(4, 8) = -14.49137674618944;
        trafo(4, 15) = 0.90571104663684;
        trafo(4, 17) = -14.49137674618944;
        trafo(4, 19) = 14.49137674618944;
        trafo(5, 2) = -2.71713313991052;
        trafo(5, 7) = 5.43426627982104;
        trafo(5, 9) = 7.24568837309472;
        trafo(5, 16) = 8.15139941973156;
        trafo(5, 18) = -21.73706511928416;
        trafo(6, 4) = -8.15139941973156;
        trafo(6, 11) = -5.43426627982104;
        trafo(6, 13) = 21.73706511928416;
        trafo(6, 22) = 2.71713313991052;
        trafo(6, 24) = -7.24568837309472;
        trafo(7, 0) = -0.49607837082461;
        trafo(7, 3) = 2.48039185412305;
        trafo(7, 5) = 4.96078370824611;
        trafo(7, 10) = 2.48039185412305;
        trafo(7, 12) = -29.76470224947665;
        trafo(7, 21) = -0.49607837082461;
        trafo(7, 23) = 4.96078370824611;
        trafo(8, 1) = -1.98431348329844;
        trafo(8, 8) = 19.84313483298443;
        trafo(8, 15) = 1.98431348329844;
        trafo(8, 17) = -19.84313483298443;
        trafo(9, 2) = 2.32681380862329;
        trafo(9, 7) = -23.26813808623286;
        trafo(9, 16) = 11.63406904311643;
        trafo(10, 4) = 11.63406904311643;
        trafo(10, 11) = -23.26813808623286;
        trafo(10, 22) = 2.32681380862329;
        trafo(11, 0) = 0.67169328938140;
        trafo(11, 3) = -10.07539934072094;
        trafo(11, 10) = 10.07539934072094;
        trafo(11, 21) = -0.67169328938140;
        trafo(12, 1) = 4.03015973628838;
        trafo(12, 6) = -13.43386578762792;
        trafo(12, 15) = 4.03015973628838;

        return trafo;
    }
    case (7):
    {
        trafo(0, 2) = -2.18750000000000;
        trafo(0, 7) = -6.56250000000000;
        trafo(0, 9) = 13.12500000000000;
        trafo(0, 16) = -6.56250000000000;
        trafo(0, 18) = 26.25000000000000;
        trafo(0, 20) = -10.50000000000000;
        trafo(0, 29) = -2.18750000000000;
        trafo(0, 31) = 13.12500000000000;
        trafo(0, 33) = -10.50000000000000;
        trafo(0, 35) = 1.00000000000000;
        trafo(1, 0) = -0.41339864235384;
        trafo(1, 3) = -1.24019592706153;
        trafo(1, 5) = 9.92156741649221;
        trafo(1, 10) = -1.24019592706153;
        trafo(1, 12) = 19.84313483298443;
        trafo(1, 14) = -19.84313483298443;
        trafo(1, 21) = -0.41339864235384;
        trafo(1, 23) = 9.92156741649221;
        trafo(1, 25) = -19.84313483298443;
        trafo(1, 27) = 5.29150262212918;
        trafo(2, 1) = -0.41339864235384;
        trafo(2, 6) = -1.24019592706153;
        trafo(2, 8) = 9.92156741649221;
        trafo(2, 15) = -1.24019592706153;
        trafo(2, 17) = 19.84313483298443;
        trafo(2, 19) = -19.84313483298443;
        trafo(2, 28) = -0.41339864235384;
        trafo(2, 30) = 9.92156741649221;
        trafo(2, 32) = -19.84313483298443;
        trafo(2, 34) = 5.29150262212918;
        trafo(3, 2) = 3.03784720237868;
        trafo(3, 7) = 3.03784720237868;
        trafo(3, 9) = -16.20185174601965;
        trafo(3, 16) = -3.03784720237868;
        trafo(3, 20) = 9.72111104761179;
        trafo(3, 29) = -3.03784720237868;
        trafo(3, 31) = 16.20185174601965;
        trafo(3, 33) = -9.72111104761179;
        trafo(4, 4) = 6.07569440475737;
        trafo(4, 11) = 12.15138880951474;
        trafo(4, 13) = -32.40370349203930;
        trafo(4, 22) = 6.07569440475737;
        trafo(4, 24) = -32.40370349203930;
        trafo(4, 26) = 19.44222209522358;
        trafo(5, 0) = 0.42961647140211;
        trafo(5, 3) = -0.42961647140211;
        trafo(5, 5) = -8.59232942804220;
        trafo(5, 10) = -2.14808235701055;
        trafo(5, 12) = 17.18465885608440;
        trafo(5, 14) = 11.45643923738960;
        trafo(5, 21) = -1.28884941420633;
        trafo(5, 23) = 25.77698828412660;
        trafo(5, 25) = -34.36931771216879;
        trafo(6, 1) = 1.28884941420633;
        trafo(6, 6) = 2.14808235701055;
        trafo(6, 8) = -25.77698828412660;
        trafo(6, 15) = 0.42961647140211;
        trafo(6, 17) = -17.18465885608440;
        trafo(6, 19) = 34.36931771216879;
        trafo(6, 28) = -0.42961647140211;
        trafo(6, 30) = 8.59232942804220;
        trafo(6, 32) = -11.45643923738960;
        trafo(7, 2) = -2.84975327879450;
        trafo(7, 7) = 14.24876639397250;
        trafo(7, 9) = 9.49917759598167;
        trafo(7, 16) = 14.24876639397250;
        trafo(7, 18) = -56.99506557588999;
        trafo(7, 29) = -2.84975327879450;
        trafo(7, 31) = 9.49917759598167;
        trafo(8, 4) = -11.39901311517800;
        trafo(8, 13) = 37.99671038392666;
        trafo(8, 22) = 11.39901311517800;
        trafo(8, 24) = -37.99671038392666;
        trafo(9, 0) = -0.47495887979908;
        trafo(9, 3) = 4.27462991819175;
        trafo(9, 5) = 5.69950655758900;
        trafo(9, 10) = 2.37479439899542;
        trafo(9, 12) = -56.99506557588999;
        trafo(9, 21) = -2.37479439899542;
        trafo(9, 23) = 28.49753278794499;
        trafo(10, 1) = -2.37479439899542;
        trafo(10, 6) = 2.37479439899542;
        trafo(10, 8) = 28.49753278794499;
        trafo(10, 15) = 4.27462991819175;
        trafo(10, 17) = -56.99506557588999;
        trafo(10, 28) = -0.47495887979908;
        trafo(10, 30) = 5.69950655758900;
        trafo(11, 2) = 2.42182459624970;
        trafo(11, 7) = -36.32736894374543;
        trafo(11, 16) = 36.32736894374543;
        trafo(11, 29) = -2.42182459624970;
        trafo(12, 4) = 14.53094757749817;
        trafo(12, 11) = -48.43649192499390;
        trafo(12, 22) = 14.53094757749817;
        trafo(13, 0) = 0.64725984928775;
        trafo(13, 3) = -13.59245683504274;
        trafo(13, 10) = 22.65409472507123;
        trafo(13, 21) = -4.53081894501425;
        trafo(14, 1) = 4.53081894501425;
        trafo(14, 6) = -22.65409472507123;
        trafo(14, 15) = 13.59245683504274;
        trafo(14, 28) = -0.64725984928775;

        return trafo;
    }
    case (8):
    {
        trafo(0, 0) = 0.27343750000000;
        trafo(0, 3) = 1.09375000000000;
        trafo(0, 5) = -8.75000000000000;
        trafo(0, 10) = 1.64062500000000;
        trafo(0, 12) = -26.25000000000000;
        trafo(0, 14) = 26.25000000000000;
        trafo(0, 21) = 1.09375000000000;
        trafo(0, 23) = -26.25000000000000;
        trafo(0, 25) = 52.50000000000000;
        trafo(0, 27) = -14.00000000000000;
        trafo(0, 36) = 0.27343750000000;
        trafo(0, 38) = -8.75000000000000;
        trafo(0, 40) = 26.25000000000000;
        trafo(0, 42) = -14.00000000000000;
        trafo(0, 44) = 1.00000000000000;
        trafo(1, 2) = -3.28125000000000;
        trafo(1, 7) = -9.84375000000000;
        trafo(1, 9) = 26.25000000000000;
        trafo(1, 16) = -9.84375000000000;
        trafo(1, 18) = 52.50000000000000;
        trafo(1, 20) = -31.50000000000000;
        trafo(1, 29) = -3.28125000000000;
        trafo(1, 31) = 26.25000000000000;
        trafo(1, 33) = -31.50000000000000;
        trafo(1, 35) = 6.00000000000000;
        trafo(2, 4) = -3.28125000000000;
        trafo(2, 11) = -9.84375000000000;
        trafo(2, 13) = 26.25000000000000;
        trafo(2, 22) = -9.84375000000000;
        trafo(2, 24) = 52.50000000000000;
        trafo(2, 26) = -31.50000000000000;
        trafo(2, 37) = -3.28125000000000;
        trafo(2, 39) = 26.25000000000000;
        trafo(2, 41) = -31.50000000000000;
        trafo(2, 43) = 6.00000000000000;
        trafo(3, 0) = -0.39218438743785;
        trafo(3, 3) = -0.78436877487570;
        trafo(3, 5) = 11.76553162313544;
        trafo(3, 12) = 11.76553162313544;
        trafo(3, 14) = -31.37475099502783;
        trafo(3, 21) = 0.78436877487570;
        trafo(3, 23) = -11.76553162313544;
        trafo(3, 27) = 12.54990039801113;
        trafo(3, 36) = 0.39218438743785;
        trafo(3, 38) = -11.76553162313544;
        trafo(3, 40) = 31.37475099502783;
        trafo(3, 42) = -12.54990039801113;
        trafo(4, 1) = -0.78436877487570;
        trafo(4, 6) = -2.35310632462709;
        trafo(4, 8) = 23.53106324627088;
        trafo(4, 15) = -2.35310632462709;
        trafo(4, 17) = 47.06212649254175;
        trafo(4, 19) = -62.74950199005566;
        trafo(4, 28) = -0.78436877487570;
        trafo(4, 30) = 23.53106324627088;
        trafo(4, 32) = -62.74950199005566;
        trafo(4, 34) = 25.09980079602227;
        trafo(5, 2) = 3.18612102524371;
        trafo(5, 7) = -3.18612102524370;
        trafo(5, 9) = -21.24080683495804;
        trafo(5, 16) = -15.93060512621853;
        trafo(5, 18) = 42.48161366991607;
        trafo(5, 20) = 16.99264546796643;
        trafo(5, 29) = -9.55836307573112;
        trafo(5, 31) = 63.72242050487411;
        trafo(5, 33) = -50.97793640389929;
        trafo(6, 4) = 9.55836307573112;
        trafo(6, 11) = 15.93060512621853;
        trafo(6, 13) = -63.72242050487411;
        trafo(6, 22) = 3.18612102524370;
        trafo(6, 24) = -42.48161366991607;
        trafo(6, 26) = 50.97793640389929;
        trafo(6, 37) = -3.18612102524371;
        trafo(6, 39) = 21.24080683495804;
        trafo(6, 41) = -16.99264546796643;
        trafo(7, 0) = 0.41132645565901;
        trafo(7, 3) = -1.64530582263602;
        trafo(7, 5) = -9.87183493581614;
        trafo(7, 10) = -4.11326455659006;
        trafo(7, 12) = 49.35917467908068;
        trafo(7, 14) = 16.45305822636023;
        trafo(7, 21) = -1.64530582263602;
        trafo(7, 23) = 49.35917467908068;
        trafo(7, 25) = -98.71834935816138;
        trafo(7, 36) = 0.41132645565901;
        trafo(7, 38) = -9.87183493581614;
        trafo(7, 40) = 16.45305822636023;
        trafo(8, 1) = 1.64530582263602;
        trafo(8, 6) = 1.64530582263602;
        trafo(8, 8) = -39.48733974326455;
        trafo(8, 15) = -1.64530582263602;
        trafo(8, 19) = 65.81223290544091;
        trafo(8, 28) = -1.64530582263602;
        trafo(8, 30) = 39.48733974326455;
        trafo(8, 32) = -65.81223290544091;
        trafo(9, 2) = -2.96611725366682;
        trafo(9, 7) = 26.69505528300138;
        trafo(9, 9) = 11.86446901466728;
        trafo(9, 16) = 14.83058626833410;
        trafo(9, 18) = -118.64469014667280;
        trafo(9, 29) = -14.83058626833410;
        trafo(9, 31) = 59.32234507333640;
        trafo(10, 4) = -14.83058626833410;
        trafo(10, 11) = 14.83058626833410;
        trafo(10, 13) = 59.32234507333640;
        trafo(10, 22) = 26.69505528300138;
        trafo(10, 24) = -118.64469014667280;
        trafo(10, 37) = -2.96611725366682;
        trafo(10, 39) = 11.86446901466728;
        trafo(11, 0) = -0.45768182862115;
        trafo(11, 3) = 6.40754560069610;
        trafo(11, 5) = 6.40754560069610;
        trafo(11, 12) = -96.11318401044156;
        trafo(11, 21) = -6.40754560069610;
        trafo(11, 23) = 96.11318401044156;
        trafo(11, 36) = 0.45768182862115;
        trafo(11, 38) = -6.40754560069610;
        trafo(12, 1) = -2.74609097172690;
        trafo(12, 6) = 6.40754560069610;
        trafo(12, 8) = 38.44527360417662;
        trafo(12, 15) = 6.40754560069610;
        trafo(12, 17) = -128.15091201392207;
        trafo(12, 28) = -2.74609097172690;
        trafo(12, 30) = 38.44527360417662;
        trafo(13, 2) = 2.50682661696018;
        trafo(13, 7) = -52.64335895616369;
        trafo(13, 16) = 87.73893159360614;
        trafo(13, 29) = -17.54778631872123;
        trafo(14, 4) = 17.54778631872123;
        trafo(14, 11) = -87.73893159360614;
        trafo(14, 22) = 52.64335895616369;
        trafo(14, 37) = -2.50682661696018;
        trafo(15, 0) = 0.62670665424004;
        trafo(15, 3) = -17.54778631872123;
        trafo(15, 10) = 43.86946579680307;
        trafo(15, 21) = -17.54778631872123;
        trafo(15, 36) = 0.62670665424004;
        trafo(16, 1) = 5.01365323392035;
        trafo(16, 6) = -35.09557263744246;
        trafo(16, 15) = 35.09557263744246;
        trafo(16, 28) = -5.01365323392035;

        return trafo;
    }
    case (9):
    {
        trafo(0, 2) = 2.46093750000000;
        trafo(0, 7) = 9.84375000000000;
        trafo(0, 9) = -26.25000000000000;
        trafo(0, 16) = 14.76562500000000;
        trafo(0, 18) = -78.75000000000000;
        trafo(0, 20) = 47.25000000000000;
        trafo(0, 29) = 9.84375000000000;
        trafo(0, 31) = -78.75000000000000;
        trafo(0, 33) = 94.50000000000000;
        trafo(0, 35) = -18.00000000000000;
        trafo(0, 46) = 2.46093750000000;
        trafo(0, 48) = -26.25000000000000;
        trafo(0, 50) = 47.25000000000000;
        trafo(0, 52) = -18.00000000000000;
        trafo(0, 54) = 1.00000000000000;
        trafo(1, 0) = 0.36685490255856;
        trafo(1, 3) = 1.46741961023424;
        trafo(1, 5) = -14.67419610234237;
        trafo(1, 10) = 2.20112941535136;
        trafo(1, 12) = -44.02258830702711;
        trafo(1, 14) = 58.69678440936949;
        trafo(1, 21) = 1.46741961023424;
        trafo(1, 23) = -44.02258830702711;
        trafo(1, 25) = 117.39356881873897;
        trafo(1, 27) = -46.95742752749559;
        trafo(1, 36) = 0.36685490255856;
        trafo(1, 38) = -14.67419610234237;
        trafo(1, 40) = 58.69678440936949;
        trafo(1, 42) = -46.95742752749559;
        trafo(1, 44) = 6.70820393249937;
        trafo(2, 1) = 0.36685490255856;
        trafo(2, 6) = 1.46741961023424;
        trafo(2, 8) = -14.67419610234237;
        trafo(2, 15) = 2.20112941535136;
        trafo(2, 17) = -44.02258830702711;
        trafo(2, 19) = 58.69678440936949;
        trafo(2, 28) = 1.46741961023424;
        trafo(2, 30) = -44.02258830702711;
        trafo(2, 32) = 117.39356881873897;
        trafo(2, 34) = -46.95742752749559;
        trafo(2, 45) = 0.36685490255856;
        trafo(2, 47) = -14.67419610234237;
        trafo(2, 49) = 58.69678440936949;
        trafo(2, 51) = -46.95742752749559;
        trafo(2, 53) = 6.70820393249937;
        trafo(3, 2) = -3.44140403305831;
        trafo(3, 7) = -6.88280806611662;
        trafo(3, 9) = 34.41404033058310;
        trafo(3, 18) = 34.41404033058310;
        trafo(3, 20) = -55.06246452893296;
        trafo(3, 29) = 6.88280806611662;
        trafo(3, 31) = -34.41404033058310;
        trafo(3, 35) = 15.73213272255227;
        trafo(3, 46) = 3.44140403305831;
        trafo(3, 48) = -34.41404033058310;
        trafo(3, 50) = 55.06246452893296;
        trafo(3, 52) = -15.73213272255227;
        trafo(4, 4) = -6.88280806611662;
        trafo(4, 11) = -20.64842419834986;
        trafo(4, 13) = 68.82808066116620;
        trafo(4, 22) = -20.64842419834986;
        trafo(4, 24) = 137.65616132233239;
        trafo(4, 26) = -110.12492905786593;
        trafo(4, 37) = -6.88280806611662;
        trafo(4, 39) = 68.82808066116620;
        trafo(4, 41) = -110.12492905786593;
        trafo(4, 43) = 31.46426544510455;
        trafo(5, 0) = -0.37548796377181;
        trafo(5, 5) = 13.51756669578516;
        trafo(5, 10) = 2.25292778263086;
        trafo(5, 12) = -13.51756669578516;
        trafo(5, 14) = -45.05855565261719;
        trafo(5, 21) = 3.00390371017448;
        trafo(5, 23) = -67.58783347892577;
        trafo(5, 25) = 90.11711130523439;
        trafo(5, 27) = 24.03122968139584;
        trafo(5, 36) = 1.12646389131543;
        trafo(5, 38) = -40.55270008735547;
        trafo(5, 40) = 135.17566695785158;
        trafo(5, 42) = -72.09368904418750;
        trafo(6, 1) = -1.12646389131543;
        trafo(6, 6) = -3.00390371017448;
        trafo(6, 8) = 40.55270008735547;
        trafo(6, 15) = -2.25292778263086;
        trafo(6, 17) = 67.58783347892577;
        trafo(6, 19) = -135.17566695785158;
        trafo(6, 30) = 13.51756669578516;
        trafo(6, 32) = -90.11711130523439;
        trafo(6, 34) = 72.09368904418750;
        trafo(6, 45) = 0.37548796377181;
        trafo(6, 47) = -13.51756669578516;
        trafo(6, 49) = 45.05855565261719;
        trafo(6, 51) = -24.03122968139584;
        trafo(7, 2) = 3.31621990421700;
        trafo(7, 7) = -13.26487961686800;
        trafo(7, 9) = -26.52975923373599;
        trafo(7, 16) = -33.16219904216999;
        trafo(7, 18) = 132.64879616867995;
        trafo(7, 20) = 26.52975923373599;
        trafo(7, 29) = -13.26487961686800;
        trafo(7, 31) = 132.64879616867995;
        trafo(7, 33) = -159.17855540241595;
        trafo(7, 46) = 3.31621990421700;
        trafo(7, 48) = -26.52975923373599;
        trafo(7, 50) = 26.52975923373599;
        trafo(8, 4) = 13.26487961686800;
        trafo(8, 11) = 13.26487961686800;
        trafo(8, 13) = -106.11903693494396;
        trafo(8, 22) = -13.26487961686800;
        trafo(8, 26) = 106.11903693494396;
        trafo(8, 37) = -13.26487961686800;
        trafo(8, 39) = 106.11903693494396;
        trafo(8, 41) = -106.11903693494396;
        trafo(9, 0) = 0.39636409043643;
        trafo(9, 3) = -3.17091272349146;
        trafo(9, 5) = -11.09819453222009;
        trafo(9, 10) = -5.54909726611005;
        trafo(9, 12) = 99.88375078998087;
        trafo(9, 14) = 22.19638906444019;
        trafo(9, 23) = 55.49097266110048;
        trafo(9, 25) = -221.96389064440191;
        trafo(9, 36) = 1.98182045218216;
        trafo(9, 38) = -55.49097266110048;
        trafo(9, 40) = 110.98194532220096;
        trafo(10, 1) = 1.98182045218216;
        trafo(10, 8) = -55.49097266110048;
        trafo(10, 15) = -5.54909726611005;
        trafo(10, 17) = 55.49097266110048;
        trafo(10, 19) = 110.98194532220096;
        trafo(10, 28) = -3.17091272349146;
        trafo(10, 30) = 99.88375078998087;
        trafo(10, 32) = -221.96389064440191;
        trafo(10, 45) = 0.39636409043643;
        trafo(10, 47) = -11.09819453222009;
        trafo(10, 49) = 22.19638906444019;
        trafo(11, 2) = -3.07022304258990;
        trafo(11, 7) = 42.98312259625865;
        trafo(11, 9) = 14.32770753208622;
        trafo(11, 18) = -214.91561298129324;
        trafo(11, 29) = -42.98312259625865;
        trafo(11, 31) = 214.91561298129324;
        trafo(11, 46) = 3.07022304258990;
        trafo(11, 48) = -14.32770753208622;
        trafo(12, 4) = -18.42133825553942;
        trafo(12, 11) = 42.98312259625864;
        trafo(12, 13) = 85.96624519251729;
        trafo(12, 22) = 42.98312259625864;
        trafo(12, 24) = -286.55415064172428;
        trafo(12, 37) = -18.42133825553942;
        trafo(12, 39) = 85.96624519251729;
        trafo(13, 0) = -0.44314852502787;
        trafo(13, 3) = 8.86297050055736;
        trafo(13, 5) = 7.09037640044589;
        trafo(13, 10) = -6.20407935039015;
        trafo(13, 12) = -148.89790440936369;
        trafo(13, 21) = -12.40815870078031;
        trafo(13, 23) = 248.16317401560616;
        trafo(13, 36) = 3.10203967519508;
        trafo(13, 38) = -49.63263480312123;
        trafo(14, 1) = -3.10203967519508;
        trafo(14, 6) = 12.40815870078031;
        trafo(14, 8) = 49.63263480312123;
        trafo(14, 15) = 6.20407935039015;
        trafo(14, 17) = -248.16317401560616;
        trafo(14, 28) = -8.86297050055736;
        trafo(14, 30) = 148.89790440936369;
        trafo(14, 45) = 0.44314852502787;
        trafo(14, 47) = -7.09037640044589;
        trafo(15, 2) = 2.58397773170915;
        trafo(15, 7) = -72.35137648785613;
        trafo(15, 16) = 180.87844121964034;
        trafo(15, 29) = -72.35137648785613;
        trafo(15, 46) = 2.58397773170915;
        trafo(16, 4) = 20.67182185367318;
        trafo(16, 11) = -144.70275297571226;
        trafo(16, 22) = 144.70275297571226;
        trafo(16, 37) = -20.67182185367318;
        trafo(17, 0) = 0.60904939217552;
        trafo(17, 3) = -21.92577811831886;
        trafo(17, 10) = 76.74022341411600;
        trafo(17, 21) = -51.16014894274400;
        trafo(17, 36) = 5.48144452957971;
        trafo(18, 1) = 5.48144452957971;
        trafo(18, 6) = -51.16014894274400;
        trafo(18, 15) = 76.74022341411600;
        trafo(18, 28) = -21.92577811831886;
        trafo(18, 45) = 0.60904939217552;

        return trafo;
    }
    default:
    {
        throw std::runtime_error("Inappropriate angular momentum given!\n");
    }
    }
}

vector<tuple<int, int, double>> LI::sphericalTrafo(const int l)
{
    arma::dmat sph_trafo_mat = returnSphericalTrafo(l);

    vector<tuple<int, int, double>> sph_trafo;
    for (size_t i = 0; i < sph_trafo_mat.n_rows; i++)
        for (size_t j = 0; j < sph_trafo_mat.n_cols; j++)
        {
            double val = sph_trafo_mat(i, j);
            if (std::fabs(val) != 0)
                sph_trafo.push_back(std::make_tuple(i, j, val));
        }

    return sph_trafo;
}

void LI::transferIntegrals(const int ipair, const ShellPairData &sp_data,
                           const arma::dmat &ints_sph, vec2d &ints)
{
    int pos_a = sp_data.offsets_sph[2 * ipair];
    int pos_b = sp_data.offsets_sph[2 * ipair + 1];
    int pos_norm_a = sp_data.offsets_norms[2 * ipair];
    int pos_norm_b = sp_data.offsets_norms[2 * ipair + 1];

    for (size_t mu = 0; mu < ints_sph.n_rows; mu++)
        for (size_t nu = 0; nu < ints_sph.n_cols; nu++)
        {
            double norm_a = sp_data.norms[pos_norm_a + mu];
            double norm_b = sp_data.norms[pos_norm_b + nu];
            double normalized_int = ints_sph(mu, nu) * norm_a * norm_b;

            ints(pos_a + mu, pos_b + nu) = normalized_int;
            ints(pos_b + nu, pos_a + mu) = normalized_int;
        }
}

void LI::transferIntegrals(const int ipair_ab, const ShellPairData &sp_data_ab,
                           const vector<double> &eri4_shells_sph, vec2d &eri4_diagonal)
{
    int dim_a = dimSphericals(sp_data_ab.la);
    int dim_b = dimSphericals(sp_data_ab.lb);
    int dim_ab = dim_a * dim_b;
    
    int pos_a = sp_data_ab.offsets_sph[2 * ipair_ab];
    int pos_b = sp_data_ab.offsets_sph[2 * ipair_ab + 1];
    int pos_norm_a = sp_data_ab.offsets_norms[2 * ipair_ab];
    int pos_norm_b = sp_data_ab.offsets_norms[2 * ipair_ab + 1];

    for (int mu = 0, munu = 0; mu < dim_a; mu++)
        for (int nu = 0; nu < dim_b; nu++, munu++)
        {
            int munumunu = munu * dim_ab + munu;

            double norm_a = sp_data_ab.norms[pos_norm_a + mu];
            double norm_b = sp_data_ab.norms[pos_norm_b + nu];

            double normalized_int = eri4_shells_sph[munumunu] * norm_a * norm_b * norm_a * norm_b;

            int a = pos_a + mu;
            int b = pos_b + nu;

            eri4_diagonal(a, b) = normalized_int;
            eri4_diagonal(b, a) = normalized_int;
        }        
}

void LI::transferIntegrals(const int ishell_a, const int ishell_b,
                           const ShellData &sh_data_a, const ShellData &sh_data_b,
                           const vector<double> &eri2_shells_sph, vec2d &eri2)
{
    int dim_a = dimSphericals(sh_data_a.l);
    int dim_b = dimSphericals(sh_data_b.l);

    int pos_a = sh_data_a.offsets_sph[ishell_a];
    int pos_b = sh_data_b.offsets_sph[ishell_b];
    int pos_norm_a = sh_data_a.offsets_norms[ishell_a];
    int pos_norm_b = sh_data_b.offsets_norms[ishell_b];

    for (int mu = 0; mu < dim_a; mu++)
        for (int nu = 0; nu < dim_b; nu++)
        {
            int munu = mu * dim_b + nu;

            double norm_a = sh_data_a.norms[pos_norm_a + mu];
            double norm_b = sh_data_b.norms[pos_norm_b + nu];

            double normalized_int = norm_a * norm_b * eri2_shells_sph[munu];

            int a = pos_a + mu;
            int b = pos_b + nu;
            
            // eri2(a, b) = normalized_int;
            // eri2(b, a) = normalized_int;

            eri2(a, b) = norm_a * norm_b; // tmp
            eri2(b, a) = norm_a * norm_b; // tmp
        }
}

void LI::transferIntegrals(const int ipair_ab, const int ishell_c,
                           const ShellData &sh_data_c, const ShellPairData &sp_data_ab,
                           const vector<double> &eri3_shells_sph, vec3d &eri3)
{
    int dim_a = dimSphericals(sp_data_ab.la);
    int dim_b = dimSphericals(sp_data_ab.lb);
    int dim_c = dimSphericals(sh_data_c.l);

    int pos_a = sp_data_ab.offsets_sph[2 * ipair_ab];
    int pos_b = sp_data_ab.offsets_sph[2 * ipair_ab + 1];
    int pos_c = sh_data_c.offsets_sph[ishell_c];

    int pos_norm_a = sp_data_ab.offsets_norms[2 * ipair_ab];
    int pos_norm_b = sp_data_ab.offsets_norms[2 * ipair_ab + 1];
    int pos_norm_c = sh_data_c.offsets_norms[ishell_c];

    for (int mu = 0, munu = 0; mu < dim_a; mu++)
        for (int nu = 0; nu < dim_b; nu++, munu++)
            for (int ka = 0; ka < dim_c; ka++)
            {
                int munuka = munu * dim_c + ka;

                double norm_a = sp_data_ab.norms[pos_norm_a + mu];
                double norm_b = sp_data_ab.norms[pos_norm_b + nu];
                double norm_c = sh_data_c.norms[pos_norm_c + ka];

                double normalized_int = eri3_shells_sph[munuka] * norm_a * norm_b * norm_c;

                int a = pos_a + mu;
                int b = pos_b + nu;
                int c = pos_c + ka;

                eri3(a, b, c) = normalized_int;
                eri3(b, a, c) = normalized_int;
            }
}

void LI::transferIntegrals(const int ipair_ab, const int ipair_cd,
                           const ShellPairData &sp_data_ab,
                           const ShellPairData &sp_data_cd,
                           const vector<double> &eri4_shells_sph, vec4d &eri4)
{
    int dim_a = dimSphericals(sp_data_ab.la);
    int dim_b = dimSphericals(sp_data_ab.lb);
    int dim_c = dimSphericals(sp_data_cd.la);
    int dim_d = dimSphericals(sp_data_cd.lb);

    int pos_a = sp_data_ab.offsets_sph[2 * ipair_ab];
    int pos_b = sp_data_ab.offsets_sph[2 * ipair_ab + 1];
    int pos_c = sp_data_cd.offsets_sph[2 * ipair_cd];
    int pos_d = sp_data_cd.offsets_sph[2 * ipair_cd + 1];

    int pos_norm_a = sp_data_ab.offsets_norms[2 * ipair_ab];
    int pos_norm_b = sp_data_ab.offsets_norms[2 * ipair_ab + 1];
    int pos_norm_c = sp_data_cd.offsets_norms[2 * ipair_cd];
    int pos_norm_d = sp_data_cd.offsets_norms[2 * ipair_cd + 1];

    int dim_cd = dim_c * dim_d;

    for (int mu = 0, munu = 0; mu < dim_a; mu++)
        for (int nu = 0; nu < dim_b; nu++, munu++)
            for (int ka = 0, kata = 0; ka < dim_c; ka++)
                for (int ta = 0; ta < dim_d; ta++, kata++)
                {
                    int munukata = munu * dim_cd + kata;

                    double norm_a = sp_data_ab.norms[pos_norm_a + mu];
                    double norm_b = sp_data_ab.norms[pos_norm_b + nu];
                    double norm_c = sp_data_cd.norms[pos_norm_c + ka];
                    double norm_d = sp_data_cd.norms[pos_norm_d + ta];
                    double normalized_int = eri4_shells_sph[munukata] * norm_a * norm_b * norm_c *
                                            norm_d;

                    int a = pos_a + mu;
                    int b = pos_b + nu;
                    int c = pos_c + ka;
                    int d = pos_d + ta;

                    eri4(a, b, c, d) = normalized_int;
                    eri4(a, b, d, c) = normalized_int;
                    eri4(b, a, c, d) = normalized_int;
                    eri4(b, a, d, c) = normalized_int;
                    eri4(c, d, a, b) = normalized_int;
                    eri4(c, d, b, a) = normalized_int;
                    eri4(d, c, a, b) = normalized_int;
                    eri4(d, c, b, a) = normalized_int;
                }
}