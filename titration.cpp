$PARAM TVCL = 1.3, TVVC=28, TVKA=0.6, WT=70, START_DOSE = 15

$SET delta= 1

$CMT GUT CENT

$PLUGIN Rcpp mrgx

$ENV
possibleDoses <- c(10, 15, 20)
  VISITT <- seq(48, 280, 48)
  
$GLOBAL
using namespace Rcpp;
NumericVector possibleDoses;
NumericVector VISITT;

bool within(Rcpp::NumericVector x, double val) {
  int n = x.size();
  for (int i = 0; i < n; ++i) {
    if (x[i] == val) {
      return true;
    }
  }
  return false;
}
double titrateDose(Rcpp::NumericVector possibleDoses, double currentDose, bool up){
  if (up) {
    possibleDoses = possibleDoses[possibleDoses >= currentDose];
    if (possibleDoses.size() > 1) {
      return possibleDoses[1]; // 2nd element - one dose higher
    }
    return possibleDoses[0]; // at max dose since only one dose remaining that is >= so keep the same
  } else {
    possibleDoses = possibleDoses[possibleDoses <= currentDose];
    if (possibleDoses.size() > 1) {
      return possibleDoses[possibleDoses.size()-2]; // 2nd to last element - one dose lower
    }
    return possibleDoses[0]; // at min dose since only one dose remaining that is <= so keep the same
  }
}

$PREAMBLE
possibleDoses = mrgx::get<Rcpp::NumericVector>("possibleDoses", self);
VISITT = mrgx::get<Rcpp::NumericVector>("VISITT", self);

$MAIN
if (NEWIND <= 1) {
  // titration dose to start on, right now not explicitly checking
  // if in possible doses, probably should do that
  F_GUT = START_DOSE;
}
if (within(VISITT, TIME)) {
  // only adjust dose on EVID == 1 or also during observation time can trigger a dose
  // adjustment if both dosing and observing at the same time and not
  // also checking EVID == 1
  if (CENT < 10 && EVID == 1) {
    F_GUT = titrateDose(possibleDoses, F_GUT, true);
  }
  if (CENT > 15 && EVID == 1) {
    F_GUT = titrateDose(possibleDoses, F_GUT, false);
  }
}
double CLi = exp(log(TVCL) + 0.75*log(WT/70) + ETA(1));
double VCi = exp(log(TVVC) + ETA(2));
double KAi = exp(log(TVKA) + ETA(3));

$OMEGA name="IIV"
0.1 0 0

$ODE
  dxdt_GUT = -KAi*GUT;
dxdt_CENT = KAi*GUT - (CLi/VCi)*CENT;

$TABLE
  double CP = CENT/VCi;
double ETA1 = ETA(1);
double ETA2 = ETA(2);

$CAPTURE ETA(1) ETA(2) F_GUT