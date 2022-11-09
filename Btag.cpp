#include <iostream>
#include <random>

using namespace std;

bool Btag(double rate) {
  random_device rd;
  uniform_int_distribution<int> irn(0, 1000);
  double rrn;
  rrn = irn(rd)/1000.0;
  if (rrn<rate) return true;
  return false;
}
