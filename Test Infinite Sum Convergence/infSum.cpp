#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <cmath>
using namespace std;

const long double PI = 3.141592653589793238463;

void infSum(){
    const size_t N = 20;
    const long double x = 0.5, y = 0.5;
    const long double L = 1, W = 1;
    vector<long double> fact(N, 0.0);
    std::ofstream output_file("infSum.txt");
    for (size_t i = 1; i < N; ++i){
        fact[i] = fact[i-1] + (pow(-1, i+1) + 1)/i * sin(i*PI*x/L) * sinh(i*PI*y/L) / sinh(i*PI*W/L);
        output_file << std::fixed << std::setprecision(std::numeric_limits<double>::digits10+2) << fact[i]<<"\n";
    }
    output_file << endl;   
}



int main() {
  infSum();
  return 0;
}