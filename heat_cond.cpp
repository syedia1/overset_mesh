#include <iostream>
#include <vector>
#include <fstream>
#include <iterator>
#include <cmath>
using namespace std;

typedef vector<vector<double> > vec2d;
typedef vector<double> vec1d;

void vecToText(const vec2d& T, string filename) {
  // saving results to file
  std::ofstream output_file(filename);
  size_t Ny = T.size()-1-2, Nx = T[0].size()-1-2; // 1 based indexing, 2 boundary cells
  for(size_t i = 2; i <= Nx+1; i++){
    for(size_t j = 2; j <= Ny+1; j++){
			output_file<<T[j][i]<<" ";
    }
		output_file<<endl;
  }
}

void meshTransfer(const vec2d &T, vec2d &T2, const size_t x2, const size_t y2) {
  // size_t Ny = T.size()-1-2, Nx = T[0].size()-1-2; // 1 based indexing, 2 boundary cells
  size_t Ny2 = T2.size()-1-2, Nx2 = T2[0].size()-1-2; // 1 based indexing, 2 boundary cells
  
  // bottom edge
  for (size_t i = 2; i <= Nx2+1; i++) {
    T2[1][i] =  T[2+y2-1][i+x2];
  } 
  // right edge
  for (size_t j = 2; j <= Ny2+1; j++) {
    T2[j][Nx2+2] =  T[j+y2][2+x2+Nx2];
  } 
  // top edge
  for (size_t i = 2; i <= Nx2+1; i++) {
    T2[Ny2+2][i] =  T[2+y2+Ny2][i+x2];
  } 
  // left edge
  for (size_t j = 2; j <= Ny2+1; j++) {
    T2[j][1] =  T[j+y2][2+x2-1];
  } 
}

void BC(vec2d &T, const vec1d T_b) {
  size_t Ny = T.size()-1-2, Nx = T[0].size()-1-2; // 1 based indexing, 2 boundary cells
  // T_b is counter clockwise from bottom edge
  
  // bottom edge
  for (size_t i = 2; i <= Nx+1; i++) {
    T[1][i] = 2 * T_b[0] - T[2][i];
  } 
  // right edge
  for (size_t j = 2; j <= Ny+1; j++) {
    T[j][Nx+2] = 2 * T_b[1] - T[j][Nx+1];
  } 
  // top edge
  for (size_t i = 2; i <= Nx+1; i++) {
    T[Ny+2][i] = 2 * T_b[2] - T[Ny+1][i];
  } 
  // left edge
  for (size_t j = 2; j <= Ny+1; j++) {
    T[j][1] = 2 * T_b[3] - T[j][2];
  } 
}

void gridSolve(vec2d &Tn1, const vec2d &Tn, const double del_t, const double alpha, const double del_x, const double del_y){
  size_t Ny = Tn.size()-1-2, Nx = Tn[0].size()-1-2; // 1 based indexing, 2 boundary cells
  for (size_t i = 2; i <= Nx+1; i++) {
    for (size_t j = 2; j <= Ny+1; j++) {
      Tn1[j][i] = Tn[j][i] + del_t * alpha * ((Tn[j-1][i] - 2 * Tn[j][i] + Tn[j+1][i])/pow(del_y, 2) + (Tn[j][i-1] - 2 * Tn[j][i] + Tn[j][i+1])/pow(del_x, 2));
    }
  }
}

double rmsT(const vec2d &Tn, const vec2d &Tn1) {
  size_t Ny = Tn.size()-1-2, Nx = Tn[0].size()-1-2; // 1 based indexing, 2 boundary cells
  double rms_T = 0.0;
  for (size_t i = 2; i <= Nx+1; i++) {
    for (size_t j = 2; j <= Ny+1; j++) {
      rms_T += pow(Tn1[j][i] - Tn[j][i], 2);
    }
  }
  rms_T = rms_T / (Nx * Ny);
  return rms_T;
}

void solver() {
  // Grid discretization AND Domain 
  size_t x2 = 1, y2 = 1;
  size_t Nx = 100, Ny = 100;
  size_t Nx2 = 40, Ny2 = 98;
  double L = 1, W = 1;
  double L2 = L * Nx2/Nx, W2 = W * Ny2/Ny;

  double del_x = L / Nx, del_y = W / Ny;
  double del_x2 = L2 / Nx2, del_y2 = W2 / Ny2;
 
  // constant
  double k = 16.2;
  double rho = 7750;
  double cp = 500.0;
  double alpha = k / (rho * cp);
  double del_t = (1 / (2 * alpha * (1/pow(del_x, 2) + 1/pow(del_y,2)))) * 0.9; // taking 90 % of the upper limit as per von Neumann stability analysis
  double del_t2 = (1 / (2 * alpha * (1/pow(del_x2, 2) + 1/pow(del_y2,2)))) * 0.9; // taking 90 % of the upper limit as per von Neumann stability analysis
  
  // 1 based indexing and 2 boundary cells
	vec2d Tn(Ny+2 + 1, vec1d(Nx+2 + 1, 0));
	vec2d Tn1(Ny+2 + 1, vec1d(Nx+2 + 1, 0));

	vec2d T2n(Ny2+2 + 1, vec1d(Nx2+2 + 1, 0));
	vec2d T2n1(Ny2+2 + 1, vec1d(Nx2+2 + 1, 0));

  // boundary condition values; counter clockwise from bottom edge
  const vec1d T_bound = {100.0, 200.0, 200.0, 200.0};
  void BC(vec2d &T, const vec1d T_b);

  // function declarations
  void meshTransfer(const vec2d &T, vec2d &T2, const size_t x2, const size_t y2);
  void gridSolve(vec2d &Tn1, const vec2d &Tn, const double del_t, const double alpha, const double del_x, const double del_y);
  double rmsT(const vec2d &Tn, const vec2d &Tn1);

  BC(Tn, T_bound);
  meshTransfer(Tn, T2n, x2, y2);

  // convergence constants
  double rms_T = 1.0;
  double rms_T2 = 1.0;
  const double convergence_factor = 1e-5;

  int n = 0;
  while (rms_T > convergence_factor || rms_T2 > convergence_factor){
    gridSolve(Tn1, Tn, del_t, alpha, del_x, del_y);
    BC(Tn1, T_bound);
    meshTransfer(Tn1, T2n, x2, y2);
    gridSolve(T2n1, T2n, del_t2, alpha, del_x2, del_y2);

    rms_T = rmsT(Tn1, Tn);
    rms_T2 = rmsT(T2n1, T2n);
    n += 1;
    Tn = Tn1;
    T2n = T2n1;
    if (n % 1000 == 0){
      cout << "Iterations = " << n << " RMS_T = " <<  rms_T << " RMS_T2 = " <<  rms_T2 << "\n";
    }
  }
  cout << "Total Iterations = " << n << '\n';
  
  vecToText(Tn1, "./mesh1.txt");
  vecToText(T2n1, "./mesh2.txt");
}

int main() {
  solver();
  return 0;
}