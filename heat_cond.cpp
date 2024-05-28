#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <iomanip>
using namespace std;

typedef vector<vector<pair<vector<size_t>, vector<double> > > > vec3id;
typedef vector<vector<double> > vec2d;
typedef vector<pair<vector<size_t>, vector<double> > > vec2id;
typedef vector<vector<size_t> > vec2i;
typedef vector<double> vec1d;
typedef vector<size_t> vec1i;

const double PI = 3.141592653589793;
const double TOLERANCE = 1024 * std::numeric_limits<double>::epsilon();
// const double TOLERANCE = 1e-6;
const size_t NS_MAXITER = 20;


// function declarations
void BC(vec2d &T, const vec1d T_b);
void meshTransfer(const vec2d &T, vec2d &T2, const vec3id &interpolation_pts);
void gridSolve(vec2d &Tn1, const vec2d &Tn, const double del_t, const double alpha, const vec1d del_mesh);
double rmsT(const vec2d &Tn, const vec2d &Tn1);
vec1d mesh2to1(const vec1d mesh2_origin, const vec1d del_mesh2, const size_t i, const size_t j);
vec3id interpolation_pts(const vec2d &T, const vec2d &T2, const vec1d del_mesh, const vec1d del_mesh2, const vec1d mesh2_origin);
vec1d interpolation_coeff(const vec2d points, const vec1d xy);

vec2d matMul(const vec2d &A, const vec2d &B);
vec1d matAdd(const vec1d &A, const vec1d &B);
vec1d solveAxB(const vec2d &A, const vec1d &B);


void vecToText(const vec2d& T, string filename) {
  // saving results to file
  std::ofstream output_file(filename);
  size_t Ny = T.size()-1-2, Nx = T[0].size()-1-2; // 1 based indexing, 2 boundary cells
  for(size_t j = 2; j <= Ny+1; j++){
    for(size_t i = 2; i <= Nx+1; i++){
			output_file << std::fixed << std::setprecision(std::numeric_limits<double>::digits10+2) << T[j][i]<<" ";
    }
		output_file << endl;
  }
}

vec2d matMul(const vec2d &A, const vec2d &B) {
  size_t m = A.size(), n = B[0].size(), r = A[0].size();
  vec2d C (m, vec1d(n, 0.0));
  // order of j and k loops is swapped for contigous access of C and B, cache friendly
  for (size_t i = 0 ; i < m; i++) {
    for (size_t k = 0 ; k < r; k++) {
      for (size_t j = 0 ; j < n; j++) {
        C[i][j] += A[i][k] * B[k][j];
      }
    }
  }
  return C;
}

vec1d matAdd(const vec1d &A, const vec1d &B) {
  size_t m = A.size();
  vec1d C (m, 0.0);
  for (size_t i = 0 ; i < m; i++) {
    C[i] = A[i] + B[i];
  }
  return C;
}

vec1d solveAxB(const vec2d &A, const vec1d &B) {
  size_t n = A.size(); 
  vec1d x(n, 0.0);

  // Pre-calculate inverses of diagonal elements (avoid division within the loop)
  vec1d inv_diag(n, 0.0);
  for (size_t i = 0; i < n; ++i) {
    inv_diag[i] = 1.0 / A[i][i];
  }

  // Loop for iterations, Guass Seidel
  size_t GS_itr = 0;
  for (; GS_itr < NS_MAXITER; ++GS_itr) {

    // Update each element of the solution vector
    double sum = 0.0;
    bool converged = true;
    for (size_t j = 0; j < n; ++j) {
      sum = 0.0;
      for (size_t k = 0; k < n; ++k) {
        if (j != k) {
          sum += A[j][k] * x[k];
        }
      }
      x[j] = (B[j] - sum) * inv_diag[j];
      converged &= (abs(x[j]*A[j][j] + sum - B[j]) <= TOLERANCE);
    }

    if (converged) {
      break;
    }
  }
  // cout << "GS=" << GS_itr << ", ";
  return x;
}

vec1d interpolation_coeff(const vec2d points, const vec1d xy){
  auto transform = [](const vec1d &coeff, const vec1d &points) -> double {
    const double phi = coeff[0], chi = coeff[1];
    return (1.0-phi)*(1.0-chi) * points[0] + (phi)*(1.0-chi) * points[1] + (phi)*(chi) * points[2] + (1.0-phi)*(chi) * points[3];
  };
  auto partialDerivative_phi = [](const vec1d &coeff, const vec1d &points) {
    const double chi = coeff[1];
    return (-1.0)*(1.0-chi) * points[0] + (1.0-chi) * points[1] + (chi) * points[2] + (-1.0)*(chi) * points[3];  
  };
  auto partialDerivative_chi = [](const vec1d &coeff, const vec1d &points) {
    const double phi = coeff[0];
    return (1.0-phi)*(-1.0) * points[0] + (phi) * (-1.0) * points[1] + (phi) * points[2] + (1.0-phi) * points[3];
  };

  vec1d coeff = {0.5, 0.5}; // initial guess value of phi,
  
  // NR, Newton raphson to minimze B, Ax = B
  size_t NR_itr = 0;
  for (; NR_itr < NS_MAXITER; NR_itr++) {
    vec2d A = {{partialDerivative_phi(coeff, points[0]), partialDerivative_chi(coeff, points[0])}, 
      {partialDerivative_phi(coeff, points[1]), partialDerivative_chi(coeff, points[1])}};
    vec1d B = {-1.0*transform(coeff, points[0])+xy[0], -1.0*transform(coeff, points[1])+xy[1]};
    
    coeff = matAdd(coeff, solveAxB(A, B));
    // check for convergence
    bool converged = true;
    for (int eq = 0; eq < 2; ++eq) {
      converged &= (abs(B[eq]) <= TOLERANCE);
    }
    if (converged) {
      break;
    }
  }
  // cout << "NR=" << NR_itr << ", ";
  return vec1d {(1.0-coeff[0])*(1.0-coeff[1]), (coeff[0])*(1.0-coeff[1]), (coeff[0])*(coeff[1]), (1.0-coeff[0])*(coeff[1])};
}

vec1d mesh2to1(const vec1d mesh2_origin, const vec1d del_mesh2, const size_t i, const size_t j){
  // coordinate transforms between 2 cartesian grids shifted by x2, y2 and rotated counter-clock by theta
  vec1d xy(2, 0.0);
  const double x2 = mesh2_origin[0], y2 = mesh2_origin[1];
  const double sinTheta = mesh2_origin[2], cosTheta = mesh2_origin[3];
  const double del_x2 = del_mesh2[0], del_y2 = del_mesh2[1];
  xy[0] = x2 + ((i-2.0)*del_x2*cosTheta - (j-2.0)*del_y2*sinTheta);
  xy[1] = y2 + ((i-2.0)*del_x2*sinTheta + (j-2.0)*del_y2*cosTheta);
  return xy;
} 

vec3id interpolation_pts(const vec2d &T, const vec2d &T2, const vec1d del_mesh, const vec1d del_mesh2, const vec1d mesh2_origin){
  // 3d vector of interpolation points for boundaries of mesh 2 counter clockwise from bottom edge
  // innermost vector has i, j i.e. first interpolation donor cells of mesh 1 and coeffs i.e phi and chi 

  // size_t Ny = T.size()-1-2, Nx = T[0].size()-1-2; // 1 based indexing, 2 boundary cells  
  size_t Ny2 = T2.size()-1-2, Nx2 = T2[0].size()-1-2; // 1 based indexing, 2 boundary cells
  vec3id inter_points;
  vec2id temp_line;
  //bottom edge
  for (size_t i = 2; i <= Nx2+1; i++) {
    vec1d xy = mesh2to1(mesh2_origin, del_mesh2, i, 2);
    size_t _i = floor(xy[0]/del_mesh[0]);
    size_t _j = floor(xy[1]/del_mesh[1]);
    vec1d coeff = interpolation_coeff(vec2d {{_i * del_mesh[0], (_i+1)*del_mesh[0], (_i+1) * del_mesh[0], _i*del_mesh[0]}, {_j*del_mesh[1], _j*del_mesh[1], (_j+1)*del_mesh[1], (_j+1)*del_mesh[1]}}, xy);
    temp_line.push_back(make_pair(vec1i {_i+2, _j+2}, coeff));
  } 
  inter_points.push_back(temp_line);
  temp_line.clear();

  //right edge
  for (size_t j = 2; j <= Ny2+1; j++) {
    vec1d xy = mesh2to1(mesh2_origin, del_mesh2, Nx2+1, j);
    size_t _i = floor(xy[0]/del_mesh[0]);
    size_t _j = floor(xy[1]/del_mesh[1]);
    vec1d coeff = interpolation_coeff(vec2d {{_i * del_mesh[0], (_i+1)*del_mesh[0], (_i+1) * del_mesh[0], (_i)*del_mesh[0]}, {_j*del_mesh[1], _j*del_mesh[1], (_j+1)*del_mesh[1], (_j+1)*del_mesh[1]}}, xy);
    temp_line.push_back(make_pair(vec1i {_i+2, _j+2}, coeff));
  } 
  inter_points.push_back(temp_line);
  temp_line.clear();

  //top edge
  for (size_t i = 2; i <= Nx2+1; i++) {
    vec1d xy = mesh2to1(mesh2_origin, del_mesh2, i, Ny2+1);
    size_t _i = floor(xy[0]/del_mesh[0]);
    size_t _j = floor(xy[1]/del_mesh[1]);
    vec1d coeff = interpolation_coeff(vec2d {{_i * del_mesh[0], (_i+1)*del_mesh[0], (_i+1) * del_mesh[0], (_i)*del_mesh[0]}, {_j*del_mesh[1], _j*del_mesh[1], (_j+1)*del_mesh[1], (_j+1)*del_mesh[1]}}, xy);
    temp_line.push_back(make_pair(vec1i {_i+2, _j+2}, coeff));
  } 
  inter_points.push_back(temp_line);
  temp_line.clear();


  //left edge
  for (size_t j = 2; j <= Ny2+1; j++) {
    vec1d xy = mesh2to1(mesh2_origin, del_mesh2, 2, j);
    size_t _i = floor(xy[0]/del_mesh[0]);
    size_t _j = floor(xy[1]/del_mesh[1]);
    vec1d coeff = interpolation_coeff(vec2d {{_i * del_mesh[0], (_i+1)*del_mesh[0], (_i+1) * del_mesh[0], (_i)*del_mesh[0]}, {_j*del_mesh[1], _j*del_mesh[1], (_j+1)*del_mesh[1], (_j+1)*del_mesh[1]}}, xy);
    temp_line.push_back(make_pair(vec1i {_i+2, _j+2}, coeff));
  } 
  inter_points.push_back(temp_line);
  return inter_points;
}

void meshTransfer(const vec2d &T, vec2d &T2, const vec3id &interpolation_pts) {
  // size_t Ny = T.size()-1-2, Nx = T[0].size()-1-2; // 1 based indexing, 2 boundary cells
  size_t Ny2 = T2.size()-1-2, Nx2 = T2[0].size()-1-2; // 1 based indexing, 2 boundary cells
  
  size_t bIdx = 0;
  // bottom edge
  for (size_t i = 2; i <= Nx2+1; i++) {
    size_t _i = interpolation_pts[bIdx][i-2].first[0], _j = interpolation_pts[bIdx][i-2].first[1];
    vec2i pts = {{_i, _j}, {_i+1, _j}, {_i, _j+1}, {_i+1, _j+1}};
    T2[1][i] = 0.0;
    for (size_t idx = 0; idx < 4; idx++){
      T2[1][i] += T[pts[idx][1]][pts[idx][0]] * interpolation_pts[bIdx][i-2].second[idx];
    }
  } 
  bIdx = 1;
  // right edge
  for (size_t j = 2; j <= Ny2+1; j++) {
    size_t _i = interpolation_pts[bIdx][j-2].first[0], _j = interpolation_pts[bIdx][j-2].first[1];
    vec2i pts = {{_i, _j}, {_i+1, _j}, {_i, _j+1}, {_i+1, _j+1}};
    T2[j][Nx2+2] = 0.0;
    for (size_t idx = 0; idx < 4; idx++){
      T2[j][Nx2+2] += T[pts[idx][1]][pts[idx][0]] * interpolation_pts[bIdx][j-2].second[idx];
    }
  } 
  bIdx = 2;
  // top edge
  for (size_t i = 2; i <= Nx2+1; i++) {
    size_t _i = interpolation_pts[bIdx][i-2].first[0], _j = interpolation_pts[bIdx][i-2].first[1];
    vec2i pts = {{_i, _j}, {_i+1, _j}, {_i, _j+1}, {_i+1, _j+1}};
    T2[Ny2+2][i] = 0.0;
    for (size_t idx = 0; idx < 4; idx++){
      T2[Ny2+2][i] += T[pts[idx][1]][pts[idx][0]] * interpolation_pts[bIdx][i-2].second[idx];
    }
  } 
  bIdx = 3;
  // left edge
  for (size_t j = 2; j <= Ny2+1; j++) {
    size_t _i = interpolation_pts[bIdx][j-2].first[0], _j = interpolation_pts[bIdx][j-2].first[1];
    vec2i pts = {{_i, _j}, {_i+1, _j}, {_i, _j+1}, {_i+1, _j+1}};
    T2[j][1] = 0.0;
    for (size_t idx = 0; idx < 4; idx++){
      T2[j][1] += T[pts[idx][1]][pts[idx][0]] * interpolation_pts[bIdx][j-2].second[idx];
    }
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

void gridSolve(vec2d &Tn1, const vec2d &Tn, const double del_t, const double alpha, const vec1d del_mesh){
  size_t Ny = Tn.size()-1-2, Nx = Tn[0].size()-1-2; // 1 based indexing, 2 boundary cells
  for (size_t i = 2; i <= Nx+1; i++) {
    for (size_t j = 2; j <= Ny+1; j++) {
      Tn1[j][i] = Tn[j][i] + del_t * alpha * ((Tn[j-1][i] - 2 * Tn[j][i] + Tn[j+1][i])/pow(del_mesh[1], 2) + (Tn[j][i-1] - 2 * Tn[j][i] + Tn[j][i+1])/pow(del_mesh[0], 2));
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
  const size_t Nx = 160, Ny = 160;
  const double L = 5.0, W = 5.0;

  const double x2 = L/3, y2 = W/3, theta = -25.0 * PI/180.0;
  const size_t Nx2 = 48, Ny2 = 48;
  const double L2 = 1.5, W2 = 1.5;

  const double del_x = L / Nx, del_y = W / Ny;
  const double del_x2 = L2 / Nx2, del_y2 = W2 / Ny2;
  const double cosTheta = cos(theta), sinTheta = sin(theta);
  
  const vec1d mesh2_origin = {x2, y2, sinTheta, cosTheta};
  const vec1d del_mesh = {del_x, del_y};
  const vec1d del_mesh2 = {del_x2, del_y2};

  // constants
  const double k = 16.2;
  const double rho = 7750;
  const double cp = 500.0;
  const double alpha = k / (rho * cp);
  const double del_t = (1 / (2 * alpha * (1/pow(del_x, 2) + 1/pow(del_y,2)))) * 0.9; // taking 90 % of the upper limit as per von Neumann stability analysis
  const double del_t2 = (1 / (2 * alpha * (1/pow(del_x2, 2) + 1/pow(del_y2,2)))) * 0.9; 
  
  // 1 based indexing and 2 boundary cells
	vec2d Tn(Ny+2 + 1, vec1d(Nx+2 + 1, 0.0));
	vec2d Tn1(Ny+2 + 1, vec1d(Nx+2 + 1, 0.0));

	vec2d T2n(Ny2+2 + 1, vec1d(Nx2+2 + 1, 0.0));
	vec2d T2n1(Ny2+2 + 1, vec1d(Nx2+2 + 1, 0.0));

  // boundary condition values; counter clockwise from bottom edge
  const vec1d T_bound = {200.0, 200.0, 100.0, 200.0};

  vec3id interpolated_pts = interpolation_pts(Tn, T2n, del_mesh, del_mesh2, mesh2_origin);

  BC(Tn, T_bound);
  meshTransfer(Tn1, T2n1, interpolated_pts);

  // convergence constants
  double rms_T = 1.0;
  double rms_T2 = 1.0;

  int n = 0;
  while (rms_T > TOLERANCE || rms_T2 > TOLERANCE){
    gridSolve(Tn1, Tn, del_t, alpha, del_mesh);
    BC(Tn1, T_bound);
    meshTransfer(Tn1, T2n1, interpolated_pts);
    gridSolve(T2n1, T2n, del_t2, alpha, del_mesh2);

    rms_T = rmsT(Tn1, Tn);
    rms_T2 = rmsT(T2n1, T2n);
    n += 1;
    Tn = Tn1;
    T2n = T2n1;
    if (n % 5000 == 0){
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