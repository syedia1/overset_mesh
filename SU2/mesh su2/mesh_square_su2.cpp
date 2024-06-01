/*!
 * \file square.cpp
 * \brief C++ script for creating a .su2 mesh of a simple square domain.
 * \author Thomas D. Economon
 * \version 1.0.
 */

#include <fstream>

using namespace std;

void genMesh(const int nNode, const int mNode, const double L, const double W, const double x_origin, const double y_origin){
  /*--- Local variables ---*/
  int KindElem, KindBound;
  int iElem;
  // int nNode, mNode;
  int iNode, jNode, nPoint, iPoint, jPoint, kPoint;
  ofstream Mesh_File;
  
  /*--- Set the VTK type for the interior elements and the boundary elements ---*/
  KindElem  = 5; // Triangle
  KindBound = 3; // Line
  
  // /*--- Store the number of nodes and output mesh filename ---*/
  // nNode     = 4;
  // mNode     = 4;

  // /* Shifting origin*/
  // double x_origin = 0.0, y_origin = 0.0; 

  // /*Dimensions in x and y i.e L and W*/
  // double L = 1.0, W = 1.0;
  
  /*--- Open .su2 grid file ---*/
	Mesh_File.precision(15);
	Mesh_File.open("square_"+to_string(nNode-1)+"x"+to_string(mNode-1)+".su2", ios::out);
  
  /*--- Write the dimension of the problem and the number of interior elements ---*/
  Mesh_File << "%" << endl;
  Mesh_File << "% Problem dimension" << endl;
  Mesh_File << "%" << endl;
  Mesh_File << "NDIME= 2" << endl;
  Mesh_File << "%" << endl;
  Mesh_File << "% Inner element connectivity" << endl;
  Mesh_File << "%" << endl;
  Mesh_File << "NELEM= " <<  2*(nNode-1)*(mNode-1) << endl;

  /*--- Write the element connectivity ---*/
  iElem = 0;
  for (jNode = 0; jNode < mNode-1; jNode++) {
    for (iNode = 0; iNode < nNode-1; iNode++) {
      iPoint = jNode*nNode + iNode;
      jPoint = jNode*nNode + iNode + 1;
      kPoint = (jNode + 1)*nNode + iNode;
      Mesh_File << KindElem << "\t" << iPoint << "\t" << jPoint << "\t" << kPoint << "\t" << iElem << endl;
      iElem ++;
      iPoint = jNode*nNode + (iNode + 1);
      jPoint = (jNode + 1)*nNode + (iNode + 1);
      kPoint = (jNode + 1)*nNode + iNode;
      Mesh_File << KindElem << "\t" << iPoint << "\t" << jPoint << "\t" << kPoint << "\t" << iElem << endl;
      iElem++;
    }
  }
  
  /*--- Compute the number of nodes and write the node coordinates ---*/
  nPoint = (nNode)*(mNode);
  Mesh_File << "%" << endl;
  Mesh_File << "% Node coordinates" << endl;
  Mesh_File << "%" << endl;
  Mesh_File << "NPOIN= " << nNode*mNode << endl;
  iPoint = 0;
  for (jNode = 0; jNode < mNode; jNode++) {
    for (iNode = 0; iNode < nNode; iNode++) {
      Mesh_File << x_origin + ((double)iNode)/((double)(nNode-1)) << "\t" << y_origin + ((double)jNode)/((double)(mNode-1)) << "\t" << iPoint << endl;
      iPoint++;
    }
  }
  
  /*--- Write the header information for the boundary markers ---*/
  Mesh_File << "%" << endl;
  Mesh_File << "% Boundary elements" << endl;
  Mesh_File << "%" << endl; 
  Mesh_File << "NMARK= 4" << endl;
  
  /*--- Write the boundary information for each marker ---*/
  Mesh_File << "MARKER_TAG= lower" << endl;
  Mesh_File << "MARKER_ELEMS= "<< (nNode-1) << endl;
  for (iNode = 0; iNode < nNode-1; iNode++) {
    Mesh_File << KindBound << "\t" << iNode << "\t" << (iNode + 1) << endl;
  }
  Mesh_File << "MARKER_TAG= right" << endl;
  Mesh_File << "MARKER_ELEMS= "<< (mNode-1) << endl;
  for (jNode = 0; jNode < mNode-1; jNode++) {
    Mesh_File << KindBound << "\t" << jNode*nNode + (nNode - 1) << "\t" << (jNode + 1)*nNode + (nNode - 1) << endl;
  }
  Mesh_File << "MARKER_TAG= upper" << endl;
  Mesh_File << "MARKER_ELEMS= "<< (nNode-1) << endl;
  for (iNode = 0; iNode < nNode-1; iNode++) {
    Mesh_File << KindBound << "\t" << (nNode*mNode - 1) - iNode << "\t" << (nNode*mNode - 1) - (iNode + 1) << endl;
  }
  Mesh_File << "MARKER_TAG= left" << endl;
  Mesh_File << "MARKER_ELEMS= "<< (mNode-1) << endl;
  for (jNode = mNode-2; jNode > mNode-4; jNode--) {
    Mesh_File << KindBound << "\t" << (jNode + 1)*nNode << "\t" << jNode*nNode << endl;
  }
  
  /*--- Close the mesh file and exit ---*/
  Mesh_File.close();
}

int main() {
  int n1 = 11, m1 = 11;
  double l1 = 1.0, w1 = 1.0;
  int n2 = 4, m2 = 4;

  genMesh(n1, m1, l1, w1, 0.0, 0.0);
  genMesh(n2, m2, ((double)(n2-1)/double(n1-1)), ((double)(m2-1)/double(m1-1)), 0.23, 0.36);
  return 0;
}
