#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <array>
#include <stack>
#include <queue>
#include <random>
#include <algorithm>
using namespace std;

using su2double = double;
const size_t SU2_BBOX_SIZE   = 4;  /*!< \brief Size of the bounding box array that is allocated for each element*/
constexpr su2double su2double_lowest = std::numeric_limits<su2double>::lowest();
constexpr su2double su2double_highest = std::numeric_limits<su2double>::max();

typedef vector<vector<pair<vector<size_t>, vector<double> > > > vec3id;
typedef vector<vector<double> > vec2d;
typedef vector<pair<vector<size_t>, vector<double> > > vec2id;
typedef vector<vector<size_t> > vec2i;
typedef vector<double> vec1d;
typedef vector<size_t> vec1i;

const double PI = 3.141592653589793238463;
const double TOLERANCE = 1024 * std::numeric_limits<double>::epsilon();
// const double TOLERANCE = 1e-6;
const size_t NS_MAXITER = 20;
const size_t BBOXSIZE = 4;

enum POINT_TYPE {
  UNUNSED = 0,
  CALCULATED = 1,          
  INTERPOLATION_DONOR = 2,
  INTERPOLATION_RECIEVER = 3,
  DONOR_BUFFER = 4,
  BC_SPECIFIED = 5,
};

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

class node_adt {
  public:
    size_t elementIndex;
    array<su2double, SU2_BBOX_SIZE> bBoxCoordinates; 
    node_adt * left;
    node_adt * right;
    
    node_adt(const size_t elementIdx, const array<su2double, SU2_BBOX_SIZE> &bBoxCoords) : elementIndex(elementIdx), bBoxCoordinates(bBoxCoords), left(nullptr), right(nullptr) {}
    virtual ~node_adt() = default;
};

class ADT{
  public:
    node_adt * root;
    size_t treeHeirarchy;
    ADT() : root(nullptr), treeHeirarchy(0) {}
    virtual ~ADT() = default;
  
    void insertNode(node_adt * node) {
      node_adt * current = root;
      node_adt * parent = nullptr;
      size_t elementHeirarchy = 0, i = 0;

      /* Traverse the tree to find the insertion point */
      while (current != nullptr) {
        parent = current;
        /* Heirarchy cycling -> (x_min, y_min, z_min, x_max, y_max, z_max)*/
        i = elementHeirarchy % 4; /* For 2d cycle with 4 values*/
        /* node with equal values are put on left branch*/
        // if (floatEqual(current->bBoxCoordinates[i], node->bBoxCoordinates[i])) {
        //   current = current->left;
        // } else 
        if (current->bBoxCoordinates[i] < node->bBoxCoordinates[i]) {
          current = current->right;
        } else {
          current = current->left;
        }
        elementHeirarchy = elementHeirarchy + 1;
      }
      if (parent == nullptr) {
        root = node;
      // } else if (floatEqual(parent->bBoxCoordinates[i], node->bBoxCoordinates[i])) {
      //   parent->left = node;
      // } else 
      } else if (parent->bBoxCoordinates[i] < node->bBoxCoordinates[i]) {
        parent->right = node;
      } else {
        parent->left = node;
      }

      current = node;
      treeHeirarchy = max(treeHeirarchy, elementHeirarchy);
      // cout << "Element added " << current->elementIndex << " at h = " << heirarchy << endl;
    }
    
    vector<size_t> searchADT(array<su2double, SU2_BBOX_SIZE> &testBBox) const {
      /*return a vector of indices of elements which intersect with the test bounding box*/
      size_t currHeirarchy = 0, i = 0;
      vector<size_t> intersectingBBox;
      node_adt * current = nullptr;
      bool intersect = true;

      stack<pair<node_adt *, size_t> > searchQ;
      searchQ.push(make_pair(root, 0));

      while (searchQ.empty() == false) {
        current = searchQ.top().first;
        currHeirarchy = searchQ.top().second;
        searchQ.pop();

        while (current != nullptr) {
          // if (current->elementIndex == );
          /* check intersection of current node*/
          intersect = true;
          for (unsigned short iDim = 0; iDim < SU2_BBOX_SIZE/2; iDim++) {
            // intersect = intersect && (testBBox[iDim] < current->bBoxCoordinates[iDim+SU2_BBOX_SIZE/2] || floatEqual(testBBox[iDim], current->bBoxCoordinates[iDim+SU2_BBOX_SIZE/2]));
            // intersect = intersect && (testBBox[iDim+SU2_BBOX_SIZE/2] > current->bBoxCoordinates[iDim] || floatEqual(testBBox[iDim+SU2_BBOX_SIZE/2], current->bBoxCoordinates[iDim]));
            intersect = intersect && (testBBox[iDim] <= current->bBoxCoordinates[iDim+SU2_BBOX_SIZE/2]);
            intersect = intersect && (testBBox[iDim+SU2_BBOX_SIZE/2] >= current->bBoxCoordinates[iDim]);
          }
          // intersect = intersect && testBBox[0] <= current->bBoxCoordinates[0+SU2_BBOX_SIZE/2];
          // intersect = intersect && testBBox[1] <= current->bBoxCoordinates[1+SU2_BBOX_SIZE/2];
          // intersect = intersect && testBBox[2] >= current->bBoxCoordinates[2-SU2_BBOX_SIZE/2];
          // intersect = intersect && testBBox[3] >= current->bBoxCoordinates[3-SU2_BBOX_SIZE/2];
          if (intersect) {
            intersectingBBox.push_back(current->elementIndex);
            // cout << "intersection found with element: " << current->elementIndex << endl;
          }
          i = currHeirarchy % 4;

          /*branching based on minimum coordinate*/
          if (i < SU2_BBOX_SIZE/2) {
            /*curr->left is always searched as the left has min coords lower than current which gives no info on intersection */
            searchQ.push(make_pair(current->left, currHeirarchy+1)); 
            // if (floatEqual(testBBox[i+SU2_BBOX_SIZE/2], current->bBoxCoordinates[i])) {
            //   current = current->right;
            // }
            /*Test _max < Current _min*/
            if (testBBox[i+SU2_BBOX_SIZE/2] < current->bBoxCoordinates[i]) {
              current = nullptr;
            }
            else {
              current = current->right;
            }
            // searchQ.push(make_pair(current->right, currHeirarchy)); 
          }
          /*branching based on maximum coordinate*/
          else{
            /*curr->right is always searched as the right has max coords higher than current which gives no info on intersection */
            searchQ.push(make_pair(current->right, currHeirarchy+1)); 
            // if (floatEqual(testBBox[i-SU2_BBOX_SIZE/2], current->bBoxCoordinates[i])) {
            //   current = current->left;
            // }
            /*Test _min > Current _max*/
            if (testBBox[i-SU2_BBOX_SIZE/2] > current->bBoxCoordinates[i] ) {
              current = nullptr;
            }
            else {
              current = current->left;
            }
            // searchQ.push(make_pair(current->left, currHeirarchy)); 
          }
          // current = nullptr;
          currHeirarchy = currHeirarchy + 1;
        }
      }
      if (intersectingBBox.size() == 0) {
        // cout << "No intersecting element BBox found." << endl;
      }
      return intersectingBBox;
    }

    /* Level order output ADT */ 
    void printLevelOrder() const {
      if (root == nullptr) {
        cerr << "Empty ADT" << endl;
        return;
      }
      queue< node_adt * > q;
      q.push(root);

      while (q.empty() == false) {
        // node_adt * node = q.front();
        // cout << node->elementIndex << " ";
        // q.pop();
        // if (node->left != nullptr)
        //     q.push(node->left);
        // if (node->right != nullptr)
        //     q.push(node->right);
        int count = q.size();
        while (count > 0){
          node_adt *node = q.front();
          cout << node->elementIndex << " ";
          q.pop();
          if (node->left != NULL)
              q.push(node->left);
          if (node->right != NULL)
              q.push(node->right);
          count--;
        }
        cout << endl;
      }
    }

    // Function to generate DOT syntax for the BST (helper function)
    std::string generateDot(node_adt* node) const {
      if (node == nullptr) {
        return "";
      }
      std::string dot;
      std::string nodeName = std::to_string(node->elementIndex);
      dot += "\"" + nodeName + "\"";
      dot += "[label=\"" + nodeName + "\"]\n";

      if (node->left != nullptr) {
        dot += "\"" + nodeName + "\"" + " -> \"" + std::to_string(node->left->elementIndex) + "\"";
        dot += "[color=red]\n";
        dot += generateDot(node->left);
      }
      else if (node->right != nullptr){
        dot += "\"" + nodeName + "_NULL\" ";
        dot += "[shape=point]\n";
        dot += "\"" + nodeName + "\"" + " -> \"" + nodeName + "_NULL\"";
        dot += "[color=red]\n";
      }
      if (node->right != nullptr) {
        dot += "\"" + nodeName + "\"" + " -> \"" + std::to_string(node->right->elementIndex) + "\"";
        dot += "[color=blue]\n";
        dot += generateDot(node->right);
      }
      else if (node->left != nullptr){
        dot += "\"" + nodeName + "_NULL\" ";
        dot += "[shape=point]\n";
        dot += "\"" + nodeName + "\"" + " -> \"" + nodeName + "_NULL\"";
        dot += "[color=blue]\n";
      }
      return dot;
    }

    // Generate DOT syntax and write to a file
    void writeDotToFile(const std::string& filename) const {
      std::ofstream file(filename);
      if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
      }
      file << "digraph BST {" << std::endl;
      file << generateDot(root);
      file << "}" << std::endl;
      file.close();
    }    
};

class StructuredMesh{
  /* 
  j, i based indexing of mesh so point coordinates are not required
  Points are numberered row wise aligned with local x coordinate i.e. 0, 1 ... Nx-1 is the bottom-most row followed by Nx, Nx+1 ... 2*Nx-1
  Elements are also numbered in this way i.e. teh bottom-most row is  0, 1, ... Nx-2
  */
  public:
    /*defines length (x), width(y), anti-clockwise rotation, origin*/
    double L, W, theta, xO, yO;
    double dx, dy;
    /* number of divisions in x and y directions*/
    size_t Nx, Ny, numberOfElements, numberOfPoints,dimension;
    
    /* defines the type of point: unused = 0, calculated = 1, interpolation donor = 2, interpolation reciever = 3 */
    vector<unsigned short> pointType;

    /* stores the donor interpolation stencil (point numbers are stored) from the other overset mesh for all the elements if it exists*/
    vector<vector<size_t> > interpolationStencil;
    /* stores the interpolation coefficients corresponding to the donor points in interpolationStencil */
    vector<vector<double> > interpolationCoeffs;

    vector<vector<double> > pointCoordinates;
    vector<vector<size_t> > neighborPointsOfPoint;
    vector<vector<size_t> > elementConnectivity;
    vector<vector<double> > elementBBox;

    /* bottom, right, top, left*/
    vector<vector<size_t> > boundaryPoints;
    /* bottom, right, top, left*/
    vector<double> boundaryCondition;
  
    ADT adtBoundingBox;

    /* variables at nodes */
    vector<double> T;
    vector<double> Tn1;

    StructuredMesh(double xOrigin, double yOrigin, double axisTheta, double length, double width, size_t _Nx, size_t _Ny) {
      dimension = 2;
      L = length, W = width, xO = xOrigin, yO = yOrigin, theta = axisTheta;
      Nx = _Nx, Ny = _Ny;
      numberOfPoints = Nx * Ny;
      numberOfElements = (Nx-1) * (Ny-1);
      dx = L / Nx;
      dy = W / Ny;

      /* Point coordinates */
      pointCoordinates.resize(dimension);
      for (unsigned short iDim = 0; iDim < dimension; ++iDim) {
        pointCoordinates[iDim].resize(numberOfPoints);
      }
      pointType.resize(numberOfPoints, CALCULATED);
      size_t iPoint = numberOfPoints;
      for (size_t j = 0; j < Ny; ++j) {
        for (size_t i = 0; i < Nx; i++) {
          iPoint = GetPointNumber(i, j);
          // pointCoordinates[0][iPoint] = i * dx; /* local x axis*/
          // pointCoordinates[1][iPoint] = j * dy; /* local y axis*/
          pointCoordinates[0][iPoint] = xO + (i*dx) * cos(theta) - (j*dy)* sin(theta); /* global x axis*/
          pointCoordinates[1][iPoint] = yO + (i*dx) * sin(theta) + (j*dy)* cos(theta); /* global y axis*/
        }
      }

      /* Element connectivity and it's bounding box */
      elementConnectivity.resize(numberOfElements);
      elementBBox.resize(numberOfElements);
      size_t iElement = numberOfElements;
      for(size_t j = 0; j < Ny-1; ++j) {
        for(size_t i = 0; i < Nx-1; ++i) {
          iElement = GetElementNumber(i, j);
          // elementConnectivity[iElement].assign({iElement, iElement+1, iElement+1+Nx, iElement+Nx});
          elementConnectivity[iElement].assign({j*Nx + i, j*Nx + i+1, (j+1)*Nx + i+1, (j+1)*Nx + i});
          
          /* Element cartesian aligned Bounding Box */
          // elementBBox[iElement].assign({i*dx, j*dy, (i+1)*dx, (j+1)*dy}); /* {Xmin, Ymin, Xmax, Ymax}*/
          elementBBox[iElement].assign({su2double_highest, su2double_highest, su2double_lowest, su2double_lowest}); /* {Xmin, Ymin, Xmax, Ymax}*/
          for (size_t iPoint = 0; iPoint < elementConnectivity[iElement].size(); ++iPoint){
              double x = pointCoordinates[0][elementConnectivity[iElement][iPoint]];
              double y = pointCoordinates[1][elementConnectivity[iElement][iPoint]];
              elementBBox[iElement][0] = min(elementBBox[iElement][0], x);
              elementBBox[iElement][1] = min(elementBBox[iElement][1], y);
              elementBBox[iElement][2] = max(elementBBox[iElement][2], x);
              elementBBox[iElement][3] = max(elementBBox[iElement][3], y);
          }
        }
      }

      /* Boundary Surfaces*/
      boundaryPoints.resize(4); /* bottom, right, top, left*/
      for(size_t iBoundPoint = 0; iBoundPoint < Nx; iBoundPoint++){
        boundaryPoints[0].push_back(iBoundPoint); /*bottom*/
        boundaryPoints[2].push_back((Ny-1)*Nx + iBoundPoint); /*top*/
      }      
      for(size_t iBoundPoint = 0; iBoundPoint < Ny; iBoundPoint++){
        boundaryPoints[3].push_back(iBoundPoint*Nx); /*left*/
        boundaryPoints[1].push_back(Nx*iBoundPoint + (Nx-1)); /*right*/
      }
      reverse(boundaryPoints[2].begin(), boundaryPoints[2].end());
      reverse(boundaryPoints[3].begin(), boundaryPoints[3].end());

      /* Neighbour points of point */
      neighborPointsOfPoint.resize(numberOfPoints);
      for (size_t j = 0; j < Ny-1; ++j) {
        for (size_t i = 0; i < Nx-1; ++i) {
          iElement = GetElementNumber(i, j);
          size_t pointA, pointB;
          for(size_t iPoint = 0; iPoint < elementConnectivity[iElement].size(); iPoint++) {
            pointA = elementConnectivity[iElement][iPoint];
            pointB = elementConnectivity[iElement][(iPoint+1)%4];
            neighborPointsOfPoint[pointA].push_back(pointB);
            neighborPointsOfPoint[pointB].push_back(pointA);
          }
        }
      }
      /* remove duplicates from the neighboring point lists*/
      iPoint = numberOfPoints;
      vector<size_t>::iterator vecIt;
      for (size_t j = 0; j < Ny; ++j) {
        for (size_t i = 0; i < Nx; ++i) {
          iPoint = GetPointNumber(i, j);
          /* sort neighboring points for each point */
          sort(neighborPointsOfPoint[iPoint].begin(), neighborPointsOfPoint[iPoint].end());
          
          /* uniquify list of neighboring points */
          vecIt = unique(neighborPointsOfPoint[iPoint].begin(), neighborPointsOfPoint[iPoint].end());

          /* adjust size of vector */
          neighborPointsOfPoint[iPoint].resize(vecIt - neighborPointsOfPoint[iPoint].begin());
        }
      }
      
      /*resize/initialize all vectors*/
      T.resize(numberOfPoints, -1);
      Tn1.resize(numberOfPoints, -1);

      interpolationStencil.resize(numberOfPoints);
      interpolationCoeffs.resize(numberOfPoints);
      GenerateADT();
    };
    
    inline size_t GetElementNumber (size_t i, size_t j) const {return j*(Nx-1) + i;}
    inline size_t GetPointNumber (size_t i, size_t j) const {return j*Nx + i;}

    virtual ~StructuredMesh() = default;

    void GenerateADT() {
      array<su2double, SU2_BBOX_SIZE> bboxCoords{};
      
      /* modifying node insertion to balance the tree */
      vector<size_t> elementOrderADT;
      elementOrderADT.resize(numberOfElements);
      iota(elementOrderADT.begin(), elementOrderADT.end(), 0);
      // shuffle(elementOrderADT.begin(), elementOrderADT.end(), std::mt19937{std::random_device{}()});
      shuffle(elementOrderADT.begin(), elementOrderADT.end(), std::mt19937{12337});

      for (size_t randomIndex = 0; randomIndex < numberOfElements; ++randomIndex) {
        size_t LocalIndex = elementOrderADT[randomIndex];
        for (unsigned short i = 0; i < SU2_BBOX_SIZE; ++i) {
          bboxCoords[i] = elementBBox[LocalIndex][i];
        }
        node_adt * tempNode = new node_adt(LocalIndex, bboxCoords);
        adtBoundingBox.insertNode(tempNode);
      }
      cout << "Max heirarchy : " << adtBoundingBox.treeHeirarchy << endl;
    }

    array<su2double, SU2_BBOX_SIZE> GetPointBBox (const size_t pointNumber) const {
      /* small value to generate a bounding box centered around a point */
      su2double epsilon = std::numeric_limits<su2double>::epsilon();
      array<su2double, 3> Coords {0.0, 0.0, 0.0};
      array<su2double, SU2_BBOX_SIZE> pointBBoxCoords{};
      for (unsigned short iDim = 0; iDim < dimension; ++iDim) {
          Coords[iDim] = pointCoordinates[iDim][pointNumber];
          pointBBoxCoords[iDim] = Coords[iDim] - epsilon; 
          pointBBoxCoords[iDim+2] = Coords[iDim] + epsilon; 
        }
      return pointBBoxCoords;
    }

    void MarkOversetBoundaryDonors(const StructuredMesh& oversetMesh) {
      /* does a type of holecutting / removal of coarse cells where overlapping finer cells are present*/

      /* iterate over all the boundaries of overset mesh to mark donors */
      for (size_t iBound = 0; iBound < oversetMesh.boundaryPoints.size(); ++iBound) {
        /* iterate over all points in iBound*/
        for (size_t iBoundPoint = 0; iBoundPoint < oversetMesh.boundaryPoints[iBound].size(); ++iBoundPoint){
          size_t point = oversetMesh.boundaryPoints[iBound][iBoundPoint];
          auto pointBBox = oversetMesh.GetPointBBox(point);
          const auto BBox = adtBoundingBox.searchADT(pointBBox);
          if (BBox.size() == 0) {
            cerr << "Point: "<< point << " No interpolation stencil found. What to do ??" << endl;
          }
          else {
            if (BBox.size() != 1) {
            cerr << "Point : "<< point << " Multiple interpolation stencil found [" << BBox.size() << "]. What to do ??" << endl;
            }
            for (auto iPoint : elementConnectivity[BBox[0]]) {
              // cout << "Marking iPoint = " << iPoint << " as donor " << endl;
              pointType[iPoint] = 2;
              /* mark all neighbors as required (for complete stencil of donor)*/
              for (auto neighborPoint : neighborPointsOfPoint[iPoint]) {
                if (pointType[neighborPoint] == 2) {continue;}
                pointType[neighborPoint] = 4; /*donor buffer*/
              }
            }
          }
        }
      }

      /*check removability of each point -> if neighbours can be interpolation cells and cell is bigger than overlapping cell*/
      for (size_t iPoint = 0; iPoint < numberOfPoints; ++iPoint) {
        /* skip interpolation donors and it's buffer */
        if (pointType[iPoint] == INTERPOLATION_DONOR || pointType[iPoint] == DONOR_BUFFER) {
          continue;
        }
        
        /* TODO: check if donors/overlapping cells are finer/smaller than iPoint. Assumed for now.*/
        bool neighborsCanBeInterpolated = true;
        for (auto neighborPoint : neighborPointsOfPoint[iPoint]) {
          if (pointType[neighborPoint] == INTERPOLATION_DONOR || pointType[neighborPoint] == DONOR_BUFFER) {neighborsCanBeInterpolated = false; break;}
          auto pointBBox = GetPointBBox(neighborPoint);
          const auto BBox = oversetMesh.adtBoundingBox.searchADT(pointBBox);
          if (BBox.size() == 0) {
            neighborsCanBeInterpolated = false;
          } 
        }
        if (neighborsCanBeInterpolated) {
          pointType[iPoint] = UNUNSED; // iPoint becomes unused
          for (auto neighborPoint : neighborPointsOfPoint[iPoint]) {
            if (pointType[neighborPoint] != CALCULATED) {continue;}
            pointType[neighborPoint] = INTERPOLATION_RECIEVER; // neighbors are marked as interpolated
          }
        }
      }

    }

    /* pType : 3 = interpolated, 5 = boundary condition specified */
    void MarkBoundaryPointType(const unsigned short ptType) {

      /* marks all the boundaries (for component meshes marked as intepolated. Assumed component meshes dont overlap each other nor do they intersect boundaries) */
      for (size_t iBound = 0; iBound < boundaryPoints.size(); ++iBound) {
        /* iterate over all points in iBound*/
        for (size_t iBoundPoint = 0; iBoundPoint < boundaryPoints[iBound].size(); ++iBoundPoint){
          size_t point = boundaryPoints[iBound][iBoundPoint];
          pointType[point] = ptType;
        }
      }
    }

    void StoreInterpolationStencil(const StructuredMesh& oversetMesh) {
      for (size_t iPoint = 0; iPoint < numberOfPoints; ++iPoint) {
        /* skip points which are not interpolated*/
        if (pointType[iPoint] != INTERPOLATION_RECIEVER ) {
          continue;
        }
        auto pointBBox = GetPointBBox(iPoint);
        const auto BBox = oversetMesh.adtBoundingBox.searchADT(pointBBox);
        if (BBox.size() == 0) {
          cerr << " No donor found for interpolation marked point." << endl;
        } else {
          // interpolationStencil[iPoint].push_back(BBox[0]);
          auto donorStencil = oversetMesh.elementConnectivity[BBox[0]];
          vector<vector<double> > donorCoords (dimension);
          for (size_t iDonor = 0; iDonor < donorStencil.size(); ++iDonor){
            interpolationStencil[iPoint].push_back(donorStencil[iDonor]);
            for (unsigned short iDim = 0; iDim < dimension; ++iDim){
              donorCoords[iDim].push_back(oversetMesh.pointCoordinates[iDim][donorStencil[iDonor]]);
            }
          }
          vector<double> receiverCoords;
          for(unsigned short iDim = 0; iDim < dimension; ++iDim) {
            receiverCoords.push_back(pointCoordinates[iDim][iPoint]);
          }
          interpolationCoeffs[iPoint] = interpolation_coeff(donorCoords, receiverCoords);
        }
      }
    }

    void WriteTxtPointType(string filename = "pointType.txt") const{
      ofstream pointTypeFile(filename);
      if (!pointTypeFile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
      }
      for (size_t iPoint = 0; iPoint < numberOfPoints; iPoint++) {
        pointTypeFile << pointType[iPoint] << " ";
        // if (iPoint % 100 == 0 && iPoint != 0) {pointType << endl;}
      }
      pointTypeFile.close();
    }

    void MeshSolve(const double del_t, const double alpha, const StructuredMesh& oversetMesh) {
      /* loop over all points even interpolated and boundary and solve for calculated and interpolated for receiver*/
      for (size_t j = 0; j < Ny; j++) {
        for (size_t  i = 0; i < Nx; ++i) {
          size_t point = GetPointNumber(i ,j), left = GetPointNumber(i-1 ,j), right = GetPointNumber(i+1 ,j), top = GetPointNumber(i ,j+1), bottom = GetPointNumber(i ,j-1);
          /* solve for calculated type i.e pointType = 1 */
          if (pointType[point] == CALCULATED || pointType[point] == INTERPOLATION_DONOR || pointType[point] == DONOR_BUFFER) {
            Tn1[point] = T[point] + del_t*alpha* ((T[left] - 2*T[point] + T[right])/pow(dx, 2) + (T[bottom] - 2*T[point] + T[top])/pow(dy, 2));
          }
          /* interpolate values for receiver type i.e pointType = 3 */
          else if (pointType[point] == INTERPOLATION_RECIEVER) {
            Tn1[point] = 0.0;
            for(unsigned short iDonor = 0; iDonor < interpolationStencil[point].size(); ++iDonor) {
              auto donorPointNumber = interpolationStencil[point][iDonor];
              Tn1[point] = Tn1[point] + oversetMesh.Tn1[donorPointNumber] * interpolationCoeffs[point][iDonor];
            }
          }
          else {
            continue;
          }
        }
      }
    }

    void ApplySpecifiedBoundaryCondition (vector<double> BC) {
      boundaryCondition = BC;
      for (unsigned short iBound = 0; iBound < boundaryPoints.size(); ++iBound) {
        for (auto iBoundPoint : boundaryPoints[iBound]) {
          T[iBoundPoint] = boundaryCondition[iBound];
          Tn1[iBoundPoint] = boundaryCondition[iBound];
        }
      }
    }

    double RMS() const {
      double rms = 0.0;
      for (size_t j = 0; j < Ny; ++j) {
        for (size_t i = 0; i < Nx; ++i) {
          size_t point = GetPointNumber(i ,j);
          if (pointType[point] == UNUNSED) {continue;}
          rms += pow(Tn1[point] - T[point], 2);
        }
      }
      rms = rms / (Nx * Ny);
      rms = pow(rms, 0.5);
      return rms;
    }

    void SaveToTxt(string filename) const {
      // saving results to file
      std::ofstream output_file(filename);
      for(size_t j = 0; j < Ny; j++){
        for(size_t i = 0; i < Nx; i++){
          size_t point = GetPointNumber(i ,j);
          output_file << std::fixed << std::setprecision(std::numeric_limits<double>::digits10+2) << Tn1[point]<<" ";
        }
        output_file << endl;
      }
    }
};

void oversetSolver(){
  /* xOrigin, yOrigin, theta, lenght, width, Nx, Ny*/
  StructuredMesh bgMesh(0.0, 0.0, 0.0, 1.0, 1.0, 50, 50);
  StructuredMesh compMesh(0.445, 0.245, 45 * PI/180.0, 0.4, 0.4, 20, 20);
  
  // bgMesh.adtBoundingBox.writeDotToFile("bgADT.txt");

  bgMesh.MarkOversetBoundaryDonors(compMesh);
  bgMesh.MarkBoundaryPointType(BC_SPECIFIED);
  compMesh.MarkBoundaryPointType(INTERPOLATION_RECIEVER);

  bgMesh.StoreInterpolationStencil(compMesh);
  compMesh.StoreInterpolationStencil(bgMesh);

  bgMesh.WriteTxtPointType("bgPtType.txt");
  compMesh.WriteTxtPointType("compPtType.txt");

  /* numerical constants for PDE */
  const double k = 16.2, rho = 7750, cp = 500.0;
  const double alpha = k / (rho * cp);
  const double del_t = (1 / (2 * alpha * (1/pow(min(bgMesh.dx, compMesh.dx), 2) + 1/pow(min(bgMesh.dy, compMesh.dy),2)))) * 0.9; // taking 90 % of the upper limit as per von Neumann stability analysis

  /* counter clockwise from bottom edge : bottom, right, top, left */
  bgMesh.ApplySpecifiedBoundaryCondition(vector<double> {200.0, 200.0, 100.0, 200.0});

  double rmsBg = 1.0, rmsComp = 1.0;
  int n = 0;
  while ((rmsBg > TOLERANCE || rmsComp > TOLERANCE)){

    bgMesh.MeshSolve(del_t, alpha, compMesh);
    compMesh.MeshSolve(del_t, alpha, bgMesh);

    rmsBg = bgMesh.RMS();
    rmsComp = compMesh.RMS();

    n = n+1;

    bgMesh.T.swap(bgMesh.Tn1);
    compMesh.T.swap(compMesh.Tn1);

    if (n % 5000 == 0) {
      cout << "Iterations = " << n << " RMS_Bg = " <<  rmsBg << " RMS_Comp = " <<  rmsComp << "\n";
    }
  }
  cout << "Final Errors  = " << " RMS_Bg = " <<  rmsBg << " RMS_Comp = " <<  rmsComp << "\n";
  cout << "Total Iterations = " << n << endl;

  bgMesh.SaveToTxt("./num_mesh_bg.txt");
  compMesh.SaveToTxt("./num_mesh_comp.txt");
}

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

double analyticalError(const vec2d &T, const double L, const double W, const vec1d del_mesh, const vec1d mesh_origin, string filename) {
  const size_t Nx = T.size()-1-2, Ny = T[0].size()-1-2;
  const double del_x = del_mesh[0], del_y = del_mesh[1];
  const double x0 = mesh_origin[0], y0 = mesh_origin[1];
  const double sinTheta = mesh_origin[2], cosTheta = mesh_origin[3];
  vec2d xCoord (Ny+2+1, vec1d(Nx+2+1, 0.0));
  vec2d yCoord (Ny+2+1, vec1d(Nx+2+1, 0.0));
  for (size_t j = 2; j <= Ny+1; ++j) {
    for (size_t i = 2; i <= Nx+1; ++i) {
      double x_local = (i-2)*del_x + del_x/2;
      double y_local = (j-2)*del_y + del_y/2;
      xCoord[j][i] = x0 + (x_local) * cosTheta - (y_local) * sinTheta;
      yCoord[j][i] = y0 + (x_local) * sinTheta + (y_local) * cosTheta;
    }
  }

  double error = 0.0;
  double fact = 0.0;
  const double T2 = 100.0, T1 = 200.0;

  vec2d T_ana(Ny+2 + 1, vec1d(Nx+2 + 1, 0.0));
  for (size_t j = 2; j <= Ny+1; ++j) {
    for (size_t i = 2; i <= Nx+1; ++i) {
      fact = 0.0;
      for (size_t k = 1; k <= 100; ++k) {
        fact = fact + (pow(-1, k+1) + 1)/k * sin(k*PI*xCoord[j][i]/L) * sinh(k*PI*yCoord[j][i]/L) / sinh(k*PI*W/L);
      }
      T_ana[j][i] = T1 + (T2-T1) * 2/PI * fact;
      // double temp_ana = T1 + (T2-T1) * 2/PI * fact;
      error = error + pow(T_ana[j][i] - T[j][i], 2.0);
    }
  }
  vecToText(T_ana, filename);
  // vecToText(xCoord, "./xCoord.txt");
  // vecToText(yCoord, "./yCoord.txt");
  return pow(error, 0.5);
  // return error;
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
  return pow(rms_T, 0.5);
}

void solver() {
  // Grid discretization AND Domain 
  const size_t Nx = 10, Ny = 10;
  const double L = 5.0, W = 5.0;

  const double x2 = L/3, y2 = W/3, theta = -25.0 * PI/180.0;
  const double L2 = 1.5, W2 = 1.5;
  const size_t Nx2 = L2/L*Nx, Ny2 = W2/W*Ny;

  const double del_x = L / Nx, del_y = W / Ny;
  const double del_x2 = L2 / Nx2, del_y2 = W2 / Ny2;
  const double cosTheta = cos(theta), sinTheta = sin(theta);
  
  const vec1d mesh1_origin = {0.0, 0.0, 0.0, 1.0};
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
    if (n % 10000 == 0){
      cout << "Iterations = " << n << " RMS_T = " <<  rms_T << " RMS_T2 = " <<  rms_T2 << "\n";
    }
  }

  cout << "Mesh 1 error - >" << analyticalError(Tn1, L, W, del_mesh, mesh1_origin, "./ana_mesh1.txt") << endl;
  cout << "Mesh 2 error - >" << analyticalError(T2n1, L, W, del_mesh2, mesh2_origin, "./ana_mesh2.txt") << endl;
  cout << "Total Iterations = " << n << endl;
  
  vecToText(Tn1, "./num_mesh1.txt");
  vecToText(T2n1, "./num_mesh2.txt");
}

int main() {
  // solver();
  oversetSolver();
  return 0;
}