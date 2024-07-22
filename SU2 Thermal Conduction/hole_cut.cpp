#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <array>
#include <cassert>
#include <iomanip> 
#include <queue> 
#include <numeric> 
#include <random>
#include <algorithm>
#include <limits>
#include <unordered_set>
#include <set>
#include <stack>


using std::vector, std::array, std::queue, std::stack;
using std::unordered_set, std::set;
using std::pair, std::make_pair;
using std::iota, std::shuffle;
using std::string, std::to_string;
using std::ifstream, std::istringstream, std::ofstream;
using std::cout, std::cin, std::cerr, std::endl;
using std::min, std::max, std::sort;
using su2double = double;

/**
 * @brief Tests if the 2 numbers are approximately equal, with the floating point rounding error taken into account.
 *
 * @param a Compared value.
 * @param b Reference value.
 * @return 1: Approximately equal. 0: Not approximately equal.
 */
inline bool floatEqual(double a, double b)
{
    double _a = fabs(a), _b = fabs(b);
    return fabs(a - b) <= ((_a < _b ? _b : _a) * std::numeric_limits<su2double>::epsilon());
}
constexpr su2double su2double_lowest = std::numeric_limits<su2double>::lowest();
constexpr su2double su2double_highest = std::numeric_limits<su2double>::max();

const size_t SU2_CONN_SIZE   = 10;  /*!< \brief Size of the connectivity array that is allocated for each element*/
const size_t SU2_BBOX_SIZE   = 4;  /*!< \brief Size of the bounding box array that is allocated for each element*/

enum GEO_TYPE {
  VERTEX = 1,         /*!< \brief VTK nomenclature for defining a vertex element. */
  LINE = 3,           /*!< \brief VTK nomenclature for defining a line element. */
  TRIANGLE = 5,       /*!< \brief VTK nomenclature for defining a triangle element. */
  QUADRILATERAL = 9,  /*!< \brief VTK nomenclature for defining a quadrilateral element. */
  TETRAHEDRON = 10,   /*!< \brief VTK nomenclature for defining a tetrahedron element. */
  HEXAHEDRON = 12,    /*!< \brief VTK nomenclature for defining a hexahedron element. */
  PRISM = 13,         /*!< \brief VTK nomenclature for defining a prism element. */
  PYRAMID = 14        /*!< \brief VTK nomenclature for defining a pyramid element. */
};
constexpr unsigned short N_ELEM_TYPES = 7;           /*!< \brief General output & CGNS defines. */

constexpr unsigned short N_POINTS_LINE = 2;          /*!< \brief General output & CGNS defines. */
constexpr unsigned short N_POINTS_TRIANGLE = 3;      /*!< \brief General output & CGNS defines. */
constexpr unsigned short N_POINTS_QUADRILATERAL = 4; /*!< \brief General output & CGNS defines. */
constexpr unsigned short N_POINTS_TETRAHEDRON = 4;   /*!< \brief General output & CGNS defines. */
constexpr unsigned short N_POINTS_HEXAHEDRON = 8;    /*!< \brief General output & CGNS defines. */
constexpr unsigned short N_POINTS_PYRAMID = 5;       /*!< \brief General output & CGNS defines. */
constexpr unsigned short N_POINTS_PRISM = 6;         /*!< \brief General output & CGNS defines. */
constexpr unsigned short N_POINTS_MAXIMUM = 8;       /*!< \brief Max. out of the above, used for static arrays, keep it up to date. */

constexpr unsigned short N_FACES_LINE = 1;           /*!< \brief General output & CGNS defines. */
constexpr unsigned short N_FACES_TRIANGLE = 3;       /*!< \brief General output & CGNS defines. */
constexpr unsigned short N_FACES_QUADRILATERAL = 4;  /*!< \brief General output & CGNS defines. */
constexpr unsigned short N_FACES_TETRAHEDRON = 4;    /*!< \brief General output & CGNS defines. */
constexpr unsigned short N_FACES_PYRAMID = 5;        /*!< \brief General output & CGNS defines. */
constexpr unsigned short N_FACES_PRISM = 5;          /*!< \brief General output & CGNS defines. */
constexpr unsigned short N_FACES_HEXAHEDRON = 6;     /*!< \brief General output & CGNS defines. */
constexpr unsigned short N_FACES_MAXIMUM = 6;        /*!< \brief Max. out of the above, used for static arrays, keep it up to date. */

/*!
 * \brief Get the number of faces of the element.
 * \param[in] elementType - element type
 * \return number of faces
 */
inline unsigned short nFacesOfElementType(unsigned short elementType) {
  switch (elementType) {
    case LINE: return N_FACES_LINE;
    case TRIANGLE: return N_FACES_TRIANGLE;
    case QUADRILATERAL: return N_FACES_QUADRILATERAL;
    case TETRAHEDRON: return N_FACES_TETRAHEDRON;
    case HEXAHEDRON: return N_FACES_HEXAHEDRON;
    case PYRAMID: return N_FACES_PYRAMID;
    case PRISM: return N_FACES_PRISM;
    default: assert(false && "Invalid element type."); return 0;
  }
}

/*!
 * \brief Get the number of points of the element.
 * \param[in] elementType - element type
 * \return number of points
 */
inline unsigned short nPointsOfElementType(unsigned short elementType) {
  switch (elementType) {
    case LINE: return N_POINTS_LINE;
    case TRIANGLE: return N_POINTS_TRIANGLE;
    case QUADRILATERAL: return N_POINTS_QUADRILATERAL;
    case TETRAHEDRON: return N_POINTS_TETRAHEDRON;
    case HEXAHEDRON: return N_POINTS_HEXAHEDRON;
    case PYRAMID: return N_POINTS_PYRAMID;
    case PRISM: return N_POINTS_PRISM;
    default: assert(false && "Invalid element type."); return 0;
  }
}

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

class SU2Mesh {
  // protected:
  public:
    size_t dimension = 0;

    size_t numberOfLocalPoints = 0;
    vector<vector<su2double> > localPointCoordinates;
    
    size_t numberOfLocalElements = 0;
    vector<size_t> localVolumeElementConnectivity;
    
    size_t numberOfMarkers = 0;
    vector<string> markerNames;
    vector<vector<size_t> > surfaceElementConnectivity;
  
    /* layout is [[elementIndex, VTK type, xmin, ymin, zmin, xmax, ymax, zmax]]*/
    vector<su2double> localVolumeElementBoundingBox;
    ADT adtBoundingBox;

    vector<vector<size_t> > neighborPointsOfPoint;
    /* defines the type of point: unused = 0, calculated = 1, interpolation donor = 2, interpolation reciever = 3 */
    vector<unsigned short> localPointType;

  // public:
    SU2Mesh(const string& filename) {
      bool foundNDIME = false, foundNPOIN = false;
      bool foundNELEM = false, foundNMARK = false;
      string text_line;  

      ifstream mesh_file(filename);
      if (mesh_file.is_open()) { 
        while(getline(mesh_file, text_line)) {
          if (text_line.size() == 0){            /* skip empty line */
            continue;  
          }   
          else if (text_line[0] == '%') {       /* skip comments */ 
            continue;
          }
          else {
            /*--- Read the dimension of the problem ---*/
            if (!foundNDIME && text_line.find("NDIME=", 0) != string::npos) {
              text_line.erase(0, 6);
              dimension = atoi(text_line.c_str());
              foundNDIME = true;
              continue;
            }

            /*--- Read the points/nodes ---*/
            if (!foundNPOIN && text_line.find("NPOIN=", 0) != string::npos) {
              text_line.erase(0, 6);
              numberOfLocalPoints = atoi(text_line.c_str());

              localPointCoordinates.resize(dimension);
              for (unsigned short k = 0; k < dimension; k++) localPointCoordinates[k].reserve(numberOfLocalPoints);
              
              array<su2double, 3> Coords {0.0, 0.0, 0.0};
              for (size_t iPoint = 0; iPoint < numberOfLocalPoints; iPoint++) {
                getline(mesh_file, text_line);
                istringstream point_line(text_line);
                /* Store the coordinates more clearly. */
                point_line >> Coords[0];
                point_line >> Coords[1];
                if (dimension == 3) {
                  point_line >> Coords[2];
                }
                for (unsigned short iDim = 0; iDim < dimension; iDim++) {
                  localPointCoordinates[iDim].push_back(Coords[iDim]);
                }
              }

              foundNPOIN = true;
              continue;
            }

            /*--- Read the internal connectivity ---*/
            if (!foundNELEM && text_line.find("NELEM=", 0) != string::npos) {
              text_line.erase(0, 6);
              numberOfLocalElements = atoi(text_line.c_str());
              
              array<size_t, N_POINTS_HEXAHEDRON> connectivity{};

              for (size_t LocalIndex = 0; LocalIndex < numberOfLocalElements; ++LocalIndex) {
                getline(mesh_file, text_line);
                istringstream elem_line(text_line);

                unsigned short VTK_Type;
                elem_line >> VTK_Type;

                const auto nPointsElem = nPointsOfElementType(VTK_Type);

                for (unsigned short i = 0; i < nPointsElem; i++) {
                  elem_line >> connectivity[i];
                }
                localVolumeElementConnectivity.push_back(LocalIndex);
                localVolumeElementConnectivity.push_back(VTK_Type);
                for (unsigned short i = 0; i < N_POINTS_HEXAHEDRON; i++) {
                  localVolumeElementConnectivity.push_back(connectivity[i]);
                }
              }

              foundNELEM = true;
              continue;
            }

            /*--- Read the boundary/markers ---*/
            if (!foundNMARK && text_line.find("NMARK=", 0) != string::npos) {
              text_line.erase(0, 6);
              numberOfMarkers = atoi(text_line.c_str());;
              
              surfaceElementConnectivity.resize(numberOfMarkers);
              markerNames.resize(numberOfMarkers);
              array<unsigned long, N_POINTS_HEXAHEDRON> connectivity{};
              for (unsigned short iMarker = 0; iMarker < numberOfMarkers; ++iMarker) {
                getline(mesh_file, text_line);
                text_line.erase(0, 11);
                string::size_type position;
                
                /* removes white spaces in/around marker names */
                for (unsigned short iChar = 0; iChar < 20; iChar++) {
                  position = text_line.find(' ', 0);
                  if (position != string::npos) text_line.erase(position, 1);
                  position = text_line.find('\r', 0);
                  if (position != string::npos) text_line.erase(position, 1);
                  position = text_line.find('\n', 0);
                  if (position != string::npos) text_line.erase(position, 1);
                }
                
                markerNames[iMarker] = text_line;
                
                getline(mesh_file, text_line);
                text_line.erase(0, 13);
                size_t nElem_Bound = atoi(text_line.c_str());
                
                for (size_t iElem_Bound = 0; iElem_Bound < nElem_Bound; iElem_Bound++) {
                  getline(mesh_file, text_line);
                  istringstream bound_line(text_line);

                  unsigned short VTK_Type;
                  bound_line >> VTK_Type;
                  const auto nPointsElem = nPointsOfElementType(VTK_Type);
                  for (unsigned short i = 0; i < nPointsElem; i++) {
                    bound_line >> connectivity[i];
                  }

                  surfaceElementConnectivity[iMarker].push_back(0); /* magic number 0 ??*/
                  surfaceElementConnectivity[iMarker].push_back(VTK_Type);
                  for (unsigned short i = 0; i < N_POINTS_HEXAHEDRON; i++) {
                    surfaceElementConnectivity[iMarker].push_back(connectivity[i]);
                  }
                }
              }
              foundNMARK = true;
              continue;
            }

          }
        mesh_file.close();
        }
      }
      else { 
        cerr << "Unable to open file -> "<< filename << endl; 
      } 

      /*create vector of neighbors of each point */
      neighborPointsOfPoint.resize(numberOfLocalPoints);
      for (size_t iElement = 0; iElement < numberOfLocalElements; ++iElement) {
        unsigned short VTK_Type = localVolumeElementConnectivity[iElement*SU2_CONN_SIZE+1];
        const auto nPointsElem = nPointsOfElementType(VTK_Type);
        size_t pointA, pointB;
        for (unsigned short iPoint = 0; iPoint < nPointsElem-1; ++iPoint) {
          pointA = localVolumeElementConnectivity[iElement*SU2_CONN_SIZE+2+ iPoint];
          pointB = localVolumeElementConnectivity[iElement*SU2_CONN_SIZE+2+ (iPoint+1)];
          neighborPointsOfPoint[pointA].push_back(pointB);
          neighborPointsOfPoint[pointB].push_back(pointA);
        }
        pointA = localVolumeElementConnectivity[iElement*SU2_CONN_SIZE+2+ 0];
        pointB = localVolumeElementConnectivity[iElement*SU2_CONN_SIZE+2+ nPointsElem-1];
        neighborPointsOfPoint[pointA].push_back(pointB);
        neighborPointsOfPoint[pointB].push_back(pointA);
      }

      /* remove duplicates from the neighboring point lists*/
      vector<size_t>::iterator vecIt;
      for (size_t iPoint = 0; iPoint < numberOfLocalPoints; iPoint++) {
        /* sort neighboring points for each point */
        sort(neighborPointsOfPoint[iPoint].begin(), neighborPointsOfPoint[iPoint].end());

        /* uniquify list of neighboring points */
        vecIt = unique(neighborPointsOfPoint[iPoint].begin(), neighborPointsOfPoint[iPoint].end());

        /* adjust size of vector */
        neighborPointsOfPoint[iPoint].resize(vecIt - neighborPointsOfPoint[iPoint].begin());
      }

      SU2Mesh::GenerateElementBoundingBox();
      SU2Mesh::GenerateADT();
    };
    virtual ~SU2Mesh() = default;
    
    inline unsigned short GetDimension() const { return dimension; }
    
    inline size_t GetNumberOfLocalPoints() const { return numberOfLocalPoints; }
    inline const vector<vector<su2double> >& GetLocalPointCoordinates() const { return localPointCoordinates; }
    
    inline size_t GetNumberOfLocalElements() const { return numberOfLocalElements; }
    inline const vector<size_t>& GetLocalVolumeElementConnectivity() const {return localVolumeElementConnectivity;}

    vector<size_t> GetPointsOfElement(size_t iElement) {
      vector<size_t> elementPoints;
      unsigned short VTK_Type = localVolumeElementConnectivity[iElement * SU2_CONN_SIZE + 1];
      const auto nPointsElem = nPointsOfElementType(VTK_Type);
      for (unsigned short i = 0; i < nPointsElem; i++) {
        elementPoints.push_back(localVolumeElementConnectivity[iElement * SU2_CONN_SIZE + 2 + i]);
      }
      return elementPoints;
    }
    
    inline size_t GetNumberOfMarkers() const { return numberOfMarkers; }
    inline const vector<string>& GetMarkerNames() const { return markerNames; }
    inline size_t GetNumberOfSurfaceElementsForMarker(int val_iMarker) const {return (size_t)surfaceElementConnectivity[val_iMarker].size() / SU2_CONN_SIZE;}
    inline const vector<size_t>& GetSurfaceElementConnectivityForMarker(int val_iMarker) const {return surfaceElementConnectivity[val_iMarker];}

    inline const vector<su2double>& GetLocalVolumeElementBoundingBox() const {return localVolumeElementBoundingBox;}
    
    void PrintMeshDetails() {
      cout << "Mesh details" << endl;
      cout << "NDIME = " << SU2Mesh::GetDimension() << endl;
      cout << "NPOIN = " << SU2Mesh::GetNumberOfLocalPoints() << endl;
      cout << "NELEM = " << SU2Mesh::GetNumberOfLocalElements() << endl;
      cout << "NMARK = "<< SU2Mesh::GetNumberOfMarkers() << endl;
      cout << "BBOX = "<< SU2Mesh::GetLocalVolumeElementBoundingBox().size() / SU2_BBOX_SIZE << endl;
      cout << endl;
    }

    void GenerateElementBoundingBox() {
      array<size_t, N_POINTS_HEXAHEDRON> connectivity{};
      
      for (size_t LocalIndex = 0; LocalIndex < numberOfLocalElements; ++LocalIndex) {  
        /*[xmin, ymin, zmin, xmax, ymax, zmax]*/
        array<su2double, SU2_BBOX_SIZE> bboxCoordinates{su2double_highest, su2double_highest, su2double_lowest, su2double_lowest};
        
        unsigned short VTK_Type = localVolumeElementConnectivity[LocalIndex*SU2_CONN_SIZE + 1];
        const auto nPointsElem = nPointsOfElementType(VTK_Type);
        for (unsigned short i = 0; i < nPointsElem; i++) {
          connectivity[i] = localVolumeElementConnectivity[LocalIndex*SU2_CONN_SIZE + 2 + i];
        }
        
        for (unsigned short iDim = 0; iDim < dimension; iDim++) {
          for (unsigned short i = 0; i < nPointsElem; i++) {
            bboxCoordinates[iDim] = min(bboxCoordinates[iDim], localPointCoordinates[iDim][connectivity[i]]);
            bboxCoordinates[iDim+SU2_BBOX_SIZE/2] = max(bboxCoordinates[iDim+SU2_BBOX_SIZE/2], localPointCoordinates[iDim][connectivity[i]]); /* +2 for 2D */
          }
        }

        // localVolumeElementBoundingBox.push_back(LocalIndex);
        // localVolumeElementBoundingBox.push_back(VTK_Type);
        for (unsigned short i = 0; i < SU2_BBOX_SIZE; i++){
          localVolumeElementBoundingBox.push_back(bboxCoordinates[i]);
        }
      }
    }

    void TraverseADT() {
      adtBoundingBox.printLevelOrder();
    }

    void GenerateADT() {
      array<su2double, SU2_BBOX_SIZE> bboxCoords{};
      
      /* modifying node insertion to balance the tree */
      vector<size_t> elementOrderADT;
      elementOrderADT.resize(numberOfLocalElements);
      iota(elementOrderADT.begin(), elementOrderADT.end(), 0);
      shuffle(elementOrderADT.begin(), elementOrderADT.end(), std::mt19937{std::random_device{}()});

      for (size_t randomIndex = 0; randomIndex < numberOfLocalElements; ++randomIndex) {
        size_t LocalIndex = elementOrderADT[randomIndex];
        for (unsigned short i = 0; i < SU2_BBOX_SIZE; ++i) {
          bboxCoords[i] = localVolumeElementBoundingBox[LocalIndex*SU2_BBOX_SIZE + i];
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
          Coords[iDim] = localPointCoordinates[iDim][pointNumber];
          pointBBoxCoords[iDim] = Coords[iDim] - epsilon; 
          pointBBoxCoords[iDim+2] = Coords[iDim] + epsilon; 
        }
      return pointBBoxCoords;
    }
    
    void MarkInterpolationDonor(const vector<SU2Mesh> &overest_mesh) {
      /* input is vector of overlapping mesh*/
      
      /*by default all points are calculated type i.e = 1*/
      localPointType.resize(numberOfLocalPoints, 1);

      /* mark interpolation donors required by other subgrids*/
      for (const auto& iMesh : overest_mesh) {
        const auto meshNumberOfMarkers = iMesh.GetNumberOfMarkers();
        
        /*iterate over all markers of overset meshes to mark interpolation donor*/
        for (unsigned short iMarker = 0; iMarker < meshNumberOfMarkers; iMarker++){
          // cout << "iMarker = " << iMarker << endl;
          const auto& iMeshSurfaceElementConnectivity = iMesh.GetSurfaceElementConnectivityForMarker(iMarker);
          set<size_t> meshSurfacePoints;
          const auto nElem_Bound = iMesh.GetNumberOfSurfaceElementsForMarker(iMarker);
          
          /*iterate over elements in iMarker and store as points[in a set to avoid duplicates]*/  
          for (size_t iElem_Bound = 0; iElem_Bound < nElem_Bound*SU2_CONN_SIZE; iElem_Bound+=SU2_CONN_SIZE){
            unsigned short VTK_Type = iMeshSurfaceElementConnectivity[iElem_Bound+1]; 
            const auto nPointsElem = nPointsOfElementType(VTK_Type);
            for (unsigned short i = 0; i < nPointsElem; ++i) {
              meshSurfacePoints.insert(iMeshSurfaceElementConnectivity[iElem_Bound+2+i]);
            }
          }
          
          auto numberSurfacePoints = meshSurfacePoints.size();
         
         /*iterate over points in surface of iMarker of iMesh */
          for (const auto iPoint : meshSurfacePoints){
            auto pointBBox = iMesh.GetPointBBox(iPoint);
            const auto BBox = adtBoundingBox.searchADT(pointBBox);
            if (BBox.size() == 0) {
              cerr << "No interpolation stencil found. What to do ??" << endl;
            }
            else {
              if (BBox.size() != 1) {
              cerr << "Multiple interpolation stencil found. What to do ??" << endl;
              }
              for (auto iPoint : GetPointsOfElement(BBox[0])) {
                // cout << "Marking iPoint = " << iPoint << " as donor " << endl;
                localPointType[iPoint] = 2;
                /* mark all neighbors as required (for complete stencil of donor)*/
                for (auto neighborPoint : neighborPointsOfPoint[iPoint]) {
                  if (localPointType[neighborPoint] == 2) {continue;}
                  localPointType[neighborPoint] = 4; /*donor buffer*/
                }
              }
            }
          }
        }

        /*check removability of each point -> if neighbours can be interpolation cells and cell is bigger than overlapping cell*/
        for (size_t iPoint = 0; iPoint < numberOfLocalPoints; iPoint++) {
          /* skip interpolation donors and it's buffer */
          if (localPointType[iPoint] == 2 || localPointType[iPoint] == 4) {
            continue;
          }

          /* TODO: check if donors/overlapping cells are finer/smaller than iPoint. Assumed for now.*/
          bool neighborsCanBeInterpolated = true;
          for (auto neighborPoint : neighborPointsOfPoint[iPoint]) {
            if (localPointType[neighborPoint] == 2 || localPointType[neighborPoint] == 4) {neighborsCanBeInterpolated = false; continue;}
            auto pointBBox = GetPointBBox(neighborPoint);
            const auto BBox = iMesh.adtBoundingBox.searchADT(pointBBox);
            if (BBox.size() == 0) {
              neighborsCanBeInterpolated = false;
            }
          }
          if (neighborsCanBeInterpolated) {
            localPointType[iPoint] = 0; // iPoint becomes unused
            for (auto neighborPoint : neighborPointsOfPoint[iPoint]) {
              if (localPointType[neighborPoint] != 1) {continue;}
              localPointType[neighborPoint] = 3; // neighbors are marked as interpolated
            }
          }

        }
      }
    }

    void WriteTxtPointType(string filename = "pointType.txt") const{
      ofstream pointType(filename);
      if (!pointType.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
      }
      for (size_t iPoint = 0; iPoint < numberOfLocalPoints; iPoint++) {
        pointType << localPointType[iPoint] << " ";
        // if (iPoint % 100 == 0 && iPoint != 0) {pointType << endl;}
      }
      pointType.close();
    }
};
 
int main() {
  SU2Mesh bg_mesh("/Users/zlatangg/Documents/Overset/overset_mesh/SU2/mesh su2/rect-0000-100x100.su2");
  SU2Mesh comp_mesh("/Users/zlatangg/Documents/Overset/overset_mesh/SU2/mesh su2/rect-0505-30x30.su2");
  cout.width(std::numeric_limits<double>::digits10+2);
  cout.precision(std::numeric_limits<double>::digits10+2);
  // bg_mesh.PrintMeshDetails();
  // comp_mesh.PrintMeshDetails();

  // bg_mesh.TraverseADT();
  
  bg_mesh.MarkInterpolationDonor(vector<SU2Mesh> {comp_mesh});

  // auto pointBBox = comp_mesh.GetPointBBox(2);
  // // array<su2double, 4> pointBBox = {0.5-std::numeric_limits<su2double>::epsilon(), 0.5-std::numeric_limits<su2double>::epsilon(), 0.5+std::numeric_limits<su2double>::epsilon(), 0.5+std::numeric_limits<su2double>::epsilon()};
  // auto BBox = bg_mesh.adtBoundingBox.searchADT(pointBBox);
  // cout << "Number of Intersecting BBox found  = " << BBox.size()<<  endl;
  // if (BBox.size() != 0) {
  //   for (auto box : BBox){
  //     cout << box << " ";
  //   }
  //   cout << endl;
  // }
  
  cout << " ----------- " << endl;

  bg_mesh.adtBoundingBox.writeDotToFile("bgADT.txt");
  comp_mesh.adtBoundingBox.writeDotToFile("compADT.txt");
  
  bg_mesh.WriteTxtPointType("bgPtType.txt");
  // comp_mesh.WriteTxtPointType("compPtTypep.txt");
  // comp_mesh.TraverseADT();
  // cout << bg_mesh.adtBoundingBox.root->right->right->right->left->right->elementIndex << endl;
  // auto interpolationStencils = bg_mesh.adtBoundingBox.searchADT(pointBBox);

  // cout << bg_mesh.adtBoundingBox.root->bBoxCoordinates << endl;
  // cout << comp_mesh.adtBoundingBox.root->elementIndex << endl;

  return 0;
}