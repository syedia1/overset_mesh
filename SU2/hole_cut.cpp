#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <array>
#include <cassert>

#include <limits>

using std::vector, std::array;
using std::string;
using std::ifstream, std::istringstream;
using std::cout, std::cin, std::cerr, std::endl;
using std::min, std::max;
using su2double = double;

constexpr su2double su2double_lowest = std::numeric_limits<su2double>::lowest();
constexpr su2double su2double_highest = std::numeric_limits<su2double>::max();

const int SU2_CONN_SIZE   = 10;  /*!< \brief Size of the connectivity array that is allocated for each element*/
const int SU2_BBOX_SIZE   = 8;  /*!< \brief Size of the bounding box array that is allocated for each element*/

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

class SU2Mesh {
  protected:
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
  public:
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
              
              double Coords[3] = {0.0, 0.0, 0.0};
              for (auto iPoint = 0ul; iPoint < numberOfLocalPoints; iPoint++) {
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

                  surfaceElementConnectivity[iMarker].push_back(0);
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
        // cout << "NDIM = " << dimension << endl;
        // cout << "NPOIN = " << localPointCoordinates[0].size() << endl;
        // cout << "NELEM = " << localVolumeElementConnectivity.size() << endl;
        // cout << "NMARK = " << markerNames.size() << endl;
      }
      else { 
        cerr << "Unable to open file!" << endl; 
      } 

      SU2Mesh::GenerateElementBoundingBox();
    };
    virtual ~SU2Mesh() = default;
    
    inline unsigned short GetDimension() const { return dimension; }
    
    inline size_t GetNumberOfLocalPoints() const { return numberOfLocalPoints; }
    inline const vector<vector<su2double> >& GetLocalPointCoordinates() const { return localPointCoordinates; }
    
    inline size_t GetNumberOfLocalElements() const { return numberOfLocalElements; }
    inline const vector<size_t>& GetLocalVolumeElementConnectivity() const {return localVolumeElementConnectivity;}
    
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
      
      // array<su2double, N_POINTS_PRISM> bboxCoordinates{-1.0, -1.0, -1.0, -1.0, -1.0, -1.0};
      for (size_t LocalIndex = 0; LocalIndex < numberOfLocalElements; ++LocalIndex) {  
        /*[xmin, ymin, zmin, xmax, ymax, zmax]*/
        array<su2double, N_POINTS_PRISM> bboxCoordinates{su2double_highest, su2double_highest, su2double_highest, su2double_lowest, su2double_lowest, su2double_lowest};
        
        unsigned short VTK_Type = localVolumeElementConnectivity[LocalIndex*SU2_CONN_SIZE + 1];
        const auto nPointsElem = nPointsOfElementType(VTK_Type);
        for (unsigned short i = 0; i < nPointsElem; i++) {
          connectivity[i] = localVolumeElementConnectivity[LocalIndex*SU2_CONN_SIZE + 2 + i];
        }
        
        for (unsigned short iDim = 0; iDim < dimension; iDim++) {
          for (unsigned short i = 0; i < nPointsElem; i++) {
            bboxCoordinates[iDim] = min(bboxCoordinates[iDim], localPointCoordinates[iDim][connectivity[i]]);
            bboxCoordinates[iDim+3] = max(bboxCoordinates[iDim+3], localPointCoordinates[iDim][connectivity[i]]);
          }
        }

        localVolumeElementBoundingBox.push_back(LocalIndex);
        localVolumeElementBoundingBox.push_back(VTK_Type);
        for (unsigned short i = 0; i < N_POINTS_PRISM; i++){
          localVolumeElementBoundingBox.push_back(bboxCoordinates[i]);
        }
      }
    }
};
 
int main() {
  SU2Mesh bg_mesh("/Users/zlatangg/Documents/Overset/heat_cond SU2/square_10x10.su2");
  SU2Mesh comp_mesh("/Users/zlatangg/Documents/Overset/heat_cond SU2/square_3x3.su2");
  
  bg_mesh.PrintMeshDetails();
  auto bg_eleConn = bg_mesh.GetLocalVolumeElementConnectivity();
  auto bg_bBox = bg_mesh.GetLocalVolumeElementBoundingBox();

  auto check_num_ele = 4;
  cout << "Element Connectivity -> " << endl;
  for(unsigned short i = 0; i < check_num_ele; ++i){
    cout << "Element " << i << " : ";
    for(unsigned short j = 0; j < SU2_CONN_SIZE; ++j) {
      cout << bg_eleConn[i*SU2_CONN_SIZE + j] << '\t';
    }
    cout << endl;
  }
  cout << endl << "Bounding Box -> " << endl;
  for(unsigned short i = 0; i < check_num_ele; ++i){
    cout << "Element " << i << " : ";
    for(unsigned short j = 0; j < SU2_BBOX_SIZE; ++j) {
      cout << bg_bBox[i*SU2_BBOX_SIZE + j] << '\t';
    }
    cout << endl;
  }
  // comp_mesh.PrintMeshDetails();

  return 0;
}