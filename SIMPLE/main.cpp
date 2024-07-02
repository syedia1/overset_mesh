#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <queue>
#include <random>
#include <stack>
#include <string>
#include <vector>

using namespace std;

const size_t ITER_PRINT_FREQ = 5000;
const double PI = 3.141592653589793238463;
const size_t NS_MAXITER = 20;
const double TOLERANCE = 1024 * std::numeric_limits<double>::epsilon();

using su2double = double;
const size_t SU2_BBOX_SIZE = 4; /*!< \brief Size of the bounding box array that is allocated for each element*/
constexpr su2double su2double_lowest = std::numeric_limits<su2double>::lowest();
constexpr su2double su2double_highest = std::numeric_limits<su2double>::max();

enum POINT_TYPE {
    UNUNSED = 0,
    CALCULATED = 1,
    INTERPOLATION_DONOR = 2,
    INTERPOLATION_RECIEVER = 3,
    DONOR_BUFFER = 4,
    BC_SPECIFIED = 5,
};

vector<double> solveAxB(const vector<vector<double> > &A, const vector<double> &B) {
    size_t n = A.size();
    vector<double> x(n, 0.0);

    // Pre-calculate inverses of diagonal elements (avoid division within the loop)
    vector<double> inv_diag(n, 0.0);
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
            converged &= (abs(x[j] * A[j][j] + sum - B[j]) <= TOLERANCE);
        }

        if (converged) {
            break;
        }
    }
    // cout << "GS=" << GS_itr << ", ";
    return x;
}

vector<double> matAdd(const vector<double> &A, const vector<double> &B) {
    size_t m = A.size();
    vector<double> C(m, 0.0);
    for (size_t i = 0; i < m; i++) {
        C[i] = A[i] + B[i];
    }
    return C;
}

vector<double> interpolation_coeff(const vector<vector<double> > points, const vector<double> xy) {
    auto transform = [](const vector<double> &coeff, const vector<double> &points) -> double {
        const double phi = coeff[0], chi = coeff[1];
        return (1.0 - phi) * (1.0 - chi) * points[0] + (phi) * (1.0 - chi) * points[1] + (phi) * (chi)*points[2] + (1.0 - phi) * (chi)*points[3];
    };
    auto partialDerivative_phi = [](const vector<double> &coeff, const vector<double> &points) {
        const double chi = coeff[1];
        return (-1.0) * (1.0 - chi) * points[0] + (1.0 - chi) * points[1] + (chi)*points[2] + (-1.0) * (chi)*points[3];
    };
    auto partialDerivative_chi = [](const vector<double> &coeff, const vector<double> &points) {
        const double phi = coeff[0];
        return (1.0 - phi) * (-1.0) * points[0] + (phi) * (-1.0) * points[1] + (phi)*points[2] + (1.0 - phi) * points[3];
    };

    vector<double> coeff = {0.5, 0.5};  // initial guess value of phi,

    // NR, Newton raphson to minimze B, Ax = B
    size_t NR_itr = 0;
    for (; NR_itr < NS_MAXITER; NR_itr++) {
        vector<vector<double> > A = {{partialDerivative_phi(coeff, points[0]), partialDerivative_chi(coeff, points[0])},
                                     {partialDerivative_phi(coeff, points[1]), partialDerivative_chi(coeff, points[1])}};
        vector<double> B = {-1.0 * transform(coeff, points[0]) + xy[0], -1.0 * transform(coeff, points[1]) + xy[1]};

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
    return vector<double>{(1.0 - coeff[0]) * (1.0 - coeff[1]), (coeff[0]) * (1.0 - coeff[1]), (coeff[0]) * (coeff[1]), (1.0 - coeff[0]) * (coeff[1])};
}

class node_adt {
   public:
    size_t elementIndex;
    array<su2double, SU2_BBOX_SIZE> bBoxCoordinates;
    node_adt *left;
    node_adt *right;

    node_adt(const size_t elementIdx, const array<su2double, SU2_BBOX_SIZE> &bBoxCoords) : elementIndex(elementIdx), bBoxCoordinates(bBoxCoords), left(nullptr), right(nullptr) {}
    virtual ~node_adt() = default;
};

class ADT {
   public:
    node_adt *root;
    size_t treeHeirarchy;
    ADT() : root(nullptr), treeHeirarchy(0) {}
    virtual ~ADT() = default;

    void insertNode(node_adt *node) {
        node_adt *current = root;
        node_adt *parent = nullptr;
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
        node_adt *current = nullptr;
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
                for (unsigned short iDim = 0; iDim < SU2_BBOX_SIZE / 2; iDim++) {
                    // intersect = intersect && (testBBox[iDim] < current->bBoxCoordinates[iDim+SU2_BBOX_SIZE/2] || floatEqual(testBBox[iDim], current->bBoxCoordinates[iDim+SU2_BBOX_SIZE/2]));
                    // intersect = intersect && (testBBox[iDim+SU2_BBOX_SIZE/2] > current->bBoxCoordinates[iDim] || floatEqual(testBBox[iDim+SU2_BBOX_SIZE/2], current->bBoxCoordinates[iDim]));
                    intersect = intersect && (testBBox[iDim] <= current->bBoxCoordinates[iDim + SU2_BBOX_SIZE / 2]);
                    intersect = intersect && (testBBox[iDim + SU2_BBOX_SIZE / 2] >= current->bBoxCoordinates[iDim]);
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
                if (i < SU2_BBOX_SIZE / 2) {
                    /*curr->left is always searched as the left has min coords lower than current which gives no info on intersection */
                    searchQ.push(make_pair(current->left, currHeirarchy + 1));
                    // if (floatEqual(testBBox[i+SU2_BBOX_SIZE/2], current->bBoxCoordinates[i])) {
                    //   current = current->right;
                    // }
                    /*Test _max < Current _min*/
                    if (testBBox[i + SU2_BBOX_SIZE / 2] < current->bBoxCoordinates[i]) {
                        current = nullptr;
                    } else {
                        current = current->right;
                    }
                    // searchQ.push(make_pair(current->right, currHeirarchy));
                }
                /*branching based on maximum coordinate*/
                else {
                    /*curr->right is always searched as the right has max coords higher than current which gives no info on intersection */
                    searchQ.push(make_pair(current->right, currHeirarchy + 1));
                    // if (floatEqual(testBBox[i-SU2_BBOX_SIZE/2], current->bBoxCoordinates[i])) {
                    //   current = current->left;
                    // }
                    /*Test _min > Current _max*/
                    if (testBBox[i - SU2_BBOX_SIZE / 2] > current->bBoxCoordinates[i]) {
                        current = nullptr;
                    } else {
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
        queue<node_adt *> q;
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
            while (count > 0) {
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
    std::string generateDot(node_adt *node) const {
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
        } else if (node->right != nullptr) {
            dot += "\"" + nodeName + "_NULL\" ";
            dot += "[shape=point]\n";
            dot += "\"" + nodeName + "\"" + " -> \"" + nodeName + "_NULL\"";
            dot += "[color=red]\n";
        }
        if (node->right != nullptr) {
            dot += "\"" + nodeName + "\"" + " -> \"" + std::to_string(node->right->elementIndex) + "\"";
            dot += "[color=blue]\n";
            dot += generateDot(node->right);
        } else if (node->left != nullptr) {
            dot += "\"" + nodeName + "_NULL\" ";
            dot += "[shape=point]\n";
            dot += "\"" + nodeName + "\"" + " -> \"" + nodeName + "_NULL\"";
            dot += "[color=blue]\n";
        }
        return dot;
    }

    // Generate DOT syntax and write to a file
    void writeDotToFile(const std::string &filename) const {
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

class DataStructure {
   public:
    int Dim;
    int nVar;

    double xO, yO, theta;
    double dx, dy, Lx, Ly;
    int Nx, Ny, numberOfPoints, numberOfElements;
    double dt;
    double volp;

    double ulid;
    double nu;
    double rho;

    vector<vector<vector<double> > > Var;
    vector<vector<vector<double> > > VarOld;
    vector<vector<vector<double> > > Ff;
    vector<double> Init;
    vector<double> residual;

    vector<unsigned short> pointType;
    vector<vector<double> > pointCoordinates;
    vector<vector<size_t> > neighborPointsOfPoint;
    vector<vector<size_t> > elementConnectivity;
    vector<vector<double> > elementBBox;

    vector<vector<size_t> > interpolationStencil;
    vector<vector<double> > interpolationCoeffs;

    /* bottom, right, top, left*/
    vector<vector<size_t> > boundaryPoints;
    /* bottom, right, top, left*/
    vector<double> boundaryCondition;

    ADT adtBoundingBox;

    DataStructure(double xOrigin, double yOrigin, double axisTheta, double lengthX, double lengthY, size_t _Nx, size_t _Ny) {
        Dim = 2;
        nVar = 3;

        xO = xOrigin, yO = yOrigin, theta = axisTheta;
        Lx = lengthX, Ly = lengthY;
        Nx = _Nx, Ny = _Ny;
        dx = Lx / Nx;
        dy = Ly / Ny;

        numberOfPoints = (Nx + 2) * (Ny + 2);
        numberOfElements = (Nx + 1) * (Ny + 1);
        dt = 0.001;
        volp = dx * dy;

        /* Physical constants*/
        ulid = 1.0;
        nu = 0.01;
        rho = 1.0;

        Init.resize(nVar, 0.0);
        Var.resize(nVar);
        VarOld.resize(nVar);

        interpolationStencil.resize(numberOfPoints);
        interpolationCoeffs.resize(numberOfPoints);

        for (int i = 0; i < nVar; i++) {
            Var[i].resize(Nx + 2);
            VarOld[i].resize(Nx + 2);
            for (int j = 0; j < (Nx + 2); j++) {
                Var[i][j].resize(Ny + 2, 0.0);
                VarOld[i][j].resize(Ny + 2, 0.0);
            }
        }
        residual.resize(nVar);
        Ff.resize(2 * Dim);  // No of faces per cell = 4  i=0,1,2,3 -> E,N,W,S
        for (int k = 0; k < (2 * Dim); k++) {
            Ff[k].resize(Nx + 2);
            for (int j = 0; j < (Nx + 2); j++) {
                Ff[k][j].resize(Ny + 2, 0.0);
            }
        }

        /* requirements for overset mesh */
        /* Point coordinates */
        pointCoordinates.resize(Dim);
        for (unsigned short iDim = 0; iDim < Dim; ++iDim) {
            pointCoordinates[iDim].resize(numberOfPoints);
        }
        pointType.resize(numberOfPoints, CALCULATED);
        size_t iPoint = numberOfPoints;
        for (size_t j = 0; j < (Ny + 2); ++j) {
            for (size_t i = 0; i < (Nx + 2); i++) {
                iPoint = GetPointNumber(i, j);
                /* cell centered so adding dx/2 and dy/2 to point coordinates*/
                pointCoordinates[0][iPoint] = dx / 2 + xO + (i * dx) * cos(theta) - (j * dy) * sin(theta); /* global x axis*/
                pointCoordinates[1][iPoint] = dy / 2 + yO + (i * dx) * sin(theta) + (j * dy) * cos(theta); /* global y axis*/
            }
        }

        /* Element connectivity and it's bounding box */
        elementConnectivity.resize(numberOfElements);
        elementBBox.resize(numberOfElements);
        size_t iElement = numberOfElements;
        for (size_t j = 0; j < Ny + 1; ++j) {
            for (size_t i = 0; i < Nx + 1; ++i) {
                iElement = GetElementNumber(i, j);
                elementConnectivity[iElement].assign({j * (Nx + 2) + i, j * (Nx + 2) + i + 1, (j + 1) * (Nx + 2) + i + 1, (j + 1) * (Nx + 2) + i});

                /* Element cartesian aligned Bounding Box */
                elementBBox[iElement].assign({su2double_highest, su2double_highest, su2double_lowest, su2double_lowest}); /* {Xmin, Ymin, Xmax, Ymax}*/
                for (size_t iPoint = 0; iPoint < elementConnectivity[iElement].size(); ++iPoint) {
                    double x = pointCoordinates[0][elementConnectivity[iElement][iPoint]];
                    double y = pointCoordinates[1][elementConnectivity[iElement][iPoint]];
                    elementBBox[iElement][0] = min(elementBBox[iElement][0], x);
                    elementBBox[iElement][1] = min(elementBBox[iElement][1], y);
                    elementBBox[iElement][2] = max(elementBBox[iElement][2], x);
                    elementBBox[iElement][3] = max(elementBBox[iElement][3], y);
                }
            }
        }

        /* Neighbour points of point */
        neighborPointsOfPoint.resize(numberOfPoints);
        for (size_t j = 0; j < Ny + 1; ++j) {
            for (size_t i = 0; i < Nx + 1; ++i) {
                iElement = GetElementNumber(i, j);
                size_t pointA, pointB;
                for (size_t iPoint = 0; iPoint < elementConnectivity[iElement].size(); iPoint++) {
                    pointA = elementConnectivity[iElement][iPoint];
                    pointB = elementConnectivity[iElement][(iPoint + 1) % 4];
                    neighborPointsOfPoint[pointA].push_back(pointB);
                    neighborPointsOfPoint[pointB].push_back(pointA);
                }
            }
        }
        /* remove duplicates from the neighboring point lists*/
        iPoint = numberOfPoints;
        vector<size_t>::iterator vecIt;
        for (size_t j = 0; j < Ny + 2; ++j) {
            for (size_t i = 0; i < Nx + 2; ++i) {
                iPoint = GetPointNumber(i, j);
                /* sort neighboring points for each point */
                sort(neighborPointsOfPoint[iPoint].begin(), neighborPointsOfPoint[iPoint].end());

                /* uniquify list of neighboring points */
                vecIt = unique(neighborPointsOfPoint[iPoint].begin(), neighborPointsOfPoint[iPoint].end());

                /* adjust size of vector */
                neighborPointsOfPoint[iPoint].resize(vecIt - neighborPointsOfPoint[iPoint].begin());
            }
        }

        /* Boundary Surfaces*/
        boundaryPoints.resize(4); /* bottom, right, top, left*/
        for (size_t iBoundPoint = 0; iBoundPoint < (Nx + 2); iBoundPoint++) {
            boundaryPoints[0].push_back(iBoundPoint);                       /*bottom*/
            boundaryPoints[2].push_back((Ny + 1) * (Nx + 2) + iBoundPoint); /*top*/
        }
        for (size_t iBoundPoint = 0; iBoundPoint < (Ny + 2); iBoundPoint++) {
            boundaryPoints[3].push_back(iBoundPoint * (Nx + 2));            /*left*/
            boundaryPoints[1].push_back((Nx + 2) * iBoundPoint + (Nx + 1)); /*right*/
        }
        reverse(boundaryPoints[2].begin(), boundaryPoints[2].end());
        reverse(boundaryPoints[3].begin(), boundaryPoints[3].end());

        GenerateADT();
    }

    inline size_t GetElementNumber(size_t i, size_t j) const { return j * (Nx + 1) + i; }
    inline size_t GetPointNumber(size_t i, size_t j) const { return j * (Nx + 2) + i; }
    inline std::tuple<size_t, size_t> GetijFromPointNumber(size_t iPoint) const { return {iPoint % (Nx + 2), iPoint / (Nx + 2)}; }

    virtual ~DataStructure() = default;

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
            node_adt *tempNode = new node_adt(LocalIndex, bboxCoords);
            adtBoundingBox.insertNode(tempNode);
        }
        cout << "Max heirarchy : " << adtBoundingBox.treeHeirarchy << endl;
    }

    array<su2double, SU2_BBOX_SIZE> GetPointBBox(const size_t pointNumber) const {
        /* small value to generate a bounding box centered around a point */
        su2double epsilon = std::numeric_limits<su2double>::epsilon();
        array<su2double, 3> Coords{0.0, 0.0, 0.0};
        array<su2double, SU2_BBOX_SIZE> pointBBoxCoords{};
        for (unsigned short iDim = 0; iDim < Dim; ++iDim) {
            Coords[iDim] = pointCoordinates[iDim][pointNumber];
            pointBBoxCoords[iDim] = Coords[iDim] - epsilon;
            pointBBoxCoords[iDim + 2] = Coords[iDim] + epsilon;
        }
        return pointBBoxCoords;
    }

    void MarkOversetBoundaryDonors(const DataStructure &oversetMesh) {
        /* does a type of holecutting / removal of coarse cells where overlapping finer cells are present*/

        /* iterate over all the boundaries of overset mesh to mark donors */
        for (size_t iBound = 0; iBound < oversetMesh.boundaryPoints.size(); ++iBound) {
            /* iterate over all points in iBound*/
            for (size_t iBoundPoint = 0; iBoundPoint < oversetMesh.boundaryPoints[iBound].size(); ++iBoundPoint) {
                size_t point = oversetMesh.boundaryPoints[iBound][iBoundPoint];
                auto pointBBox = oversetMesh.GetPointBBox(point);
                const auto BBox = adtBoundingBox.searchADT(pointBBox);
                if (BBox.size() == 0) {
                    cerr << "Point: " << point << " No interpolation stencil found. What to do ??" << endl;
                } else {
                    if (BBox.size() != 1) {
                        // cerr << "Point : "<< point << " Multiple interpolation stencil found [" << BBox.size() << "]. What to do ??" << endl;
                    }
                    for (auto iPoint : elementConnectivity[BBox[0]]) {
                        // cout << "Marking iPoint = " << iPoint << " as donor " << endl;
                        pointType[iPoint] = 2;
                        /* mark all neighbors as required (for complete stencil of donor)*/
                        for (auto neighborPoint : neighborPointsOfPoint[iPoint]) {
                            if (pointType[neighborPoint] == 2) {
                                continue;
                            }
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
                if (pointType[neighborPoint] == INTERPOLATION_DONOR || pointType[neighborPoint] == DONOR_BUFFER) {
                    neighborsCanBeInterpolated = false;
                    break;
                }
                auto pointBBox = GetPointBBox(neighborPoint);
                const auto BBox = oversetMesh.adtBoundingBox.searchADT(pointBBox);
                if (BBox.size() == 0) {
                    neighborsCanBeInterpolated = false;
                }
            }
            if (neighborsCanBeInterpolated) {
                pointType[iPoint] = UNUNSED;  // iPoint becomes unused
                for (auto neighborPoint : neighborPointsOfPoint[iPoint]) {
                    if (pointType[neighborPoint] != CALCULATED) {
                        continue;
                    }
                    pointType[neighborPoint] = INTERPOLATION_RECIEVER;  // neighbors are marked as interpolated
                }
            }
        }
    }

    /* pType : 3 = interpolated, 5 = boundary condition specified */
    void MarkBoundaryPointType(const unsigned short ptType) {
        /* marks all the boundaries (for component meshes marked as intepolated. Assumed component meshes dont overlap each other nor do they intersect boundaries) */
        for (size_t iBound = 0; iBound < boundaryPoints.size(); ++iBound) {
            /* iterate over all points in iBound*/
            for (size_t iBoundPoint = 0; iBoundPoint < boundaryPoints[iBound].size(); ++iBoundPoint) {
                size_t point = boundaryPoints[iBound][iBoundPoint];
                pointType[point] = ptType;
            }
        }
    }

    void StoreInterpolationStencil(const DataStructure &oversetMesh) {
        for (size_t iPoint = 0; iPoint < numberOfPoints; ++iPoint) {
            /* skip points which are not interpolated*/
            if (pointType[iPoint] != INTERPOLATION_RECIEVER) {
                continue;
            }
            auto pointBBox = GetPointBBox(iPoint);
            const auto BBox = oversetMesh.adtBoundingBox.searchADT(pointBBox);
            if (BBox.size() == 0) {
                cerr << " No donor found for interpolation marked point." << endl;
            } else {
                // interpolationStencil[iPoint].push_back(BBox[0]);
                auto donorStencil = oversetMesh.elementConnectivity[BBox[0]];
                vector<vector<double> > donorCoords(Dim);
                for (size_t iDonor = 0; iDonor < donorStencil.size(); ++iDonor) {
                    interpolationStencil[iPoint].push_back(donorStencil[iDonor]);
                    for (unsigned short iDim = 0; iDim < Dim; ++iDim) {
                        donorCoords[iDim].push_back(oversetMesh.pointCoordinates[iDim][donorStencil[iDonor]]);
                    }
                }
                vector<double> receiverCoords;
                for (unsigned short iDim = 0; iDim < Dim; ++iDim) {
                    receiverCoords.push_back(pointCoordinates[iDim][iPoint]);
                }
                interpolationCoeffs[iPoint] = interpolation_coeff(donorCoords, receiverCoords);
            }
        }
    }

    void WriteTxtPointType(string filename = "pointType.txt") const {
        ofstream pointTypeFile(filename);
        if (!pointTypeFile.is_open()) {
            std::cerr << "Error opening file: " << filename << std::endl;
        } else {
            for (size_t iPoint = 0; iPoint < numberOfPoints; iPoint++) {
                // if (iPoint % (Nx) == 0 && iPoint !=0) {pointTypeFile << endl;}
                pointTypeFile << pointType[iPoint] << " ";
            }
            pointTypeFile.close();
        }
    }

    void AdjustNxNy() {
        /* Mesh related data structures are using 0 to Nx indexing (including halo cells)
        but Solver is using 0 to Nx+1 based indexing to acocunt to 2 halo cells in each direction*/
        Nx = Nx - 2;
        Ny = Ny - 2;
    }
};

void CopyNewtoOld(DataStructure *rect) {
    std::copy(rect->Var.begin(), rect->Var.end(), rect->VarOld.begin());
}

void ApplyBC(DataStructure *rect, int k, DataStructure *oversetMesh) {
    if (rect->pointType[0] == INTERPOLATION_RECIEVER) {
        // cout << "Interpolating at boundary for k = " << k << endl;
        for (size_t iBound = 0; iBound < rect->boundaryPoints.size(); ++iBound) {
            for (size_t iBoundPoint = 0; iBoundPoint < rect->boundaryPoints[iBound].size(); ++iBoundPoint) {
                size_t point = rect->boundaryPoints[iBound][iBoundPoint];
                auto [i, j] = rect->GetijFromPointNumber(point);
                rect->Var[k][i][j] = 0.0;
                for (unsigned short iDonor = 0; iDonor < rect->interpolationStencil[point].size(); ++iDonor) {
                    auto donorPointNumber = rect->interpolationStencil[point][iDonor];
                    auto [_iDonor, _jDonor] = oversetMesh->GetijFromPointNumber(donorPointNumber);
                    rect->Var[k][i][j] = rect->Var[k][i][j] + oversetMesh->Var[k][_iDonor][_jDonor] * rect->interpolationCoeffs[point][iDonor];
                }
            }
        }
    } else {
        // cout << "Applying specified BC for k = " << k << endl;
        switch (k) {
            case 0:  // U

                for (int j = 1; j < rect->Ny + 1; j++) {
                    rect->Var[k][0][j] = 2 * (0.0) - rect->Var[k][1][j];                    // Left
                    rect->Var[k][rect->Nx + 1][j] = 2 * (0.0) - rect->Var[k][rect->Nx][j];  // Right
                }
                for (int i = 1; i < rect->Nx + 1; i++) {
                    rect->Var[k][i][rect->Ny + 1] = 2 * rect->ulid - rect->Var[k][i][rect->Ny];  // Top
                    rect->Var[k][i][0] = 2 * (0.0) - rect->Var[k][i][1];                         // Bottom
                }
                break;

            case 1:  // V

                for (int j = 1; j < rect->Ny + 1; j++) {
                    rect->Var[k][0][j] = 2 * (0.0) - rect->Var[k][1][j];                    // Left
                    rect->Var[k][rect->Nx + 1][j] = 2 * (0.0) - rect->Var[k][rect->Nx][j];  // Right
                }
                for (int i = 1; i < rect->Nx + 1; i++) {
                    rect->Var[k][i][rect->Ny + 1] = 2 * (0.0) - rect->Var[k][i][rect->Ny];  // Top
                    rect->Var[k][i][0] = 2 * (0.0) - rect->Var[k][i][1];                    // Bottom
                }
                break;

            case 2:  // P (All Neumann Condition)

                for (int j = 1; j < rect->Ny + 1; j++) {
                    rect->Var[k][0][j] = rect->Var[k][1][j];                    // Left
                    rect->Var[k][rect->Nx + 1][j] = rect->Var[k][rect->Nx][j];  // Right
                }
                for (int i = 1; i < rect->Nx + 1; i++) {
                    rect->Var[k][i][rect->Ny + 1] = rect->Var[k][i][rect->Ny];  // Top
                    rect->Var[k][i][0] = rect->Var[k][i][1];                    // Bottom
                }
                break;
        }
    }
}

void LinearInterpolation(DataStructure *rect) {
    for (int i = 1; i < rect->Nx + 1; i++) {
        for (int j = 1; j < rect->Ny + 1; j++) {
            size_t point = rect->GetPointNumber(i, j);
            if (rect->pointType[point] == UNUNSED || rect->pointType[point] == INTERPOLATION_RECIEVER) {
                continue;
            }
            rect->Ff[0][i][j] = (rect->Var[0][i][j] + rect->Var[0][i + 1][j]) * rect->dy * 0.5;   // East Face
            rect->Ff[1][i][j] = (rect->Var[1][i][j] + rect->Var[1][i][j + 1]) * rect->dx * 0.5;   // North Face
            rect->Ff[2][i][j] = -(rect->Var[0][i][j] + rect->Var[0][i - 1][j]) * rect->dy * 0.5;  // West Face
            rect->Ff[3][i][j] = -(rect->Var[1][i][j] + rect->Var[1][i][j - 1]) * rect->dx * 0.5;  // South Face
        }
    }
}

void initialize(DataStructure *rect, DataStructure *oversetMesh) {
    // initializing interior values
    for (int k = 0; k < rect->nVar; k++) {
        // for (int i = 1; i < rect->Nx + 1; i++) {
        //     for (int j = 1; j < rect->Ny + 1; j++) {
        //         rect->Var[k][i][j] = rect->Init[k];
        //     }
        // }
        ApplyBC(rect, k, oversetMesh);
    }

    CopyNewtoOld(rect);
    LinearInterpolation(rect);
}

void GetOutput(DataStructure *rect, std::string input) {
    FILE *fp;
    fp = fopen(input.c_str(), "w");

    for (int k = 0; k < rect->nVar; k++) {
        fprintf(fp, "\n ########## Data for k = %d ############ \n", k);
        for (int j = 0; j < rect->Ny + 2; j++) {
            for (int i = 0; i < rect->Nx + 2; i++) {
                fprintf(fp, "%lf \t", rect->Var[k][i][j]);
            }
            fprintf(fp, "\n");
        }
    }
    fclose(fp);
}

void SimpleUpwind(DataStructure *rect, double *Fc, double *ap_c, int i, int j, int k) {
    double ue, uw, un, us;
    double sum = 0;
    if (rect->Ff[0][i][j] >= 0) {
        ue = rect->Var[k][i][j];
        sum += rect->Ff[0][i][j];
    } else
        ue = rect->Var[k][i + 1][j];

    if (rect->Ff[2][i][j] >= 0) {
        uw = rect->Var[k][i][j];
        sum += rect->Ff[2][i][j];
    } else
        uw = rect->Var[k][i - 1][j];

    if (rect->Ff[1][i][j] >= 0) {
        un = rect->Var[k][i][j];
        sum += rect->Ff[1][i][j];
    } else
        un = rect->Var[k][i][j + 1];

    if (rect->Ff[3][i][j] >= 0) {
        us = rect->Var[k][i][j];
        sum += rect->Ff[3][i][j];
    } else
        us = rect->Var[k][i][j - 1];

    *Fc = (ue * rect->Ff[0][i][j] + uw * rect->Ff[2][i][j] + un * rect->Ff[1][i][j] + us * rect->Ff[3][i][j]);
    *ap_c = sum * rect->volp;
}

void Quick(DataStructure *rect, double *Fc, double *ap_c, int i, int j, int k) {
    double ue, uw, un, us;
    double sum = 0;
    // East
    if (rect->Ff[0][i][j] >= 0) {
        ue = 0.75 * rect->Var[k][i][j] + 0.375 * rect->Var[k][i + 1][j] - 0.125 * rect->Var[k][i - 1][j];
        sum += 0.75 * rect->Ff[0][i][j];
    } else {
        ue = 0.75 * rect->Var[k][i + 1][j] + 0.375 * rect->Var[k][i][j] - 0.125 * rect->Var[k][i + 2][j];
        sum += 0.375 * rect->Ff[0][i][j];
    }
    // West
    if (rect->Ff[2][i][j] >= 0) {
        uw = 0.75 * rect->Var[k][i][j] + 0.375 * rect->Var[k][i - 1][j] - 0.125 * rect->Var[k][i + 1][j];
        sum += 0.75 * rect->Ff[2][i][j];
    } else {
        uw = 0.75 * rect->Var[k][i - 1][j] + 0.375 * rect->Var[k][i][j] - 0.125 * rect->Var[k][i - 2][j];
        sum += 0.375 * rect->Ff[2][i][j];
    }
    // North
    if (rect->Ff[1][i][j] >= 0) {
        un = 0.75 * rect->Var[k][i][j] + 0.375 * rect->Var[k][i][j + 1] - 0.125 * rect->Var[k][i][j - 1];
        sum += 0.75 * rect->Ff[1][i][j];
    } else {
        un = 0.75 * rect->Var[k][i][j + 1] + 0.375 * rect->Var[k][i][j] - 0.125 * rect->Var[k][i][j + 2];
        sum += 0.375 * rect->Ff[1][i][j];
    }
    if (rect->Ff[3][i][j] >= 0) {
        us = 0.75 * rect->Var[k][i][j] + 0.375 * rect->Var[k][i][j - 1] - 0.125 * rect->Var[k][i][j + 1];
        sum += 0.75 * rect->Ff[3][i][j];
    } else {
        us = 0.75 * rect->Var[k][i][j - 1] + 0.375 * rect->Var[k][i][j] - 0.125 * rect->Var[k][i][j - 2];
        sum += 0.375 * rect->Ff[3][i][j];
    }
    *Fc = (ue * rect->Ff[0][i][j] + uw * rect->Ff[2][i][j] + un * rect->Ff[1][i][j] + us * rect->Ff[3][i][j]);
    *ap_c = sum * rect->volp;
}

void DiffusiveFlux(DataStructure *rect, double *Fd, double *ap_d, int i, int j, int k) {
    *Fd = rect->volp * ((rect->Var[k][i + 1][j] - 2.0 * rect->Var[k][i][j] + rect->Var[k][i - 1][j]) / (rect->dx * rect->dx) + (rect->Var[k][i][j + 1] - 2.0 * rect->Var[k][i][j] + rect->Var[k][i][j - 1]) / (rect->dy * rect->dy));
    *ap_d = -rect->volp * (2.0 / (rect->dx * rect->dx) + 2.0 / (rect->dy * rect->dy));
}

void UpdateFlux(DataStructure *rect) {
    for (int i = 1; i < rect->Nx + 1; i++) {
        for (int j = 1; j < rect->Ny + 1; j++) {
            size_t point = rect->GetPointNumber(i, j);
            if (rect->pointType[point] == UNUNSED || rect->pointType[point] == INTERPOLATION_RECIEVER) {
                continue;
            }
            rect->Ff[0][i][j] += -rect->dt / rect->rho * (rect->Var[2][i + 1][j] - rect->Var[2][i][j]) * rect->dy / rect->dx;  // East Face
            rect->Ff[1][i][j] += -rect->dt / rect->rho * (rect->Var[2][i][j + 1] - rect->Var[2][i][j]) * rect->dx / rect->dy;  // North Face
            rect->Ff[2][i][j] += -rect->dt / rect->rho * (rect->Var[2][i - 1][j] - rect->Var[2][i][j]) * rect->dy / rect->dx;  // West Face
            rect->Ff[3][i][j] += -rect->dt / rect->rho * (rect->Var[2][i][j - 1] - rect->Var[2][i][j]) * rect->dx / rect->dy;  // South Face
        }
    }
}

double SolveUV(DataStructure *rect, DataStructure *oversetMesh, int k) {
    /* k determines U or V. k = 0 -> U*/
    int count = 0;
    double Fc, ap_c, Fd, ap_d, ap, R, rms;
    rms = 0.0;
    for (int i = 1; i < rect->Nx + 1; i++) {
        for (int j = 1; j < rect->Ny + 1; j++) {
            size_t point = rect->GetPointNumber(i, j);
            if (rect->pointType[point] == UNUNSED) {
                continue;
            }
            if (rect->pointType[point] == CALCULATED || rect->pointType[point] == INTERPOLATION_DONOR || rect->pointType[point] == DONOR_BUFFER) {
                SimpleUpwind(rect, &Fc, &ap_c, i, j, k);
                // Quick(rect, &Fc, &ap_c, i, j, k);
                DiffusiveFlux(rect, &Fd, &ap_d, i, j, k);
                R = -(rect->volp / rect->dt * (rect->Var[k][i][j] - rect->VarOld[k][i][j]) + Fc + (-rect->nu) * Fd);
                ap = rect->volp / rect->dt + ap_c + (-rect->nu) * ap_d;
                count++;
                rect->Var[k][i][j] = rect->Var[k][i][j] + R / ap;
            } else if (rect->pointType[point] == INTERPOLATION_RECIEVER) {
                double temp = rect->Var[k][i][j];
                rect->Var[k][i][j] = 0.0;
                for (unsigned short iDonor = 0; iDonor < rect->interpolationStencil[point].size(); ++iDonor) {
                    auto donorPointNumber = rect->interpolationStencil[point][iDonor];
                    auto [_iDonor, _jDonor] = oversetMesh->GetijFromPointNumber(donorPointNumber);
                    rect->Var[k][i][j] = rect->Var[k][i][j] + oversetMesh->Var[k][_iDonor][_jDonor] * rect->interpolationCoeffs[point][iDonor];
                }
                // R = temp - rect->Var[k][i][j];
                R = 0.0;
            }
            rms = rms + R * R;
        }
    }
    ApplyBC(rect, k, oversetMesh);
    rms = sqrt(rms / (rect->Nx * rect->Ny));
    return rms;
}

double SolveP(DataStructure *rect, DataStructure *oversetMesh, int k = 2) {
    /* for P k = 2 (in Var array)*/
    int count = 0;
    double Fd, ap_d, ap, R, rms;
    rms = 0.0;
    for (int i = 1; i < rect->Nx + 1; i++) {
        for (int j = 1; j < rect->Ny + 1; j++) {
            size_t point = rect->GetPointNumber(i, j);
            if (rect->pointType[point] == UNUNSED) {
                continue;
            }
            if (rect->pointType[point] == CALCULATED || rect->pointType[point] == INTERPOLATION_DONOR || rect->pointType[point] == DONOR_BUFFER) {
                DiffusiveFlux(rect, &Fd, &ap_d, i, j, k);
                double LHS = Fd;
                double RHS = rect->rho / rect->dt * (rect->Ff[0][i][j] + rect->Ff[1][i][j] + rect->Ff[2][i][j] + rect->Ff[3][i][j]);
                R = RHS - LHS;
                ap = ap_d;
                count++;
                rect->Var[k][i][j] = rect->Var[k][i][j] + R / ap;
            } else if (rect->pointType[point] == INTERPOLATION_RECIEVER) {
                double temp = rect->Var[k][i][j];
                rect->Var[k][i][j] = 0.0;
                for (unsigned short iDonor = 0; iDonor < rect->interpolationStencil[point].size(); ++iDonor) {
                    auto donorPointNumber = rect->interpolationStencil[point][iDonor];
                    auto [_iDonor, _jDonor] = oversetMesh->GetijFromPointNumber(donorPointNumber);
                    rect->Var[k][i][j] = rect->Var[k][i][j] + oversetMesh->Var[k][_iDonor][_jDonor] * rect->interpolationCoeffs[point][iDonor];
                }
                // R = temp - rect->Var[k][i][j];
                R = 0.0;
            }
            rms = rms + R * R;
        }
    }
    ApplyBC(rect, k, oversetMesh);
    rms = sqrt(rms / (rect->Nx * rect->Ny));
    return rms;
}

void CorrectVelocity(DataStructure *rect) {
    for (int k = 0; k < rect->nVar; k++) {
        rect->residual[k] = 0.0;
    }
    // Correcting Velocity
    int k;
    for (int i = 1; i < rect->Nx + 1; i++) {
        for (int j = 1; j < rect->Ny + 1; j++) {
            size_t point = rect->GetPointNumber(i, j);
            if (rect->pointType[point] == UNUNSED || rect->pointType[point] == INTERPOLATION_RECIEVER) {
                continue;
            }
            k = 0;
            rect->Var[k][i][j] = rect->Var[k][i][j] - rect->dt / rect->rho * (rect->Var[2][i + 1][j] - rect->Var[2][i - 1][j]) / (2 * rect->dx);
            k = 1;
            rect->Var[k][i][j] = rect->Var[k][i][j] - rect->dt / rect->rho * (rect->Var[2][i][j + 1] - rect->Var[2][i][j - 1]) / (2 * rect->dy);

            // calculating residuals for u, v, p
            for (int k = 0; k < rect->nVar; ++k) {
                rect->residual[k] += (rect->Var[k][i][j] - rect->VarOld[k][i][j]) * (rect->Var[k][i][j] - rect->VarOld[k][i][j]);
            }
            // rect->residual[0] += (rect->Var[0][i][j] - rect->VarOld[0][i][j]) * (rect->Var[0][i][j] - rect->VarOld[0][i][j]);
            // rect->residual[1] += (rect->Var[1][i][j] - rect->VarOld[1][i][j]) * (rect->Var[1][i][j] - rect->VarOld[1][i][j]);
            // rect->residual[2] += (rect->Var[2][i][j] - rect->VarOld[2][i][j]) * (rect->Var[2][i][j] - rect->VarOld[2][i][j]);
        }
    }
}

void ImplicitSolve(DataStructure *rect, DataStructure *oversetMesh, int outerIter) {
    // double Fc, ap_c, Fd, ap_d, ap, R, rms;  // No of faces per cell = 4  i=0,1,2,3 -> E,N,W,S

    double rmsOldBg = 1.0, rmsOldComp = 1.0;
    // Solving for U and V
    int iter = 0;
    double rmsBg = 1.0, rmsComp = 1.0;
    const int MAXITER = 50;
    const double UVTolerance = 1e-6;
    for (int k = 0; k < 2; k++) {
        do {
            rmsOldBg = rmsBg;
            rmsOldComp = rmsComp;
            rmsBg = SolveUV(rect, oversetMesh, k);
            rmsComp = SolveUV(oversetMesh, rect, k);
            // cout << "RMS U,V: " << rmsBg << " " << rmsComp << " Iter: " << iter << endl;
            iter++;
            if (iter > MAXITER) {
                break;
            }
        } while ((rmsBg > UVTolerance || rmsComp > UVTolerance) && (rmsOldBg > rmsBg || rmsOldComp > rmsComp));
    }
    // if (outerIter % 100 == 0) {
    // 	cout << "RMS U,V: " << rmsBg << " " << rmsComp << " Iter: " << iter; // << endl;
    // }
    LinearInterpolation(rect);
    LinearInterpolation(oversetMesh);

    // Solving for P
    int k = 2;
    iter = 0;
    rmsBg = 1.0, rmsComp = 1.0;
    rmsOldBg = 1.0, rmsOldComp = 1.0;
    const double PTolerance = 1e-3;
    do {
        rmsOldBg = rmsBg;
        rmsOldComp = rmsComp;
        rmsBg = SolveP(rect, oversetMesh, k);
        rmsComp = SolveP(oversetMesh, rect, k);
        // if (iter % 1000 == 0){
        // 	cout << "RMS P: " << rmsBg << " " << rmsComp << " Iter: " << iter << endl;
        // }
        iter++;
        if (iter > MAXITER) {
            break;
        }
    } while ((rmsBg > PTolerance || rmsComp > PTolerance) && (rmsOldBg > rmsBg || rmsOldComp > rmsComp));
    // if (outerIter % (ITER_PRINT_FREQ/2) == 0) {
    // 	cout << " RMS P: " << rmsBg << " " << rmsComp << " Iter: " << iter << endl;
    // }

    CorrectVelocity(rect);
    CorrectVelocity(oversetMesh);
    for (int k = 0; k < 2; ++k) {
        ApplyBC(rect, k, oversetMesh);
        ApplyBC(oversetMesh, k, rect);
    }

    UpdateFlux(rect);
    UpdateFlux(oversetMesh);
}

bool ConvergenceCheck(DataStructure *rect, int count) {
    double rms[rect->nVar];
    const double TOLERANCE = 1e-10;
    for (int k = 0; k < rect->nVar; k++) {
        rms[k] = sqrt(rect->residual[k] / (rect->Nx * rect->Ny));
        rms[k] = rms[k] / rect->dt;
        if (count % ITER_PRINT_FREQ == 0) {
            cout << std::setprecision(3) << std::scientific << rms[k] << " ";
        }
    }

    if (rms[0] > TOLERANCE || rms[1] > TOLERANCE) {
        CopyNewtoOld(rect);
        return false;
    } else {
        return true;
    }
}

void Solve(DataStructure *rect, DataStructure *oversetMesh) {
    int iter = 0;
    bool mesh1 = true, mesh2 = false;
    do {
        ImplicitSolve(rect, oversetMesh, iter);
        // ImplicitSolve(oversetMesh, rect);
        if (iter % ITER_PRINT_FREQ == 0) {
            cout << "OuterIter: " << std::setw(5) << iter;
        }
        if (iter % ITER_PRINT_FREQ == 0) {
            cout << " NS-RMS: ";
        }
        mesh1 = ConvergenceCheck(rect, iter);
        if (iter % ITER_PRINT_FREQ == 0) {
            cout << " - ";
        }
        mesh2 = ConvergenceCheck(oversetMesh, iter);
        if (iter % ITER_PRINT_FREQ == 0) {
            cout << endl;
        }
        iter++;
    } while ((mesh1 == false || mesh2 == false) && iter < int(1e5));
    cout << "Solved in " << iter << " iterations." << endl;
}

int main() {
    DataStructure bgMesh(0.0, 0.0, 0.0 * PI / 180.0, 1.0, 1.0, 100, 100);
    DataStructure compMesh(0.425, 0.425, 30.0 * PI / 180.0, 0.4, 0.4, 40, 40);
    bgMesh.MarkOversetBoundaryDonors(compMesh);
    bgMesh.MarkBoundaryPointType(BC_SPECIFIED);
    compMesh.MarkBoundaryPointType(INTERPOLATION_RECIEVER);

    bgMesh.StoreInterpolationStencil(compMesh);
    compMesh.StoreInterpolationStencil(bgMesh);

    bgMesh.WriteTxtPointType("bgPtType.txt");
    compMesh.WriteTxtPointType("compPtType.txt");

    initialize(&bgMesh, &compMesh);
    initialize(&compMesh, &bgMesh);
    Solve(&bgMesh, &compMesh);
    GetOutput(&bgMesh, "output_Upwind_bgMesh.dat");
    GetOutput(&compMesh, "output_Upwind_compMesh.dat");
    return 0;
}