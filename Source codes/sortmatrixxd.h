#ifndef SORTMATRIXXD_H
#define SORTMATRIXXD_H
#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <QVector3D>
#include "math.h"

using namespace std;
using namespace Eigen;

struct element
{
    int rowIndex;
    int colIndex;
    double Val;
};

bool sortByVal(const element &lhs,const element &rhs);


class SortMatrixXd
{
public:
    SortMatrixXd();
    void sortByCol(Ref<MatrixXd> matrixA,int col);
    void sortByRow(Ref<MatrixXd> matrixA,int row);
    void findXY(const Ref<const MatrixXd> coordinates,double x0, double x1, double y0, double y1, MatrixXd &result);    
    void findX(const Ref<const MatrixXd> coordinates,double x0, double x1, MatrixXd &result);
    void findY(const Ref<const MatrixXd> coordinates,double y0, double y1, MatrixXd &result);
    void addMatrix(MatrixXd &sourceMatrix, const Ref<const MatrixXd> addedMatrix, int col);
    void addMatrixReplace(MatrixXd &sourceMatrix, const Ref<const MatrixXd> addedMatrix, int col);
    void findElementXY(const Ref<const MatrixXd> coordinates, const Ref <const MatrixXd> elements, double x0, double x1, double y0, double y1, MatrixXd &result);
    void swapValue(double &x, double &y);
    void PlaneLineIntersect(Vector3d PointP0, Vector3d PointP1, Vector3d PointV0, Vector3d nVector,Vector3d IntersectionVector);

private:
    vector<element> vectorA;
};

#endif // SORTMATRIXXD_H
