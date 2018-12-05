#ifndef GEOMETRY_H
#define GEOMETRY_H
#include <Eigen/Dense>
#include <QDialog>
#include <QMessageBox>
#include <iostream>
#include <QDebug>
#include <QTableWidgetItem>
#include <vector>
#include <algorithm>
#include "writetofile.h"
#include "stdio.h"
#include "math.h"
#include "BaseClass/geometrybase.h"

using namespace std;
using namespace Eigen;

namespace Ui {
class Geometry;
}

class Geometry : public QDialog
{
    Q_OBJECT

public:
    explicit Geometry(QWidget *parent = 0);
    ~Geometry();
    void updateInformation();
    void createMeshSmear();
    void createMeshNoSmear();
    void shownInformation();
    void saveDataToObject();

    bool checkInformation();
    void replot();
    bool isFound(pair<double,double> findPair, pair<double,double> soucePair);
    int findNodeIndex(double x, double y, vector<pair<double,double> > &XYCoordinates);
    void pauseSystem();
    void shownPair(vector<pair<double,double> > &XYCoordinates);
    void createCrsMesh(double H0,double R);

private slots:
    void on_updateButton_clicked();
    void on_meshButton_clicked();
    void on_analysisBox_currentIndexChanged(int index);
    void on_setSubLayer_clicked();
    void on_replotButton_clicked();
    void getGeometryBaseObject(GeometryBase geometryObject);

signals:
    void sendGeometryData(Ref<MatrixXd> elements, Ref<MatrixXd>coordinates,int analysisType, double qw);
    void sendGeometryObject(GeometryBase geometryObject);

private:
    Ui::Geometry *ui;
    double re,rw,rs;
    double length;
    int numberOfLayer;
    double qw;
    int numberOfElementSmear;
    int numberOfElementSoil;
    double surfaceElevation;
    int analysisType;
    MatrixXd coordinates;
    MatrixXd elements;
    MatrixXd layerInfo;
    bool checkValid;
    bool clickFlag=false;
    int defaultSubLayer;
    WriteToFile exportFile;

    vector<pair<double,double> > XYCoordinates;
    vector<pair<int,int> > node1D;
    GeometryBase geometryObject;
};

#endif // GEOMETRY_H
