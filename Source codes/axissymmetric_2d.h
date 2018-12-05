#ifndef AXISSYMMETRIC_2D_H
#define AXISSYMMETRIC_2D_H

#include <QObject>
#include <QFile>
#include <QString>
#include <QElapsedTimer>
#include <QDebug>
#include <QDir>

#include <iostream>
#include "math.h"
#include <stdio.h>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/PardisoSupport>
#include <Eigen/Sparse>

#include "getfile.h"
#include "writetofile.h"
#include "gauss2d.h"

#include "pvdlib.h"
#include "sortmatrixxd.h"

#include "projectsetting.h"
#include "material.h"
#include "assignboundarycondition.h"
#include "stage.h"
#include "watchlist.h"

#include "boundarycondition.h"
#include "BaseClass/boundaryconditionbase.h"
#include "BaseClass/stageboundarybase.h"
#include "BaseClass/stagebase.h"
#include "BaseClass/materialbase.h"
#include "BaseClass/baseelement.h"
#include "BaseClass/watchlistbase.h"


class AxisSymmetric_2D : public QObject
{
    Q_OBJECT
public:
    explicit AxisSymmetric_2D(QObject *parent = nullptr);
    //functions
    void prepareData();
    void initializeData();
    void assemblyGlobalMatrix();
    void solveDirect();
    void calculateStress();
    void calculateStress(int calculatedStep);    

    //commons functions
    void compareMatrix(Ref<MatrixXi> matrixA,Ref<MatrixXi> matrixB,Ref<MatrixXi> matrixC);
    void calculateNodfat();
    int countPoreBoundary(Ref<MatrixXd> matrixA,Ref<MatrixXi> nodfmt);
    void createNewPoreBoundary(Ref<MatrixXd>matrixOld, Ref<MatrixXd>matrixNew,Ref<MatrixXi> nodfmt);
    void createGlobalBoundary(Ref<MatrixXd>LocalMatrix, Ref<MatrixXd>GlobalMatrix);
    double binarySearch(const Ref<const MatrixXd> searchMatrix,int colSearch,int colResult, double inVal);

    //run functions
    void runAllAnalysis();
    void exportResults();
    void exportStress();
    void averageStress(); //Export Average stress
    void exportMaterialsParameters();
    void getPreviousStep();

    //PVDs function
    void runPVDs();
    void setPVDparameters();
    void calculateAnalyticalPVD();
    double calculatePVDError();
    double PVDError(double Cd);
    double goldenSearch(double a, double b);
    void resetSolveData();

    //Crs Test function
    void runCrsTest();
    void runBackAnalysis();
    void PlaneLineIntersect();
    void getBackAnalysisResults();
    void toVectorSolutions();
    void assemblyBackAnalysis(double K, double k);
    double averageEffectiveStress(int currentStep);

    //boundary condition functino       
    void createBoundaryConditionEachSubStep(int currentStage, double realTime);
    void createDichletBoundaryCondition();
    void createDirichletBoundaryConditionIncrementalForm();
    void createLoadEachStep(int gravityCheck);

    void createGravityLoad();
    void gravityLoadTri6p(int &eleNum);
    void gravityLoadQuad8p(int &eleNum);

    void transferLineLoadToNodalLoad(Ref<MatrixXd> nodeList, double value);
    void transferLineLoadToNodalLoadXDirection(Ref<MatrixXd> nodeList, double value);
    int countPoreDegreeOfFreedom(const Ref<const MatrixXd> boundaryMatrix, const Ref<const MatrixXi> nodfmt);

    //commons functions
    double getDispTop();
    double getPorePressure();
    void setCdFactor(double Cd);
    double getError();
    void pauseSystem();

    //element functions
    void tri6pMatrix(int &eleNum);
    void quad8pMatrix(int &eleNum);
    void line2pMatrix(int &eleNum);

    //Stress functions
    void calculateStressTri6p(int &eleNum);
    void calculateStressQuad8p(int &eleNum);

    //watchList
    double calculateWatchList(int index, int step); //index is 1-based, step is 1-based    
    void exportWatchList();

signals:
    void sendSolverInformation(Ref<MatrixXd> U, Ref<MatrixXd> V, Ref<MatrixXd> P);
    void sendPVDResults(double Cd, double error);
    //post-processing


public slots:
    void setFolder(QString folderName);    
    void getMesh(Ref<MatrixXd> coordinates,Ref<MatrixXd> elements,QString folderName);    
    void getStageBoundary(vector<StageBoundaryBase> stageBoundary);
    void getMaterial(vector<MaterialBase> material);
    void getStage(vector<StageBase> stage);
    void getBoundary(vector<BoundaryConditionBase> boundary);
    void getProjectSetting(vector<double> projectParameters);
    void getWatchList(vector<WatchListBase> watch);
    void getAndExportWatchList(vector<WatchListBase> watch);
    void getPVDsParameters(double requi,double rw,double rs,double ratioKs,double Cd0,double p0, bool NoSmear,bool goldenSearchFlag);
    void getAdditionalInformation(int analysisType,double A1D, double k1D);
    void getCrsData(Ref<MatrixXd> crsData, bool crsFlag,int crsType);

private:
    //-----------------------------------------------
    //Input Data
    vector<StageBoundaryBase> stageBoundary;
    vector<MaterialBase> material;
    vector<StageBase> stage;
    vector<BoundaryConditionBase> boundary;
    vector<double> projectParameters;
    vector<WatchListBase> watch;

    //----------------------------------------------
    QElapsedTimer timerGlobal, timerAssembly, timer;
    QString fileName, folderName;
    double Ke, Ge, ve, re, kve, khe, Se, dense, voidRatioe;    
    int non,noe,dof,doff,totalDof,step,ns;
    int nos; //number of stages;
    int nob; //nubmer of boundaries
    double fy,fx,dt,dt0,press0;
    int eleNum, eleType;    
    int analysisType;
    double k1D, A1D;

    //--------------------------
    //For PVDs analysis
    double acel,gf;
    double Cf, Cs, BiotCoeff,SkemptonCoeff;
    double error;
    double Cd0, Cd;
    double requi,rw,rs,ratioKs,Cvr,ks;
    bool NoSmear;
    MatrixXd Panalytical, Perror, PerrorEachStep;
    bool PVDFlag=false;
    bool goldenSearchFlag=false;
    MatrixXd calculationTime;
    bool watchListFlag=false;

    //for CRS simulation
    MatrixXd crsData;
    bool crsFlag=false;
    MatrixXd resultCrs;
    int crsType;

    //for back analysis
    double realStress, realPore, modelStress, modelPore;
    Vector3d n_in,V0_in,P0_in,P1_in,I_in,u_in, w_in;
    double y0Pore, y0Stress, y1Pore, y1Stress;
    double tol; //tol for stop back analysis
    double Js, Jp; //Error of stress, error of pore pressure
    MatrixXd backAnalysisResult;
    //--------------------------

    MatrixXd coordinates,elements,fixx,fixy,fixh;
    MatrixXd DirichletAll;
    MatrixXd Dirichlet, Dirichlet0, dDirichlet;

    MatrixXd FyGravity;
    bool checkOK=true;
    MatrixXd Fx, Fy;

    vector<int> nodeFind;
    MatrixXd findList;
    vector<double> xCoordFind;
    vector<double> yCoordFind;
    vector<MatrixXd> watchListResult;

    MatrixXi nodfat, nodfmt;
    MatrixXd X0, XX, X, F,F0,dF;
    MatrixXd U, V, P, PizoResult;
    MatrixXd Uinsitu, Vinsitu,Pinsitu;    
    MatrixXd HydraulicVertial, HydraulicHorizontal, BulkModulus;

    MatrixXd Sxx, Syy, Sxy, Pore;
    MatrixXd SyyInterpolation, SyyGravity;
    int step_i;
    int currentStep;
    double realCalculationTime;

    typedef Eigen::Triplet<double> Trip;
    std::vector<Trip> trip_total;
    SparseMatrix<double,RowMajor> KK;

    GetFile getfile;
    WriteToFile exportFile;
    BaseElement baseEle;
    Gauss2D gauss;
    PvdLib analyTest;
    SortMatrixXd matrixLib;

    //Watch list result
    vector<double> foundNodeResult;
};

#endif // AXISSYMMETRIC_2D_H
