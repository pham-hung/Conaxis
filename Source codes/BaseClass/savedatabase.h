#ifndef SAVEDATABASE_H
#define SAVEDATABASE_H

#include <QObject>
#include <Eigen/Dense>
#include <iostream>
#include <QFile>
#include <QDir>
#include <QFileDialog>
#include <QTextStream>
#include <QString>
#include <QInputDialog>
#include <QStringList>
#include <QDebug>

#include "projectsetting.h"
#include "material.h"
#include "stage.h"
#include "assignboundarycondition.h"
#include "axissymmetric_2d.h"
#include "boundarycondition.h"
#include "watchlist.h"

#include "BaseClass/materialbase.h"
#include "BaseClass/stagebase.h"
#include "BaseClass/boundaryconditionbase.h"
#include "BaseClass/stageboundarybase.h"
#include "BaseClass/watchlistbase.h"
#include "BaseClass/geometrybase.h"

using namespace std;


class SaveDataBase : public QObject,private BoundaryConditionBase
{
    Q_OBJECT
public:
    explicit SaveDataBase(QObject *parent = nullptr);
    void saveModel();
    void saveAsModel();
    void readModel();
    void pauseSystem();
    void listMaterial();
    void listStage();
    void listBoundary();
    void listStageBoundary();
    void listWatch();
    void sendSIGNAL();
    bool checkReadingState(){return success;}
    void resetFolderName();

signals:
    void sendProjectSetting(vector<double> projectParameters);
    void sendMesh(Ref<MatrixXd> coordinates,Ref<MatrixXd> elements);
    void sendMaterial(vector<MaterialBase> material);
    void sendStage(vector<StageBase> stage);
    void sendBoundary(vector<BoundaryConditionBase> boundary);
    void sendStageBoundary(vector<StageBoundaryBase> stageBoundary);
    void sendWatchListBase(vector<WatchListBase> watch);
    void sendGeometry(GeometryBase geometry);

public slots:
    void getMesh(Ref<MatrixXd> coordinates,Ref<MatrixXd> elements,QString folderName);
    void getStageBoundary(vector<StageBoundaryBase> stageBoundary);
    void getMaterial(vector<MaterialBase> material);
    void getStage(vector<StageBase> stage);
    void getBoundary(vector<BoundaryConditionBase> boundary);
    void getProjectSetting(vector<double> projectParameters);
    void getWatchListBase(vector<WatchListBase> watch);
    void getGeometryBase(GeometryBase geometry);

private:
    bool saveClicked=false;
    QString folderName,fileName;
    MatrixXd coordinates, elements;
    vector<MaterialBase> material;
    vector<StageBase> stage;
    vector<BoundaryConditionBase> boundary;
    vector<StageBoundaryBase> stageBoundary;
    vector<double> projectParameters;
    MatrixXd tempMatrix;
    vector<MatrixXd> vectorNode;
    vector<WatchListBase> watch;
    GeometryBase geometry;

    QFile mFile;
    QDir mDir;
    bool success;
};

#endif // SAVEDATABASE_H
