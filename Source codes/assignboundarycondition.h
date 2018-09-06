#ifndef ASSIGNBOUNDARYCONDITION_H
#define ASSIGNBOUNDARYCONDITION_H

#include <QDialog>
#include <Eigen/Dense>
#include <iostream>
#include <Eigen/Dense>
#include <QInputDialog>
#include <QFileDialog>
#include <QString>
#include <vector>
#include <QDebug>

#include "stage.h"
#include "BaseClass/stagebase.h"
#include "boundarycondition.h"
#include "BaseClass/boundaryconditionbase.h"
#include "getfile.h"
#include "BaseClass/stageboundarybase.h"

using namespace std;
using namespace Eigen;

namespace Ui {
class AssignBoundaryCondition;
}

class AssignBoundaryCondition : public QDialog,private StageBase, private BoundaryConditionBase
{
    Q_OBJECT

public:
    explicit AssignBoundaryCondition(QWidget *parent = 0);
    ~AssignBoundaryCondition();
    void updateData();
    void resetData();
    void getUserData();
    void assignUI();
    void defaultData();
    void sendSIGNAL();
    void getDatafromVector(int i,int j);
    void createCrsStageBoundary(double H, double R, int controlType);
    bool checkSuccess=true;

public slots:
    void getParametersStage(vector<StageBase> stage);
    void getParametersBoundary(vector<BoundaryConditionBase> boundary);
    void getStageBoundary(vector<StageBoundaryBase> stageBoundary);

private slots:
    void on_closeButton_clicked();
    void on_addButton_clicked();
    void on_assignTypeCombo_activated(int index);
    void on_stageCombo_activated(int index);
    void on_boundaryCombo_activated(int index);
    void on_stageBoundaryIndexLine_textChanged(const QString &arg1);
    void on_listButton_clicked();

    void on_resetButton_clicked();

signals:
    void sendParametersStageBoundary(vector<StageBoundaryBase>);

private:
    Ui::AssignBoundaryCondition *ui;
    //----------------
    vector<StageBoundaryBase> stageBoundary;

    //---------------
    double tol=1e-3;
    MatrixXd coordinates;
    vector<StageBase> stage;
    vector<BoundaryConditionBase> boundary;
    MatrixXd nodeListTemp=MatrixXd::Zero(0,0);
    //-------------
    QString fileName;
    int assignType;
    int stageBoundaryIndex;
    double x0, x1, y0, y1;

};

#endif // ASSIGNBOUNDARYCONDITION_H
