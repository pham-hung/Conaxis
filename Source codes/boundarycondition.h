#ifndef BOUNDARYCONDITION_H
#define BOUNDARYCONDITION_H

#include <QDialog>
#include <QFileDialog>
#include <QInputDialog>
#include <QMessageBox>
#include <QString>
#include <Eigen/Dense>
#include <vector>
#include "BaseClass/boundaryconditionbase.h"
#include "getfile.h"
#include <QFile>

using namespace std;
using namespace Eigen;

namespace Ui {
class BoundaryCondition;
}

class BoundaryCondition : public QDialog,private BoundaryConditionBase
{
    Q_OBJECT

public:
    explicit BoundaryCondition(QWidget *parent = 0);
    ~BoundaryCondition();
    void getUserData();
    void getDatafromVector(int i);
    void showData();
    void defaultData();
    void sendSignal();
    void createCrsBoundary(Ref<MatrixXd> testData);

private slots:
    void on_closeButton_clicked();
    void on_typeCombo_activated(int boundaryIndex);
    void on_loadTypeCombo_activated(int boundaryIndex);
    void on_indexLine_cursorPositionChanged(int arg1, int arg2);
    void on_listButton_clicked();
    void on_okButton_clicked();
    void on_indexLine_textChanged(const QString &arg1);

public slots:
    void getBoundary(vector<BoundaryConditionBase> boundary);

signals:
    void sendParameters(vector<BoundaryConditionBase>);

private:
    Ui::BoundaryCondition *ui;
    vector<BoundaryConditionBase> boundary;
    bool successCheck=true;
    GetFile getfile;
};

#endif // BOUNDARYCONDITION_H
