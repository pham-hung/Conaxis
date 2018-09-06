#ifndef CRSBACKANALYSIS_H
#define CRSBACKANALYSIS_H

#include <QDialog>
#include "getfile.h"
#include "writetofile.h"
#include <Eigen/Dense>
#include "QFileDialog"
#include "QFile"
#include "QDir"
#include "QMessageBox"
#include "QInputDialog"
#include <QString>
#include <QStringList>
#include "math.h"
#include "BaseClass/crs_base.h"

using namespace std;
using namespace Eigen;

namespace Ui {
class CRSBackAnalysis;
}

class CRSBackAnalysis : public QDialog
{
    Q_OBJECT

public:
    explicit CRSBackAnalysis(QWidget *parent = 0);
    ~CRSBackAnalysis();
    void getUserData();
    void linearTheory();
    void nonLinearTheory();
    void exportData();
    void pauseSystem();
    void prepareSimulationData();

private slots:
    void on_ASTMButton_clicked();
    void on_closeButton_clicked();
    void on_browseLine_clicked();
    void on_runCrsButton_clicked();
    void on_runBackAnalysis_clicked();

signals:
    void sendCrsData(CRS_Base crsObject,bool backAnalysisFlag);

private:
    Ui::CRSBackAnalysis *ui;

    QString fileName;
    QString folderName;
    MatrixXd testData, resultData;
    int timeCol=1, stressCol=3, poreCol=4, strainCol=2;
    double H0=0.0254;
    double Hs;
    double R=0.03175;
    double e0=1.84;
    double poission=0.4;
    int outputType=0;
    int theoryType=0;
    bool sucessFlag=false;
    CRS_Base crsObject;
    bool backAnalysisFlag=false;
};

#endif // CRSBACKANALYSIS_H
