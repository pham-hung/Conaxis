#ifndef STAGE_H
#define STAGE_H

#include <QDialog>
#include <QString>
#include <iostream>
#include <vector>
#include <QMessageBox>
#include <QFileDialog>
#include <QInputDialog>
#include "getfile.h"
#include "BaseClass/stagebase.h"

using namespace std;

namespace Ui {
class Stage;
}

class Stage : public QDialog,private StageBase
{
    Q_OBJECT

public:
    explicit Stage(QWidget *parent = 0);
    ~Stage();


    void getUserData();
    void getDataFromVector(int i);
    void showData();
    void defaultData();
    void updateTime(int i);
    void sendSignal();
    void createCrsStage(Ref<MatrixXd> testTime);

private slots:
    void on_cancelButton_clicked();
    void on_listButton_clicked();
    void on_addButton_clicked();
    void on_stageIndexLine_textChanged(const QString &arg1);
    void on_typeCombo_activated(int index);
    void on_dtLine_textChanged(const QString &arg1);
    void getStage(vector<StageBase> stage);    

    void on_timeStepCombo_activated(int index);



signals:
    void sendParameters(vector<StageBase>);

private:
    Ui::Stage *ui;
    bool successCheck=true;
    vector<StageBase> stage;
};

#endif // STAGE_H
