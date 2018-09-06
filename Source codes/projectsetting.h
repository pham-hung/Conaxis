#ifndef PROJECTSETTING_H
#define PROJECTSETTING_H

#include <QDialog>
#include <vector>
#include <iostream>
using namespace std;

namespace Ui {
class ProjectSetting;
}

class ProjectSetting : public QDialog
{
    Q_OBJECT

public:
    explicit ProjectSetting(QWidget *parent = 0);
    ~ProjectSetting();
    void getUserData();
    void assignData();
    void sendSIGNAL();

private slots:
    void on_pushButton_2_clicked();
    void on_okButton_clicked();
    void getProjectParameters(vector<double> projectParameters);
signals:
    void sendParameters(vector <double>);

private:
    Ui::ProjectSetting *ui;
    int analysisType=0;
    int solverType=0;
    int elementType=0;
    double acel=9.81;
    vector <double> projectParameters;
};

#endif // PROJECTSETTING_H
