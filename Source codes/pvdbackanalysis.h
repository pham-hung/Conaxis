#ifndef PVDBACKANALYSIS_H
#define PVDBACKANALYSIS_H

#include <QDialog>

namespace Ui {
class PVDBackAnalysis;
}

class PVDBackAnalysis : public QDialog
{
    Q_OBJECT

public:
    explicit PVDBackAnalysis(QWidget *parent = 0);
    ~PVDBackAnalysis();
    void getUserData();
signals:
    void sendPVDsParameters(double requi,double rw,double rs,double ratioKs,double Cd0,double p0, bool NoSmear,bool goldenSearchFlag);

public slots:
    void getResult(double Cd,double error);

private slots:
    void on_pushButton_clicked();
    void on_pushButton_2_clicked();
    void on_runCd0Button_clicked();
    void on_comboBox_activated(int index);
private:
    Ui::PVDBackAnalysis *ui;
    double requi, rw, rs, ratioKs, Cd0, p0;
    double Cd, currentError;
    bool NoSmear;
    bool goldenSearchFlag=false;
};

#endif // PVDBACKANALYSIS_H
