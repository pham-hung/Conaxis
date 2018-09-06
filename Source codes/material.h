#ifndef MATERIAL_H
#define MATERIAL_H

#include <QDialog>
#include <Eigen/Dense>
#include <iostream>
#include <QFileDialog>
#include <QString>
#include <vector>
#include <QMessageBox>
#include <QInputDialog>

#include "getfile.h"
#include "BaseClass/materialbase.h"

using namespace std;
using namespace Eigen;

namespace Ui {
class Material;
}

class Material : public QDialog,private MaterialBase
{
    Q_OBJECT

public:
    explicit Material(QWidget *parent = 0);
    ~Material();
    vector<MaterialBase> mat;
    void getUserData();
    void sendSIGNAL();
    int NumberOfMat();
    void createCrsMaterial(double poissionRatio,double voidRatio);

private slots:
    void on_cancelButton_clicked();
    void on_addButton_clicked();
    void on_KCombo_activated(int index);
    void on_listInfor_clicked();
    void on_matNum_textChanged(const QString &arg1);
    void on_kCombo_activated(int index);
    void listMaterial();
    void getMaterial(vector<MaterialBase> mat);
    void on_setCd_clicked();
    void on_setRatio_clicked();

signals:
    void sendMaterial(vector<MaterialBase>);

private:
    Ui::Material *ui;
    bool checkSuccess=false;
};

#endif // MATERIAL_H
