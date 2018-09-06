#ifndef SOLUTIONCONTROL_H
#define SOLUTIONCONTROL_H

#include <QDialog>
#include "BaseClass/outputexportbase.h"
#include <iostream>

using namespace std;

namespace Ui {
class SolutionControl;
}

class SolutionControl : public QDialog,private OutputExportBase
{
    Q_OBJECT

public:
    explicit SolutionControl(QWidget *parent = 0);
    ~SolutionControl();
    void getUserData();
    void assignData();
    void updateData();

private slots:
    void on_cancelButton_clicked();
    void on_okButton_clicked();

signals:
    void sendSoltuonParameters(OutputExportBase solutionParameters);

private:
    Ui::SolutionControl *ui;
    OutputExportBase solutionParameters;
};

#endif // SOLUTIONCONTROL_H
