#ifndef SCALESETTING_H
#define SCALESETTING_H

#include <QDialog>
#include "BaseClass/scalesettingbase.h"

namespace Ui {
class ScaleSetting;
}

class ScaleSetting : public QDialog,private ScaleSettingBase
{
    Q_OBJECT

public:
    explicit ScaleSetting(QWidget *parent = 0);
    ~ScaleSetting();
    void getUserData();

private slots:
    void on_okButton_clicked();
    void on_cancelButton_clicked();
signals:
    void sendSignal(ScaleSettingBase scale);

private:
    Ui::ScaleSetting *ui;
    ScaleSettingBase scale;
};

#endif // SCALESETTING_H
