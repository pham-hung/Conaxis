#ifndef COLORFORM_H
#define COLORFORM_H

#include <QDialog>
#include <QString>
#include <QLineEdit>
#include <QComboBox>
#include <QPushButton>
#include <QObject>

using namespace std;

namespace Ui {
class colorForm;
}

class colorForm : public QDialog
{
    Q_OBJECT

public:
    explicit colorForm(QWidget *parent = 0);
    ~colorForm();
    void setCurrentValue();
    void defaultValue();
signals:
    void sendColorBandInfor(int noc,int numType,int fontSize,int colorPosition,
                            int nodePlot,int meshPlot,int resultPlot,int axePlot,int datePlot,int titlePlot,int valPlot,QString title);

private slots:
    void on_pushButton_clicked();
    void on_pushButton_2_clicked();

private:
    Ui::colorForm *ui;
    int noc=10;
    int numType=1;
    int fontSize=10;
    int colorPosition=0;
    int nodePlot,  meshPlot, resultPlot, axePlot, datePlot, titlePlot, valPlot;
    QString title="PROJECT";
};

#endif // COLORFORM_H
