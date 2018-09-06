#include "scalesetting.h"
#include "ui_scalesetting.h"

ScaleSetting::ScaleSetting(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::ScaleSetting)
{
    ui->setupUi(this);
    ui->xScaleLine->setText("1");
    ui->yScaleLine->setText("1");
    ui->resultScaleLine->setText("1");
    ui->stepLine->setText("1");
    ui->autoMaxMinBox->setChecked(true);
    ui->maxValLine->setText("10");
    ui->minValLine->setText("-10");
}

ScaleSetting::~ScaleSetting()
{
    delete ui;
}

void ScaleSetting::getUserData()
{
    xScale=ui->xScaleLine->text().toDouble();
    yScale=ui->yScaleLine->text().toDouble();
    resultScale=ui->resultScaleLine->text().toDouble();
    step=ui->stepLine->text().toInt();
    if(step<1)
    {
        step=1;
    }
    if(ui->autoMaxMinBox->isChecked())
    {
        autoMaxMinVal=true;
    }
    else
    {
        autoMaxMinVal=false;
    }
    minVal=ui->minValLine->text().toDouble();
    maxVal=ui->maxValLine->text().toDouble();
    if(ui->lockViewCheckBox->isChecked())
    {
        lockView=true;
    }
    else
    {
        lockView=false;
    }
}

void ScaleSetting::on_okButton_clicked()
{
    getUserData();
    scale.xScale=xScale;
    scale.yScale=yScale;
    scale.resultScale=resultScale;
    scale.step=step;
    scale.autoMaxMinVal=autoMaxMinVal;
    scale.minVal=minVal;
    scale.maxVal=maxVal;
    scale.lockView=lockView;
    emit sendSignal(scale);
    this->close();
}

void ScaleSetting::on_cancelButton_clicked()
{
    this->close();
}
