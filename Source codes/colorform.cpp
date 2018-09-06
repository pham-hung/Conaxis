#include "colorform.h"
#include "ui_colorform.h"

colorForm::colorForm(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::colorForm)
{
    ui->setupUi(this);
    defaultValue();
    setCurrentValue();
}

colorForm::~colorForm()
{
    delete ui;
}

void colorForm::setCurrentValue()
{
    ui->nocLine->setText(QString::number(noc));
    ui->sizeLine->setText(QString::number(fontSize));
    ui->numericBox->setCurrentIndex(numType);
    ui->positionBox->setCurrentIndex(colorPosition);
    if(nodePlot==1){ui->nodeBox->setChecked(true);}
    if(meshPlot==1){ui->meshBox->setChecked(true);}
    if(resultPlot==1){ui->resultBox->setChecked(true);}
    if(axePlot==1){ui->axeBox->setChecked(true);}
    if(datePlot==1){ui->dateBox->setChecked(true);}
    if(titlePlot==1){ui->titleBox->setChecked(true);}
    if(valPlot==1){ui->valBox->setChecked(true);}
    ui->titleLine->setText(title);
}

void colorForm::defaultValue()
{
    nodePlot=1;
    meshPlot=1;
    resultPlot=1;
    axePlot=1;
    datePlot=0;
    titlePlot=1;
    valPlot=1;
}

void colorForm::on_pushButton_clicked()
{
    noc=QString(ui->nocLine->text()).toInt();
    numType=ui->numericBox->currentIndex();
    fontSize=QString(ui->sizeLine->text()).toInt();
    colorPosition=ui->positionBox->currentIndex();
    if(ui->nodeBox->isChecked())
    {
        nodePlot=1;
    }
    else
    {
        nodePlot=0;
    }

    if(ui->meshBox->isChecked())
    {
        meshPlot=1;
    }
    else
    {
        meshPlot=0;
    }

    if(ui->resultBox->isChecked())
    {
        resultPlot=1;
    }
    else
    {
        resultPlot=0;
    }

    if(ui->dateBox->isChecked())
    {
        datePlot=1;
    }
    else
    {
        datePlot=0;
    }

    if(ui->titleBox->isChecked())
    {
        titlePlot=1;
    }
    else
    {
        titlePlot=0;
    }

    if(ui->valBox->isChecked())
    {
        valPlot=1;
    }
    else
    {
        valPlot=0;
    }

    if(ui->axeBox->isChecked())
    {
        axePlot=1;
    }
    else
    {
        axePlot=0;
    }

    title=ui->titleLine->text();

    emit sendColorBandInfor(noc,numType,fontSize,colorPosition,nodePlot,meshPlot,resultPlot,axePlot,datePlot,titlePlot,valPlot,title);
    this->close();
}

void colorForm::on_pushButton_2_clicked()
{
    this->close();
}
