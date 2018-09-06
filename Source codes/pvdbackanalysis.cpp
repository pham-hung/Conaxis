#include "pvdbackanalysis.h"
#include "ui_pvdbackanalysis.h"

PVDBackAnalysis::PVDBackAnalysis(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::PVDBackAnalysis)
{
    ui->setupUi(this);
}

PVDBackAnalysis::~PVDBackAnalysis()
{
    delete ui;
}

void PVDBackAnalysis::getUserData()
{
    requi=ui->requiLine->text().toDouble();
    rw=ui->rwLine->text().toDouble();
    rs=ui->rsLine->text().toDouble();
    ratioKs=ui->ratioSmearLine->text().toDouble();
    if(ui->comboBox->currentIndex()==0)
    {
        NoSmear=false;
    }
    else
    {
        NoSmear=true;
    }
    Cd0=ui->Cd0Line->text().toDouble();
    p0=ui->p0Line->text().toDouble();
}

void PVDBackAnalysis::getResult(double Cd, double error)
{
    this->Cd=Cd;
    this->currentError=error;
    ui->CdLine->setText(QString::number(Cd));
    ui->errorLine->setText(QString::number(currentError));
    QCoreApplication::processEvents();
}

//Find Cd
void PVDBackAnalysis::on_pushButton_clicked()
{
    getUserData();
    goldenSearchFlag=true;
    emit sendPVDsParameters(requi,rw,rs,ratioKs,Cd0,p0,NoSmear,goldenSearchFlag);
}

//Close button
void PVDBackAnalysis::on_pushButton_2_clicked()
{
    this->close();
}

//Run with Cd0
void PVDBackAnalysis::on_runCd0Button_clicked()
{
    getUserData();
    goldenSearchFlag=false;
    emit sendPVDsParameters(requi,rw,rs,ratioKs,Cd0,p0,NoSmear,goldenSearchFlag);
}

void PVDBackAnalysis::on_comboBox_activated(int index)
{
    if(index==1)
    {
        ui->rsLine->setText(ui->rwLine->text());
        ui->ratioSmearLine->setText("1");
    }
}
