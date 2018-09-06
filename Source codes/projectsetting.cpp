#include "projectsetting.h"
#include "ui_projectsetting.h"

ProjectSetting::ProjectSetting(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::ProjectSetting)
{
    ui->setupUi(this);
    projectParameters.resize(4);
    projectParameters[0]=analysisType;
    projectParameters[1]=solverType;
    projectParameters[2]=elementType;
    projectParameters[3]=acel;
}

ProjectSetting::~ProjectSetting()
{
    delete ui;
}

void ProjectSetting::assignData()
{
    ui->analysisCombo->setCurrentIndex(analysisType);
    ui->solverCombo->setCurrentIndex(solverType);
    ui->elementCombo->setCurrentIndex(elementType);
    ui->accelLine->setText(QString::number(acel));
}

void ProjectSetting::sendSIGNAL()
{
    if(projectParameters.size()==0)
    {
        projectParameters.resize(4);
        projectParameters[0]=analysisType;
        projectParameters[1]=solverType;
        projectParameters[2]=elementType;
        projectParameters[3]=acel;
    }
    emit sendParameters(projectParameters);
}

void ProjectSetting::on_pushButton_2_clicked()
{
    this->close();
}

void ProjectSetting::on_okButton_clicked()
{
    projectParameters.resize(4);
    analysisType=ui->analysisCombo->currentIndex();
    solverType=ui->analysisCombo->currentIndex();
    elementType=ui->analysisCombo->currentIndex();
    acel=ui->accelLine->text().toDouble();
    projectParameters[0]=analysisType;
    projectParameters[1]=solverType;
    projectParameters[2]=elementType;
    projectParameters[3]=acel;
    assignData();
    this->close();
}

void ProjectSetting::getProjectParameters(vector<double> projectParameters)
{
    this->projectParameters=projectParameters;
}
