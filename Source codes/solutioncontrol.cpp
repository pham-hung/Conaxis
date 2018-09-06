#include "solutioncontrol.h"
#include "ui_solutioncontrol.h"

SolutionControl::SolutionControl(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::SolutionControl)
{
    ui->setupUi(this);    
}

SolutionControl::~SolutionControl()
{
    delete ui;
}

void SolutionControl::getUserData()
{
    solveType=ui->solveCombo->currentIndex();
    tolerance=ui->toleranceLine->text().toDouble();
    maxIteration=ui->maximumIterationLine->text().toInt();
    if(ui->nodalSolutionBox->isChecked())
    {
        nodalSolutionCheck=1;
    }
    else
    {
        nodalSolutionCheck=0;
    }

    if(ui->elementSolutionBox->isChecked())
    {
        elemetStressCheck=1;
    }
    else
    {
        elemetStressCheck=0;
    }

    if(ui->averageElemenBox->isChecked())
    {
        averageStressCheck=1;
    }
    else
    {
        averageStressCheck=0;
    }

    if(ui->materialParameterBox->isChecked())
    {
        materialParameterCheck=1;
    }
    else
    {
        materialParameterCheck=0;
    }
}

void SolutionControl::on_cancelButton_clicked()
{
    this->close();
}

void SolutionControl::assignData()
{
    solutionParameters.solveType=solveType;
    solutionParameters.tolerance=tolerance;
    solutionParameters.maxIteration=maxIteration;
    solutionParameters.nodalSolutionCheck=nodalSolutionCheck;
    solutionParameters.elemetStressCheck=elemetStressCheck;
    solutionParameters.averageStressCheck=averageStressCheck;
    solutionParameters.materialParameterCheck=materialParameterCheck;
}

void SolutionControl::updateData()
{
    ui->solveCombo->setCurrentIndex(solveType);
    ui->toleranceLine->setText(QString::number(tolerance));
    ui->maximumIterationLine->setText(QString::number(maxIteration,'f',0));
    if(nodalSolutionCheck==1){ui->nodalSolutionBox->setChecked(true);}
    if(elemetStressCheck==1){ui->elementSolutionBox->setChecked(true);}
    if(averageStressCheck==1){ui->averageElemenBox->setChecked(true);}
    if(materialParameterCheck==1){ui->materialParameterBox->setChecked(true);}
}

void SolutionControl::on_okButton_clicked()
{
    getUserData();
    assignData();
    emit sendSoltuonParameters(solutionParameters);
    this->close();
}
