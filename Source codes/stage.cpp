#include "stage.h"
#include "ui_stage.h"

Stage::Stage(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::Stage)
{
    ui->setupUi(this);
}

Stage::~Stage()
{
    delete ui;
}

void Stage::getUserData()
{
    successCheck=true;
    stageIndex=ui->stageIndexLine->text().toInt();
    stageName=ui->nameLine->text();
    stageType=ui->typeCombo->currentIndex();

    if(stageType<0)
    {
        stageType=2;
        ui->typeCombo->setCurrentIndex(stageType);

    }

    if(stageType==0||stageType==1)
    {
        subStep=1;
        ui->subLine->setText("1");
        timeStepType=0;
        ui->timeStepCombo->setCurrentIndex(timeStepType);
        ui->fileNameLine->setText("");
        timeStep.resize(1,1);
        timeStep.setZero();
    }
    else
    {
        subStep=ui->subLine->text().toInt();
    }

    if(stageType==0&&stageIndex!=1)
    {
        QMessageBox::warning(this,"ERROR","In-Situ analysis must be first step");
        successCheck=false;
        return;
    }
    dt=ui->dtLine->text().toDouble();
    t0=ui->t0Line->text().toDouble();
    t1=dt+t0;
    ui->t1Line->setText(QString::number(t1));

    if(t1<t0)
    {
        QMessageBox::warning(this,"ERROR","time t1 must be greater than t0");
        successCheck=false;
        return;
    }
    if(subStep<1)
    {
        QMessageBox::warning(this,"ERROR","Sub-Step must be > 0");
        successCheck=false;
        return;
    }
    if(stageIndex>1)
    {
        if(t0!=stage[stageIndex-2].t1)
        {
            t0=stage[stageIndex-2].t1;
            ui->t0Line->setText(QString::number(t0));
            t1=dt+t0;
            ui->t1Line->setText(QString::number(t1));
        }
    }

    //Gravity Check
    if(ui->gravityCheck->isChecked())
    {
        gravityLoad=1;
    }
    else
    {
        gravityLoad=0;
    }

    timeStepType=ui->timeStepCombo->currentIndex();

    if(timeStepType==1)
    {
        fileName=ui->fileNameLine->text();
        t1=timeStep(timeStep.rows()-1,0)+t0;
        ui->t1Line->setText(QString::number(t1));
        dt=t1-t0;
        ui->dtLine->setText(QString::number(dt));
        subStep=timeStep.rows()-1;
        ui->subLine->setText(QString::number(subStep));
    }
    if(timeStepType==0)
    {
        fileName="";
        timeStep.resize(1,1);
        timeStep(0,0)=0;
    }

}

void Stage::getDataFromVector(int i)
{
    stageIndex=stage[i].stageIndex;
    stageName=stage[i].stageName;
    ui->nameLine->setText(stageName);
    stageType=stage[i].stageType;
    if(stageType==0||stageType==1)
    {
        subStep=1;
        ui->subLine->setText("1");
    }
    else
    {
        subStep=stage[i].subStep;
    }
    ui->subLine->setText(QString::number(subStep));
    ui->typeCombo->setCurrentIndex(stageType);
    dt=stage[i].dt;
    ui->dtLine->setText(QString::number(dt));
    t0=stage[i].t0;
    ui->t0Line->setText(QString::number(t0));
    t1=stage[i].t1;
    ui->t1Line->setText(QString::number(t1));

    //Gravity check
    gravityLoad=stage[i].gravityLoad;
    if(gravityLoad==1)
    {
        ui->gravityCheck->setChecked(true);
    }
    else
    {
        ui->gravityCheck->setChecked(false);
    }

    timeStep=stage[i].timeStep;
    timeStepType=stage[i].timeStepType;
    ui->timeStepCombo->setCurrentIndex(timeStepType);
    ui->fileNameLine->setText(stage[i].fileName);
}

void Stage::showData()
{
    if(stage.size()==0)
    {
        cout<<"Stage is not defined"<<endl;
    }
    for (auto i=0;i<stage.size();i++)
    {
        cout<<"----------------------------"<<endl;
        cout<<"Stage number      : "<<stage[i].stageIndex<<endl;
        cout<<"stage name        : "<<stage[i].stageName.toStdString()<<endl;
        switch (stage[i].stageType) {
        case 0:
            cout<<"stage Type        : In-situ analysis"<<endl;
            break;
        case 1:
            cout<<"stage Type        : Undrained analysis"<<endl;
            break;
        case 2:
            cout<<"stage Type        : Consolidation analysis"<<endl;
            break;
        default:
            break;
        }

        cout<<"Time of analysis  : "<<stage[i].dt<<endl;
        cout<<"Start day         : "<<stage[i].t0<<" End day: "<<stage[i].t1<<endl;
        cout<<"Number of sub-step: "<<stage[i].subStep<<endl;
        cout<<"Number of sub-step: "<<stage[i].timeStep<<endl;
        cout<<"Time Step Type    : "<<stage[i].timeStepType<<endl;
    }
}

void Stage::defaultData()
{
    stageIndex=ui->stageIndexLine->text().toInt();
    ui->nameLine->setText("");
    ui->typeCombo->setCurrentIndex(-1);
    ui->dtLine->setText("");

    if(stageIndex>1)
    {
        t1=stage[stageIndex-2].t1;
        ui->t0Line->setText(QString::number(t1));
    }
    else
    {
        ui->t0Line->setText("");
    }

    ui->t1Line->setText("");
    ui->subLine->setText("");
    ui->gravityCheck->setChecked(false);
    ui->fileNameLine->setText("");
    ui->timeStepCombo->setCurrentIndex(0);
}

void Stage::updateTime(int i)
{
    for (auto j=i;j<stage.size();j++)
    {
        stage[j].t0=stage[j-1].t1;
        stage[j].t1=stage[j].t0+stage[j].dt;
    }
}

void Stage::sendSignal()
{
    emit sendParameters(stage);
}

void Stage::createCrsStage(Ref<MatrixXd> testTime)
{
    //get information from testTime
    testTime.col(0)=testTime.col(0)/1440.0f;
    double testingTime=testTime.col(0).maxCoeff();
    int numberOfStep=testTime.rows()-1;
    stage.resize(2);
    //initial stage
    int index=0;
    stage[index].stageIndex=index+1;
    stage[index].stageName="Undrained stage";
    stage[index].stageType=1;
    stage[index].dt=0;
    stage[index].t0=0;
    stage[index].t1=0;
    stage[index].subStep=1;
    stage[index].gravityLoad=0;
    stage[index].timeStepType=0;
    stage[index].timeStep.resize(1,1);
    stage[index].fileName="";

    //consolidation stage
    index=1;
    stage[index].stageIndex=index+1;
    stage[index].stageName="Testing stage";
    stage[index].stageType=2;
    stage[index].dt=0;
    stage[index].t0=0;
    stage[index].t1=testingTime;
    stage[index].subStep=numberOfStep;
    stage[index].gravityLoad=0;
    stage[index].timeStepType=1;
    stage[index].timeStep=testTime;
    stage[index].fileName="From CRS Test";
    emit sendParameters(stage);
}

void Stage::on_cancelButton_clicked()
{
    this->close();
}

void Stage::on_listButton_clicked()
{
    showData();
}

void Stage::on_addButton_clicked()
{
    stageIndex=ui->stageIndexLine->text().toInt();
    ui->stageIndexLine->setText(QString::number(stageIndex));
    getUserData();
    if(successCheck==false)
    {
        return;
    }
    else
    {
        int index=stageIndex-1;
        if(stageIndex<=stage.size())
        {
            stage[index].stageIndex=stageIndex;
            stage[index].stageName=stageName;
            stage[index].stageType=stageType;
            stage[index].dt=dt;
            stage[index].t0=t0;
            stage[index].t1=t1;
            stage[index].subStep=subStep;
            stage[index].gravityLoad=gravityLoad;
            stage[index].timeStepType=timeStepType;
            stage[index].timeStep=timeStep;
            stage[index].fileName=fileName;
            updateTime(stageIndex);
            cout<<"Stage: "<<stageIndex<<" is modified"<<endl;
        }
        else
        {
            StageBase newStage;
            newStage.stageIndex=stageIndex;
            newStage.stageName=stageName;
            newStage.stageType=stageType;
            newStage.dt=dt;
            newStage.t0=t0;
            newStage.t1=t1;
            newStage.subStep=subStep;
            newStage.gravityLoad=gravityLoad;
            newStage.timeStep=timeStep;
            newStage.timeStepType=timeStepType;
            newStage.fileName=fileName;
            stage.push_back(newStage);
            cout<<"Stage: "<<stageIndex<<" is created"<<endl;
        }
    }

}

void Stage::on_stageIndexLine_textChanged(const QString &arg1)
{
    stageIndex=ui->stageIndexLine->text().toInt();
    if(stageIndex<1)
    {
        stageIndex=1;
        ui->stageIndexLine->setText(QString::number(stageIndex));
    }
    if(stageIndex>stage.size())
    {
        stageIndex=stage.size()+1;
        ui->stageIndexLine->setText(QString::number(stageIndex));
    }

    int index=stageIndex-1;
    if(stageIndex<=stage.size())
    {
        //Modify material
        getDataFromVector(index);
    }
    else
    {
        defaultData();
    }

}

void Stage::on_typeCombo_activated(int index)
{
    if(index==0||index==1)
    {
        ui->subLine->setText("1");
        ui->dtLine->setText("0");
    }
}

void Stage::on_dtLine_textChanged(const QString &arg1)
{
    dt=ui->dtLine->text().toDouble();
    t0=ui->t0Line->text().toDouble();
    t1=dt+t0;
    ui->t1Line->setText(QString::number(t1));
}

void Stage::getStage(vector<StageBase> stage)
{
    this->stage=stage;
}

void Stage::on_timeStepCombo_activated(int index)
{
    if(index==1)
    {
        fileName=QFileDialog::getOpenFileName(Q_NULLPTR,"Choose time step file");
        ui->fileNameLine->setText(fileName);
        bool ok;
        QString timeCol=QInputDialog::getText(Q_NULLPTR,"Enter Value","Colume of time step: =",QLineEdit::Normal,"1",&ok);
        QString factor=QInputDialog::getText(Q_NULLPTR,"Enter Value","Scale value with factor (for convert unit): =",QLineEdit::Normal,"6.944444e-4",&ok);

        GetFile getFile;
        getFile.fileName=fileName;
        getFile.DoGetFile();
        if(timeCol.toInt()>getFile.data_file.cols())
        {
            QMessageBox::warning(Q_NULLPTR,"ERROR","Input value is too large");
            return;
        }
        else
        {
            timeStep.resize(getFile.data_file.rows(),1);
            timeStep.col(0)=factor.toDouble()*getFile.data_file.col(timeCol.toInt()-1);
            ui->subLine->setText(QString::number(timeStep.rows()-1));
        }
    }
}
