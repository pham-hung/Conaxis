#include "assignboundarycondition.h"
#include "ui_assignboundarycondition.h"

AssignBoundaryCondition::AssignBoundaryCondition(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::AssignBoundaryCondition)
{
    ui->setupUi(this);
}

AssignBoundaryCondition::~AssignBoundaryCondition()
{
    delete ui;
}

void AssignBoundaryCondition::resetData()
{    
    ui->stageCombo->clear();
    ui->boundaryCombo->clear();

    for(auto i=0;i<stage.size();i++)
    {
        ui->stageCombo->addItem(QString::number(stage[i].stageIndex));
    }
    for(auto i=0;i<boundary.size();i++)
    {
        ui->boundaryCombo->addItem(QString::number(boundary[i].boundaryIndex));
    }
}

void AssignBoundaryCondition::getUserData()
{
    stageIndex=ui->stageCombo->currentText().toInt();
    stageName=ui->stageNameLine->text();
    stageBoundaryIndex=ui->stageBoundaryIndexLine->text().toInt();
    boundaryIndex=ui->boundaryCombo->currentText().toInt();
    boundaryName=ui->boundaryNameLine->text();

    assignType=ui->assignTypeCombo->currentIndex();
    if(assignType<0){assignType=0;}
    x0=ui->x0Line->text().toDouble();
    x1=ui->x1Line->text().toDouble();
    y0=ui->y0Line->text().toDouble();
    y1=ui->y1Line->text().toDouble();
    fileName=ui->listFileLine->text();
    if(assignType==0)
    {
        fileName="";
        nodeListTemp.resize(1,1);
        nodeListTemp(0,0)=0;
    }

    if(ui->gradientCheck->isChecked())
    {
        gradientBool=true;
    }
    else
    {
        gradientBool=false;
    }
    aFactor=ui->aFactorLine->text().toDouble();
    bFactor=ui->bFactorLine->text().toDouble();
    cFactor=ui->cFactorLine->text().toDouble();

}

void AssignBoundaryCondition::assignUI()
{
    ui->stageNameLine->setText(stageName);
    ui->assignTypeCombo->setCurrentIndex(assignType);
    ui->boundaryNameLine->setText(boundaryName);
    ui->boundaryCombo->setCurrentText(QString::number(boundaryIndex,'f',0));
    ui->x0Line->setText(QString::number(x0,'f',5));
    ui->x1Line->setText(QString::number(x1,'f',5));
    ui->y0Line->setText(QString::number(y0,'f',5));
    ui->y1Line->setText(QString::number(y1,'f',5));
    ui->aFactorLine->setText(QString::number(aFactor,'f',5));
    ui->bFactorLine->setText(QString::number(bFactor,'f',5));
    ui->cFactorLine->setText(QString::number(cFactor,'f',5));
    ui->listFileLine->setText(fileName);
    if(gradientBool==true)
    {
        ui->gradientCheck->setChecked(true);
    }
    else
    {
        ui->gradientCheck->setChecked(false);
    }
}

void AssignBoundaryCondition::defaultData()
{
    stageBoundaryIndex=1;
    boundaryIndex=-1;
    boundaryName="";
    assignType=0;
    x0=0;
    x1=0;
    y0=0;
    y1=0;
    aFactor=0;
    bFactor=0;
    cFactor=0;
    gradientBool=false;
    fileName="";
}

void AssignBoundaryCondition::sendSIGNAL()
{
    emit sendParametersStageBoundary(stageBoundary);
}

void AssignBoundaryCondition::getDatafromVector(int i, int j)
{
    //i is stageIndex
    //j is stageBoundaryIndex
    stageIndex=i+1;
    stageName=stage[i].stageName;
    stageBoundaryIndex=j+1;
    boundaryIndex=stageBoundary[i].v_boundaryIndex[j];
    boundaryName=boundary[boundaryIndex-1].boundaryName;
    assignType=stageBoundary[i].v_assignType[j];
    x0=stageBoundary[i].v_x0[j];
    x1=stageBoundary[i].v_x1[j];
    y0=stageBoundary[i].v_y0[j];
    y1=stageBoundary[i].v_y1[j];
    gradientBool=stageBoundary[i].v_gradientBool[j];
    aFactor=stageBoundary[i].v_aFactor[j];
    bFactor=stageBoundary[i].v_bFactor[j];
    cFactor=stageBoundary[i].v_cFactor[j];
    fileName=stageBoundary[i].v_nodeFileName[j];
}

void AssignBoundaryCondition::createCrsStageBoundary(double H, double R, int controlType)
{
    //initial
    stageBoundary[0].v_stageBoundaryIndex.resize(0);
    stageBoundary[0].v_stageBoundaryIndex.push_back(1); //fix x
    stageBoundary[0].v_stageBoundaryIndex.push_back(2); //fix x
    stageBoundary[0].v_stageBoundaryIndex.push_back(3); //fix y
    stageBoundary[0].v_stageBoundaryIndex.push_back(4); //move y

    nodeListTemp=MatrixXd::Zero(0,0);
    stageBoundary[0].v_nodeList.resize(0);
    stageBoundary[0].v_nodeList.push_back(nodeListTemp); //fix x
    stageBoundary[0].v_nodeList.push_back(nodeListTemp); //fix x
    stageBoundary[0].v_nodeList.push_back(nodeListTemp); //fix y
    stageBoundary[0].v_nodeList.push_back(nodeListTemp); //move y

    fileName="";
    stageBoundary[0].v_nodeFileName.resize(0);
    stageBoundary[0].v_nodeFileName.push_back(fileName); //fix x
    stageBoundary[0].v_nodeFileName.push_back(fileName); //fix x
    stageBoundary[0].v_nodeFileName.push_back(fileName); //fix y
    stageBoundary[0].v_nodeFileName.push_back(fileName); //move y

    stageBoundary[0].v_assignType.resize(0);
    stageBoundary[0].v_assignType.push_back(0); //fix x
    stageBoundary[0].v_assignType.push_back(0); //fix x
    stageBoundary[0].v_assignType.push_back(0); //fix y
    stageBoundary[0].v_assignType.push_back(0); //move y

    stageBoundary[0].v_boundaryIndex.resize(0);
    if(controlType==0) //load control
    {
        stageBoundary[0].v_boundaryIndex.push_back(1); //fix x
        stageBoundary[0].v_boundaryIndex.push_back(1); //fix x
        stageBoundary[0].v_boundaryIndex.push_back(2); //fix y
        stageBoundary[0].v_boundaryIndex.push_back(5); //applied load
    }
    else
    {
        stageBoundary[0].v_boundaryIndex.push_back(1); //fix x
        stageBoundary[0].v_boundaryIndex.push_back(1); //fix x
        stageBoundary[0].v_boundaryIndex.push_back(2); //fix y
        stageBoundary[0].v_boundaryIndex.push_back(4); //move y
    }

    stageBoundary[0].v_x0.resize(0);
    stageBoundary[0].v_x0.push_back(0); //fix x
    stageBoundary[0].v_x0.push_back(R); //fix x
    stageBoundary[0].v_x0.push_back(0); //fix y
    stageBoundary[0].v_x0.push_back(0); //move y

    stageBoundary[0].v_x1.resize(0);
    stageBoundary[0].v_x1.push_back(0); //fix x
    stageBoundary[0].v_x1.push_back(R); //fix x
    stageBoundary[0].v_x1.push_back(R); //fix y
    stageBoundary[0].v_x1.push_back(R); //move y

    stageBoundary[0].v_y0.resize(0);
    stageBoundary[0].v_y0.push_back(0); //fix x
    stageBoundary[0].v_y0.push_back(0); //fix x
    stageBoundary[0].v_y0.push_back(0); //fix y
    stageBoundary[0].v_y0.push_back(H); //move y

    stageBoundary[0].v_y1.resize(0);
    stageBoundary[0].v_y1.push_back(H); //fix x
    stageBoundary[0].v_y1.push_back(H); //fix x
    stageBoundary[0].v_y1.push_back(0); //fix y
    stageBoundary[0].v_y1.push_back(H); //move y

    stageBoundary[0].v_gradientBool.resize(0);
    stageBoundary[0].v_gradientBool.push_back(false); //fix x
    stageBoundary[0].v_gradientBool.push_back(false); //fix x
    stageBoundary[0].v_gradientBool.push_back(false); //fix y
    stageBoundary[0].v_gradientBool.push_back(false); //move y

    stageBoundary[0].v_aFactor.resize(0);
    stageBoundary[0].v_aFactor.push_back(0); //fix x
    stageBoundary[0].v_aFactor.push_back(0); //fix x
    stageBoundary[0].v_aFactor.push_back(0); //fix y
    stageBoundary[0].v_aFactor.push_back(0); //move y

    stageBoundary[0].v_bFactor.resize(0);
    stageBoundary[0].v_bFactor.push_back(0); //fix x
    stageBoundary[0].v_bFactor.push_back(0); //fix x
    stageBoundary[0].v_bFactor.push_back(0); //fix y
    stageBoundary[0].v_bFactor.push_back(0); //move y


    stageBoundary[0].v_cFactor.resize(0);
    stageBoundary[0].v_cFactor.push_back(0); //fix x
    stageBoundary[0].v_cFactor.push_back(0); //fix x
    stageBoundary[0].v_cFactor.push_back(0); //fix y
    stageBoundary[0].v_cFactor.push_back(0); //move y

    stageBoundary[1]=stageBoundary[0];
    stageBoundary[1].v_stageBoundaryIndex.push_back(5); //Pore pressure
    stageBoundary[1].v_nodeList.push_back(nodeListTemp); //Pore pressure
    stageBoundary[1].v_nodeFileName.push_back(fileName); //Pore pressure
    stageBoundary[1].v_assignType.push_back(0); //Pore pressure
    stageBoundary[1].v_boundaryIndex.push_back(3); //Pore pressure
    stageBoundary[1].v_x0.push_back(0); //Pore pressure
    stageBoundary[1].v_x1.push_back(R); //Pore pressure
    stageBoundary[1].v_y0.push_back(H); //Pore pressure
    stageBoundary[1].v_y1.push_back(H); //Pore pressure
}

void AssignBoundaryCondition::getParametersStage(vector<StageBase> stage)
{
    this->stage=stage;
    stageBoundary.resize(stage.size());
    if(stage.size()<1)
    {
        QMessageBox::warning(this,"ERROR","Define stage analysis first");
        checkSuccess=false;
    }

}

void AssignBoundaryCondition::getParametersBoundary(vector<BoundaryConditionBase> boundary)
{
    this->boundary=boundary;
    if(boundary.size()<1)
    {
        QMessageBox::warning(this,"ERROR","Define boundary condition first");
        checkSuccess=false;
    }

}

void AssignBoundaryCondition::getStageBoundary(vector<StageBoundaryBase> stageBoundary)
{
    this->stageBoundary=stageBoundary;
}

void AssignBoundaryCondition::on_closeButton_clicked()
{
    this->close();
}

void AssignBoundaryCondition::on_addButton_clicked()
{
    getUserData();
    stageIndex=ui->stageCombo->currentText().toInt();
    int i=stageIndex-1;

    stageBoundaryIndex=ui->stageBoundaryIndexLine->text().toInt();
    int j=stageBoundaryIndex-1;

    if(j<stageBoundary[i].v_boundaryIndex.size())
    {
        stageBoundary[i].v_stageBoundaryIndex[j]=stageBoundaryIndex;
        stageBoundary[i].v_boundaryIndex[j]=boundaryIndex;
        stageBoundary[i].v_assignType[j]=assignType;
        stageBoundary[i].v_x0[j]=x0;
        stageBoundary[i].v_x1[j]=x1;
        stageBoundary[i].v_y0[j]=y0;
        stageBoundary[i].v_y1[j]=y1;
        stageBoundary[i].v_gradientBool[j]=gradientBool;
        stageBoundary[i].v_aFactor[j]=aFactor;
        stageBoundary[i].v_bFactor[j]=bFactor;
        stageBoundary[i].v_cFactor[j]=cFactor;
        stageBoundary[i].v_nodeFileName[j]=fileName;
        stageBoundary[i].v_nodeList[j]=nodeListTemp;
        cout<<"Boundary for stage"<<stageName.toStdString()<<" and boundary number "<<stageBoundaryIndex<<" is modified"<<endl;
    }
    else
    {
        stageBoundary[i].v_stageBoundaryIndex.push_back(stageBoundaryIndex);
        stageBoundary[i].v_boundaryIndex.push_back(boundaryIndex);
        stageBoundary[i].v_assignType.push_back(assignType);
        stageBoundary[i].v_x0.push_back(x0);
        stageBoundary[i].v_x1.push_back(x1);
        stageBoundary[i].v_y0.push_back(y0);
        stageBoundary[i].v_y1.push_back(y1);
        stageBoundary[i].v_gradientBool.push_back(gradientBool);
        stageBoundary[i].v_aFactor.push_back(aFactor);
        stageBoundary[i].v_bFactor.push_back(bFactor);
        stageBoundary[i].v_cFactor.push_back(cFactor);
        stageBoundary[i].v_nodeFileName.push_back(fileName);
        stageBoundary[i].v_nodeList.push_back(nodeListTemp);
        cout<<"Boundary for stage"<<stageName.toStdString()<<" and boundary number "<<stageBoundaryIndex<<" is added"<<endl;
    }

}

void AssignBoundaryCondition::on_assignTypeCombo_activated(int index)
{
    if(index==1)
    {
        GetFile getfile;
        fileName=QFileDialog::getOpenFileName();
        ui->listFileLine->setText(fileName);
        getfile.fileName=fileName;
        getfile.DoGetFile();
        if(getfile.checkFinish==true)
        {
            nodeListTemp.resize(getfile.data_file.rows(),1);
            bool ok;
            QString nodeCol=QInputDialog::getText(this,"Input Value","Colum of node list",QLineEdit::Normal,"1",&ok);
            if(ok)
            {
                nodeListTemp.resize(getfile.data_file.rows(),1);
                nodeListTemp.col(0)=getfile.data_file.col(nodeCol.toInt()-1);
                cout<<nodeListTemp<<endl;
            }
            else {return;}
        }
        else
        {
            return;
        }
    }
}

void AssignBoundaryCondition::on_stageCombo_activated(int index)
{
    bool ok=false;
    if(index>0&&stageBoundary[index].v_boundaryIndex.size()==0)
    {
        QString lastStage=QString::number(index);
        QString copyStage=QInputDialog::getText(this,"Copy Boundary","Copy From Stage: ",QLineEdit::Normal,lastStage,&ok);
        if(ok)
        {
            stageBoundary[index]=stageBoundary[index-1];
        }
    }
    else
    {
        stageIndex=index+1;
        int i=stageIndex-1;
        stageName=stage[i].stageName;
        ui->stageNameLine->setText(stageName);
        ui->boundaryCombo->setCurrentIndex(-1);
        ui->boundaryNameLine->setText("");
        ui->stageBoundaryIndexLine->setText("");
        ui->assignTypeCombo->setCurrentIndex(0);
        ui->x0Line->setText("");
        ui->x1Line->setText("");
        ui->y0Line->setText("");
        ui->y1Line->setText("");
        ui->aFactorLine->setText("");
        ui->bFactorLine->setText("");
        ui->cFactorLine->setText("");
        ui->gradientCheck->setChecked(false);
    }
}

void AssignBoundaryCondition::on_boundaryCombo_activated(int index)
{
    boundaryIndex=index+1;
    int i=boundaryIndex-1;
    ui->boundaryNameLine->setText(boundary[i].boundaryName);
}

void AssignBoundaryCondition::on_stageBoundaryIndexLine_textChanged(const QString &arg1)
{

    stageIndex=ui->stageCombo->currentText().toInt();
    int i=stageIndex-1;

    stageBoundaryIndex=ui->stageBoundaryIndexLine->text().toInt();
    int j=stageBoundaryIndex-1;

    if(stageBoundaryIndex<1)
    {
        stageBoundaryIndex=1;
    }

    if(stageBoundaryIndex>stageBoundary[i].v_boundaryIndex.size())
    {
        stageBoundaryIndex=stageBoundary[i].v_boundaryIndex.size()+1;
    }
    ui->stageBoundaryIndexLine->setText(QString::number(stageBoundaryIndex));

    if(j<stageBoundary[i].v_boundaryIndex.size())
    {
        boundaryIndex=stageBoundary[i].v_boundaryIndex[j];
        ui->boundaryCombo->setCurrentIndex(boundaryIndex-1);
        ui->boundaryNameLine->setText(boundary[boundaryIndex-1].boundaryName);
        ui->assignTypeCombo->setCurrentIndex(stageBoundary[i].v_assignType[j]);
        ui->x0Line->setText(QString::number(stageBoundary[i].v_x0[j]));
        ui->x1Line->setText(QString::number(stageBoundary[i].v_x1[j]));
        ui->y0Line->setText(QString::number(stageBoundary[i].v_y0[j]));
        ui->y1Line->setText(QString::number(stageBoundary[i].v_y1[j]));
        ui->aFactorLine->setText(QString::number(stageBoundary[i].v_aFactor[j]));
        ui->bFactorLine->setText(QString::number(stageBoundary[i].v_bFactor[j]));
        ui->cFactorLine->setText(QString::number(stageBoundary[i].v_cFactor[j]));
        if(stageBoundary[i].v_gradientBool[j]==true)
        {
            ui->gradientCheck->setChecked(true);
        }
        else
        {
            ui->gradientCheck->setChecked(false);
        }
        ui->listFileLine->setText(stageBoundary[i].v_nodeFileName[j]);
    }
    else
    {
        defaultData();
        assignUI();
    }
}

void AssignBoundaryCondition::on_listButton_clicked()
{
    for (int i=0;i<stageBoundary.size();i++)
    {

        for (int j=0;j<stageBoundary[i].v_stageBoundaryIndex.size();j++)
        {
            getDatafromVector(i,j);
            cout<<"-----------------------------------"<<endl;
            cout<<"Analysis: "<<stageName.toStdString()<<" Boundary condition: "<<boundaryName.toStdString()<<endl;
            if(assignType==0)
            {
                cout<<"Assiged to node which has X, Y Coordinates "<<endl;
                cout<<"X =("<<x0<<" : "<<x1<<" )"<<endl;
                cout<<"Y =("<<y0<<" : "<<y1<<" )"<<endl;
            }
            else
            {
                cout<<"Node from file: "<<fileName.toStdString()<<endl;
            }
        }
    }
}

void AssignBoundaryCondition::on_resetButton_clicked()
{
    stageIndex=ui->stageCombo->currentText().toInt();
    int i=stageIndex-1;
    StageBoundaryBase newBoundary;
    stageBoundary[i]=newBoundary;
}
