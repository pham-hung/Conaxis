#include "boundarycondition.h"
#include "ui_boundarycondition.h"

BoundaryCondition::BoundaryCondition(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::BoundaryCondition)
{
    ui->setupUi(this);
}

BoundaryCondition::~BoundaryCondition()
{
    delete ui;
}

void BoundaryCondition::getUserData()
{
    successCheck=true;
    boundaryIndex=ui->indexLine->text().toInt();
    boundaryName=ui->nameLine->text();
    boundaryType=ui->typeCombo->currentIndex();
    loadType=ui->loadTypeCombo->currentIndex();
    val1=ui->t1Line->text().toDouble();
    val0=ui->t0Line->text().toDouble();
    if(loadType==0)
    {
        val1=val0;
        ui->t1Line->setText(QString::number(val1));
    }

    if(loadType!=2)
    {
        fileName="";
        loadCurve=MatrixXd::Zero(1,1);
    }
    else if(loadType==2)
    {
        fileName=ui->curveFileNameLine->text();
        loadCurve=loadCurve;
        val1=0;
        val0=0;
        ui->t1Line->setText("");
        ui->t0Line->setText("");
    }
}

void BoundaryCondition::getDatafromVector(int i)
{
    boundaryName=boundary[i].boundaryName;
    boundaryType=boundary[i].boundaryType;
    loadType=boundary[i].loadType;
    val0=boundary[i].val0;
    val1=boundary[i].val1;
    fileName=boundary[i].fileName;
    loadCurve=boundary[i].loadCurve;
    ui->nameLine->setText(boundaryName);
    ui->typeCombo->setCurrentIndex(boundaryType);
    ui->loadTypeCombo->setCurrentIndex(loadType);
    ui->t0Line->setText(QString::number(val0));
    ui->t1Line->setText(QString::number(val1));
    ui->curveFileNameLine->setText(fileName);
}

void BoundaryCondition::showData()
{
    cout<<"----------------------"<<endl;
    cout<<"Boundary condition type : 1-X-disp, 2-Y-disp, 3-Pressure, 4-Line PressureX, 5-Line Pressure Y, 6-Point Load"<<endl;
    cout<<"Load type : 1-Constant Load, 2- Ramp Load"<<endl;
    for (auto i=0;i<boundary.size();i++)
    {
        cout<<"----------------------"<<endl;
        cout<<"Boundary condition number: "<<boundary[i].boundaryIndex<<endl;
        cout<<"Boundary condition name  : "<<boundary[i].boundaryName.toStdString()<<endl;
        cout<<"Boundary condition type  : "<<boundary[i].boundaryType<<endl;
        cout<<"Load   type              : "<<boundary[i].loadType<<endl;
        cout<<"Constant or Start Value  : "<<boundary[i].val0<<endl;
        cout<<"End Value                : "<<boundary[i].val1<<endl;
        cout<<"Data file                : "<<boundary[i].fileName.toStdString()<<endl;
        cout<<boundary[i].loadCurve<<endl;
    }
}

void BoundaryCondition::defaultData()
{
    ui->nameLine->setText("");
    ui->typeCombo->setCurrentIndex(-1);
    ui->loadTypeCombo->setCurrentIndex(-1);
    ui->t0Line->setText("");
    ui->t1Line->setText("");
    ui->curveFileNameLine->setText("");
    loadCurve.resize(1,1);
    fileName="";
}

void BoundaryCondition::sendSignal()
{
    emit sendParameters(boundary);
}

void BoundaryCondition::createCrsBoundary(Ref<MatrixXd> testData)
{
    boundary.resize(5);
    //fixx, fixy, pore pressure, moveY, applied load
    int j=0;
    boundary[j].boundaryIndex=j+1;
    boundary[j].boundaryName="Fix X";
    boundary[j].boundaryType=0;
    boundary[j].loadType=0;
    boundary[j].val0=0;
    boundary[j].val1=0;
    boundary[j].fileName="";
    boundary[j].loadCurve.resize(1,1);
    boundary[j].loadCurve<<0;

    j=1;
    boundary[j].boundaryIndex=j+1;
    boundary[j].boundaryName="Fix Y";
    boundary[j].boundaryType=1;
    boundary[j].loadType=0;
    boundary[j].val0=0;
    boundary[j].val1=0;
    boundary[j].fileName="";
    boundary[j].loadCurve.resize(1,1);
    boundary[j].loadCurve<<0;

    j=2;
    boundary[j].boundaryIndex=j+1;
    boundary[j].boundaryName="Fix pore pressure";
    boundary[j].boundaryType=2;
    boundary[j].loadType=0;
    boundary[j].val0=0;
    boundary[j].val1=0;
    boundary[j].fileName="";
    boundary[j].loadCurve.resize(1,1);
    boundary[j].loadCurve<<0;

    j=3;
    boundary[j].boundaryIndex=j+1;
    boundary[j].boundaryName="Move Y";
    boundary[j].boundaryType=1;
    boundary[j].loadType=2;
    boundary[j].val0=0;
    boundary[j].val1=0;
    boundary[j].fileName="From CRS data";
    boundary[j].loadCurve.resize(testData.rows(),2);
    boundary[j].loadCurve.col(0)=testData.col(0);
    boundary[j].loadCurve.col(1)=testData.col(1);

    j=4;
    boundary[j].boundaryIndex=j+1;
    boundary[j].boundaryName="Applied Load";
    boundary[j].boundaryType=3;
    boundary[j].loadType=2;
    boundary[j].val0=0;
    boundary[j].val1=0;
    boundary[j].fileName="From CRS data";
    boundary[j].loadCurve.resize(testData.rows(),2);
    boundary[j].loadCurve.col(0)=testData.col(0);
    boundary[j].loadCurve.col(1)=-testData.col(2);
    emit sendParameters(boundary);
}

void BoundaryCondition::on_closeButton_clicked()
{
    this->close();
}

void BoundaryCondition::on_typeCombo_activated(int index)
{

}

void BoundaryCondition::on_loadTypeCombo_activated(int index)
{
    if(index==0)
    {
       ui->t1Line->setDisabled(true);       
    }
    else if(index==1)
    {
      ui->t1Line->setDisabled(false);
    }
    else if(index==2)
    {
        fileName=QFileDialog::getOpenFileName();
        getfile.fileName=fileName;
        getfile.DoGetFile();
        if(getfile.data_file.rows()!=1)
        {
            bool ok;
            QString timeCol=QInputDialog::getText(this,"Input Value","Colum of time",QLineEdit::Normal,"1",&ok);
            QString valCol=QInputDialog::getText(this,"Input Value","Colum of value",QLineEdit::Normal,"2",&ok);

            QString factor=QInputDialog::getText(this,"Input Value","Scale value with factor",QLineEdit::Normal,"1",&ok);
            QString timeFactor=QInputDialog::getText(this,"Input Value","Scale Time with factor",QLineEdit::Normal,"6.944444e-4",&ok);

            if(timeCol.toInt()>getfile.data_file.cols()||valCol.toInt()>getfile.data_file.cols())
            {
              QMessageBox::warning(this,"ERROR","Input colum is larger than data file");
              return;
            }

            loadCurve.resize(getfile.data_file.rows(),2);
            loadCurve.col(0)=timeFactor.toDouble()*getfile.data_file.col(timeCol.toInt()-1);
            loadCurve.col(1)=factor.toDouble()*getfile.data_file.col(valCol.toInt()-1);
            ui->curveFileNameLine->setText(fileName);
        }
    }
}

void BoundaryCondition::on_indexLine_cursorPositionChanged(int arg1, int arg2)
{
    boundaryIndex=ui->indexLine->text().toInt();
    if(boundaryIndex<1)
    {
        boundaryIndex=1;
        ui->indexLine->setText("1");
    }
    if(boundaryIndex>(boundary.size()+1))
    {
        boundaryIndex=boundary.size()+1;
        ui->indexLine->setText(QString::number(boundaryIndex));
    }

    if(boundaryIndex<=boundary.size())
    {
        int j=boundaryIndex-1;
        getDatafromVector(j);
    }
    else
    {
        defaultData();
    }
}

void BoundaryCondition::on_listButton_clicked()
{
    showData();
}

void BoundaryCondition::on_okButton_clicked()
{

    getUserData();
    boundaryIndex=ui->indexLine->text().toInt();
    int j=boundaryIndex-1;
    if(boundaryIndex<=boundary.size())
    {
        boundary[j].boundaryIndex=boundaryIndex;
        boundary[j].boundaryName=boundaryName;
        boundary[j].boundaryType=boundaryType;
        boundary[j].loadType=loadType;
        boundary[j].val0=val0;
        boundary[j].val1=val1;
        boundary[j].fileName=fileName;
        boundary[j].loadCurve=loadCurve;
        cout<<"Boundry condition "<<boundaryIndex<<" is modified"<<endl;
    }
    else
    {
        BoundaryConditionBase newBase;
        newBase.boundaryIndex=boundaryIndex;
        newBase.boundaryName=boundaryName;
        newBase.boundaryType=boundaryType;
        newBase.loadType=loadType;
        newBase.val0=val0;
        newBase.val1=val1;
        newBase.fileName=fileName;
        newBase.loadCurve=loadCurve;
        boundary.push_back(newBase);
        cout<<"Boundry condition "<<boundaryIndex<<" is created"<<endl;
    }
}

void BoundaryCondition::getBoundary(vector<BoundaryConditionBase> boundary)
{
    this->boundary=boundary;
}

void BoundaryCondition::on_indexLine_textChanged(const QString &arg1)
{
    boundaryIndex=ui->indexLine->text().toInt();
    if(boundaryIndex<1)
    {
        boundaryIndex=1;
        ui->indexLine->setText("1");
    }
    if(boundaryIndex>(boundary.size()+1))
    {
        boundaryIndex=boundary.size()+1;
        ui->indexLine->setText(QString::number(boundaryIndex));
    }

    if(boundaryIndex<=boundary.size())
    {
        int j=boundaryIndex-1;
        getDatafromVector(j);
    }
    else
    {
        defaultData();
    }
}

