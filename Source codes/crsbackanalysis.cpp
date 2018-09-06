#include "crsbackanalysis.h"
#include "ui_crsbackanalysis.h"

CRSBackAnalysis::CRSBackAnalysis(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::CRSBackAnalysis)
{
    ui->setupUi(this);
}

CRSBackAnalysis::~CRSBackAnalysis()
{
    delete ui;
}

void CRSBackAnalysis::getUserData()
{
    testData.resize(0,0);
    H0=ui->HeightLine->text().toDouble();
    R=ui->RadiusLine->text().toDouble();
    e0=ui->e0Line->text().toDouble();
    fileName=ui->fileNameLine->text();
    poission=ui->vLine->text().toDouble();
    outputType=ui->outputCombo->currentIndex();
    Hs=H0/(1.0f+e0);

    if(ui->linearCheck->isChecked())
    {
        theoryType=0;
    }
    else
    {
        theoryType=1;
    }

    GetFile getfile;
    getfile.fileName=fileName;
    getfile.DoGetFile();
    sucessFlag=getfile.sucessFlag;
    if(sucessFlag==true)
    {
        testData.resize(getfile.data_file.rows(),4);
        testData.col(0)=getfile.data_file.col(timeCol-1);
        testData.col(1)=getfile.data_file.col(strainCol-1);
        testData.col(2)=getfile.data_file.col(stressCol-1);
        testData.col(3)=getfile.data_file.col(poreCol-1);
    }
}

void CRSBackAnalysis::linearTheory()
{
    for (int i=1;i<(testData.rows()-1);i++) //Skip first line, and last line
    {
        int i1=i+1;
        int i0=i-1;

        double dH1, dH0, dH, strain1, strain0,strain, time1, time0, kn, cv, mv, K, gf;
        double stress1, stress0, stress;
        double effStress1, effStress0, effStress;
        double pore1, pore0, pore;
        double strainRate, Hn;
        double dt;

        gf=9.81;
        time1=testData(i1,0);
        time0=testData(i0,0);
        dt=60*(time1-time0); //convert from minutes to seconds

        strain1=testData(i1,1);
        strain0=testData(i0,1);
        strain=testData(i,1);
        dH1=strain1*0.01*H0;
        dH0=strain0*0.01*H0;
        dH=strain*0.01*H0;

        strainRate=(dH1-dH0)/(H0*dt);
        Hn=H0-dH;

        stress=testData(i,2);
        stress1=testData(i1,2);
        stress0=testData(i0,2);

        pore=testData(i,3);
        pore1=testData(i1,3);
        pore0=testData(i0,3);;

        effStress=stress-2.0f*pore/3.0f;
        effStress1=stress1-2.0f*pore1/3.0f;
        effStress0=stress0-2.0f*pore0/3.0f;

        kn=strainRate*Hn*H0*gf/(2.0f*pore);
        mv=(strain1-strain0)*0.01/(effStress1-effStress0);
        cv=kn/(mv*gf);

        K=(1.0f+poission)/(3.0f*mv*(1.0f-poission));
        resultData(i,0)=effStress;
        resultData(i,1)=K;
        resultData(i,2)=kn;

        if(outputType==1)
        {
            resultData(i,3)=mv;
            resultData(i,4)=cv;
            resultData(i,5)=(Hn-Hs)/Hs;
        }

        //first line
        if(i==1)
        {
            resultData.row(0)=resultData.row(i);
            resultData(0,0)=effStress0;

        }

        //last line
        if(i==resultData.rows()-2)
        {
            resultData.row(resultData.rows()-1)=resultData.row(i);
            resultData(resultData.rows()-1,0)=effStress1;
        }
    }
}

void CRSBackAnalysis::nonLinearTheory()
{
    qDebug()<<"Using non-linear theory"<<endl;
    for (int i=1;i<(testData.rows()-1);i++) //Skip first line, and last line
    {
        int i1=i+1;
        int i0=i-1;

        double dH1, dH0, dH, strain1, strain0,strain, time1, time0, kn, cv, mv, K, gf;
        double stress1, stress0, stress;
        double effStress1, effStress0, effStress;
        double pore1, pore0, pore;
        double strainRate, Hn;
        double dt;

        gf=9.81;
        time1=testData(i1,0);
        time0=testData(i0,0);
        dt=60*(time1-time0); //convert from minutes to seconds

        strain1=testData(i1,1);
        strain0=testData(i0,1);
        strain=testData(i,1);
        dH1=strain1*0.01*H0;
        dH0=strain0*0.01*H0;
        dH=strain*0.01*H0;

        strainRate=(dH1-dH0)/(H0*dt);
        Hn=H0-dH;

        stress=testData(i,2);
        stress1=testData(i1,2);
        stress0=testData(i0,2);

        pore=testData(i,3);
        pore1=testData(i1,3);
        pore0=testData(i0,3);;

        effStress=pow((stress*stress*stress-2.0f*stress*stress*pore+stress*pore*pore),0.3333333);
        effStress1=pow((stress1*stress1*stress1-2.0f*stress1*stress1*pore1+stress1*pore1*pore1),0.3333333);
        effStress0=pow((stress0*stress0*stress0-2.0f*stress0*stress0*pore0+stress0*pore0*pore0),0.3333333);

        cv=-H0*Hn*log10(stress1/stress0)/(2.0f*dt*log10(1.0f-pore/stress));
        kn=-0.434*strainRate*H0*Hn*gf/(2.0f*effStress*log10(1.0f-pore/stress));
        mv=kn/(gf*cv);

        K=(1.0f+poission)/(3.0f*mv*(1.0f-poission));
        resultData(i,0)=effStress;
        resultData(i,1)=K;
        resultData(i,2)=kn;

        if(outputType==1)
        {
            resultData(i,3)=mv;
            resultData(i,4)=cv;
            resultData(i,5)=(Hn-Hs)/Hs;
        }

        //first line
        if(i==1)
        {
            resultData.row(0)=resultData.row(i);
            resultData(0,0)=effStress0;
        }

        //last line
        if(i==resultData.rows()-2)
        {
            resultData.row(resultData.rows()-1)=resultData.row(i);
            resultData(resultData.rows()-1,0)=effStress1;
        }
    }
}

void CRSBackAnalysis::exportData()
{
    WriteToFile exportFile;
    folderName=QFileDialog::getExistingDirectory(Q_NULLPTR,"Choose Save folder");
    bool ok;
    QString saveName=QInputDialog::getText(Q_NULLPTR,"Set file name","Set File Name",QLineEdit::Normal,"test",&ok);
    saveName=folderName+"/"+saveName+".dat";
    exportFile.fileName=saveName.toStdString();
    exportFile.ToFile(resultData);
    QMessageBox::information(Q_NULLPTR,"Done","Done");
}

void CRSBackAnalysis::pauseSystem()
{
    do {
        cout << '\n' << "Press the Enter key to continue.";
    } while (cin.get() != '\n');
}

void CRSBackAnalysis::prepareSimulationData()
{
    getUserData();
    if(sucessFlag==true)
    {
        qDebug()<<"Using linear theory"<<endl;
        if(outputType==0)
        {
            resultData.resize(testData.rows(),3); //Effective stress, K, kv
        }
        else if (outputType==1)
        {
            resultData.resize(testData.rows(),6); //Effective stress, K, kv, mv, cv, voidRatio
        }

        if(theoryType==0)
        {
            linearTheory();
        }
        else
        {
            nonLinearTheory();
        }
    }
    else
    {
        QMessageBox::warning(Q_NULLPTR,"FAILED","Check test file");
    }

    crsObject.H0=H0;
    crsObject.R=R;
    crsObject.e0=e0;
    crsObject.possionRatio=poission;
    crsObject.analysisType=ui->simulationType->currentIndex();
    crsObject.crsFlag=true;
    crsObject.crsData.resize(resultData.rows(),6); //Time, displacement, load, pore pressure, K, k;
    crsObject.crsData.col(0)=testData.col(0);
    crsObject.crsData.col(1)=-0.01*H0*testData.col(1);
    crsObject.crsData.col(2)=testData.col(2);
    crsObject.crsData.col(3)=testData.col(3);
    crsObject.crsData.col(4)=resultData.col(1);
    crsObject.crsData.col(5)=resultData.col(2);

}

void CRSBackAnalysis::on_ASTMButton_clicked()
{
    getUserData();
    if(sucessFlag==true)    {

        qDebug()<<"Using linear theory"<<endl;
        if(outputType==0)
        {
            resultData.resize(testData.rows(),3); //Effective stress, K, kv
        }
        else if (outputType==1)
        {
            resultData.resize(testData.rows(),6); //Effective stress, K, kv, mv, cv, voidRatio
        }

        if(theoryType==0)
        {
            linearTheory();
        }
        else
        {
            nonLinearTheory();
        }
        exportData();
    }
    else
    {
        QMessageBox::warning(Q_NULLPTR,"FAILED","Check test file");
    }
}

void CRSBackAnalysis::on_closeButton_clicked()
{
    this->close();
}

void CRSBackAnalysis::on_browseLine_clicked()
{
    fileName=QFileDialog::getOpenFileName(Q_NULLPTR,"Open Test Data File");
    QStringList list;
    ui->fileNameLine->setText(fileName);

    while(list.count()!=4)
    {
        bool ok;
        QString inputCol=QInputDialog::getText(Q_NULLPTR,"Column of parameters","Column of test time, stress, stress, pore pressure:",QLineEdit::Normal,"1,2,3,4",&ok);
        list=inputCol.split(",");
        if(list.count()==4)
        {
            list=inputCol.split(",");
            break;
        }
        else
        {
            QMessageBox::warning(this,"ERROR","Need 4 columns");
        }
    }

    timeCol=1, stressCol=3, poreCol=4, strainCol=2;
    QString temp;
    temp=list[0];
    timeCol=temp.toInt();

    temp=list[1];
    strainCol=temp.toInt();

    temp=list[2];
    stressCol=temp.toInt();

    temp=list[3];
    poreCol=temp.toInt();
}

void CRSBackAnalysis::on_runCrsButton_clicked()
{
    qDebug()<<"Run simulation of CRS test"<<endl;
    prepareSimulationData();
    emit sendCrsData(crsObject,false);
}

void CRSBackAnalysis::on_runBackAnalysis_clicked()
{
    qDebug()<<"Run simulation of CRS test"<<endl;
    prepareSimulationData();
    emit sendCrsData(crsObject,true);
}
