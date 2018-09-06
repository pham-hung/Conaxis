#include "terzaghi1d.h"
#include "ui_terzaghi1d.h"

Terzaghi1D::Terzaghi1D(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::Terzaghi1D)
{
    ui->setupUi(this);
}

Terzaghi1D::~Terzaghi1D()
{
    delete ui;
}

void Terzaghi1D::getUserData()
{
    k=ui->kLine->text().toDouble();
    K=ui->KLine->text().toDouble();
    v=ui->vLine->text().toDouble();
    totalTime=ui->tLine->text().toDouble();
    p0=ui->poLine->text().toDouble();
    z=ui->zLine->text().toDouble();
    h=ui->HLine->text().toDouble();
    ns=ui->nsLine->text().toInt();

    double gf=9.81;
    result.resize(ns+1,2); //Time vs pore
    G=3.0*K*(1.0-2.0*v)/2.0/(1.0+v);
    double mv=1.0/(K+4.0*G/3.0);
    Cv=k/gf/mv;
    ui->CvLine->setText(QString::number(Cv));
}

void Terzaghi1D::calculateSolution()
{
    cout<<"z, h, p0, Cv, ns "<<z<<" "<<h<<" "<<p0<<" "<<Cv<<" "<<ns<<endl;
    double PI=3.1415926535897932385;
    double H=h;
    double z_loc=z;
    double u0_ana=p0;
    for (int i=0;i<=ns;i++)
    {
        double ratio=totalTime/ns;
        double t=i*ratio;
        double sa=0;
        for (int j=1;j<1000;j++)
        {
            double mu=exp(-(2*j-1)*(2*j-1)*PI*PI*0.25*Cv*t/H/H);
            double c=cos((2*j-1)*PI*z_loc*0.5/H);
            double f=pow(-1,j-1)/(2*j-1);
            sa=sa+f*c*mu;
            result(i,0)=t;
        }
      result(i,1)=4*u0_ana*sa/PI;
    }
}

void Terzaghi1D::on_closeButton_clicked()
{
    this->close();
}

void Terzaghi1D::on_okButton_clicked()
{
    getUserData();
    calculateSolution();
    QString folderName=QFileDialog::getExistingDirectory(Q_NULLPTR,"Save as");
    bool ok;
    QString fileName=QInputDialog::getText(Q_NULLPTR,"Chose file name","Save as: ",QLineEdit::Normal,"analytical",&ok);
    fileName=folderName+"/"+fileName+".txt";
    WriteToFile exportFile;
    exportFile.fileName=fileName.toStdString();
    exportFile.ToFile(result);
}
