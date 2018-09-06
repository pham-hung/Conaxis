#include "material.h"
#include "ui_material.h"

Material::Material(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::Material)
{
    ui->setupUi(this);
}

Material::~Material()
{
    delete ui;
}

void Material::getUserData()
{
    checkSuccess=true;
    matIndex=ui->matNum->text().toInt();
    KFunction=ui->KCombo->currentIndex();
    if(KFunction==0)
    {
        ui->KLine->setText("Constant value");
        KCurve.resize(1,1);
        KCurve(0,0)=ui->bulkValLine->text().toDouble();
        if(KCurve(0,0)<=0)
        {
            QMessageBox::warning(this,"ERROR","Value must be >0");
            checkSuccess=false;
            return;
        }
    }
    else
    {
        ui->KLine->setText("Input From File");
        fileNameK=ui->bulkValLine->text();

    }
    poission=ui->vLine->text().toDouble();
    if(poission<=0)
    {
        QMessageBox::warning(this,"ERROR","Value must be >0");
        checkSuccess=false;
        return;
    }

    kFunction=ui->kCombo->currentIndex();
    if(kFunction==0)
    {
        ui->kLine->setText("Constant Value");
        kCurve.resize(1,1);
        kCurve(0,0)=ui->kValLine->text().toDouble();
        if(kCurve(0,0)<=0)
        {
            QMessageBox::warning(this,"ERROR","Value must be >0");
            checkSuccess=false;
            return;
        }

    }
    else
    {
        ui->kLine->setText("Input From File");
        fileNamek=ui->kValLine->text();
    }

    gf=ui->densLine->text().toDouble();
    if(gf<=0)
    {
        QMessageBox::warning(this,"ERROR","Value must be >0");
        checkSuccess=false;
        return;
    }
    e0=ui->voidRatioLine->text().toDouble();
    if(e0<=0)
    {
        QMessageBox::warning(this,"ERROR","Value must be >0");
        checkSuccess=false;
        return;
    }

    ratio=ui->ratioLine->text().toDouble();
    if(ratio<=0)
    {
        QMessageBox::warning(this,"ERROR","Value must be >0");
        checkSuccess=false;
        return;
    }
    Cd=ui->CdLine->text().toDouble();
    if(Cd<=0)
    {
        Cd=1.0f;
        ui->CdLine->setText("1");
    }

}

void Material::sendSIGNAL()
{
    emit sendMaterial(mat);
}

int Material::NumberOfMat()
{
    return mat.size();
}

void Material::createCrsMaterial(double poissionRatio, double voidRatio)
{
    mat.resize(1);
    int index=0;
    mat[index].KFunction=0; //constant
    mat[index].KCurve.resize(1,1);
    mat[index].KCurve<<500;
    mat[index].poission=poissionRatio;
    mat[index].kFunction=0;
    mat[index].kCurve.resize(1,1);
    mat[index].kCurve<<1e-9;
    mat[index].gf=9.81;
    mat[index].e0=voidRatio;
    mat[index].fileNameK="";
    mat[index].fileNamek="";
    mat[index].ratio=1;
    mat[index].Cd=1;
}

void Material::on_cancelButton_clicked()
{
    this->close();
}

void Material::on_addButton_clicked()
{
    getUserData();
    if(checkSuccess==false)
    {
        return;
    }
    else
    {

        if(matIndex<=mat.size())
        {
            //Edit material
            int index=matIndex-1;
            mat[index].KFunction=KFunction;
            mat[index].KCurve=KCurve;
            mat[index].poission=poission;
            mat[index].kFunction=kFunction;
            mat[index].kCurve=kCurve;
            mat[index].gf=gf;
            mat[index].e0=e0;
            mat[index].fileNameK=fileNameK;
            mat[index].fileNamek=fileNamek;
            mat[index].ratio=ratio;
            mat[index].Cd=Cd;
            cout<<"Material: "<<matIndex<<" is modified"<<endl;
        }
        else
        {
            MaterialBase newMat;
            newMat.matIndex=matIndex;
            newMat.KFunction=KFunction;
            newMat.KCurve=KCurve;
            newMat.poission=poission;
            newMat.kFunction=kFunction;
            newMat.kCurve=kCurve;
            newMat.gf=gf;
            newMat.e0=e0;
            newMat.fileNameK=fileNameK;
            newMat.fileNamek=fileNamek;
            newMat.ratio=ratio;
            newMat.Cd=Cd;
            mat.push_back(newMat);
            cout<<"Material index: "<<matIndex<<" is created"<<endl;
        }
    }
}

void Material::on_KCombo_activated(int index)
{
    if(index==1)
    {
        QString fileName=QFileDialog::getOpenFileName();
        fileNameK=fileName;
        GetFile getFile;
        getFile.fileName=fileName;
        getFile.DoGetFile();
        if(getFile.data_file.rows()!=1)
        {
            bool ok;
            QString stressCol=QInputDialog::getText(this,"Input Value","Colum of stress",QLineEdit::Normal,"1",&ok);
            QString KCol=QInputDialog::getText(this,"Input Value","Colum of Bulk Modulus",QLineEdit::Normal,"2",&ok);
            if(stressCol.toInt()>getFile.data_file.cols()||KCol.toInt()>getFile.data_file.rows())
            {
              QMessageBox::warning(this,"ERROR","Input colum is larger than data file");
              return;
            }
            KCurve.resize(getFile.data_file.rows(),2);
            KCurve.col(0)=getFile.data_file.col(stressCol.toInt()-1);
            KCurve.col(1)=getFile.data_file.col(KCol.toInt()-1);
            ui->KLine->setText("From File:");
            ui->bulkValLine->setText(fileName);
        }
        else
        {
            ui->KCombo->setCurrentIndex(0);
            return;
        }
    }
    else if(index==0)
    {
        ui->KLine->setText("Constant Value");
        ui->bulkValLine->setText("");
    }
    else if (index==2)
    {
        QString fileName=QFileDialog::getOpenFileName();
        fileNameK=fileName;
        GetFile getFile;
        getFile.fileName=fileName;
        getFile.DoGetFile();
        if(getFile.data_file.rows()!=1)
        {
            bool ok;
            QString stressCol=QInputDialog::getText(this,"Input Value","Column of time",QLineEdit::Normal,"1",&ok);
            QString KCol=QInputDialog::getText(this,"Input Value","Colum of Bulk Modulus",QLineEdit::Normal,"2",&ok);
            QString factor=QInputDialog::getText(this,"Input Value","Scale time column by factor (convert minutes to days) ",QLineEdit::Normal,"6.94444444e-4",&ok);
            if(stressCol.toInt()>getFile.data_file.cols()||KCol.toInt()>getFile.data_file.rows())
            {
              QMessageBox::warning(this,"ERROR","Input colum is larger than data file");
              return;
            }
            KCurve.resize(getFile.data_file.rows(),2);
            KCurve.col(0)=factor.toDouble()*getFile.data_file.col(stressCol.toInt()-1);
            KCurve.col(1)=getFile.data_file.col(KCol.toInt()-1);
            ui->KLine->setText("From File:");
            ui->bulkValLine->setText(fileName);
        }
        else
        {
            ui->KCombo->setCurrentIndex(0);
            return;
        }
    }
}

void Material::on_listInfor_clicked()
{
    listMaterial();
}

void Material::on_matNum_textChanged(const QString &arg1)
{
    matIndex=ui->matNum->text().toInt();
    if(matIndex<1)
    {
        matIndex=1;
        ui->matNum->setText("1");
    }
    if(matIndex>mat.size())
    {
        matIndex=mat.size()+1;
        ui->matNum->setText(QString::number(matIndex));
    }
    if(matIndex<=mat.size())
    {
        int index=matIndex-1;        
        if(mat[index].KFunction==0)
        {
            ui->KLine->setText("Constant Value");
            ui->bulkValLine->setText(QString::number(mat[index].KCurve(0,0)));
            ui->KCombo->setCurrentIndex(0);
        }
        else
        {
            ui->KLine->setText("Input From File");
            ui->bulkValLine->setText(mat[index].fileNameK);
            KCurve=mat[index].KCurve;
            ui->KCombo->setCurrentIndex(mat[index].KFunction);
        }

        if(mat[index].kFunction==0)
        {
            ui->kLine->setText("Constant Value");
            ui->kValLine->setText(QString::number(mat[index].kCurve(0,0)));
            ui->kCombo->setCurrentIndex(0);
        }
        else
        {
            ui->kLine->setText("Input From File");
            ui->kValLine->setText(mat[index].fileNamek);
            kCurve=mat[index].kCurve;
            ui->kCombo->setCurrentIndex(mat[index].kFunction);
        }

        ui->vLine->setText(QString::number(mat[index].poission));
        ui->ratioLine->setText(QString::number(mat[index].ratio));
        ui->densLine->setText(QString::number(mat[index].gf));
        ui->voidRatioLine->setText(QString::number(mat[index].e0));
        ui->CdLine->setText(QString::number(mat[index].Cd));
    }
    else
    {
        ui->KCombo->setCurrentIndex(0);
        ui->bulkValLine->setText("");
        ui->kCombo->setCurrentIndex(0);
        ui->kValLine->setText("");
        ui->vLine->setText("");
        ui->ratioLine->setText("1");
        ui->densLine->setText("1");
        ui->voidRatioLine->setText("1.8");
        ui->CdLine->setText("1.0");
    }
}

void Material::on_kCombo_activated(int index)
{
    if(index==1)
    {
        QString fileName=QFileDialog::getOpenFileName();
        fileNamek=fileName;
        GetFile getFile;
        getFile.fileName=fileName;
        getFile.DoGetFile();
        if(getFile.data_file.rows()!=1)
        {
            bool ok;
            QString stressCol=QInputDialog::getText(this,"Input Value","Colum of stress",QLineEdit::Normal,"1",&ok);
            QString KCol=QInputDialog::getText(this,"Input Value","Colum of Hydraulic conductivity",QLineEdit::Normal,"3",&ok);
            kCurve.resize(getFile.data_file.rows(),2);
            kCurve.col(0)=getFile.data_file.col(stressCol.toInt()-1);
            kCurve.col(1)=getFile.data_file.col(KCol.toInt()-1);
            ui->kLine->setText("From File:");
            ui->kValLine->setText(fileName);
        }
        else
        {
            ui->kCombo->setCurrentIndex(0);
            return;
        }
    }
    else if(index==0)
    {
        ui->kLine->setText("Constant Value");
        ui->kValLine->setText("");
    }
    else if(index==2)
    {
        QString fileName=QFileDialog::getOpenFileName();
        fileNamek=fileName;
        GetFile getFile;
        getFile.fileName=fileName;
        getFile.DoGetFile();
        if(getFile.data_file.rows()!=1)
        {
            bool ok;
            QString stressCol=QInputDialog::getText(this,"Input Value","Colums of time: ",QLineEdit::Normal,"1",&ok);
            QString KCol=QInputDialog::getText(this,"Input Value","Colum of Hydraulic conductivity",QLineEdit::Normal,"3",&ok);
            QString factor=QInputDialog::getText(this,"Input Value","Scale time column by factor (convert minutes to days) ",QLineEdit::Normal,"6.94444444e-4",&ok);
            kCurve.resize(getFile.data_file.rows(),2);
            kCurve.col(0)=factor.toDouble()*getFile.data_file.col(stressCol.toInt()-1);
            kCurve.col(1)=getFile.data_file.col(KCol.toInt()-1);
            ui->kLine->setText("From File:");
            ui->kValLine->setText(fileName);
        }
        else
        {
            ui->kCombo->setCurrentIndex(0);
            return;
        }
    }
}

void Material::listMaterial()
{
    for (int i=0;i<mat.size();i++)
    {
        cout<<"--------------"<<endl;
        cout<<"Material number : "<<mat[i].matIndex<<endl;
        cout<<"Poission ratio  : "<<mat[i].poission<<endl;
        cout<<"Ratio kh/kv     : "<<mat[i].ratio<<endl;
        cout<<"Unit Weight     : "<<mat[i].gf<<endl;
        cout<<"Void ratio      : "<<mat[i].e0<<endl;
        cout<<"Cd factor       : "<<mat[i].Cd<<endl;
        if(mat[i].KFunction==0)
        {
            cout<<"Constanst K     : "<<mat[i].KCurve(0,0)<<endl;
        }
        else
        {
            cout<<"K is function"<<endl;
            cout<<KCurve<<endl;
            cout<<"-----"<<endl;
        }

        if(mat[i].kFunction==0)
        {
            cout<<"Constanst k     : "<<mat[i].kCurve(0,0)<<endl;
        }
        else
        {
            cout<<"k is function"<<endl;
            cout<<kCurve<<endl;
            cout<<"-----"<<endl;
        }
    }
}

void Material::getMaterial(vector<MaterialBase> mat)
{
    this->mat=mat;
}

void Material::on_setCd_clicked()
{
    bool ok;
    QString value=QInputDialog::getText(this,"Input Value","Cd=",QLineEdit::Normal,"1",&ok);

    for (int i=0;i<mat.size();i++)
    {
        mat[i].Cd=value.toDouble();
    }
}

void Material::on_setRatio_clicked()
{
    bool ok;
    QString value=QInputDialog::getText(this,"Input Value","Ratio(kh/kv)=",QLineEdit::Normal,"1",&ok);

    for (int i=0;i<mat.size();i++)
    {
        mat[i].ratio=value.toDouble();
    }
}
