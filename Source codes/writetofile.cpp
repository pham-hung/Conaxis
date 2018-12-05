#include "writetofile.h"
#include<iostream>
#include<string>
#include<fstream>
#include<iomanip>
#include<QString>

void WriteToFile::ToFile(const Ref<const MatrixXd> matrixA)
{
    int row=matrixA.rows();
    int col=matrixA.cols();
    double val=0;
    ofstream file_;
    file_.open(fileName);

    for (int i=0;i<row;i=i+1)
    {
        for (int j=0;j<col;j=j+1)
        {
          val=matrixA(i,j);
          file_<<std::setw(13)<<std::setprecision(5)<<std::scientific<<val<<'\t';
        }
        file_<<'\n';
    }
    file_.close();
    cout<<"Exported successfully: "<<fileName<<endl;
}

void WriteToFile::ToFile(QString fileName, Ref<MatrixXd> matrixA)
{
    this->fileName=fileName.toStdString();
    ToFile(matrixA);
}

void WriteToFile::ToFile(QStringList header, Ref<MatrixXd> matrixA)
{
    int row=matrixA.rows();
    int col=matrixA.cols();
    double val=0;
    ofstream streamFile_;
    streamFile_.open(fileName);
    QString headerString;

    for (auto i=0;i<header.size();i++)
    {
        headerString=header[i];
        streamFile_<<std::setw(13)<<std::internal<<headerString.toStdString()<<'\t';
    }
    streamFile_<<'\n';

    for (int i=0;i<row;i=i+1)
    {
        for (int j=0;j<col;j=j+1)
        {
          val=matrixA(i,j);
          streamFile_<<std::setw(13)<<std::right<<std::setprecision(5)<<std::scientific<<val<<'\t';
        }
        streamFile_<<'\n';
    }
    streamFile_.close();
}

void WriteToFile::ToFile(QString fileName, QStringList header, Ref<MatrixXd> matrixA)
{
    this->fileName=fileName.toStdString();
    ToFile(header,matrixA);
}

WriteToFile::WriteToFile()
{

}
