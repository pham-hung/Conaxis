#ifndef WRITETOFILE_H
#define WRITETOFILE_H
#include <string>
#include <iostream>
#include <ostream>
#include <Eigen/Dense>
#include <QString>
#include <QStringList>
using namespace std;
using namespace Eigen;
using Eigen::MatrixXd;

class WriteToFile
{
public:    
    string fileName;    
    void ToFile(const Ref<const MatrixXd>matrixA);
    void ToFile(QString fileName, Ref<MatrixXd> matrixA);
    void ToFile(QStringList header, Ref<MatrixXd> matrixA);
    void ToFile(QString fileName,QStringList header,Ref<MatrixXd> matrixA);
    WriteToFile();
};

#endif // WRITETOFILE_H
