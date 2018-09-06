#ifndef GETFILE_H
#define GETFILE_H
#include <QFile>
#include <QString>
#include <Eigen/Dense>
#include <QDebug>
#include <QStringList>
#include <QTextStream>
#include <iostream>
#include <QObject>
#include <QThread>

using namespace Eigen;
using namespace std;

class GetFile:public QObject
{
    Q_OBJECT

public slots:
    void DoGetFile();

signals:
    void finishGetFile();
    void startGetFile();
    void failToOpen();

public:
    bool sucessFlag=false;
    GetFile();
    int col;
    int row;
    bool checkFinish=false;
    bool doneCheck=false;
    MatrixXd data_file=MatrixXd::Zero(1,1);
    QString fileName,folderName;
    void get_dimension(QString fileName);
    void import_file(QString fileName);
private:


};

#endif // GETFILE_H
