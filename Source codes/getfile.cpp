#include "getfile.h"

void GetFile::get_dimension(QString fileName)
{
    row=0;
    col=0;
    QFile File(fileName);
    QString Line,Temp;
    QTextStream in(&File);

    if(!File.open(QFile::ReadOnly|QFile::Text))
    {
        cout<<"Cannot open the file"<<endl;
        doneCheck=false;
        return;
    }
    else
    {
        while(!in.atEnd())
        {
            Line=in.readLine();
            row=row+1;
            //qDebug()<<Line<<endl;
            if (row==1)
            {
               QStringList list=Line.split(QRegExp("\\s+"),QString::SkipEmptyParts);
               col=list.count();
            }
        }
    }
    File.close();
    doneCheck=true;
    //cout<<"The numer of row is:"<<row<<endl;
    //cout<<"The numer of col is:"<<col<<endl;
}

void GetFile::import_file(QString fileName)
{
    //qDebug()<<"Write file to Matrix"<<endl;
    //qDebug()<<"Row and col is:"<<row<<col;
    data_file.resize(row,col);
    double Value=0;
    QFile File(fileName);
    QString Line,Temp;
    QTextStream in(&File);

    if(!File.open(QFile::ReadOnly|QFile::Text))
    {
        qDebug()<<"cannot open the file"<<endl;
    }
    else
    {
        for(int i=0;i<row;i++)
        {
            Line=in.readLine();
            QStringList list=Line.split(QRegExp("\\s+"),QString::SkipEmptyParts);
            for (int j=0;j<col;j++)
            {
                Temp=list[j];
                Value=Temp.toDouble();
                data_file(i,j)=Value;
            }
        }

    }

}

void GetFile::DoGetFile()
{
    sucessFlag=false;
    checkFinish=false;
    data_file.setZero();
    data_file.resize(1,1);
    GetFile::get_dimension(fileName);
    if(doneCheck==true)
    {
        GetFile::import_file(fileName);
        cout<<"Finished Importing data from File:"<<fileName.toStdString()<<endl;
        //cout<<"Finished, Thread ID is called: "<<QThread::currentThreadId()<<endl;
        checkFinish=true;
        emit finishGetFile();
        sucessFlag=true;
    }
    else
    {
        emit failToOpen();
        checkFinish=true;
        emit finishGetFile();

    }
    doneCheck=false;
}

GetFile::GetFile()
{
}
