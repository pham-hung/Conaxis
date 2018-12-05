#include "savedatabase.h"

SaveDataBase::SaveDataBase(QObject *parent) : QObject(parent)
{

}

void SaveDataBase::saveModel()
{
    //Get folder name and file name
    if(saveClicked==false)
    {
        QString selfilter;
        fileName=QFileDialog::getSaveFileName(Q_NULLPTR,"Chose save file name","",
                 tr("coa (*.coa)"),&selfilter);
        saveClicked=true;
    }

    //Delete file and create new file
    mFile.setFileName(fileName);
    if(mFile.exists()==true){mFile.remove();}
    qDebug()<<"File name: "<<fileName<<endl;

    mFile.setFileName(fileName);
    if (mFile.open(QIODevice::ReadWrite)) {
        QTextStream stream(&mFile);
        //Save PROJECT
        stream<<"#PRJOJECT"<<endl;
        for(int i=0;i<projectParameters.size();i++)
        {
            stream<<projectParameters[i]<<" ";
        }
        stream<<endl;
        //Save dimension

        stream<<"#GEOMETRY"<<endl;
        stream<<geometry.re<<" "<<geometry.rw<<" "<<geometry.rs<<" "<<geometry.length<<" "<<geometry.surfaceElevation<<endl;
        stream<<geometry.numberOfLayer<<" "<<geometry.numberOfElementSmear<<" "<<geometry.numberOfElementSoil<<" "<<geometry.defaultSubLayer<<" "<<geometry.qw<<" "<<geometry.analysisType<<endl;
        stream<<"#LAYERINFO"<<endl;
        stream<<geometry.layerInfo.rows()<<" "<<geometry.layerInfo.cols()<<endl;
        stream.setRealNumberNotation(QTextStream::ScientificNotation);
        stream.setRealNumberPrecision(5);
        for(int i=0;i<geometry.layerInfo.rows();i++)
        {
            for (int j=0;j<geometry.layerInfo.cols();j++)
            {
                stream<<geometry.layerInfo(i,j)<<" ";
            }
            stream<<endl;
        }

        stream<<"#DIMENSION"<<endl;
        stream<<coordinates.rows()<<" "<<coordinates.cols()<<" "<<elements.rows()<<" "<<elements.cols()<<endl;

        //SAVE COORDINATES
        stream<<"#COORD"<<endl;
        stream.setRealNumberNotation(QTextStream::ScientificNotation);
        stream.setRealNumberPrecision(5);
        for(int i=0;i<coordinates.rows();i++)
        {
            for (int j=0;j<coordinates.cols();j++)
            {
                stream<<coordinates(i,j)<<" ";
            }
            stream<<endl;
        }

        //SAVE ELEMENTS
        stream.setRealNumberNotation(QTextStream::FixedNotation);
        stream.setRealNumberPrecision(0);
        stream<<"#ELEMENT"<<endl;
        for(int i=0;i<elements.rows();i++)
        {
            for (int j=0;j<elements.cols();j++)
            {
                stream<<elements(i,j)<<" ";
            }
            stream<<endl;
        }

        //Save Material
        stream<<"#MATNUM"<<endl;
        stream.setRealNumberNotation(QTextStream::SmartNotation);
        stream.setRealNumberPrecision(5);
        stream<<material.size()<<endl;
        for(int i=0;i<material.size();i++)
        {
            stream<<"#MAT"<<" "<<material[i].matIndex<<endl;
            stream<<material[i].KFunction<<" "
                 <<material[i].poission<<" "<<material[i].kFunction<<" "
                <<material[i].ratio<<" "<<material[i].gf<<" "<<material[i].e0<<" "<<material[i].Cd<<endl;

            stream.setRealNumberNotation(QTextStream::ScientificNotation);
            stream.setRealNumberPrecision(5);
            stream<<"#BULK"<<endl;
            stream<<material[i].fileNameK<<endl;
            for(int ii=0;ii<material[i].KCurve.rows();ii++)
            {
                for (int jj=0;jj<material[i].KCurve.cols();jj++)
                {
                    stream<<material[i].KCurve(ii,jj)<<" ";
                }
                stream<<endl;
            }
            stream<<"#HYDRAU"<<endl;
            stream<<material[i].fileNameK<<endl;
            for(int ii=0;ii<material[i].kCurve.rows();ii++)
            {
                for (int jj=0;jj<material[i].kCurve.cols();jj++)
                {
                    stream<<material[i].kCurve(ii,jj)<<" ";
                }
                stream<<endl;
            }
        }

        //Save stage, add graviload at the end
        stream.setRealNumberNotation(QTextStream::SmartNotation);
        stream.setRealNumberPrecision(5);
        stream<<"#STAGENUM"<<endl;
        stream<<stage.size()<<endl;
        for(int i=0;i<stage.size();i++)
        {
            stream<<"#STAGE"<<" "<<stage[i].stageIndex<<endl;
            stream<<stage[i].stageName<<endl;
            stream<<stage[i].stageType<<" "<<stage[i].dt<<" "<<stage[i].t0<<" "<<stage[i].t1<<" "<<stage[i].subStep<<" "<<stage[i].gravityLoad<<" "<<stage[i].timeStepType<<endl;
            stream<<stage[i].fileName<<endl;
            stream<<"#TIMESTEP"<<endl;
            for (int ii=0;ii<stage[i].timeStep.rows();ii++)
            {
                stream<<stage[i].timeStep(ii,0)<<endl;
            }
        }

        //Save boundary
        stream<<"#BOUNDARYNUM"<<endl;
        stream<<boundary.size()<<endl;
        for(int i=0;i<boundary.size();i++)
        {
            stream<<"#BOUNDARY"<<" "<<boundary[i].boundaryIndex<<endl;
            stream<<boundary[i].boundaryName<<endl;
            stream<<boundary[i].boundaryType<<" "<<boundary[i].loadType<<" "<<boundary[i].val0<<" "<<boundary[i].val1<<endl;
            stream.setRealNumberNotation(QTextStream::ScientificNotation);
            stream.setRealNumberPrecision(5);
            stream<<"#TIMECURVE"<<endl;
            stream<<boundary[i].fileName<<endl;
            for(int ii=0;ii<boundary[i].loadCurve.rows();ii++)
            {
                for (int jj=0;jj<boundary[i].loadCurve.cols();jj++)
                {
                    stream<<boundary[i].loadCurve(ii,jj)<<" ";
                }
                stream<<endl;
            }
        }

        //Save boundary to stage
        stream<<"#STAGEBOUNDARYNUM"<<endl;
        stream<<stageBoundary.size()<<endl;
        for (int i=0;i<stageBoundary.size();i++)
        {
            stream.setRealNumberNotation(QTextStream::FixedNotation);
            stream.setRealNumberPrecision(0);
            int NoS=stageBoundary[i].v_stageBoundaryIndex.size();
            stream<<"#STAGEBOUNDARY"<<" "<<i+1<<endl;
            stream<<"#BCNUM"<<" "<<stageBoundary[i].v_stageBoundaryIndex.size()<<endl;

            stream<<"#BCINDEX"<<endl;
            for(int j=0;j<NoS;j++)
            {
                stream<<stageBoundary[i].v_stageBoundaryIndex[j]<<" ";
            }
            stream<<endl;


            stream<<"#NODELIST"<<endl;
            for(int j=0;j<NoS;j++)
            {
                stream<<"#NODELIST "<<j+1<<endl;
                for(int jj=0;jj<stageBoundary[i].v_nodeList[j].rows();jj++)
                {
                    stream<<stageBoundary[i].v_nodeList[j](jj,0)<<endl;
                }
            }

            stream<<"#FILENAME"<<endl;
            for(int j=0;j<NoS;j++)
            {
                stream<<stageBoundary[i].v_nodeFileName[j]<<endl;
            }

            stream<<"#ASSIGNTYPE"<<endl;
            for(int j=0;j<NoS;j++)
            {
                stream<<stageBoundary[i].v_assignType[j]<<" ";
            }
            stream<<endl;

            stream<<"#BOUNDARYINDEX"<<endl;
            for(int j=0;j<NoS;j++)
            {
                stream<<stageBoundary[i].v_boundaryIndex[j]<<" ";
            }
            stream<<endl;

            stream.setRealNumberNotation(QTextStream::SmartNotation);
            stream.setRealNumberPrecision(0);
            stream<<"#GRADIENT"<<endl;
            for(int j=0;j<NoS;j++)
            {
                stream<<stageBoundary[i].v_gradientBool[j]<<" ";
            }
            stream<<endl;

            stream.setRealNumberPrecision(5);
            stream<<"#AFACTOR"<<endl;
            for(int j=0;j<NoS;j++)
            {
                stream<<stageBoundary[i].v_aFactor[j]<<" ";
            }
            stream<<endl;

            stream<<"#BFACTOR"<<endl;
            for(int j=0;j<NoS;j++)
            {
                stream<<stageBoundary[i].v_bFactor[j]<<" ";
            }
            stream<<endl;

            stream<<"#CFACTOR"<<endl;
            for(int j=0;j<NoS;j++)
            {
                stream<<stageBoundary[i].v_cFactor[j]<<" ";
            }
            stream<<endl;

            stream<<"#X0"<<endl;
            for(int j=0;j<NoS;j++)
            {
                stream<<stageBoundary[i].v_x0[j]<<" ";
            }
            stream<<endl;

            stream<<"#X1"<<endl;
            for(int j=0;j<NoS;j++)
            {
                stream<<stageBoundary[i].v_x1[j]<<" ";
            }
            stream<<endl;

            stream<<"#Y0"<<endl;
            for(int j=0;j<NoS;j++)
            {
                stream<<stageBoundary[i].v_y0[j]<<" ";
            }
            stream<<endl;

            stream<<"#Y1"<<endl;
            for(int j=0;j<NoS;j++)
            {
                stream<<stageBoundary[i].v_y1[j]<<" ";
            }
            stream<<endl;
        }

        //Save watchlist
        stream<<"#WATCHLIST"<<endl;
        stream<<"#WATCHNUM"<<" "<<watch.size()<<endl;
        for (int j=0;j<watch.size();j++)
        {
            stream<<"#WATCH"<<" "<<watch[j].watchIndex<<endl;
            stream<<watch[j].title<<endl;
            stream.setRealNumberNotation(QTextStream::ScientificNotation);
            stream.setRealNumberPrecision(5);
            stream<<watch[j].watchType<<" "<<watch[j].x0<<" "<<watch[j].x1<<" "<<watch[j].y0<<" "<<watch[j].y1<<" "<<watch[j].beginStep<<" "<<watch[j].endStep<<" "<<watch[j].averageBool<<endl;
        }
    }
    mFile.close();
    mFile.flush();
    qDebug()<<"Save sucessfully"<<endl;
}

void SaveDataBase::saveAsModel()
{
    saveClicked=false;
    saveModel();
}

void SaveDataBase::readModel()
{    
    projectParameters.resize(0);
    QString selfilter;
    fileName=QFileDialog::getOpenFileName(Q_NULLPTR,"Open data file","",
             tr("coa (*.coa);;text (*.txt *.dat);;All files (*.*)" ),&selfilter);

    int nol=0;
    int non=0, noe=0, nom=0, nob=0, nos=0, nosb=0;
    int nobb;
    bool ok=false;
    QString Line, Temp;
    QStringList list;
    qDebug()<<"Resume from file name: "<<fileName<<endl;
    mFile.setFileName(fileName);
    if(!mFile.exists())
    {
        qDebug()<<"File is not exists, try again"<<endl;
        success=false;
        return;
    }

    //Read from top to end
    mFile.setFileName(fileName);
    if (mFile.open(QIODevice::ReadWrite))
    {
        QTextStream stream(&mFile);
        QString Line;
        while(!stream.atEnd())
        {
            Line=stream.readLine();
            list=Line.split(QRegExp("\\s+"),QString::SkipEmptyParts);
            if(Line=="#PRJOJECT") //Project
            {
                Line=stream.readLine();
                list=Line.split(QRegExp("\\s+"),QString::SkipEmptyParts);
                if(list.count()>0)
                {
                    for(int j=0;j<list.count();j++)
                    {
                        Temp=list[j];
                        projectParameters.push_back(Temp.toDouble());
                    }
                }

            }

            if(Line=="#GEOMETRY")
            {
                Line=stream.readLine();
                list=Line.split(QRegExp("\\s+"),QString::SkipEmptyParts);
                Temp=list[0];
                geometry.re=Temp.toDouble();
                Temp=list[1];
                geometry.rw=Temp.toDouble();
                Temp=list[2];
                geometry.rs=Temp.toDouble();
                Temp=list[3];
                geometry.length=Temp.toDouble();
                Temp=list[4];
                geometry.surfaceElevation=Temp.toDouble();

                Line=stream.readLine();
                list=Line.split(QRegExp("\\s+"),QString::SkipEmptyParts);
                Temp=list[0];
                geometry.numberOfLayer=Temp.toInt();
                Temp=list[1];
                geometry.numberOfElementSmear=Temp.toInt();
                Temp=list[2];
                geometry.numberOfElementSoil=Temp.toInt();
                Temp=list[3];
                geometry.defaultSubLayer=Temp.toInt();
                Temp=list[4];
                geometry.qw=Temp.toDouble();
                Temp=list[5];
                geometry.analysisType=Temp.toInt();
            }

            if(Line=="#LAYERINFO")
            {
                Line=stream.readLine();
                list=Line.split(QRegExp("\\s+"),QString::SkipEmptyParts);
                Temp=list[0];
                int layerInforRow=Temp.toInt();
                Temp=list[1];
                int layerInfoCol=Temp.toInt();
                geometry.layerInfo.resize(layerInforRow,layerInfoCol);
                for(int j=0;j<layerInforRow;j++)
                {
                    Line=stream.readLine();
                    list=Line.split(QRegExp("\\s+"),QString::SkipEmptyParts);
                    for (int jj=0;jj<layerInfoCol;jj++)
                    {
                        Temp=list[jj];
                        geometry.layerInfo(j,jj)=Temp.toDouble();
                    }
                }
            }

            if(Line=="#DIMENSION") //Dimension
            {
                Line=stream.readLine();
                list=Line.split(QRegExp("\\s+"),QString::SkipEmptyParts);
                Temp=list[0];
                non=Temp.toDouble();
                coordinates.resize(non,4);
                Temp=list[2];
                noe=Temp.toDouble();
                elements.resize(noe,12);

            }

            if(Line=="#COORD") //Coordinates
            {
                if(non>1)
                {
                    for(int j=0;j<non;j++)
                    {
                        Line=stream.readLine();
                        list=Line.split(QRegExp("\\s+"),QString::SkipEmptyParts);
                        for (int jj=0;jj<4;jj++)
                        {
                            Temp=list[jj];
                            coordinates(j,jj)=Temp.toDouble();
                        }
                    }
                }

            }

            if(Line=="#ELEMENT") //Elements
            {
                if(noe>1)
                {
                    for(int j=0;j<noe;j++)
                    {
                        Line=stream.readLine();
                        list=Line.split(QRegExp("\\s+"),QString::SkipEmptyParts);
                        for (int jj=0;jj<12;jj++)
                        {
                            Temp=list[jj];
                            elements(j,jj)=Temp.toDouble();
                        }
                    }
                }

            }

            if(Line=="#MATNUM")
            {
                Line=stream.readLine();
                list=Line.split(QRegExp("\\s+"),QString::SkipEmptyParts);
                Temp=list[0];
                nom=Temp.toDouble();
                material.resize(nom);
            }

            int ii=0;
            while(ii<nom)
            {
                QString LineCheck, LineStop1, LineStop2;
                LineCheck="#MAT "+QString::number(ii+1,'f',0);
                LineStop1="#MAT "+QString::number(ii+2,'f',0);
                LineStop2="#STAGENUM";
                list=Line.split(QRegExp("\\s+"),QString::SkipEmptyParts);
                if(Line==LineCheck)
                {
                    Line=stream.readLine();
                    list=Line.split(QRegExp("\\s+"),QString::SkipEmptyParts);

                    int matIndex=ii+1;
                    Temp=list[0];
                    int KFunction=Temp.toInt();

                    Temp=list[1];
                    double poission=Temp.toDouble();

                    Temp=list[2];
                    int kFunction=Temp.toInt();

                    Temp=list[3];
                    double ratio=Temp.toDouble();

                    Temp=list[4];
                    double gf=Temp.toDouble();

                    Temp=list[5];
                    double e0=Temp.toDouble();

                    Temp=list[6];
                    double Cd=Temp.toDouble();

                    material[ii].matIndex=matIndex;
                    material[ii].KFunction=KFunction;
                    material[ii].poission=poission;
                    material[ii].kFunction=kFunction;
                    material[ii].ratio=ratio;
                    material[ii].gf=gf;
                    material[ii].e0=e0;
                    material[ii].Cd=Cd;

                }
                Line=stream.readLine();

                if(Line=="#BULK")
                {
                    Line=stream.readLine();
                    material[ii].fileNameK=Line;
                    vector<double> vectorStress;
                    vector<double> vectorVal;
                    vectorStress.resize(0);
                    vectorVal.resize(0);
                    while(Line!="#HYDRAU")
                    {
                        Line=stream.readLine();
                        if(Line=="#HYDRAU")
                        {
                            break;
                        }
                        list=Line.split(QRegExp("\\s+"),QString::SkipEmptyParts);


                        if(list.count()==1)
                        {
                            material[ii].KCurve.resize(1,1);
                            material[ii].KCurve(0,0)=Line.toDouble();
                        }
                        else
                        {
                            Temp=list[0];
                            vectorStress.push_back(Temp.toDouble());
                            Temp=list[1];
                            vectorVal.push_back(Temp.toDouble());
                        }

                    }//End while

                    if(vectorVal.size()>0)
                    {
                        material[ii].KCurve.resize(vectorVal.size(),2);
                        for(int kk=0;kk<vectorVal.size();kk++)
                        {
                            material[ii].KCurve(kk,0)=vectorStress[kk];
                            material[ii].KCurve(kk,1)=vectorVal[kk];
                        }
                    }

                }//End if bulk

                if(Line=="#HYDRAU")
                {

                    Line=stream.readLine();
                    material[ii].fileNamek=Line;

                    vector<double> vectorStress;
                    vector<double> vectorVal;
                    vectorStress.resize(0);
                    vectorVal.resize(0);
                    while(Line!=LineStop1||Line!=LineStop2)
                    {
                        Line=stream.readLine();

                        if(Line==LineStop1||Line==LineStop2)
                        {
                            break;
                        }
                        list=Line.split(QRegExp("\\s+"),QString::SkipEmptyParts);
                        if(list.count()==1)
                        {
                            material[ii].kCurve.resize(1,1);
                            material[ii].kCurve(0,0)=Line.toDouble();
                        }
                        else
                        {
                            Temp=list[0];
                            vectorStress.push_back(Temp.toDouble());
                            Temp=list[1];
                            vectorVal.push_back(Temp.toDouble());
                        }
                    }//End while

                    if(vectorVal.size()>0)
                    {
                        material[ii].kCurve.resize(vectorVal.size(),2);
                        for(int kk=0;kk<vectorVal.size();kk++)
                        {
                            material[ii].kCurve(kk,0)=vectorStress[kk];
                            material[ii].kCurve(kk,1)=vectorVal[kk];
                        }
                    }
                    ii=ii+1;
                }//End if Hydrau
            }//En loop over Material


            if(Line=="#STAGENUM")
            {
                Line=stream.readLine();
                nos=Line.toInt();
                stage.resize(nos);

                int ii=0;
                while(ii<nos)
                {
                    QString LineCheck, LineStop1, LineStop2;
                    LineCheck="#STAGE "+QString::number(ii+1,'f',0);
                    LineStop1="#STAGE "+QString::number(ii+2,'f',0);
                    LineStop2="#BOUNDARYNUM";

                    if(Line==LineStop1||Line==LineStop2)
                    {
                        break;
                    }
                    if(ii==nos)
                    {
                        break;
                    }

                    if(Line==LineCheck)
                    {
                        stage[ii].stageIndex=ii+1;
                        Line=stream.readLine();
                        stage[ii].stageName=Line;

                        Line=stream.readLine();
                        list=Line.split(QRegExp("\\s+"),QString::SkipEmptyParts);

                        Temp=list[0];
                        stage[ii].stageType=Temp.toInt();
                        Temp=list[1];
                        stage[ii].dt=Temp.toDouble();
                        Temp=list[2];
                        stage[ii].t0=Temp.toDouble();
                        Temp=list[3];
                        stage[ii].t1=Temp.toDouble();
                        Temp=list[4];
                        stage[ii].subStep=Temp.toInt();
                        Temp=list[5];
                        stage[ii].gravityLoad=Temp.toInt();
                        Temp=list[6];
                        stage[ii].timeStepType=Temp.toInt();

                        Line=stream.readLine(); //Name of time Step
                        stage[ii].fileName=Line;
                    }

                    Line=stream.readLine(); //#TIMESTEP

                    if(Line=="#TIMESTEP")
                    {
                        vector<double> vectorValue;
                        vectorValue.resize(0);
                        while(Line!=LineStop1 && Line!=LineStop2)
                        {
                            Line=stream.readLine();
                            if(Line==LineStop1||Line==LineStop2)
                            {
                                break;
                            }
                            else
                            {
                                vectorValue.push_back(Line.toDouble());
                            }
                        }

                        if(vectorValue.size()>1)
                        {
                            stage[ii].timeStep.resize(vectorValue.size(),1);
                            for(int kk=0;kk<vectorValue.size();kk++)
                            {
                                stage[ii].timeStep(kk,0)=vectorValue[kk];
                            }
                        }
                        else
                        {
                            stage[ii].timeStep.resize(1,1);
                        }
                        ii=ii+1;
                    }
                }
            }

            if(Line=="#BOUNDARYNUM")
            {
                Line=stream.readLine();
                nob=Line.toInt();
                boundary.resize(nob);

                int ii=0;
                while(ii<nob)
                {
                    QString LineCheck, LineStop1, LineStop2;
                    LineCheck="#BOUNDARY "+QString::number(ii+1,'f',0);
                    LineStop1="#BOUNDARY "+QString::number(ii+2,'f',0);
                    LineStop2="#STAGEBOUNDARYNUM";

                    if(Line==LineStop1||Line==LineStop2)
                    {
                        cout<<"break is called"<<endl;
                        break;
                    }
                    if(ii==nob)
                    {
                        break;
                    }

                    if(Line==LineCheck)
                    {
                        boundary[ii].boundaryIndex=ii+1;
                        Line=stream.readLine();
                        boundary[ii].boundaryName=Line;
                        Line=stream.readLine();
                        list=Line.split(QRegExp("\\s+"),QString::SkipEmptyParts);
                        Temp=list[0];
                        boundary[ii].boundaryType=Temp.toInt();
                        Temp=list[1];
                        boundary[ii].loadType=Temp.toInt();
                        Temp=list[2];
                        boundary[ii].val0=Temp.toDouble();
                        Temp=list[3];
                        boundary[ii].val1=Temp.toDouble();
                    }
                    Line=stream.readLine();
                    if(Line=="#TIMECURVE")
                    {
                        Line=stream.readLine();
                        boundary[ii].fileName=Line;
                        vector<double> vectorTime;
                        vector<double> vectorValue;
                        while(Line!=LineStop1&&Line!=LineStop2)
                        {
                            Line=stream.readLine();
                            if(Line==LineStop1||Line==LineStop2)
                            {
                                break;
                            }
                            list=Line.split(QRegExp("\\s+"),QString::SkipEmptyParts);
                            if(list.count()==1)
                            {
                                boundary[ii].loadCurve=MatrixXd::Zero(1,1);
                            }
                            else
                            {
                                Temp=list[0];
                                vectorTime.push_back(Temp.toDouble());
                                Temp=list[1];
                                vectorValue.push_back(Temp.toDouble());
                            }
                        } //read matrix
                        if(vectorTime.size()>0)
                        {
                            boundary[ii].loadCurve.resize(vectorTime.size(),2);
                            for(int kk=0;kk<vectorTime.size();kk++)
                            {
                                boundary[ii].loadCurve(kk,0)=vectorTime[kk];
                                boundary[ii].loadCurve(kk,1)=vectorValue[kk];
                            }
                        }

                        ii=ii+1;
                    }

                } //end loop over eachh boundary
            }

            if(Line=="#STAGEBOUNDARYNUM")
            {

                Line=stream.readLine();
                nosb=Line.toInt();
                stageBoundary.resize(nosb);
                int ii=0;
                while(ii<nosb)
                {
                    QString LineCheck, LineStop1, LineStop2;
                    LineCheck="#STAGEBOUNDARY "+QString::number(ii+1,'f',0);
                    LineStop1="#STAGEBOUNDARY "+QString::number(ii+2,'f',0);
                    LineStop2="#WATCHLIST";
                    if(Line==LineStop1 || Line==LineStop2)
                    {
                        break;
                    }
                    if(Line==LineCheck)
                    {
                        Line=stream.readLine();
                        list=Line.split(QRegExp("\\s+"),QString::SkipEmptyParts);
                        Temp=list[1];
                        nobb=Temp.toInt();
                        stageBoundary[ii].v_stageBoundaryIndex.resize(nobb);
                        stageBoundary[ii].v_nodeFileName.resize(nobb);
                        stageBoundary[ii].v_assignType.resize(nobb);
                        stageBoundary[ii].v_boundaryIndex.resize(nobb);
                        stageBoundary[ii].v_x0.resize(nobb);
                        stageBoundary[ii].v_x1.resize(nobb);
                        stageBoundary[ii].v_y0.resize(nobb);
                        stageBoundary[ii].v_y1.resize(nobb);
                        stageBoundary[ii].v_nodeList.resize(nobb);
                        stageBoundary[ii].v_gradientBool.resize(nobb);
                        stageBoundary[ii].v_aFactor.resize(nobb);
                        stageBoundary[ii].v_bFactor.resize(nobb);
                        stageBoundary[ii].v_cFactor.resize(nobb);
                    }
                    if(Line=="#BCINDEX")
                    {
                        Line=stream.readLine();
                        list=Line.split(QRegExp("\\s+"),QString::SkipEmptyParts);
                        for(int kk=0;kk<nobb;kk++)
                        {
                            Temp=list[kk];
                            stageBoundary[ii].v_stageBoundaryIndex[kk]=Temp.toInt();
                        }

                    }

                    if(Line=="#NODELIST")
                    {
                        int jj=0;
                        vector<int> nodeVector;
                        nodeVector.resize(0);
                        int nn=0;

                        Line=stream.readLine();

                        while(nn<nobb)
                        {
                            QString LineCheck, LineStop1, LineStop2;
                            LineCheck="#NODELIST "+QString::number(nn+1,'f',0);
                            LineStop1="#NODELIST "+QString::number(nn+2,'f',0);
                            LineStop2="#FILENAME";
                            nodeVector.resize(0);
                            if(Line==LineStop2)
                            {
                                break;
                            }
                            if(Line==LineCheck)
                            {
                                while(1>0)
                                {
                                    Line=stream.readLine();
                                    if(Line==LineStop1||Line==LineStop2)
                                    {
                                        break;
                                    }
                                    nodeVector.push_back(Line.toInt());
                                }

                                int nodeNumber=nodeVector.size();
                                stageBoundary[ii].v_nodeList[nn].resize(nodeNumber,1);
                                for(int kk=0;kk<nodeVector.size();kk++)
                                {
                                    stageBoundary[ii].v_nodeList[nn](kk,0)=nodeVector[kk];
                                }
                                nn=nn+1;
                            }
                        }
                    }

                    if(Line=="#FILENAME")
                    {
                        for(int kk=0;kk<nobb;kk++)
                        {
                            Line=stream.readLine();
                            stageBoundary[ii].v_nodeFileName[kk]=Line;
                        }
                    }

                    if(Line=="#ASSIGNTYPE")
                    {
                        Line=stream.readLine();
                        list=Line.split(QRegExp("\\s+"),QString::SkipEmptyParts);
                        for(int kk=0;kk<nobb;kk++)
                        {
                            Temp=list[kk];
                            stageBoundary[ii].v_assignType[kk]=Temp.toInt();
                        }
                    }

                    if(Line=="#BOUNDARYINDEX")
                    {
                        Line=stream.readLine();
                        list=Line.split(QRegExp("\\s+"),QString::SkipEmptyParts);
                        for(int kk=0;kk<nobb;kk++)
                        {
                            Temp=list[kk];
                            stageBoundary[ii].v_boundaryIndex[kk]=Temp.toInt();
                        }

                    }

                    if(Line=="#GRADIENT")
                    {
                        Line=stream.readLine();
                        list=Line.split(QRegExp("\\s+"),QString::SkipEmptyParts);
                        for(int kk=0;kk<nobb;kk++)
                        {
                            Temp=list[kk];
                            stageBoundary[ii].v_gradientBool[kk]=Temp.toInt();
                        }

                    }

                    if(Line=="#AFACTOR")
                    {
                        Line=stream.readLine();
                        list=Line.split(QRegExp("\\s+"),QString::SkipEmptyParts);
                        for(int kk=0;kk<nobb;kk++)
                        {
                            Temp=list[kk];
                            stageBoundary[ii].v_aFactor[kk]=Temp.toDouble();
                        }

                    }

                    if(Line=="#BFACTOR")
                    {
                        Line=stream.readLine();
                        list=Line.split(QRegExp("\\s+"),QString::SkipEmptyParts);
                        for(int kk=0;kk<nobb;kk++)
                        {
                            Temp=list[kk];
                            stageBoundary[ii].v_bFactor[kk]=Temp.toDouble();
                        }

                    }

                    if(Line=="#CFACTOR")
                    {
                        Line=stream.readLine();
                        list=Line.split(QRegExp("\\s+"),QString::SkipEmptyParts);
                        for(int kk=0;kk<nobb;kk++)
                        {
                            Temp=list[kk];
                            stageBoundary[ii].v_cFactor[kk]=Temp.toDouble();
                        }

                    }

                    if(Line=="#X0")
                    {
                        Line=stream.readLine();
                        list=Line.split(QRegExp("\\s+"),QString::SkipEmptyParts);
                        for(int kk=0;kk<nobb;kk++)
                        {
                            Temp=list[kk];
                            stageBoundary[ii].v_x0[kk]=Temp.toDouble();
                        }

                    }

                    if(Line=="#X1")
                    {
                        Line=stream.readLine();
                        list=Line.split(QRegExp("\\s+"),QString::SkipEmptyParts);
                        for(int kk=0;kk<nobb;kk++)
                        {
                            Temp=list[kk];
                            stageBoundary[ii].v_x1[kk]=Temp.toDouble();
                        }
                    }

                    if(Line=="#Y0")
                    {
                        Line=stream.readLine();
                        list=Line.split(QRegExp("\\s+"),QString::SkipEmptyParts);
                        for(int kk=0;kk<nobb;kk++)
                        {
                            Temp=list[kk];
                            stageBoundary[ii].v_y0[kk]=Temp.toDouble();
                        }
                    }

                    if(Line=="#Y1")
                    {
                        Line=stream.readLine();
                        list=Line.split(QRegExp("\\s+"),QString::SkipEmptyParts);
                        for(int kk=0;kk<nobb;kk++)
                        {
                            Temp=list[kk];
                            stageBoundary[ii].v_y1[kk]=Temp.toDouble();
                        }
                        ii=ii+1;
                    }
                    Line=stream.readLine();
                }
            }
            if(Line=="#WATCHLIST")
            {
                Line=stream.readLine();
                list=Line.split(QRegExp("\\s+"),QString::SkipEmptyParts);
                Temp=list[1];
                watch.resize(Temp.toInt());
                int now=watch.size();
                int ii=0;
                while(ii<now)
                {
                    QString LineCheck, LineStop1, LineStop2;
                    LineCheck="#WATCH "+QString::number(ii+1,'f',0);
                    LineStop1="#WATCH "+QString::number(ii+2,'f',0);
                    if(Line==LineStop1)
                    {
                        break;
                    }
                    if(ii==now)
                    {
                        break;
                    }
                    if(stream.atEnd())
                    {
                        break;
                    }

                    if(Line==LineCheck)
                    {
                        watch[ii].watchIndex=ii+1;
                        Line=stream.readLine();
                        watch[ii].title=Line;
                        Line=stream.readLine();
                        list=Line.split(QRegExp("\\s+"),QString::SkipEmptyParts);
                        Temp=list[0];
                        watch[ii].watchType=Temp.toInt();

                        Temp=list[1];
                        watch[ii].x0=Temp.toDouble();

                        Temp=list[2];
                        watch[ii].x1=Temp.toDouble();

                        Temp=list[3];
                        watch[ii].y0=Temp.toDouble();

                        Temp=list[4];
                        watch[ii].y1=Temp.toDouble();

                        Temp=list[5];
                        watch[ii].beginStep=Temp.toInt();

                        Temp=list[6];
                        watch[ii].endStep=Temp.toInt();

                        if(list.size()>7)
                        {
                            Temp=list[7];
                            int averageBool=Temp.toInt();
                            if(averageBool>0)
                            {
                                watch[ii].averageBool=true;
                            }
                            else
                            {
                                watch[ii].averageBool=false;
                            }
                        }
                        ii=ii+1;
                    }
                    Line=stream.readLine();
                }
                ok=true;
            }

            if(ok==true)
            {
                break;
            }

        }
        mFile.close();
    }
    success=true;
    qDebug()<<"Import sucessfully"<<endl;
    //------------------

}

void SaveDataBase::pauseSystem()
{
    do {
        cout << '\n' << "Press the Enter key to continue.";
    } while (cin.get() != '\n');
}

void SaveDataBase::listMaterial()
{
    for (int i=0;i<material.size();i++)
    {
        cout<<"--------------"<<endl;
        cout<<"Material number : "<<material[i].matIndex<<endl;
        cout<<"Poission ratio  : "<<material[i].poission<<endl;
        cout<<"Ratio kh/kv     : "<<material[i].ratio<<endl;
        cout<<"Unit Weight     : "<<material[i].gf<<endl;
        cout<<"Void ratio      : "<<material[i].e0<<endl;
        cout<<"Cd factor       : "<<material[i].Cd<<endl;
        if(material[i].KFunction==0)
        {
            cout<<"Constanst K     : "<<material[i].KCurve(0,0)<<endl;
        }
        else
        {
            cout<<"K is function, imported from file :"<<endl;
            cout<<material[i].fileNameK.toStdString()<<endl;
            cout<<material[i].KCurve<<endl;
            cout<<"-----"<<endl;
        }

        if(material[i].kFunction==0)
        {
            cout<<"Constanst k     : "<<material[i].kCurve(0,0)<<endl;
        }
        else
        {
            cout<<"k is function"<<endl;
            cout<<material[i].fileNamek.toStdString()<<endl;
            cout<<material[i].kCurve<<endl;
            cout<<"-----"<<endl;
        }
    }
}

void SaveDataBase::listStage()
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
        cout<<"Star day          : "<<stage[i].t0<<" End day: "<<stage[i].t1<<endl;
        cout<<"Number of sub-step: "<<stage[i].subStep<<endl;
    }
}

void SaveDataBase::listBoundary()
{
    for (auto i=0;i<boundary.size();i++)
    {
        cout<<"----------------------"<<endl;
        cout<<"Boundary condition number: "<<boundary[i].boundaryIndex<<endl;
        cout<<"Boundary condition name  : "<<boundary[i].boundaryName.toStdString()<<endl;
        cout<<"Boundary condition type  : "<<boundary[i].boundaryType<<endl;
        cout<<"Load   type              : "<<boundary[i].loadType<<endl;
        cout<<"Constant or Start Value  : "<<boundary[i].val0<<endl;
        cout<<"End Value                : "<<boundary[i].val1<<endl;
    }
}

void SaveDataBase::listStageBoundary()
{
    for (int i=0;i<stageBoundary.size();i++)
    {
        for (int j=0;j<stageBoundary[i].v_stageBoundaryIndex.size();j++)
        {
            qDebug()<<"Stage: "<<stageBoundary[i].v_stageBoundaryIndex[j]<<endl;
            qDebug()<<"name: "<<stageBoundary[i].v_nodeFileName[j]<<endl;
            qDebug()<<"assingType: "<<stageBoundary[i].v_assignType[j]<<endl;
            qDebug()<<"boundary Index: "<<stageBoundary[i].v_boundaryIndex[j]<<endl;
            qDebug()<<"X0: "<<stageBoundary[i].v_x0[j]<<endl;
            qDebug()<<"X1: "<<stageBoundary[i].v_x1[j]<<endl;
            qDebug()<<"Y0: "<<stageBoundary[i].v_y0[j]<<endl;
            qDebug()<<"Y1: "<<stageBoundary[i].v_y1[j]<<endl;

        }
    }
}

void SaveDataBase::listWatch()
{
    for (int i=0;i<watch.size();i++)
    {
        cout<<"title: "<<endl;
        qDebug()<<watch[i].title<<endl;
    }
}

void SaveDataBase::sendSIGNAL()
{
    emit sendProjectSetting(projectParameters);
    emit sendMesh(coordinates,elements);
    emit sendMaterial(material);
    emit sendStage(stage);
    emit sendBoundary(boundary);
    emit sendStageBoundary(stageBoundary);
    emit sendWatchListBase(watch);
    emit sendGeometry(geometry);
    qDebug()<<"Data base is sent"<<endl;
}

void SaveDataBase::resetFolderName()
{
    saveClicked=false;
}

void SaveDataBase::getMesh(Ref<MatrixXd> coordinates, Ref<MatrixXd> elements, QString folderName)
{
    this->folderName=folderName;
    this->coordinates=coordinates;
    this->elements=elements;
    qDebug()<<"Mesh is sent"<<endl;
}

void SaveDataBase::getStageBoundary(vector<StageBoundaryBase> stageBoundary)
{
    this->stageBoundary=stageBoundary;
}

void SaveDataBase::getMaterial(vector<MaterialBase> material)
{
    this->material=material;
}

void SaveDataBase::getStage(vector<StageBase> stage)
{
    this->stage=stage;
}

void SaveDataBase::getBoundary(vector<BoundaryConditionBase> boundary)
{
    this->boundary=boundary;
}

void SaveDataBase::getProjectSetting(vector<double> projectParameters)
{
    this->projectParameters=projectParameters;
}

void SaveDataBase::getWatchListBase(vector<WatchListBase> watch)
{
    this->watch=watch;
}

void SaveDataBase::getGeometryBase(GeometryBase geometry)
{
    this->geometry=geometry;
}
