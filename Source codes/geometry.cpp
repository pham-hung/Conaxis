#include "geometry.h"
#include "ui_geometry.h"

Geometry::Geometry(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::Geometry)
{
    ui->setupUi(this);
}

Geometry::~Geometry()
{
    delete ui;
}

void Geometry::updateInformation()
{
    checkValid=false;
    re=ui->reLine->text().toDouble();
    rw=ui->rwLine->text().toDouble();
    rs=ui->rsLine->text().toDouble();

    length=ui->lengthLine->text().toDouble();
    qw=ui->qwLine->text().toDouble();
    numberOfLayer=ui->layerLine->text().toInt();
    surfaceElevation=ui->elevationLine->text().toDouble();
    numberOfElementSmear=ui->npsLine->text().toInt();
    numberOfElementSoil=ui->nplLine->text().toInt();
    analysisType=ui->analysisBox->currentIndex();
    defaultSubLayer=ui->subLayerLine->text().toInt();

    checkValid=true;
    ui->inforTable->setRowCount(numberOfLayer);
    ui->inforTable->setColumnCount(4);
    int col=4;
    layerInfo.resize(numberOfLayer,col);

    for (int i=0;i<numberOfLayer;i++)
    {
        for (int j=0;j<col;j++)
        {
            QTableWidgetItem *item=ui->inforTable->item(i,j);
            double value;
            if(item==NULL)
            {
                value=1;
                QTableWidgetItem *newItem=new QTableWidgetItem("1");
                if(j==0)
                {
                    newItem->setText(QString::number(i+1));
                    value=i+1;
                }
                if(j==3)
                {
                    newItem->setText(QString::number(defaultSubLayer));
                    value=defaultSubLayer;
                }
                ui->inforTable->setItem(i,j,newItem);
            }
            else
            {
                value=item->text().toDouble();
            }
            layerInfo(i,j)=value;
        }
    }
}

void Geometry::createMeshSmear()
{
    //With smearzones
    int noe,non, numberOfSubLayer;
    numberOfSubLayer=0;

    for (int i=0;i<layerInfo.rows();i++)
    {
        numberOfSubLayer=layerInfo(i,3)+numberOfSubLayer;
    }
    noe=numberOfElementSmear+numberOfElementSoil;
    noe=noe*numberOfSubLayer;
    elements.resize(noe,12);

    XYCoordinates.clear();
    XYCoordinates.reserve(8*noe);

    int eleCount=0;
    double x1,x2,x3,x4,x5,x6,x7,x8;
    double y1,y2,y3,y4,y5,y6,y7,y8;
    int node1,node2,node3,node4,node5,node6,node7,node8;
    double layerElevation=surfaceElevation;

    for (int i=0;i<numberOfLayer;i++) //Loop over layers
    {
        int subLayer=layerInfo(i,3);
        double thickness=layerInfo(i,1);
        double dy=thickness/subLayer;

        for (int j=0;j<subLayer;j++) //Loop over sublayers
        {
            double dxSmear=(rs-rw)/numberOfElementSmear;
            double dxSoil=(re-rs)/numberOfElementSoil;
            double yTop=layerElevation-j*dy;
            double yBot=yTop-dy;

            for (int ii=0;ii<numberOfElementSmear;ii++) //Loop over smear element
            {
                elements(eleCount,0)=eleCount+1;
                elements(eleCount,9)=8;
                elements(eleCount,10)=layerInfo(i,2); //material
                elements(eleCount,11)=2; //Type smear

                //From node1-node8
                x1=rw+ii*dxSmear;
                x3=x1+dxSmear;
                x2=x1;
                x5=x1;
                x4=x3;
                x7=x3;
                x6=0.5*(x2+x3);
                x8=0.5*(x1+x4);
                y1=yTop;
                y4=y1;
                y8=y1;
                y2=yBot;
                y3=y2;
                y6=y2;
                y5=0.5*(y1+y2);
                y7=0.5*(y3+y4); 

                elements(eleCount,1)=findNodeIndex(x2,y2,XYCoordinates);
                elements(eleCount,2)=findNodeIndex(x3,y3,XYCoordinates);
                elements(eleCount,3)=findNodeIndex(x4,y4,XYCoordinates);
                elements(eleCount,4)=findNodeIndex(x1,y1,XYCoordinates);
                elements(eleCount,5)=findNodeIndex(x6,y6,XYCoordinates);
                elements(eleCount,6)=findNodeIndex(x7,y7,XYCoordinates);
                elements(eleCount,7)=findNodeIndex(x8,y8,XYCoordinates);
                elements(eleCount,8)=findNodeIndex(x5,y5,XYCoordinates);

                eleCount++;
            }

            for (int ii=0;ii<numberOfElementSoil;ii++) //Loop over soi zone
            {
                elements(eleCount,0)=eleCount+1;
                elements(eleCount,9)=8;
                elements(eleCount,10)=layerInfo(i,2); //material
                elements(eleCount,11)=1; //Type soil

                //From node1-node8
                x1=rs+ii*dxSoil;
                x3=x1+dxSoil;
                x2=x1;
                x5=x1;
                x4=x3;
                x7=x3;
                x6=0.5*(x2+x3);
                x8=0.5*(x1+x4);
                y1=yTop;
                y4=y1;
                y8=y1;
                y2=yBot;
                y3=y2;
                y6=y2;
                y5=0.5*(y1+y2);
                y7=0.5*(y3+y4);

                elements(eleCount,1)=findNodeIndex(x2,y2,XYCoordinates);
                elements(eleCount,2)=findNodeIndex(x3,y3,XYCoordinates);
                elements(eleCount,3)=findNodeIndex(x4,y4,XYCoordinates);
                elements(eleCount,4)=findNodeIndex(x1,y1,XYCoordinates);
                elements(eleCount,5)=findNodeIndex(x6,y6,XYCoordinates);
                elements(eleCount,6)=findNodeIndex(x7,y7,XYCoordinates);
                elements(eleCount,7)=findNodeIndex(x8,y8,XYCoordinates);
                elements(eleCount,8)=findNodeIndex(x5,y5,XYCoordinates);

                eleCount++;
            }
        }
        layerElevation=layerElevation-thickness;
    }

    coordinates.resize(XYCoordinates.size(),4);
    coordinates.setZero();
    for (int i=0;i<XYCoordinates.size();i++)
    {
        coordinates(i,0)=i+1;
        double value1,value2;
        value1=XYCoordinates[i].first;
        value2=XYCoordinates[i].second;
        double error=1e-12;
        if(fabs (value1) <error)
        {
            value1=0;
        }
        if(fabs (value2) <error)
        {
            value2=0;
        }
        coordinates(i,1)=value1;
        coordinates(i,2)=value2;
    }

    //create 1D elements
    if(analysisType==3)
    {
        node1D.resize(0);
        for (int i=0;i<numberOfSubLayer;i++)
        {
            int eleIndex=1+i*(numberOfElementSmear+numberOfElementSoil);
            int node1=elements(eleIndex-1,4);
            int node2=elements(eleIndex-1,1);
            double y1=coordinates(node1-1,2);
            double y2=coordinates(node2-1,2);

            if((-y1+surfaceElevation)>length || (-y2+surfaceElevation)>length)
            {
                break;
            }
            else
            {
                auto pair=std::make_pair(node1,node2);
                node1D.push_back(pair);
            }
        }
        //Add 1D elements to elemnts
        int numberOf1DElements=node1D.size();
        int elementsSize=elements.rows();
        elements.conservativeResize(elementsSize+numberOf1DElements,elements.cols());

        for (int j=0;j<numberOf1DElements;j++)
        {
            elements(elementsSize+j,0)=elementsSize+j+1;
            elements(elementsSize+j,1)=node1D[j].first;
            elements(elementsSize+j,2)=node1D[j].second;
            elements(elementsSize+j,9)=2;
            elements(elementsSize+j,10)=1;
            elements(elementsSize+j,11)=1;

            elements(elementsSize+j,3)=0;
            elements(elementsSize+j,4)=0;
            elements(elementsSize+j,5)=0;
            elements(elementsSize+j,6)=0;
            elements(elementsSize+j,7)=0;
            elements(elementsSize+j,8)=0;
        }
    }
}

void Geometry::createMeshNoSmear()
{
    //No smear zone
    int noe,non, numberOfSubLayer;
    numberOfSubLayer=0;
    for (int i=0;i<layerInfo.rows();i++)
    {
        numberOfSubLayer=layerInfo(i,3)+numberOfSubLayer;
    }
    noe=numberOfElementSoil;
    noe=noe*numberOfSubLayer;
    elements.resize(noe,12);

    XYCoordinates.clear();
    XYCoordinates.reserve(8*noe);

    int eleCount=0;
    double x1,x2,x3,x4,x5,x6,x7,x8;
    double y1,y2,y3,y4,y5,y6,y7,y8;
    int node1,node2,node3,node4,node5,node6,node7,node8;
    double layerElevation=surfaceElevation;

    for (int i=0;i<numberOfLayer;i++) //Loop over layers
    {
        int subLayer=layerInfo(i,3);
        double thickness=layerInfo(i,1);
        double dy=thickness/subLayer;
        for (int j=0;j<subLayer;j++) //Loop over sublayers
        {
            double dxSoil=(re-rw)/numberOfElementSoil;
            double yTop=layerElevation-j*dy;
            double yBot=yTop-dy;
            for (int ii=0;ii<numberOfElementSoil;ii++) //Loop over soi zone
            {
                elements(eleCount,0)=eleCount+1;
                elements(eleCount,9)=8;
                elements(eleCount,10)=layerInfo(i,2); //material
                elements(eleCount,11)=1; //Type soil

                //From node1-node8
                x1=rw+ii*dxSoil;
                x3=x1+dxSoil;
                x2=x1;
                x5=x1;
                x4=x3;
                x7=x3;
                x6=0.5*(x2+x3);
                x8=0.5*(x1+x4);
                y1=yTop;
                y4=y1;
                y8=y1;
                y2=yBot;
                y3=y2;
                y6=y2;
                y5=0.5*(y1+y2);
                y7=0.5*(y3+y4);

                elements(eleCount,1)=findNodeIndex(x1,y1,XYCoordinates);
                elements(eleCount,2)=findNodeIndex(x2,y2,XYCoordinates);
                elements(eleCount,3)=findNodeIndex(x3,y3,XYCoordinates);
                elements(eleCount,4)=findNodeIndex(x4,y4,XYCoordinates);
                elements(eleCount,5)=findNodeIndex(x5,y5,XYCoordinates);
                elements(eleCount,6)=findNodeIndex(x6,y6,XYCoordinates);
                elements(eleCount,7)=findNodeIndex(x7,y7,XYCoordinates);
                elements(eleCount,8)=findNodeIndex(x8,y8,XYCoordinates);

                eleCount++;
            }

        }
        layerElevation=layerElevation-thickness;
    }

    //convert pair vector to MatrixXd
    coordinates.resize(XYCoordinates.size(),4);
    coordinates.setZero();
    for (int i=0;i<XYCoordinates.size();i++)
    {
        coordinates(i,0)=i+1;
        double value1,value2;
        value1=XYCoordinates[i].first;
        value2=XYCoordinates[i].second;
        double error=1e-12;
        if(fabs (value1) <error)
        {
            value1=0;
        }
        if(fabs (value2) <error)
        {
            value2=0;
        }
        coordinates(i,1)=value1;
        coordinates(i,2)=value2;
    }

    //create 1D elements
    if(analysisType==1)
    {
        node1D.resize(0);
        for (int i=0;i<numberOfSubLayer;i++)
        {
            int eleIndex=1+i*numberOfElementSoil;
            int node1=elements(eleIndex-1,4);
            int node2=elements(eleIndex-1,1);
            double y1=coordinates(node1-1,2);
            double y2=coordinates(node2-1,2);

            if((-y1+surfaceElevation)>length || (-y2+surfaceElevation)>length)
            {
                break;
            }
            else
            {
                auto pair=std::make_pair(node1,node2);
                node1D.push_back(pair);
            }
        }
        //Add 1D elements to elemnts
        int numberOf1DElements=node1D.size();
        int elementsSize=elements.rows();
        elements.conservativeResize(elementsSize+numberOf1DElements,elements.cols());

        for (int j=0;j<numberOf1DElements;j++)
        {
            elements(elementsSize+j,0)=elementsSize+j+1;
            elements(elementsSize+j,1)=node1D[j].first;
            elements(elementsSize+j,2)=node1D[j].second;
            elements(elementsSize+j,9)=2;
            elements(elementsSize+j,10)=1;
            elements(elementsSize+j,11)=1;

            elements(elementsSize+j,3)=0;
            elements(elementsSize+j,4)=0;
            elements(elementsSize+j,5)=0;
            elements(elementsSize+j,6)=0;
            elements(elementsSize+j,7)=0;
            elements(elementsSize+j,8)=0;
        }
    }
    //Without smearzones
}

bool Geometry::checkInformation()
{
    return true;
}

bool Geometry::isFound(pair<double, double> findPair, pair<double, double> soucePair)
{
    double error=1e-10;
    if(fabs(findPair.first-soucePair.first)<error && fabs(findPair.second-soucePair.second)<error)
    {
        return true;
    }
    else
    {
        return false;
    }
}

int Geometry::findNodeIndex(double x, double y, vector<pair<double, double> > &XYCoordinates)
{
    vector<pair<double,double> >::iterator it;
    pair<double,double> p=std::make_pair(x,y);
    it=std::find_if(XYCoordinates.begin(),XYCoordinates.end(),[=](pair<double,double> sourePair)->bool
    {
        return isFound(p,sourePair);
    }
    );

    if(it != XYCoordinates.end())
    {
        int index=distance(XYCoordinates.begin(),it);
        return index+1;
        //return node (already created)
    }
    else
    {
        XYCoordinates.push_back(p);
        return XYCoordinates.size();
        //Return new node index
    }
}

void Geometry::pauseSystem()
{
    do {
        cout << '\n' << "Press the Enter key to continue.";
    } while (cin.get() != '\n');
}

void Geometry::shownPair(vector<pair<double, double> > &XYCoordinates)
{
    for (int i=0;i<XYCoordinates.size();i++)
    {
        qDebug()<<XYCoordinates[i].first<<'\t'<<XYCoordinates[i].second<<endl;
    }
}

void Geometry::createCrsMesh(double H0, double R)
{
    re=0;
    rw=R;
    rs=0;
    numberOfLayer=1;
    surfaceElevation=H0;
    analysisType=0;
    numberOfElementSmear=0;
    numberOfElementSoil=10;

    layerInfo.resize(1,4);
    layerInfo<<1,H0,1,10;

    createMeshNoSmear();
    emit sendGeometryData(elements,coordinates,analysisType,qw);
}

void Geometry::on_updateButton_clicked()
{
    updateInformation();
}

void Geometry::on_meshButton_clicked()
{
    updateInformation();
    if(analysisType==0 || analysisType == 1)
    {
        createMeshNoSmear();
    }
    else
    {
        createMeshSmear();
    }
    QString fileName="C:/elements.dat";
    exportFile.fileName=fileName.toStdString();
    exportFile.ToFile(elements);

    fileName="C:/coordinates.dat";
    exportFile.fileName=fileName.toStdString();
    exportFile.ToFile(coordinates);

    emit sendGeometryData(elements,coordinates,analysisType,qw);

}

void Geometry::on_analysisBox_currentIndexChanged(int index)
{
    if(index==0)
    {
        ui->rsLine->setDisabled(true);
        ui->rsLine->setText("NO SMEAR");
        ui->qwLine->setDisabled(true);
        ui->qwLine->setText("NO WELL RESISTANCE");
    }
    else if(index==1)
    {
        ui->rsLine->setDisabled(true);
        ui->rsLine->setText("NO SMEAR");
        ui->qwLine->setEnabled(true);
        ui->qwLine->setText("1E-5");
    }
    else if(index==2)
    {
        ui->rsLine->setEnabled(true);
        ui->rsLine->setText("");
        ui->qwLine->setDisabled(true);
        ui->qwLine->setText("NO WELL RESISTANCE");
    }
    else if(index==3)
    {
        ui->rsLine->setEnabled(true);
        ui->rsLine->setText("0.1");
        ui->qwLine->setEnabled(true);
        ui->qwLine->setText("1E-5");
    }
}

void Geometry::on_setSubLayer_clicked()
{
    defaultSubLayer=ui->subLayerLine->text().toInt();
    layerInfo.col(3).fill(defaultSubLayer);
    int col=layerInfo.cols();

    for (int i=0;i<numberOfLayer;i++)
    {
        for (int j=0;j<col;j++)
        {
            QTableWidgetItem *item=ui->inforTable->item(i,j);
            double value;
            if(item==NULL)
            {
                value=1;
                QTableWidgetItem *newItem=new QTableWidgetItem("1");
                if(j==0)
                {
                    newItem->setText(QString::number(i+1));
                    value=i+1;
                }
                if(j==3)
                {
                    newItem->setText(QString::number(defaultSubLayer));
                    value=defaultSubLayer;
                }
                ui->inforTable->setItem(i,j,newItem);
            }
            else
            {
                if(j==3)
                {
                    value=layerInfo(i,j);
                    QTableWidgetItem *newItem=new QTableWidgetItem(QString::number(value));
                    ui->inforTable->setItem(i,j,newItem);
                }
            }
        }
    }
}

void Geometry::on_replotButton_clicked()
{
    this->close();
}
