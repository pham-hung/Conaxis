#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    //onpen Console
    ui->setupUi(this);
    setUpLabel();
    setUpLayout();
    connect(replotButton,SIGNAL(clicked(bool)),this,SLOT(replotButton_clicked()));
    connect(colorBandButton,SIGNAL(clicked(bool)),this,SLOT(colorBandButton_clicked()));
    connect(colorBand,SIGNAL(sendColorBandInfor(int,int,int,int,int,int,int,int,int,int,int,QString)),
            this,SLOT(getColorBandInfor(int,int,int,int,int,int,int,int,int,int,int,QString)));
    connect(resultBox,SIGNAL(currentIndexChanged(int)),this,SLOT(changePlotData()));

    connect(boundary,SIGNAL(sendParameters(vector<BoundaryConditionBase>)),stageBoundary,SLOT(getParametersBoundary(vector<BoundaryConditionBase>)));
    connect(stage,SIGNAL(sendParameters(vector<StageBase>)),stageBoundary,SLOT(getParametersStage(vector<StageBase>)));
    connect(solution,SIGNAL(sendSoltuonParameters(OutputExportBase)),this,SLOT(getSolutionParameters(OutputExportBase)));
    connect(nextStep,SIGNAL(clicked(bool)),this,SLOT(nextStep_clicked()));
    connect(geo,SIGNAL(sendGeometryData(Ref<MatrixXd>,Ref<MatrixXd>,int,double)),this,SLOT(getMeshFromGeometry(Ref<MatrixXd>,Ref<MatrixXd>,int,double)));
    this->setWindowTitle("CONAXIS");

    initilizeData();
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::setUpLayout()
{
    this->setGeometry(200,200,1366,768);
    //    this->setFixedSize(1024,768);
    ui->centralWidget->setLayout(mainLayout);
    mainLayout->addLayout(glLayout,0,0,18,10);
    mainLayout->addLayout(buttonLayout,19,0,1,10);

    buttonLayout->addWidget(resultBox);
    buttonLayout->addWidget(replotButton);
    buttonLayout->addWidget(colorBandButton);
    buttonLayout->addWidget(nextStep);
    buttonLayout->addItem(spacer1);
}

void MainWindow::setUpLabel()
{    
    resultBox->addItem("Horizontal Displacement");
    resultBox->addItem("Vertical Displacement");
    resultBox->addItem("Excess Pore Pressure");
    resultBox->addItem("Material number");
    resultBox->addItem("Element Type");
    resultBox->addItem("Custom Plot");
    resultBox->addItem("Element Result Plot");

    replotButton->setText("Replot");
    colorBandButton->setText("View Setting");
    nextStep->setText("Next Step");

    resultBox->setFixedWidth(150);
    replotButton->setFixedWidth(70);
    colorBandButton->setFixedWidth(70);
    nextStep->setFixedWidth(70);
}

void MainWindow::getFolder()
{
    QString folderName_new=QFileDialog::getExistingDirectory();
    if(folderName_new=="")
    {
        return;
    }
    else
    {
        folderName=folderName_new;
        cout<<"Set working folder to: "<<folderName.toStdString()<<endl;
    }

}

void MainWindow::importData()
{
    GetFile getfile;
    timer.start();
    QObject::connect(threadU,SIGNAL(started()),getU,SLOT(DoGetFile()));
    QObject::connect(getU,SIGNAL(finishGetFile()),threadU,SLOT(quit()));
    QObject::connect(getU,SIGNAL(finishGetFile()),threadU,SLOT(deleteLater()));
    QObject::connect(threadU,SIGNAL(finished()),threadU,SLOT(deleteLater()));

    QObject::connect(threadV,SIGNAL(started()),getV,SLOT(DoGetFile()));
    QObject::connect(getV,SIGNAL(finishGetFile()),threadV,SLOT(quit()));
    QObject::connect(getV,SIGNAL(finishGetFile()),threadV,SLOT(deleteLater()));
    QObject::connect(threadV,SIGNAL(finished()),threadV,SLOT(deleteLater()));

    QObject::connect(threadP,SIGNAL(started()),getP,SLOT(DoGetFile()));
    QObject::connect(getP,SIGNAL(finishGetFile()),threadP,SLOT(quit()));
    QObject::connect(getP,SIGNAL(finishGetFile()),threadP,SLOT(deleteLater()));
    QObject::connect(threadP,SIGNAL(finished()),threadP,SLOT(deleteLater()));

    fileName=folderName+"/"+"coordinates.dat";
    file.setFileName(fileName);
    getfile.fileName=fileName;
    getfile.DoGetFile();
    coordinates=MatrixXd::Zero(getfile.row,getfile.col);
    coordinates=getfile.data_file;

    fileName=folderName+"/"+"elements.dat";
    file.setFileName(fileName);
    getfile.fileName=fileName;
    getfile.DoGetFile();
    elements=MatrixXd::Zero(getfile.row,getfile.col);
    elements=getfile.data_file;

    getU->fileName=folderName+"/"+"U.txt";
    getV->fileName=folderName+"/"+"V.txt";
    getP->fileName=folderName+"/"+"P.txt";

    //getU thread
    getU->moveToThread(threadU);
    getV->moveToThread(threadV);
    getP->moveToThread(threadP);

    threadU->start();
    threadV->start();
    threadP->start();

    bool checkThread;
    checkThread=(getU->checkFinish&&getV->checkFinish);
    checkThread=(checkThread&&getP->checkFinish);

    while(checkThread==false)
    {
        QThread::msleep(500);
        checkThread=(getU->checkFinish&&getV->checkFinish);
        checkThread=(checkThread&&getP->checkFinish);
        cout<<"Running ";
    }

    U=MatrixXd::Zero(getU->row,getU->col);
    U=getU->data_file;

    V=MatrixXd::Zero(getV->row,getV->col);
    V=getV->data_file;

    P=MatrixXd::Zero(getP->row,getP->col);
    P=getP->data_file;

    delete getU;
    delete getV;
    delete getP;

    XX=MatrixXd::Zero(0,0);
    coordScale=MatrixXd::Zero(coordinates.rows(),coordinates.cols());
    coord=MatrixXd::Zero(coordinates.rows(),coordinates.cols());
    cout<<"Import time: "<<timer.elapsed()/1000<<" seconds"<<endl;
    non=coordinates.rows();
    noe=elements.rows();
    ns=U.cols();

    loadDataClicked=true;
    resultBox->setCurrentIndex(-1);
}

void MainWindow::getUserData()
{    
    comboPosition=resultBox->currentIndex();
    resultTitle=resultBox->currentText();
    resultTitle=resultTitle+"-Step: "+QString::number(step,'f',0);
}

void MainWindow::systemPause()
{
    do {
        cout << '\n' << "Press the Enter key to continue.";
    } while (cin.get() != '\n');
}

void MainWindow::initilizeData()
{
    U=MatrixXd::Zero(1,1);
    V=MatrixXd::Zero(1,1);
    P=MatrixXd::Zero(1,1);
    coordinates=MatrixXd::Zero(0,0);
    elements=MatrixXd::Zero(0,0);
    ns=1;
    non=0;
    noe=0;
}

void MainWindow::closeEvent(QCloseEvent *event)
{
    QMessageBox::StandardButton resBtn = QMessageBox::question( this,"EXIT PROGRAM",
                                                                tr("Are you sure?\n"),
                                                                QMessageBox::Cancel | QMessageBox::No | QMessageBox::Yes,
                                                                QMessageBox::Yes);
    if (resBtn != QMessageBox::Yes) {
        event->ignore();
    } else {
        event->accept();
    }

    material->close();
    stage->close();
    boundary->close();
    stageBoundary->close();
}

void MainWindow::runSimulationCRS()
{
    AxisSymmetric_2D *test=new AxisSymmetric_2D;
    //Connect to another class
    connect(this,SIGNAL(sendMesh(Ref<MatrixXd>,Ref<MatrixXd>,QString)),test,SLOT(getMesh(Ref<MatrixXd>,Ref<MatrixXd>,QString)));
    connect(this,SIGNAL(sendAdditionalInformation(int,double,double)),test,SLOT(getAdditionalInformation(int,double,double)));
    connect(material,SIGNAL(sendMaterial(vector<MaterialBase>)),test,SLOT(getMaterial(vector<MaterialBase>)));
    connect(stage,SIGNAL(sendParameters(vector<StageBase>)),test,SLOT(getStage(vector<StageBase>)));
    connect(boundary,SIGNAL(sendParameters(vector<BoundaryConditionBase>)),test,SLOT(getBoundary(vector<BoundaryConditionBase>)));
    connect(stageBoundary,SIGNAL(sendParametersStageBoundary(vector<StageBoundaryBase>)),test,SLOT(getStageBoundary(vector<StageBoundaryBase>)));
    connect(project,SIGNAL(sendParameters(vector<double>)),test,SLOT(getProjectSetting(vector<double>)));
    connect(test,SIGNAL(sendSolverInformation(Ref<MatrixXd>,Ref<MatrixXd>,Ref<MatrixXd>)),this,SLOT(getSolverInformation(Ref<MatrixXd>,Ref<MatrixXd>,Ref<MatrixXd>)));
    connect(this,SIGNAL(sendFolder(QString)),test,SLOT(setFolder(QString)));
    connect(watch,SIGNAL(sendSignal(vector<WatchListBase>)),test,SLOT(getWatchList(vector<WatchListBase>)));
    connect(watch,SIGNAL(sendSIGNALNow(vector<WatchListBase>)),test,SLOT(getAndExportWatchList(vector<WatchListBase>)));
    connect(this,SIGNAL(sendCrsData(Ref<MatrixXd>,bool,int)),test,SLOT(getCrsData(Ref<MatrixXd>,bool,int)));

    emit this->sendMesh(coordinates,elements,folderName);
    emit this->sendAdditionalInformation(analysisType,A1D,k1D);
    emit this->sendCrsData(crsObject.crsData,true,crsObject.analysisType);
    material->sendSIGNAL();
    boundary->sendSignal();
    stage->sendSignal();
    stageBoundary->sendSIGNAL();
    project->sendSIGNAL();
    watch->sendSIGNAL();

    test->prepareData();
    if(backAnalysisFlag==false)
    {
        test->runCrsTest();
        test->exportWatchList();
    }
    else
    {
        test->runBackAnalysis();
    }
}

void MainWindow::getSignalFromClass()
{
    emit sendParas(minVal,maxVal,meshBool,nodeBool,resultBool,autoColor);
    emit sendDataToGL(coord,elements,XX,coordScale);
    emit sendViewportToWidget(xleft0,xright0,ybot0,ytop0,lockView);
    emit sendColorInfor(noc,numType,fontSize,colorPosition,nodePlot,meshPlot,resultPlot,axePlot,datePlot,titlePlot,valPlot,QString(title+" - "+resultTitle));
    disconnect(this,SIGNAL(sendParas(double,double,int,int,int,bool)),widget,SLOT(getPars(double,double,int,int,int,bool)));
    disconnect(this,SIGNAL(sendDataToGL(Ref<MatrixXd>,Ref<MatrixXd>,Ref<MatrixXd>,Ref<MatrixXd>)),widget,SLOT(getDataGL(Ref<MatrixXd>,Ref<MatrixXd>,Ref<MatrixXd>,Ref<MatrixXd>)));
    disconnect(widget,SIGNAL(sendSignaltoGUI()),this,SLOT(getSignalFromClass()));
    disconnect(this,SIGNAL(sendColorInfor(int,int,int,int,int,int,int,int,int,int,int,QString)),
               widget,SLOT(getColorBandInfor(int,int,int,int,int,int,int,int,int,int,int,QString)));
}

void MainWindow::getCurrentValSignal(vector<double> ScreenValue)
{
    this->ScreenValue=ScreenValue;
    QString shownValue="Value at mouse position: "+QString::number(ScreenValue[0]);
    double xCoordinates=ScreenValue[1]/xScale;
    double yCoordinates=ScreenValue[2]/yScale;
    QString xLocation= "X-coordinates: "+QString::number(xCoordinates,'e',5);
    QString yLocation= "Y-coordinates: "+QString::number(yCoordinates,'e',5);
    shownValue=shownValue+"         "+xLocation+"          "+yLocation;

    ui->statusBar->showMessage(shownValue);
    //    qDebug()<<"Value at mouse position: "<<currentVal<<endl;
}

void MainWindow::getViewPort(float xleft, float xright, float ybot, float ytop)
{
    if(lockView==false)
    {
        this->xleft=xleft;
        this->xright=xright;
        this->ybot=ybot;
        this->ytop=ytop;
    }
}

void MainWindow::getColorBandInfor(int noc,int numType,int fontSize,int colorPosition,
                                   int nodePlot,int meshPlot,int resultPlot,int axePlot,int datePlot,int titlePlot,int valPlot,QString title)
{
    this->noc=noc;
    this->numType=numType;
    this->fontSize=fontSize;
    this->colorPosition=colorPosition;
    this->nodePlot=nodePlot;
    this->meshPlot=meshPlot;
    this->resultPlot=resultPlot;
    this->axePlot=axePlot;
    this->datePlot=datePlot;
    this->titlePlot=titlePlot;
    this->title=title;
    this->valPlot=valPlot;
}

void MainWindow::getMesh(Ref<MatrixXd> coordinates, Ref<MatrixXd> elements)
{
    this->coordinates=coordinates;
    this->elements=elements;
}

void MainWindow::getMeshFromGeometry(Ref<MatrixXd> elements, Ref<MatrixXd> coordinates, int analysisType, double qw)
{
    this->coordinates=coordinates;
    this->elements=elements;
    this->analysisType=analysisType;
    this->qw=qw;

    rw=coordinates.col(1).minCoeff();
    A1D=PI*rw*rw;
    k1D=qw/(2*PI*A1D);

    XX=MatrixXd::Zero(0,0);
    coordScale=MatrixXd::Zero(coordinates.rows(),coordinates.cols());
    coord=MatrixXd::Zero(coordinates.rows(),coordinates.cols());
    non=coordinates.rows();
    noe=elements.rows();
    ns=U.cols();

    //    if(loadDataClicked==false)
    //    {
    //        colorFormSetup();
    //    }

    loadDataClicked=true;
    resultBox->setCurrentIndex(3); //Material plot
    replotButton_clicked();
}

void MainWindow::getSolverInformation(Ref<MatrixXd> U, Ref<MatrixXd> V, Ref<MatrixXd> P)
{
    this->U=U;
    this->V=V;
    this->P=P;
    ns=U.cols();
}

void MainWindow::getSolutionParameters(OutputExportBase solutionParameters)
{
    this->solutionParameters=solutionParameters;
}

void MainWindow::getCrsData(CRS_Base crsObject, bool backAnalysisFlag)
{
    this->crsObject=crsObject;
    this->backAnalysisFlag=backAnalysisFlag;
    geo->createCrsMesh(crsObject.H0,crsObject.R);
    material->createCrsMaterial(crsObject.possionRatio,crsObject.e0);

    MatrixXd timeStep=MatrixXd::Zero(crsObject.crsData.rows(),1);
    timeStep=crsObject.crsData.col(0);
    stage->createCrsStage(timeStep);

    MatrixXd testData=MatrixXd::Zero(crsObject.crsData.rows(),3);
    testData.col(0)=crsObject.crsData.col(0)/1440.0f;
    testData.col(1)=crsObject.crsData.col(1);
    testData.col(2)=crsObject.crsData.col(2);
    boundary->createCrsBoundary(testData);

    stageBoundary->createCrsStageBoundary(crsObject.H0,crsObject.R,crsObject.analysisType);
    watch->createCrsWatchList(crsObject.H0,crsObject.R);

    //Run
    runSimulationCRS();
}

void MainWindow::runAnimation()
{    
    if(endStep>U.cols())
    {
        endStep=U.cols();
    }
    for (int j=beginStep;j<=endStep;j++)
    {
        if (glLayout->count()!=0)
        {

            glLayout->removeWidget(widget);
            widget->setParent(NULL);
            delete widget;
            widget=NULL;
        }
        getUserData();
        step=j;
        XX.setZero();

        if(comboPosition==0)
        {
            XX.resize(U.rows(),1);
            XX.col(0)=U.col(step-1);
        }
        else if (comboPosition==1)
        {
            XX.resize(V.rows(),1);
            XX.col(0)=V.col(step-1);

        }
        else if(comboPosition==2)
        {
            XX.resize(P.rows(),1);
            XX.col(0)=P.col(step-1);
        }
        else if (comboPosition==3)
        {
            XX.resize(1,1);
            XX(0,0)=0;
        }
        else if(comboPosition==4)
        {
            XX.resize(1,1);
            XX(0,0)=1;
        }
        else if(comboPosition==5)
        {
            XX.resize(customPlot.rows(),1);
            XX.col(0)=customPlot.col(step-1);
        }
        else if(comboPosition==6)
        {
            XX.resize(elemPlot.rows(),1);
            XX.col(0)=elemPlot.col(step-1);
        }
        coord.col(0)=coordinates.col(0);
        coord.col(1)=xScale*coordinates.col(1);
        coord.col(2)=yScale*coordinates.col(2);
        coord.col(3)=coordinates.col(3);

        if(XX.rows()!=1)
        {
            coordScale.col(0)=coordinates.col(0);
            coordScale.col(1)=xScale*coordinates.col(1)+xScale*dScale*U.col(step-1);
            coordScale.col(2)=yScale*coordinates.col(2)+yScale*dScale*V.col(step-1);
            coordScale.col(3)=coordinates.col(3);
        }
        else
        {
            coordScale=coord;
        }
        //------------------
        widget =new GLWidget; //openGL Widget

        QSurfaceFormat format;
        format.setRenderableType(QSurfaceFormat::OpenGL);
        format.setProfile(QSurfaceFormat::OpenGLContextProfile::CompatibilityProfile);
        format.setVersion(3,3);
        widget->setFormat(format);

        connect(this,SIGNAL(sendParas(double,double,int,int,int,bool)),widget,SLOT(getPars(double,double,int,int,int,bool)),Qt::DirectConnection);
        connect(this,SIGNAL(sendDataToGL(Ref<MatrixXd>,Ref<MatrixXd>,Ref<MatrixXd>,Ref<MatrixXd>)),widget,SLOT(getDataGL(Ref<MatrixXd>,Ref<MatrixXd>,Ref<MatrixXd>,Ref<MatrixXd>)),Qt::DirectConnection);

        connect(widget,SIGNAL(sendSignaltoGUI()),this,SLOT(getSignalFromClass()));
        connect(widget,SIGNAL(sendCurrentValToGUI(vector<double>)),this,SLOT(getCurrentValSignal(vector<double>)));
        connect(widget,SIGNAL(sendViewport(float,float,float,float)),this,SLOT(getViewPort(float,float,float,float)));
        connect(this,SIGNAL(sendViewportToWidget(float,float,float,float,bool)),widget,SLOT(getViewportInfor(float,float,float,float,bool)));
        connect(this,SIGNAL(sendColorInfor(int,int,int,int,int,int,int,int,int,int,int,QString))
                ,widget,SLOT(getColorBandInfor(int,int,int,int,int,int,int,int,int,int,int,QString)));

        //glLayout
        glLayout->addWidget(widget);
        QCoreApplication::processEvents();
        QThread::msleep(delay);
    }
}

void MainWindow::autoMaxMinValueCalculate()
{
    if(loadDataClicked==false)
    {
        QMessageBox::warning(this,"ERROR: NODATA","Load Data first");
        return;
    }
    else
    {
        XX.setZero();
        if(comboPosition==0 && U.rows()==coordinates.rows()) //Horizontal
        {
            XX.resize(U.rows(),1);
            XX.col(0)=U.col(step-1);
            minVal=XX.col(0).minCoeff();
            maxVal=XX.col(0).maxCoeff();
        }
        else if(comboPosition==0 && U.rows()!=coordinates.rows())
        {
            resultBox->setCurrentIndex(3);
            QMessageBox::warning(this,"ERROR: NODATA","Load results data first");
        }
        else if (comboPosition==1 && V.rows()==coordinates.rows())
        {
            XX.resize(V.rows(),1);
            XX.col(0)=V.col(step-1);
            minVal=XX.col(0).minCoeff();
            maxVal=XX.col(0).maxCoeff();

        }
        else if (comboPosition==1 && V.rows()!=coordinates.rows())
        {
            resultBox->setCurrentIndex(3);
            QMessageBox::warning(this,"ERROR: NODATA","Load results data first");
        }
        else if(comboPosition==2 && P.rows()==coordinates.rows())
        {
            XX.resize(P.rows(),1);
            XX.col(0)=P.col(step-1);
            minVal=XX.col(0).minCoeff();
            maxVal=XX.col(0).maxCoeff();
        }
        else if(comboPosition==2 && P.rows()!=coordinates.rows())
        {
            resultBox->setCurrentIndex(3);
            QMessageBox::warning(this,"ERROR: NODATA","Load results data first");
        }
        else if (comboPosition==3)
        {
            minVal=elements.col(10).minCoeff();
            maxVal=elements.col(10).maxCoeff();
        }
        else if(comboPosition==4)
        {
            minVal=elements.col(11).minCoeff();
            maxVal=elements.col(11).maxCoeff();
        }
        else if(comboPosition==5)
        {
            XX.resize(customPlot.rows(),1);
            XX.col(0)=customPlot.col(step-1);
            minVal=XX.col(0).minCoeff();
            maxVal=XX.col(0).maxCoeff();
        }
        else if(comboPosition==6)
        {
            XX.resize(elemPlot.rows(),1);
            XX.col(0)=elemPlot.col(step-1);
            minVal=XX.col(0).minCoeff();
            maxVal=XX.col(0).maxCoeff();
        }
        resultTitle=resultBox->currentText();
    }
}

void MainWindow::loadButton_clicked()
{
    getFolder();
    if(folderName=="")
    {
        return;
    }
    else
    {
        importData();
        colorFormSetup();
    }
}

void MainWindow::replotButton_clicked()
{

    if(loadDataClicked==false)
    {
        QMessageBox::warning(this,"ERROR: NODATA","Load Data first");
        return;
    }
    else
    {
        //Delete old widget
        if (glLayout->count()!=0)
        {
            glLayout->removeWidget(widget);
            widget->setParent(NULL);
            delete widget;
            widget=NULL;
        }

        //Input data
        getUserData();
        XX.setZero();
        if(comboPosition<0)
        {
            resultBox->setCurrentIndex(3);
        }
        if(comboPosition==0 && U.rows()==coordinates.rows())
        {
            XX.resize(U.rows(),1);
            XX.col(0)=U.col(step-1);
        }
        else if(comboPosition==0 && U.rows()!=coordinates.rows())
        {
            resultBox->setCurrentIndex(3);
            XX.resize(1,1);
            XX(0,0)=0;
            QMessageBox::warning(this,"ERROR","Results are not valid");
        }
        else if (comboPosition==1 && V.rows()==coordinates.rows())
        {
            XX.resize(V.rows(),1);
            XX.col(0)=V.col(step-1);

        }
        else if(comboPosition==1 && V.rows()!=coordinates.rows())
        {
            resultBox->setCurrentIndex(3);
            XX.resize(1,1);
            XX(0,0)=0;
            QMessageBox::warning(this,"ERROR","Results are not valid");
        }
        else if(comboPosition==2 && P.rows()==coordinates.rows())
        {
            XX.resize(P.rows(),1);
            XX.col(0)=P.col(step-1);
        }
        else if(comboPosition==2 && P.rows()!=coordinates.rows())
        {
            resultBox->setCurrentIndex(3);
            XX.resize(1,1);
            XX(0,0)=0;
            QMessageBox::warning(this,"ERROR","Results are not valid");
        }
        else if (comboPosition==3)
        {
            XX.resize(1,1);
            XX(0,0)=0;
        }
        else if(comboPosition==4)
        {
            XX.resize(1,1);
            XX(0,0)=1;
        }
        else if(comboPosition==5)
        {
            XX.resize(customPlot.rows(),1);
            XX.col(0)=customPlot.col(step-1);
        }
        else if(comboPosition==6)
        {
            XX.resize(elemPlot.rows(),1);
            XX.col(0)=elemPlot.col(step-1);
        }
        coord.col(0)=coordinates.col(0);
        coord.col(1)=xScale*coordinates.col(1);
        coord.col(2)=yScale*coordinates.col(2);
        coord.col(3)=coordinates.col(3);

        if(XX.rows()!=1&&U.rows()==non)
        {
            coordScale.col(0)=coordinates.col(0);
            coordScale.col(1)=xScale*coordinates.col(1)+xScale*dScale*U.col(step-1);
            coordScale.col(2)=yScale*coordinates.col(2)+yScale*dScale*V.col(step-1);
            coordScale.col(3)=coordinates.col(3);
        }
        else
        {
            coordScale=coord;
        }

        //------------------
        widget =new GLWidget; //openGL Widget

        QSurfaceFormat format;
        format.setRenderableType(QSurfaceFormat::OpenGL);
        format.setProfile(QSurfaceFormat::OpenGLContextProfile::CompatibilityProfile);
        format.setVersion(3,3);
        widget->setFormat(format);

        connect(this,SIGNAL(sendParas(double,double,int,int,int,bool)),widget,SLOT(getPars(double,double,int,int,int,bool)),Qt::DirectConnection);
        connect(this,SIGNAL(sendDataToGL(Ref<MatrixXd>,Ref<MatrixXd>,Ref<MatrixXd>,Ref<MatrixXd>)),widget,SLOT(getDataGL(Ref<MatrixXd>,Ref<MatrixXd>,Ref<MatrixXd>,Ref<MatrixXd>)),Qt::DirectConnection);

        connect(widget,SIGNAL(sendSignaltoGUI()),this,SLOT(getSignalFromClass()));
        connect(widget,SIGNAL(sendCurrentValToGUI(vector<double>)),this,SLOT(getCurrentValSignal(vector<double>)));
        connect(widget,SIGNAL(sendViewport(float,float,float,float)),this,SLOT(getViewPort(float,float,float,float)));
        connect(this,SIGNAL(sendViewportToWidget(float,float,float,float,bool)),widget,SLOT(getViewportInfor(float,float,float,float,bool)));
        connect(this,SIGNAL(sendColorInfor(int,int,int,int,int,int,int,int,int,int,int,QString))
                ,widget,SLOT(getColorBandInfor(int,int,int,int,int,int,int,int,int,int,int,QString)));

        //glLayout
        glLayout->addWidget(widget);
        widget->update();
    }

}

void MainWindow::lockViewport_clicked()
{
    lockView=true;
    xleft0=xleft;
    xright0=xright;
    ybot0=ybot;
    ytop0=ytop;
    emit sendViewportToWidget(xleft0,xright0,ybot0,ytop0,lockView);
}

void MainWindow::unlockViewport_clicked()
{
    lockView=false;
    xleft0=0;
    xright0=0;
    ybot0=0;
    ytop0=0;
    emit sendViewportToWidget(xleft0,xright0,ybot0,ytop0,lockView);
}

void MainWindow::colorBandButton_clicked()
{
    colorBand->setCurrentValue();
    colorBand->show();
}

void MainWindow::changePlotData()
{
    if(loadDataClicked==false)
    {
        QMessageBox::warning(this,"ERROR: NODATA","Load Data first");
        return;
    }
    else
    {
        getUserData();
        if(resultBox->currentIndex()==5) //Get custom nodal results
        {
            QString fileNameCustom=QFileDialog::getOpenFileName();
            GetFile getfile;
            file.setFileName(fileNameCustom);
            getfile.fileName=fileNameCustom;
            getfile.DoGetFile();
            if(getfile.data_file.rows()!=1)
            {
                customPlot=getfile.data_file;
                ns=customPlot.cols();
                step=1;
                getUserData();
            }
            else
            {
                if(customPlot.rows()==1){resultBox->setCurrentIndex(3);}
                else{autoMaxMinValueCalculate();}
                return;
            }
        }
        else if(resultBox->currentIndex()==6) //Get custom element results
        {
            QString fileNameCustom=QFileDialog::getOpenFileName();
            GetFile getfile;
            file.setFileName(fileNameCustom);
            getfile.fileName=fileNameCustom;
            getfile.DoGetFile();
            if(getfile.data_file.rows()!=1)
            {
                elemPlot=getfile.data_file;
            }
            else
            {
                if(elemPlot.rows()==1){resultBox->setCurrentIndex(3);}
                else{autoMaxMinValueCalculate();}
                return;
            }
        }
        autoMaxMinValueCalculate();
        resultTitle=resultBox->currentText();
    }
}

void MainWindow::colorFormSetup()
{
    colorBand->show();
}

void MainWindow::on_createMat_triggered()
{
    material->show();
}

void MainWindow::on_addStage_triggered()
{
    stage->show();
}

void MainWindow::on_listLoadStage_triggered()
{
    stage->showData();
}

void MainWindow::on_projectSetting_triggered()
{
    project->show();
}

void MainWindow::on_addBoundary_triggered()
{
    boundary->show();
}

void MainWindow::on_assignBoundary_triggered()
{    
    emit boundary->sendSignal();
    emit stage->sendSignal();
    stageBoundary->resetData();
    if(stageBoundary->checkSuccess==true)
    {
        stageBoundary->show();
    }
    else
    {
        stageBoundary->checkSuccess=true;
    }

}

void MainWindow::on_exit_triggered()
{
    this->close();
}

void MainWindow::on_runAnalysis_triggered()
{
    bool ok=true;
    AxisSymmetric_2D *test=new AxisSymmetric_2D;
    //Connect to another class
    connect(this,SIGNAL(sendMesh(Ref<MatrixXd>,Ref<MatrixXd>,QString)),test,SLOT(getMesh(Ref<MatrixXd>,Ref<MatrixXd>,QString)));
    connect(this,SIGNAL(sendAdditionalInformation(int,double,double)),test,SLOT(getAdditionalInformation(int,double,double)));
    connect(material,SIGNAL(sendMaterial(vector<MaterialBase>)),test,SLOT(getMaterial(vector<MaterialBase>)));
    connect(stage,SIGNAL(sendParameters(vector<StageBase>)),test,SLOT(getStage(vector<StageBase>)));
    connect(boundary,SIGNAL(sendParameters(vector<BoundaryConditionBase>)),test,SLOT(getBoundary(vector<BoundaryConditionBase>)));
    connect(stageBoundary,SIGNAL(sendParametersStageBoundary(vector<StageBoundaryBase>)),test,SLOT(getStageBoundary(vector<StageBoundaryBase>)));
    connect(project,SIGNAL(sendParameters(vector<double>)),test,SLOT(getProjectSetting(vector<double>)));
    connect(test,SIGNAL(sendSolverInformation(Ref<MatrixXd>,Ref<MatrixXd>,Ref<MatrixXd>)),this,SLOT(getSolverInformation(Ref<MatrixXd>,Ref<MatrixXd>,Ref<MatrixXd>)));
    connect(this,SIGNAL(sendFolder(QString)),test,SLOT(setFolder(QString)));
    connect(watch,SIGNAL(sendSignal(vector<WatchListBase>)),test,SLOT(getWatchList(vector<WatchListBase>)));
    connect(watch,SIGNAL(sendSIGNALNow(vector<WatchListBase>)),test,SLOT(getAndExportWatchList(vector<WatchListBase>)));

    //Emit signal
    if(coordinates.rows()==1)
    {
        QMessageBox::warning(this,"ERROR","Mesh is not imported");
        return;
    }
    else
    {
        emit this->sendMesh(coordinates,elements,folderName);
        emit this->sendAdditionalInformation(analysisType,A1D,k1D);
        material->sendSIGNAL();
        boundary->sendSignal();
        stage->sendSignal();
        stageBoundary->sendSIGNAL();
        project->sendSIGNAL();
        watch->sendSIGNAL();
    }

    int NumberOfMaterial=elements.col(10).maxCoeff();
    int realMat=material->NumberOfMat();

    if(NumberOfMaterial>realMat)
    {
        ok=false;
        QString m_String="Material from" +QString::number(realMat+1,'f',0)+" are not defined";
        QMessageBox::critical(this,"ERROR",m_String);
    }

    if(ok==true)
    {
        test->prepareData();
        test->runAllAnalysis();
        if(solutionParameters.nodalSolutionCheck==1){test->exportResults();}
        if(solutionParameters.elemetStressCheck==1){test->exportStress();}
        if(solutionParameters.averageStressCheck==1){test->averageStress();}
        if(solutionParameters.materialParameterCheck==1){test->exportMaterialsParameters();}
        test->exportWatchList();
    }
}

void MainWindow::on_saveData_triggered()
{
    connect(project,SIGNAL(sendParameters(vector<double>)),save,SLOT(getProjectSetting(vector<double>)));
    connect(this,SIGNAL(sendMesh(Ref<MatrixXd>,Ref<MatrixXd>,QString)),save,SLOT(getMesh(Ref<MatrixXd>,Ref<MatrixXd>,QString)));
    connect(material,SIGNAL(sendMaterial(vector<MaterialBase>)),save,SLOT(getMaterial(vector<MaterialBase>)));
    connect(stage,SIGNAL(sendParameters(vector<StageBase>)),save,SLOT(getStage(vector<StageBase>)));
    connect(boundary,SIGNAL(sendParameters(vector<BoundaryConditionBase>)),save,SLOT(getBoundary(vector<BoundaryConditionBase>)));
    connect(stageBoundary,SIGNAL(sendParametersStageBoundary(vector<StageBoundaryBase>)),save,SLOT(getStageBoundary(vector<StageBoundaryBase>)));
    connect(watch,SIGNAL(sendSignal(vector<WatchListBase>)),save,SLOT(getWatchListBase(vector<WatchListBase>)));
    emit this->sendMesh(coordinates,elements,folderName);

    project->sendSIGNAL();
    material->sendSIGNAL();
    stage->sendSignal();
    boundary->sendSignal();
    stageBoundary->sendSIGNAL();
    watch->sendSIGNAL();
    save->saveModel();
}

void MainWindow::on_openFile_triggered()
{
    system("cls");
    save->readModel();
    if(save->checkReadingState()==true)
    {
        connect(save,SIGNAL(sendMesh(Ref<MatrixXd>,Ref<MatrixXd>)),this,SLOT(getMesh(Ref<MatrixXd>,Ref<MatrixXd>)));
        connect(save,SIGNAL(sendBoundary(vector<BoundaryConditionBase>)),boundary,SLOT(getBoundary(vector<BoundaryConditionBase>)));
        connect(save,SIGNAL(sendMaterial(vector<MaterialBase>)),material,SLOT(getMaterial(vector<MaterialBase>)));
        connect(save,SIGNAL(sendProjectSetting(vector<double>)),project,SLOT(getProjectParameters(vector<double>)));
        connect(save,SIGNAL(sendStage(vector<StageBase>)),stage,SLOT(getStage(vector<StageBase>)));
        connect(save,SIGNAL(sendStageBoundary(vector<StageBoundaryBase>)),stageBoundary,SLOT(getStageBoundary(vector<StageBoundaryBase>)));
        connect(save,SIGNAL(sendWatchListBase(vector<WatchListBase>)),watch,SLOT(getWatchListBase(vector<WatchListBase>)));

        save->sendSIGNAL();
        emit boundary->sendSignal();
        emit stage->sendSignal();

        XX=MatrixXd::Zero(0,0);
        coordScale=MatrixXd::Zero(coordinates.rows(),coordinates.cols());
        coord=MatrixXd::Zero(coordinates.rows(),coordinates.cols());
        non=coordinates.rows();
        noe=elements.rows();
        loadDataClicked=true;
        resultBox->setCurrentIndex(-1);
        colorBand->show();
    }
}

void MainWindow::getQuitAction(bool quit, bool noSave)
{
    this->quit=quit;
    this->noSave=noSave;
}

void MainWindow::on_nodeList_triggered()
{
    ns=P.cols();
    qDebug()<<"List of nodes---------------"<<endl;
    qDebug()<<"Maximum X      : "<<coordinates.col(1).maxCoeff()<<endl;
    qDebug()<<"Minimum X      : "<<coordinates.col(1).minCoeff()<<endl;
    qDebug()<<"Maximum Y      : "<<coordinates.col(2).maxCoeff()<<endl;
    qDebug()<<"Minimum Y      : "<<coordinates.col(2).minCoeff()<<endl;
    qDebug()<<"Number of nodes: "<<coordinates.rows()<<endl;
    qDebug()<<"List of elements---------------"<<endl;
    qDebug()<<"Number of element :"<<noe<<endl;
    qDebug()<<"Working folder    :"<<folderName<<endl;
    qDebug()<<"Number of step    :"<<ns<<endl;
}

void MainWindow::on_setWorkingFolder_triggered()
{
    folderName=QFileDialog::getExistingDirectory();
    emit sendFolder(folderName);
}

void MainWindow::on_solutionControl_triggered()
{
    solution->updateData();
    solution->show();
}

void MainWindow::on_loadResults_triggered()
{
    connect(loadSolution,SIGNAL(sendFileName(QString,QString,QString,bool)),this,SLOT(getResultName(QString,QString,QString,bool)));
    loadSolution->show();
}

void MainWindow::getResultName(QString UfileName, QString VfileName, QString PfileName,bool ok)
{
    qDebug()<<"Load result files "<<endl;
    if(ok==true)
    {
        this->UfileName=UfileName;
        this->VfileName=VfileName;
        this->PfileName=PfileName;

        GetFile getfile;
        timer.start();
        QObject::connect(threadU,SIGNAL(started()),getU,SLOT(DoGetFile()));
        QObject::connect(getU,SIGNAL(finishGetFile()),threadU,SLOT(quit()));
        QObject::connect(getU,SIGNAL(finishGetFile()),threadU,SLOT(deleteLater()));
        QObject::connect(threadU,SIGNAL(finished()),threadU,SLOT(deleteLater()));

        QObject::connect(threadV,SIGNAL(started()),getV,SLOT(DoGetFile()));
        QObject::connect(getV,SIGNAL(finishGetFile()),threadV,SLOT(quit()));
        QObject::connect(getV,SIGNAL(finishGetFile()),threadV,SLOT(deleteLater()));
        QObject::connect(threadV,SIGNAL(finished()),threadV,SLOT(deleteLater()));

        QObject::connect(threadP,SIGNAL(started()),getP,SLOT(DoGetFile()));
        QObject::connect(getP,SIGNAL(finishGetFile()),threadP,SLOT(quit()));
        QObject::connect(getP,SIGNAL(finishGetFile()),threadP,SLOT(deleteLater()));
        QObject::connect(threadP,SIGNAL(finished()),threadP,SLOT(deleteLater()));

        getU->fileName=this->UfileName;
        getV->fileName=this->VfileName;
        getP->fileName=this->PfileName;

        //getU thread
        getU->moveToThread(threadU);
        getV->moveToThread(threadV);
        getP->moveToThread(threadP);

        threadU->start();
        threadV->start();
        threadP->start();

        bool checkThread;
        checkThread=(getU->checkFinish&&getV->checkFinish);
        checkThread=(checkThread&&getP->checkFinish);

        while(checkThread==false)
        {
            QThread::msleep(500);
            checkThread=(getU->checkFinish&&getV->checkFinish);
            checkThread=(checkThread&&getP->checkFinish);
            cout<<"Importing Data ";
        }

        U=MatrixXd::Zero(getU->row,getU->col);
        U=getU->data_file;

        V=MatrixXd::Zero(getV->row,getV->col);
        V=getV->data_file;

        P=MatrixXd::Zero(getP->row,getP->col);
        P=getP->data_file;

        delete getU;
        delete getV;
        delete getP;
        loadSolution->close();
        ns=U.cols();
        loadDataClicked=true;
    }
}

void MainWindow::on_animationSetting_triggered()
{
    connect(animationControl,SIGNAL(sendSignal(AnimationBaseClasee,bool)),this,SLOT(getAnimationParameters(AnimationBaseClasee,bool)));
    animationControl->show();
}

void MainWindow::getAnimationParameters(AnimationBaseClasee animationObject,bool ok)
{
    if(ok==true)
    {
        animationControl->close();
        lockView=animationObject.lockView;
        if(lockView==true)
        {
            xleft0=xleft;
            xright0=xright;
            ybot0=ybot;
            ytop0=ytop;
        }
        else
        {

        }
        beginStep=animationObject.beginStep;
        endStep=animationObject.endStep;
        delay=animationObject.delayTime;
        autoColor=animationObject.autoChangeColor;
        disconnect(animationControl,SIGNAL(sendSignal(AnimationBaseClasee,bool)),this,SLOT(getAnimationParameters(AnimationBaseClasee,bool)));
        this->runAnimation();

        //Finish run anitmation
        lockView=false;
        xleft0=0;
        xright0=0;
        ybot0=0;
        ytop0=0;
    }

}

void MainWindow::on_scaleSetting_triggered()
{
    connect(scaleControl,SIGNAL(sendSignal(ScaleSettingBase)),this,SLOT(getScaleSetting(ScaleSettingBase)));
    scaleControl->show();
}

void MainWindow::getScaleSetting(ScaleSettingBase scale)
{
    disconnect(scaleControl,SIGNAL(sendSignal(ScaleSettingBase)),this,SLOT(getScaleSetting(ScaleSettingBase)));
    scaleControl->close();
    if(scale.step>U.cols())
    {
        QMessageBox::warning(this,"ERROR","The input step is greater than maximun step");
        step=1;
    }
    step=scale.step;
    dScale=scale.resultScale;
    xScale=scale.xScale;
    yScale=scale.yScale;
    minVal=scale.minVal;
    maxVal=scale.maxVal;
    lockView=scale.lockView;

    if(lockView==true)
    {
        xleft0=xleft;
        xright0=xright;
        ybot0=ybot;
        ytop0=ytop;
    }

    if(scale.autoMaxMinVal==1)
    {
        autoMaxMinValueCalculate();
    }
    comboPosition=resultBox->currentIndex();
    resultTitle=resultBox->currentText();
    this->replotButton_clicked();
}

void MainWindow::getFolderName(QString coorFileName, QString elementFileName)
{
    GetFile getfile;
    fileName=coorFileName;
    file.setFileName(fileName);
    getfile.fileName=fileName;
    getfile.DoGetFile();
    coordinates=MatrixXd::Zero(getfile.row,getfile.col);
    coordinates=getfile.data_file;

    fileName=elementFileName;
    file.setFileName(fileName);
    getfile.fileName=fileName;
    getfile.DoGetFile();
    elements=MatrixXd::Zero(getfile.row,getfile.col);
    elements=getfile.data_file;

    XX=MatrixXd::Zero(0,0);
    coordScale=MatrixXd::Zero(coordinates.rows(),coordinates.cols());
    coord=MatrixXd::Zero(coordinates.rows(),coordinates.cols());
    non=coordinates.rows();
    noe=elements.rows();
    ns=U.cols();

    colorFormSetup();
    loadDataClicked=true;
    disconnect(loadMesh,SIGNAL(sendSignal(QString,QString)),this,SLOT(getFolderName(QString,QString)));
}

void MainWindow::on_loadMeshFile_triggered()
{
    connect(loadMesh,SIGNAL(sendSignal(QString,QString)),this,SLOT(getFolderName(QString,QString)));
    loadMesh->show();
}

void MainWindow::nextStep_clicked()
{
    autoColor=true;
    step=step+1;
    ns=U.cols();
    if(step>ns){step=ns;}
    lockView=true;
    xleft0=xleft;
    xright0=xright;
    ybot0=ybot;
    ytop0=ytop;
    emit sendViewportToWidget(xleft0,xright0,ybot0,ytop0,lockView);
    replotButton_clicked();
    lockView=false;
    emit sendViewportToWidget(xleft0,xright0,ybot0,ytop0,lockView);
}

void MainWindow::on_watchList_triggered()
{
    watch->show();
}

void MainWindow::on_CRSTest_triggered()
{    
    connect(crs,SIGNAL(sendCrsData(CRS_Base,bool)),this,SLOT(getCrsData(CRS_Base,bool)));
    crs->show();
}

void MainWindow::on_PVDbackAnalysis_triggered()
{
    PVDBackAnalysis *pvd=new PVDBackAnalysis;
    AxisSymmetric_2D *test=new AxisSymmetric_2D;
    bool ok=true;

    //Connect to another class
    connect(this,SIGNAL(sendMesh(Ref<MatrixXd>,Ref<MatrixXd>,QString)),test,SLOT(getMesh(Ref<MatrixXd>,Ref<MatrixXd>,QString)));
    connect(material,SIGNAL(sendMaterial(vector<MaterialBase>)),test,SLOT(getMaterial(vector<MaterialBase>)));
    connect(stage,SIGNAL(sendParameters(vector<StageBase>)),test,SLOT(getStage(vector<StageBase>)));
    connect(boundary,SIGNAL(sendParameters(vector<BoundaryConditionBase>)),test,SLOT(getBoundary(vector<BoundaryConditionBase>)));
    connect(stageBoundary,SIGNAL(sendParametersStageBoundary(vector<StageBoundaryBase>)),test,SLOT(getStageBoundary(vector<StageBoundaryBase>)));
    connect(project,SIGNAL(sendParameters(vector<double>)),test,SLOT(getProjectSetting(vector<double>)));
    connect(test,SIGNAL(sendSolverInformation(Ref<MatrixXd>,Ref<MatrixXd>,Ref<MatrixXd>)),this,SLOT(getSolverInformation(Ref<MatrixXd>,Ref<MatrixXd>,Ref<MatrixXd>)));
    connect(this,SIGNAL(sendFolder(QString)),test,SLOT(setFolder(QString)));
    connect(watch,SIGNAL(sendSignal(vector<WatchListBase>)),test,SLOT(getWatchList(vector<WatchListBase>)));
    connect(watch,SIGNAL(sendSIGNALNow(vector<WatchListBase>)),test,SLOT(getAndExportWatchList(vector<WatchListBase>)));
    connect(pvd,SIGNAL(sendPVDsParameters(double,double,double,double,double,double,bool,bool)),
            test,SLOT(getPVDsParameters(double,double,double,double,double,double,bool,bool)));
    connect(test,SIGNAL(sendPVDResults(double,double)),pvd,SLOT(getResult(double,double)));

    //Emit signal
    if(coordinates.rows()==1)
    {
        QMessageBox::warning(this,"ERROR","Mesh is not imported");
        return;
    }
    else
    {
        emit this->sendMesh(coordinates,elements,folderName);
        material->sendSIGNAL();
        boundary->sendSignal();
        stage->sendSignal();
        stageBoundary->sendSIGNAL();
        project->sendSIGNAL();
        watch->sendSIGNAL();
    }

    int NumberOfMaterial=elements.col(10).maxCoeff();
    int realMat=material->NumberOfMat();

    if(NumberOfMaterial>realMat)
    {
        ok=false;
        QString m_String="Material from" +QString::number(realMat+1,'f',0)+" are not defined";
        QMessageBox::critical(this,"ERROR",m_String);
    }

    if(ok==true)
    {
        pvd->show();
    }
}

void MainWindow::on_actionInput_Geometry_triggered()
{
    geo->show();
}

void MainWindow::on_saveAs_triggered()
{
    connect(project,SIGNAL(sendParameters(vector<double>)),save,SLOT(getProjectSetting(vector<double>)));
    connect(this,SIGNAL(sendMesh(Ref<MatrixXd>,Ref<MatrixXd>,QString)),save,SLOT(getMesh(Ref<MatrixXd>,Ref<MatrixXd>,QString)));
    connect(material,SIGNAL(sendMaterial(vector<MaterialBase>)),save,SLOT(getMaterial(vector<MaterialBase>)));
    connect(stage,SIGNAL(sendParameters(vector<StageBase>)),save,SLOT(getStage(vector<StageBase>)));
    connect(boundary,SIGNAL(sendParameters(vector<BoundaryConditionBase>)),save,SLOT(getBoundary(vector<BoundaryConditionBase>)));
    connect(stageBoundary,SIGNAL(sendParametersStageBoundary(vector<StageBoundaryBase>)),save,SLOT(getStageBoundary(vector<StageBoundaryBase>)));
    connect(watch,SIGNAL(sendSignal(vector<WatchListBase>)),save,SLOT(getWatchListBase(vector<WatchListBase>)));
    emit this->sendMesh(coordinates,elements,folderName);

    project->sendSIGNAL();
    material->sendSIGNAL();
    stage->sendSignal();
    boundary->sendSignal();
    stageBoundary->sendSIGNAL();
    watch->sendSIGNAL();
    save->resetFolderName();
    save->saveModel();
}
