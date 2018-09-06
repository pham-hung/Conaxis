#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QBoxLayout>
#include <QPushButton>
#include <QLineEdit>
#include <QLabel>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QCheckBox>
#include <QFormLayout>
#include <QTextEdit>
#include <QFile>
#include <QFileDialog>
#include <QElapsedTimer>
#include <QMessageBox>
#include <QComboBox>
#include <QRgb>
#include <QMouseEvent>
#include <QFrame>
#include <QThread>
#include <QProcess>
#include <QRadioButton>
#include <QSpacerItem>
#include <QDebug>

#include <Eigen/Dense>
#include <iostream>
#include <windows.h>

#include "colorform.h"
#include "glwidget.h"
#include "getfile.h"
#include "loadmesh.h"
#include "crsbackanalysis.h"
#include "pvdbackanalysis.h"
#include "geometry.h"

#include "material.h"
#include "stage.h"
#include "projectsetting.h"
#include "boundarycondition.h"
#include "assignboundarycondition.h"
#include "axissymmetric_2d.h"
#include "BaseClass/savedatabase.h"
#include "solutioncontrol.h"
#include "loadresult.h"
#include "animation.h"
#include "BaseClass/animationbaseclasee.h"
#include "scalesetting.h"
#include "watchlist.h"
#include "BaseClass/watchlistbase.h"
#include "BaseClass/crs_base.h"

#define PI 3.14159265359

using namespace std;
using namespace Eigen;

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
    void setUpLayout();
    void setUpLabel();
    void getFolder();
    void importData();
    void getUserData();
    void systemPause();
    void initilizeData();
    void closeEvent(QCloseEvent *event);

    void runSimulationCRS();

public slots:
    void getSignalFromClass();
    void getCurrentValSignal(vector<double> ScreenValue);
    void getViewPort(float xleft,float xright,float ybot,float ytop);
    void getColorBandInfor(int noc,int numType,int fontSize,int colorPosition,
                           int nodePlot,int meshPlot,int resultPlot,int axePlot,int datePlot,int titlePlot,int valPlot,QString title);
    void getMesh(Ref<MatrixXd> coordinates, Ref<MatrixXd> elements);
    void getMeshFromGeometry(Ref<MatrixXd> elements, Ref<MatrixXd> coordinates,int analysisType,double qw);
    void getSolverInformation(Ref<MatrixXd> U, Ref<MatrixXd> V, Ref<MatrixXd> P);
    void getSolutionParameters(OutputExportBase solutionParameters);
    void getCrsData(CRS_Base crsObject,bool backAnalysisFlag);
    void runAnimation();
    void autoMaxMinValueCalculate();

private slots:
    void loadButton_clicked();
    void replotButton_clicked(); //replot
    void lockViewport_clicked(); //Lock view
    void unlockViewport_clicked(); //Unlock view
    void colorBandButton_clicked();  //Color Band infor
    void changePlotData(); //Change plot Data
    void colorFormSetup(); //Color Form Setup
    void on_createMat_triggered(); //create material
    void on_addStage_triggered();  //add stage analysis
    void on_listLoadStage_triggered();  //list stage
    void on_projectSetting_triggered(); //project setting
    void on_addBoundary_triggered();  //add boundary condition
    void on_assignBoundary_triggered();  //assign boundary conditions
    void on_exit_triggered();  //exits
    void on_runAnalysis_triggered(); //Run Analysis
    void on_saveAs_triggered(); //save as
    void on_saveData_triggered(); //save
    void on_openFile_triggered(); //open file
    void getQuitAction(bool quit,bool noSave);
    void on_nodeList_triggered(); //list nodes
    void on_setWorkingFolder_triggered(); //set working folder
    void on_solutionControl_triggered();  //solutions control
    void on_loadResults_triggered();      //load results
    void getResultName(QString UfileName, QString VfileName,QString PfileName,bool ok);
    void on_animationSetting_triggered();  //animation
    void getAnimationParameters(AnimationBaseClasee animationObject,bool ok);
    void on_scaleSetting_triggered();  //scale setting
    void getScaleSetting(ScaleSettingBase scale);
    void getFolderName(QString coorFileName, QString elementFileName); //folderName
    void on_loadMeshFile_triggered(); //load mesh
    void nextStep_clicked();  //next step
    void on_watchList_triggered(); //watch list
    void on_CRSTest_triggered();   //calculate CRS parameters
    void on_PVDbackAnalysis_triggered();    
    void on_actionInput_Geometry_triggered();

signals:
    void sendParas(double,double,int,int,int,bool); //minVal, maxVal, meshBool, nodeBool, resultBool;
    void sendFolder(QString);
    void sendDataToGL(Ref<MatrixXd>, Ref<MatrixXd>, Ref<MatrixXd>, Ref<MatrixXd>); //send coordinates, elements, results;
    void sendViewportToWidget(float,float,float,float,bool);
    void sendColorInfor(int noc,int numType,int fontSize,int colorPosition,
                       int nodePlot,int meshPlot,int resultPlot,int axePlot,int datePlot,int titlePlot,int valPlot,QString title); //send plot information
    void sendMesh(Ref<MatrixXd> coordinates, Ref<MatrixXd> elements,QString folderName);  //send Mesh information
    void sendAdditionalInformation(int analysisType,double A1D, double k1D); //send 1D parameters
    void sendCrsData(Ref<MatrixXd> crsData, bool crsFlag, int crsType);

protected:

private:
    Ui::MainWindow *ui;
    int non,noe,ns,step=1;
    double minVal, maxVal;
    double xScale=1, yScale=1, dScale=1;
    int meshBool, nodeBool, resultBool;
    bool autoColorCheck;
    bool loadDataClicked=false;
    bool lockView=false;
    bool autoColor;
    bool quit, noSave;

    QString folderName, fileName;
    QFile file;

    MatrixXd coordinates, elements, U, V, P;
    MatrixXd customPlot=MatrixXd::Zero(1,1);
    MatrixXd elemPlot=MatrixXd::Zero(1,1);
    MatrixXd XX;
    MatrixXd coordScale, coord;
    QElapsedTimer timer;
    QSurfaceFormat format;
    int comboPosition=1;

    float xleft, xright, ybot, ytop;
    float xleft0=0;
    float xright0=0;
    float ybot0=0;
    float ytop0=0;
    vector<double> ScreenValue;

    int beginStep, endStep;
    float delay;
    int currentStage;    
    int analysisType;
    double qw, A1D, k1D, rw;
    //------------------------------------------------------
    QGridLayout *mainLayout=new QGridLayout;
    QVBoxLayout *glLayout = new QVBoxLayout;
    QBoxLayout *buttonLayout= new QBoxLayout(QBoxLayout::LeftToRight);

    //QPushbutton
    QPushButton *replotButton= new QPushButton;
    QPushButton *colorBandButton= new QPushButton;
    QPushButton *nextStep= new QPushButton;
    QPushButton *previousStep= new QPushButton;

    //QComboBox
    QComboBox *resultBox =new QComboBox;

    //Spacer
    QSpacerItem *spacer1 = new QSpacerItem(20,10,QSizePolicy::Expanding,QSizePolicy::Minimum);

    //Pixmap
    QPixmap *pixmap= new QPixmap;

    //GLwidget
    GLWidget *widget;
    //-------------------------------------------------------------------
    //Multi Thread variables
    GetFile *getU=new GetFile;
    GetFile *getV=new GetFile;
    GetFile *getP=new GetFile;
    QThread *threadU=new QThread;
    QThread *threadV=new QThread;
    QThread *threadP=new QThread;
    QString UfileName,VfileName,PfileName;
    //-------------------------------------------------------------------
    //Color Form
    colorForm *colorBand=new colorForm(this);
    int noc=10;
    int numType=1;
    int fontSize=10;
    int colorPosition=0;    
    int nodePlot=1;
    int meshPlot=1;
    int resultPlot=1;
    int axePlot=1;
    int datePlot=1;
    int titlePlot=1;
    int valPlot=1;
    QString title="PROJECT";
    QString resultTitle;
    QString fullTitle;
    //-------------------------------------------------------------------
    Material *material=new Material;
    Stage *stage=new Stage;
    ProjectSetting *project=new ProjectSetting;
    BoundaryCondition *boundary=new BoundaryCondition;
    AssignBoundaryCondition *stageBoundary=new AssignBoundaryCondition;
    SaveDataBase *save=new SaveDataBase;
    SolutionControl *solution=new SolutionControl;
    OutputExportBase solutionParameters;    
    LoadResult *loadSolution=new LoadResult;
    AnimationBaseClasee *animationObject=new AnimationBaseClasee;
    Animation *animationControl=new Animation;
    ScaleSetting *scaleControl=new ScaleSetting;
    LoadMesh *loadMesh=new LoadMesh;
    WatchList *watch =new WatchList;
    Geometry *geo=new Geometry;    
    //-------------------------------------------------------------------
    CRSBackAnalysis *crs=new CRSBackAnalysis;
    CRS_Base crsObject;
    bool backAnalysisFlag=false;
};

#endif // MAINWINDOW_H
