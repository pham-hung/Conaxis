#ifndef GLWIDGET_H
#define GLWIDGET_H
#include <QDebug>
#include <iostream>
#include <QOpenGLWidget>
#include <QOpenGLFunctions>
#include <QOpenGLBuffer>
#include <QOpenGLShaderProgram>
#include <QMatrix4x4>
#include <QVector4D>
#include <QOpenGLVertexArrayObject>
#include <vector>
#include <QFileDialog>
#include <QObject>
#include <QOpenGLWidget>
#include <QString>
#include <QFont>
#include "input.h"
#include <QWheelEvent>
#include <QMessageBox>
#include <QElapsedTimer>
#include <Eigen/Dense>
#include <QDate>
#include <QPoint>
#include <QCursor>
#include <QPixmap>
#include <QRgb>
#include <QTimer>

using namespace std;
using namespace Eigen;

class GLWidget:public QOpenGLWidget, protected QOpenGLFunctions
{
     Q_OBJECT
public:
   //Default constructor
    GLWidget();  
    void initializeGL();
    void paintGL();
    void resizeGL(int w,int h);

    void prepareData();

    void nodeData();
    void initializeNode();
    void paintNode();    

    void tri6EleData();
    void initializeTri6Ele();
    void paintTri6Ele();

    void tri6LineData();
    void initializeTri6Line();
    void paintTri6Line();

    void lineElementData();
    void initializeLineElement();
    void paintLineElement();

    void axisLineData();
    void initializeAxisLineData();
    void paintAxisLineData();

    //Function
    void assignColor(double inVal);
    void vectorColorBand();

public slots:
    void getPars(double minValUI, double maxValUI, int meshCheckUI, int nodeCheckUI, int resultCheckUI,bool autoColor);
    void getDataGL(Ref<MatrixXd> coordUI, Ref<MatrixXd>elementsUI, Ref<MatrixXd> XXUI, Ref<MatrixXd> coordScaleUI);
    void getViewportInfor(float xleft0,float xright0,float ybot0,float ytop0,bool lockView);
    void getColorBandInfor(int noc,int numType,int fontSize,int colorPosition,
                           int nodePlot,int meshPlot,int resultPlot,int axePlot,int datePlot,int titlePlot,int valPlot,QString title);
    void getMousePosition();
//    void updateWidget();

signals:
    void sendSignaltoGUI();
    void sendCurrentValToGUI(vector<double> ScreenValue);
    void sendViewport(float,float,float,float);

protected slots:
    void teardownGL();
    void updateWidget();

protected:
    void keyPressEvent(QKeyEvent *event);
    void keyReleaseEvent(QKeyEvent *event);
    void mousePressEvent(QMouseEvent *event);
    void mouseReleaseEvent(QMouseEvent *event);
    void wheelEvent(QWheelEvent *event);

private:
    void updateWheel();
    void updateLeftMouse();    
    //---------------------
    //Color bufffer
    QOpenGLBuffer BC_node; //For color  
    QOpenGLBuffer BC_tri6Ele;
    QOpenGLBuffer BC_tri6Line;
    QOpenGLBuffer BC_axisLine;
    QOpenGLBuffer BC_lineElement;

    //Vertex buffer
    QOpenGLBuffer B_node; //For node array    
    QOpenGLBuffer B_tri6Ele;
    QOpenGLBuffer B_tri6Line; 
    QOpenGLBuffer B_axisLine;
    QOpenGLBuffer B_lineElement;

    //Vertex array object
    QOpenGLVertexArrayObject O_node; //tet10 element object
    QOpenGLVertexArrayObject O_tri6Ele;
    QOpenGLVertexArrayObject O_tri6Line;
    QOpenGLVertexArrayObject O_axisLine;
    QOpenGLVertexArrayObject O_lineElement;

    //Vertex array shader
    QOpenGLShaderProgram *S_node; //tet10 shaders
    QOpenGLShaderProgram *S_tri6Ele;
    QOpenGLShaderProgram *S_tri6Line;
    QOpenGLShaderProgram *S_axisLine;
    QOpenGLShaderProgram *S_lineElement;

    int u_modelToWorld;
    int u_worldToCamera;
    int u_cameraToView;

    QMatrix4x4 m_projection; //Projection matrix, always set Ortho type
    QMatrix4x4 m_camera; //From camera class, for translation + rotation
    QMatrix4x4 m_transform;  //Transform class, should replaced by transform3D class
    QMatrix4x4 m_Matrix;
    //---------------------
    int deltaWheel;
    int buttonWheel;
    //---------------------
    float redVal, greenVal, blueVal;
    vector <double> ScreenValue;

    QPoint mousePoint;
    QPixmap pixMap;
    QRgb pixVal;

    //---------------------------------------------------------------------------------------------
    //Prepare input matrix
    MatrixXd coord, coordScale; //node

    vector <Vector3f> V_node; //Node plot
    vector <Vector3f> V_node_c; //node_color

    vector <Vector3f> V_tri6Ele;
    vector <Vector3f> V_tri6Ele_c;

    vector <Vector3f> V_tri6Line;
    vector <Vector3f> V_tri6Line_c;

    vector<Vector3f> V_axisLine;
    vector<Vector3f> V_axisLine_c;

    vector <Vector3f> V_lineElement;
    vector <Vector3f> V_lineElement_c;

    int quad4Count, tri3Count, line2Count;

    //For display with QPainter
    QString text;
    QFont font;

    //Prepare global variable
    double xmax, xmin, ymax, ymin, z_max, z_min;
    double xcent, ycent, z_center;
    float xleft0, xright0, ybot0, ytop0; //Defaul value of ortho projection;
    bool lockView;
    float xleft, xright, ybot, ytop;
    float scale_wheel=1;    

    int non,noe;
    bool updateCheck;
    bool autoColor;

    //---------------------------------------------------------------------------------------------
    double minVal,maxVal;
    int meshCheck, resultCheck, nodeCheck;
    MatrixXd  elements, XX;
    float red,green,blue;
    int noc,numType,fontSize,colorPosition;
    int nodePlot,meshPlot,resultPlot,axePlot,datePlot,titlePlot,valPlot;
    QString title;
    QTimer *aTimer = new QTimer;
    QTimer *bTimer= new QTimer;
    //---------------------------------------------------------------------------------------------

};

#endif // GLWIDGET_H
