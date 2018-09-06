#include "glwidget.h"
#include <QPainter>
#include <QKeyEvent>
#include <QDebug>
#include <math.h>
#include <QElapsedTimer>
#define PI 3.14159265

GLWidget::GLWidget()
{
    //    cout<<"OpenGL Widget is re-created"<<endl;
}

void GLWidget::prepareData()
{
    //Calculate center of model
    xmax=coord.col(1).maxCoeff();
    xmin=coord.col(1).minCoeff();
    ymax=coord.col(2).maxCoeff();
    ymin=coord.col(2).minCoeff();
    xcent=(xmax+xmin)*0.5;
    ycent=(ymax+ymin)*0.5;

    non=coord.rows();
    quad4Count=0;
    tri3Count=0;
    line2Count=0;

    for (int j=0;j<elements.rows();j++)
    {
        int nodeCount=elements(j,9);
        if(nodeCount==8){quad4Count++;}
        else if(nodeCount==6){tri3Count++;}
        else if(nodeCount==2){line2Count++;}
    }

    deltaWheel=0;
    m_transform.setToIdentity();
    m_transform.translate(0.0f,0.0f,0.0f);
    m_camera.setToIdentity();
    m_projection.setToIdentity();

    if(lockView==false)
    {
        xleft0=xmin;
        xright0=xmax;
        ybot0=ymin;
        ytop0=ymax;
    }

    scale_wheel=1;
    xleft=xleft0;
    xright=xright0;
    ybot=ybot0;
    ytop=ytop0;

    updateCheck=false;
    nodeData();
    tri6LineData();
    tri6EleData();
    axisLineData();
    lineElementData();
}

void GLWidget::nodeData()
{
    V_node.resize(non);
    V_node_c.resize(non);
    //Data for node-----------------------------------------------------------
    for (int i=0;i<non;i++)
    {
        V_node[i]=Vector3f(coord(i,1),coord(i,2),coord(i,3));
    }

    //Node color
    for (int i=0;i<non;i++)
    {
        V_node_c[i]=Vector3f(1.0f,0.0f,0.5f);//Red color
    }

}

void GLWidget::initializeGL()
{
    emit sendSignaltoGUI();

    initializeOpenGLFunctions();
    connect(context(), SIGNAL(aboutToBeDestroyed()), this, SLOT(teardownGL()), Qt::DirectConnection);
    connect(aTimer,SIGNAL(timeout()),this,SLOT(updateWidget()));
    aTimer->start(10);

    glClearColor(1.0f,1.0f,1.0f,1.0f);
    prepareData();
    initializeNode();
    initializeTri6Line();
    initializeTri6Ele();
    initializeAxisLineData();
    if(line2Count>0)
    {
        initializeLineElement();
    }
}

void GLWidget::initializeNode()
{
    //----------------------------------------------
    //Setup shaders for node
    S_node=new QOpenGLShaderProgram();
    S_node->addShaderFromSourceFile(QOpenGLShader::Vertex,":/shaders/node.vert");
    S_node->addShaderFromSourceFile(QOpenGLShader::Fragment,":/shaders/node.frag");
    S_node->link();
    S_node->bind();

    // Cache Uniform Locations
    u_modelToWorld = S_node->uniformLocation("modelToWorld");
    u_worldToCamera = S_node->uniformLocation("worldToCamera");
    u_cameraToView = S_node->uniformLocation("cameraToView");

    B_node.create();
    B_node.bind();
    B_node.setUsagePattern(QOpenGLBuffer::StaticDraw);
    B_node.allocate(&V_node[0],non*sizeof(V_node[0]));

    O_node.create();
    O_node.bind();
    S_node->enableAttributeArray(0); //Location 0, node coordinates
    S_node->setAttributeBuffer(0,GL_FLOAT,0,3,12); //0 offset, 3 tuplesize, 12 bytes per vector node

    BC_node.create();
    BC_node.bind();
    BC_node.setUsagePattern(QOpenGLBuffer::StaticDraw);
    BC_node.allocate(&V_node_c[0],non*sizeof(V_node_c[0]));
    S_node->enableAttributeArray(1);
    S_node->setAttributeBuffer(1,GL_FLOAT,0,3,12);

    O_node.release();
    B_node.release();
    BC_node.release();
    S_node->release();
    //---------------------------------------------------------------
}

void GLWidget::paintGL()
{
    QPainter painter(this);
    painter.beginNativePainting();
    glEnable(GL_DEPTH_TEST);
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT); //clear color
    glEnable(GL_PROGRAM_POINT_SIZE);
    glEnable(GL_LINE_WIDTH);
    glLineWidth(0.5f);
    glViewport(0,0,this->width(),this->height());
    if(nodePlot==1){paintNode();}
    if(meshPlot==1)
    {
        paintTri6Line();
        if(line2Count>0){paintLineElement();}
    }
    if(resultPlot==1){paintTri6Ele();}
    glDisable(GL_DEPTH_TEST);
    if(axePlot==1){paintAxisLineData();}
    //-------------------------------------------------------------------------------------
    painter.endNativePainting();
    //-------------------------------------------------------------------------------------
    //Draw text on background
    //Cannot bind with nodal location
    //    int fontSize=10;
    //    int noc=11; //Nunmber of contour
    MatrixXd contour=MatrixXd::Zero(noc,1);
    for (int i=0;i<noc;i++)
    {
        contour(i,0)=minVal+(maxVal-minVal)*i/(noc-1);
    }

    painter.setPen(Qt::black);
    painter.setFont(QFont("Arial",fontSize));

    QString maxString, minString, dateString, titleString,tempString;
    dateString=QDate::currentDate().toString("dd.MM.yyyy");

    if(numType==0)
    {
        maxString="Max Value (Blue Color): "+QString::number(maxVal,'f',0);
        minString="Min Value (Red Color): "+QString::number(minVal,'f',0);
    }
    else if(numType==1)
    {
        maxString="Max Value (Blue Color): "+QString::number(maxVal,'f',3);
        minString="Min Value (Red Color): "+QString::number(minVal,'f',3);
    }
    else if(numType==2)
    {
        maxString="Max Value (Blue Color): "+QString::number(maxVal,'e',3);
        minString="Min Value (Red Color): "+QString::number(minVal,'e',3);
    }

    int maxY=2*fontSize;
    if(titlePlot==1)
    {
        painter.drawText(0,maxY,title);
        maxY=maxY+2*fontSize;
    }
    if(datePlot==1)
    {
        painter.drawText(0,maxY,dateString);
        maxY=maxY+2*fontSize;
    }
    if(valPlot==1)
    {
        painter.drawText(0,maxY,maxString);
        maxY=maxY+2*fontSize;
        painter.drawText(0,maxY,minString);
        maxY=maxY+2*fontSize;
    }

    int widthSize=this->width();
    int heightSize=this->height();
    int GradientH=heightSize-maxY-fontSize;
    int GradientW=widthSize-200;
    QLinearGradient m_gradient(QPoint(0,maxY),QPoint(0,maxY+GradientH));
    //from min to max
    m_gradient.setColorAt(0.0,Qt::red);
    m_gradient.setColorAt(0.5,Qt::green);
    m_gradient.setColorAt(1.0,Qt::blue);

    if(colorPosition==0){painter.fillRect(QRect(QPoint(10,maxY),QSize(20,GradientH)),m_gradient);}
    else if(colorPosition==1)
    {
        maxY=20;
        GradientH=heightSize-maxY-fontSize;
        QLinearGradient m_gradient(QPoint(0,maxY),QPoint(0,maxY+GradientH));
        m_gradient.setColorAt(0.0,Qt::red);
        m_gradient.setColorAt(0.5,Qt::green);
        m_gradient.setColorAt(1.0,Qt::blue);
        painter.fillRect(QRect(QPoint(widthSize-100,maxY),QSize(20,GradientH)),m_gradient);
    }
    else if(colorPosition==2)
    {
        //m_gradient.sets(QPoint(0,100),QPoint(0,600));
        QLinearGradient m_gradient(QPoint(100,0),QPoint(GradientW,0));
        m_gradient.setColorAt(0.0,Qt::red);
        m_gradient.setColorAt(0.5,Qt::green);
        m_gradient.setColorAt(1.0,Qt::blue);
        painter.fillRect(QRect(QPoint(100,heightSize-25-2*fontSize),QSize(GradientW,20)),m_gradient);
    }

    for (int i=0;i<noc;i++)
    {
        if(numType==0){tempString=QString::number(contour(i,0),'f',0);}
        else if(numType==1){tempString=QString::number(contour(i,0),'f',3);}
        else if(numType==2){tempString=QString::number(contour(i,0),'e',3);}
        if(colorPosition==0){painter.drawText(30+0.5*fontSize,maxY+0.5*fontSize+i*GradientH/(noc-1),tempString);}
        else if(colorPosition==1)
        {
            maxY=20;
            GradientH=heightSize-maxY-fontSize;
            painter.drawText(widthSize-78,maxY+0.5*fontSize+i*GradientH/(noc-1),tempString);
        }
        else if(colorPosition==2)
        {
            painter.drawText(100+0.5*fontSize+i*GradientW/(noc-1),heightSize-20+0.5*fontSize,tempString);
        }
    }
    //-------------------------------------------------------------------------------------
}

void GLWidget::paintNode()
{
    //-------------------------------------------------------------------------------------
    S_node->bind();
    S_node->setUniformValue(u_worldToCamera,m_camera);
    S_node->setUniformValue(u_cameraToView, m_projection);
    {
        glEnable(GL_LINE_WIDTH);
        glLineWidth(1.0f);
        O_node.bind();
        S_node->setUniformValue(u_modelToWorld, m_transform);
        glDrawArrays(GL_POINTS,0,non); //Draw nodes
        O_node.release();
    }
    S_node->release();
}

void GLWidget::tri6EleData()
{
    noe=elements.rows();
    V_tri6Ele.resize(0);
    V_tri6Ele.reserve(3*tri3Count+6*quad4Count);
    V_tri6Ele_c.resize(0);
    V_tri6Ele_c.reserve(3*tri3Count+6*quad4Count);

    for (int i=0;i<noe;i++)
    {
        unsigned int node1=elements(i,1)-1;
        unsigned int node2=elements(i,2)-1;
        unsigned int node3=elements(i,3)-1;
        unsigned int node4=elements(i,4)-1;
        int nodeCount=elements(i,9);
        if(nodeCount==8)
        {
            V_tri6Ele.push_back(Vector3f(coordScale(node1,1),coordScale(node1,2),0.0f));
            V_tri6Ele.push_back(Vector3f(coordScale(node2,1),coordScale(node2,2),0.0f));
            V_tri6Ele.push_back(Vector3f(coordScale(node3,1),coordScale(node3,2),0.0f));
            V_tri6Ele.push_back(Vector3f(coordScale(node1,1),coordScale(node1,2),0.0f));
            V_tri6Ele.push_back(Vector3f(coordScale(node3,1),coordScale(node3,2),0.0f));
            V_tri6Ele.push_back(Vector3f(coordScale(node4,1),coordScale(node4,2),0.0f));
            if(XX.rows()!=1)
            {
                double val1,val2,val3,val4;
                if(XX.rows()==noe)
                {
                    val1=XX(i,0);
                    val2=XX(i,0);
                    val3=XX(i,0);
                    val4=XX(i,0);
                }
                else
                {
                    val1=XX(node1,0);
                    val2=XX(node2,0);
                    val3=XX(node3,0);
                    val4=XX(node4,0);
                }
                assignColor(val1);
                V_tri6Ele_c.push_back(Vector3f(redVal,greenVal,blueVal));
                assignColor(val2);
                V_tri6Ele_c.push_back(Vector3f(redVal,greenVal,blueVal));
                assignColor(val3);
                V_tri6Ele_c.push_back(Vector3f(redVal,greenVal,blueVal));
                assignColor(val1);
                V_tri6Ele_c.push_back(Vector3f(redVal,greenVal,blueVal));
                assignColor(val3);
                V_tri6Ele_c.push_back(Vector3f(redVal,greenVal,blueVal));
                assignColor(val4);
                V_tri6Ele_c.push_back(Vector3f(redVal,greenVal,blueVal));
            }
            else if(XX.rows()==1)
            {
                if(XX(0,0)==0) //Material show
                {
                    maxVal=elements.col(10).maxCoeff();
                    minVal=elements.col(10).minCoeff();
                    noc=maxVal-minVal+1;
                    if(maxVal==minVal){
                        redVal=1.0f;
                        greenVal=0.0f;
                        blueVal=0.0f;
                        noc=2;
                    }
                    else
                    {
                        double val1=elements(i,10);
                        assignColor(val1);
                    }
                }
                else if(XX(0,0)==1)
                {
                    maxVal=elements.col(11).maxCoeff();
                    minVal=elements.col(11).minCoeff();
                    noc=maxVal-minVal+1;
                    if(maxVal==minVal){
                        redVal=1.0f;
                        greenVal=0.0f;
                        blueVal=0.0f;
                        noc=2;
                    }
                    else
                    {
                        double val1=elements(i,11);
                        assignColor(val1);
                    }

                }

                V_tri6Ele_c.push_back(Vector3f(redVal,greenVal,blueVal));
                V_tri6Ele_c.push_back(Vector3f(redVal,greenVal,blueVal));
                V_tri6Ele_c.push_back(Vector3f(redVal,greenVal,blueVal));
                V_tri6Ele_c.push_back(Vector3f(redVal,greenVal,blueVal));
                V_tri6Ele_c.push_back(Vector3f(redVal,greenVal,blueVal));
                V_tri6Ele_c.push_back(Vector3f(redVal,greenVal,blueVal));
            }

        }
        else if(nodeCount==6)
        {
            V_tri6Ele.push_back(Vector3f(coordScale(node1,1),coordScale(node1,2),0.0f));
            V_tri6Ele.push_back(Vector3f(coordScale(node2,1),coordScale(node2,2),0.0f));
            V_tri6Ele.push_back(Vector3f(coordScale(node3,1),coordScale(node3,2),0.0f));
            if(XX.rows()!=1)
            {
                double val1, val2, val3;
                if(XX.rows()==noe)
                {
                    val1=XX(i,0);
                    val2=XX(i,0);
                    val3=XX(i,0);
                }
                else
                {
                    val1=XX(node1,0);
                    val2=XX(node2,0);
                    val3=XX(node3,0);
                }

                assignColor(val1);
                V_tri6Ele_c.push_back(Vector3f(redVal,greenVal,blueVal));
                assignColor(val2);
                V_tri6Ele_c.push_back(Vector3f(redVal,greenVal,blueVal));
                assignColor(val3);
                V_tri6Ele_c.push_back(Vector3f(redVal,greenVal,blueVal));
            }
            else if(XX.rows()==1)
            {
                if(XX(0,0)==0) //Material show
                {
                    maxVal=elements.col(10).maxCoeff();
                    minVal=elements.col(10).minCoeff();
                    noc=maxVal-minVal+1;
                    if(maxVal==minVal){
                        redVal=1.0f;
                        greenVal=0.0f;
                        blueVal=0.0f;
                        noc=2;
                    }
                    else
                    {
                        double val1=elements(i,10);
                        assignColor(val1);
                    }
                }
                else if(XX(0,0)==1)
                {
                    maxVal=elements.col(11).maxCoeff();
                    minVal=elements.col(11).minCoeff();
                    noc=maxVal-minVal+1;
                    if(maxVal==minVal){
                        redVal=1.0f;
                        greenVal=0.0f;
                        blueVal=0.0f;
                        noc=2;
                    }
                    else
                    {
                        double val1=elements(i,11);
                        assignColor(val1);
                    }
                }

                V_tri6Ele_c.push_back(Vector3f(redVal,greenVal,blueVal));
                V_tri6Ele_c.push_back(Vector3f(redVal,greenVal,blueVal));
                V_tri6Ele_c.push_back(Vector3f(redVal,greenVal,blueVal));
            }
        }

    }

}

void GLWidget::initializeTri6Ele()
{
    //----------------------------------------------
    //Setup shaders for node
    S_tri6Ele=new QOpenGLShaderProgram();
    S_tri6Ele->addShaderFromSourceFile(QOpenGLShader::Vertex,":/shaders/t10element.vert");
    S_tri6Ele->addShaderFromSourceFile(QOpenGLShader::Fragment,":/shaders/t10element.frag");
    S_tri6Ele->link();
    S_tri6Ele->bind();

    // Cache Uniform Locations
    u_modelToWorld = S_tri6Ele->uniformLocation("modelToWorld");
    u_worldToCamera = S_tri6Ele->uniformLocation("worldToCamera");
    u_cameraToView = S_tri6Ele->uniformLocation("cameraToView");

    B_tri6Ele.create();
    B_tri6Ele.bind();
    B_tri6Ele.setUsagePattern(QOpenGLBuffer::StaticDraw);
    B_tri6Ele.allocate(&V_tri6Ele[0],V_tri6Ele.size()*sizeof(V_tri6Ele[0]));

    O_tri6Ele.create();
    O_tri6Ele.bind();
    S_tri6Ele->enableAttributeArray(4); //Location 0, node coordinates
    S_tri6Ele->setAttributeBuffer(4,GL_FLOAT,0,3,12); //0 offset, 3 tuplesize, 12 bytes per vector node

    BC_tri6Ele.create();
    BC_tri6Ele.bind();
    BC_tri6Ele.setUsagePattern(QOpenGLBuffer::StaticDraw);
    BC_tri6Ele.allocate(&V_tri6Ele_c[0],V_tri6Ele_c.size()*sizeof(V_tri6Ele_c[0]));
    S_tri6Ele->enableAttributeArray(5);
    S_tri6Ele->setAttributeBuffer(5,GL_FLOAT,0,3,12);

    O_tri6Ele.release();
    B_tri6Ele.release();
    BC_tri6Ele.release();
    S_tri6Ele->release();
    //---------------------------------------------------------------
}

void GLWidget::paintTri6Ele()
{
    S_tri6Ele->bind();
    S_tri6Ele->setUniformValue(u_worldToCamera,m_camera);
    S_tri6Ele->setUniformValue(u_cameraToView, m_projection);

    {
        O_tri6Ele.bind();
        S_tri6Ele->setUniformValue(u_modelToWorld, m_transform);
        glDrawArrays(GL_TRIANGLES,0,3*tri3Count+6*quad4Count);
        O_tri6Ele.release();
    }
    S_tri6Ele->release();
}

void GLWidget::tri6LineData()
{
    noe=elements.rows();
    V_tri6Line.resize(0);
    V_tri6Line.reserve(6*tri3Count+8*quad4Count);
    for (int i=0;i<noe;i++)
    {
        unsigned int node1=elements(i,1)-1;
        unsigned int node2=elements(i,2)-1;
        unsigned int node3=elements(i,3)-1;
        unsigned int node4=elements(i,4)-1;
        int nodeCount=elements(i,9);
        if(nodeCount==8)
        {
            V_tri6Line.push_back(Vector3f(coord(node1,1),coord(node1,2),0.0f));
            V_tri6Line.push_back(Vector3f(coord(node2,1),coord(node2,2),0.0f));
            V_tri6Line.push_back(Vector3f(coord(node2,1),coord(node2,2),0.0f));
            V_tri6Line.push_back(Vector3f(coord(node3,1),coord(node3,2),0.0f));
            V_tri6Line.push_back(Vector3f(coord(node3,1),coord(node3,2),0.0f));
            V_tri6Line.push_back(Vector3f(coord(node4,1),coord(node4,2),0.0f));
            V_tri6Line.push_back(Vector3f(coord(node4,1),coord(node4,2),0.0f));
            V_tri6Line.push_back(Vector3f(coord(node1,1),coord(node1,2),0.0f));

        }
        else if(nodeCount==6)
        {
            V_tri6Line.push_back(Vector3f(coord(node1,1),coord(node1,2),0.0f));
            V_tri6Line.push_back(Vector3f(coord(node2,1),coord(node2,2),0.0f));
            V_tri6Line.push_back(Vector3f(coord(node2,1),coord(node2,2),0.0f));
            V_tri6Line.push_back(Vector3f(coord(node3,1),coord(node3,2),0.0f));
            V_tri6Line.push_back(Vector3f(coord(node3,1),coord(node3,2),0.0f));
            V_tri6Line.push_back(Vector3f(coord(node1,1),coord(node1,2),0.0f));
        }
    }

}

void GLWidget::initializeTri6Line()
{
    S_tri6Line=new QOpenGLShaderProgram();
    S_tri6Line->addShaderFromSourceFile(QOpenGLShader::Vertex,":/shaders/line.vert");
    S_tri6Line->addShaderFromSourceFile(QOpenGLShader::Fragment,":/shaders/line.frag");
    S_tri6Line->link();
    S_tri6Line->bind();

    // Cache Uniform Locations
    u_modelToWorld = S_tri6Line->uniformLocation("modelToWorld");
    u_worldToCamera = S_tri6Line->uniformLocation("worldToCamera");
    u_cameraToView = S_tri6Line->uniformLocation("cameraToView");

    B_tri6Line.create();
    B_tri6Line.bind();
    B_tri6Line.setUsagePattern(QOpenGLBuffer::StaticDraw);
    B_tri6Line.allocate(&V_tri6Line[0],V_tri6Line.size()*sizeof(V_tri6Line[0]));

    O_tri6Line.create();
    O_tri6Line.bind();
    S_tri6Line->enableAttributeArray(2);
    S_tri6Line->setAttributeBuffer(2,GL_FLOAT,0,3,12); //0 offset, 3 tuplesize, 12 bytes per vector node
    O_tri6Line.release();
    B_tri6Line.release();
    S_tri6Line->release();

}

void GLWidget::paintTri6Line()
{
    S_tri6Line->bind();
    S_tri6Line->setUniformValue(u_worldToCamera,m_camera);
    S_tri6Line->setUniformValue(u_cameraToView, m_projection);
    {
        glLineWidth(0.5f);
        O_tri6Line.bind();
        S_tri6Line->setUniformValue(u_modelToWorld, m_transform);
        glDrawArrays(GL_LINES,0,6*tri3Count+8*quad4Count);
        O_tri6Line.release();
    }
    S_tri6Line->release();
}

void GLWidget::lineElementData()
{
    noe=elements.rows();
    V_lineElement.resize(0);
    V_lineElement.reserve(2*line2Count);
    for (int i=0;i<noe;i++)
    {
        unsigned int node1=elements(i,1)-1;
        unsigned int node2=elements(i,2)-1;
        int nodeCount=elements(i,9);
        if(nodeCount==2)
        {
            V_lineElement.push_back(Vector3f(coord(node1,1),coord(node1,2),0.0f));
            V_lineElement.push_back(Vector3f(coord(node2,1),coord(node2,2),0.0f));
        }
    }
}

void GLWidget::initializeLineElement()
{
    S_lineElement=new QOpenGLShaderProgram();
    S_lineElement->addShaderFromSourceFile(QOpenGLShader::Vertex,":/shaders/line1D.vert");
    S_lineElement->addShaderFromSourceFile(QOpenGLShader::Fragment,":/shaders/line1D.frag");
    S_lineElement->link();
    S_lineElement->bind();

    // Cache Uniform Locations
    u_modelToWorld = S_lineElement->uniformLocation("modelToWorld");
    u_worldToCamera = S_lineElement->uniformLocation("worldToCamera");
    u_cameraToView = S_lineElement->uniformLocation("cameraToView");

    B_lineElement.create();
    B_lineElement.bind();
    B_lineElement.setUsagePattern(QOpenGLBuffer::StaticDraw);
    B_lineElement.allocate(&V_lineElement[0],V_lineElement.size()*sizeof(V_lineElement[0]));

    O_lineElement.create();
    O_lineElement.bind();

    S_lineElement->enableAttributeArray(8);
    S_lineElement->setAttributeBuffer(8,GL_FLOAT,0,3,12); //0 offset, 3 tuplesize, 12 bytes per vector node

    O_lineElement.release();
    B_lineElement.release();
    S_lineElement->release();
}

void GLWidget::paintLineElement()
{
    S_lineElement->bind();
    S_lineElement->setUniformValue(u_worldToCamera,m_camera);
    S_lineElement->setUniformValue(u_cameraToView, m_projection);
    {
        glLineWidth(3.0f);
        O_lineElement.bind();
        S_lineElement->setUniformValue(u_modelToWorld, m_transform);
        glDrawArrays(GL_LINES,0,2*line2Count);
        O_lineElement.release();
    }
    S_lineElement->release();
}

//Axis Data--------------
void GLWidget::axisLineData()
{
    V_axisLine.resize(0);
    V_axisLine_c.resize(0);
    double length1=(xmax-xmin);
    double length2=(ymax-ymin);
    double axisLength=0;
    if(length1>length2)
    {
        axisLength=0.3*length2;
    }
    else
    {
        axisLength=0.3*length1;
    }
    axisLength=abs(axisLength);

    //X Axis
    V_axisLine.push_back(Vector3f(0.0f,0.0f,0.0f));
    V_axisLine.push_back(Vector3f(axisLength,0.0f,0.0f));
    V_axisLine.push_back(Vector3f(axisLength,0.0f,0.0f));
    V_axisLine.push_back(Vector3f(3*axisLength/4.0,axisLength/4.0f,0.0f));
    V_axisLine.push_back(Vector3f(axisLength,0.0f,0.0f));
    V_axisLine.push_back(Vector3f(3*axisLength/4.0,-axisLength/4.0f,0.0f));
    //Y Axis
    V_axisLine.push_back(Vector3f(0.0f,0.0f,0.0f));
    V_axisLine.push_back(Vector3f(0.0f,axisLength,0.0f));
    V_axisLine.push_back(Vector3f(0.0f,axisLength,0.0f));
    V_axisLine.push_back(Vector3f(axisLength/4.0f,3*axisLength/4,0.0f));
    V_axisLine.push_back(Vector3f(0.0f,axisLength,0.0f));
    V_axisLine.push_back(Vector3f(-axisLength/4.0f,3*axisLength/4,0.0f));

    //X color is Red
    //Y color is Blue
    V_axisLine_c.push_back(Vector3f(0.0f,0.0f,0.0f));
    V_axisLine_c.push_back(Vector3f(0.0f,0.0f,0.0f));
    V_axisLine_c.push_back(Vector3f(0.0f,0.0f,0.0f));
    V_axisLine_c.push_back(Vector3f(0.0f,0.0f,0.0f));
    V_axisLine_c.push_back(Vector3f(0.0f,0.0f,0.0f));
    V_axisLine_c.push_back(Vector3f(0.0f,0.0f,0.0f));

    V_axisLine_c.push_back(Vector3f(0.0f,0.0f,0.0f));
    V_axisLine_c.push_back(Vector3f(0.0f,0.0f,0.0f));
    V_axisLine_c.push_back(Vector3f(0.0f,0.0f,0.0f));
    V_axisLine_c.push_back(Vector3f(0.0f,0.0f,0.0f));
    V_axisLine_c.push_back(Vector3f(0.0f,0.0f,0.0f));
    V_axisLine_c.push_back(Vector3f(0.0f,0.0f,0.0f));


}

void GLWidget::initializeAxisLineData()
{
    //Setup shaders for element
    S_axisLine=new QOpenGLShaderProgram();
    S_axisLine->addShaderFromSourceFile(QOpenGLShader::Vertex,":/shaders/axis.vert");
    S_axisLine->addShaderFromSourceFile(QOpenGLShader::Fragment,":/shaders/axis.frag");
    S_axisLine->link();
    S_axisLine->bind();

    // Cache Uniform Locations
    u_modelToWorld = S_axisLine->uniformLocation("modelToWorld");
    u_worldToCamera = S_axisLine->uniformLocation("worldToCamera");
    u_cameraToView = S_axisLine->uniformLocation("cameraToView");

    B_axisLine.create();
    B_axisLine.bind();
    B_axisLine.setUsagePattern(QOpenGLBuffer::StaticDraw);
    B_axisLine.allocate(&V_axisLine[0],V_axisLine.size()*sizeof(V_axisLine[0]));

    O_axisLine.create();
    O_axisLine.bind();
    S_axisLine->enableAttributeArray(6); //Location 0, node coordinates
    S_axisLine->setAttributeBuffer(6,GL_FLOAT,0,3,12); //0 offset, 3 tuplesize, 12 bytes per vector node

    BC_axisLine.create();
    BC_axisLine.bind();
    BC_axisLine.setUsagePattern(QOpenGLBuffer::StaticDraw);
    BC_axisLine.allocate(&V_axisLine_c[0],V_axisLine_c.size()*sizeof(V_axisLine_c[0]));

    S_axisLine->enableAttributeArray(7);
    S_axisLine->setAttributeBuffer(7,GL_FLOAT,0,3,12);

    O_axisLine.release();
    BC_axisLine.release();
    BC_axisLine.release();
    S_axisLine->release();
}

void GLWidget::paintAxisLineData()
{
    S_axisLine->bind();
    S_axisLine->setUniformValue(u_worldToCamera,m_camera);
    S_axisLine->setUniformValue(u_cameraToView, m_projection);
    {
        glLineWidth(3.0f);
        O_axisLine.bind();
        S_axisLine->setUniformValue(u_modelToWorld, m_transform);
        glDrawArrays(GL_LINES,0,12);
        O_axisLine.release();
    }
    S_axisLine->release();
}

void GLWidget::assignColor(double inVal)
{
    double midVal=0.5*(maxVal+minVal);
    redVal=0.0f;
    greenVal=0.0f;
    blueVal=0.0f;
    //minColor is red
    //maxColor is blue
    if(maxVal==minVal)
    {
        redVal=1.0f;
        greenVal=0.0f;
        blueVal=0.0f;
    }

    else if (inVal>maxVal || inVal< minVal)
    {
        redVal=0.5f;
        greenVal=0.5f;
        blueVal=0.5f;
    }
    else if(inVal>=minVal && inVal<=midVal)
    {
        blueVal=0.0f;
        redVal=float((inVal-midVal)/(minVal-midVal));
        greenVal=float((minVal-inVal)/(minVal-midVal));

    }
    else if(inVal<=maxVal && inVal>midVal)
    {
        redVal=0.0f;
        greenVal=float((inVal-maxVal)/(midVal-maxVal));
        blueVal=float((midVal-inVal)/(midVal-maxVal));

    }
    else
    {
        QMessageBox::critical(this,"ERROR","OUT OF RANGE COLOR");

    }

}

void GLWidget::getPars(double minValUI, double maxValUI, int meshCheckUI, int nodeCheckUI, int resultCheckUI,bool autoColor)
{
    this->autoColor=autoColor;
    minVal=minValUI;
    maxVal=maxValUI;
    meshCheck=meshCheckUI;
    nodeCheck=nodeCheckUI;
    resultCheck=resultCheckUI;
    //    if(lockView==false)
    //    {

    //    }
}

void GLWidget::getDataGL(Ref<MatrixXd> coordUI, Ref<MatrixXd> elementsUI, Ref<MatrixXd> XXUI, Ref<MatrixXd> coordScaleUI)
{
    coord=MatrixXd::Zero(coordUI.rows(),coordUI.cols());
    coord=coordUI;

    elements=MatrixXd::Zero(elementsUI.rows(),elementsUI.cols());
    elements=elementsUI;

    coordScale=MatrixXd::Zero(coordScaleUI.rows(),coordScaleUI.cols());
    coordScale=coordScaleUI;

    XX=MatrixXd::Zero(XXUI.rows(),XXUI.cols());
    XX=XXUI;

    if(autoColor==true)
    {
        minVal=XX.col(0).minCoeff();
        maxVal=XX.col(0).maxCoeff();
    }

}

void GLWidget::getViewportInfor(float xleft0, float xright0, float ybot0, float ytop0, bool lockView)
{
    this->lockView=lockView;
    if(lockView==true)
    {
        this->xleft0=xleft0;
        this->xright0=xright0;
        this->ybot0=ybot0;
        this->ytop0=ytop0;
    }
}

void GLWidget::getColorBandInfor(int noc,int numType,int fontSize,int colorPosition,
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

void GLWidget::getMousePosition()
{    
    mousePoint=QCursor::pos();
    mousePoint=this->mapFromGlobal(mousePoint);
}

void GLWidget::updateWidget()
{
    Input::update();
    updateWheel();
    updateLeftMouse();
    if(updateCheck==true)
    {
        float aspect;
        aspect=(float)(this->width()/this->height());
        xright=xright+(ytop-ybot)*aspect;
        m_projection.setToIdentity();
        m_projection.ortho(xleft,xright,ybot,ytop,0,0);
        resizeGL(this->width(),this->height());
        QOpenGLWidget::update();
        updateCheck=false;
        emit sendViewport(xleft,xright,ybot,ytop);
    }
}

void GLWidget::resizeGL(int w, int h)
{
    float aspect;
    aspect=(float)w/h;
    if(lockView==false)
    {
        xright=xleft+(ytop-ybot)*aspect;

    }
    else
    {
        xleft=xleft0;
        xright=xright0;
        ybot=ybot0;
        ytop=ytop0;
    }

    m_projection.setToIdentity();
    m_projection.ortho(xleft,xright,ybot,ytop,-500000,500000);
    QOpenGLWidget::update();
}

void GLWidget::teardownGL()
{
    //destroy object
    O_node.destroy();
    O_tri6Ele.destroy();
    O_tri6Line.destroy();

    //destroy buffer
    B_node.destroy();
    BC_node.destroy();
    BC_tri6Ele.destroy();
    B_tri6Ele.destroy();

    //delete shadares
    S_node->deleteLater();
    S_tri6Ele->deleteLater();
    S_tri6Line->deleteLater();
    O_node.deleteLater();
    O_tri6Ele.deleteLater();
    O_tri6Line.deleteLater();
}

//Mouse event-----------------
void GLWidget::mousePressEvent(QMouseEvent *event)
{
    double inVal;
    ScreenValue.resize(3);
    Input::registerMousePress(event->button());
    if(event->button()==1)
    {
        mousePoint=QCursor::pos();
        mousePoint=this->mapFromGlobal(mousePoint);
        QImage image = this->grabFramebuffer();
        QRgb rgb=image.pixel(mousePoint.x(),mousePoint.y());
        redVal=qRed(rgb);
        greenVal=qGreen(rgb);
        blueVal=qBlue(rgb);

        double midVal;
        midVal=0.5*(maxVal+minVal);

        if ((redVal==255 && greenVal==255) && blueVal==255)
        {
            inVal=9999999;
        }
        else if (redVal==0)
        {
            inVal=midVal-(blueVal*(midVal-maxVal)/255);

        }
        else
        {
            inVal=(redVal*(minVal-midVal)/255)+midVal;
        }
        m_Matrix=m_projection*m_camera*m_transform;
        MatrixXd myMatrix=MatrixXd::Zero(4,4);
        MatrixXd Location=MatrixXd::Zero(4,1);
        MatrixXd realLocation=MatrixXd::Zero(4,1);

        myMatrix(0,0)=m_Matrix.data()[0];
        myMatrix(1,0)=m_Matrix.data()[1];
        myMatrix(2,0)=m_Matrix.data()[2];
        myMatrix(3,0)=m_Matrix.data()[3];

        myMatrix(0,1)=m_Matrix.data()[4];
        myMatrix(1,1)=m_Matrix.data()[5];
        myMatrix(2,1)=m_Matrix.data()[6];
        myMatrix(3,1)=m_Matrix.data()[7];

        myMatrix(0,2)=m_Matrix.data()[8];
        myMatrix(1,2)=m_Matrix.data()[9];
        myMatrix(2,2)=m_Matrix.data()[10];
        myMatrix(3,2)=m_Matrix.data()[11];

        myMatrix(0,3)=m_Matrix.data()[12];
        myMatrix(1,3)=m_Matrix.data()[13];
        myMatrix(2,3)=m_Matrix.data()[14];
        myMatrix(3,3)=m_Matrix.data()[15];

        QPoint NewMousePoint=QCursor::pos();
        NewMousePoint=this->mapFromGlobal(NewMousePoint);
        int w=this->width();
        int h=this->height();
        double centerx=w/2;
        double centery=h/2;

        double xCoordinates=(NewMousePoint.x()-centerx)/(0.5f*w);
        double yCoordinates=-(NewMousePoint.y()-centery)/(0.5f*h);
        Location<<xCoordinates,yCoordinates,0,1;
        realLocation=myMatrix.inverse()*Location;
        ScreenValue[0]=inVal;
        ScreenValue[1]=realLocation(0,0);
        ScreenValue[2]=realLocation(1,0);
        emit sendCurrentValToGUI(ScreenValue);
    }
}

void GLWidget::mouseReleaseEvent(QMouseEvent *event)
{
    Input::registerMouseRelease(event->button());
}

//Key events---------------
void GLWidget::keyPressEvent(QKeyEvent *event)
{
    if (event->isAutoRepeat())
    {
        event->ignore();
    }
    else
    {
        Input::registerKeyPress(event->key());
    }
}

void GLWidget::keyReleaseEvent(QKeyEvent *event)
{
    if (event->isAutoRepeat())
    {
        event->ignore();
    }
    else
    {
        Input::registerKeyRelease(event->key());
    }
}

void  GLWidget::wheelEvent(QWheelEvent *event)
{
    deltaWheel=event->delta();
}

void GLWidget::updateWheel()
{
    if(deltaWheel!=0)
    {
        updateCheck=true;
        double sensitive=1.05;
        scale_wheel=deltaWheel/120;

        double length_x=0.5*(-xleft+xright);
        double length_y=0.5*(-ybot+ytop);
        float center_point_x=0.5*(xleft+xright);
        float center_point_y=0.5*(ybot+ytop);

        if(deltaWheel>0){
            scale_wheel=sensitive*deltaWheel/120;
            xleft=center_point_x-length_x*scale_wheel;
            xright=center_point_x+length_x*scale_wheel;
            ybot=center_point_y-length_y*scale_wheel;
            ytop=center_point_y+length_y*scale_wheel;
        }
        else

        {
            scale_wheel=-1/(sensitive*deltaWheel/120);
            xleft=center_point_x-length_x*scale_wheel;
            xright=center_point_x+length_x*scale_wheel;
            ybot=center_point_y-length_y*scale_wheel;
            ytop=center_point_y+length_y*scale_wheel;

        }
    }
    deltaWheel=0;
    scale_wheel=1;
}

void GLWidget::updateLeftMouse()
{
    //Giữ chuột trái, di chuyển
    if (Input::buttonPressed(Qt::LeftButton))
    {
        updateCheck=true;
        float dx=xright-xleft;
        float dy=ytop-ybot;
        float dd=max(dy,dx);
        static const float transSpeed = 0.002*dd; //set transpeed is 5% of dimension

        xleft=xleft-transSpeed * Input::mouseDelta().x();
        xright=xright-transSpeed * Input::mouseDelta().x();
        ybot=ybot+transSpeed * Input::mouseDelta().y();
        ytop=ytop+transSpeed * Input::mouseDelta().y();
    }
}

