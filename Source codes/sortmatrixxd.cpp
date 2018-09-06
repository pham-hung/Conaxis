#include "sortmatrixxd.h"

SortMatrixXd::SortMatrixXd()
{

}

void SortMatrixXd::sortByCol(Ref<MatrixXd> matrixA, int col)
{
    vectorA.resize(matrixA.rows());
    for(int i=0;i<matrixA.rows();i++)
    {
        vectorA[i].Val=matrixA(i,col);
        vectorA[i].rowIndex=i;
    }
    MatrixXd matrixCopy=matrixA;
    sort(vectorA.begin(),vectorA.end(),sortByVal);
    for(int i=0;i<vectorA.size();i++)
    {
        matrixA.row(i)=matrixCopy.row(vectorA[i].rowIndex);
    }
}

void SortMatrixXd::sortByRow(Ref<MatrixXd> matrixA, int row)
{
    vectorA.resize(matrixA.cols());
    for(int i=0;i<matrixA.cols();i++)
    {
        vectorA[i].Val=matrixA(row,i);
        vectorA[i].colIndex=i;
    }
    MatrixXd matrixCopy=matrixA;
    std::sort(vectorA.begin(),vectorA.end(),sortByVal);
    for(int i=0;i<vectorA.size();i++)
    {
        matrixA.col(i)=matrixCopy.col(vectorA[i].colIndex);
    }
}

void SortMatrixXd::findXY(const Ref<const MatrixXd> coordinates, double x0, double x1, double y0, double y1, MatrixXd &result)
{
    vector<int> nodeList;
    for (int i=0;i<coordinates.rows();i++)
    {
        double tol=1e-3;
        double x=coordinates(i,1);
        double y=coordinates(i,2);
        int index=coordinates(i,0);
        double temp;
        if(x0>x1)
        {
            temp=x0;
            x0=x1;
            x1=temp;
        }

        if(y0>y1)
        {
            temp=y0;
            y0=y1;
            y1=temp;
        }

        if(x>(x0-tol)&& x<(x1+tol))
        {
            if(y>(y0-tol)&& y<(y1+tol))
            {
               nodeList.push_back(index);
            }
        }
    }
    result.resize(nodeList.size(),1);
    sort(nodeList.begin(),nodeList.end());

    for(int i=0;i<nodeList.size();i++)
    {
        result(i,0)=nodeList[i];
    }    

}

void SortMatrixXd::findX(const Ref<const MatrixXd> coordinates, double x0, double x1, MatrixXd &result)
{
    vector<int> nodeList;
    for (int i=0;i<coordinates.rows();i++)
    {
        double tol=1e-3;
        double x=coordinates(i,1);
        double y=coordinates(i,2);
        int index=coordinates(i,0);
        double temp;
        if(x0>x1)
        {
            temp=x0;
            x0=x1;
            x1=temp;
        }



        if(x>(x0-tol)&&x<(x1+tol))
        {
            nodeList.push_back(index);
        }
    }
    result.resize(nodeList.size(),1);
    sort(nodeList.begin(),nodeList.end());

    for(int i=0;i<nodeList.size();i++)
    {
        result(i,0)=nodeList[i];
    }
}

void SortMatrixXd::findY(const Ref<const MatrixXd> coordinates,double y0, double y1, MatrixXd &result)
{
    vector<int> nodeList;
    for (int i=0;i<coordinates.rows();i++)
    {
        double tol=1e-3;
        double x=coordinates(i,1);
        double y=coordinates(i,2);
        int index=coordinates(i,0);
        double temp;

        if(y0>y1)
        {
            temp=y0;
            y0=y1;
            y1=temp;
        }

        if(y>(y0-tol)&&y<(y1+tol))
        {
            nodeList.push_back(index);
        }
    }
    result.resize(nodeList.size(),1);
    sort(nodeList.begin(),nodeList.end());

    for(int i=0;i<nodeList.size();i++)
    {
        result(i,0)=nodeList[i];
    }
}

void SortMatrixXd::addMatrix(MatrixXd &sourceMatrix, const Ref<const MatrixXd> addedMatrix, int col)
{
    vector<int> sourceIndex;
    sourceIndex.resize(sourceMatrix.rows());
    vector<int>::iterator it;
    vector<int> notFoundIndex;
    notFoundIndex.resize(0);

    //create sourceIndex
    for(int i=0;i<sourceMatrix.rows();i++)
    {
        sourceIndex[i]=int(sourceMatrix(i,col));
    }

    //Find elment from addedMatrix to source Matrix;
    //In founded, added
    for (int j=0;j<addedMatrix.rows();j++)
    {
        int Val=int(addedMatrix(j,col));
        it=find(sourceIndex.begin(),sourceIndex.end(),Val);
        if(it!=sourceIndex.end())
        {
            int index=distance(sourceIndex.begin(),it);
            sourceMatrix.row(index)=sourceMatrix.row(index)+addedMatrix.row(j);
            sourceMatrix(index,col)=Val;
        }
        else
        {
            notFoundIndex.push_back(j);
        }

    }

    //Added notFoundIdex to source Matrix
    MatrixXd copyMatrix=sourceMatrix;
    sourceMatrix.resize(sourceMatrix.rows()+notFoundIndex.size(),sourceMatrix.cols());
    {
        for(int i=0;i<copyMatrix.rows();i++)
        {
           sourceMatrix.row(i)= copyMatrix.row(i);
        }
        for (int i=0;i<notFoundIndex.size();i++)
        {
            int index=notFoundIndex[i];
            sourceMatrix.row(copyMatrix.rows()+i)=addedMatrix.row(index);
        }
    }
}

void SortMatrixXd::addMatrixReplace(MatrixXd &sourceMatrix, const Ref<const MatrixXd> addedMatrix, int col)
{
    vector<int> sourceIndex;
    sourceIndex.resize(sourceMatrix.rows());
    vector<int>::iterator it;
    vector<int> notFoundIndex;
    notFoundIndex.resize(0);

    //create sourceIndex
    for(int i=0;i<sourceMatrix.rows();i++)
    {
        sourceIndex[i]=int(sourceMatrix(i,col));
    }

    //Find elment from addedMatrix to source Matrix;
    //if founded matrix, replace values
    for (int j=0;j<addedMatrix.rows();j++)
    {
        int Val=int(addedMatrix(j,col));
        it=find(sourceIndex.begin(),sourceIndex.end(),Val);
        if(it!=sourceIndex.end())
        {
            int index=distance(sourceIndex.begin(),it);
            sourceMatrix.row(index)=addedMatrix.row(j);
            sourceMatrix(index,col)=Val;
        }
        else
        {
            notFoundIndex.push_back(j);
        }

    }

    //Added notFoundIdex to source Matrix
    MatrixXd copyMatrix=sourceMatrix;
    sourceMatrix.resize(sourceMatrix.rows()+notFoundIndex.size(),sourceMatrix.cols());
    {
        for(int i=0;i<copyMatrix.rows();i++)
        {
           sourceMatrix.row(i)= copyMatrix.row(i);
        }
        for (int i=0;i<notFoundIndex.size();i++)
        {
            int index=notFoundIndex[i];
            sourceMatrix.row(copyMatrix.rows()+i)=addedMatrix.row(index);
        }
    }
}

void SortMatrixXd::findElementXY(const Ref<const MatrixXd> coordinates, const Ref<const MatrixXd> elements, double x0, double x1, double y0, double y1, MatrixXd &result)
{
    swapValue(x0,x1);
    swapValue(y0,y1);
    vector<int> resultElement;
    double tol=1e-3;
    for (int i=0;i<elements.rows();i++)
    {
        int node1=elements(i,1)-1;
        int node2=elements(i,2)-1;
        int node3=elements(i,3)-1;
        int node4=elements(i,4)-1;

        double xCenter, yCenter;

        int nodeCount=elements(i,9);
        if(nodeCount==6)
        {
            xCenter=(coordinates(node1,1)+coordinates(node2,1)+coordinates(node3,1))/3;
            yCenter=(coordinates(node1,2)+coordinates(node2,2)+coordinates(node3,2))/3;
        }
        else if (nodeCount==8)
        {
            xCenter=(coordinates(node1,1)+coordinates(node2,1)+coordinates(node3,1)+coordinates(node4,1))/4;
            yCenter=(coordinates(node1,2)+coordinates(node2,2)+coordinates(node3,2)+coordinates(node4,2))/4;
        }        

        if(xCenter<x1+tol && xCenter>x0-tol)
        {
            if(yCenter < y1+tol && yCenter>y0-tol)
            {
                resultElement.push_back(i);
            }
        }
    }
    result.resize(resultElement.size(),1);
    for(int i=0;i<resultElement.size();i++)
    {
        result(i,0)=resultElement[i]+1; //increase base 1;
    }    
}

void SortMatrixXd::swapValue(double &x, double &y)
{
    double temp;
    if(x>y)
    {
        x=temp;
        x=y;
        y=temp;
    }
}

void SortMatrixXd::PlaneLineIntersect(Vector3d PointP0, Vector3d PointP1, Vector3d PointV0, Vector3d nVector,Vector3d IntersectionVector)
{
    //P0 is P0 point
    //P1 is P1 point
    //V0 is point belongs to the plane
    //nVector is the normal vector of plane
    //IntersectionVector is vector of Itersection point

    double D,N,sI;
    int checkInt=0;
    Vector3d uVector, wVector;
    checkInt=0;
    IntersectionVector<<0,0,0;
    uVector=PointP1-PointP0;
    wVector=PointP0-PointV0;
    D=nVector.dot(uVector);
    N=-nVector.dot(wVector);

    if(abs(D)<1e-15)
    {
        if(N==0)
        {
            checkInt=2;
        }
        else
        {
            checkInt=0;
        }

    }
    sI=N/D;
    IntersectionVector=PointP0+sI*uVector;
}


bool sortByVal(const element &lhs, const element &rhs)
{

    return lhs.Val<rhs.Val;
}

