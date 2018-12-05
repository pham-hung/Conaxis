#include "axissymmetric_2d.h"

AxisSymmetric_2D::AxisSymmetric_2D(QObject *parent) : QObject(parent)
{
    cout<<"----------------------"<<endl;
    cout<<"Create new analysis   "<<endl;
}

//-------------------------------------------------
//Prepare Data from Classes
void AxisSymmetric_2D::prepareData()
{
    initializeData();
    //Calculate number of step
    nos=stage.size();
    nob=boundary.size();
    ns=0;
    for (int i=0;i<nos;i++)
    {
        ns=ns+stage[i].subStep;
    }

    //Initial Load vector stress vector

    {
        F=MatrixXd::Zero(totalDof,1);
        F0=MatrixXd::Zero(totalDof,1);
        dF=MatrixXd::Zero(totalDof,1);
        FyGravity=MatrixXd::Zero(non,2);
        for(int i=0;i<non;i++)
        {
            FyGravity(i,0)=coordinates(i,0);
        }
        XX=MatrixXd::Zero(totalDof,1);
        X=MatrixXd::Zero(totalDof,ns);
        X0=MatrixXd::Zero(totalDof,1);
        U=MatrixXd::Zero(non,ns);
        V=MatrixXd::Zero(non,ns);
        P=MatrixXd::Zero(non,ns);
        PizoResult=MatrixXd(ns,1);
        Panalytical=MatrixXd::Zero(non,ns);
        Perror=MatrixXd::Zero(non,ns);
    }
    {
        Sxx=MatrixXd::Zero(noe,ns);
        Syy=MatrixXd::Zero(noe,ns);
        Sxy=MatrixXd::Zero(noe,ns);
        Pore=MatrixXd::Zero(noe,ns);
        SyyInterpolation=MatrixXd::Zero(noe,1);
        SyyGravity=MatrixXd::Zero(noe,1);
    }
    {
        HydraulicVertial=MatrixXd::Zero(noe,ns);
        HydraulicHorizontal=MatrixXd::Zero(noe,ns);
        BulkModulus=MatrixXd::Zero(noe,ns);
    }

    //Debug
    cout<<"Number of total steps     : "<<ns<<endl;
    cout<<"Number of total nodes     : "<<non<<endl;
    cout<<"Number of total elems     : "<<noe<<endl;
    cout<<"Number of total stages    : "<<nos<<endl;
    cout<<"Number of total boundaries: "<<nob<<endl;
    cout<<"Total degree of freedom   : "<<totalDof<<endl;
    createGravityLoad();
}

void AxisSymmetric_2D::initializeData()
{
    {
        Cf=1e-7;
        Cs=0;
        acel=9.81;
        gf=projectParameters[3];
        BiotCoeff=1;
        SkemptonCoeff=1;
    }
    //-----------------------------------------
    {
        non=coordinates.rows();
        noe=elements.rows();
        nodfat=MatrixXi::Zero(non,3);
        nodfmt=MatrixXi::Zero(non,2);
        calculateNodfat();
        totalDof=0;
        for (int j=0;j<non;j++)
        {
            nodfmt(j,1)=nodfat.row(j).sum();
            totalDof=totalDof+nodfmt(j,1);
            if(j<(non-1)){nodfmt(j+1,0)=nodfmt(j,0)+nodfat.row(j).sum();}
        }
    }
    //---------
    {
        Dirichlet=MatrixXd::Zero(0,2);
        dDirichlet=MatrixXd::Zero(0,2);
        Dirichlet0=MatrixXd::Zero(0,2);
        DirichletAll=MatrixXd::Zero(totalDof,3);
    }
}

void AxisSymmetric_2D::assemblyGlobalMatrix()
{
    trip_total.clear();
    trip_total.reserve(15*15*noe);
    KK.setZero();
    KK.resize(totalDof,totalDof);

    for (int ele=0;ele<noe;ele++)
    {
        eleNum=ele;
        int eleMat=elements(eleNum,10);
        int eleType=elements(eleNum,11);
        int nodeCount=elements(eleNum,9);
        int KFunction=material[eleMat-1].KFunction;
        int kFunction=material[eleMat-1].kFunction;
        double ratio=material[eleMat-1].ratio;
        double Cd=material[eleMat-1].Cd;
        double inVal;

        if(eleType==1) //soil elements
        {
            Cd=1;
        }

        if(PVDFlag==true && eleType==2) //set Cd
        {
            Cd=Cd0;
        }

        if(PVDFlag==true && eleMat==2) //set Cd
        {
            Cd=Cd0;
        }

        double e0=material[eleMat-1].e0;
        double porosity=e0/(1+e0);
        Se=porosity*Cf+(BiotCoeff-porosity)*Cs;

        //Interpolation Bulk modulus Ke
        if(KFunction==0) //Constant Bulks modulus
        {
            Ke=material[eleMat-1].KCurve(0,0);

        }
        else if(KFunction==1) //Iterpolation
        {
            inVal=SyyInterpolation(eleNum,0);
            Ke=binarySearch(material[eleMat-1].KCurve,0,1,inVal);
        }

        if(crsFlag==true)
        {
            Ke=binarySearch(crsData,0,4,realCalculationTime*1440.0f); //from minute to day
        }

        ve=material[eleMat-1].poission;
        Ge=3*Ke*(1-2*ve)/2/(1+ve);
        BulkModulus(eleNum,currentStep)=Ke;

        //Interpolation vertical hydraulic conductivity k
        if(kFunction==0)
        {
            kve=material[eleMat-1].kCurve(0,0);
        }
        else if(kFunction==1)
        {
            inVal=SyyInterpolation(eleNum,0);
            kve=binarySearch(material[eleMat-1].kCurve,0,1,inVal);
        }

        if(crsFlag==true)
        {
            kve=binarySearch(crsData,0,5,realCalculationTime*1440.0f); //from minute to day
        }

        HydraulicVertial(eleNum,currentStep)=kve;
        khe=Cd*ratio*kve;
        HydraulicHorizontal(eleNum,currentStep)=khe;

        //Assembly global matrix
        if(nodeCount==6)
        {
            tri6pMatrix(eleNum);
        }
        else if(nodeCount==8)
        {
            quad8pMatrix(eleNum);
        }
        else if(nodeCount==2)
        {
            if(analysisType==1 || analysisType==3)
            {
                line2pMatrix(eleNum);
            }
        }
    }
    KK.setFromTriplets(trip_total.begin(),trip_total.end());

    //Apply boundary condition, global level
    for (int j=0;j<dDirichlet.rows();j++)
    {
        int jj=dDirichlet(j,0);
        dF(jj,0)=dDirichlet(j,1);
        KK.coeffRef(jj,jj)=1;
    }
    KK.prune(0.0);
    KK.makeCompressed();
    trip_total.clear(); //Release memory
}

void AxisSymmetric_2D::solveDirect()
{
    PardisoLU<SparseMatrix<double,RowMajor> > solver;
    solver.compute(KK);
    XX.setZero();
    XX.col(0)=solver.solve(dF);
}

void AxisSymmetric_2D::calculateStress()
{
    //Calculate stress for all steps
    for(int jj=0;jj<ns;jj++)
    {
        step_i=jj;
        X0=X.col(step_i);
        for (int ele=0;ele<noe;ele++)
        {
            eleNum=ele;
            int nodeCount=elements(eleNum,9);
            if(nodeCount==6)
            {
                calculateStressTri6p(eleNum);
            }
            else if(nodeCount==8)
            {
                calculateStressQuad8p(eleNum);
            }
        }
    }

    exportStress();
}

void AxisSymmetric_2D::calculateStress(int calculatedStep)
{
    //calculate incremantl stress
    step_i=calculatedStep;
    for (int ele=0;ele<noe;ele++)
    {
        eleNum=ele;
        int nodeCount=elements(eleNum,9);
        int eleMat=elements(eleNum,10);
        Ke=BulkModulus(eleNum,step_i);
        ve=material[eleMat-1].poission;
        Ge=3*Ke*(1-2*ve)/2/(1+ve);
        if(nodeCount==6)
        {
            calculateStressTri6p(eleNum);
        }
        else if(nodeCount==8)
        {
            calculateStressQuad8p(eleNum);
        }
    }

    //Add to previous step
    if(step_i>0)
    {
        Pore.col(step_i)=Pore.col(step_i)+Pore.col(step_i-1);
        Sxx.col(step_i)=Sxx.col(step_i)+Sxx.col(step_i-1);
        Syy.col(step_i)=Syy.col(step_i)+Syy.col(step_i-1);
        Sxy.col(step_i)=Sxy.col(step_i)+Sxy.col(step_i-1);
    }
}

void AxisSymmetric_2D::averageStress()
{
    QString stressFolder=folderName+"/Average Nodal Stress";
    QDir mDir;
    mDir.remove(stressFolder);
    mDir.mkdir(stressFolder);
    MatrixXd Sxx_node=MatrixXd::Zero(non,ns);
    MatrixXd Syy_node=MatrixXd::Zero(non,ns);
    MatrixXd Sxy_node=MatrixXd::Zero(non,ns);
    MatrixXd nodeStressCount=MatrixXd::Zero(non,1);
    for (int i=0;i<ns;i++)
    {
        for (int j=0;j<noe;j++)
        {
            int nodeCount=elements(j,9);
            for (int jj=0;jj<nodeCount;jj++)
            {
                int nodeIndex=elements(j,jj+1)-1; //Base index is 1
                if(i==0){nodeStressCount(nodeIndex,0)=nodeStressCount(nodeIndex,0)+1;}

                Sxx_node(nodeIndex,i)=Sxx_node(nodeIndex,i)+Sxx(j,i);
                Syy_node(nodeIndex,i)=Syy_node(nodeIndex,i)+Syy(j,i);
                Sxy_node(nodeIndex,i)=Sxy_node(nodeIndex,i)+Sxy(j,i);
            }
        }
    }

    for(int i=0;i<non;i++)
    {
        Sxx_node.row(i)=Sxx_node.row(i)/nodeStressCount(i,0);
        Syy_node.row(i)=Syy_node.row(i)/nodeStressCount(i,0);
        Sxy_node.row(i)=Sxy_node.row(i)/nodeStressCount(i,0);
    }

    fileName=stressFolder+"/"+"Sxx_node.txt";
    exportFile.fileName=fileName.toStdString();
    exportFile.ToFile(Sxx_node);

    fileName=stressFolder+"/"+"Syy_node.txt";
    exportFile.fileName=fileName.toStdString();
    exportFile.ToFile(Syy_node);

    fileName=stressFolder+"/"+"Sxy_node.txt";
    exportFile.fileName=fileName.toStdString();
    exportFile.ToFile(Sxy_node);

}

void AxisSymmetric_2D::exportMaterialsParameters()
{
    QString materialFolder=folderName+"/Material Parameters";
    QDir mDir;
    mDir.remove(materialFolder);
    mDir.mkdir(materialFolder);

    fileName=materialFolder+"/"+"BulkModulus.txt";
    exportFile.fileName=fileName.toStdString();
    exportFile.ToFile(BulkModulus);

    fileName=materialFolder+"/"+"HydraulicVertical.txt";
    exportFile.fileName=fileName.toStdString();
    exportFile.ToFile(HydraulicVertial);

    fileName=materialFolder+"/"+"HydraulicHorizontal.txt";
    exportFile.fileName=fileName.toStdString();
    exportFile.ToFile(HydraulicHorizontal);
}

void AxisSymmetric_2D::getPreviousStep()
{
    X0.setZero();
    for (int jj=0;jj<non;jj++)
    {
        if(nodfmt(jj,1)==3) //u0=0; v0=0;
        {
            X0(nodfmt(jj,0)+2,0)=P(jj,currentStep-1);
        }
    }

    if(currentStep==1)
    {

        for (int jj=0;jj<Dirichlet.rows();jj++)
        {
            int dofIndex=Dirichlet(jj,0);
            double bcValue=Dirichlet(jj,1);
            X0(dofIndex,0)=bcValue;
        }
    }
}

void AxisSymmetric_2D::runPVDs()
{ 
    //Loop over all stage
    resetSolveData();
    timer.start();
    currentStep=0;
    realCalculationTime=0;
    calculationTime.resize(ns,1);

    for (int i=0;i<nos;i++) //loop all steps
    {
        int subStep=stage[i].subStep;
        int stageType=stage[i].stageType;
        int timeStepType=stage[i].timeStepType;
        int gravityCheck=stage[i].gravityLoad;

        double t0=stage[i].t0;
        double t1=stage[i].t1;

        for(int j=0;j<subStep;j++)
        {
            //Calculate real time step
            if(timeStepType==0) //constant time step
            {
                dt=(t1-t0)/subStep;
                realCalculationTime=t0+dt*(j+1);
                dt=dt*86400; //convert from day to seconds
            }
            else
            {
                t1=stage[i].timeStep(j+1,0);
                t0=stage[i].timeStep(j,0);
                dt=(t1-t0);
                dt=dt*86400;
                realCalculationTime=t1;
            }

            //Create boundary condition from v_fix, v_fixy, v_fixh
            Dirichlet0=Dirichlet;
            createBoundaryConditionEachSubStep(i,realCalculationTime);
            createDichletBoundaryCondition();
            createDirichletBoundaryConditionIncrementalForm();

            //First Step, everything is Zero
            X0.setZero();
            if(currentStep==0)
            {
                X0.setZero();
                SyyInterpolation.setZero();

            }
            else if(currentStep==1 && stage[0].stageType==0) //Ignore displacement because of gravity
            {
                X0.setZero();
            }
            else if(currentStep==1 && stage[0].stageType!=0)
            {
                getPreviousStep();
            }
            else if(currentStep>1 && stage[0].stageType!=0)
            {
                getPreviousStep();
            }
            else if(currentStep>1 && stage[0].stageType==0)
            {
                getPreviousStep();
            }
            else
            {
                qDebug()<<"No case"<<endl;
                pauseSystem();
            }

            //Create load vector
            F0=F;
            dF.setZero();
            createLoadEachStep(gravityCheck);
            dF=F-F0;

            //Assembly matrix
            SyyInterpolation.setZero();
            assemblyGlobalMatrix(); //result is KK matrix
            solveDirect(); //solve KK*XX=dF

            if(currentStep==0) //Insitu, undrained
            {
                X.col(currentStep)=XX;
            }
            else //add to incremental
            {
                X.col(currentStep)=XX+X.col(currentStep-1);
                for (int jj=0;jj<Dirichlet.rows();jj++)
                {
                    int dofIndex=Dirichlet(jj,0);
                    double bcValue=Dirichlet(jj,1);
                    X(dofIndex,currentStep)=bcValue;
                }
            }

            //assign to U, V, P arrays
            for (int jj=0;jj<non;jj++)
            {
                U(jj,currentStep)=X(nodfmt(jj,0)+0,currentStep);
                V(jj,currentStep)=X(nodfmt(jj,0)+1,currentStep);
                if(nodfmt(jj,1)==3)
                {
                    P(jj,currentStep)=X(nodfmt(jj,0)+2,currentStep);
                }
            }

            cout<<"Calculation step: "<<currentStep<<" Maximum vertical displacement: "<<V.col(currentStep).minCoeff()<<endl;
            cout<<"Running time    : "<<timer.elapsed()/1000<< " seconds"<<endl;
            calculationTime(currentStep,0)=realCalculationTime*86400;
            currentStep=currentStep+1;
        }
    }
}

void AxisSymmetric_2D::createBoundaryConditionEachSubStep(int currentStage, double realTime)
{
    fixx.resize(0,2);
    fixy.resize(0,2);
    fixh.resize(0,2);
    Fy.resize(0,2);
    Fx.resize(0,2);

    int i=currentStage;
    int nosb=stageBoundary[i].v_boundaryIndex.size();
    double totalTime=stage[i].dt;
    int stageStep=stage[i].subStep;
    int stageType=stage[i].stageType;

    double ratio;
    if(totalTime==0 || stageStep==1)
    {
        ratio=1;
    }
    else
    {
        ratio=double(realTime-stage[i].t0)/double(totalTime);
    }

    for (int j=0;j<nosb;j++)
    {
        int boundaryIndex=stageBoundary[i].v_boundaryIndex[j];
        int boundaryType=boundary[boundaryIndex-1].boundaryType; //Type fixx(0), fixy(1), fixh(2), pressureY(3)
        int assignType=stageBoundary[i].v_assignType[j]; //to location, to node list
        int loadType=boundary[boundaryIndex-1].loadType; //Constant, ramp, and time fucntion

        double x0=stageBoundary[i].v_x0[j];
        double x1=stageBoundary[i].v_x1[j];
        double y0=stageBoundary[i].v_y0[j];
        double y1=stageBoundary[i].v_y1[j];
        double val0=boundary[boundaryIndex-1].val0;
        double val1=boundary[boundaryIndex-1].val1;
        double value;

        bool gradientBool=stageBoundary[i].v_gradientBool[j];
        double aFactor=stageBoundary[i].v_aFactor[j];
        double bFactor=stageBoundary[i].v_bFactor[j];
        double cFactor=stageBoundary[i].v_cFactor[j];

        if(gradientBool==false)
        {
            aFactor=0;
            bFactor=0;
            cFactor=0;
        }

        MatrixXd nodeList=stageBoundary[i].v_nodeList[j];
        MatrixXd foundNode;
        MatrixXd tempMatrix;

        //Loadtype==0 is constant, ==1 is ramp, ==2 is interpolationl
        if(loadType==0)
        {
            value=val1;
        }
        else if (loadType==1)
        {
            if(ratio==1){value=val1;}
            else if(ratio!=1)
            {
                value=ratio*(val1-val0)+val0;
            }
        }
        else if(loadType==2)
        {
            value=binarySearch(boundary[boundaryIndex-1].loadCurve,0,1,realTime);
        }
        //        qDebug()<<"Calculation time: "<<realTime<<endl;
        //        qDebug()<<"Value           : "<<value<<endl;

        //find Node
        matrixLib.findXY(coordinates,x0,x1,y0,y1,foundNode);

        //Assign Type==0 ==> Apply via coordinates
        if(assignType==0)
        {
            tempMatrix.resize(foundNode.rows(),2);
            tempMatrix.setZero();
            tempMatrix.col(0)=foundNode.col(0);
        }
        else if (assignType==1) //Apply via node List
        {
            tempMatrix.resize(nodeList.rows(),2);
            tempMatrix.setZero();
            tempMatrix.col(0)=nodeList.col(0);
        }

        //boundaryType fixX=0, fixY=1, FixH=2, Pressure Y=3; Pressure X=4;
        if(boundaryType==0)
        {
            for (int ii=0;ii<tempMatrix.rows();ii++)
            {
                int index=tempMatrix(ii,0)-1;
                double xCoord=coordinates(index,1);
                double yCoord=coordinates(index,2);
                double valueBoundary=aFactor*xCoord+bFactor*yCoord+cFactor+value;
                tempMatrix(ii,1)=valueBoundary;
            }   
            matrixLib.addMatrixReplace(fixx,tempMatrix,0);
        }
        if(boundaryType==1)
        {
            for (int ii=0;ii<tempMatrix.rows();ii++)
            {
                int index=tempMatrix(ii,0)-1;
                double xCoord=coordinates(index,1);
                double yCoord=coordinates(index,2);
                double valueBoundary=aFactor*xCoord+bFactor*yCoord+cFactor+value;
                tempMatrix(ii,1)=valueBoundary;
            }
            matrixLib.addMatrixReplace(fixy,tempMatrix,0);
        }
        if(boundaryType==2)
        {
            for (int ii=0;ii<tempMatrix.rows();ii++)
            {
                int index=tempMatrix(ii,0)-1;
                double xCoord=coordinates(index,1);
                double yCoord=coordinates(index,2);
                double valueBoundary=aFactor*xCoord+bFactor*yCoord+cFactor+value;
                tempMatrix(ii,1)=valueBoundary;
            }
            matrixLib.addMatrixReplace(fixh,tempMatrix,0);
        }
        if(boundaryType==3) //Not support gradient yet
        {
            transferLineLoadToNodalLoad(tempMatrix,value);
            matrixLib.addMatrix(Fy,tempMatrix,0);
        }
        if(boundaryType==4) //Not support gradient yet
        {
            transferLineLoadToNodalLoadXDirection(tempMatrix,value);
            matrixLib.addMatrix(Fx,tempMatrix,0);
        }
    }

    if(stageType==0) //In-situ, all pore pressure is zero
    {
        fixh.resize(non,2);
        fixh.setZero();
        fixh.col(0)=coordinates.col(0);
    }
}

void AxisSymmetric_2D::createDichletBoundaryCondition()
{
    MatrixXd fixhNew;
    int poreCount=0;
    //Fixx
    for (int i=0;i<fixx.rows();i++)
    {
        int index=fixx(i,0)-1;
        fixx(i,0)=nodfmt(index,0)+0;
    }

    //Fixy
    for (int i=0;i<fixy.rows();i++)
    {
        int index=fixy(i,0)-1;
        fixy(i,0)=nodfmt(index,0)+1;
    }

    //Fixh
    poreCount=countPoreDegreeOfFreedom(fixh,nodfmt);
    fixhNew.resize(poreCount,2);
    poreCount=0;
    for (int j=0;j<fixh.rows();j++)
    {
        int index=fixh(j,0)-1;
        if(nodfmt(index,1)==3)
        {
            fixhNew(poreCount,0)=nodfmt(index,0)+2;
            fixhNew(poreCount,1)=fixh(j,1);
            poreCount=poreCount+1;
        }
    }

    //Local Dirichlet
    Dirichlet.resize(fixx.rows()+fixy.rows()+fixhNew.rows(),2);
    Dirichlet<<fixx,fixy,fixhNew;
}

void AxisSymmetric_2D::createDirichletBoundaryConditionIncrementalForm()
{
    //move from matrix to vector
    dDirichlet=Dirichlet;
    if(currentStep==0)
    {
        dDirichlet=Dirichlet;
    }
    else
    {
        vector<double> oldDirichet;
        vector<double> newDirichlet;
        vector<double>::iterator it;
        oldDirichet.resize(Dirichlet0.rows());
        newDirichlet.resize(Dirichlet.rows());

        VectorXd::Map(&oldDirichet[0],Dirichlet0.rows())=Dirichlet0.col(0);
        VectorXd::Map(&newDirichlet[0],Dirichlet.rows())=Dirichlet.col(0);

        //compare between Dirichlet and Dirichlet 0
        for (int i=0;i<Dirichlet.rows();i++)
        {
            double equationIndex=Dirichlet(i,0);
            it=find(oldDirichet.begin(),oldDirichet.end(),equationIndex);
            if(it!=oldDirichet.end())
            {
                int index=distance(oldDirichet.begin(),it);
                double oldValue=Dirichlet0(index,1);
                dDirichlet(i,1)=Dirichlet(i,1)-oldValue;
            }
        }
    }

    //Create global boundary condition
    DirichletAll=MatrixXd::Zero(totalDof,3);
    createGlobalBoundary(dDirichlet,DirichletAll);
}

void AxisSymmetric_2D::createLoadEachStep(int gravityCheck)
{
    F.setZero();
    for (int i=0;i<Fy.rows();i++)
    {
        int index=Fy(i,0)-1;
        F(nodfmt(index,0)+1,0)=F(nodfmt(index,0)+1,0)+Fy(i,1);
    }

    for (int i=0;i<Fx.rows();i++)
    {
        int index=Fx(i,0)-1;
        F(nodfmt(index,0)+0,0)=F(nodfmt(index,0)+0,0)+Fx(i,1);
    }

    if(gravityCheck==1)
    {
        for(int i=0;i<FyGravity.rows();i++)
        {
            int index=FyGravity(i,0)-1;
            F(nodfmt(index,0)+1,0)=F(nodfmt(index,0)+1,0)+FyGravity(i,1);
        }
    }
}

void AxisSymmetric_2D::compareMatrix(Ref<MatrixXi> matrixA, Ref<MatrixXi> matrixB, Ref<MatrixXi> matrixC)
{
    int cols=matrixA.cols();
    for (int i=0;i<cols;i++)
    {
        if (matrixA(0,i)==matrixB(0,i)){matrixC(0,i)=matrixA(0,i);}
        else if(matrixA(0,i)>matrixB(0,i)){matrixC(0,i)=matrixA(0,i);}
        else{matrixC(0,i)=matrixB(0,i);}
    }
}

void AxisSymmetric_2D::calculateNodfat()
{
    MatrixXi tri6p,tri6d,quad8p,quad8d, line2p;
    tri6p=baseEle.tri6p;
    tri6d=baseEle.tri6d;
    quad8p=baseEle.quad8p;
    quad8d=baseEle.quad8d;
    line2p=baseEle.line2p;

    for (int i=0;i<noe;i++)
    {
        int nodeCount=elements(i,9);
        if (nodeCount==6)
        {
            for (int j=0;j<6;j++)
            {
                MatrixXi matrixA= MatrixXi::Zero(1,3);
                MatrixXi matrixB= MatrixXi::Zero(1,3);
                MatrixXi matrixC= MatrixXi::Zero(1,3);
                int node=elements(i,j+1)-1;
                matrixA=nodfat.row(node);
                matrixB=tri6p.row(j);
                compareMatrix(matrixA,matrixB,matrixC);
                nodfat.row(node)=matrixC;
            }
        }
        else if(nodeCount==8)
        {
            for (int j=0;j<8;j++)
            {
                MatrixXi matrixA= MatrixXi::Zero(1,3);
                MatrixXi matrixB= MatrixXi::Zero(1,3);
                MatrixXi matrixC= MatrixXi::Zero(1,3);
                int node=elements(i,j+1)-1;
                matrixA=nodfat.row(node);
                matrixB=quad8p.row(j);
                compareMatrix(matrixA,matrixB,matrixC);
                nodfat.row(node)=matrixC;
            }
        }
        else if(nodeCount==2)
        {
            for (int j=0;j<2;j++)
            {
                MatrixXi matrixA= MatrixXi::Zero(1,3);
                MatrixXi matrixB= MatrixXi::Zero(1,3);
                MatrixXi matrixC= MatrixXi::Zero(1,3);
                int node=elements(i,j+1)-1;
                matrixA=nodfat.row(node);
                matrixB=line2p.row(j);
                compareMatrix(matrixA,matrixB,matrixC);
                nodfat.row(node)=matrixC;
            }
        }
        else
        {
            cout<<"Error: "<<endl;
            pauseSystem();
        }
    }

}

int AxisSymmetric_2D::countPoreBoundary(Ref<MatrixXd> matrixA, Ref<MatrixXi> nodfmt)
{
    int count=0;
    for (int j=0;j<matrixA.rows();j++)

    {
        int index=matrixA(j,0)-1;
        if(nodfmt(index,1)==3)
        {
            count++;
        }
    }
    return count;
}

void AxisSymmetric_2D::createNewPoreBoundary(Ref<MatrixXd> matrixOld, Ref<MatrixXd> matrixNew, Ref<MatrixXi> nodfmt)
{
    int count=0;
    for(int j=0;j<matrixOld.rows();j++)
    {
        int index=matrixOld(j,0)-1;
        if(nodfmt(index,1)==3)
        {
            matrixNew(count,0)=nodfmt(index,0)+2; //equation index
            matrixNew(count,1)=matrixOld(j,1); //value
            count++;
        }
    }
}

void AxisSymmetric_2D::createGlobalBoundary(Ref<MatrixXd> LocalMatrix, Ref<MatrixXd> GlobalMatrix)
{
    for (int j=0;j<LocalMatrix.rows();j++)
    {
        int equaNum=LocalMatrix(j,0);
        GlobalMatrix(equaNum,0)=equaNum;
        GlobalMatrix(equaNum,1)=1; //Mark as boundary condition
        GlobalMatrix(equaNum,2)=LocalMatrix(j,1);
    }
}

double AxisSymmetric_2D::binarySearch(const Ref<const MatrixXd> searchMatrix, int colSearch,int colResult, double inVal)
{
    //colSearch is searching column
    //colResult is getting value column
    //inVal is input value
    //Return output value

    int left=0;
    int right=0;
    int mid=0;
    double x0=0;
    double x1=0;
    double y0=0;
    double y1=0;
    double outVal;
    right = searchMatrix.rows()-1;
    double max_val;
    double min_val;

    max_val=searchMatrix.col(colSearch).maxCoeff();
    min_val=searchMatrix.col(colSearch).minCoeff();

    if(inVal>=max_val)
    {
        outVal=searchMatrix(searchMatrix.rows()-1,colResult);

    }

    else if(inVal<=min_val)
    {
        outVal=searchMatrix(0,colResult);
    }

    else
    {
        //first step: perform binary search
        while(left<=right)
        {
            mid=(int)((left+right)/2);
            if(inVal<=searchMatrix(mid,colSearch)&&inVal>searchMatrix(mid-1,colSearch))
            {
                break;
            }
            else if (inVal>searchMatrix(mid,colSearch))
            {
                left=mid+1;

            }
            else
            {
                right = mid - 1;
            }

        }
        //second step: linear interpolation
        x0=searchMatrix(mid-1,colSearch);
        x1=searchMatrix(mid,colSearch);
        y0=searchMatrix(mid-1,colResult);
        y1=searchMatrix(mid,colResult);
        outVal=y0+(inVal-x0)*(y1-y0)/(x1-x0);

    }
    return outVal;
}

void AxisSymmetric_2D::runAllAnalysis()
{
    qDebug()<<"Start Running Anlysis"<<endl;
    calculationTime.resize(ns,1);
    //Loop over all stage
    timer.start();
    currentStep=0;
    realCalculationTime=0;

    for (int i=0;i<nos;i++) //loop all steps
    {
        int subStep=stage[i].subStep;
        int stageType=stage[i].stageType;
        int timeStepType=stage[i].timeStepType;
        int gravityCheck=stage[i].gravityLoad;

        double t0=stage[i].t0;
        double t1=stage[i].t1;

        for(int j=0;j<subStep;j++)
        {
            //Calculate real time step
            if(timeStepType==0) //constant time step
            {
                dt=(t1-t0)/subStep;
                realCalculationTime=t0+dt*(j+1);
                dt=dt*86400; //convert from day to seconds
            }
            else
            {
                t1=stage[i].timeStep(j+1,0);
                t0=stage[i].timeStep(j,0);
                dt=(t1-t0);
                dt=dt*86400;
                realCalculationTime=t1;
            }

            if(crsFlag==true)
            {
                realCalculationTime=realCalculationTime*1440.0f;
            }

            qDebug()<<"Time step dt (s) "<<dt<<endl;

            //Create boundary condition from v_fix, v_fixy, v_fixh
            if(currentStep==1 && stage[0].stageType==0) //Ignore insitu step
            {
                //Reset
                Dirichlet.resize(0,2);
                F.setZero();
            }

            Dirichlet0=Dirichlet;
            createBoundaryConditionEachSubStep(i,realCalculationTime);
            createDichletBoundaryCondition();
            createDirichletBoundaryConditionIncrementalForm();

            //First Step, everything is Zero
            X0.setZero();
            if(currentStep==0)
            {
                X0.setZero();
                SyyInterpolation.setZero();

            }
            else if(currentStep==1 && stage[0].stageType==0) //Ignore displacement because of gravity
            {
                X0.setZero();
            }
            else if(currentStep==1 && stage[0].stageType!=0)
            {
                getPreviousStep();
            }
            else if(currentStep>1 && stage[0].stageType!=0)
            {
                getPreviousStep();
            }
            else if(currentStep>1 && stage[0].stageType==0)
            {
                getPreviousStep();
            }
            else
            {
                qDebug()<<"No case"<<endl;
                pauseSystem();
            }

            //Create load vector
            F0=F;
            dF.setZero();
            createLoadEachStep(gravityCheck);
            dF=F-F0;

            //Assembly matrix
            SyyInterpolation.setZero();
            if(currentStep==0)
            {
                SyyInterpolation.setZero(); //step 0, Syy'=0
            }
            else if(currentStep == 1 && stage[0].stageType == 0)
            {
                SyyInterpolation=SyyGravity; //step 1+insitu
            }
            else if(currentStep == 1 && stage[0].stageType != 0)
            {
                SyyInterpolation=Syy.col(currentStep-1); //step1 + without Insitute
            }
            else if(currentStep>1 && stage[0].stageType == 0)
            {
                SyyInterpolation=Syy.col(currentStep-1)+SyyGravity;  //Insitute
            }
            else if(currentStep>1 && stage[0].stageType != 0)
            {
                SyyInterpolation=Syy.col(currentStep-1); //Ignore gravity load
            }
            else
            {
                qDebug()<<" No case for stress interpolation"<<endl;
                pauseSystem();
            }

            assemblyGlobalMatrix(); //result is KK matrix
            solveDirect(); //solve KK*XX=dF

            if(currentStep==0 && stage[0].stageType!=0) //Insitu, undrained
            {
                X.col(currentStep)=XX;
            }
            else if(currentStep==0 && stage[0].stageType==0)
            {
                X.col(currentStep).setZero();
            }
            else //add to incremental
            {
                X.col(currentStep)=XX+X.col(currentStep-1);
                for (int jj=0;jj<Dirichlet.rows();jj++)
                {
                    int dofIndex=Dirichlet(jj,0);
                    double bcValue=Dirichlet(jj,1);
                    X(dofIndex,currentStep)=bcValue;
                }
            }

            //assign to U, V, P arrays
            for (int jj=0;jj<non;jj++)
            {
                U(jj,currentStep)=X(nodfmt(jj,0)+0,currentStep);
                V(jj,currentStep)=X(nodfmt(jj,0)+1,currentStep);
                if(nodfmt(jj,1)==3)
                {
                    P(jj,currentStep)=X(nodfmt(jj,0)+2,currentStep);
                }
            }

            //Calculate stress
            X0=XX;
            calculateStress(currentStep);
            if(currentStep == 0 && stage[0].stageType == 0) //Institu stress
            {
                SyyGravity=Syy.col(currentStep);
                Syy.col(currentStep).setZero();
                Sxx.col(currentStep).setZero();
                Sxy.col(currentStep).setZero();
                Pore.col(currentStep).setZero();
            }
            cout<<"Calculation step: "<<currentStep<<" Maximum vertical displacement: "<<V.col(currentStep).minCoeff()<<endl;
            cout<<"Running time    : "<<timer.elapsed()/1000<< " seconds"<<endl;
            calculationTime(currentStep,0)=realCalculationTime;
            currentStep=currentStep+1;
        }
    }
    emit sendSolverInformation(U,V,P);
}

void AxisSymmetric_2D::exportResults()
{
    if(folderName=="")
    {
        folderName=QFileDialog::getExistingDirectory(Q_NULLPTR,"Chose Save Folder");
    }

    fileName=folderName+"/"+"U.txt";
    exportFile.fileName=fileName.toStdString();
    exportFile.ToFile(U);

    fileName=folderName+"/"+"V.txt";
    exportFile.fileName=fileName.toStdString();
    exportFile.ToFile(V);

    fileName=folderName+"/"+"P.txt";
    exportFile.fileName=fileName.toStdString();
    exportFile.ToFile(P);

    fileName=folderName+"/"+"calculationTime.txt";
    exportFile.fileName=fileName.toStdString();
    exportFile.ToFile(calculationTime);

    if(PVDFlag==true)
    {
        fileName=folderName+"/"+"Panalytical.txt";
        exportFile.fileName=fileName.toStdString();
        exportFile.ToFile(Panalytical);

        fileName=folderName+"/"+"Perror.txt";
        exportFile.fileName=fileName.toStdString();
        exportFile.ToFile(Perror);

        fileName=folderName+"/"+"PerrorEachStep.txt";
        exportFile.fileName=fileName.toStdString();
        exportFile.ToFile(PerrorEachStep);
    }
}

void AxisSymmetric_2D::exportStress()
{
    QDir mDir;
    QString stressFolder=folderName+"/Stress Result";
    mDir.mkdir(stressFolder);

    fileName=stressFolder+"/"+"Sxx.txt";
    exportFile.fileName=fileName.toStdString();
    exportFile.ToFile(Sxx);

    fileName=stressFolder+"/"+"Syy.txt";
    exportFile.fileName=fileName.toStdString();
    exportFile.ToFile(Syy);

    fileName=stressFolder+"/"+"Sxy.txt";
    exportFile.fileName=fileName.toStdString();
    exportFile.ToFile(Sxy);

    fileName=stressFolder+"/"+"Pore.txt";
    exportFile.fileName=fileName.toStdString();
    exportFile.ToFile(Pore);

    fileName=stressFolder+"/"+"SyyGravity.txt";
    exportFile.fileName=fileName.toStdString();
    exportFile.ToFile(SyyGravity);

}

void AxisSymmetric_2D::createGravityLoad()
{
    for (int i=0;i<noe;i++)
    {
        eleNum=i;
        int nodeCount=elements(i,9);
        int eleMat=elements(i,10);
        dense=material[eleMat-1].gf;
        if(nodeCount==6)
        {
            gravityLoadTri6p(eleNum);
        }
        else if(nodeCount==8)
        {
            gravityLoadQuad8p(eleNum);
        }
    }
    qDebug()<<"Gravity Load is created"<<endl;
}

void AxisSymmetric_2D::gravityLoadTri6p(int &eleNum)
{
    int ii=eleNum;
    MatrixXi nodeIndex=MatrixXi::Zero(6,1); //node of elements
    MatrixXd X_coor=MatrixXd::Zero(6,1); //X coordinates
    MatrixXd Y_coor=MatrixXd::Zero(6,1); //X coordinates
    MatrixXi index_v=MatrixXi::Zero(6,1);
    //-----------------------------------------
    //get node coordinates
    for (int j=0;j<6;j++)
    {
        nodeIndex(j,0)=elements(ii,j+1)-1;
        X_coor(j,0)=coordinates(nodeIndex(j,0),1);
        Y_coor(j,0)=coordinates(nodeIndex(j,0),2);
    }

    double x1=X_coor(0,0);
    double x2=X_coor(1,0);
    double x3=X_coor(2,0);
    double x4=X_coor(3,0);
    double x5=X_coor(4,0);
    double x6=X_coor(5,0);

    double y1=Y_coor(0,0);
    double y2=Y_coor(1,0);
    double y3=Y_coor(2,0);

    double x13=x1-x3; double x21=x2-x1; double y31=y3-y1;
    double y12=y1-y2; double y23=y2-y3; double x32=x3-x2;
    double area=0.5*((x2-x1)*(y3-y1)-(y2-y1)*(x3-x1));
    double djay=2*area;
    MatrixXd X_gravity=MatrixXd::Zero(6,1);
    MatrixXd fyEle=MatrixXd::Zero(6,1);
    X_gravity<<0,0,0,x4,x5,x6;
    fyEle=-area*(dense*1.0f/3.0f)*X_gravity;

    FyGravity(nodeIndex(3,0),1)=FyGravity(nodeIndex(3,0),1)+fyEle(3,1);
    FyGravity(nodeIndex(4,0),1)=FyGravity(nodeIndex(4,0),1)+fyEle(4,1);
    FyGravity(nodeIndex(5,0),1)=FyGravity(nodeIndex(5,0),1)+fyEle(5,1);
}

void AxisSymmetric_2D::gravityLoadQuad8p(int &eleNum)
{
    int ii=eleNum;
    double g1, g2, g3, wi;
    dof=8;
    doff=4;
    MatrixXd gaussPoint=MatrixXd::Zero(1,1);
    gaussPoint=gauss.gauss4;
    //-------------------------
    MatrixXi nodeIndex=MatrixXi::Zero(dof,1); //node of elements
    MatrixXd X_coor=MatrixXd::Zero(dof,1); //X coordinates
    MatrixXd Y_coor=MatrixXd::Zero(dof,1); //X coordinates

    //-------------------------------------
    for (int j=0;j<dof;j++)
    {
        nodeIndex(j,0)=elements(ii,j+1)-1;
        X_coor(j,0)=coordinates(nodeIndex(j,0),1);
        Y_coor(j,0)=coordinates(nodeIndex(j,0),2);
    }

    MatrixXd N=MatrixXd::Zero(1,dof);
    MatrixXd dNL1=MatrixXd::Zero(1,dof);
    MatrixXd dNL2=MatrixXd::Zero(1,dof);

    MatrixXd jacobi=MatrixXd::Zero(2,2);
    MatrixXd jacobiLeft=MatrixXd::Zero(2,dof);
    MatrixXd jacobiRight=MatrixXd::Zero(dof,2);

    MatrixXd dN=MatrixXd::Zero(2,dof);
    MatrixXd dNx=MatrixXd::Zero(1,dof);
    MatrixXd dNy=MatrixXd::Zero(1,dof);

    MatrixXd radius=MatrixXd::Zero(1,1);
    MatrixXd Fyelement=MatrixXd::Zero(dof,1);

    //------------------
    //Loop over gauss point
    for (int jj=0;jj<gaussPoint.rows();jj++)
    {
        g1=gaussPoint(jj,0);
        g2=gaussPoint(jj,1);
        g3=1.0f-g1-g2;
        wi=gaussPoint(jj,2);

        N(0,0)=(1.0f/4.0f)*(1.0f-g1)*(1.0f-g2)*(-g1-g2-1.0f);
        N(0,1)=(1.0f/4.0f)*(1.0f+g1)*(1.0f-g2)*(+g1-g2-1.0f);
        N(0,2)=(1.0f/4.0f)*(1.0f+g1)*(1.0f+g2)*(+g1+g2-1.0f);
        N(0,3)=(1.0f/4.0f)*(1.0f-g1)*(1.0f+g2)*(-g1+g2-1.0f);
        N(0,4)=(1.0f/2.0f)*(1.0f-g1*g1)*(1.0f-g2);
        N(0,5)=(1.0f/2.0f)*(1.0f+g1)*(1.0f-g2*g2);
        N(0,6)=(1.0f/2.0f)*(1.0f-g1*g1)*(1.0f+g2);
        N(0,7)=(1.0f/2.0f)*(1.0f-g1)*(1.0f-g2*g2);

        dNL1(0,0)=(1.0f/4.0f)*(1.0f-g2)*(-1.0f*(-g1-g2-1.0f)-1.0f*(1.0f-g1));
        dNL1(0,1)=(1.0f/4.0f)*(1.0f-g2)*(+1.0f*(+g1-g2-1.0f)+1.0f*(1.0f+g1));
        dNL1(0,2)=(1.0f/4.0f)*(1.0f+g2)*(+1.0f*(+g1+g2-1.0f)+1.0f*(1.0f+g1));
        dNL1(0,3)=(1.0f/4.0f)*(1.0f+g2)*(-1.0f*(-g1+g2-1.0f)-1.0f*(1.0f-g1));
        dNL1(0,4)=-(1.0f-g2)*g1;
        dNL1(0,5)=+(1.0f/2.0f)*(1.0f-g2*g2);
        dNL1(0,6)=-(1.0f+g2)*g1;
        dNL1(0,7)=-(1.0f/2.0f)*(1.0f-g2*g2);

        dNL2(0,0)=(1.0f/4.0f)*(1.0f-g1)*(-1.0f*(-g1-g2-1.0f)-1.0f*(1.0f-g2));
        dNL2(0,1)=(1.0f/4.0f)*(1.0f+g1)*(-1.0f*(+g1-g2-1.0f)-1.0f*(1.0f-g2));
        dNL2(0,2)=(1.0f/4.0f)*(1.0f+g1)*(+1.0f*(+g1+g2-1.0f)+1.0f*(1.0f+g2));
        dNL2(0,3)=(1.0f/4.0f)*(1.0f-g1)*(+1.0f*(-g1+g2-1.0f)+1.0f*(1.0f+g2));
        dNL2(0,4)=-(1.0f/2.0f)*(1.0f-g1*g1);
        dNL2(0,5)=-(1.0f+g1)*g2;
        dNL2(0,6)=+(1.0f/2.0f)*(1.0f-g1*g1);
        dNL2(0,7)=-(1.0f-g1)*g2;

        radius=N*X_coor;
        double rade=radius(0,0);

        jacobiLeft<<dNL1,dNL2;
        jacobiRight<<X_coor,Y_coor;
        jacobi=jacobiLeft*jacobiRight;
        double detj;
        detj=abs(jacobi.determinant());

        Fyelement=Fyelement-detj*wi*rade*dense*N.transpose();
    }

    for(int i=0;i<dof;i++)
    {
        FyGravity(nodeIndex(i,0),1)=FyGravity(nodeIndex(i,0),1)+Fyelement(i,0);
    }
}

void AxisSymmetric_2D::transferLineLoadToNodalLoad(Ref<MatrixXd> nodeList, double value)
{
    MatrixXd nodeCoord=MatrixXd::Zero(nodeList.rows(),4);
    nodeCoord.col(0)=nodeList.col(0);
    for(int i=0;i<nodeList.rows();i++)
    {
        int index=nodeList(i,0)-1;
        double xCoord=coordinates(index,1);
        double yCoord=coordinates(index,2);
        nodeCoord(i,1)=xCoord;
        nodeCoord(i,2)=yCoord;
    }
    matrixLib.sortByCol(nodeCoord,1);

    int NumberOfElement=(nodeList.rows()-1)*0.5;
    for (int i=0;i<NumberOfElement;i++)
    {
        int firstNode=i*2;
        int secondNode=i*2+1;
        int thirdNode=i*2+2;

        double r0=nodeCoord(firstNode,1);
        double r1=nodeCoord(thirdNode,1);

        nodeCoord(firstNode,3)=value/6.0f*(r0*r1-r0*r0)+nodeCoord(firstNode,3);
        nodeCoord(secondNode,3)=value/3.0f*(r1*r1-r0*r0)+nodeCoord(secondNode,3);
        nodeCoord(thirdNode,3)=value/6.0f*(r1*r1-r0*r1)+nodeCoord(thirdNode,3);
    }
    matrixLib.sortByCol(nodeCoord,0);
    nodeList.col(1)=nodeCoord.col(3);
}

void AxisSymmetric_2D::transferLineLoadToNodalLoadXDirection(Ref<MatrixXd> nodeList, double value)
{
    MatrixXd nodeCoord=MatrixXd::Zero(nodeList.rows(),4);
    nodeCoord.col(0)=nodeList.col(0);
    for (int i=0;i<nodeList.rows();i++)
    {
        int index=nodeList(i,0)-1;
        double xCoord=coordinates(index,1);
        double yCoord=coordinates(index,2);
        nodeCoord(i,1)=xCoord;
        nodeCoord(i,2)=yCoord;
    }
    matrixLib.sortByCol(nodeCoord,2); //sort by Y direction
    int NumberOfElement=(nodeList.rows()-1)*0.5;
    for (int i=0;i<NumberOfElement;i++)
    {
        int firstNode=i*2;
        int secondNode=i*2+1;
        int thirdNode=i*2+2;

        double r0=nodeCoord(secondNode,1); //middle node
        double y0=nodeCoord(firstNode,2);
        double y1=nodeCoord(thirdNode,2);
        double length=abs(y1-y0);
        nodeCoord(firstNode,3)=r0*value/6.0f*length+nodeCoord(firstNode,3);
        nodeCoord(secondNode,3)=2*r0*value/3.0f*length+nodeCoord(secondNode,3);
        nodeCoord(thirdNode,3)=r0*value/6.0f*length+nodeCoord(thirdNode,3);
    }
    matrixLib.sortByCol(nodeCoord,0);
    nodeList.col(1)=nodeCoord.col(3); //value
}

int AxisSymmetric_2D::countPoreDegreeOfFreedom(const Ref<const MatrixXd> boundaryMatrix, const Ref<const MatrixXi> nodfmt)
{
    int count=0;

    for(int i=0;i<boundaryMatrix.rows();i++)
    {
        int index=boundaryMatrix(i,0)-1;
        if(nodfmt(index,1)==3)
        {
            count=count+1;
        }
    }
    return count;
    cout<<"Number of count :"<<endl;
}

void AxisSymmetric_2D::pauseSystem()
{
    do {
        cout << '\n' << "Press the Enter key to continue.";
    } while (cin.get() != '\n');
}

void AxisSymmetric_2D::tri6pMatrix(int &eleNum)
{
    int ii=eleNum;
    double g1, g2, g3, wi;
    dof=6;
    doff=3;
    MatrixXd gaussPoint=MatrixXd::Zero(1,1);
    gaussPoint=gauss.gauss3;
    //-------------------------
    MatrixXi nodeIndex=MatrixXi::Zero(dof,1); //node of elements
    MatrixXd X_coor=MatrixXd::Zero(dof,1); //X coordinates
    MatrixXd Y_coor=MatrixXd::Zero(dof,1); //X coordinates
    MatrixXd X_coorPore=MatrixXd::Zero(doff,1); //X coordinates
    MatrixXd Y_coorPore=MatrixXd::Zero(doff,1); //X coordinates

    MatrixXd u0l=MatrixXd::Zero(dof,1); //initial x-displacment of elements
    MatrixXd v0l=MatrixXd::Zero(dof,1); //initial y-displacment of elements
    MatrixXd p0l=MatrixXd::Zero(doff,1); //initial pore-pressure of elements

    MatrixXi index_u=MatrixXi::Zero(dof,1);
    MatrixXi index_v=MatrixXi::Zero(dof,1);
    MatrixXi index_p=MatrixXi::Zero(doff,1);
    MatrixXi index_total=MatrixXi::Zero(2*dof+doff,1);

    //element matrices
    MatrixXd Pijl = MatrixXd::Zero(dof,dof); Pijl.setZero();
    MatrixXd Qijl = MatrixXd::Zero(dof,dof); Qijl.setZero();
    MatrixXd Sijl = MatrixXd::Zero(dof,doff); Sijl.setZero();

    MatrixXd Qjil = MatrixXd::Zero(dof,dof); Qijl.setZero();
    MatrixXd Rijl = MatrixXd::Zero(dof,dof); Rijl.setZero();
    MatrixXd Tijl = MatrixXd::Zero(dof,doff); Tijl.setZero();

    MatrixXd Vijl = MatrixXd::Zero(doff,dof); Vijl.setZero();
    MatrixXd Wijl = MatrixXd::Zero(doff,dof); Wijl.setZero();
    MatrixXd Dijl = MatrixXd::Zero(doff,doff); Dijl.setZero();
    MatrixXd Cijl = MatrixXd::Zero(doff,doff); Cijl.setZero();
    MatrixXd Uijl = MatrixXd::Zero(doff,doff); Uijl.setZero();

    MatrixXd QQil = MatrixXd::Zero(doff,1); QQil.setZero();
    MatrixXd Kl=MatrixXd::Zero(2*dof+doff,2*dof+doff); Kl.setZero();
    //-------------------------------------
    for (int j=0;j<dof;j++)
    {
        nodeIndex(j,0)=elements(ii,j+1)-1;
        index_u(j,0)=nodfmt(nodeIndex(j,0),0)+0;
        index_v(j,0)=nodfmt(nodeIndex(j,0),0)+1;
        if(j<doff){index_p(j,0)=nodfmt(nodeIndex(j,0),0)+2;}

        X_coor(j,0)=coordinates(nodeIndex(j,0),1);
        Y_coor(j,0)=coordinates(nodeIndex(j,0),2);

        u0l(j,0)=X0(index_u(j,0),0);
        v0l(j,0)=X0(index_v(j,0),0);
        if(j<doff){p0l(j,0)=X0(index_p(j,0),0);}
    }
    for (int j=0;j<doff;j++)
    {
        X_coorPore(j,0)=X_coor(j,0);
        Y_coorPore(j,0)=Y_coor(j,0);
    }
    index_total<<index_u,index_v,index_p;
    MatrixXd N=MatrixXd::Zero(1,dof);
    MatrixXd dNL1=MatrixXd::Zero(1,dof);
    MatrixXd dNL2=MatrixXd::Zero(1,dof);

    MatrixXd jacobi=MatrixXd::Zero(2,2);
    MatrixXd jacobiLeft=MatrixXd::Zero(2,dof);
    MatrixXd jacobiRight=MatrixXd::Zero(dof,2);

    MatrixXd dN=MatrixXd::Zero(2,dof);
    MatrixXd dNx=MatrixXd::Zero(1,dof);
    MatrixXd dNy=MatrixXd::Zero(1,dof);

    MatrixXd Npore=MatrixXd::Zero(1,doff);
    MatrixXd dNporeL1=MatrixXd::Zero(1,doff);
    MatrixXd dNporeL2=MatrixXd::Zero(1,doff);

    MatrixXd jacobipore=MatrixXd::Zero(2,2);
    MatrixXd jacobiporeLeft=MatrixXd::Zero(2,doff);
    MatrixXd jacobiporeRight=MatrixXd::Zero(doff,2);

    MatrixXd dNpore=MatrixXd::Zero(2,doff);
    MatrixXd dNporex=MatrixXd::Zero(1,doff);
    MatrixXd dNporey=MatrixXd::Zero(1,doff);

    MatrixXd radius=MatrixXd::Zero(1,1);


    //------------------
    //Loop over gauss point
    for (int jj=0;jj<gaussPoint.rows();jj++)
    {
        g1=gaussPoint(jj,0);
        g2=gaussPoint(jj,1);
        g3=1.0f-g1-g2;
        wi=gaussPoint(jj,2);
        N<<(2.0f*g3-1.0f)*g3,(2.0f*g1-1.0f)*g1,(2.0f*g2-1.0f)*g2,4.0f*g1*g3,4.0f*g1*g2,4.0f*g2*g3;
        dNL1<<(1.0f-4.0f*g3),(4.0f*g1-1.0f),0,4.0f*(g3-g1),4.0f*g2,-4.0f*g2;
        dNL2<<(1.0f-4.0f*g3),0,(4.0f*g2-1.0f),-4.0f*g1,4.0f*g1,4.0f*(g3-g2);
        jacobiLeft<<dNL1,dNL2;
        jacobiRight<<X_coor,Y_coor;
        jacobi=jacobiLeft*jacobiRight;
        double detj;
        detj=abs(jacobi.determinant());
        dN=jacobi.inverse()*jacobiLeft;
        dNx=dN.row(0);
        dNy=dN.row(1);

        Npore<<g3,g1,g2;
        dNporeL1<<-1,1,0;
        dNporeL2<<-1,0,1;
        jacobiporeLeft<<dNporeL1,dNporeL2;
        jacobiporeRight<<X_coorPore,Y_coorPore;

        jacobipore=jacobiporeLeft*jacobiporeRight;
        double detjpore;
        detjpore=abs(jacobipore.determinant());
        dNpore=jacobipore.inverse()*jacobiporeLeft;
        dNporex=dNpore.row(0);
        dNporey=dNpore.row(1);
        radius=N*X_coor;
        double rade=radius(0,0);

        Pijl=Pijl+wi*detj*rade*(Ke+4.0f*Ge/3.0f)*(dNx.transpose()*dNx+N.transpose()*N/rade/rade)+
                wi*detj*rade*(Ke-2.0f*Ge/3.0f)*(dNx.transpose()*N/rade+N.transpose()*dNx/rade)+
                wi*detj*rade*Ge*dNy.transpose()*dNy;

        Qijl=Qijl+wi*detj*rade*(Ke-2.0f*Ge/3.0f)*(dNx.transpose()*dNy+N.transpose()*dNy/rade)+
                wi*detj*rade*Ge*dNy.transpose()*dNx;

        Sijl=Sijl-BiotCoeff*detj*wi*rade*(dNx.transpose()*Npore+N.transpose()*Npore/rade);

        Rijl=Rijl+wi*detj*rade*(Ke+4.0f*Ge/3.0f)*(dNy.transpose()*dNy)+
                wi*detj*rade*Ge*(dNx.transpose()*dNx);

        Tijl=Tijl-BiotCoeff*wi*detj*rade*dNy.transpose()*Npore;

        Cijl=Cijl+Se*wi*detj*rade*Npore.transpose()*Npore;
        Dijl=Dijl+wi*detj*rade*(khe/gf)*(dNporex.transpose()*dNporex)+
                wi*detj*rade*(kve/gf)*(dNporey.transpose()*dNporey);

    }
    Qjil=Qijl.transpose();
    Vijl=-Sijl.transpose();
    Wijl=-Tijl.transpose();

    Uijl=Cijl+dt*Dijl;
    QQil=-dt*Dijl*p0l;

    Kl<<Pijl,Qijl,Sijl,
            Qjil,Rijl,Tijl,
            Vijl,Wijl,Uijl;

    for (int j=0;j<doff;j++)
    {
        dF(index_p(j,0),0)=dF(index_p(j,0),0)+QQil(j,0);
    }

    for (int j=0;j<Kl.rows();j=j+1)
    {
        int jj=DirichletAll(index_total(j,0),1);
        double boundaryValue=DirichletAll(index_total(j,0),2);
        if(jj==1)
        {
            for (int k1=0;k1<Kl.cols();k1=k1+1)
            {
                int kk=index_total(k1,0);
                dF(kk,0)=dF(kk,0)-Kl(k1,j)*boundaryValue;
            }
            Kl.row(j).setZero();
            Kl.col(j).setZero();
            Kl(j,j)=1;
        }
    }

    //Push back Kl to trip_total
    for (int j=0;j<Kl.rows();j++)
    {
        for (int jj=0;jj<Kl.cols();jj++)
        {
            int row_i=index_total(j,0);
            int col_j=index_total(jj,0);
            double val_ij=Kl(j,jj);
            trip_total.push_back(Trip(row_i,col_j,val_ij));
        }
    }

}

void AxisSymmetric_2D::quad8pMatrix(int &eleNum)
{
    int ii=eleNum;
    double g1, g2, g3, wi;
    dof=8;
    doff=4;
    MatrixXd gaussPoint=MatrixXd::Zero(1,1);
    gaussPoint=gauss.gauss4;
    //-------------------------
    MatrixXi nodeIndex=MatrixXi::Zero(dof,1); //node of elements
    MatrixXd X_coor=MatrixXd::Zero(dof,1); //X coordinates
    MatrixXd Y_coor=MatrixXd::Zero(dof,1); //X coordinates
    MatrixXd X_coorPore=MatrixXd::Zero(doff,1); //X coordinates
    MatrixXd Y_coorPore=MatrixXd::Zero(doff,1); //X coordinates


    MatrixXd u0l=MatrixXd::Zero(dof,1); //initial x-displacment of elements
    MatrixXd v0l=MatrixXd::Zero(dof,1); //initial y-displacment of elements
    MatrixXd p0l=MatrixXd::Zero(doff,1); //initial pore-pressure of elements

    MatrixXi index_u=MatrixXi::Zero(dof,1);
    MatrixXi index_v=MatrixXi::Zero(dof,1);
    MatrixXi index_p=MatrixXi::Zero(doff,1);
    MatrixXi index_total=MatrixXi::Zero(2*dof+doff,1);

    //element matrices
    MatrixXd Pijl = MatrixXd::Zero(dof,dof); Pijl.setZero();
    MatrixXd Qijl = MatrixXd::Zero(dof,dof); Qijl.setZero();
    MatrixXd Sijl = MatrixXd::Zero(dof,doff); Sijl.setZero();

    MatrixXd Qjil = MatrixXd::Zero(dof,dof); Qijl.setZero();
    MatrixXd Rijl = MatrixXd::Zero(dof,dof); Rijl.setZero();
    MatrixXd Tijl = MatrixXd::Zero(dof,doff); Tijl.setZero();

    MatrixXd Vijl = MatrixXd::Zero(doff,dof); Vijl.setZero();
    MatrixXd Wijl = MatrixXd::Zero(doff,dof); Wijl.setZero();
    MatrixXd Dijl = MatrixXd::Zero(doff,doff); Dijl.setZero();
    MatrixXd Cijl = MatrixXd::Zero(doff,doff); Cijl.setZero();
    MatrixXd Uijl = MatrixXd::Zero(doff,doff); Uijl.setZero();

    MatrixXd QQil = MatrixXd::Zero(doff,1); QQil.setZero();
    MatrixXd Kl=MatrixXd::Zero(2*dof+doff,2*dof+doff); Kl.setZero();
    //-------------------------------------
    for (int j=0;j<dof;j++)
    {
        nodeIndex(j,0)=elements(ii,j+1)-1;
        index_u(j,0)=nodfmt(nodeIndex(j,0),0)+0;
        index_v(j,0)=nodfmt(nodeIndex(j,0),0)+1;
        if(j<doff){index_p(j,0)=nodfmt(nodeIndex(j,0),0)+2;}

        X_coor(j,0)=coordinates(nodeIndex(j,0),1);
        Y_coor(j,0)=coordinates(nodeIndex(j,0),2);

        u0l(j,0)=X0(index_u(j,0),0);
        v0l(j,0)=X0(index_v(j,0),0);
        if(j<doff){p0l(j,0)=X0(index_p(j,0),0);}
    }
    for (int j=0;j<doff;j++)
    {
        X_coorPore(j,0)=X_coor(j,0);
        Y_coorPore(j,0)=Y_coor(j,0);
    }
    index_total<<index_u,index_v,index_p;
    MatrixXd N=MatrixXd::Zero(1,dof);
    MatrixXd dNL1=MatrixXd::Zero(1,dof);
    MatrixXd dNL2=MatrixXd::Zero(1,dof);

    MatrixXd jacobi=MatrixXd::Zero(2,2);
    MatrixXd jacobiLeft=MatrixXd::Zero(2,dof);
    MatrixXd jacobiRight=MatrixXd::Zero(dof,2);

    MatrixXd dN=MatrixXd::Zero(2,dof);
    MatrixXd dNx=MatrixXd::Zero(1,dof);
    MatrixXd dNy=MatrixXd::Zero(1,dof);

    MatrixXd Npore=MatrixXd::Zero(1,doff);
    MatrixXd dNporeL1=MatrixXd::Zero(1,doff);
    MatrixXd dNporeL2=MatrixXd::Zero(1,doff);

    MatrixXd jacobipore=MatrixXd::Zero(2,2);
    MatrixXd jacobiporeLeft=MatrixXd::Zero(2,doff);
    MatrixXd jacobiporeRight=MatrixXd::Zero(doff,2);

    MatrixXd dNpore=MatrixXd::Zero(2,doff);
    MatrixXd dNporex=MatrixXd::Zero(1,doff);
    MatrixXd dNporey=MatrixXd::Zero(1,doff);

    MatrixXd radius=MatrixXd::Zero(1,1);

    //------------------
    //Loop over gauss point
    for (int jj=0;jj<gaussPoint.rows();jj++)
    {
        g1=gaussPoint(jj,0);
        g2=gaussPoint(jj,1);
        g3=1.0f-g1-g2;
        wi=gaussPoint(jj,2);

        N(0,0)=(1.0f/4.0f)*(1.0f-g1)*(1.0f-g2)*(-g1-g2-1.0f);
        N(0,1)=(1.0f/4.0f)*(1.0f+g1)*(1.0f-g2)*(+g1-g2-1.0f);
        N(0,2)=(1.0f/4.0f)*(1.0f+g1)*(1.0f+g2)*(+g1+g2-1.0f);
        N(0,3)=(1.0f/4.0f)*(1.0f-g1)*(1.0f+g2)*(-g1+g2-1.0f);
        N(0,4)=(1.0f/2.0f)*(1.0f-g1*g1)*(1.0f-g2);
        N(0,5)=(1.0f/2.0f)*(1.0f+g1)*(1.0f-g2*g2);
        N(0,6)=(1.0f/2.0f)*(1.0f-g1*g1)*(1.0f+g2);
        N(0,7)=(1.0f/2.0f)*(1.0f-g1)*(1.0f-g2*g2);

        dNL1(0,0)=(1.0f/4.0f)*(1.0f-g2)*(-1.0f*(-g1-g2-1.0f)-1.0f*(1.0f-g1));
        dNL1(0,1)=(1.0f/4.0f)*(1.0f-g2)*(+1.0f*(+g1-g2-1.0f)+1.0f*(1.0f+g1));
        dNL1(0,2)=(1.0f/4.0f)*(1.0f+g2)*(+1.0f*(+g1+g2-1.0f)+1.0f*(1.0f+g1));
        dNL1(0,3)=(1.0f/4.0f)*(1.0f+g2)*(-1.0f*(-g1+g2-1.0f)-1.0f*(1.0f-g1));
        dNL1(0,4)=-(1.0f-g2)*g1;
        dNL1(0,5)=+(1.0f/2.0f)*(1.0f-g2*g2);
        dNL1(0,6)=-(1.0f+g2)*g1;
        dNL1(0,7)=-(1.0f/2.0f)*(1.0f-g2*g2);

        dNL2(0,0)=(1.0f/4.0f)*(1.0f-g1)*(-1.0f*(-g1-g2-1.0f)-1.0f*(1.0f-g2));
        dNL2(0,1)=(1.0f/4.0f)*(1.0f+g1)*(-1.0f*(+g1-g2-1.0f)-1.0f*(1.0f-g2));
        dNL2(0,2)=(1.0f/4.0f)*(1.0f+g1)*(+1.0f*(+g1+g2-1.0f)+1.0f*(1.0f+g2));
        dNL2(0,3)=(1.0f/4.0f)*(1.0f-g1)*(+1.0f*(-g1+g2-1.0f)+1.0f*(1.0f+g2));
        dNL2(0,4)=-(1.0f/2.0f)*(1.0f-g1*g1);
        dNL2(0,5)=-(1.0f+g1)*g2;
        dNL2(0,6)=+(1.0f/2.0f)*(1.0f-g1*g1);
        dNL2(0,7)=-(1.0f-g1)*g2;

        radius=N*X_coor;
        double rade=radius(0,0);

        jacobiLeft<<dNL1,dNL2;
        jacobiRight<<X_coor,Y_coor;
        jacobi=jacobiLeft*jacobiRight;
        double detj;
        detj=abs(jacobi.determinant());
        dN=jacobi.inverse()*jacobiLeft;
        dNx=dN.row(0);
        dNy=dN.row(1);

        Npore(0,0)=(1.0f/4.0f)*(1.0f-g1)*(1.0f-g2);
        Npore(0,1)=(1.0f/4.0f)*(1.0f+g1)*(1.0f-g2);
        Npore(0,2)=(1.0f/4.0f)*(1.0f+g1)*(1.0f+g2);
        Npore(0,3)=(1.0f/4.0f)*(1.0f-g1)*(1.0f+g2);

        dNporeL1(0,0)=-(1.0f/4.0f)*(1.0f-g2);
        dNporeL1(0,1)=+(1.0f/4.0f)*(1.0f-g2);
        dNporeL1(0,2)=+(1.0f/4.0f)*(1.0f+g2);
        dNporeL1(0,3)=-(1.0f/4.0f)*(1.0f+g2);

        dNporeL2(0,0)=-(1.0f/4.0f)*(1.0f-g1);
        dNporeL2(0,1)=-(1.0f/4.0f)*(1.0f+g1);
        dNporeL2(0,2)=+(1.0f/4.0f)*(1.0f+g1);
        dNporeL2(0,3)=+(1.0f/4.0f)*(1.0f-g1);

        jacobiporeLeft<<dNporeL1,dNporeL2;
        jacobiporeRight<<X_coorPore,Y_coorPore;
        jacobipore=jacobiporeLeft*jacobiporeRight;
        double detjpore;
        detjpore=abs(jacobipore.determinant());
        dNpore=jacobipore.inverse()*jacobiporeLeft;
        dNporex=dNpore.row(0);
        dNporey=dNpore.row(1);

        Pijl=Pijl+wi*detj*rade*(Ke+4.0f*Ge/3.0f)*(dNx.transpose()*dNx+N.transpose()*N/rade/rade)+
                wi*detj*rade*(Ke-2.0f*Ge/3.0f)*(dNx.transpose()*N/rade+N.transpose()*dNx/rade)+
                wi*detj*rade*Ge*dNy.transpose()*dNy;
        Qijl=Qijl+wi*detj*rade*(Ke-2.0f*Ge/3.0f)*(dNx.transpose()*dNy+N.transpose()*dNy/rade)+
                wi*detj*rade*Ge*dNy.transpose()*dNx;
        Sijl=Sijl-BiotCoeff*detj*wi*rade*(dNx.transpose()*Npore+N.transpose()*Npore/rade);

        Rijl=Rijl+wi*detj*rade*(Ke+4.0f*Ge/3.0f)*(dNy.transpose()*dNy)+
                wi*detj*rade*Ge*(dNx.transpose()*dNx);
        Tijl=Tijl-BiotCoeff*wi*detj*rade*dNy.transpose()*Npore;

        Cijl=Cijl+Se*wi*detj*rade*Npore.transpose()*Npore;
        Dijl=Dijl+wi*detj*rade*(khe/gf)*(dNporex.transpose()*dNporex)+
                wi*detj*rade*(kve/gf)*(dNporey.transpose()*dNporey);
    }
    Qjil=Qijl.transpose();
    Vijl=-Sijl.transpose();
    Wijl=-Tijl.transpose();

    Uijl=Cijl+dt*Dijl;
    QQil=-dt*Dijl*p0l;

    Kl<<Pijl,Qijl,Sijl,
            Qjil,Rijl,Tijl,
            Vijl,Wijl,Uijl;


    for (int j=0;j<doff;j++)
    {
        dF(index_p(j,0),0)=dF(index_p(j,0),0)+QQil(j,0);
    }

    for (int j=0;j<Kl.rows();j=j+1)
    {
        int jj=DirichletAll(index_total(j,0),1);
        double boundaryValue=DirichletAll(index_total(j,0),2);
        if(jj==1)
        {
            for (int k1=0;k1<Kl.cols();k1=k1+1)
            {
                int kk=index_total(k1,0);
                dF(kk,0)=dF(kk,0)-Kl(k1,j)*boundaryValue;
            }
            Kl.row(j).setZero();
            Kl.col(j).setZero();
            Kl(j,j)=1;
        }
    }

    //Push back Kl to trip_total
    for (int j=0;j<Kl.rows();j++)
    {
        for (int jj=0;jj<Kl.cols();jj++)
        {
            int row_i=index_total(j,0);
            int col_j=index_total(jj,0);
            double val_ij=Kl(j,jj);
            trip_total.push_back(Trip(row_i,col_j,val_ij));
        }
    }

}

void AxisSymmetric_2D::line2pMatrix(int &eleNum)
{
    int ii=eleNum;
    doff=2;
    MatrixXi nodeIndex=MatrixXi::Zero(doff,1); //node of elements
    MatrixXd X_coorPore=MatrixXd::Zero(doff,1); //X coordinates
    MatrixXd Y_coorPore=MatrixXd::Zero(doff,1); //X coordinates
    MatrixXd p0l=MatrixXd::Zero(doff,1); //initial pore-pressure of elements
    double length=0;

    MatrixXi index_p=MatrixXi::Zero(doff,1);
    MatrixXi index_total=MatrixXi::Zero(doff,1);

    //element matrices
    MatrixXd Dijl = MatrixXd::Zero(doff,doff);
    MatrixXd Uijl = MatrixXd::Zero(doff,doff);
    MatrixXd QQil = MatrixXd::Zero(doff,1);
    MatrixXd Kl=MatrixXd::Zero(doff,doff);

    for (int j=0;j<doff;j++)
    {
        nodeIndex(j,0)=elements(ii,j+1)-1;
        index_p(j,0)=nodfmt(nodeIndex(j,0),0)+2;
        X_coorPore(j,0)=coordinates(nodeIndex(j,0),1);
        Y_coorPore(j,0)=coordinates(nodeIndex(j,0),2);
        p0l(j,0)=X0(index_p(j,0),0);
    }

    length=(X_coorPore(1,0)-X_coorPore(0,0))*(X_coorPore(1,0)-X_coorPore(0,0))+(Y_coorPore(1,0)-Y_coorPore(0,0))*(Y_coorPore(1,0)-Y_coorPore(0,0));
    length=sqrt(length);
    index_total=index_p;

    MatrixXd Ae=MatrixXd::Zero(2,2);
    Ae<<1.0f,-1.0f,
            -1.0f,1.0f;
    double tempValue=(A1D*k1D)/(gf*length);
    Dijl=tempValue*Ae;

    //Time incremental form
    Uijl=dt*Dijl;
    QQil=-dt*Dijl*p0l;
    Kl<<Uijl;

    for (int j=0;j<doff;j++)
    {
        dF(index_p(j,0),0)=dF(index_p(j,0),0)+QQil(j,0);
    }

    for (int j=0;j<Kl.rows();j=j+1)
    {
        int jj=DirichletAll(index_total(j,0),1);
        double boundaryValue=DirichletAll(index_total(j,0),2);
        if(jj==1)
        {
            for (int k1=0;k1<Kl.cols();k1=k1+1)
            {
                int kk=index_total(k1,0);
                dF(kk,0)=dF(kk,0)-Kl(k1,j)*boundaryValue;
            }
            Kl.row(j).setZero();
            Kl.col(j).setZero();
            Kl(j,j)=1;
        }
    }

    //Push back Kl to trip_total
    for (int j=0;j<Kl.rows();j++)
    {
        for (int jj=0;jj<Kl.cols();jj++)
        {
            int row_i=index_total(j,0);
            int col_j=index_total(jj,0);
            double val_ij=Kl(j,jj);
            trip_total.push_back(Trip(row_i,col_j,val_ij));
        }
    }
}

void AxisSymmetric_2D::calculateStressTri6p(int &eleNum)
{
    int ii=eleNum;
    double g1, g2, g3, wi;
    dof=6;
    doff=3;
    MatrixXd gaussPoint=MatrixXd::Zero(1,1);
    gaussPoint=gauss.gauss3;
    //-------------------------
    MatrixXi nodeIndex=MatrixXi::Zero(dof,1); //node of elements
    MatrixXd X_coor=MatrixXd::Zero(dof,1); //X coordinates
    MatrixXd Y_coor=MatrixXd::Zero(dof,1); //X coordinates
    MatrixXd X_coorPore=MatrixXd::Zero(doff,1); //X coordinates
    MatrixXd Y_coorPore=MatrixXd::Zero(doff,1); //X coordinates

    MatrixXd u0l=MatrixXd::Zero(dof,1); //initial x-displacment of elements
    MatrixXd v0l=MatrixXd::Zero(dof,1); //initial y-displacment of elements
    MatrixXd p0l=MatrixXd::Zero(doff,1); //initial pore-pressure of elements

    MatrixXi index_u=MatrixXi::Zero(dof,1);
    MatrixXi index_v=MatrixXi::Zero(dof,1);
    MatrixXi index_p=MatrixXi::Zero(doff,1);
    MatrixXi index_total=MatrixXi::Zero(2*dof+doff,1);

    //-------------------------------------
    MatrixXd Pore_e=MatrixXd::Zero(1,1);
    MatrixXd Sxx_e=MatrixXd::Zero(1,1);
    MatrixXd Syy_e=MatrixXd::Zero(1,1);
    MatrixXd Sxy_e=MatrixXd::Zero(1,1);
    //-------------------------

    for (int j=0;j<dof;j++)
    {
        nodeIndex(j,0)=elements(ii,j+1)-1;
        index_u(j,0)=nodfmt(nodeIndex(j,0),0)+0;
        index_v(j,0)=nodfmt(nodeIndex(j,0),0)+1;
        if(j<doff){index_p(j,0)=nodfmt(nodeIndex(j,0),0)+2;}

        X_coor(j,0)=coordinates(nodeIndex(j,0),1);
        Y_coor(j,0)=coordinates(nodeIndex(j,0),2);

        u0l(j,0)=X0(index_u(j,0),0);
        v0l(j,0)=X0(index_v(j,0),0);
        if(j<doff){p0l(j,0)=X0(index_p(j,0),0);}
    }
    for (int j=0;j<doff;j++)
    {
        X_coorPore(j,0)=X_coor(j,0);
        Y_coorPore(j,0)=Y_coor(j,0);
    }
    index_total<<index_u,index_v,index_p;
    MatrixXd N=MatrixXd::Zero(1,dof);
    MatrixXd dNL1=MatrixXd::Zero(1,dof);
    MatrixXd dNL2=MatrixXd::Zero(1,dof);

    MatrixXd jacobi=MatrixXd::Zero(2,2);
    MatrixXd jacobiLeft=MatrixXd::Zero(2,dof);
    MatrixXd jacobiRight=MatrixXd::Zero(dof,2);

    MatrixXd dN=MatrixXd::Zero(2,dof);
    MatrixXd dNx=MatrixXd::Zero(1,dof);
    MatrixXd dNy=MatrixXd::Zero(1,dof);

    MatrixXd Npore=MatrixXd::Zero(1,doff);
    MatrixXd radius=MatrixXd::Zero(1,1);
    //------------------
    //Loop over gauss point
    for (int jj=0;jj<gaussPoint.rows();jj++)
    {
        g1=gaussPoint(jj,0);
        g2=gaussPoint(jj,1);
        g3=1.0f-g1-g2;
        wi=gaussPoint(jj,2);
        N<<(2.0f*g3-1.0f)*g3,(2.0f*g1-1.0f)*g1,(2.0f*g2-1.0f)*g2,4.0f*g1*g3,4.0f*g1*g2,4.0f*g2*g3;
        dNL1<<(1.0f-4.0f*g3),(4.0f*g1-1.0f),0,4.0f*(g3-g1),4.0f*g2,-4.0f*g2;
        dNL2<<(1.0f-4.0f*g3),0,(4.0f*g2-1.0f),-4.0f*g1,4.0f*g1,4.0f*(g3-g2);
        jacobiLeft<<dNL1,dNL2;
        jacobiRight<<X_coor,Y_coor;
        jacobi=jacobiLeft*jacobiRight;
        double detj;
        detj=abs(jacobi.determinant());
        dN=jacobi.inverse()*jacobiLeft;
        dNx=dN.row(0);
        dNy=dN.row(1);

        Npore<<g3,g1,g2;

        radius=N*X_coor;
        double rade=radius(0,0);

        Pore_e=Pore_e+Npore*p0l;
        Sxx_e=Sxx_e-(Ke+4.0*Ge/3.0)*dNx*u0l-(Ke-2.0*Ge/3.0)*(N*u0l/rade+dNy*v0l);
        Syy_e=Syy_e-(Ke+4.0*Ge/3.0)*dNy*v0l-(Ke-2.0*Ge/3.0)*(N*u0l/rade+dNx*u0l);
        Sxy_e=Sxy_e-Ge*dNy*u0l-Ge*dNx*v0l;
    }
    if(step_i==0)
    {
        Pore(ii,step_i)=Pore_e(0,0)/gaussPoint.rows();
        Sxx(ii,step_i)=Sxx_e(0,0)/gaussPoint.rows();
        Syy(ii,step_i)=Syy_e(0,0)/gaussPoint.rows();
        Sxy(ii,step_i)=Sxy_e(0,0)/gaussPoint.rows();
    }
    else
    {
        Pore(ii,step_i)=Pore_e(0,0)/gaussPoint.rows();
        Sxx(ii,step_i)=Sxx_e(0,0)/gaussPoint.rows();
        Syy(ii,step_i)=Syy_e(0,0)/gaussPoint.rows();
        Sxy(ii,step_i)=Sxy_e(0,0)/gaussPoint.rows();
    }

}

void AxisSymmetric_2D::calculateStressQuad8p(int &eleNum)
{
    int ii=eleNum;
    double g1, g2, g3, wi;
    dof=8;
    doff=4;
    MatrixXd gaussPoint=MatrixXd::Zero(1,1);
    gaussPoint=gauss.gauss4;
    //-------------------------
    MatrixXi nodeIndex=MatrixXi::Zero(dof,1); //node of elements
    MatrixXd X_coor=MatrixXd::Zero(dof,1); //X coordinates
    MatrixXd Y_coor=MatrixXd::Zero(dof,1); //X coordinates
    MatrixXd X_coorPore=MatrixXd::Zero(doff,1); //X coordinates
    MatrixXd Y_coorPore=MatrixXd::Zero(doff,1); //X coordinates

    MatrixXd u0l=MatrixXd::Zero(dof,1); //initial x-displacment of elements
    MatrixXd v0l=MatrixXd::Zero(dof,1); //initial y-displacment of elements
    MatrixXd p0l=MatrixXd::Zero(doff,1); //initial pore-pressure of elements

    MatrixXi index_u=MatrixXi::Zero(dof,1);
    MatrixXi index_v=MatrixXi::Zero(dof,1);
    MatrixXi index_p=MatrixXi::Zero(doff,1);
    MatrixXi index_total=MatrixXi::Zero(2*dof+doff,1);

    //-------------------------------------
    MatrixXd Pore_e=MatrixXd::Zero(1,1);
    MatrixXd Sxx_e=MatrixXd::Zero(1,1);
    MatrixXd Syy_e=MatrixXd::Zero(1,1);
    MatrixXd Sxy_e=MatrixXd::Zero(1,1);
    //-------------------------------------
    for (int j=0;j<dof;j++)
    {
        nodeIndex(j,0)=elements(ii,j+1)-1;
        index_u(j,0)=nodfmt(nodeIndex(j,0),0)+0;
        index_v(j,0)=nodfmt(nodeIndex(j,0),0)+1;
        if(j<doff){index_p(j,0)=nodfmt(nodeIndex(j,0),0)+2;}

        X_coor(j,0)=coordinates(nodeIndex(j,0),1);
        Y_coor(j,0)=coordinates(nodeIndex(j,0),2);

        u0l(j,0)=X0(index_u(j,0),0);
        v0l(j,0)=X0(index_v(j,0),0);
        if(j<doff){p0l(j,0)=X0(index_p(j,0),0);}
    }
    for (int j=0;j<doff;j++)
    {
        X_coorPore(j,0)=X_coor(j,0);
        Y_coorPore(j,0)=Y_coor(j,0);
    }
    index_total<<index_u,index_v,index_p;
    MatrixXd N=MatrixXd::Zero(1,dof);
    MatrixXd dNL1=MatrixXd::Zero(1,dof);
    MatrixXd dNL2=MatrixXd::Zero(1,dof);

    MatrixXd jacobi=MatrixXd::Zero(2,2);
    MatrixXd jacobiLeft=MatrixXd::Zero(2,dof);
    MatrixXd jacobiRight=MatrixXd::Zero(dof,2);

    MatrixXd dN=MatrixXd::Zero(2,dof);
    MatrixXd dNx=MatrixXd::Zero(1,dof);
    MatrixXd dNy=MatrixXd::Zero(1,dof);

    MatrixXd Npore=MatrixXd::Zero(1,doff);
    MatrixXd radius=MatrixXd::Zero(1,1);

    //------------------
    //Loop over gauss point
    for (int jj=0;jj<gaussPoint.rows();jj++)
    {
        g1=gaussPoint(jj,0);
        g2=gaussPoint(jj,1);
        g3=1.0f-g1-g2;
        wi=gaussPoint(jj,2);

        N(0,0)=(1.0f/4.0f)*(1.0f-g1)*(1.0f-g2)*(-g1-g2-1.0f);
        N(0,1)=(1.0f/4.0f)*(1.0f+g1)*(1.0f-g2)*(+g1-g2-1.0f);
        N(0,2)=(1.0f/4.0f)*(1.0f+g1)*(1.0f+g2)*(+g1+g2-1.0f);
        N(0,3)=(1.0f/4.0f)*(1.0f-g1)*(1.0f+g2)*(-g1+g2-1.0f);
        N(0,4)=(1.0f/2.0f)*(1.0f-g1*g1)*(1.0f-g2);
        N(0,5)=(1.0f/2.0f)*(1.0f+g1)*(1.0f-g2*g2);
        N(0,6)=(1.0f/2.0f)*(1.0f-g1*g1)*(1.0f+g2);
        N(0,7)=(1.0f/2.0f)*(1.0f-g1)*(1.0f-g2*g2);

        dNL1(0,0)=(1.0f/4.0f)*(1.0f-g2)*(-1.0f*(-g1-g2-1.0f)-1.0f*(1.0f-g1));
        dNL1(0,1)=(1.0f/4.0f)*(1.0f-g2)*(+1.0f*(+g1-g2-1.0f)+1.0f*(1.0f+g1));
        dNL1(0,2)=(1.0f/4.0f)*(1.0f+g2)*(+1.0f*(+g1+g2-1.0f)+1.0f*(1.0f+g1));
        dNL1(0,3)=(1.0f/4.0f)*(1.0f+g2)*(-1.0f*(-g1+g2-1.0f)-1.0f*(1.0f-g1));
        dNL1(0,4)=-(1.0f-g2)*g1;
        dNL1(0,5)=+(1.0f/2.0f)*(1.0f-g2*g2);
        dNL1(0,6)=-(1.0f+g2)*g1;
        dNL1(0,7)=-(1.0f/2.0f)*(1.0f-g2*g2);

        dNL2(0,0)=(1.0f/4.0f)*(1.0f-g1)*(-1.0f*(-g1-g2-1.0f)-1.0f*(1.0f-g2));
        dNL2(0,1)=(1.0f/4.0f)*(1.0f+g1)*(-1.0f*(+g1-g2-1.0f)-1.0f*(1.0f-g2));
        dNL2(0,2)=(1.0f/4.0f)*(1.0f+g1)*(+1.0f*(+g1+g2-1.0f)+1.0f*(1.0f+g2));
        dNL2(0,3)=(1.0f/4.0f)*(1.0f-g1)*(+1.0f*(-g1+g2-1.0f)+1.0f*(1.0f+g2));
        dNL2(0,4)=-(1.0f/2.0f)*(1.0f-g1*g1);
        dNL2(0,5)=-(1.0f+g1)*g2;
        dNL2(0,6)=+(1.0f/2.0f)*(1.0f-g1*g1);
        dNL2(0,7)=-(1.0f-g1)*g2;

        radius=N*X_coor;
        double rade=radius(0,0);

        jacobiLeft<<dNL1,dNL2;
        jacobiRight<<X_coor,Y_coor;
        jacobi=jacobiLeft*jacobiRight;
        double detj;
        detj=abs(jacobi.determinant());
        dN=jacobi.inverse()*jacobiLeft;
        dNx=dN.row(0);
        dNy=dN.row(1);

        Npore(0,0)=(1.0f/4.0f)*(1.0f-g1)*(1.0f-g2);
        Npore(0,1)=(1.0f/4.0f)*(1.0f+g1)*(1.0f-g2);
        Npore(0,2)=(1.0f/4.0f)*(1.0f+g1)*(1.0f+g2);
        Npore(0,3)=(1.0f/4.0f)*(1.0f-g1)*(1.0f+g2);

        Pore_e=Pore_e+Npore*p0l;
        Sxx_e=Sxx_e-(Ke+4.0*Ge/3.0)*dNx*u0l-(Ke-2.0*Ge/3.0)*(N*u0l/rade+dNy*v0l);
        Syy_e=Syy_e-(Ke+4.0*Ge/3.0)*dNy*v0l-(Ke-2.0*Ge/3.0)*(N*u0l/rade+dNx*u0l);
        Sxy_e=Sxy_e-Ge*dNy*u0l-Ge*dNx*v0l;

    }
    Pore(ii,step_i)=Pore_e(0,0)/gaussPoint.rows();
    Sxx(ii,step_i)=Sxx_e(0,0)/gaussPoint.rows();
    Syy(ii,step_i)=Syy_e(0,0)/gaussPoint.rows();
    Sxy(ii,step_i)=Sxy_e(0,0)/gaussPoint.rows();

}

double AxisSymmetric_2D::calculateWatchList(int index, int step) //base 1, and base 1
{
    int i=index-1; //index of watchList
    int kk=step-1;
    int watchType=watch[i].watchType;
    double x0=watch[i].x0;
    double x1=watch[i].x1;
    double y0=watch[i].y0;
    double y1=watch[i].y1;
    bool averageBool=watch[i].averageBool;

    MatrixXd foundNode;

    if(watchType<=2)
    {
        matrixLib.findXY(coordinates,x0,x1,y0,y1,foundNode);
    }
    else
    {
        matrixLib.findElementXY(coordinates,elements,x0,x1,y0,y1,foundNode);
    }

    double totalResults=0;
    int nodeCount=0;

    for (int jj=0;jj<foundNode.rows();jj++)
    {
        int nodeIndex=foundNode(jj,0)-1;
        if(watchType==0) //X direction
        {
            totalResults=totalResults+U(nodeIndex,kk);
            nodeCount=nodeCount+1;
        }
        else if(watchType==1) //Y displacement
        {
            totalResults=totalResults+V(nodeIndex,kk);
            nodeCount=nodeCount+1;
        }
        else if(watchType==2) //Pore Pressure
        {
            if(nodfmt(nodeIndex,1)==3)
            {
                totalResults=totalResults+P(nodeIndex,kk);
                nodeCount=nodeCount+1;
            }
        }
        else if (watchType==3) //Effective Sxx
        {
            totalResults=totalResults+Sxx(nodeIndex,kk);
            nodeCount=nodeCount+1;
        }
        else if (watchType==4) //Effective Syy
        {
            totalResults=totalResults+Syy(nodeIndex,kk);
            nodeCount=nodeCount+1;
        }
        else if (watchType==5) //Effectiive Sxy
        {
            totalResults=totalResults+Sxy(nodeIndex,kk);
            nodeCount=nodeCount+1;
        }
        else if (watchType==6) //Pore pressure
        {
            totalResults=totalResults+Pore(nodeIndex,kk);
            nodeCount=nodeCount+1;
        }
        else if (watchType==7) //total Sxx
        {
            totalResults=totalResults+Sxx(nodeIndex,kk)+Pore(nodeIndex,kk);
            nodeCount=nodeCount+1;
        }
        else if (watchType==8) //Total Syy
        {
            totalResults=totalResults+Syy(nodeIndex,kk)+Pore(nodeIndex,kk);
            nodeCount=nodeCount+1;
        }
        else if (watchType==9) //Total Sxy
        {
            totalResults=totalResults+Sxy(nodeIndex,kk)+Pore(nodeIndex,kk);
            nodeCount=nodeCount+1;
        }
    }
    totalResults=totalResults/double(nodeCount);
    return totalResults;
}

void AxisSymmetric_2D::exportWatchList()
{
    watchListResult.resize(watch.size());
    for (int i=0;i<watch.size();i++)
    {
        //Get information of watch lists
        QString title=watch[i].title;
        int watchType=watch[i].watchType;
        double x0=watch[i].x0;
        double x1=watch[i].x1;
        double y0=watch[i].y0;
        double y1=watch[i].y1;
        bool averageBool=watch[i].averageBool;

        int firstStep=watch[i].beginStep;
        int endStep=watch[i].endStep;
        if(endStep>ns){endStep=ns;}
        if(firstStep<1){firstStep=1;}

        MatrixXd foundNode;
        vector<int> nodeList; //for coordinates output
        nodeList.resize(0);
        bool exportNodeList=false; //node list only needs to be export once only

        if(watchType<=2)
        {
            matrixLib.findXY(coordinates,x0,x1,y0,y1,foundNode);
        }
        else
        {
            matrixLib.findElementXY(coordinates,elements,x0,x1,y0,y1,foundNode);
        }

        MatrixXd outputMatrix;

        for (int kk=firstStep-1;kk<endStep;kk++)
        {
            //Loop over each step
            //If the average bool is false, format result like this
            //Step Time Result1 Result2 Result3...
            //If average bool is true, format table
            //Step Time Result_Averag

            double totalResults=0;
            int nodeCount=0;

            vector<double> nodeListResult;
            nodeListResult.resize(0);

            for (int jj=0;jj<foundNode.rows();jj++)
            {
                int index=foundNode(jj,0)-1;
                if(watchType==0) //X direction
                {
                    double nodalResult=U(index,kk);
                    if(exportNodeList==false){nodeList.push_back(index);}
                    nodeListResult.push_back(nodalResult);

                    totalResults=totalResults+nodalResult;
                    nodeCount=nodeCount+1;
                }
                else if(watchType==1)
                {
                    double nodalResult=V(index,kk);
                    nodeListResult.push_back(nodalResult);
                    if(exportNodeList==false){nodeList.push_back(index);}

                    totalResults=totalResults+nodalResult;
                    nodeCount=nodeCount+1;
                }
                else if(watchType==2)
                {
                    if(nodfmt(index,1)==3)
                    {
                        double nodalResult=P(index,kk);
                        nodeListResult.push_back(nodalResult);
                        if(exportNodeList==false){nodeList.push_back(index);}

                        totalResults=totalResults+nodalResult;
                        nodeCount=nodeCount+1;
                    }
                }
                else if (watchType==3)
                {
                    double nodalResult=Sxx(index,kk);
                    nodeListResult.push_back(nodalResult);
                    if(exportNodeList==false){nodeList.push_back(index);}

                    totalResults=totalResults+nodalResult;
                    nodeCount=nodeCount+1;
                }
                else if (watchType==4)
                {
                    double nodalResult=Syy(index,kk);
                    nodeListResult.push_back(nodalResult);
                    if(exportNodeList==false){nodeList.push_back(index);}

                    totalResults=totalResults+nodalResult;
                    nodeCount=nodeCount+1;
                }
                else if (watchType==5)
                {
                    double nodalResult=Sxy(index,kk);
                    nodeListResult.push_back(nodalResult);
                    if(exportNodeList==false){nodeList.push_back(index);}

                    totalResults=totalResults+nodalResult;
                    nodeCount=nodeCount+1;
                }
                else if (watchType==6)
                {
                    double nodalResult=Pore(index,kk);
                    nodeListResult.push_back(nodalResult);
                    if(exportNodeList==false){nodeList.push_back(index);}

                    totalResults=totalResults+nodalResult;
                    nodeCount=nodeCount+1;
                }
                else if (watchType==7)
                {
                    double nodalResult=Sxx(index,kk)+Pore(index,kk);
                    nodeListResult.push_back(nodalResult);
                    if(exportNodeList==false){nodeList.push_back(index);}

                    totalResults=totalResults+nodalResult;
                    nodeCount=nodeCount+1;
                }
                else if (watchType==8)
                {
                    double nodalResult=Syy(index,kk)+Pore(index,kk);
                    nodeListResult.push_back(nodalResult);
                    if(exportNodeList==false){nodeList.push_back(index);}

                    totalResults=totalResults+nodalResult;
                    nodeCount=nodeCount+1;
                }
                else if (watchType==9)
                {
                    double nodalResult=Sxy(index,kk);
                    nodeListResult.push_back(nodalResult);
                    if(exportNodeList==false){nodeList.push_back(index);}

                    totalResults=totalResults+nodalResult;
                    nodeCount=nodeCount+1;
                }
            }
            totalResults=totalResults/double(nodeCount);

            //for first step
            if(exportNodeList==false)
            {
                //resize output Matrix results
                if(averageBool==false) //
                {
                    outputMatrix.resize(endStep-firstStep+1,nodeListResult.size()+2);
                }
                else
                {
                    outputMatrix.resize(endStep-firstStep+1,3);
                }
                MatrixXd nodeFoundCoordinates;
                nodeFoundCoordinates.resize(nodeList.size(),3);
                for (int nn=0;nn<nodeList.size();nn++)
                {
                    int index=nodeList[nn];
                    nodeFoundCoordinates(nn,0)=index+1;
                    nodeFoundCoordinates(nn,1)=coordinates(index,1);
                    nodeFoundCoordinates(nn,2)=coordinates(index,2);
                }

                if(folderName=="")
                {
                    folderName=QFileDialog::getExistingDirectory(Q_NULLPTR,"Chose save folder");
                }
                QDir mDir;
                QString watchListFolder=folderName+"/WatchList";
                mDir.mkdir(watchListFolder);
                fileName=watchListFolder+"/"+"nodeCoordinate_"+title+".txt";
                exportFile.fileName=fileName.toStdString();
                QStringList header={"Node#","X-coordinate","Y-coordinate"};
                exportFile.ToFile(header,nodeFoundCoordinates);
            }
            exportNodeList=true;

            //Save to ouput Matrix
            if(averageBool==true)
            {
                outputMatrix(kk,0)=firstStep+kk;
                outputMatrix(kk,1)=calculationTime(kk,0);
                outputMatrix(kk,2)=totalResults;
            }
            else
            {
                outputMatrix(kk,0)=firstStep+kk;
                outputMatrix(kk,1)=calculationTime(kk,0);
                for (int nn=0;nn<nodeListResult.size();nn++)
                {
                    outputMatrix(kk,2+nn)=nodeListResult[nn];
                }
            }
        } //End loop over each step

        //Write results
        if(folderName=="")
        {
            folderName=QFileDialog::getExistingDirectory(Q_NULLPTR,"Chose save folder");
        }

        QDir mDir;
        QString watchListFolder=folderName+"/WatchList";
        mDir.mkdir(watchListFolder);

        fileName=watchListFolder+"/"+title+".txt";
        exportFile.fileName=fileName.toStdString();
        if(averageBool==true)
        {
            QStringList header={"Step#","Time(days)","Average_result"};
            exportFile.ToFile(header,outputMatrix);
        }
        else
        {
            QStringList header={"Step#","Time(days)"};
            for(int nn=0;nn<nodeList.size();nn++)
            {
                int index=nodeList[nn]+1;
                QString newString="Node_"+QString::number(index,'f',0);
                header.append(newString);
            }
            exportFile.ToFile(header,outputMatrix);
        }
        watchListResult[i]=outputMatrix;
    }
}

void AxisSymmetric_2D::setCdFactor(double Cd)
{
    this->Cd=Cd;
}

double AxisSymmetric_2D::getError()
{
    return error;
}

void AxisSymmetric_2D::setFolder(QString folderName)
{
    this->folderName=folderName;
    qDebug()<<"Setting folder to: "<<folderName<<endl;

}

//-----------------------------------------------------------------------------------------------------------
//Get Input Data
void AxisSymmetric_2D::getMesh(Ref<MatrixXd> coordinates, Ref<MatrixXd> elements,QString folderName)
{
    this->coordinates=coordinates;
    this->elements=elements;
    this->folderName=folderName;
}

void AxisSymmetric_2D::getStageBoundary(vector<StageBoundaryBase> stageBoundary)
{
    this->stageBoundary=stageBoundary;
}

void AxisSymmetric_2D::getMaterial(vector<MaterialBase> material)
{
    this->material=material;
}

void AxisSymmetric_2D::getStage(vector<StageBase> stage)
{
    this->stage=stage;
}

void AxisSymmetric_2D::getBoundary(vector<BoundaryConditionBase> boundary)
{
    this->boundary=boundary;
}

void AxisSymmetric_2D::getProjectSetting(vector<double> projectParameters)
{
    this->projectParameters=projectParameters;
    if(projectParameters.size()==0)
    {
        this->projectParameters.resize(4);
        this->projectParameters[0]=0;
        this->projectParameters[1]=0;
        this->projectParameters[2]=0;
        this->projectParameters[3]=9.81;
    }
}

void AxisSymmetric_2D::getWatchList(vector<WatchListBase> watch)
{
    this->watch=watch;
}

void AxisSymmetric_2D::getAndExportWatchList(vector<WatchListBase> watch)
{
    this->watch=watch;
    exportWatchList();
}

void AxisSymmetric_2D::getPVDsParameters(double requi, double rw, double rs, double ratioKs, double Cd0, double p0, bool NoSmear,bool goldenSearchFlag)
{
    //------------------------------
    this->requi=requi;
    this->rw=rw;
    this->rs=rs;
    this->ratioKs=ratioKs;
    this->press0=p0;
    this->NoSmear=NoSmear;
    this->Cd0=Cd0;
    this->goldenSearchFlag=goldenSearchFlag;
    //------------------------------
    PVDFlag=true;
    prepareData();
    setPVDparameters();
    runPVDs();
    calculateAnalyticalPVD();
    error= calculatePVDError();
    emit sendPVDResults(Cd0,error);
    if(goldenSearchFlag==false)
    {
        emit sendSolverInformation(U,V,P);
        exportResults();
        exportStress();
        exportWatchList();
        PVDFlag=false;
    }
    else
    {
        //------------------------------
        //Start golden section method
        goldenSearch(0.1,2);
        emit sendSolverInformation(U,V,P);
        exportResults();
        exportStress();
        exportWatchList();
        PVDFlag=false;
    }
}

void AxisSymmetric_2D::getAdditionalInformation(int analysisType, double A1D, double k1D)
{
    this->analysisType=analysisType;
    this->k1D=k1D;
    this->A1D=A1D;
}

void AxisSymmetric_2D::getCrsData(Ref<MatrixXd> crsData, bool crsFlag, int crsType)
{
    this->crsData=crsData;
    this->crsFlag=crsFlag;
    this->crsType=crsType;
}

void AxisSymmetric_2D::setPVDparameters()
{
    Panalytical.resize(non,ns);
    Perror.resize(non,ns);
    calculationTime.resize(ns,1);

    //Calculate Cvr
    double K=material[0].KCurve(0,0); //constant K
    double kv=material[0].kCurve(0,0); //constant kv
    double ratio=material[0].ratio;
    double v=material[0].poission;
    double kh=kv*ratio;
    double ks=ratioKs*kh;
    double G=3*K*(1-2*v)/2/(1+v);
    double mv=1.0f/(K+4.0f*G/3.0f);
    Cvr=kh/gf/mv;

    if(NoSmear==true)
    {
        analyTest.setParameterNoSmear(requi,rw,Cvr,press0);
    }
    else
    {
        analyTest.setParameterSmear(requi,rw,rs,Cvr,kh,ks,press0);
    }
}

void AxisSymmetric_2D::calculateAnalyticalPVD()
{
    for (int kk=0;kk<ns;kk++)
    {
        double realTime=calculationTime(kk,0);
        for (int nn=0;nn<non;nn++)
        {
            double r=coordinates(nn,1); //get r-coordinates
            if(nodfmt(nn,1)==3)
            {
                if(NoSmear==true)
                {
                    Panalytical(nn,kk)=analyTest.analyticalNoSmear(r,realTime);
                }
                else
                {
                    Panalytical(nn,kk)=analyTest.analyticalSmear(r,realTime);
                }
            }
        }
    }
}

double AxisSymmetric_2D::calculatePVDError()
{
    PerrorEachStep.resize(P.cols(),1);
    int count=0;
    double errorTotal=0;
    double errorEle=0;

    for(int i=0;i<P.cols();i++)//Loop over step
    {
        double totalErrorAverageStep=0;
        int countStep=0;

        for (int j=0;j<P.rows();j++) //Loop over node
        {

            if(nodfmt(j,1)==3)
            {
                double poreModel=P(j,i);
                double poreAnalytical=Panalytical(j,i);
                if(poreAnalytical==0)
                {
                    Perror(j,i)=0;
                }
                else
                {
                    errorEle=abs(100*(poreAnalytical-poreModel)/poreAnalytical);
                    Perror(j,i)=errorEle;
                    count=count+1;
                    errorTotal=errorTotal+errorEle;

                    totalErrorAverageStep=totalErrorAverageStep+errorEle;
                    countStep++;
                }
            }
        }
        PerrorEachStep(i,0)=totalErrorAverageStep/countStep;
    }
    error=errorTotal/count;
    cout<<"errorPVD is: "<<error<<endl;
    return error;
}

double AxisSymmetric_2D::PVDError(double Cd)
{
    Cd0=Cd;
    runPVDs();
    error=calculatePVDError();
    return error;
}

double AxisSymmetric_2D::goldenSearch(double a, double b)
{
    double gr=(1+sqrt(5.0))/2.0;
    double tol=1e-3;
    double c=b-(b-a)/gr;
    double d=a+(b-a)/gr;
    double fc=PVDError(c);
    double fd=PVDError(d);
    while(abs(c-d)>tol)
    {
        if(fc<fd)
        {
            b=d;
            d=c;
            fd=fc;
            c=b-(b-a)/gr;
            fc=PVDError(c);
        }
        else
        {
            a=c;
            c=d;
            fc=fd;
            d=a+(b-a)/gr;
            fd=PVDError(d);
        }
        emit sendPVDResults(Cd0,abs(c-d));
    }
    double Cdmin=(a+b)/2.0f;
    emit sendPVDResults(Cdmin,error);
    return Cdmin;
}

void AxisSymmetric_2D::resetSolveData()
{
    //Initial Load vector stress vector
    F=MatrixXd::Zero(totalDof,1);
    F0=MatrixXd::Zero(totalDof,1);
    dF=MatrixXd::Zero(totalDof,1);

    XX=MatrixXd::Zero(totalDof,1);
    X=MatrixXd::Zero(totalDof,ns);
    X0=MatrixXd::Zero(totalDof,1);
    U=MatrixXd::Zero(non,ns);
    V=MatrixXd::Zero(non,ns);
    P=MatrixXd::Zero(non,ns);
    Perror=MatrixXd::Zero(non,ns);

    Sxx=MatrixXd::Zero(noe,ns);
    Syy=MatrixXd::Zero(noe,ns);
    Sxy=MatrixXd::Zero(noe,ns);
    Pore=MatrixXd::Zero(noe,ns);
    SyyInterpolation=MatrixXd::Zero(noe,1);
    SyyGravity=MatrixXd::Zero(noe,1);

    HydraulicVertial=MatrixXd::Zero(noe,ns);
    HydraulicHorizontal=MatrixXd::Zero(noe,ns);
    BulkModulus=MatrixXd::Zero(noe,ns);

    Dirichlet=MatrixXd::Zero(0,2);
    dDirichlet=MatrixXd::Zero(0,2);
    Dirichlet0=MatrixXd::Zero(0,2);
    DirichletAll=MatrixXd::Zero(totalDof,3);
}

void AxisSymmetric_2D::runCrsTest()
{
    qDebug()<<"Start Running Anlysis"<<endl;
    calculationTime.resize(ns,1);
    timer.start();
    currentStep=0;
    realCalculationTime=0;
    resultCrs=MatrixXd::Zero(crsData.rows()-1,7);

    for (int i=0;i<nos;i++) //loop all steps
    {
        int subStep=stage[i].subStep;
        int stageType=stage[i].stageType;
        int timeStepType=stage[i].timeStepType;
        int gravityCheck=stage[i].gravityLoad;

        double t0=stage[i].t0;
        double t1=stage[i].t1;

        for(int j=0;j<subStep;j++)
        {
            //Calculate real time step
            if(timeStepType==0) //constant time step
            {
                dt=(t1-t0)/subStep;
                realCalculationTime=t0+dt*(j+1);
                dt=dt*86400; //convert from day to seconds
            }
            else
            {
                t1=stage[i].timeStep(j+1,0);
                t0=stage[i].timeStep(j,0);
                dt=(t1-t0);
                dt=dt*86400;
                realCalculationTime=t1;
            }

            if(currentStep>0 && crsType==1)
            {
                realStress=crsData(currentStep,2);
                realPore=crsData(currentStep,3);
            }

            qDebug()<<"Time step dt (s) "<<dt<<endl;

            //Create boundary condition from v_fix, v_fixy, v_fixh
            Dirichlet0=Dirichlet;
            createBoundaryConditionEachSubStep(i,realCalculationTime);
            createDichletBoundaryCondition();
            createDirichletBoundaryConditionIncrementalForm();

            //First Step, everything is Zero
            X0.setZero();
            if(currentStep==0)
            {
                X0.setZero();
                SyyInterpolation.setZero();

            }
            else if(currentStep==1 && stage[0].stageType==0) //Ignore displacement because of gravity
            {
                X0.setZero();
            }
            else if(currentStep==1 && stage[0].stageType!=0)
            {
                getPreviousStep();
            }
            else if(currentStep>1 && stage[0].stageType!=0)
            {
                getPreviousStep();
            }
            else if(currentStep>1 && stage[0].stageType==0)
            {
                getPreviousStep();
            }
            else
            {
                qDebug()<<"No case"<<endl;
                pauseSystem();
            }

            //Create load vector
            F0=F;
            dF.setZero();
            createLoadEachStep(gravityCheck);
            dF=F-F0;

            //Assembly matrix
            SyyInterpolation.setZero();
            assemblyGlobalMatrix(); //result is KK matrix
            solveDirect(); //solve KK*XX=dF

            if(currentStep==0) //Insitu, undrained
            {
                X.col(currentStep)=XX;
            }
            else //add to incremental
            {
                X.col(currentStep)=XX+X.col(currentStep-1);
                for (int jj=0;jj<Dirichlet.rows();jj++)
                {
                    int dofIndex=Dirichlet(jj,0);
                    double bcValue=Dirichlet(jj,1);
                    X(dofIndex,currentStep)=bcValue;
                }
            }

            //assign to U, V, P arrays
            for (int jj=0;jj<non;jj++)
            {
                U(jj,currentStep)=X(nodfmt(jj,0)+0,currentStep);
                V(jj,currentStep)=X(nodfmt(jj,0)+1,currentStep);
                if(nodfmt(jj,1)==3)
                {
                    P(jj,currentStep)=X(nodfmt(jj,0)+2,currentStep);
                }
            }

            //Calculate stress
            X0=XX;
            calculateStress(currentStep);
            cout<<"Calculation step: "<<currentStep<<" Maximum vertical displacement: "<<V.col(currentStep).minCoeff()<<endl;
            cout<<"Running time    : "<<timer.elapsed()/1000<< " seconds"<<endl;



            if(currentStep>0) //ignored first step
            {
                resultCrs(currentStep-1,4)=calculateWatchList(1,currentStep+1); //base 1 index
                resultCrs(currentStep-1,5)=calculateWatchList(3,currentStep+1);
                resultCrs(currentStep-1,6)=calculateWatchList(2,currentStep+1);
                resultCrs(currentStep-1,0)=crsData(currentStep,0);
                resultCrs(currentStep-1,1)=crsData(currentStep,1);
                resultCrs(currentStep-1,2)=crsData(currentStep,2);
                resultCrs(currentStep-1,3)=crsData(currentStep,3);
            }

            calculationTime(currentStep,0)=realCalculationTime;
            currentStep=currentStep+1;
        }
    }
    emit sendSolverInformation(U,V,P);

    //export results
    if(folderName=="")
    {
        folderName=QFileDialog::getExistingDirectory(Q_NULLPTR,"Chose save folder");
    }

    fileName=folderName+"/"+"crsResult.txt";
    QStringList header={"Time","Deform","Applied_Stress","Pore_Pressure","Model_Deform","Model_Stress","Model_Pore"};
    exportFile.ToFile(fileName,header,resultCrs);

}

void AxisSymmetric_2D::runBackAnalysis()
{
    qDebug()<<"Start Back-Anlysis for strain control"<<endl;
    calculationTime.resize(ns,1);
    timer.start();
    currentStep=0;
    realCalculationTime=0;
    backAnalysisResult=MatrixXd::Zero(crsData.rows()-1,5);

    for (int i=0;i<nos;i++) //loop all steps
    {
        int subStep=stage[i].subStep;
        int stageType=stage[i].stageType;
        int timeStepType=stage[i].timeStepType;
        int gravityCheck=stage[i].gravityLoad;
        double t0=stage[i].t0;
        double t1=stage[i].t1;

        for(int j=0;j<subStep;j++)
        {
            //Calculate real time step
            if(timeStepType==0) //constant time step
            {
                dt=(t1-t0)/subStep;
                realCalculationTime=t0+dt*(j+1);
                dt=dt*86400; //convert from day to seconds
            }
            else
            {
                t1=stage[i].timeStep(j+1,0);
                t0=stage[i].timeStep(j,0);
                dt=(t1-t0);
                dt=dt*86400;
                realCalculationTime=t1;
            }

            //Create boundary condition from v_fix, v_fixy, v_fixh
            Dirichlet0=Dirichlet;
            createBoundaryConditionEachSubStep(i,realCalculationTime);
            createDichletBoundaryCondition();
            createDirichletBoundaryConditionIncrementalForm();

            //First Step, everything is Zero
            X0.setZero();
            if(currentStep==0)
            {
                X0.setZero();
                SyyInterpolation.setZero();

            }
            else if(currentStep==1 && stage[0].stageType==0) //Ignore displacement because of gravity
            {
                X0.setZero();
            }
            else if(currentStep==1 && stage[0].stageType!=0)
            {
                getPreviousStep();
            }
            else if(currentStep>1 && stage[0].stageType!=0)
            {
                getPreviousStep();
            }
            else if(currentStep>1 && stage[0].stageType==0)
            {
                getPreviousStep();
            }
            else
            {
                qDebug()<<"No case"<<endl;
                pauseSystem();
            }
            MatrixXd X0_initial;
            X0_initial=X0;

            //Create load vector
            F0=F;
            dF.setZero();
            createLoadEachStep(gravityCheck);
            dF=F-F0;
            MatrixXd dF0;
            dF0=dF;

            //Assembly matrix
            SyyInterpolation.setZero();
            if(currentStep<1)
            {
                assemblyGlobalMatrix(); //result is KK matrix
            }
            else
            {
                assemblyBackAnalysis(Ke,kve);
            }
            solveDirect(); //solve KK*XX=dF

            if(currentStep==0) //Insitu, undrained
            {
                X.col(currentStep)=XX;
            }
            else //add to incremental
            {
                X.col(currentStep)=XX+X.col(currentStep-1);
                for (int jj=0;jj<Dirichlet.rows();jj++)
                {
                    int dofIndex=Dirichlet(jj,0);
                    double bcValue=Dirichlet(jj,1);
                    X(dofIndex,currentStep)=bcValue;
                }
            }

            //assign to U, V, P arrays
            for (int jj=0;jj<non;jj++)
            {
                U(jj,currentStep)=X(nodfmt(jj,0)+0,currentStep);
                V(jj,currentStep)=X(nodfmt(jj,0)+1,currentStep);
                if(nodfmt(jj,1)==3)
                {
                    P(jj,currentStep)=X(nodfmt(jj,0)+2,currentStep);
                }
            }

            //Calculate stress
            X0=XX;
            calculateStress(currentStep);

            //star back analysis
            if(currentStep>0 && crsType==1) //start back analysis
            {
                qDebug()<<"Start back analsysis, step: "<<currentStep<<endl;
                getBackAnalysisResults();
                error=max(abs(Jp),abs(Js));
                y0Pore=Jp;
                y0Stress=Js;
                double f=1;
                double k0, k1, K0, K1;
                k0=kve;
                K0=Ke;
                tol=1e-4;
                int loop=0;

                //Calculate error for first try
                if(error>tol)
                {
                    if(abs(Jp)>abs(Js))
                    {
                        k0=kve;
                        K0=Ke;
                        kve=(1+f*Jp)*kve;
                        k1=kve;
                        K1=Ke;
                    }
                    else
                    {
                        k0=kve;
                        K0=Ke;
                        kve=(1+f*Js)*kve;
                        Ke=(1-f*Js)*Ke;
                        k1=kve;
                        K1=Ke;
                    }
                }
                //start back analysis profcess
                while(error>tol)
                {
                    bool loopPore=false;
                    bool loopStress=false;
                    while(abs(Jp)>tol)
                    {
                        loop++;
                        qDebug()<<"Loop to reduce error of pore pressure of step: "<<currentStep<<endl;
                        qDebug()<<"Loop times: "<<loop<<endl;

                        k1=kve;
                        K1=Ke;
                        y0Pore=Jp;
                        y0Stress=Js;
                        X0=X0_initial;
                        dF=dF0;
                        assemblyBackAnalysis(Ke,kve);
                        solveDirect();
                        toVectorSolutions();
                        getBackAnalysisResults();
                        y1Pore=Jp;
                        y1Stress=Js;
                        P0_in<<k0,K1,y0Pore;
                        P1_in<<k1,K1,y1Pore;
                        PlaneLineIntersect();

                        k0=kve;
                        K0=Ke;
                        kve=I_in(0);
                        Ke=I_in(1);
                        if(kve<0)
                        {
                            kve=1e-10;
                            Ke=500;
                        }
                        loopPore=true;
                        loopStress=false;
                    }//end loop pore

                    if(abs(Js)>tol && loopPore==true)
                    {
                        k0=kve;
                        K0=Ke;
                        kve=(1+f*Js)*kve;
                        Ke=(1-f*Js)*Ke;
                    }

                    while(abs(Js)>tol)
                    {
                        loop++;
                        qDebug()<<"Loop to reduce error of applied stress, step: "<<currentStep<<endl;
                        qDebug()<<"Loop times: "<<loop<<endl;

                        k1=kve;
                        K1=Ke;

                        y0Pore=Jp;
                        y0Stress=Js;
                        X0=X0_initial;
                        dF=dF0;
                        assemblyBackAnalysis(Ke,kve);
                        solveDirect();
                        toVectorSolutions();
                        getBackAnalysisResults();
                        y1Pore=Jp;
                        y1Stress=Js;

                        P0_in<<k0,K0,y0Stress;
                        P1_in<<k1,K1,y1Stress;
                        PlaneLineIntersect();

                        k0=kve;
                        K0=Ke;

                        kve=I_in(0);
                        Ke=I_in(1);

                        if(kve<0)
                        {
                            kve=1e-10;
                            Ke=500;
                        }
                        loopStress=true;
                        loopPore=false;
                    }//end loop stress

                    if(abs(Jp)>tol && loopStress==true)
                    {
                        k0=kve;
                        K0=Ke;
                        kve=(1+f*Jp)*kve;
                    }
                    error=max(abs(Jp),abs(Js));
                }//end loop pore and stress
            }//end back analysis process

            if(currentStep>0)
            {
                backAnalysisResult(currentStep-1,0)=averageEffectiveStress(currentStep);
                backAnalysisResult(currentStep-1,1)=crsData(currentStep,4);
                backAnalysisResult(currentStep-1,2)=crsData(currentStep,5);
                backAnalysisResult(currentStep-1,3)=Ke;
                backAnalysisResult(currentStep-1,4)=kve;
            }

            cout<<"Calculation step: "<<currentStep<<" Maximum vertical displacement: "<<V.col(currentStep).minCoeff()<<endl;
            cout<<"Running time    : "<<timer.elapsed()/1000<< " seconds"<<endl;
            calculationTime(currentStep,0)=realCalculationTime;
            currentStep=currentStep+1;
        }
    }
    emit sendSolverInformation(U,V,P);
    //export results
    if(folderName=="")
    {
        folderName=QFileDialog::getExistingDirectory(Q_NULLPTR,"Chose save folder");
    }
    fileName=folderName+"/"+"backAnalysisResult.txt";
    QStringList header={"Eff_Stress","K_ASTM","kv_ASTM","K_BackAnalysis","kv_BackAnalysis"};
    exportFile.ToFile(fileName,header,backAnalysisResult);

    MatrixXd mat=MatrixXd::Zero(backAnalysisResult.rows(),3);
    mat.col(0)=backAnalysisResult.col(0);
    mat.col(1)=backAnalysisResult.col(3);
    mat.col(2)=backAnalysisResult.col(4);
    fileName=folderName+"/"+"mat.dat";
    exportFile.ToFile(fileName,mat);
}

void AxisSymmetric_2D::PlaneLineIntersect()
{
    double D,N,sI;
    int checkInt=0;

    V0_in<<0,0,0; //Point belongs to plane
    I_in<<0,0,0;
    n_in<<0,0,1; //Normal vector;

    u_in=P1_in-P0_in;
    w_in=P0_in-V0_in;
    D=n_in.dot(u_in);
    N=-n_in.dot(w_in);
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
    I_in=P0_in+sI*u_in;
}

void AxisSymmetric_2D::getBackAnalysisResults()
{
    if(currentStep>0 && crsType==1)
    {
        realStress=crsData(currentStep,2);
        realPore=crsData(currentStep,3);
        modelStress=calculateWatchList(3,currentStep+1);
        modelPore=calculateWatchList(2,currentStep+1);
        Jp=(modelPore-realPore)/realPore;
        Js=(modelStress-realStress)/realStress;
    }
}

void AxisSymmetric_2D::toVectorSolutions()
{
    if(currentStep==0) //Insitu, undrained
    {
        X.col(currentStep)=XX;
    }
    else //add to incremental
    {
        X.col(currentStep)=XX+X.col(currentStep-1);
        for (int jj=0;jj<Dirichlet.rows();jj++)
        {
            int dofIndex=Dirichlet(jj,0);
            double bcValue=Dirichlet(jj,1);
            X(dofIndex,currentStep)=bcValue;
        }
    }

    //assign to U, V, P arrays
    for (int jj=0;jj<non;jj++)
    {
        U(jj,currentStep)=X(nodfmt(jj,0)+0,currentStep);
        V(jj,currentStep)=X(nodfmt(jj,0)+1,currentStep);
        if(nodfmt(jj,1)==3)
        {
            P(jj,currentStep)=X(nodfmt(jj,0)+2,currentStep);
        }
    }

    //Calculate stress
    X0=XX;
    calculateStress(currentStep);
}

void AxisSymmetric_2D::assemblyBackAnalysis(double K, double k)
{
    trip_total.clear();
    trip_total.reserve(15*15*noe);
    KK.setZero();
    KK.resize(totalDof,totalDof);

    for (int ele=0;ele<noe;ele++)
    {
        eleNum=ele;
        int nodeCount=elements(eleNum,9);
        double ratio=material[0].ratio;
        Ke=K;
        kve=k;
        ve=material[0].poission;
        Ge=3*Ke*(1-2*ve)/2/(1+ve);
        khe=ratio*kve;
        BulkModulus(eleNum,currentStep)=Ke;
        //Assembly global matrix
        if(nodeCount==6)
        {
            tri6pMatrix(eleNum);
        }
        else if(nodeCount==8)
        {
            quad8pMatrix(eleNum);
        }
        else if(nodeCount==2)
        {
            if(analysisType==1 || analysisType==3)
            {
                line2pMatrix(eleNum);
            }
        }
    }

    KK.setFromTriplets(trip_total.begin(),trip_total.end());
    //Apply boundary condition, global level
    for (int j=0;j<dDirichlet.rows();j++)
    {
        int jj=dDirichlet(j,0);
        dF(jj,0)=dDirichlet(j,1);
        KK.coeffRef(jj,jj)=1;
    }
    KK.prune(0.0);
    KK.makeCompressed();
    trip_total.clear(); //Release memory
}

double AxisSymmetric_2D::averageEffectiveStress(int currentStep)
{
    double totalStress=0;
    for (int i=0;i<noe;i++)
    {
        totalStress=Syy(i,currentStep)+totalStress;
    }
    return totalStress/noe;
}
