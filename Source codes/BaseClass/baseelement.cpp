#include "baseelement.h"

BaseElement::BaseElement(QObject *parent) : QObject(parent)
{
    //3 DOF per node
    tri6p=MatrixXi::Zero(6,3);
    tri6d=MatrixXi::Zero(6,3);
    quad8p=MatrixXi::Zero(8,3);
    quad8d=MatrixXi::Zero(8,3);
    line2p=MatrixXi::Zero(2,3);

    tri6p<<1,1,1,
           1,1,1,
           1,1,1,
           1,1,0,
           1,1,0,
           1,1,0;

    tri6d<<1,1,0,
           1,1,0,
           1,1,0,
           1,1,0,
           1,1,0,
           1,1,0;

    quad8p<<1,1,1,
           1,1,1,
           1,1,1,
           1,1,1,
           1,1,0,
           1,1,0,
           1,1,0,
           1,1,0;

    quad8d<<1,1,0,
           1,1,0,
           1,1,0,
           1,1,0,
           1,1,0,
           1,1,0,
           1,1,0,
           1,1,0;

    line2p<<0,0,1,
            0,0,1;
}
