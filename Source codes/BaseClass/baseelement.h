#ifndef BASEELEMENT_H
#define BASEELEMENT_H
#include <iostream>
#include <QObject>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;

class BaseElement : public QObject
{
    Q_OBJECT
public:
    explicit BaseElement(QObject *parent = nullptr);
    MatrixXi tri6p, tri6d, quad8p, quad8d, line2p;
signals:

public slots:
};

#endif // BASEELEMENT_H
