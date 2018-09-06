#ifndef LOADMESH_H
#define LOADMESH_H

#include <QDialog>
#include <QFileDialog>
#include <QDebug>
#include <QString>
using namespace std;

namespace Ui {
class LoadMesh;
}

class LoadMesh : public QDialog
{
    Q_OBJECT

public:
    explicit LoadMesh(QWidget *parent = 0);
    ~LoadMesh();
    void getUserData();

private slots:
    void on_okButton_clicked();
    void on_closeButton_clicked();
    void on_coordBrowse_clicked();
    void on_elementBrowse_clicked();

signals:
    void sendSignal(QString,QString);

private:
    Ui::LoadMesh *ui;
    QString coorFileName, elementFileName;
};

#endif // LOADMESH_H
