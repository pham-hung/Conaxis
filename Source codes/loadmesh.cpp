#include "loadmesh.h"
#include "ui_loadmesh.h"

LoadMesh::LoadMesh(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::LoadMesh)
{
    ui->setupUi(this);
}

LoadMesh::~LoadMesh()
{
    delete ui;
}

void LoadMesh::getUserData()
{

}

void LoadMesh::on_okButton_clicked()
{
    emit sendSignal(coorFileName,elementFileName);
    this->close();
}

void LoadMesh::on_closeButton_clicked()
{
    this->close();
}

void LoadMesh::on_coordBrowse_clicked()
{
    coorFileName=QFileDialog::getOpenFileName();
    ui->coordLine->setText(coorFileName);
}

void LoadMesh::on_elementBrowse_clicked()
{
    elementFileName=QFileDialog::getOpenFileName();
    ui->elementLine->setText(elementFileName);
}
