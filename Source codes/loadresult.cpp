#include "loadresult.h"
#include "ui_loadresult.h"

LoadResult::LoadResult(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::LoadResult)
{
    ui->setupUi(this);
}

LoadResult::~LoadResult()
{
    delete ui;
}

void LoadResult::on_ULoadButton_clicked()
{
    UfileName=QFileDialog::getOpenFileName();
    ui->UfileNameLine->setText(UfileName);
}

void LoadResult::on_VLoadButton_clicked()
{
    VfileName=QFileDialog::getOpenFileName();
    ui->VfileNameLine->setText(VfileName);
}

void LoadResult::on_PLoadButton_clicked()
{
    PfileName=QFileDialog::getOpenFileName();
    ui->PfileNameLine->setText(PfileName);
}

void LoadResult::on_okButton_clicked()
{
    ok=true;
    emit sendFileName(UfileName,VfileName,PfileName,ok);
}

void LoadResult::on_cancelButton_clicked()
{
    ok=false;
    emit sendFileName(UfileName,VfileName,PfileName,ok);
    this->close();

}
