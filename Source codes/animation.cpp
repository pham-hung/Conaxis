#include "animation.h"
#include "ui_animation.h"

Animation::Animation(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::Animation)
{
    ui->setupUi(this);
}

Animation::~Animation()
{
    delete ui;
}

void Animation::getUserData()
{
    beginStep=ui->beginLine->text().toInt();
    endStep=ui->endLine->text().toInt();
    delayTime=ui->delayTime->text().toInt();
    if(ui->radioLockView->isChecked())
    {
        lockView=true;
        unlockView=false;
    }
    else
    {
        lockView=false;
        unlockView=true;
    }

    if(ui->autoChangeColorCheckBox->isChecked())
    {
        autoChangeColor=true;
    }
    else
    {
        autoChangeColor=false;
    }
}

void Animation::on_okButton_clicked()
{
    getUserData();
    animation.beginStep=beginStep;
    animation.endStep=endStep;
    animation.delayTime=delayTime;
    animation.lockView=lockView;
    animation.autoChangeColor=autoChangeColor;
    ok=true;
    emit sendSignal(animation,ok);
}

void Animation::on_cancelButton_clicked()
{
    ok=false;
    this->close();
}
