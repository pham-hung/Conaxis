#include "watchlist.h"
#include "ui_watchlist.h"

WatchList::WatchList(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::WatchList)
{
    ui->setupUi(this);
}

WatchList::~WatchList()
{
    delete ui;
}

void WatchList::getUserData()
{
    watchIndex=ui->watchIndexLine->text().toInt();
    title=ui->titleLine->text();
    if(title=="")
    {
        title="Please enter name";
        ui->titleLine->setText(title);
    }
    watchType=ui->watchCombo->currentIndex();
    x0=ui->x0Line->text().toDouble();
    x1=ui->x1Line->text().toDouble();
    y0=ui->y0Line->text().toDouble();
    y1=ui->y1Line->text().toDouble();
    beginStep=ui->beginStepLine->text().toInt();
    endStep=ui->endStepLine->text().toInt();
    if(endStep<beginStep)
    {
        endStep=beginStep;
        ui->endStepLine->setText(QString::number(endStep));
    }
    if(ui->averageCheck->isChecked())
    {
        averageBool=true;
    }
    else
    {
        averageBool=false;
    }
}

void WatchList::getData(int index)
{
    int i = index-1;
    x0=watchList[i].x0;
    x1=watchList[i].x1;
    y0=watchList[i].y0;
    y1=watchList[i].y1;
    watchType=watchList[i].watchType;
    beginStep=watchList[i].beginStep;
    endStep=watchList[i].endStep;
    title=watchList[i].title;
    averageBool=watchList[i].averageBool;
    showData();
}

void WatchList::showData()
{
    ui->x0Line->setText(QString::number(x0,'f',5));
    ui->x1Line->setText(QString::number(x1,'f',5));
    ui->y0Line->setText(QString::number(y0,'f',5));
    ui->y1Line->setText(QString::number(y1,'f',5));
    ui->beginStepLine->setText(QString::number(beginStep));
    ui->endStepLine->setText(QString::number(endStep));
    ui->titleLine->setText(title);
    ui->watchCombo->setCurrentIndex(watchType);
    if(averageBool==true)
    {
        ui->averageCheck->setChecked(true);
    }
    else
    {
        ui->averageCheck->setChecked(false);
    }
}

void WatchList::defaultData()
{
    ui->titleLine->setText("");
    ui->watchCombo->setCurrentIndex(-1);
    ui->x0Line->setText("");
    ui->x1Line->setText("");
    ui->y0Line->setText("");
    ui->y1Line->setText("");
    ui->beginStepLine->setText("");
    ui->endStepLine->setText("");
    ui->averageCheck->setChecked(true);
}

void WatchList::sendSIGNAL()
{
    emit sendSignal(watchList);
}

void WatchList::createCrsWatchList(double H, double R)
{
    WatchListBase topMovement;
    WatchListBase porePressure;
    WatchListBase totalStress;

    topMovement.watchIndex=1;
    topMovement.title="Deformation";
    topMovement.watchType=1;
    topMovement.x0=0;
    topMovement.x1=R;
    topMovement.y0=H;
    topMovement.y1=H;
    topMovement.beginStep=1;
    topMovement.endStep=99999;
    topMovement.averageBool=true;

    porePressure.watchIndex=2;
    porePressure.title="Pore Pressure";
    porePressure.watchType=2;
    porePressure.x0=0;
    porePressure.x1=R;
    porePressure.y0=0;
    porePressure.y1=0;
    porePressure.beginStep=1;
    porePressure.endStep=99999;
    porePressure.averageBool=true;


    totalStress.watchIndex=3;
    totalStress.title="Total Stress";
    totalStress.watchType=8;
    totalStress.x0=0;
    totalStress.x1=R;
    totalStress.y0=0;
    totalStress.y1=H;
    totalStress.beginStep=1;
    totalStress.endStep=99999;
    totalStress.averageBool=true;

    watchList.resize(0);
    watchList.push_back(topMovement);
    watchList.push_back(porePressure);
    watchList.push_back(totalStress);
}

void WatchList::getWatchListBase(vector<WatchListBase> watchList)
{
    this->watchList=watchList;
}

void WatchList::on_okButton_clicked()
{

    getUserData();
    watchIndex=ui->watchIndexLine->text().toInt();

    //Add new index
    if(watchIndex>(watchList.size()+2))
    {
        watchIndex=watchList.size()+1;
        ui->watchIndexLine->setText(QString::number(watchIndex));
    }

    //Correct index value
    if(watchIndex<1)
    {
        watchIndex=1;
        ui->watchIndexLine->setText(QString::number(watchIndex));
    }

    //Normal index value
    if(watchIndex<=watchList.size())
    {
        int i=watchIndex-1;
        watchList[i].title=title;
        watchList[i].watchType=watchType;
        watchList[i].x0=x0;
        watchList[i].x1=x1;
        watchList[i].y0=y0;
        watchList[i].y1=y1;
        watchList[i].beginStep=beginStep;
        watchList[i].endStep=endStep;
        watchList[i].averageBool=averageBool;
    }
    else
    {
       WatchListBase newBase;
       newBase.title=title;
       newBase.watchIndex=watchIndex;
       newBase.watchType=watchType;
       newBase.x0=x0;
       newBase.x1=x1;
       newBase.y0=y0;
       newBase.y1=y1;
       newBase.beginStep=beginStep;
       newBase.endStep=endStep;
       newBase.averageBool=averageBool;
       watchList.push_back(newBase);
    }

}

void WatchList::on_watchIndexLine_textChanged(const QString &arg1)
{
    watchIndex=ui->watchIndexLine->text().toInt();
    if(watchIndex>(watchList.size()+2))
    {
        watchIndex=watchList.size()+1;
        ui->watchIndexLine->setText(QString::number(watchIndex));
    }

    if(watchIndex<1)
    {
        watchIndex=1;
        ui->watchIndexLine->setText(QString::number(watchIndex));
    }

    if(watchIndex>watchList.size())
    {
        defaultData();
    }
    else
    {
        getData(watchIndex);
        showData();
    }
}

void WatchList::on_closeButton_clicked()
{
    this->close();
}

void WatchList::on_exportButton_clicked()
{
    emit sendSIGNALNow(watchList);
}

void WatchList::on_delButton_clicked()
{
    watchIndex=ui->watchIndexLine->text().toInt();
    vector<WatchListBase> watchList_new;
    watchList_new.resize(0);

    int indexCount=1;
    for(int i=0;i<watchList.size();i++)
    {
        int oldIndex=watchList[i].watchIndex;

        if(watchIndex!=oldIndex)
        {
            watchList_new.push_back(watchList[i]);
            watchList_new[indexCount-1].watchIndex=indexCount;
            indexCount=indexCount+1;
        }
    }

    watchList.resize(0);
    watchList=watchList_new;

    ui->watchIndexLine->setText("1");
}

