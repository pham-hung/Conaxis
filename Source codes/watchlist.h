#ifndef WATCHLIST_H
#define WATCHLIST_H

#include <QDialog>
#include "BaseClass/watchlistbase.h"
#include <vector>
#include <QString>

using namespace std;

namespace Ui {
class WatchList;
}

class WatchList : public QDialog,private WatchListBase
{
    Q_OBJECT

public:
    explicit WatchList(QWidget *parent = 0);
    ~WatchList();

    void getUserData();
    void getData(int index);
    void showData();
    void defaultData();
    void sendSIGNAL();
    void createCrsWatchList(double H, double R);

public slots:
    void getWatchListBase(vector<WatchListBase> watchList);

private slots:
    void on_okButton_clicked(); //add,modify watch list
    void on_watchIndexLine_textChanged(const QString &arg1);
    void on_closeButton_clicked();
    void on_exportButton_clicked();
    void on_delButton_clicked();

signals:
    void sendSignal(vector<WatchListBase> watchList);
    void sendSIGNALNow(vector<WatchListBase> watchList);
    private:
    Ui::WatchList *ui;
    vector<WatchListBase> watchList;
};

#endif // WATCHLIST_H
