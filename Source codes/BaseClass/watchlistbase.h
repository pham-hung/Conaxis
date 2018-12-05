#ifndef WATCHLISTBASE_H
#define WATCHLISTBASE_H
#include <QString>
#include <QDebug>

class WatchListBase
{
public:
    int watchIndex;
    QString title;
    int watchType;
    double x0=0;
    double x1=0;
    double y0=0;
    double y1=0;
    int beginStep=0;
    int endStep=1000;
    bool averageBool=true;
    WatchListBase();
};

#endif // WATCHLISTBASE_H
