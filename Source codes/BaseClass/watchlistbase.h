#ifndef WATCHLISTBASE_H
#define WATCHLISTBASE_H
#include <QString>

class WatchListBase
{
public:
    int watchIndex;
    QString title;
    int watchType;
    double x0, x1, y0, y1;
    int beginStep, endStep;
    WatchListBase();
};

#endif // WATCHLISTBASE_H
