#ifndef SCALESETTINGBASE_H
#define SCALESETTINGBASE_H


class ScaleSettingBase
{
public:
    ScaleSettingBase();
    double xScale=1;
    double yScale=1;
    double resultScale=1;
    bool autoMaxMinVal=true;
    double maxVal=100;
    double minVal=-100;
    int step=1;
    bool lockView;
};

#endif // SCALESETTINGBASE_H
