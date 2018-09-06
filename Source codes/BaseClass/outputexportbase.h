#ifndef OUTPUTEXPORTBASE_H
#define OUTPUTEXPORTBASE_H


class OutputExportBase
{
public:
    OutputExportBase();
    int solveType=0;
    double tolerance=1e-3;
    int maxIteration=1;
    int nodalSolutionCheck=1;
    int elemetStressCheck=1;
    int averageStressCheck=0;
    int materialParameterCheck=0;
};

#endif // OUTPUTEXPORTBASE_H
