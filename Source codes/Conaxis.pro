#-------------------------------------------------
#
# Project created by QtCreator 2017-05-20T23:16:30
#
#-------------------------------------------------

QT       += core gui opengl
CONFIG +=console

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = Conaxis
TEMPLATE = app
INCLUDEPATH += "C:/Eigen/Eigen"
INCLUDEPATH +=C:/BoostLib
INCLUDEPATH +=C:/BoostLib/stage/lib


INCLUDEPATH += "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries_2016.4.246/windows/mkl/include"
INCLUDEPATH += "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries_2016.4.246/windows/mkl/include/intel64"
INCLUDEPATH += "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries_2016.4.246/windows/mkl/lib/intel64"
INCLUDEPATH += "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries_2016.4.246/windows/mkl/lib/intel64_win"

LIBS += "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries_2016.4.246/windows/mkl/lib/intel64_win/mkl_core.lib"
LIBS += "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries_2016.4.246/windows/mkl/lib/intel64_win/mkl_intel_lp64.lib"
LIBS += "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries_2016.4.246/windows/mkl/lib/intel64_win/mkl_intel_thread.lib"
LIBS += "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries_2016.4.246/windows/compiler/lib/intel64_win/libiomp5md.lib"

win32:RC_ICONS += mainLogo.ico

SOURCES += main.cpp\
        mainwindow.cpp \
    glwidget.cpp \
    input.cpp \
    getfile.cpp \
    colorform.cpp \
    colorform.cpp \
    material.cpp \       
    stage.cpp \
    projectsetting.cpp \
    boundarycondition.cpp \    
    assignboundarycondition.cpp \   
    axissymmetric_2d.cpp \
    gauss2d.cpp \
    pvdlib.cpp \
    writetofile.cpp \
    sortmatrixxd.cpp \
    BaseClass/outputexportbase.cpp \
    solutioncontrol.cpp \
    BaseClass/baseelement.cpp \
    BaseClass/boundaryconditionbase.cpp \
    BaseClass/materialbase.cpp \
    BaseClass/savedatabase.cpp \
    BaseClass/stagebase.cpp \
    BaseClass/stageboundarybase.cpp \
    loadresult.cpp \
    BaseClass/animationbaseclasee.cpp \
    animation.cpp \
    BaseClass/scalesettingbase.cpp \
    scalesetting.cpp \
    loadmesh.cpp \
    watchlist.cpp \
    BaseClass/watchlistbase.cpp \
    crsbackanalysis.cpp \
    pvdbackanalysis.cpp \
    BaseClass/terzaghi1d.cpp \
    geometry.cpp \
    BaseClass/crs_base.cpp

HEADERS  += mainwindow.h \
    glwidget.h \
    input.h \
    getfile.h \
    colorform.h \
    material.h \    
    stage.h \
    projectsetting.h \
    boundarycondition.h \    
    assignboundarycondition.h \    
    axissymmetric_2d.h \
    gauss2d.h \
    pvdlib.h \
    writetofile.h \    
    sortmatrixxd.h \
    BaseClass/outputexportbase.h \
    solutioncontrol.h \
    BaseClass/baseelement.h \
    BaseClass/boundaryconditionbase.h \
    BaseClass/materialbase.h \
    BaseClass/savedatabase.h \
    BaseClass/stagebase.h \
    BaseClass/stageboundarybase.h \
    loadresult.h \
    BaseClass/animationbaseclasee.h \
    animation.h \
    BaseClass/scalesettingbase.h \
    scalesetting.h \
    loadmesh.h \
    watchlist.h \
    BaseClass/watchlistbase.h \
    crsbackanalysis.h \
    pvdbackanalysis.h \
    BaseClass/terzaghi1d.h \
    geometry.h \
    BaseClass/crs_base.h

FORMS    += mainwindow.ui \
    colorform.ui \
    material.ui \
    stage.ui \
    projectsetting.ui \
    boundarycondition.ui \
    assignboundarycondition.ui \
    solutioncontrol.ui \
    loadresult.ui \
    animation.ui \
    scalesetting.ui \
    loadmesh.ui \
    watchlist.ui \
    crsbackanalysis.ui \
    pvdbackanalysis.ui \
    BaseClass/terzaghi1d.ui \
    geometry.ui

RESOURCES += \
    resources.qrc

DISTFILES +=
