QT       += core gui webchannel
unix: QT += svg opengl concurrent
win32: QT += opengl printsupport

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += c++11

include(../SimCenterCommon_fork/Common/Common.pri)
include(../SimCenterCommon_fork/SimFigure/SimFigure.pri)

win32: include(C:\qwt-6.1.5\features\qwt.prf)
unix: include(/usr/local/qwt-6.1.5/features/qwt.prf)

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

SOURCES += \
    UI\BonzaTableModel.cpp \
    UI\BonzaTableView.cpp \
    main.cpp \
    mainwindow.cpp
HEADERS += \
    UI\BonzaTableModel.h \
    UI\BonzaTableView.h \
    GlobalConstances.h \
    mainwindow.h

FORMS += \
    mainwindow.ui
