#include "mainwindow.h"

#include <QApplication>
#include <QFile>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w;
    w.show();

    // set the global stylesheet
    QFile file(":/resources/stylesheet.css");
    if(file.open(QFile::ReadOnly)) {
      QString styleSheet = QLatin1String(file.readAll());
      a.setStyleSheet(styleSheet);
    }

    return a.exec();
}
