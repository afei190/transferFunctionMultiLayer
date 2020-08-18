#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <complex>
#include <QColor>
#include <qwt_plot.h>


QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

struct PLOTOBJECT {
    QwtPlotItem *itemPtr;   // pointer to the PlotItem
    int          index;     // layer index
};

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    PLOTOBJECT itemAt( const QPoint& pos ) const;
    ~MainWindow();

private slots:
    void on_addRowBtn_clicked();
    void on_delRowBtn_clicked();
    void onTableViewUpdated();
    void updateLayerID();
    void setActiveLayer(const QModelIndex &index);
    void setActiveLayer(int index);
    void on_picker_appended(const QPoint &pos);

    void on_thicknessDoubleSpinBox_valueChanged(double arg1);

    void on_shearWaveVelocityDoubleSpinBox_valueChanged(double arg1);

    void on_dampingDoubleSpinBox_valueChanged(double arg1);

    void on_densitydoubleSpinBox_valueChanged(double arg1);

private:
    Ui::MainWindow *ui;

    void updatePlots();
    void plotLayers();
    void calculateNatFreq(int& nPt, double& maxfreq);
    void updateSoilTF();

    int m_activeLayer = 1;
    QwtPlot *plot;
    QVector<double> m_freq;
    QVector<double> m_soilTF;
    QVector<double> m_Nfreq;
    QList<PLOTOBJECT> plotItemList;
};


static QVector<QColor> BRUSH_COLOR({QColor(255, 127, 0, 127),
                                       QColor(191,  95, 0, 127),
                                       QColor(127,  63, 0, 127),
                                       QColor(127, 255, 0, 127),
                                       QColor( 95, 191, 0, 127),
                                       QColor( 63, 127, 0, 127),
                                       QColor(127, 255, 255, 255),
                                       QColor( 95, 191, 255, 255),
                                       QColor( 63, 127, 255, 255),
                                       QColor(127, 255, 255, 127),
                                       QColor( 95, 191, 255, 127),
                                       QColor( 63, 127, 255, 127),
                                       QColor(255, 127, 0, 63),
                                       QColor(191,  95, 0, 63),
                                       QColor(127,  63, 0, 63),
                                       QColor(127, 255, 0, 63),
                                       QColor( 95, 191, 0, 63),
                                       QColor( 63, 127, 0, 63)
                                    });
#endif // MAINWINDOW_H
