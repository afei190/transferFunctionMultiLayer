#define _USE_MATH_DEFINES

#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <cmath>
#include <QDebug>
#include <QVector>
#include <QGridLayout>
#include <QModelIndex>
#include <QMessageBox>
#include <QGuiApplication>
#include <QScreen>

#include "Utils/dialogabout.h"
#include "qwt_picker.h"
#include "qwt_plot_picker.h"
#include "qwt_picker_machine.h"
#include "qwt_plot_shapeitem.h"
#include "qwt_plot_curve.h"


MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    ui->tableView->setSelectionBehavior(QAbstractItemView::SelectRows);
    // set up plot for layers
    plot = new QwtPlot(this);
    QGridLayout *plotLayout = new QGridLayout(this);
    plotLayout->addWidget(plot, 0, 0);
    plotLayout->setMargin(0);
    ui->groupBox->setLayout(plotLayout);

    plot->setCanvasBackground(QBrush(Qt::white));
    plot->setAxisScale(QwtPlot::xBottom, 0, 1, 1);
    plot->setAxisScale(QwtPlot::yLeft, 0, DefaultThickness, 0.5);
    plot->enableAxis(QwtPlot::xBottom, false);

    QwtText text;
    text.setText("Thickness (m)");
    QFont font;
    font.setPointSize(10);
    text.setFont(font);
    plot->setAxisTitle(QwtPlot::yLeft, text);

    //Picker
    QwtPicker *picker = new QwtPicker(plot->canvas());
    picker->setStateMachine(new QwtPickerClickPointMachine);
    picker->setTrackerMode(QwtPicker::AlwaysOff);
    picker->setRubberBand(QwtPicker::RectRubberBand);
    connect(picker, SIGNAL(appended(const QPoint &)), this, SLOT(on_picker_appended(const QPoint &)));

    // add some connections
    connect(ui->tableView->m_sqlModel, SIGNAL(dataChanged(const QModelIndex&,const QModelIndex&)), this, SLOT(onTableViewUpdated()));
    connect(ui->tableView, SIGNAL(cellClicked(const QModelIndex&)), this, SLOT(setActiveLayer(const QModelIndex&)));

    // set range for spin boxes
    ui->shearWaveVelocityDoubleSpinBox->setRange(10, 1000);
    ui->dampingDoubleSpinBox->setRange(0,30);
    ui->thicknessDoubleSpinBox->setRange(0.5,100.0);
    ui->densityDoubleSpinBox->setRange(0.1,5.0);

    // initialize transfer function plot
    ui->TransferFunctionFig->showAxisControls(false);
    ui->TransferFunctionFig->setXLabel("Freq. [Hz]");
    ui->TransferFunctionFig->setYLabel("[H]");
    ui->TransferFunctionFig->setLabelFontSize(8);

    QList<QVariant> valueListLayer;
    valueListLayer << QString ("Layer %1").arg(ui->tableView->m_sqlModel->rowCount() + 1) << DefaultThickness << DefaultDensity << DefaultVs << DefaultEType << "2";
    ui->tableView->insertAt(valueListLayer,0);
    ui->tableView->setTotalHeight(3);
    ui->totalLayerLineEdit->setText("1");
    ui->activeLayerlineEdit->setText(QString::number(m_activeLayer));

    updateSpinBox();
    createActions();
    plotLayers();
    updateSoilTF();
    updatePlots();
}

MainWindow::~MainWindow()
{
    delete ui;
}


void MainWindow::on_addRowBtn_clicked()
{
    QList<QVariant> valueListLayer;
    valueListLayer << QString ("Layer %1").arg(ui->tableView->m_sqlModel->rowCount() + 1) << DefaultThickness << DefaultDensity << DefaultVs << DefaultEType << "2";
    emit ui->tableView->insertAbove(valueListLayer);
    ui->totalLayerLineEdit->setText(QString::number(ui->tableView->m_sqlModel->rowCount()));
    updateLayerID();
    plotLayers();
}

void MainWindow::on_delRowBtn_clicked()
{
    if (ui->tableView->m_sqlModel->rowCount() > 1) {
        emit ui->tableView->removeOneRow(1);
        ui->totalLayerLineEdit->setText(QString::number(ui->tableView->m_sqlModel->rowCount()));
        ui->totalHeight->setText(QString::number(ui->tableView->m_sqlModel->getTotalHeight()));
        updateLayerID();
        plotLayers();
    }
}

void MainWindow::updateLayerID()
{
    for (int i = 1; i <= ui->tableView->m_sqlModel->rowCount();i++) {
        ui->tableView->m_sqlModel->editData(i-1, LAYERNAME, QString ("Layer %1").arg(i));
    }
    ui->tableView->update();
}

void MainWindow::onTableViewUpdated()
{
    ui->totalHeight->setText(QString::number(ui->tableView->totalHeight()));
    updateSpinBox();
    plotLayers();
    updateSoilTF();
    updatePlots();
}

void MainWindow::updateSpinBox()
{
    ui->shearWaveVelocityDoubleSpinBox->setValue(ui->tableView->m_sqlModel->getVS(m_activeLayer-1));
    ui->densityDoubleSpinBox->setValue(ui->tableView->m_sqlModel->getDensity(m_activeLayer-1));
    ui->dampingDoubleSpinBox->setValue(ui->tableView->m_sqlModel->getESize(m_activeLayer-1));
    ui->thicknessDoubleSpinBox->setValue(ui->tableView->m_sqlModel->getThickness(m_activeLayer-1));
}

void MainWindow::updatePlots()
{
    ui->TransferFunctionFig->clear();
    ui->TransferFunctionFig->plot(m_freq, m_soilTF, SimFigure::LineType::Solid, Qt::blue);
    ui->TransferFunctionFig->setLabelFontSize(8);
    ui->TransferFunctionFig->setXLim(0, 40.0);
}

// Following codes are adopted from work by Pedro Arduino and Alborz Ghorfrani of University of Washington

void MainWindow::calculateNatFreq(int& nPt, double& maxfreq)
{
    QVector<double> Nfreq(nPt);
    double dfreq = maxfreq/nPt;

    double TFtan=1.0;
    double TFtan1;
    int ii=0;
    for (int i=1; i<m_freq.size(); i++){
        TFtan1 = (m_soilTF[i] - m_soilTF[i-1])/dfreq;
        if (TFtan1*TFtan <= 0 && TFtan > 0){
            ii = ii+1;
            Nfreq[ii] = m_freq[i];
        }
        TFtan = TFtan1;
    }
    m_Nfreq = Nfreq;
}

void MainWindow::updateSoilTF()
{
    double Psi, Psi1, HH, Vs, Vs1, Rho, Rho1, aux;
    int nPt = 8000;
    double maxfreq = 40.0;
    double dfreq = maxfreq/nPt;
    QVector<double> freq(nPt), omega(nPt);
    QVector<std::complex<double> > mAA(nPt), mBB(nPt);
    QVector<std::complex<double> > mAA1(nPt), mBB1(nPt);
    QVector<std::complex<double> > mAAaux(nPt), mBBaux(nPt);
    QVector<double> Nfreq(nPt);
    QTextStream out(stdout);

    m_freq.resize(nPt);
    m_soilTF.resize(nPt);

    std::complex<double> kstar;
    std::complex<double> Vsstar;
    std::complex<double> kHstar;
    std::complex<double> alphastar;
    std::complex<double> mExpA;
    std::complex<double> mExpB;
    std::complex<double> One;
    std::complex<double> AAA;

    One = std::complex<double>(1,0);
    freq[0]=0.0;
    for (int i=1; i<m_freq.size();i++){
        freq[i]  = freq[i-1]+dfreq;
        omega[i] = 2.0*M_PI*freq[i];
        mAAaux[i]=std::complex<double>(1,0);
        mBBaux[i]=std::complex<double>(1,0);
    }
    m_freq = freq;
    mAA = mAAaux; mBB = mBBaux;

    // set dummy value for rock layer
    int numLayers = ui->tableView->m_sqlModel->rowCount() + 1;
    QVector<QString> vsVec = ui->tableView->m_sqlModel->getvsVec();
    vsVec.append(QString::number(DefaultVs));
    QVector<QString> dampingVec = ui->tableView->m_sqlModel->getesizeVec();  // damping
    dampingVec.append("2");
    QVector<QString> thicknessVec = ui->tableView->m_sqlModel->getthicknessVec();
    thicknessVec.append(QString::number(DefaultThickness));
    QVector<QString> rhoVec = ui->tableView->m_sqlModel->getdensityVec();
    rhoVec.append(QString::number(DefaultDensity));

    if (numLayers == 0)
        return;

    for(int i = 0; i<numLayers-1; i++){
        Psi = dampingVec[i].toDouble() /100.0;
        Psi1 = dampingVec[i+1].toDouble() /100.0;
        HH  = thicknessVec[i].toDouble();
        Vs = vsVec[i].toDouble();
        Vs1 = vsVec[i+1].toDouble();
        Rho = rhoVec[i].toDouble();
        Rho1 = rhoVec[i+1].toDouble();
        mAAaux = mAA; mBBaux = mBB;

        aux = Rho*Vs/(Rho1*Vs1)/(1+Psi1*Psi1);
        alphastar = std::complex<double>(aux*(1+Psi1*Psi1), -aux*(Psi1-Psi));

        for (int ii = 0; ii < freq.size(); ii++ ){
            double km = 2.0* M_PI *freq[ii]/Vs/(1.0+Psi*Psi);
            double kmHH = km*HH;
            kstar = std::complex<double>(km, -km*Psi);
            Vsstar = std::complex<double>(Vs, Vs*Psi);
            kHstar = std::complex<double>(km*HH*Psi, km*HH);
            mExpA = std::exp(kHstar);
            mExpB = std::exp(-kHstar);
            mAA[ii] = 0.5*mAAaux[ii]*(One+alphastar)*mExpA + 0.5*mBBaux[ii]*(One-alphastar)*mExpB;
            mBB[ii] = 0.5*mAAaux[ii]*(One-alphastar)*mExpA + 0.5*mBBaux[ii]*(One+alphastar)*mExpB;
        }
        if (i == 0){
            mAA1 = mAA;
            mBB1 = mBB;
        }
    }

    m_soilTF[0] = 1.0;

    for (int i=1; i<m_freq.size();i++){
        AAA = One/(mAA[i]+mBB[i]);
        m_soilTF[i]= 2*abs(AAA);
    }
    calculateNatFreq(nPt, maxfreq);
}

void MainWindow::on_thicknessDoubleSpinBox_valueChanged(double arg1)
{
    if (ui->tableView->m_sqlModel->rowCount() >= m_activeLayer) {
        ui->tableView->m_sqlModel->editData(m_activeLayer - 1, THICKNESS, arg1);
    }
}

void MainWindow::on_shearWaveVelocityDoubleSpinBox_valueChanged(double arg1)
{
    if (ui->tableView->m_sqlModel->rowCount() >= m_activeLayer) {
        ui->tableView->m_sqlModel->editData(m_activeLayer - 1, VS, arg1);
    }
}

void MainWindow::on_dampingDoubleSpinBox_valueChanged(double arg1)
{
    if (ui->tableView->m_sqlModel->rowCount() >= m_activeLayer) {
        ui->tableView->m_sqlModel->editData(m_activeLayer - 1, ESIZE, arg1); // element size for damping
        ui->tableView->update();
    }
}


void MainWindow::on_densityDoubleSpinBox_valueChanged(double arg1)
{
    if (ui->tableView->m_sqlModel->rowCount() >= m_activeLayer) {
        ui->tableView->m_sqlModel->editData(m_activeLayer - 1, DENSITY, arg1);
        // ui->tableView->update();
    }
}


void MainWindow::setActiveLayer(const QModelIndex &index)
{
    m_activeLayer = index.row() + 1;
    ui->activeLayerlineEdit->setText(QString::number(m_activeLayer));
    updateSpinBox();
}

void MainWindow::setActiveLayer(int index)
{
    m_activeLayer = index;
    // QModelIndex qindex = ui->tableView->model()->index(m_activeLayer - 1, 0);
    emit ui->tableView->onCellSingleClicked(ui->tableView->model()->index(m_activeLayer - 1, 0));
    ui->activeLayerlineEdit->setText(QString::number(m_activeLayer));
    updateSpinBox();
}

// The following part of code is modified from PGT
void MainWindow::plotLayers()
{
    foreach (PLOTOBJECT item, plotItemList) {
        item.itemPtr->detach();
        delete item.itemPtr;
    }
    plotItemList.clear();
    //
    // Plot Ground Layers
    //
    plot->detachItems();
    double currentBase = 0;
    for (int iLayer = ui->tableView->m_sqlModel->rowCount(); iLayer >= 1; iLayer--) {
        QPolygonF groundCorners;
        groundCorners << QPointF(0, currentBase)
                      << QPointF(1, currentBase)
                      << QPointF(1, currentBase + ui->tableView->m_sqlModel->getThickness(iLayer-1))
                      << QPointF(0, currentBase + ui->tableView->m_sqlModel->getThickness(iLayer-1))
                      << QPointF(0, currentBase);
        currentBase += ui->tableView->m_sqlModel->getThickness(iLayer-1);
        QwtPlotShapeItem *layerII = new QwtPlotShapeItem();
        if (iLayer == m_activeLayer) {
            layerII->setPen(QPen(Qt::red, 2));
            layerII->setBrush(QBrush(BRUSH_COLOR[3+iLayer]));
        }
        else {
            layerII->setPen(QPen(BRUSH_COLOR[iLayer], 1));
            layerII->setBrush(QBrush(BRUSH_COLOR[iLayer]));
        }

        layerII->setPolygon(groundCorners);
        layerII->attach(plot);

        PLOTOBJECT var;
        var.itemPtr = layerII;
        var.index   = iLayer;
        plotItemList.append(var);
    }
    plot->setAxisScale(QwtPlot::yLeft, 0, currentBase);
    plot->replot();
}


void MainWindow::on_picker_appended (const QPoint &pos)
{
    PLOTOBJECT    obj = itemAt(pos);
    if (obj.index >= 1) {
        setActiveLayer(obj.index);
        ui->tableView->m_sqlModel->setActive(obj.index-1);
    }
}


PLOTOBJECT MainWindow::itemAt( const QPoint& pos ) const
{
    PLOTOBJECT emptyObj;
    emptyObj.itemPtr = nullptr;
    emptyObj.index   = -1;

    if ( plot == nullptr )
        return emptyObj;

    // translate pos into the plot coordinates
    double coords[ QwtPlot::axisCnt ];
    coords[ QwtPlot::xBottom ] = plot->canvasMap( QwtPlot::xBottom ).invTransform( pos.x() );
    coords[ QwtPlot::xTop ]    = plot->canvasMap( QwtPlot::xTop ).invTransform( pos.x() );
    coords[ QwtPlot::yLeft ]   = plot->canvasMap( QwtPlot::yLeft ).invTransform( pos.y() );
    coords[ QwtPlot::yRight ]  = plot->canvasMap( QwtPlot::yRight ).invTransform( pos.y() );

    for ( int i = plotItemList.size() - 1; i >= 0; i-- )
    {
        PLOTOBJECT obj = plotItemList[i];
        QwtPlotItem *item = obj.itemPtr;
        if ( item->isVisible() && item->rtti() == QwtPlotItem::Rtti_PlotCurve )
        {
            double dist;

            QwtPlotCurve *curveItem = static_cast<QwtPlotCurve *>( item );
            const QPointF p( coords[ item->xAxis() ], coords[ item->yAxis() ] );

            if ( curveItem->boundingRect().contains( p ) )
            {
                // trace curves ...
                dist = 1000.;
                for (size_t line=0; line < curveItem->dataSize() - 1; line++)
                {
                    QPointF pnt;
                    double x, y;

                    pnt = curveItem->sample(line);
                    x = plot->canvasMap( QwtPlot::xBottom ).transform( pnt.x() );
                    y = plot->canvasMap( QwtPlot::yLeft ).transform( pnt.y() );
                    QPointF x0(x,y);

                    pnt = curveItem->sample(line+1);
                    x = plot->canvasMap( QwtPlot::xBottom ).transform( pnt.x() );
                    y = plot->canvasMap( QwtPlot::yLeft ).transform( pnt.y() );
                    QPointF x1(x,y);

                    QPointF r  = pos - x0;
                    QPointF s  = x1 - x0;
                    double s2  = QPointF::dotProduct(s,s);

                    if (s2 > 1e-6)
                    {
                        double xi  = QPointF::dotProduct(r,s) / s2;

                        if ( 0.0 <= xi && xi <= 1.0 )
                        {
                            QPointF t(-s.y()/sqrt(s2), s.x()/sqrt(s2));
                            double d1 = QPointF::dotProduct(r,t);
                            if ( d1 < 0.0 )  { d1 = -d1; }
                            if ( d1 < dist ) { dist = d1;}
                        }
                    }
                    else
                    {
                        dist = sqrt(QPointF::dotProduct(r,r));
                        //QPointF r2 = pos - x1;
                        double d2  = sqrt(QPointF::dotProduct(r,r));
                        if ( d2 < dist ) { dist = d2; }
                    }
                }

                qWarning() << "curve dist =" << dist;

                if ( dist <= 5 ) return obj;
            }
        }
        if ( item->isVisible() && item->rtti() == QwtPlotItem::Rtti_PlotShape )
        {
            QwtPlotShapeItem *shapeItem = static_cast<QwtPlotShapeItem *>( item );
            const QPointF p( coords[ item->xAxis() ], coords[ item->yAxis() ] );

            if ( shapeItem->boundingRect().contains( p ) && shapeItem->shape().contains( p ) )
            {
                return obj;
            }
        }
    }

    return emptyObj;
}


void MainWindow::version()
{
    QString versionText("Version 0.9");
    QMessageBox msgBox;
    QSpacerItem *theSpacer = new QSpacerItem(700, 0, QSizePolicy::Minimum, QSizePolicy::Expanding);
    msgBox.setText(versionText);
    QGridLayout *layout = (QGridLayout*)msgBox.layout();
    layout->addItem(theSpacer, layout->rowCount(),0,1,layout->columnCount());
    msgBox.exec();
}

void MainWindow::about()
{
    DialogAbout *dlg = new DialogAbout();
    QString aboutTitle("A SimCenter Tool for Calculating Transfer Function of Multiple Layers");
    QString aboutSource = ":/resources/textAboutthis.html";
    dlg->setTitle(aboutTitle);
    dlg->setTextSource(aboutSource);

    //
    // adjust size of application window to the available display
    //

    QRect rec = QGuiApplication::primaryScreen()->geometry();
    int height = 0.50*rec.height();
    int width  = 0.50*rec.width();
    dlg->resize(width, height);

    dlg->exec();
    delete dlg;
    /*
    QString aboutText("A SimCenter Tool For Transfer Function Calculation");
    QMessageBox msgBox;
    QSpacerItem *theSpacer = new QSpacerItem(700, 0, QSizePolicy::Minimum, QSizePolicy::Expanding);
    msgBox.setText(aboutText);
    QGridLayout *layout = (QGridLayout*)msgBox.layout();
    layout->addItem(theSpacer, layout->rowCount(),0,1,layout->columnCount());
    msgBox.exec(); */
}

void MainWindow::copyright()
{
    QMessageBox msgBox;
    QString copyrightText = QString("\
                                    <p>\
                                    The source code is licensed under a BSD 2-Clause License:<p>\
                                    \"Copyright (c) 2017-2019, The Regents of the University of California (Regents).\"\
                                    All rights reserved.<p>\
                                    <p>\
                                    Redistribution and use in source and binary forms, with or without \
                                    modification, are permitted provided that the following conditions are met:\
                                    <p>\
                                    1. Redistributions of source code must retain the above copyright notice, this\
                                    list of conditions and the following disclaimer.\
                                    \
                                    \
                                    2. Redistributions in binary form must reproduce the above copyright notice,\
                                    this list of conditions and the following disclaimer in the documentation\
                                    and/or other materials provided with the distribution.\
                                    <p>\
                                    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS \"AS IS\" AND\
                                    ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED\
                                    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE\
                                    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR\
                                    ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES\
                                    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;\
                                    LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND\
            ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT\
            (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS\
            SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.\
            <p>\
            The views and conclusions contained in the software and documentation are those\
            of the authors and should not be interpreted as representing official policies,\
            either expressed or implied, of the FreeBSD Project.\
            <p>\
            REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, \
            THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.\
            THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY, PROVIDED HEREUNDER IS \
            PROVIDED \"AS IS\". REGENTS HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT,\
            UPDATES, ENHANCEMENTS, OR MODIFICATIONS.\
            <p>\
            ------------------------------------------------------------------------------------\
            <p>\
            The compiled binary form of this application is licensed under a GPL Version 3 license.\
            The licenses are as published by the Free Software Foundation and appearing in the LICENSE file\
            included in the packaging of this application. \
            <p>\
            ------------------------------------------------------------------------------------\
            <p>\
            This software makes use of the QT packages (unmodified): core, gui, widgets and network\
                                                                     <p>\
                                                                     QT is copyright \"The Qt Company Ltd&quot; and licensed under the GNU Lesser General \
                                                                     Public License (version 3) which references the GNU General Public License (version 3)\
      <p>\
      The licenses are as published by the Free Software Foundation and appearing in the LICENSE file\
      included in the packaging of this application. \
      <p>\
      ");

    QSpacerItem *theSpacer = new QSpacerItem(700, 0, QSizePolicy::Minimum, QSizePolicy::Expanding);
    msgBox.setText(copyrightText);
    QGridLayout *layout = (QGridLayout*)msgBox.layout();
    layout->addItem(theSpacer, layout->rowCount(),0,1,layout->columnCount());
    msgBox.exec();

}


void MainWindow::createActions() {
    QMenu *helpMenu = menuBar()->addMenu(tr("&Help"));
    QAction *versionAct = helpMenu->addAction(tr("&Version"), this, &MainWindow::version);
    QAction *aboutAct = helpMenu->addAction(tr("&About"), this, &MainWindow::about);
    QAction *copyrightAct = helpMenu->addAction(tr("&License"), this, &MainWindow::copyright);
}
