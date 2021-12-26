#include <QtGui>
#include <QtWidgets/qmessagebox.h>
#include <QtWidgets/qlistwidget.h>
#include <QtWidgets/qstatusbar.h>
#include <QtWidgets/qtoolbar.h>
#include <QtWidgets/qfiledialog.h>
#include <QtWidgets/qaction.h>
#include <QtWidgets/qtablewidget.h>
#include <QtWidgets/qtreewidget.h>
#include <QtWidgets/qmenu.h>
#include <QtWidgets/qmenubar.h>
#include <QtWidgets/qpushbutton.h>
#include <QtWidgets/qlineedit.h>
#include <QtWidgets/qtextedit.h>
#include <QtWidgets/qboxlayout.h>
#include <QtWidgets/qcompleter.h>
#include <qapplication.h>

#undef UNICODE
#include <log4cplus/logger.h>
#include <log4cplus/loggingmacros.h>
#include <log4cplus/fileappender.h>
#include <log4cplus/consoleappender.h>
#include <log4cplus/layout.h>
#include <log4cplus/spi/loggingevent.h>
using namespace log4cplus;
using namespace log4cplus::helpers;
#define UNICODE

#include "mainwindow.h"
#include "ModelGlWidget.h"

#include "MeshData.h"
#include "MeshContainer3d.h"
#include "MeshContainer3dSurface.h"
#include "MeshDomainVolume.h"
#include "MeshDomainSurface.h"
#include "MeshContainer2d.h"
#include "MeshViewSet.h"
#include "MeshGenerator2d.h"

//#define USE_AUTHENTICATION 
   

//#define LOG_MUTEX(txt) LOG4CPLUS_INFO(MeshLog::logger_mesh, txt)
#define LOG_MUTEX(txt) 

MainWindow::MainWindow(int argc, char* argv[]) 
	: tId(0)
{
	//QGLFormat glFormat;
	//glFormat.setVersion( 3, 3 );
	//glFormat.setProfile( QGLFormat::CoreProfile ); // Requires >=Qt-4.8.0
	//glFormat.setSampleBuffers( true );

	modelGLView = new ModelGLWidget();
    setCentralWidget(modelGLView);

	viewDialog = new MeshViewDialog(this);
//	viewDialog->setVisible(false);

    createActions();
    createMenus();
    createToolBars();
    createStatusBar();
    createDockWindows();

    setWindowTitle(tr("QMeshGen"));

	MeshViewSet::m_view = this;

	checkParams(argc, argv);

	readSettings();

//#ifdef USE_AUTHENTICATION
//	mhttp = new MeshHttp(this);
//	validation_result = mhttp->checkStart(user_id);
//#else
//	validation_result = true;
	MeshGenerator2d::param_triangulation_type  = 0;
	MeshGenerator2d::param_quality_improvement = 1;
//#endif

	updateModelInfo();

	propDialog = new MeshPropertiesDialog(this, &(console->meshKernelThread));

//	if(validation_result)
		for(int i = 1; i < argc; i++)
			loadScript(tr(argv[i]));
}

ConsoleQtAppender::ConsoleQtAppender(MainWindow* mw) : m_mw(mw) { }

void ConsoleQtAppender::close()
{
	closed = true;
}

//void ConsoleQtAppender::registerAppender()
//{
//	log4cplus::spi::AppenderFactoryRegistry & reg
//		= log4cplus::spi::getAppenderFactoryRegistry();
//	LOG4CPLUS_REG_PRODUCT(reg, "log4cplus::", ExampleCustomAppender,
//		yournamespace::, log4cplus::spi::AppenderFactory);
//
//}

void ConsoleQtAppender::append(const spi::InternalLoggingEvent& event)
{
	if (!closed) {
		std::stringstream log_out;
		layout->formatAndAppend(log_out, event);
		m_mw->message(log_out.str(), event.getLogLevel());

	}
}

void MainWindow::initLog4WithQt(const string& log_name)
{
	SharedObjectPtr<Appender> append_console(new ConsoleQtAppender(this));
	append_console->setName(LOG4CPLUS_TEXT("QtConsole"));
	auto pattern_console = LOG4CPLUS_TEXT("%m");
	append_console->setLayout(std::make_unique<PatternLayout>(pattern_console));

	MeshLog::initLog4(append_console, log_name);
}

void MainWindow::checkParams(int & argc, char* argv[])
{
	int arg_start = 0;
	if (argc > 2 && strcmp(argv[1], "-logname") == 0) {
		initLog4WithQt(argv[2]);
		arg_start += 2;
	}
	else {
		initLog4WithQt("qtmesh");
	}

	int i = arg_start;
	while(++i < argc){
		if(strcmp(argv[i], "-userid") == 0){
			if(i+1 == argc){
				QMessageBox::critical(this, tr("QMeshGen"),
					tr("The -userid option has to be followed by an actual userid value"));
				--argc;
			}else{
				argc -= 2;
				// set user_id
				user_id = argv[i+1];
				// shift params
				for(int j = i; j < argc; j++)
					argv[j] = argv[j+2];
				--i;
				// store user_id value
				QSettings settings("KI AGH", "QMeshGen");
				settings.setValue("userid", tr(user_id.c_str()));
			}
		}
	}
}

void MainWindow::showEvent(QShowEvent * /* event */)
{
	tId = startTimer(200);
}

void MainWindow::hideEvent(QHideEvent * /* event */)
{
	killTimer(tId);
}

void MainWindow::timerEvent(QTimerEvent *event)
{
	if (event->timerId() == tId) {

		if(viewDialog->closing){
			killTimer(tId);
			return;
		}

		// *** mesh based ***
		bool need_update_mesh = false;
		bool need_update_prop = false;

		QMutex* mutex = console->meshKernelThread.getLockMesh();
		if(mutex->tryLock()){
			LOG_MUTEX("+lm MainWindow::timerEvent");

			if(console->meshKernelThread.getMeshModel().modifiedMesh()){
				need_update_mesh = true;
				console->meshKernelThread.getMeshModel().setModifiedMesh(false);
			}
			if(console->meshKernelThread.getMeshModel().modifiedProperties()){
				need_update_prop = true;
				console->meshKernelThread.getMeshModel().setModifiedProperties(false);
			}
			mutex->unlock();
			LOG_MUTEX("-lm MainWindow::timerEvent");
		}

		if(need_update_mesh){
			// update statistics & view
			updateModelInfo();
		}
		if(need_update_prop){
			prefDock->fillPrefList();
		}
		// *** status/log based ***
		updateLogView();
		console->updateCmdView();
	} else {
		QWidget::timerEvent(event);
	}
}

void MainWindow::addLogItem(const QString& text, log4cplus::LogLevel type, bool update_position)
{
	static const QIcon iconLogInfo(":/images/circ_blue.png");
	static const QIcon iconLogWarning(":/images/circ_yellow.png");
	static const QIcon iconLogError(":/images/circ_red.png");
	static const QIcon iconLogTime(":/images/clock.png");
	static const QIcon iconLogRun(":/images/right_arrow.png");

	switch(type){
	case log4cplus::INFO_LOG_LEVEL:
		logList->addItem(new QListWidgetItem(iconLogInfo, text));
		break;
	case log4cplus::WARN_LOG_LEVEL:
		logList->addItem(new QListWidgetItem(iconLogWarning, text));
		break;
	case log4cplus::ERROR_LOG_LEVEL:
		logList->addItem(new QListWidgetItem(iconLogError, text));
		break;
	case log4cplus::TRACE_LOG_LEVEL:
		logList->addItem(new QListWidgetItem(iconLogTime, text+"s"));
		break;
	case log4cplus::DEBUG_LOG_LEVEL:
		logList->addItem(new QListWidgetItem(iconLogRun, text));
		break;
	}
	if(update_position) 
		logList->setCurrentRow(logList->count()-1);
}

void MainWindow::updateLogView()
{
//	if(!validation_result)
//		QApplication::exit(-1);

	// log
	QList<MainWindow::LogMsg> temp_list;
	{
		LOG_MUTEX("+la MainWindow::updateLogView");
		QMutexLocker lock(&mutexApp);
		temp_list = toshowMsg;
		toshowMsg.clear();
		LOG_MUTEX("-la MainWindow::updateLogView");
	}
	while(!temp_list.empty()){
		LogMsg msg = temp_list.takeFirst();
		addLogItem(msg.text, msg.type, temp_list.empty());
	}

	// status + mesh_view_dlg
	QString status_text;
	ViewMesh view;
	{
		LOG_MUTEX("+la MainWindow::updateLogView");
		QMutexLocker lock(&mutexApp);
		status_text = toshowStatus;
		toshowStatus.clear();
		view = toshowView;
		toshowView.setEmpty();
		LOG_MUTEX("-la MainWindow::updateLogView");
	}

	if(status_text != "")
		statusBar()->showMessage(status_text);

	if(view.set){
		if(MeshViewSet::param_show_visualization != 0){
			QString info;
			if(view.set->m_blocks.count() > 0)
				view.set->addInfo("view blocks", view.set->m_blocks.count() );
			if(view.set->m_faces.count() > 0)
				view.set->addInfo("view faces", view.set->m_faces.count() );
			if(view.set->m_edges.count() > 0)
				view.set->addInfo("view edges", view.set->m_edges.count() );
			if(view.set->m_points.count() > 0)
				view.set->addInfo("view nodes", view.set->m_points.count() );

			viewDialog->setViewSet(view.desc, view.set, view.t, view.with_reset);
		}else{
			delete view.set;
		}
	}

}

void MainWindow::showViewSetNoReset(const string& desc, MeshViewSet* set, int t)
{
	if(viewDialog->closing || viewDialog->skipAll || 
		MeshViewSet::param_show_visualization == 0)
	{
		delete set;
		return;
	}

	{
		LOG_MUTEX("+la MainWindow::showViewSet");
		QMutexLocker lock(&mutexApp);
		if(toshowView.set != nullptr) delete toshowView.set; // only one view-to-show allowed
		toshowView.desc = tr(desc.c_str());
		toshowView.set = set;
		toshowView.t = t;
		toshowView.with_reset = false;
		viewDialog->waitingForUser = true;
		LOG_MUTEX("-la MainWindow::showViewSet");
	}
	if(t < 0) // wait for user-interaction
		viewDialog->viewWaitSemaphore.acquire();
	else if(t > 0) // wait for user-interaction, but not more than t milliseconds
		viewDialog->viewWaitSemaphore.tryAcquire(1, t);
	// else don't wait at all ...

	viewDialog->waitingForUser = false;

	if(viewDialog->closing) return;
	int av = viewDialog->viewWaitSemaphore.available();
	if(av > 0) // shouldn't really happen...
		viewDialog->viewWaitSemaphore.release(av);
}


void MainWindow::showViewSet(const string& desc, MeshViewSet* set, int t)
{
	if(viewDialog->closing || viewDialog->skipAll || 
		MeshViewSet::param_show_visualization == 0)
	{
		delete set;
		return;
	}

	{
		LOG_MUTEX("+la MainWindow::showViewSet");
		QMutexLocker lock(&mutexApp);
		if(toshowView.set != nullptr) delete toshowView.set; // only one view-to-show allowed
		toshowView.desc = tr(desc.c_str());
		toshowView.set = set;
		toshowView.t = t;
		toshowView.with_reset = true;
		viewDialog->waitingForUser = true;
		LOG_MUTEX("-la MainWindow::showViewSet");
	}
	if(t < 0) // wait for user-interaction
		viewDialog->viewWaitSemaphore.acquire();
	else if(t > 0) // wait for user-interaction, but not more than t milliseconds
		viewDialog->viewWaitSemaphore.tryAcquire(1, t);
	// else don't wait at all ...

	viewDialog->waitingForUser = false;

	if(viewDialog->closing) return;
	int av = viewDialog->viewWaitSemaphore.available();
	if(av > 0) // shouldn't really happen...
		viewDialog->viewWaitSemaphore.release(av);
}

void MainWindow::open()
{
	static QString dir = ".";
	QString fullFileName = QFileDialog::getOpenFileName(this,
                        tr("Choose a file name"), dir,
                        tr("XML Mesh file (*.xml);;Script file (*.txt)"));
    if (fullFileName.isEmpty())
  		return;

	dir = fullFileName;

	if(QFileInfo(fullFileName).suffix() == "xml")
		loadFile(fullFileName);
	else
		loadScript(fullFileName);
}

void MainWindow::loadScript(const QString& fname)
{
     QFile file(fname);
	 if (!file.open(QIODevice::ReadOnly | QIODevice::Text)){
		 addLogItem(tr("Error opening ")+fname, log4cplus::ERROR_LOG_LEVEL);
         return;
	 }

     QTextStream in(&file);
     while (!in.atEnd()) {
         QString line = in.readLine();
		 line.remove(QRegExp("#.*"));
		 line.remove(QRegExp("^\\s+"));
		 if(!line.isEmpty())
			 console->addMeshCommand(line);
     }
}

void MainWindow::saveGRD()
{
	static QString dir = ".";
	QString fullFileName = QFileDialog::getSaveFileName(this,
						tr("Choose a file name"), dir,
						tr("GRD file (*.grd)"));
	if(fullFileName.isEmpty()) return;

	console->addMeshCommand(tr("store_grd ") + fullFileName);
}

void MainWindow::saveTXT()
{
	static QString dir = ".";
	QString fullFileName = QFileDialog::getSaveFileName(this,
						tr("Choose a file name"), dir,
						tr("TXT file (*.txt)"));
	if(fullFileName.isEmpty()) return;

	if(fullFileName.endsWith(".txt"))
		fullFileName.truncate(fullFileName.size()-4);
	console->addMeshCommand(tr("store_txt ") + fullFileName);
}

void MainWindow::loadFile(const QString& fname)
{
	console->addMeshCommand(tr("load ") + fname);

	recentFiles.removeAll(fname);
	recentFiles.prepend(fname);
	updateRecentFileActions();
}

void MainWindow::updateRecentFileActions()
{
	QMutableStringListIterator it(recentFiles);
	while (it.hasNext()) {
		if (!QFile::exists(it.next()))
			it.remove();
	}

	for (int j = 0; j < MaxRecentFiles; j++) {
		if (j < recentFiles.count()) {
			QString text = tr("&%1 %2").arg(j + 1).arg(QFileInfo(recentFiles[j]).fileName());
			recentFileActions[j]->setText(text);
			recentFileActions[j]->setData(recentFiles[j]);
			recentFileActions[j]->setVisible(true);
		} else {
			recentFileActions[j]->setVisible(false);
		}
	}
	separatorAction->setVisible(!recentFiles.isEmpty());
}

void MainWindow::openRecentFile()
{
	QAction *action = qobject_cast<QAction *>(sender());
	if (action)
		loadFile(action->data().toString());
}

void MainWindow::writeSettings()
{
	QSettings settings("KI AGH", "QMeshGen");
	settings.setValue("fgeometry", frameGeometry());
	settings.setValue("geometry", geometry());
	settings.setValue("recentFiles", recentFiles);
	settings.setValue("cmdHistory", console->cmdHistory);
//	settings.setValue("prefDockSize", prefDock->size());
	settings.setValue("prefColumnWidth0", prefDock->prefTable->columnWidth(0));
	settings.setValue("prefColumnWidth1", prefDock->prefTable->columnWidth(1));

	settings.setValue("userid", tr(user_id.c_str()));
}

void MainWindow::readSettings()
{
	QSettings settings("KI AGH", "QMeshGen");

	QRect frect = settings.value("fgeometry",
		QRect(200, 200, 800, 600)).toRect();
	move(frect.topLeft());
	QRect rect = settings.value("geometry",
		QRect(200, 200, 800, 600)).toRect();
	resize(rect.size());

	recentFiles = settings.value("recentFiles").toStringList();
	updateRecentFileActions();

	console->setHistory(settings.value("cmdHistory").toStringList());

//	QSize size = settings.value("prefDockSize", QSize(0, 0)).toSize();
//	if(size.width() > 0) prefDock->resize(size);

	int width = settings.value("prefColumnWidth0", 0).toInt();
	if(width > 0) prefDock->prefTable->setColumnWidth(0, width);
	width = settings.value("prefColumnWidth1", 0).toInt();
	if(width > 0) prefDock->prefTable->setColumnWidth(1, width);

	user_id = settings.value("userid", "").toString().toStdString();
}

void MainWindow::triangulateAction()
{
	console->addMeshCommand(tr("clear_cs"));
	console->addMeshCommand(tr("triangulate"));
}

void MainWindow::quadAction()
{
	console->addMeshCommand(tr("make_quads"));
}

void MainWindow::tetrahedraAction()
{
	console->addMeshCommand(tr("clear_cs"));
	console->addMeshCommand(tr("triangulate3"));
}

void MainWindow::parseGrainAction()
{
	static QString dir = ".";
	QString fullFileName = QFileDialog::getOpenFileName(this,
                        tr("Choose a grain file"), dir,
                        tr("Grain lines file (LinesShort.txt)"));
    if (fullFileName.isEmpty())
  		return;
	dir = fullFileName.left(fullFileName.length() - 14);

	console->addMeshCommand(tr("parse_grain ") + dir);
}

void MainWindow::meshingPropAction()
{
	propDialog->show();
}

void MainWindow::about()
{
	QString text;
	text.sprintf("Mesh Generator (%s), Tomasz Jurczyk. "
		"Automated generation of unstructured anisotropic surface and volume meshes.", 
		mesh_data.version().c_str());

	QMessageBox::about(this, tr("About Mesh Generator"), text);
}

void MainWindow::createActions()
{
    openAct = new QAction(QIcon(":/images/open.png"), tr("&Open..."), this);
    openAct->setShortcut(tr("Ctrl+O"));
    openAct->setStatusTip(tr("Open mesh model or script file"));
    connect(openAct, SIGNAL(triggered()), this, SLOT(open()));

	for (int i = 0; i < MaxRecentFiles; ++i) {
		recentFileActions[i] = new QAction(this);
		recentFileActions[i]->setVisible(false);
		connect(recentFileActions[i], SIGNAL(triggered()),
			this, SLOT(openRecentFile()));
	}

	saveGRDAct = new QAction(QIcon(":/images/save_grd.png"), tr("Export to &GRD..."), this);
    saveGRDAct->setStatusTip(tr("Export mesh to GRD file"));
    connect(saveGRDAct, SIGNAL(triggered()), this, SLOT(saveGRD()));

	saveTXTAct = new QAction(QIcon(":/images/save_txt.png"), tr("Export to &TXT..."), this);
    saveTXTAct->setStatusTip(tr("Export mesh to TXT file"));
    connect(saveTXTAct, SIGNAL(triggered()), this, SLOT(saveTXT()));

    genTriangleAct = new QAction(QIcon(":/images/triangle.png"), tr("Generate triangles"), this);
    genTriangleAct->setStatusTip(tr("Automatic generation of triangular meshes for surfaces"));
    connect(genTriangleAct, SIGNAL(triggered()), this, SLOT(triangulateAction()));

    genQuadAct = new QAction(QIcon(":/images/quad.png"), tr("Generate quads"), this);
    genQuadAct->setStatusTip(tr("Automatic conversion of triangular surface meshes into quadrilateral"));
    connect(genQuadAct, SIGNAL(triggered()), this, SLOT(quadAction()));

    genTetraAct = new QAction(QIcon(":/images/tetrahedron.png"), tr("Generate tetrahedra"), this);
    genTetraAct->setStatusTip(tr("Automatic generation of tetrahedral meshes for volume blocks"));
    connect(genTetraAct, SIGNAL(triggered()), this, SLOT(tetrahedraAction()));

    parseGrainAct = new QAction(tr("Parse grain file"), this);
    parseGrainAct->setStatusTip(tr("Parse text files for multi-grain description"));
    connect(parseGrainAct, SIGNAL(triggered()), this, SLOT(parseGrainAction()));

    propDlgAct = new QAction(tr("Meshing properties"), this);
    propDlgAct->setStatusTip(tr("Show dialog for most important meshing properties"));
    connect(propDlgAct, SIGNAL(triggered()), this, SLOT(meshingPropAction()));

    aboutAct = new QAction(QIcon(":/images/spider_web.png"), tr("&About"), this);
    aboutAct->setStatusTip(tr("Show the application's About box"));
    connect(aboutAct, SIGNAL(triggered()), this, SLOT(about()));

    quitAct = new QAction(tr("&Quit"), this);
    quitAct->setShortcut(tr("Ctrl+Q"));
    quitAct->setStatusTip(tr("Quit the application"));
    connect(quitAct, SIGNAL(triggered()), this, SLOT(close()));

}

void MainWindow::createMenus()
{
    QMenu *fileMenu = menuBar()->addMenu(tr("&File"));
    fileMenu->addAction(openAct);
	separatorAction = fileMenu->addSeparator();
	for (int i = 0; i < MaxRecentFiles; i++)
		fileMenu->addAction(recentFileActions[i]);
    fileMenu->addSeparator();
	fileMenu->addAction(saveGRDAct);
	fileMenu->addAction(saveTXTAct);
    fileMenu->addSeparator();
    fileMenu->addAction(quitAct);

    viewMenu = menuBar()->addMenu(tr("&View"));

	QMenu *meshingMenu = menuBar()->addMenu(tr("&Meshing"));
	meshingMenu->addAction(genTriangleAct);
	meshingMenu->addAction(genQuadAct);
	meshingMenu->addSeparator();
	meshingMenu->addAction(genTetraAct);

	QMenu *specialMenu = menuBar()->addMenu(tr("&Special"));
	specialMenu->addAction(parseGrainAct);
	specialMenu->addAction(propDlgAct);

    menuBar()->addSeparator();

    QMenu *helpMenu = menuBar()->addMenu(tr("&Help"));
    helpMenu->addAction(aboutAct);
}

void MainWindow::createToolBars()
{
    QToolBar *fileToolBar = addToolBar(tr("File"));
    fileToolBar->addAction(openAct);
	fileToolBar->addAction(saveGRDAct);
	fileToolBar->addAction(saveTXTAct);

    QToolBar *meshing2dToolBar = addToolBar(tr("Meshing 2D"));
    meshing2dToolBar->addAction(genTriangleAct);
    meshing2dToolBar->addAction(genQuadAct);

    QToolBar *meshing3dToolBar = addToolBar(tr("Meshing 3D"));
    meshing3dToolBar->addAction(genTetraAct);
}

void MainWindow::createStatusBar()
{
    statusBar()->showMessage(tr("Ready"));
}

void MainWindow::createDockWindows()
{
	QDockWidget *dock = new QDockWidget(tr("Log"), this);
    dock->setAllowedAreas(Qt::TopDockWidgetArea | Qt::BottomDockWidgetArea);
	logList = new QListWidget(dock);
	logList->addItem(tr(mesh_data.version().c_str()));
	dock->setWidget(logList);
    addDockWidget(Qt::BottomDockWidgetArea, dock);
    viewMenu->addAction(dock->toggleViewAction());

	dock = new QDockWidget(tr("Command list"), this);
    dock->setAllowedAreas(Qt::TopDockWidgetArea | Qt::BottomDockWidgetArea);
	console = new ConsoleWidget(viewDialog, dock);
	dock->setWidget(console);
    addDockWidget(Qt::BottomDockWidgetArea, dock);
    viewMenu->addAction(dock->toggleViewAction());

    dock = new QDockWidget(tr("Mesh model"), this);
    dock->setAllowedAreas(Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea);
    statTree = new QTreeWidget(dock);

	statTree->setColumnCount(2);
	statTree->setHeaderHidden(true);

    dock->setWidget(statTree);
    addDockWidget(Qt::LeftDockWidgetArea, dock);
    viewMenu->addAction(dock->toggleViewAction());

    itemViewAct = new QAction(tr("&View"), this);
    itemViewAct->setStatusTip(tr("View selected mesh"));
    connect(itemViewAct, SIGNAL(triggered()), this, SLOT(showMeshItemView()));
	statTree->addAction(itemViewAct);

	itemStatAct = new QAction(tr("&Statistics"), this);
    itemStatAct->setStatusTip(tr("Show detailed statistics"));
    connect(itemStatAct, SIGNAL(triggered()), this, SLOT(showMeshItemStat()));
	statTree->addAction(itemStatAct);

//	itemClearAct = new QAction(tr("&Clear"), this);	
//	itemClearAct->setStatusTip(tr("Remove selected mesh"));
//	connect(itemClearAct, SIGNAL(triggered()), this, SLOT(clearMeshItem()));
//	statTree->addAction(itemClearAct);

	connect(statTree, SIGNAL(itemSelectionChanged()), this, SLOT(statTreeSelectionChanged()));
	statTree->setContextMenuPolicy(Qt::ActionsContextMenu);

	prefDock = new PropertiesWidget(&(console->meshKernelThread), this);
    addDockWidget(Qt::LeftDockWidgetArea, prefDock);
    viewMenu->addAction(prefDock->toggleViewAction());
}

MeshPropertiesDialog::MeshPropertiesDialog(QWidget *parent, MeshKernelThread* mkt)
    : QDialog(parent), meshKernelThread(mkt)
{
	prop_labels  
		<< tr("Sharp edge threshold")
		<< tr("Local shape tolerance")
		<< tr("Max relative length 3D")
		<< tr("Max relative length 2D")
		<< tr("Max size stretching")
		<< tr("Max size gradation")
		<< tr("Curvature sizing ratio")
		<< tr("Tolerance for grain surface/curves")
		<< tr("Tolerance for grain nodes");
	prop_names
		<< tr("MG3S_SHARP_EDGE_THRESHOLD")
		<< tr("MG3S_LOCAL_SHAPE_TOLERANCE")
		<< tr("ACS3D_DIAMETER_MAX_RATIO")
		<< tr("ACS2D_DIAMETER_MAX_RATIO")
		<< tr("ACS_STRETCH_MAX_RATIO")
		<< tr("ACS_GRADATION_RATIO")
		<< tr("ACS_CURVATURE_RATIO")
		<< tr("GRAIN_TOLERANCE")
		<< tr("GRAIN_NODE_TOLERANCE");

	int prop_ct = prop_labels.count();

	{
		LOG_MUTEX("+lm MeshPropertiesDialog::MeshPropertiesDialog");
		QMutexLocker lm(meshKernelThread->getLockMesh());

		for(int i = 0; i < prop_ct; i++){
			auto p = mesh_data.getProperty(prop_names[i].toStdString());
			if(!p) prop_values << tr("err");
			else{
				if(p->isInt()){
					prop_values << tr("%1").arg(p->getIntValue());
				}else{
					prop_values << tr("%1").arg(p->getDoubleValue());
				}
			}
		}
		LOG_MUTEX("-lm MeshPropertiesDialog::MeshPropertiesDialog");
	}

	QVBoxLayout *vbox = new QVBoxLayout;

	propTable = new QTableWidget(this);
	propTable->setColumnCount(2);
	propTable->setHorizontalHeaderLabels(
		QStringList() << tr("Property") << tr("Value") );

	propTable->setRowCount(prop_ct);
	for(int i = 0; i < prop_ct; i++){
		QTableWidgetItem *item = new QTableWidgetItem(prop_labels[i]);
		item->setFlags(Qt::ItemIsEnabled);
		propTable->setItem(i, 0, item);

		item = new QTableWidgetItem(prop_values[i]);
		item->setFlags(Qt::ItemIsEditable | Qt::ItemIsEnabled);
		propTable->setItem(i, 1, item);
	}

	vbox->addWidget(propTable);

	for(int i = 0; i < 2; i++)
			propTable->resizeColumnToContents(i);

	QHBoxLayout *clipBtnLayout = new QHBoxLayout;
	QPushButton *applyButton = new QPushButton(tr("Apply"));
	connect(applyButton, SIGNAL(clicked()), this, SLOT(apply()));
	QPushButton *closeButton = new QPushButton(tr("Close"));
	connect(closeButton, SIGNAL(clicked()), this, SLOT(close()));
	clipBtnLayout->addWidget(applyButton);
	clipBtnLayout->addWidget(closeButton);

	vbox->addLayout(clipBtnLayout);

	// apply layout
	setLayout(vbox);

	setWindowTitle(tr("Meshing properties"));
}

void MeshPropertiesDialog::apply()
{
	LOG_MUTEX("+lm MeshPropertiesDialog::apply");
	QMutexLocker lm(meshKernelThread->getLockMesh());

	int prop_ct = prop_labels.count();
	bool ok;
	QString strValue, oldStrValue, newStrValue;

	for(int i = 0; i < prop_ct; i++){
		strValue = propTable->item(i, 1)->text();

		auto p = mesh_data.getProperty(prop_names[i].toStdString());
		if(!p) continue;

		if(p->isInt()){
			int old_value = p->getIntValue();
			oldStrValue = tr("%1").arg(old_value);
			int value = strValue.toInt(&ok);
			if(ok){
				if(value == old_value) continue;
				p->setIntValue(value);
				newStrValue = tr("%1").arg(value);
			}
		}else{
			double old_value = p->getDoubleValue();
			oldStrValue = tr("%1").arg(old_value);
			double value = strValue.toDouble(&ok);
			if(ok){
				if(value == old_value) continue;
				p->setDoubleValue(value);
				newStrValue = tr("%1").arg(value);
			}
		}
		if(ok){
			QString info = tr("Changed %1 to %2").arg(prop_names[i]).arg(newStrValue);
			LOG4CPLUS_INFO(MeshLog::logger_console, info.toStdString());
		}else{
			QString info = tr("Parse error for %1 -> %2").arg(prop_names[i]).arg(strValue);
			LOG4CPLUS_ERROR(MeshLog::logger_console, info.toStdString());
		}
	}
	LOG_MUTEX("-lm MeshPropertiesDialog::apply");
}

void PropertiesWidget::updateMeshProperty(int row, int column)
{
	if(column != 1) return; // only if value was edited

	bool ok = false;
	QString strValue;
	QString newStrValue;
	QString oldStrValue;
	string name;

	{
		LOG_MUTEX("+lm PropertiesWidget::updateMeshProperty");
		QMutexLocker lm(meshKernelThread->getLockMesh());

		QTableWidgetItem *item = prefTable->item(row, 0);
		name = item->text().toUpper().toStdString();
		item = prefTable->item(row, column);
		strValue = item->text();

		auto p = mesh_data.getProperty(name);
//		assert(p);
		if(!p) return;

		if(p->isInt()){
			int old_value = p->getIntValue();
			oldStrValue = tr("%1").arg(old_value);
			int value = strValue.toInt(&ok);
			if(ok){
				if(value == old_value) return;
				p->setIntValue(value);
				newStrValue = tr("%1").arg(value);
			}
		}else{
			double old_value = p->getDoubleValue();
			oldStrValue = tr("%1").arg(old_value);
			double value = strValue.toDouble(&ok);
			if(ok){
				if(value == old_value) return;
				p->setDoubleValue(value);
				newStrValue = tr("%1").arg(value);
			}
		}
		if(ok){
			QString info = tr("Changed %1 to %2").arg(name.c_str()).arg(newStrValue);
			LOG4CPLUS_INFO(MeshLog::logger_console, info.toStdString());
		}else{
			QString info = tr("Parse error for %1 -> %2").arg(name.c_str()).arg(strValue);
			LOG4CPLUS_ERROR(MeshLog::logger_console, info.toStdString());
		}
		LOG_MUTEX("-lm PropertiesWidget::updateMeshProperty");
	}

	if(!ok){
		QTableWidgetItem *item = prefTable->item(row, column);
		prefTable->blockSignals(true);
		item->setText(oldStrValue);
		prefTable->blockSignals(false);
	}
}

void PropertiesWidget::fillPrefList()
{
	QStringList names;
	QStringList values;
	QStringList descriptions;

	int ct = 0;
	{
		LOG_MUTEX("+lm PropertiesWidget::fillPrefList");
		QMutexLocker lm(meshKernelThread->getLockMesh());

		auto pvalues = mesh_data.m_parameters.values();
		ct = pvalues.countInt();
		for(auto p : pvalues){
			names.append(tr(p->name.c_str()));
			descriptions.append(tr(p->description.c_str()));
			if(p->isInt()){
				values.append(tr("%1").arg(p->getIntValue()));
			}else{
				values.append(tr("%1").arg(p->getDoubleValue()));
			}
		}
		LOG_MUTEX("-lm PropertiesWidget::fillPrefList");
	}

	prefTable->clearContents();
	prefTable->setRowCount(ct);
	for(int i = 0; i < ct; i++){
		QTableWidgetItem *item = new QTableWidgetItem(names.at(i));
		item->setFlags(Qt::ItemIsEnabled);
		item->setToolTip(descriptions.at(i));
		prefTable->setItem(i, 0, item);

		item = new QTableWidgetItem(values.at(i));
		item->setFlags(Qt::ItemIsEditable | Qt::ItemIsEnabled);
		prefTable->setItem(i, 1, item);
	}
	prefTable->sortItems(0);
	prefTable->resizeColumnToContents(1);
	prefTable->resizeRowsToContents();
}

void PropertiesWidget::resizeEvent(QResizeEvent * /* event */)
{
//	prefTable->resizeColumnToContents(1);
	int diff = width() - prefTable->width();
	prefTable->setColumnWidth(0, prefTable->columnWidth(0) + diff); 
//	prefTable->resizeRowsToContents();
}

void MainWindow::message(const string& text, log4cplus::LogLevel msg_type)
{
	LOG_MUTEX("+la MainWindow::message");
	QMutexLocker lock(&mutexApp);
	toshowMsg.append(LogMsg(tr(text.c_str()), msg_type));
	LOG_MUTEX("-la MainWindow::message");
}

void MainWindow::status(const string& text)
{
	LOG_MUTEX("+la MainWindow::status");
	QMutexLocker lock(&mutexApp);
	toshowStatus = tr(text.c_str());
	LOG_MUTEX("-la MainWindow::status");
}

void MainWindow::updateModelInfo()
{
	static const QIcon iconView(":/images/blocks.png");

	QList<QTreeWidgetItem*> items;

	{
		LOG_MUTEX("+lm MainWindow::updateModelInfo");
		QMutexLocker lm(console->meshKernelThread.getLockMesh());
		MeshContainer3d* mesh = console->meshKernelThread.getMeshModel().getMesh();
		if(mesh){
			MeshContainer3dSurface* total_surface_mesh = mesh->getTotalMDV()->getSurfaceMesh();
			MeshContainer3d* total_volume_mesh = mesh->getTotalMDV()->getMesh();

			QTreeWidgetItem* model_item = nullptr;
			int mbct = mesh->getBlocksCount();
			if(mbct > 0 || total_surface_mesh || total_volume_mesh){
				model_item = new QTreeWidgetItem((QTreeWidget*)0, 
					QStringList( (mbct > 1) ? QString("model (domains: %1)").arg(mbct) : QString("model") ));
				model_item->setData(0, Qt::UserRole + MVIEW_DATA, 
					qVariantFromValue(MeshVariant(MeshVariant::DOMAIN_MODEL, mesh)));
				items.append(model_item);
			}

			if(total_surface_mesh){
				int smesh_ct = total_surface_mesh->getFacesCount();
				QTreeWidgetItem* mesh_item = 
					new QTreeWidgetItem(model_item, 
							QStringList(tr("surface mesh (NF=%1)").arg(smesh_ct)));
				mesh_item->setData(0, Qt::UserRole + MVIEW_DATA, 
					qVariantFromValue(MeshVariant(MeshVariant::MESH_MODEL_SURF, mesh)));
				mesh_item->setIcon(0, iconView);
				model_item->addChild(mesh_item);
			}
			if(total_volume_mesh){
				int vmesh_ct = total_volume_mesh->getBlocksCount();
				QTreeWidgetItem* mesh_item = 
					new QTreeWidgetItem(model_item, 
							QStringList(tr("volume mesh (NB=%1)").arg(vmesh_ct)));
				mesh_item->setData(0, Qt::UserRole + MVIEW_DATA, 
					qVariantFromValue(MeshVariant(MeshVariant::MESH_MODEL_3D, mesh)));
				mesh_item->setIcon(0, iconView);
				model_item->addChild(mesh_item);
			}

			for(int i = 0; i < mbct; i++){
				MeshDomainVolume* volume = (MeshDomainVolume*)mesh->getBlockAt(i);

				QTreeWidgetItem* volume_item;
				if(volume->getAreaID() >= 0)
					volume_item = new QTreeWidgetItem((QTreeWidget*)0, QStringList(QString("block #%1").arg(
							volume->getAreaID())));
				else
					volume_item = new QTreeWidgetItem((QTreeWidget*)0, QStringList(QString("block")));

				volume_item->setData(0, Qt::UserRole + MVIEW_DATA, 
					qVariantFromValue(MeshVariant(MeshVariant::DOMAIN_3D, volume)));
				model_item->addChild(volume_item);

				CS3dPtr ucs3d = volume->getUserControlSpace();
				if(ucs3d){
					QTreeWidgetItem* cs_item = ucs3d->isAdaptive() ?
						new QTreeWidgetItem(volume_item, QStringList(tr("user ACS3d (NP=%1)").arg(
							ucs3d->getControlNodesCount()))) :
						new QTreeWidgetItem(volume_item, QStringList(tr("user CS3d")));
					if(ucs3d->isAdaptive())
						cs_item->setData(0, Qt::UserRole + MVIEW_DATA, 
							qVariantFromValue(MeshVariant(MeshVariant::CONTROL_3D, ucs3d.get())));
					volume_item->addChild(cs_item);
				}

				CS3dPtr cs3d = volume->getControlSpace();
				if(cs3d){
					QTreeWidgetItem* cs_item = cs3d->isAdaptive() ? 
						new QTreeWidgetItem(volume_item, QStringList(tr("ACS3d (NP=%1)").arg(
							cs3d->getControlNodesCount()))) :
						new QTreeWidgetItem(volume_item, QStringList(tr("CS3d")));
					if(cs3d->isAdaptive())
						cs_item->setData(0, Qt::UserRole + MVIEW_DATA, 
							qVariantFromValue(MeshVariant(MeshVariant::CONTROL_3D, cs3d.get())));
					volume_item->addChild(cs_item);
				}

				int vntq = 0;
				for(int i = 0; i < volume->getFaceCount(); i++){
					MeshDomainSurface* ds = (MeshDomainSurface*)volume->getFace(i);
					MeshContainer2d* mesh2d = ds->getMesh();
					if(mesh2d) vntq += mesh2d->getElementsCount();
				}
				if(vntq > 0){
					QTreeWidgetItem* mesh_item = 
						new QTreeWidgetItem(volume_item, QStringList(tr("mesh2d (NTQ=%1)").arg(vntq)));
					mesh_item->setData(0, Qt::UserRole + MVIEW_DATA, 
						qVariantFromValue(MeshVariant(MeshVariant::MESH_SURF, volume)));
					mesh_item->setIcon(0, iconView);

					volume_item->addChild(mesh_item);
				}

				MeshContainer3dSurface* smesh = volume->getSurfaceMesh();
				if(smesh){
					QTreeWidgetItem* mesh_item = 
						new QTreeWidgetItem(volume_item, QStringList(tr("mesh3dsurf (NT=%1)").arg(
								smesh->getFacesCount())));
					mesh_item->setData(0, Qt::UserRole + MVIEW_DATA, 
						qVariantFromValue(MeshVariant(MeshVariant::MESH_3D_SURF, volume)));
					mesh_item->setIcon(0, iconView);

					volume_item->addChild(mesh_item);
				}

				MeshContainer3d* mesh3d = volume->getMesh();
				if(mesh3d){
					QTreeWidgetItem* mesh_item = 
						new QTreeWidgetItem(volume_item, QStringList(tr("mesh3d (NT=%1)").arg(
								mesh3d->getBlocksCount())));
					mesh_item->setData(0, Qt::UserRole + MVIEW_DATA, 
						qVariantFromValue(MeshVariant(MeshVariant::MESH_3D, mesh3d)));
					mesh_item->setIcon(0, iconView);

					volume_item->addChild(mesh_item);
				}

				QTreeWidgetItem* facegroup_item = 
					new QTreeWidgetItem(volume_item, QStringList(QString("faces (%1)").arg(
							volume->getFaceCount())));
				volume_item->addChild(facegroup_item);

				for(int j = 0; j < volume->getFaceCount(); j++){
					MeshDomainSurface* face = (MeshDomainSurface*)volume->getFace(j);

					QTreeWidgetItem* face_item = 
						new QTreeWidgetItem(facegroup_item, QStringList(QString("face #%1").arg(
								face->getIntTag(TagExtended::TAG_ID))));
					face_item->setData(0, Qt::UserRole + MVIEW_DATA, 
						qVariantFromValue(MeshVariant(MeshVariant::DOMAIN_SURF, face)));
					facegroup_item->addChild(face_item);

					CS2dPtr ucs = face->getUserControlSpace();
					if(ucs){
						QTreeWidgetItem* cs_item = ucs->isAdaptive() ? 
							new QTreeWidgetItem(face_item, QStringList(tr("user ACS2d (NP=%1)").arg(
								ucs->getControlNodesCount()))) :
							new QTreeWidgetItem(face_item, QStringList(tr("user CS2d")));
						if(ucs->isAdaptive()) 
							cs_item->setData(0, Qt::UserRole + MVIEW_DATA, 
								qVariantFromValue(MeshVariant(MeshVariant::CONTROL_2D, ucs.get())));
						face_item->addChild(cs_item);
					}

					MeshContainer2d *face_boundary = face->getBoundary();
					CS2dPtr cs = face_boundary ? face_boundary->getControlSpace() : nullptr;
					if(cs){
						QTreeWidgetItem* cs_item = cs->isAdaptive() ?
							new QTreeWidgetItem(face_item, QStringList(tr("ACS2d (NP=%1)").arg(
								cs->getControlNodesCount()))) :
							new QTreeWidgetItem(face_item, QStringList(tr("CS2d")));
						if(cs->isAdaptive())
							cs_item->setData(0, Qt::UserRole + MVIEW_DATA, 
								qVariantFromValue(MeshVariant(MeshVariant::CONTROL_2D, cs.get())));
						face_item->addChild(cs_item);
					} 

					MeshContainer2d* mesh = face->getMesh(); 
					if(mesh){
						QTreeWidgetItem* mesh_item = 
							new QTreeWidgetItem(face_item, QStringList(tr("mesh2d (NTQ=%1)").arg(
								mesh->getElementsCount())));
						mesh_item->setData(0, Qt::UserRole + MVIEW_DATA, 
							qVariantFromValue(MeshVariant(MeshVariant::MESH_FACE, mesh)));
						mesh_item->setIcon(0, iconView);
						face_item->addChild(mesh_item);
					}
				} 
			}
		}else{
			items.append(new QTreeWidgetItem((QTreeWidget*)0, QStringList() << tr("empty")));
		} 

		LOG_MUTEX("-lm MainWindow::updateModelInfo");
	}

	{
		LOG_MUTEX("+la MainWindow::updateModelInfo");
		QMutexLocker lock(&mutexApp);

		statTree->clear();
		statTree->insertTopLevelItems(0, items);
		statTree->expandAll();
		statTree->resizeColumnToContents(0);
		statTree->collapseAll();
		statTree->setCurrentItem(items.front());
		LOG_MUTEX("-la MainWindow::updateModelInfo");
	}
}

MeshViewSet* MainWindow::genViewSet(const MeshVariant& mv, int part_id)
{
	LOG_MUTEX("+lm MainWindow::genViewSet");
	QMutexLocker lm(console->meshKernelThread.getLockMesh());

	MeshViewSet* set = nullptr;

	if(mv.type == MeshVariant::DOMAIN_MODEL){
		// geometry overview-grid for all blocks
		MeshContainer3d* mesh = (MeshContainer3d*)mv.ptr;
		for(IteratorFace it = mesh->getFirstFace(); it.isValid(); it.nextFace()){
			MeshDomainSurface* domain_face = (MeshDomainSurface*)it.getFace();
			set = domain_face->getViewSet(set);
		}
		for(int i = 0; i < mesh->getBlocksCount(); i++){
			MeshDomainVolume* mdv = (MeshDomainVolume*)mesh->getBlockAt(i);
			if(mdv->getFaceCount() == 0){
				set = mdv->getViewSet(set);
			}
		}
		if(set){
			for(int i = 0; i < mesh->getPointsCount(); i++){
				MeshPoint3d* mp = mesh->getPointAt(i);
				set->addPoint(mp);
			}
			for(int i = 0; i < mesh->getBlocksCount(); i++){
				auto fps = ((MeshDomainVolume*)mesh->getBlockAt(i))->getFreePoints();
				if(fps){
					for(int i = 0; i < fps->count(); i++){
						set->addPoint(fps->get(i)->getCoordinates(), 4);
					}
				}
			}
		}
	}else if(mv.type == MeshVariant::DOMAIN_3D){
		// geometry overview-grid for the selected domain-block
		MeshDomainVolume* mdv = (MeshDomainVolume*)mv.ptr; 
		int dfct = mdv->getFaceCount();
		for(int i = 0; i < dfct; i++){
			MeshDomainSurface* domain_face = (MeshDomainSurface*)mdv->getFace(i);
			set = domain_face->getViewSet(set);
		}
		for(int i = 0; i < mdv->getPointCount(); i++){
			set->addPoint(mdv->getPoint(i));
		}
		if(dfct == 0){
			set = mdv->getViewSet(set);
		}
		auto fps = mdv->getFreePoints();
		if(fps){
			for(int i = 0; i < fps->count(); i++){
				set->addPoint(fps->get(i)->getCoordinates(), 4);
			}
		}
	}else if(mv.type == MeshVariant::DOMAIN_SURF){
		// geometry overview-grid for the selected domain-surface
		MeshDomainSurface* ds = (MeshDomainSurface*)mv.ptr;
		set = ds->getViewSet(set);
		if(set){
			for(int i = 0; i < ds->getEdgeCount(); i++){
				MeshPoint3d* mp = ds->getPoint(i);
				set->addPoint(mp);
			}
		}
	}else if(mv.type == MeshVariant::MESH_MODEL_3D){
		// mesh3d for all domain-blocks
		LOG4CPLUS_INFO(MeshLog::logger_console, "genViewSet MeshVariant::MESH_MODEL_3D");
		MeshContainer3d* mesh = (MeshContainer3d*)mv.ptr;
		MeshContainer3d* vmesh = mesh->getTotalMDV()->getMesh();
		if(vmesh){
			set = vmesh->getViewSet(set, part_id);
			LOG4CPLUS_INFO(MeshLog::logger_console, "set for vmesh created");
		}
	}else if(mv.type == MeshVariant::MESH_MODEL_SURF){
		// mesh3d for all domain-blocks
		LOG4CPLUS_INFO(MeshLog::logger_console, "genViewSet MeshVariant::MESH_MODEL_SURF");
		MeshContainer3d* mesh = (MeshContainer3d*)mv.ptr;
		MeshContainer3dSurface* smesh = mesh->getTotalMDV()->getSurfaceMesh();
		if(smesh){
			set = smesh->getViewSet(set);
			LOG4CPLUS_INFO(MeshLog::logger_console, "set for smesh created");
		}
	}else if(mv.type == MeshVariant::MESH_3D){
		// mesh3d for the selected domain-block
		MeshContainer3d* mesh3d = (MeshContainer3d*)mv.ptr;
		set = mesh3d->getViewSet(set, part_id);
	}else if(mv.type == MeshVariant::MESH_3D_SURF){
		// surface mesh 3d for the selected domain-block
		MeshDomainVolume* block = (MeshDomainVolume*)mv.ptr; 
		MeshContainer3dSurface* smesh = block->getSurfaceMesh();
		if(smesh)
			set = smesh->getViewSet(set);
	}else if(mv.type == MeshVariant::MESH_SURF){
		// surface mesh for the selected domain-block
		MeshDomainVolume* block = (MeshDomainVolume*)mv.ptr; 
		for(int i = 0; i < block->getFaceCount(); i++){
			MeshDomainSurface* ds = (MeshDomainSurface*)block->getFace(i);
			MeshContainer2d* mesh2d = ds->getMesh();
			if(mesh2d){
				bool proper_orientation = (ds->getBlock(0) == block);
				set = mesh2d->getViewSet(set, true, proper_orientation);
			}
		}
	}else if(mv.type == MeshVariant::MESH_FACE){
		// surface mesh for the selected domain-face
		MeshContainer2d* mesh2d = (MeshContainer2d*)mv.ptr;
		set = mesh2d->getViewSet(set);
	}

	LOG_MUTEX("-lm MainWindow::genViewSet");
	return set;
}

void MainWindow::statTreeSelectionChanged()
{
	QList<QTreeWidgetItem*> items = statTree->selectedItems();
	bool available = false;
	MeshViewSet* view_set = nullptr;
	MeshContainer3d* view_mesh = nullptr;

	if(items.count() > 0){
		QTreeWidgetItem* item = items.front();
		QVariant v = item->data(0, Qt::UserRole + MVIEW_DATA);
		available = !v.isNull() && v.isValid();
		if(available){
			MeshVariant mv = v.value<MeshVariant>();
			if(MeshViewSet::param_show_visualization != 0){
				if(mv.type == MeshVariant::MESH_MODEL_3D || mv.type == MeshVariant::MESH_3D ){ // show directly - without snapshot
					// mesh3d for all domain-blocks
					view_mesh = (MeshContainer3d*)mv.ptr;
				}else{
					// get part_id, if any
					int part_id = -2;
					QVariant vpart = item->data(0, Qt::UserRole + MVIEW_PART);
					if(vpart.isValid() && !vpart.isNull())
						part_id = vpart.toInt();
					// create view set
					view_set = genViewSet(mv, part_id);
				}
			}
		}
	}
	
	itemViewAct->setEnabled(available);
	itemStatAct->setEnabled(available);
//	itemClearAct->setEnabled(available);
	if(view_mesh)
		modelGLView->setViewMesh(view_mesh);
	else
		modelGLView->setViewSet(view_set);
}

void MainWindow::showMeshItemView()
{
	QList<QTreeWidgetItem*> items = statTree->selectedItems();
	if(items.count() < 1) return;

	QTreeWidgetItem* item = items.front();
	QVariant v = item->data(0, Qt::UserRole + MVIEW_DATA);
	if(v.isNull() || !v.isValid()) return;

	int part_id = -2;
	QVariant vpart = item->data(0, Qt::UserRole + MVIEW_PART);
	if(vpart.isValid() && !vpart.isNull())
		part_id = vpart.toInt();

	MeshVariant mv = v.value<MeshVariant>();

	if(mv.type == MeshVariant::MESH_MODEL_3D || mv.type == MeshVariant::MESH_3D ){ // show directly - without snapshot
		// mesh3d for all domain-blocks
		MeshContainer3d* mesh = (MeshContainer3d*)mv.ptr;
		int bct = 0;
		int pct = 0;
		MeshContainer3d* mesh3d = mesh->getTotalMDV()->getMesh();
		if(mesh3d){
			bct = mesh3d->getBlocksCount();
			pct = mesh3d->getPointsCount();
		}else{
			for(int i = 0; i < mesh->getBlocksCount(); i++){
				MeshBlock* block = mesh->getBlockAt(i);
				if(block->getType() == BLOCK_DOMAIN){
					mesh3d = ((MeshDomainVolume*)block)->getMesh();
					if(mesh3d){
						bct += mesh3d->getBlocksCount();
						pct += mesh3d->getPointsCount();
					}
				}else{
					bct++;
				}
			}
			if(pct == 0) pct = mesh->getPointsCount();
		}

		QString info = item->text(0) + "\n -> ";
		info += tr("%1 blocks ").arg(bct);
		info += tr("%1 nodes ").arg(pct);

		viewDialog->setViewMesh(info, mesh, 0);

	}else{

		MeshViewSet* set = genViewSet(mv, part_id);

		if(set){
			QString info;
			if(set->m_blocks.count() > 0)
				info += tr("%1 blocks ").arg(set->m_blocks.count());
			if(set->m_faces.count() > 0)
				info += tr("%1 faces ").arg(set->m_faces.count());
			if(set->m_edges.count() > 0)
				info += tr("%1 edges ").arg(set->m_edges.count());
			if(set->m_points.count() > 0)
				info += tr("%1 nodes ").arg(set->m_points.count());

			if(info.isEmpty())
				info = item->text(0);
			else
				info = item->text(0) + "\n -> " + info;
			viewDialog->setViewSet(info, set, 0);
		}else{
			addLogItem("View mesh item - nothing to show", log4cplus::ERROR_LOG_LEVEL);
		}
	}
}

void MainWindow::clearMeshItem()
{
	addLogItem("MainWindow::clearMeshItem", log4cplus::INFO_LOG_LEVEL);
}

void MainWindow::showMeshItemStat()
{
	addLogItem("MainWindow::showMeshItemStat", log4cplus::INFO_LOG_LEVEL);
}

PropertiesWidget::PropertiesWidget(MeshKernelThread* mkt, QWidget * parent, Qt::WindowFlags flags)
	: QDockWidget(tr("Properties"), parent, flags), meshKernelThread(mkt)
{
    prefTable = new QTableWidget(this);
	prefTable->setColumnCount(2);
	prefTable->setHorizontalHeaderLabels(QStringList() << tr("Property") << tr("Value"));
    setWidget(prefTable);

	fillPrefList();
	connect(prefTable, SIGNAL(cellChanged(int, int)),
            this, SLOT(updateMeshProperty(int, int)));
}


ConsoleWidget::ConsoleWidget(MeshViewDialog* meshViewDlg, QWidget *parent)
    : viewDlg(meshViewDlg), QWidget(parent), meshKernelThread(this)
{
	// create cmdEdit
	cmdEdit = new QLineEdit;
	QStringList preCmdList;
	preCmdList << "exit" << "load" << "load_and_prepare" 
		<< "parse_grain" << "quit" << "store_grd" << "store_txt" 
		<< "test" << "triangulate" << "triangulate3";
	QCompleter *completer = new QCompleter(preCmdList, this);
	completer->setCaseSensitivity(Qt::CaseInsensitive);
	completer->setCompletionMode(QCompleter::InlineCompletion);
	cmdEdit->setCompleter(completer);
	// create applyButton
	applyButton = new QPushButton(tr("Run"));
	applyButton->setDefault(true);
	applyButton->setEnabled(false);
	// connect the two above
	connect(cmdEdit, SIGNAL(textChanged(const QString &)),
			this, SLOT(enableApplyButton(const QString &)));
	connect(applyButton, SIGNAL(clicked()),
			this, SLOT(applyClicked()));
	connect(cmdEdit, SIGNAL(returnPressed()),
			this, SLOT(applyClicked()));
	// create chLayout for the two above
	QHBoxLayout *chLayout = new QHBoxLayout;
	chLayout->addWidget(cmdEdit);
	chLayout->addWidget(applyButton);
	// create cmdList
	cmdList = new QListWidget;
	// connect some signals
	connect(cmdList, SIGNAL(itemClicked(QListWidgetItem *)),
			this, SLOT(cmdListClicked(QListWidgetItem *)));
	connect(cmdList, SIGNAL(itemDoubleClicked(QListWidgetItem *)),
			this, SLOT(cmdListDoubleClicked(QListWidgetItem *)));
	// create cvLayout for the two above
	QVBoxLayout *cvLayout = new QVBoxLayout;
	cvLayout->addWidget(cmdList);
	cvLayout->addLayout(chLayout);
	// apply layout
	setLayout(cvLayout);

	// start kernel 
	meshKernelThread.start();
}

void ConsoleWidget::applyClicked()
{
	QString text = cmdEdit->text();
	if(!text.isEmpty()){
		cmdEdit->clear();
		addMeshCommand(text);
	}
}

void ConsoleWidget::updateCmdState(int which, int state)
{
	LOG_MUTEX("+lc ConsoleWidget::updateCmdState");
	QMutexLocker lc(meshKernelThread.getLockCmd());
	toshowCmdWhich.append(which);
	toshowCmdState.append(state);
	if(state == 1){ // running new command
		viewDlg->skipAll = false;
	}
	LOG_MUTEX("-lc ConsoleWidget::updateCmdState");
}

void ConsoleWidget::updateCmdView()
{
	static const QIcon icons[] = {
		QIcon(":/images/rect_green.png"),
		QIcon(":/images/right_arrow.png"),
		QIcon(":/images/rect_blue.png"),
		QIcon(":/images/rect_red.png"),
		QIcon(":/images/rect_gray.png")};

	int which = -1;
	int state = -1;

	while(true){
		{
			QMutexLocker lc(meshKernelThread.getLockCmd());
			if(toshowCmdWhich.empty()) break;

			which = toshowCmdWhich.takeFirst();
			state = toshowCmdState.takeFirst();
		}
		if(which < cmdList->count()){
			QListWidgetItem* item = cmdList->item(which);
			item->setIcon(icons[state]);

			//if(state > 1){
			//	TSTATUS("Ready");
			//}else{
			//	TSTATUS(string("Running: ") + item->text().toStdString());
			//}
		}
	}
	if(which > -1 && which < cmdList->count())
		cmdList->setCurrentRow(which);
}

void ConsoleWidget::enableApplyButton(const QString &text)
{
	applyButton->setEnabled(!text.isEmpty());
}

void ConsoleWidget::addMeshCommand(const QString& cmd)
{
	cmdList->addItem(cmd);
	meshKernelThread.addCommand(cmd);

	cmdHistory.removeAll(cmd);
	cmdHistory.append(cmd);
	if(cmdHistory.count() > 100)	cmdHistory.removeFirst();
}

void ConsoleWidget::cmdListClicked(QListWidgetItem * item)
{
	cmdEdit->setText(item->text());
}

void ConsoleWidget::cmdListDoubleClicked(QListWidgetItem * item)
{
	addMeshCommand(item->text());
}

void ConsoleWidget::setHistory(const QStringList& history)
{
	cmdHistory = history;
	meshKernelThread.addHistoryCommands(history);
	cmdList->addItems(history);
	for(int i = 0; i < history.count(); i++)
		updateCmdState(i, 4);
}

MeshKernelThread::MeshKernelThread(ConsoleWidget* cw)
	: console(cw), stopped(false), processedCmd(-1), lastCmd(-1), availSemaphore(0)
{
}

void MeshKernelThread::run()
{
//	emit modelChanged(meshModel);
	while(true) {
		availSemaphore.acquire();
		if(stopped) break;
		++processedCmd;

		LOG_MUTEX("+lc MeshKernelThread::run");
		lockCmd.lock();
		QString cmd = commandList.takeFirst();
		LOG_MUTEX("-lc MeshKernelThread::run");
		lockCmd.unlock();

		if(!cmd.isEmpty()){
			console->updateCmdState(processedCmd, 1);
			LOG4CPLUS_DEBUG(MeshLog::logger_console, cmd.toStdString());

			int res;
			double sec;

			{
				LOG_MUTEX("+lm MeshKernelThread::run");
				QMutexLocker lm(getLockMesh());
				clock_t clock_start = clock();
				res = model.execute(cmd.toStdString());
				sec = (clock() - clock_start)/(double)CLOCKS_PER_SEC;
				LOG_MUTEX("-lm MeshKernelThread::run");
			}

			if(res == MeshModel::CM_OK){
				auto info = tr("%1\t%2").arg(sec).arg(cmd);
				LOG4CPLUS_INFO(MeshLog::logger_console, info.toStdString());
				console->updateCmdState(processedCmd, 2);
			}else if(res == MeshModel::CM_QUIT){
				// ???
				LOG4CPLUS_INFO(MeshLog::logger_console, "Quitting...");
				QApplication::exit();
			}else{
				console->updateCmdState(processedCmd, 3);
			}
		}
	}
}

void MeshKernelThread::addHistoryCommands(const QStringList& history)
{
	processedCmd += history.count();
	lastCmd += history.count();
}

void MeshKernelThread::stop()
{
	LOG4CPLUS_WARN(MeshLog::logger_console, 
		"Waiting for meshing kernel to stop...");
	LOG_MUTEX("+lc MeshKernelThread::stop");
	lockCmd.lock();
	stopped = true;
	LOG_MUTEX("-lc MeshKernelThread::stop");
	lockCmd.unlock();
	availSemaphore.release();
}

void MainWindow::closeEvent(QCloseEvent *event)
{
	viewDialog->closing = true;

	console->meshKernelThread.stop();
	if(viewDialog->waitingForUser)
		viewDialog->viewWaitSemaphore.release(100);
	console->meshKernelThread.wait();

	writeSettings();
	event->accept();
}

int MeshKernelThread::addCommand(const QString& cmd)
{
	LOG_MUTEX("+lc MeshKernelThread::addCommand");
	lockCmd.lock();
	commandList.append(cmd);
	LOG_MUTEX("-lc MeshKernelThread::addCommand");
	lockCmd.unlock();
	console->updateCmdState(++lastCmd, 0);
	availSemaphore.release();
	return lastCmd;
}

MeshViewDialog::MeshViewDialog(QWidget *parent)
	: QDialog(parent), viewWaitSemaphore(0), waitingForUser(false), closing(false), skipAll(false)
{
	meshGLView = new ModelGLWidget(parent);

	closeButton = new QPushButton(tr("Close"));
	//closeButton->setDefault(true);
	connect(closeButton, SIGNAL(clicked()),
			this, SLOT(close()));

	continueButton = new QPushButton(tr("Continue"));
	continueButton->hide();
	connect(continueButton, SIGNAL(clicked()),
			this, SLOT(continuePressed()));

	skipButton = new QPushButton(tr("Skip all"));
	connect(skipButton, SIGNAL(clicked()),
			this, SLOT(skipAllPressed()));

	logButton = new QPushButton(tr("Log"));
	connect(logButton, SIGNAL(clicked()),
			this, SLOT(logPressed()));

	description = new QTextEdit;
	description->setReadOnly(true);
	description->setFixedHeight(100);

	QHBoxLayout *chLayout = new QHBoxLayout;
	chLayout->addWidget(description);

	infoTable = new QTableWidget(this);
	infoTable->setColumnCount(2);
	infoTable->setHorizontalHeaderLabels( QStringList() << tr("Property") << tr("Value") );
	infoTable->setFixedHeight(100);

	chLayout->addWidget(infoTable);
	

	QVBoxLayout *cvbLayout = new QVBoxLayout;
	cvbLayout->addWidget(continueButton);
	cvbLayout->addWidget(closeButton);
	cvbLayout->addWidget(skipButton);
	cvbLayout->addWidget(logButton);
	cvbLayout->addStretch();
	chLayout->addLayout(cvbLayout);

	QVBoxLayout *cvLayout = new QVBoxLayout;
	cvLayout->addWidget(meshGLView);
	cvLayout->setStretchFactor(meshGLView, 1);
	cvLayout->addLayout(chLayout);

	// apply layout
	setLayout(cvLayout);

	setWindowTitle(tr("Mesh View GL"));

	QSettings settings("KI AGH", "QMeshGen");

	QRect qr(300, 300, 400, 400);
	QRect frect = settings.value("mvd-fgeometry", qr).toRect();
	move(frect.topLeft());
	QRect rect = settings.value("mvd-geometry", qr).toRect();
	resize(rect.size());

}

void MeshViewDialog::fillInfoTable(const DataVector<string> & info)
{
	infoTable->clearContents();
	int info_ct = info.count() / 2;
	infoTable->setRowCount(info_ct);
	for(int i = 0; i < info_ct; i++){
		QTableWidgetItem *item = new QTableWidgetItem(tr(info[2*i].c_str()));
		//item->setFlags(Qt::ItemIsEnabled);
		infoTable->setItem(i, 0, item);

		item = new QTableWidgetItem(tr(info[2*i+1].c_str()));
		//item->setFlags(Qt::ItemIsEditable | Qt::ItemIsEnabled);
		infoTable->setItem(i, 1, item);
		infoTable->setRowHeight(i, 18);
	}
	for(int i = 0; i < 2; i++) 
		infoTable->resizeColumnToContents(i);
}

void MeshViewDialog::setViewSet(const QString& desc, MeshViewSet* set, long t, bool with_reset)
{
	//if(skipAll){
	//	delete set;
	//	return;
	//}

	assert(set);
	fillInfoTable(set->m_info);
	description->setPlainText(desc);
	meshGLView->setViewSet(set, with_reset);
	if(t != 0){
		closeButton->hide();
		continueButton->show();
		continueButton->setDefault(true);
	}
	this->show();
}

void MeshViewDialog::setViewMesh(const QString& desc, MeshContainer3d* mesh, long t)
{
	//if(skipAll){
	//	delete set;
	//	return;
	//}

	DataVector<string> info;
	// mesh->getInfo(info); ???
	fillInfoTable(info);
	description->setPlainText(desc);
	meshGLView->setViewMesh(mesh);
	if(t != 0){
		closeButton->hide();
		continueButton->show();
		continueButton->setDefault(true);
	}
	this->show();
}

void MeshViewDialog::closeEvent(QCloseEvent *event)
{
	if(waitingForUser)
		viewWaitSemaphore.release();
	continueButton->hide();
	closeButton->show();
	closeButton->setDefault(true);

	QSettings settings("KI AGH", "QMeshGen");
	settings.setValue("mvd-fgeometry", frameGeometry());
	settings.setValue("mvd-geometry", geometry());

	event->accept();
}

void MeshViewDialog::continuePressed()
{
	if(waitingForUser)
		viewWaitSemaphore.release();
	continueButton->hide();
	closeButton->show();
	closeButton->setDefault(true);
}

void MeshViewDialog::skipAllPressed()
{
	continuePressed();
	skipAll = true;
	this->hide();
}

void MeshViewDialog::logPressed()
{
	LOG4CPLUS_INFO(MeshLog::logger_mesh, 
		description->toPlainText().toStdString());
	int info_ct = infoTable->rowCount();
	for(int i = 0; i < info_ct; i++){
		auto info = tr("\t%1\t%2").arg(infoTable->item(i, 0)->text()).arg(infoTable->item(i, 1)->text());
		LOG4CPLUS_INFO(MeshLog::logger_mesh, info.toStdString());
	}
}