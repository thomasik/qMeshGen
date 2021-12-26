
#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QThread>
#include <QWidget>
#include <QMainWindow>
#include <QString>
#include <QSemaphore>
#include <QMutex>
#include <QList>
#include <QStringList>
#include <QDialog>
#include <QVariant>
#include <QDockWidget>
#include <QDialog>
#include <QGLFormat>

// from meshlib
#include "MeshLog.h"
#include "MeshViewSet.h"
#include "MeshModel.h"

QT_BEGIN_NAMESPACE
class QAction;
class QListWidget;
class QListWidgetItem;
class QTreeWidget;
class QTableWidget;
class QLineEdit;
class QPushButton;
class QMenu;
class QTextEdit;
class QCheckBox;
QT_END_NAMESPACE

class ConsoleWidget;
class ModelGLWidget;
class MainWindow;

class MeshKernelThread : public QThread
{
	Q_OBJECT

public:
	MeshKernelThread(ConsoleWidget* cw);

public:
	void stop();
	int addCommand(const QString& cmd);
	void addHistoryCommands(const QStringList& history);
	MeshModel& getMeshModel() { return model; }
	QMutex* getLockCmd() { return &lockCmd; }
	QMutex* getLockMesh() { return &lockMesh; }
protected:
	void run();
private:
	ConsoleWidget* console;
	volatile bool stopped;
	int processedCmd;
	int lastCmd;
	QSemaphore availSemaphore;
	QMutex lockCmd;
	QMutex lockMesh;
	QStringList commandList;
	MeshModel model;
};

class ConsoleQtAppender : public Appender
{
public:
	ConsoleQtAppender(MainWindow* mw);
	virtual void close();
//	static void registerAppender();
protected:
	virtual void append(const spi::InternalLoggingEvent& event);
private:
	MainWindow * m_mw;
};

class MeshViewDialog : public QDialog
{
    Q_OBJECT

public:
    MeshViewDialog(QWidget *parent = 0);
	void fillInfoTable(const DataVector<string> & info);
	void setViewSet(const QString& desc, MeshViewSet* set, long t = -1, bool with_reset = true);
	void setViewMesh(const QString& desc, MeshContainer3d* mesh, long t = -1);

protected:
	void closeEvent(QCloseEvent *event);

private slots:
    void continuePressed();
	void skipAllPressed();
	void logPressed();

public:
	QSemaphore viewWaitSemaphore;
	bool waitingForUser;
	bool closing;
	bool skipAll;

private:
	QPushButton *closeButton;
	QPushButton *continueButton;
	QPushButton *skipButton;
	QPushButton *logButton;
    ModelGLWidget *meshGLView;
	QTextEdit* description;
	QTableWidget *infoTable;
};

class MeshPropertiesDialog : public QDialog
{
    Q_OBJECT

public:
    MeshPropertiesDialog(QWidget *parent, MeshKernelThread* mkt);

private slots:
	void apply();

private:
	QTableWidget *propTable;
	MeshKernelThread* meshKernelThread;
	QStringList prop_labels;
	QStringList prop_names;
	QStringList prop_values;
};

class ConsoleWidget : public QWidget
{
    Q_OBJECT

public:
    ConsoleWidget(MeshViewDialog* meshViewDlg, QWidget *parent = 0);
	void updateCmdState(int which, int state);
	void updateCmdView();
	void setHistory(const QStringList& history);

public slots:
	void addMeshCommand(const QString& cmd);

private slots:
	void enableApplyButton(const QString &text);
	void applyClicked();
	void cmdListClicked(QListWidgetItem * item);
	void cmdListDoubleClicked(QListWidgetItem * item);

public:
	MeshKernelThread meshKernelThread;
	QStringList cmdHistory;

private:
	QList<int> toshowCmdWhich;
	QList<int> toshowCmdState;
private:
	QLineEdit   *cmdEdit;
	QListWidget *cmdList;
	QPushButton *applyButton;
	MeshViewDialog* viewDlg;
};

struct MeshVariant {
	enum PtrType { UNKNOWN, 
		MESH_FACE, MESH_SURF, MESH_3D, MESH_3D_SURF, MESH_MODEL_SURF, MESH_MODEL_3D,
		DOMAIN_SURF, DOMAIN_3D, DOMAIN_MODEL, 
		CONTROL_2D, CONTROL_3D };
	MeshVariant(PtrType t = UNKNOWN, void* p = 0) : type(t), ptr(p) {}
	PtrType type;
	void* ptr;
};
Q_DECLARE_METATYPE(MeshVariant);

class PropertiesWidget : public QDockWidget
{
    Q_OBJECT

public:
	PropertiesWidget(MeshKernelThread* mkt, QWidget * parent = 0, Qt::WindowFlags flags = 0);

	void fillPrefList();

private slots:
	void updateMeshProperty(int row, int column);

protected:
	void resizeEvent ( QResizeEvent * event );
public:
    QTableWidget *prefTable;
private:
	MeshKernelThread* meshKernelThread;
};

class MainWindow : public QMainWindow, public MeshViewExt
{
    Q_OBJECT

public:
    MainWindow(int argc, char* argv[]);
public:
	// from MeshLogExt interface
	virtual void message(const string& text, log4cplus::LogLevel msg_type);
	// from MeshLogExt interface
	virtual void status(const string& text);
	// from MeshViewExt interface
	virtual void showViewSet(const string& desc, MeshViewSet* set, int t);
	virtual void showViewSetNoReset(const string& desc, MeshViewSet* set, int t);

//	void invalidateAuthentication() { validation_result = false; }

protected:
	void closeEvent(QCloseEvent *event);

	void timerEvent(QTimerEvent *event);
	void showEvent(QShowEvent *event);
	void hideEvent(QHideEvent *event);

	void updateLogView();

private slots:
    void open();
	void openRecentFile();
	void saveGRD();
	void saveTXT();
    void about();
	void updateModelInfo();
	// for statTree
    void showMeshItemView();
    void showMeshItemStat();
    void clearMeshItem();
	void statTreeSelectionChanged();

	void triangulateAction();
	void quadAction();
	void tetrahedraAction();

	void parseGrainAction();
	void meshingPropAction();

private:
    void createActions();
    void createMenus();
    void createToolBars();
    void createStatusBar();
    void createDockWindows();

	void readSettings();
	void writeSettings();
	void updateRecentFileActions();
	void loadFile(const QString& fname);
	void loadScript(const QString& fname);
	void checkParams(int& argc, char* argv[]);
	void initLog4WithQt(const string& log_name = "qtmesh");

	MeshViewSet* genViewSet(const MeshVariant& mv, int part_id = -2);

	QStringList recentFiles;
	enum { MaxRecentFiles = 5 };
	QAction *recentFileActions[MaxRecentFiles];
	QAction *separatorAction;

    ModelGLWidget *modelGLView;
    QTreeWidget *statTree;
	PropertiesWidget* prefDock;
	QListWidget *logList;
	MeshViewDialog* viewDialog;
	MeshPropertiesDialog* propDialog;

	ConsoleWidget *console;
	int tId;
	string user_id;

	QMutex mutexApp;

private:
	QString     toshowStatus;
	struct LogMsg {
		LogMsg(const QString& _text, log4cplus::LogLevel _type) : text(_text), type(_type) {}
		LogMsg() : text(""), type(log4cplus::INFO_LOG_LEVEL) {}
		QString text;
		log4cplus::LogLevel type;
	};
	void addLogItem(const QString& _text, log4cplus::LogLevel _type, bool update_position = true);

	QList<LogMsg>  toshowMsg;
	struct ViewMesh {
		ViewMesh() : desc("Mesh set"), set(nullptr), t(-1), with_reset(true) {}
		void setEmpty() { set = nullptr; }
		QString desc;
		MeshViewSet* set;
		int t;
		bool with_reset;
	};
	ViewMesh toshowView;

private:  // actions for statTree
	QAction* itemViewAct;
	QAction* itemStatAct;
//	QAction* itemClearAct;
	enum StatTreeRole { MVIEW_DATA, MVIEW_PART }; 

private:
    QAction *openAct;
	QAction *saveGRDAct;
	QAction *saveTXTAct;
    QAction *aboutAct;
    QAction *quitAct;
	QAction *genTriangleAct;
	QAction *genQuadAct;
	QAction *genTetraAct;
	QAction *parseGrainAct;
	QAction *propDlgAct;

	QMenu *viewMenu;
};

#endif
