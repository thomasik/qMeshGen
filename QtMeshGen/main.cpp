#include <QApplication>
#include <iostream>
#include <QtOpenGL>

#include "mainwindow.h"

// QFileSystemWatcher

int main(int argc, char* argv[])
{
	log4cplus::Initializer initializer;

	QApplication app(argc, argv);

	if (!QGLFormat::hasOpenGL()) {
		std::cerr << "This system has no OpenGL support" << endl;
		return 1;
	}

    Q_INIT_RESOURCE(qtmeshgen);

	app.setWindowIcon(QIcon(":/images/green_web.png"));
    MainWindow mainWin(argc, argv);
    mainWin.show();
    auto ret = app.exec();
	return ret;
}