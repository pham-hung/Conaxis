#ifndef LOADRESULT_H
#define LOADRESULT_H

#include <QDialog>
#include <QFileDialog>

namespace Ui {
class LoadResult;
}

class LoadResult : public QDialog
{
    Q_OBJECT

public:
    explicit LoadResult(QWidget *parent = 0);
    ~LoadResult();    


private slots:

    void on_ULoadButton_clicked();
    void on_VLoadButton_clicked();
    void on_PLoadButton_clicked();
    void on_okButton_clicked();

    void on_cancelButton_clicked();

signals:
    void sendFileName(QString,QString,QString,bool);

private:
    Ui::LoadResult *ui;
    QString UfileName;
    QString VfileName;
    QString PfileName;
    bool ok=false;
};

#endif // LOADRESULT_H
