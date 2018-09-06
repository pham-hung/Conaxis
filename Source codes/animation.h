#ifndef ANIMATION_H
#define ANIMATION_H

#include <QDialog>
#include "BaseClass/animationbaseclasee.h"

namespace Ui {
class Animation;
}

class Animation : public QDialog,private AnimationBaseClasee
{
    Q_OBJECT

public:
    explicit Animation(QWidget *parent = 0);
    ~Animation();
    void getUserData();

signals:
    void sendSignal(AnimationBaseClasee animation,bool ok);
private slots:
    void on_okButton_clicked();
    void on_cancelButton_clicked();

private:
    Ui::Animation *ui;
    AnimationBaseClasee animation;
    bool ok=false;
};

#endif // ANIMATION_H
