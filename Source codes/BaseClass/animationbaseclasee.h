#ifndef ANIMATIONBASECLASEE_H
#define ANIMATIONBASECLASEE_H


class AnimationBaseClasee
{
public:
    AnimationBaseClasee();
    int beginStep=1;
    int endStep=1;
    int delayTime=0;
    bool lockView=false;
    bool unlockView=true;
    bool autoChangeColor=true;
    bool run=false;
};

#endif // ANIMATIONBASECLASEE_H
