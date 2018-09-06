#version 330
in highp vec4 vColor;
out highp vec4 fColor;
uniform vec3 colorBand;
uniform vec3 colorInfo;
float noc;
float max;
float min;

void main()
{
    fColor = vColor;
}
