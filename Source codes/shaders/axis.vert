#version 330
layout(location = 6) in vec3 position;
layout(location = 7) in vec3 color;
out vec4 vColor;
float gl_PointSize;

uniform mat4 modelToWorld;
uniform mat4 worldToCamera;
uniform mat4 cameraToView;

void main()
{
    gl_Position = cameraToView * worldToCamera * modelToWorld * vec4(position, 1.0);
    vColor = vec4(color,1.0f);    
}
