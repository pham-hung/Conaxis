#version 330
layout(location = 8) in vec3 position;
layout(location = 9) in vec3 color;
out vec4 vColor;
uniform mat4 modelToWorld;
uniform mat4 worldToCamera;
uniform mat4 cameraToView;

void main()
{
  gl_Position = cameraToView * worldToCamera * modelToWorld * vec4(position, 1.0);
  vColor = vec4(0.0f,0.5f,0.0f,1.0f);
}
