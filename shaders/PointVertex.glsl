#version 330 core
/// @brief the vertex passed in
layout(location=0)in vec3 inVert;
/// @brief the colour of the point
layout(location=1)in vec3 inColour;
uniform mat4 MVP;

out vec3 pointColour;

void main()
{
  pointColour=inColour;
  gl_Position = MVP * vec4(inVert,1);
}
