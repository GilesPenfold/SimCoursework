#version 330 core

layout(location=0) out vec4 fragColour;
in vec3 pointColour;
void main ()
{
  // it is important to set the Alpha here as
  // points use it (show what .rgb does in lecture!)
  fragColour=vec4(pointColour,1.0);
}

