#ifndef SPHPARTICLE_H
#define SPHPARTICLE_H
#pragma once
#include <vector>
#include <ngl/Camera.h>
#include <ngl/AbstractVAO.h>
#include <memory>

#pragma pack(push,1)
struct GLParticle
{
  GLfloat px;
  GLfloat py;
  GLfloat pz;
  GLfloat pr;
  GLfloat pg;
  GLfloat pb;

};
#pragma pack(pop)

struct int3
{
    int3(){x=0;y=0;z=0;}
    int3(int _x, int _y, int _z){x=_x; y=_y; z=_z;}
    int x,y,z;
};

class SPHParticle
{
public:
    SPHParticle();

    void Draw();
    void Update();

    void UpdatePressure();

    // Testing placement within a grid
    float RandPoint(float a, float b);
    ngl::Vec3 RandPosSphere(float _r, float _centre);

    std::vector<int> m_Neighbours;

    // Neighbourhood search radius
    float m_H;

    int m_id;
    int m_CellID;
    int3 m_cell;

    ngl::Vec3 m_pos;
    ngl::Vec3 m_vel;
    ngl::Vec3 m_acc;
    ngl::Vec3 m_adv;

    ngl::Vec3 m_col;

    float m_Gravity;
    float m_Mass = 0.02f;

    float m_Density;
    float m_Pressure;
    ngl::Vec3 m_PressureForce;
    float m_Viscosity;
    float m_SurfaceTension;

};


#endif // SPHPARTICLE_H
