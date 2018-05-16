#ifndef SPHGRID_H
#define SPHGRID_H
#pragma once
#include <SPHParticle.h>
#include <ngl/Vec3.h>
#include <unordered_map>
#include <vector>
#include <string>
#include <ngl/AbstractVAO.h>
#include <map>
#include <math.h>

struct SPHVariables
{
public:

    float PI = 3.141592653589f;

    // Set Variables
    float m_Kernel = 0.0125f;
    float m_Kernel2 = m_Kernel*m_Kernel;

    float m_WallDampening = 0.5f; // Halves the velocity when striking a wall
    float m_RestDensity = 998.2071f; // Resting density of water
    float m_GasConstant = 1.0f;//8.314f; // Basic gas constant
    float m_Viscosity = 6.0f; // Basic viscosity scalar
    float m_TimeStep = 0.03f; // Our dt time step
    float m_Mass = 0.02f; // Basic mass of water (can be variable, but not for this sim)

    // Kernels
    float m_Poly6Value = 315.0f/(64.0f * PI * pow(m_Kernel, 9));
    float m_Poly6Grad = -945.0f/(32.0f * PI * pow(m_Kernel, 9));
    float m_Poly6Lap = -945.0f/(8.0f * PI * pow(m_Kernel, 9));

    float m_SpikyValue = -45.0f/(PI * pow(m_Kernel, 6));
    float m_ViscosityValue = 45.0f/(PI * pow(m_Kernel, 6));
    float m_TensionLap = m_Kernel2 * (0-((3/4)*m_Kernel2)) * m_Poly6Lap * m_Mass;

    float m_RDens = m_Mass * m_Poly6Value * pow(m_Kernel, 6);
    float m_SurfCoef = 0.5f;
    float m_SurfTenValue = 6.0f; // Default for water

    ngl::Vec3 m_Gravity = ngl::Vec3(0.0f,-9.8f,0.0f);



};

class SPHGrid
{
public:
    SPHGrid(int _numParticles, ngl::Vec3 _dim, float _cellSize);

    void Update();
    void Draw(const ngl::Mat4 &_rot);
    void AddToMap(int _pID);
    void AddToMap(int _pID, int _cID);

    const std::string GetShaderName(){return m_shaderName;}
    void setShaderName(const std::string &_n){m_shaderName=_n;}
    void setCam(ngl::Camera *_cam){m_cam=_cam;}

    void FindNeighbours();

    bool IsPrime();
    int NextPrime();

    void MoveByCellID(int _cellID, ngl::Vec3 _dir);
    void MoveByCell(int3 _cell, ngl::Vec3 _dir);
    void MoveNeighbours(int _ParticleID, ngl::Vec3 _dir);
    void CalcDensity(int _ParticleID);
    void CalcPressure(int _ParticleID);
    void CalcViscosity(int _ParticleID);
    float CreatePoly6(float _h, float _r);

    void WriteToFile();

    bool m_Update=true;

private:
    ngl::Vec3 m_pos;
    ngl::Vec3 m_dim;
    ///@brief this will be our mapping that will
    /// take the id of a particle and assign it to
    /// a cell id on the grid
    std::unordered_map<int, std::vector<int>> m_map;
    std::vector<int> m_CellIDs;

    std::unique_ptr<SPHParticle []>  m_Particles;
    std::unique_ptr<GLParticle []> m_Glparticles;

    std::unique_ptr<ngl::AbstractVAO> m_vao;

    std::string m_shaderName;

    int m_NumParticles;
    float m_CellSize;

    ngl::Camera *m_cam;

    SPHVariables SVar;




    int frame, maxFrame;
    float k=1;

    int m_FrameTick = 0;

    int m_Prime1;
    int m_Prime2;
    int m_Prime3;


    float ZVAL = 0.0000000001f;




};

// https://gamedev.stackexchange.com/questions/56590/c-how-to-implement-a-spatial-hash-for-a-2d-game

#endif // SPHGRID_H

