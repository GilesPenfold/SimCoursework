#include <SPHGrid.h>
#include <ngl/VAOFactory.h>
#include <ngl/SimpleVAO.h>
#include <ngl/ShaderLib.h>
#include <ngl/Random.h>
#include <iostream>
#include <ngl/NGLStream.h>
#include <QElapsedTimer>
#include <ngl/Logger.h>
#include <fstream>
#include <sstream>

/*
 * Based on the SPH implementations of Dongli Zhang
 * https://github.com/finallyjustice/sphfluid
 * And Paul Kennedy
 * https://github.com/PaulKennedyDIT/SPH
 * Using Müllers 2003 paper: Particle-Based Fluid Simulation for Interactive Applications
 */





///
/// \brief SPHGrid::SPHGrid
/// \param _numParticles The number of particles in the simulation
/// \param _dim The dimensions of the containing box
/// \param _cellSize The size of the cells within the box that SPH will divide into
///
SPHGrid::SPHGrid(int _numParticles, ngl::Vec3 _dim, float _cellSize)
{
    // Allocate our input variables
    m_NumParticles = _numParticles;
    m_CellSize = _cellSize;
    m_dim = _dim;

    m_Update=false;

    // Create placeholder variables for particle creation
    SPHParticle p;
    GLParticle g;
    // Allocate space for our particles
    m_Particles.reset( new SPHParticle[m_NumParticles]);
    m_Glparticles.reset( new GLParticle[m_NumParticles]);
    // Init our random function
    ngl::Random *rando=ngl::Random::instance();

    m_CellIDs = std::vector<int>();

    // Spatial Hashing Prime Calculation
    // Uses 3 high numbered primes
    m_Prime1 = 73856093;
    m_Prime2 = 19349663;
    m_Prime3 = 83492791;

    // Set the position of our simulation to the origin
    m_pos = ngl::Vec3(0,0,0);

    frame = 1;
    maxFrame = 700;

    std::cout << "Didn't crash yet!\n";



    m_vao.reset( ngl::VAOFactory::createVAO(ngl::simpleVAO,GL_POINTS));

    // Setup our map for neighbourhood searching
    for(int x = 0; x < _dim.m_x/_cellSize; x++)
    {
        for(int y = 0; y < _dim.m_y/_cellSize; y++)
        {
            for (int z = 0; z < _dim.m_z/_cellSize; z++)
            {
                int CELLID = (m_Prime1 * x) + (m_Prime2 * y) + (m_Prime3 * z);
                m_map[CELLID].reserve(sizeof(SPHParticle));
            }
        }
    }



    // Initialise all of our particles
    for(size_t i=0; i<m_NumParticles; ++i)
    {
        // Initial Setup of Particles

        float boxMod = 4.0f;

        //p.m_pos.m_x=p.RandPoint(0,_dim.m_x/boxMod)+m_pos.m_x;
        //p.m_pos.m_y=p.RandPoint(0,_dim.m_y/boxMod)+m_pos.m_y;
        //p.m_pos.m_z=p.RandPoint(0,_dim.m_z/boxMod)+m_pos.m_z;
        p.m_pos = p.RandPosSphere(_dim.m_x/boxMod, _dim.m_x/2);
        g.px=i+m_pos.m_x;
        g.py=i+m_pos.m_y;
        g.pz=i+m_pos.m_z;
        ngl::Colour c=rando->getRandomColour();
        p.m_col.m_r=g.pr=c.m_r;
        p.m_col.m_g=g.pg=c.m_g;
        p.m_col.m_b=g.pb=c.m_b;


        p.m_vel = ngl::Vec3(p.RandPoint(-10.0f,10.0f),p.RandPoint(-1.0f,5.0f),p.RandPoint(-10.0f,10.0f));
        p.m_acc = ngl::Vec3(0.0f,0.0f,0.0f);
        p.m_adv = ngl::Vec3(0.0f,0.0f,0.0f);
        p.m_Mass=p.RandPoint(0.01f,0.03f);

        p.m_id = i;

        // Search distance variable
        p.m_H = SVar.m_Kernel;

        m_Particles[i]=p;
        m_Glparticles[i]=g;

        // Give some starting velocity
        //p.m_vel.m_x = 25.0f;
        //p.m_vel.m_z = 0.9f;

        // Spatial Hash to find all CELLIDs

        int X = (int)(m_Particles[i].m_pos.m_x/m_CellSize);
        int Y = (int)(m_Particles[i].m_pos.m_y/m_CellSize);
        int Z = (int)(m_Particles[i].m_pos.m_z/m_CellSize);

        int CELLID = m_Particles[i].m_CellID = (m_Prime1 * X) + (m_Prime2 * Y) + (m_Prime3 * Z);
        m_map[CELLID].push_back(i);
        m_Particles[i].m_cell = int3(X,Y,Z);


        m_Particles[i].m_col.m_r = m_CellSize*(float)X/m_dim.m_x;
        m_Particles[i].m_col.m_g = m_CellSize*(float)Y/m_dim.m_y;
        m_Particles[i].m_col.m_b = m_CellSize*(float)Z/m_dim.m_z;

    }

    // Bind our VAO
    m_vao->bind();
    // Input our data
    m_vao->setData(ngl::SimpleVAO::VertexData(m_NumParticles*sizeof(GLParticle),m_Particles[0].m_pos.m_x));
    m_vao->setVertexAttributePointer(0,3,GL_FLOAT,sizeof(GLParticle),0);
    m_vao->setVertexAttributePointer(1,3,GL_FLOAT,sizeof(GLParticle),3);
    m_vao->setNumIndices(m_NumParticles);
    m_vao->unbind();

    // Calculate the local neighbours of each particle - Currently unused as using brute force instead of SPH
    FindNeighbours();


}

///
/// \brief SPHGrid::Update
/// Updates the particles within this grid
void SPHGrid::Update()
{
    m_vao->bind();
    ngl::Real *glPtr=m_vao->mapBuffer();
    unsigned int glIndex=0;

    FindNeighbours();

//    std::vector<SPHParticle*> m_DebugList = std::vector<SPHParticle*>();
//    for(unsigned int i = 0; i < m_NumParticles; i++)
//    {
//        m_DebugList.push_back(&m_Particles[i]);
//    }

    #pragma omp parallel for
    for(unsigned int i=0; i<m_NumParticles; ++i)
    {


        // Cap Velocities to stop the particles having a rave
        float maxVel = 100.0f;

        // Compute Density
        CalcDensity(i);

        // Compute Pressure Force, Viscosity and Surface Tension
        CalcPressure(i);

        // Calculate Leapfrog variables (normal acceleration, just halved to get V + 1/2)

        float leapX = 0.5f*(m_Particles[i].m_acc.m_x * SVar.m_TimeStep / m_Particles[i].m_Density + SVar.m_Gravity.m_x * SVar.m_TimeStep);
        float leapY = 0.5f*(m_Particles[i].m_acc.m_y * SVar.m_TimeStep / m_Particles[i].m_Density + SVar.m_Gravity.m_y * SVar.m_TimeStep);
        float leapZ = 0.5f*(m_Particles[i].m_acc.m_z * SVar.m_TimeStep / m_Particles[i].m_Density + SVar.m_Gravity.m_z * SVar.m_TimeStep);


        // Update Velocities using leapfrog variables
        m_Particles[i].m_vel.m_x += leapX;
        m_Particles[i].m_vel.m_y += leapY;
        m_Particles[i].m_vel.m_z += leapZ;

        // Give some directional force
        //m_Particles[i].m_vel.m_x += 0.1f;

        // Cap velocities

        if(m_Particles[i].m_vel.m_y > maxVel)
            m_Particles[i].m_vel.m_y = maxVel;
        if(m_Particles[i].m_vel.m_y < -maxVel)
            m_Particles[i].m_vel.m_y = -maxVel;
        if(m_Particles[i].m_vel.m_x > maxVel)
            m_Particles[i].m_vel.m_x = maxVel;
        if(m_Particles[i].m_vel.m_x < -maxVel)
            m_Particles[i].m_vel.m_x = -maxVel;
        if(m_Particles[i].m_vel.m_z > maxVel)
            m_Particles[i].m_vel.m_z = maxVel;
        if(m_Particles[i].m_vel.m_z < -maxVel)
            m_Particles[i].m_vel.m_z = -maxVel;

        // Update Positions with leapfrog velocities (v 1/2)
        m_Particles[i].m_pos.m_x += m_Particles[i].m_vel.m_x * SVar.m_TimeStep;
        m_Particles[i].m_pos.m_y += m_Particles[i].m_vel.m_y * SVar.m_TimeStep;
        m_Particles[i].m_pos.m_z += m_Particles[i].m_vel.m_z * SVar.m_TimeStep;

        // Update Velocities using leapfrog variables again to bring v to 1
        m_Particles[i].m_vel.m_x += leapX;
        m_Particles[i].m_vel.m_y += leapY;
        m_Particles[i].m_vel.m_z += leapZ;

        // Cap velocities again

        if(m_Particles[i].m_vel.m_y > maxVel)
            m_Particles[i].m_vel.m_y = maxVel;
        if(m_Particles[i].m_vel.m_y < -maxVel)
            m_Particles[i].m_vel.m_y = -maxVel;
        if(m_Particles[i].m_vel.m_x > maxVel)
            m_Particles[i].m_vel.m_x = maxVel;
        if(m_Particles[i].m_vel.m_x < -maxVel)
            m_Particles[i].m_vel.m_x = -maxVel;
        if(m_Particles[i].m_vel.m_z > maxVel)
            m_Particles[i].m_vel.m_z = maxVel;
        if(m_Particles[i].m_vel.m_z < -maxVel)
            m_Particles[i].m_vel.m_z = -maxVel;

        //m_Particles[i].m_vel.m_x += 0.0f;
        //m_Particles[i].m_vel.m_y += -9.8f;
        //m_Particles[i].m_vel.m_z += 0.0f;

        float pushBack = 0.00001f;

        // Keep within bounds, apply dampening
        if(m_Particles[i].m_pos.m_y<=0)
        {
            m_Particles[i].m_vel.m_y = -m_Particles[i].m_vel.m_y * SVar.m_WallDampening;
            m_Particles[i].m_pos.m_y=0;
        }
        if(m_Particles[i].m_pos.m_y>=m_dim.m_y)
        {
            m_Particles[i].m_vel.m_y = -m_Particles[i].m_vel.m_y * SVar.m_WallDampening;
            m_Particles[i].m_pos.m_y=m_dim.m_y-pushBack;
        }
        if(m_Particles[i].m_pos.m_z<=0)
        {
            m_Particles[i].m_vel.m_z = -m_Particles[i].m_vel.m_z * SVar.m_WallDampening;
            m_Particles[i].m_pos.m_z=0+pushBack;
        }
        if(m_Particles[i].m_pos.m_z>=m_dim.m_z)
        {
            m_Particles[i].m_vel.m_z = -m_Particles[i].m_vel.m_z * SVar.m_WallDampening;
            m_Particles[i].m_pos.m_z=m_dim.m_z-pushBack;
        }
        if(m_Particles[i].m_pos.m_x<=0)
        {
            m_Particles[i].m_vel.m_x = -m_Particles[i].m_vel.m_x * SVar.m_WallDampening;
            m_Particles[i].m_pos.m_x=0+pushBack;
        }
        if(m_Particles[i].m_pos.m_x>=m_dim.m_x)
        {
            m_Particles[i].m_vel.m_x = -m_Particles[i].m_vel.m_x * SVar.m_WallDampening;
            m_Particles[i].m_pos.m_x=m_dim.m_x-pushBack;
        }

        // Debug to see what the hell is going on
//        if(i==100)
//            std::cout<<"ID:"<<i<<", Vel:"<<m_Particles[i].m_vel.m_x<<","<<m_Particles[i].m_vel.m_y<<","<<m_Particles[i].m_vel.m_z<<","
//                    << "Density: " <<m_Particles[i].m_Density << " Pressure: " << m_Particles[i].m_Pressure
//                    << "Pressure Force:"  <<m_Particles[i].m_PressureForce.m_x << "," <<m_Particles[i].m_PressureForce.m_y << ","<<m_Particles[i].m_PressureForce.m_z  <<'\n';

        // Apply transforms
        glPtr[glIndex]=m_Particles[i].m_pos.m_x=m_Particles[i].m_pos.m_x;
        glPtr[glIndex+1]=m_Particles[i].m_pos.m_y=m_Particles[i].m_pos.m_y;
        glPtr[glIndex+2]=m_Particles[i].m_pos.m_z=m_Particles[i].m_pos.m_z;


        // Update colour info
        glPtr[glIndex+3]=m_Particles[i].m_col.m_x;
        glPtr[glIndex+4]=m_Particles[i].m_col.m_y;
        glPtr[glIndex+5]=m_Particles[i].m_col.m_z;

      #pragma omp atomic
      glIndex+=6;

    }

    m_vao->unmapBuffer();

    m_vao->unbind();

    //WriteToFile();

//    if(frame < maxFrame)
//    {
//        frame++;
//    }
//    else
//    {
//        exit(0);
//    }
}

///
/// \brief SPHGrid::Draw
/// \param _rot Rotation of the camera
/// Draws the particles
void SPHGrid::Draw(const ngl::Mat4 &_rot)
{
    glPointSize(10.0f);
    ngl::ShaderLib *shader=ngl::ShaderLib::instance();
    shader->use(GetShaderName());
    ngl::Mat4 vp=m_cam->getVPMatrix();

    shader->setUniform("MVP",vp*_rot);

    m_vao->bind();
    m_vao->draw();
    m_vao->unbind();

}

///
/// \brief SPHGrid::FindNeighbours
/// Our SPH neighbourhood search
/// Probably in dire need of optimization
void SPHGrid::FindNeighbours()
{
    ngl::Logger *log = ngl::Logger::instance();
    log->logMessage("Starting particle neighbour allocation.\n");
    QElapsedTimer timer;
    timer.start();



    // Clear out cell maps
    for (int i = 0; i < m_CellIDs.size(); i++)
    {
        m_map[m_CellIDs[i]].clear();
    }

    // Update cell maps to current particles
    for (unsigned int i = 0; i < m_NumParticles; i++)
    {

        int X = (int)(m_Particles[i].m_pos.m_x/m_CellSize);
        int Y = (int)(m_Particles[i].m_pos.m_y/m_CellSize);
        int Z = (int)(m_Particles[i].m_pos.m_z/m_CellSize);

        int CELLID = m_Particles[i].m_CellID = (m_Prime1 * X) + (m_Prime2 * Y) + (m_Prime3 * Z);
        m_map[CELLID].push_back(i);
        m_Particles[i].m_cell = int3(X,Y,Z);

        m_Particles[i].m_col.m_r = m_CellSize*(float)X/m_dim.m_x;
        m_Particles[i].m_col.m_g = m_CellSize*(float)Y/m_dim.m_y;
        m_Particles[i].m_col.m_b = m_CellSize*(float)Z/m_dim.m_z;
    }

    // Perform neighbourhood search
    for(unsigned int i = 0; i < m_NumParticles; i++)
    {

        m_Particles[i].m_Neighbours.clear();


        for(int x = -1; x <= 1; x++)
        {
            for(int y = -1; y <= 1; y++)
            {
                for(int z = -1; z <= 1; z++)
                {
                   int3 CellCheck = m_Particles[i].m_cell;
                   CellCheck.x += x;
                   CellCheck.y += y;
                   CellCheck.z += z;
                   //float check = m_dim.m_x / m_CellSize;
                   // Make sure the current cell is valid to be checked
                   if(        CellCheck.x < 0
                           || CellCheck.y < 0
                           || CellCheck.z < 0
                           || CellCheck.x >= m_dim.m_x/m_CellSize
                           || CellCheck.y >= m_dim.m_y/m_CellSize
                           || CellCheck.z >= m_dim.m_z/m_CellSize)
                   {
                       // Invalid cell
                       break;
                   }
                   else
                   {
                       SPHParticle P = m_Particles[i];

                       // Obtain all the particles in the local cell
                       int CELLID = (CellCheck.x * m_Prime1) + (CellCheck.y * m_Prime2) + (CellCheck.z * m_Prime3);
                       std::vector<int> A = m_map.at(CELLID);

                       for(unsigned int j = 0; j < A.size(); ++j)
                       {
                           SPHParticle N = m_Particles[A[j]];

                            // Check distance against both particles
                           ngl::Vec3 P1 = m_Particles[i].m_pos;
                           ngl::Vec3 P2 = m_Particles[A[j]].m_pos;

                           // Check distance against particle smoothing kernel, to adapt for the squared distance, we square H
                           ngl::Vec3 D = ngl::Vec3(P2.m_x-P1.m_x,P2.m_y-P1.m_y,P2.m_z-P1.m_z);
                           float distance = D.lengthSquared();
                           if(distance < SVar.m_Kernel2*SVar.m_Kernel)
                               m_Particles[i].m_Neighbours.push_back(m_Particles[A[j]].m_id);

                       }
                   }

                }
            }
        }
    }


    log->logMessage("Finished neighbour allocation. Took %d milliseconds\n",timer.elapsed());
}

///
/// \brief SPHGrid::MoveByCellID
/// \param _cellID Input cell ID
/// \param _dir Direction of the movement desired
/// Takes all the particles within a cell by given cellID and moves them
void SPHGrid::MoveByCellID(int _cellID, ngl::Vec3 _dir)
{
    // Spatial Hash Code
    std::vector<int> A = m_map.at(_cellID);
    for(unsigned int i = 0; i < A.size(); ++i)
    {
        m_Particles[A[i]].m_pos += _dir;
    }

}

///
/// \brief SPHGrid::MoveByCell
/// \param _cell Specified cell
/// \param _dir  Direction of the movement desired
/// Takes all the particles within a cell and moves them
void SPHGrid::MoveByCell(int3 _cell, ngl::Vec3 _dir)
{
    // Spatial Hash Code
    int CELLID = (_cell.x * m_Prime1) + (_cell.y * m_Prime2) + (_cell.z * m_Prime3);
    std::vector<int> A = m_map.at(CELLID);
    for(unsigned int i = 0; i < A.size(); ++i)
    {
        m_Particles[A[i]].m_pos += _dir;
    }


}

///
/// \brief SPHGrid::MoveNeighbours
/// \param _ParticleID ID of particle to centre movement around
/// \param _dir Direction of the movement desired
/// Takes all neighbours of a particle and moves them
void SPHGrid::MoveNeighbours(int _ParticleID, ngl::Vec3 _dir)
{
    for(unsigned int i = 0; i < m_Particles[_ParticleID].m_Neighbours.size(); ++i)
    {
        m_Particles[m_Particles[_ParticleID].m_Neighbours[i]].m_pos += _dir;
    }
}


///
/// \brief SPHGrid::CalcDensity
/// \param _ParticleID
/// Calculate the density for a given particle
void SPHGrid::CalcDensity(int _ParticleID)
{
    SPHParticle* P = &m_Particles[_ParticleID];
    SPHParticle* N;
    P->m_Density = 0.0f;
    P->m_Pressure = 0.0f;

    ngl::Vec3 relativePosition;
    float relativeDistance2;

    for(unsigned int i = 0; i < P->m_Neighbours.size(); ++i)
    //for(unsigned int i = 0; i < m_NumParticles; ++i)
    {
        N = &m_Particles[m_Particles[_ParticleID].m_Neighbours[i]];
        if(N->m_id == _ParticleID)
            continue;
        //N = &m_Particles[i];
        relativePosition = ngl::Vec3(N->m_pos.m_x-P->m_pos.m_x, N->m_pos.m_y-P->m_pos.m_y, N->m_pos.m_z-P->m_pos.m_z);
        relativeDistance2 = (relativePosition.m_x*relativePosition.m_x) + (relativePosition.m_y*relativePosition.m_y) + (relativePosition.m_z*relativePosition.m_z);

        if(relativeDistance2 < ZVAL || relativeDistance2 > SVar.m_Kernel2)
        {
            P->m_Density = P->m_Density + P->m_Mass * SVar.m_Poly6Value * pow(SVar.m_Kernel2-relativeDistance2,3);
        }
    }

    P->m_Density = P->m_Density - SVar.m_RDens;
    P->m_Pressure = (pow(P->m_Density/SVar.m_RestDensity,7)-1)*SVar.m_GasConstant;


}

///
/// \brief SPHGrid::CalcPressure
/// \param _ParticleID
/// Calculate the pressure force, viscosity and surface tension for each particle of a given id
void SPHGrid::CalcPressure(int _ParticleID)
{

    SPHParticle* P = &m_Particles[_ParticleID];
    SPHParticle* N;

    ngl::Vec3 relativePosition;
    float relativeDistance, relativeDistance2;
    float vis, kernelRD, pressureKernel, viscKernel, temp;
    float tension = 0.0f;
    ngl::Vec3 gradTen = ngl::Vec3(0.0f,0.0f,0.0f);

    P->m_acc = ngl::Vec3(0.0f,0.0f,0.0f);

    for(unsigned int i = 0; i < P->m_Neighbours.size(); ++i)
    //for(unsigned int i = 0; i < m_NumParticles; ++i)
    {
        N = &m_Particles[m_Particles[_ParticleID].m_Neighbours[i]];
        //N = &m_Particles[i];
        relativePosition = ngl::Vec3(N->m_pos.m_x-P->m_pos.m_x, N->m_pos.m_y-P->m_pos.m_y, N->m_pos.m_z-P->m_pos.m_z);
        relativeDistance2 = (relativePosition.m_x*relativePosition.m_x) + (relativePosition.m_y*relativePosition.m_y) + (relativePosition.m_z*relativePosition.m_z);

        if(relativeDistance2 > ZVAL && relativeDistance2 < SVar.m_Kernel2)
        {
            // Pressure Gradient
            relativeDistance = sqrt(relativeDistance2);
            vis = P->m_Mass / N->m_Density / 2;
            kernelRD = SVar.m_Kernel-relativeDistance;
            pressureKernel = SVar.m_SpikyValue * kernelRD * kernelRD;
            temp = vis * (P->m_Pressure) * pressureKernel;

            P->m_acc.m_x -= relativePosition.m_x * temp / relativeDistance;
            P->m_acc.m_y -= relativePosition.m_y * temp / relativeDistance;
            P->m_acc.m_z -= relativePosition.m_z * temp / relativeDistance;

            // Viscosity Kernel

            viscKernel = SVar.m_ViscosityValue * (SVar.m_Kernel-relativeDistance);
            temp = vis * SVar.m_Viscosity * viscKernel;

            P->m_acc.m_x = P->m_acc.m_x * temp;
            P->m_acc.m_y = P->m_acc.m_y * temp;
            P->m_acc.m_z = P->m_acc.m_z * temp;

            // Surface Tension

            float k2r2 = SVar.m_Kernel2-relativeDistance2;
            float colourAdv = relativeDistance2 - (3/4) * k2r2;

            temp = -1 * SVar.m_Poly6Grad * vis * k2r2*k2r2;
            tension += SVar.m_Poly6Lap * vis * k2r2 * colourAdv;
            gradTen.m_x += temp * relativePosition.m_x;
            gradTen.m_y += temp * relativePosition.m_y;
            gradTen.m_z += temp * relativePosition.m_z;

        }
    }

    // Apply surface tension

    tension += SVar.m_TensionLap/P->m_Density;
    P->m_SurfaceTension = sqrt(gradTen.m_x*gradTen.m_x+gradTen.m_y*gradTen.m_y+gradTen.m_z*gradTen.m_z);
    if(P->m_SurfaceTension > SVar.m_SurfTenValue)
    {
        P->m_acc.m_x += SVar.m_SurfCoef * tension * gradTen.m_x / P->m_SurfaceTension;
        P->m_acc.m_y += SVar.m_SurfCoef * tension * gradTen.m_y / P->m_SurfaceTension;
        P->m_acc.m_z += SVar.m_SurfCoef * tension * gradTen.m_z / P->m_SurfaceTension;
    }

}


///@brief Checks if a number is prime
/// Ref: https://stackoverflow.com/questions/30052316/find-next-prime-number-algorithm
bool IsPrime(int number)
{

    if (number == 2 || number == 3)
        return true;

    if (number % 2 == 0 || number % 3 == 0)
        return false;

    int divisor = 6;
    while (divisor * divisor - 2 * divisor + 1 <= number)
    {

        if (number % (divisor - 1) == 0)
            return false;

        if (number % (divisor + 1) == 0)
            return false;

        divisor += 6;

    }

    return true;

}

///@brief Finds the next prime number of a given number
/// This is used to create our hash table size
/// Ref: https://stackoverflow.com/questions/30052316/find-next-prime-number-algorithm
int NextPrime(int a)
{
    while (!IsPrime(++a))
    { }
    return a;
}


void SPHGrid::WriteToFile()
{

    static int s_frame=0;
    char fname[50];
    //std::sprintf(fname,"geo/SPH_MILLION.%03d.geo",s_frame++);
    std::sprintf(fname,"../../../../../../run/media/s5005745/Backups/SPH/SPH_HTHOUSAND.%03d.geo",s_frame++);
    // we will use a stringstream as it may be more efficient
    std::stringstream ss;
    std::ofstream file;
    file.open(fname);
    if (!file.is_open())
    {
      //std::cout << "failed to Open file "<<fname<<'\n';
      std::cerr << "*** error: could not open output file\n" ;
      perror("Output file failed to open because: ");
      exit(EXIT_FAILURE);
    }
    // write header see here http://www.sidefx.com/docs/houdini15.0/io/formats/geo
    ss << "PGEOMETRY V5\n";
    ss << "NPoints " << m_NumParticles << " NPrims 1\n";
    ss << "NPointGroups 0 NPrimGroups 1\n";
    // this is hard coded but could be flexible we have 1 attrib which is Colour
    ss << "NPointAttrib 1  NVertexAttrib 0 NPrimAttrib 2 NAttrib 0\n";
    // now write out our point attrib this case Cd for diffuse colour
    ss <<"PointAttrib \n";
    // default the colour to white
    //ss <<"v 3 float 0 0 0\n";
    ss <<"Cd 3 float 1 1 1\n";

    // now we write out the particle data in the format
    // x y z 1 (attrib so in this case colour)
    for(unsigned int i=0; i<m_NumParticles; ++i)
    {
    ss<<m_Particles[i].m_pos.m_x<<" "<<m_Particles[i].m_pos.m_y<<" "<<m_Particles[i].m_pos.m_z << " 1 ";
    /*ss<<"("<<m_Particles[i].m_vel.m_x<<" "<<m_Particles[i].m_vel.m_y<<" "<< m_Particles[i].m_vel.m_z<<"*/ss<<"("<<m_Particles[i].m_col.m_r<<" "<< m_Particles[i].m_col.m_g<<" "<< m_Particles[i].m_col.m_b<<")\n";
    }

    // now write out the index values
    ss<<"PrimitiveAttrib\n";
    ss<<"generator 1 index 1 location1\n";
    ss<<"dopobject 1 index 1 /obj/AutoDopNetwork:1\n";
    ss<<"Part "<<m_NumParticles<<" ";
    for(size_t i=0; i<m_NumParticles; ++i)
    {
    ss<<i<<" ";
    }
    ss<<" [0	0]\n";
    ss<<"box_object1 unordered\n";
    ss<<"1 1\n";
    ss<<"beginExtra\n";
    ss<<"endExtra\n";
    // dump string stream to disk;
    file<<ss.rdbuf();
    file.close();


}
