#include <SPHParticle.h>
#include <math.h>


SPHParticle::SPHParticle()
{
    // Set base density to the density of water at rest
    m_Density = 998.2071f;
}

///
/// \brief SPHParticle::RandPoint
/// \param a First bounds
/// \param b Second bounds
/// \return Pseudorandom point between these two bounds
/// Returns a pseudorandom point between two give bounds _a and _b
float SPHParticle::RandPoint(float _a, float _b)
{
    float random = ((float) rand()) / (float) RAND_MAX;
    float diff = _b - _a;
    float r = random * diff;
    return _a + r;
}

///
/// \brief SPHParticle::RandPosSphere
/// \param _r
/// \param _centre
/// \return
/// Generates a random point within a given sphere
ngl::Vec3 SPHParticle::RandPosSphere(float _r, float _centre)
{

    ngl::Vec3 retVal = ngl::Vec3(0,0,0);
    // Attempt at most 10 times
    for (int i = 0; i < 10; i ++)
    {
        float random = ((float) rand()) / (float) RAND_MAX;
        float diff = _r + _r;
        float r = random * diff;
        float x = r - _r;

        random = ((float) rand()) / (float) RAND_MAX;
        r = random * diff;
        float y = r - _r;

        random = ((float) rand()) / (float) RAND_MAX;
        r = random * diff;
        float z = r - _r;

        retVal = ngl::Vec3(x,y,z);

        if(x*x + y*y + z*z < _r*_r)
            return ngl::Vec3(retVal.m_x + _centre,retVal.m_y + _centre,retVal.m_z + _centre);

    }

    return ngl::Vec3(retVal.m_x + _centre,retVal.m_y + _centre,retVal.m_z + _centre);
}

void SPHParticle::UpdatePressure()
{
    // Gas Constant multiplied by density minus rest density
    // Ref: https://github.com/PaulKennedyDIT/SPH/blob/master/Assets/Code/SPH/FluidParticle.cs
    m_Pressure = 3.0f * (((m_Density - 998.2071f)));

}
