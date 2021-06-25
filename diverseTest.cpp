#include <openvdb/openvdb.h>
#include <openvdb/tools/LevelSetSphere.h>
#include <openvdb/tools/LevelSetPlatonic.h>
#include <vector>
#include <openvdb/points/PointConversion.h>
#include <openvdb/tools/ParticlesToLevelSet.h>
#include <openvdb/points/PointAdvect.h>
#include <openvdb/points/PointMove.h>
#include <openvdb/tools/PotentialFlow.h>
#include <openvdb/tools/ValueTransformer.h>
#include <openvdb/math/Math.h>
#include <openvdb/tools/PoissonSolver.h>
#include <openvdb/tools/GridOperators.h>
#include <openvdb/Grid.h>
#include <openvdb/Types.h>
#include <openvdb/tools/Interpolation.h>
#include <cmath>
#include <math.h>
#include <stdlib.h>
#include "ppm.h"

// Used to scale every entry in a grid by the time step
struct Scale {
    Scale(const float s): dt(s) {}
    inline void operator()(const openvdb::Vec3fGrid::ValueOnIter& iter) const {
        iter.setValue(-(*iter * dt));
    }
  float dt;
};

// Used to multiply 2 grids together and divide by the time step
// result is stored in another grid
struct Divided
{
  Divided(float div): dt(div) {}
  inline void operator()(const float& a, const float& b,
			 float& result) const
  {
    result = a*b/dt;
  }
  float dt;
};

// Performs a 
struct mExchange
{
  mExchange(float v, float p): volume(v), density(p) {}
  inline void operator()(const openvdb::Vec3f& a, const openvdb::Vec3f& b,
						 openvdb::Vec3f& result) const
  {
	float mass = density*volume;
	result = a + (b/mass);
  }
  float volume;
  float density;
};

struct Scaled
{
  Scaled(float s): scale(s) {}
  inline void operator()(const openvdb::Vec3f& a, const float& b,
			 openvdb::Vec3f& result) const
  {
    result = a*b*scale;
  }
  inline void operator()(const float& a, const float& b,
			 float& result) const
  {
    result = a*b*scale;
  }
  float scale;
};



struct Corrector
{
  static inline void op(const openvdb::Vec3f& a, const float& b,
		 openvdb::Vec3f& result) 
  {
    result = -a/b;
  }
};

struct Multiply
{
  static inline void op(const openvdb::Vec3f& a, const float& b,
		 openvdb::Vec3f& result)
  {
    result = a*b;
  }
};


struct Local {
    static inline void diff(const openvdb::Vec3f& a, const openvdb::Vec3f& b,
      openvdb::Vec3f& result) {
        result = a - b;
    }
    static inline void getX(const openvdb::Vec3fGrid::ValueOnCIter& iter,
      openvdb::tree::ValueAccessor<openvdb::FloatTree>& accessor)
    {     
        if (iter.isVoxelValue()) { // set a single voxel
            accessor.setValue(iter.getCoord(), (*iter)[0]);
        } else { // fill an entire tile
            openvdb::CoordBBox bbox;
            iter.getBoundingBox(bbox);
            accessor.getTree()->fill(bbox, (*iter)[0]); 
        }
    }
    static inline void getY(const openvdb::Vec3fGrid::ValueOnCIter& iter,
      openvdb::tree::ValueAccessor<openvdb::FloatTree>& accessor)
    {     
        if (iter.isVoxelValue()) { // set a single voxel
            accessor.setValue(iter.getCoord(), (*iter)[1]);
        } else { // fill an entire tile
            openvdb::CoordBBox bbox;
            iter.getBoundingBox(bbox);
            accessor.getTree()->fill(bbox, (*iter)[1]); 
        }
    }
    static inline void getZ(const openvdb::Vec3fGrid::ValueOnCIter& iter,
      openvdb::tree::ValueAccessor<openvdb::FloatTree>& accessor)
    {     
        if (iter.isVoxelValue()) { // set a single voxel
            accessor.setValue(iter.getCoord(), (*iter)[2]);
        } else { // fill an entire tile
            openvdb::CoordBBox bbox;
            iter.getBoundingBox(bbox);
            accessor.getTree()->fill(bbox, (*iter)[2]); 
        }
    }
    static inline void bounds(const openvdb::Vec3fGrid::ValueAllIter& iter)
    {
      if (iter.getCoord().z() <= 0)
      {
        iter.setValue(openvdb::Vec3f(0));
      }
    }
};  


struct MatMul {
    openvdb::math::Mat3s M;
    MatMul(const openvdb::math::Mat3s& mat): M(mat) {}
    inline void operator()(const openvdb::Vec3SGrid::ValueOnIter& iter) const {
        iter.setValue(M.transform(*iter));
    }
};


openvdb::Vec3fGrid::Ptr laplacian(openvdb::Vec3fGrid::Ptr velocity)
{
  openvdb::Vec3fGrid::Ptr temp1 = openvdb::tools::gradient(
    *(openvdb::tools::divergence(*velocity)));
  openvdb::Vec3fGrid::Ptr temp2 = openvdb::tools::curl(
    *(openvdb::tools::curl(*velocity)));
  openvdb::Vec3fGrid::Ptr result = openvdb::Vec3fGrid::create();
  result->tree().combine2(temp1->tree(), temp2->tree(), Local::diff);
  return result;
}

openvdb::Vec3fGrid::Ptr convection(openvdb::Vec3fGrid::Ptr velocity)
{
  openvdb::FloatGrid::Ptr velX = openvdb::FloatGrid::create();
  openvdb::FloatGrid::Ptr velY = openvdb::FloatGrid::create();
  openvdb::FloatGrid::Ptr velZ = openvdb::FloatGrid::create();
  openvdb::tools::transformValues(velocity->cbeginValueOn(), *velX, Local::getX);
  openvdb::tools::transformValues(velocity->cbeginValueOn(), *velY, Local::getY);
  openvdb::tools::transformValues(velocity->cbeginValueOn(), *velZ, Local::getZ);

  openvdb::Vec3fGrid::Ptr gradX = openvdb::tools::gradient(*velX);
  openvdb::Vec3fGrid::Ptr gradY = openvdb::tools::gradient(*velY);
  openvdb::Vec3fGrid::Ptr gradZ = openvdb::tools::gradient(*velZ);

  gradX->tree().combine2(gradX->tree(), velX->tree(), Multiply::op);
  gradY->tree().combine2(gradY->tree(), velY->tree(), Multiply::op);
  gradZ->tree().combine2(gradZ->tree(), velZ->tree(), Multiply::op);
  
  openvdb::Vec3fGrid::Ptr result = openvdb::Vec3fGrid::create(openvdb::Vec3f(0));
  openvdb::tools::compSum(*result, *gradX);
  openvdb::tools::compSum(*result, *gradY);
  openvdb::tools::compSum(*result, *gradZ);
  return result;
}

// velocity + dt*(viscocity*laplacian(velocity) - velocity.dot(gradient(velocity)))
void predictor(openvdb::Vec3fGrid::Ptr velocity,
  openvdb::FloatGrid::Ptr viscocity, float dt)
{
  openvdb::Vec3fGrid::Ptr laplaceGrid = laplacian(velocity);
  openvdb::Vec3fGrid::Ptr temp;

  laplaceGrid->tree().combine2(laplaceGrid->tree(), viscocity->tree(), Scaled(dt));

  openvdb::Vec3fGrid::Ptr convectionGrid = convection(velocity);
  openvdb::tools::foreach(convectionGrid->beginValueOn(), Scale(dt));
  
  openvdb::tools::compSum(*laplaceGrid, *convectionGrid);

  openvdb::tools::compSum(*velocity, *laplaceGrid);
}

// standard boundary condition for a rectangular container with a floor at 0
// taken from the poisson unit test example
struct BoundaryOp {
    void operator()(const openvdb::Coord& ijk, const openvdb::Coord& neighbor,
        double& source, double& diagonal) const
    {
        if (neighbor.x() == ijk.x() && neighbor.z() == ijk.z()) {
            // Workaround for spurious GCC 4.8 -Wstrict-overflow warning:
            const openvdb::Coord::ValueType dy = (ijk.y() - neighbor.y());
            if (dy > 0) source -= 1.0;
            else diagonal -= 1.0;
        }
    }
};

// Taken from unit test for particle rendering
// Use this class for visualizing particles
class MyParticleList
{
protected:
    struct MyParticle {
        openvdb::Vec3R p, v;
        openvdb::Real  r;
    };
    openvdb::Real           mRadiusScale;
    openvdb::Real           mVelocityScale;
    std::vector<MyParticle> mParticleList;
public:

    typedef openvdb::Vec3R  PosType;

    MyParticleList(openvdb::Real rScale=1, openvdb::Real vScale=1)
        : mRadiusScale(rScale), mVelocityScale(vScale) {}
    void add(const openvdb::Vec3R &p, const openvdb::Real &r,
             const openvdb::Vec3R &v=openvdb::Vec3R(0,0,0))
    {
        MyParticle pa;
        pa.p = p;
        pa.r = r;
        pa.v = v;
        mParticleList.push_back(pa);
    }
    /// @return coordinate bbox in the space of the specified transfrom
    openvdb::CoordBBox getBBox(const openvdb::GridBase& grid) {
        openvdb::CoordBBox bbox;
        openvdb::Coord &min= bbox.min(), &max = bbox.max();
        openvdb::Vec3R pos;
        openvdb::Real rad, invDx = 1/grid.voxelSize()[0];
        for (size_t n=0, e=this->size(); n<e; ++n) {
            this->getPosRad(n, pos, rad);
            const openvdb::Vec3d xyz = grid.worldToIndex(pos);
            const openvdb::Real   r  = rad * invDx;
            for (int i=0; i<3; ++i) {
                min[i] = openvdb::math::Min(min[i], openvdb::math::Floor(xyz[i] - r));
                max[i] = openvdb::math::Max(max[i], openvdb::math::Ceil( xyz[i] + r));
            }
        }
        return bbox;
    }

    /// Return the total number of particles in list.
    ///  Always required!
    size_t size() const { return mParticleList.size(); }

    /// Get the world space position of n'th particle.
    /// Required by ParticledToLevelSet::rasterizeSphere(*this,radius).
    void getPos(size_t n,  openvdb::Vec3R&pos) const { pos = mParticleList[n].p; }


    void getPosRad(size_t n,  openvdb::Vec3R& pos, openvdb::Real& rad) const {
        pos = mParticleList[n].p;
        rad = mRadiusScale*mParticleList[n].r;
    }
    void getPosRadVel(size_t n,  openvdb::Vec3R& pos, openvdb::Real& rad, openvdb::Vec3R& vel) const {
        pos = mParticleList[n].p;
        rad = mRadiusScale*mParticleList[n].r;
        vel = mVelocityScale*mParticleList[n].v;
    }
    // The method below is only required for attribute transfer
    void getAtt(size_t n, openvdb::Index32& att) const { att = openvdb::Index32(n); }
};

// Probability density function for particle size
double PDF(double d)
{
  double n = 0.1;
  double d_avg = 10;
  double y = (n/d_avg)*pow(log(d)/d_avg,n-1)*pow(d,-n/d_avg);
  return y;
}

// generates a random particle size using the probability density function
double getParticleSize()
{
  double d_min = 2; // minimum diameter (microns)
  double d_max = 500; // maximum diameter (microns)
  double y_max = PDF(d_min); // maximum probability of the PDF

  double d_guess;
  double y_guess;
  bool valid = false;
  // keep generating values until within PDF
  while (!valid)
  {
    d_guess = double(rand() % int(d_max - d_min + 1)) + d_min;
    y_guess = (double(rand()) / RAND_MAX)*y_max;
    if (PDF(d_guess) > y_guess)
      return d_guess * 0.000001; // convert to meters
  }
  return d_min;
}

// velocity function for cough based on time
double getVelocity(double t)
{
	double Vel;
	
	if (t < 0.066)

	  Vel=22.06*sin(8.417508418*M_PI*t); 
	
	else if (0.066 <= t && t < 1.0)

		Vel=42.2468512*exp(-10.07685819*t); 
		
	else 
		
		Vel=0.0;	

	return Vel;
}

// Alexander and Morsi equations for drag coefficient
float dragCoeff(float density_air, float diameter_p, openvdb::Vec3f velocity_air,
		 openvdb::Vec3f velocity_p, float viscocity_air)
{
  // compute Reynolds number
  float Rep =
    density_air*diameter_p*(velocity_p - velocity_air).length() / viscocity_air;
  float Cd;
  float a1, a2, a3;
  // coefficients depend on Reynolds number
  if (Rep < 0.1)
    {a1 = 0; a2 = 24; a3 = 0;}
  else if (Rep < 1)
    {a1 = 3.69; a2 = 22.73; a3 = 0.0903;}
  else if (Rep < 10)
    {a1 = 1.222; a2 = 29.1667; a3 = -3.8889;}
  else if (Rep < 100)
    {a1 = 0.6167; a2 = 46.5; a3 = -116.67;}
  else if (Rep < 1000)
    {a1 = 0.3644; a2 = 98.33; a3 = -2778;}
  else if (Rep < 5000)
    {a1 = 0.357; a2 = 148.62; a3 = -47500;}
  else if (Rep < 10000)
    {a1 = 0.46; a2 = -490.546; a3 = 578700;}
  else
    {a1 = 0.5191; a2 = -1662.5; a3 = 5416700;}
  // calculate drag coefficient using coefficients and Re
  Cd = a1 + a2/Rep + a3/(Rep*Rep);
  return Cd;
}

// struct used to move particles using movepoints
struct ParticleDeformer
{
  // update position based on time step and particle velocity
  ParticleDeformer(const float& _dt,
		   openvdb::points::AttributeHandle<openvdb::Vec3f> _vHandle):
    dt(_dt), vHandle(_vHandle) {}

  // initialize the velocity handle for each leaf iteration
  template <typename LeafT>
  void reset(const LeafT& leafIter, size_t)
  {
    const openvdb::points::AttributeArray& varray =
      leafIter.constAttributeArray("velocity");
    openvdb::points::AttributeHandle<openvdb::Vec3f> temp(varray);
    vHandle = temp;
  }

  // perform position update based on each particle velocity and time step
  template <typename IndexIterT>
  void apply(openvdb::Vec3d& position, const IndexIterT& iter) const
  {
    position = position + vHandle.get(*iter)*dt;
  }

  // member variables
  float dt;
  openvdb::points::AttributeHandle<openvdb::Vec3f> vHandle;
};

// Updates the velocity of each particle
// takes into account current position,
// the fluid density, fluid velocity, fluid viscocity,
// time step and particle mass

// returns vector of min and max particle velocity
std::vector<float> updateParticleVelocity(
			    openvdb::points::PointDataGrid::Ptr pointGrid,
			    openvdb::tools::GridSampler<openvdb::FloatGrid,
			    openvdb::tools::BoxSampler> pSampler,
			    openvdb::tools::GridSampler<openvdb::Vec3fGrid,
			     openvdb::tools::BoxSampler> vSampler,
			    openvdb::tools::GridSampler<openvdb::FloatGrid,
			    openvdb::tools::BoxSampler> uSampler,
			    double dt,
				openvdb::Vec3fGrid::Accessor mAccessor)
{
  // keep track of min and max velocity to inform time step choice later
  std::vector<float> maxValues(2, 0);
  float maxVel = 0.0;
  float maxX = 0.0;
  
  openvdb::Vec3f gravity(0, 0, -9.81); // define acceleration of gravity
  openvdb::Vec3f vel;
  // loop through all particle leaf nodes
  for (auto leafIter = pointGrid->tree().beginLeaf(); leafIter; ++leafIter)
    {
      // get an accessor to particle mass
      const openvdb::points::AttributeArray& marray =
	leafIter->constAttributeArray("mass");
      openvdb::points::AttributeHandle<float> mHandle(marray);
      // get an accessor to particle diameter
      const openvdb::points::AttributeArray& darray =
	leafIter->constAttributeArray("diameter");
      openvdb::points::AttributeHandle<float> dHandle(darray);
      // get an accessor to particle density
      const openvdb::points::AttributeArray& parray =
	leafIter->constAttributeArray("density");
      openvdb::points::AttributeHandle<float> pHandle(parray);
      // get and accessor to particle velocity
      openvdb::points::AttributeArray& varray =
	leafIter->attributeArray("velocity");
      openvdb::points::AttributeWriteHandle<openvdb::Vec3f,
	openvdb::points::NullCodec> vHandle(varray);

      // get an accessor to particle position
      openvdb::points::AttributeArray& positionArray =
	leafIter->attributeArray("P");
      // Create read-only handles for position and radius.
      openvdb::points::AttributeWriteHandle<openvdb::Vec3f,
	openvdb::points::NullCodec> positionHandle(positionArray);
      // Iterate over the point indices in the leaf.
      for (auto iter = leafIter->beginIndexOn(); iter; ++iter)
	{
	  
	  // Extract the voxel-space position of the point.
	  openvdb::Vec3f voxelPosition = positionHandle.get(*iter);
	  // Extract the world-space position of the voxel.
	  openvdb::Vec3d xyz = iter.getCoord().asVec3d();
	  // Compute the world-space position of the point.
	  openvdb::Vec3f position =
	    pointGrid->transform().indexToWorld(voxelPosition + xyz);

	  float mp = mHandle.get(*iter);  // particle mass
	  float dp = dHandle.get(*iter);  // particle diameter
	  float pp = pHandle.get(*iter);  // particle density
	  openvdb::Vec3f vp = vHandle.get(*iter);  // particle velocity

	  float p = pSampler.wsSample(position);  // air density
	  openvdb::Vec3f v = vSampler.wsSample(position); // air speed
	  float u = uSampler.wsSample(position);  // air viscocity
	  float Cd = dragCoeff(p, dp, v, vp, u);  // compute drag coefficient
	  // compute drag force using the computed drag coefficient
	  // this is the same equation used by ANSYS
	  openvdb::Vec3f dragForce = 3*p*mp*Cd*(v-vp)*((v-vp).length()) / (4*pp*dp);
	  // change drag force if it would result in particle oscillating
	  // this should not be necessary, needs more investigation
	  if ((dragForce*dt/mp).length() > (v-vp).length())
	    dragForce = (v-vp)*mp/dt;

	  // calculate bouyancy force from fluid density and gravity
	  openvdb::Vec3f bouyancy = (pp-p) * M_PI * dp*dp*dp * gravity / 6;

	  // compute total forces on particle
	  openvdb::Vec3f totalForce = (dragForce + bouyancy);

	  // compute momentum change of particle and update density of fluid
	  openvdb::Vec3f momentum = totalForce*dt;
	  openvdb::Vec3f dv = momentum/mp;
	  openvdb::Coord xyzCoord = iter.getCoord();
		mAccessor.setValue(xyzCoord, mAccessor.getValue(xyzCoord)+momentum);

	  // particle is on the ground, stop moving
	  if (position(2) <= 0)
	    {
	      //positionHandle.set(*iter, openvdb::Vec3f(position(0),position(1),0));
	      vel =  openvdb::Vec3f(0);
	    }
	  // update particle velocity
	  else
	    {
	      vel = vp+dv;
	    }
	  vHandle.set(*iter, vel);

	  // update min/max velocity if applicable
	  if (vel.length() > maxVel)
	    {
	    maxVel = vel.length();
	    }
	  if (position(0) > maxX)
	    maxX = position(0);
	}
    }
  // return max velocity
  maxValues[0] = maxVel;
  maxValues[1] = maxX;
  return maxValues;
}

int main()
{
    openvdb::initialize();
    float voxelSize = 0.001,/*grid resolution*/
      halfWidth = 2.0f; // width in voxels of level-set surface
 
    int numPoints = 500; // number of points to simulate
    float cube_size = 0.01; // size of cube to distribute particles in

    // vectors to store particle parameters
    // will be added as properties to point grid points later
    std::vector<openvdb::Vec3R> positions(numPoints, openvdb::Vec3R(0, 0, 0));
    std::vector<float> pDensity(numPoints, 0.0);
    std::vector<float> pDiameter(numPoints, 0.0);
    std::vector<float> pMass(numPoints, 0.0);
    std::vector<openvdb::Vec3f> pVelocity(numPoints, openvdb::Vec3f(0));
// Randomize the point positions.
    std::mt19937 generator(/*seed=*/0);
    std::uniform_real_distribution<> distribution(-0.5, 0.5);
    
    double waterDensity = 997.0; // kg / m^3

    double maxD = 500*0.000001; // maximum water droplet diameter in meters
    
    double i, j, k, d; // variables for 3D position (ijk) and diameter (d) of particles
    float maxVel; // for maximum velocity of initial particles
    float vel; // variable for particle velocity
    float h = 1.75; // height of center of particle cload
    // angles determining initial direction of particle trajectories
    double theta1;
    double theta2;
    // initalize all the particles and their properties
    for (int n = 0; n < numPoints; n++)
    {
      // get a random time at which particle spawns
      double time = double(rand()) / RAND_MAX;
      // determine 3D position of particle
      i = 0; // distribution(generator)*cube_size;
      j = distribution(generator)*cube_size;
      k = distribution(generator)*cube_size + h;
      // sotre these parameters for later
      positions[n] = openvdb::Vec3R(i,j,k);
      pDensity[n] = waterDensity;
      // get a random diameter using distribution given
      d = getParticleSize(); 
      pDiameter[n] = d;
      // determine mass based on diameter and density
      pMass[n] = waterDensity * M_PI * d*d*d / 6;
      // get particle velocity based on time of spawning
      vel = getVelocity(time);
      if (abs(vel) > maxVel)
	maxVel = abs(vel);
      // randomly modify trajectory within a cone of trajectories
      theta1 = 0.5*distribution(generator);
      theta2 = 0.5*distribution(generator);
      // determine x,y,z velocity and set it
      pVelocity[n] =
	vel*openvdb::Vec3f(pow(cos(theta1),2)-pow(sin(theta2),2),
			   sin(theta1), sin(theta2));
    }

    
    openvdb::points::PointAttributeVector<openvdb::Vec3R> positionsWrapper(positions);
    // Create a transform using this voxel-size.
    openvdb::math::Transform::Ptr transform =
        openvdb::math::Transform::createLinearTransform(voxelSize);
    // Create a PointDataGrid containing these four points and using the
    // transform given. This function has two template parameters, (1) the codec
    // to use for storing the position, (2) the grid we want to create
    // (ie a PointDataGrid).
    // We use no compression here for the positions.
    openvdb::tools::PointIndexGrid::Ptr pointIndexGrid =
      openvdb::tools::createPointIndexGrid<openvdb::tools::PointIndexGrid>(
				       positionsWrapper, *transform);
    openvdb::points::PointDataGrid::Ptr pointGrid =
        openvdb::points::createPointDataGrid<openvdb::points::NullCodec,
	     openvdb::points::PointDataGrid>(*pointIndexGrid,
					     positionsWrapper, *transform);

    
    using Codec = openvdb::points::NullCodec;
    openvdb::points::TypedAttributeArray<float, Codec>::registerType();
    openvdb::NamePair densityAttribute =
        openvdb::points::TypedAttributeArray<float, Codec>::attributeType();
    openvdb::points::appendAttribute(pointGrid->tree(),
				     "density", densityAttribute);
    // Create a wrapper around the density vector.
    openvdb::points::PointAttributeVector<float> densityWrapper(pDensity);
    // Populate the attribute on the points
    openvdb::points::populateAttribute<openvdb::points::PointDataTree,
        openvdb::tools::PointIndexTree,
	openvdb::points::PointAttributeVector<float>>(
        pointGrid->tree(), pointIndexGrid->tree(), "density", densityWrapper);
    
    openvdb::NamePair diameterAttribute =
        openvdb::points::TypedAttributeArray<float, Codec>::attributeType();
    openvdb::points::appendAttribute(pointGrid->tree(),
				     "diameter", diameterAttribute);
    // Create a wrapper around the diameter vector.
    openvdb::points::PointAttributeVector<float> diameterWrapper(pDiameter);
    // Populate the attribute on the points
    openvdb::points::populateAttribute<openvdb::points::PointDataTree,
        openvdb::tools::PointIndexTree,
	openvdb::points::PointAttributeVector<float>>(
        pointGrid->tree(), pointIndexGrid->tree(), "diameter", diameterWrapper);
    
    openvdb::NamePair massAttribute =
        openvdb::points::TypedAttributeArray<float, Codec>::attributeType();
    openvdb::points::appendAttribute(pointGrid->tree(),
				     "mass", massAttribute);
    // Create a wrapper around the mass vector.
    openvdb::points::PointAttributeVector<float> massWrapper(pMass);
    // Populate the attribute on the points
    openvdb::points::populateAttribute<openvdb::points::PointDataTree,
        openvdb::tools::PointIndexTree,
	openvdb::points::PointAttributeVector<float>>(
        pointGrid->tree(), pointIndexGrid->tree(), "mass", massWrapper);
    
    openvdb::NamePair velocityAttribute =
      openvdb::points::TypedAttributeArray<openvdb::Vec3f, Codec>::attributeType();
    openvdb::points::appendAttribute(pointGrid->tree(),
				     "velocity", velocityAttribute);
    openvdb::points::TypedAttributeArray<openvdb::Vec3f, Codec>::registerType();
    // Create a wrapper around the velocity vector.
    openvdb::points::PointAttributeVector<openvdb::Vec3f> velocityWrapper(pVelocity);
    // Populate the attribute on the points
    openvdb::points::populateAttribute<openvdb::points::PointDataTree,
        openvdb::tools::PointIndexTree,
        openvdb::points::PointAttributeVector<openvdb::Vec3f>>(
        pointGrid->tree(), pointIndexGrid->tree(), "velocity", velocityWrapper);
    
    // Set the name of the grid
    pointGrid->setName("Points");

    // file names for output and intermediate names for rendering scenes
    std::string vdbFile = "grids.vdb";
    std::string ppmFile = "out.ppm";

    // size of window to view simulation
    int width = 1240;
    int height = 720;

    // variables used to change camera position
    // ignore any warnings if compiler says they are unused
    float minX = 0;
    float maxX = 0.5*cube_size;
    float minZ;
    float maxZ;
    float avgX;
    float avgZ;
    bool init;

    float blue;
    float red;

    float camDist;
    
    // create an SDL context
    SDL_Init(SDL_INIT_VIDEO);
    SDL_Window* window = SDL_CreateWindow("PPM", 100, 100, width, height, SDL_WINDOW_OPENGL);
    SDL_GLContext context = SDL_GL_CreateContext(window);
    SDL_Event event;
    
    std::cout << "Simulating for 5 seconds" << std::endl;
    
    int timeSteps = 1;
    float t = 0.0;
    float dt;

	float minfVel;
	float maxfVel;
	float ft = 0.0;

	// Initialize fluid velocity as 0
    openvdb::Vec3fGrid::Ptr velocity =
      openvdb::Vec3fGrid::create(openvdb::Vec3f(0));

    // initalize fluid density and viscocity
    float airDensity = 1.225;
    openvdb::FloatGrid::Ptr density = openvdb::FloatGrid::create(airDensity);
    openvdb::FloatGrid::Ptr viscocity = openvdb::FloatGrid::create(0.0000181);

    // create samplers for accessing values of fluid properties at particle points
    openvdb::tools::GridSampler<openvdb::FloatGrid,
				openvdb::tools::BoxSampler> densitySampler(*density);
    openvdb::tools::GridSampler<openvdb::Vec3fGrid,
				openvdb::tools::BoxSampler> velocitySampler(*velocity);
    openvdb::tools::GridSampler<openvdb::FloatGrid,
				openvdb::tools::BoxSampler> viscocitySampler(*viscocity);

    // grids used for navier stokes equations
    openvdb::FloatGrid::Ptr divGrid;
    openvdb::FloatGrid::Ptr source = openvdb::FloatGrid::create();
    openvdb::FloatGrid::Ptr pressure = openvdb::FloatGrid::create();

    openvdb::Vec3fGrid::Ptr pressureGradient;
    openvdb::FloatGrid::Ptr magnitude;
	
    openvdb::util::NullInterrupter interrupter;
	
    std::string colorGridName = "Color";

    std::vector<float> maxValues;

    // set simulation times to explicityly save snapshots
    int ct = 7;
    float compTimes[ct] = {0.1, 0.5, 1, 5, 30, 60, 300};
    int compCount = 0;

    // simulate until need no more pictures
    for (int frame = 0; t < compTimes[ct-1]; frame++) {
        if (SDL_PollEvent(&event)) {
            // break if esc key is pressed or window is closed
            if (event.type == SDL_QUIT) 
                break;
            if (
                event.type == SDL_KEYUP && 
                event.key.keysym.sym == SDLK_ESCAPE
            )
                break;
        }
	// CFL criterion	
	dt = 0.5*voxelSize / maxVel;           ;
	std::cout << "Time: " << t << "\tdt: " << dt << "\tMax Distance: " <<
	maxX << std::endl;
    
	t += dt*timeSteps;

	// create a velocity accessor
	auto leafIter = pointGrid->tree().cbeginLeaf();
	const openvdb::points::AttributeArray& varray =
	  leafIter->constAttributeArray("velocity");
	openvdb::points::AttributeHandle<openvdb::Vec3f> vHandle(varray);
	// move all the points based on their velocity attribute and time step
	ParticleDeformer deformer(dt, vHandle);
	openvdb::points::movePoints(*pointGrid, deformer);

	// update the particle and fluid properties including new particle velocity
	openvdb::Vec3fGrid::Ptr momentum = openvdb::Vec3fGrid::create(openvdb::Vec3f(0));
	openvdb::Vec3fGrid::Accessor mAccessor = momentum->getAccessor();
	maxValues = updateParticleVelocity(pointGrid, densitySampler,
									   velocitySampler, viscocitySampler, dt,
									   mAccessor);
	/*
	velocity->tree().combine(momentum->tree(), mExchange(voxelSize*voxelSize*voxelSize, airDensity));
	*/
	maxVel = maxValues[0];
	maxX = maxValues[1];

	// these are calculations to perfrom Navier Stokes equations on the fluid
	// They should work, but have little effect on the particle behavior
	/*
	if (t >= ft)
	{ 
	  magnitude = openvdb::tools::magnitude(*velocity);
      magnitude->tree().evalMinMax(minfVel, maxfVel);
	  float dft = voxelSize / maxfVel;
	  ft += dft*timeSteps;
	  std::cout << "Time: " << t << "\tdft: " << dft << std::endl;
	  
	  openvdb::tools::foreach(velocity->beginValueAll(), Local::bounds);

	  predictor(velocity, viscocity, dt);
        
	  divGrid = openvdb::tools::divergence(*velocity);
	  source->tree().combine2(density->tree(), divGrid->tree(), Scaled(dt));
	
	  openvdb::math::pcg::State state =
		openvdb::math::pcg::terminationDefaults<float>();
	  state.iterations = 100;
	  const float epsilon = openvdb::math::Delta<float>::value();
	  state.relativeError = state.absoluteError = epsilon;

	  pressure->treePtr() = openvdb::tools::poisson::solveWithBoundaryConditions(
		  source->tree(), BoundaryOp(), state, interrupter, staggered=true);

	  // velocity = -gradient(pressure) / density
	  pressureGradient = openvdb::tools::gradient(*pressure);
	  velocity->tree().combine2(pressureGradient->tree(),
								density->tree(), Corrector::op);
	}
*/

	// only show every 100 time steps to reduce execution time
	// since rendering is a significant time cost
	if (frame % 100 == 0 || t > compTimes[compCount])
	{
	  // clean up particle and velocity grids for faster execution
	  pointGrid->tree().prune();
	  velocity->tree().prune();

	  // object to hold particles in correct format for rendering
	  MyParticleList pList;
	  // define size of particles to render
	  float vSize = 0.01;
	  int scale = 3;
	  float radius = vSize*scale ;
	  openvdb::CoordBBox bbox;
	  int s = int(radius/vSize) + 1;
	  openvdb::Coord offset = openvdb::Coord(s,s,s);
	  init = false;

	  // create a grid to hold level-sets of particles
	  openvdb::FloatGrid::Ptr renderGrid = 
	    openvdb::createLevelSet<openvdb::FloatGrid>(vSize, halfWidth);
	  // create a color grid to store level-set surface color
	  openvdb::Vec3SGrid::Ptr colorGrid =
	    openvdb::Vec3SGrid::create(openvdb::Vec3s(1,1,1));
	  colorGrid->setTransform(renderGrid->transform().copy());
	  colorGrid->setName(colorGridName);
	  openvdb::Vec3SGrid::Accessor acc = colorGrid->getAccessor();
	  
	  openvdb::tools::ParticlesToLevelSet<openvdb::FloatGrid> raster(*renderGrid);
	  for (auto leafIter = pointGrid->tree().cbeginLeaf(); leafIter; ++leafIter)
	    {
	      // Extract the position attribute from the leaf by name (P is position).
	      const openvdb::points::AttributeArray& positionArray =
		leafIter->constAttributeArray("P");
	      // Create read-only handles for position and radius.
		openvdb::points::AttributeHandle<openvdb::Vec3f> positionHandle(positionArray);
	      const openvdb::points::AttributeArray& darray =
		leafIter->constAttributeArray("diameter");
	      // Create read-only handles for position and radius.
		openvdb::points::AttributeHandle<float> dHandle(darray);
	      // Iterate over the point indices in the leaf.
	      for (auto indexIter = leafIter->beginIndexOn(); indexIter; ++indexIter)
		{
		  // Extract the voxel-space position of the point.
		  openvdb::Vec3f voxelPosition = positionHandle.get(*indexIter);
		  // Extract the world-space position of the voxel.
		  openvdb::Vec3d xyz = indexIter.getCoord().asVec3d();
		  // Compute the world-space position of the point.
		  openvdb::Vec3f worldPosition =
		    pointGrid->transform().indexToWorld(voxelPosition + xyz);
		  pList.add(worldPosition, openvdb::Real(radius));

		  // keep track of min/max coordinates particles are
		  if (!init)
		    {
		      minX = worldPosition[0];
		      maxX = worldPosition[0];
		      minZ = worldPosition[2];
		      maxZ = worldPosition[2];
		      init = true;
		    }
		  else
		    {
		      if (worldPosition[0] < minX) minX = worldPosition[0];
		      if (worldPosition[0] > maxX) maxX = worldPosition[0];
		      if (worldPosition[2] < minZ) minZ = worldPosition[2];
		      if (worldPosition[2] > maxZ) maxZ = worldPosition[2];
		    }
		  // define color of particle based on particle diameter
		  float diameter = dHandle.get(*indexIter);
		  blue = (1 - diameter/maxD);
		  red = (diameter/maxD);
		  // set the color grid values in the neighborhood of the particle
		  openvdb::Coord ijk(openvdb::Vec3I(
						    colorGrid->transform().worldToIndex(worldPosition)));
		  bbox = openvdb::CoordBBox(ijk-offset, ijk+offset);
		  acc.getTree()->fill(bbox, openvdb::Vec3s(red, 0, blue));
		}
	    }
	  
  
	  raster.setGrainSize(1);//a value of zero disables threading
	  raster.rasterizeSpheres(pList);
	  raster.finalize();
	  //openvdb::tools::csgUnion(*renderGrid, *temp);
	  // Create a VDB file object and write out the grid.
	  openvdb::io::File("grids.vdb").write({renderGrid, colorGrid});
	  openvdb::io::File("pointGrid.vdb").write({pointGrid});
	  // also save the vdb file if this is an important time step
	  if (t > compTimes[compCount])
	  {
	    std::string pgFileName = "pData" + std::to_string(compCount) + ".vdb";
	    openvdb::io::File(pgFileName).write({pointGrid});
	    compCount++;
	  }
	  avgX = (minX + maxX) / 2;
	  avgZ = (minZ + maxZ) / 2;
	  camDist = (maxX-minX) * 1.75;
	  float camHeight = (h) / 2;

	  // render the level-set particles with the color grid specifications
	  std::string str = "vdb_render ";
			str = str + vdbFile + " " + ppmFile + " -translate "
			  + "1.5,-5,"
			  + std::to_string(camHeight)
	    + " -lookat "+ "1.5,0," + std::to_string(camHeight)
	    + " -up 0,0,1"
	    + " -res " + std::to_string(width) + "x" + std::to_string(height)
	    + " -color " + colorGridName
	    + " -shader normal"
	    + " -frame 60";
	  const char* command = str.c_str();
	  system(command);

	  // open file stream in binary mode
	  FILE* stream = fopen(ppmFile.c_str(), "rb");
  
	  // pixel map containing data from ppm
	  int depth;
	  unsigned int* pixmap = decodePPM(stream, &width, &height, &depth);
    
	  // clean up file stream
	  fclose(stream);

	  // draw calls here
	  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	  glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	  glRasterPos2i(-1.0f, -1.0f);
	  glDrawPixels(width, height, GL_RGBA, GL_UNSIGNED_BYTE, pixmap);
	  SDL_GL_SwapWindow(window);
	  usleep(10000);
	}
	
    }
    

    // clean up SDL context
    SDL_GL_DeleteContext(context);
    SDL_Quit();
}
