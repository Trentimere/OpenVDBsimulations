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


struct Scale {
    Scale(const float s): dt(s) {}
    inline void operator()(const openvdb::Vec3fGrid::ValueOnIter& iter) const {
        iter.setValue(-(*iter * dt));
    }
  float dt;
};

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
  static inline void op(const openvdb::BoolGrid::ValueOnCIter& iter,
						openvdb::Vec3SGrid::Accessor &accessor)
  {
    accessor.setValue(iter.getCoord(), openvdb::Vec3s(100,100,0));
  }
  
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

double PDF(double d)
{
  double n = 0.1;
  double d_avg = 10;
  double y = (n/d_avg)*pow(log(d)/d_avg,n-1)*pow(d,-n/d_avg);
  return y;
}

double getParticleSize()
{
  double d_min = 2;
  double d_max = 500;
  double y_max = PDF(d_min);

  double d_guess;
  double y_guess;
  bool valid = false;
  while (!valid)
  {
    d_guess = double(rand() % int(d_max - d_min + 1)) + d_min;
    y_guess = (double(rand()) / RAND_MAX)*y_max;
    if (PDF(d_guess) > y_guess)
      return d_guess * 0.000001;
  }
  return d_min;
}

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

float dragCoeff(float density_air, float diameter_p, openvdb::Vec3f velocity_air,
		 openvdb::Vec3f velocity_p, float viscocity_air)
{
  float Rep =
    density_air*diameter_p*(velocity_p - velocity_air).length() / viscocity_air;
  float Cd;
  float a1, a2, a3;
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
  Cd = a1 + a2/Rep + a3/(Rep*Rep);
  return Cd;
}

struct ParticleDeformer
{
  ParticleDeformer(const float& _dt,
		   openvdb::points::AttributeHandle<openvdb::Vec3f> _vHandle):
    dt(_dt), vHandle(_vHandle) {}

  template <typename LeafT>
  void reset(const LeafT& leafIter, size_t)
  {
    const openvdb::points::AttributeArray& varray =
      leafIter.constAttributeArray("velocity");
    openvdb::points::AttributeHandle<openvdb::Vec3f> temp(varray);
    vHandle = temp;
  }

  template <typename IndexIterT>
  void apply(openvdb::Vec3d& position, const IndexIterT& iter) const
  {
    position = position + vHandle.get(*iter)*dt;
  }
  float dt;
  openvdb::points::AttributeHandle<openvdb::Vec3f> vHandle;
};

// Make changes to particle velocities if they encounter an object
void handleObjectCollisions(
			    openvdb::points::PointDataGrid::Ptr pointGrid,
			    openvdb::tools::GridSampler<openvdb::Vec3fGrid,
			     openvdb::tools::BoxSampler> surfaceSampler,
				openvdb::Vec3f objectVelocity)
{
  for (auto leafIter = pointGrid->tree().beginLeaf(); leafIter; ++leafIter)
  {
    // create an accessor for particle velocity
      openvdb::points::AttributeArray& varray =
	leafIter->attributeArray("velocity");
      openvdb::points::AttributeWriteHandle<openvdb::Vec3f,
	openvdb::points::NullCodec> vHandle(varray);

      //create an accessor for particle position
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
	  openvdb::Vec3f grad = surfaceSampler.wsSample(position);

	  // If the particle is encountering the object surface
	  if (grad.length() != 0)
	  {
	    // Update particle velocity to have no component into the surface
	    // This needs to be reevaluated as it can cause crashes
	    // if the inital conditions are not correct leading to infinite velocity
	    openvdb::Vec3f n = grad.unit();
	    openvdb::Vec3f vp = vHandle.get(*iter);
	    openvdb::Vec3f vn = objectVelocity.unit();
	    vp = vp + n*(vp.dot(n));
	    openvdb::Vec3f vperp = (n - vn*(vn.dot(n)));
	    if (vp.dot(vn) < objectVelocity.length())
	      vp = vp + objectVelocity;
//	    if ((objectVelocity.length())*vperp.length() > vp.dot(vperp));
	      vp = vp + vperp*objectVelocity.length();
	    vHandle.set(*iter, vp);
	    
	  }
	  
	}
  }
}

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
  std::vector<float> maxValues(2, 0);
  float maxVel = 0.0;
  float maxX = 0.0;
  openvdb::Vec3f gravity(0, 0, -9.81);
  openvdb::Vec3f vel;
  for (auto leafIter = pointGrid->tree().beginLeaf(); leafIter; ++leafIter)
    {
      const openvdb::points::AttributeArray& marray =
	leafIter->constAttributeArray("mass");
      openvdb::points::AttributeHandle<float> mHandle(marray);
      const openvdb::points::AttributeArray& darray =
	leafIter->constAttributeArray("diameter");
      openvdb::points::AttributeHandle<float> dHandle(darray);
      const openvdb::points::AttributeArray& parray =
	leafIter->constAttributeArray("density");
      openvdb::points::AttributeHandle<float> pHandle(parray);
      openvdb::points::AttributeArray& varray =
	leafIter->attributeArray("velocity");
      openvdb::points::AttributeWriteHandle<openvdb::Vec3f,
	openvdb::points::NullCodec> vHandle(varray);
	
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
	  float Cd = dragCoeff(p, dp, v, vp, u);
	  openvdb::Vec3f dragForce = 3*p*mp*Cd*(v-vp)*((v-vp).length()) / (4*pp*dp);
	  if ((dragForce*dt/mp).length() > (v-vp).length())
	    dragForce = (v-vp)*mp/dt;
	  openvdb::Vec3f bouyancy = (pp-p) * M_PI * dp*dp*dp * gravity / 6;

	  openvdb::Vec3f totalForce = (dragForce + bouyancy);
	  openvdb::Vec3f momentum = totalForce*dt;
	  openvdb::Vec3f dv = momentum/mp;
	  openvdb::Coord xyzCoord = iter.getCoord();
		mAccessor.setValue(xyzCoord, mAccessor.getValue(xyzCoord)+momentum);

	  if (position(2) <= 0)
	    {
	      //positionHandle.set(*iter, openvdb::Vec3f(position(0),position(1),0));
	      vel =  openvdb::Vec3f(0);
	    }
	  else
	    {
	      vel = vp+dv;
	    }
	  vHandle.set(*iter, vel);
	  
	  if (vel.length() > maxVel)
	    {
	    maxVel = vel.length();
	    }
	  if (position(0) > maxX)
	    maxX = position(0);
	}
    }
  maxValues[0] = maxVel;
  maxValues[1] = maxX;
  return maxValues;
}

int main()
{
    openvdb::initialize();
    float voxelSize = 0.001, halfWidth = 4.0f;

	float size = 0.2;
	float vSize = 0.01;
	// create a level-set cube to interact with particles
    openvdb::FloatGrid::Ptr cubeGrid =
        openvdb::tools::createLevelSetPlatonic<openvdb::FloatGrid>(6,
            /*radius=*/size, /*center=*/openvdb::Vec3f(0, 0, 0),
            /*voxel size=*/vSize, /*width=*/halfWidth);
    cubeGrid->setName("Cube");

	std::string colorGridName = "Color";
	
    int numPoints = 500;
    float cube_size = 0.01;
    std::vector<openvdb::Vec3R> positions(numPoints, openvdb::Vec3R(0, 0, 0));
    std::vector<float> pDensity(numPoints, 0.0);
    std::vector<float> pDiameter(numPoints, 0.0);
    std::vector<float> pMass(numPoints, 0.0);
    std::vector<openvdb::Vec3f> pVelocity(numPoints, openvdb::Vec3f(0));
// Randomize the point positions.
    std::mt19937 generator(/*seed=*/0);
    std::uniform_real_distribution<> distribution(-0.5, 0.5);
    double waterDensity = 997.0; 

    double maxD = 500*0.000001;
    
    double i, j, k, d;
    float maxVel;
    float vel;
    float h = 1.75;
    double theta1;
    double theta2;
    for (int n = 0; n < numPoints; n++)
    {
      double time = double(rand()) / RAND_MAX;
      i = 0; // distribution(generator)*cube_size;
      j = distribution(generator)*cube_size;
      k = distribution(generator)*cube_size + h;
      positions[n] = openvdb::Vec3R(i,j,k);
      pDensity[n] = waterDensity;
      d = getParticleSize();
      pDiameter[n] = d;
      pMass[n] = waterDensity * M_PI * d*d*d / 6;
      vel = getVelocity(time);
      if (abs(vel) > maxVel)
	maxVel = abs(vel);
      theta1 = 0.5*distribution(generator);
      theta2 = 0.5*distribution(generator);
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

    std::string vdbFile = "grids.vdb";
    std::string ppmFile = "out.ppm";

    int width = 1240;
    int height = 720;

    float maxX = 0.5*cube_size;
    bool init;

    float blue;
    float red;    
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

    openvdb::Vec3fGrid::Ptr velocity =
      openvdb::Vec3fGrid::create(openvdb::Vec3f(0));
	float airDensity = 1.225;
    openvdb::FloatGrid::Ptr density = openvdb::FloatGrid::create(airDensity);
    openvdb::FloatGrid::Ptr viscocity = openvdb::FloatGrid::create(0.0000181);

    
    openvdb::tools::GridSampler<openvdb::FloatGrid,
				openvdb::tools::BoxSampler> densitySampler(*density);
    openvdb::tools::GridSampler<openvdb::Vec3fGrid,
				openvdb::tools::BoxSampler> velocitySampler(*velocity);
    openvdb::tools::GridSampler<openvdb::FloatGrid,
				openvdb::tools::BoxSampler> viscocitySampler(*viscocity);
    
    openvdb::FloatGrid::Ptr divGrid;
    openvdb::FloatGrid::Ptr source = openvdb::FloatGrid::create();;
    openvdb::FloatGrid::Ptr pressure = openvdb::FloatGrid::create();;

    openvdb::Vec3fGrid::Ptr pressureGradient;
    openvdb::FloatGrid::Ptr magnitude;
	
    openvdb::util::NullInterrupter interrupter;

    std::vector<float> maxValues;

    int ct = 7;
    float compTimes[ct] = {0.1, 0.5, 1, 5, 30, 60, 300};
    int compCount = 0;

    float theta = M_PI_4; // initial orientation of cube object

    // compute the gradient of the surface of the cube
    openvdb::FloatGrid::Ptr objectGrid = cubeGrid->deepCopy();
    //openvdb::tools::erodeVoxels(objectGrid->tree(), int(halfWidth/2)+1);
    openvdb::Vec3fGrid::Ptr surfaceGradient =
      openvdb::tools::gradient(*objectGrid);

    // make a copy of the cube for rendering purposes
    openvdb::FloatGrid::Ptr cubeCopy = cubeGrid->deepCopy();
    openvdb::math::Transform &cubeXForm = cubeCopy->transform();
	openvdb::math::Transform &gradXForm = surfaceGradient->transform();

	// rotate the cube about the y-axis by the specified amount
    const openvdb::math::Axis yaxis = openvdb::math::Y_AXIS;

    cubeXForm.preRotate(theta, yaxis);
    gradXForm.preRotate(theta, yaxis);
	openvdb::math::Mat3s rotation =
	  openvdb::math::rotation<openvdb::math::Mat3s>(yaxis, theta);
	openvdb::tools::foreach(surfaceGradient->beginValueOn(), MatMul(rotation));
	
	float speed = 2.0; // speed of cube (m/s)
    openvdb::Vec3s objectVelocity = openvdb::Vec3s(-speed,0,0);

    // define initial position of the cube
	openvdb::Vec3s objectOffset = openvdb::Vec3s(1.5, 0, 1.5);
	cubeXForm.postTranslate(objectOffset);
	gradXForm.postTranslate(objectOffset);
	
	
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

	// Update cube and gradient postiion based on cube speed and time step
	cubeXForm.postTranslate(dt*timeSteps*objectVelocity);
	gradXForm.postTranslate(dt*timeSteps*objectVelocity);
	// create a sampler of the gradient for the cube for boundary conditions
    openvdb::tools::GridSampler<openvdb::Vec3fGrid,
				openvdb::tools::BoxSampler> surfaceSampler(*surfaceGradient);

	auto leafIter = pointGrid->tree().cbeginLeaf();
	const openvdb::points::AttributeArray& varray =
	  leafIter->constAttributeArray("velocity");
	openvdb::points::AttributeHandle<openvdb::Vec3f> vHandle(varray);
	ParticleDeformer deformer(dt, vHandle);
	openvdb::points::movePoints(*pointGrid, deformer);

	openvdb::Vec3fGrid::Ptr momentum = openvdb::Vec3fGrid::create(openvdb::Vec3f(0));
	openvdb::Vec3fGrid::Accessor mAccessor = momentum->getAccessor();
	maxValues = updateParticleVelocity(pointGrid, densitySampler,
									   velocitySampler, viscocitySampler, dt,
									   mAccessor);

	// Change particle velocities based on cube posiiton and speed
	handleObjectCollisions(pointGrid, surfaceSampler, objectVelocity);
	/*
	velocity->tree().combine(momentum->tree(), mExchange(voxelSize*voxelSize*voxelSize, airDensity));
	*/
	maxVel = maxValues[0];
	maxX = maxValues[1];

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
	if (frame % 100 == 0 || t > compTimes[compCount])
	{
	  pointGrid->tree().prune();
	  velocity->tree().prune();
	  MyParticleList pList;
	  int scale = 2;
	  float radius = vSize*scale ;
	  openvdb::CoordBBox bbox;
	  int s = int(radius/vSize) + 1;
	  openvdb::Coord offset = openvdb::Coord(s,s,s);
	  init = false;
	  openvdb::FloatGrid::Ptr renderGrid = 
	    openvdb::createLevelSet<openvdb::FloatGrid>(vSize, 3);
	  openvdb::Vec3SGrid::Ptr colorGrid =
	    openvdb::Vec3SGrid::create(openvdb::Vec3s(255,255,255));
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

		  if (!init)
		    {
		      maxX = worldPosition[0];
		      init = true;
		    }
		  else
		    {
		      if (worldPosition[0] > maxX) maxX = worldPosition[0];
		    }
		  float diameter = dHandle.get(*indexIter);
		  blue = (1 - diameter/maxD);
		  red = (diameter/maxD);
		  openvdb::Coord ijk(openvdb::Vec3I(
						    colorGrid->transform().worldToIndex(worldPosition)));
		  bbox = openvdb::CoordBBox(ijk-offset, ijk+offset);
		  acc.getTree()->fill(bbox, openvdb::Vec3s(red, 0, blue));
		}
	    }
	  
  
	  raster.setGrainSize(1);//a value of zero disables threading
	  raster.rasterizeSpheres(pList);
	  raster.finalize();

	  // temporary grid to render cube into the particle level-set space
	  openvdb::FloatGrid::Ptr temp =
		openvdb::createLevelSet<openvdb::FloatGrid>(vSize, halfWidth);

	  // determine transformation from cube coordinates to partilce world cooordinates
	  const openvdb::math::Transform &targetXform = temp->transform();
	  // Compute a source grid to target grid transform.
	  // (For this example, we assume that both grids' transforms are linear,
	  // so that they can be represented as 4 x 4 matrices.)
	  openvdb::Mat4R xform =
	    cubeXForm.baseMap()->getAffineMap()->getMat4() *
	    targetXform.baseMap()->getAffineMap()->getMat4().inverse();
	  // Create the transformer.
	  openvdb::tools::GridTransformer transformer(xform);
	  // Resample using nearest-neighbor interpolation.
	  transformer.transformGrid<openvdb::tools::PointSampler,
								openvdb::FloatGrid>(*cubeCopy, *temp);
	  temp->tree().prune();

	  // combine grid with level-set particles and transformed cube into 1 grid
	  openvdb::tools::csgUnion(*renderGrid, *temp);
	  
	  // Create a VDB file object and write out the grid.
	  openvdb::io::File("grids.vdb").write({renderGrid, colorGrid});
	  openvdb::io::File("pointGrid.vdb").write({pointGrid});
	  if (t > compTimes[compCount])
	  {
	    std::string pgFileName = "pData" + std::to_string(compCount) + ".vdb";
	    openvdb::io::File(pgFileName).write({pointGrid});
	    compCount++;
	  }
	  float camHeight = (h) / 2;
	  
	  std::string str = "vdb_render ";
			str = str + vdbFile + " " + ppmFile + " -translate "
			  + "1,-5,"
			  + std::to_string(camHeight)
	    + " -lookat "+ "1,0," + std::to_string(camHeight)
	    + " -up 0,0,1"
	    + " -res " + std::to_string(width) + "x" + std::to_string(height)
	    + " -color " + colorGridName
			  //+ " -shader normal"
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
