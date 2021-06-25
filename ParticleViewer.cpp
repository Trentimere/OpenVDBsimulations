#include <iostream>
#include <openvdb/openvdb.h>
#include <openvdb/points/PointConversion.h>
#include <openvdb/tools/ParticlesToLevelSet.h>

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


int main(int argc, char* argv[]) {
    // must have a filename argument
    if (argc < 5) {
        std::cout << "Usage: pointFile.vdb out.vdb grid_name particle_size\n";
        return 0;
    }

    // get path to file
    std::string vdbin = "";
    vdbin.append(argv[1]);
    std::string vdbout = "";
    vdbout.append(argv[2]);
    std::string gridName = "";
    gridName.append(argv[3]);
    float p_size = atof(argv[4]);
    std::string radiusLabel = "";
    bool customRadii = false;
    if (argc == 6)
    {
      customRadii = true;
      radiusLabel.append(argv[5]);
    }
    
	 
    // Initialize grid types and point attributes types.
    openvdb::initialize();
    openvdb::io::File newFile(vdbin);
    // Open the file. This reads the file header, but not any grids.
    newFile.open();
    // Read the grid by name.
    openvdb::GridBase::Ptr baseGrid = newFile.readGrid(gridName);
    newFile.close();
    
    // From the example above, "Points" is known to be a PointDataGrid,
    // so cast the generic grid pointer to a PointDataGrid pointer.
    openvdb::points::PointDataGrid::Ptr grid =
      openvdb::gridPtrCast<openvdb::points::PointDataGrid>(baseGrid);
    // Iterate over all the leaf nodes in the grid.
    MyParticleList pList;
    float minRadius = 0;
    float radius;
    std::cout << "Reading points ..." << std::endl;
    for (auto leafIter = grid->tree().cbeginLeaf(); leafIter; ++leafIter)
    {
        // Extract the position attribute from the leaf by name (P is position).
        const openvdb::points::AttributeArray& positionArray =
            leafIter->constAttributeArray("P");
        // Create read-only handles for position and radius.
        openvdb::points::AttributeHandle<openvdb::Vec3f> positionHandle(positionArray);
        // Extract the radius attribute from the leaf by name (pscale is radius).
       
       
        // Iterate over the point indices in the leaf.
        for (auto indexIter = leafIter->beginIndexOn(); indexIter; ++indexIter)
	{
            // Extract the voxel-space position of the point.
            openvdb::Vec3f voxelPosition = positionHandle.get(*indexIter);
            // Extract the world-space position of the voxel.
            openvdb::Vec3d xyz = indexIter.getCoord().asVec3d();
            // Compute the world-space position of the point.
            openvdb::Vec3f worldPosition =
                grid->transform().indexToWorld(voxelPosition + xyz);
            // Extract the radius of the point.
	    if (customRadii)
	    {
	      const openvdb::points::AttributeArray& radiusArray =
		leafIter->constAttributeArray(radiusLabel);
	      openvdb::points::AttributeHandle<float> radiusHandle(radiusArray);
	      radius = p_size*radiusHandle.get(*indexIter);
	    }
	    else
	      radius = p_size;
	    pList.add(worldPosition, openvdb::Real(radius));
	    if (minRadius == 0 || radius < minRadius)
	      minRadius = radius;
        }
    }

    const float voxelSize = minRadius/2, halfWidth = 2.0f;
    std::cout << "#Points: " << pList.size() << std::endl;
    std::cout << "Voxel Size: " << voxelSize << std::endl; 
    openvdb::FloatGrid::Ptr pGrid =
      openvdb::createLevelSet<openvdb::FloatGrid>(voxelSize, halfWidth);
    openvdb::tools::ParticlesToLevelSet<openvdb::FloatGrid> raster(*pGrid);

    raster.setGrainSize(1);//a value of zero disables threading
    raster.rasterizeSpheres(pList);
    raster.finalize();
    pGrid->setName("Points");
    openvdb::io::File(vdbout).write({pGrid});
}
