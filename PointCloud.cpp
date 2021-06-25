#include <iostream>
#include <openvdb/openvdb.h>
#include <openvdb/tools/LevelSetSphere.h>
#include <openvdb/points/PointCount.h>
#include <vector>
#include <openvdb/points/PointConversion.h>
int main()
{
    // Initialize grid types and point attributes types.
    openvdb::initialize();
    int numPoints = 10000;
    float cube_size = 20.0;
    std::vector<openvdb::Vec3R> positions(numPoints, openvdb::Vec3R(0, 0, 0));
// Randomize the point positions.
    std::mt19937 generator(/*seed=*/0);
    std::uniform_real_distribution<> distribution(-0.5, 0.5);
    double i, j, k;
      for (int n = 0; n < numPoints; n++)
      {
	i = distribution(generator)*cube_size;
	j = distribution(generator)*cube_size;
	k = distribution(generator)*cube_size;
	positions[n] = openvdb::Vec3R(i,j,k);
      }
	  
    
    openvdb::points::PointAttributeVector<openvdb::Vec3R> positionsWrapper(positions);
    // This method computes a voxel-size to match the number of
    // points / voxel requested. Although it won't be exact, it typically offers
    // a good balance of memory against performance.
    int pointsPerVoxel = 8;
    float voxelSize =
        openvdb::points::computeVoxelSize(positionsWrapper, pointsPerVoxel);
    // Print the voxel-size to cout
    std::cout << "VoxelSize=" << voxelSize << std::endl;
    // Create a transform using this voxel-size.
    openvdb::math::Transform::Ptr transform =
        openvdb::math::Transform::createLinearTransform(voxelSize);
    // Create a PointDataGrid containing these four points and using the
    // transform given. This function has two template parameters, (1) the codec
    // to use for storing the position, (2) the grid we want to create
    // (ie a PointDataGrid).
    // We use no compression here for the positions.
    openvdb::points::PointDataGrid::Ptr grid =
        openvdb::points::createPointDataGrid<openvdb::points::NullCodec,
                        openvdb::points::PointDataGrid>(positions, *transform);
    openvdb::points::appendAttribute(grid->tree(), "radius",
      openvdb::points::TypedAttributeArray<float>::attributeType());
    std::uniform_real_distribution<> radDist(0.01, 0.1);
    float size;
    // Iterate over leaf nodes.
    for (auto leafIter = grid->tree().beginLeaf(); leafIter; ++leafIter)
    {
    // Create a read-write samples handle.
      openvdb::points::AttributeArray& array(
        leafIter->attributeArray("radius"));
      openvdb::points::AttributeWriteHandle<float> handle(array);
    // Iterate over the point indices in the leaf.
      for (auto indexIter = leafIter->beginIndexOn(); indexIter; ++indexIter)
      {
	size = radDist(generator);
        handle.set(*indexIter, 0,  size);
      }
    }
    // Set the name of the grid
    grid->setName("Points");
    // Create a VDB file object and write out the grid.
    openvdb::io::File("mypoints.vdb").write({grid});
}
