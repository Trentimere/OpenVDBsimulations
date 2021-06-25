// file to parse an STL file
// IMPORTANT: file msut be binary format (not ASCII)

#include <fstream>
#include <iostream>
#include <sstream>
#include <streambuf>
#include <map>
#include <vector>
#include <iostream>
#include <cassert>

struct Point
{
  float x;
  float y;
  float z;

  Point(float x, float y, float z): x(x), y(y), z(z) {}
  
  bool operator<(const struct Point &p2) const
  {
    return (x < p2.x) || (x == p2.x && y < p2.y) ||
	   (x == p2.x && y == p2.y && z < p2.z);
  }
};

std::ostream& operator<<(std::ostream& out, const struct Point p)
{
  out << "(" << p.x << ", " << p.y << ", " << p.z << ")";
  return out;
}

struct TrianglePoints
{
  int v1;
  int v2;
  int v3;

  TrianglePoints(int i, int j, int k): v1(i), v2(j), v3(k) {}
};

std::ostream& operator<<(std::ostream& out, const struct TrianglePoints tp)
{
  out << "[" << tp.v1 << ", " << tp.v2 << ", " << tp.v3 << "]";
  return out;
}

struct Mesh
{
  int getVertex(const struct Point p)
  {
    if (vertices.count(p) == 0)
    {
      vertices[p] = numPoints;
      sortedVertices.push_back(p);
      numPoints++;
    }
    return vertices[p];
  }

  void addTriangle(const struct Point p1, const struct Point p2,
		   const struct Point p3)
  {
    struct TrianglePoints triangle(getVertex(p1), getVertex(p2), getVertex(p3));
    triangles.push_back(triangle);
  }

  std::vector<struct TrianglePoints> triangles;
  std::map<struct Point, int> vertices;
  std::vector<struct Point> sortedVertices;
  int numPoints = 0;
};

float parse_float(std::ifstream& s)
{
  char buf[sizeof(float)];
  s.read(buf, 4);
  float* fptr = (float*) buf;
  return *fptr;
}

Point parse_point(std::ifstream& s) {
  float x = parse_float(s);
  float y = parse_float(s);
  float z = parse_float(s);
  return Point(x, y, z);
}

struct Mesh parse_stl(const std::string& stl_path) {
  std::ifstream stl_file(stl_path.c_str(), std::ios::in | std::ios::binary);
  if (!stl_file) {
    std::cout << "ERROR: COULD NOT READ FILE" << std::endl;
    assert(false);
  }

  struct Mesh object;
  
  char header[80] = "";
  char numTriangles[4];
  stl_file.read(header, 80);
  stl_file.read(numTriangles, 4);
  unsigned int* r = (unsigned int*) numTriangles;
  unsigned int n = *r;
  for (unsigned int i = 0; i < n; i++) {
    parse_point(stl_file); // read in the normal vector (not used)
    auto p1 = parse_point(stl_file);
    auto p2 = parse_point(stl_file);
    auto p3 = parse_point(stl_file);
    object.addTriangle(p1, p2, p3);
    char dummy[2];
    stl_file.read(dummy, 2);
  }
  return object;
}
