#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Plane_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Exact_predicates_tag Tag;
struct FaceInfo {
  bool interior, processed;
  FaceInfo() { // useful, use it!
    processed = false;
    interior = false;
  }
};
typedef CGAL::Triangulation_vertex_base_2<Kernel> VertexBase;
typedef CGAL::Constrained_triangulation_face_base_2<Kernel> FaceBase;
typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo, Kernel, FaceBase> FaceBaseWithInfo;
typedef CGAL::Triangulation_data_structure_2<VertexBase, FaceBaseWithInfo> TriangulationDataStructure;
typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel, TriangulationDataStructure, Tag> Triangulation;
typedef Kernel::Point_2 Point2;
typedef Kernel::Point_3 Point3;


//const std::string input_file = "/Users/ken/Downloads/hw1/NL.IMBAG.Pand.0503100000032914-0.obj"; // the faculty building
const std::string input_file = "/mnt/c/Users/LarsB/OneDrive/Documenten/GitHub/GEO1004_ass_1/hw1_obj_files/NL.IMBAG.Pand.0503100000000138-0.obj"; // simple test file
//const std::string output_file = "/Users/ken/Downloads/faculty.obj";
const std::string output_file = "/mnt/c/Users/LarsB/OneDrive/Documenten/GitHub/GEO1004_ass_1/hw1_output_files/test.obj";

// struct are like public classes

struct Vertex {
  double x, y, z;
};

struct Face {
  std::list<int> boundary; // ids of vertices vector
  Kernel::Plane_3 best_plane; // place to store best fitting plane
  Triangulation triangulation; // place to store triangulation
};

int main(int argc, const char * argv[]) {
    std::vector<Vertex> vertices;
    std::vector<Face> faces;

    // Step 1: Import the OBJ file and load it into a data structure
    std::ifstream input_stream;
    input_stream.open(input_file);
    if (input_stream.is_open()) {
        std::string line;

        // Parse line by line
        while (getline(input_stream, line)) {
            //std::cout << line << std::endl;

            std::istringstream line_stream(line);
            char line_type;
            line_stream >> line_type;

            // Vertex
            if (line_type == 'v') {
                vertices.emplace_back();
                double x, y, z;
                line_stream >> x >> y >> z;
                vertices.back().x = x;
                vertices.back().y = y;
                vertices.back().z = z;
            }

            // Face
            if (line_type == 'f') {
                faces.emplace_back();
                int v;
                while (!line_stream.eof()) {
                    line_stream >> v;
                    faces.back().boundary.emplace_back(v-1);
            }
          }

        }
  }

  for (auto vertex : vertices) {
      std::cout << vertex.x << " " << vertex.y << " " << vertex.z << std::endl;
  }

  // Step 2: For each of the faces compute the best-fitting plane
  for (auto &face : faces) {
      std::vector<Kernel::Point_3> points; // face.boundary -> just the ids of the boundary vertices
      for (const int &index : face.boundary) {
          points.emplace_back(vertices[index].x, vertices[index].y, vertices[index].z); // using emplace_back here because we are constructing it
      }
      // Compute the best-fitting plane
      Kernel::Plane_3 plane; // initialize a plane for the set of points of the face
      CGAL::linear_least_squares_fitting_3(points.begin(), points.end(), plane, CGAL::Dimension_tag<0>()); // fit the plane using iterators from the points, CGAL::Dimension_tag<0>() is used to specify to cgal that it deals with points.
      face.best_plane = plane; // store the best fitting plane to the face
  }

  // Step 3: Triangulate faces
  for (auto &face : faces) {
      // For every face fit points onto 2d
      std::vector<Point2> facePoints2D;

      for (int vertex : face.boundary) {
          Point3 point3d(vertices[vertex].x, vertices[vertex].y, vertices[vertex].z);
          Point2 pt = face.best_plane.to_2d(point3d);
          facePoints2D.push_back(pt);
      }

      // create a triangulation for the facePoints2D.
      Triangulation triangulation;

      // insert all the points into the triangulation
      for (const Point2& pt : facePoints2D) {
          triangulation.insert(pt);
      }

      // add begin point again to facePoints2D, so we can make a closed loop with a for loop
      facePoints2D.push_back(facePoints2D[0]);

      // iterate over the points and add pairs of them as constraints
      for (int i = 0; i < facePoints2D.size() - 1; ++i) {
          triangulation.insert_constraint(facePoints2D[i], facePoints2D[i + 1]);
      }

      // when the triangulation is finished assign it to the face
      face.triangulation = triangulation;

      // Step 4: Label triangulation (per face)
      std::list<Triangulation::Face_handle> to_check;
      triangulation.infinite_face()->info().processed = true;
      CGAL_assertion(triangulation.infinite_face()->info().processed == true); // check that the infinite face is always marked as processed
      CGAL_assertion(triangulation.infinite_face()->info().interior == false); // check that the infinite face is always marked as exterior
      to_check.push_back(triangulation.infinite_face());
      while (!to_check.empty()) {
          CGAL_assertion(to_check.front()->info().processed == true);
          for (int neighbour = 0; neighbour < 3; ++neighbour) {
              if (to_check.front()->neighbor(neighbour)->info().processed) {

              } else {
                  to_check.front()->neighbor(neighbour)->info().processed = true;
                  CGAL_assertion(to_check.front()->neighbor(neighbour)->info().processed == true);
                  if (triangulation.is_constrained(Triangulation::Edge(to_check.front(), neighbour))) {
                      to_check.front()->neighbor(neighbour)->info().interior = !to_check.front()->info().interior;
                      to_check.push_back(to_check.front()->neighbor(neighbour));
                  } else {
                      to_check.front()->neighbor(neighbour)->info().interior = to_check.front()->info().interior;
                      to_check.push_back(to_check.front()->neighbor(neighbour));
                  }
              }
          } to_check.pop_front();
      }
      
  }

  // Label triangulation (to do)

  // Export triangles (to do)

  // test

  return 0;
}
