#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>

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
  FaceInfo() {
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


//const std::string input_file = "/mnt/c/Users/LarsB/OneDrive/Documenten/GitHub/GEO1004_ass_1/hw1_obj_files/NL.IMBAG.Pand.0503100000000138-0.obj";
const std::string input_file = "/mnt/c/Users/LarsB/OneDrive/Documenten/GitHub/GEO1004_ass_1/hw1_obj_files/NL.IMBAG.Pand.0503100000025027-0.obj";
//const std::string output_file = "/mnt/c/Users/LarsB/OneDrive/Documenten/GitHub/GEO1004_ass_1/hw1_output_files/simple.obj";
const std::string output_file = "/mnt/c/Users/LarsB/OneDrive/Documenten/GitHub/GEO1004_ass_1/hw1_output_files/complex.obj";

// struct are like public classes
struct Vertex {
  double x, y, z;
};

struct Face {
  std::list<int> boundary;
  Kernel::Plane_3 best_plane;
  Triangulation triangulation;
};

// struct to maintain a strict weak ordering of Point3 objects
// to check for duplicate Point3 values within epsilon tolerance later
struct epsilonCompare {
    bool operator() (const Point3& a, const Point3& b) const {
        const double EPSILON = 0.0001;

        if (std::abs(a.x() - b.x()) < EPSILON) {
            if (std::abs(a.y() - b.y()) < EPSILON) {
                return (a.z() + EPSILON) < b.z(); // Compare z if x and y are "close enough"
            }
            return (a.y() + EPSILON) < b.y(); // Compare y if x is "close enough"
        }
        return (a.x() + EPSILON) < b.x(); // Compare x coordinates
    }
};

int main(int argc, const char * argv[]) {
    // initialize vectors to read in OBJ data
    std::vector<Vertex> vertices;
    std::vector<Face> faces;

    // initialize map to check for duplicate vertices using epsilonCompare
    std::map<Point3, int, epsilonCompare> vertex_indices;

    // initialize vertex and face vectors to store output vertices and faces before writing to OBJ file
    std::vector<Point3> output_vertices;
    std::vector<std::vector<int>> output_faces;

    // Step 1: Import the OBJ file and load it into a data structure
    std::ifstream input_stream;
    input_stream.open(input_file);
    if (input_stream.is_open()) {
        std::string line;

        // Parse line by line
        while (getline(input_stream, line)) {
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

  // Step 2: For each of the faces compute the best-fitting plane
  for (auto &face : faces) {
      // initialize vector for storing vertex coordinates of the face
      std::vector<Kernel::Point_3> vertex_coord;

      // iterate over the indices from the face and use them to construct coordinates
      for (const int &index : face.boundary) {
          vertex_coord.emplace_back(vertices[index].x, vertices[index].y, vertices[index].z); // using emplace_back here because we are constructing the point
      }

      // compute the best-fitting plane
      Kernel::Plane_3 plane; // initialize a plane for the set of points of the face
      CGAL::linear_least_squares_fitting_3(vertex_coord.begin(), vertex_coord.end(), plane, CGAL::Dimension_tag<0>()); // fit the plane using iterators from the points, CGAL::Dimension_tag<0>() is used to specify to cgal that it deals with points.
      face.best_plane = plane; // store the best fitting plane to the face
  }

  // Step 3: Triangulate faces
  for (auto &face : faces) {

      // vector to store 2D points of face
      std::vector<Point2> facePoints2D;

      // iterate over 3D vertices of face and store the 2D vertices using best-fitting plane
      for (int vertex : face.boundary) {
          Point3 point3d(vertices[vertex].x, vertices[vertex].y, vertices[vertex].z);
          Point2 pt = face.best_plane.to_2d(point3d);
          facePoints2D.push_back(pt);
      }

      // initialize a triangulation object for the facePoints2D.
      Triangulation triangulation;

      // insert all the 2D points into the triangulation
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

      // Step 4: Label face triangulation with interior or exterior
      std::list<Triangulation::Face_handle> to_check;
      triangulation.infinite_face()->info().processed = true;
      CGAL_assertion(triangulation.infinite_face()->info().processed == true);
      CGAL_assertion(triangulation.infinite_face()->info().interior == false);
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

      // Step 5: export the interior faces back to 3d using the best fitting plane
      // for storing point indices
      int index = 0;

      // iterate over faces of the triangulation
      for (auto it = triangulation.finite_faces_begin(); it != triangulation.finite_faces_end(); ++it) {

          if (it->info().interior) {
              std::vector<int> face_vertex_indices;
              for (int i = 0; i < 3; ++i) {
                  Point2 pt2 = it->vertex(i)->point();
                  Point3 pt3 = face.best_plane.to_3d(pt2);

                  // Check if this vertex is already added
                  if (vertex_indices.find(pt3) == vertex_indices.end()) {
                      // new vertex, add to vector and map
                      output_vertices.push_back(pt3);
                      vertex_indices[pt3] = ++index;
                  }

                  // use the vertex_indices to find the vertex number to store for the face
                  face_vertex_indices.push_back(vertex_indices[pt3]);
              }

              // store the vertex indices for the face
              output_faces.push_back(face_vertex_indices);
          }
      }

  }

  // Step 6 : Output to OBJ file
  std::ofstream output_stream(output_file);
  if (output_stream.is_open()) {
      // write vertices
      for (const auto &vertex : output_vertices) {
          output_stream << "v " << vertex.x() << " " << vertex.y() << " " << vertex.z() << std::endl;
      }

      // write faces
      for (const auto &face : output_faces) {
          output_stream << "f";
          for (int vertex_index : face) {
              output_stream << " " << vertex_index;
          }
          output_stream << std::endl;
      }
      output_stream.close();
  } else {
      std::cerr << "Could not open output file." << std::endl;
      return -1;
  }

  return 0;
}
