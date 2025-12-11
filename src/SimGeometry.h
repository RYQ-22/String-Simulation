#ifndef __SIM_GEOMETRY_H__
#define __SIM_GEOMETRY_H__

// std
#include <vector>
// Eigen
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Sparse>
// polyscope
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

class OctopusHead {
private:
    std::string name;
    std::vector<Eigen::Vector3d> vertices;
    std::vector<std::vector<size_t>> faces;
    polyscope::SurfaceMesh* mesh_ptr = nullptr;

public:
    OctopusHead(std::string name, double radius, int stacks = 10, int slices = 20) : name(name) {
        generateGeometry(radius, stacks, slices);
    }

    void registerMesh() {
        mesh_ptr = polyscope::registerSurfaceMesh(name, vertices, faces);
        mesh_ptr->setSmoothShade(true);
        mesh_ptr->setSurfaceColor({0.8, 0.3, 0.3});
    }

private:
    void generateGeometry(double r, int stacks, int slices) {
        // 1. generate vertices
        vertices.clear();

        vertices.push_back(Eigen::Vector3d(0, r, 0)); 

        for (int i = 0; i < stacks; ++i) {
            double phi = (M_PI / 2.0) * (double)(i + 1) / (double)stacks;
            
            for (int j = 0; j < slices; ++j) {
                double theta = 2.0 * M_PI * (double)j / (double)slices;

                double x = r * std::sin(phi) * std::cos(theta);
                double z = r * std::sin(phi) * std::sin(theta);
                double y = r * std::cos(phi);

                vertices.push_back(Eigen::Vector3d(x, y, z));
            }
        }

        vertices.push_back(Eigen::Vector3d(0, 0, 0));
        size_t center_idx = vertices.size() - 1;

        // generate faces
        faces.clear();
        
        for (int j = 0; j < slices; ++j) {
            size_t p0 = 0;
            size_t p1 = 1 + j;
            size_t p2 = 1 + (j + 1) % slices;
            faces.push_back({p0, p1, p2});
        }

        for (int i = 0; i < stacks - 1; ++i) {
            size_t ring_start = 1 + i * slices;
            size_t next_ring_start = 1 + (i + 1) * slices;

            for (int j = 0; j < slices; ++j) {
                size_t next_j = (j + 1) % slices;
                
                size_t v0 = ring_start + j;
                size_t v1 = next_ring_start + j;
                size_t v2 = next_ring_start + next_j;
                size_t v3 = ring_start + next_j;

                faces.push_back({v0, v1, v2});
                faces.push_back({v0, v2, v3});
            }
        }

        size_t last_ring_start = 1 + (stacks - 1) * slices;
        for (int j = 0; j < slices; ++j) {
            size_t v_curr = last_ring_start + j;
            size_t v_next = last_ring_start + (j + 1) % slices;
            
            faces.push_back({center_idx, v_next, v_curr});
        }
    }
};

#endif
