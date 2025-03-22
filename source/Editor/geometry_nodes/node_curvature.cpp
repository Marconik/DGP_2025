#include <pxr/base/vt/array.h>

#include <vector>

#include "GCore/Components/MeshOperand.h"
#include "OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh"
#include "geom_node_base.h"
#include "nodes/core/def/node_def.hpp"

typedef OpenMesh::TriMesh_ArrayKernelT<> MyMesh;

namespace Surf
{
    float TriAngle(MyMesh::Point A, MyMesh::Point B, MyMesh::Point C);
    float TriArea(MyMesh::Point A,MyMesh::Point B,MyMesh::Point C);
    float BaryCen(MyMesh::Point A,MyMesh::Point B,MyMesh::Point C);

    float TriAngle(MyMesh::Point A, MyMesh::Point B, MyMesh::Point C)
    {
        auto ea = (B - C).norm();
        auto eb = (A - B).norm();
        auto ec = (A - C).norm();
        auto r = (eb * eb + ec * ec - ea * ea) / (2 * eb * ec);
        if (std::abs(r) <= 1)
            return std::acos(r);
        else
            return 0;
    }
    float TriArea(MyMesh::Point A,MyMesh::Point B,MyMesh::Point C)
    {
        auto ea = (B - C).norm();
        auto eb = (A - B).norm();
        auto ec = (A - C).norm();
        auto r = (eb * eb + ec * ec - ea * ea) / (2 * eb * ec);
        if(std::abs(r) <= 1)
        {
            auto s = std::sqrt(1 - r * r);
            return s * eb * ec / 2;
        }
        else
            return 0;
    }
    float BaryCen(MyMesh::Point A,MyMesh::Point B,MyMesh::Point C)
    {
        auto G=(A+B+C)/3;
        auto M1=(A+C)/2;
        auto M2=(A+B)/2;
        return TriArea(A,G,M1)+TriArea(A,G,M2);
    }
}


void compute_mean_curvature(
    const MyMesh& omesh,
    pxr::VtArray<float>& mean_curvature)
{
    // TODO: Implement the mean curvature computation
    //  You need to fill in `mean_curvature`
    for(auto v_it = omesh.vertices_begin(); v_it != omesh.vertices_end(); ++v_it)
    {
        MyMesh::Point v = omesh.point(*v_it);
        std::vector<MyMesh::Point> node_neighbors;
        for(auto vv_it = omesh.cvv_iter(*v_it); vv_it.is_valid(); ++vv_it)
        {
            node_neighbors.push_back(omesh.point(*vv_it));
        }
        auto n = node_neighbors.size();     
        auto area = Surf::BaryCen(v, node_neighbors[0], node_neighbors[n-1]);
        auto laplacian=(node_neighbors[0]-v)/std::tan(Surf::TriAngle(node_neighbors[n-1],v,node_neighbors[0]));
        laplacian += (node_neighbors[n-1] - v) /
                     std::tan(Surf::TriAngle(node_neighbors[0], v, node_neighbors[n-1]));
        for(int i = 0;i<n-1;++i)
        {
            area += Surf::BaryCen(v, node_neighbors[i], node_neighbors[i+1]);
            laplacian += (node_neighbors[i] - v) /
                         std::tan(Surf::TriAngle(
                             node_neighbors[i + 1], node_neighbors[i], v));
            laplacian += (node_neighbors[i + 1] - v) /
                         std::tan(Surf::TriAngle(
                             node_neighbors[i], node_neighbors[i + 1], v));
        }
        auto r=laplacian.norm();
        r=r/(4*area);
        mean_curvature.push_back(r);
    }
}

void compute_gaussian_curvature(
    const MyMesh& omesh,
    pxr::VtArray<float>& gaussian_curvature)
{ 
    // TODO: Implement the Gaussian curvature computation
    //  You need to fill in `gaussian_curvature`
    for(auto v_it = omesh.vertices_begin(); v_it != omesh.vertices_end(); ++v_it)
    {
        MyMesh::Point v = omesh.point(*v_it);
        float omega = 0;
        float area = 0;
        std::vector<MyMesh::Point> node_neighbors;
        for (auto vv_it = omesh.cvv_iter(*v_it); vv_it.is_valid(); ++vv_it)
        {
            node_neighbors.push_back(omesh.point(*vv_it));
        }
        auto n = node_neighbors.size();
        for(int i=0;i<n-1;++i)
        {
            omega +=
                Surf::TriAngle(v, node_neighbors[i], node_neighbors[i + 1]);
            area += Surf::TriArea(v, node_neighbors[i], node_neighbors[i + 1]);
        }
        omega += Surf::TriAngle(v, node_neighbors[0], node_neighbors[n - 1]);
        area += Surf::TriArea(v, node_neighbors[0], node_neighbors[n - 1]);
        auto K=(2*std::acos(-1)-omega)/area;
        if(K>=0)
            gaussian_curvature.push_back(K);
        else
            gaussian_curvature.push_back(0);
    }
}

NODE_DEF_OPEN_SCOPE

NODE_DECLARATION_FUNCTION(mean_curvature)
{
    b.add_input<Geometry>("Mesh");
    b.add_output<pxr::VtArray<float>>("Mean Curvature");
}

NODE_EXECUTION_FUNCTION(mean_curvature)
{
    auto geometry = params.get_input<Geometry>("Mesh");
    auto mesh = geometry.get_component<MeshComponent>();
    auto vertices = mesh->get_vertices();
    auto face_vertex_indices = mesh->get_face_vertex_indices();
    auto face_vertex_counts = mesh->get_face_vertex_counts();

    // Convert the mesh to OpenMesh
    MyMesh omesh;

    // Add vertices
    std::vector<OpenMesh::VertexHandle> vhandles;
    vhandles.reserve(vertices.size());

    for (auto vertex : vertices) {
        omesh.add_vertex(OpenMesh::Vec3f(vertex[0], vertex[1], vertex[2]));
    }

    // Add faces
    size_t start = 0;
    for (int face_vertex_count : face_vertex_counts) {
        std::vector<OpenMesh::VertexHandle> face;
        face.reserve(face_vertex_count);
        for (int j = 0; j < face_vertex_count; j++) {
            face.push_back(
                OpenMesh::VertexHandle(face_vertex_indices[start + j]));
        }
        omesh.add_face(face);
        start += face_vertex_count;
    }

    // Compute mean curvature
    pxr::VtArray<float> mean_curvature;
    mean_curvature.reserve(omesh.n_vertices());

    compute_mean_curvature(omesh, mean_curvature);

    params.set_output("Mean Curvature", mean_curvature);

    return true;
}

NODE_DECLARATION_UI(mean_curvature);

NODE_DECLARATION_FUNCTION(gaussian_curvature)
{
    b.add_input<Geometry>("Mesh");
    b.add_output<pxr::VtArray<float>>("Gaussian Curvature");
}

NODE_EXECUTION_FUNCTION(gaussian_curvature)
{
    auto geometry = params.get_input<Geometry>("Mesh");
    auto mesh = geometry.get_component<MeshComponent>();
    auto vertices = mesh->get_vertices();
    auto face_vertex_indices = mesh->get_face_vertex_indices();
    auto face_vertex_counts = mesh->get_face_vertex_counts();

    // Convert the mesh to OpenMesh
    MyMesh omesh;

    // Add vertices
    std::vector<OpenMesh::VertexHandle> vhandles;
    vhandles.reserve(vertices.size());

    for (auto vertex : vertices) {
        omesh.add_vertex(OpenMesh::Vec3f(vertex[0], vertex[1], vertex[2]));
    }

    // Add faces
    size_t start = 0;
    for (int face_vertex_count : face_vertex_counts) {
        std::vector<OpenMesh::VertexHandle> face;
        face.reserve(face_vertex_count);
        for (int j = 0; j < face_vertex_count; j++) {
            face.push_back(
                OpenMesh::VertexHandle(face_vertex_indices[start + j]));
        }
        omesh.add_face(face);
        start += face_vertex_count;
    }

    // Compute Gaussian curvature
    pxr::VtArray<float> gaussian_curvature;
    gaussian_curvature.reserve(omesh.n_vertices());

    compute_gaussian_curvature(omesh, gaussian_curvature);

    params.set_output("Gaussian Curvature", gaussian_curvature);

    return true;
}

NODE_DECLARATION_UI(gaussian_curvature);

NODE_DEF_CLOSE_SCOPE
