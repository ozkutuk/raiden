#include "parser.h"
#include "tinyxml2.h"
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>

void parser::Scene::loadFromXml(const std::string &filepath) {
    tinyxml2::XMLDocument file;
    std::stringstream stream;

    auto res = file.LoadFile(filepath.c_str());
    if (res) {
        throw std::runtime_error("Error: The xml file cannot be loaded.");
    }

    auto root = file.FirstChild();
    if (!root) {
        throw std::runtime_error("Error: Root is not found.");
    }

    // Get BackgroundColor
    auto element = root->FirstChildElement("BackgroundColor");
    if (element) {
        stream << element->GetText() << std::endl;
    } else {
        stream << "0 0 0" << std::endl;
    }
    stream >> background_color.x >> background_color.y >> background_color.z;

    // Get ShadowRayEpsilon
    element = root->FirstChildElement("ShadowRayEpsilon");
    if (element) {
        stream << element->GetText() << std::endl;
    } else {
        stream << "0.001" << std::endl;
    }
    stream >> shadow_ray_epsilon;

    // Get MaxRecursionDepth
    element = root->FirstChildElement("MaxRecursionDepth");
    if (element) {
        stream << element->GetText() << std::endl;
    } else {
        stream << "0" << std::endl;
    }
    stream >> max_recursion_depth;

    // Get Cameras
    element = root->FirstChildElement("Cameras");
    element = element->FirstChildElement("Camera");
    Camera camera;
    while (element) {

        std::string camera_type;
        if (element->Attribute("type"))
            camera_type = element->Attribute("type");
        else
            camera_type = "default";

        auto child = element->FirstChildElement("Position");
        stream << child->GetText() << std::endl;
        if (camera_type == "lookAt") {
            child = element->FirstChildElement("GazePoint");
            stream << child->GetText() << std::endl;
            child = element->FirstChildElement("Up");
            stream << child->GetText() << std::endl;
            child = element->FirstChildElement("FovY");
        } else {
            child = element->FirstChildElement("Gaze");
            stream << child->GetText() << std::endl;
            child = element->FirstChildElement("Up");
            stream << child->GetText() << std::endl;
            child = element->FirstChildElement("NearPlane");
            stream << child->GetText() << std::endl;
        }

        child = element->FirstChildElement("NearDistance");
        stream << child->GetText() << std::endl;
        child = element->FirstChildElement("ImageResolution");
        stream << child->GetText() << std::endl;
        child = element->FirstChildElement("ImageName");
        stream << child->GetText() << std::endl;
        child = element->FirstChildElement("NumSamples");
        if (child) {
            stream << child->GetText() << std::endl;
        } else {
            stream << "1" << std::endl;
        }
        child = element->FirstChildElement("FocusDistance");
        if (child) {
            stream << child->GetText() << std::endl;
        } else {
            stream << "0" << std::endl; // TODO: 0 may be a valid number for focus
        }
        child = element->FirstChildElement("ApertureSize");
        if (child) {
            stream << child->GetText() << std::endl;
        } else {
            stream << "0" << std::endl;
        }

        stream >> camera.position.x >> camera.position.y >> camera.position.z;
        if (camera_type == "lookAt") {
            Vec3f gaze_point;
            stream >> gaze_point.x >> gaze_point.y >> gaze_point.z;
            Vec3f gaze_vector = gaze_point - camera.position;
            camera.gaze = gaze_vector;
        } else {
            stream >> camera.gaze.x >> camera.gaze.y >> camera.gaze.z;
        }

        stream >> camera.up.x >> camera.up.y >> camera.up.z;

        if (camera_type == "lookAt") {
            double fov_y;
            stream >> fov_y;
            stream >> camera.near_distance;
            stream >> camera.image_width >> camera.image_height;
            const double pi = std::acos(-1);
            double fov_radians = (pi / 180) * fov_y;

            float plane_x, plane_y;
            plane_y = camera.near_distance * std::tan(fov_radians / 2);
            float ratio = static_cast<float>(camera.image_width) / camera.image_height;
            plane_x = ratio * plane_y;

            camera.near_plane.x = -1 * plane_x;
            camera.near_plane.y = plane_x;
            camera.near_plane.z = -1 * plane_y;
            camera.near_plane.w = plane_y;
        } else {
            stream >> camera.near_plane.x >> camera.near_plane.y >> camera.near_plane.z >>
                camera.near_plane.w;
            stream >> camera.near_distance;
            stream >> camera.image_width >> camera.image_height;
        }
        stream >> camera.image_name;
        stream >> camera.n_samples;
        stream >> camera.focus_distance;
        stream >> camera.aperture_size;

        cameras.push_back(camera);

        element = element->NextSiblingElement("Camera");
    }

    // Get Lights
    element = root->FirstChildElement("Lights");
    auto child = element->FirstChildElement("AmbientLight");
    stream << child->GetText() << std::endl;
    stream >> ambient_light.x >> ambient_light.y >> ambient_light.z;
    auto points = element->FirstChildElement("PointLight");
    PointLight point_light;
    while (points) {
        child = points->FirstChildElement("Position");
        stream << child->GetText() << std::endl;
        child = points->FirstChildElement("Intensity");
        stream << child->GetText() << std::endl;

        stream >> point_light.position.x >> point_light.position.y >> point_light.position.z;
        stream >> point_light.intensity.x >> point_light.intensity.y >> point_light.intensity.z;

        point_lights.push_back(point_light);
        points = points->NextSiblingElement("PointLight");
    }
    auto areas = element->FirstChildElement("AreaLight");
    AreaLight area_light;
    while (areas) {
        child = areas->FirstChildElement("Position");
        stream << child->GetText() << std::endl;
        child = areas->FirstChildElement("Radiance");
        stream << child->GetText() << std::endl;
        child = areas->FirstChildElement("Normal");
        stream << child->GetText() << std::endl;
        child = areas->FirstChildElement("Size");
        stream << child->GetText() << std::endl;

        stream >> area_light.position.x >> area_light.position.y >> area_light.position.z;
        stream >> area_light.intensity.x >> area_light.intensity.y >> area_light.intensity.z;
        stream >> area_light.normal.x >> area_light.normal.y >> area_light.normal.z;
        stream >> area_light.size;

        area_lights.push_back(area_light);
        areas = areas->NextSiblingElement("AreaLight");
    }

    // Get Materials
    element = root->FirstChildElement("Materials");
    element = element->FirstChildElement("Material");
    Material material;
    while (element) {
        child = element->FirstChildElement("AmbientReflectance");
        stream << child->GetText() << std::endl;
        child = element->FirstChildElement("DiffuseReflectance");
        stream << child->GetText() << std::endl;
        child = element->FirstChildElement("SpecularReflectance");
        stream << child->GetText() << std::endl;
        child = element->FirstChildElement("MirrorReflectance");
        if (child)
            stream << child->GetText() << std::endl;
        else
            stream << "0 0 0" << std::endl;
        child = element->FirstChildElement("PhongExponent");
        if (child)
            stream << child->GetText() << std::endl;
        else
            stream << "1.0" << std::endl;
        child = element->FirstChildElement("RefractionIndex");
        if (child)
            stream << child->GetText() << std::endl;
        else
            stream << "0.0" << std::endl;
        child = element->FirstChildElement("Transparency");
        if (child)
            stream << child->GetText() << std::endl;
        else
            stream << "1.0 1.0 1.0" << std::endl;
        child = element->FirstChildElement("Roughness");
        if (child)
            stream << child->GetText() << std::endl;
        else
            stream << "0.0" << std::endl;

        stream >> material.ambient.x >> material.ambient.y >> material.ambient.z;
        stream >> material.diffuse.x >> material.diffuse.y >> material.diffuse.z;
        stream >> material.specular.x >> material.specular.y >> material.specular.z;
        stream >> material.mirror.x >> material.mirror.y >> material.mirror.z;
        stream >> material.phong_exponent;
        stream >> material.refractive_index;
        stream >> material.transparency.x >> material.transparency.y >> material.transparency.z;
        stream >> material.roughness;

        materials.push_back(material);
        element = element->NextSiblingElement("Material");
    }

    // Get VertexData
    element = root->FirstChildElement("VertexData");
    if (element) {
        stream << element->GetText() << std::endl;
        Vec3f vertex;
        while (!(stream >> vertex.x).eof()) {
            stream >> vertex.y >> vertex.z;
            vertex_data.push_back(vertex);
        }
        stream.clear();
    }

    // Get Meshes
    element = root->FirstChildElement("Objects");
    element = element->FirstChildElement("Mesh");
    Mesh mesh;
    while (element) {
        child = element->FirstChildElement("Material");
        stream << child->GetText() << std::endl;
        stream >> mesh.material_id;

        child = element->FirstChildElement("Faces");

        auto ply_file = child->Attribute("plyFile");

        if (ply_file) {
            // TODO: import ply
            std::cerr << "Ply support not yet implemented" << std::endl;
        } else {
            stream << child->GetText() << std::endl;
            Face face;
            while (!(stream >> face.v0_id).eof()) {
                stream >> face.v1_id >> face.v2_id;
                mesh.faces.push_back(face);
            }
            stream.clear();

            meshes.push_back(mesh);
            mesh.faces.clear();
        }
        element = element->NextSiblingElement("Mesh");
    }
    stream.clear();

    // Get Triangles
    element = root->FirstChildElement("Objects");
    element = element->FirstChildElement("Triangle");
    Triangle triangle;
    while (element) {
        child = element->FirstChildElement("Material");
        stream << child->GetText() << std::endl;
        stream >> triangle.material_id;

        child = element->FirstChildElement("Indices");
        stream << child->GetText() << std::endl;
        stream >> triangle.indices.v0_id >> triangle.indices.v1_id >> triangle.indices.v2_id;

        triangles.push_back(triangle);
        element = element->NextSiblingElement("Triangle");
    }

    // Get Spheres
    element = root->FirstChildElement("Objects");
    element = element->FirstChildElement("Sphere");
    Sphere sphere;
    while (element) {
        child = element->FirstChildElement("Material");
        stream << child->GetText() << std::endl;
        stream >> sphere.material_id;

        child = element->FirstChildElement("Center");
        stream << child->GetText() << std::endl;
        stream >> sphere.center_vertex_id;

        child = element->FirstChildElement("Radius");
        stream << child->GetText() << std::endl;
        stream >> sphere.radius;

        spheres.push_back(sphere);
        element = element->NextSiblingElement("Sphere");
    }
}
