#pragma once

#include <vector>

class Mat4;

class Vec2 {
    public:
    float x = 0.0f, y = 0.0f;
    Vec2(float a, float b) {
        x = a;
        y = b;
    }
};

class Vec3
{
public:
    float x = 0.0f, y = 0.0f, z = 0.0f;
    Vec3() {}
    Vec3(float a, float b, float c)
    {
        x = a;
        y = b;
        z = c;
    }
    Vec3 operator+(Vec3 &v);
    Vec3 operator-(Vec3 &v);
    float &operator[](int idx);

    float norm();
    Vec3 normalize();
};

Vec3 cross(Vec3 &v1, Vec3 &v2);
float dot(Vec3 &v1, Vec3 &v2);

class Vec4
{
public:
    float x = 0.0f, y = 0.0f, z = 0.0f, w = 0.0f;
    Vec4() {}
    Vec4(Vec3 &v)
    {
        x = v.x;
        y = v.y;
        z = v.z;
        w = 1;
    }
    float &operator[](int idx);
};

class Mat4
{
    std::vector<std::vector<float>> data;

public:
    Mat4();
    std::vector<float> &operator[](const int row);
    Vec4 operator*(Vec4 &v);
};

Mat4 lookat(Vec3 &eye, Vec3 &center, Vec3 &up);
/**
 * Perspective projection
 */
Mat4 proj(float fov, float near, float far, float aspect);