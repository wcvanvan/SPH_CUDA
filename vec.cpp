#include <assert.h>
#include <cmath>
#include "const.h"
#include "vec.h"

Vec3 Vec3::operator+(Vec3 &v)
{
    Vec3 result;
    result.x = x + v.x;
    result.y = y + v.y;
    result.z = z + v.z;
    return result;
}

Vec3 Vec3::operator-(Vec3 &v)
{
    Vec3 result;
    result.x = x - v.x;
    result.y = y - v.y;
    result.z = z - v.z;
    return result;
}

float &Vec3::operator[](int idx)
{
    assert(idx >= 0 && idx < 3);
    if (idx == 0)
        return x;
    if (idx == 1)
        return y;
    return z;
}

float &Vec4::operator[](int idx)
{
    assert(idx >= 0 && idx < 4);
    if (idx == 0)
        return x;
    if (idx == 1)
        return y;
    if (idx == 2)
        return z;
    return w;
}

Vec3 cross(Vec3 &v1, Vec3 &v2)
{
    Vec3 result;
    result.x = v1.y * v2.z - v1.z * v2.y;
    result.y = v1.z * v2.x - v1.x * v2.z;
    result.z = v1.x * v2.y - v1.y * v2.x;
    return result;
}

float dot(Vec3 &v1, Vec3 &v2)
{
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

float Vec3::norm()
{
    return std::sqrt(x * x + y * y + z * z);
}

Vec3 Vec3::normalize()
{
    Vec3 result = *this;
    float norm = this->norm();
    result.x /= norm;
    result.y /= norm;
    result.z /= norm;
    return result;
}

Mat4::Mat4()
{
    this->data = std::vector<std::vector<float>>(4, std::vector<float>(4, 0.0f));
}

std::vector<float> &Mat4::operator[](const int row)
{
    assert(row >= 0 && row < 4);
    return data[row];
}

Mat4 lookat(Vec3 &eye, Vec3 &center, Vec3 &up)
{
    Vec3 camForward = Vec3(center - eye).normalize();
    Vec3 camRight = cross(camForward, up).normalize();
    Vec3 camUp = cross(camRight, camForward);

    float translation_x = dot(camRight, eye) * (-1);
    float translation_y = dot(camUp, eye) * (-1);
    float translation_z = dot(camForward, eye);

    Mat4 viewMat;
    for (int i = 0; i < 3; i++)
    {
        viewMat[0][i] = camRight[i];
        viewMat[1][i] = camUp[i];
        viewMat[2][i] = -camForward[i];
    }
    viewMat[0][3] = translation_x;
    viewMat[1][3] = translation_y;
    viewMat[2][3] = translation_z;
    viewMat[3][3] = 1;
    return viewMat;
}

Mat4 proj(float fov, float near, float far, float aspectRatio)
{
    Mat4 projMat;
    float tanHalfFov = tan(fov / 2.0f);
    projMat[0][0] = 1.0f / (aspectRatio * tanHalfFov);
    projMat[1][1] = 1.0f / tanHalfFov;
    projMat[2][2] = -(far + near) / (far - near);
    projMat[2][3] = -(2.0f * far * near) / (far - near);
    projMat[3][2] = -1.0f;
    return projMat;
}

Vec4 Mat4::operator*(Vec4 &v)
{
    Vec4 result;
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            result[i] += this->data[i][j] * v[j];
        }
    }
    return result;
}