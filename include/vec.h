#pragma once

#include <vector>

#ifndef __CUDACC__
#define __host__
#define __device__
#endif

class Mat4;

class Vec2 {
 public:
  float x = 0.0f, y = 0.0f;
  __host__ __device__ Vec2() : x(0.0f), y(0.0f) {}
  __host__ __device__ Vec2(float a, float b) : x(a), y(b) {}
};

class Vec3 {
 public:
  float x = 0.0f, y = 0.0f, z = 0.0f;
  __host__ __device__ Vec3() {}
  __host__ __device__ Vec3(float a, float b, float c) : x(a), y(b), z(c) {}
  __host__ __device__ Vec3 operator+(Vec3 v);
  __host__ __device__ Vec3 operator+(float f);
  __host__ __device__ Vec3 operator-(Vec3 v);
  __host__ __device__ Vec3 operator*(float f);
  __host__ __device__ Vec3 operator/(float f);
  __host__ __device__ float &operator[](int idx);

  __host__ __device__ float norm();
  __host__ __device__ Vec3 normalize();
};

__host__ __device__ Vec3 cross(Vec3 &v1, Vec3 &v2);
__host__ __device__ float dot(Vec3 &v1, Vec3 &v2);

class Vec4 {
 public:
  float x = 0.0f, y = 0.0f, z = 0.0f, w = 0.0f;
  __host__ __device__ Vec4() {}
  __host__ __device__ Vec4(Vec3 &v) {
    x = v.x;
    y = v.y;
    z = v.z;
    w = 1;
  }
  __host__ __device__ float &operator[](int idx);
};

class Mat4 {
  std::vector<std::vector<float>> data;

 public:
  Mat4();
  std::vector<float> &operator[](int row);
  Vec4 operator*(Vec4 &v);
  Mat4 operator*(Mat4 &m);
};

Mat4 lookat(Vec3 &eye, Vec3 &center, Vec3 &up);
/**
 * Perspective projection
 */
Mat4 proj(float fov, float near, float far, float aspect);