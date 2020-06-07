
#pragma once

#include <cmath>
#include <cassert>
#include <iostream>
#include <iomanip>

namespace mafs {

using i32 = int32_t;
using i64 = int64_t;

using Scalar = double;
const inline Scalar pi = 3.14159265;

#define MAFS_FOR(N)  for (i64 i = 0; i < (N); i++)
#define MAFS_FOR2(N, M)  MAFS_FOR(N) for (i64 j = 0; j < (M); j++)
#define vec_fn  template <i64 n> auto
#define mat_fn  template <i64 n, i64 m> auto


template <i64 n>
union Vector {
    static_assert(n > 1, "Vector size must be at least 2");

    Scalar data[n];
    struct { Scalar x, y, z, w; };

    auto operator[](i64 i) -> Scalar& {
        return data[i];
    }
};

#define SPECIALIZE_VECTOR(N, components)  \
template<> union Vector<N> {              \
    Scalar data[N];                       \
    struct { Scalar components; };        \
    auto operator[](i64 i) -> Scalar& {   \
        return data[i];                   \
    }                                     \
}                                         \

#define COMMA ,
SPECIALIZE_VECTOR(2, x COMMA y);
SPECIALIZE_VECTOR(3, x COMMA y COMMA z);
SPECIALIZE_VECTOR(4, x COMMA y COMMA z COMMA w);
#undef COMMA
#undef SPECIALIZE_VECTOR


template <i64 n, i64 m>
struct Matrix {
    Scalar data[n][m];

    auto operator[](i64 i) -> Scalar* {
        return data[i];
    }
};


union Quaternion {
    Vector<4> vec;
    struct { Scalar r, i, j, k; };
};
const inline Quaternion identity_quat = { 1, 0, 0, 0 };



inline auto sign(Scalar x) -> Scalar {
    if (x < 0.) return -1.;
    return 1.;
}


#define OPERATION_ON_VECTORS(op)                                               \
mat_fn operator op(Vector<n> lhs, Vector<m> rhs) -> Vector<std::max(n, m)> {   \
    const i64 greater  = std::max(n, m);                                       \
    const i64 smaller = std::min(n, m);                                        \
                                                                               \
    Vector<greater> out;                                                       \
                                                                               \
    MAFS_FOR(greater) {                                                        \
        if (i < smaller) out[i] = rhs[i] op lhs[i];                            \
        else out[i] = n > m ? lhs[i] : rhs[i];                                 \
    }                                                                          \
                                                                               \
    return out;                                                                \
}                                                                              \

OPERATION_ON_VECTORS(+)
OPERATION_ON_VECTORS(-)
#undef OPERATION_ON_VECTOR


#define OPERATION_ON_VECTOR_AND_SCALAR(op)                                     \
vec_fn operator op(Vector<n> lhs, Scalar rhs) -> Vector<n> {                   \
    Vector<n> out;                                                             \
    MAFS_FOR(n) {                                                              \
        out[i] = lhs[i] op rhs;                                                \
    }                                                                          \
    return out;                                                                \
}                                                                              \

OPERATION_ON_VECTOR_AND_SCALAR(+)
OPERATION_ON_VECTOR_AND_SCALAR(-)
OPERATION_ON_VECTOR_AND_SCALAR(*)
OPERATION_ON_VECTOR_AND_SCALAR(/)
#undef OPERATION_ON_VECTOR_AND_SCALAR


#define OPERATION_ON_SCALAR_AND_VECTOR(op)                                     \
vec_fn operator op(Scalar lhs, Vector<n> rhs) -> Vector<n> {                   \
    return rhs op lhs;                                                         \
}                                                                              \

OPERATION_ON_SCALAR_AND_VECTOR(+)
OPERATION_ON_SCALAR_AND_VECTOR(*)
#undef OPERATION_ON_SCALAR_AND_VECTOR


#define OPERATION_ON_MATRICES(op)                                              \
mat_fn operator op(Matrix<n, m> lhs, Matrix<n, m> rhs) -> Matrix<n, m> {       \
    Matrix<n, m> out;                                                          \
    MAFS_FOR2(n, m) {                                                          \
        out[i][j] = lhs[i][j] op rhs[i][j];                                    \
    }                                                                          \
    return out;                                                                \
}                                                                              \

OPERATION_ON_MATRICES(+)
OPERATION_ON_MATRICES(-)
#undef OPERATION_ON_MATRICES

#define OPERATION_ON_MATRIX_AND_SCALAR(op)                                     \
mat_fn operator op(Matrix<n, m> lhs, Scalar rhs) -> Matrix<n, m> {             \
    Matrix<n, m> out;                                                          \
    MAFS_FOR2(n, m) {                                                          \
        out[i][j] = lhs[i][j] op rhs;                                          \
    }                                                                          \
    return out;                                                                \
}                                                                              \

OPERATION_ON_MATRIX_AND_SCALAR(+)
OPERATION_ON_MATRIX_AND_SCALAR(-)
OPERATION_ON_MATRIX_AND_SCALAR(*)
OPERATION_ON_MATRIX_AND_SCALAR(/)
#undef OPERATION_ON_MATRIX_AND_SCALAR

#define OPERATION_ON_SCALAR_AND_MATRIX(op)                                     \
mat_fn operator op(Scalar lhs, Matrix<n, m> rhs) -> Matrix<n, m> {             \
    return rhs op lhs;                                                         \
}                                                                              \

OPERATION_ON_SCALAR_AND_MATRIX(+)
OPERATION_ON_SCALAR_AND_MATRIX(*)
#undef OPERATION_ON_SCALAR_AND_MATRIX




vec_fn operator==(Vector<n> a, Vector<n> b) -> bool {
    MAFS_FOR(n) if (a[i] != b[i]) return false;
    return true;
}

vec_fn operator-(Vector<n> vec) -> Vector<n> {
    Vector<n> neg = { 0. };
    MAFS_FOR(n) neg[i] = -vec[i];
    return neg;
}

vec_fn dot(Vector<n> lhs, Vector<n> rhs) -> Scalar {
    Scalar sum = 0.;
    MAFS_FOR(n) sum += lhs[i] * rhs[i];
    return sum;
}

inline auto cross(Vector<3> lhs, Vector<3> rhs) -> Vector<3> {
    const Scalar x = { lhs.y * rhs.z - lhs.z * rhs.y };
    const Scalar y = { lhs.z * rhs.x - lhs.x * rhs.z };
    const Scalar z = { lhs.x * rhs.y - lhs.y * rhs.x };
    return { x, y, z };
}

vec_fn magnitude(Vector<n> vec) -> Scalar {
    Scalar abs = 0.;
    MAFS_FOR(n) abs += vec[i] * vec[i];
    return std::sqrt(abs);
}

vec_fn normalise(Vector<n> vec) -> Vector<n> {
    Scalar abs = magnitude(vec);
    return vec / abs;
}

vec_fn are_linearly_dependant(Vector<n> a, Vector<n> b) -> bool {
    if (normalise(a) == normalise(b) || normalise(a) == -normalise(b)) return true;
    return false;
}

mat_fn vector_from_column(Matrix<n, m> mat, i64 c) -> Vector<m> {
    Vector<m> out;
    MAFS_FOR(m) {
        out[i] = mat[c][i];
    }
    return out;
}



mat_fn transform(Matrix<n, m> mat, Vector<n> vec) -> Vector<m> {
    Vector<m> out = { 0 };
    MAFS_FOR(n) {
        out = out + vector_from_column(mat, i) * vec[i];
    }
    return out;
}

mat_fn operator*(Matrix<n, m> mat, Vector<n> vec) -> Vector<m> {
    return transform(mat, vec);
}

mat_fn identity_matrix() -> Matrix<n, m> {
    Matrix<n, m> mat = { 0. };
    MAFS_FOR(std::min(n, m)) mat[i][i] = 1.f;
    return mat;
}

mat_fn transpose(Matrix<n, m> mat) -> Matrix<m, n> {
    Matrix<m, n> out;
    MAFS_FOR(n) {
        for (i64 j = 0; j < m; j++) out[i][j] = mat[j][i];
    }
    return out;
}

vec_fn compose(Matrix<n, n> a, Matrix<n, n> b) -> Matrix<n, n> {
    Matrix<n, n> out;
    MAFS_FOR(n) {
        Vector<n> col = vector_from_column(b, i);
        col = transform(a, col);
        for (i64 j = 0; j < n; j++) out[i][j] = col[j];
    }
    return out;
}

vec_fn operator*(Matrix<n, n> lhs, Matrix<n, n> rhs) -> Matrix<n, n> {
    return compose(lhs, rhs);
}

vec_fn determinant(Matrix<n, n> mat) -> Scalar {
    if (n == 2)  return mat[0][0] *  mat[1][1] - mat[1][0] * mat[0][1];
    if (n == 3)  return mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1])
                      - mat[1][0] * (mat[0][1] * mat[2][2] - mat[0][2] * mat[2][1])
                      + mat[2][0] * (mat[0][1] * mat[1][2] - mat[1][1] * mat[0][2]);
    assert(false && "Determinant of arbitrary square matrices will be added in the future");
    return 0.;
}


enum class Direction: i64 {
    x, y, z, w
};

vec_fn unit_vector(Direction dir) -> Vector<n> {
    assert((i64) dir < n && "Vector is too small to have a component in this direction");

    Vector<n> unit;
    MAFS_FOR(n) unit[i] = 0.;
    unit[(i64) dir] = 1.;

    return unit;
}



inline auto conjugate(Quaternion quat) -> Quaternion {
    return Quaternion { quat.r, -quat.i, -quat.j, -quat.k };
}

inline auto inverse(Quaternion quat) -> Quaternion {
    Quaternion out;
    out.vec = { conjugate(quat).vec / magnitude(quat.vec) };
    return out;
}

inline auto operator*(Quaternion lhs, Quaternion rhs) -> Quaternion {
    Vector<3> lhs_vector = { lhs.i, lhs.j, lhs.k };
    Vector<3> rhs_vector = { rhs.i, rhs.j, rhs.k };
    Vector<3> product = cross(lhs_vector, rhs_vector);

    Vector<3> vector_part = product + rhs_vector * lhs.r + lhs_vector * rhs.r;
    Scalar scalar_part = lhs.r * rhs.r - dot(lhs_vector, rhs_vector);

    return Quaternion { scalar_part, vector_part[0], vector_part[1], vector_part[2] };
}

inline auto euler_to_quat(Vector<3> euler_angles) -> Quaternion {
    const Scalar roll  = euler_angles.x / 2.;
    const Scalar pitch = euler_angles.y / 2.;
    const Scalar yaw   = euler_angles.z / 2.;

    Quaternion quat;

    quat.r = std::cos(roll) * std::cos(pitch) * std::cos(yaw) + std::sin(roll) * std::sin(pitch) * std::sin(yaw);
    quat.i = std::sin(roll) * std::cos(pitch) * std::cos(yaw) - std::cos(roll) * std::sin(pitch) * std::sin(yaw);
    quat.j = std::cos(roll) * std::sin(pitch) * std::cos(yaw) + std::sin(roll) * std::cos(pitch) * std::sin(yaw);
    quat.k = std::cos(roll) * std::cos(pitch) * std::sin(yaw) - std::sin(roll) * std::sin(pitch) * std::cos(yaw);

    quat.vec = normalise(quat.vec);
    return quat;
}

inline auto quat_to_euler(Quaternion quat) -> Vector<3> {
    Vector<3> euler_angles = { 0 };
    quat.vec = normalise(quat.vec);
    {
        const Scalar a =      2. * (quat.r * quat.i + quat.j * quat.k);
        const Scalar b = 1. - 2. * (quat.i * quat.i + quat.j * quat.j);
        euler_angles.x = std::atan2(a, b);
    }

    const Scalar intermediate   = 2. * (quat.r * quat.j - quat.k * quat.i);
    if (std::abs(intermediate) >= 1.) euler_angles.y = sign(intermediate) * pi / 2.;
    else euler_angles.y = std::asin(intermediate);

    {
        const Scalar a =      2. * (quat.r * quat.k + quat.i * quat.j);
        const Scalar b = 1. - 2. * (quat.j * quat.j + quat.k * quat.k);
        euler_angles.z = std::atan2(a, b);
    }
    return euler_angles;
}

inline auto rotation_matrix_of(Quaternion quat) -> Matrix<3, 3> {
    quat.vec = normalise(quat.vec);
    Scalar x = quat.i, y = quat.j, z = quat.k, w = quat.r;
    return transpose(Matrix<3, 3> { 1. - 2*y*y - 2*z*z, 2*x*y - 2*x*z, 2*x*z + 2*y*w,
                                    2*x*y + 2*z*w, 1. - 2*x*x - 2*z*z, 2*y*z - 2*x*w,
                                    2*x*z - 2*y*w, 2*y*z + 2*x*w, 1. - 2*x*x - 2*y*y });
}

inline auto nlerp_quaternions(Quaternion start, Quaternion end, Scalar t) -> Quaternion {
    Vector<4> vec = normalise(start.vec * (1. - t) + end.vec * t);
    return Quaternion { vec };
}



inline auto print(Scalar scalar) -> std::string {
    std::string str = std::to_string(scalar);
    if (str.size() > 4) return str.substr(0, str.size() - 3);
    return str;
}

vec_fn print(Vector<n> vec) -> std::string {
    std::string out = "[ ";
    MAFS_FOR(n-1) out += print(vec[i])   + ", ";
    return        out +  print(vec[n-1]) + " ]";
}

mat_fn print(Matrix<n, m> mat) -> std::string {
    
    std::string out = "[";
    MAFS_FOR(n) {
        out += print(vector_from_column(mat, i));
    }
    return out + "]";
}

inline auto print(Quaternion quat) -> std::string {
    return print(quat.vec);
}

}
