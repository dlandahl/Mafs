
#include <cmath>
#include <cassert>
#include <iostream>
#include <iomanip>

namespace mafs {

using i32 = int32_t;
using i64 = int64_t;

using Scalar = float;

#define MAFS_FOR(N)  for (i64 i = 0; i < (N); i++)
#define vec_fn  template <i64 n> auto
#define mat_fn  template <i64 n, i64 m> auto


template <i64 n>
union Vector {
    static_assert(n > 1, "Vector size must be at least 2");

    Scalar data[n] = { 0.f };
    struct { Scalar x, y, z, w; };

    Scalar& operator[](i64 i);
};

// Specializing the vector type allows us to prevent use of
// vector components that aren't in a vector of that size.
// So, access specifying the `z` component of a `Vector<2>` is a compile-time error
#define SPECIALIZE_VECTOR(N, components)  \
template<> union Vector<N> {              \
    Scalar data[N] = { 0.f };             \
    struct { Scalar components ; };       \
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
    Scalar* operator[](i64 i);
};



mat_fn Matrix<n, m>::operator[](i64 i) -> Scalar* {
    return data[i];
}

vec_fn Vector<n>::operator[](i64 i) -> Scalar& {
    return data[i];
}



#define OPERATION_ON_VECTORS(op)                                               \
template <i64 a, i64 b>                                                        \
auto operator op(Vector<a> lhs, Vector<b> rhs) -> Vector<std::max(a, b)> {     \
    const i64 larger  = std::max(a, b);                                        \
    const i64 smaller = std::min(a, b);                                        \
                                                                               \
    Vector<larger> out;                                                        \
                                                                               \
    MAFS_FOR(larger) {                                                         \
        if (i < smaller) out[i] = rhs[i] op lhs[i];                            \
        else out[i] = a > b ? lhs[i] : rhs[i];                                 \
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


vec_fn negate(Vector<n> vec) -> Vector<n> {
    Vector<n> neg = 0.f;
    MAFS_FOR(n) neg = -vec[i];
    return neg;
}

vec_fn dot(Vector<n> lhs, Vector<n> rhs) -> Scalar {
    Scalar sum = 0.f;
    MAFS_FOR(n) sum += lhs[i] * rhs[i];
    return sum;
}

vec_fn magnitude(Vector<n> vec) -> Scalar {
    Scalar abs = 0.f;
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
    MAFS_FOR(m) out[i] = mat[c][i];
    return out;
}



mat_fn transform(Matrix<n, m> mat, Vector<n> vec) -> Vector<m> {
    Vector<m> out = { 0 };
    MAFS_FOR(n) out = out + vector_from_column(mat, i) * vec[i];
    return out;
}

mat_fn operator*(Matrix<n, m> mat, Vector<n> vec) -> Vector<m> {
    return transform(mat, vec);
}




mat_fn transpose(Matrix<n, m> mat) -> Matrix<m, n> {
    Matrix<m, n> out;
    MAFS_FOR(n) {
        for (i64 j = 0; j < m; j++) mat[i][j] = out[j][i];
    }
    return out;
}

vec_fn compose(Matrix<n, n> a, Matrix<n, n> b) -> Matrix<n, n> {
    Matrix<n, n> out;
    MAFS_FOR(n) {
        Vector<n> col = vector_from_column(b, i);
        col =  transform(a, col);
        for (i64 j = 0; j < n; j++) out[i][j] = col[j];
    }
    return out;
}

vec_fn determinant(Matrix<n, n> mat) -> Scalar {
    if (n == 2)  return mat[0][0] *  mat[1][1] - mat[1][0] * mat[0][1];
    if (n == 3)  return mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1])
                      - mat[1][0] * (mat[0][1] * mat[2][2] - mat[0][2] * mat[2][1])
                      + mat[2][0] * (mat[0][1] * mat[1][2] - mat[1][1] * mat[0][2]);
    assert(false && "Determinant of arbitrary square matrices will be added in the future");
    return 0.f;
}


enum class Direction: i64 {
    x, y, z, w
};

vec_fn unit_vector(Direction dir) -> Vector<n> {
    assert((i64) dir < n && "Vector is too small to have a component in this direction");

    Vector<n> unit;
    MAFS_FOR(n) unit[i] = 0.f;
    unit[(i64) dir] = 1.f;

    return unit;
}



auto print(Scalar scalar) -> std::string {
    std::string str = std::to_string(scalar);
    return str.substr(0, str.size() - 3);
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

}
