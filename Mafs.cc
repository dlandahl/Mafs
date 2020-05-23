
#include "Mafs.hh"

using namespace mafs;

auto main() -> i32 {

    std::cout << R"(
    //
    // Transform a vector by a matrix
    //
    // [ 1  3  2 ] [-1 ]   [ 9 ]
    // [-2  0  0 ] [ 2 ] = [ 2 ]
    //             [ 2 ]
    //
    // Matrices are mat[c][r], which may be less cache friendly than [r][c],
    // but it's more intuitive in terms of the maths that's happening.
    // It's an experiment, I may change it in the future for cache and data locality niceness.
    //
    )";

    {
        Vector<3>    vec = { -1, 2, 2 };
        Matrix<3, 2> mat = { 1, -2, 3, 0, 2, 0 };

        std::cout << "\tTransformed Result: ";
        std::cout << print(mat * vec) << "\n\n";
    }


    std::cout << R"(
    //
    // Compose two matrices
    //
    // [ 0  2 ] [ 1 -2 ]   [ 2  0 ]
    // [ 1  0 ] [ 1  0 ] = [ 1 -2 ]
    //
    )";

    {
        Matrix<2, 2> a = { 0, 1, 2, 0 };
        Matrix<2, 2> b = { 1, 1,-2, 0 };

        std::cout << "\tComposition Result: ";
        std::cout << print(compose(a, b)) << "\n\n";
    }


    std::cout << R"(
    //
    // Find the determinant of a matrix
    //
    //      [ 1  2 ]
    // det ([ 1 -1 ]) = -3
    //
    )";

    {
        Matrix<2, 2> mat = { 1, 1, 2, -1 };
        std::cout << "\tDeterminant: ";
        std::cout << print(determinant(mat)) << "\n\n";
    }

    std::cout << std::flush;
}
