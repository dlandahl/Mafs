
#include "Mafs.hh"

using namespace mafs;

auto demo() -> void {

    std::cout << R"(
    //
    // Transform a vector by a matrix
    //
    // [ 1  3  2 ] [-1 ]   [ 9 ]
    // [-2  0  0 ] [ 2 ] = [ 2 ]
    //             [ 2 ]
    //
    // Matrices are mat[c][r], which may be less cache friendly than [r][c], but
    // it's more intuitive for me in terms of the maths that's happening.
    // I may fix it in the future for data locality.
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


    std::cout << R"(
    //
    // Go between Vector<3> euler angles and quaternions
    //
    // q ([ pi, 0, -pi/2 ])
    //
    // We use z-up
    //
    )";

    {
        Vector<3> euler_angles = { pi, 0, -pi / 2. };
        Quaternion q = euler_to_quat(euler_angles);

        std::cout << "\tQuaternion from Eulers: ";
        std::cout << print(q.vec) << "\n";
        std::cout << "\tEulers from Quaternion: ";
        std::cout << print(quat_to_euler(q)) << "\n\n";
    }


    std::cout << R"(
    //
    // Interpolate Quaternions using normalised linear interpolation
    //
    // We take two Vector<3>s, convert them to quaternions, nlerp, and convert them back
    // NLERP([ 0, 0, pi], [pi, 0, 0], 0.5)
    //
    )";

    {
        Vector<3> start = { 0, 0, pi };
        Vector<3> end = { pi, 0, 0 };

        std::cout << "\tInterpolated: ";
        std::cout << print(quat_to_euler(nlerp_quaternions(euler_to_quat(start), euler_to_quat(end), 0.5))) << "\n\n";
    }


    std::cout << R"(
    //
    // Get the rotation matrix of a quaternion
    //
    // q = [ 0, 0.707, 0.707, 0 ]
    //     q is xyzw     
    //
    //     [ 0, 1, 0 ]
    // R = [ 1, 0, 0 ]
    //     [ 0, 0,-1 ]
    //
    )";

    {
        Quaternion quat = { 0, 0.707, 0.707, 0 };

        std::cout << "\tMatrix of this quaternion: ";
        std::cout << print(rotation_matrix_of(quat)) << "\n\n";
    }

    std::cout << "\n\nMore things to do:\n";
    std::cout << " - Cross and dot products:\n\t"     << print(cross(Vector<3> { 3, -1, 0.3 }, Vector<3> { -6, 2, 9.2 })) << "\n\n";
    std::cout << " - Generate identity matrices:\n\t" << print(identity_matrix<3, 3>()) << "\n\n";

    std::cout << " - Check linear dependance:\n\t{ 3,-1 } and {-6, 2 } are linearly dependant: "
              << (are_linearly_dependant(Vector<2> { 3, -1 }, Vector<2> { -6, 2 }) ? "true" : "false") << "\n\n";

    std::cout << " - Get magnitude and normalise:\n\t{ 3, 6, 5 } has magnitude: "
              << print(magnitude(Vector<3> {6, 3, 2})) << ". Normed: " << print(normalise(Vector<3> {6, 3, 2})) << "\n\n";

    std::cout << " - Multiply, invert, or conjugate Quaternions:\n\t" << print(Quaternion { 0.3, 0.2, 0.2, 0.6 } *
                                                                               Quaternion { 0.3, 0.2, 0.2, 0.6}) << "\n\n";
}

auto main() -> i32 {
    demo();
}
