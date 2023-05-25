/**
 * Copyright (C) 2015 by Liangliang Nan (liangliang.nan@gmail.com)
 * https://3d.bk.tudelft.nl/liangliang/
 *
 * This file is part of Easy3D. If it is useful in your research/work,
 * I would be grateful if you show your appreciation by citing it:
 * ------------------------------------------------------------------
 *      Liangliang Nan.
 *      Easy3D: a lightweight, easy-to-use, and efficient C++
 *      library for processing and rendering 3D data. 2018.
 * ------------------------------------------------------------------
 * Easy3D is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License Version 3
 * as published by the Free Software Foundation.
 *
 * Easy3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "triangulation.h"
#include "matrix_algo.h"
#include <easy3d/optimizer/optimizer_lm.h>
#include <iostream>
#include <fstream>


using namespace easy3d;


/**
 * TODO: Finish this function for reconstructing 3D geometry from corresponding image points.
 * @return True on success, otherwise false. On success, the reconstructed 3D points must be written to 'points_3d'
 *      and the recovered relative pose must be written to R and t.
 */
double distance(Vector2D pt, Vector2D centroid){
    return sqrt(pow(pt[0] - centroid[0], 2) + pow(pt[1] - centroid[1], 2));
}

bool Triangulation::triangulation(
        double fx, double fy,     /// input: the focal lengths (same for both cameras)
        double cx, double cy,     /// input: the principal point (same for both cameras)
        const std::vector<Vector2D> &points_0,  /// input: 2D image points in the 1st image.
        const std::vector<Vector2D> &points_1,  /// input: 2D image points in the 2nd image.
        std::vector<Vector3D> &points_3d,       /// output: reconstructed 3D points
        Matrix33 &R,   /// output: 3 by 3 matrix, which is the recovered rotation of the 2nd camera
        Vector3D &t    /// output: 3D vector, which is the recovered translation of the 2nd camera
) const
{
    /// NOTE: there might be multiple workflows for reconstructing 3D geometry from corresponding image points.
    ///       This assignment uses the commonly used one explained in our lecture.
    ///       It is advised to define a function for the sub-tasks. This way you have a clean and well-structured
    ///       implementation, which also makes testing and debugging easier. You can put your other functions above
    ///       triangulation(), or put them in one or multiple separate files.

    std::cout << "\nTODO: I am going to implement the triangulation() function in the following file:" << std::endl
              << "\t    - triangulation_method.cpp\n\n";

    std::cout << "[Liangliang]:\n"
                 "\tFeel free to use any provided data structures and functions. For your convenience, the\n"
                 "\tfollowing three files implement basic linear algebra data structures and operations:\n"
                 "\t    - Triangulation/matrix.h  Matrices of arbitrary dimensions and related functions.\n"
                 "\t    - Triangulation/vector.h  Vectors of arbitrary dimensions and related functions.\n"
                 "\t    - Triangulation/matrix_algo.h  Determinant, inverse, SVD, linear least-squares...\n"
                 "\tPlease refer to the above files for a complete list of useful functions and their usage.\n\n"
                 "\tIf you choose to implement the non-linear method for triangulation (optional task). Please\n"
                 "\trefer to 'Tutorial_NonlinearLeastSquares/main.cpp' for an example and some explanations.\n\n"
                 "\tIn your final submission, please\n"
                 "\t    - delete ALL unrelated test or debug code and avoid unnecessary output.\n"
                 "\t    - include all the source code (and please do NOT modify the structure of the directories).\n"
                 "\t    - do NOT include the 'build' directory (which contains the intermediate files in a build step).\n"
                 "\t    - make sure your code compiles and can reproduce your results without ANY modification.\n\n" << std::flush;

    /// Below are a few examples showing some useful data structures and APIs.

    /// define a 2D vector/point
    Vector2D b(1.1, 2.2);

    /// define a 3D vector/point
    Vector3D a(1.1, 2.2, 3.3);

    /// get the Cartesian coordinates of a (a is treated as Homogeneous coordinates)
    Vector2D p = a.cartesian();

    /// get the Homogeneous coordinates of p
    Vector3D q = p.homogeneous();

    /// define a 3 by 3 matrix (and all elements initialized to 0.0)
    Matrix33 A;

    /// define and initialize a 3 by 3 matrix
    Matrix33 T(1.1, 2.2, 3.3,
               0, 2.2, 3.3,
               0, 0, 1);

    /// define and initialize a 3 by 4 matrix
    Matrix34 M(1.1, 2.2, 3.3, 0,
               0, 2.2, 3.3, 1,
               0, 0, 1, 1);

    /// set first row by a vector
    M.set_row(0, Vector4D(1.1, 2.2, 3.3, 4.4));

    /// set second column by a vector
    M.set_column(1, Vector3D(5.5, 5.5, 5.5));

    /// define a 15 by 9 matrix (and all elements initialized to 0.0)
    Matrix W(15, 9, 0.0);
    /// set the first row by a 9-dimensional vector
    W.set_row(0, {0, 1, 2, 3, 4, 5, 6, 7, 8}); // {....} is equivalent to a std::vector<double>

    /// get the number of rows.
    int num_rows = W.rows();

    /// get the number of columns.
    int num_cols = W.cols();

    /// get the the element at row 1 and column 2
    double value = W(1, 2);

    /// get the last column of a matrix
    Vector last_column = W.get_column(W.cols() - 1);

    /// define a 3 by 3 identity matrix
    Matrix33 I = Matrix::identity(3, 3, 1.0);

    /// matrix-vector product
    Vector3D v = M * Vector4D(1, 2, 3, 4); // M is 3 by 4

    ///For more functions of Matrix and Vector, please refer to 'matrix.h' and 'vector.h'

    // TODO: delete all above example code in your final submission

    //--------------------------------------------------------------------------------------------------------------
    // implementation starts ...

    // TODO: check if the input is valid (always good because you never known how others will call your function).

    // TODO: Estimate relative pose of two views. This can be subdivided into
    //      - estimate the fundamental matrix F;
    //      - compute the essential matrix E;
    //      - recover rotation R and t.
    //NORMALIZE
    //compute centroid of points
    //centroid points_0
    double allx_0 = 0.0;
    double ally_0 = 0.0;
    for(auto &pt : points_0){
        allx_0 = allx_0 + pt[0];
        ally_0 = ally_0 + pt[1];
    }
    Vector2D centroid0(allx_0 / points_0.size(), ally_0 / points_0.size());
    //centroid points_1
    double allx_1 = 0.0;
    double ally_1 = 0.0;
    for(auto &pt : points_1){
        allx_1 = allx_1 + pt[0];
        ally_1 = ally_1 + pt[1];
    }
    Vector2D centroid1(allx_1 / points_1.size(), ally_1 / points_1.size());
    std::cout<< "centroids " << centroid0 << " & " << centroid1 << std::endl;

    //compute scaling factor 0
    double distances_0 = 0.0;
    for (auto &pt: points_0) {
        double dis = distance(pt, centroid0);
        distances_0 = distances_0 + dis;
    }
    double mean_dis0 = distances_0 / points_0.size();
    std::cout << " mean dist  " << mean_dis0 << std::endl;

    double scaling_fact0 = sqrt(2) / mean_dis0;
    std::cout << " total distances 0 is " << distances_0 << std::endl;

    //compute scaling factor 1
    double distances_1 = 0.0;
    for (auto &pt: points_1) {
        double dis = distance(pt, centroid1);
        distances_1 = distances_1 + dis;
    }
    double mean_dis1 = distances_1 / points_1.size();
    double scaling_fact1 = sqrt(2) / mean_dis1;
    std::cout << "scaling factors " << scaling_fact0 << " and " <<scaling_fact1 << std::endl;

    //make translation matrices
    Matrix33 T0 (scaling_fact0, 0, centroid0[0], 0, scaling_fact0, centroid0[1], 0, 0, 1 );
    Matrix33 T1 (scaling_fact1, 0, centroid1[0], 0, scaling_fact1, centroid1[1], 0, 0, 1 );
    std::cout << "translation/scaling matrices " << T0 << " and " <<T1 << std::endl;

 // apply translation/scaling: T * p. Gives the normalized points.
    std::vector<Vector3D> norm_points_0;
    for (auto &pt : points_0) {
//        std::cout << " point " << pt << std::endl;
        Vector3D homo_pt(pt[0], pt[1], 1);
//        std::cout << "homo point " << homo_pt << std::endl;
        Vector3D translated_point0 = mult(T0 , homo_pt);
//        std::cout << "norm point " << translated_point0 << std::endl;
        norm_points_0.push_back(translated_point0);
    }
//    std::cout << " norm points " << norm_points_0[0] << " and " << norm_points_0[95] << " and " << norm_points_0[72] << std::endl;



    std::vector<Vector3D> norm_points_1;
    for (auto &pt : points_1) {
        Vector3D homo_pt(pt[0], pt[1], 1);
        Vector3D translated_point1 = mult(T1, homo_pt);
        norm_points_1.push_back(translated_point1);
    }

//    std::cout << " 1 normalized point is " << norm_points_0[15] << " which before was " << points_0[15] << std::endl;
//    std::cout << " another normalized point is " << norm_points_0[34] << " which before was " << points_0[34] << std::endl;



    //CONSTRUCTING THE W MATRIX
    int nrrows = points_0.size();
    Matrix W_matrix(nrrows, 9, 0.0) ;
    for (int i = 0; i < points_0.size(); ++i) {
        W_matrix.set_row(i, {norm_points_0[i][0] * norm_points_1[i][0], norm_points_0[i][1]*norm_points_1[i][0], norm_points_1[i][0], norm_points_0[i][0]*norm_points_1[i][1], norm_points_0[i][1]*norm_points_1[i][1], norm_points_1[i][1], norm_points_0[i][0], norm_points_0[i][1], 1});
    }
//    std::cout << " W matrix: \n" << W_matrix <<std::endl;

    //SVD DECOMPOSE W
    Matrix U_matrix(points_0.size(), points_0.size(), 0.0);   // initialized with 0s
    Matrix D_matrix(points_0.size(), 9, 0.0);   // initialized with 0s
    Matrix V_matrix(9, 9, 0.0);   // initialized with 0s
    svd_decompose(W_matrix, U_matrix, D_matrix, V_matrix);
//    std::cout << " U " << U_matrix << std::endl;
//    std::cout << " D " << D_matrix << std::endl;
//    std::cout << " V " << V_matrix << std::endl; //is this transpose? NO! W = U * D * V^T

    //last column of V is F hat vec
    Vector F_hat_vec = V_matrix.get_column(V_matrix.cols() - 1);
    std::cout << "F hat vec " << F_hat_vec << std::endl;
    //make F hat matrix
    Matrix33 F_hat_mat(F_hat_vec[0], F_hat_vec[1], F_hat_vec[2], F_hat_vec[3], F_hat_vec[4], F_hat_vec[5], F_hat_vec[6], F_hat_vec[7], F_hat_vec[8]);

    //decompose F hat
    Matrix U(3, 3, 3);
    Matrix D(3, 3, 3);
    Matrix V(3, 3, 3);
    svd_decompose(F_hat_mat, U, D, V);
    std::cout << " U " << U << std::endl;
    std::cout << " D " << D << std::endl;
    std::cout << " V " << V<< std::endl;

    //get F
    std::cout << " D " << D << std::endl;
    double d1 = D.get(0,0);
    double d2 = D.get(1,1);
    Matrix33 D_r2(d1, 0.0, 0.0, 0.0, d2, 0.0, 0.0, 0.0, 0.0);
    std::cout << " D r2" << D_r2 << std::endl;
    Matrix33 F = mult(mult(U, D_r2), transpose(V));
    std::cout << " F " << F << std::endl;

    // TODO: Reconstruct 3D points. The main task is
    //      - triangulate a pair of image points (i.e., compute the 3D coordinates for each corresponding point pair)

    // TODO: Don't forget to
    //          - write your recovered 3D points into 'points_3d' (so the viewer can visualize the 3D points for you);
    //          - write the recovered relative pose into R and t (the view will be updated as seen from the 2nd camera,
    //            which can help you check if R and t are correct).
    //       You must return either 'true' or 'false' to indicate whether the triangulation was successful (so the
    //       viewer will be notified to visualize the 3D points and update the view).
    //       There are a few cases you should return 'false' instead, for example:
    //          - function not implemented yet;
    //          - input not valid (e.g., not enough points, point numbers don't match);
    //          - encountered failure in any step.
    return points_3d.size() > 0;

}