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
//    Matrix33 A;

    /// define and initialize a 3 by 3 matrix
    Matrix33 T(1.1, 2.2, 3.3,
               0, 2.2, 3.3,
               0, 0, 1);

//    /// define and initialize a 3 by 4 matrix
//    Matrix34 M(1.1, 2.2, 3.3, 0,
//               0, 2.2, 3.3, 1,
//               0, 0, 1, 1);

    /// set first row by a vector
//    M.set_row(0, Vector4D(1.1, 2.2, 3.3, 4.4));

    /// set second column by a vector
//    M.set_column(1, Vector3D(5.5, 5.5, 5.5));

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
//    Vector3D v = M * Vector4D(1, 2, 3, 4); // M is 3 by 4

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

    //scaling factor 0
    double dis0x = 0.0;
    double dis0y = 0.0;
    for (auto &pt : points_0) {
        double disx = abs(pt.x() - centroid0.x());
        dis0x += disx;
        double disy = abs(pt.y() - centroid0.y());
        dis0y += disy;
    }
    double mean_0x = dis0x / points_0.size();
    double mean_0y = dis0y / points_0.size();
    double scaling0x = sqrt(2) / mean_0x;
    double scaling0y = sqrt(2) / mean_0y;

    //scaling factor 1
    double dis1x = 0.0;
    double dis1y = 0.0;
    for (auto &pt : points_1) {
        double disx = abs(pt.x() - centroid1.x());
        dis1x += disx;
        double disy = abs(pt.y() - centroid1.y());
        dis1y += disy;
    }
    double mean_1x = dis1x / points_1.size();
    double mean_1y = dis1y / points_1.size();
    double scaling1x = sqrt(2) / mean_1x;
    double scaling1y = sqrt(2) / mean_1y;


    //make translation/scaling matrices
    Matrix33 T0 (scaling0x, 0, centroid0[0], 0, scaling0y, centroid0[1], 0, 0, 1 );
    Matrix33 T1 (scaling1x, 0, centroid1[0], 0, scaling1y, centroid1[1], 0, 0, 1 );
    std::cout << "translation/scaling matrices " << T0 << " and " <<T1 << std::endl;

    // apply translation/scaling: T * p. Gives the normalized points.
    std::vector<Vector3D> norm_points_0;
    for (auto &pt : points_0) {
        Vector3D homo_pt(pt[0], pt[1], 1);
        Vector3D translated_point0 = mult(T0 , homo_pt);
        norm_points_0.push_back(translated_point0);
    }
    std::vector<Vector3D> norm_points_1;
    for (auto &pt : points_1) {
        Vector3D homo_pt(pt[0], pt[1], 1);
        Vector3D translated_point1 = mult(T1, homo_pt);
        norm_points_1.push_back(translated_point1);
    }


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
//    std::cout << " V " << V_matrix << std::endl; //is this transpose? YA?! noo?  W = U * D * V^T

    //last column of V is F hat vec
    Vector F_hat_vec = transpose(V_matrix).get_column(transpose(V_matrix).cols() - 1);
//    std::cout << "F hat vec " << F_hat_vec << std::endl;
    //make F hat matrix
    Matrix33 F_hat_mat(F_hat_vec[0], F_hat_vec[1], F_hat_vec[2], F_hat_vec[3], F_hat_vec[4], F_hat_vec[5], F_hat_vec[6], F_hat_vec[7], F_hat_vec[8]);
//    std::cout << "f hat mat " << F_hat_mat << std::endl;
    //decompose F hat
    Matrix U(3, 3, 0.0);
    Matrix D(3, 3, 0.0);
    Matrix V(3, 3, 0.0);
    svd_decompose(F_hat_mat, U, D, V);
    std::cout << " U " << U << std::endl;
    std::cout << " D " << D << std::endl;
    std::cout << " V " << V<< std::endl;

    //get F
    double d1 = D.get(0,0);
    double d2 = D.get(1,1);
    Matrix33 D_r2(d1, 0.0, 0.0, 0.0, d2, 0.0, 0.0, 0.0, 0.0);
    std::cout << " D r2" << D_r2 << std::endl;
    Matrix33 F_q = mult(mult(U, D_r2), V);
    std::cout << " F " << F_q << std::endl;

    //denormalize F
    Matrix33 F = mult(mult(transpose(T1), F_q), T0);
    std::cout << "Denormalized F " << F << std::endl;


    
    
    // TODO: Estimate relative pose of two views. This can be subdivided into
    //      - estimate the fundamental matrix F;
    //      - compute the essential matrix E;
    //      - recover rotation R and t.

    Matrix33 K (fx, 0, cx,
                0, fy, cy,
                0, 0, 1);

    //Define matrix E
    Matrix E = K.transpose() * F * K;
    
    Matrix E_U (3, 3, 0.0);
    Matrix E_D (3, 3, 0.0);
    Matrix E_V (3, 3, 0.0);

    svd_decompose(E, E_U, E_D, E_V);

    //Obtain t (last column of U)
    Vector t1 = E_U.get_column(E_U.cols()-1);
    Vector t2 = -1 * E_U.get_column()-1);

    //Compute rotation matrix R
    Matrix33 E_W (0, -1, 0,
                  1, 0, 0,
                  0, 0, 1);

    double determ1 = determinant(E_U * E_W * E_V.transpose());
    double determ2 = determinant(E_U * E_W.transpose() * E_V.transpose());

    Matrix R1 = determ1 * E_U * E_W * E_V.transpose();
    Matrix R2 = determ2 * E_U * E_W.transpose() * E_V.transpose();

    int calcNumberBehindOrigin(Matrix33 K, Matrix33 R, Vector3D t, const std::vector<Vector2D> &points_0, const std::vector<Vector2D> &points_1);
    //The four candidate relative poses
    int pose0 = calcNumberBehindOrigin(K, R1, t1, points_0, points_1);
    int pose1 = calcNumberBehindOrigin(K, R1, t2, points_0, points_1);
    int pose2 = calcNumberBehindOrigin(K, R2, t1, points_0, points_1);
    int pose3 = calcNumberBehindOrigin(K, R2, t2, points_0, points_1);

    if(pose0 < pose1 && pose0 < pose2 && pose0 < pose3) {
        std::cout << "relative pose 0" << std::endl;
        std::cout << determinant(R1) << std::endl;
        t = t1;
        R = R1;
        
    } else if (pose1 < pose2 && pose1 < pose3) {
        std::cout << "relative pose 1" << std::endl;
        std::cout << determinant(R1) << std::endl;
        t = t2;
        R = R1;
        
    } else if(pose2 < pose3){
        std::cout << "relative pose 2" << std::endl;
        std::cout << determinant(R2) << std::endl;
        t = t1;
        R = R2;
        
    } else {
        std::cout << "relative pose 3" << std::endl;
        std::cout << determinant(R2) << std::endl;
        t = t2;
        R = R2;
        
    }



    
    
    // TODO: Reconstruct 3D points. The main task is
    //      - triangulate a pair of image points (i.e., compute the 3D coordinates for each corresponding point pair)

    // Compute PROJECTION MATRIX (M_matrix)
    Matrix33 K_prime, K;
    Matrix34 RT;
    Matrix34 I0 (1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0);
    RT.set_row(0, {R(0,0), R(0, 1), R(0,2), t[0]});
    RT.set_row(1, {R(1,0), R(1, 1), R(1,2), t[1]});
    RT.set_row(2, {R(2,0), R(2, 1), R(2,2), t[2]});
    Matrix34 M, M_prime;
    M_prime = mult(K_prime, RT);
    M = mult(K, I0);
    std::cout << "M\n" << M << std::endl;

    int id = 0;
    for (auto &pt : points_0) {
        Vector2D pt_1 = points_1[id];
        // construct A matrix
        Matrix A(4, 4, 0.0);
        std:: cout << "try: " << (pt.x() * M.get_row(2)) - M.get_row(0) << std::endl;
        A.set_row(0, (pt.x() * M.get_row(2)) - M.get_row(0));
        A.set_row(1, (pt.y() * M.get_row(2)) - M.get_row(1));
        A.set_row(2, (pt_1.x() * M_prime.get_row(2)) - M_prime.get_row(0));
        A.set_row(3, (pt_1.y() * M_prime.get_row(2)) - M_prime.get_row(1));

        // get P using SVD?
        Matrix U_mat(A.rows(), A.rows(), 0.0);
        Matrix D_mat(A.rows(), 4, 0.0);
        Matrix V_mat(4, 4, 0.0);
        svd_decompose(A, U_mat, D_mat, V_mat);
        std::cout << "decomposed\n" << std::endl;
        Vector P_vec = V_mat.get_column(V_mat.cols() - 1);
        std::cout << "p vec \n" << P_vec << std::endl;

        // assign 3d point result to points_3d
        points_3d[id].x() = P_vec[0];
        points_3d[id].y() = P_vec[1];
        points_3d[id].z() = P_vec[2];

        id++;
    }


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
