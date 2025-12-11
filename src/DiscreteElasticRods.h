#ifndef __DISCRETE_ELASTIC_RODS_H__
#define __DISCRETE_ELASTIC_RODS_H__

// std
#include <vector>
// Eigen
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Sparse>
// project
#include "SimParameters.h"
#define M_PI 3.14159265359

class DiscreteElasticRods {
public:
    double t = 0.;
    Eigen::VectorXd x;
    Eigen::VectorXd x_iter;
    Eigen::VectorXd x_delta;
    std::vector<bool> is_fixed;
    std::vector<bool> is_connected;
    std::vector<int> rod_id;
    std::vector<bool> is_collision;
    Eigen::VectorXd v;
    Eigen::VectorXd e;
    Eigen::VectorXd length_rest;
    Eigen::VectorXd l_ref;
    Eigen::VectorXd length;
    Eigen::MatrixX3d kb;
    Eigen::MatrixX2d kappa_ref;
    Eigen::MatrixX2d kappa_rest;
    Eigen::VectorXd twist_rest;
    int nv;

    // time-parallel frame
    Eigen::MatrixX3d d1_ref;
    Eigen::MatrixX3d d2_ref;
    Eigen::MatrixX3d d1;
    Eigen::MatrixX3d d2;
    Eigen::MatrixX3d d1_vis;
    Eigen::MatrixX3d d2_vis;
    Eigen::VectorXd theta;

    // simulation parameters
    SimParameters params;
    bool verbose = false;

    // visualization
    double stretching_energy = 0.;
    double bending_energy = 0.;
    double twisting_energy = 0.;
    std::vector<Eigen::Vector3d> vis_gradient;
    std::vector<Eigen::Vector3d> vis_stretching_force;
    std::vector<Eigen::Vector3d> vis_bending_force;
    std::vector<Eigen::Vector3d> vis_twisting_force;
    std::vector<Eigen::Vector3d> vis_actuation_force;

    // collision
    double weight = 0.;
    Eigen::SparseMatrix<double> md;
    Eigen::SparseMatrix<double> n;
    Eigen::SparseMatrix<double> w;

    // dragging
    double angle = 0.;
    Eigen::Matrix3d A_angle;
    Eigen::Matrix3d B_angle;

    // control
    Eigen::VectorXd u; // muscle_excitation, in [-1, 1]
    Eigen::VectorXd u_new;
    Eigen::Vector3d actuation_force;
    double total_impulse = 0.;
    double muscle_axis_angle = 0.;

    DiscreteElasticRods();

    void initSimulation(int nv_, Eigen::VectorXd x_, Eigen::VectorXd theta_, std::vector<bool> is_fixed_,
            std::vector<bool> is_connected_, std::vector<int> rod_id_, SimParameters params_);

    void simulateOneStep();

    void updateCenterlinePosition(void);

    void updateMaterialFrame(Eigen::MatrixX3d prev_d3, Eigen::MatrixX3d d3);

    std::tuple<Eigen::MatrixXd, Eigen::SparseMatrix<double> > createZeroGradientAndHessian();

    void computeGradientAndHessian(Eigen::VectorXd& gradient,
            Eigen::SparseMatrix<double>& hessian,
            Eigen::MatrixX3d& d3,
            Eigen::VectorXd& twist);

    void updateCenterlineVelocity(Eigen::VectorXd& gradient);

    void updateFrameTheta(Eigen::VectorXd& gradient);

    void updateEdge();

    void updateLength();

    void updateCurvatureBinormal(Eigen::MatrixX3d d3);

    Eigen::MatrixX3d unitTangents(Eigen::VectorXd& x_);

    Eigen::Vector3d parallelTransport(Eigen::Vector3d v, Eigen::Vector3d r1, Eigen::Vector3d r2);

    void getMaterialCurvature(Eigen::MatrixX2d& kappa);

    Eigen::VectorXd getTwist(Eigen::MatrixX3d& d2, Eigen::MatrixX3d& d3);

    void getVoronoiLength(Eigen::VectorXd& l);

    Eigen::Matrix3d getCrossMatrix(Eigen::Vector3d v);

    double applyStretchingForce(Eigen::VectorXd& gradient,
            std::vector<Eigen::Triplet<double>>& hessian,
            Eigen::MatrixX3d& d3,
            Eigen::VectorXd& stretching_force);

    double applyBendingForce(Eigen::VectorXd& gradient,
            std::vector<Eigen::Triplet<double>>& hessian,
            Eigen::MatrixX3d& d3,
            Eigen::VectorXd& bending_force);

    double applyTwistingForce(Eigen::VectorXd& gradient,
            std::vector<Eigen::Triplet<double>>& hessian,
            Eigen::VectorXd& twist,
            Eigen::VectorXd& twisting_force);

    double applyGravity(Eigen::VectorXd& gradient);

    void buildVisualization(Eigen::VectorXd& gradient);

    void buildForceVisualization(Eigen::VectorXd& stretching_force,
            Eigen::VectorXd& bending_force,
            Eigen::VectorXd& twisting_force);

    std::tuple<double, Eigen::Vector3d, Eigen::Vector2d> getMinimumDistance(
            const Eigen::Vector3d& ri1,
            const Eigen::Vector3d& ri2,
            const Eigen::Vector3d& rj1,
            const Eigen::Vector3d& rj2);

    void computeDisplacements();

    double applyCollision(Eigen::VectorXd& gradient);

    double applyDraggingForce(Eigen::VectorXd& gradient);

    double computeCPGSignal();

    void updateExcitation();

    void updateRestShape();

protected:

private:
};

#endif
