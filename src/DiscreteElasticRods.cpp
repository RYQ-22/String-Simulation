#include "DiscreteElasticRods.h"
#include <iostream>
#include <utility>

DiscreteElasticRods::DiscreteElasticRods() { }

void DiscreteElasticRods::initSimulation(int nv_, Eigen::VectorXd x_, Eigen::VectorXd theta_,
        std::vector<bool> is_fixed_, std::vector<bool> is_connected_, SimParameters params_)
{
    nv = nv_;
    x = std::move(x_);
    is_fixed = std::move(is_fixed_);
    is_connected = std::move(is_connected_);
    x_iter.resize(nv*3);
    x_delta.resize(nv*3);
    v.resize(nv*3);
    v.setZero();

    e.resize((nv-1)*3);
    e.setZero();
    length_rest.resize(nv-1);
    length_rest.setZero();
    length.resize(nv-1);
    length.setZero();
    updateEdge();
    updateLength();

    int ne_vis = 0;

    for (int i = 0; i<nv-1; i++) {
        length_rest(i) = length(i);
        if (is_connected[i]) ne_vis++; 
    }
    Eigen::MatrixX3d d3 = unitTangents(x);

    // set up reference frame
    theta = std::move(theta_);
    // theta.resize(nv-1);
    // theta.setZero();
    d1_ref.resize(nv-1, 3);
    d2_ref.resize(nv-1, 3);
    d1_vis.resize(ne_vis, 3);
    d2_vis.resize(ne_vis, 3);

    for (int i = 0, i_vis = 0; i<nv-1; i++) {
        Eigen::Vector3d d3_i = d3.row(i).transpose();
        if (d3_i(2) == 1.) d1_ref.row(i) = Eigen::RowVector3d(1, 0, 0);
        else d1_ref.row(i) = Eigen::RowVector3d(-d3_i(1), d3_i(0), 0);
        Eigen::Vector3d d1_i = d1_ref.row(i).transpose();
        d2_ref.row(i) = (d1_i.cross(d3_i)).transpose();
        if (is_connected[i]) {
            d1_vis.row(i_vis) = d1_ref.row(i);
            d2_vis.row(i_vis) = d2_ref.row(i);
            i_vis++;
        }
    }
    d1.resize(nv-1, 3);
    d1 = d1_ref;
    d2.resize(nv-1, 3);
    d2 = d2_ref;

    kb.resize(nv-2, 3);
    updateCurvatureBinormal(d3);
    getMaterialCurvature(kappa_ref);
    twist_rest.resize(nv-2);
    twist_rest.setZero();
    getVoronoiLength(l_ref);

    params = params_;

    vis_gradient.resize(nv);
    vis_stretching_force.resize(nv);
    vis_bending_force.resize(nv);
    vis_twisting_force.resize(nv);
}

void DiscreteElasticRods::simulateOneStep()
{
    Eigen::MatrixX3d prev_d3 = unitTangents(x);
    updateCenterlinePosition();
    Eigen::MatrixX3d d3 = unitTangents(x);
    updateMaterialFrame(prev_d3, d3);

    updateEdge();
    updateLength();
    updateCurvatureBinormal(d3);
    Eigen::VectorXd twist = getTwist(d2_ref, d3);

    Eigen::VectorXd gradient;
    Eigen::SparseMatrix<double> hessian;
    if (verbose) std::cout << "- Computing Energy" << std::endl;
    computeGradientAndHessian(gradient, hessian, d3, twist);
    updateCenterlineVelocity(gradient);
    updateFrameTheta(gradient);
    buildVisualization(gradient);
}

void DiscreteElasticRods::updateCenterlinePosition(void)
{
    const double h = params.time_step;
    #pragma omp parallel for
    for (int i = 0; i<nv; i++) {
        if (!is_fixed[i]) x.segment<3>(3*i) += h*v.segment<3>(3*i);
    }
}

void DiscreteElasticRods::updateMaterialFrame(Eigen::MatrixX3d prev_d3, Eigen::MatrixX3d d3)
{
    for (int i = 0, i_vis = 0; i<nv-1; i++) {
        double frame_theta = theta(i);
        d1_ref.row(i) = parallelTransport(d1_ref.row(i), prev_d3.row(i).transpose(), d3.row(i).transpose());
        d2_ref.row(i) = parallelTransport(d2_ref.row(i), prev_d3.row(i).transpose(), d3.row(i).transpose());
        d1.row(i) = std::cos(frame_theta)*d1_ref.row(i)+std::sin(frame_theta)*d2_ref.row(i);
        d2.row(i) = -std::sin(frame_theta)*d1_ref.row(i)+std::cos(frame_theta)*d2_ref.row(i);
        if (is_connected[i]) {
            d1_vis.row(i_vis) = d1.row(i);
            d2_vis.row(i_vis) = d2.row(i);
            i_vis++;
        }
    }
}

std::tuple<Eigen::MatrixXd, Eigen::SparseMatrix<double> > DiscreteElasticRods::createZeroGradientAndHessian()
{
    int ndof = 3*nv+nv-1;
    Eigen::VectorXd gradient(ndof);
    gradient.setZero();
    Eigen::SparseMatrix<double> hessian(ndof, ndof);
    return std::make_tuple(gradient, hessian);
}

void DiscreteElasticRods::computeGradientAndHessian(Eigen::VectorXd& gradient,
        Eigen::SparseMatrix<double>& hessian,
        Eigen::MatrixX3d& d3,
        Eigen::VectorXd& twist)
{
    Eigen::VectorXd stretching_force, bending_force, twisting_force;
    stretching_force.resize(3*nv);
    stretching_force.setZero();
    bending_force.resize(3*nv);
    bending_force.setZero();
    twisting_force.resize(3*nv);
    twisting_force.setZero();

    std::tie(gradient, hessian) = createZeroGradientAndHessian();
    std::vector<Eigen::Triplet<double>>
    hessian_triplets;

    if (params.stretching_energy_enabled) {
        if (verbose) std::cout << "    - Computing Stretching Energy" << std::endl;
        stretching_energy = applyStretchingForce(gradient, hessian_triplets, d3, stretching_force);
    }

    if (params.bending_energy_enabled) {
        if (verbose) std::cout << "    - Computing Bending Energy" << std::endl;
        bending_energy = applyBendingForce(gradient, hessian_triplets, d3, bending_force);
    }

    if (params.twisting_energy_enabled) {
        if (verbose) std::cout << "    - Computing Twisting Energy" << std::endl;
        twisting_energy = applyTwistingForce(gradient, hessian_triplets, twist, twisting_force);
    }

    if (params.gravity_enabled) {
        if (verbose) std::cout << "    - Adding Gravity" << std::endl;
        applyGravity(gradient);
    }

    if (params.collision_enabled) {
        if (verbose) std::cout << "    - Applying Collision" << std::endl;
        applyCollision(gradient);
    }

    hessian.setFromTriplets(hessian_triplets.begin(), hessian_triplets.end());

    buildForceVisualization(stretching_force, bending_force, twisting_force);
}

void DiscreteElasticRods::updateCenterlineVelocity(Eigen::VectorXd& gradient)
{
    const double h = params.time_step;
    #pragma omp parallel for
    for (int i = 0; i<nv; i++) {
        // TODO: add mass
        double m = 1;
        if (!is_fixed[i]) v.segment<3>(3*i) -= h*gradient.segment<3>(3*i)/m;
    }

}

void DiscreteElasticRods::updateFrameTheta(Eigen::VectorXd& gradient)
{
    const double h = params.time_step;
    #pragma omp parallel for
    for (int i = 0; i<nv-1; i++) {
        if (!is_fixed[i] || !is_fixed[i+1]) theta(i) -= h*gradient(3*nv+i);
    }
}

void DiscreteElasticRods::updateEdge()
{
    #pragma omp parallel for
    for (int i = 0; i<nv-1; i++) {
        e.segment<3>(3*i) = x.segment<3>(3*(i+1))-x.segment<3>(3*i);
    }
}

void DiscreteElasticRods::updateLength()
{
    #pragma omp parallel for
    for (int i = 0; i<nv-1; i++) {
        length(i) = e.segment<3>(3*i).norm();
    }
}

void DiscreteElasticRods::updateCurvatureBinormal(Eigen::MatrixX3d d3)
{
    #pragma omp parallel for
    for (int i = 0; i<nv-2; i++) {
        Eigen::Vector3d d3_i = d3.row(i).transpose();
        Eigen::Vector3d d3_ip1 = d3.row(i+1).transpose();
        Eigen::Vector3d kb_i = 2*(d3_i).cross(d3_ip1)/(1+d3_i.dot(d3_ip1));
        kb.row(i) = kb_i.transpose();
    }
}

Eigen::MatrixX3d DiscreteElasticRods::unitTangents(Eigen::VectorXd& x_)
{
    Eigen::MatrixX3d unit_tangents(nv-1, 3);
    #pragma omp parallel for
    for (int i = 0; i<nv-1; i++)
        unit_tangents.row(i) = ((x_.segment<3>(3*(i+1))-x_.segment<3>(3*i)).normalized()).transpose();
    return unit_tangents;
}

Eigen::Vector3d DiscreteElasticRods::parallelTransport(Eigen::Vector3d v, Eigen::Vector3d r1, Eigen::Vector3d r2)
{
    Eigen::Vector3d k = r1.cross(r2).normalized();
    double theta = std::atan2((r1.cross(r2)).dot(k), r1.dot(r2));
    return v*std::cos(theta)+(k.cross(v))*std::sin(theta)+k*(k.dot(v))*(1-std::cos(theta));
}

void DiscreteElasticRods::getMaterialCurvature(Eigen::MatrixX2d& kappa)
{
    kappa.resize(nv-2, 2);
    #pragma omp parallel for
    for (int i = 0; i<nv-2; i++) {
        Eigen::VectorXd kb_i = kb.row(i).transpose();
        kappa(i, 0) = kb_i.dot((d2.row(i)+d2.row(i+1)).transpose())/2;
        kappa(i, 1) = -kb_i.dot((d1.row(i)+d1.row(i+1)).transpose())/2;
    }
}

Eigen::VectorXd DiscreteElasticRods::getTwist(Eigen::MatrixX3d& d2, Eigen::MatrixX3d& d3)
{
    Eigen::VectorXd twist(nv-2);
    #pragma omp parallel for
    for (int i = 1; i<nv-1; i++) {
        Eigen::Vector3d vec1 = parallelTransport(d2.row(i-1).transpose(), d3.row(i-1).transpose(),
                d3.row(i).transpose());
        Eigen::Vector3d vec2 = d2.row(i).transpose();
        Eigen::Vector3d n = (vec1.cross(vec2)).normalized();
        double reference_twist = std::atan2((vec1.cross(vec2)).dot(d3.row(i).transpose()), vec1.dot(vec2));
        twist(i-1) = theta(i)-theta(i-1)+reference_twist;
    }
    return twist;
}

void DiscreteElasticRods::getVoronoiLength(Eigen::VectorXd& l)
{
    l.resize(nv-2);
    #pragma omp parallel for
    for (int i = 0; i<nv-2; i++) {
        l(i) = (length_rest(i)+length_rest(i+1))/2;
    }
}

Eigen::Matrix3d DiscreteElasticRods::getCrossMatrix(Eigen::Vector3d v)
{
    Eigen::Matrix3d cross_matrix;
    cross_matrix << 0, -v(2), v(1),
            v(2), 0, -v(0),
            -v(1), v(0), 0;
    return cross_matrix;
}

double DiscreteElasticRods::applyStretchingForce(Eigen::VectorXd& gradient,
        std::vector<Eigen::Triplet<double>>& hessian,
        Eigen::MatrixX3d& d3,
        Eigen::VectorXd& stretching_force)
{
    std::vector<Eigen::Triplet<double>> dF_stretching_triplets;
    Eigen::MatrixX2d kappa;
    getMaterialCurvature(kappa);

    const double r = params.segment_radius;
    const double k = params.stretching_modulus;
    double E_s = 0.0;
    #pragma omp parallel for
    for (int i = 0; i<nv-1; i++) {
        if (!is_connected[i]) continue;
        // ======== compute stretching energy E_s ========
        const double a_i = r;
        const double b_i = r;
        // TODO: fix A_i = pi * a_j * b_j
        const double A_i = M_PI*a_i*b_i;
        E_s += k*A_i*(std::pow(length(i)-length_rest(i), 2)-1)*length(i);

        // ======== compute stretching force derivative ========
        double e_j_rest = length_rest(i);
        double e_j = length(i);
        Eigen::Vector3d t = d3.row(i).transpose();
        gradient.segment<3>(3*i) += k*(e_j/e_j_rest-1)*-t;
        gradient.segment<3>(3*(i+1)) += k*(e_j/e_j_rest-1)*t;
        stretching_force.segment<3>(3*i) -= k*(e_j/e_j_rest-1)*-t;
        stretching_force.segment<3>(3*(i+1)) -= k*(e_j/e_j_rest-1)*t;

        // TODO: check stretching hessian
        // ======== compute bending force hessian ========
    }
    E_s /= 2;

    hessian.insert(hessian.end(), dF_stretching_triplets.begin(), dF_stretching_triplets.end());

    return E_s;
}

double DiscreteElasticRods::applyBendingForce(Eigen::VectorXd& gradient,
        std::vector<Eigen::Triplet<double>>& hessian,
        Eigen::MatrixX3d& d3,
        Eigen::VectorXd& bending_force)
{
    std::vector<Eigen::Triplet<double>> hess_bending_triplets;
    const double r = params.segment_radius;
    const double E = params.bending_modulus;
    double E_b = 0.0;
    Eigen::MatrixX2d kappa;
    getMaterialCurvature(kappa);

    #pragma omp parallel for
    for (int i = 1; i<nv-1; i++) {
        if (!is_connected[i-1] || !is_connected[i]) continue;
        const double a_i = r;
        const double b_i = r;
        const double A_i = M_PI*a_i*b_i;
        const double B11 = E*A_i*pow(a_i, 2)/4;
        const double B22 = E*A_i*pow(b_i, 2)/4;
        const double l_i = l_ref(i-1);
        const double kappa1_ref_i = kappa_ref(i-1, 0);
        const double kappa2_ref_i = kappa_ref(i-1, 1);
        const double kappa1_i = kappa(i-1, 0);
        const double kappa2_i = kappa(i-1, 1);
        // ======== compute bending energy E_b =========
        E_b += (B11*std::pow((kappa1_i-kappa1_ref_i), 2)+
                B22*std::pow((kappa2_i-kappa2_ref_i), 2))/l_i;

        // ======== compute bending force gradient ========
        Eigen::VectorXd gradient_kappa1_i(11); // 3*3 (x DOF) + 2 (theta DOF)
        Eigen::VectorXd gradient_kappa2_i(11); // 3*3 (x DOF) + 2 (theta DOF)
        const double chi = 1+d3.row(i-1).dot(d3.row(i));
        Eigen::Vector3d kb_i = kb.row(i-1).transpose();
        Eigen::Vector3d d1_im1 = d1.row(i-1).transpose();
        Eigen::Vector3d d1_i = d1.row(i).transpose();
        Eigen::Vector3d d1_tilde = (d1_i+d1_im1)/chi;
        Eigen::Vector3d d2_im1 = d2.row(i-1).transpose();
        Eigen::Vector3d d2_i = d2.row(i).transpose();
        Eigen::Vector3d d2_tilde = (d2_i+d2_im1)/chi;
        Eigen::Vector3d t_im1 = d3.row(i-1).transpose();
        Eigen::Vector3d t_i = d3.row(i).transpose();
        Eigen::Vector3d t_tilde = (t_i+t_im1)/chi;

        Eigen::Matrix3d de_dx_im1 = -Eigen::Matrix3d::Identity();
        Eigen::Matrix3d de_dx_i = Eigen::Matrix3d::Identity();

        Eigen::Vector3d dkappa1_i_de_im1 = d2_tilde.cross(-t_i/length_rest(i-1))-kappa1_i*t_tilde/length_rest(i-1);
        Eigen::Vector3d dkappa1_i_de_i = d2_tilde.cross(t_im1/length_rest(i))-kappa1_i*t_tilde/length_rest(i);
        Eigen::Vector3d dkappa2_i_de_im1 = d1_tilde.cross(t_i/length_rest(i-1))-kappa2_i*t_tilde/length_rest(i-1);
        Eigen::Vector3d dkappa2_i_de_i = d1_tilde.cross(-t_im1/length_rest(i))-kappa2_i*t_tilde/length_rest(i);
        gradient_kappa1_i.segment<3>(0) = dkappa1_i_de_im1.transpose()*de_dx_im1;
        gradient_kappa1_i.segment<3>(3) = dkappa1_i_de_im1.transpose()*de_dx_i+dkappa1_i_de_i.transpose()*de_dx_im1;
        gradient_kappa1_i.segment<3>(6) = dkappa1_i_de_i.transpose()*de_dx_i;
        gradient_kappa2_i.segment<3>(0) = dkappa2_i_de_im1.transpose()*de_dx_im1;
        gradient_kappa2_i.segment<3>(3) = dkappa2_i_de_im1.transpose()*de_dx_i+dkappa2_i_de_i.transpose()*de_dx_im1;
        gradient_kappa2_i.segment<3>(6) = dkappa2_i_de_i.transpose()*de_dx_i;

        gradient_kappa1_i(9) = kb_i.dot(d1_im1)/2;
        gradient_kappa1_i(10) = kb_i.dot(d1_i)/2;
        gradient_kappa2_i(9) = -kb_i.dot(d2_im1)/2;
        gradient_kappa2_i(10) = -kb_i.dot(d2_i)/2;

        // update dx
        Eigen::Matrix<double, 9, 1> gradient_dx = (B11*(kappa1_i-kappa1_ref_i)*gradient_kappa1_i.segment<9>(0)
                +B22*(kappa2_i-kappa2_ref_i)*gradient_kappa2_i.segment<9>(0))/l_i;
        gradient.segment<9>(3*(i-1)) += gradient_dx;
        bending_force.segment<9>(3*(i-1)) -= gradient_dx.segment<9>(0);

        // update dtheta
        gradient.segment<2>(3*nv+i-1) += (B11*(kappa1_i-kappa1_ref_i)*gradient_kappa1_i.segment<2>(9)
                +B22*(kappa2_i-kappa2_ref_i)*gradient_kappa2_i.segment<2>(9))/l_i;

        // TODO: compute bending hessian
        // ======== compute bending force hessian ========
    }
    E_b /= 2;

    hessian.insert(hessian.end(), hess_bending_triplets.begin(), hess_bending_triplets.end());

    return E_b;
}

double DiscreteElasticRods::applyTwistingForce(Eigen::VectorXd& gradient,
        std::vector<Eigen::Triplet<double>>& hessian,
        Eigen::VectorXd& twist,
        Eigen::VectorXd& twisting_force)
{
    std::vector<Eigen::Triplet<double>> hess_twisting_triplets;
    const double r = params.segment_radius;
    const double G = params.twisting_modulus;
    double E_t = 0.0;
    Eigen::MatrixX2d kappa;
    getMaterialCurvature(kappa);
    #pragma omp parallel for
    for (int i = 1; i<nv-1; i++) {
        if (!is_connected[i-1] || !is_connected[i]) continue;
        const double a_i = r;
        const double b_i = r;
        //TODO: fix A_i = pi * a_j * b_j
        const double A_i = M_PI*a_i*b_i;
        double beta_i = G*A_i*(std::pow(a_i, 2)+std::pow(b_i, 2))/4;
        const double l_i = l_ref(i-1);
        const double m_i = twist(i-1);
        const double m_i_rest = twist_rest(i-1);
        // ======== compute twisting energy E_t ========
        E_t += beta_i/l_i*std::pow((m_i-m_i_rest), 2);

        // ======== compute twisting force gradient ========
        Eigen::VectorXd gradient_m_i(11); // 3*3 (x DOF) + 2 (theta DOF)

        Eigen::Vector3d kb_i = kb.row(i-1).transpose();
        double e_im1 = length(i-1);
        double e_i = length(i);
        // update dx
        gradient_m_i.segment<3>(0) = -kb_i/(2*e_im1);
        gradient_m_i.segment<3>(3) = -kb_i/(2*e_i)+kb_i/(2*e_im1);
        gradient_m_i.segment<3>(6) = kb_i/(2*e_i);
        // update dtheta
        gradient_m_i(9) = -1.;
        gradient_m_i(10) = 1.;

        Eigen::Matrix<double, 9, 1> gradient_dx = beta_i/l_i*(m_i-m_i_rest)*gradient_m_i.segment<9>(0);
        gradient.segment<9>(3*(i-1)) += gradient_dx;
        twisting_force.segment<9>(3*(i-1)) -= gradient_dx;
        gradient.segment<2>(3*nv+i-1) += beta_i/l_i*(m_i-m_i_rest)*gradient_m_i.segment<2>(9);

        // TODO: compute twisting Hessian dF
        // ======== compute bending force hessian ========
    }
    E_t /= 2;

    hessian.insert(hessian.end(), hess_twisting_triplets.begin(), hess_twisting_triplets.end());

    return E_t;
}

double DiscreteElasticRods::applyGravity(Eigen::VectorXd& gradient)
{
    const double g = params.gravity_G;
    #pragma omp parallel for
    for (int i = 0; i<nv; i++) {
        // TODO: add mass
        double m = 1;
        gradient(3*i+1) -= m*g;
    }
    return 0;
}

void DiscreteElasticRods::buildVisualization(Eigen::VectorXd& gradient)
{
    #pragma omp parallel for
    for (int i = 0; i<nv; i++) {
        vis_gradient[i] = gradient.segment<3>(3*i);
    }
}

void DiscreteElasticRods::buildForceVisualization(Eigen::VectorXd& stretching_force,
        Eigen::VectorXd& bending_force,
        Eigen::VectorXd& twisting_force)
{
    #pragma omp parallel for
    for (int i = 0; i<nv; i++) {
        vis_stretching_force[i] = stretching_force.segment<3>(3*i);
        vis_bending_force[i] = bending_force.segment<3>(3*i);
        vis_twisting_force[i] = twisting_force.segment<3>(3*i);
    }
}

std::tuple<double, Eigen::Vector3d, Eigen::Vector2d> DiscreteElasticRods::getMinimumDistance(
        const Eigen::Vector3d& ri1,
        const Eigen::Vector3d& ri2,
        const Eigen::Vector3d& rj1,
        const Eigen::Vector3d& rj2)
{
    Eigen::Vector3d ei = ri2 - ri1;
    Eigen::Vector3d ej = rj2 - rj1;
    Eigen::Vector3d w = ri1 - rj1;

    const double a1 = ei.dot(ei);
    const double a2 = ei.dot(ej);
    const double a3 = ej.dot(ej);
    const double a4 = ei.dot(w);
    const double a5 = ej.dot(w);
    const double det = a1*a3-a2*a2;
    double s, t;

    if (det > 0.) { // nonparallel segmets
        const double b1 = a2*a5;
        const double b2 = a3*a4;
        if (b1 <= b2) { // s <= 0
            if (a5 <= 0.) { // t <= 0 (region 6)
                s = -a4>=a1 ? 1. : (-a4>=0. ? -a4/a1 : 0.);
                t = 0.;
            }
            else if (a5 < a3) { // 0 < t < 1 (region 5)
                s = 0.;
                t = a5/a3;
            }
            else { // t >= 1 (region 4)
                s = a2-a4>=a1 ? 1. : (a2-a4>=0. ? (a2-a4)/a1 : 0.);
                t = 1.;
            }
        }
        else { // s > 0
            s = b1-b2;
            if (s >= det) { // s >= 1
                if (a2+a5 <= 0.) { // t <= 0 (region 8)
                    s = -a4<=0. ? 0. : (-a4<a1 ? -a4/a1 : 1.);
                    t = 0.;
                }
                else if (a2+a5 < a3) { // 0 < t < 1 (region 1)
                    s = 1.;
                    t = (a2+a5)/a3;
                }
                else { // t >= 1(region 2)
                    s = a2-a4<=0. ? 0. : (a2-a4<a1 ? (a2-a4)/a1 : 1.);
                    t = 1.;
                }
            }
            else { // 0 < s < 1
                const double b3 = a1*a5;
                const double b4 = a2*a4;
                if (b3 <= b4) { // t <= 0 (region 7)
                    s = -a4<=0. ? 0. : (-a4>=a1 ? 1. : -a4/a1);
                    t = 0.;
                }
                else { // t > 0
                    t = b3-b4;
                    if (t >= det) { // t >= 1 (region 3)
                        s = a2-a4<=0. ? 0. : (a2-a4>=a1 ? 1. : (a2-a4)/a1);
                        t = 1.;
                    }
                    else { // 0 < t < 1 (region 0)
                        s /= det;
                        t /= det;
                    }
                }
            }
        }
    }
    else { // parallel segments
        if (a5 <= 0.) {
            s = -a4<=0. ? 0.: (-a4>=a1 ? 1. : -a4/a1);
            t = 0.;
        }
        else if (a5 >= a3) {
            s = a2-a4<=0. ? 0.: (a2-a4>=a1 ? 1. : (a2-a4)/a1);
            t = 1.;
        }
        else {
            s = 0.;
            t = a5/a3;
        }
    }

    Eigen::Vector2d wij(1.-s, 1.-t);
    Eigen::Vector3d nij = (rj1+t*(rj2-rj1))-(ri1+s*(ri2-ri1));
    double mdij = nij.norm();

    return std::make_tuple(mdij, nij, wij);
}

void DiscreteElasticRods::computeDisplacements()
{
    const double d = 2.*params.segment_radius;
    Eigen::Vector3d ri1, ri2, rj1, rj2;
    double mdij;
    Eigen::Vector2d wij;
    Eigen::Vector3d nij;
    std::vector<Eigen::Triplet<double>> md_triplets;
    std::vector<Eigen::Triplet<double>> n_triplets;
    std::vector<Eigen::Triplet<double>> w_triplets;
    md.resize(nv-1,nv-1);
    n.resize(nv-1,3*(nv-1));
    w.resize(nv-1,2*(nv-1));
    for (int i = 0; i<nv-1; i++) {
        if (!is_connected[i]) continue;
        for (int j = i+2; j<nv-1; j++) {
            if (!is_connected[j]) continue;
            ri1 = x_iter.segment<3>(3*i);
            ri2 = x_iter.segment<3>(3*(i+1));
            rj1 = x_iter.segment<3>(3*j);
            rj2 = x_iter.segment<3>(3*(j+1));
            std::tie(mdij, nij, wij) = getMinimumDistance(ri1, ri2, rj1, rj2);
            if (mdij<d) {
                md_triplets.push_back(Eigen::Triplet<double>(i,j,mdij));
                n_triplets.push_back(Eigen::Triplet<double>(i,3*j,nij(0)));
                n_triplets.push_back(Eigen::Triplet<double>(i,3*j+1,nij(1)));
                n_triplets.push_back(Eigen::Triplet<double>(i,3*j+2,nij(2)));
                w_triplets.push_back(Eigen::Triplet<double>(i,2*j,wij(0)));
                w_triplets.push_back(Eigen::Triplet<double>(i,2*j+1,wij(1)));
            }
        }
    } 
    md.setFromTriplets(md_triplets.begin(),md_triplets.end());
    n.setFromTriplets(n_triplets.begin(),n_triplets.end());
    w.setFromTriplets(w_triplets.begin(),w_triplets.end());
}

double DiscreteElasticRods::applyCollision(Eigen::VectorXd& gradient)
{
    const double d = 2.*params.segment_radius;
    const double h = params.time_step;
    const int iters = params.collision_max_iters;
    double mdij;
    Eigen::Vector2d wij;
    Eigen::Vector3d nij;
    x_iter = x;
    x_delta.setZero();

    // solve x_delta iteratively
    double err = -1.;
    for (int iter = 0; iter<iters; iter++) {
        computeDisplacements();
        // if (err!=-1.) std::cout << "err: " << err << std::endl;
        for (int k = 0; k < md.outerSize(); k++) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(md,k); it; ++it) {
                const int i = it.row();
                const int j = it.col();
                // std::cout << "(i,j): (" << i << "," << j << ")" << std::endl; 
                mdij = md.coeff(i,j);
                nij(0) = n.coeff(i,3*j);
                nij(1) = n.coeff(i,3*j+1);
                nij(2) = n.coeff(i,3*j+2);
                wij(0) = w.coeff(i,2*j);
                wij(1) = w.coeff(i,2*j+1);
                err = std::max(err,std::max(d-mdij,0.));
                // std::cout << "nij: " << nij(0) << ", " << nij(1) << ", " << nij(2) << std::endl;
                // std::cout << "mdij: " << mdij << std::endl;
                // std::cout << "wij: " << wij(0) << ", " << wij(1) << "\n" << std::endl;
                // m = 1 for all points
                weight = 0.5;
                x_delta.segment<3>(3*i) += nij*weight*(mdij-d)*wij(0);
                x_delta.segment<3>(3*(i+1)) += nij*weight*(mdij-d)*(1.-wij(0));
                x_delta.segment<3>(3*j) += nij*(1.-weight)*(d-mdij)*wij(1);
                x_delta.segment<3>(3*(j+1)) += nij*(1.-weight)*(d-mdij)*(1.-wij(1));

                x_iter.segment<3>(3*i) = x.segment<3>(3*i) + x_delta.segment<3>(3*i);
                x_iter.segment<3>(3*(i+1)) = x.segment<3>(3*(i+1)) + x_delta.segment<3>(3*(i+1));
                x_iter.segment<3>(3*j) = x.segment<3>(3*j) + x_delta.segment<3>(3*j);
                x_iter.segment<3>(3*(j+1)) = x.segment<3>(3*(j+1)) + x_delta.segment<3>(3*(j+1));
            }
        }
    }

    #pragma omp parallel for
    for (int i = 0; i<nv; i++) {
        // m = 1 for all points
        const double m = 1.;
        gradient.segment<3>(3*i) += -x_delta.segment<3>(3*i)*m/(h*h);
        // std::cout << (x_delta.segment<3>(3*i)*m/(h*h)).norm() << std::endl;
    }

    return 0;
}
