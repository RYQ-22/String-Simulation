#ifndef SIMPARAMETERS_H
#define SIMPARAMETERS_H
#include <cstddef>

struct NerveParameters {
    float theta1; // duration
    float theta2; // amplitude
    float theta3; // period
    float alpha; // delivery factor
    float beta;  // attenuation factor
    int kappa; // propagation speed
};

struct SimParameters {
    double time_step;

    int newton_max_iters;
    double newton_tolerance;

    bool gravity_enabled;
    double gravity_G;

    double stretching_modulus;
    double bending_modulus;
    double twisting_modulus;
    double segment_radius;

    bool stretching_energy_enabled;
    bool bending_energy_enabled;
    bool twisting_energy_enabled;

    // collision
    int collision_max_iters;
    bool collision_enabled;

    // dragging
    bool dragging_enabled;
    float water_velocity;
    double volume_fraction;
    double d_parallel, d_perp;

    // control
    bool control_enabled;
    double muscle_max_angle;
    NerveParameters nerve_params;

  SimParameters()
    {
        time_step = 5e-4;

        newton_max_iters = 20;
        newton_tolerance = 1e-8;

        gravity_enabled = true;
        gravity_G = -9.8;

        stretching_modulus = 1e4;
        bending_modulus = 1.e7;
        twisting_modulus = 1e5;
        segment_radius = 0.12;

        stretching_energy_enabled = true;
        bending_energy_enabled = true;
        twisting_energy_enabled = true;

        collision_max_iters = 50;
        collision_enabled = false;

        dragging_enabled = true;
        water_velocity = 0.0f;

        control_enabled = false;
        muscle_max_angle = 0.2;
        nerve_params.theta1 = 2.0f;
        nerve_params.theta2 = 1.f;
        nerve_params.theta3 = 4.0f;
        nerve_params.alpha = 0.3f;
        nerve_params.beta = 0.7f;
        nerve_params.kappa = 1;
    }
};

#endif
