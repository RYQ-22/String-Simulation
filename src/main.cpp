// std
#include <iostream>
#include <omp.h>
// polyscope
#include "polyscope/polyscope.h"
// project
#include "InputParser.h"
#include "DiscreteElasticRodsDriver.h"

DiscreteElasticRodsDriver discrete_elastic_rods_driver;

void discreteElasticRodsDriverCallback()
{
    ImGui::Text("Status: %s", (discrete_elastic_rods_driver.running) ? "Running" : "Stopped");
    ImGui::Text("Stretching Energy = %f", discrete_elastic_rods_driver.discrete_elastic_rods.stretching_energy);
    ImGui::Text("Bending Energy = %f", discrete_elastic_rods_driver.discrete_elastic_rods.bending_energy);
    ImGui::Text("Twisting Energy = %f", discrete_elastic_rods_driver.discrete_elastic_rods.twisting_energy);
    ImGui::Text("Total Impules = %f", discrete_elastic_rods_driver.discrete_elastic_rods.total_impulse);
    ImGui::SliderFloat("Water Velocity", &discrete_elastic_rods_driver.discrete_elastic_rods.params.water_velocity, -15.0f, 15.0f);
    ImGui::SliderFloat("theta1", &discrete_elastic_rods_driver.discrete_elastic_rods.params.nerve_params.theta1, 1.0f, 4.0f);
    ImGui::SliderFloat("theta2", &discrete_elastic_rods_driver.discrete_elastic_rods.params.nerve_params.theta2, -1.0f, 1.0f);
    ImGui::SliderFloat("theta3", &discrete_elastic_rods_driver.discrete_elastic_rods.params.nerve_params.theta3, 4.0f, 10.0f);
    if (ImGui::Button("Run/Stop Simulation"))
        discrete_elastic_rods_driver.running = !discrete_elastic_rods_driver.running;
    ImGui::SameLine();
    if (ImGui::Button("Run One Time Step Simulation")) discrete_elastic_rods_driver.simulateOneStep();
    if (discrete_elastic_rods_driver.running) discrete_elastic_rods_driver.simulateOneStep();
}

int main(int argc, char* argv[])
{
    omp_set_num_threads(1);
    InputParser input(argc, argv);
    int driver = input.getIntegerCmdOption("-driver");
    int test = input.getIntegerCmdOption("-test");

    // polyscope setup
    polyscope::options::autocenterStructures = true;
    polyscope::options::autoscaleStructures = true;
    switch (driver) {
    case 0: {
        polyscope::options::programName = "Discrete Elastic Rods Simulation";
        polyscope::options::printPrefix = "[DER] ";
    }
        break;
    case 1: {
        polyscope::options::programName = "Discrete Elastic Rods Based Additive Manufactured Parts Simulation";
        polyscope::options::printPrefix = "[DER-AM] ";
    }
        break;
    default: {
        std::cerr << "[ERROR] Please Select A Valid Driver.";
    }
    }

    // driver setup
    polyscope::view::setWindowSize(1920,1080);
    polyscope::init();
    switch (driver) {
    case 0: {
        discrete_elastic_rods_driver.test = test;
        discrete_elastic_rods_driver.initialize();
        polyscope::state::userCallback = discreteElasticRodsDriverCallback;
    }
        break;
    case 1: {
        // TODO: discrete elastic rods based AM simulation
    }
        break;
    default: {
        std::cerr << "[ERROR] Please Select A Valid Driver.";
    }
    }
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;
    polyscope::options::enableVSync = false;
    polyscope::options::maxFPS = 1000;
    polyscope::show();
}
