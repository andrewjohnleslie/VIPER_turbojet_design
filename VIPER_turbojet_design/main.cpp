// This test case is the preliminary test case for a simple turbojet (single-spool) engine.
// The engine being considered and modelled is the Rolls-Royce Viper [Mark 601] turbojet engine.
// This program places the previous code into an inherited class structure.

// The data provided for this engine is at a DESIGN POINT for a sea-level test bed.
// F_G = 15167 (N)
// sfc = 0.993 (kg/h/kg)
// m_a = 23.81 (kg/s)
// m_f = 0.4267 (kg/s)
// OPR [P_03/P_02] = 5.5

#include <iostream>
#include "Design.h"

int main() {
    std::cout << "DESIGN POINT CALCULATION " << std::endl;
    std::cout << "============================== " << std::endl;
    std::cout << "      " << std::endl;

    std::cout << "Ta = " << Design().T_a << std::endl;
    std::cout << "Pa = " << Design().P_a << std::endl;
    std::cout << "Rho_a = " << Design().rho << std::endl;
    std::cout << "T_02 = " << Design().T_02 << std::endl;
    std::cout << "P_02 = " << Design().P_02 << std::endl;
    std::cout << "Theta = " << Design().theta << std::endl;
    std::cout << "Delta = " << Design().delta << std::endl;
    std::cout << "T_03 = " << Design().T_03 << std::endl;
    std::cout << "P_03 = " << Design().P_03 << std::endl;
    std::cout << "T_04 = " << Design().T_04 << std::endl;
    std::cout << "P_04 = " << Design().P_04 << std::endl;
    std::cout << "T_05 = " << Design().T_05 << std::endl;
    std::cout << "P_05 = " << Design().P_05 << std::endl;
    std::cout << "T_09 = " << Design().T_09 << std::endl;
    //std::cout << "T_9 = " << Design().T_9 << std::endl;
    //std::cout << "P_9 = " << Design().P_9 << std::endl;
    //std::cout << "rho_9 = " << Design().rho_9 << std::endl;
    //std::cout << "V_9 (old) = " << Design().V_9_calc << std::endl;          // Calculated jet velocity
    //std::cout << "V_9 (new) = " << Design().V_9 << std::endl;
    std::cout << "mbar_choke_ge = " << Design().m_norm_ge << std::endl;
    std::cout << "K_H = " << Design().k_H << std::endl;
    std::cout << "Turbine nozzle area  (m^2) = " << Design().A_nozz_t << std::endl;
    std::cout << "Propulsive nozzle area  (m^2) = " << Design().A_nozz_p << std::endl;
    //std::cout << "Turbine nozzle area  (m^2) = " << Design().A_4 << std::endl;
    //std::cout << "Nozzle exit area  (m^2) = " << Design().A_9 << std::endl;
    std::cout << "Critical pressure ratio = " << Design().Pcrit << std::endl;
    //std::cout << "Critical pressure = " << Design().Pc << std::endl;
    std::cout << "Turbine pressure ratio = " << Design().x << std::endl;
    std::cout << "Nozzle pressure ratio = " << Design().y << std::endl;
    //std::cout << "Gross thrust (F_T) = " << Design().F_g << std::endl;
    //std::cout << "Net thrust (F_T) = " << Design().Ft << std::endl;
    //std::cout << "Specific thrust (F_S) = " << Design().Fs << std::endl;
    //std::cout << "SFC = " << Design().sfc << std::endl;

    Design nozzle;
    nozzle.callmyCustomFunction();


    //Design preform;
    //preform.callmyCustomFunction1();
}
