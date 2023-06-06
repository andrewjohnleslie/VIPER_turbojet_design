//
// Created by Andrew on 05/06/2023.
//

#ifndef VIPER_TURBOJET_DESIGN_DESIGN_H
#define VIPER_TURBOJET_DESIGN_DESIGN_H

#include <cmath>

class Design {

    // Public inheritance: Public member variables in the base class become protected in the derived class; therefore, not
    //                      accessible.
public:

// INLET
// =========================
// =========================

    //Member variables
    double gama = 1.4;
    double gamma_e = 1.3;
    double cp = 1005;
    double cpe = 1244;
    double R = 287;

// Ambient conditions - SEA-LEVEL STATIC TEST BED (0 km, 0 m/s)
    double h = 0;
    double T = 15 - (0.0065 * h);
    double T_a = T + 273.15;
    double T_sls = 288.15;                                               // At h=0 this is equivalent to T_sls_std
    double P_sls = 101325;                                               // At h=0 this is equivalent to P_sls_std
    double P_a = P_sls * (pow( (1 - (0.0065 * (h/T_a))), 5.2561));

    double delta = P_a/P_sls;
    double theta = T_a/T_sls;


    double M_f = 0;
    double a = pow((gama * R * T_sls), 0.5);
    double V_f = a * M_f;

    //Member variables
    double T_02 = T_a * (1 + ((gama - 1) / 2) * pow(M_f, 2));
    double P_02 = P_a * pow((1 + ((gama - 1) / 2) * pow(M_f, 2)),gama/(gama - 1));
    double rho = P_02 / (R * T_02);



// COMPRESSOR
// =======================
// =======================

    double pi_c = 5.5;              // OPR = P_03/P_02
    double e_c = 0.9;

    double T_03 = T_02 * pow(pi_c, (gama - 1)/(e_c * gama));
    double P_03 = pi_c * P_02;
    double tau_c = T_03 / T_02;

    double DelT_c = T_03 - T_02;
    double DelT_t = (cp/cpe) * DelT_c;



// TURBINE
// =======================
// =======================
/*
    void myCustomFunction() {
         if ( y < Pcrit){                   // If P05/Pa < P05/Pcrit (Nozzle PR < Critical PR)
             std::cout << std::endl;
             std::cout << "Nozzle is not choked! \n" << std::endl;
         }else{
             std::cout << std::endl;
             std::cout << "Nozzle is choked! \n " << std::endl;
             // double m_norm_ge = pow(gamma_e, (gamma_e - 1)) * pow(2*(pow(P_a/P_05, 2/gamma_e) - pow(P_a/P_05, (gamma_e + 1)/gamma_e)) , 0.5);     pg 77, 177
         }

    }
*/

    double e_t = 0.85;
    double mdot_a = 23.81;
    double mdot_f = 0.4267;
    double mdot_tot = mdot_a + mdot_f;
    double LCV = 43100000;

    double energy_rel = mdot_f * LCV;               //[Energy_rel = Qr]
    // double sp_energy_rel = energy_rel/mdot_a                             [per kg of air-flow]

    double T_04 = ((energy_rel + (mdot_a * cp * (T_03 - 298))) / ((mdot_a + mdot_f) * cpe)) + 298;
    double P_04 = P_03;                          // Pressure loss in the combustor is neglected: P_04 = (1 - delP_c) * P_03

    double T_05 = T_04 - DelT_t;
    double P_05 = P_04 * pow((T_05/T_04), (gamma_e / (e_t * (gamma_e-1))));

    double x = P_04/P_05;                       // Turbine pressure ratio
    double y = P_05/P_a;                        // Nozzle pressure ratio

    //double Pcrit = 1 / pow((1 - ((gamma_e - 1) / (gamma_e + 1))), (gamma_e / (gamma_e - 1)));
    double Pcrit = pow(((gamma_e + 1)/2), ((gamma_e)/(gamma_e - 1)));                  // Equivalent to P05/Pcrit
    //double Pc = P_05 * Pcrit;

    double m_norm_ge = (gamma_e / pow((gamma_e - 1), 0.5)) * pow(((gamma_e + 1) / 2), ((-1 * (gamma_e + 1)) / (2 * (gamma_e - 1))));      // norm_mf = 1.389 for gamma = 1.3
    // m_norm_4 == m_norm_9

    double A_nozz_t =  ((mdot_tot * pow((cpe * T_04), 0.5)) / (P_04 * m_norm_ge));
    double A_nozz_p =  ((mdot_tot * pow((cpe * T_05), 0.5)) / (P_05 * m_norm_ge));     // THIS COULD BE USED TO WORK OUT K_H

    double k_H = (T_04 - T_05) / T_04;

/*
    void callmyCustomFunction(){
        myCustomFunction();
    }
*/

    //double choke_t = P_05/Pcrit;        // If P_04/P* >= 1.832, the engine is choked for gamma values = 1.3!  pg 177

    double SigT = (2 * e_t * (gamma_e - 1)) / ((2 * gamma_e) - (e_t * (gamma_e - 1)));        // pg 180
    double SigP = (2 * gamma_e) / ((2 * gamma_e) - (e_t * (gamma_e - 1)));                    // pg 180

    double pi_t = pow((A_nozz_t/A_nozz_p), SigP);
    double tau_t = pow((A_nozz_t/A_nozz_p), SigT);



// NOZZLE
// =======================
// =======================

    double T_09 = T_05;

    void myCustomFunction() {
        if ( y < Pcrit){                   // If P05/Pa < P05/Pcrit (Nozzle PR < Critical PR)
            std::cout << std::endl;
            std::cout << "Nozzle is not choked! \n" << std::endl;

            double P_9 = P_a;                    // For NOT CHOCKED: P9 = P09/(P05/Pa) = Pa   [Pa == P9]
            double T_9 = T_05 * pow((P_9/P_05), ((gamma_e - 1)/gamma_e));
            double V_9 = pow((2 * cpe * (T_09 - T_9)), 0.5);

            double rho_9 = P_9 / (R * T_9);
            double A_9 = mdot_tot / (rho_9 * V_9);

            double F_g = V_9 * mdot_tot;

            std::cout << "T_9 = " << T_9 << std::endl;
            std::cout << "P_9 = " << P_9 << std::endl;
            std::cout << "rho_9 = " << rho_9 << std::endl;
            std::cout << "V_9 = " << V_9 << std::endl;          // Calculated jet velocity
            std::cout << "Nozzle exit area  (m^2) = " << A_9 << std::endl;


            if (V_f == 0) {
                double Ft = F_g;
                double sfc = mdot_f / Ft;     // [kg/h/N]
                double Fs = Ft / mdot_a;

                std::cout << "\nNet thrust (F_T) = " << Ft << std::endl;
                std::cout << "Specific thrust (F_S) = " << Fs << std::endl;
                std::cout << "SFC = " << sfc << std::endl;

            } else {
                double Ft = (mdot_tot * V_9) - (mdot_a * V_f) + (A_nozz_p * (P_05 - P_a));
                double sfc = mdot_f / Ft;     // [kg/h/N]
                double Fs = Ft / mdot_a;

                std::cout << "\nNet thrust (F_T) = " << Ft << std::endl;
                std::cout << "Specific thrust (F_S) = " << Fs << std::endl;
                std::cout << "SFC = " << sfc << std::endl;
            }



        }else{
            std::cout << std::endl;
            std::cout << "Nozzle is choked! \n " << std::endl;
            // double m_norm_ge = pow(gamma_e, (gamma_e - 1)) * pow(2*(pow(P_a/P_05, 2/gamma_e) - pow(P_a/P_05, (gamma_e + 1)/gamma_e)) , 0.5);     pg 77, 177

            double P_9 = P_05 / Pcrit;                     // For CHOKED: P9 = P09/Pcrit
            double T_9 = T_05 * pow((P_9/P_05), ((gamma_e - 1)/gamma_e));
            double V_9 = 1 * pow((gamma_e * R * T_9), 0.5);

            double rho_9 = P_9 / (R * T_9);
            double A_9 = mdot_tot / (rho_9 * V_9);

            double F_g = V_9 * mdot_tot;

            std::cout << "T_9 = " << T_9 << std::endl;
            std::cout << "P_9 = " << P_9 << std::endl;
            std::cout << "rho_9 = " << rho_9 << std::endl;
            std::cout << "V_9 = " << V_9 << std::endl;          // Calculated jet velocity
            std::cout << "Nozzle exit area  (m^2) = " << A_9 << std::endl;

             if (V_f == 0) {
                double Ft = F_g;
                double sfc = mdot_f / Ft;     // [kg/h/N]
                double Fs = Ft / mdot_a;

                std::cout << "\nNet thrust (F_T) [GROSS] = " << Ft << std::endl;
                std::cout << "Specific thrust (F_S) = " << Fs << std::endl;
                std::cout << "SFC = " << sfc << std::endl;

            } else {
                double Ft = (mdot_tot * V_9) - (mdot_a * V_f) + (A_nozz_p * (P_05 - P_a));
                double sfc = mdot_f / Ft;     // [kg/h/N]
                double Fs = Ft / mdot_a;

                std::cout << "\nNet thrust (F_T) [NET] = " << Ft << std::endl;
                std::cout << "Specific thrust (F_S) = " << Fs << std::endl;
                std::cout << "SFC = " << sfc << std::endl;
            }

            
        }
    }

    void callmyCustomFunction(){
        myCustomFunction();
    }


/*
    //double P_9 = P_a;     //Assume the nozzle fully expands the flow so that the core jet stream exits at ambient pressure
    double P_9 = P_05 / Pcrit;                     // For CHOKED: P9 = P09/Pcrit
                                                   // For NOT CHOCKED: P9 = P09/(P05/Pa) = Pa   [Pa == P9]
    double T_9 = T_05 * pow((P_9/P_05), ((gamma_e - 1)/gamma_e));
    //double T_09 = T_05;
    double rho_9 = P_9 / (R * T_9);
    // P_09 = P_05;     pg 182


// Calculated thermodynamic cycle
    double V_9_calc = pow((2 * cpe * (T_09 - T_9)), 0.5);        // V_9_calc = V_j
    double V_9 = 1 * pow((gamma_e * R * T_9), 0.5);

    double A_9 = mdot_tot / (rho_9 * V_9);
    //double A_ratio = (m_norm_ge / (y * m_norm_ge)) / ((P_05/P_04) * pow(T_05/T_04, 0.5));
    //double A_4 = A_9 / A_ratio;

    //double F_g = V_9_calc * mdot_tot;       // F_G = 15167 (N)
     double F_g = V_9 * mdot_tot;

// Discrepancy of V_j values is likely due to neglecting the combustor pressure loss, and the simplifying assumptions
    // made for the gas properties (i.e. constant gamma).


    void myCustomFunction1() {
        if (V_f == 0) {
            double Ft = F_g;
            double sfc = mdot_f / Ft;     // [kg/h/N]
            double Fs = Ft / mdot_a;

            std::cout << "Net thrust (F_T) = " << Ft << std::endl;
            std::cout << "Specific thrust (F_S) = " << Fs << std::endl;
            std::cout << "SFC = " << sfc << std::endl;

        } else {
            double Ft = (mdot_tot * V_9_calc) - (mdot_a * V_f) + (A_nozz_p * (P_05 - P_a));
            double sfc = mdot_f / Ft;     // [kg/h/N]
            double Fs = Ft / mdot_a;

            std::cout << "Net thrust (F_T) = " << Ft << std::endl;
            std::cout << "Specific thrust (F_S) = " << Fs << std::endl;
            std::cout << "SFC = " << sfc << std::endl;
        }
    };

    void callmyCustomFunction1(){
        myCustomFunction1();
    }
*/

    Design() = default;          //Constructor
    ~Design() = default;         //Destructor

};
#endif //VIPER_TURBOJET_DESIGN_DESIGN_H
