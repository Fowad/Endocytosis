#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <ctime>
#include <new>
#include <vector>
#include <random>

using namespace std;

// const int n_filaments = 64, n_central = 16; // n_central = 4, 16, 36

const int n_filaments = 100, n_central = 16; // n_central = 4, 16, 36

const int m = 10, m_p = 4; // m = sqrt(n_filaments); m_p = sqrt(n_central)

double A_push = 20, B_push = 0, A_pull = 131, B_pull = 65.5, sh_pull = 1;
double kapa1_push = 0.9, kapa2_push = 0.3, kapa1_pull = 0.9, kapa2_pull = 0.3, sh_push = 1;

int bend_on = 1; // 1 means tip bending is on and 0 means it's off
int elas_on = 1; // 1 means elasticity is on and 0 means it's off

const unsigned long long int Ntstep = 2*pow(10,9); // Number of time steps.
const unsigned long long int Ntstep_to_store_data = pow(10,6); // Number of time steps that I'm storing data in the text files.
const int limit = Ntstep/(2*pow(10,9));
double dt = pow(10,-8); //[sec]
double t_equilibrium = 0.1; //[sec]

double ss = 10000; // Speed-up Technique

double D_memb = 5*pow(10,3); // [nm^2/sec] obstacle diffusion constant, HINT: this parameter shows the response to the thermal forces.
double mu_memb = 1.22*pow(10,3); // [sec/g] the mobility parameter (1 over drag coefficient) shows the response to the other physical forces. D = mu*kT
double z = 0; //[nm] membrane coordinate.
double opengl_window_scaling_metric = 100; // meaning 1 unit in the opengl window equals that number in nanometer. And the point (0,0) is at the center!
double delta = 2.7; // nm subunit(monomer) size
double beta = 1/4.1;

double D_bend = 5*pow(10,3); // filaments bending diffusion constant
double mu_bend = 1.22*pow(10,3);
double k_bend = 0.1; //(4.1*10000)/(50*50*50) ~ 0.3 [g/sec^2] stiffness of the bending restoring force on the filaments.

double D_elas = 5*pow(10,3);
double mu_elas = 1.22*pow(10,3);
double k_elas = 0.5; // [g/sec^2] spring constant for elastic force.

double dz_bend [m][m];
double dz_bend_avg [m][m], dz_bend_sqr [m][m], dz_bend_sqr_avg [m][m];
double dz_bend_total = 0, dz_bend_sqr_total = 0;
double xb [m][m], xb_prime [m][m];

double dz_elas [m][m];
double dz_elas_avg [m][m], dz_elas_sqr [m][m], dz_elas_sqr_avg [m][m];
double dz_elas_total = 0, dz_elas_sqr_total = 0;
double S [m][m];
double xe [m][m], xe_prime [m][m];

int Num_subunits [m][m];
int Num_subunits_total_pushers = 0, Num_subunits_avg_pushers = 0, Num_subunits_total_pullers = 0, Num_subunits_avg_pullers = 0;
int Num_poly_events [m][m], Num_depoly_events [m][m];
int Num_poly_events_total = 0, Num_depoly_events_total = 0;

double k_on [m][m], k_off [m][m]; // monomer on and off-rate constants [1/microMolar * 1/sec] 11.6 and 1 are the standard values respectively.
double C_0 = 1; // [microMolar] Bulk Monomer Concentration
double Poly_rate [m][m], Depoly_rate [m][m];
double Poly_prob [m][m], Depoly_prob [m][m];

double uni_rand_membrane = 0;
double uni_rand_bend [m][m], uni_rand_elas [m][m];
double uni_rand_poly [m][m], uni_rand_depoly [m][m];

double Filament_Height [m][m];
double R [m][m];
double U [m][m];
double F [m][m];
double F_total_fil_n [m][m];

double F_total_pull = 0, F_total_push = 0, F_total = 0;
double F_total_pull_store = 0, F_total_push_store = 0, F_total_store = 0;
double F_pull_avg = 0, F_push_avg = 0, F_total_avg = 0;

double v_memb = 0; // membrane velocity
double v_memb_total = 0;
double v_memb_avg = 0;

double sc = sqrt (24*dt);
double g_memb = 0; // g has the dimension of [sec]^1/2 derived from <g(t)g(t')> = 2*dt*delta(t,t')
double g_bend [m][m];
double g_elas [m][m];
double F_ext = 0; //pN

double Potential (double R_tip){
    
    for (int r = 0; r < m; r++){for (int c = 0; c < m; c++){

        if ( (r >= (m_p - 1) && r <= (m - m_p)) && (c >= (m_p - 1) && c <= (m - m_p)) )
            
        {
            U[r][c] = (A_pull/kapa1_pull)*exp(-kapa1_pull*(R_tip - sh_pull)) - (B_pull/kapa2_pull)*exp(-kapa2_pull*(R_tip - sh_pull));
        }
        
        else
            
        {
            U[r][c] = (A_push/kapa1_push)*exp(-kapa1_push*(R_tip - sh_push)) - (B_push/kapa2_push)*exp(-kapa2_push*(R_tip - sh_push));
        }
        
        return U[r][c];
    }}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//========================================================================================================================================//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main () {
    
    //    ifstream F_ext_data ("F_ext_data.txt");
    //    F_ext_data.is_open();
    //    F_ext_data >> F_ext;
    
    ofstream Sc_zdata ("Scaled_zdata.txt", ofstream::app); // Scaled zdata is for -1 to 1 normalization of GL window
    ofstream Time_Course ("Time_Course.txt", ofstream::app);

    ofstream Force_Dist ("Force_Dist.txt");
//    ofstream F_V_Plot ("F_V_Plot.txt", ofstream::app);
//    ofstream F_V_Plot_per_fil ("F_V_Plot_per_fil.txt", ofstream::app);
//    ofstream Power ("Power_Load.txt", ofstream::app);
    ofstream Data ("Data.txt", ofstream::app);
    
    random_device seed;
    mt19937 rgen(seed());
    
    /////////////////////////////////////////////
    //=========================================//
    /////////////////////////////////////////////
    

    for (int r = 0; r < m; r++){for (int c = 0; c < m; c++){

        if ( (r >= (m_p - 1) && r <= (m - m_p)) && (c >= (m_p - 1) && c <= (m - m_p)) )

        {k_on [r][c] = 11.6; k_off [r][c] = 1; Poly_rate [r][c] = k_on[r][c]*C_0; Depoly_rate [r][c] = k_off[r][c]*C_0; Poly_prob [r][c] = Poly_rate[r][c]*dt; Depoly_prob [r][c] = Depoly_rate[r][c]*dt; Num_subunits [r][c] = 1; Filament_Height [r][c] = Num_subunits[r][c]*delta + (double(r)/double(m))*delta; dz_bend[r][c] = 0; dz_elas[r][c] = 0;}
        
        else {k_on [r][c] = 11.6; k_off [r][c] = 1; Poly_rate [r][c] = k_on[r][c]*C_0; Depoly_rate [r][c] = k_off[r][c]*C_0; Poly_prob [r][c] = Poly_rate[r][c]*dt; Depoly_prob [r][c] = Depoly_rate[r][c]*dt; Num_subunits [r][c] = 1; Filament_Height [r][c] = Num_subunits[r][c]*delta + (double(r)/double(m))*delta; dz_bend[r][c] = 0; dz_elas[r][c] = 0;}
    
    }}
    
    
    for (int r = 0; r < m; r++){for (int c = 0; c < m; c++){ofstream Num_subunits_data (("Num_subunits_" + to_string(r) + "_" + to_string(c) + ".txt"), ofstream::trunc); ofstream dz_bend_data (("dz_bend_" + to_string(r) + "_" + to_string(c) + ".txt"), ofstream::trunc); ofstream dz_elas_data (("dz_elas_" + to_string(r) + "_" + to_string(c) + ".txt"), ofstream::trunc);}}
    
    
    uniform_int_distribution<unsigned long long int> dist_membrane(0,pow(10,9));
    uniform_int_distribution<unsigned long long int> dist_poly(0,pow(10,9));
    uniform_int_distribution<unsigned long long int> dist_depoly(0,pow(10,9));
    uniform_int_distribution<unsigned long long int> dist_bend(0,pow(10,9));
    uniform_int_distribution<unsigned long long int> dist_elas(0,pow(10,9));

    
    /////////////////////////////////////////////
    //=========================================//
    /////////////////////////////////////////////
    
    for (int j=0; j<limit; j++){
        
        for (int i=0; i<(Ntstep/limit); i++){
            
            double t = i*dt + j*(Ntstep/limit)*dt;
            
            uni_rand_membrane = (double(dist_membrane(rgen))/double(pow(10,9))) - 0.5;
            g_memb = sc*uni_rand_membrane;
            
            /////////////////////////////////////////////
            //=========================================//
            /////////////////////////////////////////////
            
            
            if (i % 10000 == 0){
                
                for (int r = 0; r < m; r++){for (int c = 0; c < m; c++){
    
                    uni_rand_poly [r][c] = double(dist_poly(rgen))/(double(pow(10,9)));
                    
                    uni_rand_depoly [r][c] = double(dist_depoly(rgen))/(double(pow(10,9)));
                    
                    double R_tip = R[r][c];
                    
                    if ( Potential ((R_tip - delta)) > Potential (R_tip) )
                    {
                        k_on[r][c] = 11.6*exp(-beta*( Potential((R_tip - delta)) - Potential(R_tip) ));
                    }
                        
                    else
                    {
                        k_on[r][c] = 11.6;
                    }
                       
                        
                    if ( Potential ((R_tip + delta)) > Potential (R_tip) )
                    {
                        k_off[r][c] = 1*exp(-beta*( Potential((R_tip + delta)) - Potential(R_tip) ));
                    }
                    
                    else
                    {
                        k_off[r][c] = 1;
                    }
                    
                    Poly_prob[r][c] = k_on[r][c]*C_0*dt;
                    Depoly_prob[r][c] = k_off[r][c]*C_0*dt;
                    
                    if (ss*Poly_prob[r][c] > uni_rand_poly[r][c]){Num_subunits[r][c]++; Num_poly_events[r][c]++;}
                    if (ss*Depoly_prob[r][c] > uni_rand_depoly[r][c]){Num_subunits[r][c]--; Num_depoly_events[r][c]++;}
                }}
                
            }
            
            
            
//            if (i % Ntstep_to_store_data*10 == 0){
//
//            for (int r = 0; r < m; r++){for (int c = 0; c < m; c++){
//
//                double R_tip = R[r][c];
//
//                cout << Potential (R_tip) << ' ' << Potential (R_tip-delta) << endl;
//
//            }}}
            
            
            
            for (int r = 0; r < m; r++){for (int c = 0; c < m; c++){
                
                if (bend_on == 1){uni_rand_bend [r][c] = (double(dist_bend(rgen))/double(pow(10,9))) - 0.5; g_bend [r][c] = sc*uni_rand_bend [r][c];}
                if (elas_on == 1){uni_rand_elas [r][c] = (double(dist_elas(rgen))/double(pow(10,9))) - 0.5; g_elas [r][c] = sc*uni_rand_elas [r][c];}
                
                Filament_Height [r][c] = Num_subunits[r][c]*delta + (double(r)/double(m))*delta + dz_bend[r][c] + dz_elas[r][c];
                
                R [r][c] = z - Filament_Height [r][c];
                
                if ( (r >= (m_p - 1) && r <= (m - m_p)) && (c >= (m_p - 1) && c <= (m - m_p)) )
                {F[r][c] = A_pull*exp(-kapa1_pull*(R[r][c] - sh_pull)) - B_pull*exp(-kapa2_pull*(R[r][c] - sh_pull)); F_total_pull += F[r][c];}
                else {F[r][c] = A_push*exp(-kapa1_push*(R[r][c] - sh_push)) - B_push*exp(-kapa2_push*(R[r][c] - sh_push)); F_total_push += F[r][c];}
                
                F_total_fil_n [r][c] += F [r][c];
                
                if (i % Ntstep_to_store_data == 0){ofstream Num_subunits_data (("Num_subunits_" + to_string(r) + "_" + to_string(c) + ".txt"), ofstream::app); Num_subunits_data << Num_subunits[r][c] << endl;}
                
            }}
            
            
            
            /////////////////////////////////////////////
            //=========================================//
            /////////////////////////////////////////////
            
            if (i % Ntstep_to_store_data == 0){Sc_zdata << z/opengl_window_scaling_metric << endl; Time_Course << t << ' ' << z << endl;}
            
            /////////////////////////////////////////////
            //=========================================//
            /////////////////////////////////////////////
            
            if (elas_on==1){
                
                for (int r = 0; r < m; r++){for (int c = 0; c < m; c++){
                    
                    if ( (r = 0) && (c = 0) ) { S[r][c] = (dz_elas[r][c+1] - dz_elas[r][c]) + (dz_elas[r+1][c] - dz_elas[r][c]); }
                    
                    if ( (r = 0) && (c = m-1) ) { S[r][c] = ( dz_elas[r][c-1] - dz_elas[r][c]) + (dz_elas[r+1][c] - dz_elas[r][c]); }
                    
                    if ( (r = m-1) && (c = 0) ) { S[r][c] = ( dz_elas[r-1][c] - dz_elas[r][c]) + (dz_elas[r][c+1] - dz_elas[r][c]); }
                    
                    if ( (r = m-1) && (c = m-1) ) { S[r][c] = (dz_elas[r][c-1] - dz_elas[r][c]) + (dz_elas[r][c] - dz_elas[r-1][c]); }

                    if ( (r = 0) && (c > 0 && c < m-1) ) { S[r][c] = (dz_elas[r][c+1] - dz_elas[r][c]) + (dz_elas[r][c-1] - dz_elas[r][c]) + (dz_elas[r+1][c] - dz_elas[r][c]); }
                    
                    if ( (c = 0) && (r > 0 && r < m-1) ) { S[r][c] = (dz_elas[r-1][c] - dz_elas[r][c]) + (dz_elas[r+1][c] - dz_elas[r][c]) + (dz_elas[r][c+1] - dz_elas[r][c]); }
                    
                    if ( (r = m-1) && (c > 1 && c < m-1) ) { S[r][c] = (dz_elas[r][c+1] - dz_elas[r][c]) + (dz_elas[r][c-1] - dz_elas[r][c]) + (dz_elas[r-1][c] - dz_elas[r][c]); }
                
                    if ( (c = m-1) && (r > 1 && r < m-1) ) { S[r][c] = (dz_elas[r-1][c] - dz_elas[r][c]) + (dz_elas[r+1][c] - dz_elas[r][c]) + (dz_elas[r][c-1] - dz_elas[r][c]); }
                    
                    else { S[r][c] = (dz_elas[r-1][c] - dz_elas[r][c]) + (dz_elas[r][c+1] - dz_elas[r][c]) + (dz_elas[r+1][c] - dz_elas[r][c]) + (dz_elas[r][c-1] - dz_elas[r][c]); }
                
                }}
                
            }
            
            /////////////////////////////////////////////
            //=========================================//
            /////////////////////////////////////////////
            
            for (int r = 0; r < m; r++){for (int c = 0; c < m; c++){
                
                if (bend_on == 1){
                    if ( (r >= (m_p - 1) && r <= (m - m_p)) && (c >= (m_p - 1) && c <= (m - m_p)) ){dz_bend [r][c] = 0;}else{ dz_bend [r][c] = dz_bend [r][c] + g_bend[r][c]*sqrt(D_bend) + mu_bend*dt*(-F[r][c] - k_bend*(dz_bend[r][c]));}}
                
//                if (bend_on == 1 && (r == 0 && c == 0)){
//                    dz_bend [r][c] = dz_bend [r][c] + g_bend[r][c]*sqrt(D_bend) + mu_bend*dt*(-F[r][c] - k_bend*(dz_bend[r][c]));
//                }
                
                if (elas_on == 1){dz_elas [r][c] = dz_elas [r][c] + g_elas[r][c]*sqrt(D_elas) + mu_elas*dt*(-F[r][c] + k_elas*S[r][c] - k_elas*dz_elas[r][c]);}
                
                xb[r][c] += dz_bend[r][c]; dz_bend_sqr[r][c] = dz_bend[r][c]*dz_bend[r][c]; xb_prime[r][c] += dz_bend_sqr[r][c];
                xe[r][c] += dz_elas[r][c]; dz_elas_sqr[r][c] = dz_elas[r][c]*dz_elas[r][c]; xe_prime[r][c] += dz_elas_sqr[r][c];
                
            }}
            
            
            if(i % Ntstep_to_store_data == 0){for (int r = 0; r < m; r++){for (int c = 0; c < m; c++){ofstream dz_bend_data (("dz_bend_" + to_string(r) + "_" + to_string(c) + ".txt"), ofstream::app); dz_bend_data << dz_bend[r][c] << endl; ofstream dz_elas_data (("dz_elas_" + to_string(r) + "_" + to_string(c) + ".txt"), ofstream::app); dz_elas_data << dz_elas[r][c] << endl;}}}
            

            /////////////////////////////////////////////
            //=========================================//
            /////////////////////////////////////////////
            
            F_total = F_total_pull + F_total_push;
            
            z = z + g_memb*sqrt(D_memb) + mu_memb*dt*(-F_ext + F_total);
            
            v_memb = (g_memb*sqrt(D_memb) + mu_memb*dt*(-F_ext + F_total))/dt;
            
            
            if (t>t_equilibrium){
                v_memb_total += v_memb;
                F_total_pull_store += F_total_pull;
                F_total_push_store += F_total_push;
                F_total_store += F_total;}
            
            
            F_total_pull = 0;
            F_total_push = 0;
            F_total = 0;
            
        } // closing the time iteration loop
        
        cout << "j = " << j << ' ' << " limit = " << limit << ' ' << " Num_subunits[0][0] = " << Num_subunits[0][0] << endl;
    }
    
    /////////////////////////////////////////////
    //=========================================//
    /////////////////////////////////////////////
    
    const unsigned long long int Ntstep_after_equ = Ntstep*(1-(t_equilibrium/(Ntstep*dt)));
    
    F_pull_avg = F_total_pull_store/Ntstep_after_equ;
    F_push_avg = F_total_push_store/Ntstep_after_equ;
    F_total_avg = F_total_store/Ntstep_after_equ;
    v_memb_avg = v_memb_total/Ntstep_after_equ;
    
    for (int r = 0; r < m; r++){for (int c = 0; c < m; c++){
        
        Force_Dist << r+1 << ' ' << c+1 << ' ' << F_total_fil_n[r][c]/Ntstep_after_equ << endl;

        dz_bend_avg[r][c] = (xb[r][c]/Ntstep_after_equ);
        dz_bend_sqr_avg[r][c] = (xb_prime[r][c]/Ntstep_after_equ);
        dz_elas_avg[r][c] = (xe[r][c]/Ntstep_after_equ);
        dz_elas_sqr_avg[r][c] = (xe_prime[r][c]/Ntstep_after_equ);
        
        if ( (r >= (m_p - 1) && r <= (m - m_p)) && (c >= (m_p - 1) && c <= (m - m_p)) ){
            
            Num_subunits_total_pullers += Num_subunits [r][c]; Num_subunits_avg_pullers = Num_subunits_total_pullers/(n_central);}
        
        else {Num_subunits_total_pushers += Num_subunits [r][c]; Num_subunits_avg_pushers = Num_subunits_total_pushers/(n_filaments - n_central);}
        
        
    }}
    
    /////////////////////////////////////////////
    //=========================================//
    /////////////////////////////////////////////
    
    Data << endl << "/////////////////////////////////////////////////////////////////" << endl << "n_filaments is:" << n_filaments << " and n_central is:" << n_central << " and the run time is:" << Ntstep*dt << " sec" << endl << endl <<
    "A_push:" << A_push << " B_push:" << B_push << " A_pull:" << A_pull << " B_pull:" << B_pull << " sh_pull:" << sh_pull << endl << endl <<
    "kapa1_push:" << kapa1_push << " kapa2_push:" << kapa2_push << " kapa1_pull:" << kapa1_pull << " kapa2_pull:" << kapa2_pull << " sh_push:" <<
    sh_push << endl << endl << "sqrt(2DT)/(mu*<F_total>*T) = " << (sqrt(2*D_memb*dt*Ntstep))/(mu_memb*F_total_avg*dt*Ntstep) << endl <<
    "<v_membrane>/(mu*<F_total>) = " << v_memb_avg/(mu_memb*F_total_avg) << endl << "F_ext is: " << F_ext << endl << endl << "<F_pull> is: " <<
    F_pull_avg << "pN and <F_push> is: " << F_push_avg << "pN and <F_total> is: " << F_total_avg << "pN" << endl << "<v_membrane> is: " << v_memb_avg
    << " nm/sec and <v_puller> is: " << (Num_subunits_avg_pullers*delta)/(Ntstep*dt) << " nm/sec and <v_pusher> is: " <<
    (Num_subunits_avg_pushers*delta)/(Ntstep*dt) << " nm/sec" << endl << endl << endl;
    
    for (int r = 0; r < m; r++){for (int c = 0; c < m; c++){Data << "Filament[" << r << "][" << c << "] has " << Num_subunits[r][c] << " subunits, therefore " << Num_poly_events[r][c] << " polymerization events and " << Num_depoly_events[r][c] << " depolymerization events." << endl;}}
    
    Data << endl << endl;
    
    for (int r = 0; r < m; r++){for (int c = 0; c < m; c++){Data << "Filament[" << r << "][" << c << "] has <dz_bend>: " << dz_bend_avg[r][c] << " nm and <dz_bend^2>: " << dz_bend_sqr_avg[r][c] << "nm^2 and <dz_elas>: " << dz_elas_avg[r][c] << " nm and <dz_elas^2>: " << dz_elas_sqr_avg[r][c] << endl;}}
    
    
    /////////////////////////////////////////////
    //=========================================//
    /////////////////////////////////////////////
    
    
    cout << endl;
    
    cout << "n_filaments is:" << n_filaments << " and n_central is:" << n_central << " and the run time is: " << Ntstep*dt << " sec" << endl;
    
    cout << endl;
    
    if (bend_on == 1){cout << "bending is on!" << endl;}else{cout << "bending is off!" << endl;}
    if (elas_on == 1){cout << "elasticity is on!" << endl;}else{cout << "elasticity is off!" << endl;}
    
    cout << endl;
    
//    for (int r = 0; r < m; r++){for (int c = 0; c < m; c++){cout << "Filament[" << r << "][" << c << "] has " << Num_subunits[r][c] << " subunits, therefore " << Num_poly_events[r][c] << " polymerization events and " << Num_depoly_events[r][c] << " depolymerization events." << endl;}}
//
//    cout << endl;
//
//    for (int r = 0; r < m; r++){for (int c = 0; c < m; c++){cout << "Filament[" << r << "][" << c << "] has <dz_bend>: " << dz_bend_avg[r][c] << " nm and <dz_bend^2>: " << dz_bend_sqr_avg[r][c] << "nm^2 and <dz_elas>: " << dz_elas_avg[r][c] << " nm and <dz_elas^2>: " << dz_elas_sqr_avg[r][c] << endl;}}
    
    cout << endl;
    
    cout << "sqrt(2DT)/(mu*<F_total>*T) = " << (sqrt(2*D_memb*dt*Ntstep))/(mu_memb*F_total_avg*dt*Ntstep) << endl;
    cout << "<v_membrane>/(mu*<F_total>) = " << v_memb_avg/(mu_memb*F_total_avg) << endl;
    cout << "F_ext is: " << F_ext << endl;
    
    cout << "<F_pull> is: " << F_pull_avg << "pN and <F_push> is: " << F_push_avg << "pN and <F_total> is: " << F_total_avg << "pN" << endl;
    cout << "<v_membrane> is: " << v_memb_avg << " nm/sec and <v_puller> is: " << (Num_subunits_avg_pullers*delta)/(Ntstep*dt) << " nm/sec and <v_pusher> is: " << (Num_subunits_avg_pushers*delta)/(Ntstep*dt) << " nm/sec" << endl;
    
    cout << endl;
    
    //      system ("xmgrace Scaled_zdata.txt");
    
    return 0;
}
