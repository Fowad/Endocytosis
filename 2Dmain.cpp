#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <ctime>
#include <new>
#include <random>
#include <vector>
#include <mpi.h>

using namespace std;

const int n_filaments = 144, n_central = 36; // n_central = 4, 16, 36
const int n_inner_central = 16, n_outer_central = 20;

const int m = 12, m_p = 6; // m = sqrt(n_filaments); m_p = sqrt(n_central)

int seed_1 = 0;

double A_push = 20, B_push = 0, kappa1_push = 0.9, kappa2_push = 0.3, sh_push = 1;
double A_pull = 26, B_pull = 13, kappa1_pull = 0.9, kappa2_pull = 0.3, sh_pull = 1;

const unsigned long long int Ntstep = 20*pow(10,9); // Number of time steps.
const unsigned long long int Ntstep_to_store_data = pow(10,6); // Number of time steps that I'm storing data in the text files.
const int limit = 10; // limit = Ntstep/(1*10^9) in denominator is the size of nested "time" for loop.

double dt = 5*pow(10,-10); //[sec]
double t_equilibrium = 3; //[sec] j = 0, 1, and 2 won't be included in error bar calculations then.
double t_break = 0;
int Break = 0;

const int ss = 100; // Speed-up Technique

double D_memb = pow(10,4); // [nm^2/sec] obstacle diffusion constant, HINT: this parameter shows the response to the thermal forces.
double mu_memb = 0.2433*pow(10,4); // [sec/g] the mobility parameter (1 over drag coefficient) shows the response to the other physical forces. D = mu*kT
double z = 4.5; //[nm] membrane coordinate.
double opengl_window_scaling_metric = 100; // meaning 1 unit in the opengl window equals that number in nanometer. And the point (0,0) is at the center!
double delta = 2.21; // nm subunit(monomer) size. Theta is 35 degree.
double beta = 0.2433;

double D_bend = pow(10,5); // filaments bending diffusion constant
double mu_bend = 0.2433*pow(10,5);
double k_bend = 4.17; //// Lp = 17.5 micron; L = 54 nm; Theta = 35 degree.
double dz_bend [m][m];
double Error_bar_dz_bend [m][m];
double dz_bend_total_each_fil [m][m];
double dz_bend_total_each_fil_each_j [m][m];
double dz_bend_avg_each_fil [m][m];

double D_elas = pow(10,5);
double mu_elas = 0.2433*pow(10,5);
double k_elas = 0.53; // [g/sec^2] spring constant for elastic force.
double dz_elas [m][m];
double Error_bar_dz_elas [m][m];
double dz_elas_total_each_fil [m][m];
double dz_elas_total_each_fil_each_j [m][m];
double dz_elas_avg_each_fil [m][m];
double S [m][m];

double dz_bend_avg_total_pullers = 0, dz_elas_avg_total_pullers = 0;
double dz_bend_avg_total_pushers = 0, dz_elas_avg_total_pushers = 0;
double dz_bend_avg_array_pullers = 0, dz_elas_avg_array_pullers = 0;
double dz_bend_avg_array_pushers = 0, dz_elas_avg_array_pushers = 0;

double dz_elas_avg_pullers_at_each_timestep = 0, dummy1 = 0;
double dz_elas_avg_inner_pullers_at_each_timestep = 0, dummy2 = 0;
double dz_elas_avg_outer_pullers_at_each_timestep = 0, dummy3 = 0;
double dz_elas_avg_pushers_at_each_timestep = 0, dummy4 = 0;

int Num_subunits [m][m];
int Num_subunits_total_pullers = 0, Num_subunits_avg_pullers = 0;
int Num_subunits_total_pushers = 0, Num_subunits_avg_pushers = 0;
int Num_poly_events [m][m], Num_depoly_events [m][m];

double k_on [m][m], k_off [m][m]; // monomer on and off-rate constants [1/microMolar * 1/sec] 11.6 and 1 are the standard free values respectively.
double C_0 = 5.3; // [microMolar] Bulk Actin Monomer Concentration
double Poly_rate [m][m], Depoly_rate [m][m];
double Poly_prob [m][m], Depoly_prob [m][m];

double uni_rand_membrane = 0;
double uni_rand_elas [m][m], uni_rand_bend [m][m];
double uni_rand_poly [m][m], uni_rand_depoly [m][m];

double filament_height [m][m];
double initial_staggered_length [m][m];
double r [m][m];
double U [m][m];
double F [m][m];

double Error_bar_Force [m][m];
double F_total_each_fil [m][m];
double F_total_each_fil_each_j [m][m];
double F_avg_each_fil [m][m];

double Error_bar_symmetrized_Force [m][m];
double F_symmetrized_each_fil [m][m];
double F_symmetrized_each_fil_each_j [m][m];
double F_symmetrized_avg_each_fil [m][m];

double Error_bar_Pulling_Force = 0, Error_bar_Pushing_Force = 0;
double Error_bar_Inner_Pulling_Force = 0, Error_bar_Outer_Pulling_Force = 0;
double Error_bar_Total_Force = 0, Error_bar_v_memb = 0;

double F_total_pull_to_draw = 0, F_total_push_to_draw = 0, F_total_to_draw = 0;
double F_total_pull = 0, F_total_push = 0, F_total_inner_pull = 0, F_total_outer_pull = 0, F_total = 0;
double F_total_pull_each_j = 0, F_total_push_each_j = 0, F_total_inner_pull_each_j = 0, F_total_outer_pull_each_j = 0;

double F_pull_avg = 0, F_push_avg = 0, F_total_avg = 0;
double F_inner_pull_avg = 0, F_outer_pull_avg = 0;

double v_memb = 0; // membrane velocity
double v_memb_total = 0;
double v_memb_total_each_j = 0;
double v_memb_avg = 0;
//double v_fil_avg [m][m];
//double Error_bar_v_fil [m][m];

double sc = sqrt (24*dt);
double g_memb = 0; // g has the dimension of [sec]^1/2 derived from <g(t)g(t')> = 2*dt*delta(t,t')
double g_bend [m][m];
double g_elas [m][m];
double F_ext = 0; //pN

double sum_Pr = 0;

double Potential (double r, int row, int col){
    
        if ( (row >= ((m - m_p)/2) && row <= (m + m_p)/2 - 1) && (col >= ((m - m_p)/2) && col <= (m + m_p)/2 - 1) )
            
        {
            U[row][col] = (A_pull/kappa1_pull)*exp(-kappa1_pull*(r - sh_pull)) - (B_pull/kappa2_pull)*exp(-kappa2_pull*(r - sh_pull));
        }
        
        else
            
        {
            U[row][col] = (A_push/kappa1_push)*exp(-kappa1_push*(r - sh_push)) - (B_push/kappa2_push)*exp(-kappa2_push*(r - sh_push));
        }
        
        return U[row][col];
}

double Force (double r, int row, int col){
    
    if ( (row >= ((m - m_p)/2) && row <= (m + m_p)/2 - 1) && (col >= ((m - m_p)/2) && col <= (m + m_p)/2 - 1) )
        
    {
        F[row][col] = A_pull*exp(-kappa1_pull*(r - sh_pull)) - B_pull*exp(-kappa2_pull*(r - sh_pull));
    }
    
    else
        
    {
        F[row][col] = A_push*exp(-kappa1_push*(r - sh_push)) - B_push*exp(-kappa2_push*(r - sh_push));
    }
    
    return F[row][col];
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//========================================================================================================================================//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main (int argc, char** argv ) {

   MPI_Init(NULL,NULL); //initialize MPI
   
   //Get number of processes
   int world_size;
   MPI_Comm_size(MPI_COMM_WORLD, &world_size);
   
   //Get number of processes
   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   // Get the name of the processor
   char processor_name[MPI_MAX_PROCESSOR_NAME];
   int name_len;
   MPI_Get_processor_name(processor_name, &name_len);
   printf("Processor %s, rank %d out of %d processors\n",processor_name,rank, world_size);

   //set parameters for run
   const int MAXPROC = 40;
   vector<vector<double>> params; //parameter table for MPI parallel run
   vector<vector<double>>::size_type ppos = 0; //iterator for vector of vectors
   vector<double> empty_dbls; //empty dbl vector for initialization

   //create vector of empty vectors
   for( ppos=0; ppos!=MAXPROC; ppos++)
   {
     params.push_back(empty_dbls);
   }

   //fill in parameter table

   //construct set 0
   ppos = 0;
   params[ppos].push_back(0.53);          //add k_elas
   params[ppos].push_back(0.1);         //add t_equilibrium
   params[ppos].push_back(26.0);        //add A_pull
   params[ppos].push_back(13.0);         //add B_pull
   params[ppos].push_back(20.0);        //add A_push
   params[ppos].push_back(0.0);         //add B_push
   params[ppos].push_back(0.0);         //add F_ext
   params[ppos].push_back(0);         //add seed_1

    
   //construct set 1
   ppos = 1;
   params[ppos].push_back(0.53);          //add k_elas
   params[ppos].push_back(3);         //add t_equilibrium
   params[ppos].push_back(131.0);        //add A_pull
   params[ppos].push_back(65.5);         //add B_pull
   params[ppos].push_back(20.0);        //add A_push
   params[ppos].push_back(0.0);         //add B_push
   params[ppos].push_back(0.0);         //add F_ext
   params[ppos].push_back(1);         //add seed_1


   //construct set 2
   ppos = 2;
   params[ppos].push_back(0.53);        //add k_elas
   params[ppos].push_back(3);          //add t_equilibrium
   params[ppos].push_back(262.0);        //add A_pull
   params[ppos].push_back(131.0);         //add B_pull
   params[ppos].push_back(20.0);        //add A_push
   params[ppos].push_back(0.0);         //add B_push
   params[ppos].push_back(0.0);         //add F_ext
   params[ppos].push_back(2);         //add seed_1


   //construct set 3
   ppos = 3;
   params[ppos].push_back(0.53);          //add k_elas
   params[ppos].push_back(3);          //add t_equilibrium
   params[ppos].push_back(520.0);        //add A_pull
   params[ppos].push_back(260.0);         //add B_pull
   params[ppos].push_back(20.0);        //add A_push
   params[ppos].push_back(0.0);         //add B_push
   params[ppos].push_back(0.0);         //add F_ext
   params[ppos].push_back(3);         //add seed_1

    
   //construct set 4
   ppos = 4;
   params[ppos].push_back(0.265);          //add k_elas
   params[ppos].push_back(3);          //add t_equilibrium
   params[ppos].push_back(131.0);        //add A_pull
   params[ppos].push_back(65.5);         //add B_pull
   params[ppos].push_back(20.0);        //add A_push
   params[ppos].push_back(0.0);         //add B_push
   params[ppos].push_back(0.0);         //add F_ext
   params[ppos].push_back(4);         //add seed_1


   //construct set 5
   ppos = 5;
   params[ppos].push_back(1.06);          //add k_elas
   params[ppos].push_back(3);          //add t_equilibrium
   params[ppos].push_back(131.0);        //add A_pull
   params[ppos].push_back(65.5);         //add B_pull
   params[ppos].push_back(20.0);        //add A_push
   params[ppos].push_back(0.0);         //add B_push
   params[ppos].push_back(0.0);         //add F_ext
   params[ppos].push_back(5);         //add seed_1


   //construct set 6
   ppos = 6;
   params[ppos].push_back(0.265);          //add k_elas
   params[ppos].push_back(3);          //add t_equilibrium
   params[ppos].push_back(262.0);        //add A_pull
   params[ppos].push_back(131.0);         //add B_pull
   params[ppos].push_back(20.0);        //add A_push
   params[ppos].push_back(0.0);         //add B_push
   params[ppos].push_back(0.0);         //add F_ext
   params[ppos].push_back(6);         //add seed_1


   //construct set 7
   ppos = 7;
   params[ppos].push_back(1.06);          //add k_elas
   params[ppos].push_back(3);          //add t_equilibrium
   params[ppos].push_back(262.0);        //add A_pull
   params[ppos].push_back(131.0);         //add B_pull
   params[ppos].push_back(20.0);        //add A_push
   params[ppos].push_back(0.0);         //add B_push
   params[ppos].push_back(0.0);         //add F_ext
   params[ppos].push_back(7);         //add seed_1

    
   //construct set 8
   ppos = 8;
   params[ppos].push_back(0.53);          //add k_elas
   params[ppos].push_back(3);          //add t_equilibrium
   params[ppos].push_back(131.0);        //add A_pull
   params[ppos].push_back(65.5);         //add B_pull
   params[ppos].push_back(13.0);        //add A_push
   params[ppos].push_back(6.5);         //add B_push
   params[ppos].push_back(0.0);         //add F_ext
   params[ppos].push_back(8);         //add seed_1

    
    //construct set 9
    ppos = 9;
    params[ppos].push_back(0.53);          //add k_elas
    params[ppos].push_back(3);          //add t_equilibrium
    params[ppos].push_back(131.0);        //add A_pull
    params[ppos].push_back(65.5);         //add B_pull
    params[ppos].push_back(26.0);        //add A_push
    params[ppos].push_back(13.0);         //add B_push
    params[ppos].push_back(0.0);         //add F_ext
    params[ppos].push_back(9);         //add seed_1

    
    //construct set 10
    ppos = 10;
    params[ppos].push_back(0.53);          //add k_elas
    params[ppos].push_back(3);          //add t_equilibrium
    params[ppos].push_back(131.0);        //add A_pull
    params[ppos].push_back(65.5);         //add B_pull
    params[ppos].push_back(52.0);        //add A_push
    params[ppos].push_back(26.0);         //add B_push
    params[ppos].push_back(0.0);         //add F_ext
    params[ppos].push_back(10);         //add seed_1

    
    //construct set 11
    ppos = 11;
    params[ppos].push_back(0.53);          //add k_elas
    params[ppos].push_back(3);          //add t_equilibrium
    params[ppos].push_back(262.0);        //add A_pull
    params[ppos].push_back(131.0);         //add B_pull
    params[ppos].push_back(13.0);        //add A_push
    params[ppos].push_back(6.5);         //add B_push
    params[ppos].push_back(0.0);         //add F_ext
    params[ppos].push_back(11);         //add seed_1

    
    //construct set 12
    ppos = 12;
    params[ppos].push_back(0.53);          //add k_elas
    params[ppos].push_back(3);          //add t_equilibrium
    params[ppos].push_back(262.0);        //add A_pull
    params[ppos].push_back(131.0);         //add B_pull
    params[ppos].push_back(26.0);        //add A_push
    params[ppos].push_back(13.0);         //add B_push
    params[ppos].push_back(0.0);         //add F_ext
    params[ppos].push_back(12);         //add seed_1

    
    //construct set 13
    ppos = 13;
    params[ppos].push_back(0.53);          //add k_elas
    params[ppos].push_back(3);          //add t_equilibrium
    params[ppos].push_back(262.0);        //add A_pull
    params[ppos].push_back(131.0);         //add B_pull
    params[ppos].push_back(52.0);        //add A_push
    params[ppos].push_back(26.0);         //add B_push
    params[ppos].push_back(0.0);         //add F_ext
    params[ppos].push_back(13);         //add seed_1

    
    if( rank == 0 ){
       if(world_size > MAXPROC){
         cout << "Max number of processes is %i !" << MAXPROC <<  endl;
         return 0;
       } 
    }

    //Now load the parameters for each process
    vector<vector<double>>::size_type my_params = static_cast<vector<vector<double>>::size_type>(rank);
    k_elas = params[my_params][0]; //add k_elas
    t_equilibrium = params[my_params][1]; //add t_equilibrium
    A_pull = params[my_params][2];        //add A_pull
    B_pull = params[my_params][3];         //add B_pull
    A_push = params[my_params][4];        //add A_push
    B_push = params[my_params][5];         //add B_push
    F_ext = params[my_params][6];         //add F_ext
    seed_1 = params[my_params][7];         //add seed_1

    
    //next, let each process create its own directory
    //they will be labeld run00, run01, and so on
    string command="mkdir -p ";//string for command
    string dirname="Run"; //part of directory name
    string dirnamefull ; // directory name
    if( rank < 10) dirname+="0"; //0 padded directory names
    dirname += to_string(rank); //construct directory name
    dirnamefull = dirname + "/";
    string path = dirnamefull;
    //const char* path = dirnamefull.c_str(); //store directory name
    command += dirname; //add directory name to mkdir command
    const char* create_dir = command.c_str(); //convert string to const char* as required by system()
    //cout << dirname << endl;
    const int dir_err = system(create_dir); //create directory
    if (dir_err == -1 )
    {
       printf("Can't create directory. Exiting...\n");
       exit(1);
    }
    
    /////////////////////////////////////////////
    //=========================================//
    /////////////////////////////////////////////

    ofstream Paramfile (path + "Parameters.txt", ofstream::app); // parameters for this run
    ofstream Sc_zdata (path + "Scaled_zdata.txt", ofstream::app); // Scaled zdata is for -1 to 1 normalization of GL window
    ofstream Time_Course (path + "Time_Course.txt", ofstream::app);
    
    ofstream F_pulling_vs_time (path + "F_pulling_vs_time.txt", ofstream::app);
    ofstream F_pushing_vs_time (path + "F_pushing_vs_time.txt", ofstream::app);
    ofstream F_total_vs_time (path + "F_total_vs_time.txt", ofstream::app);
    
    ofstream Gel_deformation_avg_in_time (path + "Gel_deformation_avg_in_time.txt", ofstream::app);
    ofstream Edge_deformation_avg_in_time (path + "Edge_deformation_avg_in_time.txt", ofstream::app);
    ofstream dz_elas_avg_pullers_vs_time (path + "dz_elas_avg_pullers_vs_time.txt", ofstream::app);
    ofstream dz_elas_avg_inner_pullers_vs_time (path + "dz_elas_avg_inner_pullers_vs_time.txt", ofstream::app);
    ofstream dz_elas_avg_outer_pullers_vs_time (path + "dz_elas_avg_outer_pullers_vs_time.txt", ofstream::app);
    ofstream dz_elas_avg_pushers_vs_time (path + "dz_elas_avg_pushers_vs_time.txt", ofstream::app);
    
    ofstream k_off_data (path + "k_off.txt", ofstream::app);
    ofstream k_on_data (path + "k_on.txt", ofstream::app);

    ofstream Force_Dist (path + "Force_Dist.txt");
 //   ofstream Velocity_Dist (path + "Velocity_Dist.txt");
    ofstream Gel_Deformation_Dist (path + "Gel_Deformation_Dist.txt");
    ofstream Force_Dist_Symmetrized (path + "Force_Dist_Symmetrized.txt");
    ofstream Force_Dist_Row_5 (path + "Force_Dist_Row_5.txt");
    ofstream Force_Dist_Row_3 (path + "Force_Dist_Row_3.txt");
    ofstream Data (path + "Data.txt", ofstream::app);
    ofstream F_ext_vs_v_memb (path + "F_ext_vs_v_memb.txt", ofstream::app);
    ofstream F_ext_vs_v_fil (path + "F_ext_vs_v_fil.txt", ofstream::app);
    
    Paramfile << " PID: " << rank << ", k_bend: " <<  k_bend << ", k_elas: " <<  k_elas << ", t_eq: " << t_equilibrium
    << ", A_pull: " << A_pull << ", B_pull: " << B_pull << ", A_push: " << A_push << ", B_push: " << B_push << ", F_ext: " << F_ext << endl;
   
    //random_device seed;
    //mt19937 rgen(seed());
    mt19937 rgen(seed_1);
    
    //srand(time(NULL));
    srand(seed_1);

    /////////////////////////////////////////////
    //=========================================//
    /////////////////////////////////////////////
    
    double Ntstep_after_equilibrium = Ntstep - t_equilibrium/dt;
    const int n = limit - int(t_equilibrium); // Number of means to use in error bar calculations.
    
    double dz_bend_avg_each_fil_this_j [m][m][n];
    double dz_elas_avg_each_fil_this_j [m][m][n];
    double F_avg_each_fil_this_j [m][m][n];
    double F_symmetrized_avg_each_fil_this_j [m][m][n];
    double F_avg_pull_this_j [n], F_avg_push_this_j [n], F_avg_total_this_j [n];
    double F_avg_inner_pull_this_j [n], F_avg_outer_pull_this_j [n];
    double v_memb_avg_this_j [n];
//    double v_fil_this_j [m][m][n];
//    double Num_subunits_this_j [m][m][n];
//    double DeltaN_in_this_j [m][m][n];
    
    for (int row = 0; row < m; row++)
    {
        for (int col = 0; col < m; col++)
        {
            for (int j = 0; j < n; j++)
            {
                dz_bend_avg_each_fil_this_j [m][m][n] = 0;
                dz_elas_avg_each_fil_this_j [m][m][n] = 0;
                F_avg_each_fil_this_j [m][m][n] = 0;
                F_symmetrized_avg_each_fil_this_j [m][m][n] = 0;
//                v_fil_this_j [m][m][n] = 0;
//                Num_subunits_this_j [m][m][n] = 0;
//                DeltaN_in_this_j [m][m][n] = 0;
            }
        }
    }
    
    for (int j = 0; j < n; j++)
    {
        F_avg_pull_this_j [n] = 0;
        F_avg_push_this_j [n] = 0;
        F_avg_total_this_j [n] = 0;
        F_avg_inner_pull_this_j [n] = 0;
        F_avg_outer_pull_this_j [n] = 0;
        v_memb_avg_this_j [n] = 0;
    }


    for (int row = 0; row < m; row++)
    {
      for (int col = 0; col < m; col++)
      {
          
         double randh = ((double) rand() / (RAND_MAX));
         initial_staggered_length [row][col] = randh*delta;
          
         if ( (row >= ((m - m_p)/2) && row <= (m + m_p)/2 - 1) && (col >= ((m - m_p)/2) && col <= (m + m_p)/2 - 1) )
         {
            k_on [row][col] = 11.6; 
            k_off [row][col] = 1.4; 
            Poly_rate [row][col] = k_on[row][col]*C_0; 
            Depoly_rate [row][col] = k_off[row][col]; 
            Poly_prob [row][col] = Poly_rate[row][col]*dt; 
            Depoly_prob [row][col] = Depoly_rate[row][col]*dt; 
            Num_subunits [row][col] = 1;
            dz_elas[row][col] = 0; dz_bend[row][col] = 0;
            filament_height [row][col] = Num_subunits[row][col]*delta + initial_staggered_length[row][col]; 
         }
         else 
         {
            k_on [row][col] = 11.6; 
            k_off [row][col] = 1.4; 
            Poly_rate [row][col] = k_on[row][col]*C_0; 
            Depoly_rate [row][col] = k_off[row][col]; 
            Poly_prob [row][col] = Poly_rate[row][col]*dt; 
            Depoly_prob [row][col] = Depoly_rate[row][col]*dt; 
            Num_subunits [row][col] = 1;
            dz_elas[row][col] = 0; dz_bend[row][col] = 0;
            filament_height [row][col] = Num_subunits[row][col]*delta + initial_staggered_length[row][col];
         }
       }
    
    }
    
    for (int row = 0; row < m; row++)
    {
       for (int col = 0; col < m; col++)
       {
          ofstream Num_subunits_data ((path + "Num_subunits_" + to_string(row) + "_" + to_string(col) + ".txt"), ofstream::trunc); 
          ofstream dz_elas_data ((path + "dz_elas_" + to_string(row) + "_" + to_string(col) + ".txt"), ofstream::trunc);
          ofstream dz_bend_data ((path + "dz_bend_" + to_string(row) + "_" + to_string(col) + ".txt"), ofstream::trunc);
          ofstream Force_vs_time_data ((path + "Force_vs_time_" + to_string(row) + "_" + to_string(col) + ".txt"), ofstream::trunc);
          ofstream r_data ((path + "r_" + to_string(row) + "_" + to_string(col) + ".txt"), ofstream::trunc);
       }
     }
    
    uniform_int_distribution<unsigned long long int> dist_membrane(0,pow(10,9));
    uniform_int_distribution<unsigned long long int> dist_poly(0,pow(10,9));
    uniform_int_distribution<unsigned long long int> dist_depoly(0,pow(10,9));
    uniform_int_distribution<unsigned long long int> dist_bend(0,pow(10,9));
    uniform_int_distribution<unsigned long long int> dist_elas(0,pow(10,9));
    
    /////////////////////////////////////////////
    //=========================================//
    /////////////////////////////////////////////
    
    for (size_t j=0; j<limit; j++) // This is the "error bar" for loop!
    {
        
       for (size_t i=0; i<(Ntstep/limit); i++)
       {
           
          double t = i*dt + j*(Ntstep/limit)*dt;
          uni_rand_membrane = (double(dist_membrane(rgen))/double(pow(10,9))) - 0.5;
          g_memb = sc*uni_rand_membrane;
           
    /////////////////////////////////////////////
    //=========================================//
    /////////////////////////////////////////////
         
          if (i % ss == 0)
          {
             for (int row = 0; row < m; row++)
             {
                for (int col = 0; col < m; col++)
                {
                   uni_rand_poly [row][col] = double(dist_poly(rgen))/(double(pow(10,9)));
                   uni_rand_depoly [row][col] = double(dist_depoly(rgen))/(double(pow(10,9)));
                   
                   double r_tip = r [row][col];
                   double pot = Potential(r_tip,row,col); 
                   double potmi = Potential(r_tip - delta,row,col); 
                   double potpl = Potential(r_tip + delta,row,col); 
                   
                   k_on [row][col] = 11.6*exp(-beta*( potmi - pot ));
                   k_off [row][col] = 1.4*exp(-beta*( potpl - pot ));

                   if ( potmi < pot )
                   {
                       k_on [row][col] = 11.6;
                   }

                   if ( potpl < pot )
                   {
                       k_off [row][col] = 1.4;
                   }
                   
                   Poly_prob [row][col] = k_on[row][col]*C_0*dt;
                   Depoly_prob [row][col] = k_off[row][col]*dt;
                   
                   if (ss*Poly_prob[row][col] > uni_rand_poly[row][col])
                   {
                      Num_subunits[row][col]++; Num_poly_events[row][col]++;
                   }
                   if (ss*Depoly_prob[row][col] > uni_rand_depoly[row][col])
                   {
                      Num_subunits[row][col]--; Num_depoly_events[row][col]++;
                   }
                    
                }
             }
          }
            
            
          for (int row = 0; row < m; row++)
          {
             for (int col = 0; col < m; col++)
             {
                
                uni_rand_bend [row][col] = (double(dist_bend(rgen))/double(pow(10,9))) - 0.5;
                g_bend [row][col] = sc*uni_rand_bend [row][col];

                uni_rand_elas [row][col] = (double(dist_elas(rgen))/double(pow(10,9))) - 0.5;
                g_elas [row][col] = sc*uni_rand_elas [row][col];
                 
                filament_height [row][col] = Num_subunits[row][col]*delta + dz_elas[row][col] + dz_bend[row][col] + initial_staggered_length[row][col];
                r [row][col] = z - filament_height [row][col];

                double r_tip = r [row][col];
                F [row][col] = Force(r_tip,row,col);
                 
                if ((row >= ((m - m_p)/2) && row <= (m + m_p)/2 - 1) && (col >= ((m - m_p)/2) && col <= (m + m_p)/2 - 1)) {

                     F_total_pull_to_draw += F[row][col];
                 }
                 
                else { F_total_push_to_draw += F[row][col]; }
                 
                 
                if (t>t_equilibrium && Break == 0)
                {
                    
                    if ((row >= ((m - m_p)/2) && row <= (m + m_p)/2 - 1) && (col >= ((m - m_p)/2) && col <= (m + m_p)/2 - 1)) {
                        
                        F_total_pull += F[row][col];
                        F_total_pull_each_j += F[row][col];
                        
                        if ((row >= ((m - m_p)/2) + 1 && row <= (m + m_p)/2 - 2) && (col >= ((m - m_p)/2) + 1 && col <= (m + m_p)/2 - 2)) {
                                F_total_inner_pull += F[row][col];
                                F_total_inner_pull_each_j += F[row][col];
                             }
                        else { F_total_outer_pull += F[row][col]; F_total_outer_pull_each_j += F[row][col]; }
                        
                    }
                    
                    else { F_total_push += F[row][col]; F_total_push_each_j += F[row][col]; }
                    
                    F_total = F_total_pull + F_total_push; // This is zero until t_eq!
                 
                    F_total_each_fil [row][col] += F [row][col];
                   
                    F_total_each_fil_each_j [row][col] += F [row][col];

                    F_symmetrized_each_fil [row][col] += (0.125)*( F[row][col] + F[11-row][col] + F[row][11-col] + F[11-row][11-col] + F[col][row] + F[11-col][row] + F[col][11-row] + F[11-col][11-row] );
                 
                    F_symmetrized_each_fil_each_j [row][col] += (0.125)*( F[row][col] + F[11-row][col] + F[row][11-col] + F[11-row][11-col] + F[col][row] + F[11-col][row] + F[col][11-row] + F[11-col][11-row] );
                    
                }
                 
                
                if (i % Ntstep_to_store_data == 0)
                    {
                        ofstream Num_subunits_data ((path + "Num_subunits_" + to_string(row) + "_" + to_string(col) + ".txt"), ofstream::app);
                        Num_subunits_data << Num_subunits[row][col] << endl; // keep them like this for opengl resource directory!
                        ofstream r_data ((path + "r_" + to_string(row) + "_" + to_string(col) + ".txt"), ofstream::app);
                        r_data << r [row][col] << endl;
                        ofstream Force_vs_time_data ((path + "Force_vs_time_" + to_string(row) + "_" + to_string(col) + ".txt"), ofstream::app);
                        Force_vs_time_data << F[row][col] << endl;
                    }
                 
               } 
            }

        /////////////////////////////////////////////
        //=========================================//
        /////////////////////////////////////////////

        for (int row = 0; row < m; row++)
            {
                for (int col = 0; col < m; col++)
                {
                     if ( (row == 0) && (col == 0) ) //upper left corner 
                     { 
                        S[row][col] = (dz_elas[row][col+1] - dz_elas[row][col]) + (dz_elas[row+1][col] - dz_elas[row][col]); 
                     }
                     else if ( (row == 0) && (col == m-1) ) //upper right corner
                     { 
                        S[row][col] = ( dz_elas[row][col-1] - dz_elas[row][col]) + (dz_elas[row+1][col] - dz_elas[row][col]); 
                     }
                     else if ( (row == m-1) && (col == 0) ) //lower left corner
                     { 
                        S[row][col] = ( dz_elas[row-1][col] - dz_elas[row][col]) + (dz_elas[row][col+1] - dz_elas[row][col]); 
                     }
                     else if ( (row == m-1) && (col == m-1) ) //lower right corner
                     { 
                        S[row][col] = (dz_elas[row][col-1] - dz_elas[row][col]) + (dz_elas[row-1][col] - dz_elas[row][col]); 
                     }
                     else if ( (row == 0) && (col > 0 && col < m-1) ) 
                     { 
                        S[row][col] = (dz_elas[row][col+1] - dz_elas[row][col]) + (dz_elas[row][col-1] - dz_elas[row][col]) + (dz_elas[row+1][col] - dz_elas[row][col]); 
                     }
                     else if ( (col == 0) && (row > 0 && row < m-1) ) 
                     { 
                        S[row][col] = (dz_elas[row-1][col] - dz_elas[row][col]) + (dz_elas[row+1][col] - dz_elas[row][col]) + (dz_elas[row][col+1] - dz_elas[row][col]); 
                     }
                     else if ( (row == m-1) && (col > 0 && col < m-1) ) 
                     { 
                        S[row][col] = (dz_elas[row][col+1] - dz_elas[row][col]) + (dz_elas[row][col-1] - dz_elas[row][col]) + (dz_elas[row-1][col] - dz_elas[row][col]); 
                     }
                     else if ( (col == m-1) && (row > 0 && row < m-1) ) 
                     { 
                        S[row][col] = (dz_elas[row-1][col] - dz_elas[row][col]) + (dz_elas[row+1][col] - dz_elas[row][col]) + (dz_elas[row][col-1] - dz_elas[row][col]); 
                     }
                     else 
                     { 
                        S[row][col] = (dz_elas[row-1][col] - dz_elas[row][col]) + (dz_elas[row][col+1] - dz_elas[row][col]) + (dz_elas[row+1][col] - dz_elas[row][col]) + (dz_elas[row][col-1] - dz_elas[row][col]); 
                     }
                
                }
            }
           
            /////////////////////////////////////////////
            //=========================================//
            /////////////////////////////////////////////
            
            for (int row = 0; row < m; row++)
            {
               for (int col = 0; col < m; col++)
               {

                  dz_bend [row][col] = dz_bend [row][col] + g_bend[row][col]*sqrt(D_bend) + mu_bend*dt*(-F[row][col] - k_bend*(dz_bend[row][col]));
                   
                  if (dz_bend [row][col] > 9.77) // L = 54nm; Theta = 35 degree.
                  {
                     dz_bend [row][col] = 9.77;
                  }
                  if (dz_bend [row][col] < -44.23)
                  {   
                     dz_bend [row][col] = -44.23;
                  }
                   
                   dz_elas [row][col] = dz_elas [row][col] + g_elas[row][col]*sqrt(D_elas) + mu_elas*dt*(-F[row][col] + k_elas*S[row][col] - k_elas*dz_elas[row][col]); // As an item on the list of further improvements on model: we can distinguish between these two k_elas values, base to neighbor base, and base to actin gel body connections.
                   
                  if (t>t_equilibrium && Break == 0)
                  {
                      dz_bend_total_each_fil [row][col] += dz_bend [row][col];
                      dz_bend_total_each_fil_each_j [row][col] += dz_bend [row][col];
                      
                      dz_elas_total_each_fil [row][col] += dz_elas [row][col];
                      dz_elas_total_each_fil_each_j [row][col] += dz_elas [row][col];
                  }
                   
               }
            }
           
           
            if(i % Ntstep_to_store_data == 0)
            {
               for (int row = 0; row < m; row++)
               {
                  for (int col = 0; col < m; col++)
                  {
                      if ( (row >= ((m - m_p)/2) && row <= (m + m_p)/2 - 1) && (col >= ((m - m_p)/2) && col <= (m + m_p)/2 - 1) )
                      {
                          dummy1 += dz_elas [row][col];
                          
                          if ((row >= ((m - m_p)/2) + 1 && row <= (m + m_p)/2 - 2) && (col >= ((m - m_p)/2) + 1 && col <= (m + m_p)/2 - 2)) {dummy2 += dz_elas [row][col];}
                          else {dummy3 += dz_elas [row][col];}
                      }
                      
                      else {dummy4 += dz_elas [row][col];}
                      
                     ofstream dz_elas_data ((path + "dz_elas_" + to_string(row) + "_" + to_string(col) + ".txt"), ofstream::app); 
                     dz_elas_data << dz_elas[row][col] << endl; 
                     ofstream dz_bend_data ((path + "dz_bend_" + to_string(row) + "_" + to_string(col) + ".txt"), ofstream::app); 
                     dz_bend_data << dz_bend[row][col] << endl;
                   }
                }
                
                dz_elas_avg_pullers_at_each_timestep = dummy1/(n_central);
                dz_elas_avg_inner_pullers_at_each_timestep = dummy2/(n_inner_central);
                dz_elas_avg_outer_pullers_at_each_timestep = dummy3/(n_outer_central);
                dz_elas_avg_pushers_at_each_timestep = dummy4/(n_filaments - n_central);

                dummy1 = 0;
                dummy2 = 0;
                dummy3 = 0;
                dummy4 = 0;
            }
           
        /////////////////////////////////////////////
        //=========================================//
        /////////////////////////////////////////////
           
           F_total_to_draw = F_total_pull_to_draw + F_total_push_to_draw; // This is non-zero right away from the beginning!
           
           z = z + g_memb*sqrt(D_memb) + mu_memb*dt*(F_ext + F_total_to_draw);
           
           v_memb = (g_memb*sqrt(D_memb) + mu_memb*dt*(F_ext + F_total_to_draw))/dt;
           
           if (t>t_equilibrium && Break == 0)
           {
               v_memb_total += v_memb;
               v_memb_total_each_j += v_memb;
           }
           
           sum_Pr += (C_0*k_on[5][5]-k_off[5][5]);
           
        /////////////////////////////////////////////
        //=========================================//
        /////////////////////////////////////////////
           
           if (i % Ntstep_to_store_data == 0){
               
               Sc_zdata << z/opengl_window_scaling_metric << endl;
               Time_Course << t << ' ' << z << endl;
               
               k_off_data << k_off[5][5] << endl;
               k_on_data << k_on[5][5] << endl;
               
               F_pulling_vs_time << t << ' ' << F_total_pull_to_draw << endl;
               F_pushing_vs_time << t << ' ' << F_total_push_to_draw << endl;
               F_total_vs_time << t << ' ' << F_total_to_draw << endl;
               
               Gel_deformation_avg_in_time << t << ' ' << dz_elas_avg_pullers_at_each_timestep - dz_elas_avg_pushers_at_each_timestep << endl;
               Edge_deformation_avg_in_time << t << ' ' << dz_elas_avg_outer_pullers_at_each_timestep - dz_elas_avg_inner_pullers_at_each_timestep << endl;
               
               dz_elas_avg_pullers_vs_time << t << ' ' << dz_elas_avg_pullers_at_each_timestep << endl;
               dz_elas_avg_inner_pullers_vs_time << t << ' ' << dz_elas_avg_inner_pullers_at_each_timestep << endl;
               dz_elas_avg_outer_pullers_vs_time << t << ' ' << dz_elas_avg_outer_pullers_at_each_timestep << endl;

               dz_elas_avg_pushers_vs_time << t << ' ' << dz_elas_avg_pushers_at_each_timestep << endl;
               
               if ( (t >= 2 && dz_elas_avg_pullers_at_each_timestep < 0) && Break == 0 ) {Break = 1; t_break = t;}
               
           }
        
            F_total_pull_to_draw = 0;
            F_total_push_to_draw = 0;
            F_total_to_draw = 0;
            
        } // closing the time iteration loop
        
    ////////////////////////////////////////////////
    //=====// closing the time iteration loop=====//
    ////////////////////////////////////////////////
        
    if (j >= int(t_equilibrium)){
        
        int jj = j - int(t_equilibrium);
        
        for (int row = 0; row < m; row++)
        {
            for (int col = 0; col < m; col++)
            {
                dz_bend_avg_each_fil_this_j [row][col][jj] = dz_bend_total_each_fil_each_j [row][col]/(Ntstep/limit);
                dz_bend_total_each_fil_each_j [row][col] = 0;
            
                dz_elas_avg_each_fil_this_j [row][col][jj] = dz_elas_total_each_fil_each_j [row][col]/(Ntstep/limit);
                dz_elas_total_each_fil_each_j [row][col] = 0;
            
                F_avg_each_fil_this_j [row][col][jj] = F_total_each_fil_each_j [row][col]/(Ntstep/limit);
                F_total_each_fil_each_j [row][col] = 0;
            
                F_symmetrized_avg_each_fil_this_j [row][col][jj] = F_symmetrized_each_fil_each_j [row][col]/(Ntstep/limit);
                F_symmetrized_each_fil_each_j [row][col] = 0;
                
//                DeltaN_in_this_j [row][col][jj] = Num_subunits[row][col] - Num_subunits_this_j [row][col][jj];
//                v_fil_this_j [row][col][jj] = (DeltaN_in_this_j [row][col][jj]*delta)/((Ntstep/limit)*dt);
//                Num_subunits_this_j [row][col][jj] = Num_subunits[row][col];

            }
        }
        
            F_avg_total_this_j [jj] = (F_total_pull_each_j+F_total_push_each_j)/(Ntstep/limit);

            F_avg_pull_this_j [jj] = F_total_pull_each_j/(Ntstep/limit);
            F_total_pull_each_j = 0;
        
            F_avg_inner_pull_this_j [jj] = F_total_inner_pull_each_j/(Ntstep/limit);
            F_total_inner_pull_each_j = 0;
        
            F_avg_outer_pull_this_j [jj] = F_total_outer_pull_each_j/(Ntstep/limit);
            F_total_outer_pull_each_j = 0;
        
            F_avg_push_this_j [jj] = F_total_push_each_j/(Ntstep/limit);
            F_total_push_each_j = 0;
        
            v_memb_avg_this_j [jj] = v_memb_total_each_j/(Ntstep/limit);
            v_memb_total_each_j = 0;
        
    }
        
        cout << "Current j is: " << j << ' ' << " out of: " << limit - 1 << ' ' << ", and Num_subunits[5][5] = " << Num_subunits[5][5] << ' ' << " and z = " << z << endl;
        
} // closing the entire time loop
    
    /////////////////////////////////////////////
    //==============TIME LOOP ENDS=============//
    /////////////////////////////////////////////

    if (t_break > t_equilibrium){ Ntstep_after_equilibrium = (t_break-t_equilibrium)/dt; }
    
    for (int row = 0; row < m; row++)
    {
        for (int col = 0; col < m; col++)
        {
            
            dz_bend_avg_each_fil [row][col] = dz_bend_total_each_fil [row][col]/(Ntstep_after_equilibrium);
            
            dz_elas_avg_each_fil [row][col] = dz_elas_total_each_fil [row][col]/(Ntstep_after_equilibrium);

            F_avg_each_fil [row][col] = F_total_each_fil [row][col]/(Ntstep_after_equilibrium);

            F_symmetrized_avg_each_fil [row][col] = F_symmetrized_each_fil [row][col]/(Ntstep_after_equilibrium);
            
//            v_fil_avg [row][col] = (Num_subunits [row][col]*delta)/(Ntstep*dt);
            
            if ( (row >= ((m - m_p)/2) && row <= (m + m_p)/2 - 1) && (col >= ((m - m_p)/2) && col <= (m + m_p)/2 - 1) )
            {
                Num_subunits_total_pullers += Num_subunits [row][col];
                dz_bend_avg_total_pullers += dz_bend_avg_each_fil [row][col];
                dz_elas_avg_total_pullers += dz_elas_avg_each_fil [row][col];
            }
            else
            {
                Num_subunits_total_pushers += Num_subunits [row][col];
                dz_bend_avg_total_pushers += dz_bend_avg_each_fil [row][col];
                dz_elas_avg_total_pushers += dz_elas_avg_each_fil [row][col];
            }
            
        }
    }
    
    F_pull_avg = F_total_pull/Ntstep_after_equilibrium;
    F_inner_pull_avg = F_total_inner_pull/Ntstep_after_equilibrium;
    F_outer_pull_avg = F_total_outer_pull/Ntstep_after_equilibrium;
    
    F_push_avg = F_total_push/Ntstep_after_equilibrium;
    F_total_avg = F_total/Ntstep_after_equilibrium;
    
    v_memb_avg = v_memb_total/Ntstep_after_equilibrium;
    
    Num_subunits_avg_pullers = Num_subunits_total_pullers/(n_central);
    Num_subunits_avg_pushers = Num_subunits_total_pushers/(n_filaments - n_central);
    
    dz_bend_avg_array_pullers = dz_bend_avg_total_pullers/(n_central);
    dz_elas_avg_array_pullers = dz_elas_avg_total_pullers/(n_central);
    
    dz_bend_avg_array_pushers = dz_bend_avg_total_pushers/(n_filaments - n_central);
    dz_elas_avg_array_pushers = dz_elas_avg_total_pushers/(n_filaments - n_central);
    
    ///////////////////////////////////////////////
    //==============ERROR BAR LOOP===============//
    ///////////////////////////////////////////////
    
    double dummy_sum_1 = 0;
    double dummy_sum_2 = 0;
    double dummy_sum_3 = 0;
    double dummy_sum_4 = 0;
    double dummy_sum_5 = 0;
    double dummy_sum_6 = 0;
    double dummy_sum_7 = 0;
    double dummy_sum_8 = 0;
    double dummy_sum_9 = 0;
    double dummy_sum_10 = 0;
//    double dummy_sum_11 = 0;

    
    for (int row = 0; row < m; row++)
    {
        for (int col = 0; col < m; col++)
        {
                
            for (int j =0; j < n; j++)
            {
                dummy_sum_1 += (F_avg_each_fil_this_j [row][col][j] - F_avg_each_fil [row][col])*(F_avg_each_fil_this_j [row][col][j] - F_avg_each_fil [row][col]);
                dummy_sum_2 += (F_symmetrized_avg_each_fil_this_j [row][col][j] - F_symmetrized_avg_each_fil [row][col])*(F_symmetrized_avg_each_fil_this_j [row][col][j] - F_symmetrized_avg_each_fil [row][col]);
                dummy_sum_3 += (dz_elas_avg_each_fil_this_j [row][col][j] - dz_elas_avg_each_fil [row][col])*(dz_elas_avg_each_fil_this_j [row][col][j] - dz_elas_avg_each_fil [row][col]);
                dummy_sum_4 += (dz_bend_avg_each_fil_this_j [row][col][j] - dz_bend_avg_each_fil [row][col])*(dz_bend_avg_each_fil_this_j [row][col][j] - dz_bend_avg_each_fil [row][col]);
      //          dummy_sum_11 += (v_fil_this_j [row][col][j] - v_fil_avg [row][col])*((v_fil_this_j [row][col][j] - v_fil_avg [row][col]);

            }
                
                Error_bar_Force [row][col] = sqrt( dummy_sum_1/(n*(n-1)) );
                dummy_sum_1 = 0;
                
                Error_bar_symmetrized_Force [row][col] = sqrt( dummy_sum_2/(n*(n-1)) );
                dummy_sum_2 = 0;
                
                Error_bar_dz_elas [row][col] = sqrt( dummy_sum_3/(n*(n-1)) );
                dummy_sum_3 = 0;
                
                Error_bar_dz_bend [row][col] = sqrt( dummy_sum_4/(n*(n-1)) );
                dummy_sum_4 = 0;
                                                                                     
//                Error_bar_v_fil [row][col] = sqrt( dummy_sum_11/(n*(n-1)) );
//                dummy_sum_11 = 0;
        }
    }
    
    for (int j =0; j < n; j++)
    {
        dummy_sum_5 += (F_avg_pull_this_j [j] - F_pull_avg)*(F_avg_pull_this_j [j] - F_pull_avg);
        dummy_sum_6 += (F_avg_inner_pull_this_j [j] - F_inner_pull_avg)*(F_avg_inner_pull_this_j [j] - F_inner_pull_avg);
        dummy_sum_7 += (F_avg_outer_pull_this_j [j] - F_outer_pull_avg)*(F_avg_outer_pull_this_j [j] - F_outer_pull_avg);
        dummy_sum_8 += (F_avg_push_this_j [j] - F_push_avg)*(F_avg_push_this_j [j] - F_push_avg);
        dummy_sum_9 += (F_avg_total_this_j [j] - F_total_avg)*(F_avg_total_this_j [j] - F_total_avg);
        dummy_sum_10 += (v_memb_avg_this_j [j] - v_memb_avg)*(v_memb_avg_this_j [j] - v_memb_avg);
    }
    
    Error_bar_Pulling_Force = sqrt( dummy_sum_5/(n*(n-1)) );
    Error_bar_Inner_Pulling_Force = sqrt( dummy_sum_6/(n*(n-1)) );
    Error_bar_Outer_Pulling_Force = sqrt( dummy_sum_7/(n*(n-1)) );
    Error_bar_Pushing_Force = sqrt( dummy_sum_8/(n*(n-1)) );
    Error_bar_Total_Force = sqrt( dummy_sum_9/(n*(n-1)) );
    Error_bar_v_memb = sqrt( dummy_sum_10/(n*(n-1)) );
    
    /////////////////////////////////////////////
    //=========================================//
    /////////////////////////////////////////////
    
    for (int row = 0; row < m; row++)
    {
       for (int col = 0; col < m; col++)
       {

           Force_Dist << row+1 << ' ' << col+1 << ' ' << F_avg_each_fil [row][col] << ' ' << Error_bar_Force [row][col] << endl;
           Force_Dist_Symmetrized << row+1 << ' ' << col+1 << ' ' << F_symmetrized_avg_each_fil [row][col] << ' ' << Error_bar_symmetrized_Force [row][col] << endl;
 //          Velocity_Dist << row+1 << ' ' << col+1 << ' ' << v_fil_avg [row][col] << ' ' << Error_bar_v_fil [row][col] << endl;

           Gel_Deformation_Dist << row+1 << ' ' << col+1 << ' ' << dz_elas_avg_each_fil [row][col] << ' ' << Error_bar_dz_elas [row][col] << endl;
           
           if (row==5){ Force_Dist_Row_5 << col+1 << ' ' << F_symmetrized_avg_each_fil [row][col] << ' ' << Error_bar_symmetrized_Force [row][col] << endl; }
           if (row==3){ Force_Dist_Row_3 << col+1 << ' ' << F_symmetrized_avg_each_fil [row][col] << ' ' << Error_bar_symmetrized_Force [row][col] << endl; }
       }
    }
    
    F_ext_vs_v_memb << F_ext << ' ' << v_memb_avg << ' ' << Error_bar_v_memb << endl;
    
    /////////////////////////////////////////////
    //=========================================//
    /////////////////////////////////////////////
    
    Data << "Total simulation time: " << Ntstep*dt << " sec." << endl << "Equilibrium time: " << t_equilibrium << " sec." << endl
    << "Obstacle breaks apart at: " << t_break << " sec." << endl;
    
    Data << "A_pull: " << A_pull << " pN." << endl << "External load: " << F_ext << " pN." << endl << "Pulling Force: " << F_pull_avg << " +- " <<
    Error_bar_Pulling_Force << " pN." << endl << "Pushing Force: " << F_push_avg << " +- " << Error_bar_Pushing_Force << " pN." << endl << "Total Force: " <<
    F_total_avg << " +- " << Error_bar_Total_Force << " pN." << endl << "Inner Pulling Force: " << F_inner_pull_avg << " +- " << Error_bar_Inner_Pulling_Force <<
    " pN." << endl << "Outer Pulling Force: " << F_outer_pull_avg << " +- " << Error_bar_Outer_Pulling_Force << " pN." << endl << "Membrane Average Velocity: " <<
    v_memb_avg << " +- " << Error_bar_v_memb << " nm/sec." << endl << "<Pullers Polymerization Speed>: " << (Num_subunits_avg_pullers*delta)/(Ntstep*dt) << " nm/sec."
    << endl << "<Pushers Polymerization Speed>: " << (Num_subunits_avg_pushers*delta)/(Ntstep*dt) << " nm/sec." << endl;
    
    Data << "Actin gel surface deformation in pulling region: " << dz_elas_avg_array_pullers << " nm, and in pushing region: " << dz_elas_avg_array_pushers << " nm."
    << endl << "Filament tip region stretch and bend respectively: " << dz_bend_avg_array_pullers << " and " << dz_bend_avg_array_pushers << " nm." << endl;
    
    Data << "Safety check ratio (must be about 1): " << (double(Num_subunits [5][5])/double(Ntstep*dt))/(sum_Pr/Ntstep) << endl;

    for (int row = 0; row < m; row++)
    {
       for (int col = 0; col < m; col++)
       {
           Data << "Filament[" << row << "][" << col << "]: "
                << Num_subunits[row][col] << " subunits, with "
                << Num_poly_events[row][col] << " polymerization events and " 
                << Num_depoly_events[row][col] << " depolymerization events." << endl;
        }
    }
    
    Data << endl;
    
    for (int row = 0; row < m; row++)
    {   
       for (int col = 0; col < m; col++)
       {
           Data << "Filament[" << row << "][" << col << "] <dz_elas>: " <<
           dz_elas_avg_each_fil [row][col] << " +- " << Error_bar_dz_elas [row][col] << " nm." << endl;
        }
     }
    
    Data << endl;

    for (int row = 0; row < m; row++)
    {
        for (int col = 0; col < m; col++)
        {
            Data << "Filament[" << row << "][" << col << "] <dz_bend>: " <<
            dz_bend_avg_each_fil [row][col] << " +- " << Error_bar_dz_bend [row][col] << " nm." << endl;
        }
    }
    
    /////////////////////////////////////////////
    //=========================================//
    /////////////////////////////////////////////
    
    cout << "F_ext is: " << F_ext << endl;
    
    MPI_Finalize();
    return 0;
}
