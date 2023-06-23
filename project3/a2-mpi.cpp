#include <iostream>
#include <iomanip>
#include <chrono>
#include <cmath>
#include <mpi.h>
#include "a2-helpers.hpp"

using namespace std;

int main(int argc, char **argv)
{    
    int max_iterations = 1000;
    double epsilon = 1.0e-3;
    bool verify = true, print_config = false;

    // default values for M rows and N columns
    int N = 12;
    int M = 12;
    process_input(argc, argv, N, M, max_iterations, epsilon, verify, print_config);

    if ( print_config )
        std::cout << "Configuration: m: " << M << ", n: " << N << ", max-iterations: " << max_iterations << ", epsilon: " << epsilon << std::endl;


    // START: the main part of the code that needs to use MPI/OpenMP timing routines 
    //  MPI hint: remember to initialize MPI first 
    MPI_Init(&argc, &argv);

    int numprocs, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double t1 = MPI_Wtime();

    int i, j;
    double diffnorm;
    int iteration_count = 0;
    int internal_M = M / numprocs;

   

    int internal_M_start, internal_M_end;
    // initialzie the internal_M for  M%numprocs != 0 
    if (M% numprocs != 0){
        if (rank == numprocs - 1){

                internal_M = internal_M + (M% numprocs);
        }
        else { // this is normal process
            internal_M_start = rank * internal_M; // start is the same place
            internal_M_end = internal_M_start + internal_M; // we shift  

        }
    }
    else{ // its easily divided by amount of threads
        internal_M_start = rank * internal_M;
        internal_M_end = internal_M_start + internal_M;
    }
        
   
    Mat U(internal_M + 2, N);
    Mat W(internal_M + 2, N);
    // Init & Boundary
    
    for (i = 0; i < internal_M + 2; ++i) {
        for (j = 0; j < N; ++j) {
            W[i][j] = U[i][j] = 0.0;
        }

        W[i][0] = U[i][0] = 0.05; // left side
        W[i][N-1] = U[i][N-1] = 0.1; // right side
    }
    // had to be added ! top init
    if (rank == 0){
          for (j = 0; j < N; ++j) {
            W[1][j] = U[1][j] = 0.02; // top 
          }
        
    }

    else if (rank == numprocs - 1){
        for (j = 0; j < N; ++j) {
            W[internal_M][j] = U[internal_M ][j] = 0.2; // bottom
        }
    }
   
 

    iteration_count = 0;
    do
    {
        iteration_count++;
        diffnorm = 0.0;
        MPI_Request secondReq[2], firstReq[2];
        MPI_Status statusFirst[2], statusSecond[2];
        int firstCount = 0;
        int secondCount = 0;
        double local_diffnorm = 0;
        // SENDINT AND RECEIVING MESSEGES - HALO REGION EXCHANGE ///
        // if  not last
        if (rank < numprocs - 1) {
            // Send the last line
            MPI_Isend(&U[internal_M][0], N, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &firstReq[firstCount]); // Send to next process
            firstCount++;
            MPI_Irecv(&U[internal_M + 1][0], N, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &firstReq[firstCount]); // Non-blocking receive from the next process
            firstCount++;
        }
        // if not 0 
        if (rank > 0) {
            // Send the first line

            MPI_Irecv(&U[0][0], N, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD,  &secondReq[secondCount]); // Receive in-between row
            secondCount++;
            MPI_Isend(&U[1][0], N, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &secondReq[secondCount]); // Send first row
            secondCount++;
        }

        //  different waits for different ranks for safety
        if (rank == 0){
            MPI_Waitall(firstCount, firstReq, statusFirst);
        }
        else if (rank == numprocs - 1){
            MPI_Waitall(secondCount, secondReq, statusSecond);
        }
        else{
            MPI_Waitall(secondCount, secondReq, statusSecond);
            MPI_Waitall(firstCount, firstReq, statusFirst);
        }
        
        // rank 0 
        if (rank == 0){
            for (i = 2; i < internal_M +1  ; ++i) {   
                for (j = 1; j < N - 1; ++j)
                {
                    W[i][j] = (U[i][j + 1] + U[i][j - 1] + U[i + 1][j] + U[i - 1][j]) * 0.25;
                    local_diffnorm += (W[i][j] - U[i][j]) * (W[i][j] - U[i][j]);
                }
            }
            for (i = 2; i < internal_M +1 ; ++i){
                for (j = 1; j < N - 1; ++j){
                    U[i][j] = W[i][j];
                }
            }
        }
        // rank last
        else if (rank == numprocs - 1){
            for (i = 1; i < internal_M   ; ++i) {   
              for (j = 1; j < N - 1; ++j)
                {
                    W[i][j] = (U[i][j + 1] + U[i][j - 1] + U[i + 1][j] + U[i - 1][j]) * 0.25;
                    local_diffnorm += (W[i][j] - U[i][j]) * (W[i][j] - U[i][j]);
                }
            }
            for (i = 1; i < internal_M  ; ++i){
                for (j = 1; j < N - 1; ++j){
                    U[i][j] = W[i][j];
                }
            }

        }
        // rest
        else {
            for (i = 1; i < internal_M +1  ; ++i) {   
                for (j = 1; j < N - 1; ++j){
                    W[i][j] = (U[i][j + 1] + U[i][j - 1] + U[i + 1][j] + U[i - 1][j]) * 0.25;
                    local_diffnorm += (W[i][j] - U[i][j]) * (W[i][j] - U[i][j]);
                }
            }
            for (i = 1; i < internal_M + 1 ; ++i){
                for (j = 1; j < N - 1; ++j){
                    U[i][j] = W[i][j];
                }
            }
        }
        // WE GOTTA ALL REDUCE THE DIFFNORM !!
        MPI_Allreduce(&local_diffnorm, &diffnorm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        diffnorm = sqrt(diffnorm); // all processes need to know when to stop
        
    } while (epsilon <= diffnorm && iteration_count < max_iterations);

    Mat bigU(M, N); // create bigU on rank 0 to store the gathered data
    
    // Pieces of data can be different in size (last packet can be different)
    // lets adjut the data size
    int receive_displs[numprocs];  
    int receive_counts[numprocs]; 
    int sendCount = internal_M * N;
    for (int i = 0; i < numprocs; ++i) {
            receive_counts[i] = internal_M * N;
            receive_displs[i] = internal_M * N *i;
    }
    receive_counts[numprocs-1] += ((M% numprocs)*N);

    if (rank == 0){
        // time of processing data ( should be nearly the samge for all rands)
        double t2 = MPI_Wtime();
        double elapsed_time = t2 -t1;
    
        cout << "Elapsed time: "; 
        cout << std::fixed << std::setprecision(4) << elapsed_time << "Rank:"<<rank; // modify accordingly for MPI/OpenMP
        cout << " seconds, iterations: " << iteration_count << endl; 
        double t3 = MPI_Wtime();

        MPI_Gatherv(&U[1][0], sendCount, MPI_DOUBLE, &bigU[0][0], receive_counts, receive_displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        double t4 = MPI_Wtime();
        // no need for reducing time, because gatherv waits 
        double gathering_time = t4-t3;
        cout << "Gathering time: "; 
        cout << std::fixed << std::setprecision(4) << gathering_time <<std::endl;;

        if ( verify ) {
            Mat U_sequential(M, N); // init another matrix for the verification

            int iteration_count_seq = 0;
            heat2d_sequential(U_sequential, max_iterations, epsilon, iteration_count_seq); 
            // Here we need both results - from the sequential (U_sequential) and also from the OpenMP/MPI version, then we compare them with the compare(...) function 
            cout << "Verification: " << ( bigU.compare(U_sequential) && iteration_count == iteration_count_seq ? "OK" : "NOT OK") << std::endl;
                        // U_sequential.print();

        }
    }
  
    else{
        // send different data
        MPI_Gatherv(&U[1][0], sendCount, MPI_DOUBLE, &bigU[1][0], receive_counts, receive_displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    }
    MPI_Finalize();
    return 0;
}