#include <iostream>
#include <mpi.h>
#include <string>
#include <stdlib.h>
#include <vector>

#define MASTER 0

using namespace std;

void generate_rand_arr(int C, int rank, vector<double> &rand_arr){
	srand48(rank+C);
	for(int i=0;i<rand_arr.size();i++){
		rand_arr[i] = drand48();
	}
}

void sum_local_rand_arr(const vector<double> &rand_arr, double &sum_local){
	for(int i=0;i<rand_arr.size();i++){
		sum_local+=rand_arr[i];
	}
}

void sum_parallel(double &sum_local, const MPI_Comm &comm, int rank, int P){
	double recv_sum;
	MPI_Status stat;
	for(int i=1;i<P;i=i<<1){
		if((rank & i)!=0){
			MPI_Send(&sum_local, 1, MPI_DOUBLE, (rank^i), 111, comm);
			break;
		}
		else{
			MPI_Recv(&recv_sum, 1, MPI_DOUBLE, (rank^i), 111, comm, &stat);
			sum_local+=recv_sum;
		}
	}
}

int main(int argc, char* argv[]){
	int rank;
	int N = 0;
	int C = 0;
	int P;
	double sum_local = 0;

	MPI_Init(&argc, &argv);
	MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &P);
	
	//broadcast relevant info from master processor
	if(rank==MASTER){
		if(argc==3){
			N=stoi(argv[1]);
			C=stoi(argv[2]);
		}
		else if(argc==2){
			cout<<"warning missing C, please check if its not intentional"<<endl;
			cout<<"default value of C = 111 being used"<<endl;
			N=stoi(argv[1]);
			C=111;
			}
		else{
			cout<<"problem with the argument expected, please check N and C provided"<<endl;
			cout<<"default value of N=500000 and C=111 is used"<<endl;
			N=500000;
			C=111;
		}
	}
        MPI_Bcast(&N, 1, MPI_INT, MASTER, comm);
	MPI_Bcast(&C, 1, MPI_INT, MASTER, comm);
         
	MPI_Barrier(comm);
	double t0 = MPI_Wtime();

	vector<double> rand_arr;
	rand_arr.resize(N/P);
	//rand arr generated
	generate_rand_arr(C, rank, rand_arr);	
	//local sum
	sum_local_rand_arr(rand_arr, sum_local);	
        //sum parallely using point to point comm
	sum_parallel(sum_local, comm, rank, P);
	
	MPI_Barrier(comm);
	double t1 = MPI_Wtime();
	
	//printing the final sum by process 0
	if(rank==MASTER){
		cout<<"N = "<<N<<", P = "<<P<<", C = "<<C<<", S = "<<fixed<<sum_local<<endl;
		cout<<"Time = "<<fixed<< (t1-t0)<<endl;
	}

	MPI_Finalize();
	return 0;
}
