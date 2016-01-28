#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
using namespace std;
#include "mpi.h"

extern "C" {
#include "xtcio.h"
}
#include "theobald_rmsd.h"

#define PATHSIZE 200
#define TRAJECTORYLIST "trajlistshort"
#define DATADIR "/home/server.171.64.65.58/gboxer/kcenters/data/"
#define TRAJECTORYDIR "/home/server.171.64.65.58/gboxer/kcenters/trajectories/"
#define GENERATORDIR "/home/server.171.64.65.58/gboxer/kcenters/generators/"
#define INDEXFILE "atom_indices"
#define SNAPSHOTSPERFRAME 200

#define atoms 576

int main(int argc,char** argv){
	int process;
	int numprocess;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&process);
	MPI_Comm_size(MPI_COMM_WORLD,&numprocess);

	int* atomindicies;
	int natomindicies;
	ifstream indiciesfile(INDEXFILE);
	indiciesfile >> natomindicies;
	atomindicies=new int[natomindicies];
	for(int i=0;i<natomindicies;i++){
		indiciesfile >> atomindicies[i];
		atomindicies[i]--;
	}
	indiciesfile.close();

	int numtrajs;
	int* lengths;
	real** distance;
	int** state;
	char** trajnames;

	//only process 0 needs this array
	real* farthestdists;
	if(process==0)
		farthestdists=new real[numprocess];

	int trajlistsize;
	char tempstr[PATHSIZE];
	ifstream trajlist(TRAJECTORYLIST);
	trajlist >> trajlistsize;
	numtrajs=trajlistsize/numprocess+(trajlistsize%numprocess>process);

	lengths=new int[numtrajs];
	distance=new real*[numtrajs];
	state=new int*[numtrajs];
	trajnames=new char*[numtrajs];
	
	//Read in my trajectory assignments.  I take every one where myprocess=trajectory%numprocesses
	//Store the name of each ones trajectory file and data file
	for(int i=0;i<trajlistsize;i++){
		if(i%numprocess==process){
			trajnames[i/numprocess]=new char[PATHSIZE];
			trajlist >> trajnames[i/numprocess];
		}else
			trajlist >> tempstr;
	}
	trajlist.close();

	//Initialize the data (distance and state) for each snapshot.
	//For now we abuse the fact that each "frame" file has 200 snapshots but I probably shouldn't do this
	for(int i=0;i<numtrajs;i++){
		char trajfilename[PATHSIZE];
		sprintf(trajfilename,"%s%s",TRAJECTORYDIR,trajnames[i]);
		ifstream trajlist(trajfilename);
		char filename[PATHSIZE];
		int frames=0;
		while(trajlist>>filename)
			frames++;
		lengths[i]=frames*SNAPSHOTSPERFRAME+1; //This is shady
		distance[i]=new real[lengths[i]];
		state[i]=new int[lengths[i]];
		trajlist.close();
	}

	rvec generator[atoms];
	//choose the first generator, currently we just take the first snapshot from the first trajectory.  is there a better way to do this?
	if(process==0){
		char trajfilename[PATHSIZE];
		sprintf(trajfilename,"%s%s",TRAJECTORYDIR,trajnames[0]);
		ifstream trajlist(trajfilename);
		char filename[PATHSIZE];
		trajlist >> filename;
		//Stuff for reading in xtc files
		int fp = open_xtc(filename,"r");
		int natoms,step;
		real time;
		matrix box;
		rvec *x;
		real prec;
		bool bOK;
		read_first_xtc(fp,&natoms,&step,&time,box,&x,&prec,&bOK);

		//DEBUG:
		if(natoms!=atoms)
			cout << "DRAMATIC ERROR ENCOUNTERED: input frame has the wrong number of atoms" << endl;

		//write out the first generator to file
		//make the generator filename and open the file
		sprintf(filename,"%sgenerator.%i",GENERATORDIR,0);
		ofstream genfile(filename);
		//write out the trajectory, snapshot, generator number, and dist
		genfile << "Generator# " << 0 << endl;
		genfile << "Trajectory: " << trajnames[0] << endl;
		genfile << "Snapshot: " << 0 << endl;
		genfile << "Distance: " << "N/A" << endl;
		//write out the configuration
		genfile << "Configuration: " << endl;
		for(int atom=0;atom<atoms;atom++)
			genfile << x[atom][0] << " " << x[atom][1] << " " << x[atom][2] << endl;
		genfile.close();

		for(int atom=0;atom<atoms;atom++){
			generator[atom][0]=x[atom][0];
			generator[atom][1]=x[atom][1];
			generator[atom][2]=x[atom][2];
		}
		close_xtc(fp);
		delete[] x;
		trajlist.close();
	}
	//broadcast the first generator to everyone
	MPI_Bcast(generator,sizeof(rvec)*atoms,MPI_BYTE,0,MPI_COMM_WORLD);

	//The main loop
	bool done=0;
	int iteration=0;
	while(!done){
		real farthestdist=-1.0;
		rvec farthestconfig[atoms];
		int farthesttraj;
		int farthestsnapshot;
		//Iterate over trajectories
		for(int t=0;t<numtrajs;t++){
			char trajfilename[PATHSIZE];
			sprintf(trajfilename,"%s%s",TRAJECTORYDIR,trajnames[t]);
			ifstream trajlist(trajfilename);
			char filename[PATHSIZE];
			real prevtime=-50.0;
			int snapshotindex=0;
			//Iterate over frames in the trajectory
			while(trajlist >> filename){
				//Stuff for reading in xtc files
				int fp = open_xtc(filename,"r");
				int natoms,step;
				real time;
				matrix box;
				rvec *x;
				real prec;
				bool bOK;
				//Process each snapshot
				read_first_xtc(fp,&natoms,&step,&time,box,&x,&prec,&bOK);
				do{
					//DEBUG:
					if(natoms!=atoms)
						cout << "DRAMATIC ERROR ENCOUNTERED: input frame has the wrong number of atoms" << endl;
					
					//Check if this is not a repeated frame (could use a bit of fine tuning)
					if(time-prevtime>1.0){
						prevtime=time;
						real dist=ls_rmsd(atoms,x,generator,atomindicies,natomindicies);
						//Possibly update this snapshots state and distance
						if(iteration!=0&&dist<distance[t][snapshotindex]||iteration==0){
							distance[t][snapshotindex]=dist;
							state[t][snapshotindex]=iteration;
						}
						//Possibly update this processes farthest snapshot distance and its configuration
						if(distance[t][snapshotindex]>farthestdist){
							farthestdist=distance[t][snapshotindex];
							farthesttraj=t;
							farthestsnapshot=snapshotindex;
							//There might be a slightly better way to do this
							for(int atom=0;atom<atoms;atom++){
								farthestconfig[atom][0]=x[atom][0];
								farthestconfig[atom][1]=x[atom][1];
								farthestconfig[atom][2]=x[atom][2];
							}
						}
						snapshotindex++;
					}
				}while(read_next_xtc(fp,natoms,&step,&time,box,x,&prec,&bOK));
				close_xtc(fp);
				delete[] x;
			}
			trajlist.close();
			//DEBUG:
			if(snapshotindex!=lengths[t]){
				cout << "SLIGHTLY LESS DRAMATIC ERROR ENCOUNTERED: wrong number of snapshots in trajectory" << endl;
				cout << "In trajectory: " << trajnames[t] << endl;
				cout << "Expected: " << lengths[t] << endl;
				cout << "Got: " << snapshotindex << endl;
				if(snapshotindex>lengths[t])
					cout << "There were too many snapshots (how the heck did that happen????? ) so you are basically screwed." << endl;
				else{
					cout << "There were not enough snapshots so this will just be ignored, but you should really find out why this is happening since there could potentially be jumps in your data now" << endl;
					lengths[t]=snapshotindex;
				}
			}
		}

		real globalfarthest;
		int farthestpros;

		//process 0 gathers the distances and finds the largest
		MPI_Gather(&farthestdist,sizeof(real),MPI_BYTE,farthestdists,sizeof(real),MPI_BYTE,0,MPI_COMM_WORLD);
		if(process==0){
			globalfarthest=-1.0;
			for(int i=0;i<numprocess;i++){
				if(farthestdists[i]>globalfarthest){
					globalfarthest=farthestdists[i];
					farthestpros=i;
				}
			}
			cout << globalfarthest << endl;
		}
		//process 0 broadcasts the process that has the longest dist
		MPI_Bcast(&farthestpros,1,MPI_INT,0,MPI_COMM_WORLD);

		//the process that has the new generator writes it to a file and then broadcasts its configuration to everyone
		//note that this means that the last generator file is not actually a generator!!!!!
		if(process==farthestpros){
			//make the generator filename and open the file
			char filename[PATHSIZE];
			sprintf(filename,"%sgenerator.%i",GENERATORDIR,iteration+1);
			ofstream genfile(filename);
			//write out the trajectory, snapshot, generator number, and dist
			genfile << "Generator# " << iteration+1 << endl;
			genfile << "Trajectory: " << trajnames[farthesttraj] << endl;
			genfile << "Snapshot: " << farthestsnapshot << endl;
			genfile << "Distance: " << farthestdist << endl;
			//write out the configuration
			genfile << "Configuration: " << endl;
			for(int atom=0;atom<atoms;atom++)
				genfile << farthestconfig[atom][0] << " " << farthestconfig[atom][1] << " " << farthestconfig[atom][2] << endl;
			genfile.close();

			//copy over the next generator
			for(int atom=0;atom<atoms;atom++){
				generator[atom][0]=farthestconfig[atom][0];
				generator[atom][1]=farthestconfig[atom][1];
				generator[atom][2]=farthestconfig[atom][2];
			}
		}
		MPI_Bcast(generator,sizeof(rvec)*atoms,MPI_BYTE,farthestpros,MPI_COMM_WORLD);

		iteration++;
		if(iteration==10||iteration==20||iteration==50||iteration==100||iteration==200||iteration==500||iteration==1000||iteration==2000||iteration==3000){
			//write the data to the data files
			for(int t=0;t<numtrajs;t++){
				char datafilename[PATHSIZE];
				sprintf(datafilename,"%s%s.%i",DATADIR,trajnames[t],iteration);
				ofstream datafile(datafilename);
				for(int ss=0;ss<lengths[t];ss++)
					datafile << state[t][ss] << " " << distance[t][ss] << endl;
				datafile.close();
			}
		}

		if(iteration==3000)
			done=1;
	}

	/*
	//write the data to the data files
	for(int t=0;t<numtrajs;t++){
		char datafilename[PATHSIZE];
		sprintf(datafilename,"%s%s",DATADIR,trajnames[t]);
		ofstream datafile(datafilename);
		for(int ss=0;ss<lengths[t];ss++)
			datafile << state[t][ss] << " " << distance[t][ss] << endl;
		datafile.close();
	}*/
	
	for(int i=0;i<numtrajs;i++){
		delete[] distance[i];
		delete[] state[i];
		delete[] trajnames[i];
	}

	if(process==0){
		delete[] farthestdists;
	}

	delete[] trajnames;
	delete[] lengths;
	delete[] distance;
	delete[] state;
	delete[] atomindicies;

	MPI_Finalize();
}
