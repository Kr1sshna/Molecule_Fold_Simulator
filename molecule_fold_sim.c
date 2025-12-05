#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

typedef struct atom{ //defining a struct atom which stores the atom name, its type(hydrophobic or polar), coordinates, velocity vector components, force vector components and mass of the atom
	char type;
	char name[5];
	double x;
	double y;
	double z;
	double vx;
	double vy;
	double vz;
	double fx;
	double fy;
	double fz;
	double mass;
} ATOM;

typedef struct bond{  //defining a struct for linear bonds, stores the index of atom 1 and atom 2 which are bonded, the spring stiffness(Hooke's constant) of the bond and the optimal length of the bond
    int a1; 
    int a2; 
    double k;  
    double len; 
} BOND;


double calc_rg(ATOM a[], int size) {   //to calculate radius of gyration of the molecule
    double com_x = 0, com_y = 0, com_z = 0;   //initialises the centre of mass coordinates at origin
    
    for(int i = 0; i < size; i++) {   
        com_x += a[i].x;
        com_y += a[i].y;
        com_z += a[i].z;
    }
    com_x /= size;
    com_y /= size;
    com_z /= size;  //finds the centre of mass coordinates for a step

    double sum_sq_dist = 0;
    for(int i = 0; i < size; i++) {
        double dx = a[i].x - com_x;
        double dy = a[i].y - com_y;
        double dz = a[i].z - com_z;
        sum_sq_dist += (dx*dx + dy*dy + dz*dz);   //calculates the root mean square difference of all atom position to that of the centre of mass of the molecule, giving us the radius of gyration
    }

    return sqrt(sum_sq_dist / size);
}

double gaussian_rand() {    //function to produce random thermal noises which occur in real
    double u1 = ((double) rand() / RAND_MAX);   
    double u2 = ((double) rand() / RAND_MAX);   
    
    if (u1 < 1e-10) u1 = 1e-10;
    
    double z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * 3.1415926535 * u2);   //uses the box-muller Transform equation to give us a gaussian rather than a random number with equal probability
    return z0;
}

double force(double dist,double epsilon){ //function to define the Lennard-Jones Force, uses normalised value
	if (dist < 0.000001){   //checks for division by zero error
		dist = 0.000001;
	}
	double force = epsilon * 24 * (2*pow(1/dist,13) - pow(1/dist,7));  //epsilon determine the force increase and decrease due to the type of atom1 and atom2(neighbouring atoms)
	return force;
}

double dist_x(ATOM a, ATOM b){  //calculates the x axis distance
	double x_dist = a.x - b.x;
	return x_dist;
}

double dist_y(ATOM a, ATOM b){   //calculates the y axis distance
	double y_dist = a.y - b.y;
	return y_dist;
}

double dist_z(ATOM a, ATOM b){   //calculates the z axis distance
	double z_dist = a.z - b.z;
	return z_dist;
}

double distance(ATOM a, ATOM b){   //calculates the euclidean distance
	double x_dist = a.x - b.x;
	double y_dist = a.y-b.y;
	double z_dist = a.z-b.z;
	double distance = sqrt(pow(x_dist,2) + pow(y_dist,2) + pow(z_dist,2));
	return distance;
}

double new_vel(double old_vel,double force,double time,double mass){   //calculates the new velocity from newton 2nd law of motion (F = m*a)
	double new_vel = old_vel + ((force*time)/mass);
	return new_vel;
}

double new_pos(double old_pos, double new_vel, double time){   //calculates new position from the instantaneus velocity formula
	double new_pos = old_pos + (new_vel * time);
	return new_pos;
}

int simulation(ATOM a[], BOND b[]){     //the function which carries all th simulation task
	
	int size = 20;     //20 atoms are used 
	int num_bond = size-1;   //linear bonds will always be equal to no. of atoms - 1
	int steps = 100000;   //the amount of steps the simulation will run for, which is about 1000 seconds(steps * dt)
	double time = 0;    //time initialised at zero
	double dt = 0.01;   //time increment of 0.01 seconds
	double dist;   
	double dx;
	double dy;
	double dz;
	double force_lj = 0;
	double temperature = 5.0;  //temprture of 5 unit
	double kb = 1.0;   //normalised boltsman constant
	double gamma = 2;   //gamma to stimulate the friction in water
	double epsilon = 0;  //initialises epsilon to 0
	
	FILE *fp = fopen("trajectory.xyz", "w");   //opens a file pointer for storing trajectory data
	if (fp == NULL) {
    		printf("Error opening file!\n");
    		return 1;
	}
	
	FILE *fp_rg = fopen("analysis.csv", "w");   //opens a file pointer for storing the radius of gyration data
	if (fp_rg == NULL) {
    		printf("Error opening file!\n");
    		return 1;
	}
	
	for(int i = 0; i<num_bond; i++){  //assigns bonds to all neighbouring atoms
			b[i].a1 = i;
			b[i].a2 = i+1;
			b[i].k = 50.00;
			b[i].len = 1.5;
		}
	
	for(int step = 0; step < steps; step++){  
		time = step * dt;    //time of the step
		
		for(int i = 0; i < size; i++){  //zeroes the force vector after each iteration
			a[i].fx = 0;
			a[i].fy = 0;
			a[i].fz = 0;
		}
		
		for(int j = 0; j < size; j++){  //runs a loop of all combination possible for 2 atom to interct that is, the time complexity is O(size!/(size-2)!*2!)
			for(int k = j+1; k < size; k++){
				
				dist = distance(a[j],a[k]);
				dx = dist_x(a[j],a[k]);
				dy = dist_y(a[j],a[k]);
				dz = dist_z(a[j],a[k]);
				
				if (a[j].type == 'H' && a[k].type == 'H'){  //assigns the epsilon based on the types of the 2 atom, for ex:- hydrophobic atoms(denoted by 'H') will try to stay closer together
					epsilon = 10;
				} else if (a[j].type == 'P' && a[k].type == 'P'){
					epsilon = 1;
				} else{
					epsilon = 0.5;
				}
				force_lj = force(dist,epsilon);  //finds the Lennard-Jones Force for the 2 atom
		
				double fx = force_lj * (dx/dist);  //finds the vector components of the force
				double fy = force_lj * (dy/dist);
				double fz = force_lj * (dz/dist);
				
				a[j].fx += fx;  //by newton third law of motion, every action has a equal and oposite reaction and thus, one atom will have a positive force vector and the other will have a negative force vector
				a[k].fx -= fx;
				a[j].fy += fy;
				a[k].fy -= fy;
				a[j].fz += fz;
				a[k].fz -= fz;
			}
		}
		
		for(int i = 0; i < num_bond; i++){   //runs the loop the size of num_bonds
		
			int ind1 = b[i].a1;   //finds the atom1 and atom2 indexes for bond[i]
			int ind2 = b[i].a2;
			
			dist = distance(a[ind1],a[ind2]);
			dx = dist_x(a[ind1],a[ind2]);
			dy = dist_y(a[ind1],a[ind2]);
			dz = dist_z(a[ind1],a[ind2]);
			
			double f_bond = -(b[i].k)*(dist-b[i].len);  //calculates the bonding force between 2 neighbouring atoms
			
			double fx = f_bond * (dx/dist);  //finds the vector components of the force
			double fy = f_bond * (dy/dist);
			double fz = f_bond * (dz/dist);
				
			a[ind1].fx += fx;  // A positive LJ force denotes repulsion and a negative bond force denotes a stretched condition(i.e. attraction)
			a[ind2].fx -= fx;
			a[ind1].fy += fy;
			a[ind2].fy -= fy;
			a[ind1].fz += fz;
			a[ind2].fz -= fz;
		}
		
		for (int i = 0; i < size; i++){ //subtracts the friction force for a damped oscillator and adds the random thermal noise for the wiggly, conformation changing nature of proteins(this molecule is just a replica of protein, it is not a protein, it just simulates the same folding of protein due to hydrophobic and hydrophilic interaction)
		
			a[i].fx -= (gamma * a[i].vx);
			a[i].fy -= (gamma * a[i].vy);
			a[i].fz -= (gamma * a[i].vz);
		
			double sigma = sqrt((2.0 * a[i].mass * gamma * kb * temperature) / dt);
			
			a[i].fx += sigma * gaussian_rand();
    			a[i].fy += sigma * gaussian_rand();
    			a[i].fz += sigma * gaussian_rand();
		
			a[i].vx = new_vel(a[i].vx,a[i].fx,dt,a[i].mass);  //finally after the total force is calculated( F_total = F_lj + F_bond - F_friction + F_thermal), we calculate the velocity at that step
			a[i].vy = new_vel(a[i].vy,a[i].fy,dt,a[i].mass);
			a[i].vz = new_vel(a[i].vz,a[i].fz,dt,a[i].mass);
		
			a[i].x = new_pos(a[i].x,a[i].vx,dt);  //and then calculates the new position at that step
			a[i].y = new_pos(a[i].y,a[i].vy,dt);
			a[i].z = new_pos(a[i].z,a[i].vz,dt);
		}
		
		if (step % 100 == 0) {   //writes the radius of gyration values in analysis.csv at every 100 steps
        		double current_rg = calc_rg(a, size);
        
        		fprintf(fp_rg, "%lf,%lf\n", time, current_rg);
    		}
		
		if (step % 10 == 0) {    //writes the position value in trajectory.xyz at every 10 steps, so that we can use VMD to virtually observe the folding and fluctuations
    			fprintf(fp, "%d\n", size); 
    
    			fprintf(fp, "Step %d, Time %.3f\n", step, time);
    
    			for (int i = 0; i < size; i++) {
        			fprintf(fp, "%s %lf %lf %lf\n",a[i].name, a[i].x, a[i].y, a[i].z);
    			}
		}
	}
	fclose(fp); //closes the file pointer for both the files
	fclose(fp_rg);
	printf("Simulation Complete. Data saved to trajectory.xyz\n");  //prints to the terminal to let the user know
	printf("Analysis complete. Data saved to analysis.csv\n");
}
	

int main(){
	srand(time(NULL));  //to get a new random value every time the program runs
	ATOM a[20]= {  //an array of 20 atoms
	{'H',"C",0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 12.0},
	{'H',"C",2.0, 0.1, -0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 12.0},
	{'H',"C",4.0, -0.1, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 12.0},
	{'H',"C",6.0, 0.1, -0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 12.0},
	{'H',"C",8.0, -0.1, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 12.0},
	{'H',"C",10.0, 0.1, -0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 12.0},
	{'H',"C",12.0, -0.1, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 12.0},
	{'P',"O",14.0, 0.1, -0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 12.0},
	{'P',"O",16.0, -0.1, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 12.0},
	{'P',"O",18.0, 0.1, -0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 12.0},
	{'P',"O",20.0, -0.1, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 12.0},
	{'P',"O",22.0, 0.1, -0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 12.0},
	{'P',"O",24.0, -0.1, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 12.0},
	{'P',"O",26.0, 0.1, -0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 12.0},
	{'P',"O",28.0, -0.1, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 12.0},
	{'H',"C",30.0, 0.1, -0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 12.0},
	{'H',"C",32.0, -0.1, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 12.0},
	{'H',"C",34.0, 0.1, -0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 12.0},
	{'H',"C",36.0, -0.1, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 12.0},
	{'H',"C",38.0, 0.1, -0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 12.0}
	};
	
	BOND b[19]; //an array of 19 bonds
	
	int result = simulation(a,b);  //runs the simulation
	if(result == 1){  //if error in opening file
		return 1;  //ends the program with return 1
	}
	
	return 0;
}
