#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
//#include <conio.h>

// The units used in this program are reduced units.
// 1 unit of mass = 39.948 amu. (Mass of 1 argon).
// 1 unit of length = 3.4 angstroms (Sigma for argon).
// 1 unit of energy = 165.677856 (g(angs)^2)/sec. (Epsilon for argon in non SI units).
double pos[7000][3], vel[7000][3], acc[7000][3], force[7000][3], new_pos[7000][3], pre_pos[7000][3], mom[3];
float t_step = 0.02; // End time and time step
double mass_ar = 1, box = 20.0, box2 = 10.0;
float t_min = 0.0, t_max = 400.0;
FILE *kinetic, *velocities, *in_pos, *ini_pos, *ini_vel, *distance, *forces, *acceleration, *energy, *final_pos, *total_energy, *potential;
double sigma = 1.0, epsilon = 1.0; 
int n;
double sig12 = pow(sigma, 12.0);
double sig6 = pow(sigma, 6.0);

void Initial_positions(double *a1, double *a2, double *a3)
{
    int j, i;
    double r1, r2, r3, x, y ,z;
    r1 = (double)rand()/(double)RAND_MAX; // Generating random number bet. [0,1]
    r2 = (double)rand()/(double)RAND_MAX;
    r3 = (double)rand()/(double)RAND_MAX;
    x = r1*box; // Boxlength = 30.08 angs. = 8 * 1.88 * 2
    y = r2*box; // 1.88 angs. is Van der Waal's radius for argon.
    z = r3*box;
    *a1 = x;
    *a2 = y;
    *a3 = z;
}

void Distances()
{
    int i, j;
    double dxij, dyij, dzij, d1, d2, dx1, dy1, dz1, d3, d4;
    for ( i = 0; i < n ; i++ )
    {
        fprintf(distance, "For atom %d :\n", i+1);
        for ( j = 0; j < n && j != i; j++ )
        {
            dxij = &pos[i][0] - &pos[j][0];
            dyij = &pos[i][1] - &pos[j][1];
            dzij = &pos[i][2] - &pos[j][2];
            dx1 = dxij - box*floor(dxij/box);
            dy1 = dyij - box*floor(dyij/box);
            dz1 = dzij - box*floor(dzij/box);
            d1 = dxij*dxij + dyij*dyij + dzij*dzij;
            d2 = sqrt(d1);
	    d3 = dx1*dx1 + dy1*dy1 + dz1*dz1;
	    d4 = sqrt(d3);
            fprintf(distance, "\tDistance (d) with atom %d is : %f\n", j+1, d2);
	    fprintf(distance, "\tDistance (d1) with atom %d is : %f\n", j+1, d4);
        }
        fprintf(distance, "\tEnd of distance calculation for atom %d !!\n", i+1);
    }
}

void Initial_velocities(double *b1, double *b2, double *b3)
{
    int j, i;
    double r1, r2, r3, r4, vx, vy, vz, v1, v2, v3, v4, z1, z2;
    do
    {
        r1 = (double)rand()/(double)RAND_MAX; // Generating random number bet. [-1,1]
        r2 = (double)rand()/(double)RAND_MAX;
        r3 = (double)rand()/(double)RAND_MAX;
        r4 = (double)rand()/(double)RAND_MAX;
        v1 = 2.0*r1 - 1.0;
        v2 = 2.0*r2 - 1.0;
        v3 = 2.0*r3 - 1.0;
        v4 = 2.0*r4 - 1.0;
        z1 = v1*v1 + v2*v2;
        z2 = v3*v3 + v4*v4;
    } while ( z1 >= 1.0 || z2 >= 1.0);
    z1 = sqrt((-2.0*log(z1)) / z1);
    z2 = sqrt((-2.0*log(z2)) / z2);
    vx = v1*z1;
    vy = v2*z1;
    vz = v3*z2;
    *b1 = vx;
    *b2 = vy;
    *b3 = vz;
}

void Momentum()
{
    int i;
    double tot_mass = mass_ar*n;
    mom[0] = 0.0;
    mom[1] = 0.0;
    mom[2] = 0.0;
    for ( i=0 ; i<n ; i++ )
    {
        mom[0] = mom[0] + mass_ar*vel[i][0];
        mom[1] = mom[1] + mass_ar*vel[i][1];
        mom[2] = mom[2] + mass_ar*vel[i][2];
    }
    mom[0] = mom[0] / (mass_ar*n);
    mom[1] = mom[1] / (mass_ar*n);
    mom[2] = mom[2] / (mass_ar*n);
    for ( i=0 ; i<n ; i++ )
    {
        vel[i][0] = vel[i][0] - mom[0]/tot_mass;
        vel[i][1] = vel[i][1] - mom[1]/tot_mass;
        vel[i][2] = vel[i][2] - mom[2]/tot_mass;
    }
                          // previous value (needed in Verlet step)
    for ( i=0 ; i<n ; i++ )
    {
        pre_pos[i][0] = pos[i][0] - t_step*vel[i][0];
        pre_pos[i][1] = pos[i][1] - t_step*vel[i][1];
        pre_pos[i][2] = pos[i][2] - t_step*vel[i][2];
        pre_pos[i][0] = pre_pos[i][0] - box*floor(pre_pos[i][0]/box);
        pre_pos[i][1] = pre_pos[i][1] - box*floor(pre_pos[i][1]/box);
        pre_pos[i][2] = pre_pos[i][2] - box*floor(pre_pos[i][2]/box);
    }
}

void Force(float step)
{
    int i, j;
    double xij, yij, zij, fx, fy, fz, r2, r8, r14, coeff;
    for ( i = 0; i < n ; i++ )
    {
        force[i][0] = 0.0;
        force[i][1] = 0.0;
        force[i][2] = 0.0;
    }

    for ( i = 0; i < n-1 ; i++ )
    {
        for ( j = i+1; j < n ; j++ )
        {
            xij = pos[i][0] - pos[j][0]; // Cal. diff. bet. x, y & z coordinates of atoms.
            yij = pos[i][1] - pos[j][1];
            zij = pos[i][2] - pos[j][2];
            xij = xij - box*floor(xij/box);
            yij = yij - box*floor(yij/box);
            zij = zij - box*floor(zij/box);
            r2 = xij*xij + yij*yij + zij*zij; // Cal. r^8 and r^14.
            r8 = r2*r2*r2*r2;
            r14 = r8*r2*r2*r2;
            coeff = 4.0*epsilon*((12.0*sig12)/r14 - (6.0*sig6)/r8); // Cal. coefficient for force.
            fx = coeff*xij;
            fy = coeff*yij;
            fz = coeff*zij;
            force[i][0] = force[i][0] + fx;
            force[i][1] = force[i][1] + fy;
            force[i][2] = force[i][2] + fz;
            force[j][0] = force[j][0] - fx;
            force[j][1] = force[j][1] - fy;
            force[j][2] = force[j][2] - fz;
        }
    }
    fprintf(forces, "\tFor step = %6.1f\n", step*50.0);
    fprintf(forces, "\tx, y, z components of forces on each particle are: \n");
    for ( i = 0; i < n ; i++ )
    {
        fprintf(forces, "\t%f\t\t%f\t\t%f\n", force[i][0], force[i][1], force[i][2]);
    }
    fprintf(forces, "\n\tEnd of forces for step %6.1f !!\n", step*50.0);
    for ( i = 0; i < n; i++ )
    {
        acc[i][0] = force[i][0] / mass_ar;
        acc[i][1] = force[i][1] / mass_ar;
        acc[i][2] = force[i][2] / mass_ar;
    }
}

void Potentials(float step)
{
    double v2, e_kin, e_pot, e_total, factor, P[3];
    double xij, yij, zij, r2, r6, r12, r1;
    int i, j;

    e_kin = 0.0; // Cal. kinetic energy.
    for ( i = 0; i < n ; i++)
    {
        v2 = vel[i][0]*vel[i][0] + vel[i][1]*vel[i][1] + vel[i][2]*vel[i][2];
        e_kin = e_kin + v2;
    }
    e_kin = e_kin*0.5*mass_ar;
    e_pot = 0.0; // Cal. potential energy.
    for ( i = 0; i < n-1 ; i++ )
    {
        for ( j = i+1; j < n ; j++ )
        {
            xij = pos[i][0] - pos[j][0];
            yij = pos[i][1] - pos[j][1];
            zij = pos[i][2] - pos[j][2];
            xij = xij - box*floor(xij/box);
            yij = yij - box*floor(yij/box);
            zij = zij - box*floor(zij/box);
            r2 = xij*xij + yij*yij + zij*zij;
	    r1 = sqrt(r2);
            r6 = r2*r2*r2;
            r12 = r6*r6;
            factor = 4*epsilon*( sig12/r12 - sig6/r6 );
            e_pot = e_pot + factor;
	    fprintf(potential, "%f\t%f\n", r1, factor);
        }
    }
    e_total = e_kin + e_pot;

    P[0] = 0.0;
    P[1] = 0.0;
    P[2] = 0.0;
    for ( i = 0; i < n ; i++ ) // For centre of mass.
    {
        P[0] = P[0] + mass_ar*vel[i][0];
        P[1] = P[1] + mass_ar*vel[i][1];
        P[2] = P[2] + mass_ar*vel[i][2];
    }

    P[0] = P[0] / (mass_ar*n);
    P[1] = P[1] / (mass_ar*n);
    P[2] = P[2] / (mass_ar*n);

    fprintf (energy, "\n\tKinetic Energy = %f,\tPotential Energy = %f,\tTotal Energy = %f\nCentre of Mass momentum = %f\t%f\t%f\t\n", e_kin, e_pot, e_total, P[0], P[1], P[2]);
    fprintf (total_energy, "%f\t%f\n", step*50.0, e_total);
    fprintf (energy, "\tEnd of energies for step %6.1f \n", step*50.0);
    fprintf (kinetic, "%f\t%f\n", step*50.0, e_kin);
}

void Verlet(float step) // Velocity Verlet.
{
    int i;
    fprintf(forces, "\n\tFor step %6.1f :\n", step*50.0);
    fprintf(energy, "\n\tFor step %6.1f :\n", step*50.0);
    Force(step);
    Potentials(step);
    for ( i = 0; i < n ; i++ )
    {
        new_pos[i][0] = 2*pos[i][0] - pre_pos[i][0] + acc[i][0]*t_step*t_step;
        new_pos[i][1] = 2*pos[i][1] - pre_pos[i][1] + acc[i][1]*t_step*t_step;
        new_pos[i][2] = 2*pos[i][2] - pre_pos[i][2] + acc[i][2]*t_step*t_step;
        new_pos[i][0] = new_pos[i][0] - box*floor(new_pos[i][0]/box);
        new_pos[i][1] = new_pos[i][1] - box*floor(new_pos[i][1]/box);
        new_pos[i][2] = new_pos[i][2] - box*floor(new_pos[i][2]/box);
        vel[i][0] = (new_pos[i][0] - pre_pos[i][0]) / (2*t_step);
        vel[i][1] = (new_pos[i][1] - pre_pos[i][1]) / (2*t_step);
        vel[i][2] = (new_pos[i][2] - pre_pos[i][2]) / (2*t_step);
    }
    for ( i = 0; i < n ; i++ )
    {
        pre_pos[i][0] = pos[i][0];
        pre_pos[i][1] = pos[i][1];
        pre_pos[i][2] = pos[i][2];
        pos[i][0] = new_pos[i][0];
        pos[i][1] = new_pos[i][1];
        pos[i][2] = new_pos[i][2];
    }
}

int main()
{
    // Declaring mass, positions, velocities and accelerations
    int i, j;  // Number of atoms
    double t, x, y, z, xij, yij, zij, r1, r2;
    double ori_vel[7000][3], ori_pos[7000][3];
    srand(time(NULL));
    n = 64; //512 // No. of argon atoms.
    //    printf ("Simulation will run for 100,000 steps. \n");
    in_pos = fopen("Initial Positions.txt","w");
    ini_pos = fopen("Initial Positions.xyz", "w");
    ini_vel = fopen("Initial Velocities.txt","w");
    distance = fopen("Initial Distances.txt","w");
    forces = fopen("Forces.txt","w");
    energy = fopen("Energy.txt","w");
    final_pos = fopen("Final Positions.xyz","w");
    total_energy = fopen("Total Energies.txt","w");
    kinetic = fopen("Kinetic Energies.txt","w");
    potential = fopen("Potential vs. r.txt","w");
    Initial_positions(&pos[0][0], &pos[0][1], &pos[0][2]);
//    printf ("%f\t%f\t%f\n", pos[0][0], pos[0][1], pos[0][2]);
    for ( i=1; i<n; i++)
    {
        Initial_positions(&pos[i][0], &pos[i][1], &pos[i][2]);
        Initial_velocities(&vel[i][0], &vel[i][1], &vel[i][2]);
        for ( j=0; j<i; j++ )
        {
            xij = &pos[i][0] - &pos[j][0];
            yij = &pos[i][1] - &pos[j][1];
            zij = &pos[i][2] - &pos[j][2];
            r2 = xij*xij + yij*yij + zij*zij;
            r1 = sqrt(r1);
            if (r1 <= 1.106)
            {
                break;
            }

        }
    }
    fprintf(ini_pos, "%d\n\n", n);
    for ( i = 0; i < n ; i++ )
    {
        ori_pos[i][0] = pos[i][0];
        ori_pos[i][1] = pos[i][1];
        ori_pos[i][2] = pos[i][2];
        ori_vel[i][0] = vel[i][0];
        ori_vel[i][1] = vel[i][1];
        ori_vel[i][2] = vel[i][2];
        fprintf(in_pos, "Ar\t\t%f\t\t%f\t\t%f\n", ori_pos[i][0], ori_pos[i][1], ori_pos[i][2]);
        fprintf(ini_pos, "Ar\t\t%f\t\t%f\t\t%f\n", ori_pos[i][0], ori_pos[i][1], ori_pos[i][2]);
        fprintf(ini_vel, "\t%3d\t%f\t\t%f\t\t%f\n", i+1, ori_vel[i][0], ori_vel[i][1], ori_vel[i][2]);
    }
    fclose(ini_pos);
    fclose(ini_vel);
    fclose(in_pos);
    Distances();
    fclose(distance);
    Momentum();
    t = t_min;
    while (t < t_max)
    {
        t = t + t_step;
        Verlet(t);
        fprintf(final_pos, "%d\n\n", n);
        for ( i = 0; i < n ; i++ )
        {
            fprintf(final_pos, "Ar\t\t%f\t\t%f\t\t%f\n", pos[i][0], pos[i][1], pos[i][2]);
        }
	//        fprintf(final_pos, "\n\tTrajectory ended for step %6.1f !!\n", t * 10000.0);
    }
    fclose(forces);
    fclose(energy);
    fclose(final_pos);
    fclose(kinetic);
    fclose(total_energy);
    fclose(potential);
//    getch();
}
