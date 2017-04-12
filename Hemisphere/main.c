#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <libiomp/omp.h>

#include "calculate.h"
#include "data.h"
#include "file_io.h"
#include "init.h"

/*THIS CODE SIMULATES JONATHAN'S EXPERIMENTS
 
 Rat walks straight up south wall, turns 90 to the left, then walks over to the west wall. The rat then turns 180 degrees
 and returns to the south wall. It continues over to the east wall.
 
 */

int main(int argc, char **argv)
{
	char *fname;
	int timestep;   //keeps track of timestep
	params *pptr;   //initialises parameters
	data *dptr;     //intialises data
	connect *cptr;  //initialises connections
	cArray **bptr;  //initialises the conduction array
	int counter;
	int connection;
    float increment;
	double start, end;  //keeps track of runtime for openMP speed testing
	
	if(argc!=3)
	{
		fprintf(stderr, "Usage: gravity_vector -fname <parameters_file>\n");  //warning in case of missing parameter inputs
		exit(1);
	}
	else
	{
		pptr = (params *)malloc(sizeof(params));    //allocating memory for pptr
		cptr = (connect *)malloc(sizeof(connect));  //allocating memory for cptr
		fname = argv[2];
	}

    read_parameters(fname, pptr, cptr); //reading in parameters from 3Dparams.dat
	
	dptr = (data *)malloc(sizeof(data));    //allocating memory for dptr
   
	bptr = (cArray **)malloc(sizeof(cArray *) * pptr->num_HD_cells);    //allocating memory for bptr (indexing conduction array)
    
	for(counter=0; counter<pptr->num_HD_cells; counter++)
	{
		bptr[counter] = malloc(sizeof(cArray)); //allocating individual bptr arrays
	}
		
	build_network(dptr, pptr, cptr, bptr);  //allocating memory for everything to do with network
	
   	set_connectivity(cptr, pptr);   //setting up connectivity
    set_favoured_view(dptr, pptr);  //setting up favoured view
	
	set_RC_weights(pptr, dptr, cptr); //setting up the recurrent weight profile
	
		
	//load_full_RC(cptr, pptr); ONLY NEEDED IF CONNECTIVITY IS INCOMPLETE
	
	save_RCweights(cptr, pptr); //recording the initialised weight profile
    	
	/* Simulate Network */
	
	zero_states(dptr, pptr, bptr, cptr);    //this and zero_input making sure everything's initiliased to zero
	zero_input(dptr, pptr);
		
	printf("\nNetwork now being simulated...\n");
	fflush(stdout);
	
    dptr->combined_yaw_input_loc = pptr->input_loc;   //making sure egocentric/gravity input starts at packet location

   	
	printf("Applying input...\n");
	fflush(stdout);

	/* START OF PARALLEL SECTION */
	omp_set_num_threads(pptr->threads); //sets number of threads
	start = omp_get_wtime();            //records start time for speed testing
	
#pragma omp parallel private (timestep) //Beginning of the parallel region
{
	
	for(timestep=0; timestep<pptr->input_timesteps; timestep++)
	{
#pragma omp single nowait
	{
		printf("\rTimestep %d", timestep+1);    //This code just displays current timestep in terminal window
		fflush(stdout);
	}

		record_previous_states(dptr, pptr, cptr);   //Recording network state
		
		set_input(dptr, pptr);                      //Applying input packet
        
		calculate_activations_excitatory(dptr, pptr, cptr, bptr, timestep); //This function calculates HD activation
		calculate_rates_excitatory(dptr, pptr);     //Transfer function from acts to rates
		

	
#pragma omp for private (connection)
		for(connection=0; connection<pptr->num_HD_cells; connection++) 
		{
			fill_RC_buffer(dptr, pptr, bptr, connection);   //Loading conduction delay buffer with presynaptic weights
		}

				
#pragma omp single
		{
				record_activations_rates(dptr, pptr, (timestep));   //Recording HD activations and rates
				record_input(dptr,pptr,(timestep)); //Recording inputs to network
		
			
		}			
	}
	
	
#pragma omp single nowait
	{
		printf("\nInput removed...");
		fflush(stdout);
	}
	
		zero_input(dptr, pptr);
    
	
#pragma omp single nowait
    {
        printf("\nRat exploring hemisphere...\n");
        fflush(stdout);
    }
    
    //NOW READING IN THE HEADING AND POSITION DATA FOR THE RAT'S PATH
    
#pragma omp single
    {
        read_in_path(pptr,dptr);
    }
    
    
	for(timestep=pptr->input_timesteps; timestep<pptr->input_timesteps + pptr->path_timesteps; timestep++)
	{
#pragma omp single nowait
	{
		printf("\rTimestep %d", timestep+1);
		fflush(stdout);
	}

#pragma omp single
        {
            
            calculate_pvector(dptr, pptr);
            dptr->pvector_time[timestep - pptr->input_timesteps] = dptr->pvector;
            
            set_egocentric_yaw_increment(dptr, pptr, timestep);
            //ALT_egocentric_yaw_increment(dptr, pptr, timestep);
            //set_gravity_yaw_increment(dptr, pptr, timestep);
        }
        
#pragma omp single
            {
                
                
                increment = dptr->gravity_increment + dptr->egocentric_increment;
                
                printf("\n\nG:%f E:%f T:%f\n\n",dptr->gravity_increment,dptr->egocentric_increment,increment);
                fflush(stdout);
                
                dptr->actual_increment[timestep-pptr->input_timesteps] = increment;
                
                dptr->combined_yaw_input_loc += increment;
                
                if(dptr->combined_yaw_input_loc >= 360.0) //if loc goes above 360
                {
                    dptr->combined_yaw_input_loc -= 360.0;
                }
                if(dptr->combined_yaw_input_loc < 0.0) //if loc goes below 0
                {
                    dptr->combined_yaw_input_loc += 360.0;
                }
                
                dptr->combined_input_pvector_time[timestep - pptr->input_timesteps] = dptr->combined_yaw_input_loc;
            }
            
        set_combined_yaw_input(dptr, pptr);
            
        		
		record_previous_states(dptr, pptr, cptr);
		calculate_activations_excitatory(dptr, pptr, cptr, bptr, timestep);
		calculate_rates_excitatory(dptr, pptr);
				
#pragma omp for private (connection)
		for(connection=0; connection<pptr->num_HD_cells; connection++)		
		{
			fill_RC_buffer(dptr, pptr, bptr, connection);
		}
		
	
#pragma omp single
		{
		
				record_activations_rates(dptr, pptr, (timestep));
				record_input(dptr,pptr,(timestep));
            
            
		}		
		
	}
    

}
	
	end = omp_get_wtime();
	
	printf("\n\nRuntime = %f seconds\n", end-start);
	
	
	printf("\nSaving results...");
	fflush(stdout);

	
	
	save_input(dptr,pptr);                  //This saves the input over time
	save_rates_activations(dptr, pptr);     //This saves rates and activations over time
    save_pvector(dptr, pptr);               //This save the pvector over time, as well input pvector over time
	free_memory(dptr, cptr, bptr, pptr);
	
	free(dptr);
	free(cptr);
	free(pptr);
	
	for(connection=0;connection<pptr->num_HD_cells; connection++)
	{
		free(bptr[connection]);
	}
	
	free(bptr);
	
	
	printf("\nProgram successfully finished\n\n");
	fflush(stdout);
	
	return 0;
	
}



