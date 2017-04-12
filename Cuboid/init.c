#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_randist.h>

#include "init.h"
#include "data.h"
#include "file_io.h"
#include "calculate.h"

void read_parameters(const char *fname, params *pptr, connect *cptr)
{
	FILE *fptr;
	char str[200], *tokptr;
	
	fptr = file_open(fname, "r");
	
	while(!feof(fptr))
	{
		fgets(str, 200, fptr);
		
		tokptr = strtok(str, "=");
		
		while(tokptr!=NULL)
		{
			trim_space(tokptr);
			
			if(strcmp(tokptr, "time")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->time = atof(tokptr);
			}
            if(strcmp(tokptr, "input time")==0)
            {
                tokptr = strtok(NULL, "=");
                pptr->input_time = atof(tokptr);
            }
            if(strcmp(tokptr, "path time")==0)
            {
                tokptr = strtok(NULL, "=");
                pptr->path_time = atof(tokptr);
            }
            else if(strcmp(tokptr, "total HD cells")==0)
            {
                tokptr = strtok(NULL, "=");
                pptr->num_HD_cells = atoi(tokptr);
            }
			else if(strcmp(tokptr, "RC connections")==0)
			{
				tokptr = strtok(NULL, "=");
				cptr->num_RC_connections = atoi(tokptr);
			}
			else if(strcmp(tokptr, "phi RC")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->phi_RC = atof(tokptr);
			}
            else if(strcmp(tokptr, "sigma RC")==0)
            {
                tokptr = strtok(NULL, "=");
                pptr->sigma_RC = atof(tokptr);
            }
            else if(strcmp(tokptr, "tau HD")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->tau_HD = atof(tokptr);
			}
			else if(strcmp(tokptr, "alpha HD")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->alpha_HD = atof(tokptr);
			}
			else if(strcmp(tokptr, "beta HD")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->beta_HD = atof(tokptr);
			}
			else if(strcmp(tokptr, "timestep size")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->timestep_size = atof(tokptr);
			}
            else if(strcmp(tokptr, "external inhibition")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->extern_inhib = atof(tokptr);
			}
			else if(strcmp(tokptr, "global inhibition")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->global_inhibition = atof(tokptr);
			}
			else if(strcmp(tokptr, "input strength")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->input_str = atof(tokptr);
			}
			else if(strcmp(tokptr, "input location")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->input_loc = atof(tokptr);
			}
            else if(strcmp(tokptr, "input sigma")==0)
            {
                tokptr = strtok(NULL, "=");
                pptr->input_sigma = atof(tokptr);
            }
            else if(strcmp(tokptr, "combined yaw input strength")==0)
            {
                tokptr = strtok(NULL, "=");
                pptr->combined_yaw_input_str = atof(tokptr);
            }
            else if(strcmp(tokptr, "combined yaw input sigma")==0)
            {
                tokptr = strtok(NULL, "=");
                pptr->combined_yaw_input_sigma = atof(tokptr);
            }
            else if(strcmp(tokptr, "normalise")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->normalise = atoi(tokptr);
			}
            else if(strcmp(tokptr, "conduction delay")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->conduction_delay = atof(tokptr);
			}			
			else if(strcmp(tokptr, "parallel threads")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->threads = atoi(tokptr);
			}
            else if(strcmp(tokptr, "sigmoid")==0)
            {
                tokptr = strtok(NULL, "=");
                pptr->sigmoid = atoi(tokptr);
            }
            else if(strcmp(tokptr, "gravity increment switch")==0)
            {
                tokptr = strtok(NULL, "=");
                pptr->gravity_increment_switch = atoi(tokptr);
            }

            else
				tokptr = strtok(NULL, "=");
			
		}
	}
	
	fclose(fptr);
	
	return;
}


void build_network(data *dptr, params *pptr, connect *cptr, cArray **bptr)
{
	int idx;
	int connection;
	
	pptr->timesteps = (int)(pptr->time/pptr->timestep_size);
    pptr->input_timesteps = (int)(pptr->input_time/pptr->timestep_size);
    pptr->path_timesteps = (int)(pptr->path_time/pptr->timestep_size);
    
    //Here's the conduction buffer allocation stuff
	
	pptr->conduction_buffer_size = (int)(pptr->conduction_delay/pptr->timestep_size);
	
	for (connection=0; connection<pptr->num_HD_cells; connection++)						
	{
		bptr[connection]->array = (float *)malloc(pptr->conduction_buffer_size *sizeof(float));
	}
    
    //Here's HD rates and activations
    
	dptr->activations_HD = (float *)malloc(pptr->num_HD_cells*sizeof(float));
	dptr->prev_activations_HD = (float *)malloc(pptr->num_HD_cells*sizeof(float));
   
	dptr->rates_HD = (float *)malloc(pptr->num_HD_cells*sizeof(float));
	dptr->prev_rates_HD = (float *)malloc(pptr->num_HD_cells*sizeof(float));
    
	dptr->rates_HD_time = (float **)malloc(pptr->num_HD_cells*sizeof(float *));
	dptr->rates_HD_time[0] = (float *)malloc((pptr->num_HD_cells*(pptr->timesteps))*sizeof(float));
	
	for(idx=1; idx<pptr->num_HD_cells; idx++)
	{
		dptr->rates_HD_time[idx] = dptr->rates_HD_time[0] + (idx * (pptr->timesteps));
	}
  
	dptr->activations_HD_time = (float **)malloc(pptr->num_HD_cells*sizeof(float *));
	dptr->activations_HD_time[0] = (float *)malloc((pptr->num_HD_cells*(pptr->timesteps))*sizeof(float));
	
	for(idx=1; idx<pptr->num_HD_cells; idx++)
	{
		dptr->activations_HD_time[idx] = dptr->activations_HD_time[0] + (idx * (pptr->timesteps));
	}
	
    //Here's RC connectivity
    
	cptr->RC_connections = (int **)malloc(pptr->num_HD_cells*sizeof(int *));
	cptr->RC_connections[0] = (int *)malloc((pptr->num_HD_cells*cptr->num_RC_connections)*sizeof(int));
	
	for(idx=1; idx<pptr->num_HD_cells; idx++)
	{
		cptr->RC_connections[idx] = cptr->RC_connections[0] + (idx * cptr->num_RC_connections);
	}
		
	cptr->RC_weights =(float **)malloc(pptr->num_HD_cells*sizeof(float *));
	cptr->RC_weights[0] = (float *)malloc((pptr->num_HD_cells*cptr->num_RC_connections)*sizeof(float));
	
	for(idx=1; idx<pptr->num_HD_cells; idx++)
	{
		cptr->RC_weights[idx] = cptr->RC_weights[0] + (idx * cptr->num_RC_connections);
	}
    
	//Here's a matrix for the full possible RC connectivity, not the sample from sparse connectivity (if needed)
	
	cptr->full_RC =(float **)malloc(pptr->num_HD_cells*sizeof(float *));
	cptr->full_RC[0] = (float *)malloc((pptr->num_HD_cells*pptr->num_HD_cells)*sizeof(float));
	
	for(idx=1; idx<pptr->num_HD_cells; idx++)
	{
		cptr->full_RC[idx] = cptr->full_RC[0] + (idx * pptr->num_HD_cells);
	}

	cptr->prev_weights_RC =(float **)malloc(pptr->num_HD_cells*sizeof(float *));
	cptr->prev_weights_RC[0] = (float *)malloc((pptr->num_HD_cells*cptr->num_RC_connections)*sizeof(float));
	
	for(idx=1; idx<pptr->num_HD_cells; idx++)
	{
		cptr->prev_weights_RC[idx] = cptr->prev_weights_RC[0] + (idx * cptr->num_RC_connections);
	}	
 
    //Here's misc extra stuff - inputs and favoured views
    
    dptr->favoured_view = (float *)malloc(pptr->num_HD_cells*sizeof(float));

    dptr->input = (float *)malloc(pptr->num_HD_cells*sizeof(float));
	
	dptr->input_location_time = (float **)malloc(pptr->num_HD_cells*sizeof(float *));
	dptr->input_location_time[0] = (float *)malloc((pptr->num_HD_cells*(pptr->timesteps))*sizeof(float));
	
	for(idx=1; idx<pptr->num_HD_cells; idx++)
	{
		dptr->input_location_time[idx] = dptr->input_location_time[0] + (idx * (pptr->timesteps));
	}

    dptr->combined_yaw_input = (float *)malloc(pptr->num_HD_cells*sizeof(float));
    
    dptr->combined_yaw_input_location_time = (float **)malloc(pptr->num_HD_cells*sizeof(float *));
    dptr->combined_yaw_input_location_time[0] = (float *)malloc((pptr->num_HD_cells*(pptr->timesteps))*sizeof(float));
    
    for(idx=1; idx<pptr->num_HD_cells; idx++)
    {
        dptr->combined_yaw_input_location_time[idx] = dptr->combined_yaw_input_location_time[0] + (idx * (pptr->timesteps));
    }

    
    //Here's path related stuff read in from MATLAB - it's timesteps+1 as that's how MATLAB outputs
    
    dptr->positions_x = (float *)malloc((pptr->path_timesteps+1)*sizeof(float));
    dptr->positions_y = (float *)malloc((pptr->path_timesteps+1)*sizeof(float));
    dptr->positions_z = (float *)malloc((pptr->path_timesteps+1)*sizeof(float));
    
    dptr->headings_x = (float *)malloc((pptr->path_timesteps+1)*sizeof(float));
    dptr->headings_y = (float *)malloc((pptr->path_timesteps+1)*sizeof(float));
    dptr->headings_z = (float *)malloc((pptr->path_timesteps+1)*sizeof(float));

    dptr->pvector_time = (float *)malloc(pptr->path_timesteps*sizeof(float));

    dptr->combined_input_pvector_time = (float *)malloc(pptr->path_timesteps*sizeof(float));
    
    dptr->ideal_increment = (float *)malloc(pptr->path_timesteps*sizeof(float));
    dptr->actual_increment = (float *)malloc(pptr->path_timesteps*sizeof(float));
    
    
    //Saving the reference vector
    dptr->refVec_time = (double **)malloc(3*sizeof(double *));
    dptr->refVec_time[0] = (double *)malloc((3*(pptr->path_timesteps))*sizeof(double));
    
    for(idx=1; idx<3; idx++)
    {
        dptr->refVec_time[idx] = dptr->refVec_time[0] + (idx * (pptr->path_timesteps));
    }

    //Saving rotated reference vector
    
    dptr->rotated_refVec_time = (double **)malloc(3*sizeof(double *));
    dptr->rotated_refVec_time[0] = (double *)malloc((3*(pptr->path_timesteps))*sizeof(double));
    
    for(idx=1; idx<3; idx++)
    {
        dptr->rotated_refVec_time[idx] = dptr->rotated_refVec_time[0] + (idx * (pptr->path_timesteps));
    }
    
    //Saving rotated heading vector
    
    dptr->rotated_heading_time = (double **)malloc(3*sizeof(double *));
    dptr->rotated_heading_time[0] = (double *)malloc((3*(pptr->path_timesteps))*sizeof(double));
    
    for(idx=1; idx<3; idx++)
    {
        dptr->rotated_heading_time[idx] = dptr->rotated_heading_time[0] + (idx * (pptr->path_timesteps));
    }

    //Saving start angle
    dptr->start_angle = (double *)malloc(pptr->path_timesteps*sizeof(double));

    
	return;
}


void set_connectivity(connect *cptr, params *pptr)
{
    
    /* This function randomly sets connectivity (i.e. which cells are connected to which). This actually doesn't do anything when the network is fully interconnected....*/
    
	int cell;
	const gsl_rng_type *t=NULL;
    gsl_rng *r=NULL;
	
	gsl_rng_env_setup();
    
    t=gsl_rng_default;
    r=gsl_rng_alloc(t);
    
    gsl_rng_set(r, 1);
	
	int a[pptr->num_HD_cells];
	
	for(cell=0; cell<pptr->num_HD_cells; cell++)
	{
		a[cell] = cell;
	}
	
	
	
	for(cell=0; cell<pptr->num_HD_cells; cell++)
	{
		gsl_ran_choose(r, &(cptr->RC_connections[cell][0]), cptr->num_RC_connections, a, pptr->num_HD_cells, sizeof(int)); 

	}
	
	gsl_rng_free(r);
	
	return;
}


void set_favoured_view(data *dptr, params *pptr)
{
    /*  This function generates arrays of favoured view in the yaw plane for each cell
     Cells are simply indexed by num_HD_cells*/
     
     int cell;
     float increment;
     
     increment = 360.0/(float)pptr->num_HD_cells;
     
     for(cell=0; cell<pptr->num_HD_cells; cell++)
     {
     dptr->favoured_view[cell] = (float)cell * increment;
     }
     

    
    save_prefdirs(dptr, pptr);
    
    return;
}


void set_RC_weights(params *pptr, data *dptr, connect *cptr)
{
    
    /* This function sets recurrent collateral weights, from the perspective of a 
     postsynaptic HD cell (i.e. setting all the weights that a given cell receives) */
	
    int cell, connection, presynaptic;
    float distance1, distance2, distance;

	for(cell=0; cell<pptr->num_HD_cells; cell++)
	{
		
		for(connection=0; connection<cptr->num_RC_connections; connection++)
		{
			presynaptic = cptr->RC_connections[cell][connection];
			
            distance1 = fabs(dptr->favoured_view[presynaptic] - dptr->favoured_view[cell]);
            distance2 = fabs(360.0 - distance1);
			
			if(distance1<=distance2)
				distance = distance1;
			else
			{
				distance = distance2;
			}
			
            cptr->RC_weights[cell][connection] = exp(-0.5 * (distance/pptr->sigma_RC) * (distance/pptr->sigma_RC));
			
		}
	}
	
		
	normalise_RC_weights(cptr, pptr);
	
	return;
}




void free_memory(data *dptr, connect *cptr, cArray **bptr, params *pptr)
{
	int connection;
	
	free(dptr->activations_HD);
	free(dptr->prev_activations_HD);

	free(dptr->rates_HD);
	free(dptr->prev_rates_HD);
			
	free(dptr->rates_HD_time[0]);
	free(dptr->rates_HD_time);
	
	free(dptr->activations_HD_time[0]);
	free(dptr->activations_HD_time);
	
    free(dptr->favoured_view);
    free(dptr->input);
    
    free(dptr->combined_yaw_input);
		
	free(cptr->RC_connections[0]);
	free(cptr->RC_connections);
	
	free(cptr->RC_weights[0]);
	free(cptr->RC_weights);
		
	free(cptr->full_RC[0]);
	free(cptr->full_RC);
		
	free(cptr->prev_weights_RC[0]);
	free(cptr->prev_weights_RC);
		
	free(dptr->input_location_time[0]);
	free(dptr->input_location_time);
    
    free(dptr->combined_yaw_input_location_time[0]);
    free(dptr->combined_yaw_input_location_time);
    
	for(connection=0; connection<pptr->num_HD_cells; connection++)
	{
		free(bptr[connection]->array);
	}
    
    free(dptr->positions_x);
    free(dptr->positions_y);
    free(dptr->positions_z);
    
    free(dptr->headings_x);
    free(dptr->headings_y);
    free(dptr->headings_z);
    
    free(dptr->pvector_time);
    free(dptr->start_angle);
    free(dptr->combined_input_pvector_time);
    free(dptr->ideal_increment);
    free(dptr->actual_increment);
    
    free(dptr->refVec_time[0]);
    free(dptr->refVec_time);
    
    free(dptr->rotated_refVec_time[0]);
    free(dptr->rotated_refVec_time);
    
    free(dptr->rotated_heading_time[0]);
    free(dptr->rotated_heading_time);
    
  
        
    
	return;
}
