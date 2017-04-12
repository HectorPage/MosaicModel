#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "calculate.h"
#include "data.h"

void zero_states(data *dptr, params *pptr, cArray **bptr, connect *cptr)
{
    
    /* This function makes sure everything is initially set to zero*/
    
	int cell;
	int connection;
	
	for(cell=0; cell<pptr->num_HD_cells; cell++)
	{
		dptr->activations_HD[cell] = 0.0;
		dptr->prev_activations_HD[cell] = 0.0;
		dptr->rates_HD[cell] = 0.0;
		dptr->prev_rates_HD[cell] = 0.0;
        dptr->combined_yaw_input[cell] = 0.0;
    }
	
	for(cell=0; cell<pptr->conduction_buffer_size; cell++)
	{
		for(connection=0; connection<pptr->num_HD_cells; connection++) 
		{
			bptr[connection]->array[cell] = 0.0;			
		}
	}    
	
		return;
}


void set_input(data *dptr, params *pptr)
{
    /* This function generates an input packet in the HD layer, specified by an input_str and input_loc parameters*/
    
    int cell;
    float distance1, distance2, distance;


        
#pragma omp for private(cell, distance, distance1, distance2)
        for(cell=0; cell<pptr->num_HD_cells; cell++)
        {
            distance1 = fabs(pptr->input_loc - dptr->favoured_view[cell]);
            distance2 = fabs(360.0 - distance1);
            
            if(distance1<=distance2)
                distance = distance1;
            else
                distance = distance2;
            
            dptr->input[cell] = pptr->input_str * exp(-0.5 * (distance/pptr->input_sigma) * (distance/pptr->input_sigma));
        }
    
}

void set_combined_yaw_input(data *dptr, params *pptr)
{
    /* This function generates an input packet in the HD layer, resulting from egocentric yaw rotations*/
    
    int cell;
    float distance1, distance2, distance;
    

#pragma omp for private(cell, distance, distance1, distance2)
    for(cell=0; cell<pptr->num_HD_cells; cell++)
    {
        distance1 = fabs(dptr->combined_yaw_input_loc - dptr->favoured_view[cell]);
        distance2 = fabs(360.0 - distance1);
        
        if(distance1<=distance2)
            distance = distance1;
        else
            distance = distance2;
        
        dptr->combined_yaw_input[cell] = pptr->combined_yaw_input_str * exp(-0.5 * (distance/pptr->combined_yaw_input_sigma) * (distance/pptr->combined_yaw_input_sigma));
      
    }

}


void zero_input(data *dptr, params *pptr)
{
    /* This function ensures that input is zero unless it is explicity initialised */
	int cell;
#pragma omp for private(cell)
	for(cell=0; cell<pptr->num_HD_cells; cell++)
	{
		dptr->input[cell] = 0.0;
	}
	
	return;
}



void record_previous_states(data *dptr, params *pptr, connect *cptr)
{
    /* This function stores HD activations, rates, and weights - to be used at the next timestep */
    
#pragma omp sections
	{
#pragma omp section
		{
			memcpy(dptr->prev_activations_HD, dptr->activations_HD, pptr->num_HD_cells*sizeof(float));
		}
#pragma omp section
		{
			memcpy(dptr->prev_rates_HD, dptr->rates_HD, pptr->num_HD_cells*sizeof(float));
		}
#pragma omp section
		{
			memcpy(&(cptr->prev_weights_RC[0][0]), &(cptr->RC_weights[0][0]), pptr->num_HD_cells*cptr->num_RC_connections*sizeof(float));
		}
	}
	return;
}



void calculate_activations_excitatory(data *dptr, params *pptr, connect *cptr, cArray **bptr, int timestep)
{
	int post_cell, synapse, connection;
	float HD_coefficient;
	float inhibitory, recurrent;
	float RC_scale;

	
	int cArray_index;
	
	cArray_index = 0;
	
	HD_coefficient = pptr->timestep_size/pptr->tau_HD;
		
    RC_scale = pptr->phi_RC/(float)cptr->num_RC_connections;

	
	inhibitory = 0.0;
	
	for(synapse = 0; synapse<pptr->num_HD_cells; synapse++)
	{
		inhibitory += pptr->global_inhibition * dptr->prev_rates_HD[synapse];   //Calculating global inhibition
	}
	

	
#pragma omp for private(post_cell, recurrent, connection, synapse)
	for(post_cell=0; post_cell<pptr->num_HD_cells; post_cell++)
	{
		recurrent = 0.0;
				
		for(connection=0; connection<cptr->num_RC_connections; connection++)
		{
			synapse = cptr->RC_connections[post_cell][connection];  //This is for random incomplete connectivity - tells you identity of presynaptic cell for each connection onto current postsynaptic cell
			
			recurrent += cptr->RC_weights[post_cell][connection] * read_RC_buffer(bptr, synapse); //Reads delayed RC input
		}
		
		dptr->activations_HD[post_cell] = ((1.0 - HD_coefficient) * dptr->prev_activations_HD[post_cell])
		+ (HD_coefficient * (RC_scale * recurrent))
		+ (HD_coefficient * dptr->input[post_cell])
		- (HD_coefficient * inhibitory)
		- (HD_coefficient * pptr->extern_inhib)    //This line is older code and is typically unused.
        + (HD_coefficient * dptr->combined_yaw_input[post_cell]);
       	}
		
	return;
}


void calculate_rates_excitatory(data *dptr, params *pptr)
{
    
    /* This function converts HD activations to rates: Option of simoid or hyperbolic tangent transfer function available */
    
	int cell;
	
	if(pptr->sigmoid)
	   {

#pragma omp for private(cell)
	for(cell=0; cell<pptr->num_HD_cells; cell++)
	{
		
	   dptr->rates_HD[cell] = 1.0 / (1.0 + exp(-2.0 * pptr->beta_HD * (dptr->activations_HD[cell] - pptr->alpha_HD)));
	
	}
		}
	   
	   else
	   {
#pragma omp for private(cell)
	   for(cell=0; cell<pptr->num_HD_cells; cell++)
	   {
		   
		   dptr->rates_HD[cell] = (float)tanh(dptr->activations_HD[cell]);
		   
		   if(dptr->rates_HD[cell] < 0.0)
			   dptr->rates_HD[cell] = 0.0;
		   
	   }
		   
	   }
	
	return;
}

void record_activations_rates(data *dptr, params *pptr, int index)
{
    /* This function records HD activations and rates over the entire timecourse of simulation */
    
	int cell;
	
	for(cell=0; cell<pptr->num_HD_cells; cell++)
	{
		dptr->activations_HD_time[cell][index] = dptr->activations_HD[cell];
		dptr->rates_HD_time[cell][index] = dptr->rates_HD[cell];
	}
	
	return;
}

void record_input(data *dptr, params *pptr, int index)
{
    /* This function records network input over entire timecourse of simulation */
    
	int cell;
	
	for(cell=0; cell<pptr->num_HD_cells; cell++)
	{
		dptr->input_location_time[cell][index] = dptr->input[cell];
	}
    
    for(cell=0; cell<pptr->num_HD_cells; cell++)
    {
        dptr->combined_yaw_input_location_time[cell][index] = dptr->combined_yaw_input[cell];
    }
   
	return;
}


void normalise_RC_weights(connect *cptr, params *pptr)
{
    /*This function normalises the RC weight packet*/
    
	int post_cell, pre_cell;
	
	float sumsq;
#pragma omp for private(post_cell, pre_cell, sumsq)	
	for(post_cell = 0; post_cell<pptr->num_HD_cells; post_cell++)
	{
		sumsq = 0.0;
		
		for(pre_cell=0;pre_cell<cptr->num_RC_connections;pre_cell++)
		{
			sumsq+= cptr->RC_weights[post_cell][pre_cell]*cptr->RC_weights[post_cell][pre_cell];
		}
		
		if(sumsq==0.0)
			continue;
		
		sumsq = sqrt(sumsq);
		
		for(pre_cell=0;pre_cell<cptr->num_RC_connections;pre_cell++)
		{
			cptr->RC_weights[post_cell][pre_cell] = cptr->RC_weights[post_cell][pre_cell]/sumsq;
		}
	}
	
	return;
	
}

void fill_RC_buffer(data *dptr, params *pptr, cArray **bptr, int cell)
{
    /* This function loads a circular buffer storing presynaptic firing (for simulating conduction delay */
    
	bptr[cell]->array[bptr[cell]->index++] = dptr->rates_HD[cell];
	
	if(bptr[cell]->index == pptr->conduction_buffer_size)
	{
		bptr[cell]->index = 0;
	}
	
	return;
	
}

float calculate_azimuth_angle(float x, float y)
{
    
    /*This function gets azimuth angle, in radians, of coordinate*/
    
    float azimuth_angle;
  
    azimuth_angle = atan2f(y,x);
    
    if(azimuth_angle < 0)
       azimuth_angle = azimuth_angle + (2*PI);  //putting in 0 to 2pi range
    
    
    return azimuth_angle;

}

double dotProduct(double Vectora[3], double Vectorb[3])
{
    double dP =(Vectora[0]*Vectorb[0])+(Vectora[1]*Vectorb[1])+(Vectora[2]*Vectorb[2]);
    
    return dP;
}


void ALT_egocentric_yaw_increment(data *dptr, params *pptr, int timestep)
{
    
    
    int current_path_timestep = timestep - pptr->input_timesteps;
    double surface_normal1[3];
    double surface_normal2[3];
    int index;
    double zaxis[3] = {0, 0, 1};
    
    //DOING ROTATION FOR CURRENT TIMESTEP
    
    
    surface_normal1[0] = (2*dptr->positions_x[current_path_timestep]);
    surface_normal1[1] = (2*dptr->positions_y[current_path_timestep]);
    surface_normal1[2] = (2*dptr->positions_z[current_path_timestep]);
    
    
    //Finding point on current plane where z=1 NEW REFERENCE VECTOR METHOD
    
    double curr_y, curr_x;
    double curr_z = 1.0;
    
    if(!surface_normal1[0])
    {
        curr_y = dptr->positions_y[current_path_timestep] + (((surface_normal1[2]*dptr->positions_z[current_path_timestep]) - surface_normal1[2])/surface_normal1[1]);
        curr_x = 0.0;
    }else if(!surface_normal1[1])
    {
        curr_y = 0.0;
        curr_x = dptr->positions_x[current_path_timestep] + ((surface_normal1[2] * dptr->positions_z[current_path_timestep])/surface_normal1[0]) - (surface_normal1[2]/surface_normal1[0]);
    }
    else{
        
        curr_x = ((surface_normal1[0]*dptr->positions_x[current_path_timestep]) + (surface_normal1[1]*dptr->positions_y[current_path_timestep]) + (surface_normal1[2]*dptr->positions_z[current_path_timestep]) - surface_normal1[2])/(surface_normal1[0] + ((surface_normal1[1]*dptr->positions_y[current_path_timestep])/dptr->positions_x[current_path_timestep]));
        
        curr_y = (dptr->positions_y[current_path_timestep]/dptr->positions_x[current_path_timestep]) * curr_x;
    }
    
    double refVec1[3] = {(curr_x-dptr->positions_x[current_path_timestep]),(curr_y-dptr->positions_y[current_path_timestep]), (curr_z-dptr->positions_z[current_path_timestep])};
    double length_RV1 = sqrt((refVec1[0]*refVec1[0])+(refVec1[1]*refVec1[1])+(refVec1[2]*refVec1[2]));
    refVec1[0] = refVec1[0]/length_RV1;
    refVec1[1] = refVec1[1]/length_RV1;
    refVec1[2] = refVec1[2]/length_RV1;
    
    
    for(index=0; index<3; index++)
    {
        dptr->refVec_time[index][current_path_timestep] = refVec1[index];
    }

    
    double length_sn1 = sqrt((surface_normal1[0]*surface_normal1[0])+(surface_normal1[1]*surface_normal1[1])+(surface_normal1[2]*surface_normal1[2]));
    surface_normal1[0] = surface_normal1[0]/length_sn1;
    surface_normal1[1] = surface_normal1[1]/length_sn1;
    surface_normal1[2] = surface_normal1[2]/length_sn1;

    
    double cross_product_1[3] = {surface_normal1[1]*zaxis[2]-surface_normal1[2]*zaxis[1], surface_normal1[2]*zaxis[0]-surface_normal1[0]*zaxis[2],surface_normal1[0]*zaxis[1]-surface_normal1[1]*zaxis[0]};
    double norm_cross_product_1 = sqrt((cross_product_1[0]*cross_product_1[0])+(cross_product_1[1]*cross_product_1[1])+(cross_product_1[2]*cross_product_1[2]));

    //finding rotation to get surface_normal to be aligned with z axis
    double rot_angle = atan2(norm_cross_product_1, dotProduct(surface_normal1,zaxis));
    double rot_axis[3] = {zaxis[1]*surface_normal1[2]-zaxis[2]*surface_normal1[1], zaxis[2]*surface_normal1[0]-zaxis[0]*surface_normal1[2],zaxis[0]*surface_normal1[1]-zaxis[1]*surface_normal1[0]};

    //Now rotating normalised current heading down to xy plane (i.e. by rotation that takes surface normal into alignment with z axis)
    
    double current_heading[3] = {dptr->headings_x[current_path_timestep],dptr->headings_y[current_path_timestep],dptr->headings_z[current_path_timestep]};
    double length_current_heading = sqrt((current_heading[0]*current_heading[0])+(current_heading[1]*current_heading[1])+(current_heading[2]*current_heading[2]));
    current_heading[0] = current_heading[0]/length_current_heading;
    current_heading[1] = current_heading[1]/length_current_heading;
    current_heading[2] = current_heading[2]/length_current_heading;
    
    //normalise current heading is intial quaternion
    double initial_quat[4] = {0,current_heading[0], current_heading[1],current_heading[2]};

    //generating rotation quaternion
    double rotation_quat[4] = {0,0,0,0};
    rotation_quat[0] = cos(rot_angle/2);
    rotation_quat[1] = (sin(rot_angle/2)) * rot_axis[0];
    rotation_quat[2] = (sin(rot_angle/2)) * rot_axis[1];
    rotation_quat[3] = (sin(rot_angle/2)) * rot_axis[2];
    
    //intermediate quat is rotation * initial
    double intermediate_quat[4] = {0,0,0,0};
    
    intermediate_quat[0] = (rotation_quat[0] * initial_quat[0]) - (rotation_quat[1] * initial_quat[1])
    - (rotation_quat[2] * initial_quat[2]) - (rotation_quat[3] * initial_quat[3]);
    
    intermediate_quat[1] = (rotation_quat[0] * initial_quat[1]) + (rotation_quat[1] * initial_quat[0])
    - (rotation_quat[2] * initial_quat[3]) + (rotation_quat[3] * initial_quat[2]);
    
    intermediate_quat[2] = (rotation_quat[0] * initial_quat[2]) + (rotation_quat[1] * initial_quat[3])
    + (rotation_quat[2] * initial_quat[0]) - (rotation_quat[3] * initial_quat[1]);
    
    intermediate_quat[3] = (rotation_quat[0] * initial_quat[3]) - (rotation_quat[1] * initial_quat[2])
    + (rotation_quat[2] * initial_quat[1]) + (rotation_quat[3] * initial_quat[0]);

    
    //now finding inverse rotation quat
    
    double sum = 0;
    for(index=0;index<4;index++)
    {
        sum = sum + rotation_quat[index] * rotation_quat[index];
    }
    
    double rotation_quat_magnitude = sqrt(sum);
    
    //Input quat conjugate
    double rotation_quat_conjugate[4] = {0,0,0,0};;
    rotation_quat_conjugate[0] = rotation_quat[0];
    
    for(index=1;index<4;index++)
    {
        rotation_quat_conjugate[index] = -rotation_quat[index];
    }
    
    //Now using above to inverse rotation quaternion
    double inverse_quat[4] = {0,0,0,0};
    
    for(index=0;index<4;index++)
    {
        inverse_quat[index] = rotation_quat_conjugate[index]/rotation_quat_magnitude;
    }
    
    //now multiplying intermediate quat by inverse quat
    double final_quat[4] = {0,0,0,0};
    
    final_quat[0] = (intermediate_quat[0] * inverse_quat[0]) - (intermediate_quat[1] * inverse_quat[1])
    - (intermediate_quat[2] * inverse_quat[2]) - (intermediate_quat[3] * inverse_quat[3]);
    
    final_quat[1] = (intermediate_quat[0] * inverse_quat[1]) + (intermediate_quat[1] * inverse_quat[0])
    - (intermediate_quat[2] * inverse_quat[3]) + (intermediate_quat[3] * inverse_quat[2]);
    
    final_quat[2] = (intermediate_quat[0] * inverse_quat[2]) + (intermediate_quat[1] * inverse_quat[3])
    + (intermediate_quat[2] * inverse_quat[0]) - (intermediate_quat[3] * inverse_quat[1]);
    
    final_quat[3] = (intermediate_quat[0] * inverse_quat[3]) - (intermediate_quat[1] * inverse_quat[2])
    + (intermediate_quat[2] * inverse_quat[1]) + (intermediate_quat[3] * inverse_quat[0]);

    double h1_azi = calculate_azimuth_angle(final_quat[1], final_quat[2]);
    
    
    for(index=0; index<3; index++)
    {
        dptr->rotated_heading_time[index][current_path_timestep] = final_quat[index+1];
    }

    
    //now rotating refVec1 by the same angle about same axis//
    ///
    ///
    ///
    ///

    //normalise current heading is intial quaternion
    initial_quat[0] = 0;
    initial_quat[1] = refVec1[0];
    initial_quat[2] = refVec1[1];
    initial_quat[3] = refVec1[2];
    
    //generating rotation quaternion
    rotation_quat[0] = cos(rot_angle/2);
    rotation_quat[1] = (sin(rot_angle/2)) * rot_axis[0];
    rotation_quat[2] = (sin(rot_angle/2)) * rot_axis[1];
    rotation_quat[3] = (sin(rot_angle/2)) * rot_axis[2];
    
    //intermediate quat is rotation * initial
    
    
    intermediate_quat[0] = (rotation_quat[0] * initial_quat[0]) - (rotation_quat[1] * initial_quat[1])
    - (rotation_quat[2] * initial_quat[2]) - (rotation_quat[3] * initial_quat[3]);
    
    intermediate_quat[1] = (rotation_quat[0] * initial_quat[1]) + (rotation_quat[1] * initial_quat[0])
    - (rotation_quat[2] * initial_quat[3]) + (rotation_quat[3] * initial_quat[2]);
    
    intermediate_quat[2] = (rotation_quat[0] * initial_quat[2]) + (rotation_quat[1] * initial_quat[3])
    + (rotation_quat[2] * initial_quat[0]) - (rotation_quat[3] * initial_quat[1]);
    
    intermediate_quat[3] = (rotation_quat[0] * initial_quat[3]) - (rotation_quat[1] * initial_quat[2])
    + (rotation_quat[2] * initial_quat[1]) + (rotation_quat[3] * initial_quat[0]);
    
    
    //now finding inverse rotation quat
    
    sum = 0;
    for(index=0;index<4;index++)
    {
        sum = sum + rotation_quat[index] * rotation_quat[index];
    }
    
    rotation_quat_magnitude = sqrt(sum);
    
    //Input quat conjugate
    
    rotation_quat_conjugate[0] = rotation_quat[0];
    
    for(index=1;index<4;index++)
    {
        rotation_quat_conjugate[index] = -rotation_quat[index];
    }
    
    //Now using above to inverse rotation quaternion
    
    
    for(index=0;index<4;index++)
    {
        inverse_quat[index] = rotation_quat_conjugate[index]/rotation_quat_magnitude;
    }
    
    //now multiplying intermediate quat by inverse quat
    
    
    final_quat[0] = (intermediate_quat[0] * inverse_quat[0]) - (intermediate_quat[1] * inverse_quat[1])
    - (intermediate_quat[2] * inverse_quat[2]) - (intermediate_quat[3] * inverse_quat[3]);
    
    final_quat[1] = (intermediate_quat[0] * inverse_quat[1]) + (intermediate_quat[1] * inverse_quat[0])
    - (intermediate_quat[2] * inverse_quat[3]) + (intermediate_quat[3] * inverse_quat[2]);
    
    final_quat[2] = (intermediate_quat[0] * inverse_quat[2]) + (intermediate_quat[1] * inverse_quat[3])
    + (intermediate_quat[2] * inverse_quat[0]) - (intermediate_quat[3] * inverse_quat[1]);
    
    final_quat[3] = (intermediate_quat[0] * inverse_quat[3]) - (intermediate_quat[1] * inverse_quat[2])
    + (intermediate_quat[2] * inverse_quat[1]) + (intermediate_quat[3] * inverse_quat[0]);
    
    double R1_azi = calculate_azimuth_angle(final_quat[1], final_quat[2]);

    for(index=0; index<3; index++)
    {
        dptr->rotated_refVec_time[index][current_path_timestep] = final_quat[index+1];
    }
    
    
    printf("\n\nH1 AZI = %f R1 AZI = %f\n\n", h1_azi *(180/PI), R1_azi * (180/PI));
    fflush(stdout);
    
    
    //DOING ROTATING FOR NEXT TIMESTEP
    surface_normal2[0] = (2*dptr->positions_x[current_path_timestep+1]);
    surface_normal2[1] = (2*dptr->positions_y[current_path_timestep+1]);
    surface_normal2[2] = (2*dptr->positions_z[current_path_timestep+1]);
    
    
    //Finding point on current plane where z=1 NEW REFERENCE VECTOR METHOD
    
    double curr_y2, curr_x2;
    double curr_z2 = 1.0;
    
    if(!surface_normal2[0])
    {
        curr_y2 = dptr->positions_y[current_path_timestep+1] + (((surface_normal1[2]*dptr->positions_z[current_path_timestep+1]) - surface_normal1[2])/surface_normal1[1]);
        curr_x2 = 0.0;
    }else if(!surface_normal1[1])
    {
        curr_y2 = 0.0;
        curr_x2 = dptr->positions_x[current_path_timestep+1] + ((surface_normal1[2] * dptr->positions_z[current_path_timestep+1])/surface_normal1[0]) - (surface_normal1[2]/surface_normal1[0]);
    }
    else{
        
        curr_x2 = ((surface_normal1[0]*dptr->positions_x[current_path_timestep+1]) + (surface_normal1[1]*dptr->positions_y[current_path_timestep+1]) + (surface_normal1[2]*dptr->positions_z[current_path_timestep+1]) - surface_normal1[2])/(surface_normal1[0] + ((surface_normal1[1]*dptr->positions_y[current_path_timestep+1])/dptr->positions_x[current_path_timestep+1]));
        
        curr_y2 = (dptr->positions_y[current_path_timestep+1]/dptr->positions_x[current_path_timestep+1]) * curr_x2;
    }
    
    double refVec2[3] = {(curr_x2-dptr->positions_x[current_path_timestep+1]),(curr_y2-dptr->positions_y[current_path_timestep+1]), (curr_z2-dptr->positions_z[current_path_timestep+1])};
    double length_RV2 = sqrt((refVec2[0]*refVec2[0])+(refVec2[1]*refVec2[1])+(refVec2[2]*refVec2[2]));
    refVec2[0] = refVec2[0]/length_RV2;
    refVec2[1] = refVec2[1]/length_RV2;
    refVec2[2] = refVec2[2]/length_RV2;
    
    
    double length_sn2 = sqrt((surface_normal2[0]*surface_normal2[0])+(surface_normal2[1]*surface_normal2[1])+(surface_normal2[2]*surface_normal2[2]));
    surface_normal2[0] = surface_normal2[0]/length_sn2;
    surface_normal2[1] = surface_normal2[1]/length_sn2;
    surface_normal2[2] = surface_normal2[2]/length_sn2;
    
    
    double cross_product_2[3] = {surface_normal2[1]*zaxis[2]-surface_normal2[2]*zaxis[1], surface_normal2[2]*zaxis[0]-surface_normal2[0]*zaxis[2],surface_normal2[0]*zaxis[1]-surface_normal2[1]*zaxis[0]};
    double norm_cross_product_2 = sqrt((cross_product_2[0]*cross_product_2[0])+(cross_product_2[1]*cross_product_2[1])+(cross_product_2[2]*cross_product_2[2]));
    
    //finding rotation to get surface_normal2 to be aligned with z axis
    double rot_angle2 = atan2(norm_cross_product_2, dotProduct(surface_normal2,zaxis));
    double rot_axis2[3] = {zaxis[1]*surface_normal2[2]-zaxis[2]*surface_normal2[1], zaxis[2]*surface_normal2[0]-zaxis[0]*surface_normal2[2],zaxis[0]*surface_normal2[1]-zaxis[1]*surface_normal2[0]};
    
    //Now rotating normalised current heading down to xy plane (i.e. by rotation that takes surface normal into alignment with z axis)
    
    double next_heading[3] = {dptr->headings_x[current_path_timestep+1],dptr->headings_y[current_path_timestep+1],dptr->headings_z[current_path_timestep+1]};
    double length_next_heading = sqrt((next_heading[0]*next_heading[0])+(next_heading[1]*next_heading[1])+(next_heading[2]*next_heading[2]));
    next_heading[0] = next_heading[0]/length_next_heading;
    next_heading[1] = next_heading[1]/length_next_heading;
    next_heading[2] = next_heading[2]/length_next_heading;
    
    //normalise next heading is intial quaternion
    double initial_quat2[4] = {0,next_heading[0], next_heading[1],next_heading[2]};
    
    //generating rotation quaternion
    double rotation_quat2[4] = {0,0,0,0};
    rotation_quat2[0] = cos(rot_angle2/2);
    rotation_quat2[1] = (sin(rot_angle2/2)) * rot_axis2[0];
    rotation_quat2[2] = (sin(rot_angle2/2)) * rot_axis2[1];
    rotation_quat2[3] = (sin(rot_angle2/2)) * rot_axis2[2];
    
    //intermediate quat is rotation * initial
    double intermediate_quat2[4] = {0,0,0,0};
    
    intermediate_quat2[0] = (rotation_quat2[0] * initial_quat2[0]) - (rotation_quat2[1] * initial_quat2[1])
    - (rotation_quat2[2] * initial_quat2[2]) - (rotation_quat2[3] * initial_quat2[3]);
    
    intermediate_quat2[1] = (rotation_quat2[0] * initial_quat2[1]) + (rotation_quat2[1] * initial_quat2[0])
    - (rotation_quat2[2] * initial_quat2[3]) + (rotation_quat2[3] * initial_quat2[2]);
    
    intermediate_quat2[2] = (rotation_quat2[0] * initial_quat2[2]) + (rotation_quat2[1] * initial_quat2[3])
    + (rotation_quat2[2] * initial_quat2[0]) - (rotation_quat2[3] * initial_quat2[1]);
    
    intermediate_quat2[3] = (rotation_quat2[0] * initial_quat2[3]) - (rotation_quat2[1] * initial_quat2[2])
    + (rotation_quat2[2] * initial_quat2[1]) + (rotation_quat2[3] * initial_quat2[0]);
    
    
    //now finding inverse rotation quat2
    
    sum = 0;
    for(index=0;index<4;index++)
    {
        sum = sum + rotation_quat2[index] * rotation_quat2[index];
    }
    
    double rotation_quat2_magnitude = sqrt(sum);
    
    //Input quat2 conjugate
    double rotation_quat2_conjugate[4] = {0,0,0,0};
    rotation_quat2_conjugate[0] = rotation_quat2[0];
    
    for(index=1;index<4;index++)
    {
        rotation_quat2_conjugate[index] = -rotation_quat2[index];
    }
    
    //Now using above to inverse rotation quat2ernion
    double inverse_quat2[4] = {0,0,0,0};
    
    for(index=0;index<4;index++)
    {
        inverse_quat2[index] = rotation_quat2_conjugate[index]/rotation_quat2_magnitude;
    }
    
    //now multiplying intermediate quat2 by inverse quat2
    double final_quat2[4] = {0,0,0,0};
    
    final_quat2[0] = (intermediate_quat2[0] * inverse_quat2[0]) - (intermediate_quat2[1] * inverse_quat2[1])
    - (intermediate_quat2[2] * inverse_quat2[2]) - (intermediate_quat2[3] * inverse_quat2[3]);
    
    final_quat2[1] = (intermediate_quat2[0] * inverse_quat2[1]) + (intermediate_quat2[1] * inverse_quat2[0])
    - (intermediate_quat2[2] * inverse_quat2[3]) + (intermediate_quat2[3] * inverse_quat2[2]);
    
    final_quat2[2] = (intermediate_quat2[0] * inverse_quat2[2]) + (intermediate_quat2[1] * inverse_quat2[3])
    + (intermediate_quat2[2] * inverse_quat2[0]) - (intermediate_quat2[3] * inverse_quat2[1]);
    
    final_quat2[3] = (intermediate_quat2[0] * inverse_quat2[3]) - (intermediate_quat2[1] * inverse_quat2[2])
    + (intermediate_quat2[2] * inverse_quat2[1]) + (intermediate_quat2[3] * inverse_quat2[0]);
    
    double h2_azi = calculate_azimuth_angle(final_quat2[1], final_quat2[2]);
    
    //now rotating refVec2 by the same angle about same axis//
    ///
    ///
    ///
    ///
    
    //normalised current heading is intial quaternion
    initial_quat2[0] = 0;
    initial_quat2[1] = refVec2[0];
    initial_quat2[2] = refVec2[1];
    initial_quat2[3] = refVec2[2];
    
    //generating rotation quaternion
    rotation_quat2[0] = cos(rot_angle2/2);
    rotation_quat2[1] = (sin(rot_angle2/2)) * rot_axis2[0];
    rotation_quat2[2] = (sin(rot_angle2/2)) * rot_axis2[1];
    rotation_quat2[3] = (sin(rot_angle2/2)) * rot_axis2[2];
    
    //intermediate quat is rotation * initial
    
    
    intermediate_quat2[0] = (rotation_quat2[0] * initial_quat2[0]) - (rotation_quat2[1] * initial_quat2[1])
    - (rotation_quat2[2] * initial_quat2[2]) - (rotation_quat2[3] * initial_quat2[3]);
    
    intermediate_quat2[1] = (rotation_quat2[0] * initial_quat2[1]) + (rotation_quat2[1] * initial_quat2[0])
    - (rotation_quat2[2] * initial_quat2[3]) + (rotation_quat2[3] * initial_quat2[2]);
    
    intermediate_quat2[2] = (rotation_quat2[0] * initial_quat2[2]) + (rotation_quat2[1] * initial_quat2[3])
    + (rotation_quat2[2] * initial_quat2[0]) - (rotation_quat2[3] * initial_quat2[1]);
    
    intermediate_quat2[3] = (rotation_quat2[0] * initial_quat2[3]) - (rotation_quat2[1] * initial_quat2[2])
    + (rotation_quat2[2] * initial_quat2[1]) + (rotation_quat2[3] * initial_quat2[0]);
    
    
    //now finding inverse rotation quat2
    
    sum = 0;
    for(index=0;index<4;index++)
    {
        sum = sum + rotation_quat2[index] * rotation_quat2[index];
    }
    
    rotation_quat2_magnitude = sqrt(sum);
    
    //Input quat2 conjugate
    
    rotation_quat2_conjugate[0] = rotation_quat2[0];
    
    for(index=1;index<4;index++)
    {
        rotation_quat2_conjugate[index] = -rotation_quat2[index];
    }
    
    //Now using above to inverse rotation quat2ernion
    
    
    for(index=0;index<4;index++)
    {
        inverse_quat2[index] = rotation_quat2_conjugate[index]/rotation_quat2_magnitude;
    }
    
    //now multiplying intermediate quat2 by inverse quat2
    
    
    final_quat2[0] = (intermediate_quat2[0] * inverse_quat2[0]) - (intermediate_quat2[1] * inverse_quat2[1])
    - (intermediate_quat2[2] * inverse_quat2[2]) - (intermediate_quat2[3] * inverse_quat2[3]);
    
    final_quat2[1] = (intermediate_quat2[0] * inverse_quat2[1]) + (intermediate_quat2[1] * inverse_quat2[0])
    - (intermediate_quat2[2] * inverse_quat2[3]) + (intermediate_quat2[3] * inverse_quat2[2]);
    
    final_quat2[2] = (intermediate_quat2[0] * inverse_quat2[2]) + (intermediate_quat2[1] * inverse_quat2[3])
    + (intermediate_quat2[2] * inverse_quat2[0]) - (intermediate_quat2[3] * inverse_quat2[1]);
    
    final_quat2[3] = (intermediate_quat2[0] * inverse_quat2[3]) - (intermediate_quat2[1] * inverse_quat2[2])
    + (intermediate_quat2[2] * inverse_quat2[1]) + (intermediate_quat2[3] * inverse_quat2[0]);
    
    double R2_azi = calculate_azimuth_angle(final_quat2[1], final_quat2[2]);
    
    printf("\n\nH2 AZI = %f R2 AZI = %f\n\n", h2_azi *(180/PI), R2_azi * (180/PI));
    fflush(stdout);
    
    
    
    double difference1 = h1_azi - R1_azi;
    double difference2 = h2_azi - R2_azi;

    
    if(difference1<0 && difference2>0)
        difference2-=(2*PI);
   
    if(difference1>0 && difference2<0)
        difference2+=(2*PI);
    
    /*
    if(difference1< - PI)
        difference1 += (2*PI);
    
    if(difference1 > PI)
        difference1-= (2*PI);
  

    
    if(difference2< - PI)
        difference2 += (2*PI);
    
    if(difference2 > PI)
        difference2-= (2*PI);
     */
    
    double egoinc;
    
    egoinc = difference2 - difference1;
    
    if(fabs(fabs(difference2-difference1))>PI)
    {
        if(difference2>0)
        {
            egoinc -= (2*PI);
        }else{
            egoinc += (2*PI);
        }
    }
    
    
    dptr->egocentric_increment = egoinc * (180/PI);

    if(dptr->positions_x[current_path_timestep]==0 && dptr->positions_y[current_path_timestep]==0)
    {
        dptr->egocentric_increment = (h1_azi - h2_azi) * (180/PI);
    }
    
    
    printf("\n\nD1: %f \nD2: %f \nEGOINC = %f\n\n", difference1*(180/PI), difference2*(180/PI), egoinc*(180/PI));
    fflush(stdout);

    
    return;

}

void set_egocentric_yaw_increment(data *dptr, params *pptr, int timestep)
{
    /*This function calculates the amount of egocentric head rotation this timestep, and increments accordingly*/
    
    int current_path_timestep = timestep - pptr->input_timesteps;
    double surface_normal1[3];
    double surface_normal2[3];
    double difference;
    int index;
    
    double current_heading[3] = {dptr->headings_x[current_path_timestep],dptr->headings_y[current_path_timestep],dptr->headings_z[current_path_timestep]};
    
    //Using surface normal to define the tangent plane equation - THIS WILL NEED TO CHANGE FOR CUBOID SIMULATION, SURF NORMAL IS DIFFERENT!
    
    surface_normal1[0] = (2*dptr->positions_x[current_path_timestep]);
    surface_normal1[1] = (2*dptr->positions_y[current_path_timestep]);
    surface_normal1[2] = (2*dptr->positions_z[current_path_timestep]);
    
    
    
    
    //Finding point on current plane where z=1 NEW REFERENCE VECTOR METHOD
    
    double curr_y, curr_x;
    double curr_z = 1.0;
    
    if(!surface_normal1[0])
    {
        curr_y = dptr->positions_y[current_path_timestep] + (((surface_normal1[2]*dptr->positions_z[current_path_timestep]) - surface_normal1[2])/surface_normal1[1]);
        curr_x = 0.0;
    }else if(!surface_normal1[1])
    {
        curr_y = 0.0;
        curr_x = dptr->positions_x[current_path_timestep] + ((surface_normal1[2] * dptr->positions_z[current_path_timestep])/surface_normal1[0]) - (surface_normal1[2]/surface_normal1[0]);
    }
    else{
        
        curr_x = ((surface_normal1[0]*dptr->positions_x[current_path_timestep]) + (surface_normal1[1]*dptr->positions_y[current_path_timestep]) + (surface_normal1[2]*dptr->positions_z[current_path_timestep]) - surface_normal1[2])/(surface_normal1[0] + ((surface_normal1[1]*dptr->positions_y[current_path_timestep])/dptr->positions_x[current_path_timestep]));
        
        curr_y = (dptr->positions_y[current_path_timestep]/dptr->positions_x[current_path_timestep]) * curr_x;
    }

     double refVec1[3] = {(curr_x-dptr->positions_x[current_path_timestep]),(curr_y-dptr->positions_y[current_path_timestep]), (curr_z-dptr->positions_z[current_path_timestep])};
    
    
    printf("\n\nRefVec1: %f %f %f", refVec1[0], refVec1[1], refVec1[2]);
    fflush(stdout);
    /*
    
    //Finding place on z axis were the tangent plane intersects, based on position in plane and surface normal
   
    double z1 = ((surface_normal1[0]*dptr->positions_x[current_path_timestep])+(surface_normal1[1]*dptr->positions_y[current_path_timestep])+(surface_normal1[2]*dptr->positions_z[current_path_timestep]))/surface_normal1[2];
    
    //Reference is vector from position to point where tangent plane intersects with the z axis
    
    double refVec1[3] = {(0-dptr->positions_x[current_path_timestep]),(0-dptr->positions_y[current_path_timestep]), (z1-dptr->positions_z[current_path_timestep])};
    */
    
    
    for(index=0; index<3; index++)
    {
        dptr->refVec_time[index][current_path_timestep] = refVec1[index];
    }
    

    
    //Now find the angle between current heading, and the reference vector
    
    //normalising refVec1 to unit length
    double length_RV1 = sqrt((refVec1[0]*refVec1[0])+(refVec1[1]*refVec1[1])+(refVec1[2]*refVec1[2]));
    refVec1[0] = refVec1[0]/length_RV1;
    refVec1[1] = refVec1[1]/length_RV1;
    refVec1[2] = refVec1[2]/length_RV1;
    
    //normalising current_heading to unit length
    double length_current_heading = sqrt((current_heading[0]*current_heading[0])+(current_heading[1]*current_heading[1])+(current_heading[2]*current_heading[2]));
    current_heading[0] = current_heading[0]/length_current_heading;
    current_heading[1] = current_heading[1]/length_current_heading;
    current_heading[2] = current_heading[2]/length_current_heading;
    
    printf("\n\nHeadLength pre norm: %f\nPostNorm: %f\n\n", length_current_heading, sqrt((current_heading[0]*current_heading[0])+(current_heading[1]*current_heading[1])+(current_heading[2]*current_heading[2])));
    fflush(stdout);
    
    printf("\n\nRefLength pre norm: %f\nPostNorm: %f\n\n", length_RV1, sqrt((refVec1[0]*refVec1[0])+(refVec1[1]*refVec1[1])+(refVec1[2]*refVec1[2])));
    fflush(stdout);

    
    //Taking cross product of unit vectors
    double cross_product_1[3] = {refVec1[1]*current_heading[2]-refVec1[2]*current_heading[1], refVec1[2]*current_heading[0]-refVec1[0]*current_heading[2],refVec1[0]*current_heading[1]-refVec1[1]*current_heading[0]};
   
    
    double norm_cross_product_1 = sqrt((cross_product_1[0]*cross_product_1[0])+(cross_product_1[1]*cross_product_1[1])+(cross_product_1[2]*cross_product_1[2]));
    
/*
    //normalising the surface normal to unit length
    double length_sn1 = sqrt((surface_normal1[0]*surface_normal1[0])+(surface_normal1[1]*surface_normal1[1])+(surface_normal1[2]*surface_normal1[2]));
    surface_normal1[0] = surface_normal1[0]/length_sn1;
    surface_normal1[1] = surface_normal1[1]/length_sn1;
    surface_normal1[2] = surface_normal1[2]/length_sn1;

 
    //now angle is atan2(determinant, dot product)
    double angle1 = atan2(dotProduct(surface_normal1,cross_product_1), dotProduct(refVec1,current_heading));

*/
    //unsigned angle between refVec and current_heading
    double angle1 = atan2(norm_cross_product_1, dotProduct(refVec1,current_heading));
    

    //making it signed, and then converting to range 0 to 2pi

    if(dotProduct(surface_normal1, cross_product_1)<0)
    {
        angle1 = -angle1 + (2*PI);    //make it negative based on cross product dot with surface normal.
    }

  
    //NOW DOING ABOVE PROCESS AGAIN FOR THE NEXT TIMESTEP
    
    
    double next_heading[3] = {(dptr->headings_x[current_path_timestep+1]),(dptr->headings_y[current_path_timestep+1]),(dptr->headings_z[current_path_timestep+1])};
 
    //Using surface normal to define the tangent plane equation - THIS WILL NEED TO CHANGE FOR CUBOID SIMULATION, SURF NORMAL IS DIFFERENT!
    
    surface_normal2[0] = (2*dptr->positions_x[current_path_timestep+1]);
    surface_normal2[1] = (2*dptr->positions_y[current_path_timestep+1]);
    surface_normal2[2] = (2*dptr->positions_z[current_path_timestep+1]);
    
   
    //Finding point on current plane where z=1 NEW REFERENCE VECTOR METHOD
    
    double next_x, next_y;
    double next_z = 1.0;
    
    if(!surface_normal2[0])
    {
        next_y = dptr->positions_y[current_path_timestep+1] + (((surface_normal2[2]*dptr->positions_z[current_path_timestep+1]) - surface_normal2[2])/surface_normal2[1]);
        next_x = 0.0;
    }else if(!surface_normal2[1])
    {
        next_y = 0.0;
        next_x = dptr->positions_x[current_path_timestep+1] + ((surface_normal2[2] * dptr->positions_z[current_path_timestep+1])/surface_normal2[0]) - (surface_normal2[2]/surface_normal2[0]);
    }
    else{
        
        next_x = ((surface_normal2[0]*dptr->positions_x[current_path_timestep+1]) + (surface_normal2[1]*dptr->positions_y[current_path_timestep+1]) + (surface_normal2[2]*dptr->positions_z[current_path_timestep+1]) - surface_normal2[2])/(surface_normal2[0] + ((surface_normal2[1]*dptr->positions_y[current_path_timestep+1])/dptr->positions_x[current_path_timestep+1]));
        
        next_y = (dptr->positions_y[current_path_timestep+1]/dptr->positions_x[current_path_timestep+1]) * next_x;
    }
    
    double refVec2[3] = {(next_x-dptr->positions_x[current_path_timestep+1]),(next_y-dptr->positions_y[current_path_timestep+1]), (next_z-dptr->positions_z[current_path_timestep+1])};
    
    
    
    /*

    //Finding place on z axis were the tangent plane intersects, using surf norm and plane position at that point
    
    double z2 = ((surface_normal2[0]*dptr->positions_x[current_path_timestep+1])+(surface_normal2[1]*dptr->positions_y[current_path_timestep+1])+(surface_normal2[2]*dptr->positions_z[current_path_timestep+1]))/surface_normal2[2];
    
    //Reference is vector from position to point where tangent plane intersects with the z axis
    
    double refVec2[3] = {(0-dptr->positions_x[current_path_timestep+1]),(0-dptr->positions_y[current_path_timestep+1]), (z2-dptr->positions_z[current_path_timestep+1])};
    */
    
    //normalising refVec2 to unit length
    double length_RV2 = sqrt((refVec2[0]*refVec2[0])+(refVec2[1]*refVec2[1])+(refVec2[2]*refVec2[2]));
    refVec2[0] = refVec2[0]/length_RV2;
    refVec2[1] = refVec2[1]/length_RV2;
    refVec2[2] = refVec2[2]/length_RV2;
    
    //normalising next_heading to unit length
    double length_next_heading = sqrt((next_heading[0]*next_heading[0])+(next_heading[1]*next_heading[1])+(next_heading[2]*next_heading[2]));
    next_heading[0] = next_heading[0]/length_next_heading;
    next_heading[1] = next_heading[1]/length_next_heading;
    next_heading[2] = next_heading[2]/length_next_heading;

    double cross_product_2[3] = {refVec2[1]*next_heading[2]-refVec2[2]*next_heading[1], refVec2[2]*next_heading[0]-refVec2[0]*next_heading[2],refVec2[0]*next_heading[1]-refVec2[1]*next_heading[0]};
    //double cross_product_2[3] = {next_heading[1]*refVec2[2]-next_heading[2]*refVec2[1], next_heading[2]*refVec2[0]-next_heading[0]*refVec2[2],next_heading[0]*refVec2[1]-next_heading[1]*refVec2[0]};
    
    double norm_cross_product_2 = sqrt(((cross_product_2[0]*cross_product_2[0])+(cross_product_2[1]*cross_product_2[1])+(cross_product_2[2]*cross_product_2[2])));
    
    
    printf("\n\nHeadLength2 pre norm: %f\nPostNorm: %f\n\n", length_next_heading, sqrt((next_heading[0]*next_heading[0])+(next_heading[1]*next_heading[1])+(next_heading[2]*next_heading[2])));
    fflush(stdout);
    
    printf("\n\nRefLength2 pre norm: %f\nPostNorm: %f\n\n", length_RV2, sqrt((refVec2[0]*refVec2[0])+(refVec2[1]*refVec2[1])+(refVec2[2]*refVec2[2])));
    fflush(stdout);

/*
    //normalising the surface normal to unit length
    double length_sn2 = sqrt((surface_normal2[0]*surface_normal2[0])+(surface_normal2[1]*surface_normal2[1])+(surface_normal2[2]*surface_normal2[2]));
    surface_normal2[0] = surface_normal2[0]/length_sn2;
    surface_normal2[1] = surface_normal2[1]/length_sn2;
    surface_normal2[2] = surface_normal2[2]/length_sn2;

    cross_product_2[0] = cross_product_2[0]/norm_cross_product_2;
    cross_product_2[1] = cross_product_2[1]/norm_cross_product_2;
    cross_product_2[2] = cross_product_2[2]/norm_cross_product_2;

    //now angle is atan2(determinant, dot product)
    double angle2 = atan2(dotProduct(surface_normal2,cross_product_2), dotProduct(refVec2,next_heading));

*/
    //unsigned angle between refVec2 and next_heading
    double angle2 = atan2(norm_cross_product_2, dotProduct(refVec2,next_heading));
    

    //making it signed, and then converting to range 0 to 2pi
    if(dotProduct(surface_normal2, cross_product_2)<0)
    {
        angle2 = -angle2  + (2*PI);    //make it negative based on cross product dot with surface normal.
    }

/*
    if(dptr->positions_x[current_path_timestep] == 0 && dptr->positions_y[current_path_timestep] == 0 && dptr->positions_z[current_path_timestep] == 1)
    {
        angle1 = PI;
    }
 */
    //The egocentric increment will be the change in angle relative to reference vector at each timestep (i.e. angle2 - angle1)
    
    
    difference = angle2 - angle1;
    
    if(fabs(fabs(angle2-angle1))>PI)
    {
        if(angle2>angle1)
        {
            difference -= (2*PI);
            
        }else{
            difference += (2*PI);
        }
    }
    
    
    dptr->egocentric_increment = difference * (180/PI);
    
    if(dptr->positions_x[current_path_timestep]==0 && dptr->positions_y[current_path_timestep]==0)
    {
        
        dptr->egocentric_increment = acos(dotProduct(current_heading, next_heading)) * (180/PI);
    }
    
    //difference = atan2(sin(angle2-angle1),cos(angle2-angle1));
    
    //dptr->egocentric_increment =  difference * (180/PI);
    
    printf("\n\nAngle1: %lf Angle2: %lf\n\n", angle1 * (180/PI), angle2 * (180/PI));
    fflush(stdout);
    
    dptr->start_angle[current_path_timestep] = angle1 * (180/PI);
    
    return;
}


void set_gravity_yaw_increment(data *dptr, params *pptr, int timestep)
{
    /* This function sets the gravity yaw increment, based on changes in position and heading generated in MATLAB script*/
    /*At each timestep, a vector is caulculated from the rat's position to the north pole of the hemisphere (or the north pole of the cuboid), the change in azimuth angle of this vector between timesteps is the gravity increment.*/

    double gravity_vector_change;
    double init_axis_azimuth, fin_axis_azimuth;
    double init_head_azimuth, fin_head_azimuth;
    double surface_normal1[3], surface_normal2[3];
    double length_sn1, length_sn2;
    
    //Gravity vector change is change in azimuth of DV axis (i.e. surface normal at current position)
    
    int current_path_timestep = timestep - pptr->input_timesteps;
    
    //Surface normal at current timestep
    surface_normal1[0] = 2*dptr->positions_x[current_path_timestep];
    surface_normal1[1] = 2*dptr->positions_y[current_path_timestep];
    surface_normal1[2] = 2*dptr->positions_z[current_path_timestep];
    
    //Normalising it
    length_sn1 = sqrt((surface_normal1[0]*surface_normal1[0])+ (surface_normal1[1]*surface_normal1[1]) + (surface_normal1[2]*surface_normal1[2]));
    
    surface_normal1[0] =  surface_normal1[0]/length_sn1;
    surface_normal1[1] =  surface_normal1[1]/length_sn1;
    surface_normal1[2] =  surface_normal1[2]/length_sn1;
    
    
    //Surface Normal at next timestep
    surface_normal2[0] = 2*dptr->positions_x[current_path_timestep+1];
    surface_normal2[1] = 2*dptr->positions_y[current_path_timestep+1];
    surface_normal2[2] = 2*dptr->positions_z[current_path_timestep+1];
    
    //Normalising it
    length_sn2 = sqrt((surface_normal2[0]*surface_normal2[0])+ (surface_normal2[1]*surface_normal2[1]) + (surface_normal2[2]*surface_normal2[2]));
    
    surface_normal2[0] =  surface_normal2[0]/length_sn2;
    surface_normal2[1] =  surface_normal2[1]/length_sn2;
    surface_normal2[2] =  surface_normal2[2]/length_sn2;
    
    
    //work out azimuth angles of gravity vectors in range 0 to 2pi

    if(surface_normal1[0] == 0 && surface_normal1[1] == 0 && surface_normal1[2]*length_sn1 ==2)
    {
        dptr->gravity_increment = 0.0;
        
    }else{
    
    init_axis_azimuth = calculate_azimuth_angle(surface_normal1[0], surface_normal1[1]);
    fin_axis_azimuth = calculate_azimuth_angle(surface_normal2[0], surface_normal2[1]);

    //work out raw azimuth angles in range 0 to 2pi - this code is useful to tell us what the combined increment SHOULD be....
    init_head_azimuth = calculate_azimuth_angle(dptr->headings_x[current_path_timestep], dptr->headings_y[current_path_timestep]);
    fin_head_azimuth = calculate_azimuth_angle(dptr->headings_x[current_path_timestep+1], dptr->headings_y[current_path_timestep+1]);
     
    dptr->ideal_increment[current_path_timestep] = atan2(sin(fin_head_azimuth-init_head_azimuth), cos(fin_head_azimuth-init_head_azimuth)) * (180.0/PI);
    
        printf("\n\nIdealInc: %f\n\n", dptr->ideal_increment[current_path_timestep]);
        fflush(stdout);
        
        
    //this should give the shortest angle in range  -pi to +pi
    
    gravity_vector_change = atan2(sin(fin_axis_azimuth-init_axis_azimuth), cos(fin_axis_azimuth-init_axis_azimuth));
    dptr->gravity_increment = gravity_vector_change * (180.0/PI);
    }
        printf("\n\nGravity Increment: %f\nInit azimuth: %f\nFin Azimuth: %f\n", dptr->gravity_increment, init_axis_azimuth * (180.0/PI), fin_axis_azimuth * (180.0/PI));
        fflush(stdout);
    
}

inline float read_RC_buffer(cArray **bptr, int cell)
{
    /* This function reads out presynaptic firing stored in conduction delay buffer */
    
	return bptr[cell]->array[bptr[cell]->index];
}


void calculate_pvector(data *dptr, params *pptr)
{
    
    int cell;
    float vector1 = 0.0;
    float vector2 = 0.0;
    
    for(cell=0; cell<pptr->num_HD_cells; cell++)
    {
        vector1+= dptr->rates_HD[cell] * (sinf(dptr->favoured_view[cell]*(PI/180.0)));
        vector2+= dptr->rates_HD[cell] * (cosf(dptr->favoured_view[cell]*(PI/180.0)));
    }
    
    if (vector1 > 0.0 && vector2 > 0.0)
    {
        dptr->pvector = (atanf((vector1/vector2)) * (180.0/PI));
    }
    else if(vector2 < 0.0)
    {
        dptr->pvector = (atanf((vector1/vector2)) * (180.0/PI)) + 180.0;
    }
    else
    {
        dptr->pvector = (atanf((vector1/vector2)) * (180.0/PI))+ 360.0;
    }	
    
    
    return;
}
	