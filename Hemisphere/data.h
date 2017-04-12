#ifndef __DATA_INCLUDED__
#define __DATA_INCLUDED__
#define PI 3.141592654


typedef struct {
	
	int timesteps;
	float time;
    
    int input_timesteps;
    float input_time;
    
    int path_timesteps;
    float path_time;

    int num_HD_cells;
		
	float phi_RC;
    float sigma_RC;

	float tau_HD;
	float alpha_HD;
	float beta_HD;
	
	float timestep_size;
	
	float extern_inhib;
	float global_inhibition;
	
	float input_str;
	float input_loc;
    float input_sigma;
    
    float combined_yaw_input_str;
    float combined_yaw_input_sigma;
   
    int normalise;
	
	float conduction_delay;
	int conduction_buffer_size;
	
	int threads;
	
	int sigmoid;
    
    int gravity_increment_switch;

}params;

typedef struct {
	
	float *activations_HD, *prev_activations_HD;

	float *rates_HD, *prev_rates_HD;
	
	float *input;
		
	float **rates_HD_time;
		
	float **activations_HD_time;
	
	float *favoured_view;
  	
	float *time;
	
	float sumsq;
	
	float **input_location_time;
    
    float *combined_yaw_input;
    
    float combined_yaw_input_loc;

    float **combined_yaw_input_location_time;
    
    float egocentric_velocity;  //stores the velocity of egocentric rotation
    
    float euler_angle; //stores the angle between start_heading and end_heading (taking location into account)
    
    int rotation_start_timestep;    //stores start of rotation, used for slerp
    
    float gravity_increment;       //this stores the gravity-based increment to yaw input
    float egocentric_increment;    //this stores the ahv-based increment to yaw input

    float pvector;  //stores current HD activity packet location
    
    float *positions_x;
    float *positions_y;
    float *positions_z;
    
    float *headings_x;
    float *headings_y;
    float *headings_z;
    
    float *pvector_time;
    
    double *start_angle;
    
    float *combined_input_pvector_time;
    float *ideal_increment;
    float *actual_increment;
    
    double **refVec_time;
    
    double **rotated_refVec_time;
    double **rotated_heading_time;
 
    
	
}data;

typedef struct {
	
	int num_RC_connections;
	
	int **RC_connections;
	
	float **RC_weights;
	float **prev_weights_RC;
	
	float **full_RC;
		
	
	int inhib_random, excite_random;
	
	unsigned long int random_seed;	
	
	
}connect;

typedef struct {
	
	int index;
	float *array;
}cArray;


#endif