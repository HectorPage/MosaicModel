#ifndef __CALCULATE_INCLUDED__
#define __CALCULATE_INCLUDED__

#include "data.h"

extern void zero_states(data *dptr, params *pptr, cArray **bptr, connect *cptr);
extern void record_previous_states(data *dptr, params *pptr, connect *cptr);

extern void set_input(data *dptr, params *pptr);
extern void zero_input(data *dptr, params *pptr);

extern void calculate_activations_excitatory(data *dptr, params *pptr, connect *cptr, cArray **bptr, int timestep);
extern void calculate_rates_excitatory(data *dptr, params *pptr);

extern void record_activations_rates(data *dptr, params *pptr, int index);
extern void record_input(data *dptr, params *pptr, int index);

extern void normalise_RC_weights(connect *cptr, params *pptr);

extern void fill_RC_buffer(data *dptr, params *pptr, cArray **bptr, int cell);
extern float read_RC_buffer(cArray **bptr, int cell);

extern void set_combined_yaw_input(data *dptr, params *pptr);

extern float calculate_azimuth_angle(float x, float y);

extern void set_gravity_yaw_increment(data *dptr, params *pptr, int timestep);

extern void set_egocentric_yaw_increment(data *dptr, params *pptr, int timestep);

extern void calculate_pvector(data *dptr, params *pptr);

extern double dotProduct(double Vectora[3], double Vectorb[3]);




#endif