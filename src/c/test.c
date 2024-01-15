#include <stdio.h>
#include <stdint.h>
#include <feos.h>
// typedef struct feos_eos feos_eos_t;
// extern feos_eos_t *feos_eos_from_json(const char *str);
// extern void feos_eos_free(feos_eos_t *);

// typedef struct feos_state feos_state_t;
// extern feos_state_t *feos_state_new_npt(const feos_eos_t *, double, double, const double*, size_t, const char *);
// extern void feos_state_free(feos_state_t *);
// extern double feos_state_pressure(feos_state_t *, size_t);

int main(void)
{
    char *str = "{ \
    \"residual_model\": \"PC-SAFT\", \
    \"ideal_gas_model\": \"\", \
    \"residual_substance_parameters\": [ \
        { \
            \"identifier\": { \
                \"cas\": \"74-82-8\", \
                \"name\": \"methane\", \
                \"iupac_name\": \"methane\", \
                \"smiles\": \"C\", \
                \"inchi\": \"InChI=1/CH4/h1H4\", \
                \"formula\": \"CH4\" \
            }, \
            \"model_record\": { \
                \"m\": 1.0, \
                \"sigma\": 3.7039, \
                \"epsilon_k\": 150.03 \
            }, \
            \"molarweight\": 16.043 \
        }, \
        { \
            \"identifier\": { \
                \"cas\": \"110-82-7\", \
                \"name\": \"cyclohexane\", \
                \"iupac_name\": \"cyclohexane\", \
                \"smiles\": \"C1CCCCC1\", \
                \"inchi\": \"InChI=1/C6H12/c1-2-4-6-5-3-1/h1-6H2\", \
                \"formula\": \"C6H12\" \
            }, \
            \"model_record\": { \
                \"m\": 2.5303, \
                \"sigma\": 3.8499, \
                \"epsilon_k\": 278.11 \
            }, \
            \"molarweight\": 84.147 \
        } \
    ], \
    \"residual_binary_parameters\": [ \
        { \
            \"id1\": { \
                \"cas\": \"67-56-1\", \
                \"name\": \"methanol\", \
                \"iupac_name\": \"methanol\", \
                \"smiles\": \"CO\", \
                \"inchi\": \"InChI=1/CH4O/c1-2/h2H,1H3\", \
                \"formula\": \"CH4O\" \
            }, \
            \"id2\": { \
                \"cas\": \"110-82-7\", \
                \"name\": \"cyclohexane\", \
                \"iupac_name\": \"cyclohexane\", \
                \"smiles\": \"C1CCCCC1\", \
                \"inchi\": \"InChI=1/C6H12/c1-2-4-6-5-3-1/h1-6H2\", \
                \"formula\": \"C6H12\" \
            }, \
            \"model_record\": { \
                \"k_ij\": 0.051 \
            } \
        } \
    ] \
}";
    char *str2 = "{ \
    \"residual_model\": \"PC-SAFT\", \
    \"ideal_gas_model\": \"\", \
    \"residual_substance_parameters\": [ \
        { \
            \"identifier\": { \
                \"cas\": \"74-82-8\", \
                \"name\": \"methane\", \
                \"iupac_name\": \"methane\", \
                \"smiles\": \"C\", \
                \"inchi\": \"InChI=1/CH4/h1H4\", \
                \"formula\": \"CH4\" \
            }, \
            \"model_record\": { \
                \"m\": 1.0, \
                \"sigma\": 3.7039, \
                \"epsilon_k\": 150.03 \
            }, \
            \"molarweight\": 16.043 \
        } \
    ] \
}";
    feos_equation_of_state_t *eos = feos_eos_from_json(str2);

    // state
    double temperature = 200.0; // K
    double pressure = 15.0;     // bar
    double moles[] = {1.0};     // mol for each component
    size_t n = sizeof(moles) / sizeof(size_t); // number of components

    feos_state_t *state = feos_state_new_npt(
        eos, // equation of state
        temperature, // temperature in K
        pressure, // pressure in bar
        moles, // number of mols
        n, // number of components
        "stable" // which phase to compute (stable, liquid, vapor)
    );

    double pressure_calculated = feos_state_pressure(state, 2);
    double density = feos_state_density(state);
    
    printf("State stable? %s\n", feos_state_is_stable(state) ? "true" : "false");
    printf("pressure: %f bar\n", pressure);
    printf("density: %f mol/m3\n", density);

    // free equation of state and exit
    feos_state_free(state);
    feos_eos_free(eos);

    return 0;
}