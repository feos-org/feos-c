#include <stdio.h>
#include <stdint.h>

typedef struct eos eos_t;
extern eos_t *eos_from_json(const char *str);
extern void free_eos_ptr(eos_t *);
extern double pressure_bar(eos_t *eos_ptr, double temperature, double density, const double *molefracs, size_t len);
extern double reduced_entropy(eos_t *eos_ptr, double temperature, double density, const double *molefracs, size_t len);
extern double da_dv(eos_t *eos_ptr, double temperature, double density, const double *molefracs, size_t len);
extern int get_Arxy(const long long int uuid, const int NT, const int ND, const double T, const double rho, const double* molefrac, const int Ncomp, double *val, char* errmsg, int errmsg_length);

int main(void)
{
    char *str = "{ \
    \"model\": \"PC-SAFT\", \
    \"substance_parameters\": [ \
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
    \"binary_parameters\": [ \
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
    eos_t *eos = eos_from_json(str);

    // state
    double temperature = 300.0; // K
    double density = 10000;   // mol / m³
    double molefracs[] = {0.1, 0.9};
    double rgas = 8.3145; // J / mol / K
    
    // error message
    int errmsg_length = 300;
    char errmsg[errmsg_length];

    // status and return value    
    int status;
    double value = 0.0;
    double a00 = 0.0, a10 = 0.0;

    status = get_Arxy(eos, 0, 1, temperature, density, molefracs, 2, &value, errmsg, 100);
    printf("Status: %d, 01: %f\n", status, 1e-5 * density * rgas * temperature * (1.0 + value));
    
    // 
    double pressure = pressure_bar(eos, temperature, density, molefracs, 2);
    printf("pressure: %f bar\n", pressure);

    double s_red = reduced_entropy(eos, temperature, density, molefracs, 2);
    printf("reduced_entropy: %f\n", s_red);

    get_Arxy(eos, 0, 0, temperature, density, molefracs, 2, &a00, errmsg, 100);
    get_Arxy(eos, 1, 0, temperature, density, molefracs, 2, &a10, errmsg, 100);
    printf("reduced_entropy: %f (a00: %f, a10: %f)\n", a10 - a00, a00, a10);

    printf("address: %ld \n", (int64_t)eos);
    double a_v_cast = da_dv((int64_t)eos, temperature, density, molefracs, 2);
    printf("dA_dv (cast): %f\n", a_v_cast);
    double a_v = da_dv(eos, temperature, density, molefracs, 2);
    printf("dA_dv: %f\n", a_v);

    // print error message if there is one
    if (status == -1) {
        printf("Error: %s\n", errmsg);
    }

    // free equation of state and exit
    free_eos_ptr(eos);
    return status;
}