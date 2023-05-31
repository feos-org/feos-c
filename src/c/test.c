#include <stdio.h>
#include <stdint.h>

typedef struct eos eos_t;
extern eos_t *eos_from_json(const char *str);
extern void free_eos_ptr(eos_t *);
extern double pressure_bar(eos_t *eos_ptr, double temperature, double density, const double *molefracs, size_t len);
extern double da_dv(eos_t *eos_ptr, double temperature, double density, const double *molefracs, size_t len);

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
    double temperature = 300.0;
    double density = 10000.0;
    double molefracs[] = {0.1, 0.9};
    double pressure = pressure_bar(eos, temperature, density, molefracs, 2);
    printf("pressure: %f bar\n", pressure);
    printf("address: %ld \n", (int64_t)eos);
    double a_v = da_dv(eos, temperature, density, molefracs, 2);
    printf("dA_dv: %f\n", a_v);
    free_eos_ptr(eos);
}