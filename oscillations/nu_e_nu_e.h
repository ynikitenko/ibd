// Oscillations of \nu_e to \nu_e
#include <stdio.h>
#include <math.h>

double sin2_theta13 = 0.0214; // PDG16 best fit value for delta m2 > 0.
// double sin2_theta13 = 0.0218; // PDG16 best fit value for delta m2 < 0.
double sin2_theta13_3sigma_min = 0.0185; // PDG16, lower bound (3 sigma) for sin^2 \theta_{13}
double sin2_theta13_3sigma_max = 0.0248; // PDG16, upper bound (3 sigma) for sin^2 \theta_{13}
// 30/214. = 0.140, relative precision is about 14%.

double delta_m_atm_2 = 2.50E-3; // PDG16 best fit for m_1 < m_2 < m_3
// double delta_m_atm_2 = 2.46-3; // PDG16 best fit for m_3 < m_1 < m_2
double delta_m_atm_2_3sigma_min = 2.33E-3; // PDG16, lower bound (3 sigma) for delta m^2
double delta_m_atm_2_3sigma_max = 2.63E-3; // PDG16, lower bound (3 sigma) for delta m^2
// 0.2/2.33 = 0.0858, relative precision is about 9%.

double oscill_const = 1.267;    // constant from oscillation formula

double E_nu_min = 1.8; // minimum antineutrino energy in MeV
double E_nu_max = 9.25; // maximum antineutrino energy in MeV
double step = 0.01;

double p_surv_nue(double E_nu, double L) {
    // Probability of survival of electron (anti)neutrino of energy E_nu (MeV) at distance L (metres).
    double sin2_2_theta13 = 4. * sin2_theta13 * (1 - sin2_theta13); 
    double res = 1. - sin2_2_theta13 * pow(sin(oscill_const * delta_m_atm_2 * L / E_nu), 2);
    return res;
}

double *get_survival_prob_arr(double L, FILE *fout) {
    fprintf(fout, "# Survival probability for electron (anti)neutrino of the given energy\n");
    int nbins = (E_nu_max - E_nu_min) / step; 
    double *surv_prob = malloc((nbins + 1) * sizeof(double));
    int i = 0;
    for (double E_nu = E_nu_min; E_nu <= E_nu_max; E_nu += step) {
        double sp = p_surv_nue(E_nu, L);
        fprintf(fout, "%f %f\n", E_nu, sp);
        surv_prob[i] = sp;
        ++i; 
    }
    return surv_prob;
}

void generate_data(double L) { 
    double sin2_theta13_mean = 0.5 * (sin2_theta13_3sigma_min + sin2_theta13_3sigma_max);
    double sin2_theta13_sigma = (sin2_theta13_mean - sin2_theta13_3sigma_min) / 3.;
    double sin2_theta13_5sigma_min = sin2_theta13_mean - 5 * sin2_theta13_sigma;
    double sin2_theta13_5sigma_max = sin2_theta13_mean + 5 * sin2_theta13_sigma;

    double delta_m_atm_2_mean = 0.5 * (delta_m_atm_2_3sigma_min + delta_m_atm_2_3sigma_max); 
    double delta_m_atm_2_sigma = (delta_m_atm_2_3sigma_max - delta_m_atm_2_3sigma_min)/ 6.;
    double delta_m_atm_2_5sigma_min = delta_m_atm_2_mean - 5. * delta_m_atm_2_sigma;
    double delta_m_atm_2_5sigma_max = delta_m_atm_2_mean + 5. * delta_m_atm_2_sigma;

    FILE *fout;
    double *surv_prob;

    fout = fopen("nu_e_nu_e_sin2_theta13_dm2m0.txt", "w"); // sin^2 \theta_{13} for delta m^2 > 0.
    surv_prob = get_survival_prob_arr(L, fout);
    free(surv_prob);
    fclose(fout);

    const short NPLOTS = 4;
    double sin2_theta13s[] = { sin2_theta13_5sigma_min, sin2_theta13_3sigma_min, sin2_theta13_3sigma_max, sin2_theta13_5sigma_max};
    const char *fouts[] = { 
        "nu_e_nu_e_sin2_theta13_5sigma_min.txt",
        "nu_e_nu_e_sin2_theta13_3sigma_min.txt",
        "nu_e_nu_e_sin2_theta13_3sigma_max.txt",
        "nu_e_nu_e_sin2_theta13_5sigma_max.txt"
    };
    for(int i = 0; i < NPLOTS; ++i) {
        fout = fopen(fouts[i], "w");
        sin2_theta13 = sin2_theta13s[i];
        surv_prob = get_survival_prob_arr(L, fout);
        free(surv_prob);
        fclose(fout);
    }

    fout = fopen("nu_e_nu_e_sin2_theta13_5sigma_min_deltam2_5sigma_min.txt", "w");
    sin2_theta13 = sin2_theta13_5sigma_min;
    delta_m_atm_2 = delta_m_atm_2_5sigma_min;
    surv_prob = get_survival_prob_arr(L, fout);
    free(surv_prob);
    fclose(fout);

    fout = fopen("nu_e_nu_e_sin2_theta13_5sigma_max_deltam2_5sigma_max.txt", "w");
    sin2_theta13 = sin2_theta13_5sigma_max;
    delta_m_atm_2 = delta_m_atm_2_5sigma_max;
    surv_prob = get_survival_prob_arr(L, fout);
    free(surv_prob);
    fclose(fout);
    return;
}

