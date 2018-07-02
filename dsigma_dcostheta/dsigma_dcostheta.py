def generate_dsigma_dcostheta():
    import numpy as np
    nends = 201
    nbins = nends - 1
    costhetas = np.linspace(-1, 1, nends)

    def write_dsigma_Enu(Enu):
        def write_outfile(outfilename, data):
            with open(outfilename, "w") as outfile:
                for ct, ds in data:
                    line = "{:1.2f} {}\n".format(ct, ds)
                    outfile.write(line)

        outfilename = "dsigma_dcostheta_Enu_{}MeV.txt".format(Enu)
        dsigmas = map(lambda costheta: dsigma_dcostheta(costheta, Enu), costhetas)
        write_outfile(outfilename, zip(costhetas, dsigmas))

        norm = sum(dsigmas) / nends
        dsigmas_norm = map(lambda ct: ct / norm / 2, dsigmas)
        outfilename_norm = "dsigma_dcostheta_norm_Enu_{}MeV.txt".format(Enu)
        write_outfile(outfilename_norm, zip(costhetas, dsigmas_norm))

    Enus = [2, 4]
    for Enu in Enus:
        write_dsigma_Enu(Enu)
    return 

def dsigma_dcostheta(costheta, Enu):
    from ctypes import CDLL, c_double
    libibd = CDLL("libibd.so.1")

    # in fact, that is not needed, since Enu is given in the parameters
    # to the function. This is not a very good design though, 
    # one should get rid of the global variable if that's not needed.
    # Enu_lib = c_bool.in_dll(libibd, "neutrinoE")
    # Enu_lib.value = c_double(Enu)

    dsigma_dcostheta = libibd.find_DCS_free_cosTheta_RC
    dsigma_dcostheta.restype = c_double
    res = dsigma_dcostheta(c_double(Enu), c_double(costheta))
    return res

if __name__ == "__main__":
    generate_dsigma_dcostheta()

