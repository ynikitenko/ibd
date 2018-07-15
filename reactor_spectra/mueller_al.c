#include "mueller_al.h"

int main() {
    spectrum_vector *spec = get_spectrum_vector();
    print_spectrum("mueller_al.txt", "# Reactor spectrum for antineutrino energy", spec);
    delete spec;
    return 0;
}

