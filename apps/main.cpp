#include "ccrh.h"
#include "io.h"

int main() {
    try {
        std::cout << "[" << current_time() << "] Starting program..." << std::endl;

        CCRH ccrh;
        ccrh.evolve();
        ccrh.dump();

        std::cout << "[" << current_time() << "] Program completed successfully." << std::endl;

        return EXIT_SUCCESS;
    } catch (const std::exception &e) {
        std::cerr << "[" << current_time() << "] Exception: " << e.what() << std::endl;

        return EXIT_FAILURE;
    }
}
