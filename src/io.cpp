#include "io.h"

// Function to open file and return output stream
std::ofstream open_output_file(const std::string& filename) {
    std::ofstream outfile(std::string("output/") + filename);
    if (!outfile.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        throw std::ios_base::failure("Failed to open output file.");
    }
    return outfile;
}

// Helper function to get the current time as a string
std::string current_time() {
    std::time_t now = std::time(nullptr);

    std::tm time_info;
    localtime_r(&now, &time_info);

    char buffer[26];  // Enough to hold "Day Mon dd hh:mm:ss yyyy\n\0"
    std::strftime(buffer, sizeof(buffer), "%a %b %d %H:%M:%S %Y", &time_info);

    return std::string(buffer);
}
