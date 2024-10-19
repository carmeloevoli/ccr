#ifndef CCRH_IO_H_
#define CCRH_IO_H_

#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>

std::string current_time();

std::ofstream open_output_file(const std::string& filename);

#endif