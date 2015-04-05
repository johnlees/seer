/*
 * File: pangwasErr.cpp
 *
 * Throws error messages
 *
 */

#include "pangwas.hpp"

void badCommand(const std::string& command, const std::string& value)
{
   throw std::runtime_error("Bad " + command + " specified: " + value + "\n");
}

