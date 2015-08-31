/*
 * File: seerErr.cpp
 *
 * Throws error messages
 *
 */

#include "seercommon.hpp"

void badCommand(const std::string& command, const std::string& value)
{
   throw std::runtime_error("Bad " + command + " specified: " + value + "\n");
}

