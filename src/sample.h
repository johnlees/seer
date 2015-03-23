/*
 * sample.h
 * Header file for sample class
 */

#include <iostream>
#include <string>
#include <stdexcept>

// Samples have a name and phenotype
// Covariates will be added later
class Sample
{
   public:
      // Initialisation
      Sample(int phenotype, std::string name);
      Sample();

      // nonmodifying operations
      int pheno() const { return phenotype; }
      std::string iid() const { return name; }
      static bool compareSamples(Sample lhs, Sample rhs) { return (lhs.name < rhs.name); }

   private:
      int phenotype;
      std::string name;
};

// Overload input operator
std::istream& operator>>(std::istream &is, Sample& s);
