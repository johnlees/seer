/*
 * sample.h
 * Header file for sample class
 */

#include <iostream>
#include <string>

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


   private:
      int phenotype;
      std::string name;
};

// Overload input operator
std::istream& operator>>(std::istream &is, Sample& s);
