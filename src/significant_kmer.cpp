/*
 * File: significant_kmer.cpp
 *
 * Helper functions for the significant kmer class
 * Reads input into class
 *
 */

#include "significant_kmer.hpp"

Significant_kmer::Significant_kmer()
   :_line_nr(0), _num_covars(default_covars)
{
}

Significant_kmer::Significant_kmer(const int num_covars)
   :_line_nr(0),_num_covars(num_covars)
{
}

Significant_kmer::Significant_kmer(const std::string& word, const std::vector<std::string>& samples, const double maf, const double unadj_p, const double adj_p, const double lrt_p, const double beta, const double se, const std::string& comments)
   :_line_nr(0), _word(word), _samples(samples), _maf(maf), _unadj_p(unadj_p), _adj_p(adj_p), _adj_lrt_p(lrt_p), _beta(beta), _se(se), _comment(comments), _num_covars(default_covars)
{
}

// Fills significant kmer object from seer output file
// Sample vector is returned sorted
std::istream& operator>>(std::istream &is, Significant_kmer& sk)
{
   double maf, unadj_p, adj_p, lrt_p, beta, se;
   std::string sequence, sample, line_in, comments = "";
   std::vector<std::string> sample_list;

   // Read the line convert to stringstream, extract sequence and stats fields
   std::getline(is, line_in);
   std::stringstream line_stream(line_in);

   line_stream >> sequence;

   line_stream >> maf;
   line_stream >> unadj_p;
   line_stream >> adj_p;
   line_stream >> lrt_p;
   line_stream >> beta;
   line_stream >> se;

   // Ignore the covariate fields
   for (unsigned int i = 0; i < sk.num_covars(); ++i)
   {
      line_stream.ignore(std::numeric_limits<std::streamsize>::max(), '\t');
   }

   line_stream >> comments;

   // Read remainder of line, which is sample names, in the same way they were
   // written
   std::copy(std::istream_iterator<std::string>(line_stream), std::istream_iterator<std::string>(), std::back_inserter(sample_list));

   // Ensure vector remains sorted on sample name
   std::sort(sample_list.begin(), sample_list.end());
   sk = Significant_kmer(sequence, sample_list, maf, unadj_p, adj_p, lrt_p, beta, se, comments);

   return is;
}

// Print fields tab sep, identical to input. Doesn't print newline
std::ostream& operator<<(std::ostream &os, const Significant_kmer& k)
{
   os << std::fixed << std::setprecision(3) << k.sequence() << "\t" << k.maf()
      << "\t" << std::scientific << k.unadj() << "\t" << k.p_val() << "\t" << k.lrt_p_val()
      << "\t" << k.beta() << "\t" << k.se();

   std::vector<std::string> samples_found = k.samples_found();
   if (samples_found.size() > 0)
   {
      os << "\t";
      std::copy(samples_found.begin(), samples_found.end() - 1, std::ostream_iterator<std::string>(os, "\t"));
      os << samples_found.back();
   }

   return os;
}

unsigned int Significant_kmer::num_covars() const
{
   unsigned int covars = 0;
   if (_covar_p.size() > 0)
   {
      covars =_covar_p.size();
   }
   else
   {
      covars =_num_covars;
   }

   return covars;
}

// Returns the reverse complement of a k-mer
// Thanks to http://stackoverflow.com/questions/33074574/creating-complement-of-dna-sequence-and-reversing-it-c
// for good c++11 style
std::string Significant_kmer::rev_comp() const
{
   std::string rev_seq(_word);
   auto lambda = [](const char c)
   {
      switch (c)
      {
        case 'A':
            return 'T';
        case 'G':
            return 'C';
        case 'C':
            return 'G';
        case 'T':
            return 'A';
        default:
            throw std::runtime_error("Invalid nucleotide");
        }
    };

    std::transform(rev_seq.cbegin(), rev_seq.cend(), rev_seq.begin(), lambda);
    std::reverse(rev_seq.begin(), rev_seq.end());
    return rev_seq;
}


// Returns number of excess columns (i.e. number of covariate fields)
int parseHeader(const std::string& header_line)
{
   std::istringstream iss(header_line);

   int num_cols = 0;
   do
   {
      std::string col;
      iss >> col;
      if (col.length() > 0)
      {
         num_cols++;
      }
   }
   while(iss);

   return num_cols - standard_cols;
}

sortSigKmer::sortSigKmer()
   :_sort_field(0)
{
}

sortSigKmer::sortSigKmer(const std::string& sort_field)
   :_sort_field(0)
{
   this->addSortField(sort_field);
}


// Sets the sort field for sig kmers
void sortSigKmer::addSortField(const std::string& sort_field)
{
   if (sort_field == "maf")
   {
      this->_sort_field = 1;
   }
   else if (sort_field == "chisq")
   {
      this->_sort_field = 2;
   }
   else if (sort_field == "pval")
   {
      this->_sort_field = 3;
   }
   else if (sort_field == "beta")
   {
      this->_sort_field = 4;
   }
   else if (sort_field == "sequence")
   {
      this->_sort_field = 5;
   }
   else
   {
      this->_sort_field = 3;
   }
}

bool sortSigKmer::operator() (const Significant_kmer& sk1, const Significant_kmer& sk2) const
{
   bool compare;
   switch (_sort_field)
   {
      case 1:
         compare = sk1.maf() < sk2.maf();
         break;
      case 2:
         compare = sk1.unadj() < sk2.unadj();
         break;
      case 3:
         compare = sk1.p_val() < sk2.p_val();
         break;
      case 4:
         compare = sk1.beta() < sk2.beta();
         break;
      case 5:
         compare = sk1.sequence().length() < sk2.sequence().length();
         break;
      default:
         compare = sk1.p_val() < sk2.p_val();
   }
   return compare;
}

