#ifndef ARBITRARY_GUARD_MACRO_CODE_OF_DOOM_WHICH_IS_CALLED_FCS_UTIL
#define ARBITRARY_GUARD_MACRO_CODE_OF_DOOM_WHICH_IS_CALLED_FCS_UTIL

#include<vector>
#include<exception>

namespace FCSTools
{

  // de facto byte-length in bits
  static const std::size_t DEFACTO_BYTEL = 8;

  // ERROR MESSAGES
  struct fcs_error : std::exception {
    virtual const char* what () const throw () {
      return "FCS Error";}
  };
  struct no_keyword : fcs_error {};
  struct no_keyword_fatal : no_keyword {};
  template <typename Base>
  struct T_no_nkeyword : Base {
    T_no_nkeyword (int n=0) : n_(n) {}
    virtual int N () const { return this->n_; }
    int n_;
  };
  typedef T_no_nkeyword<no_keyword> no_nkeyword;
  typedef T_no_nkeyword<no_keyword_fatal> no_nkeyword_fatal;
  struct invalid_datatype : fcs_error {
    virtual const char* what () const throw () {
      return "Invalid datatype"; }
  };
  
  struct unsupported_fcs_format : fcs_error {
    virtual const char* what () const throw () {
      return "Unsupported FCS format";}
  };
  struct text_miscalculation_error : fcs_error {
    virtual const char* what () const throw () {
      return "Position of text portion of FCS file unknown";}
  };
  struct data_miscalculation_error : fcs_error {
    virtual const char* what () const throw () {
      return "Position of data portion of FCS file unknown";}
  };

  struct invalid_bit_length : fcs_error {
    virtual const char* what () const throw () {
      return "Invalid bit length for parameter"; }
  };
  struct uncorrelated_mode : fcs_error {
    virtual const char* what () const throw () {
      return "Uncorrelated mode is unsupported"; }
  };
  struct correlated_mode : fcs_error {
    virtual const char* what () const throw () {
      return "Correlated mode is unsupported"; }
  };
  struct float_datatype : invalid_datatype {
    virtual const char* what () const throw () {
      return "Float is an invalid datatype";}
  };
  struct double_datatype : invalid_datatype {
    virtual const char* what () const throw () {
      return "Double is an invalid datatype";}
  };
  struct fixed_ascii_datatype : invalid_datatype {
    virtual const char* what () const throw () {
      return "Fixed ASCII is an invalid datatype";}
  };
  

  // COMPLIANT KEYWORDS
  struct no_total_keyword : no_keyword {
    virtual const char* what () const throw () {
      return "TOT"; }
  };
  struct no_nextdata_keyword : no_keyword {
    virtual const char* what () const throw () {
      return "NEXTDATA"; }
  };

  // FATAL KEYWORDS
  struct no_datatype_keyword : no_keyword_fatal {
    virtual const char* what () const throw () {
      return "DATATYPE"; }
  };
  struct no_byteorder_keyword : no_keyword_fatal {
    virtual const char* what () const throw () {
      return "BYTEORDER"; }
  };
  struct no_mode_keyword : no_keyword_fatal {
    virtual const char* what () const throw () {
      return "MODE"; }
  };
  struct no_par_keyword : no_keyword_fatal {
    virtual const char* what () const throw () {
      return "PAR"; }
  };
  struct no_bitsize_keyword : no_nkeyword_fatal {
    no_bitsize_keyword (int n=0) : no_nkeyword_fatal (n) {}
    virtual const char* what () const throw () {
      return "PNB"; }
  };
  struct no_range_keyword : no_nkeyword_fatal {
    no_range_keyword (int n=0) : no_nkeyword_fatal (n) {}
    virtual const char* what () const throw () {
      return "PNR"; }
  };

  // String2T conversion
  template <typename T>
  T convert (std::string const& str)
  {
    std::stringstream ss(str);
    T t;
    ss >> t;
    return t;
  }
  
  // Converts from one type of List mode data to
  // another type of List mode data
  // The cast is UNSAFE. No promises here.
  // We do not use the Boost casting mechanism
  // because Boost is too heavy-duty to add for
  // now.
  template <typename LValueType, typename RValueType>
  void convert (std::vector< std::vector<LValueType> >& Result,
		std::vector< std::vector<RValueType> > const& Data)
  {
    const std::size_t DSize = Data.size ();
    Result = std::vector< std::vector<LValueType> >(DSize);
    for (std::size_t i=0; i<DSize; ++i)
      for (std::size_t j=0; j<Data[i].size (); ++j)
    	Result[i].push_back ((LValueType)Data[i][j]);
  }

  // DEBUGGING UTILITY TO PRINT VECTORS

#ifdef DEBUG_FCS
  template <typename T>
  std::ostream& operator << (std::ostream& ostr, std::vector<T> const& V)
  {
    for (std::size_t i=0; i<V.size(); ++i)
      ostr << V[i] << " ";
    return ostr;
  }
#endif

  // ByteOrder stores the arbitrary byte-order
  // of the FCS file. 1234 is x86; 4321 is PPC
  // 2143 is Sparc, etc.
  typedef std::vector<std::size_t> byteorder_t;
  typedef std::vector<byteorder_t> byteorders_t;
  
  // shorthand
  typedef unsigned long blessed_integral;
  typedef long double blessed_float;

  typedef std::vector<blessed_integral> integral_vector_t;
  typedef std::vector<integral_vector_t> integral_vector_list_t;

  typedef std::vector<long double> ldouble_vector_t;
  typedef std::vector<ldouble_vector_t> ldouble_vector_list_t;

  // Human.Readable
  struct fcshr_error : fcs_error {
    virtual const char* what () const throw () {
      return "FCS:H.R error";}
  };
  struct improper_state_in_file : fcshr_error {
    virtual const char* what () const throw () {
      return "Improper state in file";}
  };
  
  typedef std::vector<double> element;
  typedef std::vector<element> dataset;
  
  std::string lower (std::string const& str)
  {
    std::string result (str);
    for (std::size_t i=0; i<result.size (); ++i)
      result[i] = std::tolower (result[i]);
    return result;
  }

}

#endif//ARBITRARY_GUARD_MACRO_CODE_OF_DOOM_WHICH_IS_CALLED_FCS_UTIL
