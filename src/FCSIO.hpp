#ifndef ARBITRARY_GUARD_MACRO_CODE_OF_DOOM_WHICH_IS_CALLED_FCS_2_IO
#define ARBITRARY_GUARD_MACRO_CODE_OF_DOOM_WHICH_IS_CALLED_FCS_2_IO

#include<cmath>
#include<vector>
#include<map>
#include<string>
#include<iostream>
#include<sstream>
#include<exception>
#include<iomanip>
#include<FCSUtil.hpp>

namespace FCSTools
{

  // Some data in the FCS format is specified
  // as `N'; that is, it is defined for each parameter
  //
  // Name $PnN: name of the column (in List mode)
  // Specification $PnS: light-source
  // BitSize $PnB: number of bits (bytes) used
  //   to store this column's data
  // Range $PnR: specifies the data-range from [0,value(PnR)]
  // Scale $PnE: decade range of the data, 0,0 is linear
  // ByteOrder ---: local ordering of the bytes to ease
  //   reading the data
  struct NDatum
  {
    NDatum ()
      : Name(), Specification(), BitSize(), Range()
      , Scale(), ByteOrder() {}
    std::string Name;
    std::string Specification;
    std::size_t BitSize;
    std::size_t Range;
    std::pair<std::size_t,std::size_t> Scale;
    byteorder_t ByteOrder;
  };
  typedef std::vector<NDatum> NData;

  // Header is the Key/Value section of the FCS file
  //
  // File $FIL: original name of the file for this data
  // Mode $MODE: List, Uncorrelated, Correlated; only List
  //   is supported and used (de facto)
  // Datatype $DATATYPE: ASCII, Integer, Double, Float;
  //   we support (Variable) ASCII [when $PnB is a '*'],
  //   and Integer
  // Begin/EndTime $BTIM, $ETIM: start and stop of initial
  //   data acquisition
  // System $SYS: operating system the data was acquired
  //   under
  // Cytometer(?) $CYT: machine the data was acquired with
  // Date $DATE: day (unknown format) data was acquired on
  // Total $TOT: the total number of events in the file
  // Nextdata $NEXTDATA: offset to the next data set in the
  //   file
  // ByteOrder $BYTEORDER: permutation (ordering) of the
  //   bytes for Integer
  // Parameter: the per-parameter data
  // AllKeywords: optional/unknown/vender-specific keywords
  //   in addition to those above
  struct Header
  {
    typedef NData parameter_type;
    Header ()
      : File(), Mode(), Datatype(), BeginTime(), EndTime()
      , System(), Cytometer(), Date(), Total(), Nextdata()
      , ByteOrder(), Parameter(), AllKeywords() {}

    std::string File;
    std::string Mode;
    std::string Datatype;
    std::string BeginTime;
    std::string EndTime;
    std::string System;
    std::string Cytometer;
    std::string Date;
    std::size_t Total;
    std::size_t Nextdata;
    byteorder_t ByteOrder;
    NData Parameter;
    std::map<std::string, std::string> AllKeywords;
    bool has_keyword (std::string const& key) const {
      return this->AllKeywords.end() != this->AllKeywords.find(key);
    }
    std::string const& operator [] (std::string const& key) const {
      return this->AllKeywords.find(key)->second;
    }
    std::string& operator [] (std::string const& key) {
      return this->AllKeywords[key];
    }
  };

 
  // An FCS file is a Head and a Data section
  template <typename ValueType>
  struct FCS
  {
    typedef Header head_type;
    typedef std::vector<std::vector<ValueType> > data_type;

    FCS () : Head(), Data() {}

    head_type Head;
    data_type Data;
  };

  // Acquires data in List format in integral form
  // List data is a sequential set of vectors (positions in space)
  // representing elements. Integral data is stored in bytes
  // (although unspecified) with some arbitrary permutation
  integral_vector_list_t
  mode_L_integral (std::vector<unsigned char> const& Buffer,
       Header const& Head)
  {
    const std::size_t bordsz = Head.Parameter.size ();
    std::size_t veclength = 0;
    byteorders_t nByteOrders (bordsz);

    // This is for the algebraic permutation of the data
    // to get it to the native format.
    // The idea is that we grab a byte of data from
    // the stream and multiply it by the appropriate scale.
    // For instance, if the byte-order is: 1342 then we
    // scale the first byte by `1', the second byte by
    // 2^(DEFACTO_BYTEL*2) {e.g. 256*256}, the third byte by
    // 2^(DEFACTO_BYTEL*3) {e.g. 256^3}, and the fourth byte
    // by 2^DEFACTO_BYTEL {e.g. 256}
    // This computes those values.
    for (std::size_t i=0; i<bordsz; ++i)
      {
        for (std::size_t j=0; j<Head.Parameter[i].ByteOrder.size (); ++j)
          nByteOrders[i].push_back (1 << (Head.Parameter[i].ByteOrder[j]
              * DEFACTO_BYTEL));
        veclength += (Head.Parameter[i].BitSize / DEFACTO_BYTEL);
      }
    
    // The number of elements, this will agree with
    // $TOT if all goes well
    const std::size_t nElements = Buffer.size () / veclength;

    // Acquire the data; terminate when we have
    // all of the data
    integral_vector_list_t ivl (nElements);
    std::size_t element = 0;

    std::vector<unsigned char>::const_iterator cData = Buffer.begin ();

    while (element < nElements)
      {
        integral_vector_t iv (bordsz, 0);
        for (std::size_t i=0; i<bordsz; ++i)
          {
            blessed_integral sl = 0;
            // perform the algebraic permutation
            for (std::size_t j=0; j<nByteOrders[i].size (); ++j, ++cData)
              sl += nByteOrders[i][j] * *cData;
            iv[i] = sl;
          }
        ivl[element] = iv;
        ++element;
      }
    return ivl;
  }
  
  // ASCII data (variable length) in List mode is similar
  // to Integral List data, but the data is stored in a 
  // human-readable format
  integral_vector_list_t
  mode_L_ASCII_variable (std::vector<unsigned char> const& Buffer,
       Header const& Head)
  {
    const std::size_t Dimension = Head.Parameter.size ();
    std::string sData (Buffer.begin (), Buffer.end ());
    std::stringstream ss (sData);
    integral_vector_list_t ivl;
    integral_vector_t iv (Dimension, 0);

    // This is personal opinion: trying to guess
    // the reserve limit: sData.size () / (5*Dimension)
    // gives a rough-bound assuming (typical) 10-bit
    // precision data plus one space
    ivl.reserve (sData.size () / (5*Dimension));

    // Simplicity itself: pack the data into a stream
    // and read the numbers out one at a time; stuff
    // them into vectors and append to the list!
    bool loop = true;
    while (loop)
      {
        for (std::size_t i=0; i<Dimension; ++i)
          {
            ss >> iv[i];
            if (! ss)
              {
          loop = false;
          break;
              }
          }
        if (loop)
          ivl.push_back (iv);
      }
    return ivl;
  }
 
  // Writes an FCS file to a stream
  template <typename ValueType>
  bool Writer (std::ostream& FCSFile, FCS<ValueType> const& fcs)
  {
    // first, compute how long the data is going to be
    std::string VariableString;
    std::size_t DataLength; // <- place here

    // we use our own definition for BitSize
    std::vector<std::size_t> BitSize (fcs.Head.Parameter.size (),
              sizeof(blessed_integral) * DEFACTO_BYTEL);

    if ("*" == fcs.Head.Datatype)
      {
        std::stringstream ss;
        for (std::size_t i=0; i<fcs.Data.size (); ++i)
          for (std::size_t j=0; j<fcs.Data[i].size (); ++j)
            ss << fcs.Data[i][j] << " ";
        VariableString = ss.str ();
        DataLength = VariableString.size ();
      }
    else if ("I" == fcs.Head.Datatype)
      { // precompute the length
        // Now we figure out how to write the data in the proper
        // width; first we determine if all of the BitSizes are
        // char-aligned. Then we recreate the byte-length
        // for each one.
        for (std::size_t i=0; i<fcs.Head.Parameter.size (); ++i)
          if ((fcs.Head.Parameter[i].BitSize%DEFACTO_BYTEL) != 0)
            { // smallest width > than BitSize that -is- char-aligned
              std::size_t floor = fcs.Head.Parameter[i].BitSize / DEFACTO_BYTEL;
              BitSize[i] = (floor + 1) * DEFACTO_BYTEL; // i.e. ceil
            }
          else
            BitSize[i] = fcs.Head.Parameter[i].BitSize;
        std::size_t DataWidth = 0;
        for (std::size_t i=0; i<BitSize.size (); ++i)
          DataWidth += BitSize[i] / DEFACTO_BYTEL;
      
        DataLength = DataWidth * fcs.Data.size ();
      }
    else if ("A" == fcs.Head.Datatype)
      throw fixed_ascii_datatype ();
    else if ("F" == fcs.Head.Datatype)
      throw float_datatype ();
    else if ("D" == fcs.Head.Datatype)
      throw double_datatype ();
    else
      throw invalid_datatype ();

    // buffer for the output data
    std::vector<unsigned char> Data (DataLength);

    // Pack the data in
    if ("*" == fcs.Head.Datatype)
      for (std::size_t i=0; i<DataLength; ++i)
        Data[i] = VariableString[i];
    else if ("I" == fcs.Head.Datatype)
      {

        std::size_t Cursor = 0;
        std::vector<std::size_t> ByteLength (fcs.Head.Parameter.size (),
                     sizeof(blessed_integral));
        for (std::size_t i=0; i<fcs.Head.Parameter.size (); ++i)
          ByteLength[i] = BitSize[i] / DEFACTO_BYTEL;
      
        for (std::size_t i=0; i<fcs.Data.size (); ++i)
          for (std::size_t j=0; j<fcs.Data[i].size (); ++j)
            {
              unsigned int val = fcs.Data[i][j];
              // we're not going to fuck with clever byte-reording
              // we're going to do this algebraically
              // it's slower, but guaranteed to be platform independent
              for (std::size_t k=0; k<ByteLength[j]; ++k, ++Cursor)
                {
                  const std::size_t Modulus = (1UL<<(DEFACTO_BYTEL * (k+1)));
                  const std::size_t Divisor = (1UL<<(DEFACTO_BYTEL * k));
                  Data[Cursor] = (val % Modulus) / Divisor;
                }
            }
      }

    std::stringstream ssKeywords;
    ssKeywords << "/$TOT/" << fcs.Data.size ()
         << "/$PAR/" << fcs.Head.Parameter.size ()
         << "/$NEXTDATA/0"
         << "/$DATATYPE/" << fcs.Head.Datatype
         << "/$MODE/" << fcs.Head.Mode
         << "/$BYTEORD";
    for (std::size_t i=0; i<sizeof(blessed_integral); ++i)
      ssKeywords << (0==i?"/":",") << (i+1);

    const std::size_t MaximumDataLengthForFCS2_0 = 99999999;
    const bool Needs3_0 = DataLength > MaximumDataLengthForFCS2_0;
    const std::size_t OffsetToText = 256;
    const std::size_t DataPadding = 25;

    std::size_t sBitSize;
    if ("I" == fcs.Head.Datatype)
      sBitSize = sizeof(blessed_integral) * DEFACTO_BYTEL;
    for (std::size_t i=0; i<fcs.Head.Parameter.size (); ++i)
      {
        std::size_t N = i+1;
        ssKeywords << "/$P" << N << "B/";
        if (sBitSize > 0) ssKeywords << BitSize[i];
        else ssKeywords << "*";
        ssKeywords << "/$P" << N << "R/" << fcs.Head.Parameter[i].Range;
        if (fcs.Head.Parameter[i].Name.size () > 0)
          ssKeywords << "/$P" << N << "N/" << fcs.Head.Parameter[i].Name;
        ssKeywords << "/$P" << N << "E/"
             << fcs.Head.Parameter[i].Scale.first << ","
             << fcs.Head.Parameter[i].Scale.second;
      }
    // So here we have a conundrum: we must ascertain the length
    // of the keyword section to determine where the DataSection begins,
    // but if the DataSection is too long we have to insert more
    // keywords. The trick is to know this:
    // /$BEGINDATA/-100-chars-/$ENDDATA/-100-chars/
    // 123456789012           3456789012          3
    // 000000000111           1111111222          2
    // 23 + 200 = 223: KeyWords.size () + OffsetToText + 223 + DataPadding = DataOffset!
    std::size_t DataOffset = 
      223 + ssKeywords.str ().size ()
      + OffsetToText + DataPadding;
    std::size_t DataEnd = DataOffset + DataLength;
    ssKeywords << "/$BEGINDATA/" << std::setw (100) << DataOffset
         << "/$ENDDATA/" << std::setw (100) << DataEnd;
    ssKeywords << "/";
    std::string Keywords = ssKeywords.str ();


    if (Needs3_0)
      FCSFile << "FCS3.0    ";
    else
      FCSFile << "FCS2.0    ";
    FCSFile.fill (' ');
    FCSFile << std::setw (8) << OffsetToText
      << std::setw (8) << (OffsetToText+Keywords.size()-1);
    if (Needs3_0)
      FCSFile << std::setw (8) << 0 << std::setw (8) << 0;
    else
      FCSFile << std::setw (8) << (OffsetToText+Keywords.size()+25)
        << std::setw (8) << (OffsetToText+Keywords.size()+25+DataLength);
    FCSFile << std::setw (8) << 0
      << std::setw (8) << 0
      << std::setw (OffsetToText-58) << " "
      << Keywords
      << std::setw (DataPadding) << " ";
    // write the data in 4096-element chunks
    const std::size_t PageSize = 4096;
    std::size_t Chunk = 0;
    char Buffer[PageSize];
    while (Chunk < Data.size ())
      {
        std::size_t GetC = PageSize;
        if ((Chunk+GetC) > Data.size ())
          GetC = Data.size () - Chunk;
        for (std::size_t i=0; i<GetC; ++i)
          Buffer[i] = Data[Chunk+i];
        FCSFile.write (Buffer, GetC);
        Chunk += PageSize;
      }
    return FCSFile;
  }

  // Attempts to read FCS files with the following caveats:
  // 1) there is, deliberately, no attempt to read fixed ASCII,
  //    correlated, or uncorrelated data
  // 2) there is no interpreter for Float or Double data
  // 3) is quite liberal with the set of `necessary' keywords
  // 4) if there is a problem, it throws an exception; I know
  //    some people really hate exceptions and call a pox on
  //    my head for using them. However, they -are- standard,
  //    they -are- used throughout iostreams, and they -are-
  //    understandable. Sorry.
  template <typename ValueType>
  FCS<ValueType>
  Reader (std::istream& FCSFile, bool Compliance_P = false)
  {
    typedef std::vector<ValueType> ElementT;
    typedef std::vector<ElementT> VElementT;

    const std::size_t InitialOffset = FCSFile.tellg ();

    FCS<ValueType> fcs;

    // This is the `header' section as defined in the FCS2.0 standard;
    // There are some sort of noises about having `8' characters for this,
    // `10' characters for that, etc. etc.
    // I assume the data is white-space delimited. This is non-standard
    // but if you *do* have a crap-load of data, FCS2.0 is not the standard
    // for you.

    std::string FCSKind;
    std::size_t KeysBegin = 0, KeysEnd = 0,
                KeySection = 0, DataSection = 0,
                DataBegin = 0, DataEnd = 0;
    FCSFile >> FCSKind;

    double Kind = 0.0;

    if ("FCS2.0" == FCSKind)
      Kind = 2.0;
    else if ("FCS3.0" == FCSKind)
      Kind = 3.0;

    // Want any teeth pulled?
    const std::size_t KindLength = 10;
    const std::size_t LocationSize = 8;
    const std::size_t NumberLocations = 4;
    const std::size_t HeaderLength = KindLength + LocationSize * NumberLocations;
    char DumbBuffer[HeaderLength+10];
    FCSFile.seekg (InitialOffset);
    FCSFile.read (DumbBuffer, HeaderLength);
    DumbBuffer[HeaderLength] = 0;

    FCSKind = std::string (DumbBuffer, DumbBuffer+KindLength);
    std::stringstream ssLocations;
    for (std::size_t i=0; i<NumberLocations; ++i)
      {
        std::string Loc (DumbBuffer + KindLength + i * LocationSize,
             DumbBuffer + KindLength + (i+1) * LocationSize);
        ssLocations << Loc << " ";
      }
    ssLocations >> KeysBegin >> KeysEnd >> DataBegin >> DataEnd;

    if (2.0 != Kind && 3.0 != Kind)
      throw unsupported_fcs_format ();

    bool KeyWordDataSection = false;

    if (KeysEnd < KeysBegin)
      throw text_miscalculation_error ();
    if (DataEnd < DataBegin)
      if (DataEnd == DataBegin && 0 == DataBegin)
        {
          if ((Compliance_P && 3.0 == Kind) || !Compliance_P)
            KeyWordDataSection = true;
          else
            throw data_miscalculation_error ();
        }

    FCSFile.seekg (InitialOffset+KeysBegin);
    
    // Now we get the key/value pairs; they
    // are delimited by the first character in the
    // `Text' section. I don't know why.
    char Delimiter = (char)FCSFile.peek ();
    KeySection = KeysEnd - KeysBegin + 1;
    char cKeyWordBuffer[KeySection+2];

    FCSFile.get (cKeyWordBuffer, KeySection+1, 0);
    std::string KeyWordBuffer (cKeyWordBuffer+1, KeySection+1);

    std::stringstream ss (KeyWordBuffer);
    do
      {
        std::string key, value;
        std::getline (ss, key, Delimiter);
        std::getline (ss, value, Delimiter);
        fcs.Head[key] = value;
      }
    while (ss);
    
    if (KeyWordDataSection)
      {
        if (fcs.Head.has_keyword ("$BEGINDATA"))
          DataBegin = convert<std::size_t>(fcs.Head["$BEGINDATA"]);
        else
          throw data_miscalculation_error ();
        if (fcs.Head.has_keyword ("$ENDDATA"))
          DataEnd = convert<std::size_t>(fcs.Head["$ENDDATA"]);
        else
          throw data_miscalculation_error ();  
      }

    // Go to the Data section; we will read the data in
    // as 4096-byte-pages
    FCSFile.seekg (InitialOffset+DataBegin);
    DataSection = DataEnd - DataBegin + 1;
    std::vector<unsigned char> DataBuffer;
    if (true)
      {
        const std::size_t PageSize = 128000;
        DataBuffer.reserve(DataSection);
        std::size_t Chunk = 0;
        char Buffer[PageSize];
        while (Chunk < DataSection)
          {
            std::size_t GetC = PageSize;
            if ((Chunk+GetC) > DataSection)
              GetC = DataSection - Chunk;
            FCSFile.read (Buffer, GetC);
            DataBuffer.insert (DataBuffer.end(), Buffer, Buffer+GetC);
            Chunk += PageSize;
          }
      }
    else
      {
        char *mBuffer = new char[DataSection];
        FCSFile.read (mBuffer, DataSection);
        DataBuffer.insert (DataBuffer.end(), mBuffer, mBuffer+DataSection);
      }

    // Pick through and get a bunch of keywords real-quick
    // Keywords that "throw" are "fatal": they are -required-
    // to interpret the FCS file.
    if (fcs.Head.has_keyword ("$PAR"))
      fcs.Head.Parameter = NData(convert<std::size_t>(fcs.Head["$PAR"]));
    else
      throw no_par_keyword ();

    const std::size_t nParameters = fcs.Head.Parameter.size ();

    // some keywords throw if they are 'required' by the standard
    // but are not required to actually read the file; they
    // throw when the Compliance_P switch is set
    if (fcs.Head.has_keyword ("$NEXTDATA"))
      fcs.Head.Nextdata = convert<std::size_t>(fcs.Head["$NEXTDATA"]);
    else if (Compliance_P)
      throw no_nextdata_keyword ();

    if (fcs.Head.has_keyword ("$SYSTEM"))
      fcs.Head.System = fcs.Head["$SYS"];

    if (fcs.Head.has_keyword ("$BTIM"))
      fcs.Head.BeginTime = fcs.Head["$BTIM"];

    if (fcs.Head.has_keyword ("$ETIM"))
      fcs.Head.EndTime = fcs.Head["$ETIM"];

    if (fcs.Head.has_keyword ("$FIL"))
      fcs.Head.File = fcs.Head["$FIL"];

    if (fcs.Head.has_keyword ("$DATE"))
      fcs.Head.Date = fcs.Head["$DATE"];

    bool Variable = false;

    for (std::size_t npar = 0; npar < nParameters; ++npar)
      {
        std::stringstream par;
        par << "$P" << (npar+1) << "B";
        std::string key;
        par >> key;
        if (fcs.Head.has_keyword (key))
          {
            // Bit-Size OR Byte-Size, depending upon
            // the data type; since we ignore fixed ASCII,
            // it is -always- Bit-Size
            if ("*" == fcs.Head[key])
              {
                Variable = true;
                fcs.Head.Parameter[npar].BitSize = 0; // variable bit-size
              }
            else
              fcs.Head.Parameter[npar].BitSize = convert<std::size_t>(fcs.Head[key]);
          }
        else
          throw no_bitsize_keyword (npar);
        // must be byte-aligned to a char
        if ((fcs.Head.Parameter[npar].BitSize%DEFACTO_BYTEL) != 0)
          throw invalid_bit_length ();
        par.clear ();
        par.str ("");
        
        par << "$P" << (npar+1) << "R";
        par >> key;
        if (fcs.Head.has_keyword (key))
          fcs.Head.Parameter[npar].Range = convert<std::size_t>(fcs.Head[key]);
        else if (Compliance_P)
          // if we ever support Correlated or Uncorrelated data
          // we will need this keyword
          throw no_range_keyword (npar);
        par.clear ();
        par.str ("");
        
        // an `optional' keyword; morons
        par << "$P" << (npar+1) << "N";
        par >> key;
        if (fcs.Head.has_keyword (key))
          fcs.Head.Parameter[npar].Name = fcs.Head[key];
        par.clear ();
        par.str ("");
        
        // I have no clue what this is; it is, however,
        // ubiquitous; a number of readers fail if it is
        // missing
        par << "$P" << (npar+1) << "E";
        par >> key;
        if (fcs.Head.has_keyword (key))
          {
            std::stringstream exp (fcs.Head[key]);
            std::size_t L, R;
            std::string sL, sR;
            std::getline (exp, sL, ',');
            std::getline (exp, sR);
            L = convert<std::size_t>(sL);
            R = convert<std::size_t>(sR);
            fcs.Head.Parameter[npar].Scale = std::make_pair (L, R);
          }
        par.clear ();
        par.str ("");
        
        // Yeah. No clue.
        par << "$P" << (npar+1) << "S";
        par >> key;
        if (fcs.Head.has_keyword (key))
          fcs.Head.Parameter[npar].Specification = fcs.Head[key];
        par.clear ();
        par.str ("");
      }

    // Get the permutation of the data. The permutation
    // is delimted by a ','. I suppose, in theory, we
    // could have data that uses enough bytes to require
    // two-digits, which could be troublesome. Ugh. Now that
    // I think about it, this goddamn file format may last
    // that long.
    std::stringstream byteorder;
    if (fcs.Head.has_keyword ("$BYTEORD"))
      byteorder << fcs.Head["$BYTEORD"];
    else
      throw no_byteorder_keyword ();

    while (true)
      {
        std::string byteo;
        std::getline (byteorder, byteo, ',');
        if (! byteorder)
          break;
        fcs.Head.ByteOrder.push_back (convert<std::size_t>(byteo)-1);
      }
   
    // This creates a custom permutation for each parameter;
    // this just trivializes some later computations. What it
    // does is throw away permutation values that are "too big"
    // for our particular value, making look-up trivial.
    for (std::size_t i=0; i<nParameters; ++i)
      {
        for (std::size_t j=0; j<fcs.Head.ByteOrder.size (); ++j)
          if (fcs.Head.ByteOrder[j] < (fcs.Head.Parameter[i].BitSize/DEFACTO_BYTEL))
            fcs.Head.Parameter[i].ByteOrder.push_back (fcs.Head.ByteOrder[j]);
      }

    // Why? is this some lame attempt at redundancy checking?
    // Idiocy.
    if (fcs.Head.has_keyword ("$TOT"))
      fcs.Head.Total = convert<std::size_t>(fcs.Head["$TOT"]);
    else if (Compliance_P)
      throw no_total_keyword ();

    // legal: [L]ist, [U]ncorrelated, [C]orrelated
    if (fcs.Head.has_keyword ("$MODE"))
      fcs.Head.Mode = fcs.Head["$MODE"];
    else
      throw no_mode_keyword ();
    if ("U" == fcs.Head.Mode)
      throw uncorrelated_mode ();
    if ("C" == fcs.Head.Mode)
      throw correlated_mode ();

    // legal: [A]SCII, [I]nteger, [D]ouble, [F]loat
    if (fcs.Head.has_keyword ("$DATATYPE"))
      fcs.Head.Datatype = fcs.Head["$DATATYPE"];
    else
      throw no_datatype_keyword ();
    if ("A" == fcs.Head.Datatype && Variable)
      fcs.Head.Datatype = "*";
    if ("A" == fcs.Head.Datatype)
      throw fixed_ascii_datatype ();
    if ("F" == fcs.Head.Datatype)
      throw float_datatype ();
    if ("D" == fcs.Head.Datatype)
      throw double_datatype ();
    
    if (fcs.Head.has_keyword ("$CYT"))
      fcs.Head.Cytometer = fcs.Head["$CYT"];

    if ("L" == fcs.Head.Mode)
      { // list of vectors
        // we grab up to, but no more than Total
        integral_vector_list_t ivl;
        switch (fcs.Head.Datatype[0])
          {
          case 'I': ivl
              = mode_L_integral (DataBuffer, fcs.Head);
            break;
          case '*': ivl
              = mode_L_ASCII_variable (DataBuffer, fcs.Head);
            break;
          default : break;
          }
        convert (fcs.Data, ivl);
      }

    return fcs;
  }

}

#endif//ARBITRARY_GUARD_MACRO_CODE_OF_DOOM_WHICH_IS_CALLED_FCS_2_IO
