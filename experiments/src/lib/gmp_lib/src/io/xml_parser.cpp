#include <gmp_lib/io/xml_parser.h>
#include <cctype>

namespace as64_
{

namespace gmp_
{

bool isNumber(const std::string &str)
{
  for (int i=0;i<str.size();i++)
  {
    if (!isdigit(str[i])) return false;
  }
  return true;
}

std::string errMsgAtLine(const std::string &fname, int line_number, const std::string &msg)
{
  std::ostringstream oss;

  oss << "[XmlParser Error]: Error parsing \"" << fname << "\", line " << line_number << ": " << msg;

  return oss.str();
}

std::string errMsg(const std::string &fname, const std::string &msg)
{
  std::ostringstream oss;

  oss << "[XmlParser Error]: Error parsing \"" << fname << ": " << msg;

  return oss.str();
}

bool XmlParser::valid_key(const std::string& key)
{
    if(par_map.find(key) == par_map.end() )
    {
        std::cerr << "\033[33mXmlParser::valid_key: Parameter " + key + " not found!\033[00m" << std::endl;
        return false;
    }
    return true;
}

std::complex<double> XmlParser::parse_cx(std::string str)
{
    double re,im;
    char i_ch;
    std::stringstream iss(str);

    // Parse full ...
    if(iss >> re >> im >> i_ch && (i_ch=='i'|| i_ch=='j')) return std::complex<double>(re,im);

    // ... or only imag
    iss.clear();
    iss.seekg(0,iss.beg);
    if(iss >> im >> i_ch && (i_ch=='i'|| i_ch=='j')) return std::complex<double>(0.0,im);

    // .. or only real
    iss.clear();
    iss.seekg(0,iss.beg);
    if(iss >> re) return std::complex<double>(re,0.0);

    // ... otherwise
    throw std::ios_base::failure("XmlParser::parse_cx: Could not parse complex number!");
}

XmlParser::XmlParser(const std::string& fname)
{
  this->fname = fname;
  std::ifstream fh;
  std::string line;
  size_t mark = 0;
  line_number = 0;
  column_number = 0;

  // Clear parameter map
  par_map.clear();

  // Open file
  fh.open(fname.c_str());
  if (!fh)
  {
    throw std::ios_base::failure(errMsg(fname, "Could not find " + fname));
  }
  else
  {
    // Parse
    while (std::getline(fh, line))
    {
      column_number = 1;
      line_number++;

      std::string keyS="";
      std::string dataS="";

      // Skip empty lines
      if (line.empty())
        continue;

      // Skip lines with only whitespace
      if(line.find_first_not_of("\t ")==std::string::npos)
        continue;

      // Remove comment
      mark = line.find("#");
      if(mark!=std::string::npos)
      {
        line.erase(mark,line.length());
        if (line.empty()) continue;
      }

      // Do we have a '=' or a ':'
      mark = line.find("=");
      if(mark==std::string::npos) mark = line.find(":");
      if(mark!=std::string::npos)
      {
        // Find key
        size_t i = line.find_first_not_of("\t ");
        keyS = line.substr(i,mark-i);
        keyS = keyS.substr(0,keyS.find_last_not_of("\t ")+1);
        column_number += i;

        if (keyS.empty()) throw std::ios_base::failure(errMsgAtLine(fname, line_number,"Missing key."));
        if(keyS.find_first_of("\t ")!=std::string::npos) throw std::ios_base::failure(errMsgAtLine(fname, line_number,"Invalide key: \033[31m" +  keyS + "\033[00m. Contains whitespace characters."));

        // Find data
        dataS = line.substr(mark+1,line.length());
        i = dataS.find_first_not_of("\t ");
        if (i==std::string::npos) throw std::ios_base::failure(errMsgAtLine(fname, line_number,"Missing value."));
        size_t i2 = dataS.find_last_not_of("\t ");
        dataS = dataS.substr(i, i2-i+1);

        // Do we have a vector/matrix
        mark = 0;
        if (dataS[0] == '[')
        {
          if (dataS.back() != ']') throw std::ios_base::failure(errMsgAtLine(fname, line_number,"Missing \"]\"."));
          dataS = dataS.substr(1, dataS.length()-2);

          if (dataS.find("\"") != std::string::npos) // array with strings
          {
            std::vector<std::string> temp_value;
            getStringsFromStringVector(dataS, temp_value);
          }
        }
        else if (dataS.back() == ']')
        {
          throw std::ios_base::failure(errMsgAtLine(fname, line_number,"Missing \"[\"."));
        }
        // Do we have a string
        else if (dataS[0] == '\"')
        {
          if (dataS.back() != '\"') throw std::ios_base::failure(errMsgAtLine(fname, line_number,"Missing \" at the end."));
          dataS = dataS.substr(1, dataS.length()-2);
        }
        else if (dataS.back() == '\"')
        {
          throw std::ios_base::failure(errMsgAtLine(fname, line_number,"Missing \" at the beginning."));
        }
        else
        {
          std::istringstream in_temp(dataS);
          std::string temp;
          in_temp >> temp;
          if (!temp.compare("true")) dataS="1";
          else if (!temp.compare("false")) dataS="0";
          // else throw std::ios_base::failure(errMsgAtLine(fname, line_number, "Unrecognized sequence: " + temp + ". If it's a string put it in quotes \"\""));
        }

        // std::cout << "dataS = " << dataS << "\n";

        // Insert to map
        par_map.insert(std::pair<std::string, std::string>(keyS, dataS));
      }
    }

    // Close file
    fh.close();
  }
}

void XmlParser::validateArrayValues(std::string &dataS)
{
  // check if the array contrains strings
  if (dataS.find("\"") == std::string::npos)
  {
    // no strings, so we expect numeric data
    std::string full_str, row, col;
    std::istringstream iss_full;

    full_str = dataS;
    int R = static_cast<int>(std::count(full_str.begin(),full_str.end(),';')+1);

    iss_full.str(full_str);
    std::getline(iss_full, row, ';');
    int C = static_cast<int>(std::count(row.begin(),row.end(),',')+1);

    // arma::Mat<T> x(R,C);
    //
    // iss_full.seekg(0,iss_full.beg);
    // for(int r=0; r<R; r++)
    // {
    //     std::getline(iss_full, row, ';');
    //     std::istringstream iss_row(row);
    //     for(int c=0; c<C; c++)
    //     {
    //         std::getline(iss_row, col, ',');
    //         std::istringstream iss_col(col);
    //         iss_col >> x(r,c);
    //     }
    // }
    // value = x;

    return;
  }

  // contrains strings!
  std::string str;

  std::cout << "dataS = " << dataS << "\n";

  int i1 = 0; // index of first \" of the string
  while (i1 < dataS.size())
  {
    int i2 = -1; // index of second \" of the string

    if (dataS[i1] != '\"')
    {
      i1++;
      continue;
    }

    for (int i=i1+1;i<dataS.size();i++)
    {
      if (dataS[i] == '\"')
      {
        i2 = i;
        break;
      }
    }

    if (i2 < 0) throw std::ios_base::failure(errMsgAtLine(fname, line_number,"Missing \"."));

    i1++;
    i2--;

    if (i1==i2) str += " " ;
    else str += dataS.substr(i1,i2-i1+1) + " ";

    std::cout << "dataS.substr(i1,i2-i1+1) = " << dataS.substr(i1,i2-i1+1) << "\n";

    i1 = i2+2;
  }
  dataS = "";
  for (int i=0;i<str.size();i++) dataS += str[i] + " ";

  std::cout << "====> Array with strings!\n";
}

XmlParser::~XmlParser() {}


bool XmlParser::getString(const std::string key, std::string &value)
{
    if(!valid_key(key)) return false;

    value = par_map.find(key)->second;
    return true;
}

bool XmlParser::getVectorString(const std::string key, std::vector<std::string> &value)
{
    if(!valid_key(key)) return false;

    std::string dataS = par_map.find(key)->second;

    getStringsFromStringVector(dataS, value);

    return true;
}

void XmlParser::getStringsFromStringVector(const std::string &dataS, std::vector<std::string> &value)
{
    bool found_separator = true; // comma
    int i1 = 0; // index of first \" of the string
    while (i1 < dataS.size())
    {
      int i2 = -1; // index of second \" of the string



      if (!found_separator && dataS[i1] == '\"') throw std::ios_base::failure(errMsgAtLine(fname, line_number, "Missing separator: , ."));

      if (dataS[i1] == ',')
      {
        if (found_separator) throw std::ios_base::failure(errMsgAtLine(fname, line_number, "No string between two consecuive separators , ."));
        found_separator = true;
        i1++;
        continue;
      }

      if (dataS[i1]==' ' || dataS[i1]=='\t')
      {
        i1++;
        continue;
      }

      if (dataS[i1]!='\"') throw std::ios_base::failure(errMsgAtLine(fname, line_number, "Invalid symbol outside of string quotes."));

      // a single separator with white space characters only before and after the separator and at the end a \" character
      for (int i=i1+1;i<dataS.size();i++)
      {
        if (dataS[i] == '\"')
        {
          i2 = i;
          break;
        }
      }

      if (i2 < 0) throw std::ios_base::failure("XmlParser::getVectorString]: Missing \".");

      i1++;
      i2--;

      if (i1==i2) value.push_back("");
      else value.push_back(dataS.substr(i1,i2-i1+1));

      i1 = i2+2;

      found_separator = false;
    }

    // std::string col,str=par_map.find(key)->second;
    // std::istringstream full_str(str);
    // int K = static_cast<int>(std::count(str.begin(),str.end(),',')+1);
    // value.resize(K);
    // for(int k=0; k<K; k++)
    // {
    //     std::getline(full_str, col, ',');
    //     std::stringstream iss(col);
    //     iss >> value[k];
    // }
}


bool XmlParser::getCxCol(const std::string key, arma::cx_vec &value)
{
    if(!valid_key(key)) return false;

    std::string row,str=par_map.find(key)->second;
    std::istringstream full_str(str);
    int K = static_cast<int>(std::count(str.begin(),str.end(),';')+1);
    arma::cx_vec x(K);
    for(int k=0; k<K; k++)
    {
        std::getline(full_str, row, ';');
        x(k) = parse_cx(row);
    }
    value = x;
    return true;
}


bool XmlParser::getCxRow(const std::string key, arma::cx_rowvec &value)
{
    if(!valid_key(key)) return false;

    std::string col,str=par_map.find(key)->second;
    std::istringstream full_str(str);
    int K = static_cast<int>(std::count(str.begin(),str.end(),',')+1);
    arma::cx_rowvec x(K);
    for(int k=0; k<K; k++)
    {
        std::getline(full_str, col, ',');
        x(k) = parse_cx(col);
    }
    value = x;
    return true;
}


bool XmlParser::getCxMat(const std::string key, arma::cx_mat &value)
{
    if(!valid_key(key)) return false;

    std::string full_str,row,col;
    std::istringstream iss_full;

    full_str = par_map.find(key)->second;
    int R = static_cast<int>(std::count(full_str.begin(),full_str.end(),';')+1);

    iss_full.str(full_str);
    std::getline(iss_full, row, ';');
    int C = static_cast<int>(std::count(row.begin(),row.end(),',')+1);

    arma::cx_mat x(R,C);

    iss_full.seekg(0,iss_full.beg);
    for(int r=0; r<R; r++)
    {
        std::getline(iss_full, row, ';');
        std::istringstream iss_row(row);
        for(int c=0; c<C; c++)
        {
            std::getline(iss_row, col, ',');
            x(r,c)=parse_cx(col);
        }
    }
    value = x;
    return true;
}

bool XmlParser::getParam(const std::string key, std::vector<bool> &value)
{
  arma::Row<unsigned> v;
  bool ret = getRow(key, v);
  value.resize(v.size());
  for (int i=0;i<v.size();i++) value[i] = (bool)(v(i));

  return ret;
}

} // namespace gmp_

} // namespace as64_
