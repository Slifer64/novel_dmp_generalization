#ifndef AS64_XML_PARSER_H
#define AS64_XML_PARSER_H

#include <iostream>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <sstream>
#include <exception>
#include <algorithm>
#include <map>

#include <armadillo>

namespace as64_
{

namespace io_
{
  ///
  /// @defgroup parser Parser
  /// \brief Parameter file parser functions.
  /// @{

  ///
  /// \brief A parser class.
  ///
  /// Implements parsing from text file for different types
  ///
  class XmlParser
  {
    private:
      std::map<std::string, std::string> par_map; ///< Map structure to store parameter and value

      ////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief Check if key is valid.
      /// @returns TRUE if key is valid FALSE otherwise
      /// @param key Key string
      ////////////////////////////////////////////////////////////////////////////////////////////
      bool valid_key(const std::string& key);


      ////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief Parse complex value.
      ///
      /// String may be in format e.g. "1+1i","1+1j" or entirely real or imaginary "1", "1i" or "1j"
      /// @returns Complex value
      /// @param str Complex notation string
      ////////////////////////////////////////////////////////////////////////////////////////////
      std::complex<double> parse_cx(std::string str);

      int line_number = 0;
      int column_number = 0;
      std::string fname;

      void validateArrayValues(std::string &dataS);

      void getStringsFromStringVector(const std::string &dataS, std::vector<std::string> &value);

    public:
      ////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief Constructor.
      ///
      /// Opens parameter file and puts the key and value in the map structure
      /// @param fname Parameter file name
      ////////////////////////////////////////////////////////////////////////////////////////////
      XmlParser(const std::string& fname);


      ////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief Destructor.
      ////////////////////////////////////////////////////////////////////////////////////////////
      ~XmlParser();

      ////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief Generic type get function.
      /// @returns Value of type T, if it not exists in file the default value is returned
      /// @param key Parameter name
      /// @param def_val Default value if key was not found in file
      ////////////////////////////////////////////////////////////////////////////////////////////
      template <typename T>
      bool getParam(const std::string key, T &value)
      {
        if(!valid_key(key)) return false;

        std::istringstream iss(par_map.find(key)->second);
        iss >> value;
        return true;
      }

      bool getParam(const std::string key, std::string &value)
      {
        return getString(key, value);
      }

      template <typename T>
      bool getParam(const std::string key, std::vector<T> &value)
      {
        arma::Row<T> v;
        bool ret = getRow(key, v);
        value.resize(v.size());
        std::copy(v.begin(), v.end(), value.begin());

        return ret;
      }

      bool getParam(const std::string key, std::vector<bool> &value);

      template <typename T>
      bool getParam(const std::string key, arma::Col<T> &value)
      {
        return getCol(key, value);
      }

      template <typename T>
      bool getParam(const std::string key, arma::Row<T> &value)
      {
        return getRow(key, value);
      }

      template <typename T>
      bool getParam(const std::string key, arma::Mat<T> &value)
      {
        return getMat(key, value);
      }

      bool getParam(const std::string key, arma::cx_vec &value)
      {
        return getCxCol(key, value);
      }

      bool getParam(const std::string key, arma::cx_rowvec &value)
      {
        return getCxRow(key, value);
      }

      bool getParam(const std::string key, arma::cx_mat &value)
      {
        return getCxMat(key, value);
      }

      bool getParam(const std::string key, std::vector<std::string> &value)
      {
        return getVectorString(key, value);
      }

      ////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief String type get function.
      /// @returns String value, if it not exists in file the default value is returned
      /// @param key Parameter name
      /// @param def_val Default value if key was not found in file
      ////////////////////////////////////////////////////////////////////////////////////////////
      bool getString(const std::string key, std::string &value);


      ////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief Col type get function.
      /// @returns Col value, if it not exists in file the default value is returned
      /// @param key Parameter name
      /// @param def_val Default value if key was not found in file
      ////////////////////////////////////////////////////////////////////////////////////////////
      template <typename T>
      bool getCol(const std::string key, arma::Col<T> &value)
      {
          if(!valid_key(key)) return false;

          std::string row,str=par_map.find(key)->second;
          std::istringstream full_str(str);
          int K = static_cast<int>(std::count(str.begin(),str.end(),';')+1);
          arma::Col<T> x(K);
          for(int k=0; k<K; k++)
          {
              std::getline(full_str, row, ';');
              std::stringstream iss(row);
              iss >> x(k);
          }
          value = x;
          return true;
      }

      ////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief cx_vec type get function.
      /// @returns cx_vec value, if it not exists in file the default value is returned
      /// @param key Parameter name
      /// @param def_val Default value if key was not found in file
      ////////////////////////////////////////////////////////////////////////////////////////////
      bool getCxCol(const std::string key, arma::cx_vec &value);


      ////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief Row type get function.
      /// @returns Row value, if it not exists in file the default value is returned
      /// @param key Parameter name
      /// @param def_val Default value if key was not found in file
      ////////////////////////////////////////////////////////////////////////////////////////////
      template <typename T>
      bool getRow(const std::string key, arma::Row<T> &value)
      {
          if(!valid_key(key)) return false;

          std::string col,str=par_map.find(key)->second;
          std::istringstream full_str(str);
          int K = static_cast<int>(std::count(str.begin(),str.end(),',')+1);
          arma::Row<T> x(K);
          for(int k=0; k<K; k++)
          {
              std::getline(full_str, col, ',');
              std::stringstream iss(col);
              iss >> x(k);
          }
          value = x;
          return true;
      }

      ////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief cx_rowvec type get function.
      /// @returns cx_rowvec value, if it not exists in file the default value is returned
      /// @param key Parameter name
      /// @param def_val Default value if key was not found in file
      ////////////////////////////////////////////////////////////////////////////////////////////
      bool getCxRow(const std::string key, arma::cx_rowvec &value);

      ////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief Mat type get function.
      /// @returns Mat value, if it not exists in file the default value is returned
      /// @param key Parameter name
      /// @param def_val Default value if key was not found in file
      ////////////////////////////////////////////////////////////////////////////////////////////
      template <typename T>
      bool getMat(const std::string key, arma::Mat<T> &value)
      {
          if(!valid_key(key)) return false;

          std::string full_str,row,col;
          std::istringstream iss_full;

          full_str = par_map.find(key)->second;
          int R = static_cast<int>(std::count(full_str.begin(),full_str.end(),';')+1);

          iss_full.str(full_str);
          std::getline(iss_full, row, ';');
          int C = static_cast<int>(std::count(row.begin(),row.end(),',')+1);

          arma::Mat<T> x(R,C);

          iss_full.seekg(0,iss_full.beg);
          for(int r=0; r<R; r++)
          {
              std::getline(iss_full, row, ';');
              std::istringstream iss_row(row);
              for(int c=0; c<C; c++)
              {
                  std::getline(iss_row, col, ',');
                  std::istringstream iss_col(col);
                  iss_col >> x(r,c);
              }
          }
          value = x;
          return true;
      }

      ////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief cx_mat type get function.
      /// @returns cx_mat value, if it not exists in file the default value is returned
      /// @param key Parameter name
      /// @param def_val Default value if key was not found in file
      ////////////////////////////////////////////////////////////////////////////////////////////
      bool getCxMat(const std::string key, arma::cx_mat &value);

      bool getVectorString(const std::string key, std::vector<std::string> &value);

  }; // end Class
  /// @}

} // namespace io_

} // namespace as64_

#endif // AS64_XML_PARSER_H
