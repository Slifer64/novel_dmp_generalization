#ifndef AS64_GMP_LIB_FILE_IO_H
#define AS64_GMP_LIB_FILE_IO_H

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <vector>
#include <exception>
#include <iomanip>
#include <type_traits>
#include <map>

#include <armadillo>
#include <Eigen/Dense>

namespace as64_
{

namespace gmp_
{

typedef long long_t; // type for matrix/vector dimensions
typedef unsigned long size_t_; // size_t type

#define FILE_IO_fun_ std::string("[FileIO::") + __func__ + "]: "

// -------------------------------------------------------------

// Warning: Add new 'Type' or  new 'ScalarType' at the end of the enums, to preserve backwards compatibility!

enum Type
{
  SCALAR = 0,   // scalar, see @ScalarType
  ARMA = 1,     // arma (Mat, Col, Row)
  EIGEN = 2,    // Eigen (Matrix, Vecotr, RowVector)
  STD_VEC = 3,  // std::vector
  STD_STRING = 4 // std::string
  // CUSTOM ?
};

enum ScalarType
{

  BOOL = 0,     // bool
  INT = 1,      // int
  UINT = 2,     // unsigned int
  LONG = 3,     // long
  ULONG = 4,    // unsigned long
  LLONG = 5,    // long long
  ULLONG = 6,   // unsigned long long
  FLOAT = 7,    // float
  DOUBLE = 8,   // double
  UNKNOWN = 9,  // unsupported type
  ScType_NA = 10 // for Types where the scalar type is redundant, like strings
};

enum Error
{
  EMPTY_HEADER = 0,
  CORRUPTED_HEADER,
  ENTRY_NOT_EXIST,
  DUPLICATE_ENTRY,
  TYPE_MISMATCH,
  DIMENSIONS_MISMATCH,
  INVALID_OP_FOR_OPENMODE,
  UNKNOWN_TYPE
};


// ################################
// #######   CLASS FileIO   #######
// ################################


class FileIO
{
// ============================
// =======   PUBLIC  ==========
// ============================
public:

  enum OpenMode
  {
    in = 1<<0,     ///< Allow input operations.
    out = 1<<1,    ///< Allow output operations.
    trunc = 1<<2   ///< Discard all previous contents from the file.
    // binary = 1<<3
  };

  /* Opens a stream for input/output to a file.
   * If the file exist, it opens this files and reads it header.
   * If it doesn't exist, it creates a new file with an empty header.
   * @param[in] filename: the name of the file.
   * @param[in] open_mode: Determines the open mode. Combination of flags from @FileIO::OpenMode. (optional, default=in|out)
   */
  FileIO(const std::string &filename, int open_mode = FileIO::in|FileIO::out);

  ~FileIO();

  /* Closes the IO stream to the file.
   * No more IO operations can be performed afterwards.
   */
  void close() const { this->fid.close(); }

// -------------------------------------------------------------

  /* Writes a scalar to the file.
   * @param[in] name_: the name of the scalar.
   * @param[in] s: the scalar.
   */
  template<typename T>
  void write(const std::string &name_, T s)
  {
    if (!out_flag) throw std::runtime_error(FILE_IO_fun_ + getErrMsg(INVALID_OP_FOR_OPENMODE) + ": \"" + getOpenModeName() + "\"");

    int i = findNameIndex(name_);
    if (i>=0) throw std::runtime_error(FILE_IO_fun_ + getErrMsg(DUPLICATE_ENTRY) + ": \"" + name_ + "\"");

    current_name = name_;

    this->name.push_back(name_);
    this->type.push_back(Type::SCALAR);
    this->sc_type.push_back(findScalarType(s));
    this->name_map[name_] = this->name.size()-1;
    fid.seekp(this->header_start, fid.beg);
    this->i_pos.push_back(fid.tellp());
    FileIO::writeScalar_(s, this->fid);
    this->header_start = fid.tellp();
    this->writeHeader(); // overwrites previous header
  }

  /* Writes an std::string to the file.
   * @param[in] name_: the name of the string.
   * @param[in] s: the string.
   */
  void write(const std::string &name_, const std::string &s);
  void write(const std::string &name_, const char *s);

  /* Writes an std vector to the file.
   * @param[in] name_: the name of the vector.
   * @param[in] v: std::vector<T> of scalar type T.
   */
  template<typename T>
  void write(const std::string &name_, const std::vector<T> &v)
  {
    if (!out_flag) throw std::runtime_error(FILE_IO_fun_ + getErrMsg(INVALID_OP_FOR_OPENMODE) + ": \"" + getOpenModeName() + "\"");

    long_t n_rows = v.size();
    long_t n_cols = 1;
    writeMatDim<T>(name_, Type::STD_VEC, n_rows, n_cols);

    long_t n_elem = n_rows;
    T *buff = new T[n_elem];
    int k=0;
    for (int i=0;i<n_elem;i++) buff[k++] = v[i];
    this->fid.write(reinterpret_cast<const char *>(buff), n_elem*sizeof(T));
    delete []buff;

    this->header_start = fid.tellp();
    this->writeHeader(); // overwrites previous header
    // the new header is always at least as big as the previous one, so the previous one will be completely overwritten
  }

  /* Writes an arma matrix to the file.
   * @param[in] name_: the name of the matrix.
   * @param[in] m: arma::Mat<T> of scalar type T.
   */
  template<typename T>
  void write(const std::string &name_, const arma::Mat<T> &m)
  {
    if (!out_flag) throw std::runtime_error(FILE_IO_fun_ + getErrMsg(INVALID_OP_FOR_OPENMODE) + ": \"" + getOpenModeName() + "\"");

    long_t n_rows = m.n_rows;
    long_t n_cols = m.n_cols;
    writeMatDim<T>(name_, Type::ARMA, n_rows, n_cols);
    writeMatrix<const arma::Mat<T>, T>(m, n_rows, n_cols);
  }

  /* Writes an arma column vector to the file.
   * @param[in] name_: the name of the vector.
   * @param[in] m: arma::Col<T> of scalar type T.
   */
  template<typename T>
  void write(const std::string &name_, const arma::Col<T> &m)
  {
    if (!out_flag) throw std::runtime_error(FILE_IO_fun_ + getErrMsg(INVALID_OP_FOR_OPENMODE) + ": \"" + getOpenModeName() + "\"");

    long_t n_rows = m.n_rows;
    long_t n_cols = 1;
    writeMatDim<T>(name_, Type::ARMA, n_rows, n_cols);
    writeVector<const arma::Col<T>, T>(m, n_rows);
  }

  /* Writes an arma row vector to the file.
   * @param[in] name_: the name of the vector.
   * @param[in] m: arma::Row<T> of scalar type T.
   */
  template<typename T>
  void write(const std::string &name_, const arma::Row<T> &m)
  {
    if (!out_flag) throw std::runtime_error(FILE_IO_fun_ + getErrMsg(INVALID_OP_FOR_OPENMODE) + ": \"" + getOpenModeName() + "\"");

    long_t n_rows = 1;
    long_t n_cols = m.n_cols;
    writeMatDim<T>(name_, Type::ARMA, n_rows, n_cols);
    writeVector<const arma::Row<T>, T>(m, n_cols);
  }

  /* Writes an Eigen matrix to the file.
   * @param[in] name_: the name of the matrix.
   * @param[in] m: Eigen::Matrix<T> of scalar type T.
   */
  template<typename T>
  void write(const std::string &name_, const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> &m)
  {
    if (!out_flag) throw std::runtime_error(FILE_IO_fun_ + getErrMsg(INVALID_OP_FOR_OPENMODE) + ": \"" + getOpenModeName() + "\"");

    long_t n_rows = m.rows();
    long_t n_cols = m.cols();
    writeMatDim<T>(name_, Type::EIGEN, n_rows, n_cols);
    writeMatrix<const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>, T>(m, n_rows, n_cols);
  }

  /* Writes an Eigen column vector to the file.
   * @param[in] name_: the name of the vector.
   * @param[in] m: Eigen::Vector<T> of scalar type T.
   */
  template<typename T>
  void write(const std::string &name_, const Eigen::Matrix<T,Eigen::Dynamic,1> &m)
  {
    if (!out_flag) throw std::runtime_error(FILE_IO_fun_ + getErrMsg(INVALID_OP_FOR_OPENMODE) + ": \"" + getOpenModeName() + "\"");

    long_t n_rows = m.rows();
    long_t n_cols = 1;
    writeMatDim<T>(name_, Type::EIGEN, n_rows, n_cols);
    writeVector<const Eigen::Matrix<T,Eigen::Dynamic,1>, T>(m, n_rows);
  }

  /* Writes an Eigen row vector to the file.
   * @param[in] name_: the name of the vector.
   * @param[in] m: Eigen::RowVector<T> of scalar type T.
   */
  template<typename T>
  void write(const std::string &name_, const Eigen::Matrix<T,1,Eigen::Dynamic> &m)
  {
    if (!out_flag) throw std::runtime_error(FILE_IO_fun_ + getErrMsg(INVALID_OP_FOR_OPENMODE) + ": \"" + getOpenModeName() + "\"");

    long_t n_rows = 1;
    long_t n_cols = m.cols();
    writeMatDim<T>(name_, Type::EIGEN, n_rows, n_cols);
    writeVector<const Eigen::Matrix<T,1,Eigen::Dynamic>, T>(m, n_cols);
  }

// -------------------------------------------------------------

  /* Reads a scalar from the file.
   * @param[in] name_: the name of the scalar.
   * @param[out] s: the scalar.
   */
  template<typename T>
  void read(const std::string &name_, T &s)
  {
    if (!in_flag) throw std::runtime_error(FILE_IO_fun_ + getErrMsg(INVALID_OP_FOR_OPENMODE) + ": \"" + getOpenModeName() + "\"");

    int i = findNameIndex(name_);
    if (i<0) throw std::runtime_error(FILE_IO_fun_ + getErrMsg(ENTRY_NOT_EXIST) + "\"" + name_ + "\"");

    current_name = name_;

    checkType(i, Type::SCALAR, findScalarType(s));
    size_t_ pos = i_pos[i];
    fid.seekg(pos, fid.beg);

    FileIO::readScalar_(s, this->fid);
  }

  /* Reads an std::string from the file.
   * @param[in] name_: the name of the string.
   * @param[out] s: the string.
   */
  void read(const std::string &name_, std::string &s);

  /* Reads an std vector from the file.
   * @param[in] name_: the name of the vector.
   * @param[out] v: std::vector<T> of scalar type T.
   */
  template<typename T>
  void read(const std::string &name_, std::vector<T> &m)
  {
    if (!in_flag) throw std::runtime_error(FILE_IO_fun_ + getErrMsg(INVALID_OP_FOR_OPENMODE) + ": \"" + getOpenModeName() + "\"");

    long_t n_rows, n_cols;
    readMatDim<T>(name_, Type::STD_VEC, n_rows, n_cols);

    if (n_rows!=1 && n_cols!=1)
      throw std::runtime_error(FILE_IO_fun_ + getErrMsg(DIMENSIONS_MISMATCH) + ": entry \"" + name_
                               + "\": input argument is std::vector, file data is " + std::to_string(n_rows) + " x " + std::to_string(n_cols) + " matrix" );

    long_t n_elem = n_rows*n_cols;
    m.resize(n_elem);
    T *buff = new T[n_elem];
    fid.read((char *)(buff), n_elem*sizeof(T));
    int k=0;
    for (int i=0;i<n_elem;i++) m[i] = buff[k++];
    delete []buff;
  }

  /* Reads an arma matrix from the file.
   * @param[in] name_: the name of the matrix.
   * @param[out] m: arma::Mat<T> of scalar type T.
   */
  template<typename T>
  void read(const std::string &name_, arma::Mat<T> &m)
  {
    if (!in_flag) throw std::runtime_error(FILE_IO_fun_ + getErrMsg(INVALID_OP_FOR_OPENMODE) + ": \"" + getOpenModeName() + "\"");

    long_t n_rows, n_cols;
    readMatDim<T>(name_, Type::ARMA, n_rows, n_cols);
    readMatrix<arma::Mat<T>, T>(m, n_rows, n_cols);
  }

  /* Reads an arma column vector from the file.
   * @param[in] name_: the name of the vector.
   * @param[out] m: arma::Col<T> of scalar type T.
   */
  template<typename T>
  void read(const std::string &name_, arma::Col<T> &m)
  {
    if (!in_flag) throw std::runtime_error(FILE_IO_fun_ + getErrMsg(INVALID_OP_FOR_OPENMODE) + ": \"" + getOpenModeName() + "\"");

    long_t n_rows, n_cols;
    readMatDim<T>(name_, Type::ARMA, n_rows, n_cols);

    if (n_cols != 1)
      throw std::runtime_error(FILE_IO_fun_ + getErrMsg(DIMENSIONS_MISMATCH) + ": "
      + "entry \"" + name_ + "\": input argument is arma::Col, file data is " + std::to_string(n_rows) + " x " + std::to_string(n_cols) + " matrix" );

    readVector<arma::Col<T>, T>(m, n_rows);
  }

  /* Reads an arma row vector from the file.
   * @param[in] name_: the name of the vector.
   * @param[out] m: arma::Row<T> of scalar type T.
   */
  template<typename T>
  void read(const std::string &name_, arma::Row<T> &m)
  {
    if (!in_flag) throw std::runtime_error(FILE_IO_fun_ + getErrMsg(INVALID_OP_FOR_OPENMODE) + ": \"" + getOpenModeName() + "\"");

    long_t n_rows, n_cols;
    readMatDim<T>(name_, Type::ARMA, n_rows, n_cols);

    if (n_rows != 1)
      throw std::runtime_error(FILE_IO_fun_ + getErrMsg(DIMENSIONS_MISMATCH) + ": entry \"" + name_
      + "\": input argument is arma::Row, file data is " + std::to_string(n_rows) + " x " + std::to_string(n_cols) + " matrix" );

    readVector<arma::Row<T>, T>(m, n_cols);
  }

  /* Reads an Eigen matrix from the file.
   * @param[in] name_: the name of the matrix.
   * @param[out] m: Eigen::Matrix<T> of scalar type T.
   */
  template<typename T>
  void read(const std::string &name_, Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> &s)
  {
    if (!in_flag) throw std::runtime_error(FILE_IO_fun_ + getErrMsg(INVALID_OP_FOR_OPENMODE) + ": \"" + getOpenModeName() + "\"");

    long_t n_rows, n_cols;
    readMatDim<T>(name_, Type::EIGEN, n_rows, n_cols);
    readMatrix<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>, T>(s, n_rows, n_cols);
  }

  /* Reads an Eigen column vector from the file.
   * @param[in] name_: the name of the vector.
   * @param[out] m: Eigen::Vector<T> of scalar type T.
   */
  template<typename T>
  void read(const std::string &name_, Eigen::Matrix<T,Eigen::Dynamic,1> &s)
  {
    if (!in_flag) throw std::runtime_error(FILE_IO_fun_ + getErrMsg(INVALID_OP_FOR_OPENMODE) + ": \"" + getOpenModeName() + "\"");

    long_t n_rows, n_cols;
    readMatDim<T>(name_, Type::EIGEN, n_rows, n_cols);

    if (n_cols != 1)
      throw std::runtime_error(FILE_IO_fun_ + getErrMsg(DIMENSIONS_MISMATCH) + ": entry \"" + name_
      + "\": input argument is Eigen::Vector, file data is " + std::to_string(n_rows) + " x " + std::to_string(n_cols) + " matrix" );

    readVector<Eigen::Matrix<T,Eigen::Dynamic,1>, T>(s, n_rows);
  }

  /* Reads an Eigen row vector from the file.
   * @param[in] name_: the name of the vector.
   * @param[out] m: Eigen::RowVector<T> of scalar type T.
   */
  template<typename T>
  void read(const std::string &name_, Eigen::Matrix<T,1,Eigen::Dynamic> &s)
  {
    if (!in_flag) throw std::runtime_error(FILE_IO_fun_ + getErrMsg(INVALID_OP_FOR_OPENMODE) + ": \"" + getOpenModeName() + "\"");

    long_t n_rows, n_cols;
    readMatDim<T>(name_, Type::EIGEN, n_rows, n_cols);

    if (n_rows != 1)
      throw std::runtime_error(FILE_IO_fun_ + getErrMsg(DIMENSIONS_MISMATCH) + ": entry \"" + name_
      + "\": input argument is Eigen::RowVector, file data is " + std::to_string(n_rows) + " x " + std::to_string(n_cols) + " matrix" );

    readVector<Eigen::Matrix<T,1,Eigen::Dynamic>, T>(s, n_cols);
  }

// -------------------------------------------------------------

  /* Prints the header of the file.
   * @param[in] out: the output stream where the header will be printed (optinal, default=std::cout)
   */
  void printHeader(std::ostream &out=std::cout) const;

// =============================================
// =======   Public Static functions  ==========
// =============================================
public:
  template<typename T>
  static void writeScalar_(T s, std::ostream &out)
  {
    out.write(reinterpret_cast<const char *>(&s), sizeof(s));
  }

  template<typename T>
  static void readScalar_(T &s, std::istream &in)
  {
    in.read(reinterpret_cast<char *>(&s), sizeof(s));
  }

// ===============================
// =======   PROTECTED  ==========
// ===============================
protected:
// -------------------------------------------------------------

  int findNameIndex(const std::string &name_) const;

  void checkType(int i, enum Type t2, enum ScalarType sc_t2) const;

  template<typename T>
  ScalarType findScalarType(T s) const
  {
    ScalarType t;
    if (std::is_same<T, bool>::value) t = ScalarType::BOOL;
    else if (std::is_same<T, int>::value) t = ScalarType::INT;
    else if (std::is_same<T, unsigned int>::value) t = ScalarType::UINT;
    else if (std::is_same<T, long>::value) t = ScalarType::LONG;
    else if (std::is_same<T, unsigned long>::value) t = ScalarType::ULONG;
    else if (std::is_same<T, long long>::value) t = ScalarType::LLONG;
    else if (std::is_same<T, unsigned long long>::value) t = ScalarType::ULLONG;
    else if (std::is_same<T, float>::value) t = ScalarType::FLOAT;
    else if (std::is_same<T, double>::value) t = ScalarType::DOUBLE;
    else t = ScalarType::UNKNOWN;

    if (t == ScalarType::UNKNOWN) throw std::runtime_error(FILE_IO_fun_ + "entry \"" + current_name + "\", " + getErrMsg(UNKNOWN_TYPE));

    return t;
  }

  std::string getMatDim2str(int k) const;

// -------------------------------------------------------------

  void writeHeader() const;

  void readHeader();

// -------------------------------------------------------------

  template<typename Mat_t, typename T>
  void readMatrix(Mat_t &s, long_t n_rows, long_t n_cols)
  {
    s.resize(n_rows, n_cols);

    T *buff = new T[n_rows*n_cols];
    fid.read((char *)(buff), n_rows*n_cols*sizeof(T));
    int k=0;
    for (int j=0;j<n_cols;j++){
      for (int i=0;i<n_rows;i++) s(i,j) = buff[k++];
    }
    delete []buff;
  }

  template<typename Mat_t, typename T>
  void readVector(Mat_t &s, long_t n_elem)
  {
    s.resize(n_elem);
    T *buff = new T[n_elem];
    fid.read((char *)(buff), n_elem*sizeof(T));
    int k=0;
    for (int i=0;i<n_elem;i++) s(i) = buff[k++];
    delete []buff;
  }

  template<typename T>
  void readMatDim(const std::string &name_, Type t, long_t &n_rows, long_t &n_cols)
  {
    int i = findNameIndex(name_);
    if (i<0) throw std::runtime_error(FILE_IO_fun_ + getErrMsg(ENTRY_NOT_EXIST) + "\"" + name_ + "\"");

    current_name = name_;

    T sc_t; // temporary variable to determine the scalar type
    checkType(i, t, findScalarType(sc_t));
    size_t_ pos = i_pos[i];
    fid.seekg(pos, fid.beg);

    FileIO::readScalar_(n_rows, this->fid);
    FileIO::readScalar_(n_cols, this->fid);

    if (n_rows<0 || n_cols<0)
      throw std::runtime_error(FILE_IO_fun_ +
                               getErrMsg(CORRUPTED_HEADER) + ": Negative matrix dimensions: \"" +
                               name_ + "\" " + std::to_string(n_rows) + " x " + std::to_string(n_cols) );
  }

  template<typename Mat_t, typename T>
  void writeMatrix(const Mat_t &m, long_t n_rows, long_t n_cols)
  {
    T *buff = new T[n_rows*n_cols];
    int k=0;
    for (int j=0;j<n_cols;j++){
      for (int i=0;i<n_rows;i++) buff[k++] = m(i,j);
    }
    this->fid.write(reinterpret_cast<const char *>(buff), n_rows*n_cols*sizeof(T));
    delete []buff;

    this->header_start = fid.tellp();
    this->writeHeader(); // overwrites previous header
    // the new header is always at least as big as the previous one, so the previous one will be completely overwritten
  }

  template<typename Mat_t, typename T>
  void writeVector(Mat_t &m, long_t n_elem)
  {
    T *buff = new T[n_elem];
    int k=0;
    for (int i=0;i<n_elem;i++) buff[k++] = m(i);
    this->fid.write(reinterpret_cast<const char *>(buff), n_elem*sizeof(T));
    delete []buff;

    this->header_start = fid.tellp();
    this->writeHeader(); // overwrites previous header
    // the new header is always at least as big as the previous one, so the previous one will be completely overwritten
  }

  template<typename T>
  void writeMatDim(const std::string &name_, Type t, long_t n_rows, long_t n_cols)
  {
    int i = findNameIndex(name_);
    if (i>=0) throw std::runtime_error(FILE_IO_fun_ + getErrMsg(DUPLICATE_ENTRY) + ": \"" + name_ + "\"");

    current_name = name_;

    this->name.push_back(name_);
    this->type.push_back(t);
    T temp; // temporary variable to determine the type of T
    this->sc_type.push_back(findScalarType(temp));
    this->name_map[name_] = this->name.size()-1;
    fid.seekp(this->header_start, fid.beg);
    this->i_pos.push_back(fid.tellp());

    FileIO::writeScalar_(n_rows, this->fid);
    FileIO::writeScalar_(n_cols, this->fid);
  }

  std::string getOpenModeName() const;

// -------------------------------------------------------------

  static std::string getFullTypeName(enum Type type, enum ScalarType sc_type);

  static std::string getTypeName(enum Type type);

  static std::string getScalarTypeName(enum ScalarType type);

  static std::string getErrMsg(enum Error err_id);

// -------------------------------------------------------------

  mutable std::fstream fid; // the stream for reading/writing to the file

  long_t header_start; // the position where the header starts in the file

  std::map<std::string, int> name_map; // maps the names of the data to their index in the following vectors
  std::vector<std::string> name; // vector with the names of the data contained in the file
  std::vector<Type> type; // vector with the type (see @Type) of the data contained in the file
  std::vector<ScalarType> sc_type; // vector with the scalar type (see @ScalarType) of the data contained in the file
  std::vector<size_t_> i_pos; // vector with the position of the data in the file

  bool in_flag; // true when the open_mode is "in"
  bool out_flag; // true when the open_mode is "out"

  static const char* error_msg[];
  static const char *TypeName[];
  static const char *ScalarTypeName[];

  std::string current_name; // stores the name of the entry that is currently proccessed (used to for display purposes on errors)

};

} // namespace gmp_

} // namespace as64_

#endif // AS64_GMP_LIB_FILE_IO_H
