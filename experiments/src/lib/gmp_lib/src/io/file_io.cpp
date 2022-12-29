#include <gmp_lib/io/file_io.h>

namespace as64_
{

namespace gmp_
{

const char* FileIO::error_msg[] =
{
  "EMPTY_HEADER",
  "CORRUPTED_HEADER",
  "ENTRY DOES NOT EXIST:",
  "DUPLICATE_ENTRY",
  "TYPE_MISMATCH",
  "DIMENSIONS_MISMATCH",
  "INVALID_OP_FOR_OPENMODE",
  "UNKNOWN_TYPE"
};

const char *FileIO::TypeName[] =
{
  "scalar",
  "arma::Mat",
  "Eigen::Matrix",
  "std::vector",
  "std::string"
};

const char *FileIO::ScalarTypeName[] =
{
  "bool",
  "int",
  "unsigned int",
  "long",
  "unsigned long",
  "long long",
  "unsigned long long",
  "float",
  "double",
  "N/A",
  "",  // ScType_NA
};


FileIO::FileIO(const std::string &filename, int open_mode)
{
  this->header_start = 0;

  in_flag = open_mode & FileIO::in;
  out_flag = open_mode & FileIO::out;
  bool trunc_flag = open_mode & FileIO::trunc;

  bool file_exist;
  { // check if file exists
    std::ifstream in(filename, std::ios::in);
    file_exist = in.good();
  }

  if (!out_flag)
  {
    if (!in_flag) throw std::runtime_error(FILE_IO_fun_ + "Open-mode \"in\" and/or \"out\" must be specified!");
    if (in_flag && trunc_flag) throw std::runtime_error(FILE_IO_fun_ + "Setting only \"in\" and \"trunc\" makes no sense :P ...");
    if (in_flag && !file_exist) throw std::runtime_error(FILE_IO_fun_ + "Open-mode \"in\" was specified but file \"" + filename + "\" doesn't exist or cannot be accessed...");
  }

  if (trunc_flag && file_exist)
  {
    if (std::remove(filename.c_str()) != 0) // or fid.open(filename, ios::out | ios::trunc) ?
    // if (!std::ofstream(filename, std::ios::out | std::ios::trunc))  // have to also write empty header...
      throw std::runtime_error(FILE_IO_fun_ + "Cannot discard the previous contents of the file \"" + filename + "\"...");
    file_exist = false;
  }

  if (out_flag && !file_exist) // create the file and add empty header
  {
    fid.open(filename, std::ios::out);
    if (!fid) throw std::runtime_error(FILE_IO_fun_ + "Failed to create file \"" + filename + "\"...");
    this->writeHeader(); // write empty header
    fid.close();
  }

  std::ios::openmode op_mode = std::ios::binary;
  if (in_flag) op_mode |= std::ios::in;
  if (out_flag) op_mode |= std::ios::out | std::ios::in; // add "in" to avoid overriding previous contents

  fid.open(filename, op_mode);
  if (!fid) throw std::ios_base::failure(std::string(FILE_IO_fun_ + "Failed to open file \"") + filename + "\". The file "
                                                      "may not exist or it cannot be accessed...\n");

  this->readHeader();
}

FileIO::~FileIO()
{
  fid.close();
}

void FileIO::write(const std::string &name_, const char *s)
{
  write(name_, std::string(s) );
}

void FileIO::write(const std::string &name_, const std::string &s)
{
  if (!out_flag) throw std::runtime_error(FILE_IO_fun_ + getErrMsg(INVALID_OP_FOR_OPENMODE) + ": \"" + getOpenModeName() + "\"");

  int i = findNameIndex(name_);
  if (i>=0) throw std::runtime_error(FILE_IO_fun_ + getErrMsg(DUPLICATE_ENTRY) + ": \"" + name_ + "\"");

  current_name = name_;

  this->name.push_back(name_);
  this->type.push_back(Type::STD_STRING);
  this->sc_type.push_back(ScalarType::ScType_NA);
  this->name_map[name_] = this->name.size()-1;
  fid.seekp(this->header_start, fid.beg);
  this->i_pos.push_back(fid.tellp());

  // write string:
  long_t n_elem = s.size();
  this->fid.write(reinterpret_cast<const char *>(&n_elem), sizeof(n_elem) );
  this->fid.write(s.data(), n_elem*sizeof(char) );

  this->header_start = fid.tellp();
  this->writeHeader(); // overwrites previous header
}

void FileIO::read(const std::string &name_, std::string &s)
{
  if (!in_flag) throw std::runtime_error(FILE_IO_fun_ + getErrMsg(INVALID_OP_FOR_OPENMODE) + ": \"" + getOpenModeName() + "\"");

  int i = findNameIndex(name_);
  if (i<0) throw std::runtime_error(FILE_IO_fun_ + getErrMsg(ENTRY_NOT_EXIST) + "\"" + name_ + "\"");

  current_name = name_;

  checkType(i, Type::STD_STRING, ScalarType::ScType_NA);
  size_t_ pos = i_pos[i];
  fid.seekg(pos, fid.beg);
  // read string
  long_t n_elem = s.size();
  this->fid.read(reinterpret_cast<char *>(&n_elem), sizeof(n_elem) );
  s.resize(n_elem);
  this->fid.read(&s[0], n_elem*sizeof(char) );
}


void FileIO::printHeader(std::ostream &out) const
{
  int name_len = 0;
  for (int k=0; k<name.size(); k++)
  {
    if (name[k].length() > name_len) name_len = name[k].length();
  }
  name_len += 5;

  std::vector<std::string> type_name(type.size());

  int type_len = 0;
  for (int k=0; k<type.size(); k++)
  {
    type_name[k] = getFullTypeName(type[k], sc_type[k]);
    if (type[k]==ARMA || type[k]==EIGEN) type_name[k] += getMatDim2str(k);
    int n = type_name[k].length();
    if (n > type_len) type_len = n;
  }
  type_len += 5;

  std::string horiz_line;
  horiz_line.resize(name_len + type_len + 6);
  for (int i=0; i<horiz_line.length(); i++) horiz_line[i] = '-';

  out << horiz_line << "\n";
  out << std::setw(name_len) << std::left << "Name" << std::setw(type_len) << std::left << "Type" << "Pos\n";
  out << horiz_line << "\n";
  for (int k=0; k<name.size(); k++)
  {
    out << std::setw(name_len) << std::left << name[k]
        << std::setw(type_len) << std::left << type_name[k]
        << i_pos[k] << "\n";
  }
}

int FileIO::findNameIndex(const std::string &name_) const
{
  auto it = name_map.find(name_);
  if (it != name_map.end()) return it->second;
  else return -1;
}

void FileIO::checkType(int i, enum Type t2, enum ScalarType sc_t2) const
{
  Type t = type[i];
  ScalarType sc_t = sc_type[i];

  if (t!=t2 | sc_t!=sc_t2)
    throw std::runtime_error(FILE_IO_fun_ + "entry \"" + current_name + "\", " + getErrMsg(TYPE_MISMATCH) +
                             ": " + getFullTypeName(t,sc_t) + " != " + getFullTypeName(t2,sc_t2));
}

std::string FileIO::getMatDim2str(int k) const
{
  fid.seekg(i_pos[k], fid.beg);

  long_t n_rows;
  long_t n_cols;
  FileIO::readScalar_(n_rows, this->fid);
  FileIO::readScalar_(n_cols, this->fid);

  std::ostringstream out;
  out << "(" << n_rows << "," << n_cols << ")";
  return out.str();
}

void FileIO::writeHeader() const
  {
    this->fid.seekp(this->header_start, this->fid.beg);

    size_t_ i1 = this->fid.tellp();

    for (int k=0; k<name.size(); k++)
    {
      Type t = type[k];
      ScalarType sc_t = sc_type[k];
      size_t_ i = i_pos[k];
      long_t name_len = name[k].length();
      this->fid.write(reinterpret_cast<const char *>(&name_len), sizeof(name_len));
      this->fid.write(reinterpret_cast<const char *>(name[k].c_str()), name_len);
      this->fid.write(reinterpret_cast<const char *>(&t), sizeof(t));
      this->fid.write(reinterpret_cast<const char *>(&sc_t), sizeof(sc_t));
      this->fid.write(reinterpret_cast<const char *>(&i), sizeof(i));
    }

    size_t_ i2 = this->fid.tellp();

    long_t header_len = i2-i1 + sizeof(header_len);
    this->fid.write(reinterpret_cast<const char *>(&header_len), sizeof(header_len));
    this->fid.flush();
  }

void FileIO::readHeader()
{
  this->name_map.clear();
  this->name.clear();
  this->type.clear();
  this->sc_type.clear();
  this->i_pos.clear();

  fid.seekg (0, fid.beg);
  size_t_ i_start = fid.tellg();
  fid.seekg (0, fid.end);
  size_t_ i_end = fid.tellg();
  if (i_start == i_end) throw std::runtime_error(FILE_IO_fun_ + FileIO::getErrMsg(EMPTY_HEADER));

  long_t header_len;
  this->fid.seekg(-sizeof(header_len), this->fid.end);
  i_end = this->fid.tellg();
  this->fid.read(reinterpret_cast<char *>(&header_len), sizeof(header_len));
  if (header_len < 0) throw std::runtime_error(FILE_IO_fun_ + FileIO::getErrMsg(CORRUPTED_HEADER));

  this->fid.seekg(-header_len, this->fid.end);
  this->header_start = this->fid.tellg();
  while (this->fid.tellg() != i_end)
  {
    Type t;
    ScalarType sc_t;
    size_t_ i;
    long_t len;
    char *name_i;
    this->fid.read(reinterpret_cast<char *>(&len), sizeof(len));
    if (len < 0) throw std::runtime_error(FILE_IO_fun_ + FileIO::getErrMsg(CORRUPTED_HEADER));
    name_i = new char[len+1];
    this->fid.read(reinterpret_cast<char *>(name_i), len);
    name_i[len] = '\0';
    this->fid.read(reinterpret_cast<char *>(&t), sizeof(t));
    this->fid.read(reinterpret_cast<char *>(&sc_t), sizeof(sc_t));
    this->fid.read(reinterpret_cast<char *>(&i), sizeof(i));

    type.push_back(t);
    sc_type.push_back(sc_t);
    i_pos.push_back(i);
    name.push_back(name_i);
    name_map[name_i] = name.size()-1;

    delete name_i;
  }

}

std::string FileIO::getOpenModeName() const
{
  std::string op_mode;
  if (in_flag) op_mode = "in";
  else if (out_flag) op_mode = "out";
  else op_mode = "in|out";
  return op_mode;
}

// ======================================
// =======   STATIC Functions  ==========
// ======================================

std::string FileIO::getFullTypeName(enum Type type, enum ScalarType sc_type)
{
  std::string sc_t = getScalarTypeName(sc_type);
  if (sc_t.empty()) return getTypeName(type);
  else return getTypeName(type) + "<" + sc_t + ">";
}

std::string FileIO::getTypeName(enum Type type)
{
  return std::string(FileIO::TypeName[static_cast<int>(type)]);
}

std::string FileIO::getScalarTypeName(enum ScalarType type)
{
  return std::string(FileIO::ScalarTypeName[static_cast<int>(type)]);
}

std::string FileIO::getErrMsg(enum Error err_id)
{
  return std::string(FileIO::error_msg[err_id]);
}


} // namespace gmp_

} // namespace as64_
