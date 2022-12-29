#ifndef AS64_gmp_IO_H
#define AS64_gmp_IO_H

#include <iostream>
#include <fstream>
#include <memory>
#include <iomanip>
#include <cerrno>
#include <streambuf>
#include <armadillo>

#include <gmp_lib/GMP/GMP.h>
#include <gmp_lib/GMP/GMPo.h>
#include <gmp_lib/io/file_io.h>

namespace as64_
{

namespace gmp_
{

/** Write the GMP model to a file.
   * @param[in] gmp: Pointer to a @GMP object.
   * @param[in] filename: String with the name of the file.
   * @param[in] prefix: The prefix that will be used for writing the names of all GMP params (optional, default="").
   */
  void write(const gmp_::GMP *gmp, const std::string &filename, const std::string &prefix="");

  /** Write the GMP model to a file.
   * @param[in] gmp: Pointer to a @GMP object.
   * @param[in] fid: Object of type @FileIO associated with the file.
   * @param[in] prefix: The prefix that will be used for writing the names of all GMP params (optional, default="").
   */
  void write(const gmp_::GMP *gmp, gmp_::FileIO &fid, const std::string &prefix="");

  /** Reads the GMP model from a file.
   * @param[in] gmp: Pointer to a preallocated @GMP object.
   * @param[in] fid: Filename string or object of type @FileIO associated with the file.
   * @param[in] prefix: The prefix that will be used for reading the names of all GMP params (optional, default="").
   */
  void read(gmp_::GMP *gmp, const std::string &filename, const std::string &prefix="");

  /** Reads the GMP model from a file.
   * @param[in] gmp: Pointer to a preallocated @GMP object.
   * @param[in] fid: Filename string or object of type @FileIO associated with the file.
   * @param[in] prefix: The prefix that will be used for reading the names of all GMP params (optional, default="").
   */
  void read(gmp_::GMP *gmp, gmp_::FileIO &fid, const std::string &prefix="");


  /** Write the GMP model to a file.
   * @param[in] gmp: Pointer to a @GMPo object.
   * @param[in] filename: String with the name of the file.
   * @param[in] prefix: The prefix that will be used for writing the names of all GMP params (optional, default="").
   */
  void write(const gmp_::GMPo *gmp, const std::string &filename, const std::string &prefix="");

  /** Write the GMP model to a file.
   * @param[in] gmp: Pointer to a @GMPo object.
   * @param[in] fid: Object of type @FileIO associated with the file.
   * @param[in] prefix: The prefix that will be used for writing the names of all GMP params (optional, default="").
   */
  void write(const gmp_::GMPo *gmp, gmp_::FileIO &fid, const std::string &prefix="");

  /** Reads the GMP model from a file.
   * @param[in] gmp: Pointer to a preallocated @GMPo object.
   * @param[in] fid: Filename string or object of type @FileIO associated with the file.
   * @param[in] prefix: The prefix that will be used for reading the names of all GMP params (optional, default="").
   */
  void read(gmp_::GMPo *gmp, const std::string &filename, const std::string &prefix="");

  /** Reads the GMP model from a file.
   * @param[in] gmp: Pointer to a preallocated @GMPo object.
   * @param[in] fid: Filename string or object of type @FileIO associated with the file.
   * @param[in] prefix: The prefix that will be used for reading the names of all GMP params (optional, default="").
   */
  void read(gmp_::GMPo *gmp, gmp_::FileIO &fid, const std::string &prefix="");

} // namespace gmp_

} // namespace as64_


#endif // AS64_gmp_IO_H
