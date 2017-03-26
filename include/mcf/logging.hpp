// vim: expandtab:ts=2:sw=2
#ifndef MCF_LOGGING_HPP
#define MCF_LOGGING_HPP

#include <iostream>

namespace mcf {

/***
 * A small logging class that either prints to standard output or suppresses
 * messages.
 */
class Printer {
 public:
  /**
   * Constructor.
   *
   * @param verbose If true, print messages to standard output. Otherwise,
   *        suppress messages.
   */
  Printer(bool verbose = false);

  //! Set verbosity state.
  void set_verbose(bool verbose);

  //! Get verbosity state.
  bool is_verbose() const;

  /**
   * Log message.
   *
   * @param arg The first argument to print.
   * @param args Following arguments to print.
   */
  template <typename T, typename... Args>
  void Print(const T& arg, const Args&... args) const;

  /**
   * Log message.
   *
   * @param arg One argument to print.
   */
  template <typename T>
  void Print(const T& arg) const;

 private:
  bool verbose_;
};

//! \return Returns a global Printer object for printing log messages.
Printer& printer();

//! Print log message to global Printer object (convenience function).
template <typename... Args>
void print(const Args&... args);

}  // namespace mcf

#include <mcf/impl/logging_inl.hpp>

#endif
