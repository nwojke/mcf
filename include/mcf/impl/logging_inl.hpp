// vim: expandtab:ts=2:sw=2
#ifndef MCF_LOGGING_INL_HPP
#define MCF_LOGGING_INL_HPP

#include <mcf/logging.hpp>

namespace mcf {

inline Printer::Printer(const bool verbose) : verbose_(verbose) {}

inline void Printer::set_verbose(const bool verbose) { verbose_ = verbose; }

inline bool Printer::is_verbose() const { return verbose_; }

template <typename T, typename... Args>
void Printer::Print(const T& arg, const Args&... args) const {
  if (!verbose_) {
    return;
  }
  std::cout << arg;
  this->Print(args...);
}

template <typename T>
void Printer::Print(const T& arg) const {
  if (!verbose_) {
    return;
  }
  std::cout << arg << std::endl;
}

inline Printer& printer() {
  static Printer printer;
  return printer;
}

template <typename... Args>
void print(const Args&... args) {
  printer().Print(args...);
}

}  // namespace mcf

#endif
