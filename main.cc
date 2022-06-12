#include "fuga.h"

int main(int argc, const char* argv[]) {
  if (argc != 4) {
    cerr << "usage: a [in.R1.sam] [in.R2.sam] [in.Ref.fa]\n";
    return 1;
  } else {
    Fuga fuga(argv[1], argv[2], argv[3]);
    fuga.Main();
  }
  return 0;
}
