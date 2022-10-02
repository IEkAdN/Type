#ifndef TYPE_H_
#define TYPE_H_

#include "split.h"
#include "alignment.h"

class Type {
 public:
  Type(string R1SamNom, string R2SamNom, string RefFaNom)
      : kR1SamNom_(R1SamNom), kR2SamNom_(R2SamNom), kRefFaNom_(RefFaNom),
        kJudgeFrameshift_(0), kIndelRate_(0) {}
  ~Type() {}
  void Main();
  const string kR1SamNom_;
  const string kR2SamNom_;
  const string kRefFaNom_;
  bool kJudgeFrameshift_;
  const double kIndelRate_;

 private:
  void ReadRef();
  void ReadSam(string FNom, unordered_map<string, vector<SamL> >* Sam);
  void EditRefInfoM(unordered_map<string, vector<SamL> >* Sam);
  void SummarizeRefInfo();
  void SetRefScore();
  void PrintHit();
  // {read ID: vector<each alignment info of the read>}
  unordered_map<string, vector<SamL> > R1Sam_;
  // {read ID: vector<each alignment info of the read>}
  unordered_map<string, vector<SamL> > R2Sam_;
  // {ref ID, RefInfo}
  unordered_map<string, RefInfo> RefInfoM_;
  // {identity * coverage: vector<ref ID>}
  // identity: average identity in region covered by top hit reads
  // coverage: ratio of region covered by top hit reads
  // variable to sort ref IDs by identity * coverage
  // use score as key to sort ref IDs by it
  map<double, vector<string>, greater<double> > RefScore_;
};

#endif  // TYPE_H_
