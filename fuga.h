#ifndef FUGA_H_
#define FUGA_H_

#include "hoge.h"
#include "lib.h"

class Fuga {
 public:
  Fuga(string R1SamNom, string R2SamNom, string RefFaNom)
      : kR1SamNom_(R1SamNom), kR2SamNom_(R2SamNom), kRefFaNom_(RefFaNom),
        kJudgeFrameshift_(0), kIndelRate_(0) {}
  ~Fuga() {}
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
  // {read ID, R1 readの全alignmentのvector}
  unordered_map<string, vector<SamL> > R1Sam_;
  // {read ID, R2 readの全alignmentのvector}
  unordered_map<string, vector<SamL> > R2Sam_;
  // {ref ID, RefInfo}
  unordered_map<string, RefInfo> RefInfoM_;
  // {readがtop hitした領域の平均iden * readがtop hitした領域の割合, [ref ID]}
  // iden * covが大きい順にref IDをソートするための変数
  map<double, vector<string>, greater<double> > RefScore_;
};

#endif  // FUGA_H_
