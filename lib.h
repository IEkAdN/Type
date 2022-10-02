#ifndef LIB_H_
#define LIB_H_

#include "hoge.h"

// each ref pos has one RefPosInfo object
class RefPosInfo {
 public:
  RefPosInfo() : AvgIden_(0), Ins_(0), Dp_(0), AllDp_(0), Pos_(0) {}
  ~RefPosInfo() {}
  void EditIden(double Iden) { Iden_.emplace_back(Iden); }
  void EditDel(u32 Len) { ++Del_[Len]; }
  void EditIns() { ++Ins_; }
  void EditDp() { ++Dp_; }
  void EditAllDp() { ++AllDp_; }
  void SetAvgIden() { AvgIden_ = accumulate(Iden_.begin(), Iden_.end(), 0.0) / Dp_; }
  void SetPos(u32 Pos) { Pos_ = Pos; }
  unordered_map<u32, u32> Del() const { return Del_; }
  u32 Ins() const { return Ins_; }
  double AvgIden() const { return AvgIden_; }
  u32 Dp() const { return Dp_; }
  u32 AllDp() const { return AllDp_; }
  u32 Pos() const { return Pos_; }

 private:
  // vector<identity b/w the pos and each top hit read>
  // appended even if the read supports ref deletion at the pos
  vector<double> Iden_;
  double AvgIden_;
  // number of top hit reads supporting ref insertion
  // for experimental function
  u32 Ins_;
  // {ref deletion length: number of top hit reads supporting it}
  // for experimental function
  unordered_map<u32, u32> Del_;
  // coverage depth of top hit reads
  // counted up even if a read supports ref deletion at the pos
  u32 Dp_;
  // coverage depth of all hit reads
  u32 AllDp_;
  u32 Pos_;
};

class RefInfo {
 public:
  RefInfo(double IndelRate) : kIndelRate_(IndelRate), Iden_(0), Cov_(0),
                              AllCov_(0), Dp_(0), HasFrameshift_(false) {}
  ~RefInfo() {}
  void Resize(u32 SeqLen);
  void EditPosIden(u32 Pos, double Iden) { PosInfo_.at(Pos).EditIden(Iden); }
  void EditPosDel(u32 Pos, double Len) { PosInfo_.at(Pos).EditDel(Len); }
  void EditPosIns(u32 Pos) { PosInfo_.at(Pos).EditIns(); }
  void EditPosDp(u32 Pos) { PosInfo_.at(Pos).EditDp(); }
  void EditPosAllDp(u32 Pos) { PosInfo_.at(Pos).EditAllDp(); }
  void SummarizeInfo();
  // never called (experimental function)
  void SetHasFrameshift();
  double Iden() const { return Iden_; }
  double Cov() const { return Cov_; }
  double AllCov() const { return AllCov_; }
  double Dp() const { return Dp_; }
  // never called (experimental function)
  bool HasFrameshift() const { return HasFrameshift_; }
  // never called (experimental function)
  string Frameshift() const { return Frameshift_; }
  // never used (for experimental function)
  const double kIndelRate_;

 private:
  // whether there are frameshifts due to short deletions
  // never called (experimental function)
  bool HasShortDelFrameshift();
  // whether there are frameshifts due to short insertions
  // never called (experimental function)
  bool HasShortInsFrameshift();
  // whether there are frameshifts due to long insertions
  // never called (experimental function)
  bool HasLongInsFrameshift();

  vector<RefPosInfo> PosInfo_;
  // average of RefPosInfo::AvgIden_ in region covered by top hit reads
  double Iden_;
  // ratio of region covered by top hit reads
  double Cov_;
  // ratio of region covered by all hit reads
  double AllCov_;
  // average coverage depth of top hit reads
  double Dp_;
  bool HasFrameshift_;
  string Frameshift_;
};

// each read pos has one ReadPosInfo object
class ReadPosInfo {
 public:
  ReadPosInfo() : As_(0), Iden_(0) {}
  ~ReadPosInfo() {}
  // initialize RefMatch_, RefDel_, RefIns_, and update As_, Iden_
  void InitTopHitInfo(u32 As, double Iden, const string& RefId);
  void EditRefMatch(const string& RefId, u32 RefPos) { RefMatch_[RefId].insert(RefPos); }
  void EditRefDel(const string& RefId, u32 RefPos) { ++RefDel_[RefId][RefPos]; }
  void EditRefIns(const string& RefId, u32 RefPos) { RefIns_[RefId].insert(RefPos); }
  void EditAllRefMatch(const string& RefId, u32 RefPos) { AllRefMatch_[RefId].insert(RefPos); }
  void EditAllRefIns(const string& RefId, u32 RefPos) { AllRefIns_[RefId].insert(RefPos); }
  // edit Fuga::RefInfoM_
  void EditRefInfoM(unordered_map<string, RefInfo>* RefInfoM);
  void ParseRefMatch(unordered_map<string, RefInfo>* RefInfoM);
  void ParseRefDel(unordered_map<string, RefInfo>* RefInfoM);
  void ParseRefIns(unordered_map<string, RefInfo>* RefInfoM);
  void ParseAllRefMatch(unordered_map<string, RefInfo>* RefInfoM);
  void ParseAllRefIns(unordered_map<string, RefInfo>* RefInfoM);
  u32 As() const { return As_; }

 private:
  // RefMatch_, RefDel_, RefIns_, As_, Iden_ are set considering only top hits
  // top hits are judged based on BWA's AS tag
  // sometimes there are multiple top hits for one read pos
  // AllRefMatch_, AllRefIns_ are set considering all hits
  // these variables stores corresponding ref pos, because they are set to edit
  // Fuga::RefInfoM_

  // {top hit ref ID: set<top hit ref pos>}
  unordered_map<string, unordered_set<u32> > RefMatch_;
  // {top hit ref ID, {top hit ref pos, ref deletion length on the right of it}
  unordered_map<string, unordered_map<u32, u32> > RefDel_;
  // {top hit ref ID: set<top hit ref pos deleted>]}
  unordered_map<string, unordered_set<u32> > RefIns_;
  // same as RefMatch_, but considering all hit reads, not only top hit reads
  unordered_map<string, unordered_set<u32> > AllRefMatch_;
  // same as RefIns_, but considering all hit reads, not only top hit reads
  unordered_map<string, unordered_set<u32> > AllRefIns_;
  // BWA's AS tag of top hit
  u32 As_;
  // identity of top hit
  double Iden_;
};

// each alignment info is stored in one SamL (SAM line) object
class SamL {
 public:
  SamL() : Flag_(0), LeftRefPos_(0), RightRefPos_(0), ReadLen_(0),
           AlignedRefLen_(0), ReadIsMapped_(true), ReadIsRC_(false), Nm_(0),
           As_(0), Iden_(0) {}
  ~SamL() {}
  void ReadL(const string& L, const vector<string>& LSp);
  // edit Fuga::EditRefInfoM()::ReadInfo
  void EditReadInfo(vector<ReadPosInfo>* ReadInfo) const;
  u32 Flag() const { return Flag_; }
  u32 LeftRefPos() const { return LeftRefPos_; }
  u32 RightRefPos() const { return RightRefPos_; }
  u32 ReadLen() const { return ReadLen_; }

 private:
  void ParseFlag();
  void ParseCigar(const string& L, const string& Cigar);
  void ParseTag(const string& L, const vector<string>& LSp);
  void MoveReadPos(u32* ReadPos) const;
  string RefId_;
  vector<char> Cigar_;
  u32 Flag_;
  u32 LeftRefPos_;
  u32 RightRefPos_;
  u32 ReadLen_;
  u32 AlignedRefLen_;
  bool ReadIsMapped_;
  bool ReadIsRC_;
  u32 Nm_;
  u32 As_;
  double Iden_;
};

#endif  // LIB_H_
