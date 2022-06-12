#ifndef LIB_H_
#define LIB_H_

#include "hoge.h"

// 各ref posに対しRefPosInfo型のobject 1つが対応
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
  // [当該ref posへのtop hitのiden]
  // あるtop hitにより当該posの欠失が支持されている場合も含む
  vector<double> Iden_;
  double AvgIden_;
  // ref側での挿入を支持するtop hitの数
  u32 Ins_;
  // {当該ref posと当該ref pos + 1の間にあるref側の欠失塩基数,
  //  その欠失を支持するtop hitの数}
  unordered_map<u32, u32> Del_;
  // top hitによるcoverage depth
  // あるtop hitにより当該posの欠失が支持されている場合もカウント
  u32 Dp_;
  // 全hitによるcoverage depth
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
  // 各refへのtop hit readを出力したい場合用
  // void EditTopHitReadId(string ReadId) { TopHitReadId_.insert(ReadId); }
  // void PrintTopHitReadId() { for (auto i = TopHitReadId_.begin(); i != TopHitReadId_.end(); ++i) { cout << *i << "\n"; } }
  void SummarizeInfo();
  void SetHasFrameshift();
  double Iden() const { return Iden_; }
  double Cov() const { return Cov_; }
  double AllCov() const { return AllCov_; }
  double Dp() const { return Dp_; }
  bool HasFrameshift() const { return HasFrameshift_; }
  string Frameshift() const { return Frameshift_; }
  const double kIndelRate_;

 private:
  // alignmentがまたげる程度の長さのref側欠失(CIGARがIになる)による
  // frameshiftがあるか
  bool HasShortDelFrameshift();
  // alignmentがまたげる程度の長さのref側挿入(CIGARがDになる)による
  // frameshiftがあるか
  bool HasShortInsFrameshift();
  // alignmentがまたげない程度の長さのref側挿入(CIGARがSになる)による
  // frameshiftがあるか
  bool HasLongInsFrameshift();
  // alignmentがまたげない程度の長さのref側欠失(CIGARがSになる)による
  // frameshiftがあるかは判定困難(挿入の長さがSの長さ以上になる可能性があるため)

  vector<RefPosInfo> PosInfo_;
  // 各refへのtop hit readを出力したい場合用
  // unordered_set<string> TopHitReadId_;
  // top hitのidentityの平均の平均(RefPosInfo::AvgIden_の全pos平均)
  double Iden_;
  // top hitによるカバー率
  double Cov_;
  // 全hitによるカバー率
  double AllCov_;
  // top hitによるcoverage depthの平均
  double Dp_;
  bool HasFrameshift_;
  string Frameshift_;
};

// 各read posに対しReadPosInfo型のobject 1つが対応
class ReadPosInfo {
 public:
  ReadPosInfo() : As_(0), Iden_(0) {}
  ~ReadPosInfo() {}
  // RefMatch_, RefDel_, RefIns_を初期化, As_, Iden_を更新
  void InitTopHitInfo(u32 As, double Iden, const string& RefId);
  void EditRefMatch(const string& RefId, u32 RefPos) { RefMatch_[RefId].insert(RefPos); }
  void EditRefDel(const string& RefId, u32 RefPos) { ++RefDel_[RefId][RefPos]; }
  void EditRefIns(const string& RefId, u32 RefPos) { RefIns_[RefId].insert(RefPos); }
  void EditAllRefMatch(const string& RefId, u32 RefPos) { AllRefMatch_[RefId].insert(RefPos); }
  void EditAllRefIns(const string& RefId, u32 RefPos) { AllRefIns_[RefId].insert(RefPos); }
  // Fuga::RefInfoM_を編集
  void EditRefInfoM(unordered_map<string, RefInfo>* RefInfoM);
  // 各refへのtop hit readを出力したい場合用
  // void EditRefInfoM(unordered_map<string, RefInfo>* RefInfoM, const string& ReadId);
  void ParseRefMatch(unordered_map<string, RefInfo>* RefInfoM);
  void ParseRefDel(unordered_map<string, RefInfo>* RefInfoM);
  void ParseRefIns(unordered_map<string, RefInfo>* RefInfoM);
  void ParseAllRefMatch(unordered_map<string, RefInfo>* RefInfoM);
  void ParseAllRefIns(unordered_map<string, RefInfo>* RefInfoM);
  u32 As() const { return As_; }

 private:
  // RefMatch_, RefDel_, RefIns_, As_, Iden_にはtop hit
  // (AS基準, 複数の可能性あり)の情報を記憶
  // AllRefMatch_, AllRefIns_には全hitの情報を記憶
  // Fuga::RefInfoM_の編集が目的なのでref posを基準に記憶

  // {top hit ref ID, [top hit ref pos]}
  // split alignmentにより1 readの複数alignmentが同じマッチを支持している場合に
  // ダブルカウントしないようにvectorではなくunordered_set
  unordered_map<string, unordered_set<u32> > RefMatch_;
  // {top hit ref ID, {top hit ref pos, そのref posの右側に何文字欠失があるか}
  unordered_map<string, unordered_map<u32, u32> > RefDel_;
  // {top hit ref ID, [ref側で欠失しているtop hit ref pos]}
  // split alignmentにより1 readの複数alignmentが同じ欠失を支持している場合に
  // ダブルカウントしないようにvectorではなくunordered_set
  unordered_map<string, unordered_set<u32> > RefIns_;
  // RefMatch_のtop hitに限定しないversion
  unordered_map<string, unordered_set<u32> > AllRefMatch_;
  // RefIns_のtop hitに限定しないversion
  unordered_map<string, unordered_set<u32> > AllRefIns_;
  // top hitのAS
  u32 As_;
  // top hitのidentity
  double Iden_;
};

class SamL {
 public:
  SamL() : Flag_(0), LeftRefPos_(0), RightRefPos_(0), ReadLen_(0),
           AlignedRefLen_(0), ReadIsMapped_(true), ReadIsRC_(false), Nm_(0),
           As_(0), Iden_(0) {}
  ~SamL() {}
  void ReadL(const string& L, const vector<string>& LSp);
  // Fuga::EditRefInfoM()::ReadInfoを編集
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
