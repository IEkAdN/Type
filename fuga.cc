#include "fuga.h"

void Fuga::Main() {
  ReadRef();
  ReadSam(kR1SamNom_, &R1Sam_);
  ReadSam(kR2SamNom_, &R2Sam_);
  EditRefInfoM(&R1Sam_);
  EditRefInfoM(&R2Sam_);
  SummarizeRefInfo();
  SetRefScore();
  PrintHit();
}

void Fuga::ReadRef() {
  ifstream F(kRefFaNom_);
  string L;
  string Id;
  u32 SeqLen(0);
  bool Is1stL(true);
  while (getline(F, L)) {
    if (! L.empty()) {
      if (L.at(0) == '>') {
        if (Is1stL) {
          Is1stL = false;
        } else {
          RefInfoM_.insert(make_pair(Id, kIndelRate_));
          RefInfoM_.at(Id).Resize(SeqLen);
          SeqLen = 0;
        }
        Id = L.substr(1, L.find(" ") - 1);
      } else {
        SeqLen += L.size();
      }
    }
  }
  RefInfoM_.insert(make_pair(Id, kIndelRate_));
  RefInfoM_.at(Id).Resize(SeqLen);
}

void Fuga::ReadSam(string FNom, unordered_map<string, vector<SamL> >* Sam) {
  ifstream F(FNom);
  string L;
  while (getline(F, L)) {
    if (L.at(0) != '@') {
      vector<string> LSp;
      split(LSp, L, "\t");
      string ReadId(LSp.at(0));
      (*Sam)[ReadId].emplace_back();
      Sam->at(ReadId).back().ReadL(L, LSp);
    }
  }
}

void Fuga::EditRefInfoM(unordered_map<string, vector<SamL> >* Sam) {
  for (auto i = Sam->begin(); i != Sam->end(); ++i) {
    // 各refへのtop hit readを出力したい場合用
    // const string& ReadId(i->first);
    const vector<SamL>& SamLV(i->second);
    u32 ReadLen(SamLV.begin()->ReadLen());
    vector<ReadPosInfo> ReadInfo(ReadLen);
    for (auto j = SamLV.begin(); j != SamLV.end(); ++j) {
      j->EditReadInfo(&ReadInfo);
    }
    for (auto j = ReadInfo.begin(); j != ReadInfo.end(); ++j) {
      j->EditRefInfoM(&RefInfoM_);
      // 各refへのtop hit readを出力したい場合用
      // j->EditRefInfoM(&RefInfoM_, ReadId);
    }
  }
}

void Fuga::SummarizeRefInfo() {
  for (auto i = RefInfoM_.begin(); i != RefInfoM_.end(); ++i) {
    RefInfo* Info(&i->second);
    Info->SummarizeInfo();
    if (kJudgeFrameshift_) {
      Info->SetHasFrameshift();
    }
  }
}

void Fuga::SetRefScore() {
  for (auto i = RefInfoM_.begin(); i != RefInfoM_.end(); ++i) {
    string Id(i->first);
    RefInfo* Info(&i->second);
    RefScore_[Info->Iden() * Info->Cov()].emplace_back(Id);
  }
}

void Fuga::PrintHit() {
  for (auto i = RefScore_.begin(); i != RefScore_.end(); ++i) {
    const vector<string>& IdV(i->second);
    for (auto j = IdV.begin(); j != IdV.end(); ++j) {
      const RefInfo& _RefInfo(RefInfoM_.at(*j));
      cout << *j << "\t" <<_RefInfo.Iden() << "\t" << _RefInfo.Cov() << "\t"
           << _RefInfo.AllCov() << "\n";
    }
  }
}
