#include "type.h"

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

// initialize object storing ref info
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

// parse and store alignment info for each line in SAM
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

// link ref pos and mapped read pos
void Fuga::EditRefInfoM(unordered_map<string, vector<SamL> >* Sam) {
  // for each read ID
  for (auto i = Sam->begin(); i != Sam->end(); ++i) {
    const vector<SamL>& SamLV(i->second);
    u32 ReadLen(SamLV.begin()->ReadLen());
    vector<ReadPosInfo> ReadInfo(ReadLen);
    // for each alignment
    for (auto j = SamLV.begin(); j != SamLV.end(); ++j) {
      // edit info of each read position
      // ex) corresponding ref position
      j->EditReadInfo(&ReadInfo);
    }
    // for each read pos j
    for (auto j = ReadInfo.begin(); j != ReadInfo.end(); ++j) {
      // edit info of each ref position corresponding to j
      // ex) identity of the alignment, coverage depth
      j->EditRefInfoM(&RefInfoM_);
    }
  }
}

// calculate average identity and coverage for each ref ID
void Fuga::SummarizeRefInfo() {
  // for each ref ID
  for (auto i = RefInfoM_.begin(); i != RefInfoM_.end(); ++i) {
    RefInfo* Info(&i->second);
    // calculate average identity and coverage
    Info->SummarizeInfo();
    // always false (for experimental function)
    if (kJudgeFrameshift_) {
      Info->SetHasFrameshift();
    }
  }
}

// calculate score for each ref ID to sort ref IDs by it
void Fuga::SetRefScore() {
  // for each ref ID
  for (auto i = RefInfoM_.begin(); i != RefInfoM_.end(); ++i) {
    string Id(i->first);
    RefInfo* Info(&i->second);
    // calculate and store score (identity * coverage)
    RefScore_[Info->Iden() * Info->Cov()].emplace_back(Id);
  }
}

// print identity and coverage for each ref ID
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
