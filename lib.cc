# include "lib.h"

// read and parse a line in SAM
void SamL::ReadL(const string& L, const vector<string>& LSp) {
  Flag_ = stoul(LSp.at(1));
  ParseFlag();
  if (ReadIsMapped_) {
    RefId_ = LSp.at(2);
    LeftRefPos_ = stoul(LSp.at(3)) - 1;
    string Cigar = LSp.at(5);
    ParseCigar(L, Cigar);
    ParseTag(L, LSp);
  }
}

void SamL::ParseFlag() {
  if (Flag_ & 0x4) {
    ReadIsMapped_ = false;
  } else if (Flag_ & 0x10) {
    ReadIsRC_ = true;
  }
}

void SamL::ParseCigar(const string& L, const string& Cigar) {
  if (Cigar.find_first_not_of("0123456789ISDMH") != string::npos) {
    cerr << "ERROR: following line has unexpected CIGAR\n"
         << L << "\n";
    exit(1);
  }

  vector<string> CigarChar;
  vector<string> CigarIntStr;
  vector<u32> CigarInt;
  split(CigarChar, Cigar, "0123456789", true);
  split(CigarIntStr, Cigar, "ISDMH", true);
  u32 CigarBlockNum(CigarChar.size());
  CigarInt.resize(CigarBlockNum);
  for (u32 i = 0; i < CigarBlockNum; ++i) {
    CigarInt.at(i) = stoul(CigarIntStr.at(i));
  }

  u32 CigarLen(0);
  for (u32 i = 0; i != CigarBlockNum; ++i) {
    if (CigarChar.at(i) == "M" || CigarChar.at(i) == "D") {
      AlignedRefLen_ += CigarInt.at(i);
    }
    if (CigarChar.at(i) != "D") {
      ReadLen_ += CigarInt.at(i);
    }
    CigarLen += CigarInt.at(i);
  }
  RightRefPos_ = LeftRefPos_ + AlignedRefLen_ - 1;

  Cigar_.resize(CigarLen);
  u32 CigarIdx(0);  
  for (u32 i = 0; i != CigarBlockNum; ++i) {
    for (u32 j = 0; j != CigarInt.at(i); ++j) {
      Cigar_.at(CigarIdx) = CigarChar.at(i).at(0);
      ++CigarIdx;
    }
  }
}

void SamL::ParseTag(const string& L, const vector<string>& LSp) {
  if (LSp.at(11).substr(0, 2) != "NM" || LSp.at(13).substr(0, 2) != "AS") {
    cerr << "ERROR: following line has unexpected tag order\n"
         << L << "\n";
    exit(1);
  }
  Nm_ = stoul(LSp.at(11).substr(5));
  Iden_ = (AlignedRefLen_ - Nm_) / double(AlignedRefLen_);
  As_ = stoul(LSp.at(13).substr(5));
}

// scan a mapped read from left to right
// link each read pos to aligned ref pos
void SamL::EditReadInfo(vector<ReadPosInfo>* ReadInfo) const {
  u32 ReadPos(0);
  if (ReadIsRC_) {
    ReadPos = ReadLen_ - 1;
  }
  u32 RefPos(LeftRefPos_);
  for (auto i = Cigar_.begin(); i != Cigar_.end(); ++i) {
    if (As_ > ReadInfo->at(ReadPos).As()) {
      ReadInfo->at(ReadPos).InitTopHitInfo(As_, Iden_, RefId_);
    }
    if (*i == 'M') {
      if (As_ == ReadInfo->at(ReadPos).As()) {
        ReadInfo->at(ReadPos).EditRefMatch(RefId_, RefPos);
      }
      ReadInfo->at(ReadPos).EditAllRefMatch(RefId_, RefPos);
      MoveReadPos(&ReadPos);
      ++RefPos;
    } else if (*i == 'D') {
      if (As_ == ReadInfo->at(ReadPos).As()) {
        ReadInfo->at(ReadPos).EditRefIns(RefId_, RefPos);
      }
      ReadInfo->at(ReadPos).EditAllRefIns(RefId_, RefPos);
      ++RefPos;
    } else if (*i == 'I') {
      if (As_ == ReadInfo->at(ReadPos).As()) {
        ReadInfo->at(ReadPos).EditRefDel(RefId_, RefPos - 1);
      }
      MoveReadPos(&ReadPos);
    } else {
      MoveReadPos(&ReadPos);
    }
  }
}

// move current read pos considering the strand of the read
void SamL::MoveReadPos(u32* ReadPos) const {
  if (ReadIsRC_) {
    --(*ReadPos);
  } else {
    ++(*ReadPos);
  }
}

// called when a new hit with higher identity was found
void ReadPosInfo::InitTopHitInfo(u32 As, double Iden, const string& RefId) {
  RefMatch_.clear();
  RefDel_.clear();
  RefIns_.clear();
  As_ = As;
  Iden_ = Iden;
}

// edit info of ref pos corresponding to the read pos
// ex) identity of the alignment, coverage depth
// functions are called for each alignment type, but nothing happens in
// ParseRefDel() if the read pos does not support deletion of ref
void ReadPosInfo::EditRefInfoM(unordered_map<string, RefInfo>* RefInfoM) {
  ParseRefMatch(RefInfoM);
  // for experimental function
  ParseRefDel(RefInfoM);
  ParseRefIns(RefInfoM);
  ParseAllRefMatch(RefInfoM);
  ParseAllRefIns(RefInfoM);
}

void ReadPosInfo::ParseRefMatch(unordered_map<string, RefInfo>* RefInfoM) {
  for (auto i = RefMatch_.begin(); i != RefMatch_.end(); ++i) {
    const string& Id(i->first);
    const unordered_set<u32>& Pos(i->second);
    for (auto j = Pos.begin(); j != Pos.end(); ++j) {
      // append identity of the alignment (averaged afterwards)
      RefInfoM->at(Id).EditPosIden(*j, Iden_);
      RefInfoM->at(Id).EditPosDp(*j);
    }
  }
}

// for experimental function
void ReadPosInfo::ParseRefDel(unordered_map<string, RefInfo>* RefInfoM) {
  for (auto i = RefDel_.begin(); i != RefDel_.end(); ++i) {
    const string& Id(i->first);
    const unordered_map<u32, u32>& Del(i->second);
    for (auto j = Del.begin(); j != Del.end(); ++j) {
      u32 Pos(j->first);
      u32 Len(j->second);
      RefInfoM->at(Id).EditPosDel(Pos, Len);
      // EditPosIden() and EditPosDp() are not required here, because they are
      // called in ParseRefMatch() for this Pos
    }
  }
}

void ReadPosInfo::ParseRefIns(unordered_map<string, RefInfo>* RefInfoM) {
  for (auto i = RefIns_.begin(); i != RefIns_.end(); ++i) {
    const string& Id(i->first);
    const unordered_set<u32>& Pos(i->second);
    for (auto j = Pos.begin(); j != Pos.end(); ++j) {
      // append identity of the alignment (averaged afterwards)
      RefInfoM->at(Id).EditPosIden(*j, Iden_);
      RefInfoM->at(Id).EditPosDp(*j);
      RefInfoM->at(Id).EditPosIns(*j);
    }
  }
}

void ReadPosInfo::ParseAllRefMatch(unordered_map<string, RefInfo>* RefInfoM) {
  for (auto i = AllRefMatch_.begin(); i != AllRefMatch_.end(); ++i) {
    const string& Id(i->first);
    const unordered_set<u32>& Pos(i->second);
    for (auto j = Pos.begin(); j != Pos.end(); ++j) {
      RefInfoM->at(Id).EditPosAllDp(*j);
    }
  }
}

void ReadPosInfo::ParseAllRefIns(unordered_map<string, RefInfo>* RefInfoM) {
  for (auto i = AllRefIns_.begin(); i != AllRefIns_.end(); ++i) {
    const string& Id(i->first);
    const unordered_set<u32>& Pos(i->second);
    for (auto j = Pos.begin(); j != Pos.end(); ++j) {
      RefInfoM->at(Id).EditPosAllDp(*j);
    }
  }
}

void RefInfo::Resize(u32 SeqLen) {
  PosInfo_.resize(SeqLen);
  for (u32 i = 0; i != SeqLen; ++i) {
    PosInfo_.at(i).SetPos(i);
  }
}

// calculate average identity and coverage
void RefInfo::SummarizeInfo() {
  double IdenSum(0);
  // number of positions covered by top hit reads
  u32 CoveredPosNum(0);
  // number of positions covered by all hit reads
  u32 AllCoveredPosNum(0);
  u32 DpSum(0);
  for (auto i = PosInfo_.begin(); i != PosInfo_.end(); ++i) {
    if (i->Dp()) {
      i->SetAvgIden();
      IdenSum += i->AvgIden();
      DpSum += i->Dp();
      ++CoveredPosNum;
    }
    if (i->AllDp()) {
      ++AllCoveredPosNum;
    }
  }
  Cov_ = double(CoveredPosNum) / PosInfo_.size();
  AllCov_ = double(AllCoveredPosNum) / PosInfo_.size();
  if (CoveredPosNum) {
    Iden_ = IdenSum / CoveredPosNum;
    Dp_ = double(DpSum) / CoveredPosNum;
  } else {
    Iden_ = 0;
    Dp_ = 0;
  }
}

// never called (experimental function)
void RefInfo::SetHasFrameshift() {
  if (HasShortDelFrameshift() || HasShortInsFrameshift()) {
    HasFrameshift_ = true;
  }
}

// never called (experimental function)
bool RefInfo::HasShortDelFrameshift() {
  for (auto i = PosInfo_.begin(); i != PosInfo_.end(); ++i) {
    u32 Dp(i->Dp());
    const unordered_map<u32, u32>& Del(i->Del());
    for (auto j = Del.begin(); j != Del.end(); ++j) {
      u32 Len(j->first);
      u32 ReadNum(j->second);
      if (Len % 3 != 0 && double(ReadNum) / Dp >= kIndelRate_) {
        u32 Pos(i->Pos());
        // read側の挿入
        Frameshift_ = to_string((long long int) Len) + " bp insertion between " +
                      to_string((long long int) Pos + 1) + " and " + to_string((long long int) Pos + 2);
        return true;
      }
    }
  }
  return false;
}

// never called (experimental function)
bool RefInfo::HasShortInsFrameshift() {
  u32 InsBeg(0);
  u32 ContinuousIns(0);
  for (auto i = PosInfo_.begin(); i != PosInfo_.end(); ++i) {
    u32 Dp(i->Dp());
    u32 Ins(i->Ins());
    u32 Pos(i->Pos());
    if (double(Ins) / Dp >= kIndelRate_) {
      if (ContinuousIns == 0) {
        InsBeg = Pos;
      }
      ++ContinuousIns;
    } else {
      if (ContinuousIns % 3 != 0) {
        u32 InsEnd(Pos);
        // read側の欠失
        Frameshift_ = "deletion [" + to_string((long long int) InsBeg + 1) + "," +
                      to_string((long long int) InsEnd) + "]";
        return true;
      }
      ContinuousIns = 0;
    }
  }
  return false;
}

// never called (experimental function)
bool RefInfo::HasLongInsFrameshift() {
  u32 InsBeg(0);
  u32 ContinuousDp0(0);
  for (auto i = PosInfo_.begin(); i != PosInfo_.end(); ++i) {
    u32 Dp(i->Dp());
    u32 Pos(i->Pos());
    if (Dp == 0) {
      if (ContinuousDp0 == 0) {
        InsBeg = Pos;
      }
      ++ContinuousDp0;
    } else {
      if (ContinuousDp0 % 3 != 0) {
        u32 InsEnd(Pos);
        // read側の欠失
        Frameshift_ = "deletion [" + to_string((long long int) InsBeg + 1) + "," +
                      to_string((long long int) InsEnd) + "]";
        return true;
      } else {
        ContinuousDp0 = 0;
      }
    }
  }
  return false;
}
