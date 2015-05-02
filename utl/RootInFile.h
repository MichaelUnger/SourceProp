#ifndef _RootInFile_h_
#define _RootInFile_h_

/*
  GPLv2 and 2C-BSD
  Copyright (c) Darko Veberic, 2014
*/

#include "SaveCurrentTDirectory.h"
#include <TChain.h>
#include <TChainElement.h>
#include <TFile.h>
#include <string>
#include <stdexcept>


template<class Entry>
class RootInFile {
public:
  class Iterator {
  public:
    typedef Iterator self_type;
    typedef Entry value_type;
    typedef Entry& reference;
    typedef Entry* pointer;
    typedef int difference_type;
    typedef std::forward_iterator_tag iterator_category;

    Iterator(RootInFile& file, unsigned int index) : fFile(file), fIndex(index) { }
    Iterator& operator++() { ++fIndex; return *this; } // prefix ++it
    const Entry& operator*() { return fFile[fIndex]; }
    const Entry* operator->() { return &operator*(); }
    bool operator==(const Iterator& it) const { IsSameFile(it); return fIndex == it.fIndex; }
    bool operator!=(const Iterator& it) const { IsSameFile(it); return !operator==(it); }
  private:
    void
    IsSameFile(const Iterator& it)
      const
    {
      if (&fFile != &it.fFile)
        throw std::logic_error("mixing iterators of different RootFiles");
    }

    RootInFile& fFile;
    unsigned int fIndex;
  };

  RootInFile(const std::string& filename) :
    fChain(0),
    fEntryBuffer(0)
  {
    Open(filename);
  }

  ~RootInFile() { Close(); }

  unsigned int
  GetSize()
  {
    Check();
    return fChain->GetEntries();
  }

  const Entry&
  operator[](const unsigned int i)
  {
    Check();
    if (!fChain->GetEntry(i))
      throw std::out_of_range("no entry found");
    return *fEntryBuffer;
  }

  Iterator
  Begin()
  {
    Check();
    return Iterator(*this, 0);
  }

  Iterator
  End()
  {
    Check();
    return Iterator(*this, GetSize());
  }

  template<class T>
  const T&
  Get()
  {
    TObjArray* files;
    if (!fChain || !(files = fChain->GetListOfFiles()))
      Error("file not open");
    const T* obj = 0;
    TIter next(files);
    for (TChainElement* c = (TChainElement*)next(); !obj && c; c = (TChainElement*)next())
      TFile(c->GetTitle()).GetObject(T::Class_Name(), obj);
    if (!obj)
      Error((std::string("no object '") + T::Class_Name() + "' found in file").c_str());
    return *obj;
  }

  void
  Close()
  {
    delete fEntryBuffer;
    fEntryBuffer = 0;
    delete fChain;
    fChain = 0;
  }

private:
  // prevent copying
  RootInFile(const RootInFile&);
  RootInFile& operator=(const RootInFile&);

  void
  Error(const char* const message)
  {
    Close();
    throw std::runtime_error(message);
  }

  void
  Open(const std::string& filename)
  {
    const SaveCurrentTDirectory save;
    fChain = new TChain((std::string(Entry::Class_Name()) + "Tree").c_str());
    const int nFiles = fChain->Add(filename.c_str());
    if (!nFiles)
      Error("no files found");
    fEntryBuffer = new Entry;
    fChain->SetBranchAddress(Entry::Class_Name(), &fEntryBuffer);
    Check();
  }

  void
  Check()
  {
    if (!fChain || !fEntryBuffer)
      Error("no entry found");
  }

  TChain* fChain;
  Entry* fEntryBuffer;
};


#endif
