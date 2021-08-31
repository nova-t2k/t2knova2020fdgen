#pragma once

#include "TDirectory.h"
#include "TFile.h"
#include "TH3.h"
#include "THn.h"

#include <iostream>
#include <unordered_map>

namespace t2knova {

struct FlagBlob {
  FlagBlob()
      : Mode(0), flagCCINC(false), flagCC0pi(false), flagCC1cpi(false),
        flagCC1pi0(false), flagNCINC(false), flagNC0pi(false),
        flagNC1cpi(false), flagNC1pi0(false) {}
  int Mode;
  bool flagCCINC;
  bool flagCC0pi;
  bool flagCC1cpi;
  bool flagCC1pi0;
  bool flagNCINC;
  bool flagNC0pi;
  bool flagNC1cpi;
  bool flagNC1pi0;
};

inline FlagBlob GetFlagBlob(bool iscc, int NFSCpi, int NFSpi0, int NFSOther,
                            int Mode = 0) {

  FlagBlob fb;
  fb.Mode = Mode;

  fb.flagCCINC = iscc;
  fb.flagNCINC = !iscc;

  if (NFSOther == 0) { // Normalish event with just nucleons and pions
    if ((NFSpi0 + NFSCpi) == 0) {
      fb.flagCC0pi = iscc;
      fb.flagNC0pi = !iscc;
    } else if ((NFSpi0 == 0) && ((NFSCpi) == 1)) {
      fb.flagCC1cpi = iscc;
      fb.flagNC1cpi = !iscc;
    } else if ((NFSpi0 == 1) && ((NFSCpi) == 0)) {
      fb.flagCC1pi0 = iscc;
      fb.flagNC1pi0 = !iscc;
    } // else CCOther
  }   // else CCOther

  return fb;
}

template <typename TH> struct THTraits {};

template <> struct THTraits<TH1D> { using Base = TH1; };
template <> struct THTraits<TH2D> { using Base = TH1; };
template <> struct THTraits<TH3D> { using Base = TH1; };

template <> struct THTraits<THnD> { using Base = THnBase; };

template <typename TH> struct modeblob {
  TH TopoHist;
  std::unordered_map<int, TH> ModeHists;

  std::string fName;

  using THBase = typename THTraits<TH>::Base;

  modeblob() {}

  template <typename... Args>
  modeblob(std::string const &name, std::string const &title, Args... binning)
      : fName(name), TopoHist(name.c_str(), title.c_str(), binning...) {

    TopoHist.SetDirectory(nullptr);
    TopoHist.Sumw2(true);
  }

  void SetName(std::string const &name) {
    fName = name;

    TopoHist.SetName(fName.c_str());

    for (auto &h : ModeHists) {
      h.second.SetName((fName + "_Mode_" + (h.first < 0 ? "m" : "") +
                        std::to_string(std::abs(h.first)))
                           .c_str());
    }
  }
  void SetTitle(std::string const &title) {
    Apply([=](TH &h) { h.SetTitle(title.c_str()); });
  }

  void Write(TDirectory *f, bool scale = false) {
    Apply([=](TH &h) {
      if (!h.Integral()) {
        return;
      }

      if (scale) {
        THBase *hc = (THBase *)h.Clone();
        hc->SetDirectory(nullptr);
        hc->Scale(1, "width");
        f->WriteTObject(hc, h.GetName());
        delete hc;
      } else {
        f->WriteTObject(&h, h.GetName());
      }
    });
  }

  template <typename... XY> void Fill(double w, int mode, XY... xy) {
    TopoHist.Fill(xy..., w);

    if (mode != 0) {
      if (!ModeHists.count(mode)) {
        ModeHists.emplace(mode, TopoHist);

        ModeHists[mode].Reset();
        ModeHists[mode].SetName((fName + "_Mode_" + (mode < 0 ? "m" : "") +
                                 std::to_string(std::abs(mode)))
                                    .c_str());

        ModeHists[mode].SetDirectory(nullptr);
      }
      ModeHists[mode].Fill(xy..., w);
    }
  }

  template <typename THOut>
  modeblob<THOut> Transform(std::function<THOut(TH const &)> f) {
    modeblob<THOut> out;
    out.TopoHist = THOut(f(TopoHist));

    for (auto const &h : ModeHists) {
      out.ModeHists[h.first] = THOut(f(h.second));
    }
    return out;
  }

  void Apply(std::function<void(TH &)> f) {
    f(TopoHist);
    for (auto &h : ModeHists) {
      f(h.second);
    }
  }
};

template <typename TH> struct hblob {
  modeblob<TH> CCInc;
  modeblob<TH> CC0Pi;
  modeblob<TH> CC1CPi;
  modeblob<TH> CC1Pi0;
  modeblob<TH> CCOther;
  modeblob<TH> NCInc;
  modeblob<TH> NC0Pi;
  modeblob<TH> NC1CPi;
  modeblob<TH> NC1Pi0;
  modeblob<TH> NCOther;

  std::string fName;

  using THBase = typename THTraits<TH>::Base;

  hblob() {}

  template <typename... Args>
  hblob(std::string const &name, std::string const &title, Args... binning)
      : fName(name), CCInc(name + "_CCInc", title, binning...),
        CC0Pi(name + "_CC0Pi", title, binning...),
        CC1CPi(name + "_CC1CPi", title, binning...),
        CC1Pi0(name + "_CC1Pi0", title, binning...),
        CCOther(name + "_CCOther", title, binning...),
        NCInc(name + "_NCInc", title, binning...),
        NC0Pi(name + "_NC0Pi", title, binning...),
        NC1CPi(name + "_NC1CPi", title, binning...),
        NC1Pi0(name + "_NC1Pi0", title, binning...),
        NCOther(name + "_NCOther", title, binning...) {
  }

  void SetName(std::string const &name) {
    fName = name;

    CCInc.SetName(fName + "_CCInc");
    CC0Pi.SetName(fName + "_CC0Pi");
    CC1CPi.SetName(fName + "_CC1CPi");
    CC1Pi0.SetName(fName + "_CC1Pi0");
    CCOther.SetName(fName + "_CCOther");
    NCInc.SetName(fName + "_NCInc");
    NC0Pi.SetName(fName + "_NC0Pi");
    NC1CPi.SetName(fName + "_NC1CPi");
    NC1Pi0.SetName(fName + "_NC1Pi0");
    NCOther.SetName(fName + "_NCOther");
  }

  void SetTitle(std::string const &title) {
    Apply([=](TH &h) { h.SetTitle(title); });
  }

  void Write(TDirectory *f, bool scale = false) {
    CCInc.Write(f, scale);
    CC0Pi.Write(f, scale);
    CC1CPi.Write(f, scale);
    CC1Pi0.Write(f, scale);
    CCOther.Write(f, scale);
    NCInc.Write(f, scale);
    NC0Pi.Write(f, scale);
    NC1CPi.Write(f, scale);
    NC1Pi0.Write(f, scale);
    NCOther.Write(f, scale);
  }

  template <typename... XY> void Fill(double w, FlagBlob const &blb, XY... xy) {
    if (blb.flagCCINC) {

      CCInc.Fill(w, blb.Mode, xy...);

      if (blb.flagCC0pi) {
        CC0Pi.Fill(w, blb.Mode, xy...);
      } else if (blb.flagCC1cpi) {
        CC1CPi.Fill(w, blb.Mode, xy...);
      } else if (blb.flagCC1pi0) {
        CC1Pi0.Fill(w, blb.Mode, xy...);
      } else {
        CCOther.Fill(w, blb.Mode, xy...);
      }
    } else if (blb.flagNCINC) {

      NCInc.Fill(w, blb.Mode, xy...);

      if (blb.flagNC0pi) {
        NC0Pi.Fill(w, blb.Mode, xy...);
      } else if (blb.flagNC1cpi) {
        NC1CPi.Fill(w, blb.Mode, xy...);
      } else if (blb.flagNC1pi0) {
        NC1Pi0.Fill(w, blb.Mode, xy...);
      } else {
        NCOther.Fill(w, blb.Mode, xy...);
      }
    }
  }

  template <typename THOut>
  hblob<THOut> Transform(std::function<THOut(TH const &)> f) {
    hblob<THOut> out;
    out.CCInc = CCInc.Transform(f);
    out.CC0Pi = CC0Pi.Transform(f);
    out.CC1CPi = CC1CPi.Transform(f);
    out.CC1Pi0 = CC1Pi0.Transform(f);
    out.CCOther = CCOther.Transform(f);
    out.NCInc = NCInc.Transform(f);
    out.NC0Pi = NC0Pi.Transform(f);
    out.NC1CPi = NC1CPi.Transform(f);
    out.NC1Pi0 = NC1Pi0.Transform(f);
    out.NCOther = NCOther.Transform(f);
    return out;
  }

  void Apply(std::function<void(TH &)> f) {
    CCInc.Apply(f);
    CC0Pi.Apply(f);
    CC1CPi.Apply(f);
    CC1Pi0.Apply(f);
    CCOther.Apply(f);
    NCInc.Apply(f);
    NC0Pi.Apply(f);
    NC1CPi.Apply(f);
    NC1Pi0.Apply(f);
    NCOther.Apply(f);
  }
};

} // namespace t2knova