#pragma once

#include <iostream>
#include <sstream>
#include <type_traits>
#include <unordered_map>

#include "T2KNOvAFakeDataHelper.hxx"
#include "TDirectory.h"
#include "TFile.h"
#include "TH3.h"
#include "THn.h"
#include "toml/toml_helper.h"

namespace t2knova {

template <typename TH>
struct TrueChannelHist {
  std::unordered_map<int, TH> Hists;
  std::string fName;
  using THBase = typename THTraits<TH>::Base;

  const int kInclusive = 0;

  TrueChannelHist() {}

  template <typename... Args>
  TrueChannelHist(std::string const &name, std::string const &title,
                  Args... binning)
      : fName(name) {
    Hists.emplace(kInclusive, TH(name.c_str(), title.c_str(), binning...));
    Hists[kInclusive].SetDirectory(nullptr);
    Hists[kInclusive].Sumw2(true);
  }

  void SetName(std::string const &name) {
    fName = name;

    for (auto &h : Hists) {
      if (h.first == kInclusive) {
        h.second.SetName(fName.c_str());
      } else {
        h.second.SetName((fName + "_Mode_" + (h.first < 0 ? "m" : "") +
                          std::to_string(std::abs(h.first)))
                             .c_str());
      }
    }
  }

  void Apply(std::function<void(TH &)> f) {
    for (auto &h : Hists) {
      f(h.second);
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

  void AddModeHist(int TrueChannel) {
    if (!Hists.count(TrueChannel)) {
      Hists.emplace(TrueChannel, Hists[kInclusive]);

      Hists[TrueChannel].Reset();
      Hists[TrueChannel].SetName((fName + "_Mode_" +
                                  (TrueChannel < 0 ? "m" : "") +
                                  std::to_string(std::abs(TrueChannel)))
                                     .c_str());

      Hists[TrueChannel].SetDirectory(nullptr);
    }
  }

  template <typename... XY>
  void Fill(double w, int TrueChannel, XY... xy) {
    Hists[kInclusive].Fill(xy..., w);

    if (TrueChannel != kInclusive) {
      AddModeHist(TrueChannel);
      Hists[TrueChannel].Fill(xy..., w);
    }
  }

  template <typename THOut>
  TrueChannelHist<THOut> Transform(std::function<THOut(TH const &)> f) {
    TrueChannelHist<THOut> out;
    for (auto const &h : Hists) {
      out.Hists[h.first] = THOut(f(h.second));
    }
    return out;
  }
};

template <typename TH>
struct SelectionHists {
  std::vector<std::string> fSelections;
  std::vector<TrueChannelHist<TH>> Hists;
  std::string fName;
  using THBase = typename THTraits<TH>::Base;

  SelectionHists() {}

  template <typename... Args>
  SelectionHists(std::string const &name, std::string const &title,
                 std::vector<std::string> const &Selections, Args... binning)
      : fName(name), fSelections(Selections) {
    for (int i = 0; i < fSelections.size(); ++i) {
      Hists.emplace_back((fName + "_" + fSelections[i]).c_str(), title.c_str(),
                         binning...);
    }
  }

  void SetName(std::string const &name) {
    fName = name;

    for (int i = 0; i < fSelections.size(); ++i) {
      Hists[i].SetName((fName + "_" + fSelections[i]).c_str());
    }
  }

  void SetTitle(std::string const &title) {
    Apply([=](TH &h) { h.SetTitle(title); });
  }

  void Apply(std::function<void(TH &)> f) {
    for (auto &h : Hists) {
      h.Apply(f);
    }
  }

  void Write(TDirectory *f, bool scale = false) {
    for (auto &h : Hists) {
      h.Write(f, scale);
    }
  }

  template <typename... XY>
  void Fill(double w, std::vector<int> const &Selections, int TrueChannel,
            XY... xy) {
    for (int sel : Selections) {
      if (sel < Hists.size()) {
        Hists[sel].Fill(w, TrueChannel, xy...);
      }
    }
  }

  template <typename THOut>
  SelectionHists<THOut> Transform(std::function<THOut(TH const &)> f) {
    SelectionHists<THOut> out;
    for (auto &h : Hists) {
      out.Hists.emplace_back(h.Transform(f));
    }
    return out;
  }
};

std::vector<double> BinningFromTOML(toml::value const &binning_config) {
  std::string binning_type = toml_h::find<std::string>(binning_config, 0);
  std::vector<double> const &binning_details =
      toml_h::find<std::vector<double>>(binning_config, 1);

  std::vector<double> binning;
  if (binning_type == "uniform") {
    for (int i = 0; i < (binning_details[0] + 1); ++i) {
      binning.push_back(binning_details[1] +
                        double(i) * ((binning_details[2] - binning_details[1]) /
                                     binning_details[0]));
    }
  } else {
    binning = binning_details;
  }

  return binning;
}

template <typename TN>
SelectionHists<TN> *SelectionHistsFromTOML(
    typename std::enable_if<THTraits<TN>::NDims == 1, std::string>::type const
        &plotname,
    toml::value const &config) {
  std::vector<std::string> axis_titles =
      toml_h::find<std::vector<std::string>>(config, "axis_titles");

  while (axis_titles.size() < 2) {
    axis_titles.push_back("");
  }

  std::string titlestr = toml_h::find<std::string>(config, "title") + ";" +
                         axis_titles[0] + ";" + axis_titles[1];

  std::vector<double> x_binning =
      BinningFromTOML(toml_h::find(config, "x_binning"));

  return new SelectionHists<TN>(
      plotname, titlestr,
      toml_h::find<std::vector<std::string>>(config, "selections"),
      x_binning.size() - 1, x_binning.data());
}

template <typename TN>
SelectionHists<TN> *SelectionHistsFromTOML(
    typename std::enable_if<THTraits<TN>::NDims == 2, std::string>::type const
        &plotname,
    toml::value const &config) {
  std::vector<std::string> axis_titles =
      toml_h::find<std::vector<std::string>>(config, "axis_titles");

  while (axis_titles.size() < 3) {
    axis_titles.push_back("");
  }

  std::string titlestr = toml_h::find<std::string>(config, "title") + ";" +
                         axis_titles[0] + ";" + axis_titles[1] + ";" +
                         axis_titles[2];

  std::vector<double> x_binning =
      BinningFromTOML(toml_h::find(config, "x_binning"));
  std::vector<double> y_binning =
      BinningFromTOML(toml_h::find(config, "y_binning"));

  return new SelectionHists<TN>(
      plotname, titlestr,
      toml_h::find<std::vector<std::string>>(config, "selections"),
      x_binning.size() - 1, x_binning.data(), y_binning.size() - 1,
      y_binning.data());
}

template <typename TN>
SelectionHists<TN> *SelectionHistsFromTOML(
    typename std::enable_if<THTraits<TN>::NDims == 3, std::string>::type const
        &plotname,
    toml::value const &config) {
  std::vector<std::string> axis_titles =
      toml_h::find<std::vector<std::string>>(config, "axis_titles");

  while (axis_titles.size() < 4) {
    axis_titles.push_back("");
  }

  std::string titlestr = toml_h::find<std::string>(config, "title") + ";" +
                         axis_titles[0] + ";" + axis_titles[1] + ";" +
                         axis_titles[2] + ";" + axis_titles[3];

  std::vector<double> x_binning =
      BinningFromTOML(toml_h::find(config, "x_binning"));
  std::vector<double> y_binning =
      BinningFromTOML(toml_h::find(config, "y_binning"));
  std::vector<double> z_binning =
      BinningFromTOML(toml_h::find(config, "z_binning"));

  return new SelectionHists<TN>(
      plotname, titlestr,
      toml_h::find<std::vector<std::string>>(config, "selections"),
      x_binning.size() - 1, x_binning.data(), y_binning.size() - 1,
      y_binning.data(), z_binning.size() - 1, z_binning.data());
}

}  // namespace t2knova