#pragma once

#include "T2KNOvA/FakeDataHelper.hxx"

#include "plotutils.h"

#ifndef NO_TOML
#include "toml/toml_helper.h"
#endif

#include "TDirectory.h"
#include "TFile.h"
#include "TH3.h"
#include "THn.h"

#include <iostream>
#include <sstream>
#include <type_traits>
#include <unordered_map>

namespace t2knova {

template <typename TH> struct TrueChannelHist {
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

  void SetXAxisTitle(std::string const &title) {
    Apply([=](TH &h) { h.GetXaxis()->SetTitle(title.c_str()); });
  }

  void SetYAxisTitle(std::string const &title) {
    Apply([=](TH &h) { h.GetYaxis()->SetTitle(title.c_str()); });
  }

  void SetZAxisTitle(std::string const &title) {
    Apply([=](TH &h) { h.GetZaxis()->SetTitle(title.c_str()); });
  }

  void Write(TDirectory *f, bool width_scale = false) {
    Apply([=](TH &h) {
      if (!h.Integral()) {
        return;
      }

      if (width_scale) {
        THBase *hc = (THBase *)h.Clone();
        hc->SetDirectory(nullptr);
        Scale(hc, 1.0, width_scale);
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

  template <typename... XY> void Fill(double w, int TrueChannel, XY... xy) {
    Hists[kInclusive].Fill(xy..., w);

    if (TrueChannel != kInclusive) {
      AddModeHist(TrueChannel);
      Hists[TrueChannel].Fill(xy..., w);
    }
  }

  template <typename THOut>
  TrueChannelHist<THOut> Transform(std::function<THOut(TH const &)> f) {
    TrueChannelHist<THOut> out;
    out.fName = fName;
    for (auto const &h : Hists) {
      out.Hists[h.first] = THOut(f(h.second));
    }
    return out;
  }
};

template <typename TH> struct SelectionHists {
  std::vector<std::string> fSelections;
  std::vector<int> fSelectionsIds;
  std::unordered_map<int, TrueChannelHist<TH>> Hists;
  std::string fName;
  using THBase = typename THTraits<TH>::Base;

  SelectionHists() {}

  template <typename... Args>
  SelectionHists(std::string const &name, std::string const &title,
                 std::vector<std::string> const &Selections, Args... binning)
      : fName(name), fSelections(Selections) {
    for (int i = 0; i < fSelections.size(); ++i) {

      int sel = -1;
      for (int s = 0; s < t2knova::SelectionList.size(); ++s) {
        if (fSelections[i] == t2knova::SelectionList[s]) {
          sel = s;
          break;
        }
      }
      if (sel == -1) {
        std::cout << "{ERROR} Failed to parse selection: " << fSelections[i]
                  << std::endl;
        abort();
      }

      fSelectionsIds.push_back(sel);
      Hists.emplace(sel,
                    TrueChannelHist<TH>((fName + "_" + fSelections[i]).c_str(),
                                        title.c_str(), binning...));
    }
  }

  void SetName(std::string const &name) {
    fName = name;

    for (int i = 0; i < fSelections.size(); ++i) {
      Hists[fSelectionsIds[i]].SetName((fName + "_" + fSelections[i]).c_str());
    }
  }

  void SetTitle(std::string const &title) {
    Apply([=](TH &h) { h.SetTitle(title); });
  }

  void SetXAxisTitle(std::string const &title) {
    Apply([=](TH &h) { h.GetXaxis()->SetTitle(title.c_str()); });
  }
  void SetYAxisTitle(std::string const &title) {
    Apply([=](TH &h) { h.GetYaxis()->SetTitle(title.c_str()); });
  }
  void SetZAxisTitle(std::string const &title) {
    Apply([=](TH &h) { h.GetZaxis()->SetTitle(title.c_str()); });
  }

  void Apply(std::function<void(TH &)> f) {
    for (auto &h : Hists) {
      h.second.Apply(f);
    }
  }

  void Write(TDirectory *f, bool width_scale = false) {
    for (auto &h : Hists) {
      h.second.Write(f, width_scale);
    }
  }

  template <typename... XY>
  void Fill(double w, std::vector<int> const &Selections, int TrueChannel,
            XY... xy) {
    for (int sel : Selections) {
      if (Hists.count(sel)) {
        Hists[sel].Fill(w, TrueChannel, xy...);
      }
    }
  }

  template <typename THOut>
  SelectionHists<THOut> Transform(std::function<THOut(TH const &)> f) {
    SelectionHists<THOut> out;
    out.fSelections = fSelections;
    out.fSelectionsIds = fSelectionsIds;
    out.fName = fName;
    for (auto &h : Hists) {
      out.Hists.emplace(h.first, h.second.Transform(f));
    }
    return out;
  }
};

#ifndef NO_TOML

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
    toml::value const &plots_config) {

  if (!plots_config.contains(
          plotname)) { // If we cannot find the config just silently fail and
                       // create a histo blob that will never get
                       // filled/printed
    return new SelectionHists<TN>(plotname, "", {}, 1, 0, 1);
  }

  toml::value const &config = toml::find(plots_config, plotname);

  std::vector<std::string> axis_titles =
      toml_h::find<std::vector<std::string>>(config, "axis_titles");

  while (axis_titles.size() < 2) {
    axis_titles.push_back("");
  }

  std::string titlestr = toml_h::find<std::string>(config, "title") + ";" +
                         axis_titles[0] + ";" + axis_titles[1];

  std::vector<double> x_binning =
      BinningFromTOML(toml_h::find<toml::value>(config, "x_binning"));

  return new SelectionHists<TN>(
      plotname, titlestr,
      toml_h::find<std::vector<std::string>>(config, "selections"),
      x_binning.size() - 1, x_binning.data());
}

template <typename TN>
SelectionHists<TN> *SelectionHistsFromTOML(
    typename std::enable_if<THTraits<TN>::NDims == 2, std::string>::type const
        &plotname,
    toml::value const &plots_config) {

  if (!plots_config.contains(
          plotname)) { // If we cannot find the config just silently fail and
                       // create a histo blob that will never get
                       // filled/printed
    return new SelectionHists<TN>(plotname, "", {}, 1, 0, 1, 1, 0, 1);
  }

  toml::value const &config = toml::find(plots_config, plotname);

  std::vector<std::string> axis_titles =
      toml_h::find<std::vector<std::string>>(config, "axis_titles");

  while (axis_titles.size() < 3) {
    axis_titles.push_back("");
  }

  std::string titlestr = toml_h::find<std::string>(config, "title") + ";" +
                         axis_titles[0] + ";" + axis_titles[1] + ";" +
                         axis_titles[2];

  std::vector<double> x_binning =
      BinningFromTOML(toml_h::find<toml::value>(config, "x_binning"));
  std::vector<double> y_binning =
      BinningFromTOML(toml_h::find<toml::value>(config, "y_binning"));

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
    toml::value const &plots_config) {

  if (!plots_config.contains(
          plotname)) { // If we cannot find the config just silently fail and
                       // create a histo blob that will never get
                       // filled/printed
    return new SelectionHists<TN>(plotname, "", {}, 1, 0, 1, 1, 0, 1, 1, 0, 1);
  }

  toml::value const &config = toml::find(plots_config, plotname);

  std::vector<std::string> axis_titles =
      toml_h::find<std::vector<std::string>>(config, "axis_titles");

  while (axis_titles.size() < 4) {
    axis_titles.push_back("");
  }

  std::string titlestr = toml_h::find<std::string>(config, "title") + ";" +
                         axis_titles[0] + ";" + axis_titles[1] + ";" +
                         axis_titles[2] + ";" + axis_titles[3];

  std::vector<double> x_binning =
      BinningFromTOML(toml_h::find<toml::value>(config, "x_binning"));
  std::vector<double> y_binning =
      BinningFromTOML(toml_h::find<toml::value>(config, "y_binning"));
  std::vector<double> z_binning =
      BinningFromTOML(toml_h::find<toml::value>(config, "z_binning"));

  return new SelectionHists<TN>(
      plotname, titlestr,
      toml_h::find<std::vector<std::string>>(config, "selections"),
      x_binning.size() - 1, x_binning.data(), y_binning.size() - 1,
      y_binning.data(), z_binning.size() - 1, z_binning.data());
}
#endif

} // namespace t2knova