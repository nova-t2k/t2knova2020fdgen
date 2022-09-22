#pragma once

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

#include <iostream>
#include <memory>
#include <string>

template <typename TOut, typename TIn>
std::unique_ptr<TOut> dynamic_cast_uptr(std::unique_ptr<TIn> &&in) {
  TOut *optr = dynamic_cast<TOut *>(in.get());
  if (optr) { // take the pointer if its a good one
    in.release();
  }
  return std::unique_ptr<TOut>(optr);
}

template <typename T, typename... args>
std::unique_ptr<T> make_unique(args... a) {
  return std::unique_ptr<T>(new T(a...));
}

template <typename TH> struct THTraits {};

template <> struct THTraits<TH1> {
  using Base = TH1;
  using FloatType = TH1F;
  static size_t const NDims = 1;
};
template <> struct THTraits<TH1D> {
  using Base = TH1;
  using FloatType = TH1F;
  static size_t const NDims = 1;
};
template <> struct THTraits<TH1F> {
  using Base = TH1;
  using FloatType = TH1F;
  static size_t const NDims = 1;
};

template <> struct THTraits<TH2> {
  using Base = TH2;
  using FloatType = TH2F;
  static size_t const NDims = 2;
};
template <> struct THTraits<TH2D> {
  using Base = TH2;
  using FloatType = TH2F;
  static size_t const NDims = 2;
};
template <> struct THTraits<TH2F> {
  using Base = TH2;
  using FloatType = TH2F;
  static size_t const NDims = 2;
};

template <> struct THTraits<TH3> {
  using Base = TH3;
  using FloatType = TH3F;
  static size_t const NDims = 3;
};
template <> struct THTraits<TH3D> {
  using Base = TH3;
  using FloatType = TH3F;
  static size_t const NDims = 3;
};
template <> struct THTraits<TH3F> {
  using Base = TH3;
  using FloatType = TH3F;
  static size_t const NDims = 3;
};

inline std::vector<double> GetTAxisBinning(TAxis *const ax) {
  std::vector<double> binning;

  for (int i = 0; i < ax->GetNbins(); ++i) {
    binning.push_back(ax->GetBinLowEdge(i + 1));
  }
  binning.push_back(ax->GetBinUpEdge(ax->GetNbins()));
  return binning;
}

template <typename TN>
std::unique_ptr<typename THTraits<TN>::FloatType>
EnsureTHF(typename std::enable_if<THTraits<TN>::NDims == 1,
                                  std::unique_ptr<TN>>::type const &hin) {
  std::vector<double> xbin = GetTAxisBinning(hin->GetXaxis());

  std::unique_ptr<typename THTraits<TN>::FloatType> ret(
      new typename THTraits<TN>::FloatType(hin->GetName(), hin->GetTitle(),
                                           xbin.size() - 1, xbin.data()));

  for (int i = 0; i < ret->GetXaxis()->GetNbins(); ++i) {
    ret->SetBinContent(i + 1, hin->GetBinContent(i + 1));
    ret->SetBinError(i + 1, hin->GetBinError(i + 1));
  }

  return ret;
}

template <typename TN>
std::unique_ptr<typename THTraits<TN>::FloatType>
EnsureTHF(typename std::enable_if<THTraits<TN>::NDims == 2,
                                  std::unique_ptr<TN>>::type const &hin) {
  std::vector<double> xbin = GetTAxisBinning(hin->GetXaxis());
  std::vector<double> ybin = GetTAxisBinning(hin->GetYaxis());

  std::unique_ptr<typename THTraits<TN>::FloatType> ret(
      new typename THTraits<TN>::FloatType(hin->GetName(), hin->GetTitle(),
                                           xbin.size() - 1, xbin.data(),
                                           ybin.size() - 1, ybin.data()));

  for (int i = 0; i < ret->GetXaxis()->GetNbins(); ++i) {
    for (int j = 0; j < ret->GetXaxis()->GetNbins(); ++j) {
      ret->SetBinContent(i + 1, j + 1, hin->GetBinContent(i + 1, j + 1));
      ret->SetBinError(i + 1, j + 1, hin->GetBinError(i + 1, j + 1));
    }
  }

  return ret;
}

inline std::unique_ptr<TH1> GetTH1(std::unique_ptr<TFile> &f,
                                   std::string const &name, bool quiet = true) {
  TDirectory *odir = gDirectory;

  TH1 *h;
  f->GetObject(name.c_str(), h);
  if (!h) {
    if (!quiet) {
      std::cout << "[ERROR]: Failed to find histogram named: " << name
                << std::endl;
    }
  } else {
    h->SetDirectory(nullptr);
  }

  if (odir) {
    gDirectory->cd();
  }

  return std::unique_ptr<TH1>(h);
}

inline std::unique_ptr<TH1> GetTH1(std::string const &fname,
                                   std::string const &name, bool quiet = true) {
  TDirectory *odir = gDirectory;

  auto f = ::make_unique<TFile>(fname.c_str());
  std::unique_ptr<TH1> h = GetTH1(f, name, quiet);

  if (odir) {
    gDirectory->cd();
  }

  return h;
}

template <typename TH>
std::unique_ptr<TH>
GetTH(typename std::enable_if<std::is_base_of<TH1, TH>::value,
                              std::unique_ptr<TFile>>::type &f,
      std::string const &name, bool quiet = true) {
  std::unique_ptr<TH> h = dynamic_cast_uptr<TH>(GetTH1(f, name, quiet));
  if (!h) {
    if (!quiet) {
      std::cout << "[ERROR]: Failed to find histogram named: " << name
                << std::endl;
    }
  }
  return h;
}

template <typename TH>
inline std::unique_ptr<TH>
GetTH(typename std::enable_if<std::is_base_of<TH1, TH>::value,
                              std::string const>::type &fname,
      std::string const &name, bool quiet = true) {
  std::unique_ptr<TH> h = dynamic_cast_uptr<TH>(GetTH1(fname, name, quiet));
  if (!h) {
    if (!quiet) {
      std::cout << "[ERROR]: Failed to find histogram named: " << name
                << std::endl;
    }
  }
  return h;
}

inline double EvalHist3D(std::unique_ptr<TH1> const &h1, double x, double y,
                         double z, bool interpolate = true) {
  TH3 *h = static_cast<TH3 *>(h1.get());

  int xbin = h->GetXaxis()->FindFixBin(x);
  int ybin = h->GetYaxis()->FindFixBin(y);
  int zbin = h->GetZaxis()->FindFixBin(z);

  if ((xbin != 0) && (xbin != (h->GetXaxis()->GetNbins() + 1)) && (ybin != 0) &&
      (ybin != (h->GetYaxis()->GetNbins() + 1)) && (zbin != 0) &&
      (zbin != (h->GetZaxis()->GetNbins() + 1))) {
    if (interpolate && (xbin > 1) && (xbin < h->GetXaxis()->GetNbins()) &&
        (ybin > 1) && (ybin < h->GetYaxis()->GetNbins()) && (zbin > 1) &&
        (zbin < h->GetZaxis()->GetNbins())) {
      return h->Interpolate(x, y, z);
    }
    return h->GetBinContent(xbin, ybin, zbin);
  }

  return 1;
}

inline double EvalHist1D(std::unique_ptr<TH1> const &h, double x,
                         bool interpolate = true) {
  int xbin = h->GetXaxis()->FindFixBin(x);

  if ((xbin != 0) && (xbin != (h->GetXaxis()->GetNbins() + 1))) {
    if (interpolate && (xbin > 1) && (xbin < h->GetXaxis()->GetNbins())) {
      return h->Interpolate(x);
    }
    return h->GetBinContent(xbin);
  }

  return 1;
}

TDirectory *MakeDirectoryStructure(TDirectory *din, std::string path = "") {
  if (path.size()) {
    std::string tab = "";
    while (path.find("/") != std::string::npos) {
      std::string ndir = path.substr(0, path.find("/"));
      // std::cout << tab << "Next Dir: " << ndir << std::endl;
      if (!din->GetDirectory(ndir.c_str())) {
        // std::cout << tab << "Need to make it" << std::endl;
        din = din->mkdir(ndir.c_str());
      } else {
        // std::cout << tab << "--Have it!" << std::endl;
        din = din->GetDirectory(ndir.c_str());
      }
      path = path.substr(path.find("/") + 1);
      // std::cout << tab << "Remaining path: " << path << std::endl;
      tab = tab + "  ";
    }
    if (path.size()) {
      // std::cout << tab << "Next Dir: " << path << std::endl;
      if (!din->GetDirectory(path.c_str())) {
        // std::cout << tab << "Need to make it" << std::endl;
        din = din->mkdir(path.c_str());
      } else {
        // std::cout << tab << "--Have it!" << std::endl;
        din = din->GetDirectory(path.c_str());
      }
    }
  }
  return din;
}

double IntegralTH2(std::unique_ptr<TH2> const &H,
                   double frac_error_threshold = 1,
                   std::string const &opt = "width") {
  if (!H) {
    return 0;
  }

  bool dowidth = (opt == "width");
  double sum = 0;

  for (int i = 0; i < H->GetXaxis()->GetNbins(); ++i) {
    double bw_x = H->GetXaxis()->GetBinUpEdge(i + 1) -
                  H->GetXaxis()->GetBinLowEdge(i + 1);
    for (int j = 0; j < H->GetYaxis()->GetNbins(); ++j) {
      double bw_y = H->GetYaxis()->GetBinUpEdge(i + 1) -
                    H->GetYaxis()->GetBinLowEdge(i + 1);

      double bc = H->GetBinContent(i + 1, j + 1);
      double bin_frac_error = H->GetBinError(i + 1, j + 1) / bc;

      if ((std::isnormal(bin_frac_error)) &&
          (bin_frac_error < frac_error_threshold)) {
        sum += dowidth ? (bw_x * bw_y * bc) : (bc);
      }
    }
  }
  return sum;
}

void ScrubLowStatsBins(std::unique_ptr<TH3> &num, std::unique_ptr<TH3> &denom,
                       std::unique_ptr<TH3> &ratio,
                       double frac_error_threshold) {
  for (int i = 0; i < num->GetXaxis()->GetNbins(); ++i) {
    for (int j = 0; j < num->GetYaxis()->GetNbins(); ++j) {
      for (int k = 0; k < num->GetZaxis()->GetNbins(); ++k) {
        double num_frac_error = num->GetBinError(i + 1, j + 1, k + 1) /
                                num->GetBinContent(i + 1, j + 1, k + 1);
        double denom_frac_error = denom->GetBinError(i + 1, j + 1, k + 1) /
                                  denom->GetBinContent(i + 1, j + 1, k + 1);

        if (!std::isnormal(num_frac_error) ||
            !std::isnormal(denom_frac_error)) {
          ratio->SetBinContent(i + 1, j + 1, k + 1, 0);
        } else if (num_frac_error > frac_error_threshold) {
          ratio->SetBinContent(i + 1, j + 1, k + 1, 0);
        } else if (denom_frac_error > frac_error_threshold) {
          ratio->SetBinContent(i + 1, j + 1, k + 1, 1);
        }
      }
    }
  }
}

void ScrubLowStatsBins(std::unique_ptr<TH2> &num, std::unique_ptr<TH2> &denom,
                       std::unique_ptr<TH2> &ratio,
                       double frac_error_threshold) {
  for (int i = 0; i < num->GetXaxis()->GetNbins(); ++i) {
    for (int j = 0; j < num->GetYaxis()->GetNbins(); ++j) {
      double num_frac_error =
          num->GetBinError(i + 1, j + 1) / num->GetBinContent(i + 1, j + 1);
      double denom_frac_error =
          denom->GetBinError(i + 1, j + 1) / denom->GetBinContent(i + 1, j + 1);

      if (!std::isnormal(num_frac_error) || !std::isnormal(denom_frac_error)) {
        ratio->SetBinContent(i + 1, j + 1, 0);
      } else if (num_frac_error > frac_error_threshold) {
        ratio->SetBinContent(i + 1, j + 1, 0);
      } else if (denom_frac_error > frac_error_threshold) {
        ratio->SetBinContent(i + 1, j + 1, 1);
      }
    }
  }
}

void ScrubLowStatsBins(std::unique_ptr<TH1> &num, std::unique_ptr<TH1> &denom,
                       std::unique_ptr<TH1> &ratio,
                       double frac_error_threshold) {

  if (num->GetDimension() == 2) { // this is awful, I know
    std::unique_ptr<TH2> num_2(dynamic_cast<TH2 *>(num.get()));
    std::unique_ptr<TH2> denom_2(dynamic_cast<TH2 *>(denom.get()));
    std::unique_ptr<TH2> ratio_2(dynamic_cast<TH2 *>(ratio.get()));
    ScrubLowStatsBins(num_2, denom_2, ratio_2, frac_error_threshold);
    num_2.release();
    denom_2.release();
    ratio_2.release();
    return;
  } else if (num->GetDimension() == 3) {
    std::unique_ptr<TH3> num_3(dynamic_cast<TH3 *>(num.get()));
    std::unique_ptr<TH3> denom_3(dynamic_cast<TH3 *>(denom.get()));
    std::unique_ptr<TH3> ratio_3(dynamic_cast<TH3 *>(ratio.get()));
    ScrubLowStatsBins(num_3, denom_3, ratio_3, frac_error_threshold);
    num_3.release();
    denom_3.release();
    ratio_3.release();
    return;
  }

  for (int i = 0; i < num->GetXaxis()->GetNbins(); ++i) {
    double num_frac_error = num->GetBinError(i + 1) / num->GetBinContent(i + 1);
    double denom_frac_error =
        denom->GetBinError(i + 1) / denom->GetBinContent(i + 1);

    if (!std::isnormal(num_frac_error) || !std::isnormal(denom_frac_error)) {
      ratio->SetBinContent(i + 1, 0);
    } else if (num_frac_error > frac_error_threshold) {
      ratio->SetBinContent(i + 1, 0);
    } else if (denom_frac_error > frac_error_threshold) {
      ratio->SetBinContent(i + 1, 1);
    }
  }
}
