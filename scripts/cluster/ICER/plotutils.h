#pragma once

#include "TAxis.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TGaxis.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TPad.h"
#include "TSpline.h"

#include <iostream>
#include <limits>
#include <memory>
#include <vector>

TGraph GetGraphFromTH1(TH1 *h, double min, double max, bool close_graph = false,
                       int nsteps = 1000) {
  TSpline3 spl(h, 0, min, max);

  TGraph g = TGraph(nsteps + (close_graph ? 1 : 0));
  if (close_graph) {
    g.SetPoint(0, min, 0);
  }
  for (int i = (close_graph ? 1 : 0); i < nsteps; ++i) {
    double e = min + i * ((max - min) / double(nsteps));
    g.SetPoint(i, e, spl.Eval(e));
  }
  if (close_graph) {
    g.SetPoint(nsteps, max, 0);
  }

  return g;
}

double GetMaximumTGraph(TGraph const &g) {
  double m = -std::numeric_limits<double>::max();
  for (int i = 0; i < g.GetN(); ++i) {
    double x, y;
    g.GetPoint(i, x, y);
    m = std::max(m, y);
  }
  return m;
}

double GetMinimumTGraph(TGraph const &g) {
  double m = std::numeric_limits<double>::max();
  for (int i = 0; i < g.GetN(); ++i) {
    double x, y;
    g.GetPoint(i, x, y);
    m = std::min(m, y);
  }
  return m;
}

double GetMaximumXTGraph(TGraph const &g) {
  double m = -std::numeric_limits<double>::max();
  for (int i = 0; i < g.GetN(); ++i) {
    double x, y;
    g.GetPoint(i, x, y);
    m = std::max(m, x);
  }
  return m;
}

double GetMinimumXTGraph(TGraph const &g) {
  double m = std::numeric_limits<double>::max();
  for (int i = 0; i < g.GetN(); ++i) {
    double x, y;
    g.GetPoint(i, x, y);
    m = std::min(m, x);
  }
  return m;
}

double GetMaximumTGraphs(std::vector<TGraph> gs) {
  double max = -std::numeric_limits<double>::max();

  for (auto const &g : gs) {
    max = std::max(max, GetMaximumTGraph(g));
  }

  return max;
}

TGraph ScaleTGraph(TGraph const &g, double s) {
  TGraph o(g.GetN());

  for (int i = 0; i < g.GetN(); ++i) {
    double x, y;
    g.GetPoint(i, x, y);
    o.SetPoint(i, x, y * s);
  }
  return o;
}

TGraph CloneTGraph(TGraph const &g) { return ScaleTGraph(g, 1); }

void StyleTGraphLine(TGraph &h, int color, int width = 1, int style = 1) {
  h.SetLineStyle(style);
  h.SetLineColor(color);
  h.SetLineWidth(width);
}

void StyleTGraphFill(TGraph &h, int color, int style = 1001) {
  h.SetFillStyle(style);
  h.SetFillColor(color);
}

void DrawTGraphs(std::vector<TGraph> gs, bool first = true,
                 std::string opts = "") {
  for (auto const &g : gs) {
    if (first) {
      g.GetYaxis()->SetRangeUser(0, GetMaximumTGraphs(gs) * 1.05);
    }
    std::string dopt = std::string(first ? "A" : "") + opts;
    g.DrawClone(dopt.c_str());
    first = false;
  }
}

TGraph ConcatenateTGraphs(std::vector<TGraph> gs) {
  TGraph gout;
  int it = 0;
  for (auto const &g : gs) {
    for (int i = 0; i < g.GetN(); ++i) {
      double x, y;
      g.GetPoint(i, x, y);
      gout.SetPoint(it++, x, y);
    }
  }
  return gout;
}

TGraph SumTGraphs(std::vector<TGraph> gs) {
  TGraph gout;
  for (int i = 0; i < gs.front().GetN(); ++i) {
    double yall = 0;
    double x;
    for (auto const &g : gs) {
      double y;
      g.GetPoint(i, x, y);
      yall += y;
    }
    gout.SetPoint(i, x, yall);
  }

  gout.SetLineColor(gs.front().GetLineColor());
  gout.SetLineStyle(gs.front().GetLineStyle());
  gout.SetLineWidth(gs.front().GetLineWidth());
  return gout;
}

TGraph PointMultiplyTGraphs(std::vector<TGraph> gs) {
  TGraph gout;
  for (int i = 0; i < gs.front().GetN(); ++i) {
    double yall = 1;
    double x;
    for (auto const &g : gs) {
      double y;
      g.GetPoint(i, x, y);
      yall *= y;
    }
    gout.SetPoint(i, x, yall);
  }

  gout.SetLineColor(gs.front().GetLineColor());
  gout.SetLineStyle(gs.front().GetLineStyle());
  gout.SetLineWidth(gs.front().GetLineWidth());
  return gout;
}

TGraph MultiplyTGraphs(std::vector<TGraph> gs) {
  TGraph gout;
  for (int i = 0; i < gs.front().GetN(); ++i) {
    double yall;
    double x;
    int ctr = 0;
    for (auto const &g : gs) {

      if (ctr++ == 0) {
        g.GetPoint(i, x, yall);
      } else {
        yall *= g.Eval(x);
      }
    }
    gout.SetPoint(i, x, yall);
  }

  gout.SetLineColor(gs.front().GetLineColor());
  gout.SetLineStyle(gs.front().GetLineStyle());
  gout.SetLineWidth(gs.front().GetLineWidth());
  return gout;
}

TH1 *GetTH1(TFile *f, std::string const &name) {
  TDirectory *odir = gDirectory;

  TH1 *h;
  f->GetObject(name.c_str(), h);
  if (!h) {
    std::cout << "[ERROR]: Failed to find histogram named: " << name
              << std::endl;
  } else {
    h->SetDirectory(nullptr);
  }

  if (odir) {
    gDirectory->cd();
  }

  return h;
}

template <typename TH> double IntegralTH(TH *h) {
  if (!h) {
    return 0;
  }
  return h->Integral();
}

TH1 *GetTH1(std::string const &fname, std::string const &name) {
  TDirectory *odir = gDirectory;

  std::unique_ptr<TFile> f(TFile::Open(fname.c_str()));

  TH1 *h = GetTH1(f.get(), name);

  if (odir) {
    gDirectory->cd();
  }

  return h;
}

TH2 *GetTH2(TFile *f, std::string const &name) {
  TH2 *h = dynamic_cast<TH2 *>(GetTH1(f, name));
  if (!h) {
    std::cout << "[ERROR]: Failed to find histogram named: " << name
              << std::endl;
  }
  return h;
}

TH2 *GetTH2(std::string const &fname, std::string const &name) {
  TH2 *h = dynamic_cast<TH2 *>(GetTH1(fname, name));
  if (!h) {
    std::cout << "[ERROR]: Failed to find histogram named: " << name
              << std::endl;
  }
  return h;
}

void StyleAxis(TAxis *a, double scale = 1, double titleoffsetscale = 1,
               double labeloffsetscale = 1) {
  a->SetNdivisions(505);
  a->SetLabelSize(0.05 * scale);
  a->SetLabelOffset(0.01 * labeloffsetscale);
  a->SetTitleSize(0.05 * scale);
  a->SetTitleOffset(1.2 * titleoffsetscale);
}

void StyleAxes(TH1 *a, double scale = 1, double titleoffsetscale = 1,
               double labeloffsetscale = 1) {
  StyleAxis(a->GetYaxis(), scale, titleoffsetscale, labeloffsetscale);
  StyleAxis(a->GetXaxis(), scale, titleoffsetscale, labeloffsetscale);
}

void HideAxis(TAxis *a) {
  a->SetLabelSize(0);
  a->SetLabelOffset(10);
  a->SetTitleSize(0);
  a->SetTitleOffset(10);
}

void StyleAxis(TGaxis *g, TAxis *a) {
  g->SetLabelFont(a->GetLabelFont());
  g->SetLabelSize(a->GetLabelSize());
  g->SetLabelOffset(2 * a->GetLabelOffset());
  g->SetTitleFont(a->GetTitleFont());
  g->SetTitleSize(a->GetTitleSize());
  g->SetTitleOffset(a->GetTitleOffset());
}

void StyleTH1Line(TH1 *h, int color, int width = 1, int style = 1,
                  double alpha = 1) {
  if (!h) {
    return;
  }
  h->SetLineStyle(style);
  h->SetLineColorAlpha(color, alpha);
  h->SetLineWidth(width);
}

void StyleTH1Fill(TH1 *h, int color, int style = 1001, double alpha = 1) {
  if (!h) {
    return;
  }
  h->SetFillStyle(style);
  h->SetFillColorAlpha(color, alpha);
}

double GetMaximumTH1s(std::vector<TH1 *> hs) {
  double max = -std::numeric_limits<double>::max();

  for (auto const &h : hs) {
    if (!h) {
      continue;
    }
    max = std::max(max, h->GetMaximum());
  }

  return max;
}

void DrawTH1s(std::vector<TH1 *> hs, std::string opts = "", bool first = true) {
  for (auto const &h : hs) {
    if (!h) {
      continue;
    }
    if (first) {
      h->GetYaxis()->SetRangeUser(0, GetMaximumTH1s(hs) * 1.05);
    }
    std::string dopt = std::string(!first ? "SAME" : "") + opts;
    h->DrawClone(dopt.c_str());
    first = false;
  }
}

void ScaleTH1s(std::vector<TH1 *> hs, double s, std::string opts = "") {
  for (auto &h : hs) {
    if (!h) {
      continue;
    }
    h->Scale(s, opts.c_str());
  }
}

TCanvas *MakeCanvas(double mt = 0.03, double ml = 0.15, double mr = 0.03,
                    double mb = 0.15, std::string name = "c1") {
  TCanvas *c1 = new TCanvas(name.c_str(), "", 600, 600);
  c1->SetTopMargin(mt);
  c1->SetLeftMargin(ml);
  c1->SetRightMargin(mr);
  c1->SetBottomMargin(mb);
  return c1;
}

TPad *MakeRatioTopPad(double mt = 0.05, double ml = 0.15, double mr = 0.03,
                      double mb = 0.03, double midpoint = 0.4,
                      std::string name = "ptop") {
  TPad *c1 = new TPad(name.c_str(), "", 0, midpoint, 1, 1);
  c1->SetTopMargin(mt);
  c1->SetLeftMargin(ml);
  c1->SetRightMargin(mr);
  c1->SetBottomMargin(mb);
  return c1;
}
TPad *MakeRatioTopPadTopLegend(double mt = 0.3, double ml = 0.15,
                               double mr = 0.03, double mb = 0.03,
                               double midpoint = 0.3,
                               std::string name = "ptop") {
  return MakeRatioTopPad(mt, ml, mr, mb, midpoint, name);
}
TPad *MakeRatioBottomPad(double mt = 0.03, double ml = 0.15, double mr = 0.03,
                         double mb = 0.2, double midpoint = 0.4,
                         std::string name = "pbottom") {
  TPad *c1 = new TPad(name.c_str(), "", 0, 0, 1, midpoint);
  c1->SetTopMargin(mt);
  c1->SetLeftMargin(ml);
  c1->SetRightMargin(mr);
  c1->SetBottomMargin(mb);
  return c1;
}
TPad *MakeRatioBottomPadTopLegend(double mt = 0.03, double ml = 0.15,
                                  double mr = 0.03, double mb = 0.25,
                                  double midpoint = 0.3,
                                  std::string name = "pbottom") {
  return MakeRatioBottomPad(mt, ml, mr, mb, midpoint, name);
}

TCanvas *MakeCanvasTopLegend(std::string name = "c1") {
  return MakeCanvas(0.2, 0.15, 0.03, 0.15, name);
}

TLegend *MakeTopLegend() {
  auto leg = new TLegend(0, 0.81, 0.97, 1);
  leg->SetNColumns(2);
  leg->SetTextSize(0.07);
  leg->SetBorderSize(0);
  leg->SetFillStyle(4001);
  return leg;
}

TGaxis *MakeOtherAxis(double rightmax, double &scale) {
  scale = gPad->GetUymax() / rightmax;

  return new TGaxis(gPad->GetUxmax(), gPad->GetUymin(), gPad->GetUxmax(),
                    gPad->GetUymax(), 0, rightmax, 510, "+L");
}