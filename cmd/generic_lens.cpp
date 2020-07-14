/*

      This file is part of the <goptical/core library.

      The <goptical/core library is free software; you can redistribute it
      and/or modify it under the terms of the GNU General Public
      License as published by the Free Software Foundation; either
      version 3 of the License, or (at your option) any later version.

      The <goptical/core library is distributed in the hope that it will be
      useful, but WITHOUT ANY WARRANTY; without even the implied
      warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
      See the GNU General Public License for more details.

      You should have received a copy of the GNU General Public
      License along with the <goptical/core library; if not, write to the
      Free Software Foundation, Inc., 59 Temple Place, Suite 330,
      Boston, MA 02111-1307 USA

      Copyright (C) 2010-2011 Free Software Foundation, Inc
      Author: Alexandre Becoulet

*/

/* -*- indent-tabs-mode: nil -*- */

#include <fstream>
#include <iostream>

#include <goptical/core/math/vector.hpp>

#include <goptical/core/material/abbe.hpp>

#include <goptical/core/sys/image.hpp>
#include <goptical/core/sys/lens.hpp>
#include <goptical/core/sys/source.hpp>
#include <goptical/core/sys/source_point.hpp>
#include <goptical/core/sys/source_rays.hpp>
#include <goptical/core/sys/system.hpp>

#include <goptical/core/trace/tracer.hpp>

#include <goptical/core/analysis/focus.hpp>
#include <goptical/core/analysis/rayfan.hpp>
#include <goptical/core/analysis/spot.hpp>
#include <goptical/core/data/plot.hpp>

#include <goptical/core/io/renderer_svg.hpp>
#include <goptical/core/io/rgb.hpp>
#include <goptical/core/light/spectral_line.hpp>

#include <goptical/core/math/transform.hpp>

#include <goptical/core/io/import_bclaff.hpp>
#include <goptical/core/trace/distribution.hpp>
#include <goptical/core/trace/sequence.hpp>

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <iterator>
#include <memory>
#include <string>
#include <vector>

using namespace goptical;

struct BaseFileNames
{
  std::string layout_file;
  std::string spot_file;
  std::string spot_intensity_file;
  std::string opd_fan_file;
  std::string transverse_fan_file;
  std::string longitudinal_fan_file;
};

struct Args
{
  std::string input_file;
  bool refocus;
  unsigned scenario;
};

void analysis_fan (std::shared_ptr<sys::System> &sys,
                   const std::shared_ptr<sys::SourcePoint> &source_point,
                   const struct BaseFileNames &base_file_names);
void analysis_spot (std::shared_ptr<sys::System> &sys,
                    std::shared_ptr<sys::SourcePoint> &source_point,
                    const struct BaseFileNames &base_file_names,
                    bool skew = false);
void layout (const std::shared_ptr<sys::System> &sys,
             const std::shared_ptr<sys::SourcePoint> &source_point,
             const struct BaseFileNames &base_file_names, bool skew = false);

static std::string
get_base_name (const std::string &name)
{
  auto pos = name.find_last_of ('.');
  if (pos == std::string::npos)
    return name;
  return name.substr (0, pos);
}

std::shared_ptr<sys::SourcePoint>
setup_point_source (std::shared_ptr<sys::System> &sys, double angleOfView,
                    bool parallel)
{
  sys::SourceInfinityMode mode = sys::SourceAtInfinity;
  auto unit_vector = math::vector3_001;
  if (!parallel)
    {
      // Construct unit vector at an angle
      //      double z1 = cos (angleOfView);
      //      double y1 = sin (angleOfView);
      //      unit_vector = math::Vector3 (0, y1, z1);

      math::Matrix<3> r;
      math::get_rotation_matrix (r, 0, angleOfView);
      unit_vector = r * math::vector3_001;
    }
  auto source_point = std::make_shared<sys::SourcePoint> (mode, unit_vector);
  source_point->clear_spectrum ();
  source_point->add_spectral_line (light::SpectralLine::d);
  source_point->add_spectral_line (light::SpectralLine::C);
  source_point->add_spectral_line (light::SpectralLine::F);
  sys->add (source_point);

  return source_point;
}

// TODO
// setup_source_rays() {
//  double z = 1e5;
//  double h = tan(angleOfView) * z;
//
//  sys::SourceRays  source_rays(math::Vector3(0, h, -z));
//
// sys::SourcePoint source_point(sys::SourceAtFiniteDistance,
//                              math::Vector3(0, h, -z));
// sys::SourcePoint source_point(sys::SourceAtFiniteDistance, math::Vector3(0,
// h, -z));
//  sys::SourcePoint source_point(sys::SourceAtInfinity, math::Vector3(0, 0,
//  1));

// add sources to system
//  sys.add(source_rays);
// configure sources
//    source_rays.add_chief_rays(sys);
//    source_rays.add_marginal_rays(sys, 14);

//  source_rays.clear_spectrum();
//  source_rays.add_spectral_line(light::SpectralLine::d);
//  source_rays.add_spectral_line(light::SpectralLine::C);
//  source_rays.add_spectral_line(light::SpectralLine::F);
//}

static void
get_base_file_names (const std::string &input_file,
                     BaseFileNames *base_file_names)
{
  std::string output_base = get_base_name (input_file);
  base_file_names->layout_file = output_base + "_layout";
  base_file_names->spot_file = output_base + "_spot";
  base_file_names->spot_intensity_file = output_base + "_spot_intensity";
  base_file_names->opd_fan_file = output_base + "_opd_fan";
  base_file_names->transverse_fan_file = output_base + "_transverse_fan";
  base_file_names->longitudinal_fan_file = output_base + "_longitudinal_fan";
}

static bool
get_arguments (int argc, const char *argv[], Args *args)
{
  args->refocus = false;
  args->scenario = 0;
  if (argc < 2)
    {
      fprintf (stderr, "Please supply a data file\n");
      return false;
    }
  for (int i = 2; i < argc; i++)
    {
      if (strcmp (argv[i], "--refocus") == 0)
        args->refocus = true;
      else if (strcmp (argv[i], "--scenario") == 0 && i + 1 < argc)
        {
          i++;
          args->scenario = (unsigned)atoi (argv[i]);
        }
    }
  args->input_file = std::string (argv[1]);
  return true;
}

static void
do_system_parallel_rays (io::BClaffLensImporter *importer,
                         const BaseFileNames &base_file_names,
                         const Args &args)
{
  auto sys = importer->buildSystem (args.scenario);
  double angleOfView = importer->getAngleOfViewInRadians (args.scenario);

  /* anchor sources */
  auto source_point = setup_point_source (sys, angleOfView, true);

  /* anchor seq */
  auto seq = std::make_shared<trace::Sequence> (*sys);
  sys->get_tracer_params ().set_sequential_mode (seq);

  std::cout << "system:" << std::endl << *sys;
  std::cout << "sequence:" << std::endl << *seq;
  /* anchor end */

  if (args.refocus)
    {
      /* anchor focus */
      analysis::Focus focus (sys);
      auto best_focus = focus.get_best_focus ();
      std::cout << "Best focus found at " << best_focus << "\n";
      importer->get_image ()->set_plane (best_focus);
    }
  layout (sys, source_point, base_file_names, false);
  analysis_spot (sys, source_point, base_file_names, false);
  analysis_fan (sys, source_point, base_file_names);
}

static void
do_skew_rays (io::BClaffLensImporter *importer,
              const BaseFileNames &base_file_names, const Args &args)
{
  auto sys = importer->buildSystem (args.scenario);
  double angleOfView = importer->getAngleOfViewInRadians (args.scenario);

  /* anchor sources */
  auto source_point = setup_point_source (sys, angleOfView, false);

  /* anchor seq */
  auto seq = std::make_shared<trace::Sequence> (*sys);
  sys->get_tracer_params ().set_sequential_mode (seq);

  if (args.refocus)
    {
      /* anchor focus */
      analysis::Focus focus (sys);
      auto best_focus = focus.get_best_focus ();
      std::cout << "Best focus found at " << best_focus << "\n";
      importer->get_image ()->set_plane (best_focus);
    }
  layout (sys, source_point, base_file_names, true);
  analysis_spot (sys, source_point, base_file_names, true);
}

int
main (int argc, const char *argv[])
{
  //**********************************************************************
  // Optical system definition
  Args arguments;
  if (!get_arguments (argc, argv, &arguments))
    {
      exit (1);
    }

  io::BClaffLensImporter importer;

  BaseFileNames base_file_names;
  get_base_file_names (arguments.input_file, &base_file_names);

  if (!importer.parseFile (arguments.input_file))
    {
      std::cerr << "Failed to parse file " << arguments.input_file << "\n";
      exit (1);
    }

  do_system_parallel_rays (&importer, base_file_names, arguments);
  do_skew_rays (&importer, base_file_names, arguments);

  return 0;
}

//**********************************************************************
// Drawing rays and layout

void
layout (const std::shared_ptr<sys::System> &sys,
        const std::shared_ptr<sys::SourcePoint> &source_point,
        const struct BaseFileNames &base_file_names, bool skew)
{
  {
    /* anchor layout */
    std::string filename
        = base_file_names.layout_file + (skew ? "_skew" : "") + ".svg";
    io::RendererSvg renderer (filename.c_str (), 800, 400);

    // draw 2d system layout
    sys->draw_2d_fit (renderer);
    sys->draw_2d (renderer);

    trace::Tracer tracer (sys.get ());

#if 0
    // trace and draw rays from rays source
    sys.enable_single<sys::Source>(source_rays);
    tracer.get_trace_result().set_generated_save_state(source_rays);
    tracer.trace();
    tracer.get_trace_result().draw_2d(renderer);
    /* anchor end */
#else
    // trace and draw rays from source
    tracer.get_params ().set_default_distribution (
        trace::Distribution (trace::MeridionalDist, 50));
    tracer.get_trace_result ().set_generated_save_state (*source_point);
    tracer.trace ();
    tracer.get_trace_result ().draw_2d (renderer);
#endif
  }
}
void
analysis_spot (std::shared_ptr<sys::System> &sys,
               std::shared_ptr<sys::SourcePoint> &source_point,
               const struct BaseFileNames &base_file_names, bool skew)
{
  /* anchor spot */
  sys->enable_single<sys::Source> (*source_point);

  sys->get_tracer_params ().set_default_distribution (
      trace::Distribution (trace::HexaPolarDist, 20));

  analysis::Spot spot (sys);

  /* anchor end */
  {
    /* anchor spot */
    std::string filename
        = base_file_names.spot_file + (skew ? "_skew" : "") + ".svg";
    io::RendererSvg renderer (filename.c_str (), 300, 300, io::rgb_black);

    spot.draw_diagram (renderer);
    /* anchor end */
  }

  {
    /* anchor spot_plot */
    std::string filename
        = base_file_names.spot_intensity_file + (skew ? "_skew" : "") + ".svg";
    io::RendererSvg renderer (filename.c_str (), 640, 480);

    std::shared_ptr<data::Plot> plot = spot.get_encircled_intensity_plot (50);

    plot->draw (renderer);
    /* anchor end */
  }
}
void
analysis_fan (std::shared_ptr<sys::System> &sys,
              const std::shared_ptr<sys::SourcePoint> &source_point,
              const struct BaseFileNames &base_file_names)
{
  /* anchor opd_fan */
  sys->enable_single<sys::Source> (*source_point);

  analysis::RayFan fan (sys);

  /* anchor end */
  {
    /* anchor opd_fan */
    std::string filename = base_file_names.opd_fan_file + ".svg";
    io::RendererSvg renderer (filename.c_str (), 640, 480);

    std::shared_ptr<data::Plot> fan_plot = fan.get_plot (
        analysis::RayFan::EntranceHeight, analysis::RayFan::OpticalPathDiff);

    fan_plot->draw (renderer);

    /* anchor end */
  }

  {
    /* anchor transverse_fan */
    std::string filename = base_file_names.transverse_fan_file + ".svg";
    io::RendererSvg renderer (filename.c_str (), 640, 480);

    std::shared_ptr<data::Plot> fan_plot
        = fan.get_plot (analysis::RayFan::EntranceHeight,
                        analysis::RayFan::TransverseDistance);

    fan_plot->draw (renderer);

    /* anchor end */
  }

  {
    /* anchor longitudinal_fan */
    std::string filename = base_file_names.longitudinal_fan_file + ".svg";
    io::RendererSvg renderer (filename.c_str (), 640, 480);

    std::shared_ptr<data::Plot> fan_plot
        = fan.get_plot (analysis::RayFan::EntranceHeight,
                        analysis::RayFan::LongitudinalDistance);

    fan_plot->draw (renderer);

    /* anchor end */
  }
}
