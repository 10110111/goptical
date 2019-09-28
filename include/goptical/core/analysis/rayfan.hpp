/*

      This file is part of the <goptical/core Core library.
  
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


#ifndef GOPTICAL_ANALYSIS_RAYFAN_HH_
#define GOPTICAL_ANALYSIS_RAYFAN_HH_

#include "goptical/core/common.hpp"

#include "goptical/core/io/renderer_axes.hpp"
#include "goptical/core/math/vector.hpp"
#include "goptical/core/trace/tracer.hpp"
#include "goptical/core/trace/distribution.hpp"
#include "goptical/core/data/plot.hpp"

namespace _goptical
{

  namespace analysis
  {

    /**
       @short RayFan diagram analysis
       @header <goptical/core/analysis/RayFan
       @module {Core}
       @main

       This class is designed to compute various ray fan plots.

       @xsee {tuto_fan1, tuto_fan2}
    */
    class RayFan
    {
    public:
      /** Specify ray aberration values to plot. Angle and Distance
          aberrations values are considered in selected aberration
          plane. Entrance Height and Angle are considered in selected
          distribution plane. */
      enum rayfan_plot_type_e
        {
          /** Normalized ray height (radial distance) on entrance pupil */
          EntranceHeight,
          /** Angle of ray on entrance pupil */
          EntranceAngle,
          /** Distance on the surface from the intercept of the chief ray */
          TransverseDistance,
          /** Distance along the chief ray from the surface to the measured ray. */
          LongitudinalDistance,
          /** Angle of ray striking the target surface */
          ImageAngle,
          /** Angle of ray leaving (generated by) the target surface */
          ExitAngle,
          /** Optical path difference in waves */
          OpticalPathDiff,
        };

      /** Specify aberration analysis plane on target surface */
      enum rayfan_plane_e
        {
          SagittalAberration = 0,
          TangentialAberration = 1
        };

      RayFan(const sys::System &system, enum rayfan_plane_e plane = TangentialAberration);

      /** Set entrance pupil ray distribution plane. */
      void set_plane(enum rayfan_plane_e plane);

      /** Get internal distribution object */
      inline trace::Distribution & get_distribution();

      /** Aberrations are considered in the given plane on the
          target surface. Default is to use the same plane as entrance
          pupil ray distribution plane. */
      inline void set_aberration_plane(enum rayfan_plane_e plane);

      /** Specify entrance pupil surface to use for analysis, query
          system for entrance pupil if none defined here. */
      inline void set_entrance_surface(const sys::Surface &s);

      /** Specify target surface (image or exit pupil) to use for
          analysis, query system for image surface if none defined
          here. */
      inline void set_target_surface(const sys::Surface &s);

      /** Set longitudinal reference ray (local to target
          surface). Longitudinal aberration computes distance between
          each rays and target surface plane along this reference
          vector. Default value is along the Z axis. */
      inline void set_longitudinal_reference(const math::VectorPair3 &ref);

      /** Get aberration plot, requested x value is plotted against
          requested y value. */
      ref<data::Plot> get_plot(enum rayfan_plot_type_e x,
                               enum rayfan_plot_type_e y);

      /** Invalidate current analysis data and raytrace again on next
          plot request */
      void invalidate();

    private:
      void process_trace();

      typedef double (RayFan::*get_value_t)(const trace::Ray &r, const trace::Ray &chief) const;

      const trace::Ray & find_chief_ray(const trace::rays_queue_t &intercepts, double wavelen);

      double get_entrance_height(const trace::Ray &r, const trace::Ray &chief) const;
      double get_entrance_angle(const trace::Ray &r, const trace::Ray &chief) const;
      double get_transverse_distance(const trace::Ray &r, const trace::Ray &chief) const;
      double get_longitudinal_distance(const trace::Ray &r, const trace::Ray &chief) const;
      double get_optical_path_len(const trace::Ray &r, const trace::Ray &chief) const;
      double get_image_angle(const trace::Ray &r, const trace::Ray &chief) const;
      double get_exit_angle(const trace::Ray &r, const trace::Ray &chief) const;

      trace::Tracer     _tracer;
      bool              _processed_trace;

      const sys::Surface *_entrance;
      const sys::Surface *_exit;

      trace::Distribution _dist;
      math::VectorPair3 _ref_ray;
      enum rayfan_plane_e _dist_plane;
      enum rayfan_plane_e _ab_plane;
    };

  }
}

#endif

