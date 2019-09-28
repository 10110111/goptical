/*

      This file is part of the <goptical/core Core library.
  
      The Goptical library is free software; you can redistribute it
      and/or modify it under the terms of the GNU General Public
      License as published by the Free Software Foundation; either
      version 3 of the License, or (at your option) any later version.
  
      The Goptical library is distributed in the hope that it will be
      useful, but WITHOUT ANY WARRANTY; without even the implied
      warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
      See the GNU General Public License for more details.
  
      You should have received a copy of the GNU General Public
      License along with the Goptical library; if not, write to the
      Free Software Foundation, Inc., 59 Temple Place, Suite 330,
      Boston, MA 02111-1307 USA
  
      Copyright (C) 2010-2011 Free Software Foundation, Inc
      Author: Alexandre Becoulet

*/


#ifndef GOPTICAL_ANALYSIS_RAYFAN_HXX_
#define GOPTICAL_ANALYSIS_RAYFAN_HXX_

#include "goptical/core/io/renderer_axes.hxx"
#include "goptical/core/math/vector.hxx"
#include "goptical/core/trace/tracer.hxx"
#include "goptical/core/trace/distribution.hxx"
#include "goptical/core/data/plot.hxx"

namespace _goptical
{

  namespace analysis
  {
    void RayFan::set_entrance_surface(const sys::Surface &s)
    {
      _entrance = &s;
    }

    void RayFan::set_target_surface(const sys::Surface &s)
    {
      _exit = &s;
    }

    void RayFan::set_aberration_plane(enum rayfan_plane_e plane)
    {
      assert(plane == 0 || plane == 1);
      _ab_plane = plane;
    }

    trace::Distribution & RayFan::get_distribution()
    {
      return _dist;
    }

  }

}

#endif

