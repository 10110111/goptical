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

#include <fstream>
#include <algorithm>
#include <goptical/core/io/renderer.hpp>

#include <goptical/core/trace/params.hpp>
#include <goptical/core/trace/ray.hpp>
#include <goptical/core/trace/result.hpp>

#include <goptical/core/sys/image.hpp>

#include <goptical/core/math/triangle.hpp>
#include <goptical/core/math/vector_pair.hpp>

#include <goptical/core/light/ray.hpp>
#include <goptical/core/light/spectral_line.hpp>

namespace goptical
{

namespace io
{

enum
{
  RANDOM_MAX = 0x7fffffff
};
/* https://en.wikipedia.org/wiki/Lehmer_random_number_generator */
/* because we want deterministic sequence of random numbers */
uint32_t
lcg_parkmiller (uint32_t *state)
{
  return *state = (uint64_t)*state * 48271 % 0x7fffffff;
}

/* NOTE not threadsafe */
static uint32_t seed = 42;
static inline uint32_t
random ()
{
  return lcg_parkmiller (&seed);
}

Renderer::Renderer ()
    : _feature_size (20.), _ray_color_mode (RayColorWavelen),
      _intensity_mode (IntensityIgnore)
{
  _styles_color[StyleForeground] = Rgb (1.0, 1.0, 1.0);
  _styles_color[StyleBackground] = Rgb (0.0, 0.0, 0.0);
  _styles_color[StyleRay] = Rgb (1.0, 0.0, 0.0);
  _styles_color[StyleSurface] = Rgb (0.5, 0.5, 1.0);
  _styles_color[StyleGlass] = Rgb (0.8, 0.8, 1.0);
}

void
Renderer::group_begin (const std::string &name)
{
}

void
Renderer::group_end ()
{
}

/**********************************************************************
 * Optical elements drawing
 */

void
Renderer::draw_element_2d (const sys::Element &e, const sys::Element *ref)
{
  group_begin ("element");
  e.draw_2d_e (*this, ref);
  group_end ();
}

void
Renderer::draw_element_3d (const sys::Element &e, const sys::Element *ref)
{
  group_begin ("element");
  e.draw_3d_e (*this, ref);
  group_end ();
}

/**********************************************************************
 * light ray drawing
 */

template <unsigned D, bool draw_lost>
bool
Renderer::draw_traced_ray_recurs (const trace::Ray &ray, double lost_len,
                                  const sys::Element *ref, bool hit_image)
{

  const math::Transform<3> &t1 = ray.get_creator ()->get_transform_to (ref);
  math::VectorPair3 p;
  sys::Element *i_element = 0;

  p[0] = t1.transform (ray.origin ());

  if (ray.is_lost ())
    {
      if (!draw_lost)
        return false;

      p[1] = t1.transform (ray.origin () + ray.direction () * lost_len);
    }
  else
    {
      i_element = &ray.get_intercept_element ();
      const math::Transform<3> &t2 = i_element->get_transform_to (ref);
      p[1] = t2.transform (ray.get_intercept_point ());
    }

  bool done = false;

  for (trace::Ray *r = ray.get_first_child (); r; r = r->get_next_child ())
    done |= draw_traced_ray_recurs<D, false> (*r, lost_len, ref, hit_image);

  if (!done && hit_image && !dynamic_cast<const sys::Image *> (i_element))
    return false;

  switch (D)
    {
    case 2:
      // skip non tangential rays in 2d mode
      if (fabs (p.x1 ()) > 1e-6)
        return false;

      draw_ray_line (
          math::VectorPair2 (p[0].project_zy (), p[1].project_zy ()), ray);
      break;

    case 3:
      draw_ray_line (p, ray);
      break;
    }

  return true;
}

template <unsigned D>
void
Renderer::draw_trace_result (const trace::Result &result,
                             const sys::Element *ref, bool hit_image)
{
  const trace::Result::sources_t &sl = result.get_source_list ();
  double lost_len = result.get_params ().get_lost_ray_length ();

  if (sl.empty ())
    throw Error ("No source found in trace result");

  _max_intensity = result.get_max_ray_intensity ();

  for (auto &s : sl)
    {
      try
        {
          const trace::rays_queue_t &rl
              = result.get_generated (*(sys::Element *)s);

          for (auto &r : rl)
            {
              group_begin ("ray");
              draw_traced_ray_recurs<D, false> (*r, lost_len, ref, hit_image);
              group_end ();
            }
        }
      catch (Error &e)
        {
        }
    }
}

void
Renderer::draw_trace_result_2d (const trace::Result &result, bool hit_image,
                                const sys::Element *ref)
{
  group_begin ("rays");
  draw_trace_result<2> (result, ref, hit_image);
  group_end ();
}

void
Renderer::draw_trace_result_3d (const trace::Result &result, bool hit_image,
                                const sys::Element *ref)
{
  group_begin ("rays");
  draw_trace_result<3> (result, ref, hit_image);
  group_end ();
}

const int width=300, height=300;
double minX=INFINITY, maxX=-INFINITY, minY=INFINITY, maxY=-INFINITY;
static int ignoredPoints=0;
static std::vector<double> img;
void plotPoint(const double x, const double y, Rgb const& color)
{
//    std::cerr << "x,y: " << x << ", " << y << " => rgb: "
//              << color.r << "," << color.g << "," << color.b << "\n";
    const auto xi = std::round((x-minX)/(maxX-minX)*(width-1));
    const auto yi = std::round((y-minY)/(maxY-minY)*(height-1));
    if(xi<0 || xi>=width || yi<0 || yi>=height)
    {
        ++ignoredPoints;
        return;
    }
    img[(xi+height*yi)*3+0]+=color.r;
    img[(xi+height*yi)*3+1]+=color.g;
    img[(xi+height*yi)*3+2]+=color.b;
}

#pragma pack(push,1)
struct BitmapHeader
{
    uint16_t signature; // 0x4d42
    uint32_t fileSize;
    uint32_t zero;
    uint32_t dataOffset;
    uint32_t bitmapInfoHeaderSize; // 40
    uint32_t width;
    uint32_t height;
    uint16_t numOfPlanes; // 1
    uint16_t bpp;
    uint32_t compressionType; // none is 0
    uint32_t dataSize;
    uint32_t horizPPM;
    uint32_t vertPPM;
    uint32_t numOfColors;
    uint32_t numOfImportantColors;
};
#pragma pack(pop)

void
Renderer::draw_intercepts (const trace::Result &result, const sys::Surface &s)
{
  _max_intensity = result.get_max_ray_intensity ();

  img.clear();
  img.resize(width*height*3);
  ignoredPoints=0;

  // Find the domain extents
  minX=minY=INFINITY;
  maxX=maxY=-INFINITY;
  for (auto &i : result.get_intercepted (s))
  {
      const trace::Ray &ray = *i;
      // dont need global transform here, draw ray intercept points in
      // surface local coordinates.
      const auto pos=ray.get_intercept_point().project_xy();
      if(pos.x()<minX) minX=pos.x();
      if(pos.y()<minY) minY=pos.y();
      if(pos.x()>maxX) maxX=pos.x();
      if(pos.y()>maxY) maxY=pos.y();
  }

  // Make sure the domain being plotted is a square that contains all the points
  const auto deltaX=maxX-minX, deltaY=maxY-minY;
  if(deltaX>deltaY)
  {
      maxY+=(deltaX-deltaY)/2;
      minY-=(deltaX-deltaY)/2;
  }
  else if(deltaX<deltaY)
  {
      maxX+=(deltaY-deltaX)/2;
      minX-=(deltaY-deltaX)/2;
  }

  // Now we can do the plotting
  for (auto &i : result.get_intercepted (s))
    {
      const trace::Ray &ray = *i;
      // dont need global transform here, draw ray intercept points in
      // surface local coordinates.
      draw_point (ray.get_intercept_point ().project_xy (), ray_to_rgb (ray));
      const auto pos=ray.get_intercept_point().project_xy();
      plotPoint(pos.x(), pos.y(), ray_to_rgb(ray));
    }

  // And finally save the BMP
  {
      static int seqNum=0;
      ++seqNum;
      std::ofstream file("/tmp/plot"+std::to_string(seqNum)+".bmp");
      BitmapHeader header={};
      header.signature=0x4d42;
      header.fileSize=((width+3)&~3)*height*3+sizeof header;
      header.dataOffset=sizeof header;
      header.bitmapInfoHeaderSize=40;
      header.width=width;
      header.height=height;
      header.numOfPlanes=1;
      header.bpp=24;
      file.write(reinterpret_cast<const char*>(&header), sizeof header);

      const auto maxImgVal=*std::max_element(img.begin(), img.end());
      std::cerr << "Min x,y: " << minX << "," << minY << "; max x,y: " << maxX << "," << maxY << "\n";
      std::cerr << "Ignored points: " << ignoredPoints << "\n";
      std::cerr << "Max image color channel value: " << maxImgVal << "\n";
      for(int yi=0; yi<height; ++yi)
      {
          for(int xi=0; xi<width; ++xi)
          {
              const uint8_t red  = std::pow(std::max(0., img[(xi+height*yi)*3+0]/maxImgVal), 1/2.2)*255;
              const uint8_t green= std::pow(std::max(0., img[(xi+height*yi)*3+1]/maxImgVal), 1/2.2)*255;
              const uint8_t blue = std::pow(std::max(0., img[(xi+height*yi)*3+2]/maxImgVal), 1/2.2)*255;
              const uint8_t bgr[3]={blue, green, red};
              file.write(reinterpret_cast<const char*>(bgr), sizeof bgr);
          }
          constexpr auto scanLineAlignment=4;
          const char align[scanLineAlignment-1]={};
          const auto alignSize=(sizeof header-file.tellp())%scanLineAlignment;
          if(alignSize)
              file.write(align,alignSize);
      }
  }
}

const Rgb
Renderer::ray_to_rgb (const light::Ray &ray)
{
  // FIXMEG
  float r = static_cast<float> (random ()) / static_cast<float> (RANDOM_MAX);
  float g = static_cast<float> (random ()) / static_cast<float> (RANDOM_MAX);
  float b = static_cast<float> (random ()) / static_cast<float> (RANDOM_MAX);
  return Rgb (r, g, b);

  switch (_ray_color_mode)
    {
    case RayColorWavelen:
      return light::SpectralLine::get_wavelen_color (ray.get_wavelen ());

    default:
    case RayColorFixed:
      return get_style_color (StyleRay);
    }
}

float
Renderer::ray_to_alpha (const light::Ray &ray) const
{
  switch (_intensity_mode)
    {
    case IntensityIgnore:
      return 0.0;

    case IntensityShade:
      return 1.0 - std::min (ray.get_intensity () / _max_intensity, 1.0);

    case IntensityLogShade: // FIXME add log
      return 1.0 - std::min (ray.get_intensity () / _max_intensity, 1.0);
    }

  return 0;
}

void
Renderer::draw_ray_line (const math::VectorPair2 &l, const trace::Ray &ray)
{
  draw_segment (l, ray_to_rgb (ray));
}

void
Renderer::draw_ray_line (const math::VectorPair3 &l, const trace::Ray &ray)
{
  draw_segment (l, ray_to_rgb (ray));
}

/**********************************************************************
 * Misc shapes 2d drawing
 */

void
Renderer::draw_polygon (const math::Vector2 *array, unsigned int count,
                        const Rgb &rgb, bool filled, bool closed)
{
  unsigned int i;

  if (count < 3)
    return;

  for (i = 0; i + 1 < count; i++)
    draw_segment (math::VectorPair2 (array[i], array[i + 1]), rgb);

  if (closed)
    draw_segment (math::VectorPair2 (array[i], array[0]), rgb);
}

void
Renderer::draw_circle (const math::Vector2 &v, double r, const Rgb &rgb,
                       bool filled)
{
  unsigned int count
      = std::min (100, std::max (6, (int)(2. * M_PI * r / _feature_size)));
  DPP_VLARRAY (math::Vector2, count, p);
  double astep = 2. * M_PI / count;
  double a = astep;
  p[0] = math::Vector2 (r, 0);

  for (unsigned int i = 0; i < count; i++, a += astep)
    p[i] = v + math::Vector2 (r * cos (a), r * sin (a));

  draw_polygon (&p[0], count, rgb, filled, true);
}

void
Renderer::draw_triangle (const math::Triangle<3> &t, const Rgb &rgb)
{
  draw_polygon (&t[0], 3, rgb, false, true);
}

void
Renderer::draw_triangle (const math::Triangle<3> &t,
                         const math::Triangle<3> &gradient, const Rgb &rgb)
{
  draw_triangle (t, rgb);
}

void
Renderer::draw_polygon (const math::Vector3 *array, unsigned int count,
                        const Rgb &rgb, bool filled, bool closed)
{
  if (count < 3)
    return;

  unsigned int i;

  for (i = 0; i + 1 < count; i++)
    draw_segment (array[i], array[i + 1], rgb);

  if (closed)
    draw_segment (array[i], array[0], rgb);
}

void
Renderer::draw_box (const math::VectorPair2 &c, const Rgb &rgb)
{
  draw_segment (math::Vector2 (c[0].x (), c[0].y ()),
                math::Vector2 (c[1].x (), c[0].y ()), rgb);
  draw_segment (math::Vector2 (c[1].x (), c[1].y ()),
                math::Vector2 (c[1].x (), c[0].y ()), rgb);
  draw_segment (math::Vector2 (c[1].x (), c[1].y ()),
                math::Vector2 (c[0].x (), c[1].y ()), rgb);
  draw_segment (math::Vector2 (c[0].x (), c[0].y ()),
                math::Vector2 (c[0].x (), c[1].y ()), rgb);
}

void
Renderer::draw_triangle (const math::Triangle<2> &t, bool filled,
                         const Rgb &rgb)
{
  draw_polygon (&t[0], 3, rgb, filled, true);
}

}

}
