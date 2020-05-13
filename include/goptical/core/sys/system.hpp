/*

      This file is part of the Goptical Core library.
  
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

#ifndef GOPTICAL_SYSTEM_HH_
#define GOPTICAL_SYSTEM_HH_

#include <iostream>

#include "goptical/core/common.hpp"

#include "goptical/core/sys/element.hpp"
#include "goptical/core/sys/container.hpp"
#include "goptical/core/trace/params.hpp"
#include "goptical/core/material/proxy.hpp"

namespace _goptical {

  namespace sys {

    class SystemBuilder;
    /**
       @short Optical system
       @header <goptical/core/sys/system
       @module {Core}
       @main

       This class is used to describe an optical system. Any optical
       @ref Element {element} may be part of a system. This class handle 3d
       transformation between elements local coordinates.

       @xsee {tuto_system}
    */
    class System : public Container
    {
      friend class Element;

    public:
      /** Create a new empty system. */
      System();
      virtual ~System();

      System(const System&) = delete;
      System& operator=(const System&) = delete;

      /** Define an entrance pupil surface used to project source rays */
      inline void set_entrance_pupil(const std::shared_ptr<Surface> &entrance);
      /** Discard defined entrance pupil */
      inline void undef_entrance_pupil();
      /** Get defined entrance pupil surface or try to guess it if none defined */
      const Surface &get_entrance_pupil() const;
      /** Test if an entrance pupil has been defined */
      inline bool has_entrance_pupil() const;

      /** Define an exit pupil surface */
      inline void set_exit_pupil(const std::shared_ptr<Surface> &exit);
      /** Get exit pupil */
      inline const Surface &get_exit_pupil() const;
      /** Test if an exit pupil has been defined */
      inline bool has_exit_pupil() const;

      /** Get default tracer parameters */
      inline const trace::Params & get_tracer_params() const;
      /** Get default tracer parameters */
      inline trace::Params & get_tracer_params();

      /** Get transform between two elements local coordinates */
      inline const math::Transform<3> & get_transform(const Element &from, const Element &to) const;

      /** Get transform from element local to global coordinates */
      inline const math::Transform<3> & get_global_transform(const Element &from) const;

      /** Get transform from global to element local coordinates */
      inline const math::Transform<3> & get_local_transform(const Element &to) const;

      /** Get system version. version is updated each time system or
          associated elements properties are changed */
      inline unsigned int get_version() const;

      /** Get the number of registered elements in the system */
      inline unsigned int get_element_count() const;

      /** Get registered element. first element has index 1 */
      inline Element & get_element(unsigned int index) const;

      /** Increase current system version */
      inline void update_version();

      /** Find surface which colides with the given ray and update intersection point */
      Surface * colide_next(const trace::Params &params,
                            math::VectorPair3 &intersect,
                            const trace::Ray &ray) const;

      /** set environment material */
      void set_environment(const std::shared_ptr<material::Base> &env);

      /** get environment material */
      inline const std::shared_ptr<material::Base> & get_environment() const;

      /** @internal get environment material proxy */
      inline const std::shared_ptr<material::Base> & get_environment_proxy() const;

      /** @internal Dump 3d transforms cache */
      void transform_cache_dump(std::ostream &o) const;

      void add (std::shared_ptr<Element> e);


    private:

      /** get an new element identifier */
      unsigned int index_get(Element &element);
      /** free the identifier associated with the given element, and associated transformations */
      void index_put(const Element &element);

      /** Get a reference to cache entry for transform between 2 elements (ids) */
      inline math::Transform<3> * & transform_cache_entry(unsigned int from, unsigned int to) const;
      /* Delete transformation from / to */
      void transform_cache_delete(unsigned int from, unsigned int to);

      /** Compute and get 3d transform between element local and global coordinates */
      const math::Transform<3> & transform_l2g_cache_update(const Element &e) const;
      /** Compute and get 3d transform between element local and global coordinates */
      const math::Transform<3> & transform_g2l_cache_update(const Element &e) const;
      /** Compute and get 3d transform between two elements local coordinates */
      const math::Transform<3> & transform_cache_update(const Element &from, const Element &to) const;

      /** Flush all cached transforms associated with a given element */
      void transform_cache_flush(const Element &element);
      /** Flush all cached transforms */
      void transform_cache_flush();

      /** Resize transform cache size */
      void transform_cache_resize(unsigned int newsize);

      unsigned int              _version;

      std::shared_ptr<Surface>        _entrance;
      std::shared_ptr<Surface>        _exit;
      std::shared_ptr<material::Base>   _env_proxy;
      trace::Params             _tracer_params;
      // Count of number of elements in the system
      // This includes all child elements groups
      unsigned int              _e_count;
      // _index_map is used to assign ids to elements
      std::vector<Element *>    _index_map;
      // Holds cache entries in two dimensions _e_count by _e_count
      // transforming from element to element using element id as row and
      // column indices
      std::vector<math::Transform<3> *> _transform_cache;
    };
  }
}

#endif

