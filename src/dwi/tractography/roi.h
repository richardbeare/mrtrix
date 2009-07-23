/*
    Copyright 2008 Brain Research Institute, Melbourne, Australia

    Written by J-Donald Tournier, 27/06/08.

    This file is part of MRtrix.

    MRtrix is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MRtrix is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with MRtrix.  If not, see <http://www.gnu.org/licenses/>.

*/

#ifndef __dwi_tractography_roi_h__
#define __dwi_tractography_roi_h__

#include "point.h"
#include "ptr.h"
#include "image/interp.h"
#include "math/simulation.h"


namespace MR {
  namespace DWI {
    namespace Tractography {

      class ROI 
      {
        public:
          typedef enum {
            Seed = 0,
            Include = 1,
            Exclude = 2,
            Mask = 3,
            Undefined = G_MAXINT
          } Type;

          ROI (Type type_id, const Point& sphere_pos, float sphere_radius) : type (type_id), position (sphere_pos), radius (sphere_radius) { }
          ROI (Type type_id, RefPtr<Image::Object> mask_image) : type (type_id), radius (GSL_NAN), mask (mask_image->name()), mask_object (mask_image) { }
          ROI (Type type_id, const String& spec) : type (type_id), radius (GSL_NAN) {
            try {
              Exception::Lower s (1);
              std::vector<float> F (parse_floats (spec));
              if (F.size() != 4) throw 1;
              position.set (F[0], F[1], F[2]);
              radius = F[3];
            }
            catch (...) { 
              info ("error parsing spherical ROI specification \"" + spec + "\" - assuming mask image");
              mask = spec; 
            }
          }

          Type   type;
          Point  position;
          float  radius;
          String mask;
          RefPtr<Image::Object> mask_object;

          String  type_description () const {
            switch (type) { 
              case Seed: return ("seed"); 
              case Include: return ("include"); 
              case Exclude: return ("exclude");
              case Mask: return ("mask");
              default: return ("undefined");
            }
          }

          String shape () const { return (mask.size() ? "image" : "sphere"); }
          String parameters () const { return (mask.size() ? mask : str(position[0]) + "," + str(position[1]) + "," + str(position[2]) + "," + str(radius)); }
          String specification () const { return (type_description() + " " + parameters()); }
      };


      inline std::ostream& operator<< (std::ostream& stream, const ROI& roi)
      {
        stream << roi.type_description() << " (" << roi.parameters() << ")";
        return (stream);
      }



    }
  }
}

#endif


