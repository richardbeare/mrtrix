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

#include "app.h"
#include "image/position.h"

using namespace std; 
using namespace MR; 

SET_VERSION_DEFAULT;

DESCRIPTION = {
  "average an image along a specific axis.",
  NULL
};

ARGUMENTS = {
  Argument ("image1", "first input image", "the first input image.").type_image_in (),
  Argument ("output", "output image", "the output image.").type_image_out (),
  Argument::End
};

OPTIONS = { 
  Option ("axis", "average axis", "specify axis along which averaging should be performed. By default, the program will use the last non-singleton axis of the input image.")
    .append (Argument ("axis", "axis", "the concatenation axis").type_integer (0, INT_MAX, 0)),

  Option ("geometric", "geometric mean", "produce geometric mean. By default, the program will produce the arithmetic mean."),

  Option::End 
};




EXECUTE {
  int axis = -1;

  std::vector<OptBase> opt = get_options (0); // axis
  if (opt.size()) axis = opt[0][0].get_int();

  Image::Object &in_obj (*argument[0].get_image());
  Image::Header header (in_obj);

  bool geometric = false;
  opt = get_options (1); // geometric
  if (opt.size()) {
    if (in_obj.is_complex()) 
      throw Exception ("geometric mean not supported for complex data.");
    geometric = true;
  }

  int lastdim;
  for (lastdim = header.axes.ndim()-1; lastdim > 0 && header.axes.dim[lastdim] <= 1; lastdim--);

  if (axis < 0) axis = lastdim;
  else if (axis >= lastdim) 
    throw Exception ("averaging along singleton dimension");

  info ("averaging along axis " + str (axis));

  header.axes.set_ndim (lastdim+1);
  header.axes.dim[axis] = 1;
  if (in_obj.is_complex()) header.data_type = DataType::CFloat32;
  else header.data_type = DataType::Float32;

  Image::Position in (in_obj);
  Image::Position out (*argument[1].get_image (header));


  float norm = 1.0 / in.dim (axis);
  
  ProgressBar::init (out.voxel_count(), "averaging...");

  do {
    float re = 0.0, im = 0.0;

    for (in.set (3,0); in[3] < in.dim(3); in.inc(3)) {
      if (geometric) re += log (in.re());
      else {
        re += in.re();
        if (in.is_complex()) im += in.im();
      }
    }

    re *= norm;

    if (geometric) re = exp (re);
    out.re (re);
    if (in.is_complex()) out.im (norm*im);

    in++;
    ProgressBar::inc();
  } while (out++);

  ProgressBar::done();
}
