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


    21-07-2008 J-Donald Tournier <d.tournier@brain.org.au>
    * add definitions for all static const declarations to avoid linking errors.
  
*/

#include "data_type.h"

namespace MR {

  const guchar DataType::ComplexNumber;
  const guchar DataType::Signed;
  const guchar DataType::LittleEndian;
  const guchar DataType::BigEndian;
  const guchar DataType::Text;
  const guchar DataType::GroupStart;
  const guchar DataType::GroupEnd;
  const guchar DataType::Undefined;

  const guchar DataType::Bit;
  const guchar DataType::UInt8;
  const guchar DataType::UInt16;
  const guchar DataType::UInt32;
  const guchar DataType::Float32;
  const guchar DataType::Float64;
  const guchar DataType::Int8;
  const guchar DataType::Int16;
  const guchar DataType::Int16LE;
  const guchar DataType::UInt16LE;
  const guchar DataType::Int16BE;
  const guchar DataType::UInt16BE;
  const guchar DataType::Int32;
  const guchar DataType::Int32LE;
  const guchar DataType::UInt32LE;
  const guchar DataType::Int32BE;
  const guchar DataType::UInt32BE;
  const guchar DataType::Float32LE;
  const guchar DataType::Float32BE;
  const guchar DataType::Float64LE;
  const guchar DataType::Float64BE;
  const guchar DataType::CFloat32;
  const guchar DataType::CFloat32LE;
  const guchar DataType::CFloat32BE;
  const guchar DataType::CFloat64;
  const guchar DataType::CFloat64LE;
  const guchar DataType::CFloat64BE;
  const guchar DataType::Native;





  void DataType::parse (const String& spec)
  {
    String str (lowercase (spec));

    if (str == "float32")    { dt = Float32;     return; }
    if (str == "float32le")  { dt = Float32LE;   return; }
    if (str == "float32be")  { dt = Float32BE;   return; }
 
    if (str == "float64")    { dt = Float64;     return; }
    if (str == "float64le")  { dt = Float64LE;   return; }
    if (str == "float64be")  { dt = Float64BE;   return; }

    if (str == "int32")      { dt = Int32;       return; }
    if (str == "uint32")     { dt = UInt32;      return; }
    if (str == "int32le")    { dt = Int32LE;     return; }
    if (str == "uint32le")   { dt = UInt32LE;    return; }
    if (str == "int32be")    { dt = Int32BE;     return; }
    if (str == "uint32be")   { dt = UInt32BE;    return; }

    if (str == "int16")      { dt = Int16;       return; }
    if (str == "uint16")     { dt = UInt16;      return; }
    if (str == "int16le")    { dt = Int16LE;     return; }
    if (str == "uint16le")   { dt = UInt16LE;    return; }
    if (str == "int16be")    { dt = Int16BE;     return; }
    if (str == "uint16be")   { dt = UInt16BE;    return; }

    if (str == "cfloat32")   { dt = CFloat32;    return; }
    if (str == "cfloat32le") { dt = CFloat32LE;  return; }
    if (str == "cfloat32be") { dt = CFloat32BE;  return; }

    if (str == "cfloat64")   { dt = CFloat64;    return; }
    if (str == "cfloat64le") { dt = CFloat64LE;  return; }
    if (str == "cfloat64be") { dt = CFloat64BE;  return; }

    if (str == "int8")       { dt = Int8;        return; }
    if (str == "uint8")      { dt = UInt8;       return; }
 
    if (str == "bit")        { dt = Bit;         return; }

    throw Exception ("invalid data type \"" + spec + "\"");
  }




  guint DataType::bits () const
  {
    switch (dt) {
      case Bit:
        return (1);
      case Int8:
      case UInt8:
        return (8*sizeof (gint8));
      case Int16:
      case UInt16:
      case Int16LE:
      case UInt16LE:
      case Int16BE:
      case UInt16BE:
        return (8*sizeof (gint16));
      case Int32:
      case UInt32:
      case Int32LE:
      case UInt32LE:
      case Int32BE:
      case UInt32BE:
        return (8*sizeof (gint32));
      case Float32:
      case Float32LE:
      case Float32BE:
        return (8*sizeof (float32));
      case Float64:
      case Float64LE:
      case Float64BE:
        return (8*sizeof (float64));
      case CFloat32:
      case CFloat32LE:
      case CFloat32BE:
        return (16*sizeof (float32));
      case CFloat64:
      case CFloat64LE:
      case CFloat64BE:
        return (16*sizeof (float64));
    }
    
    error ("invalid datatype specifier");
    return(0);
  }


  const gchar* DataType::description() const
  {
    switch (dt) {
      case Bit:        return ("bitwise");

      case Int8:       return ("signed 8 bit integer");
      case UInt8:      return ("unsigned 8 bit integer");

      case Int16LE:    return ("signed 16 bit integer (little endian)");
      case UInt16LE:   return ("unsigned 16 bit integer (little endian)");
      case Int16BE:    return ("signed 16 bit integer (big endian)");
      case UInt16BE:   return ("unsigned 16 bit integer (big endian)");

      case Int32LE:    return ("signed 32 bit integer (little endian)");
      case UInt32LE:   return ("unsigned 32 bit integer (little endian)");
      case Int32BE:    return ("signed 32 bit integer (big endian)");
      case UInt32BE:   return ("unsigned 32 bit integer (big endian)");

      case Float32LE:  return ("32 bit float (little endian)");
      case Float32BE:  return ("32 bit float (big endian)");

      case Float64LE:  return ("64 bit float (little endian)");
      case Float64BE:  return ("64 bit float (big endian)");

      case CFloat32LE: return ("Complex 32 bit float (little endian)");
      case CFloat32BE: return ("Complex 32 bit float (big endian)");

      case CFloat64LE: return ("Complex 64 bit float (little endian)");
      case CFloat64BE: return ("Complex 64 bit float (big endian)");

      case Undefined:  return ("undefined");
      case Text:       return ("string");
      case GroupStart: return ("start of group");
      case GroupEnd:   return ("end of group");

      default:      return (NULL);
    }

    return (NULL);
  }



  const gchar* DataType::specifier() const
  {
    switch (dt) {
      case Bit:        return ("Bit");

      case Int8:       return ("Int8");
      case UInt8:      return ("UInt8");

      case Int16LE:    return ("Int16LE");
      case UInt16LE:   return ("UInt16LE");
      case Int16BE:    return ("Int16BE");
      case UInt16BE:   return ("UInt16BE");

      case Int32LE:    return ("Int32LE");
      case UInt32LE:   return ("UInt32LE");
      case Int32BE:    return ("Int32BE");
      case UInt32BE:   return ("UInt32BE");

      case Float32LE:  return ("Float32LE");
      case Float32BE:  return ("Float32BE");

      case Float64LE:  return ("Float64LE");
      case Float64BE:  return ("Float64BE");

      case CFloat32LE: return ("CFloat32LE");
      case CFloat32BE: return ("CFloat32BE");

      case CFloat64LE: return ("CFloat64LE");
      case CFloat64BE: return ("CFloat64BE");

      case Int16:      return ("Int16");
      case UInt16:     return ("UInt16");
      case Int32:      return ("Int32");
      case UInt32:     return ("UInt32");
      case Float32:    return ("Float32");
      case Float64:    return ("Float64");
      case CFloat32:   return ("CFloat32");
      case CFloat64:   return ("CFloat64");

      case Undefined:  return ("Undefined");
      case Text:       return ("String");
      case GroupStart: return ("GroupStart");
      case GroupEnd:   return ("GroupEnd");

      default:         return ("invalid");
    }

    return (NULL);
  }


}

