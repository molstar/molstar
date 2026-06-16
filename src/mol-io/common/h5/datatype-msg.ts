/**
 * Vendored from jsfive v0.4.0 (https://github.com/usnistgov/jsfive)
 *   "jsfive is in the public domain.
 *    It is based in large part on the pyfive library
 *    https://github.com/jjhelmus/pyfive
 *    Copyright (c) 2016 Jonathan J. Helmus. All rights reserved."
 *
 * Minimal TypeScript adaptation: relative `.js` import suffixes stripped,
 * file extensions changed to `.ts`. Behavior unchanged.
 *
 * @author Jonathan J. Helmus (pyfive original)
 * @author Brian Maranville (jsfive)
 */

/* eslint-disable */
// @ts-nocheck

import {_structure_size, _unpack_struct_from} from './core';

export class DatatypeMessage {
  //""" Representation of a HDF5 Datatype Message. """
  //# Contents and layout defined in IV.A.2.d.

  constructor(buf, offset) {
    this.buf = buf;
    this.offset = offset;
    this.dtype = this.determine_dtype();
  }
  determine_dtype() {
    //""" Return the dtype (often numpy-like) for the datatype message.  """
    let datatype_msg = _unpack_struct_from(DATATYPE_MSG, this.buf, this.offset);
    this.offset += DATATYPE_MSG_SIZE
    //# last 4 bits
    let datatype_class = datatype_msg.get('class_and_version') & 0x0F;

    if (datatype_class == DATATYPE_FIXED_POINT) {
      return this._determine_dtype_fixed_point(datatype_msg);
    }
    else if (datatype_class == DATATYPE_FLOATING_POINT) {
      return this._determine_dtype_floating_point(datatype_msg)
    }
    else if (datatype_class == DATATYPE_TIME) {
      throw "Time datatype class not supported."
    }
    else if (datatype_class == DATATYPE_STRING) {
      return this._determine_dtype_string(datatype_msg);
    }
    else if (datatype_class == DATATYPE_BITFIELD) {
      throw "Bitfield datatype class not supported."
    }
    else if (datatype_class == DATATYPE_OPAQUE) {
      throw "Opaque datatype class not supported."
    }
    else if (datatype_class == DATATYPE_COMPOUND) {
      return this._determine_dtype_compound(datatype_msg);
    }
    else if (datatype_class == DATATYPE_REFERENCE) {
      return ['REFERENCE', datatype_msg.get('size')];
    }
    else if (datatype_class == DATATYPE_ENUMERATED) {
      // enumerated base class datatype message starts at end of
      // enum datatype message, and offset is already advanced above,
      // so just run the same function again to get base class:
      return this.determine_dtype();
    }
    else if (datatype_class == DATATYPE_ARRAY) {
      throw "Array datatype class not supported."
    }
    else if (datatype_class == DATATYPE_VARIABLE_LENGTH) {
      let vlen_type = this._determine_dtype_vlen(datatype_msg);
      if (vlen_type[0] == 'VLEN_SEQUENCE') {
        let base_type = this.determine_dtype();
        vlen_type = ['VLEN_SEQUENCE', base_type];
      }
      return vlen_type
    }
    else {
      throw 'Invalid datatype class ' + datatype_class;
    }
  }

  _determine_dtype_fixed_point(datatype_msg) {
    //""" Return the NumPy dtype for a fixed point class. """
    //# fixed-point types are assumed to follow IEEE standard format
    let length_in_bytes = datatype_msg.get('size');
    if (!([1, 2, 4, 8].includes(length_in_bytes))) {
        throw "Unsupported datatype size";
    }

    let signed = datatype_msg.get('class_bit_field_0') & 0x08;
    var dtype_char;
    if (signed > 0) {
      dtype_char = 'i';
    }
    else {
      dtype_char = 'u';
    }

    let byte_order = datatype_msg.get('class_bit_field_0') & 0x01;
    var byte_order_char;
    if (byte_order == 0) {
      byte_order_char = '<';  //# little-endian
    }
    else {
      byte_order_char = '>';  //# big-endian
    }

    //# 4-byte fixed-point property description
    //# not read, assumed to be IEEE standard format
    this.offset += 4

    return byte_order_char + dtype_char + length_in_bytes.toFixed();
  }

  _determine_dtype_floating_point(datatype_msg) {
    //""" Return the NumPy dtype for a floating point class. """
    //# Floating point types are assumed to follow IEEE standard formats
    let length_in_bytes = datatype_msg.get('size');
    if (!([1, 2, 4, 8].includes(length_in_bytes))) {
      throw "Unsupported datatype size";
    }

    let dtype_char = 'f'

    let byte_order = datatype_msg.get('class_bit_field_0') & 0x01;
    var byte_order_char;
    if (byte_order == 0) {
      byte_order_char = '<';  //# little-endian
    }
    else {
      byte_order_char = '>';  //# big-endian
    }

    //# 12-bytes floating-point property description
    //# not read, assumed to be IEEE standard format
    this.offset += 12

    return byte_order_char + dtype_char + length_in_bytes.toFixed();
  }
  
  _determine_dtype_string(datatype_msg) {
    //""" Return the NumPy dtype for a string class. """
    return 'S' + datatype_msg.get('size').toFixed();
  }
  
  _determine_dtype_vlen(datatype_msg) {
    //""" Return the dtype information for a variable length class. """
    let vlen_type = datatype_msg.get('class_bit_field_0') & 0x01;
    if (vlen_type != 1) {
      return ['VLEN_SEQUENCE', 0, 0];
    }
    let padding_type = datatype_msg.get('class_bit_field_0') >> 4;  //# bits 4-7
    let character_set = datatype_msg.get('class_bit_field_1') & 0x01; 
    return ['VLEN_STRING', padding_type, character_set];
  }
  _determine_dtype_compound(datatype_msg) {
    throw "Compound type not yet implemented!";
  }
}   


var DATATYPE_MSG = new Map([
  ['class_and_version', 'B'],
  ['class_bit_field_0', 'B'],
  ['class_bit_field_1', 'B'],
  ['class_bit_field_2', 'B'],
  ['size', 'I'],
]);
var DATATYPE_MSG_SIZE = _structure_size(DATATYPE_MSG);


var COMPOUND_PROP_DESC_V1 = new Map([
    ['offset', 'I'],
    ['dimensionality', 'B'],
    ['reserved_0', 'B'],
    ['reserved_1', 'B'],
    ['reserved_2', 'B'],
    ['permutation', 'I'],
    ['reserved_3', 'I'],
    ['dim_size_1', 'I'],
    ['dim_size_2', 'I'],
    ['dim_size_3', 'I'],
    ['dim_size_4', 'I'],
]);
var COMPOUND_PROP_DESC_V1_SIZE = _structure_size(COMPOUND_PROP_DESC_V1);


//# Datatype message, datatype classes
var DATATYPE_FIXED_POINT = 0;
var DATATYPE_FLOATING_POINT = 1;
var DATATYPE_TIME = 2;
var DATATYPE_STRING = 3;
var DATATYPE_BITFIELD = 4;
var DATATYPE_OPAQUE = 5;
var DATATYPE_COMPOUND = 6;
var DATATYPE_REFERENCE = 7;
var DATATYPE_ENUMERATED = 8;
var DATATYPE_VARIABLE_LENGTH = 9;
var DATATYPE_ARRAY = 10;