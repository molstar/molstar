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

import { inflateSync } from '../../../mol-util/zip/zip';
import { struct } from './core';

const zlib_decompress = function (buf, itemsize) {
  let input_array = new Uint8Array(buf);
  let out = inflateSync(input_array);
  // Return an ArrayBuffer with the inflated bytes only.
  return out.buffer.slice(out.byteOffset, out.byteOffset + out.byteLength);
}

const unshuffle = function (buf, itemsize) {
  let buffer_size = buf.byteLength;
  let unshuffled_view = new Uint8Array(buffer_size);
  let step = Math.floor(buffer_size / itemsize);
  let shuffled_view = new DataView(buf);
  for (var j = 0; j < itemsize; j++) {
    for (var i = 0; i < step; i++) {
      unshuffled_view[j + i * itemsize] = shuffled_view.getUint8(j * step + i);
    }
  }
  return unshuffled_view.buffer;
}

const fletch32 = function (buf, itemsize) {
  _verify_fletcher32(buf);
  //# strip off 4-byte checksum from end of buffer
  return buf.slice(0, -4);
}

function _verify_fletcher32(chunk_buffer) {
  //""" Verify a chunk with a fletcher32 checksum. """
  //# calculate checksums
  var odd_chunk_buffer = ((chunk_buffer.byteLength % 2) != 0);
  var data_length = chunk_buffer.byteLength - 4;
  var view = new DataView(chunk_buffer);

  var sum1 = 0;
  var sum2 = 0;
  for (var offset=0; offset<(data_length-1); offset+=2) {
    let datum = view.getUint16(offset, true); // little-endian
    sum1 = (sum1 + datum) % 65535
    sum2 = (sum2 + sum1) % 65535
  }
  if (odd_chunk_buffer) {
    // process the last item:
    let datum = view.getUint8(data_length-1);
    sum1 = (sum1 + datum) % 65535
    sum2 = (sum2 + sum1) % 65535
  }

  //# extract stored checksums
  var [ref_sum1, ref_sum2] = struct.unpack_from('>HH', chunk_buffer, data_length); // .fromstring(chunk_buffer[-4:], '>u2')
  ref_sum1 = ref_sum1 % 65535
  ref_sum2 = ref_sum2 % 65535

  //# compare
  if (sum1 != ref_sum1 || sum2 != ref_sum2) {
    throw 'ValueError("fletcher32 checksum invalid")';
  }
  return true
}

//# IV.A.2.l The Data Storage - Filter Pipeline message
var RESERVED_FILTER = 0;
const GZIP_DEFLATE_FILTER = 1;
const SHUFFLE_FILTER = 2;
const FLETCH32_FILTER = 3;
var SZIP_FILTER = 4;
var NBIT_FILTER = 5;
var SCALEOFFSET_FILTER = 6;

// To register a new filter, add a function (ArrayBuffer) => ArrayBuffer
// the the following map, using a key that corresponds to filter_id (int)
export const Filters = new Map([
  [GZIP_DEFLATE_FILTER, zlib_decompress],
  [SHUFFLE_FILTER, unshuffle],
  [FLETCH32_FILTER, fletch32]
]);