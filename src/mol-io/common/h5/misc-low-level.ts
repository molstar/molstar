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

import { _structure_size, _padded_size, _unpack_struct_from, struct, assert, _unpack_integer, bitSize } from './core';

export class SuperBlock {
  constructor(fh, offset) {
    let version_hint = struct.unpack_from('<B', fh, offset + 8)[0];
    var contents;
    if (version_hint == 0) {
      contents = _unpack_struct_from(SUPERBLOCK_V0, fh, offset);
      this._end_of_sblock = offset + SUPERBLOCK_V0_SIZE;
    }
    else if (version_hint == 2 || version_hint == 3) {
      contents = _unpack_struct_from(SUPERBLOCK_V2_V3, fh, offset);
      this._end_of_sblock = offset + SUPERBLOCK_V2_V3_SIZE;
    } else {
      throw ("unsupported superblock version: " + version_hint.toFixed())
    }
    // verify contents
    if (contents.get('format_signature') != FORMAT_SIGNATURE) {
      throw 'Incorrect file signature: ' + contents.get('format_signature');
    }
    if (contents.get('offset_size') != 8 || contents.get('length_size') != 8) {
      throw 'File uses non-64-bit addressing';
    }
    this.version = contents.get('superblock_version');
    this._contents = contents
    this._root_symbol_table = null;
    this._fh = fh;
  }
  get offset_to_dataobjects() {
    //""" The offset to the data objects collection for the superblock. """
    if (this.version == 0) {
      var sym_table = new SymbolTable(this._fh, this._end_of_sblock, true);
      this._root_symbol_table = sym_table
      return sym_table.group_offset;
    }
    else if (this.version == 2 || this.version == 3) {
      return this._contents.get('root_group_address');
    }
    else {
      throw ("Not implemented version = " + this.version.toFixed());
    }
  }
}

export class Heap {
  /*
  """
  HDF5 local heap.
  """
  */
  constructor(fh, offset) {
    //""" initalize. """

    //fh.seek(offset)
    let local_heap = _unpack_struct_from(LOCAL_HEAP, fh, offset);
    assert(local_heap.get('signature') == 'HEAP');
    assert(local_heap.get('version') == 0);
    let data_offset = local_heap.get('address_of_data_segment');
    let heap_data = fh.slice(data_offset, data_offset + local_heap.get('data_segment_size'));
    local_heap.set('heap_data', heap_data);
    this._contents = local_heap;
    this.data = heap_data;
  }

  get_object_name(offset) {
    //""" Return the name of the object indicated by the given offset. """
    let end = new Uint8Array(this.data).indexOf(0, offset);
    let name_size = end - offset;
    let name = struct.unpack_from('<' + name_size.toFixed() + 's', this.data, offset)[0];
    return name
  }
}

export class SymbolTable {
  /*
  """
  HDF5 Symbol Table.
  """
  */
  constructor(fh, offset, root = false) {
    //""" initialize, root=True for the root group, False otherwise. """
    var node;
    if (root) {
      //# The root symbol table has no Symbol table node header
      //# and contains only a single entry
      node = new Map([['symbols', 1]]);
    }
    else {
      node = _unpack_struct_from(SYMBOL_TABLE_NODE, fh, offset);
      if (node.get('signature') != 'SNOD') { throw "incorrect node type" }
      offset += SYMBOL_TABLE_NODE_SIZE;
    }
    var entries = [];
    var n_symbols = node.get('symbols');
    for (var i = 0; i < n_symbols; i++) {
      entries.push(_unpack_struct_from(SYMBOL_TABLE_ENTRY, fh, offset))
      offset += SYMBOL_TABLE_ENTRY_SIZE;
    }
    if (root) {
      this.group_offset = entries[0].get('object_header_address');
    }
    this.entries = entries
    this._contents = node
  }

  assign_name(heap) {
    //""" Assign link names to all entries in the symbol table. """
    this.entries.forEach(function (entry) {
      let offset = entry.get('link_name_offset');
      let link_name = heap.get_object_name(offset);
      entry.set('link_name', link_name);
    });
  }

  get_links(heap) {
    //""" Return a dictionary of links (dataset/group) and offsets. """
    var links = {}
    this.entries.forEach(function (e) {
      let cache_type = e.get('cache_type');
      let link_name = e.get('link_name');
      if (cache_type == 0 || cache_type == 1) {
        links[link_name] = e.get('object_header_address');
      }
      else if (cache_type == 2) {
        let scratch = e.get('scratch');
        let buf = new ArrayBuffer(4);
        let bufView = new Uint8Array(buf);
        for (var i = 0; i < 4; i++) {
          bufView[i] = scratch.charCodeAt(i);
        }
        let offset = struct.unpack_from('<I', buf, 0)[0];
        links[link_name] = heap.get_object_name(offset);
      }
    });
    return links
  }
}

export class GlobalHeap {
  /*
  HDF5 Global Heap collection.
  */
  constructor(fh, offset) {
    let header = _unpack_struct_from(GLOBAL_HEAP_HEADER, fh, offset);
    offset += GLOBAL_HEAP_HEADER_SIZE;
    //assert(header.get('signature') == 'GCOL');
    //assert(header.get('version') == 1);
    let heap_data_size = header.get('collection_size') - GLOBAL_HEAP_HEADER_SIZE;
    let heap_data = fh.slice(offset, offset + heap_data_size);
    //assert(heap_data.byteLength == heap_data_size); //# check for early end of file

    this.heap_data = heap_data;
    this._header = header;
    this._objects = null;
  }

  get objects() {
    //""" Dictionary of objects in the heap. """
    if (this._objects == null) {
      this._objects = new Map();
      var offset = 0;
      while (offset <= this.heap_data.byteLength - GLOBAL_HEAP_OBJECT_SIZE) {
        let info = _unpack_struct_from(
          GLOBAL_HEAP_OBJECT, this.heap_data, offset);
        if (info.get('object_index') == 0) {
          break
        }
        offset += GLOBAL_HEAP_OBJECT_SIZE;
        let obj_data = this.heap_data.slice(offset, offset + info.get('object_size'));
        this._objects.set(info.get('object_index'), obj_data);
        offset += _padded_size(info.get('object_size'));
      }
    }
    return this._objects
  }
}

export class FractalHeap {
  /*
  HDF5 Fractal Heap.
  */

  constructor(fh, offset) {
    this.fh = fh;
    let header = _unpack_struct_from(FRACTAL_HEAP_HEADER, fh, offset);
    offset += _structure_size(FRACTAL_HEAP_HEADER);
    assert(header.get('signature') == 'FRHP');
    assert(header.get('version') == 0);

    if (header.get('filter_info_size') > 0) {
      throw "Filter info size not supported on FractalHeap";
    }

    if (header.get("btree_address_huge_objects") == UNDEFINED_ADDRESS) {
      header.set("btree_address_huge_objects", null);
    }
    else {
      throw "Huge objects not implemented in FractalHeap";
    }

    if (header.get("root_block_address") == UNDEFINED_ADDRESS) {
      header.set("root_block_address", null);
    }

    let nbits = header.get("log2_maximum_heap_size");
    let block_offset_size = this._min_size_nbits(nbits);
    let h = new Map([
      ['signature', '4s'],
      ['version', 'B'],
      ['heap_header_adddress', 'Q'],
      //['block_offset', `${block_offset_size}s`]
      ['block_offset', `${block_offset_size}B`]
      //this._int_format(block_offset_size)]
    ]);
    this.indirect_block_header = new Map(h); // make shallow copy;
    this.indirect_block_header_size = _structure_size(h);
    if ((header.get("flags") & 2) == 2) {
      h.set('checksum', 'I');
    }
    this.direct_block_header = h;
    this.direct_block_header_size = _structure_size(h);

    let maximum_dblock_size = header.get('maximum_direct_block_size');
    this._managed_object_offset_size = this._min_size_nbits(nbits);
    let value = Math.min(maximum_dblock_size, header.get('max_managed_object_size'));
    this._managed_object_length_size = this._min_size_integer(value);

    let start_block_size = header.get('starting_block_size');
    let table_width = header.get('table_width');
    if (!(start_block_size > 0)) {
      throw "Starting block size == 0 not implemented";
    }

    let log2_maximum_dblock_size = Number(Math.floor(Math.log2(maximum_dblock_size)));
    assert(1n << BigInt(log2_maximum_dblock_size) == maximum_dblock_size);

    let log2_start_block_size = Number(Math.floor(Math.log2(start_block_size)));
    assert(1n << BigInt(log2_start_block_size) == start_block_size);

    this._max_direct_nrows = log2_maximum_dblock_size - log2_start_block_size + 2;

    let log2_table_width = Math.floor(Math.log2(table_width)); // regular number (H, not Q format)
    assert(1 << log2_table_width == table_width);
    this._indirect_nrows_sub = log2_table_width + log2_start_block_size - 1;

    this.header = header;
    this.nobjects = header.get("managed_object_count") + header.get("huge_object_count") + header.get("tiny_object_count");

    let managed = [];
    let root_address = header.get("root_block_address");
    let nrows = 0;
    if (root_address != null) {
      nrows = header.get("indirect_current_rows_count");
    }
    if (nrows > 0) {
      for (let data of this._iter_indirect_block(fh, root_address, nrows)) {
        managed.push(data);
      }
    }
    else {
      let data = this._read_direct_block(fh, root_address, start_block_size);
      managed.push(data);
    }
    let data_size = managed.reduce((p, c) => p + c.byteLength, 0);
    let combined = new Uint8Array(data_size);
    let moffset = 0;
    managed.forEach((m) => {combined.set(new Uint8Array(m), moffset); moffset += m.byteLength});
    this.managed = combined.buffer;
  }

  _read_direct_block(fh, offset, block_size) {
    let data = fh.slice(offset, offset + block_size);
    let header = _unpack_struct_from(this.direct_block_header, data)
    assert(header.get("signature") == "FHDB");
    return data;
  }

  get_data(heapid) {
    let firstbyte = struct.unpack_from('<B', heapid, 0)[0];

    let reserved = firstbyte & 15;  // bit 0-3
    let idtype = (firstbyte >> 4) & 3;  // bit 4-5
    let version = firstbyte >> 6  // bit 6-7
    let data_offset = 1;
    if (idtype == 0) { // managed
      assert(version == 0);
      let nbytes = this._managed_object_offset_size;
      let offset = _unpack_integer(nbytes, heapid, data_offset);
      // add heap offset:
      //offset += this.offset;
      data_offset += nbytes;

      nbytes = this._managed_object_length_size;
      let size = _unpack_integer(nbytes, heapid, data_offset);

      return this.managed.slice(offset, offset + size);
    }
    else if (idtype == 1) { // tiny
      throw "tiny objectID not supported in FractalHeap"
    }
    else if (idtype == 2) { // huge
      throw "huge objectID not supported in FractalHeap"
    }
    else {
      throw "unknown objectID type in FractalHeap"
    }
  }

  _min_size_integer(integer) {
    // """ Calculate the minimal required bytes to contain an integer. """
    return this._min_size_nbits(bitSize(integer));
  }

  _min_size_nbits(nbits) {
    //""" Calculate the minimal required bytes to contain a number of bits. """
    return Math.ceil(nbits / 8);
  }

  * _iter_indirect_block(fh, offset, nrows) {
    let header = _unpack_struct_from(this.indirect_block_header, fh, offset);
    offset += this.indirect_block_header_size;
    assert(header.get("signature") == "FHIB");
    let block_offset_bytes = header.get("block_offset");
    // equivalent to python int.from_bytes with byteorder="little":
    let block_offset = block_offset_bytes.reduce((p, c, i) => p + (c << (i * 8)), 0);
    header.set("block_offset", block_offset);

    let [ndirect, nindirect] = this._indirect_info(nrows);

    let direct_blocks = [];
    for (let i = 0; i < ndirect; i++) {
      let address = struct.unpack_from('<Q', fh, offset)[0]
      offset += 8;
      if (address == UNDEFINED_ADDRESS) {
        break
      }
      let block_size = this._calc_block_size(i);
      direct_blocks.push([address, block_size]);
    }

    let indirect_blocks = [];
    for (let i = ndirect; i < ndirect + nindirect; i++) {
      let address = struct.unpack_from('<Q', fh, offset)[0];
      offset += 8;
      if (address == UNDEFINED_ADDRESS) {
        break
      }
      let block_size = this._calc_block_size(i);
      let nrows = this._iblock_nrows_from_block_size(block_size);
      indirect_blocks.push([address, nrows]);
    }
    for (let [address, block_size] of direct_blocks) {
      let obj = this._read_direct_block(fh, address, block_size);
      yield obj
    }

    for (let [address, nrows] of indirect_blocks) {
      for (let obj of this._iter_indirect_block(fh, address, nrows)) {
        yield obj;
      }
    }
  }

  _calc_block_size(iblock) {
    let row = Math.floor(iblock / this.header.get("table_width"));
    return 2 ** Math.max(row - 1, 0) * this.header.get('starting_block_size');
  }

  _iblock_nrows_from_block_size(block_size) {
    let log2_block_size = Math.floor(Math.log2(block_size));
    assert(2 ** log2_block_size == block_size);
    return log2_block_size - this._indirect_nrows_sub
  }

  _indirect_info(nrows) {
    let table_width = this.header.get('table_width');
    let nobjects = nrows * table_width;
    let ndirect_max = this._max_direct_nrows * table_width;
    let ndirect, nindirect;
    if (nrows <= ndirect_max) {
      ndirect = nobjects;
      nindirect = 0;
    }
    else {
      ndirect = ndirect_max;
      nindirect = nobjects - ndirect_max;
    }
    return [ndirect, nindirect];
  }

  _int_format(bytelength) {
    return ["B", "H", "I", "Q"][bytelength - 1];
  }

}

var FORMAT_SIGNATURE = struct.unpack_from('8s', new Uint8Array([137, 72, 68, 70, 13, 10, 26, 10]).buffer)[0];
var UNDEFINED_ADDRESS = struct.unpack_from('<Q', new Uint8Array([255, 255, 255, 255, 255, 255, 255, 255]).buffer)[0];

// Version 0 SUPERBLOCK
var SUPERBLOCK_V0 = new Map([
  ['format_signature', '8s'],

  ['superblock_version', 'B'],
  ['free_storage_version', 'B'],
  ['root_group_version', 'B'],
  ['reserved_0', 'B'],

  ['shared_header_version', 'B'],
  ['offset_size', 'B'],            // assume 8
  ['length_size', 'B'],            // assume 8
  ['reserved_1', 'B'],

  ['group_leaf_node_k', 'H'],
  ['group_internal_node_k', 'H'],

  ['file_consistency_flags', 'L'],

  ['base_address_lower', 'Q'],                  // assume 8 byte addressing
  ['free_space_address', 'Q'],            // assume 8 byte addressing
  ['end_of_file_address', 'Q'],
  ['driver_information_address', 'Q']     // assume 8 byte addressing
]);
var SUPERBLOCK_V0_SIZE = _structure_size(SUPERBLOCK_V0);

var SUPERBLOCK_V2_V3 = new Map([
  ['format_signature', '8s'],

  ['superblock_version', 'B'],
  ['offset_size', 'B'],
  ['length_size', 'B'],
  ['file_consistency_flags', 'B'],

  ['base_address', 'Q'],                  // assume 8 byte addressing
  ['superblock_extension_address', 'Q'],  // assume 8 byte addressing
  ['end_of_file_address', 'Q'],           // assume 8 byte addressing
  ['root_group_address', 'Q'],            // assume 8 byte addressing

  ['superblock_checksum', 'I']
]);
var SUPERBLOCK_V2_V3_SIZE = _structure_size(SUPERBLOCK_V2_V3);

var SYMBOL_TABLE_ENTRY = new Map([
  ['link_name_offset', 'Q'],     // 8 byte address
  ['object_header_address', 'Q'],
  ['cache_type', 'I'],
  ['reserved', 'I'],
  ['scratch', '16s'],
]);
var SYMBOL_TABLE_ENTRY_SIZE = _structure_size(SYMBOL_TABLE_ENTRY);

var SYMBOL_TABLE_NODE = new Map([
  ['signature', '4s'],
  ['version', 'B'],
  ['reserved_0', 'B'],
  ['symbols', 'H'],
]);
var SYMBOL_TABLE_NODE_SIZE = _structure_size(SYMBOL_TABLE_NODE);

// III.D Disk Format: Level 1D - Local Heaps
var LOCAL_HEAP = new Map([
  ['signature', '4s'],
  ['version', 'B'],
  ['reserved', '3s'],
  ['data_segment_size', 'Q'],         // 8 byte size of lengths
  ['offset_to_free_list', 'Q'],       // 8 bytes size of lengths
  ['address_of_data_segment', 'Q']   // 8 byte addressing
]);

// III.E Disk Format: Level 1E - Global Heap
var GLOBAL_HEAP_HEADER = new Map([
  ['signature', '4s'],
  ['version', 'B'],
  ['reserved', '3s'],
  ['collection_size', 'Q']
])
var GLOBAL_HEAP_HEADER_SIZE = _structure_size(GLOBAL_HEAP_HEADER);

var GLOBAL_HEAP_OBJECT = new Map([
  ['object_index', 'H'],
  ['reference_count', 'H'],
  ['reserved', 'I'],
  ['object_size', 'Q']   // 8 byte addressing,
])
var GLOBAL_HEAP_OBJECT_SIZE = _structure_size(GLOBAL_HEAP_OBJECT);

//# III.G. Disk Format: Level 1G - Fractal Heap
var FRACTAL_HEAP_HEADER = new Map([
  ['signature', '4s'],
  ['version', 'B'],

  ['object_index_size', 'H'],
  ['filter_info_size', 'H'],
  ['flags', 'B'],

  ['max_managed_object_size', 'I'],
  ['next_huge_object_index', 'Q'],      // 8 byte addressing
  ['btree_address_huge_objects', 'Q'],   // 8 byte addressing

  ['managed_freespace_size', 'Q'],       // 8 byte addressing
  ['freespace_manager_address', 'Q'],    // 8 byte addressing
  ['managed_space_size', 'Q'],           // 8 byte addressing
  ['managed_alloc_size', 'Q'],           // 8 byte addressing
  ['next_directblock_iterator_address', 'Q'], // 8 byte addressing

  ['managed_object_count', 'Q'],         // 8 byte addressing
  ['huge_objects_total_size', 'Q'],      // 8 byte addressing
  ['huge_object_count', 'Q'],            // 8 byte addressing
  ['tiny_objects_total_size', 'Q'],      // 8 byte addressing
  ['tiny_object_count', 'Q'],            // 8 byte addressing

  ['table_width', 'H'],
  ['starting_block_size', 'Q'],          // 8 byte addressing
  ['maximum_direct_block_size', 'Q'],    // 8 byte addressing
  ['log2_maximum_heap_size', 'H'],
  ['indirect_starting_rows_count', 'H'],
  ['root_block_address', 'Q'],           // 8 byte addressing
  ['indirect_current_rows_count', 'H']
])

