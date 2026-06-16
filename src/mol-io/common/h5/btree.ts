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

import {_unpack_struct_from, _structure_size, struct, dtype_getter, bitSize, DataView64} from './core';
import { Filters } from './filters';

class AbstractBTree {
  //B_LINK_NODE = null;
  //NODE_TYPE = null;

  constructor(fh, offset) {
    //""" initalize. """
    this.fh = fh;
    this.offset = offset;
    this.depth = null;
  }

  init() {
    this.all_nodes = new Map();
    this._read_root_node();
    this._read_children();
  }
  
  _read_children() {
    // # Leaf nodes: level 0
    // # Root node: level "depth"
    let node_level = this.depth;
    while (node_level > 0) {
      for (var parent_node of this.all_nodes.get(node_level)) {
        for (var child_addr of parent_node.get('addresses')) {
          this._add_node(this._read_node(child_addr, node_level-1));
        }
      }
      node_level--;
    }
  }
  
  _read_root_node() {
    let root_node = this._read_node(this.offset, null);
    this._add_node(root_node);
    this.depth = root_node.get('node_level');
  }
  
  _add_node(node) {
    let node_level = node.get('node_level');
    if (this.all_nodes.has(node_level)) {
      this.all_nodes.get(node_level).push(node);
    }
    else {
      this.all_nodes.set(node_level, [node]);
    }
  }
  
  _read_node(offset, node_level) {
    // """ Return a single node in the B-Tree located at a given offset. """
    node = this._read_node_header(offset, node_level);
    node.set('keys', []);
    node.set('addresses', []);
    return node
  }

  _read_node_header(offset) {
    //""" Return a single node header in the b-tree located at a give offset. """
    throw "NotImplementedError: must define _read_node_header in implementation class";
  }
}

export class BTreeV1 extends AbstractBTree {
  /*
  """
  HDF5 version 1 B-Tree.
  """
  */
  
  B_LINK_NODE = new Map([
    ['signature', '4s'],

    ['node_type', 'B'],
    ['node_level', 'B'],
    ['entries_used', 'H'],

    ['left_sibling', 'Q'],    // 8 byte addressing
    ['right_sibling', 'Q']    // 8 byte addressing
  ])
  
  _read_node_header(offset, node_level) {
    // """ Return a single node header in the b-tree located at a give offset. """
    let node = _unpack_struct_from(this.B_LINK_NODE, this.fh, offset);
    //assert node['signature'] == b'TREE'
    //assert node['node_type'] == this.NODE_TYPE
    if (node_level != null) {
      if (node.get("node_level") != node_level) {
        throw "node level does not match"
      }
    }
    return node;
  }
  
}


export class BTreeV1Groups extends BTreeV1 {
  /*
  """
  HDF5 version 1 B-Tree storing group nodes (type 0).
  """
  */
  NODE_TYPE = 0;

  constructor(fh, offset) {
    super(fh, offset);
    this.init();
  }

  _read_node(offset, node_level) {
    // """ Return a single node in the B-Tree located at a given offset. """
    let node = this._read_node_header(offset, node_level);
    offset += _structure_size(this.B_LINK_NODE);
    let keys = [];
    let addresses = [];
    let entries_used = node.get('entries_used');
    for (var i=0; i<entries_used; i++) {
      let key = struct.unpack_from('<Q', this.fh, offset)[0];
      offset += 8;
      let address = struct.unpack_from('<Q', this.fh, offset)[0];
      offset += 8;
      keys.push(key);
      addresses.push(address);
    }
    //# N+1 key
    keys.push(struct.unpack_from('<Q', this.fh, offset)[0]);
    node.set('keys', keys);
    node.set('addresses', addresses);
    return node;
  }
    
  symbol_table_addresses() {
    //""" Return a list of all symbol table address. """
    var all_address = [];
    var root_nodes = this.all_nodes.get(0);
    for (var node of root_nodes) {
      all_address = all_address.concat(node.get('addresses'));
    }
    return all_address
  }
}

export class BTreeV1RawDataChunks extends BTreeV1 {
  /*
  HDF5 version 1 B-Tree storing raw data chunk nodes (type 1).
  */
  NODE_TYPE = 1;
  
  constructor(fh, offset, dims) {
    //""" initalize. """
    super(fh, offset);
    this.dims = dims;
    this.init();
  }

  _read_node(offset, node_level) {
    //""" Return a single node in the b-tree located at a give offset. """
    //this.fh.seek(offset)
    let node = this._read_node_header(offset, node_level);
    offset += _structure_size(this.B_LINK_NODE);
    //assert node['signature'] == b'TREE'
    //assert node['node_type'] == 1

    var keys = [];
    var addresses = [];
    let entries_used = node.get('entries_used');
    for (var i=0; i<entries_used; i++) {
      let [chunk_size, filter_mask] = struct.unpack_from('<II', this.fh, offset);
      offset += 8;
      let fmt = '<' + this.dims.toFixed() + 'Q';
      let fmt_size = struct.calcsize(fmt);
      let chunk_offset = struct.unpack_from(fmt, this.fh, offset);
      //console.log(struct.unpack_from('<8B', this.fh, offset));
      offset += fmt_size;
      let chunk_address = struct.unpack_from('<Q', this.fh, offset)[0];
      offset += 8;

      keys.push(new Map([
          ['chunk_size', chunk_size],
          ['filter_mask', filter_mask],
          ['chunk_offset', chunk_offset]
      ]))
      addresses.push(chunk_address);
    }
    node.set('keys', keys);
    node.set('addresses', addresses);
    return node
  }

  construct_data_from_chunks(chunk_shape, data_shape, dtype, filter_pipeline) {
    //""" Build a complete data array from chunks. """
    var true_dtype;
    var item_getter, item_big_endian, item_size;
    if (dtype instanceof Array) {
      true_dtype = dtype;
      let dtype_class = dtype[0];
      if (dtype_class == 'REFERENCE') {
        let size = dtype[1];
        if (size != 8) {
          throw "NotImplementedError('Unsupported Reference type')";
        }
        var dtype = '<u8';
        item_getter = 'getUint64';
        item_big_endian = false;
        item_size = 8;
      }
      else if (dtype_class == 'VLEN_STRING' || dtype_class == 'VLEN_SEQUENCE') {
        item_getter = 'getVLENStruct';
        item_big_endian = false;
        item_size = 16;
      }
      else {
        throw "NotImplementedError('datatype not implemented')";
      }
    }
    else {
      true_dtype = null;
      [item_getter, item_big_endian, item_size] = dtype_getter(dtype);
    }

    //# create array to store data
    var data_size = data_shape.reduce(function(a,b) { return a * b }, 1);
    var chunk_size = chunk_shape.reduce(function(a,b) { return a * b }, 1);
    let dims = (data_shape.length);
    var current_stride = 1;
    var chunk_strides = chunk_shape.slice().map(function(d) {
      let s = current_stride; 
      current_stride *= d; 
      return s
    });
    var current_stride = 1;
    var data_strides = data_shape.slice().reverse().map(function(d) {
      let s = current_stride; 
      current_stride *= d; 
      return s
    }).reverse();
    var data = new Array(data_size);
    let chunk_buffer_size = chunk_size * item_size;
    for (var node of this.all_nodes.get(0)) {
      //console.log(node);
      let node_keys = node.get('keys');
      let node_addresses = node.get('addresses');
      let nkeys = node_keys.length;
      for (var ik=0; ik<nkeys; ik++) {
        let node_key = node_keys[ik];
        let addr = node_addresses[ik];
        var chunk_buffer;
        if (filter_pipeline == null) {
          chunk_buffer = this.fh.slice(addr, addr + chunk_buffer_size);
        }
        else {
          chunk_buffer = this.fh.slice(addr, addr + node_key.get('chunk_size'));
          let filter_mask = node_key.get('filter_mask');
          chunk_buffer = this._filter_chunk(
              chunk_buffer, filter_mask, filter_pipeline, item_size);
        }

        var chunk_offset = node_key.get('chunk_offset').slice(); //(0, -1);
        var apos = chunk_offset.slice();
        var cpos = apos.map(function() {return 0}); // start chunk pos at 0,0,0...
        var cview = new DataView64(chunk_buffer);

        for (var ci=0; ci<chunk_size; ci++) {
          for (var d=dims-1; d>=0; d--) {
            if (cpos[d] >= chunk_shape[d]) {
              cpos[d] = 0;
              apos[d] = chunk_offset[d];
              if (d > 0) {
                cpos[d-1] += 1;
                apos[d-1] += 1;
              }
            }
            else {
              break;
            }
          }
          let inbounds = apos.slice(0,-1).every(function(p, d) { return p < data_shape[d] });
          if (inbounds) {
            let cb_offset = ci * item_size;
            let datum = cview[item_getter](cb_offset, !item_big_endian, item_size);
            let ai = apos.slice(0,-1).reduce(function(prev, curr, index) { 
              return curr * data_strides[index] + prev }, 0);
            data[ai] = datum;
          }
          cpos[dims-1] += 1;
          apos[dims-1] += 1;
        }
      }
    }
    return data;
  }

  _filter_chunk(chunk_buffer, filter_mask, filter_pipeline, itemsize) {
    //""" Apply decompression filters to a chunk of data. """
    let num_filters = filter_pipeline.length;
    let buf = chunk_buffer.slice();
    for (var filter_index=num_filters-1; filter_index >=0; filter_index--) {
      //for i, pipeline_entry in enumerate(filter_pipeline[::-1]):

      //# A filter is skipped is the bit corresponding to its index in the
      //# pipeline is set in filter_mask
      if (filter_mask & (1 << filter_index)) {
        continue
      }
      let pipeline_entry = filter_pipeline[filter_index];
      let filter_id = pipeline_entry.get('filter_id');
      let client_data = pipeline_entry.get('client_data');
      if (Filters.has(filter_id)) {
        buf = Filters.get(filter_id)(buf, itemsize, client_data);
      }
      else {
        throw 'NotImplementedError("Filter with id:' + filter_id.toFixed() + ' not supported")';
      }
    }
    return buf;
  }   
}

export class BTreeV2 extends AbstractBTree {
  /*
  HDF5 version 2 B-Tree.
  */

  // III.A.2. Disk Format: Level 1A2 - Version 2 B-trees
  B_TREE_HEADER = new Map([
    ['signature', '4s'],

    ['version', 'B'],
    ['node_type', 'B'],
    ['node_size', 'I'],
    ['record_size', 'H'],
    ['depth', 'H'],
    ['split_percent', 'B'],
    ['merge_percent', 'B'],

    ['root_address', 'Q'],     // 8 byte addressing
    ['root_nrecords', 'H'],
    ['total_nrecords', 'Q'],   // 8 byte addressing
  ]);

  B_LINK_NODE = new Map([
      ['signature', '4s'],

      ['version', 'B'],
      ['node_type', 'B'],
  ])

  constructor(fh, offset) {
    super(fh, offset);
    this.init();
  }

  _read_root_node() {
    let h = this._read_tree_header(this.offset);
    this.address_formats = this._calculate_address_formats(h);
    this.header = h;
    this.depth = h.get("depth");

    let address = [h.get("root_address"), h.get("root_nrecords"), h.get("total_nrecords")];
    let root_node = this._read_node(address, this.depth);
    this._add_node(root_node);
  }
  
  _read_tree_header(offset) {
    let header = _unpack_struct_from(this.B_TREE_HEADER, this.fh, this.offset);
    //assert header['signature'] == b'BTHD'
    //assert header['node_type'] == this.NODE_TYPE
    return header;
  }
  
  _calculate_address_formats(header) {
    let node_size = header.get("node_size");
    let record_size = header.get("record_size");
    let nrecords_max = 0;
    let ntotalrecords_max = 0;
    let address_formats = new Map();
    let max_depth = header.get("depth");
    for (var node_level=0; node_level <= max_depth; node_level++) {
      let offset_fmt = "";
      let num1_fmt = "";
      let num2_fmt = "";
      let offset_size, num1_size, num2_size;
      if (node_level == 0) { // leaf node
        offset_size = 0;
        num1_size = 0;
        num2_size = 0;
      }
      else if (node_level == 1) { // internal node (twig node)
        offset_size = 8;
        offset_fmt = "<Q";
        num1_size = this._required_bytes(nrecords_max);
        num1_fmt = this._int_format(num1_size);
        num2_size = 0;
      }
      else {  // internal node
        offset_size = 8;
        offset_fmt = "<Q";
        num1_size = this._required_bytes(nrecords_max);
        num1_fmt = this._int_format(num1_size);
        num2_size = this._required_bytes(ntotalrecords_max);
        num2_fmt = this._int_format(num2_size);
      }
      address_formats.set(node_level, [
        offset_size, num1_size, num2_size,
        offset_fmt, num1_fmt, num2_fmt]);
      if (node_level < max_depth) {
        let addr_size = offset_size + num1_size + num2_size;
        nrecords_max = this._nrecords_max(node_size, record_size, addr_size);
        if (ntotalrecords_max > 0) {
          ntotalrecords_max *= nrecords_max;
        }
        else {
          ntotalrecords_max = nrecords_max;
        }
      }
    }

    return address_formats
  }
  
  _nrecords_max(node_size, record_size, addr_size) {
    // """ Calculate the maximal records a node can contain. """
    // node_size = overhead + nrecords_max*record_size + (nrecords_max+1)*addr_size
    //
    // overhead = size(B_LINK_NODE) + 4 (checksum)
    //
    // Leaf node (node_level = 0)
    //   addr_size = 0
    // Internal node (node_level = 1)
    //   addr_size = offset_size + num1_size
    // Internal node (node_level > 1)
    //   addr_size = offset_size + num1_size + num2_size
    return Math.floor((node_size - 10 - addr_size)/(record_size + addr_size));
  }
  
  _required_bytes(integer) {
    // """ Calculate the minimal required bytes to contain an integer. """
    return Math.ceil(bitSize(integer) / 8);
  }
  
  _int_format(bytelength) {
    return ["<B", "<H", "<I", "<Q"][bytelength-1];
  }
  
  _read_node(address, node_level) {
    // """ Return a single node in the B-Tree located at a given offset. """
    let [offset, nrecords, ntotalrecords] = address;
    let node = this._read_node_header(offset, node_level);
    offset += _structure_size(this.B_LINK_NODE);
    let record_size = this.header.get('record_size');
    let keys = [];
    for (let i=0; i<nrecords; i++) {
      let record = this._parse_record(this.fh, offset, record_size);
      offset += record_size;
      keys.push(record);
    }

    let addresses = [];
    let fmts = this.address_formats.get(node_level);
    if (node_level != 0) {
      let [offset_size, num1_size, num2_size, offset_fmt, num1_fmt, num2_fmt] = fmts;
      for (let j=0; j<=nrecords; j++) {
        let address_offset = struct.unpack_from(offset_fmt, this.fh, offset)[0];
        offset += offset_size;
        let num1 = struct.unpack_from(num1_fmt, this.fh, offset)[0];
        offset += num1_size;
        let num2 = num1;
        if (num2_size > 0) {
          num2 = struct.unpack_from(num2_fmt, this.fh, offset)[0];
          offset += num2_size;
        }
        addresses.push([address_offset, num1, num2]);
      }
    }

    node.set('keys', keys);
    node.set('addresses', addresses);
    return node
  }
  
  _read_node_header(offset, node_level) {
    // """ Return a single node header in the b-tree located at a give offset. """
    let node = _unpack_struct_from(this.B_LINK_NODE, this.fh, offset);
    //assert node['node_type'] == this.NODE_TYPE
    if (node_level > 0) {
      // Internal node (has children)
      // assert node['signature'] == b'BTIN'
    }
    else {
      // Leaf node (has no children)
      // assert node['signature'] == b'BTLF'
    }
    node.set("node_level", node_level);
    return node
  }

  * iter_records() {
    // """ Iterate over all records. """
    for (let nodelist of this.all_nodes.values()) {
      for (let node of nodelist) {
        for (let key of node.get('keys')) {
          yield key
        }
      }
    }
  }
     

  _parse_record(record) {
    throw "NotImplementedError"
  }
}

export class BTreeV2GroupNames extends BTreeV2 {
  /*
  HDF5 version 2 B-Tree storing group names (type 5).
  */
  NODE_TYPE = 5

  _parse_record(buf, offset, size) {
    let namehash = struct.unpack_from("<I", buf, offset)[0];
    offset += 4;
    return new Map([['namehash', namehash], ['heapid', buf.slice(offset, offset+7)]]);
  }
}


export class BTreeV2GroupOrders extends BTreeV2 {
  /*
  HDF5 version 2 B-Tree storing group creation orders (type 6).
  */    
  NODE_TYPE = 6

  _parse_record(buf, offset, size) {
    let creationorder = struct.unpack_from("<Q", buf, offset)[0];
    offset += 8;
    return new Map([['creationorder', creationorder], ['heapid', buf.slice(offset, offset+7)]]);
  }
}
