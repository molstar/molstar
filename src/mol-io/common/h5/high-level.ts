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

import {DataObjects} from './dataobjects';
import {SuperBlock} from './misc-low-level';
export { Filters } from './filters';

export class Group {
  /*
    An HDF5 Group which may hold attributes, datasets, or other groups.
    Attributes
    ----------
    attrs : dict
        Attributes for this group.
    name : str
        Full path to this group.
    file : File
        File instance where this group resides.
    parent : Group
        Group instance containing this group.
  */

  /**
   *
   *
   * @memberof Group
   * @member {Group|File} parent;
   * @member {File} file;
   * @member {string} name;
   * @member {DataObjects} _dataobjects;
   * @member {Object} _attrs;
   * @member {Array<string>} _keys;
   */
  // parent;
  // file;
  // name;
  // _links;
  // _dataobjects;
  // _attrs;
  // _keys;

  /**
   * 
   * @param {string} name 
   * @param {DataObjects} dataobjects 
   * @param {Group} [parent] 
   * @param {boolean} [getterProxy=false]
   * @returns {Group}
   */
  constructor(name, dataobjects, parent, getterProxy=false) {
    if (parent == null) {
      this.parent = this;
      this.file = this;
    }
    else {
      this.parent = parent;
      this.file = parent.file;
    }
    this.name = name;

    this._links = dataobjects.get_links();
    this._dataobjects = dataobjects;
    this._attrs = null;  // cached property
    this._keys = null;
    if (getterProxy) {
      return new Proxy(this, groupGetHandler);
    }
  }

  get keys() {
    if (this._keys == null) {
      this._keys = Object.keys(this._links);
    }
    return this._keys.slice();
  }

  get values() {
    return this.keys.map(k => this.get(k));
  }

  length() {
    return this.keys.length;
  }

  _dereference(ref) {
    //""" Deference a Reference object. """
    if (!ref) {
      throw 'cannot deference null reference';
    }
    let obj = this.file._get_object_by_address(ref);
    if (obj == null) {
      throw 'reference not found in file';
    }
    return obj
  }

  get(y) {
    //""" x.__getitem__(y) <==> x[y] """
    if (typeof(y) == 'number') {
      return this._dereference(y);
    }

    var path = normpath(y);
    if (path == '/') {
      return this.file;
    }

    if (path == '.') {
      return this
    }    
    if (/^\//.test(path)) {
      return this.file.get(path.slice(1));
    }

    if (posix_dirname(path) != '') {
      var [next_obj, additional_obj] = path.split(/\/(.*)/);
    }
    else {
      var next_obj = path;
      var additional_obj = '.'
    }
    if (!(next_obj in this._links)) {
      throw next_obj + ' not found in group';
    }

    var obj_name = normpath(this.name + '/' + next_obj);
    let link_target = this._links[next_obj];

    if (typeof(link_target) == "string") {
      try {
        return this.get(link_target)
      } catch (error) {
        return null
      } 
    }

    var dataobjs = new DataObjects(this.file._fh, link_target);
    if (dataobjs.is_dataset) {
      if (additional_obj != '.') {
        throw obj_name + ' is a dataset, not a group';
      }
      return new Dataset(obj_name, dataobjs, this);
    }
    else {
      var new_group = new Group(obj_name, dataobjs, this);
      return new_group.get(additional_obj);
    }
  }

  visit(func) {
    /*
    Recursively visit all names in the group and subgroups.
    func should be a callable with the signature:
        func(name) -> None or return value
    Returning None continues iteration, return anything else stops and
    return that value from the visit method.
    */
    return this.visititems((name, obj) => func(name));
  }

  visititems(func) {
    /*
    Recursively visit all objects in this group and subgroups.
    func should be a callable with the signature:
        func(name, object) -> None or return value
    Returning None continues iteration, return anything else stops and
    return that value from the visit method.
    */
    var root_name_length = this.name.length;
    if (!(/\/$/.test(this.name))) {
      root_name_length += 1;
    }
    //queue = deque(this.values())
    var queue = this.values.slice();
    while (queue) {
      let obj = queue.shift();
      if (queue.length == 1) console.log(obj);
      let name = obj.name.slice(root_name_length);
      let ret = func(name, obj);
      if (ret != null) {
        return ret
      }
      if (obj instanceof Group) {
        queue = queue.concat(obj.values);
      }
    }
    return null
  }

  get attrs() {
    //""" attrs attribute. """
    if (this._attrs == null) {
      this._attrs = this._dataobjects.get_attributes();
    }
    return this._attrs
  }

}

const groupGetHandler = {
  get: function(target, prop, receiver) {
    if (prop in target) {
      return target[prop];
    }
    return target.get(prop);
  }
};


export class File extends Group {
  /*
  Open a HDF5 file.
  Note in addition to having file specific methods the File object also
  inherit the full interface of **Group**.
  File is also a context manager and therefore supports the with statement.
  Files opened by the class will be closed after the with block, file-like
  object are not closed.
  Parameters
  ----------
  filename : str or file-like
      Name of file (string or unicode) or file like object which has read
      and seek methods which behaved like a Python file object.
  Attributes
  ----------
  filename : str
      Name of the file on disk, None if not available.
  mode : str
      String indicating that the file is open readonly ("r").
  userblock_size : int
      Size of the user block in bytes (currently always 0).
  */

  constructor (fh, filename) {
    //""" initalize. """
    //if hasattr(filename, 'read'):
    //    if not hasattr(filename, 'seek'):
    //        raise ValueError(
    //            'File like object must have a seek method')
    
    var superblock = new SuperBlock(fh, 0);
    var offset = superblock.offset_to_dataobjects;
    var dataobjects = new DataObjects(fh, offset);
    super('/', dataobjects, null);
    this.parent = this;

    this._fh = fh
    this.filename = filename || '';

    this.file = this;
    this.mode = 'r';
    this.userblock_size = 0;
  }

  _get_object_by_address(obj_addr) {
    //""" Return the object pointed to by a given address. """
    if (this._dataobjects.offset == obj_addr) {
      return this
    }
    return this.visititems(
      (y) => {(y._dataobjects.offset == obj_addr) ? y : null;}
    );
  }
}

export class Dataset extends Array {
  /*
  A HDF5 Dataset containing an n-dimensional array and meta-data attributes.
  Attributes
  ----------
  shape : tuple
      Dataset dimensions.
  dtype : dtype
      Dataset's type.
  size : int
      Total number of elements in the dataset.
  chunks : tuple or None
      Chunk shape, or NOne is chunked storage not used.
  compression : str or None
      Compression filter used on dataset.  None if compression is not enabled
      for this dataset.
  compression_opts : dict or None
      Options for the compression filter.
  scaleoffset : dict or None
      Setting for the HDF5 scale-offset filter, or None if scale-offset
      compression is not used for this dataset.
  shuffle : bool
      Whether the shuffle filter is applied for this dataset.
  fletcher32 : bool
      Whether the Fletcher32 checksumming is enabled for this dataset.
  fillvalue : float or None
      Value indicating uninitialized portions of the dataset. None is no fill
      values has been defined.
  dim : int
      Number of dimensions.
  dims : None
      Dimension scales.
  attrs : dict
      Attributes for this dataset.
  name : str
      Full path to this dataset.
  file : File
      File instance where this dataset resides.
  parent : Group
      Group instance containing this dataset.
  */

  /**
   *
   *
   * @memberof Dataset
   * @member {Group|File} parent;
   * @member {File} file;
   * @member {string} name;
   * @member {DataObjects} _dataobjects;
   * @member {Object} _attrs;
   * @member {string} _astype;
   */
  // parent;
  // file;
  // name;
  // _dataobjects;
  // _attrs;
  // _astype;

  constructor(name, dataobjects, parent) {
    //""" initalize. """
    super();
    this.parent = parent;
    this.file = parent.file
    this.name = name;

    this._dataobjects = dataobjects
    this._attrs = null;
    this._astype = null;
  }

  get value() {
    var data = this._dataobjects.get_data();
    if (this._astype == null) {
      return data
    }
    return data.astype(this._astype);
  }

  get shape() {
    return this._dataobjects.shape;
  }

  get attrs() {
    return this._dataobjects.get_attributes();
  }

  get dtype() {
    return this._dataobjects.dtype;
  }

  get fillvalue() {
    return this._dataobjects.fillvalue;
  }
}


function posix_dirname(p) {
  let sep = '/';
  let i = p.lastIndexOf(sep) + 1;
  let head = p.slice(0, i);
  let all_sep = new RegExp('^' + sep + '+$');
  let end_sep = new RegExp(sep + '$');
  if (head && !(all_sep.test(head))) {
    head = head.replace(end_sep, '');
  }
  return head
}

function normpath(path) {
  return path.replace(/\/(\/)+/g, '/'); 
  // path = posixpath.normpath(y)
}
