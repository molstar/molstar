/**
 * @file Bit array
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Paul Pillot <paulpillot@gmail.com>
 * @private
 */

/**
 * Compute the Hamming weight of a 32-bit unsigned integer
 * @param  {Integer} v - a 32-bit unsigned integer
 * @return {Integer} the Hamming weight
 */
export function hammingWeight(v: number) {
  // works with signed or unsigned shifts
  v -= (v >>> 1) & 0x55555555;
  v = (v & 0x33333333) + ((v >>> 2) & 0x33333333);
  return (((v + (v >>> 4)) & 0xf0f0f0f) * 0x1010101) >>> 24;
}

/**
 * Bit array
 *
 * Based heavily on https://github.com/lemire/FastBitSet.js
 * which is licensed under the Apache License, Version 2.0.
 */
export default class BitArray {
  private _words: Uint32Array;
  public length: number;

  /**
   * @param  {Integer} length - array length
   * @param  {Boolean} [setAll] - initialize with true
   */
  constructor(length: number, setAll?: boolean) {
    this.length = length;
    this._words = new Uint32Array((length + 32) >>> 5);
    if (setAll === true) {
      this.setAll();
    }
  }

  /**
   * Get value at index
   * @param  {Integer} index - the index
   * @return {Boolean} value
   */
  get(index: number) {
    return (this._words[index >>> 5] & (1 << index)) !== 0;
  }

  /**
   * Set value at index to true
   * @param  {Integer} index - the index
   * @return {undefined}
   */
  set(index: number) {
    this._words[index >>> 5] |= 1 << index;
  }

  /**
   * Set value at index to false
   * @param  {Integer} index - the index
   * @return {undefined}
   */
  clear(index: number) {
    this._words[index >>> 5] &= ~(1 << index);
  }

  /**
   * Flip value at index
   * @param  {Integer} index - the index
   * @return {undefined}
   */
  flip(index: number) {
    this._words[index >>> 5] ^= 1 << index;
  }

  _assignRange(start: number, end: number, value: boolean) {
    if (end < start) return;
    const words = this._words;
    const wordValue = value === true ? 0xffffffff : 0;
    const wordStart = start >>> 5;
    const wordEnd = end >>> 5;
    // set complete words when applicable
    for (let k = wordStart; k < wordEnd; ++k) {
      words[k] = wordValue;
    }
    // set parts of the range not spanning complete words
    const startWord = wordStart << 5;
    const endWord = wordEnd << 5;
    if (value === true) {
      if (end - start < 32) {
        for (let i = start, n = end + 1; i < n; ++i) {
          words[i >>> 5] |= 1 << i;
        }
      } else {
        for (let i = start, n = startWord; i < n; ++i) {
          words[i >>> 5] |= 1 << i;
        }
        for (let i = endWord, n = end + 1; i < n; ++i) {
          words[i >>> 5] |= 1 << i;
        }
      }
    } else {
      if (end - start < 32) {
        for (let i = start, n = end + 1; i < n; ++i) {
          words[i >>> 5] &= ~(1 << i);
        }
      } else {
        for (let i = start, n = startWord; i < n; ++i) {
          words[i >>> 5] &= ~(1 << i);
        }
        for (let i = endWord, n = end + 1; i < n; ++i) {
          words[i >>> 5] &= ~(1 << i);
        }
      }
    }
    return this;
  }

  /**
   * Set bits of the given range
   * @param {Integer} start - start index
   * @param {Integer} end - end index
   * @return {BitArray} this object
   */
  setRange(start: number, end: number) {
    return this._assignRange(start, end, true);
  }

  /**
   * Clear bits of the given range
   * @param {Integer} start - start index
   * @param {Integer} end - end index
   * @return {BitArray} this object
   */
  clearRange(start: number, end: number) {
    return this._assignRange(start, end, false);
  }

  /**
   * Set bits at all given indices
   * @param {...Integer} arguments - indices
   * @return {Boolean} this object
   */
  setBits(...indices: number[]) {
    const words = this._words;
    const n = indices.length;
    for (let i = 0; i < n; ++i) {
      const index = indices[i];
      words[index >>> 5] |= 1 << index;
    }
    return this;
  }

  /**
   * Clear bits at all given indices
   * @param {...Integer} arguments - indices
   * @return {Boolean} this object
   */
  clearBits(...indices: number[]) {
    const words = this._words;
    const n = indices.length;
    for (let i = 0; i < n; ++i) {
      const index = indices[i];
      words[index >>> 5] &= ~(1 << index);
    }
    return this;
  }

  /**
   * Set all bits of the array
   * @return {BitArray} this object
   */
  setAll() {
    return this._assignRange(0, this.length - 1, true);
  }

  /**
   * Clear all bits of the array
   * @return {BitArray} this object
   */
  clearAll() {
    return this._assignRange(0, this.length - 1, false);
  }

  /**
   * Flip all the values in the array
   * @return {BitArray} this object
   */
  flipAll() {
    const count = this._words.length;
    const words = this._words;
    const bs = 32 - (this.length % 32);
    for (let k = 0; k < count - 1; ++k) {
      words[k] = ~words[k];
    }
    words[count - 1] = ~(words[count - 1] << bs) >>> bs;
    return this;
  }

  _isRangeValue(start: number, end: number, value: boolean) {
    if (end < start) return;
    const words = this._words;
    const wordValue = value === true ? 0xffffffff : 0;
    const wordStart = start >>> 5;
    const wordEnd = end >>> 5;
    // set complete words when applicable
    for (let k = wordStart; k < wordEnd; ++k) {
      if (words[k] !== wordValue) return false;
    }
    // set parts of the range not spanning complete words
    if (end - start < 32) {
      for (let i = start, n = end + 1; i < n; ++i) {
        if (!!(words[i >>> 5] & (1 << i)) !== value) return false;
      }
    } else {
      const startWord = wordStart << 5;
      const endWord = wordEnd << 5;
      for (let i = start, n = startWord << 5; i < n; ++i) {
        if (!!(words[i >>> 5] & (1 << i)) !== value) return false;
      }
      for (let i = endWord, n = end + 1; i < n; ++i) {
        if (!!(words[i >>> 5] & (1 << i)) !== value) return false;
      }
    }
    return true;
  }

  /**
   * Test if bits in given range are set
   * @param {Integer} start - start index
   * @param {Integer} end - end index
   * @return {BitArray} this object
   */
  isRangeSet(start: number, end: number) {
    return this._isRangeValue(start, end, true);
  }

  /**
   * Test if bits in given range are clear
   * @param {Integer} start - start index
   * @param {Integer} end - end index
   * @return {BitArray} this object
   */
  isRangeClear(start: number, end: number) {
    return this._isRangeValue(start, end, false);
  }

  /**
   * Test if all bits in the array are set
   * @return {Boolean} test result
   */
  isAllSet() {
    return this._isRangeValue(0, this.length - 1, true);
  }

  /**
   * Test if all bits in the array are clear
   * @return {Boolean} test result
   */
  isAllClear() {
    return this._isRangeValue(0, this.length - 1, false);
  }

  /**
   * Test if bits at all given indices are set
   * @param {...Integer} arguments - indices
   * @return {Boolean} test result
   */
  isSet(...indices: number[]) {
    const words = this._words;
    const n = indices.length;
    for (let i = 0; i < n; ++i) {
      const index = indices[i];
      if ((words[index >>> 5] & (1 << index)) === 0) return false;
    }
    return true;
  }

  /**
   * Test if bits at all given indices are clear
   * @param {...Integer} arguments - indices
   * @return {Boolean} test result
   */
  isClear(...indices: number[]) {
    const words = this._words;
    const n = indices.length;
    for (let i = 0; i < n; ++i) {
      const index = indices[i];
      if ((words[index >>> 5] & (1 << index)) !== 0) return false;
    }
    return true;
  }

  /**
   * Test if two BitArrays are identical in all their values
   * @param {BitArray} otherBitarray - the other BitArray
   * @return {Boolean} test result
   */
  isEqualTo(otherBitarray: BitArray) {
    const words1 = this._words;
    const words2 = otherBitarray._words;
    const count = Math.min(words1.length, words2.length);
    for (let k = 0; k < count; ++k) {
      if (words1[k] !== words2[k]) {
        return false;
      }
    }
    return true;
  }

  /**
   * How many set bits?
   * @return {Integer} number of set bits
   */
  getSize() {
    const count = this._words.length;
    const words = this._words;
    let size = 0;
    for (let i = 0; i < count; ++i) {
      size += hammingWeight(words[i]);
    }
    return size;
  }

  /**
   * Calculate difference betwen this and another bit array.
   * Store result in this object.
   * @param  {BitArray} otherBitarray - the other bit array
   * @return {BitArray} this object
   */
  difference(otherBitarray: BitArray) {
    const words1 = this._words;
    const words2 = otherBitarray._words;
    const count = Math.min(words1.length, words2.length);
    for (let k = 0; k < count; ++k) {
      words1[k] = words1[k] & ~words2[k];
    }
    for (let k = words1.length; k < count; ++k) {
      words1[k] = 0;
    }
    return this;
  }

  /**
   * Calculate union betwen this and another bit array.
   * Store result in this object.
   * @param  {BitArray} otherBitarray - the other bit array
   * @return {BitArray} this object
   */
  union(otherBitarray: BitArray) {
    const words1 = this._words;
    const words2 = otherBitarray._words;
    const count = Math.min(words1.length, words2.length);
    for (let k = 0; k < count; ++k) {
      words1[k] |= words2[k];
    }
    for (let k = words1.length; k < count; ++k) {
      words1[k] = 0;
    }
    return this;
  }

  /**
   * Calculate intersection betwen this and another bit array.
   * Store result in this object.
   * @param  {BitArray} otherBitarray - the other bit array
   * @return {BitArray} this object
   */
  intersection(otherBitarray: BitArray) {
    const words1 = this._words;
    const words2 = otherBitarray._words;
    const count = Math.min(words1.length, words2.length);
    for (let k = 0; k < count; ++k) {
      words1[k] &= words2[k];
    }
    for (let k = words1.length; k < count; ++k) {
      words1[k] = 0;
    }
    return this;
  }

  /**
   * Test if there is any intersection betwen this and another bit array.
   * @param  {BitArray} otherBitarray - the other bit array
   * @return {Boolean} test result
   */
  intersects(otherBitarray: BitArray) {
    const words1 = this._words;
    const words2 = otherBitarray._words;
    const count = Math.min(words1.length, words2.length);
    for (let k = 0; k < count; ++k) {
      if ((words1[k] & words2[k]) !== 0) {
        return true;
      }
    }
    return false;
  }

  /**
   * Calculate the number of bits in common betwen this and another bit array.
   * @param  {BitArray} otherBitarray - the other bit array
   * @return {Integer} size
   */
  getIntersectionSize(otherBitarray: BitArray) {
    const words1 = this._words;
    const words2 = otherBitarray._words;
    const count = Math.min(words1.length, words2.length);
    let size = 0;
    for (let k = 0; k < count; ++k) {
      size += hammingWeight(words1[k] & words2[k]);
    }
    return size;
  }

  /**
   * Calculate intersection betwen this and another bit array.
   * Store result in a new bit array.
   * @param  {BitArray} otherBitarray - the other bit array
   * @return {BitArray} the new bit array
   */
  makeIntersection(otherBitarray: BitArray) {
    const words1 = this._words;
    const words2 = otherBitarray._words;
    const count = Math.min(words1.length, words2.length);
    const wordsA = new Uint32Array(count);
    const intersection = Object.create(BitArray.prototype);
    intersection._words = wordsA;
    intersection.length = Math.min(this.length, otherBitarray.length);
    for (let k = 0; k < count; ++k) {
      wordsA[k] = words1[k] & words2[k];
    }
    return intersection;
  }

  /**
   * Iterate over all set bits in the array
   * @param  {function( index: Integer, i: Integer )} callback - the callback
   * @return {undefined}
   */
  forEach(callback: (index: number, i: number) => any) {
    const count = this._words.length;
    const words = this._words;
    let i = 0;
    for (let k = 0; k < count; ++k) {
      let w = words[k];
      while (w !== 0) {
        const t = w & -w;
        const index = (k << 5) + hammingWeight(t - 1);
        callback(index, i);
        w ^= t;
        ++i;
      }
    }
  }

  /**
   * Get an array with the set bits
   * @return {Array} bit indices
   */
  toArray() {
    const words = this._words;
    const answer = new Array(this.getSize());
    const count = this._words.length;
    let pos = 0;
    for (let k = 0; k < count; ++k) {
      let w = words[k];
      while (w !== 0) {
        const t = w & -w;
        answer[pos++] = (k << 5) + hammingWeight(t - 1);
        w ^= t;
      }
    }
    return answer;
  }

  toString() {
    return "{" + this.toArray().join(",") + "}";
  }

  toSeleString() {
    const sele = this.toArray().join(",");
    return sele ? "@" + sele : "NONE";
  }

  /**
   * Clone this object
   * @return {BitArray} the cloned object
   */
  clone() {
    const clone = Object.create(BitArray.prototype);
    clone.length = this.length;
    clone._words = new Uint32Array(this._words);
    return clone;
  }
}
