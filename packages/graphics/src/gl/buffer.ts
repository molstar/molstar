/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { ValueCell } from '@molstar/core/util';
import type { ValueOf } from '@molstar/core/util/type-helpers';

export type UsageHint = 'static' | 'dynamic' | 'stream';
export type DataType = 'uint8' | 'int8' | 'uint16' | 'int16' | 'uint32' | 'int32' | 'float32';

export type DataTypeArrayType = {
  uint8: Uint8Array;
  int8: Int8Array;
  uint16: Uint16Array;
  int16: Int16Array;
  uint32: Uint32Array;
  int32: Int32Array;
  float32: Float32Array;
};
export type ArrayType = ValueOf<DataTypeArrayType>;
export type ArrayKind = keyof DataTypeArrayType;

export type AttributeItemSize = 1 | 2 | 3 | 4 | 16;
export type AttributeKind = 'float32';

export type AttributeDefs = {
  [k: string]: { kind: AttributeKind; itemSize: AttributeItemSize; divisor: number };
};
export type AttributeValues = { [k: string]: ValueCell<ArrayType> };

export type ElementsType = Uint16Array | Uint32Array;
export type ElementsKind = 'uint16' | 'uint32';
