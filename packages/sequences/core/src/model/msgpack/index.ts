// MessagePack (seqspec:format#values): Mol*'s encoder and decoder, vendored here (VENDORED.md), and SVS's typed-array
// extension, ours: an array of numbers as packed little-endian elements, one extension type per element kind, only the
// kinds that hold their numbers exactly. It decodes to a typed array of the kind written.

import { decodeMsgPack, type MsgPackExtDecoder } from "./decode.ts";
import { encodeMsgPack, type MsgPackExtData, type MsgPackExtEncoder } from "./encode.ts";

export { decodeMsgPack, encodeMsgPack, type MsgPackExtData, type MsgPackExtDecoder, type MsgPackExtEncoder };

type Kind =
  | typeof Int8Array
  | typeof Uint8Array
  | typeof Int16Array
  | typeof Uint16Array
  | typeof Int32Array
  | typeof Uint32Array
  | typeof Float32Array
  | typeof Float64Array;
type Typed = InstanceType<Kind>;

/** Extension type codes and their element kinds. */
export const TypedArrayExt: ReadonlyMap<number, Kind> = new Map<number, Kind>([
  [0x11, Int8Array],
  [0x12, Uint8Array],
  [0x13, Int16Array],
  [0x14, Uint16Array],
  [0x15, Int32Array],
  [0x16, Uint32Array],
  [0x17, Float32Array],
  [0x18, Float64Array],
]);
const CODE = new Map<unknown, number>([...TypedArrayExt].map(([code, kind]) => [kind, code]));
const LITTLE = new Uint8Array(Uint16Array.of(1).buffer)[0] === 1;

/** Element bytes reversed in place (a big-endian platform's elements to little-endian, and back). */
function swap(bytes: Uint8Array, size: number) {
  for (let i = 0; size > 1 && i < bytes.length; i += size) bytes.subarray(i, i + size).reverse();
  return bytes;
}

/**
 * Typed arrays as the extension: the lossless kinds as they are (a view, no copy, on a little-endian platform), a
 * `Uint8ClampedArray` as uint8.
 */
export const encodeTypedArrays: MsgPackExtEncoder = (v) => {
  if (!ArrayBuffer.isView(v)) return undefined;
  const a: unknown = v instanceof Uint8ClampedArray ? new Uint8Array(v.buffer, v.byteOffset, v.length) : v;
  const type = CODE.get((a as object).constructor);
  if (type === undefined) return undefined;
  const t = a as Typed;
  const bytes = new Uint8Array(t.buffer, t.byteOffset, t.byteLength);
  return { type, data: LITTLE ? bytes : swap(bytes.slice(), t.BYTES_PER_ELEMENT) };
};

/** The extension as a typed array of its kind (a copy, aligned); other extension types are left to the decoder. */
export const decodeTypedArrays: MsgPackExtDecoder = (type, data) => {
  const Kind = TypedArrayExt.get(type);
  if (!Kind) return undefined;
  if (data.length % Kind.BYTES_PER_ELEMENT !== 0)
    throw new Error(`a typed array (ext ${type}) of ${data.length} bytes`);
  const bytes = new Uint8Array(data);
  if (!LITTLE) swap(bytes, Kind.BYTES_PER_ELEMENT);
  return new Kind(bytes.buffer);
};
