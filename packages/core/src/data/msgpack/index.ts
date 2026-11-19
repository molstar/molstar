/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// @molstar/core/data/msgpack: MessagePack's encoder and decoder with their extension hooks; the tree's codec on them is
// `@molstar/core/repr/codecs/msgpack` (repr/_spec/README.md "Codecs", MessagePack).

export {
  createMsgPackState,
  decodeMsgPack,
  type MsgPackExtDecoder,
  type MsgPackState,
  parseMsgPackValue,
  readMsgPackArrayLength,
  readMsgPackMapLength,
  skipMsgPackValue,
} from './decode.js';
export { encodeMsgPack, type MsgPackExtData, type MsgPackExtEncoder } from './encode.js';
