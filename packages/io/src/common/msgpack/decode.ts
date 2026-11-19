/**
 * Copyright (c) 2017-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// MessagePack's decoder, moved to `@molstar/core/data/msgpack`; re-exported at its old path (_spec/README.md).

export {
  createMsgPackState,
  decodeMsgPack,
  type MsgPackExtDecoder,
  type MsgPackState,
  parseMsgPackValue,
  readMsgPackArrayLength,
  readMsgPackMapLength,
  skipMsgPackValue,
} from '@molstar/core/data/msgpack';
