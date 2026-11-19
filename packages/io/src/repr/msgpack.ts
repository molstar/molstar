/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The MessagePack codec of the tree, moved to `@molstar/core/repr/codecs/msgpack`; re-exported at its old path
// (_spec/README.md).

export {
  decodeMsgPackTree,
  encodeMsgPackTree,
  MsgPackCodec,
  MsgPackExt,
  MsgPackIssues,
  swapBytes,
} from '@molstar/core/repr/codecs/msgpack';
