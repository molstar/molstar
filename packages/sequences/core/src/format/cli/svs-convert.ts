#!/usr/bin/env node
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `svs-convert <in> <out> [--pretty] [--plain]`: a document between JSON and MessagePack (seqspec:format#values). The
// input is told by its first byte; the output by its extension (`.svsm`, `.msgpack`: MessagePack, else JSON).
// `--pretty` indents JSON; `--plain` writes MessagePack without the typed-array extension. Exit 0: written; 1: not a
// document; 2: usage, or a file that cannot be read or written. Runs with Node's type stripping.

import { readFileSync, writeFileSync } from "node:fs";
import { readDocument, toJSON, toMsgPack } from "../index.js";

const args = process.argv.slice(2);
const [from, to] = args.filter((a) => !a.startsWith("--"));
if (!from || !to || args.some((a) => a.startsWith("--") && a !== "--pretty" && a !== "--plain")) {
  console.error("usage: svs-convert <in> <out.svsj|.svsm> [--pretty] [--plain]");
  process.exit(2);
}
try {
  const { doc, issues } = readDocument(readFileSync(from));
  for (const i of issues) console.error(`${from}#${i.path}: ${i.severity} ${i.code}: ${i.message}`);
  if (!doc) process.exit(1);
  const msgpack = /\.(svsm|msgpack)$/.test(to);
  writeFileSync(
    to,
    msgpack
      ? toMsgPack(doc, { typedArrays: !args.includes("--plain") })
      : toJSON(doc, { pretty: args.includes("--pretty") }),
  );
} catch (e) {
  console.error(`svs-convert: ${e instanceof Error ? e.message : e}`);
  process.exit(2);
}
