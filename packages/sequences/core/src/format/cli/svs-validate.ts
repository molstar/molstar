#!/usr/bin/env node
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `svs-validate <document>... [--strict]`: the static validation (seqspec:format#validation) of document files, JSON or
// MessagePack (seqspec:format#values), each against the kernel's and MVS's registries, one issue per line. Exit 0: no
// errors in any (warnings and notes allowed); 1: errors; 2: usage, or a file that cannot be read (the others are still
// validated). Runs with Node's type stripping.

import { readFileSync } from "node:fs";
import { createRegistries, KernelRegistry, MvsRegistry, readDocument, register, validate } from "../index.js";

const args = process.argv.slice(2);
const files = args.filter((a) => !a.startsWith("--"));
if (files.length === 0 || args.some((a) => a.startsWith("--") && a !== "--strict")) {
  console.error("usage: svs-validate <document.svsj|.svsm>... [--strict]");
  process.exit(2);
}
const reg = createRegistries();
register(reg, [KernelRegistry, MvsRegistry]);
let code = 0;
for (const file of files) {
  let bytes: Uint8Array;
  try {
    bytes = readFileSync(file);
  } catch (e) {
    console.error(`svs-validate: ${e instanceof Error ? e.message : e}`);
    code = 2;
    continue;
  }
  const { doc, issues: envelope } = readDocument(bytes);
  const issues = doc ? validate(doc, reg, { strict: args.includes("--strict") }) : envelope;
  for (const i of issues) console.log(`${file}#${i.path}: ${i.severity} ${i.code}: ${i.message}`);
  const count = (s: string) => issues.filter((i) => i.severity === s).length;
  const errors = count("error");
  console.error(`${file}: ${errors} error(s), ${count("warning")} warning(s), ${count("info")} note(s)`);
  if (errors > 0 && code === 0) code = 1;
}
process.exit(code);
