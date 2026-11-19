/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFile, writeFile } from 'node:fs/promises';
import * as fs from 'node:fs';
import { dirname, extname, resolve } from 'node:path';
import { parseArgs } from 'node:util';
import { parseJson } from '@molstar/core/repr/codecs/json';
import { parseSexpr } from '@molstar/core/repr/codecs/sexpr';
import { hasErrors, type Issue, type Node } from '@molstar/core/repr/tree';
import { DefaultAssetManager } from '@molstar/core/util/assets';
import { setFSModule } from '@molstar/core/util/data-source';
import { headlessAssetPolicy, policyAssets } from '@molstar/seqx-core/context';
import { AlignmentRuntime, FastaRuntime, type FormatRuntime, SeqRuntime } from '@molstar/seqx-core/data';
import { AlignmentFileRuntime } from '@molstar/seqx-core/data/alignment-readers';
import { RecordsRuntime } from '@molstar/seqx-core/data/records';
import { renderPng, renderSvg } from '@molstar/seqx-core/render';
import { resvgRasterizer } from './rasterize.js';

// seqx-render: a SeqViewSpec file to an SVG or a PNG, headless (README.md). The spec's assets are read under the
// headless policy (_spec: the spec's folder only; `https:` hosts only when named by `--allow-host`).

/** The exit codes. */
export const ExitCode = {
  Ok: 0,
  /** The spec has errors, or cannot be read as a spec (its issues are printed). */
  Spec: 1,
  /** The command line is wrong. */
  Usage: 2,
  /** The data did not arrive before the timeout: nothing is written. */
  Unsettled: 3,
  /** A file could not be read or written. */
  Io: 4,
} as const;

export const Usage = `Usage: seqx-render <spec.json|spec.seq> -o <out.svg|out.png> [options]

Renders a SeqViewSpec (JSON, S-expression or MessagePack, told by its content) to an SVG or a PNG, told by the output's
extension. Data the spec names is read from files under the spec's folder only; nothing is read from the network unless
--allow-host names a host.

Options:
  -o, --output <file>     the output, .svg or .png (required)
      --width <px>        figure width in CSS px (default 960)
      --height <px|auto>  figure height in CSS px (default auto: the content's)
      --dpr <n>           device pixel ratio of the data's images and of the PNG's pixels (default 1)
      --scale <n>         PNG output pixels per CSS px, times dpr (default 1)
      --theme <light|dark>  the token table (default light)
      --view <id>         the view's address (default: the spec's first)
      --allow-host <host> an https host the spec's data may be read from (repeatable)
      --timeout <ms>      how long to wait for the data (default 10000)
  -h, --help              this text

Exit codes: 0 done; 1 the spec has errors (its issues are printed, one a line); 2 usage; 3 the data did not arrive
in time (nothing written); 4 a file could not be read or written.
`;

export interface CliIo {
  readonly out: (text: string) => void;
  readonly err: (text: string) => void;
}

const stdio: CliIo = {
  out: (t) => process.stdout.write(t),
  err: (t) => process.stderr.write(t),
};

/** One issue a line: `error seq.lang.x /data/p53: message (hint)`. */
export function formatIssue(i: Issue): string {
  const where = i.address ?? (i.path ? `/${i.path.join('/')}` : i.offset !== undefined ? `@${i.offset}` : '');
  return `${i.severity} ${i.code}${where ? ` ${where}` : ''}: ${i.message}${i.hint ? ` (${i.hint})` : ''}`;
}

/** A MessagePack map's first byte (fixmap, map 16, map 32): the only start of a spec in MessagePack. */
const isMapMarker = (b: number | undefined) => b !== undefined && ((b & 0xf0) === 0x80 || b === 0xde || b === 0xdf);

/** A spec's bytes as its tree: MessagePack by its first byte, an S-expression by `(`, else JSON. */
export async function readSpec(bytes: Uint8Array): Promise<{ node?: Node; issues: readonly Issue[] }> {
  if (isMapMarker(bytes[0])) {
    const { MsgPackCodec } = await import('@molstar/core/repr/codecs/msgpack');
    return MsgPackCodec.decode(bytes);
  }
  const text = new TextDecoder('utf-8').decode(bytes);
  return /^﻿?\s*\(/.test(text) ? parseSexpr(text) : parseJson(text);
}

const Runtimes: readonly FormatRuntime[] = [
  SeqRuntime,
  FastaRuntime,
  RecordsRuntime,
  AlignmentRuntime,
  AlignmentFileRuntime,
];

function numberOf(v: string | undefined, name: string, min: number): number | string {
  if (v === undefined) return `--${name} needs a value`;
  const n = Number(v);
  return Number.isFinite(n) && n >= min ? n : `--${name} is not a number of at least ${min}: ${v}`;
}

const Options = {
  output: { type: 'string', short: 'o' },
  width: { type: 'string' },
  height: { type: 'string' },
  dpr: { type: 'string' },
  scale: { type: 'string' },
  theme: { type: 'string' },
  view: { type: 'string' },
  'allow-host': { type: 'string', multiple: true },
  timeout: { type: 'string' },
  help: { type: 'boolean', short: 'h' },
} as const;

const parse = (argv: readonly string[]) => parseArgs({ args: [...argv], options: Options, allowPositionals: true });

/** The command: arguments (without the node and script), its streams; resolves with the exit code, never throws. */
export async function run(argv: readonly string[], io: CliIo = stdio): Promise<number> {
  let parsed: ReturnType<typeof parse>;
  try {
    parsed = parse(argv);
  } catch (e) {
    io.err(`seqx-render: ${(e as Error).message}\n\n${Usage}`);
    return ExitCode.Usage;
  }
  const { values, positionals } = parsed;
  if (values.help) {
    io.out(Usage);
    return ExitCode.Ok;
  }
  const usage = (m: string) => {
    io.err(`seqx-render: ${m}\n\n${Usage}`);
    return ExitCode.Usage;
  };
  if (positionals.length !== 1) return usage('give one spec file');
  if (!values.output) return usage('-o <out.svg|out.png> is required');
  const ext = extname(values.output).toLowerCase();
  if (ext !== '.svg' && ext !== '.png')
    return usage(`the output is .svg or .png, not ${ext || 'without an extension'}`);
  const width = values.width === undefined ? 960 : numberOf(values.width, 'width', 1);
  const height =
    values.height === undefined || values.height === 'auto' ? 'auto' : numberOf(values.height, 'height', 1);
  const dpr = values.dpr === undefined ? 1 : numberOf(values.dpr, 'dpr', 0.1);
  const scale = values.scale === undefined ? 1 : numberOf(values.scale, 'scale', 0.1);
  const timeout = values.timeout === undefined ? 10_000 : numberOf(values.timeout, 'timeout', 1);
  for (const n of [width, height, dpr, scale, timeout]) if (typeof n === 'string' && n !== 'auto') return usage(n);
  if (values.theme !== undefined && values.theme !== 'light' && values.theme !== 'dark')
    return usage(`--theme is light or dark, not ${values.theme}`);
  const theme = (values.theme as 'light' | 'dark' | undefined) ?? 'light';
  if ((dpr as number) * (scale as number) > 8) return usage('dpr × scale is at most 8');

  const specPath = resolve(positionals[0]);
  let bytes: Uint8Array;
  try {
    bytes = await readFile(specPath);
  } catch (e) {
    io.err(`seqx-render: cannot read ${positionals[0]}: ${(e as NodeJS.ErrnoException).code ?? 'error'}\n`);
    return ExitCode.Io;
  }
  const decoded = await readSpec(bytes);
  if (!decoded.node) {
    for (const i of decoded.issues) io.err(`${formatIssue(i)}\n`);
    return ExitCode.Spec;
  }

  // the spec is untrusted input: its files only under its own folder, no network unless a host is named
  setFSModule(fs);
  const assets = policyAssets(
    new DefaultAssetManager(),
    headlessAssetPolicy({ root: dirname(specPath), hosts: values['allow-host'] ?? [] }),
  );
  const common = {
    size: { width: width as number, height: height as 'auto' | number },
    theme,
    runtimes: Runtimes,
    assets,
    timeoutMs: timeout as number,
    dpr: dpr as number,
    ...(values.view ? { view: values.view } : {}),
  };
  let issues: readonly Issue[];
  let settled: boolean;
  let output: string | Uint8Array;
  if (ext === '.svg') {
    const r = await renderSvg(decoded.node, common);
    ({ issues, settled } = r);
    output = r.svg;
  } else {
    const r = await renderPng(decoded.node, {
      ...common,
      // the pixels per CSS px are scale × dpr (renderPng), the data's images are drawn at dpr
      scale: scale as number,
      rasterize: resvgRasterizer,
    });
    ({ issues, settled } = r);
    output = r.png;
  }
  for (const i of issues) io.err(`${formatIssue(i)}\n`);
  if (hasErrors(issues)) return ExitCode.Spec;
  if (!settled) return ExitCode.Unsettled;
  try {
    await writeFile(values.output, output);
  } catch (e) {
    io.err(`seqx-render: cannot write ${values.output}: ${(e as NodeJS.ErrnoException).code ?? 'error'}\n`);
    return ExitCode.Io;
  }
  return ExitCode.Ok;
}
