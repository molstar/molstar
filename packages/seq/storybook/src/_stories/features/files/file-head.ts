/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { decodeMsgPack } from '@molstar/core/data/msgpack';
import { msgPackColumnDecoder } from '@molstar/core/repr/codecs/msgpack';
import { printSexpr } from '@molstar/core/repr/codecs/sexpr';
import type { Node } from '@molstar/core/repr/tree';
import type { AssetReader } from '@molstar/seqx-core/data';
import { SeqLanguage } from '@molstar/seqx-core/lang';
import { RecordedFiles } from '@molstar/seqx-examples';
import { recordedAssets } from '@molstar/seqx-examples/recorded-pages';

// What the Files stories print beside a view: the head of the file the view reads (the text as it is, or a MessagePack
// file decoded to the JSON it holds, a compact JSON file indented), and the `source` nodes of its spec: the file and
// the transform that says what its columns mean. The story and its headless twin print the same text.

/** How many lines of a text file, and characters of a line, the head holds. */
export const HeadLines = 14;
export const HeadWidth = 120;

const clip = (line: string) => (line.length > HeadWidth ? `${line.slice(0, HeadWidth - 1)}…` : line);

/** The head of a file's bytes: its first lines; a MessagePack file (`format`) as the JSON it decodes to. */
export function fileHead(bytes: Uint8Array, format: string): string {
  let text: string;
  let prefix = '';
  if (format === 'msgpack') {
    try {
      text = JSON.stringify(decodeMsgPack(bytes, msgPackColumnDecoder), null, 2);
      prefix = `MessagePack, ${bytes.byteLength} bytes, as JSON:\n`;
    } catch {
      return `MessagePack, ${bytes.byteLength} bytes (not decoded)`;
    }
  } else {
    text = new TextDecoder('utf-8').decode(bytes);
    // a compact JSON file (the fixtures are stored so) is shown indented, as the JSON it holds
    if (format === 'json' && !text.trimEnd().includes('\n')) {
      try {
        text = JSON.stringify(JSON.parse(text), null, 2);
      } catch {
        // shown as written
      }
    }
  }
  const lines = text.split(/\r\n|\r|\n/);
  const shown = lines.slice(0, HeadLines).map(clip);
  return `${prefix}${shown.join('\n')}${lines.length > HeadLines ? '\n…' : ''}`;
}

/**
 * A file's bytes as its view reads them: from its recording (`<uri>.pages.json`) when the examples declare one
 * (`RecordedFiles`; no other file is asked for a recording), gzip inflated (a
 * genome file's text), so its head is the text the runtime parses.
 */
export async function fileBytes(reader: AssetReader, uri: string, signal: AbortSignal): Promise<Uint8Array> {
  const bytes = await recordedAssets(reader, RecordedFiles).readBytes(uri, signal);
  if (bytes[0] !== 0x1f || bytes[1] !== 0x8b) return bytes;
  const stream = new Blob([bytes as Uint8Array<ArrayBuffer>]).stream().pipeThrough(new DecompressionStream('gzip'));
  return new Uint8Array(await new Response(stream).arrayBuffer());
}

/** The `source` nodes of a spec. */
export function sourceNodes(spec: Node): Node[] {
  const out: Node[] = [];
  const walk = (n: Node) => {
    if (n.kind === 'source') out.push(n);
    for (const c of n.children ?? []) walk(c);
  };
  walk(spec);
  return out;
}

/** The source nodes of a spec as they print. */
export function printSources(spec: Node): string {
  return sourceNodes(spec)
    .map((n) => printSexpr(n, { language: SeqLanguage }))
    .join('\n');
}

/** A file of a view: where it is, and the format it is read as. */
export interface FileOf {
  readonly uri: string;
  readonly format: string;
}

/** The files the `source` nodes of a spec name by location. */
export function filesOf(spec: Node): FileOf[] {
  const out: FileOf[] = [];
  for (const n of sourceNodes(spec)) {
    const uri = n.params?.uri;
    const format = n.params?.format;
    if (typeof uri === 'string' && typeof format === 'string') out.push({ uri, format });
  }
  return out;
}
