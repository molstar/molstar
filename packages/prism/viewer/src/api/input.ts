/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { type ApiContext, type ApiFailure, type ApiModuleKey, ApiResult } from './area.js';
import { FileParam, filesOf } from './files.js';
import { OptionalParam } from './params.js';

// One loading vocabulary (api.md §1, "Inputs"): an operation that reads content takes
// exactly one of `url`, `data` (text, `'base64,…'`, or bytes) and `file`, each optional; none or several is `invalid`.
// `inputOf` checks the shape and says which one was given; `contentOf` reads it. What Mol* builds from it (a download,
// raw data, a read file) is the loader's.

/** What `data` takes: text, or bytes (scripts); `'base64,…'` text is bytes too. */
export type DataInput = string | Uint8Array | ArrayBuffer;

/**
 * `data`'s param, typed by its value for scripts (text or bytes, absent: not given); a multiline text param underneath,
 * so the command's form and the agent's tool see text, and `checkValues` lets bytes through (`takesBytes`).
 */
export type DataParamDef = PD.Value<DataInput | undefined>;

export function DataParam(description: string): DataParamDef {
  const text = OptionalParam(PD.Text('', { multiline: true, description }));
  return { ...text, takesBytes: true } as unknown as DataParamDef;
}

/** A text param that takes bytes as well (`DataParam`). */
export function takesBytes(p: PD.Any): boolean {
  return (p as { takesBytes?: unknown }).takesBytes === true;
}

export function isBytes(v: unknown): v is Uint8Array | ArrayBuffer {
  return v instanceof Uint8Array || v instanceof ArrayBuffer;
}

/** The three inputs, described for what the operation reads. */
export function InputParams(o: { url: string; data: string; file: string; accept: string }) {
  return {
    url: OptionalParam(PD.Text('', { description: o.url })),
    data: DataParam(o.data),
    file: OptionalParam(FileParam({ accept: o.accept, description: o.file })),
  };
}

/** The values `inputOf` reads (any operation's that has `InputParams`). */
export interface InputValues {
  readonly url?: string;
  readonly data?: DataInput;
  readonly file?: File | null;
}

/** The one input given; `name` is a file name (the address's last part, the file's, the label), maybe empty. */
export type Input =
  | { readonly kind: 'url'; readonly url: string; readonly name: string }
  | { readonly kind: 'data'; readonly data: string | Uint8Array; readonly name: string }
  | { readonly kind: 'file'; readonly file: File; readonly name: string };

const Base64 = /^base64,/i;

/**
 * The input given, or why not: none (`'<Verb> what? Give url, data or file.'`) or several. `label` names data (its
 * extension may tell the format); `or` is another way the operation takes (`session.open`'s `id`): given alone, the
 * result is `undefined`.
 */
export function inputOf(
  values: InputValues,
  o: { label?: string; verb?: string; or?: { readonly name: string; readonly given: boolean } } = {},
): Input | ApiFailure | undefined {
  const url = typeof values.url === 'string' ? values.url.trim() : '';
  const data = values.data;
  const hasData = typeof data === 'string' ? !!data.trim() : isBytes(data);
  const file = filesOf(values.file)[0];
  const names = ['url', 'data', 'file', ...(o.or ? [o.or.name] : [])];
  const list = `${names.slice(0, -1).join(', ')} or ${names.at(-1)}`;
  const given = [!!url, hasData, !!file, !!o.or?.given].filter(Boolean).length;
  if (given > 1) {
    return ApiResult.fail('invalid', `Give one of ${list}.`, 'Nothing is loaded: say which one is meant.');
  }
  if (!given) return ApiResult.fail('invalid', `${o.verb ?? 'Load'} what? Give ${list}.`);
  if (url) return { kind: 'url', url, name: fileNameOf(url) };
  if (file) return { kind: 'file', file, name: file.name };
  if (!hasData) return undefined;
  let raw: string | Uint8Array;
  if (typeof data === 'string') {
    try {
      raw = Base64.test(data) ? decodeBase64(data) : data;
    } catch {
      return ApiResult.fail('invalid', 'The base64 data does not decode.', "Give 'base64,' and the file's bytes.");
    }
  } else raw = data instanceof Uint8Array ? data : new Uint8Array(data as ArrayBuffer);
  return { kind: 'data', data: raw, name: o.label?.trim() ?? '' };
}

function decodeBase64(text: string): Uint8Array {
  const bin = atob(text.replace(Base64, '').replace(/\s+/g, ''));
  const out = new Uint8Array(bin.length);
  for (let i = 0; i < bin.length; i++) out[i] = bin.charCodeAt(i);
  return out;
}

/** What the host must have for the input (an address: the network; a file: access to files), or why it cannot. */
export function inputModuleProblem(input: Input | { kind: 'id' }, ctx: ApiContext, title: string) {
  const need: ApiModuleKey | undefined =
    input.kind === 'url' || input.kind === 'id' ? 'network' : input.kind === 'file' ? 'files' : undefined;
  if (!need || ctx.modules.has(need)) return undefined;
  return ApiResult.fail(
    'unavailable',
    `${title} needs ${need === 'network' ? 'the network' : 'access to files'}, which this viewer does not have.`,
    need === 'network' ? 'Give its data or a file instead of a url.' : 'Give its data instead of a file.',
  );
}

/** The input's content: an address's and a file's bytes, data as given (text or bytes). */
export async function contentOf(ctx: ApiContext, input: Input): Promise<string | Uint8Array> {
  switch (input.kind) {
    case 'url':
      return ctx.plugin.runTask(ctx.plugin.fetch({ url: input.url, type: 'binary' }));
    case 'file':
      return new Uint8Array(await input.file.arrayBuffer());
    case 'data':
      return input.data;
  }
}

export const asText = (x: string | Uint8Array) => (typeof x === 'string' ? x : new TextDecoder().decode(x));
export const asBytes = (x: string | Uint8Array) => (typeof x === 'string' ? new TextEncoder().encode(x) : x);

/** A zip's signature (`PK`): an `.molx`, `.mvsx` or `.svsx` archive. */
export const isZip = (x: string | Uint8Array) => typeof x !== 'string' && x[0] === 0x50 && x[1] === 0x4b;

/** A file's name from its address (the last part of its path, decoded); empty when it has none. */
export function fileNameOf(url: string): string {
  try {
    return decodeURIComponent(new URL(url, 'http://x/').pathname.split('/').filter(Boolean).pop() ?? '');
  } catch {
    return url.split(/[?#]/)[0].split('/').pop() ?? '';
  }
}

/**
 * A name's extension (a file name, or an address's last part), lowercase, a compression suffix skipped (`1abc.cif.gz`
 * → `cif`); empty when it has none.
 */
export function extensionOf(name: string): string {
  const parts = (name.toLowerCase().split(/[?#]/)[0].split('/').pop() ?? '').split('.');
  while (parts.length > 2 && (parts.at(-1) === 'gz' || parts.at(-1) === 'zip')) parts.pop();
  return parts.length > 1 ? parts.at(-1)! : '';
}
