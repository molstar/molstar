/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { getFileNameInfo } from '@molstar/core/util/file-info';
import { utf8Read } from '@molstar/core/util/utf8';
import { DataFormatProvider } from '@molstar/plugin/state/formats/provider';
import type { PluginStateObject } from '@molstar/plugin/state/objects';
import { type ApiContext, type ApiFailure, ApiResult } from './area.js';
import { extensionOf } from './input.js';

// The loaders' formats (api.md §1, "Inputs"): `format` is `'auto'` by default, read from the name's extension by
// Mol*'s format registry (`plugin.dataFormats`: each format's string and binary extensions, `autoOrder` and
// `isApplicable`, as `DataFormatRegistry.auto` reads them), else from the data (data without a name: a few signatures);
// whether an address or a file is binary follows from the extension or the format unless given. A format that cannot
// be told fails `invalid`, naming the formats the operation takes.

/** The formats an operation takes: some of the registry's categories (`VolumeFormatCategory`, …), or all. */
export interface FormatScope {
  readonly categories?: readonly string[];
  /** What they are, for the messages ('a volume'). */
  readonly what?: string;
}

/** A format settled before the data is read. */
export interface FormatChoice {
  /** The format named or told by the data; undefined: by the extension once the data is read (`providerFor`). */
  readonly provider: DataFormatProvider | undefined;
  /** The formats the extension may be, best first (the registry's `autoOrder`); empty when one was named. */
  readonly candidates: readonly DataFormatProvider[];
  /** Read the address or the file as binary. */
  readonly binary: boolean;
}

const inScope = (scope: FormatScope, p: DataFormatProvider) =>
  !scope.categories || scope.categories.includes(p.category ?? '');

/** The formats an operation takes, by name, for the messages. */
export function formatNames(ctx: ApiContext, scope: FormatScope = {}): string[] {
  return ctx.plugin.dataFormats.list.filter((e) => inScope(scope, e.provider)).map((e) => e.name);
}

const formatsHint = (ctx: ApiContext, scope: FormatScope) =>
  `Give its format: auto, ${formatNames(ctx, scope).join(', ')}.`;

/**
 * The format of what `name` names (a file name, an address's last part, a label): `format` by name (`bcif` is mmCIF read
 * as binary), or `'auto'` (empty too) from the extension, else from `content` when the data is at hand; `binary` given,
 * else the extension's kind, else the format's (one that reads only binary is binary).
 */
export async function chooseFormat(
  ctx: ApiContext,
  o: {
    readonly format: string | undefined;
    readonly name: string;
    readonly binary: boolean | undefined;
    readonly content?: string | Uint8Array;
    readonly scope?: FormatScope;
  },
): Promise<FormatChoice | ApiFailure> {
  const scope = o.scope ?? {};
  const registry = ctx.plugin.dataFormats;
  let format = (o.format ?? '').trim().toLowerCase() || 'auto';
  let binary = o.binary;
  if (format === 'bcif') {
    format = 'mmcif';
    binary ??= true;
  }
  const ext = extensionOf(o.name);
  if (format === 'auto') {
    // A lazy extension that reads the extension loads first (viewer-spec.md §3).
    if (ext) {
      const content = o.content;
      await ctx.extensions.ensureFor({
        fileExtensions: [ext],
        files: content !== undefined && o.name ? [{ name: o.name, content: async () => content }] : undefined,
      });
    }
    const candidates = registry.autoOrder
      .map((e) => e.provider)
      .filter((p) => inScope(scope, p) && (p.binaryExtensions?.includes(ext) || p.stringExtensions?.includes(ext)));
    if (ext && candidates.length) {
      return {
        provider: undefined,
        candidates,
        binary: binary ?? candidates.some((p) => p.binaryExtensions?.includes(ext)),
      };
    }
    const told = o.content !== undefined ? formatOfContent(o.content, scope) : undefined;
    if (told && registry.has(told)) {
      return { provider: registry.get(told), candidates: [], binary: binary ?? isBinaryContent(o.content!, told) };
    }
    const what = o.name ? `${o.name}` : 'the data';
    return ApiResult.fail(
      'invalid',
      ext ? `'.${ext}' is not ${scope.what ?? 'a'} format the viewer reads (${what}).` : `Which format is ${what}?`,
      formatsHint(ctx, scope),
    );
  }
  // A lazy extension that reads the format loads first.
  await ctx.extensions.ensureFor({ formats: [format] });
  if (!registry.has(format)) {
    return ApiResult.fail('invalid', `Unknown data format '${format}'.`, formatsHint(ctx, scope));
  }
  const provider = registry.get(format)!;
  if (!inScope(scope, provider)) {
    return ApiResult.fail('invalid', `${provider.label} is not ${scope.what ?? 'a'} format.`, formatsHint(ctx, scope));
  }
  const kind = provider.binaryExtensions?.includes(ext)
    ? true
    : provider.stringExtensions?.includes(ext)
      ? false
      : // A density server's cell is BinaryCIF unless asked otherwise (v5's default for `dscif`).
        provider.name === 'dscif' || (!!provider.binaryExtensions?.length && !provider.stringExtensions?.length);
  return { provider, candidates: [], binary: binary ?? kind };
}

/** The format once the data is read: the one chosen, else the registry's pick by the extension and the data. */
export function providerFor(
  ctx: ApiContext,
  choice: FormatChoice,
  name: string,
  data: PluginStateObject.Data.Binary | PluginStateObject.Data.String | PluginStateObject.Data.Response,
): DataFormatProvider | undefined {
  if (choice.provider) return choice.provider;
  // `DataFormatRegistry.auto` over the formats in scope (the candidates, in its order).
  return choice.candidates.find(DataFormatProvider.applicable(getFileNameInfo(name), data));
}

/**
 * The format a choice names without the data: the one chosen, else the extension's general one (the first without a
 * `priority`: those that take precedence are specialized variants a data check picks, `DataFormatProvider.priority`).
 */
export function formatNameOf(choice: FormatChoice): string {
  return (choice.provider ?? choice.candidates.find((p) => !p.priority) ?? choice.candidates[0]).name;
}

const Structures = ['Trajectory'];
const Volumes = ['Volume'];

/**
 * A format from the data itself, for data without a name: mmCIF text (`data_`), PDB records, BinaryCIF (MessagePack
 * naming its data blocks), a CCP4/MRC map (`MAP ` at byte 208). The registry reads extensions only.
 */
/** Whether content told as `format` is binary: bytes of a binary format, not text as its bytes. */
function isBinaryContent(content: string | Uint8Array, format: string): boolean {
  if (typeof content === 'string') return false;
  return format === 'ccp4' || (content[0] & 0xf0) === 0x80;
}

function formatOfContent(content: string | Uint8Array, scope: FormatScope): string | undefined {
  const takes = (categories: string[]) => !scope.categories || categories.some((c) => scope.categories!.includes(c));
  if (typeof content === 'string') {
    if (!takes(Structures)) return undefined;
    const head = content.slice(0, 4096);
    if (/^\s*(#[^\n]*\n\s*)*data_/.test(head)) return 'mmcif';
    if (/^(HEADER|ATOM {2}|HETATM|CRYST1|MODEL |COMPND|TITLE )/m.test(head)) return 'pdb';
    return undefined;
  }
  if (takes(Volumes) && content.length > 212 && String.fromCharCode(...content.subarray(208, 212)) === 'MAP ') {
    return 'ccp4';
  }
  if (takes(Structures) && (content[0] & 0xf0) === 0x80) {
    const head = String.fromCharCode(...content.subarray(0, 256));
    if (head.includes('dataBlocks')) return 'mmcif';
  }
  // text as its bytes (a response holds the bytes): told by its head, decoded
  return formatOfContent(utf8Read(content, 0, Math.min(4096, content.length)), scope);
}
