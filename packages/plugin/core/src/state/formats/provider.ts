/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { binaryCifHasCategory, binaryCifHasColumn, getBinaryCifHeader } from '@molstar/io/common/binary-cif';
import type { StringLike } from '@molstar/core/util/string-like';
import { Sniffable } from '@molstar/core/util/sniffable';
import type { PluginContext } from '@molstar/plugin/context';
import type { StateObject, StateObjectRef, StateTransformer } from '@molstar/core/state';
import { RuntimeContext, Task } from '@molstar/core/task';
import type { FileNameInfo } from '@molstar/core/util/file-info';
import { PluginStateObject } from '../objects.js';

export interface DataFormatProvider<P = any, R = any, V = any, D = any, Id extends string = string> {
  /** Registration name of the format, used in format params and `DownloadFile` options. */
  readonly name: Id;
  label: string;
  description: string;
  category?: string;
  stringExtensions?: string[];
  binaryExtensions?: string[];
  /**
   * Controls the order in which `DataFormatRegistry.auto` tries providers that share the
   * same extension, independent of registration order. Higher values are tried first.
   * Defaults to 0. Use a positive value for providers with a restrictive `isApplicable`
   * that should take precedence over more generic fallback providers (e.g. a specialized
   * variant of a shared file extension).
   */
  priority?: number;
  /**
   * Sniffs the data: text as a string, binary as bytes; Response Data as itself (a `Sniffable`) with `sniffsResponse`.
   */
  isApplicable?(info: FileNameInfo, data: StringLike | Uint8Array | Sniffable): boolean;
  /**
   * `isApplicable` takes a response as a `Sniffable` (its header, `sniffBytes`, `sniffText`): told without converting
   * it, where a provider without this gets its text, which turns the body to text (`ResponseBody`). The CIF providers,
   * whose parser reads the bytes, set it.
   */
  sniffsResponse?: boolean;
  parse(
    plugin: PluginContext,
    data: StateObjectRef<
      PluginStateObject.Data.Binary | PluginStateObject.Data.String | PluginStateObject.Data.Response
    >,
    params?: P,
  ): Promise<R>;
  /**
   * Parse the data into plain objects without creating state tree nodes. In contrast to `parse`
   * this can be called from within a state update, e.g. from a transformer.
   */
  parseRaw?(plugin: PluginContext, ctx: RuntimeContext, data: StringLike | Uint8Array, params?: P): Promise<D>;
  visuals?(plugin: PluginContext, data: R): Promise<V> | undefined;
  defaultData?: D;
}

/** Identity helper that type checks a provider and keeps its `name` a literal type. */
export function DataFormatProvider<const T extends DataFormatProvider>(provider: T): T {
  return provider;
}

const namedCopies = new WeakMap<object, Map<string, DataFormatProvider>>();

export namespace DataFormatProvider {
  /**
   * Whether a provider takes `o` by its `isApplicable` (one without takes everything). A response is sniffed as itself
   * by a provider that `sniffsResponse`, as its bytes by one that reads `info.ext` as binary, as its text by the others
   * (the body turns to text, `.v6/data-response`).
   */
  export function applicable(
    info: FileNameInfo,
    o: PluginStateObject.Data.String | PluginStateObject.Data.Binary | PluginStateObject.Data.Response,
  ): (p: DataFormatProvider) => boolean {
    const response = o instanceof PluginStateObject.Data.Response ? o.data : undefined;
    return (p) => {
      if (!p.isApplicable) return true;
      if (!response) return p.isApplicable(info, o.data as StringLike | Uint8Array);
      if (p.sniffsResponse) return p.isApplicable(info, response);
      if (p.binaryExtensions?.includes(info.ext)) return p.isApplicable(info, response.bytes);
      return p.isApplicable(info, response.textNow());
    };
  }

  /** A provider that may lack a `name`, such as the user-supplied `customFormats` of the Viewer. */
  export type Unnamed<P = any, R = any, V = any, D = any, Id extends string = string> = Omit<
    DataFormatProvider<P, R, V, D, Id>,
    'name'
  > & { name?: Id };

  /**
   * Returns `provider` when it is already called `name`, otherwise a copy that has `name`. Copies are
   * memoized, so repeated calls with the same arguments return the same object (the copy is a separate
   * identity from `provider`).
   */
  export function withName<P = any, R = any, V = any, D = any>(
    provider: Unnamed<P, R, V, D>,
    name: string,
  ): DataFormatProvider<P, R, V, D> {
    if (provider.name === name) return provider as DataFormatProvider<P, R, V, D>;
    let copies = namedCopies.get(provider);
    if (!copies) {
      copies = new Map();
      namedCopies.set(provider, copies);
    }
    let copy = copies.get(name);
    if (!copy) {
      copy = { ...provider, name };
      copies.set(name, copy);
    }
    return copy as DataFormatProvider<P, R, V, D>;
  }
}

export function rawDataObject(data: StringLike | Uint8Array) {
  return data instanceof Uint8Array
    ? new PluginStateObject.Data.Binary(data as Uint8Array<ArrayBuffer>)
    : new PluginStateObject.Data.String(data);
}

/** Applies a transformer directly, without adding nodes to the state tree. */
export function applyTransformerRaw<A extends StateObject, B extends StateObject, P extends {}>(
  plugin: PluginContext,
  ctx: RuntimeContext,
  transformer: StateTransformer<A, B, P>,
  // not checked against `A`: raw data is String or Binary Data, which a parser's `from` (String or Binary, and
  // Response) may list only one of
  a: StateObject,
  params?: Partial<P>,
): Promise<B> | B {
  const p = transformer.createDefaultParams(a as A, plugin);
  if (params) {
    for (const k of Object.keys(params) as (keyof P)[]) {
      if (params[k] !== undefined) p[k] = params[k]!;
    }
  }
  // no built-in transformer uses the spine, which is only available during a state update
  const result = transformer.definition.apply({ a: a as A, params: p, cache: {}, spine: undefined as any }, plugin);
  return Task.is(result) ? Task.resolveInContext(result, ctx) : result;
}

type CifVariants = 'dscif' | 'segcif' | 'sfcif' | 'coreCif' | -1;
/**
 * The last data sniffed (bytes, held weakly: a big file is not kept) and its variant: every CIF provider's
 * `isApplicable` asks for the same data in turn (`DataFormatRegistry.auto`), and each ask scans it whole when it is not
 * a variant. A string is not cached (it cannot be held weakly), as before.
 */
let lastSniff: { data: WeakRef<object>; ext: string; variant: CifVariants } | undefined;

export function guessCifVariant(info: FileNameInfo, data: Uint8Array | StringLike | Sniffable): CifVariants {
  if (typeof data !== 'object') return sniffCifVariant(info, data);
  const last = lastSniff;
  if (last && last.ext === info.ext && last.data.deref() === data) return last.variant;
  const variant = sniffCifVariant(info, data);
  lastSniff = { data: new WeakRef(data), ext: info.ext, variant };
  return variant;
}

function sniffCifVariant(info: FileNameInfo, data: Uint8Array | StringLike | Sniffable): CifVariants {
  const sniffable = Sniffable.of(data);
  if (info.ext === 'bcif') {
    // read whole: the header walks the MessagePack (a body holding text holds no BinaryCIF)
    const bytes = sniffable.heldBytes;
    if (!bytes) return -1;
    try {
      const header = getBinaryCifHeader(bytes);
      if (header.encoder.startsWith('VolumeServer')) return 'dscif';
      // Assumes volseg-volume-server only serves segments
      if (header.encoder.startsWith('volseg-volume-server')) return 'segcif';

      if (binaryCifHasCategory(header, '_volume_data_3d_info')) {
        if (binaryCifHasCategory(header, '_volume_data_3d')) return 'dscif';
        if (binaryCifHasCategory(header, '_segmentation_data_3d')) return 'segcif';
      }
      if (binaryCifHasCategory(header, '_refln')) {
        if (binaryCifHasColumn(header, '_refln', 'pdbx_FWT') || binaryCifHasColumn(header, '_refln', 'pdbx_DELFWT')) {
          return 'sfcif';
        }
      }
    } catch (e) {
      console.error(e);
    }
  } else if (info.ext === 'cif') {
    // a string, bytes or a response, searched in the form it is held in (`Sniffable.sniffText`)
    const has = (needle: string, where?: 'start') => sniffable.sniffText(needle, where);
    if (has('data_SERVER\n#\n_density_server_result', 'start')) return 'dscif';
    if (has('data_SERVER\n#\ndata_SEGMENTATION_DATA', 'start')) return 'segcif';

    if (has('_volume_data_3d_info.')) {
      if (has('_volume_data_3d.')) return 'dscif';
      if (has('_segmentation_data_3d.')) return 'segcif';
    }

    // Structure factor CIF: has _refln category with map coefficients
    if (has('_refln.pdbx_FWT') || has('_refln.pdbx_DELFWT')) return 'sfcif';

    if (has('atom_site_fract_x') || has('atom_site.fract_x')) return 'coreCif';
  }
  return -1;
}
