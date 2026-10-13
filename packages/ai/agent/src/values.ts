/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/** Finite JSON values. Untrusted data is decoded at the boundary that receives it. */
export type Json = null | boolean | number | string | JsonObject | readonly Json[];

export interface JsonObject {
  readonly [key: string]: Json;
}

/** An application resource. A retained reference names immutable bytes or carries their digest. */
export interface Reference {
  readonly uri: string;
  readonly mediaType?: string;
  readonly name?: string;
  readonly byteLength?: number;
  /** `sha256:<lowercase hex>` for inline bytes; resolvers may use their own schemes. */
  readonly digest?: string;
}
