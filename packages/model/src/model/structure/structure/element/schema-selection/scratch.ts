/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/** The `ScratchMemory` context (`@molstar/core/util/int-buffer`) of the buffers which the resolution of rows reuses between calls:
 * `ScratchMemory.release(SchemaSelectionScratch)` frees them. */
export const SchemaSelectionScratch = 'schema-selection';
