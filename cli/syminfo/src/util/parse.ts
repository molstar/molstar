/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/**
 * A single `begin_spacegroup...end_spacegroup` block from CCP4 `syminfo.lib`
 * (see `data/sym/syminfo.lib`), reduced to the fields needed to build
 * `src/mol-math/geometry/spacegroup/syminfo.ts` (`RawSpacegroupData`) and
 * `src/mol-math/geometry/spacegroup/_spec/syminfo.lib.ts` (this same shape,
 * emitted verbatim as the test-only regression oracle).
 */
export interface SyminfoEntry {
    /** International Tables (ITA) spacegroup number. */
    readonly number: number;
    /** CCP4-style spacegroup number (`0` if this setting has none). */
    readonly ccp4: number;
    /** Raw `symbol Hall` value (quotes/leading space stripped, any `(...)` suffix kept as-is). */
    readonly hall: string;
    /** Raw `symbol xHM` value (quotes stripped, may carry a ` :H`/`:R`/`:1`/`:2` setting-qualifier suffix). */
    readonly xHM: string;
    /** Raw `basisop` change-of-basis coordinate-triplet expression, e.g. `'z,x,y'`. */
    readonly basisop: string;
    /** Non-empty `symbol old` alias name tokens (there may be zero, one, or several). */
    readonly old: readonly string[];
    /** `symop` coordinate-triplet expressions (point-group coset representatives). */
    readonly symops: readonly string[];
    /** `cenop` coordinate-triplet expressions (lattice-centering translations). */
    readonly cenops: readonly string[];
}

function requiredField(block: string, number: string | number, re: RegExp, fieldName: string): string {
    const m = block.match(re);
    if (!m) throw new Error(`syminfo.lib block (number ${number}) missing '${fieldName}' field`);
    return m[1];
}

function quotedTokens(line: string | undefined): string[] {
    if (!line) return [];
    const tokens: string[] = [];
    const tokenRe = /'([^']*)'/g;
    let tokenMatch: RegExpExecArray | null;
    while ((tokenMatch = tokenRe.exec(line))) {
        const token = tokenMatch[1].trim();
        if (token) tokens.push(token);
    }
    return tokens;
}

/** Parses the raw CCP4 `syminfo.lib` text (`data/sym/syminfo.lib`) into `SyminfoEntry` records. */
export function parseSyminfoLib(raw: string): SyminfoEntry[] {
    const entries: SyminfoEntry[] = [];
    // Anchored to whole lines: the file's header comment also mentions
    // "begin_spacegroup / end_spacegroup records" in prose, which an
    // unanchored match would (wrongly) treat as a tiny nested block.
    const blockRe = /^begin_spacegroup$([\s\S]*?)^end_spacegroup$/gm;
    let blockMatch: RegExpExecArray | null;
    while ((blockMatch = blockRe.exec(raw))) {
        const block = blockMatch[1];

        const number = parseInt(requiredField(block, '?', /^number\s+(\d+)/m, 'number'), 10);
        const basisop = requiredField(block, number, /^basisop\s+(.+)$/m, 'basisop').trim();
        const ccp4 = parseInt(requiredField(block, number, /^symbol ccp4\s+(\d+)/m, 'symbol ccp4'), 10);
        const hall = requiredField(block, number, /^symbol Hall\s+'([^']*)'/m, 'symbol Hall').trim();
        const xHM = quotedTokens(block.match(/^symbol xHM\s+(.*)$/m)?.[1])[0] ?? '';
        const old = quotedTokens(block.match(/^symbol old\s+(.*)$/m)?.[1]);

        const symops = Array.from(block.matchAll(/^symop\s+(.+)$/gm)).map(m => m[1].trim());
        const cenops = Array.from(block.matchAll(/^cenop\s+(.+)$/gm)).map(m => m[1].trim());
        if (!symops.length) throw new Error(`syminfo.lib block (number ${number}, ccp4 ${ccp4}) has no symop lines`);
        if (!cenops.length) throw new Error(`syminfo.lib block (number ${number}, ccp4 ${ccp4}) has no cenop lines`);

        entries.push({ number, ccp4, hall, xHM, basisop, old, symops, cenops });
    }
    return entries;
}

/** Strips a trailing ` :X` CCP4/ITA setting-qualifier suffix, e.g. `'R -3 :H'` -> `'R -3'`. */
export function stripSettingQualifier(xHM: string): string {
    const colon = xHM.indexOf(':');
    return (colon < 0 ? xHM : xHM.slice(0, colon)).trim();
}
