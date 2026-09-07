#!/usr/bin/env node
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Pillot <paul.pillot@tandemai.com>
 */

/**
 * Standalone validation for bond-order perception against the CCD.
 *
 * For each entry in <data-dir>/manifest.json this loads the ligand excerpt (<pdbId>_<compId>.pdb),
 * runs perception (`computeUnitBondOrders` 'force'), and compares the perceived per-edge orders/flags with
 * the CCD-derived expectations (<pdbId>_<compId>_expected.json). The data files are gitignored — run
 * the sibling build-data CLI (build-data.ts) after `npm run build` to populate them first.
 *
 * Usage: node lib/commonjs/cli/bond-order-validation/index.js [--filter COMP] [--data-dir PATH]
 * Exits non-zero if any non-skipped bond is mis-perceived.
 */

import * as argparse from 'argparse';
import * as fs from 'fs';
import * as path from 'path';
import { Tokenizer } from '../../mol-io/reader/common/text/tokenizer';
import { PdbFile } from '../../mol-io/reader/pdb/schema';
import { pdbToMmCif } from '../../mol-model-formats/structure/pdb/to-cif';
import { trajectoryFromMmCIF } from '../../mol-model-formats/structure/mmcif';
import { Task } from '../../mol-task';
import { Structure, Unit } from '../../mol-model/structure';
import { BondType } from '../../mol-model/structure/model/types';
import { computeUnitBondOrders } from '../../mol-model-props/computed/chemistry/bond-orders';

// --- structure helpers (mirrors the jest spec) -------------------------------

function makePdb(pdbText: string): PdbFile {
    const lines = Tokenizer.readAllLines(pdbText);
    return { lines, variant: 'pdb' };
}

async function structureFromPdb(pdbText: string) {
    const cif = await pdbToMmCif(makePdb(pdbText));
    const trajectory = await trajectoryFromMmCIF(cif).run();
    const model = await Task.resolveInContext(trajectory.getFrameAtIndex(0));
    return Structure.ofModel(model);
}

/**
 * Resolve the perceived per-edge order/flags for the structure's first unit. Perception is an opt-in
 * computed property, so the bond graph no longer carries perceived orders — we run
 * `computeUnitBondOrders` and read its per-edge override arrays (parallel to the unit's bond edges).
 */
function perceivedEdges(structure: Structure) {
    const unit = structure.units[0] as Unit.Atomic;
    const { edgeCount, offset, b } = unit.bonds;
    // `force` re-derives every order from geometry (resetting file/table orders first) so the suite
    // exercises our assignment on canonical residues too, not just distance/CONECT-derived bonds.
    const ov = computeUnitBondOrders(structure, unit, 'force');
    const order = ov ? ov.order : unit.bonds.edgeProps.order;
    const flags = ov ? ov.flags : unit.bonds.edgeProps.flags;
    return { unit, edgeCount, offset, b, order, flags };
}

// --- validation model --------------------------------------------------------

interface ExpectedBond {
    a: string;
    /** Single atom name, or an array when any of the listed partners is equally valid
     *  (e.g. interchangeable nitro O, carboxylate O detected from CCD by build-data.mjs). */
    b: string | string[];
    order?: 2 | 3;
    aromatic?: true;
}
/** A bond (atom-name pair, unordered) whose perceived order the suite knowingly gets wrong; `reason`
 *  documents the open issue so it stays tracked. Exempts the pair from BOTH the forward check (an
 *  expected multiple that is missed) and the reverse check (a spurious multiple, e.g. when a missed
 *  double leaves an atom's valence to be matched elsewhere). The structure still runs and every
 *  other bond is checked. */
interface SkipBond { a: string; b: string; reason: string; }
interface ManifestEntry { pdbId: string; compId: string; skipBonds?: SkipBond[]; pass?: boolean; }

interface ComponentResult {
    label: string;
    /** Expected multiples that perception missed. */
    missed: string[];
    /** Perceived multiples the CCD doesn't have. */
    spurious: string[];
    /** skipBonds entries that are now perceived correctly and could be removed. */
    obsoleteSkips: string[];
}

/**
 * Compare the perceived bonds of `structure` with the CCD-derived `expected` set. Mirrors the
 * forward / reverse / skipBonds semantics of the former jest data suite, but accumulates results
 * instead of asserting.
 */
function checkPerceived(structure: Structure, expected: ExpectedBond[], skipBonds: SkipBond[] = []): ComponentResult {
    const { unit, offset, b, order, flags } = perceivedEdges(structure);
    const { label_atom_id, label_comp_id } = unit.model.atomicHierarchy.atoms;
    const nameToLocal = new Map<string, number>();
    const resname = label_comp_id.value(unit.elements[0]);
    for (let i = 0; i < unit.elements.length; i++) {
        nameToLocal.set(label_atom_id.value(unit.elements[i]), i);
    }

    const result: ComponentResult = { label: resname, missed: [], spurious: [], obsoleteSkips: [] };

    // Set of every atom-pair the CCD considers non-single (double / triple / aromatic),
    // as a sorted "min|max" local-index key. Interchangeable partners contribute all
    // their variants. Used by the reverse check to reject spurious perceived multiples.
    const expectedKey = (i: number, j: number) => `${Math.min(i, j)}|${Math.max(i, j)}`;
    const expectedPairs = new Set<string>();

    // Known-unperceivable bonds whose forward check is exempted (see SkipBond / manifest.json).
    const localOf = (name: string) => {
        const i = nameToLocal.get(name);
        if (i === undefined) throw new Error(`skipBonds atom not found: ${name}`);
        return i;
    };
    const skipKeys = new Set(skipBonds.map(s => expectedKey(localOf(s.a), localOf(s.b))));

    for (const exp of expected) {
        const bNames = Array.isArray(exp.b) ? exp.b : [exp.b];
        const bondLabel = `${resname} ${exp.a}-(${bNames.join('|')})`;
        const u = nameToLocal.get(exp.a);
        if (u === undefined) continue; // atom absent from PDB (e.g. OXT) — skip silently
        const vs = bNames.map(name => nameToLocal.get(name)).filter((v): v is number => v !== undefined);
        for (const v of vs) expectedPairs.add(expectedKey(u, v));
        // For interchangeable partners, succeed as soon as any partner has the expected bond.
        let found = false;
        let foundButWrongAromatic = false;
        for (const v of vs) {
            for (let t = offset[u]; t < offset[u + 1]; t++) {
                if (b[t] !== v) continue;
                if (exp.order !== undefined && order[t] !== exp.order) continue;
                foundButWrongAromatic = foundButWrongAromatic || (!!exp.aromatic && !(flags[t] & (BondType.Flag.Aromatic | BondType.Flag.AromaticHuckel)));
                if (exp.aromatic && !(flags[t] & (BondType.Flag.Aromatic | BondType.Flag.AromaticHuckel))) continue;
                found = true;
            }
            if (found) break;
        }
        // A known-unperceivable bond is exempt from the forward check; flag it if it is now
        // perceived correctly so the obsolete skip can be removed from the manifest.
        if (vs.some(v => skipKeys.has(expectedKey(u, v)))) {
            if (found) result.obsoleteSkips.push(`${bondLabel} is now perceived correctly`);
            continue;
        }
        if (!found) result.missed.push(bondLabel + (foundButWrongAromatic ? ' (aromatic flag missing)' : ''));
    }

    // Reverse check: perception must not invent multiple bonds the CCD doesn't have.
    // (A perceived order>1 inside an aromatic ring is allowed against the expected aromatic
    // pair — our Kekule may differ from the CCD's, but it stays within the same ring bonds.)
    // `activeSpuriousSkips` records skip pairs that are *currently* suppressing a real spurious
    // multiple; a spurious-type skip missing from it is obsolete (the spurious is gone).
    const nameOf = (i: number) => label_atom_id.value(unit.elements[i]);
    const activeSpuriousSkips = new Set<string>();
    for (let u = 0; u < unit.elements.length; u++) {
        for (let t = offset[u]; t < offset[u + 1]; t++) {
            const v = b[t];
            if (u >= v) continue;
            if (order[t] <= 1) continue;
            const key = expectedKey(u, v);
            if (skipKeys.has(key)) { // known-wrong pair, exempt both ways
                if (!expectedPairs.has(key)) activeSpuriousSkips.add(key);
                continue;
            }
            if (!expectedPairs.has(key)) {
                result.spurious.push(`${resname} ${nameOf(u)}=${nameOf(v)} (order ${order[t]})`);
            }
        }
    }

    // A spurious-type skip (its pair is not an expected multiple) is obsolete once perception no
    // longer invents that multiple. Missed-type skips are flagged obsolete in the forward loop above.
    for (const s of skipBonds) {
        const key = expectedKey(localOf(s.a), localOf(s.b));
        if (expectedPairs.has(key) || activeSpuriousSkips.has(key)) continue;
        result.obsoleteSkips.push(`${resname} ${s.a}=${s.b} is no longer perceived (spurious gone)`);
    }

    return result;
}

// --- driver ------------------------------------------------------------------

async function run(manifestPath: string, filter?: string) {
    const dataDir = path.dirname(manifestPath);
    if (!fs.existsSync(manifestPath)) {
        console.error(`manifest not found: ${manifestPath}`);
        process.exit(2);
    }
    const manifest: ManifestEntry[] = JSON.parse(fs.readFileSync(manifestPath, 'utf8'));

    let passed = 0;
    let failed = 0;
    let skippedMissingData = 0;
    let obsoleteSkips = 0;
    let skippedKnownBonds = 0;

    for (const { pdbId, compId, skipBonds, pass } of manifest) {
        if (filter && compId !== filter && pdbId !== filter) continue;
        const stem = `${pdbId}_${compId}`;
        const pdbPath = path.join(dataDir, `${stem}.pdb`);
        const expectedPath = path.join(dataDir, `${stem}_expected.json`);
        if (!fs.existsSync(pdbPath) || !fs.existsSync(expectedPath)) {
            console.log(`SKIP  ${pdbId} / ${compId} — data missing (run build-data.mjs)`);
            skippedMissingData++;
            continue;
        }

        const text = fs.readFileSync(pdbPath, 'utf8');
        const expected: ExpectedBond[] = JSON.parse(fs.readFileSync(expectedPath, 'utf8'));
        const structure = await structureFromPdb(text);
        const r = checkPerceived(structure, expected, skipBonds ?? []);

        const ok = r.missed.length === 0 && r.spurious.length === 0;
        if (ok) {
            if (r.obsoleteSkips.length === 0 && skipBonds?.length && !pass) {
                skippedKnownBonds++;
                console.log(`SKIP  ${pdbId} / ${compId} — known problematic bonds`);
            } else {
                passed++;
                console.log(`PASS  ${pdbId} / ${compId}`);
            }
        } else {
            failed++;
            console.log(`FAIL  ${pdbId} / ${compId}`);
            for (const m of r.missed) console.log(`        missed:   ${m}`);
            for (const s of r.spurious) console.log(`        spurious: ${s}`);
        }
        for (const o of r.obsoleteSkips) {
            obsoleteSkips++;
            console.warn(`        WARN obsolete skipBonds: ${o}`);
        }
    }

    console.log('');
    console.log(`${passed} passed, ${failed} failed, ${skippedMissingData} skipped (missing data)` +
        (skippedKnownBonds ? `, ${skippedKnownBonds} skipped (known problematic bonds)` : '') +
        (obsoleteSkips ? `, ${obsoleteSkips} obsolete skipBonds` : ''));

    if (failed > 0) process.exit(1);
}

const parser = new argparse.ArgumentParser({
    add_help: true,
    description: 'Validate bond-order perception against CCD-derived expectations'
});
parser.add_argument('--filter', '-f', {
    help: 'Only validate entries whose compId (or pdbId) matches this value'
});
parser.add_argument('--manifest-path', '-m', {
    help: 'Path to the manifest JSON file; *.pdb / *_expected.json are read from the same directory',
    default: path.resolve(process.cwd(), 'src/cli/bond-order-validation/data/manifest.json')
});
interface Args {
    filter?: string;
    manifest_path: string;
}
const args: Args = parser.parse_args();

run(args.manifest_path, args.filter).catch(e => {
    console.error(e);
    process.exit(1);
});
