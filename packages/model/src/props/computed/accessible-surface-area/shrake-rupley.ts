/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Task, type RuntimeContext } from '@molstar/core/task';
// import { BitFlags } from '../../../mol-util';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { type Structure, type StructureElement, StructureProperties, type Unit } from '@molstar/model/model/structure';
import { assignRadiusForHeavyAtoms } from './shrake-rupley/radii.js';
import { type ShrakeRupleyContext, VdWLookup, MaxAsa, DefaultMaxAsa } from './shrake-rupley/common.js';
import { computeArea } from './shrake-rupley/area.js';
import { SortedArray } from '@molstar/core/data/int';

export const ShrakeRupleyComputationParams = {
  numberOfSpherePoints: PD.Numeric(
    92,
    { min: 12, max: 360, step: 1 },
    {
      description:
        'Number of sphere points to sample per atom: 92 (original paper), 960 (BioJava), 3000 (EPPIC) - see Shrake A, Rupley JA: Environment and exposure to solvent of protein atoms. Lysozyme and insulin. J Mol Biol 1973.',
    },
  ),
  probeSize: PD.Numeric(
    1.4,
    { min: 0.1, max: 4, step: 0.01 },
    { description: 'Corresponds to the size of a water molecule: 1.4 (original paper), 1.5 (occassionally used)' },
  ),
  // buriedRasaThreshold: PD.Numeric(0.16, { min: 0.0, max: 1.0 }, { description: 'below this cutoff of relative accessible surface area a residue will be considered buried - see: Rost B, Sander C: Conservation and prediction of solvent accessibility in protein families. Proteins 1994.' }),
  nonPolymer: PD.Boolean(false, { description: 'Include non-polymer atoms as occluders.' }),
  traceOnly: PD.Boolean(false, {
    description:
      'Compute only using alpha-carbons, if true increase probeSize accordingly (e.g., 4 A). Considers only canonical amino acids.',
  }),
  symmetry: PD.Select<'all' | 'per-group'>(
    'all',
    [
      ['all', 'Every unit'],
      ['per-group', 'Once per symmetry group'],
    ],
    {
      description:
        "Once per symmetry group computes, for each group of symmetry-equivalent units, the unit nearest the structure's centre against the whole structure, and gives its residue areas to the other units of the group: much faster on large assemblies, but approximate. A copy whose surroundings differ from the computed unit's gets that unit's areas: at the edges of a crystal-packing block, in a helical segment, or in an assembly whose chains do not all share one closed group of operators. When every chain is under the same closed point group (a capsid, a symmetric homo-oligomer), every copy has the same surroundings, and the areas differ from Every unit only by the sphere points' sampling noise, which does not rotate with a copy: most residues change, by up to ~15 Å² at 92 points and ~5 Å² at 360.",
    },
  ),
};
export type ShrakeRupleyComputationParams = typeof ShrakeRupleyComputationParams;
export type ShrakeRupleyComputationProps = PD.Values<ShrakeRupleyComputationParams>;

// TODO
// - add back buried and relative asa

export { AccessibleSurfaceArea };

interface AccessibleSurfaceArea {
  readonly serialResidueIndex: ArrayLike<number>;
  readonly area: ArrayLike<number>;
}

namespace AccessibleSurfaceArea {
  /**
   * Adapts the BioJava implementation by Jose Duarte. That implementation is based on the publication by Shrake, A., and
   * J. A. Rupley. "Environment and Exposure to Solvent of Protein Atoms. Lysozyme and Insulin." JMB (1973).
   */
  export function compute(structure: Structure, props: Partial<ShrakeRupleyComputationProps> = {}) {
    const p = { ...PD.getDefaultValues(ShrakeRupleyComputationParams), ...props };
    return Task.create('Compute Accessible Surface Area', async (runtime) => {
      return await calculate(runtime, structure, p);
    });
  }

  async function calculate(
    runtime: RuntimeContext,
    structure: Structure,
    props: ShrakeRupleyComputationProps,
  ): Promise<AccessibleSurfaceArea> {
    const ctx = initialize(structure, props);

    assignRadiusForHeavyAtoms(ctx);
    if (props.symmetry === 'per-group') {
      const representatives = groupRepresentatives(structure);
      await computeArea(runtime, ctx, unitRanges(structure, representatives));
      copyToSymmetryMates(ctx, representatives);
    } else {
      await computeArea(runtime, ctx, [0, structure.elementCount]);
    }

    const { area, serialResidueIndex } = ctx;
    return { area, serialResidueIndex };
  }

  function initialize(structure: Structure, props: ShrakeRupleyComputationProps): ShrakeRupleyContext {
    const { elementCount, atomicResidueCount } = structure;
    const { probeSize, nonPolymer, traceOnly, numberOfSpherePoints } = props;

    return {
      structure,
      probeSize,
      nonPolymer,
      traceOnly,
      spherePoints: generateSpherePoints(numberOfSpherePoints),
      scalingConstant: (4.0 * Math.PI) / numberOfSpherePoints,
      maxLookupRadius: 2 * props.probeSize + 2 * VdWLookup[2], // 2x probe size + 2x largest VdW
      atomRadiusType: new Int8Array(elementCount),
      serialResidueIndex: new Int32Array(elementCount),
      area: new Float32Array(atomicResidueCount),
    };
  }

  /**
   * The unit of each symmetry group whose centre is nearest the structure's: in a crystal-packing block or an open
   * assembly the most likely to have all its neighbours (the first unit is often at a corner). Any unit would do in a
   * closed point-group assembly.
   */
  function groupRepresentatives(structure: Structure) {
    const center = structure.boundary.sphere.center;
    const unitCenter = Vec3();
    const ret: Unit[] = [];
    for (const { units } of structure.unitSymmetryGroups) {
      let best = units[0];
      let bestDistance = Infinity;
      for (const unit of units) {
        Vec3.transformMat4(unitCenter, unit.boundary.sphere.center, unit.conformation.operator.matrix);
        const d = Vec3.squaredDistance(unitCenter, center);
        if (d < bestDistance) {
          best = unit;
          bestDistance = d;
        }
      }
      ret.push(best);
    }
    return ret;
  }

  /** The serial element ranges `[begin, end)` of `units`, flattened */
  function unitRanges(structure: Structure, units: readonly Unit[]) {
    const { unitIndexMap } = structure;
    const { cumulativeUnitElementCount } = structure.serialMapping;
    const ranges: number[] = [];
    for (const unit of units) {
      const begin = cumulativeUnitElementCount[unitIndexMap.get(unit.id)];
      ranges.push(begin, begin + unit.elements.length);
    }
    return ranges;
  }

  /**
   * Gives the other units of each symmetry group the residue areas of its representative. The units of a group have
   * the same elements of the same model (`areUnitsEquivalent`), so the j-th element of each is the same atom and
   * `assignRadiusForHeavyAtoms` gave it the same residue pattern: element for element, residue for residue.
   */
  function copyToSymmetryMates(ctx: ShrakeRupleyContext, representatives: readonly Unit[]) {
    const { structure, serialResidueIndex, area } = ctx;
    const { unitIndexMap, unitSymmetryGroups } = structure;
    const { cumulativeUnitElementCount } = structure.serialMapping;
    for (let g = 0, gl = unitSymmetryGroups.length; g < gl; ++g) {
      const representative = representatives[g];
      const source = cumulativeUnitElementCount[unitIndexMap.get(representative.id)];
      const count = representative.elements.length;
      for (const unit of unitSymmetryGroups[g].units) {
        if (unit === representative) continue;
        const offset = cumulativeUnitElementCount[unitIndexMap.get(unit.id)];
        for (let j = 0; j < count; ++j) {
          const rSI = serialResidueIndex[offset + j];
          if (rSI !== -1) area[rSI] = area[serialResidueIndex[source + j]];
        }
      }
    }
  }

  /** Creates a collection of points on a sphere by the Golden Section Spiral algorithm. */
  function generateSpherePoints(numberOfSpherePoints: number): Vec3[] {
    const points: Vec3[] = [];
    const inc = Math.PI * (3.0 - Math.sqrt(5.0));
    const offset = 2.0 / numberOfSpherePoints;
    for (let k = 0; k < numberOfSpherePoints; ++k) {
      const y = k * offset - 1.0 + offset / 2.0;
      const r = Math.sqrt(1.0 - y * y);
      const phi = k * inc;
      points[points.length] = Vec3.create(Math.cos(phi) * r, y, Math.sin(phi) * r);
    }
    return points;
  }

  export const enum Flags {
    NA = 0x0,
    Buried = 0x1,
    Accessible = 0x2,
  }

  export const Flag = {
    NA: Flags.NA,
    Buried: Flags.Buried,
    Accessible: Flags.Accessible,
  } as const;
  export type Flag = (typeof Flag)[keyof typeof Flag];

  /** Get relative area for a given component id */
  export function normalize(compId: string, asa: number) {
    const maxAsa = MaxAsa[compId] || DefaultMaxAsa;
    return asa / maxAsa;
  }

  /**
   * The value is computed on the root structure: find the element in the root's unit with the location's unit id (an
   * element repeats across the units of an assembly, so the element index alone does not identify it).
   */
  export function getValue(location: StructureElement.Location, accessibleSurfaceArea: AccessibleSurfaceArea) {
    const { area, serialResidueIndex } = accessibleSurfaceArea;
    const { root } = location.structure;
    const unitIndex = root.unitIndexMap.get(location.unit.id);
    if (unitIndex === undefined) return -1;
    const index = SortedArray.indexOf(root.units[unitIndex].elements, location.element);
    if (index < 0) return -1;
    const rSI = serialResidueIndex[root.serialMapping.cumulativeUnitElementCount[unitIndex] + index];
    if (rSI === -1) return -1;
    return area[rSI];
  }

  export function getNormalizedValue(
    location: StructureElement.Location,
    accessibleSurfaceArea: AccessibleSurfaceArea,
  ) {
    const value = getValue(location, accessibleSurfaceArea);
    return value === -1 ? -1 : normalize(StructureProperties.atom.label_comp_id(location), value);
  }

  export function getFlag(location: StructureElement.Location, accessibleSurfaceArea: AccessibleSurfaceArea) {
    const value = getNormalizedValue(location, accessibleSurfaceArea);
    return value === -1 ? Flags.NA : value < 0.16 ? Flags.Buried : Flags.Accessible;
  }
}
