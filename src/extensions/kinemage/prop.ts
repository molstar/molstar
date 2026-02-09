/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Russ Taylor <russ@reliasolve.com>
 */

/** Based on the ../anvil extension. */

import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Structure, StructureProperties, Unit } from '../../mol-model/structure';
import { CustomPropertyDescriptor } from '../../mol-model/custom-property';
import { KinemageParams, KinemageProps, computeKinemage, isInMembranePlane } from './algorithm';
import { CustomStructureProperty } from '../../mol-model-props/common/custom-structure-property';
import { CustomProperty } from '../../mol-model-props/common/custom-property';
import { Vec3 } from '../../mol-math/linear-algebra';
import { QuerySymbolRuntime } from '../../mol-script/runtime/query/base';
import { CustomPropSymbol } from '../../mol-script/language/symbol';
import { Type } from '../../mol-script/language/type';
import { Task } from '../../mol-task';

import { Kinemage } from '../../mol-io/reader/kin/schema';
import { parseKin } from '../../mol-io/reader/kin/parser';


export const KinemageDataParams = {
    ...KinemageParams
};
export type KinemageDataParams = typeof KinemageDataParams
export type KinemageDataProps = PD.Values<KinemageDataParams>

export { KinemageData };

interface KinemageData {
    /// @todo Remove these
    // point in membrane boundary
    readonly planePoint1: Vec3,
    // point in opposite side of membrane boundary
    readonly planePoint2: Vec3,
    // normal vector of membrane layer
    readonly normalVector: Vec3,
    // the radius of the membrane layer
    readonly radius: number,
    readonly centroid: Vec3,

    /**
     * List of Kinemages read from one or more files.
     */
    readonly kinemages: Kinemage[],

    /**
     * Index of the active KinemageData
     */
    activeKinemage: number
}

const FileSourceParams = {
  input: PD.File({ accept: '.kin', multiple: false })
};
type FileSourceProps = PD.Values<typeof FileSourceParams>

namespace KinemageData {
    export enum Tag {
        Representation = 'kinemage-3d'
    }

    const pos = Vec3();
    export const symbols = {
        isTransmembrane: QuerySymbolRuntime.Dynamic(CustomPropSymbol('computed', 'kinemage.is-kinemage', Type.Bool),
            ctx => {
                const { unit, structure } = ctx.element;
                const { x, y, z } = StructureProperties.atom;
                if (!Unit.isAtomic(unit)) return 0;
                const KinemageData = KinemageDataProvider.get(structure).value;
                if (!KinemageData) return 0;
                Vec3.set(pos, x(ctx.element), y(ctx.element), z(ctx.element));
                const { normalVector, planePoint1, planePoint2 } = KinemageData!;
                return isInMembranePlane(pos, normalVector, planePoint1, planePoint2);
            })
    };

    async function loadKinemageData(data: string): Promise<Kinemage[]> {
      const task = parseKin(data);
      const result = await task.run();
      if (result.isError) {
        throw new Error('Failed to parse KIN data');
      }
      return result.result;
    }

    export async function open(file: FileSourceProps | File): Promise<KinemageData> {

      let fileToRead: File;

      if (file instanceof File) {
        fileToRead = file;
      } else if (file && file.input && file.input.file) {
        fileToRead = file.input.file;
      } else {
        throw new Error('No file given');
      }

      const task = Task.create('Load KIN file', async ctx => {
        const data = await fileToRead.text();
        const kinemages = await loadKinemageData(data);
        return kinemages;
      });

      const kinemages = await task.run();
      const activeKinemage = kinemages.length - 1;
      /// @todo Remove these once we no longer need them
      const planePoint1 = Vec3.create(0, 0, 0);
      const planePoint2 = Vec3.create(0, 0, 0);
      const normalVector = Vec3.create(0, 0, 1);
      const radius = 0;
      const centroid = Vec3.create(0, 0, 0);
      return { kinemages, activeKinemage, planePoint1, planePoint2, normalVector, radius, centroid };
    }
}

export const KinemageDataProvider: CustomStructureProperty.Provider<KinemageDataParams, KinemageData> = CustomStructureProperty.createProvider({
    label: 'Kinemage',
    descriptor: CustomPropertyDescriptor({
        name: 'Kinemage_loaded_data',
        symbols: KinemageData.symbols,
    }),
    type: 'root',
    defaultParams: KinemageDataParams,
    getParams: (data: Structure) => KinemageDataParams,
    isApplicable,
    obtain: async (ctx: CustomProperty.Context, data: Structure, props: Partial<KinemageDataProps>) => {
        const p = { ...PD.getDefaultValues(KinemageDataParams), ...props };
        try {
            return { value: await computeKinemageProps(ctx, data, p) };
        } catch (e) {
            // the "Residues Embedded in Membrane" symbol may bypass isApplicable() checks
            console.warn('Failed to predict membrane orientation. This happens for short peptides and entries without amino acids.');
            return { value: undefined };
        }
    }
});

function isApplicable(structure: Structure) {
    if (!structure.isAtomic) return false;

    for (const model of structure.models) {
        const { byEntityKey } = model.sequence;
        for (const key of Object.keys(byEntityKey)) {
            const { kind, length } = byEntityKey[+key].sequence;
            if (kind !== 'protein') continue; // can only process protein chains
            if (length >= 15) return true; // short peptides might fail
        }
    }
    return false;
}

async function computeKinemageProps(ctx: CustomProperty.Context, data: Structure, props: Partial<KinemageProps>): Promise<KinemageData> {
    const p = { ...PD.getDefaultValues(KinemageParams), ...props };
    return await computeKinemage(data, p).runInContext(ctx.runtime);
}