/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Russ Taylor <russ@reliasolve.com>
 */

/** Based on the ../anvil extension. */

import { Structure } from '../../mol-model/structure';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { CustomPropertyDescriptor } from '../../mol-model/custom-property';
import { CustomStructureProperty } from '../../mol-model-props/common/custom-structure-property';
import { CustomProperty } from '../../mol-model-props/common/custom-property';
import { Task } from '../../mol-task';

import { Kinemage } from './reader/schema';
import { parseKin } from './reader/parser';

export const KinemageParams = {
};
export type KinemageParams = typeof KinemageParams
export type KinemageProps = PD.Values<KinemageParams>

export const KinemageDataParams = {
  ...KinemageParams
};
export type KinemageDataParams = typeof KinemageDataParams
export type KinemageDataProps = PD.Values<KinemageDataParams>

export { KinemageData };

interface KinemageData {
    /**
     * List of Kinemages read from one or more files.
     */
    readonly kinemages: Kinemage[]
}

const FileSourceParams = {
  input: PD.File({ accept: '.kin', multiple: false })
};
type FileSourceProps = PD.Values<typeof FileSourceParams>

namespace KinemageData {
    export enum Tag {
        Representation = 'kinemage-3d'
    }

    export const symbols = {
    }

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
      return { kinemages };
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
    return false;
}

async function computeKinemageProps(ctx: CustomProperty.Context, data: Structure, props: Partial<KinemageProps>): Promise<KinemageData> {
  // Return an empty KinemageData object since the actual data will be loaded asynchronously via the `open` method.
  // This allows the property to be attached to the structure without blocking on file loading.
  return {
    kinemages: []
  };
}