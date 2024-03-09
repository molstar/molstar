import { loadCVSXFromAnything } from '.';
import { hashFnv32a } from '../../../mol-data/util';
import { DataFormatProvider } from '../../../mol-plugin-state/formats/provider';
import { PluginStateObject, PluginStateObject as SO } from '../../../mol-plugin-state/objects';
import { Download } from '../../../mol-plugin-state/transforms/data';
import { PluginContext } from '../../../mol-plugin/context';
import { StateAction, StateObjectRef, StateObjectSelector } from '../../../mol-state';
import { RuntimeContext, Task } from '../../../mol-task';
import { Asset, AssetManager } from '../../../mol-util/assets';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { unzip } from '../../../mol-util/zip/zip';

/** Data format provider for CVSX format.
 */
export const CVSXFormatProvider: DataFormatProvider<{}, StateObjectSelector<PluginStateObject.Data.Binary>, any> = DataFormatProvider({
    label: 'CVSX',
    description: 'CVSX',
    category: 'Miscellaneous',
    binaryExtensions: ['cvsx'],
    parse: async (plugin, data) => {
        return loadCVSXFromAnything(plugin, data);
    },
});