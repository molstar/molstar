import { PluginStateObject } from '../../../mol-plugin-state/objects';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';

import { Choice } from '../choice';


export const DEFAULT_MESH_SERVER = 'http://localhost:9000/v2';


export class MeshServerInfo extends PluginStateObject.Create<MeshServerInfo.Data>({ name: 'Volume Server', typeClass: 'Object' }) { }

export namespace MeshServerInfo {
    export const MeshSourceChoice = new Choice({ empiar: 'EMPIAR', emdb: 'EMDB' }, 'empiar');
    export type MeshSource = Choice.Values<typeof MeshSourceChoice>;

    export const Params = {
        serverUrl: PD.Text(DEFAULT_MESH_SERVER),
        source: MeshSourceChoice.PDSelect(),
        entryId: PD.Text(''),
    };
    export type Data = PD.Values<typeof Params>;
}
