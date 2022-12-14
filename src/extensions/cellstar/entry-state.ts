import { ParamDefinition, PluginStateObject } from '../meshes/molstar-lib-imports';

export const CellstarStateParams = {
    opacity: ParamDefinition.Numeric(1, { min: 0, max: 1, step: 0.05 }),
    selectedSegment: ParamDefinition.Numeric(-1, { step: 1 }),
    visibleSegments: ParamDefinition.ObjectList({ segmentId: ParamDefinition.Numeric(0) }, s => s.segmentId.toString()),
    visibleModels: ParamDefinition.ObjectList({ pdbId: ParamDefinition.Text('') }, s => s.pdbId.toString()),
};
export type CellstarStateData = ParamDefinition.Values<typeof CellstarStateParams>;


export class CellstarState extends PluginStateObject.Create<CellstarStateData>({ name: 'Vol & Seg Entry State', typeClass: 'Data' }) { }


export const CELLSTAR_STATE_FROM_ENTRY_TRANSFORMER_NAME = 'cellstar-state-from-entry'; // defined here to avoid cyclic dependency
