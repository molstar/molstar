import { ParamDefinition, PluginStateObject } from "../meshes/molstar-lib-imports";

export const CellStarStateParams = {
    opacity: ParamDefinition.Numeric(1, { min: 0, max: 1, step: 0.05 }),
}
export type CellStarStateData = ParamDefinition.Values<typeof CellStarStateParams>;

// interface CellStarStateData {
// }

export class CellStarState extends PluginStateObject.Create<CellStarStateData>({ name: 'CellStar Entry State', typeClass: 'Data' }) { }

export const CELLSTAR_STATE_FROM_ENTRY_TRANSFORMER_NAME = 'cellstar-state-from-entry'; // defined here to avoid cyclic dependency