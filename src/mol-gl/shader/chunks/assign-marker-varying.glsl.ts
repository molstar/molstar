

import { loadCheckpoint } from '../../../mol-util/debug';
loadCheckpoint(`mol-gl/shader/chunks/assign-marker-varying.glsl.ts::start`);
export const assign_marker_varying = `
#if defined(dNeedsMarker)
    #if defined(dMarkerType_instance)
        vMarker = readFromTexture(tMarker, aInstance, uMarkerTexDim).a;
    #elif defined(dMarkerType_groupInstance)
        vMarker = readFromTexture(tMarker, aInstance * float(uGroupCount) + group, uMarkerTexDim).a;
    #endif
#endif
`;
loadCheckpoint(`mol-gl/shader/chunks/assign-marker-varying.glsl.ts::end`);
