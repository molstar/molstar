export const assign_marker_varying = `
vMarker = readFromTexture(tMarker, aInstance * float(uGroupCount) + group, uMarkerTexDim).a;
`;