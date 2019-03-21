interface AccessibleSurfaceArea {
    readonly atomRadius: ArrayLike<number>,
    readonly accessibleSurfaceArea: ArrayLike<number>,
    readonly relativeAccessibleSurfaceArea: ArrayLike<number>,
    readonly buried: any
}

interface AccessibleSurfaceAreaComputationParameters {
    numberOfSpherePoints: number
    probeSize: number
}

export { AccessibleSurfaceArea, AccessibleSurfaceAreaComputationParameters }