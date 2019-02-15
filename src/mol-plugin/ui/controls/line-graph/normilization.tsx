import { Vec2 } from "mol-math/linear-algebra";

export function normalize(height:number, width:number, point: Vec2, padding:number=0) {
    const min = padding/2;
    const maxX = width+min;
    const maxY = height+min;
    const normalizedX = (point[0]*(maxX-min))+min;
    const normalizedY = (point[1]*(maxY-min))+min;
    
    let reverseY;
    reverseY = height+padding-normalizedY;
    const newPoint = Vec2.create(normalizedX, reverseY);
    return newPoint;
}

export function unNormalize(height:number, width: number, point: Vec2, padding:number=0) {
    const min=padding/2;
    const maxX = width+min;
    const maxY = height+min;
    const unNormalizedX = (point[0]-min)/(maxX-min);
    // We have to take into account that we reversed y when we first normalized it
    const unNormalizedY = ((height+padding)-point[1]-min)/(maxY-min);
    return Vec2.create(unNormalizedX, unNormalizedY);
}