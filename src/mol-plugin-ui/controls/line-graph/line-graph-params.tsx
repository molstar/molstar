
export interface LineGraphParams {
    height: number,
    width: number,
    padding: number,
    paddingYUnnormalized: number,
    paddingXUnnormalized: number,
    // offsetY: number,
    baseline: number
    baselineUnnormalized: number
    // real svg values (400, 600 etc.)
    roof: number,
    // [0; 1]
    roofUnnormalized: number
}


export const LineGraphParams: LineGraphParams = {
    height: 400,
    width: 600,
    padding: 70,
    // TODO: function that caclulates this based on height etc.
    // or rather contant numbers as const vars above that 
    paddingYUnnormalized: 70 / 400,
    paddingXUnnormalized: 70 / 600,
    // TODO: add offset unnormalized or
    // offsetY: 35,
    baseline: 35,
    baselineUnnormalized: 35 / 400,
    roof: 35,
    // [0; 1]
    roofUnnormalized: 35 / 400
};
