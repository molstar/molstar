
export interface LineGraphParams {
    height: number,
    width: number,
    padding: number,
    paddingYUnnormalized: number,
    paddingXUnnormalized: number,
    offsetY: number,
}


export const LineGraphParams: LineGraphParams = {
    height: 400,
    width: 600,
    padding: 70,
    // TODO: function that caclulates this based on height etc.
    // or rather contant numbers as const vars above that 
    paddingYUnnormalized: 70 / 400,
    paddingXUnnormalized: 70 / 400,
    // TODO: add offset unnormalized or
    offsetY: 35,
};
