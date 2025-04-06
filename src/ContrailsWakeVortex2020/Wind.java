package ContrailsWakeVortex2020;

/**
 * Created by s9236912 on 25.04.2019.
 */
public class Wind {


    public static double getHeadingFromDirection(double direction_deg){
        if (direction_deg>=180.0){
            return direction_deg-180.0;
        }else{
            return direction_deg+180.0;
        }
    }


    public static double getDirectionFromHeading(double heading_deg){
        if (heading_deg>=180.0){
            return heading_deg-180.0;
        }else{
            return heading_deg+180.0;
        }
    }

    /**
     * wind flows towards this direction
     *
     * @param windU
     * @param windV
     * @return in degree
     */
    public static double getWindheading(double windU, double windV) {
        /* in [-90;90] */
        double hdg = 90.0 - Math.toDegrees(Math.atan(windV / windU));

        /* turn by 180 deg if east-bound direction is negative */
        if (windU < 0) {
            hdg += 180.0;
        }

        if (hdg < 0.0) hdg += 360.0;
        if (hdg > 360.0) hdg -= 360.0;
        return hdg;
    }

    /**
     * positive: headwind, negative: tailwind
     * @param acHeading_deg
     * @param windspeed_mps
     * @param windHeading_deg
     * @return
     */
    public static double getHeadwind(double acHeading_deg, double windspeed_mps, double windHeading_deg) { //acHeading in Degree, returns Headwind Component in meter per second for the given acHeading (tailwinds are negative)

        //wind is in acHeading (direction in which wind blows)
        //convert to usual direction

        windHeading_deg=windHeading_deg-180.0;

        if (windHeading_deg < 0.0) windHeading_deg += 360.0;
        if (windHeading_deg > 360.0) windHeading_deg -= 360.0;

        double headwind=windspeed_mps*Math.cos(Math.toRadians((windHeading_deg)-acHeading_deg));
        return headwind;
    }

    public static double[] getHeadWindAndWindCorrectionAngle(double windU, double windV, double airspeed_mps, double heading) {

        Windprops windprops=getWinddirectionAndSpeed(windU,windV);
        double HeadWindWCA[]={0.0,0.0};
        HeadWindWCA[0]=getHeadwindComponentMperSec(heading, windprops.getSpeed_mps(), windprops.getDirection_deg());
        HeadWindWCA[1] = getWindCorrectionAngle(heading, windprops.getSpeed_mps(), windprops.getDirection_deg(), airspeed_mps);
        return HeadWindWCA;
    }

    private static double getHeadwindComponentMperSec(double heading, double windspeed, double windDirection) { //Heading in Degree, returns Headwind Component in meter per second for the given Heading (tailwinds are negative)

        double headwind=windspeed*Math.cos(Math.toRadians((windDirection)-heading));
        return headwind;
    }

    private static double getWindCorrectionAngle(double heading, double windspeed, double windDirection, double airspeed) { //correct Heading with this value to compensate wind drift

        double WCA=windspeed/airspeed* Math.sin(Math.toRadians((windDirection)-heading));// airspeed in meter per second

        //arc sin
        WCA=Math.toDegrees(Math.asin(WCA)); //positive values: correct by positive degrees
        return WCA;
    }

    public static Windprops getWinddirectionAndSpeed(double WindU, double WindV) { // alt in meter

        double windspeed=Math.sqrt(WindU*WindU+WindV*WindV);
        double winddirection=Math.floor(((180.0)+Math.toDegrees((Math.atan2(WindU, WindV))))*100.0)/100.0;

        Windprops windprops=new Windprops(winddirection,windspeed);

        //Direction: 90° corresponds to wind from east, 180° from south, 270° from west and 360° wind from north. 0° is used for no wind

        return windprops;
    }

    public static class Windprops{

        double mdirection_deg;
        double mspeed_mps;


        public Windprops(double direction_deg,double speed_mps){
            this.mdirection_deg=direction_deg;
            this.mspeed_mps=speed_mps;
        }


        public double getDirection_deg() {
            return mdirection_deg;
        }

        public double getSpeed_mps() {
            return mspeed_mps;
        }

    }
}
