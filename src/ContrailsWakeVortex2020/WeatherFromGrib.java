package ContrailsWakeVortex2020;
//import de.tudresden.ifl.coala.tools.interpolate;
//import de.tudresden.ifl.util.Constants;
//import Wind;

//import org.apache.commons.logging.*;
//import org.apache.commons.logging.impl.LogFactoryImpl;
import ucar.nc2.grib.grib2.*;
import ucar.unidata.io.RandomAccessFile;

import java.io.IOException;
import java.time.ZoneId;
import java.time.ZonedDateTime;
import java.util.Map;
import java.util.TreeMap;

import  org.slf4j.Logger;
import  org.slf4j.LoggerFactory;

import org.apache.commons.logging.Log;

public class WeatherFromGrib {


    public static void main(String[] args)
    {
        WeatherFromGrib TEST = new WeatherFromGrib("C:\\Users\\mluo\\eclipse-workspace\\ContrailsWakeVortex2020-master\\grib_01_2016-02-07_1200.grib2");

        //double testest = TEST.getWindUFromAltitude(50,10,1000);
        double testest = TEST.getTempFromAltitude(-50,150,10464);
        //System.out.println(testest);
        //double testest1 = TEST.getThetaFromAltitude(51,11,16500);
        //System.out.println(testest1);
        //double testest2 = TEST.getThetaFromPressure(88,1,80000);
        System.out.println(testest);
        //System.out.println(TEST.getAltitudeFromPressure(50,100,0));
    }
    
    Logger logger = LoggerFactory.getLogger(WeatherFromGrib.class);
    //public final static Log LOG= LogFactoryImpl.getLogger(WeatherFromGrib.class);

    public final double La1;
    public final double Lo1;
    public final double dLat;
    public final double dLon;
    public final int nLon;

    public final TreeMap<Double, float[]> mTemperature;// contains temperature in K
    public final TreeMap<Double, float[]> mTheta; //0,0,3
    public final TreeMap<Double, float[]> mRelativeHumidity;// contains relative Humidity in %
    public final TreeMap<Double, float[]> mWind_u;// contains u component of the wind in m/s
    public final TreeMap<Double, float[]> mWind_v;// contains v component of the wind in m/s
    public final TreeMap<Double, float[]> mGeoPotentialHeight;// contains geopotential height in m
    public final TreeMap<Double, float[]> mWindSpeed;// contains windspeed component of the wind in m/s
    public final TreeMap<Double, float[]> mWindHeading;// contains windheading ("flows to") component of the wind in degree
    public final TreeMap<Double, float[]> mvvel;// vertical velocity 0,2,8


    public final Timestamp mTimestamp;

    public final static int LEVEL_TYPE_ISOBARIC_SURFACE = 100;

    public final static int DISCIPLINE_METEOROLOGICAL = 0;

    public final static int CATEGORY_TEMPERATURE = 0;
    public final static int CATEGORY_MOISTURE = 1;
    public final static int CATEGORY_MOMENTUM = 2;
    public final static int CATEGORY_MASS = 3;
    public final static int CATEGORY_TEST = 6;

    public final static int NUMBER_TEMPERATURE = 0;
    public final static int NUMBER_RELATIVE_HUMIDITY = 1;
    public final static int NUMBER_WIND_U_COMPONENT = 2;
    public final static int NUMBER_THETA_E = 2;
    public final static int NUMBER_WIND_V_COMPONENT = 3;
    public final static int NUMBER_GEOPOTENTIAL_HEIGHT = 5;
    public final static int NUMBER_VERTICAL_VELOCITY = 8;


    /*
    public WeatherFromGrib(String gribFileName, WeatherUncertainty uncertainty) {
        this(gribFileName);

        for (WeatherUncertainty.WeatherParameter wp : WeatherUncertainty.WeatherParameter.values()) {
            TreeMap<Double, float[]> valueList = null;
            switch (wp) {
                case RH: valueList = mRelativeHumidity; break;
                case TEMP: valueList = mTemperature; break;
                case WIND_U: valueList = mWind_u; break;
                case WIND_V: valueList = mWind_v; break;
            }

            for (Map.Entry<Double, float[]> values : valueList.entrySet()) {
                float[] value = values.getValue();
                for (int i = 0; i < value.length; ++i) {
                    value[i] += uncertainty.nextUncertaintyValue(wp);
                }
            }
        }
    }
    */

    public WeatherFromGrib(TreeMap<Double, float[]> temperature,
                           TreeMap<Double, float[]> theta,
                           TreeMap<Double, float[]> relativeHumidity,
                           TreeMap<Double, float[]> wind_u,
                           TreeMap<Double, float[]> wind_v,
                           TreeMap<Double, float[]> geoPotentialHeight,
                           TreeMap<Double, float[]> windSpeed,
                           TreeMap<Double, float[]> windHeading,
                           TreeMap<Double, float[]> vvel,
                           double La1,
                           double Lo1,
                           double dLat,
                           double dLon,
                           int nLon,
                           Timestamp timestamp) {
        this.mTemperature = temperature;
        this.mTheta = theta;
        this.mRelativeHumidity = relativeHumidity;
        this.mWind_u = wind_u;
        this.mWind_v = wind_v;
        this.mGeoPotentialHeight = geoPotentialHeight;
        this.mWindSpeed = windSpeed;
        this.mWindHeading = windHeading;
        this.mvvel = vvel;
        this.La1 = La1;
        this.Lo1 = Lo1;
        this.dLat = dLat;
        this.dLon = dLon;
        this.nLon = nLon;
        this.mTimestamp = timestamp;
    }

    public WeatherFromGrib(String gribFileName) {
        RandomAccessFile gribDataFile = null;
        Grib2RecordScanner recordScanner = null;
        try {
            gribDataFile = new RandomAccessFile(gribFileName, "r");
            gribDataFile.order(RandomAccessFile.BIG_ENDIAN);
            recordScanner = new Grib2RecordScanner(gribDataFile);
        } catch (IOException e) {
            //Log.error("Could not open GRIB file: " + gribFileName, e);
        }

        // ---== read WeatherGrib ==---
        mTemperature = new TreeMap<>();
        mTheta = new TreeMap<>();
        mRelativeHumidity = new TreeMap<>();
        mWind_u = new TreeMap<>();
        mWind_v = new TreeMap<>();
        mGeoPotentialHeight = new TreeMap<>();
        mWindSpeed = new TreeMap<>();
        mWindHeading = new TreeMap<>();
        mvvel = new TreeMap<>();
        Grib2Record record = null;

        try {
            while (recordScanner != null && recordScanner.hasNext()) {
                record = recordScanner.next();
                Grib2SectionIndicator is = record.getIs();
                Grib2SectionIdentification id = record.getId();
                Grib2Pds pdsv = record.getPDS();               

                if ((is.getDiscipline() != DISCIPLINE_METEOROLOGICAL) || (pdsv.getLevelType1() != LEVEL_TYPE_ISOBARIC_SURFACE)) {
                    continue;
                }

                Map<Double,float[]> targetMap = null;
                switch (pdsv.getParameterCategory()) {
                    case CATEGORY_TEMPERATURE:
                        switch (pdsv.getParameterNumber()) {
                            case NUMBER_TEMPERATURE: targetMap = mTemperature; break;
                            case NUMBER_THETA_E: targetMap = mTheta; break;
                        }
                        break;
                    case CATEGORY_MOISTURE:
                        switch (pdsv.getParameterNumber()) {
                            case NUMBER_RELATIVE_HUMIDITY: targetMap = mRelativeHumidity; break;
                        }
                        break;
                    case CATEGORY_MOMENTUM:
                        switch (pdsv.getParameterNumber()) {
                            case NUMBER_WIND_U_COMPONENT: targetMap = mWind_u; break;
                            case NUMBER_WIND_V_COMPONENT: targetMap = mWind_v; break;
                            case NUMBER_VERTICAL_VELOCITY: targetMap = mvvel; break;
                        }
                        break;
                    //case CATEGORY_TEST:
                      //  switch (pdsv.getParameterNumber()) {
                        //    case NUMBER_THETA_E: targetMap = mTheta; break;
                        //}
                        //break;
                    case CATEGORY_MASS:
                        switch (pdsv.getParameterNumber()) {
                            case NUMBER_GEOPOTENTIAL_HEIGHT: targetMap = mGeoPotentialHeight;
                        }
                        break;
                }

                if (targetMap != null) {
                    extractValuesFromLevel(record, pdsv, gribDataFile, targetMap);
                }
            }
        } catch (IOException e) {
          //Log.error("Error while reading grib data");
            e.printStackTrace();
        }

        // ---== read grid information ==---
        if (record != null && !record.getGDS().isLatLon()) {
            //Log.error("Error: GRIB2 GDS is not of type LatLon!");
        }

        Grib2Gds.LatLon gdsv = (Grib2Gds.LatLon) record.getGDS();
        this.La1 = gdsv.la1;
        this.Lo1 = gdsv.lo1;
        double lo2 = gdsv.lo2;
        this.dLat = Math.abs(gdsv.deltaLat);
        this.dLon = Math.abs(gdsv.deltaLon);
        this.nLon = (int) ((lo2 - Lo1) / dLon) + 1;

        ZonedDateTime date = ZonedDateTime.of(record.getId().getYear(), record.getId().getMonth(), record.getId().getDay(), record.getId().getHour(), record.getId().getMinute(), record.getId().getSecond(), 0, ZoneId.of("UCT"));

        mTimestamp = new Timestamp(date, record.getId().getHour(), record.getPDS().getForecastTime());
        addwindspeed();
    }

    public void extractValuesFromLevel(Grib2Record record, Grib2Pds pdsv, RandomAccessFile gribDataFile, Map<Double, float[]> targetMap) {
        float[] values = null;
        try {
            values = record.readData(gribDataFile);
        } catch (IOException e) {
          //Log.error("Failed to read data from record.");
        }

        double level1 = pdsv.getLevelValue1();
        targetMap.put(level1, values);
    }

    public void addwindspeed() {
        if (mWind_v.size() != mWind_u.size()) {
          //Log.error("v values size (" + mWind_v.size() + ") not equal to u values size (" + mWind_u.size() + ")");
        }

        for (int i = 0; i < mWind_v.size(); ++i) {
            double layer_v = (double)mWind_v.keySet().stream().toArray()[i];
            double layer_u = (double)mWind_u.keySet().stream().toArray()[i];

            float[] wind_v = mWind_v.get(layer_v);
            float[] wind_u = mWind_u.get(layer_u);
            float[] speed = new float[wind_v.length];
            float[] heading = new float[wind_v.length];

            if (wind_v.length != wind_u.length) {
              //Log.error("v values size (" + wind_v.length + ") not equal to u values size (" + wind_u.length + ")");
            }
            for (int it = 0; it < wind_v.length; ++it) {
                float u = wind_u[it];
                float v = wind_v[it];

                double s = Math.sqrt(u * u + v * v);
                speed[it] = (float)s;

                double hdg = Wind.getWindheading(u, v);
                heading[it] = (float) hdg;
            }
            mWindSpeed.put(layer_v, speed);
            mWindHeading.put(layer_v, heading);
        }
    }

    /**
     * For a given coordinate, retrieve the index of the cell
     * within the underlying grid, which is represented as a list.
     *
     * @param lat Latitude in [deg] [-90,90]
     * @param lon longitude in [deg] [0,360]
     * @return the index within the internal list.
     */
    int getIndex(double lat, double lon) {
        if (lat < -90) {
            lon += 180;
            lat = (-180 - lat);
        }

        if (lat > 90) {
            lon += 180;
            lat = (180 - lat);
        }

        if (lon > 360) {
            lon -= 360;
        }

        if (lon < 0) {
            lon += 360;
        }

        if (lat < -90 || lat > 90) {
            throw new IllegalArgumentException("Latitude out of range (" + lat + ")");
        }

        if (lon < 0 || lon > 360) {
            throw new IllegalArgumentException("Longitude out of range (" + lon + ")");
        }

        int iLat = (int) Math.round((La1 - lat) / dLat);
        int iLon = (int) Math.round((lon - Lo1) / dLon);
        if (iLon == 360) {
            iLon = 0; // FIXME is this correct for any arbitrary resolution?
        }

        return iLat * nLon + iLon;
    }

    /**
     * Calculate height at which the given pressure prevails.
     *
     * @param lat latitude in degress [-90,90]
     * @param lon longitude in degrees [-180,180]
     * @param p   pressure in [Pa]
     * @return height in [m]
     */
    public double getAltitudeFromPressure(double lat, double lon, double p) {
        // double p must be in Pa but is caused by coordinate initilization sometimes 0m (at dep/dest airport)
        //the method below converts the airport altitude=0 to the surface press 10130pa
        if (Double.compare(p, 0.0) == 0) {
            p = 101300.0;
        }

        double gph = getValueFromPressure(lat, lon, p, mGeoPotentialHeight);
        double c = (1. - 2.6373e-3 * Math.cos(2. * lat * Math.PI / 180.) + 5.9e-6 * Math.pow(Math.cos(2. * lat * Math.PI / 180.), 2));

        return 1. / (c / gph - 1. / Constants.EARTH_RADIUS_M);
    }


    public double getWindHeadingFromAltitude(double lat, double lon, double alt) {
        return getValueFromAltitude_windheading(lat, lon, alt, mWindHeading);
    }

    public double getWindHeadingFromPressure(double lat, double lon, double p) {
        return getValueFromPressure(lat, lon, p, mWindHeading);
    }

    public double getWindspeedFromPressure(double lat, double lon, double p) {
        return getValueFromPressure(lat, lon, p, mWindSpeed);
    }

    public double getWindspeedFromAltitude(double lat, double lon, double alt) {
        return getValueFromAltitude(lat, lon, alt, mWindSpeed);
    }

    /**
     * Retrieve temperature at the given pressure altitude.
     *
     * @param lat latitude in degress [-90,90]
     * @param lon longitude in degrees [-180,180]
     * @param p   pressure in [Pa]
     * @return temperature in [K]
     */
    public double getTempFromPressure(double lat, double lon, double p) {
        return getValueFromPressure(lat, lon, p, mTemperature);
    }

    public double getThetaFromPressure(double lat, double lon, double p) {
        return getValueFromPressure(lat, lon, p, mTheta);
    }

    public double getVvelFromPressure(double lat, double lon, double p) {
        return getValueFromPressure(lat, lon, p, mvvel);
    }

    /**
     * Retrieve the temperature at the given coordinates and altitude.
     * A 3D-interpolation will be conducted.
     *
     * @param lat latitude [-90;90]
     * @param lon longitue [-180;180]
     * @param alt altitude in [m]
     * @return temperature in [K]
     */
    public double getTempFromAltitude(double lat, double lon, double alt) {
        return getValueFromAltitude(lat, lon, alt, mTemperature);
    }

    public double getThetaFromAltitude(double lat, double lon, double alt) {
        return getValueFromAltitude(lat, lon, alt, mTheta);
    }

    public double getVvelFromAltitude(double lat, double lon, double alt) {
        return getValueFromAltitude(lat, lon, alt, mvvel);
    }
    /**
     * Retrieve the relative humidity at the given coordinates and altitude.
     * A 3D-interpolation will be conducted.
     *
     * @param lat latitude [-90;90]
     * @param lon longitue [-180;180]
     * @param alt altitude in [m]
     * @return relative humidity in [%]
     */
    public double getHumidityFromAltitude(double lat, double lon, double alt) {
        return getValueFromAltitude(lat, lon, alt, mRelativeHumidity);
    }

    /**
     * Retrieve the wind's u-component at the given coordinates and altitude.
     * A 3D-interpolation will be conducted.
     *
     * @param lat latitude [-90;90]
     * @param lon longitue [-180;180]
     * @param alt altitude in [m]
     * @return u-component of wind in [m/s]
     */
    public double getWindUFromAltitude(double lat, double lon, double alt) {
        return getValueFromAltitude(lat, lon, alt, mWind_u);
    }

    /**
     * Retrieve the wind's v-component at the given coordinates and altitude.
     * A 3D-interpolation will be conducted.
     *
     * @param lat latitude [-90;90]
     * @param lon longitue [-180;180]
     * @param alt altitude in [m]
     * @return v-component of wind in [m/s]
     */
    public double getWindVFromAltitude(double lat, double lon, double alt) {
        return getValueFromAltitude(lat, lon, alt, mWind_v);
    }

    /**
     * Retrieve the 3D-interpolated pressure at the given coordinates and for the given altitude.
     *
     * @param lat latitude [-90;90]
     * @param lon longitude [-180;180]
     * @param alt the altitude at which the value should be retrieved
     * @return the pressure at coordinates lat,lon and altitude alt in [hPa]
     */
    public double getPressureFromAltitude(double lat, double lon, double alt) {
        if (lon < 0) {
            lon += 360.0;
        }

        final double lat0 = lat - Math.abs(lat % this.dLat);
        final double lat1 = Math.min(90.0, lat0 + this.dLat);
        final double lon0 = lon - Math.abs(lon % this.dLon);
        final double lon1 = Math.min(360.0, lon0 + this.dLon);
        final double w00 = this.getPressureFromAltitudeOnGrid(lat0, lon0, alt);
        final double w01 = this.getPressureFromAltitudeOnGrid(lat0, lon1, alt);
        final double w10 = this.getPressureFromAltitudeOnGrid(lat1, lon0, alt);
        final double w11 = this.getPressureFromAltitudeOnGrid(lat1, lon1, alt);
        final double wlat0 = lon1 == lon0 ? w00 : w00 + (w01 - w00) / (lon1 - lon0) * (lon - lon0);
        final double wlat1 = lon1 == lon0 ? w10 : w10 + (w11 - w10) / (lon1 - lon0) * (lon - lon0);
        return lat1 == lat0 ? wlat0 : wlat0 + (wlat1 - wlat0) / (lat1 - lat0) * (lat - lat0);
    }

    /**
     * Retrieve the 3D-interpolated value at the given coordinates and for the given altitude.
     *
     * @param lat     latitude [-90;90]
     * @param lon     longitude [-180;180]
     * @param alt     the altitude at which the value should be retrieved
     * @param valList a list of values from which one will be chosen according to the altitude
     * @return a value from valList at coordinates lat,lon and altitude alt
     */
    public double getValueFromAltitude(double lat, double lon, double alt, TreeMap<Double, float[]> valList) {
        if (lon < 0) {
            lon += 360.0;
        }

        final double lat0 = lat - Math.abs(lat % this.dLat);
        final double lat1 = Math.min(90.0, lat0 + this.dLat);
        final double lon0 = lon - Math.abs(lon % this.dLon);
        final double lon1 = Math.min(360.0, lon0 + this.dLon);
        final double w00 = this.getValueFromAltitudeOnGrid(lat0, lon0, alt, valList);
        final double w01 = this.getValueFromAltitudeOnGrid(lat0, lon1, alt, valList);
        final double w10 = this.getValueFromAltitudeOnGrid(lat1, lon0, alt, valList);
        final double w11 = this.getValueFromAltitudeOnGrid(lat1, lon1, alt, valList);

        final double wlat0 = lon1 == lon0 ? w00 : w00 + (w01 - w00) / (lon1 - lon0) * (lon - lon0);
        final double wlat1 = lon1 == lon0 ? w10 : w10 + (w11 - w10) / (lon1 - lon0) * (lon - lon0);
        return lat1 == lat0 ? wlat0 : wlat0 + (wlat1 - wlat0) / (lat1 - lat0) * (lat - lat0);
    }

    /**
     * Retrieve the 3D-interpolated wind heading at the given coordinates and for the given altitude.
     *
     * @param lat     latitude [-90;90]
     * @param lon     longitude [-180;180]
     * @param alt     the altitude at which the value should be retrieved
     * @param valList a list of values from which one will be chosen according to the altitude
     * @return a value from valList at coordinates lat,lon and altitude alt
     */
    public double getValueFromAltitude_windheading(double lat, double lon, double alt, TreeMap<Double, float[]> valList) {
        if (lon < 0) {
            lon += 360.0;
        }

        final double lat0 = lat - Math.abs(lat % this.dLat);
        final double lat1 = Math.min(90.0, lat0 + this.dLat);
        final double lon0 = lon - Math.abs(lon % this.dLon);
        final double lon1 = Math.min(360.0, lon0 + this.dLon);
        final double w00 = this.getValueFromAltitudeOnGrid_windheading(lat0, lon0, alt, valList);
        final double w01 = this.getValueFromAltitudeOnGrid_windheading(lat0, lon1, alt, valList);
        final double w10 = this.getValueFromAltitudeOnGrid_windheading(lat1, lon0, alt, valList);
        final double w11 = this.getValueFromAltitudeOnGrid_windheading(lat1, lon1, alt, valList);
        return interpolatewindheading(lat, lon, lat0, lon0, lat1, lon1, w00, w01, w10, w11);
    }

    protected static double interpolatewindheading(double lat, double lon, double lat0, double lon0, double lat1, double lon1, double w00, double w01, double w10, double w11) {
        double amount_wlon0 = (lat - lat0) / (lat1 - lat0);
        double intlon0 = interpolate.lerpDegrees(w00, w10, amount_wlon0);
        double amount_wlon1 = (lat - lat0) / (lat1 - lat0);
        double intlon1 = interpolate.lerpDegrees(w01, w11, amount_wlon1);
        double amount_wlat1 = (lon - lon0) / (lon1 - lon0);
        return interpolate.lerpDegrees(intlon0, intlon1, amount_wlat1);
    }

    /**
     * Retrieve the 3D-interpolated value at the given coordinates and for the given pressure.
     *
     * @param lat     latitude [-90;90]
     * @param lon     longitude [-180;180]
     * @param p       the pressure at which the value should be retrieved
     * @param valList a list of values from which one will be chosen according to the pressure
     * @return a value from valList at coordinates lat,lon and pressure p
     */
    public double getValueFromPressure(double lat, double lon, double p, TreeMap<Double, float[]> valList) {
        if (lon < 0) {
            lon += 360.0;
        }

        final double lat0 = lat - Math.abs(lat % this.dLat);
        final double lat1 = Math.min(90.0, lat0 + this.dLat);
        final double lon0 = lon - Math.abs(lon % this.dLon);
        final double lon1 = Math.min(360.0, lon0 + this.dLon);
        final double w00 = this.getValueFromPressureOnGrid(lat0, lon0, p, valList);
        final double w01 = this.getValueFromPressureOnGrid(lat0, lon1, p, valList);
        final double w10 = this.getValueFromPressureOnGrid(lat1, lon0, p, valList);
        final double w11 = this.getValueFromPressureOnGrid(lat1, lon1, p, valList);

        final double wlat0 = lon1 == lon0 ? w00 : w00 + (w01 - w00) / (lon1 - lon0) * (lon - lon0);
        final double wlat1 = lon1 == lon0 ? w10 : w10 + (w11 - w10) / (lon1 - lon0) * (lon - lon0);
        return lat1 == lat0 ? wlat0 : wlat0 + (wlat1 - wlat0) / (lat1 - lat0) * (lat - lat0);
    }

    /**
     * Retrieve the (interpolated) pressure at the given coordinates and for the given altitude.
     * Interpolation will be conducted over altitude scale only.
     * There is no interpolation between coordinates.
     * For this, refer to {@link WeatherFromGrib#getValueFromAltitude(double, double, double, TreeMap)}
     *
     * @param lat latitude [-90;90]
     * @param lon longitude [-180;180]
     * @param alt the altitude at which the value should be retrieved in [m]
     * @return the pressure at altitude alt in [hPa]
     */
    public double getPressureFromAltitudeOnGrid(double lat, double lon, double alt) {
        // lon is transferred to 0...360 whereas 0...180 is identical to above
        if (lon < 0) {
            lon += 360.;
        }
        int c = getIndex(lat, lon);

        // transfer data from vector to array (for interpolation)
        int nv = mGeoPotentialHeight.size();

        double[] p = new double[nv];
        double[] gph = new double[nv];

        // build list from behind, since small altitudes correspond to high pressures
        for (int i = 0; i < nv; ++i) {
            gph[nv - i - 1] = mGeoPotentialHeight.get(mGeoPotentialHeight.keySet().stream().toArray()[i])[c];
            p[nv - i - 1] = (double) mGeoPotentialHeight.keySet().stream().toArray()[i];
        }

        // calculate geopotantial height of altitude.
        double alt_gph = (1. - 2.6373e-3 * Math.cos(2 * lat * Math.PI / 180.) + 5.9e-6 * Math.pow(Math.cos(2 * lat * Math.PI / 180.), 2)) * (alt / (1 + alt / Constants.EARTH_RADIUS_M));
        return interpolate.linear(gph, p, alt_gph, true) / 100.;
    }

    /**
     * Retrieve the (interpolated) value at the given coordinates and for the given altitude.
     * Interpolation will be conducted over altitude scale only.
     * There is no interpolation between coordinates.
     * For this, refer to {@link WeatherFromGrib#getValueFromAltitude(double, double, double, TreeMap)}
     *
     * @param lat     latitude [-90;90]
     * @param lon     longitude [-180;180]
     * @param alt     the altitude at which the value should be retrieved in [m]
     * @param valList a list of values from which one will be chosen according to the altitude
     * @return a value from valList at altitude alt
     */
    public double getValueFromAltitudeOnGrid_windheading(double lat, double lon, double alt, TreeMap<Double, float[]> valList) {
        if (lon < 0) {
            lon += 360.0;
        }

        int c = getIndex(lat, lon);
        // transfer data from vector to array (for interpolation)
        int nv = valList.size(); //how many milibar layers
        double[] gph = new double[nv];
        double[] val = new double[nv];

        // for altitude, retrieve the value from valList
        // build list from behind, since small altitudes correspond to hight pressures
        for (int i = 0; i < nv; ++i) {
            gph[nv - 1 - i] = mGeoPotentialHeight.get(mGeoPotentialHeight.keySet().stream().toArray()[i])[c];
            val[nv - 1 - i] = valList.get(valList.keySet().stream().toArray()[i])[c];
        }
        // calculate geopotantial height of altitude.
        double alt_gph = (1. - 2.6373e-3 * Math.cos(2 * lat * Math.PI / 180.) + 5.9e-6 * Math.pow(Math.cos(2 * lat * Math.PI / 180.), 2)) * (alt / (1 + alt / Constants.EARTH_RADIUS_M));

        return interpolate.linear_windheading(gph, val, alt_gph, true);
    }

    /**
     * Retrieve the (interpolated) value at the given coordinates and for the given altitude.
     * Interpolation will be conducted over altitude scale only.
     * There is no interpolation between coordinates.
     * For this, refer to {@link WeatherFromGrib#getValueFromAltitude(double, double, double, TreeMap)}
     *
     * @param lat     latitude [-90;90]
     * @param lon     longitude [-180;180]
     * @param alt     the altitude at which the value should be retrieved in [m]
     * @param valList a list of values from which one will be chosen according to the altitude
     * @return a value from valList at altitude alt
     */
    public double getValueFromAltitudeOnGrid(double lat, double lon, double alt, TreeMap<Double, float[]> valList) {
        if (lon < 0) {
            lon += 360.0;
        }

        int c = getIndex(lat, lon);
        // transfer data from vector to array (for interpolation)
        int nv = valList.size(); //how many milibar layers
        double[] gph = new double[nv];
        double[] val = new double[nv];

        // for altitude, retrieve the value from valList
        // build list from behind, since small altitudes correspond to hight pressures
        for (int i = 0; i < nv; ++i) {
            gph[nv - 1 - i] = mGeoPotentialHeight.get(mGeoPotentialHeight.keySet().stream().toArray()[i])[c];
            val[nv - 1 - i] = valList.get(valList.keySet().stream().toArray()[i])[c];
        }
        // calculate geopotantial height of altitude.
        double alt_gph = (1. - 2.6373e-3 * Math.cos(2 * lat * Math.PI / 180.) + 5.9e-6 * Math.pow(Math.cos(2 * lat * Math.PI / 180.), 2)) * (alt / (1 + alt / Constants.EARTH_RADIUS_M));

        return interpolate.linear(gph, val, alt_gph, true);
    }

    /**
     * Retrieve the (interpolated) value at the given coordinates and for the given pressure.
     * Interpolation will be conducted over pressure scale only.
     * There is no interpolation between coordinates.
     * For this, refer to {@link WeatherFromGrib#getValueFromPressure(double, double, double, TreeMap)}
     *
     * @param lat      latitude [-90;90]
     * @param lon      longitude [-180;180]
     * @param pressure the pressure at which the value should be retrieved in [Pa]
     * @param valList  a list of values from which one will be chosen according to the pressure
     * @return a value from valList at pressure p
     */
    public double getValueFromPressureOnGrid(double lat, double lon, double pressure, TreeMap<Double, float[]> valList) {
        // lon is transferred to 0...360 whereas 0...180 is identical to above
        if (lon < 0) {
            lon += 360.0;
        }
        int c = getIndex(lat, lon);

        // transfer data from vector to array (for interpolation)
        int nv = valList.size();
        double[] p = new double[nv];
        double[] val = new double[nv];

        // for each pressure level, retrieve the value from valList
        for (int i = 0; i < nv; ++i) {
            p[i] = (double) valList.keySet().stream().toArray()[i]; //pList.get(i);
            val[i] = valList.get((double) valList.keySet().stream().toArray()[i])[c];
        }
        // linear interpolation to get exact value

        return interpolate.linear(p, val, pressure, true);
    }

    public double getResolution(){
        return getDLat();
    }
    /**
     * The resolution of the underlying grid in latitude direction.
     *
     * @return Size of a grid cell in [deg]
     */
    public double getDLat() {
        return this.dLat;
    }

    public double getLa1() {
        return La1;
    }

    public double getLo1() {
        return Lo1;
    }

    public double getdLat() {
        return dLat;
    }

    public double getDLon() {
        return dLon;
    }

    public int getnLon() {
        return nLon;
    }

    public TreeMap<Double, float[]> getTemperature() {
        return mTemperature;
    }

    public TreeMap<Double, float[]> getTheta() {
        return mTheta;
    }

    public TreeMap<Double, float[]> getRelativeHumidity() {
        return mRelativeHumidity;
    }

    public TreeMap<Double, float[]> getWind_u() {
        return mWind_u;
    }

    public TreeMap<Double, float[]> getWind_v() {
        return mWind_v;
    }

    public TreeMap<Double, float[]> getGeoPotentialHeight() {
        return mGeoPotentialHeight;
    }

    public TreeMap<Double, float[]> getWindSpeed() {
        return mWindSpeed;
    }

    public TreeMap<Double, float[]> getWindHeading() {
        return mWindHeading;
    }

    public TreeMap<Double, float[]> getVerticalVelocity() {
        return mvvel;
    }

    public Timestamp getTimestamp() {
        return mTimestamp;
    }

    public static class Timestamp {
        public final java.time.ZonedDateTime mIssued_dateTime;

        public int getCycle() {
            return mCycle;
        }

        public int getForecast() {
            return mForecast;
        }

        public final int mCycle;
        public final int mForecast;

        @Override
        public String toString() {
            return "Timestamp{" +
                    "Issued_dateTime=" + mIssued_dateTime +
                    ", Prediction_dateTime=" + getPrediction_dateTime() +
                    ", Cycle=" + mCycle +
                    ", Forecast=" + mForecast +
                    '}';
        }

        public Timestamp(java.time.ZonedDateTime issued_dateTime, int cycle, int forecast) {

            mIssued_dateTime = issued_dateTime;
            mCycle = cycle;
            mForecast = forecast;
        }

        public ZonedDateTime getIssued_dateTime() {
            return mIssued_dateTime;
        }

        /**
         * retrieve issued date and correspondign forecast time
         *
         * @return
         */
        public ZonedDateTime getPrediction_dateTime() {
            return mIssued_dateTime.plusHours(mForecast);
        }
    }
}
