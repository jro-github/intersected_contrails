package ContrailsWakeVortex2020TEST;

import org.ini4j.Ini;
import java.io.*;
import java.lang.Math;


//public class ACTest {



//}

public class test {

    public contrail_parameters cp;
    public String name = "Hallo";



    public static boolean testxyz(){
    	boolean intersected=false;
    	if ((1<2+3)&&(2<3+4)) {
    		intersected = true;
    		System.out.println("1");
		}else{
			intersected = false;
			System.out.println("2");
		}
		
		return intersected;
    }

    public double theta(){
        this.cp.Theta=888;
        return(this.cp.Theta);
    }



    public static void main (String args[]) {

        //ACTest aircraft = new ACTest();
        //aircraft.Aircraft();



//        int nm_to_metre = 1852;
//        double altitude = 3000;
//        //WeatherFromGrib gribdata = new WeatherFromGrib("/Users/marco/Desktop/testordner/gribtest.grb2");
//        WeatherFromGrib gribdata = new WeatherFromGrib("/Users/marco/Desktop/testordner/grib_01_2016-02-07_1200.grib2");
//        //WeatherFromGrib gribdata = new WeatherFromGrib("/Users/marco/Downloads/gfs.t06z.pgrb2full.0p50.f000");
//        //WeatherFromGrib gribdata = new WeatherFromGrib("/Users/marco/Desktop/testordner/gribtesttest");
//        //double dlat = gribdata.getDLat();
//        //double dlon = gribdata.getDLon();
//        double latitude = 37;
//        double longitude = 49;
//        //System.out.println(altitude);
//        //System.out.println(gribdata.getWindUFromAltitude(latitude,longitude,altitude));
//        System.out.println(gribdata.getVvelFromAltitude(latitude,longitude,altitude));
//        System.out.println(gribdata.getVvelFromPressure(latitude,longitude,100000));
//        System.out.println(gribdata.getTempFromPressure(latitude,longitude,100000));
//        System.out.println("\n");
//        System.out.println(gribdata.getWindUFromAltitude(latitude,longitude,altitude));
//        System.out.println(gribdata.getWindVFromAltitude(latitude,longitude,altitude));
//        System.out.println("\n");
//        System.out.println(gribdata.getVvelFromAltitude(latitude,longitude,altitude));
//        System.out.println(gribdata.getTempFromAltitude(latitude,longitude,altitude));
//        System.out.println(gribdata.getPressureFromAltitude(latitude,longitude,altitude));
//        System.out.println(-(gribdata.getVvelFromAltitude(latitude,longitude,altitude))*287.085*(gribdata.getTempFromAltitude(latitude,longitude,altitude))/(9.81*100*(gribdata.getPressureFromAltitude(latitude,longitude,altitude))));
//        test var123 = new test();
//        System.out.println(var123.name);
//        System.out.println(var123.testxyz());
//        System.out.println("TEST");
//
//        System.out.println(gribdata.getVvelFromAltitude(54, 14,10500));
//        System.out.println(Math.log(gribdata.getVvelFromAltitude(54, 14,10500)));
//        
    	boolean A = testxyz();
    	System.out.println("TEST"+ A);


    }

        /*

        // Wind Components u and v
        //double wind_u0 = gribdata.getWindUFromAltitude(50,longitude,altitude);
        double wind_v0 = gribdata.getWindVFromAltitude(51,longitude,altitude);
        double wind_u1 = gribdata.getWindUFromAltitude(latitude+dlat,longitude,altitude);
        double wind_v1= gribdata.getWindVFromAltitude(latitude+dlat,longitude+dlon,altitude);

        //double wind_u =


        double lat0 = 50;
        double lat1 = 0;
        double lon0 = 13;
        double lon1;
        double test;
        double delta_t = 1; //Zeitschritt
        double wind_u0 = 55560;// in m/s
        */

        //System.out.println(Math.cos(Math.toRadians(60)));
        //double testest = gribdata.getThetaFromAltitude(50,longitude,altitude);
        //System.out.println(testest);
        //System.out.println(gribdata.getThetaFromAltitude(50,longitude,altitude));
        //System.out.println(gribdata.getHumidityFromAltitude(lat0,lon0,altitude));
/*

        // Weather at Position lat0, lon0

        temperature = gribdata.getTempFromAltitude(lat0,lon0,altitude);
        pressure = gribdata.getPressureFromAltitude(lat0,lon0,altitude);
        humidity = gribdata.getHumidityFromAltitude(lat0,lon0,altitude); // Relative Humidity
        density = pressure/(altitude*9.81);
        gribdata.get(lat0,lon0,altitude);


        // Displaced Position (lat1 [deg], lon1 [deg]) after Timestep delta_t

        lon1 = lon0 + (wind_u0/(delta_t*60*1852*Math.cos(Math.toRadians(lat0)))); // 60nm per degree at lat = 0° and 1852 metres per nm
        //System.out.println(lon1);
        lat1 = lat0 + (wind_v0/(delta_t*60*1852));




        lon1 = lon0 + (gribdata.getWindUFromAltitude(lat0,lon0,altitude)/(delta_t*60*1852*Math.cos(Math.toRadians(lat0)))); // 60nm per degree at lat = 0° and 1852 metres per nm
        //System.out.println(lon1);
        lat1 = lat0 + (gribdata.getWindVFromAltitude(lat0,lon0,altitude)/(delta_t*60*1852));

*/


/*
        System.out.println(gribdata.getWindUFromAltitude(latitude,longitude,altitude));         //get Wind u from Grib data in m/s
        System.out.println(gribdata.getWindUFromAltitude(51,longitude,altitude));
        //System.out.println(gribdata.getWindUFromAltitude(51,longitude,altitude));

        System.out.println(0.75*wind_u0+0.25*wind_u1);
        System.out.println(gribdata.getWindUFromAltitude(50.25,longitude,altitude));


        System.out.println(gribdata.getWindUFromAltitude(latitude+dlat,longitude+dlon,altitude));
        System.out.println(gribdata.getWindUFromAltitude(51,longitude,altitude));
        System.out.println(gribdata.getWindVFromAltitude(latitude,longitude,altitude));         //get Wind v from Grib data
        System.out.println(gribdata.getTempFromAltitude(latitude,longitude,altitude)-273.15);  //get Temperature from Grib data
        System.out.println(gribdata.getPressureFromAltitude(latitude,longitude,altitude));      //get Pressure from Grib data
        System.out.println(gribdata.getHumidityFromAltitude(latitude,longitude,altitude));      //get Humidity from Grib data
        System.out.println(gribdata.getDLat());
        System.out.println(gribdata.getdLat());

        System.out.println(gribdata.getTempFromAltitude(latitude,longitude,altitude)-273.15);
        */
    }
