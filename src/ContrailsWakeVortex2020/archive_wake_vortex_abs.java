package ContrailsWakeVortex2020;

//import com.sun.xml.internal.fastinfoset.tools.FI_SAX_Or_XML_SAX_SAXEvent;

import java.lang.Math;

// for file handling
import java.io.FileWriter;
import java.io.FileReader;
import java.io.BufferedReader;
import java.io.IOException;

//**********************************************************
//This is the original wake_vortex file copy for tests.
//**********************************************************

public class archive_wake_vortex_abs
{
    private FileWriter fwDataOut;
    private radiosonde rsonde;
    private contrail_parameters cp;
    private String fn_output;
    private WeatherFromGrib gribdata;

    private Extinktionseffizienz Ext;
    private double lterr;
    private int bsol;

    private double altitude; // [m]
    private double TAS; // [m/s]
    private double dt_fp; // [s]
    private double distance; //[m]
    private double M; // plane mass [kg]
    private double fuelflow; // [kg/s]

    //private double eps; // [?]
    private double dt_init; // initial time step

    //private double lat, lon;

    // Specific Gas Constants
    final Double Rs_air = 287.058;
    final Double Rs_water = 461.523;


    public archive_wake_vortex_abs(String ConfigFileName){
        // read all possible parameters from config-file and initialise variables


        // specify Grib File for weather data

        //gribdata = new WeatherFromGrib("/Users/marco/Desktop/testordner/gribtesttest");
        //gribdata = new WeatherFromGrib("/Users/marco/Desktop/testordner/grib_01_2016-02-07_1200.grib2");


        // used for calculating extinction and scattering efficiencies as well as assymetrie
        Ext = new Extinktionseffizienz();
        lterr = 10.; // Terrestrische Band

        // micrometer 3.08 ... 99.99
        bsol = 0; // Solares Band
        //  band ... das gewuenschte spektrale Band
        //     0 = 0.55 micrometer
        //     1 = 1.35 micrometer
        //     2 = 2.25 micrometer
        //     3 = 3.0125 micrometer
        //     4 = 3.775 micrometer
        //     5 = 4.5 micrometer

        // open config-file
        BufferedReader file= null;
        try {
            file=new BufferedReader(new FileReader(ConfigFileName));
            System.out.println("Datei "+ConfigFileName+" gefunden" );
        } catch(IOException e) {System.err.println("Datei "+ConfigFileName+" nicht gefunden" );};

        // read parameters from config file
        String eineZeile; // zum Zeilenweise einlesen der Datei
        try {
            this.cp=new contrail_parameters();
            // ignore header
            eineZeile = file.readLine();
            // Read filename of Radiosondenaufstieg and initialise rsonde
            this.rsonde = new radiosonde(file.readLine());
            // Read output filename
            eineZeile = file.readLine();
            this.fn_output = file.readLine();
            // Read flight level [m]
            eineZeile = file.readLine();
            this.cp.altitude = Double.valueOf(file.readLine());
            // Read aircraft speed [m/s]
            eineZeile = file.readLine();
            this.TAS = Double.valueOf(file.readLine());
            System.out.println("TAS "+this.TAS);
            // Read time step of flight profile [s]
            eineZeile = file.readLine();
            this.dt_fp = Double.valueOf(file.readLine());
            System.out.println("dt_fp"+this.dt_fp);
            // Calculate distance [m]
            this.distance = this.TAS*this.dt_fp;
            System.out.println("distance "+this.distance);
            // Read plane mass [kg]
            eineZeile = file.readLine();
            this.M = Double.valueOf(file.readLine());
            System.out.println("plane mass"+this.M);
            // Read fuelflow [kg/s]
            eineZeile = file.readLine();
            this.fuelflow = Double.valueOf(file.readLine());
            System.out.println("fuel flow" + this.fuelflow);
            // Read epsilon turbulenz [?]
            eineZeile = file.readLine();
            this.cp.eps = Double.valueOf(file.readLine());
            System.out.println("eps " + this.cp.eps);
            // Read vertical Diffusion contant [m2/s]
            eineZeile = file.readLine();
            this.cp.Dv = Double.valueOf(file.readLine());
            System.out.println("Dv (nicht verwendet)"+ this.cp.Dv);
            // Read horizontal Diffusion contant [m2/s]
            eineZeile = file.readLine();
            this.cp.Dh = Double.valueOf(file.readLine());
            System.out.println("Dh "+ this.cp.Dh);
            // Read shearing Diffusion contant [m2/s]
            eineZeile = file.readLine();
            this.cp.Ds = Double.valueOf(file.readLine());
            System.out.println("Ds "+ this.cp.Ds);
            // Read shearing factor [1/s]
            eineZeile = file.readLine();
            this.cp.s = Double.valueOf(file.readLine());
            System.out.println("s "+ this.cp.s);
            // Read mean wind velocity [m/s]
            eineZeile = file.readLine();
            this.cp.vo = Double.valueOf(file.readLine());
            System.out.println("vo " + this.cp.vo);
            // Read vertical wind velocity [m/s]
            eineZeile = file.readLine();
            this.cp.vz = Double.valueOf(file.readLine());
            System.out.println("vz " + this.cp.vz);
            // Read horizontal wind fluctuation [m/s]
            eineZeile = file.readLine();
            this.cp.uprime = Double.valueOf(file.readLine());
            System.out.println("uprime " + this.cp.uprime);
            System.out.println("T "+this.cp.T+", RH "+this.cp.RH);
            // Read initial time step [s]
            eineZeile = file.readLine();
            this.dt_init = Double.valueOf(file.readLine());
            System.out.println("dt " + this.dt_init);
            // Read initial latitude [deg]
            eineZeile = file.readLine();
            this.cp.lat = Double.valueOf(file.readLine());
            System.out.println("lat " + this.cp.lat);
            // Read initial longitude [deg]
            eineZeile = file.readLine();
            this.cp.lon = Double.valueOf(file.readLine());
            System.out.println("lon " + this.cp.lon);
            eineZeile = file.readLine();
            // Read filename of Gribfile and initialise gribdata
            String gribfilename = file.readLine();
            this.gribdata = new WeatherFromGrib(gribfilename);

            file.close(); // Close config-file
        } catch (IOException ex) {System.err.println("Fehler beim Einlesen");};

        /*
        double lat = this.cp.lat;
        double lon = this.cp.lon;
        double altitude = this.cp.altitude;
        */

        // Read Weather Data from Gribfile

        this.cp.vx = gribdata.getWindUFromAltitude(this.cp.lat,this.cp.lon,this.cp.altitude); // U-Component of wind
        this.cp.vy = gribdata.getWindVFromAltitude(this.cp.lat,this.cp.lon,this.cp.altitude); // V-Component of wind
        this.cp.Tflight=gribdata.getTempFromAltitude(this.cp.lat,this.cp.lon,this.cp.altitude); //K
        this.cp.T=this.cp.Tflight; //K
        this.cp.RH=gribdata.getHumidityFromAltitude(this.cp.lat,this.cp.lon,this.cp.altitude); //%

        //this.cp.density = weatherdata_helper(gribdata, this.cp, lat0, lon0, altitude)
        this.cp.density = density(); //FIXME
        this.cp.vz = vz();
        this.cp.Theta = theta();
        this.cp.thetadz = dthetadz();
        this.cp.eps = eps_calculator();
        //this.cp.eps = eps_calculator(this.cp.lat, this.cp.lon, this.cp.altitude);


        // calculate initial state of Contrail with GRIB
        initial_state istate = this.Contrail_depth(this.M, this.TAS, gribdata.getPressureFromAltitude(this.cp.lat, this.cp.lon, this.cp.altitude),
                this.cp.Tflight, this.cp.Theta, this.cp.thetadz, this.cp.altitude, this.cp.eps);
        this.cp.sigma0v=(float)(istate.Height/2.2); //50m
        this.cp.sigma0h=(float)(istate.Width/2.2); //20m


        // Calculate horizontal Diffusion coefficient [m2/s]
        // According to Schuhmann, J. Geophys. Res. 100, 14147(1995)
        // this.cp.Dh = 0.1 * this.cp.uprime * this.cp.sigma0h; // c_h * u' * sigmah
        //System.out.println("Dh (errechnet) "+ this.cp.Dh);
        // Read vertical Diffusion coefficient [m2/s]
        this.cp.Dv = this.Calc_Dv(this.cp.altitude);
        System.out.println("Dv (errechnet) "+ this.cp.Dv);

        //eingefügt
        this.cp.Dh = this.Calc_Dh(this.cp.eps);
        System.out.println("Dh (interpolert) "+ this.cp.Dh);

        System.out.println("Constructor has finished \n\n");
        return;
    }




    public double density (){

        /***

         Calculating air density depending on Relative Humidity


         // From Common Sources

         density = p * M / (R * T)

         p ... Pressure
         M ... molar Mass
         R ... Universal Gas Constant 8 J / (kg * K)
         T ... Temperature

         Rs = R / M

         density = p / (Rs * T)

         Rs_air ... 287,058 J / (kg * K)
         Rs_water ... 461.523 J / (kg * K)
         pd ... Saturation vapour pressure
         p  ... Air Pressure
         RH ... Relative Humidity

         // Relative Values for air and water to get Specific Gas Constant for humid air
         // as derived from Dalton´s Law and also shown on: de.wikipedia.org/wiki/Luftdichte

         Rs = Rs_air / (1 - RF * pd/p * (1 - Rs_air / Rs_water))


         // Saturation vapour pressure at air temperature T from WMO,
         // Source: WMO GUIDE TO METEOROLOGICAL INSTRUMENTS AND METHODS OF OBSERVATION
         // WMO-No. 8 (2014 edition, Updated in 2017), page 346 )
         // https://library.wmo.int/doc_num.php?explnum_id=10179

         pd(T) = 6.107 * 10 ^ (7.5*T/(237.3+T)) ... T has to be in °C

         OR after SONTAG

         */

        // sat_vap_pres after WMO GUIDE TO METEOROLOGICAL INSTRUMENTS AND METHODS OF OBSERVATION WMO-No. 8 (2014 edition, Updated in 2017), page 346 ) https://library.wmo.int/doc_num.php?explnum_id=10179
        // NOT IN USE, replaced by SONTAG
        // this.cp.sat_vap_pres = 6.107 * Math.pow(10,(7.5*(this.cp.T-273.15)/(237.3+(this.cp.T-273.15)))); //with Temp Conversion to deg Celsius

        // sat_vap_pres after: SONTAG
        this.cp.sat_vap_pres = Math.exp(-6096.9385/cp.T+16.635794-2.711193E-2*cp.T+1.673952E-5*cp.T*cp.T+2.433502*Math.log(cp.T));
        this.cp.Rs_humidair = Rs_air / 1 - (this.cp.RH * (this.cp.sat_vap_pres/gribdata.getPressureFromAltitude(this.cp.lat,this.cp.lon,this.cp.altitude) * (1 - Rs_air / Rs_water)));

        this.cp.density = 100 * gribdata.getPressureFromAltitude(this.cp.lat,this.cp.lon,this.cp.altitude) / (this.cp.Rs_humidair * this.cp.T); //gribdata pressure value is in hPa

        return (this.cp.density);

    }


    public double vz (){

        /***

         Calculating Vertical Wind vz

         // Vertical velocity VVEL is the speed at which the air is rising or sinking.
         // Measured in negative microbars per second as a distance.
         // Instead of meters the distance measured is pressure, so a negative distance in pressure
         // is actually a positive distance in altitude, hence a negative microbar distance is used here.

         // p V = m Rs_humidair T ==> rho = p / Rs_humidair T
         // vvel = omega ... vz = -vvel Rs_humidair T / p g

         */

        //unter Einberechnung der Feuchtigkeit
        this.cp.sat_vap_pres = Math.exp(-6096.9385/cp.T+16.635794-2.711193E-2*cp.T+1.673952E-5*cp.T*cp.T+2.433502*Math.log(cp.T)); // in hPa
        this.cp.Rs_humidair = Rs_air / 1 - (this.cp.RH * (this.cp.sat_vap_pres/gribdata.getPressureFromAltitude(this.cp.lat,this.cp.lon,this.cp.altitude) * (1 - Rs_air / Rs_water))); // pressure in hPa
        this.cp.vz = -(gribdata.getVvelFromAltitude(this.cp.lat,this.cp.lon,this.cp.altitude))*this.cp.Rs_humidair*(gribdata.getTempFromAltitude(this.cp.lat,this.cp.lon,this.cp.altitude))/(9.81*100*(gribdata.getPressureFromAltitude(this.cp.lat,this.cp.lon,this.cp.altitude)));

        return this.cp.vz;


    }



    public double theta (){


        // Source: Approximation according to D. Bolton, 1980 as shown on:
        // http://www.uni-koeln.de/math-nat-fak/geomet/meteo/winfos/radiosonden/Europa/radiosondengrafiken.pdf

        Double mixture = (Rs_air / Rs_water) * (this.cp.sat_vap_pres/((gribdata.getPressureFromAltitude(this.cp.lat,this.cp.lon,this.cp.altitude))-this.cp.sat_vap_pres));
        this.cp.Theta = this.cp.T * Math.pow((1000/(gribdata.getPressureFromAltitude(this.cp.lat,this.cp.lon,this.cp.altitude))),0.2854*(1-0.28*mixture));

        return this.cp.Theta;

    }


    public double dthetadz (){

        // calculation of theta dz
        int x = 1000; // step width for interpolating theta dz in metres
        Double mixture = (Rs_air / Rs_water) * (this.cp.sat_vap_pres/((gribdata.getPressureFromAltitude(this.cp.lat,this.cp.lon,this.cp.altitude))-this.cp.sat_vap_pres));

        if (altitude-x < 0) {
            this.cp.thetadz = this.cp.T * Math.pow((1000/(gribdata.getPressureFromAltitude(this.cp.lat,this.cp.lon,this.cp.altitude+x))),0.2854*(1-0.28*mixture)) - (this.cp.T * Math.pow((1000/(gribdata.getPressureFromAltitude(this.cp.lat,this.cp.lon,this.cp.altitude))),0.2854*(1-0.28*mixture)));
        }
        else {

            this.cp.thetadz = this.cp.T * Math.pow((1000/(gribdata.getPressureFromAltitude(this.cp.lat,this.cp.lon,this.cp.altitude+x))),0.2854*(1-0.28*mixture)) - (this.cp.T * Math.pow((1000/(gribdata.getPressureFromAltitude(this.cp.lat,this.cp.lon,this.cp.altitude-x))),0.2854*(1-0.28*mixture)));
        }
        return this.cp.thetadz;

    }


    public double eps_calculator() {//(double lat, double lon, double altitude) {


        /**

         Calculation of Eddy Dissipation Rate
         after formula eps = Math.exp (C1 + C2 * (lnD - lnDe) / SD)
         as shown in: SHARMAN AND PEARSON
         "Prediction of Energy Dissipation Rates for Aviation Turbulence.
         Part I: Forecasting Nonconvective Turbulence" (2016)

         Average Values above altitude 0m for

         C1 = -2.572
         C2 = 0.5067

         **/

        double lat = this.cp.lat;
        double lon = this.cp.lon;
        double altitude = this.cp.altitude;


            // defining the full degrees range around current lat/lon which shall be considered determining the eps
            int x = 15; // 15 degrees range was tested to be sufficient
            int y = x;

            if (lat < -90 || lat > 90 || lon < -180 || lon > 180) {
                throw new IllegalArgumentException("One or more values out of range (lat = " + lat + " ,lon = " + lon);
            }
            if (y < 0 || y > 90 || x < 0 || x > 180) {
                throw new IllegalArgumentException("One or more values out of range (full degrees range around current lat/lon to be considered... lateral = " + y + " ,longitudinal = " + x);
            }


            double current_lon = -777;
            double current_lat = -777;
            double current_vvel = 0;

            double SD_lnVVEL = 0;
            double Mean_lnVVEL = 0.; //Standard Deviation and Mean Value of logarithmic VVEL values between lat/lon start and end

            double Mean_lnVVELx = 0;
            double testsjkjksjs = 0;

            // longitudinal starting location is from lon - y and calculation is always performed in eastern direction until lon +y
            // lateral starting location is from lat - x and wents to final location lat + x (even if that means to go south first and then north or vice versa in case boundary values for lat = +- 90 deg are tangered

            // calculation of Mean_lnVVEL

            for (int j = 0; j < 2 * x; j++) {
                current_lon = lon - x + j;
                if (current_lon < -180) {
                    current_lon = current_lon + 360;
                } else if (current_lon > 180) {
                    current_lon = current_lon - 360;
                }
                for (int i = 0; i < 2 * y; i++) {
                    current_lat = lat - y + i;
                    if (current_lat < -90) {
                        current_lat = -(current_lat + 180);
                    } else if (current_lat > 90) {
                        current_lat = 180 - current_lat;
                    }

                    // handles rare error that interpolated vvel value is 0
                    current_vvel = gribdata.getVvelFromAltitude(current_lat, current_lon, altitude);
                    if (current_vvel == 0) {
                        current_vvel = current_vvel + 1e-6;
                    }
                    Mean_lnVVEL = Mean_lnVVEL + Math.log(Math.abs(current_vvel));

                }
            }

            Mean_lnVVEL = Mean_lnVVEL / ((2 * y + 1) * (2 * x + 1));



            // calculation of SD_lnVVEL

            for (int j = 0; j < 2 * x; j++) {
                current_lon = lon - x + j;
                if (current_lon < -180) {
                    current_lon = current_lon + 360;
                } else if (current_lon > 180) {
                    current_lon = current_lon - 360;
                }
                for (int i = 0; i < 2 * y; i++) {
                    current_lat = lat - y + i;
                    if (current_lat < -90) {
                        current_lat = -(current_lat + 180);
                    } else if (current_lat > 90) {
                        current_lat = 180 - current_lat;
                    }

                    // handles rare error that interpolated vvel value is 0
                    current_vvel = gribdata.getVvelFromAltitude(current_lat, current_lon, altitude);
                    if (current_vvel == 0) {
                        current_vvel = current_vvel + 1e-6;
                    }
                    SD_lnVVEL = SD_lnVVEL + Math.pow(Math.log(Math.abs(current_vvel)) - Mean_lnVVEL, 2);

                }
            }

            SD_lnVVEL = Math.pow(SD_lnVVEL / ((2 * y + 1) * (2 * x + 1) - 1), 0.5);

            Double C1 = -2.572;
            Double C2 = 0.5067;
            this.cp.eps = Math.exp(C1 + C2 * (Math.log(Math.abs(gribdata.getVvelFromAltitude(lat, lon, altitude))) - Mean_lnVVEL) / SD_lnVVEL);

            return this.cp.eps;
    }




    public double Calc_Dh(double eps) {

        /*
           According to Tab 3.2, Rosenow, "Optical Properties of Contrails", pg. 68
           Eddy dissipation rate ε [m2s−3] ... column 0
           Horizontal diffusivity Dh [m2s−2] ... column 1

           Interpolation of values
           No extrapolation beyond lower and upper boundaries!

           Returns interpolated value for Dh
       */

        double[][] eps_Dh_Table = { {Math.pow(10,-6),5.00} , {5*Math.pow(10,-6),8.75} , {Math.pow(10,-5),12.50} , {5*Math.pow(10,-5),16.254} , {Math.pow(10,-4),20.00} };
        double returnvalue = 0;

        if      (eps < eps_Dh_Table [0][0]) {
            returnvalue = eps_Dh_Table [0][1]; // lower boundary
        }
        else if (eps >= eps_Dh_Table [(eps_Dh_Table.length)-1][0]) {
            returnvalue = eps_Dh_Table[(eps_Dh_Table.length)-1][1]; // upper boundary
        }
        else {
            for (int i=1; i<=((eps_Dh_Table.length)-1); i++) {
                if (eps >= eps_Dh_Table [i-1][0] && eps < eps_Dh_Table [i][0]) {
                    returnvalue = ((eps - eps_Dh_Table [i-1][0]) / (eps_Dh_Table [i][0] - eps_Dh_Table [i-1][0]) * (eps_Dh_Table [i][1] - eps_Dh_Table [i-1][1]) + eps_Dh_Table [i-1][1]);
                }
            }

        }

     return returnvalue;

    }




    public double Calc_Dv(double height)
    {

        double Nsquare=9.81/theta()*dthetadz(); //square of Brunt Vaisala frequency [a.u. 1/s^2]
        return 0.1*this.cp.eps/Nsquare;


    }


    public initial_state Contrail_depth(double M, double TAS, double p, double T, double Theta, double dThetadz, double z, double eps)
    {
        // M...aircraft mass [kg]
        // TAS... True airspeed [m/s]
        // p...pressure [hPa]
        // T...Temperature [K]
        // Theta...potential temperature [K]
        // dThetadz...derivation potential temperature with height [K/m]
        // z...Altitude [m]
        // eps...Eddy dissipation rate [m2/s3]
        //
        // calculation of initial contrail shape with the help of wake vortex characteristics
        //

        // aircraft parameters of A320
        double rc=1.1935; // [m] core radius
        double B= 34.1; //[m] span width
        double bo=Math.PI/4.*B; //[m] initial vortex spacing

        // to be calculated in this function
        double r=10; //vortex radius at the beginning of the dispersion regime [m]
        double h=20; //contrail height [m]
        double w=30; //contrail width [m]

        // characteristic scales
        // bo ... initial vortex spacing
        // wo ... [m/s] initial tangential velocity
        // Gammao ... [m2/s] initial vortex circulation
        // tprime...[s] characteristic time scale
        double rho=p*100/(287.058*T); //density of air
        //FIXME
        //double rho=density();
        //System.out.println("rho "+rho);
        double Gammao=M*9.81/(rho*bo*TAS);
        double wo=Gammao/(2.*Math.PI*bo);
        double tprime=2*Math.PI*Math.pow(bo,2.)/Gammao;

        // calculation of initial vortex radius
        // with the help of Dispersion model (tangential velocity) according to Holzaepfel 2003
        // r is the distance, where w=1/e^(1/2) * w_max
        r=(Math.exp(0.5)+Math.pow(Math.E-1,0.5))*rc;

        // calculation of contrail width
        w=2*r+bo;

        //calculation of contrail height
        double epsstar=Math.pow(eps*bo,1/3.)/wo;

        //System.out.println("epsstar "+epsstar);
        double Nstar=Math.sqrt(9.81/Theta*dThetadz)*tprime; //Brunt Vaisala frequency [a.u. 1/s]
        //System.out.println("Nstar "+Nstar);
        double Tstar=1;

        double T20star=5;
        if(epsstar>0.2535)
        {
            T20star=0.804*Math.pow(epsstar,-0.75)-1;
        }
        else if(epsstar>0.0235)
        {
            int i = 0;
            double f=1;
            double df=1;
            double dx=1;
            Tstar=1; // starting value for Newton-iteration
            do
            {
                f=Math.pow(Tstar,0.25)*Math.exp(-0.7*Tstar)-epsstar;
                df=(0.25*Math.pow(Tstar,-0.75)-0.7*Math.pow(Tstar,0.25))*Math.exp(-0.7*Tstar);
                dx=f/df;
                Tstar=Tstar-dx;
                i=i+1;
            } while(Math.abs(dx)>1E-3 && i<11);
            if (i>=11) System.out.println("Iteration not long enough for Tstar");
            T20star=Tstar-1.;
        };

        //System.out.println("T20star "+T20star);
        double T2star=T20star*Math.exp(-0.185*T20star*Nstar);

        double wstar=1.-Math.exp(-1.257*Math.pow(0.4*bo/rc,2));

        h=wstar*T2star*bo+2*r;
        System.out.println("T2star "+T2star + " wstar "+ wstar);
        System.out.println("T2 "+(T2star*tprime));
        System.out.println("r " + r + ", h " + h + ", w " + w);
        System.out.println("w* = " + wstar +", w0 = "+ wo +", w = "+(wstar*wo));
        return new initial_state(r, h, w);
    }

    public double ContrailCrossSection(double t, contrail_parameters cp)
    {
        // calculates area of contrail cross-section
        // t...time [s]
        // Dv, Dh, Ds... vertical, horizontal and sheared diffusion [m2/s]
        // vo...mean value of wind velocity at flight level [m/s]
        // s....constant wind shear 1/s

        double sigmav=2*cp.Dv*t+Math.pow(cp.sigma0v,2);
        cp.s = (gribdata.getWindspeedFromAltitude(this.cp.lat,this.cp.lon,this.cp.altitude+sigmav)-gribdata.getWindspeedFromAltitude(this.cp.lat,this.cp.lon,this.cp.altitude-sigmav))/(2*sigmav); //calculates the change of windspeed between upper and lower boundary of contrail per metre
        double sigmah=2./3*cp.s*cp.s*cp.Dv*Math.pow(t,3)+(2*cp.Ds+cp.s*Math.pow(cp.sigma0v,2))*cp.s*Math.pow(t,2)+2*cp.Dh*t+Math.pow(cp.sigma0h,2);
        double sigmas=cp.s*cp.Dv*Math.pow(t,2)+(2*cp.Ds+cp.s*Math.pow(cp.sigma0v,2))*t;

        // System.out.println("Sigmas "+sigmas);

        double detsigma=sigmah*sigmav-Math.pow(sigmas,2);
        double A=2.*Math.PI*Math.sqrt(detsigma);  //Schumann 1995

        return A;
    }

    public double AdiabaticCurve(contrail_parameters cp, double z)
    {
        //SATURATION PRESSUREs FROM SONNTAG, Temp IN K, PSAT IN hPa
        //  double eStar=Math.exp(-6096.9385/cp.T+16.635794-2.711193E-2*cp.T+1.673952E-5*cp.T*cp.T+2.433502*Math.log(cp.T));
        double eStarIce=Math.exp(-6024.5282/cp.Tflight+ 24.721994+1.0613868E-2*cp.Tflight- 1.3198825E-5*cp.Tflight*cp.Tflight- 0.49382577*Math.log(cp.Tflight));
        double deStarIce=(6024.5282/Math.pow(cp.Tflight,2)-1.0613868E-2-2*1.3198825E-5*cp.Tflight- 0.49382577/cp.Tflight)*eStarIce;

        // Saettigungsmischungsverhaeltnis
        double pressure=gribdata.getPressureFromAltitude(this.cp.lat,this.cp.lon,this.cp.altitude);
        double mStar=0.622*eStarIce/(pressure-eStarIce); // kg/kg
        // LatentHeat = Verdampfungswaerme
        double L= 2.5e6; // J/kg
        // spezifische Waerme bei konstantem Druck
        double c_p=1004; // J/(kg K)
        // spezifische Gaskonstante von Luft
        double RL=287; // J/(kg K)
        // Adiabatic Curve
        double Gammas=9.81/c_p*(1+L*mStar/(RL*cp.Tflight))/(1+L*mStar*deStarIce/(c_p*eStarIce));

        return Gammas; // [K/m]
    }

    public double kinvisc_Suth(double p, double T)
    {
        // p...pressure[hPa], T...Temperature[K],
        // calculates kinematic viscosity [m2/s] of air according to
        // Philosophical Magazine 5 (1893) page 507 W. Sutherland

        //Sutherland parameter
        double T0=273; //in K
        double mu0=1.716e-5; //in (Ns)/m2
        double C=111; // in K

        //calculation of density
        double rho=1;  //density of air
        rho=p*100/(287.058*T);

        //calculation of viscosity
        double mu=1; //dynamic viscosity
        mu=mu0*(T0+C)/(T+C)*Math.pow(T/T0, 1.5); //Sutherland formular
        return mu/rho; //kinematic viscosity
    }

    public double dynvisc_Suth(double T)
    {
        // p...pressure[hPa], T...Temperature[K],
        // calculates dynamic viscosity [kg/m s] of air according to
        // Philosophical Magazine 5 (1893) page 507 W. Sutherland

        //Sutherland parameter
        double T0=273; //in K
        double mu0=1.716e-5; //in (Ns)/m2
        double C=111; // in K

        //calculation of density
        //double rho=1;  //density of air
        //rho=p*100/(287.058*T);

        //calculation of viscosity
        double mu=1; //dynamic viscosity
        mu=mu0*(T0+C)/(T+C)*Math.pow(T/T0, 1.5); //Sutherland formular
        return mu; //dynamic viscosity
    }

    public double SedimentationSpeed(double dice, contrail_parameters cp)
    {
        // Stokessches Sedimentations-Gesetz
        double vsed=2*Math.pow(dice/2.,2)*917*9.81/(9*dynvisc_Suth(cp.T));
        return vsed;
    }

    public double lambda_e(double sigma_h, double sigma_v, double sigma_s, double r, double Nice)
    {
        double Qabs = 4.e-5;
        double Qsca = 2.;
        double le=0.;
        double detsigma = Math.pow(sigma_h*sigma_v,2)-Math.pow(sigma_s,2);
        double n = Nice/(this.distance*2*Math.PI*Math.pow(detsigma,0.5));
        le = 1./((Qabs+Qsca)*n*Math.PI*r*r);
        return le;
    }

    public double Re(double r, double vsed, double vkin)
    {
        // Berechnung der Reynolds-Zahl [einheitenlos]
        //
        // r... Partikelradius [m]
        // vsed... Sedimentationsgeschwindigkeit [m/s]
        // vkin... kinematische Viskositaet [m2/s]
        return 2*r*vsed/vkin;
    }




    public void ContrailEvolution()
    {
        contrail_parameters cp = this.cp;
        // Simulation of Contrail Evolution
        // mH2O...emitted water mass [kg]
        // distance...distance from last flight profile step [m]
        // Parameters which describe the contrail shape:
        // Dv, Dh, Ds... vertical, horizontal and sheared diffusion [m2/s]
        // vo...mean value of wind velocity at flight level [m/s]
        // s....constant wind shear

        // === 1st step: calculate initial number of ice particles ===
        // Sussmann, J. Geophys. Res., 106(D5), 4899 (2001)
        double Nice=1e15*this.fuelflow*this.distance/this.TAS; // number
        System.out.println("number of ice particles: "+Nice);

        // ========== Beginning of contrail lifetime ==========
        double t=0; // at t=0
        double z=0; // at flight level
        // FIXME
        cp.T=cp.Tflight; // contrail temperature is equal to ambient temperature at flight level
        // delta t for time integration steps
        double delta_t=this.dt_init; // [s]

        // === calculate ice water content ===
        double CCS=ContrailCrossSection(0, cp);
        //FIXME
        //SATURATION PRESSUREs FROM SONNTAG, Temp IN K, PSAT IN PA - METEOROL. Z., 3 (1994) 51-66.
        double eStar=Math.exp(-6096.9385/cp.T+16.635794-2.711193E-2*cp.T+1.673952E-5*cp.T*cp.T+2.433502*Math.log(cp.T));
        double eStarIce=Math.exp(-6024.5282/cp.T+ 24.721994+1.0613868E-2*cp.T- 1.3198825E-5*cp.T*cp.T- 0.49382577*Math.log(cp.T));

        double RHice=gribdata.getHumidityFromAltitude(this.cp.lat,this.cp.lon,this.cp.altitude+z)*eStar/eStarIce;

        if (RHice < 100.) {
            System.out.println();
            System.out.println("RHice = "+RHice+" < 100 ... Relative Humidity of ambient air is only "+Math.round((gribdata.getHumidityFromAltitude(this.cp.lat,this.cp.lon,this.cp.altitude+z)))+"% ... no contrail evolves!");
            System.exit(1);
        }
        double delta_e=eStarIce*(RHice/100.-1.);
        double IWCs=0.21667*delta_e/cp.T;
        double IWC=1.24*this.fuelflow/(this.TAS*CCS)+IWCs;
        //double testpressure=gribdata.getPressureFromAltitude(this.cp.lat,this.cp.lon,10500);

        double IWC_old=0.;
        double t_old=0.;
        double slope=0.;
        double sigmah=cp.sigma0h;
        double sigmav=cp.sigma0v;
        double sigmas=cp.s*cp.Dv*Math.pow(t,2)+(2*cp.Ds+cp.s*Math.pow(cp.sigma0v,2))*t;

        double dice=Math.pow(IWC*CCS*this.distance/(Nice*917)*6./Math.PI,1./3.);
        double rext= dice/2.*1.e6;
        try {
            this.fwDataOut =new FileWriter(fn_output);
        } catch (IOException ex) {System.err.println("Fehler beim Anlegen der Datei");};

        try {
//            this.fwDataOut.write("t, z, sigmah, sigmav, sigmas, rice, Nice, g_t, Qabs_t, Qsca_t, g_s, Qabs_s, Qsca_s, CCS, IWC, RHice, Re, v_sed, windshear_s, wind_u, wind_v, wind_w, lat, lon, altitude\n");
//            rext= dice/2.*1.e6;
//            this.fwDataOut.write(t+", "+z+", "+sigmah+", "+sigmav+", "+sigmas+", "+(dice/2)+", "+Nice+", "+
//                    Ext.Calc_g_terr(rext,lterr)+", "+Ext.Calc_Qabs_terr(rext,lterr)+", "+Ext.Calc_Qsca_terr(rext,lterr)+", "+
//                    Ext.Calc_g_sol(rext,bsol)+", "+Ext.Calc_Qabs_sol(rext,bsol)+", "+Ext.Calc_Qsca_sol(rext,bsol)+ ", "+ CCS+", "+IWC+", "+RHice+", "+
//                    Re(dice/2.,SedimentationSpeed(dice, cp),kinvisc_Suth(gribdata.getPressureFromAltitude(this.cp.lat,this.cp.lon,this.cp.altitude),cp.T))+", "+
//                    SedimentationSpeed(dice, cp)+", "+this.cp.s+", "+gribdata.getWindUFromAltitude(this.cp.lat,this.cp.lon,this.cp.altitude)+", "+
//                    gribdata.getWindVFromAltitude(this.cp.lat,this.cp.lon,this.cp.altitude)+", "+this.cp.vz+", "+((int) (this.cp.lat*10000)/10000.)+", "+
//                    ((int) (this.cp.lon*10000)/10000.)+", "+Math.round(this.cp.altitude)+"\n");
            this.fwDataOut.write("t, z, T, RH, RHice, IWC, CCS, sigmah, sigmav, sigmas, Dh, Dv, Dh/Dv, v_sedi, rice, lambda_e\n");
            this.fwDataOut.write(t+", "+z+", "+cp.T+", "+this.rsonde.get_RH(this.altitude+z)+", "+RHice+", "+IWC+", " +CCS+", "+sigmah+", "+sigmav+", "+sigmas+", "+cp.Dh+", "+cp.Dv+", "+(cp.Dh/cp.Dv)+", "+SedimentationSpeed(dice, cp)+", "+(dice/2)+", "+lambda_e(sigmah, sigmav, sigmas, (dice/2.), Nice)+"\n");
        } catch (IOException ex) {System.err.println("Fehler beim Schreiben");};
        //System.out.println(t+", "+z+", "+cp.T+", "+sigmav+", "+sigmah+", "+dice+", "+IWC);
        do {
            // calculate ice particle diameter
            dice=Math.pow(IWC*CCS*this.distance/(Nice*917)*6/Math.PI,1/3.);
            //dice = Math.cbrt(IWC*CCS*this.distance/(Nice*917)*6/Math.PI);

            // next time step
            t_old=t;
            t=t+delta_t;
            System.out.println("Timestep t: "+t);
            //cp.Dht = cp.Dht + cp.Dh*t;
            // contrail geometry
            sigmah=Math.pow(2./3.*cp.s*cp.s*cp.Dv*Math.pow(t,3)+(2*cp.Ds+cp.s*Math.pow(cp.sigma0v,2))*cp.s*Math.pow(t,2)+2*cp.Dh*t+Math.pow(cp.sigma0h,2),0.5);
            sigmav=Math.pow(2.*cp.Dv*t+Math.pow(cp.sigma0v,2),0.5);
            sigmas=cp.s*cp.Dv*Math.pow(t,2)+(2*cp.Ds+cp.s*Math.pow(cp.sigma0v,2))*t;
            //cp.Dh = 0.1 * cp.uprime * sigmah;


            // eStarIce before adiabatic heating
            double eStarIce_old=Math.exp(-6024.5282/cp.T+ 24.721994+1.0613868E-2*cp.T- 1.3198825E-5*cp.T*cp.T- 0.49382577*Math.log(cp.T));


            // Contrail horizontal movement
            // Calculates displaced contrail position (lat, lon)
            //this.cp.lat = this.cp.lat + (-16.0*(delta_t/(60000*1.852)));
            this.cp.lat = this.cp.lat + (gribdata.getWindVFromAltitude(this.cp.lat,this.cp.lon,this.cp.altitude)*(delta_t/(60000*1.852)));
            this.cp.lon = this.cp.lon + (gribdata.getWindUFromAltitude(this.cp.lat,this.cp.lon,this.cp.altitude)*(delta_t/(60000*1.852*Math.cos(Math.toRadians(this.cp.lat))))); // 60nm per degree at lat = 0° and 1852 metres per nm


            // Contrail vertical movement
            // Let contrail drift down
            cp.vz = vz();
            double delta_z=(cp.vz-SedimentationSpeed(dice, cp))*delta_t; // negative delta_z means sinking not thinking
            z=z+delta_z; //z ist die Summe der delta z je Durchlauf
            this.cp.altitude = this.cp.altitude + delta_z;


            // Temperature
            cp.T=(float)(cp.Tflight-AdiabaticCurve(cp,z)*z); // adiabatic heating


            // calculate additional ice water content

            //SATURATION PRESSUREs FROM SONNTAG, Temp IN K, PSAT IN PA - METEOROL. Z., 3 (1994) 51-66.
            eStar=Math.exp(-6096.9385/cp.T+16.635794-2.711193E-2*cp.T+1.673952E-5*cp.T*cp.T+2.433502*Math.log(cp.T));
            eStarIce=Math.exp(-6024.5282/cp.T+ 24.721994+1.0613868E-2*cp.T- 1.3198825E-5*cp.T*cp.T- 0.49382577*Math.log(cp.T));
            // FIXME
            //RHice=gribdata.getHumidityFromAltitude(this.cp.lat,this.cp.lon,this.cp.altitude)*eStar/eStarIce;
            RHice=this.rsonde.get_RH(this.altitude+z)*eStar/eStarIce;
            System.out.println(t+","+z+","+RHice+","+IWC);
            if (RHice<100)
            {
                System.out.println("RHice = "+RHice+" < 100");
            } else {
                // important for next delta_t
                IWC_old=IWC;
                // get on with it
                delta_e=eStarIce*(RHice/100.-1.);
                IWCs=0.21667*delta_e/cp.T; // Masse des Wasser ueber Saettigung in der Umgebung (pro Volumen)
                // calculate new ice water content
                // before IWCs is added, IWC is reduced because of adiabatic heating
                IWC=IWC-(eStarIce-eStarIce_old)/(462*cp.T);
                double CCSn=ContrailCrossSection(t+delta_t,cp);
                delta_t =0.1*2*cp.Dh*cp.Dv*t;
                IWC=(IWCs*(CCSn-CCS)+IWC*CCS)/CCSn;
                CCS=CCSn;

                //delta_t=0.1/Math.abs(delta_z/delta_t);
    /*if (delta_t > 60.) {
     System.out.println(delta_t);
     delta_t=60.;
     };*/
                //System.out.println(t+", "+z+", "+cp.T+", "+sigmav+", "+sigmah+", "+dice+", "+IWC);
                try {
                    rext= dice/2.*1.e6;
//                    this.fwDataOut.write(t+", "+z+", "+sigmah+", "+sigmav+", "+sigmas+", "+(dice/2)+", "+Nice+", "+
//                            Ext.Calc_g_terr(rext,lterr)+", "+Ext.Calc_Qabs_terr(rext,lterr)+", "+Ext.Calc_Qsca_terr(rext,lterr)+", "+
//                            Ext.Calc_g_sol(rext,bsol)+", "+Ext.Calc_Qabs_sol(rext,bsol)+", "+Ext.Calc_Qsca_sol(rext,bsol)+", "+CCS+", "+IWC+", "+RHice+", "+Re(dice/2.,SedimentationSpeed(dice, cp),kinvisc_Suth(gribdata.getPressureFromAltitude(this.cp.lat,this.cp.lon,this.cp.altitude),cp.T))+", "+SedimentationSpeed(dice, cp)+", "+this.cp.s+", "+gribdata.getWindUFromAltitude(this.cp.lat,this.cp.lon,this.cp.altitude)+", "+gribdata.getWindVFromAltitude(this.cp.lat,this.cp.lon,this.cp.altitude)+", "+this.cp.vz+", "+((int) (this.cp.lat*10000)/10000.)+", "+((int) (this.cp.lon*10000)/10000.)+", "+Math.round(this.cp.altitude)+"\n");
                    this.fwDataOut.write(t+", "+z+", "+cp.T+", "+this.rsonde.get_RH(this.altitude+z)+", "+RHice+", "+IWC+", " +CCS+", "+sigmah+", "+sigmav+", "+sigmas+", "+cp.Dh+", "+cp.Dv+", "+(cp.Dh/cp.Dv)+", "+SedimentationSpeed(dice, cp)+", "+(dice/2)+", "+lambda_e(sigmah, sigmav, sigmas, (dice/2.), Nice)+"\n");
                } catch (IOException ex){System.err.println("Fehler beim Schreiben");};
            };
        } while ((IWC>10e-8) && (RHice>=100));
        if (RHice>=100) {
            System.out.println("IWC = "+IWC+" <= 10e-8");
        };

        try {
            this.fwDataOut.close();
        } catch (IOException ex){System.err.println("Fehler beim Schliessen der Datei");};

        return;
    }

    public static void main(String[] args)
    {
        //wake_vortex wv=new wake_vortex(args[0]);
        archive_wake_vortex_abs wv= new archive_wake_vortex_abs ("Config_NEU_03APR20_abs.txt");
        wv.ContrailEvolution();
    }
}