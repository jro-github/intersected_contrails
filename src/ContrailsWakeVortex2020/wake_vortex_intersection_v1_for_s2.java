package ContrailsWakeVortex2020;

//import com.sun.xml.internal.fastinfoset.tools.FI_SAX_Or_XML_SAX_SAXEvent;

import java.lang.Math;

//for file handling
import java.io.FileWriter;
import java.io.FileReader;
import java.io.BufferedReader;
import java.io.IOException;

public class wake_vortex_intersection_v1_for_s2 {

	////////////////////////////////////////////////////
	/* contrail_0 parameter definition of the old one */
	////////////////////////////////////////////////////
	private FileWriter fwDataOut_0;
	private radiosonde rsonde_0;
	private contrail_parameters cp_0;
	private String fn_output_0;
	private Extinktionseffizienz Ext_0;
	private double lterr_0;
	private int bsol_0;

	private double altitude_0; // [m]
	private double TAS_0; // [m/s]
	private double dt_fp_0; // [s]
	private double distance_0; // [m]
	private double M_0; // plane mass [kg]
	private double fuelflow_0; // [kg/s]

	// private double eps; // [?]
	private double dt_init_0; // initial time step

	////////////////////////////////////////////////////
	/* contrail_1 parameter definition of the new one */
	////////////////////////////////////////////////////
	private FileWriter fwDataOut_1;
	private radiosonde rsonde_1;
	private contrail_parameters cp_1;
	private String fn_output_1;
	private Extinktionseffizienz Ext_1;
	private double lterr_1;
	private int bsol_1;

	private double altitude_1; // [m]
	private double TAS_1; // [m/s]
	private double dt_fp_1; // [s]
	private double distance_1; // [m]
	private double M_1; // plane mass [kg]
	private double fuelflow_1; // [kg/s]

	private double dt_init_1; // initial time step

	// new parameters to define the time difference and intersection angle

	private contrail_parameters cp_inter;
	private FileWriter fwDataOut_inter;
	private Extinktionseffizienz Ext_inter;

	final double time_diff = 900; // time difference of contrail_1 after contrail_0, unit [s]
	double angle = 0.; // intersection angle, limited to [0, pi/2]

	///////////////////////////////////////////////////
	/////// * common used file and parameters *////////
	///////////////////////////////////////////////////

	private WeatherFromGrib gribdata; // the weather atmosphere always same, so weather data exclude contrail_index

	// Specific Gas Constants
	final Double Rs_air = 287.058;
	final Double Rs_water = 461.523;

	///////////////////////////////////////////////////
	///////////////////////////////////////////////////
	///////////////////////////////////////////////////

	public wake_vortex_intersection_v1_for_s2(String ConfigFileName_0, String ConfigFileName_1) {
		// read all possible parameters from config-file and initialise variables

		// specify Grib File for weather data

		// gribdata = new
		// WeatherFromGrib("/Users/marco/Desktop/testordner/gribtesttest");
		// gribdata = new
		// WeatherFromGrib("/Users/marco/Desktop/testordner/grib_01_2016-02-07_1200.grib2");

		//////////////////////////////////////////
		/* contrail_0 parameters configurations */
		//////////////////////////////////////////

		// used for calculating extinction and scattering efficiencies as well as
		// assymetrie
		Ext_0 = new Extinktionseffizienz();
		lterr_0 = 10.; // Terrestrische Band

		// micrometer 3.08 ... 99.99
		bsol_0 = 0; // Solares Band
		// band ... das gewuenschte spektrale Band
		// 0 = 0.55 micrometer
		// 1 = 1.35 micrometer
		// 2 = 2.25 micrometer
		// 3 = 3.0125 micrometer
		// 4 = 3.775 micrometer
		// 5 = 4.5 micrometer

		// open config-file
		BufferedReader file_0 = null;
		try {
			file_0 = new BufferedReader(new FileReader(ConfigFileName_0));
			System.out.println("Datei " + ConfigFileName_0 + " gefunden");
		} catch (IOException e) {
			System.err.println("Datei " + ConfigFileName_0 + " nicht gefunden");
		}
		;

		// read parameters from config file
		String eineZeile_0; // zum Zeilenweise einlesen der Datei
		try {
			this.cp_0 = new contrail_parameters();
			// ignore header
			eineZeile_0 = file_0.readLine();
			// Read filename of Radiosondenaufstieg and initialise rsonde
			eineZeile_0 = file_0.readLine();
			// this.rsonde = new radiosonde(file.readLine());
			// Read output filename
			eineZeile_0 = file_0.readLine();
			// Read the line of invalid result file path;
			eineZeile_0 = file_0.readLine();
			this.fn_output_0 = "records\\results\\s2_Ausgabe_LATm50_LON150_vortex_0.txt";
			// Read flight level [m]
			eineZeile_0 = file_0.readLine();
			this.cp_0.altitude = Double.valueOf(file_0.readLine());
			System.out.println("altitude_0: " + this.cp_0.altitude);
			// Read aircraft speed [m/s]
			eineZeile_0 = file_0.readLine();
			this.TAS_0 = Double.valueOf(file_0.readLine());
			System.out.println("TAS_0: " + this.TAS_0);
			// Read time step of flight profile [s]
			eineZeile_0 = file_0.readLine();
			this.dt_fp_0 = Double.valueOf(file_0.readLine());
			System.out.println("dt_fp_0: " + this.dt_fp_0);
			// Calculate distance [m]
			this.distance_0 = this.TAS_0 * this.dt_fp_0;
			System.out.println("distance_0: " + this.distance_0);
			// Read plane mass [kg]
			eineZeile_0 = file_0.readLine();
			this.M_0 = Double.valueOf(file_0.readLine());
			System.out.println("plane mass_0: " + this.M_0);
			// Read fuelflow [kg/s]
			eineZeile_0 = file_0.readLine();
			this.fuelflow_0 = Double.valueOf(file_0.readLine());
			System.out.println("fuel flow_0: " + this.fuelflow_0);
			// Read epsilon turbulenz [?]
			eineZeile_0 = file_0.readLine();
			this.cp_0.eps = Double.valueOf(file_0.readLine());
			System.out.println("eps_0: " + this.cp_0.eps);
			// Read vertical Diffusion contant [m2/s]
			eineZeile_0 = file_0.readLine();
			this.cp_0.Dv = Double.valueOf(file_0.readLine());
			System.out.println("Dv_0 (nicht verwendet): " + this.cp_0.Dv);
			// Read horizontal Diffusion contant [m2/s]
			eineZeile_0 = file_0.readLine();
			this.cp_0.Dh = Double.valueOf(file_0.readLine());
			System.out.println("Dh_0: " + this.cp_0.Dh);
			// Read shearing Diffusion contant [m2/s]
			eineZeile_0 = file_0.readLine();
			this.cp_0.Ds = Double.valueOf(file_0.readLine());
			System.out.println("Ds_0: " + this.cp_0.Ds);
			// Read shearing factor [1/s]
			eineZeile_0 = file_0.readLine();
			this.cp_0.s = Double.valueOf(file_0.readLine());
			System.out.println("s_0: " + this.cp_0.s);
			// Read mean wind velocity [m/s]
			eineZeile_0 = file_0.readLine();
			this.cp_0.vo = Double.valueOf(file_0.readLine());
			System.out.println("vo_0: " + this.cp_0.vo);
			// Read vertical wind velocity [m/s]
			eineZeile_0 = file_0.readLine();
			this.cp_0.vz = Double.valueOf(file_0.readLine());
			System.out.println("vz_0: " + this.cp_0.vz);
			// Read horizontal wind fluctuation [m/s]
			eineZeile_0 = file_0.readLine();
			this.cp_0.uprime = Double.valueOf(file_0.readLine());
			System.out.println("uprime_0: " + this.cp_0.uprime);
			System.out.println("T_0: " + this.cp_0.T + ", RH_0: " + this.cp_0.RH);
			// Read initial time step [s]
			eineZeile_0 = file_0.readLine();
			this.dt_init_0 = Double.valueOf(file_0.readLine());
			System.out.println("dt_0: " + this.dt_init_0);
			// Read initial latitude [deg]
			eineZeile_0 = file_0.readLine();
			this.cp_0.lat = Double.valueOf(file_0.readLine());
			System.out.println("lat_0: " + this.cp_0.lat);
			// Read initial longitude [deg]
			eineZeile_0 = file_0.readLine();
			this.cp_0.lon = Double.valueOf(file_0.readLine());
			System.out.println("lon_0: " + this.cp_0.lon);
			eineZeile_0 = file_0.readLine();
			// Read filename of Gribfile and initialise gribdata
			String gribfilename_0 = file_0.readLine();
			this.gribdata = new WeatherFromGrib(gribfilename_0);

			file_0.close(); // Close config-file
		} catch (IOException ex) {
			System.err.println("Fehler beim Einlesen");
		}
		;

		/*
		 * double lat = this.cp.lat; double lon = this.cp.lon; double altitude =
		 * this.cp.altitude;
		 */

		// Read Weather Data from Gribfile

		this.cp_0.vx = gribdata.getWindUFromAltitude(this.cp_0.lat, this.cp_0.lon, this.cp_0.altitude); // U-Component
																										// of wind
		this.cp_0.vy = gribdata.getWindVFromAltitude(this.cp_0.lat, this.cp_0.lon, this.cp_0.altitude); // V-Component
																										// of wind
		this.cp_0.Tflight = gribdata.getTempFromAltitude(this.cp_0.lat, this.cp_0.lon, this.cp_0.altitude); // K
		this.cp_0.T = this.cp_0.Tflight; // K
		System.out.println("cp_0.T: " + this.cp_0.T);
		this.cp_0.RH = gribdata.getHumidityFromAltitude(this.cp_0.lat, this.cp_0.lon, this.cp_0.altitude); // %

		// this.cp.density = weatherdata_helper(gribdata, this.cp, lat0, lon0, altitude)
		this.cp_0.density = density(this.cp_0); // FIXME
		this.cp_0.vz = vz(this.cp_0);
		this.cp_0.Theta = theta(this.cp_0);
		this.cp_0.thetadz = dthetadz(this.cp_0);
		this.cp_0.eps = eps_calculator(this.cp_0);
		System.out.println("cp_0.eps: " + this.cp_0.eps);

		// this.cp.eps = eps_calculator(this.cp.lat, this.cp.lon, this.cp.altitude);

		// calculate initial state of Contrail with GRIB
		initial_state istate_0 = this.Contrail_depth(this.M_0, this.TAS_0,
				gribdata.getPressureFromAltitude(this.cp_0.lat, this.cp_0.lon, this.cp_0.altitude), this.cp_0.Tflight,
				this.cp_0.Theta, this.cp_0.thetadz, this.cp_0.altitude, this.cp_0.eps);
		this.cp_0.sigma0v = (float) (istate_0.Height / 2.2); // 50m
		this.cp_0.sigma0h = (float) (istate_0.Width / 2.2); // 20m

		// Calculate horizontal Diffusion coefficient [m2/s]
		// According to Schuhmann, J. Geophys. Res. 100, 14147(1995)
		// this.cp.Dh = 0.1 * this.cp.uprime * this.cp.sigma0h; // c_h * u' * sigmah
		// System.out.println("Dh (errechnet) "+ this.cp.Dh);
		// Read vertical Diffusion coefficient [m2/s]
		this.cp_0.Dv = this.Calc_Dv(this.cp_0.altitude, this.cp_0);
		System.out.println("Dv_0 (errechnet): " + this.cp_0.Dv);

		// eingefügt
		this.cp_0.Dh = this.Calc_Dh(this.cp_0.eps);
		System.out.println("Dh_0 (interpolert): " + this.cp_0.Dh);
		System.out.println("////////////////////////////////////////////////////////////////////////////");
		System.out.println("Above parameters belong contrail_0, constructor for contrail_0 has finished!");
		System.out.println("//////////////////////////////////////////////////////////////////////////// \n\n");

		//////////////////////////////////////////
		/* contrail_1 parameters configurations */
		//////////////////////////////////////////

		// used for calculating extinction and scattering efficiencies as well as
		// assymetrie
		Ext_1 = new Extinktionseffizienz();
		lterr_1 = 10.; // Terrestrische Band

		// micrometer 3.08 ... 99.99
		bsol_1 = 0; // Solares Band
		// band ... das gewuenschte spektrale Band
		// 0 = 0.55 micrometer
		// 1 = 1.35 micrometer
		// 2 = 2.25 micrometer
		// 3 = 3.0125 micrometer
		// 4 = 3.775 micrometer
		// 5 = 4.5 micrometer

		// open config-file
		BufferedReader file_1 = null;
		try {
			file_1 = new BufferedReader(new FileReader(ConfigFileName_1));
			System.out.println("Datei " + ConfigFileName_1 + " gefunden");
		} catch (IOException e) {
			System.err.println("Datei " + ConfigFileName_1 + " nicht gefunden");
		}
		;

		// read parameters from config file
		String eineZeile_1; // zum Zeilenweise einlesen der Datei
		try {
			this.cp_1 = new contrail_parameters();
			// ignore header
			eineZeile_1 = file_1.readLine();
			// Read filename of Radiosondenaufstieg and initialise rsonde
			eineZeile_1 = file_1.readLine();
			// this.rsonde = new radiosonde(file.readLine());
			// Read output filename
			eineZeile_1 = file_1.readLine();
			// Read the line of invalid result file path;
			eineZeile_1 = file_1.readLine();
			this.fn_output_1 = "records\\results\\s2_Ausgabe_LATm50_LON150_vortex_1.txt";
			// Read flight level [m]
			eineZeile_1 = file_1.readLine();
			this.cp_1.altitude = Double.valueOf(file_1.readLine());
			System.out.println("altitude_1: " + this.cp_1.altitude);
			// Read aircraft speed [m/s]
			eineZeile_1 = file_1.readLine();
			this.TAS_1 = Double.valueOf(file_1.readLine());
			System.out.println("TAS_1: " + this.TAS_1);
			// Read time step of flight profile [s]
			eineZeile_1 = file_1.readLine();
			this.dt_fp_1 = Double.valueOf(file_1.readLine());
			System.out.println("dt_fp_1: " + this.dt_fp_1);
			// Calculate distance [m]
			this.distance_1 = this.TAS_1 * this.dt_fp_1;
			System.out.println("distance_1: " + this.distance_1);
			// Read plane mass [kg]
			eineZeile_1 = file_1.readLine();
			this.M_1 = Double.valueOf(file_1.readLine());
			System.out.println("plane mass_1: " + this.M_1);
			// Read fuelflow [kg/s]
			eineZeile_1 = file_1.readLine();
			this.fuelflow_1 = Double.valueOf(file_1.readLine());
			System.out.println("fuel flow_1: " + this.fuelflow_1);
			// Read epsilon turbulenz [?]
			eineZeile_1 = file_1.readLine();
			this.cp_1.eps = Double.valueOf(file_1.readLine());
			System.out.println("eps_1: " + this.cp_1.eps);
			// Read vertical Diffusion contant [m2/s]
			eineZeile_1 = file_1.readLine();
			this.cp_1.Dv = Double.valueOf(file_1.readLine());
			System.out.println("Dv_1 (nicht verwendet): " + this.cp_1.Dv);
			// Read horizontal Diffusion contant [m2/s]
			eineZeile_1 = file_1.readLine();
			this.cp_1.Dh = Double.valueOf(file_1.readLine());
			System.out.println("Dh_1: " + this.cp_1.Dh);
			// Read shearing Diffusion contant [m2/s]
			eineZeile_1 = file_1.readLine();
			this.cp_1.Ds = Double.valueOf(file_1.readLine());
			System.out.println("Ds_1: " + this.cp_1.Ds);
			// Read shearing factor [1/s]
			eineZeile_1 = file_1.readLine();
			this.cp_1.s = Double.valueOf(file_1.readLine());
			System.out.println("s_1: " + this.cp_1.s);
			// Read mean wind velocity [m/s]
			eineZeile_1 = file_1.readLine();
			this.cp_1.vo = Double.valueOf(file_1.readLine());
			System.out.println("vo_1: " + this.cp_1.vo);
			// Read vertical wind velocity [m/s]
			eineZeile_1 = file_1.readLine();
			this.cp_1.vz = Double.valueOf(file_1.readLine());
			System.out.println("vz_1: " + this.cp_1.vz);
			// Read horizontal wind fluctuation [m/s]
			eineZeile_1 = file_1.readLine();
			this.cp_1.uprime = Double.valueOf(file_1.readLine());
			System.out.println("uprime_1: " + this.cp_1.uprime);
			System.out.println("T_1: " + this.cp_1.T + ", RH_1: " + this.cp_1.RH);
			// Read initial time step [s]
			eineZeile_1 = file_1.readLine();
			this.dt_init_1 = Double.valueOf(file_1.readLine());
			System.out.println("dt_1: " + this.dt_init_1);
			// Read initial latitude [deg]
			eineZeile_1 = file_1.readLine();
			this.cp_1.lat = Double.valueOf(file_1.readLine());
			System.out.println("lat_1: " + this.cp_1.lat);
			// Read initial longitude [deg]
			eineZeile_1 = file_1.readLine();
			this.cp_1.lon = Double.valueOf(file_1.readLine());
			System.out.println("lon_1: " + this.cp_1.lon);
			eineZeile_1 = file_1.readLine();
			// Read filename of Gribfile and initialise gribdata
			// gribdata is the common used weather file, so here we do not need to reload
			// it!
			// String gribfilename_1 = file_1.readLine();
			// this.gribdata = new WeatherFromGrib(gribfilename_1);

			file_1.close(); // Close config-file
		} catch (IOException ex) {
			System.err.println("Fehler beim Einlesen");
		}
		;

		/*
		 * double lat = this.cp.lat; double lon = this.cp.lon; double altitude =
		 * this.cp.altitude;
		 */

		// Read Weather Data from Gribfile

		this.cp_1.vx = gribdata.getWindUFromAltitude(this.cp_1.lat, this.cp_1.lon, this.cp_1.altitude); // U-Component
																										// of wind
		this.cp_1.vy = gribdata.getWindVFromAltitude(this.cp_1.lat, this.cp_1.lon, this.cp_1.altitude); // V-Component
																										// of wind
		this.cp_1.Tflight = gribdata.getTempFromAltitude(this.cp_1.lat, this.cp_1.lon, this.cp_1.altitude); // K
		this.cp_1.T = this.cp_1.Tflight; // K
		System.out.println("cp_1.T: " + this.cp_1.T);
		this.cp_1.RH = gribdata.getHumidityFromAltitude(this.cp_1.lat, this.cp_1.lon, this.cp_1.altitude); // %

		// this.cp.density = weatherdata_helper(gribdata, this.cp, lat0, lon0, altitude)
		this.cp_1.density = density(this.cp_1); // FIXME
		this.cp_1.vz = vz(this.cp_1);
		this.cp_1.Theta = theta(this.cp_1);
		this.cp_1.thetadz = dthetadz(this.cp_1);
		this.cp_1.eps = eps_calculator(this.cp_1);
		System.out.println("cp_1.eps: " + this.cp_1.eps);
		// this.cp.eps = eps_calculator(this.cp.lat, this.cp.lon, this.cp.altitude);

		// calculate initial state of Contrail with GRIB
		initial_state istate_1 = this.Contrail_depth(this.M_1, this.TAS_1,
				gribdata.getPressureFromAltitude(this.cp_1.lat, this.cp_1.lon, this.cp_1.altitude), this.cp_1.Tflight,
				this.cp_1.Theta, this.cp_1.thetadz, this.cp_1.altitude, this.cp_1.eps);
		this.cp_1.sigma0v = (float) (istate_1.Height / 2.2); // 50m
		this.cp_1.sigma0h = (float) (istate_1.Width / 2.2); // 20m

		// Calculate horizontal Diffusion coefficient [m2/s]
		// According to Schuhmann, J. Geophys. Res. 100, 14147(1995)
		// this.cp.Dh = 0.1 * this.cp.uprime * this.cp.sigma0h; // c_h * u' * sigmah
		// System.out.println("Dh (errechnet) "+ this.cp.Dh);
		// Read vertical Diffusion coefficient [m2/s]
		this.cp_1.Dv = this.Calc_Dv(this.cp_1.altitude, this.cp_1);
		System.out.println("Dv_1 (errechnet): " + this.cp_1.Dv);

		// eingefügt
		this.cp_1.Dh = this.Calc_Dh(this.cp_1.eps);
		System.out.println("Dh_1 (interpolert): " + this.cp_1.Dh);
		System.out.println("////////////////////////////////////////////////////////////////////////////");
		System.out.println("Above parameters belong contrail_1, constructor for contrail_1 has finished!");
		System.out.println("//////////////////////////////////////////////////////////////////////////// \n\n");

		System.out.println("////////////////////////////////////////////////////////////////////////////");
		System.out.println("/////////////contrail_0 and contrail_1 configuration completed!/////////////");
		System.out.println("//////////////////////////////////////////////////////////////////////////// \n\n");

		return;
	}

	public double density(contrail_parameters cp) {

		/***
		 * 
		 * Calculating air density depending on Relative Humidity
		 * 
		 * 
		 * // From Common Sources
		 * 
		 * density = p * M / (R * T)
		 * 
		 * p ... Pressure M ... molar Mass R ... Universal Gas Constant 8 J / (kg * K) T
		 * ... Temperature
		 * 
		 * Rs = R / M
		 * 
		 * density = p / (Rs * T)
		 * 
		 * Rs_air ... 287,058 J / (kg * K) Rs_water ... 461.523 J / (kg * K) pd ...
		 * Saturation vapour pressure p ... Air Pressure RH ... Relative Humidity
		 * 
		 * // Relative Values for air and water to get Specific Gas Constant for humid
		 * air // as derived from Dalton´s Law and also shown on:
		 * de.wikipedia.org/wiki/Luftdichte
		 * 
		 * Rs = Rs_air / (1 - RF * pd/p * (1 - Rs_air / Rs_water))
		 * 
		 * 
		 * // Saturation vapour pressure at air temperature T from WMO, // Source: WMO
		 * GUIDE TO METEOROLOGICAL INSTRUMENTS AND METHODS OF OBSERVATION // WMO-No. 8
		 * (2014 edition, Updated in 2017), page 346 ) //
		 * https://library.wmo.int/doc_num.php?explnum_id=10179
		 * 
		 * pd(T) = 6.107 * 10 ^ (7.5*T/(237.3+T)) ... T has to be in °C
		 * 
		 * OR after SONTAG
		 * 
		 */

		// sat_vap_pres after WMO GUIDE TO METEOROLOGICAL INSTRUMENTS AND METHODS OF
		// OBSERVATION WMO-No. 8 (2014 edition, Updated in 2017), page 346 )
		// https://library.wmo.int/doc_num.php?explnum_id=10179
		// NOT IN USE, replaced by SONTAG
		// this.cp.sat_vap_pres = 6.107 *
		// Math.pow(10,(7.5*(this.cp.T-273.15)/(237.3+(this.cp.T-273.15)))); //with Temp
		// Conversion to deg Celsius

		// sat_vap_pres after: SONTAG
		cp.sat_vap_pres = Math.exp(-6096.9385 / cp.T + 16.635794 - 2.711193E-2 * cp.T + 1.673952E-5 * cp.T * cp.T
				+ 2.433502 * Math.log(cp.T));
		cp.Rs_humidair = Rs_air / 1 - (cp.RH * (cp.sat_vap_pres
				/ gribdata.getPressureFromAltitude(cp.lat, cp.lon, cp.altitude) * (1 - Rs_air / Rs_water)));
		cp.density = 100 * gribdata.getPressureFromAltitude(cp.lat, cp.lon, cp.altitude) / (cp.Rs_humidair * cp.T); // gribdata
																													// pressure
																													// value
																													// is
																													// in
																													// hPa

		return cp.density;

	}

	public double vz(contrail_parameters cp) {

		/***
		 * 
		 * Calculating Vertical Wind vz
		 * 
		 * // Vertical velocity VVEL is the speed at which the air is rising or sinking.
		 * // Measured in negative microbars per second as a distance. // Instead of
		 * meters the distance measured is pressure, so a negative distance in pressure
		 * // is actually a positive distance in altitude, hence a negative microbar
		 * distance is used here.
		 * 
		 * // p V = m Rs_humidair T ==> rho = p / Rs_humidair T // vvel = omega ... vz =
		 * -vvel Rs_humidair T / p g
		 * 
		 */

		// unter Einberechnung der Feuchtigkeit
		cp.sat_vap_pres = Math.exp(-6096.9385 / cp.T + 16.635794 - 2.711193E-2 * cp.T + 1.673952E-5 * cp.T * cp.T
				+ 2.433502 * Math.log(cp.T)); // in hPa
		cp.Rs_humidair = Rs_air / 1 - (cp.RH * (cp.sat_vap_pres
				/ gribdata.getPressureFromAltitude(cp.lat, cp.lon, cp.altitude) * (1 - Rs_air / Rs_water))); // pressure
																												// in
																												// hPa
		cp.vz = -(gribdata.getVvelFromAltitude(cp.lat, cp.lon, cp.altitude)) * cp.Rs_humidair
				* (gribdata.getTempFromAltitude(cp.lat, cp.lon, cp.altitude))
				/ (9.81 * 100 * (gribdata.getPressureFromAltitude(cp.lat, cp.lon, cp.altitude)));

		return cp.vz;

	}

	public double theta(contrail_parameters cp) {

		// Source: Approximation according to D. Bolton, 1980 as shown on:
		// http://www.uni-koeln.de/math-nat-fak/geomet/meteo/winfos/radiosonden/Europa/radiosondengrafiken.pdf

		cp.sat_vap_pres = Math.exp(-6096.9385 / cp.T + 16.635794 - 2.711193E-2 * cp.T + 1.673952E-5 * cp.T * cp.T
				+ 2.433502 * Math.log(cp.T)); // in hPa
		Double mixture = (Rs_air / Rs_water) * (cp.sat_vap_pres
				/ ((gribdata.getPressureFromAltitude(cp.lat, cp.lon, cp.altitude)) - cp.sat_vap_pres));
		cp.Theta = cp.T * Math.pow((1000 / (gribdata.getPressureFromAltitude(cp.lat, cp.lon, cp.altitude))),
				0.2854 * (1 - 0.28 * mixture));

		return cp.Theta;

	}

	public double dthetadz(contrail_parameters cp) {

		// calculation of theta dz
		int x = 1000; // step width for interpolating theta dz in metres

		cp.sat_vap_pres = Math.exp(-6096.9385 / cp.T + 16.635794 - 2.711193E-2 * cp.T + 1.673952E-5 * cp.T * cp.T
				+ 2.433502 * Math.log(cp.T)); // in hPa
		Double mixture = (Rs_air / Rs_water) * (cp.sat_vap_pres
				/ ((gribdata.getPressureFromAltitude(cp.lat, cp.lon, cp.altitude)) - cp.sat_vap_pres));

		if (altitude_0 - x < 0) {
			cp.thetadz = cp.T
					* Math.pow((1000 / (gribdata.getPressureFromAltitude(cp.lat, cp.lon, cp.altitude + x))),
							0.2854 * (1 - 0.28 * mixture))
					- (cp.T * Math.pow((1000 / (gribdata.getPressureFromAltitude(cp.lat, cp.lon, cp.altitude))),
							0.2854 * (1 - 0.28 * mixture)));
		} else {
			cp.thetadz = cp.T
					* Math.pow((1000 / (gribdata.getPressureFromAltitude(cp.lat, cp.lon, cp.altitude + x))),
							0.2854 * (1 - 0.28 * mixture))
					- (cp.T * Math.pow((1000 / (gribdata.getPressureFromAltitude(cp.lat, cp.lon, cp.altitude - x))),
							0.2854 * (1 - 0.28 * mixture)));
		}
		return cp.thetadz;

	}

	public double eps_calculator(contrail_parameters cp) {// (double lat, double lon, double altitude) {

		/**
		 * 
		 * Calculation of Eddy Dissipation Rate after formula eps = Math.exp (C1 + C2 *
		 * (lnD - lnDe) / SD) as shown in: SHARMAN AND PEARSON "Prediction of Energy
		 * Dissipation Rates for Aviation Turbulence. Part I: Forecasting Nonconvective
		 * Turbulence" (2016)
		 * 
		 * Average Values above altitude 0m for
		 * 
		 * C1 = -2.572 C2 = 0.5067
		 * 
		 **/

		double lat = cp.lat;
		double lon = cp.lon;
		double altitude = cp.altitude;

		// defining the full degrees range around current lat/lon which shall be
		// considered determining the eps
		int x = 15; // 15 degrees range was tested to be sufficient
		int y = x;

		if (lat < -90 || lat > 90 || lon < -180 || lon > 180) {
			throw new IllegalArgumentException("One or more values out of range (lat = " + lat + " ,lon = " + lon);
		}
		if (y < 0 || y > 90 || x < 0 || x > 180) {
			throw new IllegalArgumentException(
					"One or more values out of range (full degrees range around current lat/lon to be considered... lateral = "
							+ y + " ,longitudinal = " + x);
		}

		double current_lon = -777;
		double current_lat = -777;
		double current_vvel = 0;

		double SD_lnVVEL = 0;
		double Mean_lnVVEL = 0.; // Standard Deviation and Mean Value of logarithmic VVEL values between lat/lon
									// start and end

		// double Mean_lnVVELx = 0;
		// double testsjkjksjs = 0;

		// longitudinal starting location is from lon - y and calculation is always
		// performed in eastern direction until lon +y
		// lateral starting location is from lat - x and wents to final location lat + x
		// (even if that means to go south first and then north or vice versa in case
		// boundary values for lat = +- 90 deg are tangered

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
		cp.eps = Math.exp(C1 + C2 * (Math.log(Math.abs(gribdata.getVvelFromAltitude(lat, lon, altitude))) - Mean_lnVVEL)
				/ SD_lnVVEL);

		return cp.eps;
	}

	public double Calc_Dh(double eps) {

		/*
		 * According to Tab 3.2, Rosenow, "Optical Properties of Contrails", pg. 68 Eddy
		 * dissipation rate ε [m2s−3] ... column 0 Horizontal diffusivity Dh [m2s−2] ...
		 * column 1
		 * 
		 * Interpolation of values No extrapolation beyond lower and upper boundaries!
		 * 
		 * Returns interpolated value for Dh
		 */

		double[][] eps_Dh_Table = { { Math.pow(10, -6), 5.00 }, { 5 * Math.pow(10, -6), 8.75 },
				{ Math.pow(10, -5), 12.50 }, { 5 * Math.pow(10, -5), 16.254 }, { Math.pow(10, -4), 20.00 } };
		double returnvalue = 0;

		if (eps < eps_Dh_Table[0][0]) {
			returnvalue = eps_Dh_Table[0][1]; // lower boundary
		} else if (eps >= eps_Dh_Table[(eps_Dh_Table.length) - 1][0]) {
			returnvalue = eps_Dh_Table[(eps_Dh_Table.length) - 1][1]; // upper boundary
		} else {
			for (int i = 1; i <= ((eps_Dh_Table.length) - 1); i++) {
				if (eps >= eps_Dh_Table[i - 1][0] && eps < eps_Dh_Table[i][0]) {
					returnvalue = ((eps - eps_Dh_Table[i - 1][0]) / (eps_Dh_Table[i][0] - eps_Dh_Table[i - 1][0])
							* (eps_Dh_Table[i][1] - eps_Dh_Table[i - 1][1]) + eps_Dh_Table[i - 1][1]);
				}
			}

		}

		return returnvalue;

	}

	public double Calc_Dv(double height, contrail_parameters cp) {

		double Nsquare = 9.81 / theta(cp) * dthetadz(cp); // square of Brunt Vaisala frequency [a.u. 1/s^2]
		return 0.1 * cp.eps / Nsquare;

	}

	public initial_state Contrail_depth(double M, double TAS, double p, double T, double Theta, double dThetadz,
			double z, double eps) {
		// M...aircraft mass [kg]
		// TAS... True airspeed [m/s]
		// p...pressure [hPa]
		// T...Temperature [K]
		// Theta...potential temperature [K]
		// dThetadz...derivation potential temperature with height [K/m]
		// z...Altitude [m]
		// eps...Eddy dissipation rate [m2/s3]
		//
		// calculation of initial contrail shape with the help of wake vortex
		// characteristics
		//

		// aircraft parameters of A320
		double rc = 1.1935; // [m] core radius
		double B = 34.1; // [m] span width
		double bo = Math.PI / 4. * B; // [m] initial vortex spacing

		// to be calculated in this function
		double r = 10; // vortex radius at the beginning of the dispersion regime [m]
		double h = 20; // contrail height [m]
		double w = 30; // contrail width [m]

		// characteristic scales
		// bo ... initial vortex spacing
		// wo ... [m/s] initial tangential velocity
		// Gammao ... [m2/s] initial vortex circulation
		// tprime...[s] characteristic time scale
		double rho = p * 100 / (287.058 * T); // density of air
		// FIXME
		// double rho=density();
		// System.out.println("rho "+rho);
		double Gammao = M * 9.81 / (rho * bo * TAS);
		double wo = Gammao / (2. * Math.PI * bo);
		double tprime = 2 * Math.PI * Math.pow(bo, 2.) / Gammao;

		// calculation of initial vortex radius
		// with the help of Dispersion model (tangential velocity) according to
		// Holzaepfel 2003
		// r is the distance, where w=1/e^(1/2) * w_max
		r = (Math.exp(0.5) + Math.pow(Math.E - 1, 0.5)) * rc;

		// calculation of contrail width
		w = 2 * r + bo;

		// calculation of contrail height
		double epsstar = Math.pow(eps * bo, 1 / 3.) / wo;

		// System.out.println("epsstar "+epsstar);
		double Nstar = Math.sqrt(9.81 / Theta * dThetadz) * tprime; // Brunt Vaisala frequency [a.u. 1/s]
		// System.out.println("Nstar "+Nstar);
		double Tstar = 1;

		double T20star = 5;
		if (epsstar > 0.2535) {
			T20star = 0.804 * Math.pow(epsstar, -0.75) - 1;
		} else if (epsstar > 0.0235) {
			int i = 0;
			double f = 1;
			double df = 1;
			double dx = 1;
			Tstar = 1; // starting value for Newton-iteration
			do {
				f = Math.pow(Tstar, 0.25) * Math.exp(-0.7 * Tstar) - epsstar;
				df = (0.25 * Math.pow(Tstar, -0.75) - 0.7 * Math.pow(Tstar, 0.25)) * Math.exp(-0.7 * Tstar);
				dx = f / df;
				Tstar = Tstar - dx;
				i = i + 1;
			} while (Math.abs(dx) > 1E-3 && i < 11);
			if (i >= 11)
				System.out.println("Iteration not long enough for Tstar");
			T20star = Tstar - 1.;
		}
		;

		// System.out.println("T20star "+T20star);
		double T2star = T20star * Math.exp(-0.185 * T20star * Nstar);

		double wstar = 1. - Math.exp(-1.257 * Math.pow(0.4 * bo / rc, 2));

		h = wstar * T2star * bo + 2 * r;
		System.out.println("T2star: " + T2star + " wstar: " + wstar);
		System.out.println("T2: " + (T2star * tprime));
		System.out.println("r: " + r + ", h: " + h + ", w: " + w);
		System.out.println("w* = " + wstar + ", w0 = " + wo + ", w = " + (wstar * wo));
		return new initial_state(r, h, w);
	}

	// calculates area of contrail cross-section
	public double ContrailCrossSection(double t, contrail_parameters cp) {
		// t...time [s]
		// Dv, Dh, Ds... vertical, horizontal and sheared diffusion [m2/s]
		// vo...mean value of wind velocity at flight level [m/s]
		// s....constant wind shear 1/s

		double sigmav = 2 * cp.Dv * t + Math.pow(cp.sigma0v, 2); // "s" was set to 2 initially?!
		cp.s = (gribdata.getWindspeedFromAltitude(cp.lat, cp.lon, cp.altitude + sigmav)
				- gribdata.getWindspeedFromAltitude(cp.lat, cp.lon, cp.altitude - sigmav)) / (2 * sigmav); // calculates
																											// the
																											// change of
																											// windspeed
																											// between
																											// upper and
																											// lower
																											// boundary
																											// of
																											// contrail
																											// per metre
		double sigmah = 2. / 3 * cp.s * cp.s * cp.Dv * Math.pow(t, 3)
				+ (2 * cp.Ds + cp.s * Math.pow(cp.sigma0v, 2)) * cp.s * Math.pow(t, 2) + 2 * cp.Dh * t
				+ Math.pow(cp.sigma0h, 2);
		double sigmas = cp.s * cp.Dv * Math.pow(t, 2) + (2 * cp.Ds + cp.s * Math.pow(cp.sigma0v, 2)) * t;

		// System.out.println("cp name: " + cp.toString());
		// System.out.println("Sigmas in func CCS: " + sigmas);
		// System.out.println("t in CCS: " + t);
		double detsigma = sigmah * sigmav - Math.pow(sigmas, 2);
		double A = 2. * Math.PI * Math.sqrt(detsigma); // Schumann 1995

		// System.out.println("Cross section size (A): " + A);
		// System.out.println("Cross section size (Eclipse cal): " +
		// Math.PI*2.2*sigmah*2.2*sigmav); this size calculation not considered!

		return A;
	}

	// compared to the same func in wake_vortex.java, the input "double z" was
	// removed because "z" unused in this func.
	public double AdiabaticCurve(contrail_parameters cp) {
		// SATURATION PRESSUREs FROM SONNTAG, Temp IN K, PSAT IN hPa
		// double
		// eStar=Math.exp(-6096.9385/cp.T+16.635794-2.711193E-2*cp.T+1.673952E-5*cp.T*cp.T+2.433502*Math.log(cp.T));
		double eStarIce = Math.exp(-6024.5282 / cp.Tflight + 24.721994 + 1.0613868E-2 * cp.Tflight
				- 1.3198825E-5 * cp.Tflight * cp.Tflight - 0.49382577 * Math.log(cp.Tflight));
		double deStarIce = (6024.5282 / Math.pow(cp.Tflight, 2) - 1.0613868E-2 - 2 * 1.3198825E-5 * cp.Tflight
				- 0.49382577 / cp.Tflight) * eStarIce;

		// Saettigungsmischungsverhaeltnis
		double pressure = gribdata.getPressureFromAltitude(cp.lat, cp.lon, cp.altitude);
		double mStar = 0.622 * eStarIce / (pressure - eStarIce); // kg/kg
		// LatentHeat = Verdampfungswaerme
		double L = 2.5e6; // J/kg
		// spezifische Waerme bei konstantem Druck
		double c_p = 1004; // J/(kg K)
		// spezifische Gaskonstante von Luft
		double RL = 287; // J/(kg K)
		// Adiabatic Curve
		double Gammas = 9.81 / c_p * (1 + L * mStar / (RL * cp.Tflight))
				/ (1 + L * mStar * deStarIce / (c_p * eStarIce));

		return Gammas; // [K/m]
	}

	public double kinvisc_Suth(double p, double T) {
		// p...pressure[hPa], T...Temperature[K],
		// calculates kinematic viscosity [m2/s] of air according to
		// Philosophical Magazine 5 (1893) page 507 W. Sutherland

		// Sutherland parameter
		double T0 = 273; // in K
		double mu0 = 1.716e-5; // in (Ns)/m2
		double C = 111; // in K

		// calculation of density
		double rho = 1; // density of air
		rho = p * 100 / (287.058 * T);

		// calculation of viscosity
		double mu = 1; // dynamic viscosity
		mu = mu0 * (T0 + C) / (T + C) * Math.pow(T / T0, 1.5); // Sutherland formular
		return mu / rho; // kinematic viscosity
	}

	public double dynvisc_Suth(double T) {
		// p...pressure[hPa], T...Temperature[K],
		// calculates dynamic viscosity [kg/m s] of air according to
		// Philosophical Magazine 5 (1893) page 507 W. Sutherland

		// Sutherland parameter
		double T0 = 273; // in K
		double mu0 = 1.716e-5; // in (Ns)/m2
		double C = 111; // in K

		// calculation of density
		// double rho=1; //density of air
		// rho=p*100/(287.058*T);

		// calculation of viscosity
		double mu = 1; // dynamic viscosity
		mu = mu0 * (T0 + C) / (T + C) * Math.pow(T / T0, 1.5); // Sutherland formular
		return mu; // dynamic viscosity
	}

	public double SedimentationSpeed(double dice, contrail_parameters cp) {
		// Stokessches Sedimentations-Gesetz
		double vsed = 2 * Math.pow(dice / 2., 2) * 917 * 9.81 / (9 * dynvisc_Suth(cp.T));
		return vsed;
	}

	public double lambda_e(double sigma_h, double sigma_v, double sigma_s, double r, double Nice) {
		double Qabs = 4.e-5;
		double Qsca = 2.;
		double le = 0.;
		double detsigma = Math.pow(sigma_h * sigma_v, 2) - Math.pow(sigma_s, 2);
		double n = Nice / (this.distance_0 * 2 * Math.PI * Math.pow(detsigma, 0.5));
		le = 1. / ((Qabs + Qsca) * n * Math.PI * r * r);
		return le;
	}

	public double Re(double r, double vsed, double vkin) {
		// Berechnung der Reynolds-Zahl [einheitenlos]
		//
		// r... Partikelradius [m]
		// vsed... Sedimentationsgeschwindigkeit [m/s]
		// vkin... kinematische Viskositaet [m2/s]
		return 2 * r * vsed / vkin;
	}

	// if intersected, the two contrails should be in the same level， but here
	// besides horizontal, still check the vertical as well
	public boolean intersection_check(contrail_parameters cp0, double sigmaV0, double sigmaH0, contrail_parameters cp1,
			double sigmaV1, double sigmaH1) {
		// Initial flag
		boolean intersected = false;
		// Earth radius
		final double earth_r = 6371000.0;

		double cp0_half_height = sigmaV0 * 2.2 / 2;
		double cp0_half_width = sigmaH0 * 2.2 / 2;
		double cp1_half_height = sigmaV1 * 2.2 / 2;
		double cp1_half_width = sigmaH1 * 2.2 / 2;

		// Calculate horizontal distance:
		// Converting latitude and longitude to radians
		double d_lat = Math.toRadians(cp1.lat - cp0.lat);
		double d_lon = Math.toRadians(cp1.lon - cp0.lon);
		double lat0 = Math.toRadians(cp0.lat);
		double lat1 = Math.toRadians(cp1.lat);

		// Applying the Haversine formula
		double a = Math.sin(d_lat / 2) * Math.sin(d_lat / 2)
				+ Math.sin(d_lon / 2) * Math.sin(d_lon / 2) * Math.cos(lat0) * Math.cos(lat1);
		double dis_h = 2 * earth_r * Math.asin(Math.sqrt(a)); // Return distance (meters)

		// Calculate vertical distance:
		double dis_v = Math.abs(cp0.altitude - cp1.altitude);

		// Check if intersected:
		if ((dis_h < cp0_half_width + cp1_half_width) && (dis_v < cp0_half_height + cp1_half_height)) {
			intersected = true;
		} else {
			intersected = false;
		}

		return intersected;
	}

	// if intersected, the two contrails should be in the same level， but here
	// besides horizontal, still check the vertical as well
	// parameter c0_full_cover_c1 and c0_full_cover_c1 only affect in the case of
	// angle = 0.
	public double intersection_volumn_cal(double ang, contrail_parameters cp0, double sigmaV0, double sigmaH0,
			double CCS0, contrail_parameters cp1, double sigmaV1, double sigmaH1, double CCS1, double distance) {

		// Earth radius
		final double earth_r = 6371000.0;

		double cp0_half_height = sigmaV0 * 2.2 / 2;
		double cp0_half_width = sigmaH0 * 2.2 / 2;
		double cp1_half_height = sigmaV1 * 2.2 / 2;
		double cp1_half_width = sigmaH1 * 2.2 / 2;

		// Calculate horizontal distance:
		// Converting latitude and longitude to radians
		double d_lat = Math.toRadians(cp1.lat - cp0.lat);
		double d_lon = Math.toRadians(cp1.lon - cp0.lon);
		double lat0 = Math.toRadians(cp0.lat);
		double lat1 = Math.toRadians(cp1.lat);

		// Applying the Haversine formula
		double a = Math.sin(d_lat / 2) * Math.sin(d_lat / 2)
				+ Math.sin(d_lon / 2) * Math.sin(d_lon / 2) * Math.cos(lat0) * Math.cos(lat1);
		double dis_h = 2 * earth_r * Math.asin(Math.sqrt(a)); // Return distance (meters)

		// Calculate vertical distance:
		// double dis_v = Math.abs(cp0.altitude - cp1.altitude);

		double inter_volumn = 0.;

		if (ang == 0) {
//			if (c0_full_cover_c1) {
//				inter_volumn = CCS1*distance;
//				System.out.println("c0_full_cover_c1 is true, parallel, c0 fully covers c1!");
//				return inter_volumn;
//			} else if (c0_partly_cover_c1) {
//				
//				double inter_h = (cp0_half_width+cp1_half_width-dis_h)/2.;
//				double inter_v =  cp0_half_height*Math.pow(1-Math.pow((cp0_half_width-inter_h)/(cp0_half_width),2), 0.5);
//				double inter_size = Math.PI*inter_h*inter_v;
//				if (inter_size >= CCS1) {
//					inter_size = CCS1;
//				}
//				inter_volumn = inter_size*distance;
//				System.out.println("c0_partly_cover_c1 is true, parallel, c0 partly covers c1!");
//				return inter_volumn;
//			} else {
			if (dis_h + cp1_half_width <= cp0_half_width) {
				// full covered
				inter_volumn = CCS1 * distance;
				System.out.println("intersection_volumn_cal: Parallel, c0 fully covers c1!");
				return inter_volumn;
			} else if ((dis_h + cp1_half_width > cp0_half_width) && (dis_h <= cp0_half_width + cp1_half_width)) {
				// partly covered
				double inter_h = (cp0_half_width + cp1_half_width - dis_h) / 2.;
				double inter_v = cp0_half_height
						* Math.pow(1 - Math.pow((cp0_half_width - inter_h) / (cp0_half_width), 2), 0.5);
				// System.out.println("inter_h: " + inter_h);
				// System.out.println("inter_v: " + inter_v);
				double inter_size = Math.PI * inter_h * inter_v;
				if (inter_size >= CCS1) {
					inter_size = CCS1;
				}
				inter_volumn = inter_size * distance;
				System.out.println("intersection_volumn_cal: Parallel, c0 partly covers c1!");
				return inter_volumn;
			} else {
				inter_volumn = 0;
				System.out.println("intersection_volumn_cal: Parallel, c0 and c1 have separated already!");
				return inter_volumn;
			}
//			}				

		} else if (ang == 90.0) {
			// intersected
			double top_steoro_length = 0;
			if (cp0_half_width >= cp1_half_width) {
				top_steoro_length = cp0_half_width
						- Math.pow(Math.pow(cp0_half_width, 2) - Math.pow(cp1_half_width, 2), 0.5);
			} else {
				top_steoro_length = 0;
			}
			double top_steoro_volumn = (1 / 3.) * CCS1 * top_steoro_length;
			double middle_cylinder_volumn = CCS1 * 2 * (cp0_half_width - top_steoro_length);
			inter_volumn = middle_cylinder_volumn + 2 * top_steoro_volumn;
			System.out.println("intersection_volumn_cal: Vertical intersected, c0 is intercepted by c1!");
			return inter_volumn;

		} else {
			// intersected
			double top_steoro_length = 2 * cp1_half_width / Math.tan(Math.toRadians(ang));
			double top_steoro_volumn = (1 / 3.) * CCS1 * top_steoro_length;
			double inter_length = 2 * cp0_half_width / Math.sin(Math.toRadians(ang));
			double middle_cylinder_volumn = CCS1 * (inter_length - top_steoro_length);
			inter_volumn = middle_cylinder_volumn + 2 * top_steoro_volumn;
			System.out.println("intersection_volumn_cal: Intersected angle: " + ang + ", c0 is intercepted by c1!");
			return inter_volumn;

		}
	}

	public double intersection_percentage_cal(double interVolumn, double CCS, double distance) {
		double inter_percent = interVolumn / (CCS * distance);
		return inter_percent;
	}

	public double intersection_2D_size_esti(double ang, contrail_parameters cp0, double sigmaV0, double sigmaH0,
			double CCS0, contrail_parameters cp1, double sigmaV1, double sigmaH1, double CCS1) {

		// Earth radius
		final double earth_r = 6371000.0;

		double cp0_half_height = sigmaV0 * 2.2 / 2;
		double cp0_half_width = sigmaH0 * 2.2 / 2;
		double cp1_half_height = sigmaV1 * 2.2 / 2;
		double cp1_half_width = sigmaH1 * 2.2 / 2;

		// Calculate horizontal distance:
		// Converting latitude and longitude to radians
		double d_lat = Math.toRadians(cp1.lat - cp0.lat);
		double d_lon = Math.toRadians(cp1.lon - cp0.lon);
		double lat0 = Math.toRadians(cp0.lat);
		double lat1 = Math.toRadians(cp1.lat);

		// Applying the Haversine formula
		double a = Math.sin(d_lat / 2) * Math.sin(d_lat / 2)
				+ Math.sin(d_lon / 2) * Math.sin(d_lon / 2) * Math.cos(lat0) * Math.cos(lat1);
		double dis_h = 2 * earth_r * Math.asin(Math.sqrt(a)); // Return distance (meters)

		// System.out.println("dis_h in 2D esti: " +dis_h);

		// Calculate vertical distance:
		// double dis_v = Math.abs(cp0.altitude - cp1.altitude);

		double inter_2D_size = 0.;

		if (ang == 0) {
//			if (c0_full_cover_c1) {
//				inter_2D_size = CCS1;
//				System.out.println("c0_full_cover_c1 is true, parallel, c0 fully covers c1!");
//				return inter_2D_size;
//			} else if (c0_partly_cover_c1) {
//				
//				double inter_h = (cp0_half_width+cp1_half_width-dis_h)/2.;
//				double inter_v =  cp0_half_height*Math.pow(1-Math.pow((cp0_half_width-inter_h)/(cp0_half_width),2), 0.5);
//				inter_2D_size = Math.PI*inter_h*inter_v;
//				if (inter_2D_size >= CCS1) {
//					inter_2D_size = CCS1;
//				}
//				
//				System.out.println("c0_partly_cover_c1 is true, parallel, c0 partly covers c1!");
//				return inter_2D_size;
//			} else {
			if (dis_h + cp1_half_width <= cp0_half_width) {
				// full covered
				inter_2D_size = CCS1;
				System.out.println("intersection_2D_size_esti: Parallel, c0 fully covers c1!");
				return inter_2D_size;
			} else if ((dis_h + cp1_half_width > cp0_half_width) && (dis_h <= cp0_half_width + cp1_half_width)) {
				// partly covered
				double inter_h = (cp0_half_width + cp1_half_width - dis_h) / 2.;
				double inter_v = cp0_half_height
						* Math.pow(1 - Math.pow((cp0_half_width - inter_h) / (cp0_half_width), 2), 0.5);
				inter_2D_size = Math.PI * inter_h * inter_v;
				if (inter_2D_size >= CCS1) {
					inter_2D_size = CCS1;
				}
				System.out.println("intersection_2D_size_esti: Parallel, c0 partly covers c1!");
				return inter_2D_size;
			} else {
				inter_2D_size = 0;
				System.out.println("intersection_2D_size_esti: Parallel, c0 and c1 have separated already!");
				return inter_2D_size;
			}
//			}				

		} else if (ang == 90.0) {
			// intersected
			double top_steoro_length = 0;
			if (cp0_half_width >= cp1_half_width) {
				top_steoro_length = cp0_half_width
						- Math.pow(Math.pow(cp0_half_width, 2) - Math.pow(cp1_half_height, 2), 0.5);
			} else {
				top_steoro_length = 0;
			}
			double top_steoro_volumn = (1 / 3.) * CCS1 * top_steoro_length;
			double middle_cylinder_volumn = CCS1 * 2 * (cp0_half_width - top_steoro_length);
			double inter_volumn = middle_cylinder_volumn + 2 * top_steoro_volumn;
			System.out.println("intersection_2D_size_esti: Vertical intersected, c0 is intercepted by c1!");
			inter_2D_size = inter_volumn / (2 * cp0_half_width);
			return inter_2D_size;

		} else {
			// intersected
			double top_steoro_length = 2 * cp1_half_width / Math.tan(Math.toRadians(ang));
			double top_steoro_volumn = (1 / 3.) * CCS1 * top_steoro_length;
			double inter_length = 2 * cp0_half_width / Math.sin(Math.toRadians(ang));
			double middle_cylinder_volumn = CCS1 * (inter_length - top_steoro_length);
			double inter_volumn = middle_cylinder_volumn + 2 * top_steoro_volumn;
			System.out.println("intersection_2D_size_esti: Intersected angle: " + ang + ", c0 is intercepted by c1!");
			inter_2D_size = inter_volumn / inter_length;
			return inter_2D_size;

		}
	}

	public double ContrailCrossSection_inter_cal_v1(double t0, double t1, double ang, double dis_between,
			contrail_parameters cp0, double sigmaV0, double sigmaH0, contrail_parameters cp1, double sigmaV1,
			double sigmaH1) {

		double A0 = ContrailCrossSection(t0, cp0);
		double A1 = ContrailCrossSection(t1, cp1);

		double cp0_half_height = sigmaV0 * 2.2 / 2;
		double cp0_half_width = sigmaH0 * 2.2 / 2;
		double cp1_half_height = sigmaV1 * 2.2 / 2;
		double cp1_half_width = sigmaH1 * 2.2 / 2;

		double ccs_inter = 0;

		if (ang == 0) {
			if (dis_between + cp1_half_width <= cp0_half_width) {
				// full covered
				ccs_inter = A1;
				System.out.println("CCS_inter_cal: Parallel, c0 fully covers c1, CCS total from c1!");

			} else if ((dis_between + cp1_half_width > cp0_half_width)
					&& (dis_between <= cp0_half_width + cp1_half_width)) {
				// partly covered
				double inter_h = (cp0_half_width + cp1_half_width - dis_between) / 2.;
				double inter_v = cp0_half_height
						* Math.pow(1 - Math.pow((cp0_half_width - inter_h) / (cp0_half_width), 2), 0.5);
				ccs_inter = Math.PI * inter_h * inter_v;
				if (ccs_inter >= A1) {
					ccs_inter = A1;
				}
				System.out.println("CCS_inter_cal: Parallel, c0 partly covers c1, CCS from c0 and c1 intersection!");

			} else {
				ccs_inter = 0;
				System.out.println("CCS_inter_cal: Parallel, c0 and c1 have separated already!");
			}

		} else {
			ccs_inter = A1;
			System.out.println("CCS_inter_cal: Intersected with angle " + ang + ", c0 is intercepted by c1!");
		}

		return ccs_inter;
	}

	public double ContrailCrossSection_inter_cal_v2(double t0, double t1, double ang, double dis_between,
			contrail_parameters cp0, double sigmaV0, double sigmaH0, contrail_parameters cp1, double sigmaV1,
			double sigmaH1) {

		double A0 = ContrailCrossSection(t0, cp0);
		double A1 = ContrailCrossSection(t1, cp1);

		double cp0_half_height = sigmaV0 * 2.2 / 2;
		double cp0_half_width = sigmaH0 * 2.2 / 2;
		double cp1_half_height = sigmaV1 * 2.2 / 2;
		double cp1_half_width = sigmaH1 * 2.2 / 2;

		double ccs_inter = 0;

		if (ang == 0) {
			if (dis_between + cp1_half_width <= cp0_half_width) {
				// full covered
				ccs_inter = A1;
				System.out.println("CCS_inter_cal: Parallel, c0 fully covers c1, CCS total from c1!");

			} else if ((dis_between + cp1_half_width > cp0_half_width)
					&& (dis_between <= cp0_half_width + cp1_half_width)) {
				// partly covered
				double inter_h = (cp0_half_width + cp1_half_width - dis_between) / 2.;
				double inter_v = cp0_half_height
						* Math.pow(1 - Math.pow((cp0_half_width - inter_h) / (cp0_half_width), 2), 0.5);
				ccs_inter = Math.PI * inter_h * inter_v;
				if (ccs_inter >= A1) {
					ccs_inter = A1;
				}
				System.out.println("CCS_inter_cal: Parallel, c0 partly covers c1, CCS from c0 and c1 intersection!");

			} else {
				ccs_inter = 0;
				System.out.println("CCS_inter_cal: Parallel, c0 and c1 have separated already!");
			}

		} else if (ang == 90.0) {
			// intersected
			double top_steoro_length = 0;
			if (cp0_half_width >= cp1_half_width) {
				top_steoro_length = cp0_half_width
						- Math.pow(Math.pow(cp0_half_width, 2) - Math.pow(cp1_half_height, 2), 0.5);
			} else {
				top_steoro_length = 0;
			}
			double top_steoro_volumn = (1 / 3.) * A1 * top_steoro_length;
			double middle_cylinder_volumn = A1 * 2 * (cp0_half_width - top_steoro_length);
			double inter_volumn = middle_cylinder_volumn + 2 * top_steoro_volumn;
			System.out.println("intersection_2D_size_esti: Vertical intersected, c0 is intercepted by c1!");
			ccs_inter = inter_volumn / (2 * cp0_half_width);

		} else {
			// intersected
			double top_steoro_length = 2 * cp1_half_width / Math.tan(Math.toRadians(ang));
			double top_steoro_volumn = (1 / 3.) * A1 * top_steoro_length;
			double inter_length = 2 * cp0_half_width / Math.sin(Math.toRadians(ang));
			double middle_cylinder_volumn = A1 * (inter_length - top_steoro_length);
			double inter_volumn = middle_cylinder_volumn + 2 * top_steoro_volumn;
			System.out.println("intersection_2D_size_esti: Intersected angle: " + ang + ", c0 is intercepted by c1!");
			ccs_inter = inter_volumn / inter_length;

		}

		return ccs_inter;
	}

	public double intersection_state_indicator_decide(double ang, double dis_between, contrail_parameters cp0,
			double sigmaV0, double sigmaH0, contrail_parameters cp1, double sigmaV1, double sigmaH1) {

		double cp0_half_height = sigmaV0 * 2.2 / 2;
		double cp0_half_width = sigmaH0 * 2.2 / 2;
		double cp1_half_height = sigmaV1 * 2.2 / 2;
		double cp1_half_width = sigmaH1 * 2.2 / 2;

		double inter_state_indicator = 0;

		if (ang == 0) {
			if (dis_between + cp1_half_width <= cp0_half_width) {
				// full covered
				inter_state_indicator = 0;
				System.out.println("intersection_state_indicator: " + inter_state_indicator
						+ ", meaning parallel that c0 fully covers c1.");

			} else if ((dis_between + cp1_half_width > cp0_half_width)
					&& (dis_between <= cp0_half_width + cp1_half_width)) {
				// partly covered
				inter_state_indicator = 1;
				System.out.println("intersection_state_indicator: " + inter_state_indicator
						+ ", meaning parallel but c0 partly covers c1.");

			} else {
				inter_state_indicator = 2;
				System.out.println("intersection_state_indicator: " + inter_state_indicator
						+ ", meaning parallel but c0 and c1 have separated totally.");
			}

		} else {
			inter_state_indicator = 3;
			System.err.println("intersection_state_indicator: " + inter_state_indicator
					+ ", Attention that angle should not be " + ang + ",");
			System.err.println("which can not be simulated in this code version.");

		}

		return inter_state_indicator;
	}

	// if intersected, the two contrails should be considered in the same level
	public double contrail_centrals_distence_cal(contrail_parameters cp0, contrail_parameters cp1) {
		// Earth radius
		final double earth_r = 6371000.0;

		// Calculate horizontal distance:
		// Converting latitude and longitude to radians
		double d_lat = Math.toRadians(cp1.lat - cp0.lat);
		double d_lon = Math.toRadians(cp1.lon - cp0.lon);
		double lat0 = Math.toRadians(cp0.lat);
		double lat1 = Math.toRadians(cp1.lat);

		// Applying the Haversine formula
		double a = Math.sin(d_lat / 2) * Math.sin(d_lat / 2)
				+ Math.sin(d_lon / 2) * Math.sin(d_lon / 2) * Math.cos(lat0) * Math.cos(lat1);
		double dis_h = 2 * earth_r * Math.asin(Math.sqrt(a)); // Return distance (meters)

		return dis_h;
	}

	public double contrail0_exposed_surface_percentage_cal(double ang, double dis_between, double sigmaV0,
			double sigmaH0, double CCS0, double sigmaV1, double sigmaH1, double CCS1, double distance) {

		double cp0_half_height = sigmaV0 * 2.2 / 2;
		double cp0_half_width = sigmaH0 * 2.2 / 2;
		double cp1_half_height = sigmaV1 * 2.2 / 2;
		double cp1_half_width = sigmaH1 * 2.2 / 2;

		double c0_exposed_surface_percent = 0;
		if (ang == 0) {
			if (dis_between + cp1_half_width <= cp0_half_width) {
				// full covered
				c0_exposed_surface_percent = 1;
				System.out.println("c0 esp cal: Parallel, c0 fully covers c1!");

			} else if ((dis_between + cp1_half_width > cp0_half_width)
					&& (dis_between <= cp0_half_width + cp1_half_width)) {
				// partly covered
				double inter_h = (cp0_half_width + cp1_half_width - dis_between) / 2.;
				double inter_v = cp0_half_height
						* Math.pow(1 - Math.pow((cp0_half_width - inter_h) / (cp0_half_width), 2), 0.5);
				double contrail0_covered_side = (2 * Math.PI * inter_v + 4 * (inter_h - inter_v)) / 2.;
				double contrail0_perimeter = 2 * Math.PI * cp0_half_height + 4 * (cp0_half_width - cp0_half_height);
				c0_exposed_surface_percent = 1 - contrail0_covered_side / contrail0_perimeter;
				System.out.println("c0 esp cal: Parallel, c0 partly covers c1!");
			} else {
				c0_exposed_surface_percent = 1;
				System.out.println("c0 esp cal: Parallel, c0 and c1 have separated already!");
			}

		} else {
			// intersected
			double contrail0_total_surface = (2 * Math.PI * cp0_half_height + 4 * (cp0_half_width - cp0_half_height))
					* distance;
			// double contrail0_intercepted_surface =
			// 2*Math.PI*cp1_half_height*cp1_half_width/Math.sin(Math.toRadians(ang));
			double contrail0_intercepted_surface = 2 * CCS1 / Math.sin(Math.toRadians(ang));
			c0_exposed_surface_percent = (contrail0_total_surface - contrail0_intercepted_surface)
					/ contrail0_total_surface;
			System.out.println("c0 esp cal: Intersected with angle " + ang + ", c0 is intercepted by c1!");
		}

		if (c0_exposed_surface_percent <= 0) {
			System.out.println("Attention from c0 esp cal: current c0_exposed_surface_percent "
					+ c0_exposed_surface_percent + " impossible!!!");
		}

		return c0_exposed_surface_percent;
	}

	public double contrail1_exposed_surface_percentage_cal(double ang, double dis_between, double sigmaV0,
			double sigmaH0, double CCS0, double sigmaV1, double sigmaH1, double CCS1, double distance) {

		double cp0_half_height = sigmaV0 * 2.2 / 2;
		double cp0_half_width = sigmaH0 * 2.2 / 2;
		double cp1_half_height = sigmaV1 * 2.2 / 2;
		double cp1_half_width = sigmaH1 * 2.2 / 2;

		double c1_exposed_surface_percent = 0;

		if (ang == 0) {
			if (dis_between + cp1_half_width <= cp0_half_width) {
				// full covered by c0
				c1_exposed_surface_percent = 0;
				System.out.println("c1 esp cal: Parallel, c0 fully covers c1!");
			} else if ((dis_between + cp1_half_width > cp0_half_width)
					&& (dis_between <= cp0_half_width + cp1_half_width)) {
				// partly covered
				double inter_h = (cp0_half_width + cp1_half_width - dis_between) / 2.;
				double inter_v = cp0_half_height
						* Math.pow(1 - Math.pow((cp0_half_width - inter_h) / (cp0_half_width), 2), 0.5);
				double contrail1_covered_side = (2 * Math.PI * inter_v + 4 * (inter_h - inter_v)) / 2.;
				double contrail1_perimeter = 2 * Math.PI * cp1_half_height + 4 * (cp1_half_width - cp1_half_height);
				c1_exposed_surface_percent = 1 - contrail1_covered_side / contrail1_perimeter;
				System.out.println("c1 esp cal: Parallel, c0 partly covers c1!");
			} else {
				c1_exposed_surface_percent = 1;
				System.out.println("c1 esp cal: Parallel, c0 and c1 have separated already!");
			}

		} else {
			double contrail1_total_surface = (2 * Math.PI * cp1_half_height + 4 * (cp1_half_width - cp1_half_height))
					* distance;
			double contrail1_intercepted_length = 2 * cp0_half_width / Math.sin(Math.toRadians(ang));
			double contrail1_intercepted_surface = (2 * Math.PI * cp1_half_height
					+ 4 * (cp1_half_width - cp1_half_height)) * contrail1_intercepted_length;
			c1_exposed_surface_percent = (contrail1_total_surface - contrail1_intercepted_surface)
					/ contrail1_total_surface;
			System.out.println("c1 esp cal: Intersected with angle " + ang + ", c0 is intercepted by c1!");
		}

		if (c1_exposed_surface_percent <= 0) {
			System.out.println("Attention from c1 esp cal, current c1_exposed_surface_percent: "
					+ c1_exposed_surface_percent + ",");
			System.out.println("which should be not <=0, then increasing contrail distance!!!!");
			c1_exposed_surface_percent = 0;
		}

		return c1_exposed_surface_percent;
	}

	public void ContrailEvolution() {

		/////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////
		// initializing the parameters of old contrail cp_0, put them to the new class
		// cp_inCE_0
		contrail_parameters cp_inCE_0 = this.cp_0;
		// Simulation of Contrail Evolution
		// mH2O...emitted water mass [kg]
		// distance...distance from last flight profile step [m]
		// Parameters which describe the contrail shape:
		// Dv, Dh, Ds... vertical, horizontal and sheared diffusion [m2/s]
		// vo...mean value of wind velocity at flight level [m/s]
		// s....constant wind shear

		// === 1st step: calculate initial number of ice particles ===
		// Sussmann, J. Geophys. Res., 106(D5), 4899 (2001)
		double Nice_0 = 1e15 * this.fuelflow_0 * this.distance_0 / this.TAS_0; // number
		double Nice_rest_0 = 0;
		System.out.println("number of ice particles: " + Nice_0);

		// Check the contrail 0 evolution state and if its done, stop data record:
		boolean evo_state_0 = false;
		boolean evo_done_0 = false;

		// ========== Beginning of contrail lifetime ==========
		double t_0 = 0; // at t=0
		double z_0 = 0; // at flight level
		// FIXME
		cp_inCE_0.T = cp_inCE_0.Tflight; // contrail temperature is equal to ambient temperature at flight level
		// delta t for time integration steps
		double delta_t_0 = this.dt_init_0; // [s]

		// === calculate ice water content ===
		double CCS_0 = ContrailCrossSection(0, cp_inCE_0);
		double CCSc_0 = 0;
		double CCS_rest_0 = CCS_0;
		double V_rest_0 = CCS_rest_0 * this.distance_0;
		// FIXME
		// SATURATION PRESSUREs FROM SONNTAG, Temp IN K, PSAT IN PA - METEOROL. Z., 3
		// (1994) 51-66.
		double eStar_0 = Math.exp(-6096.9385 / cp_inCE_0.T + 16.635794 - 2.711193E-2 * cp_inCE_0.T
				+ 1.673952E-5 * cp_inCE_0.T * cp_inCE_0.T + 2.433502 * Math.log(cp_inCE_0.T));
		double eStarIce_0 = Math.exp(-6024.5282 / cp_inCE_0.T + 24.721994 + 1.0613868E-2 * cp_inCE_0.T
				- 1.3198825E-5 * cp_inCE_0.T * cp_inCE_0.T - 0.49382577 * Math.log(cp_inCE_0.T));

		double RHice_0 = gribdata.getHumidityFromAltitude(cp_inCE_0.lat, cp_inCE_0.lon, cp_inCE_0.altitude + z_0)
				* eStar_0 / eStarIce_0;

		// inter_state_indicator initial setup, value = 100.0
		double s2_intersection_state_indicator = 100.;

		if (RHice_0 < 100.) {
			System.out.println();
			System.out.println("RHice = " + RHice_0 + " < 100 ... Relative Humidity of ambient air is only "
					+ Math.round(
							(gribdata.getHumidityFromAltitude(cp_inCE_0.lat, cp_inCE_0.lon, cp_inCE_0.altitude + z_0)))
					+ "% ... no contrail evolves!");
			System.exit(1);
		}

		double c0_esp = 1.;
		double delta_e_0 = eStarIce_0 * (RHice_0 / 100. - 1.);
		double IWCs_0 = 0.21667 * delta_e_0 / cp_inCE_0.T;
		double IWC_0 = 1.24 * this.fuelflow_0 / (this.TAS_0 * CCS_0) + IWCs_0;
		// double
		// testpressure=gribdata.getPressureFromAltitude(this.cp.lat,this.cp.lon,10500);

		double sigmah_0 = cp_inCE_0.sigma0h;
		double sigmav_0 = cp_inCE_0.sigma0v;
		double sigmas_0 = cp_inCE_0.s * cp_inCE_0.Dv * Math.pow(t_0, 2)
				+ (2 * cp_inCE_0.Ds + cp_inCE_0.s * Math.pow(cp_inCE_0.sigma0v, 2)) * t_0;

		double dice_0 = Math.pow(IWC_0 * CCS_0 * this.distance_0 / (Nice_0 * 917) * 6. / Math.PI, 1. / 3.);
		double rext_0 = dice_0 / 2. * 1.e6;
		try {
			this.fwDataOut_0 = new FileWriter(fn_output_0);
		} catch (IOException ex) {
			System.err.println("Fehler beim Anlegen der Datei");
		}
		;

		try {
			this.fwDataOut_0.write(
					"t, z, eps, Dh, Dv, sigmah, sigmav, sigmas, rice, Nice, g_t, Qabs_t, Qsca_t, g_s, Qabs_s, Qsca_s, state_indicator, c0_esp, CCS, CCS_rest, V, V_rest, IWC, RHice, Re, v_sed, windshear_s, wind_u, wind_v, wind_w, lat, lon, altitude\n");
			rext_0 = dice_0 / 2. * 1.e6;
			this.fwDataOut_0
					.write(t_0 + ", " + z_0 + ", " + cp_inCE_0.eps + ", " + cp_inCE_0.Dh + ", " + cp_inCE_0.Dv + ", "
							+ sigmah_0 + ", " + sigmav_0 + ", " + sigmas_0 + ", " + (dice_0 / 2) + ", " + Nice_0 + ", "
							+ Ext_0.Calc_g_terr(rext_0, lterr_0) + ", " + Ext_0.Calc_Qabs_terr(rext_0, lterr_0) + ", "
							+ Ext_0.Calc_Qsca_terr(rext_0, lterr_0) + ", " + Ext_0.Calc_g_sol(rext_0, bsol_0) + ", "
							+ Ext_0.Calc_Qabs_sol(rext_0, bsol_0) + ", " + Ext_0.Calc_Qsca_sol(rext_0, bsol_0) + ", "
							+ +s2_intersection_state_indicator + ", " + ((int) (c0_esp * 10000) / 10000.) + ", " + CCS_0
							+ ", " + CCS_rest_0 + ", " + (CCS_0 * this.distance_0) + ", " + V_rest_0 + ", " + IWC_0
							+ ", " + RHice_0 + ", "
							+ Re(dice_0 / 2., SedimentationSpeed(dice_0, cp_inCE_0),
									kinvisc_Suth(gribdata.getPressureFromAltitude(cp_inCE_0.lat, cp_inCE_0.lon,
											cp_inCE_0.altitude), cp_inCE_0.T))
							+ ", " + SedimentationSpeed(dice_0, cp_inCE_0) + ", " + cp_inCE_0.s + ", "
							+ gribdata.getWindUFromAltitude(cp_inCE_0.lat, cp_inCE_0.lon, cp_inCE_0.altitude) + ", "
							+ gribdata.getWindVFromAltitude(cp_inCE_0.lat, cp_inCE_0.lon, cp_inCE_0.altitude) + ", "
							+ cp_inCE_0.vz + ", " + ((int) (cp_inCE_0.lat * 10000) / 10000.) + ", "
							+ ((int) (cp_inCE_0.lon * 10000) / 10000.) + ", " + Math.round(cp_inCE_0.altitude) + "\n");

		} catch (IOException ex) {
			System.err.println("Fehler beim Schreiben");
		}
		;

		/////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////
		// initializing the parameters of new contrail cp_1, put them to the new class
		// cp_inCE_1
		contrail_parameters cp_inCE_1 = this.cp_1;
		// Simulation of Contrail Evolution
		// mH2O...emitted water mass [kg]
		// distance...distance from last flight profile step [m]
		// Parameters which describe the contrail shape:
		// Dv, Dh, Ds... vertical, horizontal and sheared diffusion [m2/s]
		// vo...mean value of wind velocity at flight level [m/s]
		// s....constant wind shear

		// === 1st step: calculate initial number of ice particles ===
		// Sussmann, J. Geophys. Res., 106(D5), 4899 (2001)
		double Nice_1 = 1e15 * this.fuelflow_1 * this.distance_1 / this.TAS_1; // number
		double Nice_rest_1 = 0;
		System.out.println("number of ice particles: " + Nice_1);

		// Check the contrail 1 evolution state and if its done, stop data record:
		boolean evo_state_1 = false;
		boolean evo_done_1 = false;

		// ========== Beginning of contrail lifetime ==========
		double t_1 = 0; // at t=0
		double z_1 = 0; // at flight level
		// FIXME
		cp_inCE_1.T = cp_inCE_1.Tflight; // contrail temperature is equal to ambient temperature at flight level
		// delta t for time integration steps
		double delta_t_1 = this.dt_init_1; // [s]

		// === calculate ice water content ===
		double CCS_1 = ContrailCrossSection(0, cp_inCE_1);
		double CCSc_1 = 0;
		double CCS_rest_1 = CCS_1;
		double V_rest_1 = CCS_rest_1 * this.distance_1;
		// FIXME
		// SATURATION PRESSUREs FROM SONNTAG, Temp IN K, PSAT IN PA - METEOROL. Z., 3
		// (1994) 51-66.
		double eStar_1 = Math.exp(-6096.9385 / cp_inCE_1.T + 16.635794 - 2.711193E-2 * cp_inCE_1.T
				+ 1.673952E-5 * cp_inCE_1.T * cp_inCE_1.T + 2.433502 * Math.log(cp_inCE_1.T));
		double eStarIce_1 = Math.exp(-6024.5282 / cp_inCE_1.T + 24.721994 + 1.0613868E-2 * cp_inCE_1.T
				- 1.3198825E-5 * cp_inCE_1.T * cp_inCE_1.T - 0.49382577 * Math.log(cp_inCE_1.T));

		double RHice_1 = gribdata.getHumidityFromAltitude(cp_inCE_1.lat, cp_inCE_1.lon, cp_inCE_1.altitude + z_1)
				* eStar_1 / eStarIce_1;

		if (RHice_1 < 100.) {
			System.out.println();
			System.out.println("RHice = " + RHice_1 + " < 100 ... Relative Humidity of ambient air is only "
					+ Math.round(
							(gribdata.getHumidityFromAltitude(cp_inCE_1.lat, cp_inCE_1.lon, cp_inCE_1.altitude + z_1)))
					+ "% ... no contrail evolves!");
			System.exit(1);
		}

		double c1_esp = 1.;
		// according to angle, setting the delta_e_1
		double delta_e_1 = eStarIce_1 * (RHice_1 / 100. - 1.);

		double IWCs_1 = 0.21667 * delta_e_1 / cp_inCE_1.T;
		double IWC_1 = 1.24 * this.fuelflow_1 / (this.TAS_1 * CCS_1) + IWCs_1;
		// double
		// testpressure=gribdata.getPressureFromAltitude(this.cp.lat,this.cp.lon,10500);

		double sigmah_1 = cp_inCE_1.sigma0h;
		double sigmav_1 = cp_inCE_1.sigma0v;
		double sigmas_1 = cp_inCE_1.s * cp_inCE_1.Dv * Math.pow(t_1, 2)
				+ (2 * cp_inCE_1.Ds + cp_inCE_1.s * Math.pow(cp_inCE_1.sigma0v, 2)) * t_1;

		double dice_1 = Math.pow(IWC_1 * CCS_1 * this.distance_1 / (Nice_1 * 917) * 6. / Math.PI, 1. / 3.);
		double rext_1 = dice_1 / 2. * 1.e6;

		double distance_between_two_contrails = contrail_centrals_distence_cal(cp_inCE_0, cp_inCE_1);

		try {
			this.fwDataOut_1 = new FileWriter(fn_output_1);
		} catch (IOException ex) {
			System.err.println("Fehler beim Anlegen der Datei");
		}
		;

		try {
			this.fwDataOut_1.write(
					"angle, time_diff, t, z, eps, Dh, Dv, sigmah, sigmav, sigmas, rice, Nice, g_t, Qabs_t, Qsca_t, g_s, Qabs_s, Qsca_s, distance_between_contrails, state_indicator, c1_esp, CCS, CCS_rest, V, V_rest, IWC, RHice, Re, v_sed, windshear_s, wind_u, wind_v, wind_w, lat, lon, altitude\n");
			rext_1 = dice_1 / 2. * 1.e6;
			this.fwDataOut_1
					.write((this.angle) + ", " + time_diff + ", " + t_1 + ", " + z_1 + ", " + cp_inCE_1.eps + ", "
							+ cp_inCE_1.Dh + ", " + cp_inCE_1.Dv + ", " + sigmah_1 + ", " + sigmav_1 + ", " + sigmas_1
							+ ", " + (dice_1 / 2) + ", " + Nice_1 + ", " + Ext_1.Calc_g_terr(rext_1, lterr_1) + ", "
							+ Ext_1.Calc_Qabs_terr(rext_1, lterr_1) + ", " + Ext_1.Calc_Qsca_terr(rext_1, lterr_1)
							+ ", " + Ext_1.Calc_g_sol(rext_1, bsol_1) + ", " + Ext_1.Calc_Qabs_sol(rext_1, bsol_1)
							+ ", " + Ext_1.Calc_Qsca_sol(rext_1, bsol_1) + ", " + distance_between_two_contrails + ", "
							+ s2_intersection_state_indicator + ", " + ((int) (c1_esp * 10000) / 10000.) + ", " + CCS_1
							+ ", " + CCS_rest_1 + ", " + (CCS_1 * this.distance_1) + ", " + V_rest_1 + ", " + +IWC_1
							+ ", " + RHice_1 + ", "
							+ Re(dice_1 / 2., SedimentationSpeed(dice_1, cp_inCE_1),
									kinvisc_Suth(gribdata.getPressureFromAltitude(cp_inCE_1.lat, cp_inCE_1.lon,
											cp_inCE_1.altitude), cp_inCE_1.T))
							+ ", " + SedimentationSpeed(dice_1, cp_inCE_1) + ", " + cp_inCE_1.s + ", "
							+ gribdata.getWindUFromAltitude(cp_inCE_1.lat, cp_inCE_1.lon, cp_inCE_1.altitude) + ", "
							+ gribdata.getWindVFromAltitude(cp_inCE_1.lat, cp_inCE_1.lon, cp_inCE_1.altitude) + ", "
							+ cp_inCE_1.vz + ", " + ((int) (cp_inCE_1.lat * 10000) / 10000.) + ", "
							+ ((int) (cp_inCE_1.lon * 10000) / 10000.) + ", " + Math.round(cp_inCE_1.altitude) + "\n");

		} catch (IOException ex) {
			System.err.println("Fehler beim Schreiben");
		}
		;

		System.out.println("contrail_0 and contrail_1 initialitation finished!");
		/////////////////////////////////////////////////////////////////////////////
		///////////// *contrail_0 and contrail_1 initialitation finished!*////////////
		/////////////////////////////////////////////////////////////////////////////

		double eStarIce_old_0 = 0.;
		double eStarIce_old_1 = 0.;

		// c_inter parameter initialization
		boolean cp_inter_initialization_flag = false;
		boolean intersection_check_flag = false;
		boolean intersection_done_flag = false;

		boolean c1_initial_full_inside_c0_flag = false;
		boolean c1_initial_partly_inside_c0_flag = false;
		boolean evo_state_inter = false;
		boolean evo_done_inter = false;
		double eStarIce_old_inter = 0.;

		this.cp_inter = new contrail_parameters();
		double dice_inter = 0;
		double t_inter = 0;
		double z_inter = 0;
		double delta_z_inter = 0;

		double sigmah_inter = 0;
		double sigmav_inter = 0;
		double sigmas_inter = 0;
		double volumn_inter = 0;
		double Nice_inter = 0;
		double CCS_inter = 0;
		double CCSc_inter = 0;
		double CCS_inter_esti = 0;
		double IWC_inter = 0;
		double lterr_inter = 0;
		int bsol_inter = 0;
		double eStar_inter = 0;
		double eStarIce_inter = 0;
		double RHice_inter = 0;
		double delta_e_inter = 0;
		double IWCs_inter = 0;
		double delta_t_inter = 0;

		// t_universal setup
		double t_universal = 1;
		double delta_t_universal = 0;

		do {

			if ((evo_state_0 == false) && (evo_done_0 == false)) {
				// calculate ice particle diameter
				if (intersection_check_flag == true) {

					System.out.println("================================================");
					if (intersection_done_flag == false) {
						Nice_rest_0 = Nice_0 * (1 - intersection_percentage_cal(volumn_inter, CCS_0, this.distance_0));
					} else {
						Nice_rest_0 = Nice_rest_0;
					}

					dice_0 = Math.pow(IWC_0 * CCS_0 * this.distance_0
							* (1 - intersection_percentage_cal(volumn_inter, CCS_0, this.distance_0))
							/ (Nice_rest_0 * 917) * 6 / Math.PI, 1 / 3.);

					System.out.println("Nice_rest_0 = " + Nice_rest_0);
					System.out.println("intersection_percentage_cal = "
							+ intersection_percentage_cal(volumn_inter, CCS_0, this.distance_0));
				} else {
					Nice_rest_0 = Nice_0;
					dice_0 = Math.pow(IWC_0 * CCS_0 * this.distance_0 / (Nice_0 * 917) * 6 / Math.PI, 1 / 3.);
				}

				// dice = Math.cbrt(IWC*CCS*this.distance/(Nice*917)*6/Math.PI);

				// next time step update:
				t_0 = t_0 + delta_t_0;

				// cp.Dht = cp.Dht + cp.Dh*t;
				// contrail geometry

				// Contrail horizontal movement
				// Calculates displaced contrail position (lat, lon)
				// this.cp.lat = this.cp.lat + (-16.0*(delta_t/(60000*1.852)));

				cp_inCE_0.lat = cp_inCE_0.lat
						+ (gribdata.getWindVFromAltitude(cp_inCE_0.lat, cp_inCE_0.lon, cp_inCE_0.altitude)
								* (delta_t_0 / (60000 * 1.852)));
				cp_inCE_0.lon = cp_inCE_0.lon
						+ (gribdata.getWindUFromAltitude(cp_inCE_0.lat, cp_inCE_0.lon, cp_inCE_0.altitude)
								* (delta_t_0 / (60000 * 1.852 * Math.cos(Math.toRadians(cp_inCE_0.lat))))); // 60nm per
																											// degree
				// at lat = 0°
				// and 1852
				// metres per nm

				// Contrail vertical movement
				// Let contrail drift down
				cp_inCE_0.vz = vz(cp_inCE_0);
				double delta_z_0 = (cp_inCE_0.vz - SedimentationSpeed(dice_0, cp_inCE_0)) * delta_t_0; // negative
																										// delta_z
																										// means sinking
																										// not
																										// thinking
				z_0 = z_0 + delta_z_0; // z ist die Summe der delta z je Durchlauf
				cp_inCE_0.altitude = cp_inCE_0.altitude + delta_z_0;

//				System.out.println("contrail 0 lat:" + ((int) (cp_inCE_0.lat * 10000) / 10000.));
//				System.out.println("contrail 0 lon:" + cp_inCE_0.lon);
//				System.out.println("contrail 0 alt:" + cp_inCE_0.altitude);

				sigmav_0 = Math.pow(2. * cp_inCE_0.Dv * t_0 + Math.pow(cp_inCE_0.sigma0v, 2), 0.5);
				cp_inCE_0.s = (gribdata.getWindspeedFromAltitude(cp_inCE_0.lat, cp_inCE_0.lon,
						cp_inCE_0.altitude + sigmav_0)
						- gribdata.getWindspeedFromAltitude(cp_inCE_0.lat, cp_inCE_0.lon,
								cp_inCE_0.altitude - sigmav_0))
						/ (2 * sigmav_0); // calculates
											// the
											// change of
											// windspeed
											// between
											// upper and
											// lower
											// boundary
											// of
											// contrail
											// per metre
				sigmah_0 = Math.pow(
						2. / 3. * cp_inCE_0.s * cp_inCE_0.s * cp_inCE_0.Dv * Math.pow(t_0, 3)
								+ (2 * cp_inCE_0.Ds + cp_inCE_0.s * Math.pow(cp_inCE_0.sigma0v, 2)) * cp_inCE_0.s
										* Math.pow(t_0, 2)
								+ 2 * cp_inCE_0.Dh * t_0 + Math.pow(cp_inCE_0.sigma0h, 2),
						0.5);
				sigmas_0 = cp_inCE_0.s * cp_inCE_0.Dv * Math.pow(t_0, 2)
						+ (2 * cp_inCE_0.Ds + cp_inCE_0.s * Math.pow(cp_inCE_0.sigma0v, 2)) * t_0;

//				System.out.println("contrail 0 sigmav_0:" + sigmav_0);
//				System.out.println("contrail 0 sigmah_0:" + sigmah_0);

				// System.out.println("Sigmas in main body: " + sigmas);
				// cp.Dh = 0.1 * cp.uprime * sigmah;

				// eStarIce before adiabatic heating
				eStarIce_old_0 = Math.exp(-6024.5282 / cp_inCE_0.T + 24.721994 + 1.0613868E-2 * cp_inCE_0.T
						- 1.3198825E-5 * cp_inCE_0.T * cp_inCE_0.T - 0.49382577 * Math.log(cp_inCE_0.T));

				// Temperature
				cp_inCE_0.T = (float) (cp_inCE_0.Tflight - AdiabaticCurve(cp_inCE_0) * z_0); // adiabatic heating

				// calculate additional ice water content

				// SATURATION PRESSUREs FROM SONNTAG, Temp IN K, PSAT IN PA - METEOROL. Z., 3
				// (1994) 51-66.
				eStar_0 = Math.exp(-6096.9385 / cp_inCE_0.T + 16.635794 - 2.711193E-2 * cp_inCE_0.T
						+ 1.673952E-5 * cp_inCE_0.T * cp_inCE_0.T + 2.433502 * Math.log(cp_inCE_0.T));
				eStarIce_0 = Math.exp(-6024.5282 / cp_inCE_0.T + 24.721994 + 1.0613868E-2 * cp_inCE_0.T
						- 1.3198825E-5 * cp_inCE_0.T * cp_inCE_0.T - 0.49382577 * Math.log(cp_inCE_0.T));
				RHice_0 = gribdata.getHumidityFromAltitude(cp_inCE_0.lat, cp_inCE_0.lon, cp_inCE_0.altitude) * eStar_0
						/ eStarIce_0;

				CCSc_0 = ContrailCrossSection(t_0, cp_inCE_0);

				if (t_0 <= time_diff) {
					if (RHice_0 < 100) {
						System.out.println(
								"CCS and IWC keep the same with the values in the second last loop, which not updated any more!");
						System.out.println("RHice = " + RHice_0 + " < 100");
						evo_state_0 = true;
					} else {

						// FIXME
						// get on with it
						delta_e_0 = eStarIce_0 * (RHice_0 / 100. - 1.);
						IWCs_0 = 0.21667 * delta_e_0 / cp_inCE_0.T; // Masse des Wasser ueber Saettigung in der Umgebung
																	// (pro
																	// Volumen)
						// calculate new ice water content
						// before IWCs is added, IWC is reduced because of adiabatic heating
						IWC_0 = IWC_0 - (eStarIce_0 - eStarIce_old_0) / (462 * cp_inCE_0.T);

						// current CCSc is the CCSn in last code version, here should not be next, it's
						// current CCS.
						delta_t_0 = 0.1 * 2 * cp_inCE_0.Dh * cp_inCE_0.Dv * t_0;
						IWC_0 = (IWCs_0 * (CCSc_0 - CCS_0) + IWC_0 * CCS_0) / CCSc_0;
						// CCS means the old state.
						CCS_0 = CCSc_0;
						/*
						 * IWC and CCS affect the next loop calculation.
						 */
						evo_state_0 = false;

						CCS_rest_0 = CCS_0;
						V_rest_0 = CCS_rest_0 * this.distance_0;
					}
				}
			}

			if (t_0 > time_diff) {
				if ((evo_state_1 == false) && (evo_done_1 == false)) {
					if (c1_initial_full_inside_c0_flag == false) {
						// start contrail 1 evoluting

						// calculate ice particle diameter
						if (intersection_check_flag == true) {

							// volumn_inter = intersection_volumn_cal(angle, false,
							// c1_partly_inside_c0_flag, cp_inCE_0, sigmav_0, sigmah_0,
							// CCS_0, cp_inCE_1, sigmav_1, sigmah_1, CCS_1, this.distance_1);
							if (intersection_done_flag == false) {
								Nice_rest_1 = Nice_1
										* (1 - intersection_percentage_cal(volumn_inter, CCS_1, this.distance_1));
							} else {
								Nice_rest_1 = Nice_rest_1;
							}

							// System.out.println("See v1 intersection percentage: " +
							// intersection_percentage_cal(volumn_inter,CCS_1,this.distance_1 ));
							dice_1 = Math.pow(IWC_1 * CCS_1 * this.distance_1
									* (1 - intersection_percentage_cal(volumn_inter, CCS_1, this.distance_1))
									/ (Nice_rest_1 * 917) * 6 / Math.PI, 1 / 3.);

						} else {
							Nice_rest_1 = Nice_1;
							dice_1 = Math.pow(IWC_1 * CCS_1 * this.distance_1 / (Nice_1 * 917) * 6 / Math.PI, 1 / 3.);
						}

						System.out.println("=======================");
						System.out.println("Nice_rest_1 = " + Nice_rest_1);
						System.out.println("intersection_percentage_cal = "
								+ intersection_percentage_cal(volumn_inter, CCS_1, this.distance_1));

						// next time step update:
						t_1 = t_1 + delta_t_1;

						// contrail geometry

						// Contrail horizontal movement
						// Calculates displaced contrail position (lat, lon)

						cp_inCE_1.lat = cp_inCE_1.lat
								+ (gribdata.getWindVFromAltitude(cp_inCE_1.lat, cp_inCE_1.lon, cp_inCE_1.altitude)
										* (delta_t_1 / (60000 * 1.852)));
						cp_inCE_1.lon = cp_inCE_1.lon
								+ (gribdata.getWindUFromAltitude(cp_inCE_1.lat, cp_inCE_1.lon, cp_inCE_1.altitude)
										* (delta_t_1 / (60000 * 1.852 * Math.cos(Math.toRadians(cp_inCE_1.lat))))); // 60nm
																													// per
																													// degree
						// at lat = 0°
						// and 1852
						// metres per nm

						// Contrail vertical movement
						// Let contrail drift down
						cp_inCE_1.vz = vz(cp_inCE_1);
						double delta_z_1 = (cp_inCE_1.vz - SedimentationSpeed(dice_1, cp_inCE_1)) * delta_t_1; // negative
																												// delta_z
																												// means
																												// sinking
																												// not
																												// thinking
						z_1 = z_1 + delta_z_1; // z ist die Summe der delta z je Durchlauf
						cp_inCE_1.altitude = cp_inCE_1.altitude + delta_z_1;

						// System.out.println("contrail 1 lat:" + cp_inCE_1.lat);
						// System.out.println("contrail 1 lon:" + cp_inCE_1.lon);
						// System.out.println("contrail 1 alt:" + cp_inCE_1.altitude);

						sigmav_1 = Math.pow(2. * cp_inCE_1.Dv * t_1 + Math.pow(cp_inCE_1.sigma0v, 2), 0.5);
						cp_inCE_1.s = (gribdata.getWindspeedFromAltitude(cp_inCE_1.lat, cp_inCE_1.lon,
								cp_inCE_1.altitude + sigmav_1)
								- gribdata.getWindspeedFromAltitude(cp_inCE_1.lat, cp_inCE_1.lon,
										cp_inCE_1.altitude - sigmav_1))
								/ (2 * sigmav_1); // calculates
													// the
													// change of
													// windspeed
													// between
													// upper and
													// lower
													// boundary
													// of
													// contrail
													// per metre
						sigmah_1 = Math
								.pow(2. / 3. * cp_inCE_1.s * cp_inCE_1.s * cp_inCE_1.Dv * Math.pow(t_1, 3)
										+ (2 * cp_inCE_1.Ds + cp_inCE_1.s * Math.pow(cp_inCE_1.sigma0v, 2))
												* cp_inCE_1.s * Math.pow(t_1, 2)
										+ 2 * cp_inCE_1.Dh * t_1 + Math.pow(cp_inCE_1.sigma0h, 2), 0.5);
						sigmas_1 = cp_inCE_1.s * cp_inCE_1.Dv * Math.pow(t_1, 2)
								+ (2 * cp_inCE_1.Ds + cp_inCE_1.s * Math.pow(cp_inCE_1.sigma0v, 2)) * t_1;

						// System.out.println("contrail 1 sigmav_1:" + sigmav_1);
						// System.out.println("contrail 1 sigmah_1:" + sigmah_1);

						// System.out.println("Sigmas in main body: " + sigmas);
						// cp.Dh = 0.1 * cp.uprime * sigmah;

						// eStarIce before adiabatic heating
						eStarIce_old_1 = Math.exp(-6024.5282 / cp_inCE_1.T + 24.721994 + 1.0613868E-2 * cp_inCE_1.T
								- 1.3198825E-5 * cp_inCE_1.T * cp_inCE_1.T - 0.49382577 * Math.log(cp_inCE_1.T));

						// Temperature
						cp_inCE_1.T = (float) (cp_inCE_1.Tflight - AdiabaticCurve(cp_inCE_1) * z_1); // adiabatic
																										// heating

						// calculate additional ice water content

						// SATURATION PRESSUREs FROM SONNTAG, Temp IN K, PSAT IN PA - METEOROL. Z., 3
						// (1994) 51-66.
						eStar_1 = Math.exp(-6096.9385 / cp_inCE_1.T + 16.635794 - 2.711193E-2 * cp_inCE_1.T
								+ 1.673952E-5 * cp_inCE_1.T * cp_inCE_1.T + 2.433502 * Math.log(cp_inCE_1.T));
						eStarIce_1 = Math.exp(-6024.5282 / cp_inCE_1.T + 24.721994 + 1.0613868E-2 * cp_inCE_1.T
								- 1.3198825E-5 * cp_inCE_1.T * cp_inCE_1.T - 0.49382577 * Math.log(cp_inCE_1.T));
						RHice_1 = gribdata.getHumidityFromAltitude(cp_inCE_1.lat, cp_inCE_1.lon, cp_inCE_1.altitude)
								* eStar_1 / eStarIce_1;
						CCSc_1 = ContrailCrossSection(t_1, cp_inCE_1);

					}
				}

				// if intersected, the two contrails should be in the same level

				distance_between_two_contrails = contrail_centrals_distence_cal(cp_inCE_0, cp_inCE_1);
				System.out.println("distance_between_two_contrails: " + distance_between_two_contrails);
				s2_intersection_state_indicator = intersection_state_indicator_decide(angle,
						distance_between_two_contrails, cp_inCE_0, sigmav_0, sigmah_0, cp_inCE_1, sigmav_1, sigmah_1);

				if ((intersection_check_flag == false)
						&& intersection_check(cp_inCE_0, sigmav_0, sigmah_0, cp_inCE_1, sigmav_1, sigmah_1)
						&& (evo_state_0 == false) && (evo_done_0 == false) && (evo_state_1 == false)
						&& (evo_done_1 == false) && (RHice_0 >= 100) && (RHice_1 >= 100)) {

					// save intersected data into a new output file

					this.cp_inter.lat = (cp_inCE_0.lat + cp_inCE_1.lat) / 2.;
					this.cp_inter.lon = (cp_inCE_0.lon + cp_inCE_1.lon) / 2.;
					this.cp_inter.altitude = (cp_inCE_0.altitude + cp_inCE_1.altitude) / 2.;
					this.cp_inter.T = (cp_inCE_0.T + cp_inCE_1.T) / 2;
					this.cp_inter.Tflight = this.cp_inter.T;
					// this.cp_inter.Dv = (cp_inCE_0.Dv + cp_inCE_1.Dv)/2;
					// this.cp_inter.Dh = (cp_inCE_0.Dh + cp_inCE_1.Dh)/2;
					// this.cp_inter.Ds = (cp_inCE_0.Ds + cp_inCE_1.Ds)/2;
					this.cp_inter.RH = gribdata.getHumidityFromAltitude(this.cp_inter.lat, this.cp_inter.lon,
							this.cp_inter.altitude);
					;
					this.cp_inter.vz = vz(this.cp_inter);
					this.cp_inter.Theta = theta(this.cp_inter);
					this.cp_inter.thetadz = dthetadz(this.cp_inter);
					this.cp_inter.eps = eps_calculator(this.cp_inter);
					System.out.println("cp_inter.eps: " + this.cp_inter.eps);

					this.cp_inter.Dv = this.Calc_Dv(this.cp_inter.altitude, this.cp_inter);
					this.cp_inter.Dh = this.Calc_Dh(this.cp_inter.eps);
					this.cp_inter.Ds = (cp_inCE_0.Ds + cp_inCE_1.Ds) / 2;

					this.cp_inter.density = density(this.cp_inter);
					this.cp_inter.vo = cp_inCE_1.vo;
					this.cp_inter.uprime = cp_inCE_1.uprime;
					this.cp_inter.vx = gribdata.getWindUFromAltitude(this.cp_inter.lat, this.cp_inter.lon,
							this.cp_inter.altitude); // U-Component
					// of wind
					this.cp_inter.vy = gribdata.getWindVFromAltitude(this.cp_inter.lat, this.cp_inter.lon,
							this.cp_inter.altitude); // V-Component
					// of wind

					initial_state istate_inter = this.Contrail_depth(this.M_1, this.TAS_1,
							gribdata.getPressureFromAltitude(this.cp_inter.lat, this.cp_inter.lon,
									this.cp_inter.altitude),
							this.cp_inter.Tflight, this.cp_inter.Theta, this.cp_inter.thetadz, this.cp_inter.altitude,
							this.cp_inter.eps);
					this.cp_inter.sigma0v = (float) (istate_inter.Height / 2.2); // 50m
					this.cp_inter.sigma0h = (float) (istate_inter.Width / 2.2); // 20m

					dice_inter = (dice_0 + dice_1) / 2;
					t_inter = t_1;
					z_inter = z_1;
					t_universal = t_inter;

					sigmah_inter = sigmah_1;
					sigmav_inter = sigmav_1;
					sigmas_inter = sigmas_1;

					this.cp_inter.s = (gribdata.getWindspeedFromAltitude(this.cp_inter.lat, this.cp_inter.lon,
							this.cp_inter.altitude + sigmav_inter)
							- gribdata.getWindspeedFromAltitude(this.cp_inter.lat, this.cp_inter.lon,
									this.cp_inter.altitude - sigmav_inter))
							/ (2 * sigmav_inter);

					volumn_inter = intersection_volumn_cal(angle, cp_inCE_0, sigmav_0, sigmah_0, CCSc_0, cp_inCE_1,
							sigmav_1, sigmah_1, CCSc_1, this.distance_1);
					Nice_inter = Nice_0 * intersection_percentage_cal(volumn_inter, CCSc_0, this.distance_0)
							+ Nice_1 * intersection_percentage_cal(volumn_inter, CCSc_1, this.distance_1);
					Nice_rest_0 = Nice_0 * (1 - intersection_percentage_cal(volumn_inter, CCSc_0, this.distance_0));
					CCS_inter = ContrailCrossSection_inter_cal_v2(t_0, t_1, angle, distance_between_two_contrails,
							cp_inCE_0, sigmav_0, sigmah_0, cp_inCE_1, sigmav_1, sigmah_1);
					CCS_inter_esti = intersection_2D_size_esti(angle, cp_inCE_0, sigmav_0, sigmah_0, CCSc_0, cp_inCE_1,
							sigmav_1, sigmah_1, CCSc_1); // FIXME
					IWC_inter = IWC_0 + IWC_1;

					if (angle == 0) { // FIXME
						if (intersection_percentage_cal(volumn_inter, CCSc_1, this.distance_1) == 1.) {
							Nice_rest_1 = 0;
							c1_initial_full_inside_c0_flag = true;
							IWC_1 = 10e-8;
							RHice_1 = 99.9999;
							System.out.println(
									"c1_full_inside_c0_flag is set to true, this scenario 1 can not be calculated in code v_1,");
							System.out.println("which is simulated in v_2, p.s., to stop the loops,");
							System.out.println("IWC_1 is set to the artifical 10e-8 and RHice_1 is 99.9999.");
						} else {
							c1_initial_partly_inside_c0_flag = true;
							Nice_rest_1 = Nice_1
									* (1 - intersection_percentage_cal(volumn_inter, CCSc_1, this.distance_1));
						}
					} else {
						Nice_rest_1 = Nice_1 * (1 - intersection_percentage_cal(volumn_inter, CCSc_1, this.distance_1));
					}

					Ext_inter = new Extinktionseffizienz();
					lterr_inter = 10.0;
					bsol_inter = 0;

					eStar_inter = Math.exp(-6096.9385 / this.cp_inter.T + 16.635794 - 2.711193E-2 * this.cp_inter.T
							+ 1.673952E-5 * this.cp_inter.T * this.cp_inter.T + 2.433502 * Math.log(this.cp_inter.T));
					eStarIce_inter = Math.exp(-6024.5282 / this.cp_inter.T + 24.721994 + 1.0613868E-2 * this.cp_inter.T
							- 1.3198825E-5 * this.cp_inter.T * this.cp_inter.T
							- 0.49382577 * Math.log(this.cp_inter.T));

					RHice_inter = gribdata.getHumidityFromAltitude(this.cp_inter.lat, this.cp_inter.lon,
							this.cp_inter.altitude + z_inter) * eStar_inter / eStarIce_inter;

					if (cp_inter_initialization_flag == false) {
						String fn_output_inter = "records\\results\\s2_Ausgabe_LATm50_LON150_vortex_inter.txt";

						try {
							fwDataOut_inter = new FileWriter(fn_output_inter);
						} catch (IOException ex) {
							System.err.println("Fehler beim Anlegen der Datei");
						}
						;

						try {
							fwDataOut_inter.write(
									"angle, time_diff, t, z, eps, Dh, Dv, sigmah, sigmav, sigmas, rice, Nice, g_t, Qabs_t, Qsca_t, g_s, Qabs_s, Qsca_s, distance_between_contrails, state_indicator, CCS, CCS_esti, V_inter, IWC, RHice, Re, v_sed, windshear_s, wind_u, wind_v, wind_w, lat, lon, altitude\n");
							Ext_inter = new Extinktionseffizienz();

						} catch (IOException ex) {
							System.err.println("Fehler beim Schreiben");
						}
						;

						cp_inter_initialization_flag = true;

//							try {
//								double rext_inter = (dice_1+dice_0) / 4. * 1.e6;  
//								fwDataOut_inter
//										.write( (this.angle) + ", " + time_diff +", "+ t_inter + ", " + z_inter + ", " + sigmah_inter + ", " + sigmav_inter
//												+ ", " + sigmas_inter + ", " + dice_inter/2. + ", " + Nice_inter + ", "
//												+ Ext_inter.Calc_g_terr(rext_inter, lterr_inter) + ", " + Ext_inter.Calc_Qabs_terr(rext_inter, lterr_inter) + ", "
//												+ Ext_inter.Calc_Qsca_terr(rext_inter, lterr_inter) + ", " + Ext_inter.Calc_g_sol(rext_inter, bsol_inter) + ", "
//												+ Ext_inter.Calc_Qabs_sol(rext_inter, bsol_inter) + ", " + Ext_inter.Calc_Qsca_sol(rext_inter, bsol_inter) + ", "
//												+ CCS_inter + ", " + CCS_inter_esti +", " + volumn_inter +", "+ IWC_inter + ", " + RHice_inter + ", "
//												+ Re(dice_inter/2. , SedimentationSpeed(dice_inter, this.cp_inter),
//														kinvisc_Suth(gribdata.getPressureFromAltitude(this.cp_inter.lat, this.cp_inter.lon,
//																this.cp_inter.altitude), this.cp_inter.T))
//												+ ", " + SedimentationSpeed(dice_inter, this.cp_inter) + ", " + this.cp_inter.s + ", "
//												+ gribdata.getWindUFromAltitude(this.cp_inter.lat, this.cp_inter.lon, this.cp_inter.altitude) + ", "
//												+ gribdata.getWindVFromAltitude(this.cp_inter.lat, this.cp_inter.lon, this.cp_inter.altitude) + ", "
//												+ this.cp_inter.vz + ", " + ((int) (this.cp_inter.lat * 10000) / 10000.) + ", "
//												+ ((int) (this.cp_inter.lon * 10000) / 10000.) + ", " + Math.round(this.cp_inter.altitude) + "\n");
						//
//							} catch (IOException ex) {
//								System.err.println("Fehler beim Schreiben");
//							}
//							;
					}
					intersection_check_flag = true;
					System.out
							.println("Contrail 0 and 1 start intersected on time stamp t_0: " + t_0 + ", t_1: " + t_1);
				}

				if ((intersection_check_flag) && (evo_state_inter == false) && (evo_done_inter == false)) {

					// calculate ice particle diameter

					Nice_inter = Nice_0 * intersection_percentage_cal(volumn_inter, CCS_0, this.distance_0)
							+ Nice_1 * intersection_percentage_cal(volumn_inter, CCS_1, this.distance_1);

					// FIXME, which one to choose!
					dice_inter = Math.pow(IWC_inter * volumn_inter / (Nice_inter * 917) * 6 / Math.PI, 1 / 3.);
					// dice_inter = (dice_0 + dice_1)/2;

					// next time step update:
					t_inter = t_inter + delta_t_inter;

					// contrail geometry

					// Contrail horizontal movement
					// Calculates displaced contrail position (lat, lon)

					this.cp_inter.lat = this.cp_inter.lat + (gribdata.getWindVFromAltitude(this.cp_inter.lat,
							this.cp_inter.lon, this.cp_inter.altitude) * (delta_t_inter / (60000 * 1.852)));
					this.cp_inter.lon = this.cp_inter.lon + (gribdata.getWindUFromAltitude(this.cp_inter.lat,
							this.cp_inter.lon, this.cp_inter.altitude)
							* (delta_t_inter / (60000 * 1.852 * Math.cos(Math.toRadians(this.cp_inter.lat))))); // 60nm
					// per
					// degree
					// at lat = 0°
					// and 1852
					// metres per nm

					// Contrail vertical movement
					// Let contrail drift down
					this.cp_inter.vz = vz(this.cp_inter);
					delta_z_inter = (this.cp_inter.vz - SedimentationSpeed(dice_inter, this.cp_inter)) * delta_t_inter; // negative
					// delta_z
					// means
					// sinking
					// not
					// thinking
					z_inter = z_inter + delta_z_inter; // z ist die Summe der delta z je Durchlauf
					this.cp_inter.altitude = this.cp_inter.altitude + delta_z_inter;

//							System.out.println("contrail 1 lat:" + cp_inCE_1.lat);
//							System.out.println("contrail 1 lon:" + cp_inCE_1.lon);
//							System.out.println("contrail 1 alt:" + cp_inCE_1.altitude);

					sigmav_inter = Math.pow(2. * this.cp_inter.Dv * t_inter + Math.pow(this.cp_inter.sigma0v, 2), 0.5);
					this.cp_inter.s = (gribdata.getWindspeedFromAltitude(this.cp_inter.lat, this.cp_inter.lon,
							this.cp_inter.altitude + sigmav_inter)
							- gribdata.getWindspeedFromAltitude(this.cp_inter.lat, this.cp_inter.lon,
									this.cp_inter.altitude - sigmav_inter))
							/ (2 * sigmav_inter); // calculates
													// the
													// change of
													// windspeed
													// between
													// upper and
													// lower
													// boundary
													// of
													// contrail
													// per metre
					sigmah_inter = Math
							.pow(2. / 3. * this.cp_inter.s * this.cp_inter.s * this.cp_inter.Dv * Math.pow(t_inter, 3)
									+ (2 * this.cp_inter.Ds + this.cp_inter.s * Math.pow(this.cp_inter.sigma0v, 2))
											* this.cp_inter.s * Math.pow(t_inter, 2)
									+ 2 * this.cp_inter.Dh * t_inter + Math.pow(this.cp_inter.sigma0h, 2), 0.5);
					sigmas_inter = this.cp_inter.s * this.cp_inter.Dv * Math.pow(t_inter, 2)
							+ (2 * this.cp_inter.Ds + this.cp_inter.s * Math.pow(this.cp_inter.sigma0v, 2)) * t_inter;

//							System.out.println("contrail 1 sigmav_1:" + sigmav_1);
//							System.out.println("contrail 1 sigmah_1:" + sigmah_1);

					// System.out.println("Sigmas in main body: " + sigmas);
					// cp.Dh = 0.1 * cp.uprime * sigmah;

					// eStarIce before adiabatic heating
					eStarIce_old_inter = Math.exp(-6024.5282 / this.cp_inter.T + 24.721994
							+ 1.0613868E-2 * this.cp_inter.T - 1.3198825E-5 * this.cp_inter.T * this.cp_inter.T
							- 0.49382577 * Math.log(this.cp_inter.T));

					// Temperature
					this.cp_inter.T = (float) (this.cp_inter.Tflight - AdiabaticCurve(this.cp_inter) * z_inter); // adiabatic
																													// heating

					// calculate additional ice water content

					// SATURATION PRESSUREs FROM SONNTAG, Temp IN K, PSAT IN PA - METEOROL. Z., 3
					// (1994) 51-66.
					eStar_inter = Math.exp(-6096.9385 / this.cp_inter.T + 16.635794 - 2.711193E-2 * this.cp_inter.T
							+ 1.673952E-5 * this.cp_inter.T * this.cp_inter.T + 2.433502 * Math.log(this.cp_inter.T));
					eStarIce_inter = Math.exp(-6024.5282 / this.cp_inter.T + 24.721994 + 1.0613868E-2 * this.cp_inter.T
							- 1.3198825E-5 * this.cp_inter.T * this.cp_inter.T
							- 0.49382577 * Math.log(this.cp_inter.T));
					RHice_inter = gribdata.getHumidityFromAltitude(this.cp_inter.lat, this.cp_inter.lon,
							this.cp_inter.altitude) * eStar_inter / eStarIce_inter;

					CCSc_inter = ContrailCrossSection_inter_cal_v2(t_0, t_inter, angle, distance_between_two_contrails,
							cp_inCE_0, sigmav_0, sigmah_0, cp_inCE_1, sigmav_1, sigmah_1);

				}

				// time synchronization:
				t_universal = t_universal + delta_t_universal;
				// if calculate c0 and c1 independently, change Dh and Dv of "this.cp_inter" to
				// "cp_inCE_1"
				delta_t_universal = 0.1 * 2 * this.cp_inter.Dh * this.cp_inter.Dv * t_universal;
				if (delta_t_universal > 60) {
					delta_t_universal = 60;
				}

				// c0 parameter update:
				if ((evo_state_0 == false) && (evo_done_0 == false)) {
					if ((RHice_0 < 100) || (IWC_0 < 10e-8)) {
						System.out.println(
								"CCS and IWC keep the same with the values in the second last loop, which not updated any more!");
						System.out.println("Current RHice_0 = " + RHice_0 + " IWC_0 = " + IWC_0 + ", at t_0 = " + t_0);
						evo_state_0 = true;
					} else {

						if (intersection_check_flag == true) {
							// distance_between_two_contrails = contrail_centrals_distence_cal(cp_inCE_0,
							// cp_inCE_1);
							// System.out.println("distance_between_two_contrails
							// C0:"+distance_between_two_contrails);
							// System.out.println("sigmah_0 in code: "+sigmah_0);
							// System.out.println("sigmah_1 in code: "+sigmah_1);

							if (intersection_done_flag == false) {
								c0_esp = contrail0_exposed_surface_percentage_cal(angle, distance_between_two_contrails,
										sigmav_0, sigmah_0, CCS_0, sigmav_1, sigmah_1, CCS_1, this.distance_0);
							} else {
								c0_esp = 1;
							}
							System.out.println("c0_esp: " + c0_esp);
							delta_e_0 = c0_esp * eStarIce_0 * (RHice_0 / 100. - 1.);
						} else {
							delta_e_0 = eStarIce_0 * (RHice_0 / 100. - 1.);
						}

						IWCs_0 = 0.21667 * delta_e_0 / cp_inCE_0.T; // Masse des Wasser ueber Saettigung in der
																	// Umgebung
																	// (pro Volumen)
						// calculate new ice water content
						// before IWCs is added, IWC is reduced because of adiabatic heating
						IWC_0 = IWC_0 - (eStarIce_0 - eStarIce_old_0) / (462 * cp_inCE_0.T);

						// current CCSc is the CCSn in last code version, here should not be next, it's
						// current CCS.
						delta_t_0 = delta_t_universal;
						IWC_0 = (IWCs_0 * (CCSc_0 - CCS_0) + IWC_0 * CCS_0) / CCSc_0;
						// CCS means the old state.
						CCS_0 = CCSc_0;
						/*
						 * IWC and CCS affect the next loop calculation.
						 */
						evo_state_0 = false;

					}

				}

				// c1 parameter update:
				if ((evo_state_1 == false) && (evo_done_1 == false)) {
					if ((RHice_1 < 100) || (IWC_1 < 10e-8)) {
						System.out.println(
								"CCS and IWC keep the same with the values in the second last loop, which not updated any more!");
						System.out.println("Current RHice_1 = " + RHice_1 + " IWC_1 = " + IWC_1 + ", at t_1 = " + t_1);
						evo_state_1 = true;
					} else {

						if (intersection_check_flag == true) {
							// distance_between_two_contrails = contrail_centrals_distence_cal(cp_inCE_0,
							// cp_inCE_1);
							// System.out.println("distance_between_two_contrails
							// C1:"+distance_between_two_contrails);
							if (intersection_done_flag == false) {
								c1_esp = contrail1_exposed_surface_percentage_cal(angle, distance_between_two_contrails,
										sigmav_0, sigmah_0, CCS_0, sigmav_1, sigmah_1, CCS_1, this.distance_1);
							} else {
								c1_esp = 1;
							}
							System.out.println("c1_esp: " + c1_esp);
							delta_e_1 = c1_esp * eStarIce_1 * (RHice_1 / 100. - 1.);
						} else {
							delta_e_1 = eStarIce_1 * (RHice_1 / 100. - 1.);
						}

						IWCs_1 = 0.21667 * delta_e_1 / cp_inCE_1.T; // Masse des Wasser ueber Saettigung in der
																	// Umgebung
																	// (pro Volumen)
						// calculate new ice water content
						// before IWCs is added, IWC is reduced because of adiabatic heating
						IWC_1 = IWC_1 - (eStarIce_1 - eStarIce_old_1) / (462 * cp_inCE_1.T);

						// current CCSc is the CCSn in last code version, here should not be next, it's
						// current CCS.

						delta_t_1 = delta_t_universal;
						IWC_1 = (IWCs_1 * (CCSc_1 - CCS_1) + IWC_1 * CCS_1) / CCSc_1;
						// CCS means the old state.
						CCS_1 = CCSc_1;
						/*
						 * IWC and CCS affect the next loop calculation.
						 */
						evo_state_1 = false;

					}

				}

				if ((intersection_check_flag == true) && (evo_state_inter == false) && (evo_done_inter == false)) {
					if ((RHice_inter < 100) || (IWC_inter < 10e-8)) {
						System.out.println(
								"CCS and IWC keep the same with the values in the second last loop, which not updated any more!");
						System.out.println("Current RHice_inter = " + RHice_inter + " IWC_inter = " + IWC_inter
								+ ", at t_inter = " + t_inter);
						evo_state_inter = true;
						intersection_done_flag = true;

						volumn_inter = intersection_volumn_cal(angle, cp_inCE_0, sigmav_0, sigmah_0, CCS_0, cp_inCE_1,
								sigmav_1, sigmah_1, CCS_1, this.distance_1);
						Nice_rest_0 = Nice_0 * (1 - intersection_percentage_cal(volumn_inter, CCS_0, this.distance_0));
						Nice_rest_1 = Nice_1 * (1 - intersection_percentage_cal(volumn_inter, CCS_1, this.distance_1));
					} else {
						// get on with it
						delta_e_inter = 0;
						IWCs_inter = 0.21667 * delta_e_inter / this.cp_inter.T; // Masse des Wasser ueber Saettigung
																				// in der
						// Umgebung
						// (pro Volumen)
						// calculate new ice water content
						// before IWCs is added, IWC is reduced because of adiabatic heating
						IWC_inter = IWC_inter - (eStarIce_inter - eStarIce_old_inter) / (462 * this.cp_inter.T);

						// current CCSc is the CCSn in last code version, here should not be next, it's
						// current CCS.

						// distance_between_two_contrails = contrail_centrals_distence_cal(cp_inCE_0,
						// cp_inCE_1);

						System.out.println("distance_between_two_contrails C_inter:" + distance_between_two_contrails);
						System.out.println("=======================");
						System.out.println("Contrail 0 and 1 intersected on time stamp t_0: " + t_0 + ", t_1: "
								+ t_inter + "; distance between two contrails (m): " + distance_between_two_contrails);

						delta_t_inter = delta_t_universal;
						IWC_inter = (IWCs_inter * (CCSc_inter - CCS_inter) + IWC_inter * CCS_inter) / CCSc_inter;
						// CCS means the old state.
						CCS_inter = CCSc_inter;
						/*
						 * IWC and CCS affect the next loop calculation.
						 */
						evo_state_inter = false;

					}

				}

				if (intersection_check_flag == true) {
					if (c1_initial_full_inside_c0_flag) {

						volumn_inter = intersection_volumn_cal(angle, cp_inCE_0, sigmav_0, sigmah_0, CCS_0,
								this.cp_inter, sigmav_inter, sigmah_inter, CCS_inter, this.distance_0);
						CCS_inter_esti = intersection_2D_size_esti(angle, cp_inCE_0, sigmav_0, sigmah_0, CCS_0,
								this.cp_inter, sigmav_inter, sigmah_inter, CCS_inter);
						CCS_rest_0 = CCS_0 - CCS_inter_esti;
						CCS_rest_1 = 0;
						V_rest_0 = CCS_rest_0 * this.distance_0;
						V_rest_1 = CCS_rest_1 * this.distance_1;
					} else if (c1_initial_partly_inside_c0_flag) {
						volumn_inter = intersection_volumn_cal(angle, cp_inCE_0, sigmav_0, sigmah_0, CCS_0, cp_inCE_1,
								sigmav_1, sigmah_1, CCS_1, this.distance_1);
						CCS_inter_esti = intersection_2D_size_esti(angle, cp_inCE_0, sigmav_0, sigmah_0, CCS_0,
								cp_inCE_1, sigmav_1, sigmah_1, CCS_1);
						CCS_rest_0 = CCS_0 - CCS_inter_esti;
						CCS_rest_1 = CCS_1 - CCS_inter_esti;
						V_rest_0 = CCS_rest_0 * this.distance_0;
						V_rest_1 = CCS_rest_1 * this.distance_1;
					} else {
						volumn_inter = intersection_volumn_cal(angle, cp_inCE_0, sigmav_0, sigmah_0, CCS_0, cp_inCE_1,
								sigmav_1, sigmah_1, CCS_1, this.distance_1);
						CCS_inter_esti = intersection_2D_size_esti(angle, cp_inCE_0, sigmav_0, sigmah_0, CCS_0,
								cp_inCE_1, sigmav_1, sigmah_1, CCS_1);
						CCS_rest_0 = CCS_0;
						CCS_rest_1 = CCS_1;
						V_rest_0 = CCS_0 * this.distance_0 - volumn_inter;
						V_rest_1 = CCS_1 * this.distance_1 - volumn_inter;
					}
				} else {
					Nice_rest_0 = Nice_0;
					CCS_rest_0 = CCS_0;
					V_rest_0 = CCS_rest_0 * this.distance_0;
					Nice_rest_1 = Nice_1;
					CCS_rest_1 = CCS_1;
					V_rest_1 = CCS_rest_1 * this.distance_1;
				}

				// contrail 1 records
				if ((c1_initial_full_inside_c0_flag == false) && (evo_state_1 == false)) {
					try {
						rext_1 = dice_1 / 2. * 1.e6;
						this.fwDataOut_1.write((this.angle) + ", " + time_diff + ", " + t_1 + ", " + z_1 + ", "
								+ cp_inCE_1.eps + ", " + cp_inCE_1.Dh + ", " + cp_inCE_1.Dv + ", " + sigmah_1 + ", "
								+ sigmav_1 + ", " + sigmas_1 + ", " + (dice_1 / 2) + ", " + Nice_rest_1 + ", "
								+ Ext_1.Calc_g_terr(rext_1, lterr_1) + ", " + Ext_1.Calc_Qabs_terr(rext_1, lterr_1)
								+ ", " + Ext_1.Calc_Qsca_terr(rext_1, lterr_1) + ", " + Ext_1.Calc_g_sol(rext_1, bsol_1)
								+ ", " + Ext_1.Calc_Qabs_sol(rext_1, bsol_1) + ", "
								+ Ext_1.Calc_Qsca_sol(rext_1, bsol_1) + ", " + distance_between_two_contrails + ", "
								+ s2_intersection_state_indicator + ", " + ((int) (c1_esp * 10000) / 10000.) + ", "
								+ CCS_1 + ", " + CCS_rest_1 + ", " + (CCS_1 * this.distance_1) + ", " + V_rest_1 + ", "
								+ IWC_1 + ", " + RHice_1 + ", "
								+ Re(dice_1 / 2., SedimentationSpeed(dice_1, cp_inCE_1),
										kinvisc_Suth(gribdata.getPressureFromAltitude(cp_inCE_1.lat, cp_inCE_1.lon,
												cp_inCE_1.altitude), cp_inCE_1.T))
								+ ", " + SedimentationSpeed(dice_1, cp_inCE_1) + ", " + cp_inCE_1.s + ", "
								+ gribdata.getWindUFromAltitude(cp_inCE_1.lat, cp_inCE_1.lon, cp_inCE_1.altitude) + ", "
								+ gribdata.getWindVFromAltitude(cp_inCE_1.lat, cp_inCE_1.lon, cp_inCE_1.altitude) + ", "
								+ cp_inCE_1.vz + ", " + ((int) (cp_inCE_1.lat * 10000) / 10000.) + ", "
								+ ((int) (cp_inCE_1.lon * 10000) / 10000.) + ", " + Math.round(cp_inCE_1.altitude)
								+ "\n");

					} catch (IOException ex) {
						System.err.println("Fehler beim Schreiben");
					}
					;
				} else if ((c1_initial_full_inside_c0_flag == false) && (evo_done_1 == false)) {
					try {
						rext_1 = dice_1 / 2. * 1.e6;
						this.fwDataOut_1.write((this.angle) + ", " + time_diff + ", " + t_1 + ", " + z_1 + ", "
								+ cp_inCE_1.eps + ", " + cp_inCE_1.Dh + ", " + cp_inCE_1.Dv + ", " + sigmah_1 + ", "
								+ sigmav_1 + ", " + sigmas_1 + ", " + (dice_1 / 2) + ", " + Nice_rest_1 + ", "
								+ Ext_1.Calc_g_terr(rext_1, lterr_1) + ", " + Ext_1.Calc_Qabs_terr(rext_1, lterr_1)
								+ ", " + Ext_1.Calc_Qsca_terr(rext_1, lterr_1) + ", " + Ext_1.Calc_g_sol(rext_1, bsol_1)
								+ ", " + Ext_1.Calc_Qabs_sol(rext_1, bsol_1) + ", "
								+ Ext_1.Calc_Qsca_sol(rext_1, bsol_1) + ", " + distance_between_two_contrails + ", "
								+ s2_intersection_state_indicator + ", " + ((int) (c1_esp * 10000) / 10000.) + ", "
								+ CCS_1 + ", " + CCS_rest_1 + ", " + (CCS_1 * this.distance_1) + ", " + V_rest_1 + ", "
								+ IWC_1 + ", " + RHice_1 + ", "
								+ Re(dice_1 / 2., SedimentationSpeed(dice_1, cp_inCE_1),
										kinvisc_Suth(gribdata.getPressureFromAltitude(cp_inCE_1.lat, cp_inCE_1.lon,
												cp_inCE_1.altitude), cp_inCE_1.T))
								+ ", " + SedimentationSpeed(dice_1, cp_inCE_1) + ", " + cp_inCE_1.s + ", "
								+ gribdata.getWindUFromAltitude(cp_inCE_1.lat, cp_inCE_1.lon, cp_inCE_1.altitude) + ", "
								+ gribdata.getWindVFromAltitude(cp_inCE_1.lat, cp_inCE_1.lon, cp_inCE_1.altitude) + ", "
								+ cp_inCE_1.vz + ", " + ((int) (cp_inCE_1.lat * 10000) / 10000.) + ", "
								+ ((int) (cp_inCE_1.lon * 10000) / 10000.) + ", " + Math.round(cp_inCE_1.altitude)
								+ "\n");

					} catch (IOException ex) {
						System.err.println("Fehler beim Schreiben");
					}
					;
					evo_done_1 = true;
				}
			}

			// contrail 0 records
			if (evo_state_0 == false) {
				try {
					rext_0 = dice_0 / 2. * 1.e6;
					this.fwDataOut_0.write(t_0 + ", " + z_0 + ", " + cp_inCE_0.eps + ", " + cp_inCE_0.Dh + ", "
							+ cp_inCE_0.Dv + ", " + sigmah_0 + ", " + sigmav_0 + ", " + sigmas_0 + ", " + (dice_0 / 2)
							+ ", " + Nice_rest_0 + ", " + Ext_0.Calc_g_terr(rext_0, lterr_0) + ", "
							+ Ext_0.Calc_Qabs_terr(rext_0, lterr_0) + ", " + Ext_0.Calc_Qsca_terr(rext_0, lterr_0)
							+ ", " + Ext_0.Calc_g_sol(rext_0, bsol_0) + ", " + Ext_0.Calc_Qabs_sol(rext_0, bsol_0)
							+ ", " + Ext_0.Calc_Qsca_sol(rext_0, bsol_0) + ", " + s2_intersection_state_indicator + ", "
							+ ((int) (c0_esp * 10000) / 10000.) + ", " + CCS_0 + ", " + CCS_rest_0 + ", "
							+ (CCS_0 * this.distance_0) + ", " + V_rest_0 + ", " + IWC_0 + ", " + RHice_0 + ", "
							+ Re(dice_0 / 2., SedimentationSpeed(dice_0, cp_inCE_0),
									kinvisc_Suth(gribdata.getPressureFromAltitude(cp_inCE_0.lat, cp_inCE_0.lon,
											cp_inCE_0.altitude), cp_inCE_0.T))
							+ ", " + SedimentationSpeed(dice_0, cp_inCE_0) + ", " + cp_inCE_0.s + ", "
							+ gribdata.getWindUFromAltitude(cp_inCE_0.lat, cp_inCE_0.lon, cp_inCE_0.altitude) + ", "
							+ gribdata.getWindVFromAltitude(cp_inCE_0.lat, cp_inCE_0.lon, cp_inCE_0.altitude) + ", "
							+ cp_inCE_0.vz + ", " + ((int) (cp_inCE_0.lat * 10000) / 10000.) + ", "
							+ ((int) (cp_inCE_0.lon * 10000) / 10000.) + ", " + Math.round(cp_inCE_0.altitude) + "\n");

				} catch (IOException ex) {
					System.err.println("Fehler beim Schreiben");
				}
				;
			} else if (evo_done_0 == false) {
				try {
					rext_0 = dice_0 / 2. * 1.e6;
					this.fwDataOut_0.write(t_0 + ", " + z_0 + ", " + cp_inCE_0.eps + ", " + cp_inCE_0.Dh + ", "
							+ cp_inCE_0.Dv + ", " + sigmah_0 + ", " + sigmav_0 + ", " + sigmas_0 + ", " + (dice_0 / 2)
							+ ", " + Nice_rest_0 + ", " + Ext_0.Calc_g_terr(rext_0, lterr_0) + ", "
							+ Ext_0.Calc_Qabs_terr(rext_0, lterr_0) + ", " + Ext_0.Calc_Qsca_terr(rext_0, lterr_0)
							+ ", " + Ext_0.Calc_g_sol(rext_0, bsol_0) + ", " + Ext_0.Calc_Qabs_sol(rext_0, bsol_0)
							+ ", " + Ext_0.Calc_Qsca_sol(rext_0, bsol_0) + ", " + s2_intersection_state_indicator + ", "
							+ ((int) (c0_esp * 10000) / 10000.) + ", " + CCS_0 + ", " + CCS_rest_0 + ", "
							+ (CCS_0 * this.distance_0) + ", " + V_rest_0 + ", " + IWC_0 + ", " + RHice_0 + ", "
							+ Re(dice_0 / 2., SedimentationSpeed(dice_0, cp_inCE_0),
									kinvisc_Suth(gribdata.getPressureFromAltitude(cp_inCE_0.lat, cp_inCE_0.lon,
											cp_inCE_0.altitude), cp_inCE_0.T))
							+ ", " + SedimentationSpeed(dice_0, cp_inCE_0) + ", " + cp_inCE_0.s + ", "
							+ gribdata.getWindUFromAltitude(cp_inCE_0.lat, cp_inCE_0.lon, cp_inCE_0.altitude) + ", "
							+ gribdata.getWindVFromAltitude(cp_inCE_0.lat, cp_inCE_0.lon, cp_inCE_0.altitude) + ", "
							+ cp_inCE_0.vz + ", " + ((int) (cp_inCE_0.lat * 10000) / 10000.) + ", "
							+ ((int) (cp_inCE_0.lon * 10000) / 10000.) + ", " + Math.round(cp_inCE_0.altitude) + "\n");

				} catch (IOException ex) {
					System.err.println("Fehler beim Schreiben");
				}
				;
				evo_done_0 = true;
			}

			// contrail intersection records
			if ((intersection_check_flag) && (evo_state_inter == false)) {
				try {
					double rext_inter = dice_inter / 2. * 1.e6;
					System.out.println("Nice_inter = " + Nice_inter);
					fwDataOut_inter.write((this.angle) + ", " + time_diff + ", " + t_inter + ", " + z_inter + ", "
							+ this.cp_inter.eps + ", " + this.cp_inter.Dh + ", " + this.cp_inter.Dv + ", "
							+ sigmah_inter + ", " + sigmav_inter + ", " + sigmas_inter + ", " + dice_inter / 2. + ", "
							+ Nice_inter + ", " + Ext_inter.Calc_g_terr(rext_inter, lterr_inter) + ", "
							+ Ext_inter.Calc_Qabs_terr(rext_inter, lterr_inter) + ", "
							+ Ext_inter.Calc_Qsca_terr(rext_inter, lterr_inter) + ", "
							+ Ext_inter.Calc_g_sol(rext_inter, bsol_inter) + ", "
							+ Ext_inter.Calc_Qabs_sol(rext_inter, bsol_inter) + ", "
							+ Ext_inter.Calc_Qsca_sol(rext_inter, bsol_inter) + ", "
							+ distance_between_two_contrails / 2 + ", " + s2_intersection_state_indicator + ", "
							+ CCS_inter + ", " + CCS_inter_esti + ", " + volumn_inter + ", " + IWC_inter + ", "
							+ RHice_inter + ", "
							+ Re(dice_inter / 2., SedimentationSpeed(dice_inter, this.cp_inter),
									kinvisc_Suth(gribdata.getPressureFromAltitude(this.cp_inter.lat, this.cp_inter.lon,
											this.cp_inter.altitude), this.cp_inter.T))
							+ ", " + SedimentationSpeed(dice_inter, this.cp_inter) + ", " + this.cp_inter.s + ", "
							+ gribdata
									.getWindUFromAltitude(this.cp_inter.lat, this.cp_inter.lon, this.cp_inter.altitude)
							+ ", "
							+ gribdata.getWindVFromAltitude(this.cp_inter.lat, this.cp_inter.lon,
									this.cp_inter.altitude)
							+ ", " + this.cp_inter.vz + ", " + ((int) (this.cp_inter.lat * 10000) / 10000.) + ", "
							+ ((int) (this.cp_inter.lon * 10000) / 10000.) + ", " + Math.round(this.cp_inter.altitude)
							+ "\n");

				} catch (IOException ex) {
					System.err.println("Fehler beim Schreiben");
				}
				;
			} else if ((intersection_check_flag) && (evo_done_inter == false)) {
				try {
					double rext_inter = dice_inter / 2. * 1.e6;
					fwDataOut_inter.write((this.angle) + ", " + time_diff + ", " + t_inter + ", " + z_inter + ", "
							+ this.cp_inter.eps + ", " + this.cp_inter.Dh + ", " + this.cp_inter.Dv + ", "
							+ sigmah_inter + ", " + sigmav_inter + ", " + sigmas_inter + ", " + dice_inter / 2. + ", "
							+ Nice_inter + ", " + Ext_inter.Calc_g_terr(rext_inter, lterr_inter) + ", "
							+ Ext_inter.Calc_Qabs_terr(rext_inter, lterr_inter) + ", "
							+ Ext_inter.Calc_Qsca_terr(rext_inter, lterr_inter) + ", "
							+ Ext_inter.Calc_g_sol(rext_inter, bsol_inter) + ", "
							+ Ext_inter.Calc_Qabs_sol(rext_inter, bsol_inter) + ", "
							+ Ext_inter.Calc_Qsca_sol(rext_inter, bsol_inter) + ", "
							+ distance_between_two_contrails / 2 + ", " + s2_intersection_state_indicator + ", "
							+ CCS_inter + ", " + CCS_inter_esti + ", " + volumn_inter + ", " + IWC_inter + ", "
							+ RHice_inter + ", "
							+ Re(dice_inter / 2., SedimentationSpeed(dice_inter, this.cp_inter),
									kinvisc_Suth(gribdata.getPressureFromAltitude(this.cp_inter.lat, this.cp_inter.lon,
											this.cp_inter.altitude), this.cp_inter.T))
							+ ", " + SedimentationSpeed(dice_inter, this.cp_inter) + ", " + this.cp_inter.s + ", "
							+ gribdata
									.getWindUFromAltitude(this.cp_inter.lat, this.cp_inter.lon, this.cp_inter.altitude)
							+ ", "
							+ gribdata.getWindVFromAltitude(this.cp_inter.lat, this.cp_inter.lon,
									this.cp_inter.altitude)
							+ ", " + this.cp_inter.vz + ", " + ((int) (this.cp_inter.lat * 10000) / 10000.) + ", "
							+ ((int) (this.cp_inter.lon * 10000) / 10000.) + ", " + Math.round(this.cp_inter.altitude)
							+ "\n");

				} catch (IOException ex) {
					System.err.println("Fehler beim Schreiben");
				}
				;
				evo_done_inter = true;
			}

		} while (((IWC_0 > 10e-8) && (RHice_0 >= 100)) || ((IWC_1 > 10e-8) && (RHice_1 >= 100))
				|| ((IWC_inter > 10e-8) && (RHice_inter >= 100)));

		if (RHice_0 >= 100) {
			System.out.println("IWC_0 = " + IWC_0 + " <= 10e-8, at t_0 = " + t_0);
		}

		if (RHice_1 >= 100) {
			System.out.println("IWC_1 = " + IWC_1 + " <= 10e-8, at t_1 = " + t_1);
		}

		if (RHice_inter >= 100) {
			System.out.println("IWC_inter = " + IWC_inter + " <= 10e-8, at t_inter = " + t_inter);
		}

		try {
			this.fwDataOut_0.close();
			this.fwDataOut_1.close();
			this.fwDataOut_inter.close();
		} catch (IOException ex) {
			System.err.println("Fehler beim Schliessen der Datei");
		}
		;

		return;
	}

	public static void run() {
		wake_vortex_intersection_v1_for_s2 wvi = new wake_vortex_intersection_v1_for_s2(
				"Config_vortex0.txt", "Config_vortex1_for_s2.txt");
		wvi.ContrailEvolution();
	}

	public static void main(String[] args) {
		run();
	}
}