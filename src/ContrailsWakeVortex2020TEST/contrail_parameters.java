package ContrailsWakeVortex2020TEST;

public class contrail_parameters
{
//  public float sabs_t, sabs_s, ssca_t, ssca_s; // keine Konstanten im Contrail
 public double Dv, Dh, Dht, Ds, vo, s, sigma0h, sigma0v;
 //meteorology
 public double T, Tflight, RH, vx, vy, vz, uprime, density, Theta, sat_vap_pres, Rs_humidair, thetadz, eps;
 //location
 public double lat, lon, altitude;
 
 public contrail_parameters()
 {
  // Dv, Dh, Ds... vertical, horizontal and sheared diffusion [m2/s]
  this.Dv=-777;
  this.Dh=-777;
  this.Ds=-777;
  this.Dht=0.;
  // vo...mean value of wind velocity at flight level [m/s]
  this.vo=-777;
  // s...constant wind shear 
  this.s=-777;
  // sigma0h/v...initial standard deviation describing contrail width and height
  this.sigma0h=-777;
  this.sigma0v=-777;  
  
  // meteorology
  this.T=-777; //temperature of contrail in K
  this.Tflight=-777; //temperature at flightlevel in K
  this.RH=-777; //relative humidity in %
  this.vx=-777; // wind velocity in m/s
  this.vy=-777; // wind velocity in m/s
  this.vz=-777; // wind velocity in m/s
  this.uprime=-777; // horizontal wind fluctuation in m/s
  //NEU
  this.density=-777; //density
  this.Theta=-779; //potential temperature
  this.sat_vap_pres=-777; //saturation vapour pressure
  this.Rs_humidair=-777; //specific gas constant of humid air
  this.thetadz=-777; //theta dz
  this.eps=-777; //Eddy Rate
  this.lat=-777;
  this.lon=-777;
  this.altitude=-777;
 }
}