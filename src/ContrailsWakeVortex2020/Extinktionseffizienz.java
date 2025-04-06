package ContrailsWakeVortex2020;
// Status 120731
//
// ---Berechnung der Extinktionseffizienten in Abh. von Partikeldurchmesser d_ice---
//
// ---== Solares Spektrum ==---
//
//Berechnet solare Extinktionsefizienten nach Yang, Journal of Geophysical Research, Vol 105, No. D4, pp. 4699-4718, 2000
//uebergeben: Koeffizienten aus Yang2000 in Abh. von Wellenlaenge Lambda
//
// for individual ice crystals
// public double Calc_Qext_sol_i(double Dmax, int band, int shape)
// public double Calc_Qabs_sol_i(double Dmax, int band, int shape)
// public double Calc_Qsca_sol_i(double Dmax, int band, int shape)
// public double Calc_g_sol_i(double Dmax, int band, int shape)
//
// for mixtures
// public double Calc_Qext_sol(double Dmax, int band)
// public double Calc_Qabs_sol(double Dmax, int band)
// public double Calc_Qsca_sol(double Dmax, int band)
// public double Calc_g_sol(double Dmax, int band)
//
// Parameters:
//  Dmax ... maximum dimension [micrometer]
//  band ... das gewuenschte spektrale Band
//     0 = 0.55 micrometer
//     1 = 1.35 micrometer
//     2 = 2.25 micrometer
//     3 = 3.0125 micrometer
//     4 = 3.775 micrometer
//     5 = 4.5 micrometer
//  shape ... 0 = Plate
//            1 = Column
//            2 = Hollow Column
//            3 = Bullet rosettes-4
//            4 = Bullet rosettes-6
//            5 = Aggregates
//
//
// ---== Terrestrisches Spektrum ==---
//
//Berechnet terrestrische Extinktionsefizienten nach Yang, Applied Optics, Vol 44, No. 26, 2005
//uebergeben: Koeffizienten aus Yang2005_terrestrSpektrum.txt eta, xi und zeta in Abh. von Wellenlaenge Lambda
//
// public double Calc_Qext_terr(double Dmax, double lambda)
// public double Calc_Qabs_terr(double Dmax, double lambda)
// public double Calc_Qsca_terr(double Dmax, double lambda)
// public double Calc_Qg_terr(double Dmax, double lambda)
//
// Parameters:
//  Dmax ... maximum dimension [micrometer] 2 < Dmax < 10'000
//  lambda ... wavelength [micrometer] 3.08 < lambda < 99.99

import java.lang.Math;

// for file handling
import java.io.FileReader;
import java.io.BufferedReader;
import java.io.IOException;
import java.lang.IndexOutOfBoundsException;
import java.util.regex.Pattern;
import java.util.*;

public class Extinktionseffizienz
{
 private String fn_output;
 private double[] bands;
  // Wellenlaenge des gewaehlten spektralen Bands
  // 0 = 0.55 micrometer
  // 1 = 1.35 micrometer
  // 2 = 2.25 micrometer
  // 3 = 2.75 micrometer
  // 4 = 3.0125 micrometer
  // 5 = 3.775 micrometer
  // 6 = 4.5 micrometer
 private double[][] table1;
  // Index 1: das gewaehlte spektrale Band [0...6]
  // Index 2: 0 = lambda
  //          1 = lambda1
  //          2 = lambda2
  //          3 = mr
  //          4 = mi
  //          5 = dS/S
 private double[][][] table3;
  // Index 1: das gewaehlte spektrale Band [0...6]
  // Index 2: shape 0 = Plate
  //           1 = Column
  //           2 = Hollow Column
  //           3 = Bullet rosettes-4
  //           4 = Bullet rosettes-6
  //           5 = Aggregates
  // Index 3: 0 = lambda
  //          1 = eta1_ext
  //          2 = eta2_ext
  //          3 = beta1
  //          4 = beta2
  //          5 = alpha
  //          6 = xi_ext
 private double[][][] table4;
  // Index 1: das gewaehlte spektrale Band [0...6]
  // Index 2: shape 0 = Plate
  //           1 = Column
  //           2 = Hollow Column
  //           3 = Bullet rosettes-4
  //           4 = Bullet rosettes-6
  //           5 = Aggregates
  // Index 3: 0 = lambda
  //          1 = eta1_abs
  //          2 = eta2_abs
  //          3 = xi1_abs
  //          4 = xi2_abs
 private double[][][] table5;
  // Index 1: das gewaehlte spektrale Band [0...6]
  // Index 2: shape 0 = Plate
  //           1 = Column
  //           2 = Hollow Column
  //           3 = Bullet rosettes-4
  //           4 = Bullet rosettes-6
  //           5 = Aggregates
  // Index 3: 0 = lambda
  //          1 = xi1_g
  //          2 = xi2_g
  //          3 = xi3_g
  //          4 = xi4_g
  //          5 = xi5_g
  //          6 = xi6_g
  //          7 = xi7_g
 private double[][] table_terr1;
  // Fuer die Berechnung im terrestrischen Wellenlaengenbereich
  // Index 1: das gewaehlte spektrale Band [0...6]
  // Index 2: 0 = lambda
  //          1 = eta1
  //          2 = eta2
  //          3 = eta3
  //          4 = xi0
  //          5 = xi1
  //          6 = xi2
  //          7 = xi3
  //          8 = zeta0
  //          9 = zeta1
  //         10 = zeta2
  //         11 = zeta3
 
 public Extinktionseffizienz(){
  // Erzeuger
  // Hier werden alle notwendigen Parameter eingelesen.
  bands = new double[] {0.55, 1.35, 2.25, 2.75, 3.0125, 3.775, 4.5}; //[micrometer]
  
  
  // ---== Read Table 1 ==---
  // Index 1: das gewaehlte spektrale Band [0...6]
  // Index 2: 0 = lambda
  //          1 = lambda1
  //          2 = lambda2
  //          3 = mr
  //          4 = mi
  //          5 = dS/S
  BufferedReader file= null;
  try {
   file=new BufferedReader(new FileReader("Parameter/Yang2000_OpticalConstants.txt"));
   System.out.println("Datei Parameter/Yang2000_OpticalConstants.txt gefunden" );
  } catch(IOException e) {System.err.println("Datei Parameter/Yang2000_OpticalConstants.txt nicht gefunden" );};
    
  //Datei ist geoeffnet und kann ueber file gelesen werden
  // Zerlegen der Eingangsdaten
  String eineZeile; // zum Zeilenweise einlesen der Datei
  Pattern p = Pattern.compile( "\\s" ); // wird zum Teilen der Zeilen verwendet
  ArrayList<String> list= new ArrayList<String>(); // wird mit den Messwerten gefuellt
  try {
   // ignore header (1 line)
   eineZeile= file.readLine();
   while (( eineZeile= file.readLine())!=null){ // Solange noch eine Zeile gelesen werden kann
    list.add(eineZeile);
   }
   file.close(); // schliessen der Datei
  } catch (IOException ex){System.err.println("Fehler beim Einlesen (Table1)");} // Falls ein Fehler auftrat
    
  // lege Daten in Array table1 ab
  table1= new double[7][6]; // Initialisieren des Arrays mit der richtigen Groesse
  for (int i=0; i<7; i++) {
   try {
    String InfoString[] = p.split( list.get(i) ); // Teile die Zeile in ein Array aus Stringvariablen
	for (int j = 0; j<6; j++){
     table1[i][j]=Double.valueOf(InfoString[j]);
    };
   } catch(IndexOutOfBoundsException e){
    System.err.println("Fehler beim Ablegen der Daten in Table1,");
    System.err.println("wahrscheinlich ist ein Formatierungsfehler in Parameter/Yang2000_OpticalConstants.txt" );
   };
  };
  
  
  // ---== Read Table 3 ==---
  // Index 1: das gewaehlte spektrale Band [0...6]
  // Index 2: shape 0 = Plate
  //           1 = Column
  //           2 = Hollow Column
  //           3 = Bullet rosettes-4
  //           4 = Bullet rosettes-6
  //           5 = Aggregates
  // Index 3: 0 = lambda
  //          1 = eta1_ext
  //          2 = eta2_ext
  //          3 = beta1
  //          4 = beta2
  //          5 = alpha
  //          6 = xi_ext
  file= null;
  try {
   file=new BufferedReader(new FileReader("Parameter/Yang2000_Extinction.txt"));
   System.out.println("Datei Parameter/Yang2000_Extinction.txt gefunden" );
  } catch(IOException e) {System.err.println("Datei Parameter/Yang2000_Extinction.txt nicht gefunden" );};
    
  //Datei ist geoeffnet und kann ueber file gelesen werden
  // Zerlegen der Eingangsdaten
  p = Pattern.compile( "\\s" ); // wird zum Teilen der Zeilen verwendet
  list= new ArrayList<String>(); // wird mit den Messwerten gefuellt
  try {
   // ignore header (2 lines)
   eineZeile= file.readLine();
   eineZeile= file.readLine();
   int i = 1;
   while (( eineZeile= file.readLine())!=null){ // Solange noch eine Zeile gelesen werden kann
    if ((i%8)!=0){
     list.add(eineZeile);
     }
    i++;
   }
   file.close(); // schliessen der Datei
  } catch (IOException ex){System.err.println("Fehler beim Einlesen (Table 3)");} // Falls ein Fehler auftrat
    
  // lege Daten in Array table3 ab
  table3= new double[7][6][7]; // Initialisieren des Arrays mit der richtigen Groesse
  for (int i=0; i<7; i++) { // das spektrale Band
   for (int j=0; j<6; j++) { // shape
    try {
     String InfoString[] = p.split( list.get(i+7*j) ); // Teile die Zeile in ein Array aus Stringvariablen
	 for (int k = 0; k<7; k++) { // Parameter
      table3[i][j][k]=Double.valueOf(InfoString[k]);
     };
    } catch(IndexOutOfBoundsException e) {
     System.err.println("Fehler beim Ablegen der Daten in Table3,");
     System.err.println("wahrscheinlich ist ein Formatierungsfehler in Parameter/Yang2000_Extinction.txt" );
    };
   };
  };
  
  
  // ---== Read Table 4 ==---
  // Index 1: das gewaehlte spektrale Band [0...6]
  // Index 2: shape 0 = Plate
  //           1 = Column
  //           2 = Hollow Column
  //           3 = Bullet rosettes-4
  //           4 = Bullet rosettes-6
  //           5 = Aggregates
  // Index 3: 0 = lambda
  //          1 = eta1_abs
  //          2 = eta2_abs
  //          3 = xi1_abs
  //          4 = xi2_abs
  file= null;
  try {
   file=new BufferedReader(new FileReader("Parameter/Yang2000_Absorption.txt"));
   System.out.println("Datei Parameter/Yang2000_Absorption.txt gefunden" );
  } catch(IOException e) {System.err.println("Datei Parameter/Yang2000_Absorption.txt nicht gefunden" );};
    
  //Datei ist geoeffnet und kann ueber file gelesen werden
  // Zerlegen der Eingangsdaten
  p = Pattern.compile( "\\s" ); // wird zum Teilen der Zeilen verwendet
  list= new ArrayList<String>(); // wird mit den Messwerten gefuellt
  try {
   // ignore header (2 lines)
   eineZeile= file.readLine();
   eineZeile= file.readLine();
   int i = 1;
   while (( eineZeile= file.readLine())!=null){ // Solange noch eine Zeile gelesen werden kann
    if ((i%8)!=0){
     list.add(eineZeile);
     }
    i++;
   }
   file.close(); // schliessen der Datei
  } catch (IOException ex){System.err.println("Fehler beim Einlesen (Table 4)");} // Falls ein Fehler auftrat
    
  // lege Daten in Array table4 ab
  table4= new double[7][6][5]; // Initialisieren des Arrays mit der richtigen Groesse
  for (int i=0; i<7; i++) { // das spektrale Band
   for (int j=0; j<6; j++) { // shape
    try {
     String InfoString[] = p.split( list.get(i+7*j) ); // Teile die Zeile in ein Array aus Stringvariablen
	 for (int k = 0; k<5; k++) { // Parameter
      table4[i][j][k]=Double.valueOf(InfoString[k]);
     };
    } catch(IndexOutOfBoundsException e) {
     System.err.println("Fehler beim Ablegen der Daten in Table4,");
     System.err.println("wahrscheinlich ist ein Formatierungsfehler in Parameter/Yang2000_Absorption.txt" );
    };
   };
  };


  // ---== Read Table 5 ==---
  // Index 1: das gewaehlte spektrale Band [0...6]
  // Index 2: shape 0 = Plate
  //           1 = Column
  //           2 = Hollow Column
  //           3 = Bullet rosettes-4
  //           4 = Bullet rosettes-6
  //           5 = Aggregates
  // Index 3: 0 = lambda
  //          1 = xi1_g
  //          2 = xi2_g
  //          3 = xi3_g
  //          4 = xi4_g
  //          5 = xi5_g
  //          6 = xi6_g
  //          7 = xi7_g
  file= null;
  try {
   file=new BufferedReader(new FileReader("Parameter/Yang2000_Asymmetry.txt"));
   System.out.println("Datei Parameter/Yang2000_Asymmetry.txt gefunden" );
  } catch(IOException e) {System.err.println("Datei Parameter/Yang2000_Asymmetry.txt nicht gefunden" );};
    
  //Datei ist geoeffnet und kann ueber file gelesen werden
  // Zerlegen der Eingangsdaten
  p = Pattern.compile( "\\s" ); // wird zum Teilen der Zeilen verwendet
  list= new ArrayList<String>(); // wird mit den Messwerten gefuellt
  try {
   // ignore header (2 lines)
   eineZeile= file.readLine();
   eineZeile= file.readLine();
   int i = 1;
   while (( eineZeile= file.readLine())!=null){ // Solange noch eine Zeile gelesen werden kann
    if ((i%8)!=0){
     list.add(eineZeile);
     }
    i++;
   }
   file.close(); // schliessen der Datei
  } catch (IOException ex){System.err.println("Fehler beim Einlesen (Table 5)");} // Falls ein Fehler auftrat
    
  // lege Daten in Array table5 ab
  table5= new double[7][6][8]; // Initialisieren des Arrays mit der richtigen Groesse
  for (int i=0; i<7; i++) { // das spektrale Band
   for (int j=0; j<6; j++) { // shape
    try {
     String InfoString[] = p.split( list.get(i+7*j) ); // Teile die Zeile in ein Array aus Stringvariablen
	 for (int k = 0; k<8; k++) { // Parameter
      table5[i][j][k]=Double.valueOf(InfoString[k]);
     };
    } catch(IndexOutOfBoundsException e) {
     System.err.println("Fehler beim Ablegen der Daten in Table5,");
     System.err.println("wahrscheinlich ist ein Formatierungsfehler in Parameter/Yang2000_Asymmetry.txt" );
    };
   };
  };
  
  
  // ---== Read table_terr1 ==---
  // Fuer die Berechnung im terrestrischen Wellenlaengenbereich
  // Index 1: das gewaehlte spektrale Band [0...6]
  // Index 2: 0 = lambda
  //          1 = eta1
  //          2 = eta2
  //          3 = eta3
  //          4 = xi0
  //          5 = xi1
  //          6 = xi2
  //          7 = xi3
  //          8 = zeta0
  //          9 = zeta1
  //         10 = zeta2
  //         11 = zeta3
  file= null;
  try {
   file=new BufferedReader(new FileReader("Parameter/Yang2005_terrestrSpektrum.txt"));
   System.out.println("Datei Parameter/Yang2005_terrestrSpektrum.txt gefunden" );
  } catch(IOException e) {System.err.println("Datei Parameter/Yang2005_terrestrSpektrum.txt nicht gefunden" );};
    
  //Datei ist geoeffnet und kann ueber file gelesen werden
  // Zerlegen der Eingangsdaten
  p = Pattern.compile( "\\s" ); // wird zum Teilen der Zeilen verwendet
  list= new ArrayList<String>(); // wird mit den Messwerten gefuellt
  try {
   // ignore header (1 line)
   eineZeile= file.readLine();
   while (( eineZeile= file.readLine())!=null){ // Solange noch eine Zeile gelesen werden kann
    list.add(eineZeile);
   }
   file.close(); // schliessen der Datei
  } catch (IOException ex){System.err.println("Fehler beim Einlesen (table_terr1)");} // Falls ein Fehler auftrat
    
  // lege Daten in Array table1 ab
  table_terr1= new double[49][12]; // Initialisieren des Arrays mit der richtigen Groesse
  for (int i=0; i<49; i++) {
   try {
    String InfoString[] = p.split( list.get(i) ); // Teile die Zeile in ein Array aus Stringvariablen
	for (int j = 0; j<12; j++){
     table_terr1[i][j]=Double.valueOf(InfoString[j]);
    };
   } catch(IndexOutOfBoundsException e){
    System.err.println("Fehler beim Ablegen der Daten in table_terr1,");
    System.err.println("wahrscheinlich ist ein Formatierungsfehler in Parameter/Yang2005_terrestrSpektrum.txt" );
   };
  };
  
  return;
 }
 
 private double Calc_Da(double Dmax, int shape){
  // spherical diameter with equivalent projected area
  //
  // Parameters:
  // Dmax ... maximum dimension [micrometer]
  // shape ... 0 = Plate
  //           1 = Column
  //           2 = Hollow Column
  //           3 = Bullet rosettes-4
  //           4 = Bullet rosettes-6
  //           5 = Aggregates
  double[][] a = {{0.43773, 0.33401, 0.33401, 0.15909, 0.14195, -0.47737},
                  {0.75497, 0.36477, 0.36477, 0.84308, 0.84394, 0.10026e1},
                  {0.19033e-1, 0.30855, 0.30855, 0.70161e-2, 0.72125e-2, -0.10030e-2},
                  {0.35191e-3, -0.55631e-1, -0.55631e-1, -0.11003e-2, -0.11219e-2, 0.15166e-3},
                  {-0.70782e-4, 0.30162e-2, 0.30162e-2, 0.45161e-4, 0.45819e-4, -0.78433e-5}};
  double Da = 0;
  double help = 0;
  for (int n = 0; n<5; n+=1) {
   help += a[n][shape] * Math.pow(Math.log(Dmax),n);
  };
  Da = Math.exp(help);
  return Da;
 }

 private double Calc_Dv(double Dmax, int shape){
  // spherical diameter with equivalent volume
  //
  // Parameters:
  // Dmax ... maximum dimension [micrometer]
  // shape ... 0 = Plate
  //           1 = Column
  //           2 = Hollow Column
  //           3 = Bullet rosettes-4
  //           4 = Bullet rosettes-6
  //           5 = Aggregates
  double[][] b = {{0.31228, 0.30581, 0.24568, -0.97940e-1, -0.10318, -0.70160},
                  {0.80874, 0.26252, 0.26202, 0.85683, 0.86290, 0.99215},
                  {0.29287e-2, 0.35458, 0.35479, 0.29483e-2, 0.70665e-3, 0.29322e-2},
                  {-0.44378e-3, -0.63202e-1, -0.63236e-1, -0.14341e-2, -0.11055e-2, -0.40492e-3},
                  {0.23109e-4, 0.33755e-2, 0.33773e-2, 0.74627e-4, 0.57906e-4, 0.18841e-4}};
  double Dv = 0;
  double help = 0;
  for (int n = 0; n<5; n+=1) {
   help += b[n][shape] * Math.pow(Math.log(Dmax),n);
  };
  Dv = Math.exp(help);
  return Dv;
 }

 private double Calc_de(double Dmax, int shape){
  // effective diameter
  //
  // Parameters:
  // Dmax ... maximum dimension [micrometer]
  // shape ... 0 = Plate
  //           1 = Column
  //           2 = Hollow Column
  //           3 = Bullet rosettes-4
  //           4 = Bullet rosettes-6
  //           5 = Aggregates
  double de = Math.pow(Calc_Dv(Dmax, shape),3)/Math.pow(Calc_Da(Dmax, shape),2);
  return de;
 }

 private double Calc_de_mix(double Dmax){
  // effective diameter
  //
  // Parameters:
  // Dmax ... maximum dimension [micrometer]
  double de = 0;
  if (Dmax<70){
    // 50% Bullet rosettes-6
    // 25% hollow columns
    // 25% plates
    de = (0.5*Math.pow(Calc_Dv(Dmax, 4),3)+0.25*Math.pow(Calc_Dv(Dmax, 2),3)+0.25*Math.pow(Calc_Dv(Dmax, 0),3))/
         (0.5*Math.pow(Calc_Da(Dmax, 4),2)+0.25*Math.pow(Calc_Da(Dmax, 2),2)+0.25*Math.pow(Calc_Da(Dmax, 0),2));
  } else {
    // 30% Aggregates
    // 30% Bullet rosettes-6
    // 20% hollow columns
    // 20% plates
    de = (0.3*Math.pow(Calc_Dv(Dmax, 4),3)+0.3*Math.pow(Calc_Dv(Dmax, 4),3)+0.2*Math.pow(Calc_Dv(Dmax, 2),3)+0.2*Math.pow(Calc_Dv(Dmax, 0),3))/
         (0.3*Math.pow(Calc_Da(Dmax, 4),2)+0.3*Math.pow(Calc_Da(Dmax, 4),2)+0.2*Math.pow(Calc_Da(Dmax, 2),2)+0.2*Math.pow(Calc_Da(Dmax, 0),2));
  };
  return de;
 }
 
 private double Calc_Qextprime(double Dmax, int band, int shape){
  // Extincion efficiency of a sphere
  // Is needed for calculating Qext
  //
  // Parameters:
  // Dmax ... maximum dimension [micrometer]
  // band ... das gewuenschte spektrale Band
  //   0 = 0.55 micrometer
  //   1 = 1.35 micrometer
  //   2 = 2.25 micrometer
  //   3 = 3.0125 micrometer
  //   4 = 3.775 micrometer
  //   5 = 4.5 micrometer
  // shape ... 0 = Plate
  //           1 = Column
  //           2 = Hollow Column
  //           3 = Bullet rosettes-4
  //           4 = Bullet rosettes-6
  //           5 = Aggregates
  double lambda = bands[band]; // wavelength [micrometer]
  double mr = table1[band][3];
  double rhoe = 2*Math.PI*Calc_de(Dmax, shape)*Math.abs(mr-1)/lambda; //effective phase delay
  double eta1 = table3[band][shape][1];  
  double rho = eta1 * rhoe;
  double beta = table3[band][shape][3];
  double Qextprime = 2. - 4.*Math.exp(-rho*Math.tan(beta)) *
                    (Math.cos(beta)/rho*Math.sin(rho-beta) + Math.pow(Math.cos(beta)/rho,2) * Math.cos(rho-2*beta)) +
                    4.*Math.pow(Math.cos(beta)/rho,2)*Math.cos(2*beta);
  return Qextprime;
 }
 
 private double Calc_Qextprimeprime(double Dmax, int band, int shape){
  // Extincion efficiency of nonspherical ice crystal without complex ray behavior
  // Is needed for calculating Qext
  //
  // Parameters:
  // Dmax ... maximum dimension [micrometer]
  // band ... das gewuenschte spektrale Band
  //   0 = 0.55 micrometer
  //   1 = 1.35 micrometer
  //   2 = 2.25 micrometer
  //   3 = 3.0125 micrometer
  //   4 = 3.775 micrometer
  //   5 = 4.5 micrometer
  // shape ... 0 = Plate
  //           1 = Column
  //           2 = Hollow Column
  //           3 = Bullet rosettes-4
  //           4 = Bullet rosettes-6
  //           5 = Aggregates
  double lambda = bands[band]; // wavelength [micrometer]
  double mr = table1[band][3];
  double rhoe = 2*Math.PI*Calc_de(Dmax, shape)*Math.abs(mr-1)/lambda; //effective phase delay
  double eta2 = table3[band][shape][2];
  double beta2 = table3[band][shape][4];
  double alpha = table3[band][shape][5];
  double Qextprimeprime = 2*(1-Math.exp(-2./3.*rhoe*eta2*Math.tan(beta2))*Math.cos(2./3.*rhoe*eta2+alpha));
  return Qextprimeprime;
 }
 
 public double Calc_Qext_sol_i(double Dmax, int band, int shape){
  // Extincion efficiency for nonspherical ice crystal under consideration of complex ray behavior
  //
  // Parameters:
  // Dmax ... maximum dimension [micrometer]
  // band ... das gewuenschte spektrale Band
  //   0 = 0.55 micrometer
  //   1 = 1.35 micrometer
  //   2 = 2.25 micrometer
  //   3 = 3.0125 micrometer
  //   4 = 3.775 micrometer
  //   5 = 4.5 micrometer
  // shape ... 0 = Plate
  //           1 = Column
  //           2 = Hollow Column
  //           3 = Bullet rosettes-4
  //           4 = Bullet rosettes-6
  //           5 = Aggregates
  double Qextprime = Calc_Qextprime(Dmax, band, shape);
  double Qextprimeprime = Calc_Qextprimeprime(Dmax, band, shape);
  double xi_ext = table3[band][shape][6];
  
  double Qext_sol =  (1.-xi_ext) * Qextprime + xi_ext * Qextprimeprime;
  return Qext_sol;
 }
 
 public double Calc_Qext_sol(double Dmax, int band){
  // Extincion efficiency for nonspherical ice crystal under consideration of complex ray behavior
  //
  // Parameters:
  // Dmax ... maximum dimension [micrometer]
  // band ... das gewuenschte spektrale Band
  //   0 = 0.55 micrometer
  //   1 = 1.35 micrometer
  //   2 = 2.25 micrometer
  //   3 = 3.0125 micrometer
  //   4 = 3.775 micrometer
  //   5 = 4.5 micrometer
  // shape ... 0 = Plate
  //           1 = Column
  //           2 = Hollow Column
  //           3 = Bullet rosettes-4
  //           4 = Bullet rosettes-6
  //           5 = Aggregates
  double[] f = {0., 0., 0., 0., 0., 1.};
  if (Dmax<70){
    // 50% Bullet rosettes-6
    // 25% hollow columns
    // 25% plates
    f = new double[] {0.25, 0., 0.25, 0., 0.5, 0.};
  } else {
    // 30% Aggregates
    // 30% Bullet rosettes-6
    // 20% hollow columns
    // 20% plates
    f = new double[] {0.2, 0., 0.2, 0., 0.3, 0.3};
  };
  double zaehler =0.;
  double nenner =0.;
  for (int i = 0; i<6; i++) {
   zaehler += f[i]*Math.pow(Calc_Da(Dmax,i),2)*Calc_Qext_sol_i(Dmax, band, i);
   nenner += f[i]*Math.pow(Calc_Da(Dmax,i),2);
  };
    
  double Qext_sol =  zaehler/nenner;
  return Qext_sol;
 }
 
 private double Calc_Qabsprime(double Dmax, int band, int shape){
  // Absorpion efficiency of a sphere
  // Is needed for calculating Qabs
  //
  // Parameters:
  // Dmax ... maximum dimension [micrometer]
  // band ... das gewuenschte spektrale Band
  //   0 = 0.55 micrometer
  //   1 = 1.35 micrometer
  //   2 = 2.25 micrometer
  //   3 = 3.0125 micrometer
  //   4 = 3.775 micrometer
  //   5 = 4.5 micrometer
  // shape ... 0 = Plate
  //           1 = Column
  //           2 = Hollow Column
  //           3 = Bullet rosettes-4
  //           4 = Bullet rosettes-6
  //           5 = Aggregates
  double lambda = bands[band]; // wavelength [micrometer]
  double mi = table1[band][4];
  double chie = 2*Math.PI*Calc_de(Dmax, shape)/2/lambda; //effective size parameter of nonspherical ice crystal
  double gammae = 4*chie*mi;
  double eta1_abs = table4[band][shape][1];
  double Qabsprime = 1+2.*Math.exp(-1.*gammae*eta1_abs)/(gammae*eta1_abs)+2.*(Math.exp(-1*gammae*eta1_abs)-1)/Math.pow(gammae*eta1_abs,2);
  return Qabsprime;
 }
 
 private double Calc_Qabsprimeprime(double Dmax, int band, int shape){
  // Absorpion efficiency of nonspherical ice crystal without complex ray behavior
  // Is needed for calculating Qabs
  //
  // Parameters:
  // Dmax ... maximum dimension [micrometer]
  // band ... das gewuenschte spektrale Band
  //   0 = 0.55 micrometer
  //   1 = 1.35 micrometer
  //   2 = 2.25 micrometer
  //   3 = 3.0125 micrometer
  //   4 = 3.775 micrometer
  //   5 = 4.5 micrometer
  // shape ... 0 = Plate
  //           1 = Column
  //           2 = Hollow Column
  //           3 = Bullet rosettes-4
  //           4 = Bullet rosettes-6
  //           5 = Aggregates
  double lambda = bands[band]; // wavelength [micrometer]
  double mi = table1[band][4];
  double chie = 2.*Math.PI*Calc_de(Dmax, shape)/2./lambda; //effective size parameter of nonspherical ice crystal
  double gammae = 4.*chie*mi;
  double eta2_abs = table4[band][shape][2];
  
  double Qabsprimeprime = 1-Math.exp(-2./3.*gammae*eta2_abs);
  return Qabsprimeprime;
 }
 
 public double Calc_Qabs_sol_i(double Dmax, int band, int shape){
  // Absorption efficiency for nonspherical ice crystal under consideration of complex ray behavior
  //
  // Parameters:
  // Dmax ... maximum dimension [micrometer]
  // band ... das gewuenschte spektrale Band
  //   0 = 0.55 micrometer
  //   1 = 1.35 micrometer
  //   2 = 2.25 micrometer
  //   3 = 3.0125 micrometer
  //   4 = 3.775 micrometer
  //   5 = 4.5 micrometer
  // shape ... 0 = Plate
  //           1 = Column
  //           2 = Hollow Column
  //           3 = Bullet rosettes-4
  //           4 = Bullet rosettes-6
  //           5 = Aggregates
  double Qabsprime = Calc_Qabsprime(Dmax, band, shape);
  double Qabsprimeprime = Calc_Qabsprimeprime(Dmax, band, shape);
  double xi1_abs = table4[band][shape][3];
  double xi2_abs = table4[band][shape][4];
  
  double Qabs_sol = (1-xi1_abs)*((1-xi2_abs)*Qabsprime + xi2_abs*Qabsprimeprime);
  return Qabs_sol;
 }
 
 public double Calc_Qabs_sol(double Dmax, int band){
  // Absorption efficiency for nonspherical ice crystal under consideration of complex ray behavior
  //
  // Parameters:
  // Dmax ... maximum dimension [micrometer]
  // band ... das gewuenschte spektrale Band
  //   0 = 0.55 micrometer
  //   1 = 1.35 micrometer
  //   2 = 2.25 micrometer
  //   3 = 3.0125 micrometer
  //   4 = 3.775 micrometer
  //   5 = 4.5 micrometer
  // shape ... 0 = Plate
  //           1 = Column
  //           2 = Hollow Column
  //           3 = Bullet rosettes-4
  //           4 = Bullet rosettes-6
  //           5 = Aggregates
  double[] f = {0., 0., 0., 0., 0., 1.};
  if (Dmax<70){
    // 50% Bullet rosettes-6
    // 25% hollow columns
    // 25% plates
    f = new double[] {0.25, 0., 0.25, 0., 0.5, 0.};
  } else {
    // 30% Aggregates
    // 30% Bullet rosettes-6
    // 20% hollow columns
    // 20% plates
    f = new double[] {0.2, 0., 0.2, 0., 0.3, 0.3};
  };
  double zaehler =0.;
  double nenner =0.;
  for (int i = 0; i<6; i++) {
   //System.out.println(Dmax+" "+i+" "+Calc_Da(Dmax,i)+" "+Calc_Qabs_sol_i(Dmax, band, i));
   zaehler += f[i]*Math.pow(Calc_Da(Dmax,i),2)*Calc_Qabs_sol_i(Dmax, band, i);
   nenner += f[i]*Math.pow(Calc_Da(Dmax,i),2);
  };
    
  double Qabs_sol =  zaehler/nenner;
  return Qabs_sol;
 }
 
 public double Calc_Qsca_sol_i(double Dmax, int band, int shape){
  // Scattering efficiency for nonspherical ice crystal under consideration of complex ray behavior
  //
  // Parameters:
  // Dmax ... maximum dimension [micrometer]
  // band ... das gewuenschte spektrale Band
  //   0 = 0.55 micrometer
  //   1 = 1.35 micrometer
  //   2 = 2.25 micrometer
  //   3 = 3.0125 micrometer
  //   4 = 3.775 micrometer
  //   5 = 4.5 micrometer
  // shape ... 0 = Plate
  //           1 = Column
  //           2 = Hollow Column
  //           3 = Bullet rosettes-4
  //           4 = Bullet rosettes-6
  //           5 = Aggregates
  double Qsca_sol = Calc_Qext_sol_i(Dmax, band, shape) - Calc_Qabs_sol_i(Dmax, band, shape);
  return Qsca_sol;
 }
 
 public double Calc_Qsca_sol(double Dmax, int band){
  // Scattering efficiency for nonspherical ice crystal under consideration of complex ray behavior
  //
  // Parameters:
  // Dmax ... maximum dimension [micrometer]
  // band ... das gewuenschte spektrale Band
  //   0 = 0.55 micrometer
  //   1 = 1.35 micrometer
  //   2 = 2.25 micrometer
  //   3 = 3.0125 micrometer
  //   4 = 3.775 micrometer
  //   5 = 4.5 micrometer
  double Qsca_sol = Calc_Qext_sol(Dmax, band) - Calc_Qabs_sol(Dmax, band);
  return Qsca_sol;
 }
 
 public double Calc_g_sol_i(double Dmax, int band, int shape){
  // solar asymmetrie factor of nonspherical ice crystal under consideration of complex ray behavior
  //
  // Parameters:
  // Dmax ... maximum dimension [micrometer]
  // band ... das gewuenschte spektrale Band
  //   0 = 0.55 micrometer
  //   1 = 1.35 micrometer
  //   2 = 2.25 micrometer
  //   3 = 3.0125 micrometer
  //   4 = 3.775 micrometer
  //   5 = 4.5 micrometer
  // shape ... 0 = Plate
  //           1 = Column
  //           2 = Hollow Column
  //           3 = Bullet rosettes-4
  //           4 = Bullet rosettes-6
  //           5 = Aggregates
  double lambda = bands[band]; // wavelength [micrometer]
  double chie = 2.*Math.PI*Calc_de(Dmax, shape)/2./lambda; //effective size parameter of nonspherical ice crystal
  double xi1_g = table5[band][shape][1];
  double xi2_g = table5[band][shape][2];
  double xi3_g = table5[band][shape][3];
  double xi4_g = table5[band][shape][4];
  double xi5_g = table5[band][shape][5];
  double xi6_g = table5[band][shape][6];
  double xi7_g = table5[band][shape][7];
  
  double f1 = (1-xi1_g) * ( 1 - (1-Math.exp(-1.*xi2_g*(chie + xi3_g))) / (xi2_g*(chie+xi3_g)) );
  double f2 = (1-xi4_g) * (1 - Math.exp(-1.*xi5_g*(chie+xi6_g)));
  double g_sol = (1-xi7_g)*f1 + xi7_g*f2; //solarer Asymmetriefaktor 
  return g_sol;
 }
 
 public double Calc_g_sol(double Dmax, int band){
  // solar asymmetrie factor of nonspherical ice crystal under consideration of complex ray behavior
  //
  // Parameters:
  // Dmax ... maximum dimension [micrometer]
  // band ... das gewuenschte spektrale Band
  //   0 = 0.55 micrometer
  //   1 = 1.35 micrometer
  //   2 = 2.25 micrometer
  //   3 = 3.0125 micrometer
  //   4 = 3.775 micrometer
  //   5 = 4.5 micrometer
  // shape ... 0 = Plate
  //           1 = Column
  //           2 = Hollow Column
  //           3 = Bullet rosettes-4
  //           4 = Bullet rosettes-6
  //           5 = Aggregates
  double[] f = {0., 0., 0., 0., 0., 1.};
  if (Dmax<70){
    // 50% Bullet rosettes-6
    // 25% hollow columns
    // 25% plates
    f = new double[] {0.25, 0., 0.25, 0., 0.5, 0.};
  } else {
    // 30% Aggregates
    // 30% Bullet rosettes-6
    // 20% hollow columns
    // 20% plates
    f = new double[] {0.2, 0., 0.2, 0., 0.3, 0.3};
  };
  double zaehler =0.;
  double nenner =0.;
  for (int i = 0; i<6; i++) {
   zaehler += f[i]*Math.pow(Calc_Da(Dmax,i),2)*Calc_Qsca_sol_i(Dmax, band, i)*Calc_g_sol_i(Dmax, band, i);
   nenner += f[i]*Math.pow(Calc_Da(Dmax,i),2)*Calc_Qsca_sol_i(Dmax, band, i);
  };
    
  double g_sol =  zaehler/nenner;
  return g_sol;
 }
 
 public double Calc_Qext_terr(double Dmax, double lambda){
  // Extinction efficiency for nonspherical ice crystal under consideration of complex ray behavior
  //
  // Parameters:
  // Dmax ... maximum dimension [micrometer]
  // lambda ... wavelength [micrometer] 3.08 < lambda < 99.99
  //
  if (lambda < 3.08) System.out.println("Error Calculating Qext_terr: Lambda "+lambda+" < 3.08\n Result should not be trusted.");
  if (lambda > 99.9) System.out.println("Error Calculating Qext_terr: Lambda "+lambda+" > 99.9\n Result should not be trusted.");
  
  double de = Calc_de_mix(Dmax);
  // Find parameters in table_terr1
  int band = 0;
  while ((lambda > (table_terr1[band][0]+table_terr1[band+1][0])/2.) && (band < 47)) {
   band++;
  };
  if (lambda > (table_terr1[47][0]+table_terr1[48][0])/2.) band = 48;
  double eta1 = table_terr1[band][1];
  double eta2 = table_terr1[band][2];
  double eta3 = table_terr1[band][3];
  
  double Qext_terr = (2.+eta1/de)/(1.+eta2/de+eta3/Math.pow(de,2.));
  return Qext_terr;
}
 
 public double Calc_Qabs_terr(double Dmax, double lambda){
  // Absorption efficiency for nonspherical ice crystal under consideration of complex ray behavior
  //
  // Parameters:
  // Dmax ... maximum dimension [micrometer]
  // lambda ... wavelength [micrometer] 3.08 < lambda < 99.99
  //
  if (lambda < 3.08) System.out.println("Error Calculating Qabs_terr: Lambda "+lambda+" < 3.08\n Result should not be trusted.");
  if (lambda > 99.9) System.out.println("Error Calculating Qabs_terr: Lambda "+lambda+" > 99.9\n Result should not be trusted.");
  
  double de = Calc_de_mix(Dmax);
  // Find parameters in table_terr1
  int band = 0;
  while ((lambda > (table_terr1[band][0]+table_terr1[band+1][0])/2.) && (band < 47)) {
   band++;
  };
  if (lambda > (table_terr1[47][0]+table_terr1[48][0])/2.) band = 48;
  
  double xi0 = table_terr1[band][4];
  double xi1 = table_terr1[band][5];
  double xi2 = table_terr1[band][6];
  double xi3 = table_terr1[band][7];
  
  double Qabs_terr = (xi0+xi1/de)/(1.+xi2/de+xi3/Math.pow(de,2.));
  return Qabs_terr;
}
 
 public double Calc_Qsca_terr(double Dmax, double lambda){
  // Scattering efficiency for nonspherical ice crystal under consideration of complex ray behavior
  //
  // Parameters:
  // Dmax ... maximum dimension [micrometer]
  // lambda ... wavelength [micrometer] 3.08 < lambda < 99.99
  //
  if (lambda < 3.08) System.out.println("Error Calculating Qsca_terr: Lambda "+lambda+" < 3.08\n Result should not be trusted.");
  if (lambda > 99.9) System.out.println("Error Calculating Qsca_terr: Lambda "+lambda+" > 99.9\n Result should not be trusted.");
  
  double Qsca_terr = Calc_Qext_terr(Dmax, lambda)-Calc_Qabs_terr(Dmax, lambda);
  return Qsca_terr;
}
 
 public double Calc_g_terr(double de, double lambda){
  // Assymetry factor for nonspherical ice crystal under consideration of complex ray behavior
  //
  // Parameters:
  // de... effective diameter [micrometer]
  // lambda ... wavelength [micrometer] 3.08 < lambda < 99.99
  //
  if (lambda < 3.08) System.out.println("Error Calculating g_terr: Lambda "+lambda+" < 3.08\n Result should not be trusted.");
  if (lambda > 99.9) System.out.println("Error Calculating g_terr: Lambda "+lambda+" > 99.9\n Result should not be trusted.");
  
  // Find parameters in table_terr1
  int band = 0;
  while ((lambda > (table_terr1[band][0]+table_terr1[band+1][0])/2.) && (band < 47)) {
   band++;
  };
  if (lambda > (table_terr1[47][0]+table_terr1[48][0])/2.) band = 48;
  
  double zeta0 = table_terr1[band][8];
  double zeta1 = table_terr1[band][9];
  double zeta2 = table_terr1[band][10];
  double zeta3 = table_terr1[band][11];
  
  double g_terr = (zeta0+zeta1/de)/(1.+zeta2/de+zeta3/Math.pow(de,2.));
  return g_terr;
}

 public static void main (String[] args) {
   Extinktionseffizienz a = new Extinktionseffizienz();
   
   System.out.println("Q_ext_sol "+a.Calc_Qext_sol_i(1500., 4, 0));
   System.out.println("Q_ext_sol "+a.Calc_Qext_sol_i(1500., 4, 1));
   System.out.println("Q_ext_sol "+a.Calc_Qext_sol_i(1500., 4, 2));
   System.out.println("Q_ext_sol "+a.Calc_Qext_sol_i(1500., 4, 3));
   System.out.println("Q_ext_sol "+a.Calc_Qext_sol_i(1500., 4, 4));
   System.out.println("Q_ext_sol "+a.Calc_Qext_sol_i(1500., 4, 5));
   System.out.println("Q_ext_sol_mix "+a.Calc_Qext_sol(1500., 4)+"\n");
   System.out.println("Q_abs_sol "+a.Calc_Qabs_sol_i(11.03, 0, 0));
   System.out.println("Q_abs_sol "+a.Calc_Qabs_sol_i(11.03, 0, 1));
   System.out.println("Q_abs_sol "+a.Calc_Qabs_sol_i(11.03, 0, 2));
   System.out.println("Q_abs_sol "+a.Calc_Qabs_sol_i(11.03, 0, 3));
   System.out.println("Q_abs_sol "+a.Calc_Qabs_sol_i(11.03, 0, 4));
   System.out.println("Q_abs_sol "+a.Calc_Qabs_sol_i(11.03, 0, 5));
   System.out.println("Q_abs_sol_mix "+a.Calc_Qabs_sol(11.03, 0)+"\n");
   System.out.println("g_sol "+a.Calc_g_sol_i(11., 4, 0));
   System.out.println("g_sol "+a.Calc_g_sol_i(11., 4, 1));
   System.out.println("g_sol "+a.Calc_g_sol_i(11., 4, 2));
   System.out.println("g_sol "+a.Calc_g_sol_i(11., 4, 3));
   System.out.println("g_sol "+a.Calc_g_sol_i(11., 4, 4));
   System.out.println("g_sol "+a.Calc_g_sol_i(11., 4, 5));
   System.out.println("g_sol_mix "+a.Calc_g_sol(11., 4)+"\n");
   
   System.out.println("Q_ext_terr "+a.Calc_Qext_terr(10.,8));
   System.out.println("Q_ext_terr "+a.Calc_Qext_terr(40.,8));
   System.out.println("Q_ext_terr "+a.Calc_Qext_terr(100.,8)+"\n");
   System.out.println("Q_abs_terr "+a.Calc_Qabs_terr(10.,8));
   System.out.println("Q_abs_terr "+a.Calc_Qabs_terr(40.,8));
   System.out.println("Q_abs_terr "+a.Calc_Qabs_terr(100.,8)+"\n");
   System.out.println("g_terr "+a.Calc_g_terr(10.,8));
   System.out.println("g_terr "+a.Calc_g_terr(40.,8));
   System.out.println("g_terr "+a.Calc_g_terr(100.,8)+"\n");
 }
}