package ContrailsWakeVortex2020TEST;

import java.io.FileReader;
import java.io.IOException;
import java.lang.IndexOutOfBoundsException;
import java.io.BufferedReader;
import java.util.regex.Pattern;
import java.lang.Float;
import java.util.*;

public class radiosonde
{
    //-------------------------------------------=======
    // Height	Temp  Pres  RHwManipuli  Theta + dTHETA/dz 
    // km      degC  hPa   %             K         K/m    
    //-------------------------------------------=======
  double[][] feuchtewerte;
  
  // Erzeuger zum Laden der Feuchtewerte des Radiosondenaufstiegs
  public radiosonde(String FileName)
  {
    // Filename steht fuer den Radiosondenaufstieg
    // Dateistruktur
    // Height	Temp	Pres	RHwManipuli	Theta
    // km degC	hPa	%	K
    // Daten sind in Spalten fester breite durch Tabulatoren getrennt.
    // Filehandler (Oeffnen der Radiosondendatei)
    BufferedReader file= null;
    try 
    {
      file=new BufferedReader(new FileReader(FileName));
      System.out.println("Datei "+FileName+" gefunden" );	
    }
    catch(IOException e)
    {
      System.err.println("Datei "+FileName+" nicht gefunden" );
    }
    
    //Datei ist geoeffnet und kann ueber file gelesen werden
    // Zerlegen der Eingangsdaten
    String eineZeile; // zum Zeilenweise einlesen der Datei
    Pattern p = Pattern.compile( "\t" ); // wird zum Teilen der Zeilen verwendet
    ArrayList<String> list= new ArrayList<String>(); // wird mit den Messwerten gefuellt
    try
    {
      // ignore header (2 lines)
      eineZeile= file.readLine();
      eineZeile= file.readLine();
      while (( eineZeile= file.readLine())!=null) // Solange noch eine Zeile gelesen werden kann
      {
       list.add(eineZeile);
      }
      file.close(); // schliessen der Datei
    }
    catch (IOException ex){System.err.println("Fehler beim Einlesen");} // Falls ein Fehler auftrat
    
    //Schreiben der eingelesenen Werte in das Array feuchtewerte
    this.feuchtewerte= new double[list.size()][12]; // Initialisieren des Arrays mit der richtigen Groesse
    int len = list.size(); // Anzahl der Messpunkte im Radiosondenaufstieg
    for (int i=0;i<len; i++)
    { 
      try 
      {
	String InfoString[] = p.split( list.get(i) ); // Teile die Zeile in ein Array aus Stringvariablen
	for (int j = 0; j<InfoString.length; j++){
	   //System.out.println("c "+InfoString[j]+ "; d " + InfoString.length);
           feuchtewerte[i][j]=Double.valueOf(InfoString[j]);
        };
      }
      catch(IndexOutOfBoundsException e)
      {
        System.err.println("Dieser Fehler kann eigentlich nicht auftreten. X)" );
      }
    };
    
    // Berechnen von dTheta/dz [K/m]
    feuchtewerte[0][5]= (feuchtewerte[1][4]-feuchtewerte[0][4])/(feuchtewerte[1][0]-feuchtewerte[0][0])/1000.;
    feuchtewerte[len-1][5]= (feuchtewerte[len-1][4]-feuchtewerte[len-2][4])/(feuchtewerte[len-1][0]-feuchtewerte[len-2][0])/1000.;
    for (int i=1; i< len-1; i++){
      feuchtewerte[i][5]= (feuchtewerte[i+1][4]-feuchtewerte[i-1][4])/(feuchtewerte[i+1][0]-feuchtewerte[i-1][0])/1000.;
    };
    
    return; // Alles ist glatt gelaufen der Erzeuger kann seine Arbeit einstellen.
  }
  
  private double interpol(double altitude, int colum)
  {
   // Gibt den gewuenschten Messwert in einer spezifizierten Hoehe zurueck.
   // Notfalls wird linear interpoliert
   altitude=altitude/1000.;
   int n=0;
   while (this.feuchtewerte[n][0]<altitude)
   {
     n=n+1;
   }
   double f=0;
   if(n==0)
   {
     f=this.feuchtewerte[0][colum];
   } else {
     double h2=this.feuchtewerte[n][0];
     double h1=this.feuchtewerte[n-1][0];
     double f2=this.feuchtewerte[n][colum];
     double f1=this.feuchtewerte[n-1][colum];
   
     double a=(f2-f1)/(h2-h1);
     f=f1+(altitude-h1)*a;
   }
   return f;
  }
  
  public double get_p(double alt){
     //pressure in hPa
     return interpol(alt, 2);
  }
  
  public double get_T(double alt){
     //temperature in degree Celsius
     return interpol(alt, 1);
  }
  
  public double get_TK(double alt){
     //temperature in Kelvin
     return interpol(alt, 1)+273.15;
  }
  
  public double get_RH(double alt){
     return interpol(alt, 3);
  }
  
  public double get_Theta(double alt){
     return interpol(alt, 4);
  }
  
  public double get_dThetadz(double alt){
     return interpol(alt, 5);
  }
}