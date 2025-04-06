package ContrailsWakeVortex2020;

import java.lang.Math;
import java.lang.System;

public class start_wake_vortex
{
 
 public static void main(String[] args)
 {
  WeatherFromGrib TEST = new WeatherFromGrib("/Users/marco/Desktop/testordner/gribtesttest");
  System.out.println(TEST.getWindHeadingFromPressure(0,0,20000));
 }
} 
